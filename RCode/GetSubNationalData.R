################
### ODDpolys ###
################



library(openxlsx)
library(rgeos)

cleanSubNatData <- function(SubNatData){
  # Clean data from xlsx file by converting data to the appropriate formats/data types
  SubNatData$source_date <- openxlsx::convertToDate(SubNatData$source_date)
  SubNatData$sdate <- openxlsx::convertToDate(SubNatData$sdate)
  SubNatData$fdate <- openxlsx::convertToDate(SubNatData$fdate)
  SubNatData$mortality <- as.integer(SubNatData$mortality)
  SubNatData$displacement <- as.integer(SubNatData$displacement)
  SubNatData$buildDam <- as.integer(SubNatData$buildDam)
  SubNatData$buildDest <- as.integer(SubNatData$buildDest)
  return(SubNatData)
}

getPolyData <- function(polygon_name, subregion, region, country, iso3){
  # Uses getGADM() or getbb() to retrieve the polygon of a region
  # Details:
  #    - Polygon_name is just used to label the polygon at the end (and isn't used to actually find the region)
  #    - Any missing values for subregion or region should be set to NA
  #    - Returns a list containing:
  #         - $polygon_name set to polygon_name
  #         - $sf_polygon set to a polygon (that works with the sf package) corresponding to this region, or NULL if polygon not found
  #    - Attempts to use ecochange::getGADM() function but if this fails will attempt the osmdata::getbb() function

  GADM_level <- 2
  GADM_array <- c(subregion, region, country)
  GADM_iso3 <- iso3
  if (is.na(GADM_array[1]) & is.na(GADM_array[2])){
    GADM_level <- 0
    GADM_array <- GADM_array[3]
  } else if (is.na(GADM_array[1])){
    GADM_level <- 1
    GADM_array <- GADM_array[2]
  } else {
    GADM_level <- 2
    if (is.na(GADM_array[2])){
      GADM_array <- GADM_array[1]
    } else {
      GADM_array <- GADM_array[c(1,2)]
    }
  }
  
  # getGADM() gives off warning 'the condition has length > 1 and only the first element will be used' if GADM_array has 
  # length greater than one (which is sometimes necessary!)
  polygon <- tryCatch(getGADM(unit.nm=GADM_array, level=GADM_level, country=GADM_iso3),error=function(e) NULL) 
  if (is.null(polygon) || NROW(polygon)==0){
    regions_found <- c(getGADM(level=1, country=GADM_iso3), getGADM(level=2, country=GADM_iso3))
    warning(paste0('No polygons (or more than one polygon) found for ', polygon_name, 
                   ' using getGADM(). Closest region found: ', regions_found[which.min(adist(GADM_array[1], regions_found ))]))
    
    #attempt getbb()
    if(GADM_level != 0) GADM_array <- c(GADM_array, country)
    polygon <- getbb(paste(GADM_array, collapse=','), format_out='sf_polygon')
    if (!is.null(polygon)){
      if(!is.null(polygon$multipolygon)) polygon <- polygon@multipolygon
      polygon <- polygon %>% as("Spatial")
      print('Polygon found using getbb() but not getGADM()') 
    }
    if(NROW(polygon)>1){
      warning(paste('Multiple Polygons for', polygon_name))
      return(list(polygon_name = polygon_name, sf_polygon = polygon[1,]))
    } 
  }
  return(list(polygon_name = polygon_name, sf_polygon = polygon))
}

getSubNatImpact <- function(SubNatEvent, subnational=TRUE){
  # Uses data from SubNatEvent to populate the 'impact' slot in an ODD object
  #  - Impact has one row per impact and region (and source)
  #  - Also creates polygons_list such that polygons_list[[i]] is the sf_polygon corresponding to the polygon with id i. 

  impact <- data.frame(iso3 = character(), polygon = numeric(), impact = character(), 
                       observed = numeric(), qualifier = character())
  
  polygons_list <- list()
  nans_polygon_id <- c()
  
  if (!subnational){ #if not interested in subnational data then just take the nationally aggregated data:
    SubNatEvent <- SubNatEvent %>% filter(is.na(Region) & is.na(Subregion))
    i = 1
    for (iso3 in unique(SubNatEvent$iso3)){
      SubNatEvent$polygon_id <- i
      polygons_list[[i]] <- list(polygon_name=countrycode(iso3, origin='iso3c', destination='country.name'), sf_polygon=NULL) #sf_polygon not required for countries
      i = i + 1
    }
  } else { #work with the subnational data
    SubNatEvent$polygon_name <- paste0(ifelse(is.na(SubNatEvent$Subregion),'',paste0(SubNatEvent$Subregion,', ')), 
                                             ifelse(is.na(SubNatEvent$Region),'',paste0(SubNatEvent$Region,', ')), 
                                             ifelse(SubNatEvent$country=='TOTAL','',SubNatEvent$country))
    
    # group by polygon name and assign each an id
    SubNatEvent %<>% group_by(polygon_name) %>% mutate(polygon_id=cur_group_id()) %>% ungroup()
    
    for (i in sort(unique(SubNatEvent$polygon_id))){
      #collect polygon information for each polygon
      polygon_name <- SubNatEvent[which(SubNatEvent$polygon_id == i)[1],]$polygon_name
      region_data <- as.character(SubNatEvent[which(SubNatEvent$polygon_id == i)[1], c('Subregion', 'Region', 'country', 'iso3')])
      polygons_list[[i]] <- getPolyData(polygon_name, region_data[1], region_data[2], region_data[3], region_data[4])

      if (is.null(polygons_list[[i]]$sf_polygon) & (!is.na(region_data[1]) | !is.na(region_data[2]))) nans_polygon_id %<>% append(i) #remove empty non-national polygons
      # should probably plot polygons here to make sure they all seem reasonable!
    }
  }
  
  SubNatEvent %<>% filter(!(polygon_id %in% nans_polygon_id))
  
  for (impact_type in Model$impacts$labels){
    
    notnans <- which(!is.na(SubNatEvent[,impact_type]))
    
    if (length(notnans)>0){
      #assign each source an id:
      SubNatEvent[notnans,'source_id'] <- notnans
      sources_selected <- c()
      SubNatEvent_by_polygon <- SubNatEvent[notnans,] %>% group_by(polygon_id) %>% group_split()
      for (i in 1:length(SubNatEvent_by_polygon)){
        if(NROW(SubNatEvent_by_polygon[[i]])==1){ #only one source for that polygon
          sources_selected %<>% append(SubNatEvent_by_polygon[[i]]$source_id) 
          next
        } 
        if(length(unique(pull(SubNatEvent_by_polygon[[i]][,impact_type])))==1){ #all matching, so just take first
          sources_selected %<>% append(SubNatEvent_by_polygon[[i]]$source_id[1])
          next
        }
        #conflicting sources:
        cat(paste('Please select between the following sources, by typing the id of the chosen source and pressing return.\n',
                  'If you would like to select more than one source, please separate each id with a comma and we will use the mean \n')) 
        if(impact_type == 'buildDam' || impact_type == 'buildDest'){
          print(cbind(id=1:NROW(SubNatEvent_by_polygon[[i]]), SubNatEvent_by_polygon[[i]][,c('event_name', 'source_type', 'source_date', 'building_type', impact_type, 'notes', paste0(impact_type,'_term'))]))
        } else {
          print(cbind(id=1:NROW(SubNatEvent_by_polygon[[i]]), SubNatEvent_by_polygon[[i]][,c('event_name', 'source_type', 'source_date', impact_type, 'notes', paste0(impact_type,'_term'))]))
        }
        ids_chosen <- as.integer(unlist(strsplit(readline(prompt='id selected: '), ",")))
        
        if (length(ids_chosen)>1){ #if multiple ids selected then choose the first and edit its value to be the mean of all
          observed_mean <- round(mean(pull(SubNatEvent_by_polygon[[i]][ids_chosen, impact_type])))
          SubNatEvent[SubNatEvent_by_polygon[[i]]$source_id[ids_chosen[1]], impact_type] <- observed_mean
          ids_chosen <- ids_chosen[1]
        } 
        sources_selected %<>% append(SubNatEvent_by_polygon[[i]]$source_id[ids_chosen])
      }
      impact %<>% rbind(SubNatEvent[sources_selected,] %>% transmute(iso3=iso3, sdate=sdate,
                                                                         impact=impact_type,
                                                                         observed=!!sym(impact_type), 
                                                                         qualifier=!!sym(paste0(impact_type, '_qualifier')), 
                                                                         polygon=polygon_id))
    }
  }
  return(list(impact=impact, polygons_list=polygons_list))
}

updateODDSubNat <- function(dir, ODDy, event_sdate, event_fdate, subnat_file='EQ_SubNational.xlsx'){
  # Update an existing ODD object (ODDy) using subnational data:
  #   - Edit the 'impact' slot to include subnational data from the provided spreadsheet (default is EQ_subnational.xlsx)
  #   - Edit the 'polygons' slot to identify the pixels belonging to each polygon for which we have impact data
  
  SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData %<>% cleanSubNatData()
  
  # Use event date and iso3 to identify the event in subnat data that correspond to the ODD object
  # (note that all rows for the same event in SubNatData should have the same name and sdate!)
  SubNatData_match <- SubNatData %>% group_by(event_name, sdate) %>% filter(
    length(intersect(iso3,unique(ODDy@data$ISO3C)))>0,
    sdate > (event_sdate - 1) &  fdate < (event_fdate + 1) #LOOSEEND: HAZDATES FOR EXISTING ODD OBJECTS ARE DODGEY
  ) %>% group_split()
  
  #if more than one matching event is found in subnat data then choose manually
  id_chosen <- 1
  if(length(SubNatData_match)>1){
    cat(paste('Please select between the following',length(SubNatData_match),'events by typing the id of the chosen source and pressing return.\n'))
    print('Desired Event:')
    print(paste('Event Date:', event_date, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
    for(i in 1:length(SubNatData_match)){
      print(paste0('Option ', i,':'))
      print(SubNatData_match[[i]][1, c('event_name', 'iso3', 'sdate', 'fdate', 'notes')])
    }
    id_chosen <- as.integer(readline(prompt='id selected: '))
  }
  SubNatEvent <- SubNatData_match[[id_chosen]]
  
    SubNatImpact <- getSubNatImpact(SubNatEvent, subnational=TRUE) #populate the 'impact' slot in ODDy using the national and subnational data
  ODDy@impact <- SubNatImpact$impact
  ODDy <- addODDPolygons(ODDy, SubNatImpact$polygons_list) #populate the 'polygons' slot in ODDy using the regions for which we have impact data
  
  return(ODDy)
}

pixelsInPoly <- function(poly, ODDy){
  
  poly_sp <- SpatialPolygons(poly)
  
  insidepoly<-rep(FALSE,nrow(ODDy@data))
  lon_cellsize <- ODDy@grid@cellsize[1]
  lat_cellsize <- ODDy@grid@cellsize[2]
  
  # Get rid of values outside the bounding box first
  coords <- ODDy@coords
  indies<- which(coords[,1] >= (poly_sp@bbox[1,1]-lon_cellsize/2) &
    coords[,1]<= (poly_sp@bbox[1,2]+ lon_cellsize/2) &
    coords[,2]>= (poly_sp@bbox[2,1]-lat_cellsize/2) &
    coords[,2]<= (poly_sp@bbox[2,2]+lon_cellsize/2))
  
  # Now we only need to calculate a few points that lie inside the polygon!
  for (j in indies){
    #this isn't very efficient but is only done once when reading in the impact data
    pixel <- SpatialPolygons(list(Polygons(list(Polygon(rbind(c(coords[j,1]-lon_cellsize/2, coords[j,2]-lat_cellsize/2),
                           c(coords[j,1]-lon_cellsize/2, coords[j,2]+lat_cellsize/2), 
                           c(coords[j,1]+lon_cellsize/2, coords[j,2]+lat_cellsize/2), 
                           c(coords[j,1]+lon_cellsize/2, coords[j,2]-lat_cellsize/2), 
                           c(coords[j,1]-lon_cellsize/2, coords[j,2]-lat_cellsize/2)
                           ))), 1)))
    insidepoly[j]<-gIntersects(poly_sp, pixel)
  }
  
  return(which(insidepoly == "TRUE"))
}

addODDPolygons <- function(ODDy, polygons_list){
  # Populates the polygons slot in ODDy with a list, where each list element is another list containing:
  #   - polygon_name: the same as the polygon name chosen in getSubNatImpact
  #   - indexes: the indexes of the pixels in ODDy which belong to the polygon
  
  polygons_indexes <- list()
  
  for (i in 1:length(polygons_list)){
    print(i)
    if(!grepl(',', polygons_list[[i]]$polygon_name, fixed = TRUE)){ #check if subnational by searching for comma in polygon name
      if (polygons_list[[i]]$polygon_name=='TOTAL'){ #if 'TOTAL' then add all pixels to the polygon
        inPolyInds <- which(!is.na(tmp$iso3))
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      } else { #if national then add pixels with appropriate ISO3C value to the polygon
        iso3 <- countrycode(polygons_list[[i]]$polygon_name, origin='country.name', destination='iso3c')
        inPolyInds <- which(ODDy@data$ISO3C == iso3)
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      }
    } else {
      if(!is.null(polygons_list[[i]]$sf_polygon)){
        #Find indexes of the spatial pixels inside the polygon and double check that the iso3 codes match. 
        iso3 <- countrycode(sub('.*,\\s*', '', polygons_list[[i]]$polygon_name), origin='country.name', destination='iso3c')
        #inPolyInds <- intersect(which(inPoly((polygons_list[[i]]$sf_polygon)@polygons[[1]], pop = ODDy@coords)$indies), which(ODDy@data$ISO3C==iso3))
        inPolyInds <- pixelsInPoly((polygons_list[[i]]$sf_polygon)@polygons, ODDy)
        if(length(inPolyInds)==0) warning(paste('Region not found in impacted area. Polygon name: ', polygons_list[[i]]$polygon_name))
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      }
    }
  }

  ODDy@polygons <- polygons_indexes
  
  return(ODDy)
  
}

# plot(ODDy@coords[inPolyInds,])
# points(ODDy@coords[which(ODDy$Population>0),], col='green', pch=16)
# points(ODDy@coords[inPolyInds,], col='blue', pch=16)
# points(poly_sp@polygons[[1]]@Polygons[[1]]@coords, type='l')
# points(poly_sp@polygons[[1]]@Polygons[[2]]@coords, type='l')

updateAllODDSubNat <- function(dir, subnat_file='EQ_SubNational.xlsx'){
  # Takes all the ODD objects currently in IIDIPUS_Input/ODDobjects and updates them using the data from the subnational data spreadsheet
  
  folderin<- paste0(dir, 'IIDIPUS_Input/ODDobjects/') #"/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/"
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
  for (ODD_file in ufiles){
    event_date <- as.Date(substr(ODD_file, 3, 10), "%Y%m%d") #seems to be some issues with ODD@hazdate
    ODDy <- readRDS(paste0(folderin, ODD_file))
    ODDy <- updateODDSubNat(dir, ODDy, event_date, subnat_file)
    #saveRDS(ODDy, paste0(folderin, ufiles[i])) 
  }
}

getSubNatData <- function(dir, haz="EQ",extractedData=T, subnat_file= 'EQ_SubNational.xlsx'){
  # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if 
  # no corresponding existing ODD object can be found, creates a new ODD object.
  
  #existingODDloc <- paste0(dir, "IIDIPUS_Input/ODDobjects/")
  existingODDloc <-"/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/"
  existingODDfiles <- na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T)) 
  
  SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData %<>% cleanSubNatData()
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% group_split()
  
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (i in 1:length(SubNatDataByEvent)){
    
    # Subset displacement and disaster database objects
    miniDam<-SubNatDataByEvent[[i]]
    
    miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3), sdate=miniDam$sdate[1], 
                          fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
    
    #CHECK IF EXISTING ODD OBJECT AND UPDATE INSTEAD:
    filename_options <- paste0( miniDam$hazard[1], rep(gsub("-", "", seq(miniDam$sdate[1]-3, miniDam$sdate[1]+3, by='days')), each=length(miniDamSimplified$iso3)), miniDamSimplified$iso3)
    
    if( any(filename_options %in% substr(existingODDfiles, 1, 13))){
      matched_index <- match(filename_options, substr(existingODDfiles, 1, 13))
      file_match <- existingODDfiles[matched_index[which(!is.na(matched_index ))]]
      ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/', file_match))
      print('Updating the event')
      event_date <- as.Date(substr(file_match, 3, 10), "%Y%m%d") 
      print(paste('Event Date:', event_date, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
      print('Using the data:')
      print(head(miniDam))
      
      ODDy <- updateODDSubNat(dir, ODDy, miniDam$sdate[1])
      #save ODDy here
      next
    }
    
    # Set some confining dates for the rest of the data to be assigned to this event
    maxdate<-miniDamSimplified$sdate[1]-5
    if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3

    # Match displacement and hazard data and extract hazard maps
    # HazSDF includes SpatialPixelDataFrame object of hazmean & hazsd per date 
    # (list of each, bbox-cropped to remove M < minmag)
    #lhazSDF<-tryCatch(GetDisaster(miniDam,miniDACS),error=function(e) NULL)
    lhazSDF<-tryCatch(GetDisaster(miniDamSimplified),error=function(e) NULL)
    if(is.null(lhazSDF)) {
      print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
                   " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
      next
    }
    
    # Create the ODD object:
    ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified),error=function(e) NULL)
    if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
    
    ODDy <- updateODDSubNat(dir, ODDy, miniDam$sdate[1])
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDamSimplified$iso3)[1],
                  "_",ODDy@eventid)
    
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_Input/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    HAZARDpath<-paste0(dir,"IIDIPUS_Input/HAZARDobjects/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
    
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDamSimplified$iso3) & 
                               sdate<mindate & sdate>maxdate)
    
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDamSimplified,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",ev, " ",unique(miniDamSimplified$iso3)[1]," ", unique(miniDamSimplified$sdate)[1])) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input/BDobjects/",namer)
      # Save it out!
      saveRDS(BDy, BDpath)
    }
    
    # Path to file and ID
    path%<>%rbind(data.frame(ODDpath=ODDpath,
                             BDpath=BDpath,
                             eventid=ODDy@eventid))
    # Save some RAM
    rm(ODDy,BDy,miniDam)
  }
  
  return(path)
  
}

#download Philippines ODD objects from scratch
getSubNatDataFiltered <- function(dir, haz="EQ",subnat_file= 'EQ_SubNational.xlsx', iso_filter='PHL'){
  # Works through EQ_Subnational.xlsx and creates a new ODD object for each event in The Philippines
  
  SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData %<>% cleanSubNatData()
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% filter(iso_filter %in% iso3) %>% group_split()
  
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (i in 1:length(SubNatDataByEvent)){
    
    # Subset displacement and disaster database objects
    miniDam<-SubNatDataByEvent[[i]]
    
    miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3), sdate=miniDam$sdate[1], 
                                    fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
    
    
    # Set some confining dates for the rest of the data to be assigned to this event
    maxdate<-miniDamSimplified$sdate[1]-5
    if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
    
    # Match displacement and hazard data and extract hazard maps
    # HazSDF includes SpatialPixelDataFrame object of hazmean & hazsd per date 
    # (list of each, bbox-cropped to remove M < minmag)
    #lhazSDF<-tryCatch(GetDisaster(miniDam,miniDACS),error=function(e) NULL)
    #lhazSDF<-tryCatch(GetDisaster(miniDamSimplified),error=function(e) NULL)
    lhazSDF <- GetDisaster(miniDamSimplified)
    if(is.null(lhazSDF)) {
      print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
                   " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
      next
    }
    
    # Create the ODD object:
    #ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified),error=function(e) NULL)
    ODDy <- new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified)
    if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
    
    ODDy_with_impact <- updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1])
    ODDy <- ODDy_with_impact
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDamSimplified$iso3)[1],
                  "_",ODDy@eventid)
    
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"PHL_IIDIPUS_Input/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    HAZARDpath<-paste0(dir,"PHL_IIDIPUS_Input/HAZARDobjects/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    #ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
    
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDamSimplified$iso3) & 
                               sdate<mindate & sdate>maxdate)
    
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDamSimplified,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",ev, " ",unique(miniDamSimplified$iso3)[1]," ", unique(miniDamSimplified$sdate)[1])) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input/BDobjects/",namer)
      # Save it out!
      saveRDS(BDy, BDpath)
    }
    
    # Path to file and ID
    path%<>%rbind(data.frame(ODDpath=ODDpath,
                             BDpath=BDpath,
                             eventid=ODDy@eventid))
    # Save some RAM
    rm(ODDy,BDy,miniDam)
  }
  
  return(path)
  
}

library(randomcoloR)
ncols <- length(ODDy@polygons)
palette <- distinctColorPalette(ncols)
plot(ODDy@coords[ODDy@polygons[[45]]$indexes,], col='black')

for (i in c(1:44)){
  Sys.sleep(1)
  points(ODDy@coords[ODDy@polygons[[i]]$indexes,], col=palette[i])
}


# -----------------------------------------------------------------------------------------

# # # Loop through and check EQ Disagg Data for errors/issues
# # 
# subnat_file='IIDIPUS_Input/EQ_SubNational.xlsx'
# SubNatData <- read.xlsx(paste0(dir, subnat_file), colNames = TRUE , na.strings = c("","NA"))
# SubNatData %<>% cleanSubNatData()
# 
# # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
# SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% group_split()
# 
# for (i in 1:length(SubNatDataByEvent)){
# 
#   miniDam<-SubNatDataByEvent[[i]]
#   print(miniDam[, c('sdate', 'iso3', 'Region', 'Subregion', 'source', 'mortality', 'displacement', 'buildDam', 'buildDest')])
#   miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3), sdate=miniDam$sdate[1],
#                                   fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
#   #miniDamSimplified$sdate <-miniDamSimplified$sdate - 3
#   #miniDamSimplified$fdate <-miniDamSimplified$fdate + 3
#   maxdate<-miniDamSimplified$sdate[1]-5
#   if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
# 
#   # Match displacement and hazard data and extract hazard maps
#   # HazSDF includes SpatialPixelDataFrame object of hazmean & hazsd per date
#   # (list of each, bbox-cropped to remove M < minmag)
#   #lhazSDF<-tryCatch(GetDisaster(miniDam,miniDACS),error=function(e) NULL)
#   lhazSDF<-tryCatch(GetDisaster(miniDamSimplified),error=function(e) NULL)
#   if(is.null(lhazSDF)) {
#     print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
#                  " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
#     next
#   }
# 
#   ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified),error=function(e) NULL)
#   if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
# 
# 
#   plot(ODDy@data$Population, ODDy@data$nBuildings)
# 
#   ODDy <- UpdateODDSubNat(dir, ODDy, miniDam$sdate[1])
# 
#   DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)
# }

# # -------------------------------------------------------------
# 
# # Loop through events and check which regions can be found using getbb()/getGADM()
# 
# subnat_file='/home/manderso/Downloads/EQ_SubNational.xlsx'
# SubNatData <- read.xlsx(subnat_file, colNames = TRUE , na.strings = c("","NA"))
# SubNatData %<>% cleanSubNatData()
# 
# for (i in 616:660){
#   if(!is.na(SubNatData$Subregion[i]) | !is.na(SubNatData$Region[i])){
#     polyname <- paste(SubNatData$Subregion[i], SubNatData$Region[i])
#     getPolyData(polyname, SubNatData$Subregion[i], SubNatData$Region[i], SubNatData$country[i], SubNatData$iso3[i])
#   }
# }



# SubNatData %>% group_by(event_name, sdate) %>% filter(
#   iso3 %in% c('BGD', 'COD', 'DZA', 'IDN', 'MOZ', 'MWI','NPL', 'PHL', 'TLS')
# ) %>% group_split()













