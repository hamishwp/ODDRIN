################
### ODDpolys ###
################

library('openxlsx')

cleanSubNatData <- function(SubNatData){
  SubNatData$source_date <- openxlsx::convertToDate(SubNatData$source_date)
  SubNatData$sdate <- openxlsx::convertToDate(SubNatData$sdate)
  SubNatData$fdate <- openxlsx::convertToDate(SubNatData$fdate)
  SubNatData$mortality <- as.integer(SubNatData$mortality)
  SubNatData$displacement <- as.integer(SubNatData$displacement)
  SubNatData$buildDam <- as.integer(SubNatData$buildDam)
  SubNatData$buildDest <- as.integer(SubNatData$buildDest)
  return(SubNatData)
}

getSubNatImpact <- function(SubNatEvent, subnational=TRUE){
  #Use data from EQ_subnational.xlsx to generate data for 'impact' slot in ODD object
  #   - impact has one row per impact and aggregation level
  #   - also creates polygons_list such that polygons_list[[i]] is the sf_polygon corresponding to the polygon with id i. 

  impact <- data.frame(iso3 = character(), polygon = numeric(), impact = character(), 
                       observed = numeric(), qualifier = character())
  
  polygons_list <- list()
  nans_polygon_id <- c()
  
  if (!subnational){
    SubNatEvent <- SubNatEvent %>% filter(is.na(Region) & is.na(Subregion))
    i = 1
    for (iso3 in unique(SubNatEvent$iso3)){
      SubNatEvent$polygon_id <- i
      polygons_list[[i]] <- list(polygon_name=countrycode(iso3, origin='iso3c', destination='country.name'), sf_polygon=NULL)
      i = i + 1
    }
  } else { 
    SubNatEvent$polygon_name <- paste0(ifelse(is.na(SubNatEvent$Subregion),'',paste0(SubNatEvent$Subregion,', ')), 
                                             ifelse(is.na(SubNatEvent$Region),'',paste0(SubNatEvent$Region,', ')), 
                                             ifelse(SubNatEvent$country=='TOTAL','',SubNatEvent$country))
    
    SubNatEvent %<>% group_by(polygon_name) %>% mutate(polygon_id=cur_group_id()) %>% ungroup()
    
    for (i in sort(unique(SubNatEvent$polygon_id))){
      polygon_name <- SubNatEvent$polygon_name[which(SubNatEvent$polygon_id == i)[1]]
      if (grepl(',', polygon_name, fixed = TRUE)){ #check if subnational by searching for comma in polygon name
        polygon <- getbb(polygon_name, format_out='sf_polygon')
        if(NROW(polygon)>1){
          print(paste('Multiple Polygons for', polygon_name))
          polygons_list[[i]] <- list(polygon_name = polygon_name, sf_polygon = polygon[1,])
        } else {
          polygons_list[[i]] <- list(polygon_name = polygon_name, sf_polygon = polygon)
        }
        if (is.null(polygons_list[[i]]$sf_polygon)) nans_polygon_id %<>% append(i)
        
        # should plot polygons here to make sure they all seem reasonable!
      } else {
        polygons_list[[i]] <- list(polygon_name = polygon_name, sf_polygon = NULL)
      }
    }
  }
  
  SubNatEvent %<>% filter(!(polygon_id %in% nans_polygon_id))
  
  for (impact_type in Model$impacts$labels){
    
    notnans <- which(!is.na(SubNatEvent[,impact_type]))
    
    if (length(notnans)>0){
      SubNatEvent[notnans,'source_id'] <- notnans
      sources_selected <- c()
      SubNatEvent_by_polygon <- SubNatEvent[notnans,] %>% group_by(polygon_id) %>% group_split()
      for (i in 1:length(SubNatEvent_by_polygon)){
        if(NROW(SubNatEvent_by_polygon[[i]])==1){#only one source for that polygon
          sources_selected %<>% append(SubNatEvent_by_polygon[[i]]$source_id) 
          next
        } 
        if(length(unique(pull(SubNatEvent_by_polygon[[i]][,impact_type])))==1){ #all matching, so just take first
          sources_selected %<>% append(SubNatEvent_by_polygon[[i]]$source_id[1])
          next
        }
        cat(paste('Please select between the following sources, by typing the id of the chosen source and pressing return.\n',
                  'If you would like to select more than one source, please separate each id with a comma and we will use the mean \n')) 
        if(impact_type == 'buildDam' || impact_type == 'buildDest'){
          print(cbind(id=1:NROW(SubNatEvent_by_polygon[[i]]), SubNatEvent_by_polygon[[i]][,c('event_name', 'source_type', 'source_date', 'building_type', impact_type, 'notes', paste0(impact_type,'_term'))]))
        } else {
          print(cbind(id=1:NROW(SubNatEvent_by_polygon[[i]]), SubNatEvent_by_polygon[[i]][,c('event_name', 'source_type', 'source_date', impact_type, 'notes', paste0(impact_type,'_term'))]))
        }
        ids_chosen <- as.integer(unlist(strsplit(readline(prompt='id selected: '), ",")))
        
        if (length(ids_chosen)>1){
          observed_mean <- round(mean(pull(SubNatEvent_by_polygon[[i]][ids_chosen, impact_type])))
          SubNatEvent[SubNatEvent_by_polygon[[i]]$source_id[ids_chosen[1]], impact_type] <- observed_mean
          ids_chosen <- ids_chosen[1]
        } 
        sources_selected %<>% append(SubNatEvent_by_polygon[[i]]$source_id[ids_chosen])
      }
      impact %<>% rbind(SubNatEvent[sources_selected,] %>% transmute(iso3=iso3,
                                                                         impact=impact_type,
                                                                         observed=!!sym(impact_type), 
                                                                         qualifier=!!sym(paste0(impact_type, '_qualifier')), 
                                                                         polygon=polygon_id))
    }
  }
  return(list(impact=impact, polygons_list=polygons_list))
}

updateODDSubNat <- function(dir, ODDy, event_date, subnat_file='IIDIPUS_Input/EQ_SubNational.xlsx'){
  #Update ODD object using subnational data:
  #   - Edit 'impact' slot to include subnational data from EQ_subnational.xlsx
  #   - Adds columns for each impact (e.g. polyMort) to ODD@data which identifies the polygon to which each pixel belongs for that impact
  SubNatData <- read.xlsx(paste0(dir, subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData %<>% cleanSubNatData()
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatData_match <- SubNatData %>% group_by(event_name, sdate) %>% filter(
    length(intersect(iso3,unique(ODDy@data$ISO3C)))>0,
    sdate > (event_date - 3) &  fdate < (event_date + 3) #NEED TO FIX HAZDATES
  ) %>% group_split()
  
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
  
  SubNatImpact <- getSubNatImpact(SubNatEvent, subnational=TRUE)
  ODDy@impact <- SubNatImpact$impact
  #ODDy@impact$observed <- as.integer(ODDy@impact$observed)
  ODDy <- addODDPolygons(ODDy, SubNatImpact$polygons_list)
  
  return(ODDy)
}

addODDPolygons <- function(ODDy, polygons_list){
  #Input: 
  #  - ODD Object containing a slot named 'impact' with rows corresponding to different impacts in different polygons
  # Output
  #  - 
  polygons_indexes <- list()
  
  for (i in 1:length(polygons_list)){
    if(!grepl(',', polygons_list[[i]]$polygon_name, fixed = TRUE)){ #check if subnational by searching for comma in polygon name
      if (polygons_list[[i]]$polygon_name=='TOTAL'){
        inPolyInds <- which(!is.na(tmp$iso3))
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      } else {
        iso3 <- countrycode(polygons_list[[i]]$polygon_name, origin='country.name', destination='iso3c')
        inPolyInds <- which(ODDy@data$ISO3C == iso3)
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      }
    } else {
      if(!is.null(polygons_list[[i]]$sf_polygon)){
        inPolyInds <- inPoly((polygons_list[[i]]$sf_polygon %>% as("Spatial"))@polygons[[1]], pop = ODDy@coords)
        if(length(inPolyInds)==0) warning(paste('Region not found in impacted area. Polygon name: ', polygons_list[[i]]$polygon_name))
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      }
    }
  }
  #what about national?
  # } else if (polygon_id == 0){
  #   inPolyInds <- which(ODDy@data$ISO3C==ODDy@impact[i,'iso3'])
  #   if(any(ODDy@data[inPolyInds,poly_label] != -1)) warning('Trying to assign a pixel to two polygons for the same impact')
  #   ODDy@data[inPolyInds,poly_label] <- ODDy@impact[i,'polygon']
  # }

  ODDy@polygons <- polygons_indexes
  
  return(ODDy)
  
}

inPoly<-function(poly, pop, iii = 1, sumFn = "sum"){
  #determines which pixels are inside a polygon
  if(any(class(pop) == "SpatialPointsDataFrame") | any(class(pop) == "SpatialPixelsDataFrame")){
    coords<-pop@coords
    data<-pop@data
  } else {
    coords<-as.data.frame(pop[,c("Longitude","Latitude")])
    data<-as.data.frame(pop)
  }
  
  insidepoly<-rep(FALSE,nrow(pop))
  
  for (i in 1:length(poly@Polygons)){
    # Get rid of values outside the bounding box first
    minipoly<-rep(FALSE,length(insidepoly))
    indies<-coords[,1]>=min(poly@Polygons[[i]]@coords[,1]) &
      coords[,1]<=max(poly@Polygons[[i]]@coords[,1]) &
      coords[,2]>=min(poly@Polygons[[i]]@coords[,2]) &
      coords[,2]<=max(poly@Polygons[[i]]@coords[,2])
    # Now we only need to calculate a few points that lie inside the polygon!
    minipoly[indies]<-sp::point.in.polygon(coords[indies,1],
                                           coords[indies,2],
                                           poly@Polygons[[i]]@coords[,1],
                                           poly@Polygons[[i]]@coords[,2])>0
    # Add to the total
    insidepoly<- insidepoly | minipoly
  }
  #outer<-match.fun(sumFn)(data[insidepoly,iii],na.rm=T)
  #return(list(vals=outer,indies=insidepoly))
  return(indies = which(insidepoly == "TRUE"))
}

updateAllODDSubNat <- function(dir){
  #updates all the current ODD objects:
  
  #folderin<-paste0(dir, "IIDIPUS_Input/ODDobjects/")
  folderin<-"/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/"
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
  for (ODD_file in ufiles[4]){
    event_date <- as.Date(substr(ODD_file, 3, 10), "%Y%m%d") #seems to be some issues with ODD@hazdate
    #ODDy <- readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/", ODD_file))
    ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/', ODD_file))
    ODDy <- updateODDSubNat(dir, ODDy, event_date)
  }
}

getSubNatData <- function(dir, haz="EQ",extractedData=T){
  #works through EQ_Subnational.xlsx and, for each event, either updates the current ODD object or creates a new ODD object
  
  #existingODDloc <- paste0(dir, "IIDIPUS_Input/ODDobjects/")
  existingODDloc <-"/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/"
  existingODDfiles <- na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T)) 
  
  subnat_file='IIDIPUS_Input/EQ_SubNational.xlsx'
  SubNatData <- read.xlsx(paste0(dir, subnat_file), colNames = TRUE , na.strings = c("","NA"))
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
      
      ODDy <- updateODDsubnat(dir, ODDy, miniDam$sdate[1])
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
    
    ODDy <- updateODDsubnat(dir, ODDy, miniDam$sdate[1])
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDam$iso3)[1],
                  "_",ODDy@eventid)
    
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_Input/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    HAZARDpath<-paste0(dir,"IIDIPUS_Input/HAZARDobjects/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
    
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDam$iso3) & 
                               sdate<mindate & sdate>maxdate)
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
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

# -----------------------------------------------------------------------------------------

# Loop through and check EQ Disagg Data for errors/issues
subnat_file='IIDIPUS_Input/EQ_SubNational.xlsx'
SubNatData <- read.xlsx(paste0(dir, subnat_file), colNames = TRUE , na.strings = c("","NA"))
SubNatData %<>% cleanSubNatData()

# Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% group_split()

for (i in 1:length(SubNatDataByEvent)){
  
  miniDam<-SubNatDataByEvent[[i]]
  print(miniDam[, c('sdate', 'iso3', 'Region', 'Subregion', 'source', 'mortality', 'displacement', 'buildDam', 'buildDest')])
  miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3), sdate=miniDam$sdate[1], 
                                  fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
  #miniDamSimplified$sdate <-miniDamSimplified$sdate - 3
  #miniDamSimplified$fdate <-miniDamSimplified$fdate + 3
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
  
  ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified),error=function(e) NULL)
  if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
  
  
  plot(ODDy@data$Population, ODDy@data$nBuildings)
  
  ODDy <- updateODDSubNat(dir, ODDy, miniDam$sdate[1])
  
  DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)
}

# -------------------------------------------------------------

