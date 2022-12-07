################
### ODDpolys ###
################

library('openxlsx')

getDisaggImpact <- function(DisaggEvent, disaggregate=FALSE){
  #Use data from EQ_subnational.xlsx to generate data for 'impact' slot in ODD object
  #   - impact has one row per impact and aggregation level
  #   - also creates polygons_list such that polygons_list[[i]] is the sf_polygon corresponding to the polygon with id i. 
  #NEED TO EDIT TO ALLOW FOR DOUBLE COUNTING (allow overlaps!)
  impact <- data.frame(iso3 = character(), polygon = numeric(), impact = character(), 
             observed = numeric(), qualifier = character())
  polygons_list <- list()
  
  DisaggEvent_national <- DisaggEvent %>% filter(is.na(Region) & is.na(Subregion))
  
  if (disaggregate){ #use getbb() to source polygons corresponding to the disaggregated regions
    DisaggEvent_local <- DisaggEvent %>% filter(!is.na(Region) | !is.na(Subregion))
    DisaggEvent_local$polygon_name <- paste0(ifelse(is.na(DisaggEvent_local$Subregion),'',paste0(DisaggEvent_local$Subregion,', ')), 
                                             ifelse(is.na(DisaggEvent_local$Region),'',paste0(DisaggEvent_local$Region,', ')), 
                                             ifelse(DisaggEvent_local$country=='TOTAL','',DisaggEvent_local$country))
    
    DisaggEvent_local %<>% group_by(polygon_name) %>% mutate(polygon_id=cur_group_id()) %>% ungroup()
    for (i in sort(unique(DisaggEvent_local$polygon_id))){
      polygon_name <- DisaggEvent_local$polygon_name[which(DisaggEvent_local$polygon_id == i)[1]]
      polygon <- getbb(polygon_name, format_out='sf_polygon')
      if(NROW(polygon)>1){
        print(paste('Multiple Polygons for', polygon_name))
        polygons_list[[i]] <- polygon[1,]
      } else {
        polygons_list[[i]] <- polygon
      }
      # plot polygons to make sure they all seem reasonable!
    }
    if (length(polygons_list)>0){
      notnans_polygon_id <- which(sapply(polygons_list, function(x) !is.null(x)))
      DisaggEvent_local %<>% filter(polygon_id %in% notnans_polygon_id) 
    } else {
      DisaggEvent_local <- DisaggEvent_local[0,]
    }
  }
  
  for (impact_type in Model$impacts$labels){
    
      notnans_national <- which(!is.na(DisaggEvent_national[,impact_type]))
      notnans_local <- which(!is.na(DisaggEvent_local[,impact_type]))
        if (NROW(DisaggEvent_national[notnans_national,])>0){
        national_source_selected <- c()
        for (i in 1:length(unique(DisaggEvent_national[notnans_national, ]$iso3))){
          notnans_national_iso3filter <- notnans_national[which(DisaggEvent_national[notnans_national, ]$iso3==unique(DisaggEvent_national[notnans_national,]$iso3)[i])]
          if (DisaggEvent_national[notnans_national_iso3filter,]$iso3[1]=='TOT'){
            if ( length(unique(DisaggEvent_national[notnans_national,]$iso3)) > 1 | length(notnans_local) > 0 ){
              next #don't mix local/national and total observations over the same impact type
            }
          }
          national_source_selected <- append(national_source_selected, notnans_national_iso3filter[1])
          if(length(notnans_national_iso3filter)>1){ #if multiple different sources
            if(!all(pull(DisaggEvent_national[notnans_national_iso3filter[1],impact_type])==pull(DisaggEvent_national[notnans_national_iso3filter,impact_type]))){ #if conflicting values, choose one
              cat(paste('Please select between the following',length(notnans_national_iso3filter),'sources, by typing the id of the chosen source and pressing return.\n',
                        'If you would like to select more than one source, please separate each id with a comma\n')) #haven't yet included this functionality
              if(impact_type == 'buildDam' || impact_type == 'buildDest'){
                print(cbind(id=1:length(notnans_national_iso3filter), DisaggEvent_national[notnans_national_iso3filter,c('event_name', 'source_type', 'source_date', 'building_type', impact_type, 'notes', paste0(impact_type,'_term'))]))
              } else {
                print(cbind(id=1:length(notnans_national_iso3filter), DisaggEvent_national[notnans_national_iso3filter,c('event_name', 'source_type', 'source_date', impact_type, 'notes', paste0(impact_type,'_term'))]))
              }
              id_chosen <- as.integer(readline(prompt='id selected: '))
              national_source_selected[length(national_source_selected)] <- notnans_national_iso3filter[id_chosen]
            }
          }
        }
        if(!disaggregate | (length(notnans_local) == 0)){
          impact %<>% rbind(DisaggEvent_national[national_source_selected,] %>% dplyr::select(iso3=iso3,observed=!!sym(impact_type), 
                                                                                              qualifier=paste0(impact_type, '_qualifier')) %>% 
                              cbind(polygon=0, impact=impact_type)) #need to match polygon = 0 with national in ODD object creation
        }
      }
      
      if (disaggregate & (length(notnans_local) > 0)){
        #find which polygons overlap
        intersections <- matrix(FALSE, nrow=length(notnans_local), ncol=length(notnans_local))
        for(i in 1:NROW(intersections)){
          for(j in 1:i){
            intersections[i,j] <- st_intersects(polygons_list[[DisaggEvent_local$polygon_id[notnans_local[i]]]], polygons_list[[DisaggEvent_local$polygon_id[notnans_local[j]]]], sparse=FALSE)
            next
          }
        }
        diag(intersections) <- FALSE
        #remove overlapping polygons
        while(!all(!intersections)){ #while there are still overlapping polygons
          intersect_ij <- which(intersections, arr.ind=TRUE)[1,]
          
          #if polygons are equal then choose manually
          if(DisaggEvent_local$polygon_id[intersect_ij[1]] == DisaggEvent_local$polygon_id[intersect_ij[2]]){
            cat(paste('Please select between the following sources, by typing the id of the chosen source and pressing return.\n',
                      'If you would like to select more than one source, please separate each id with a comma\n')) #haven't yet included this functionality
            if(impact_type == 'buildDam' || impact_type == 'buildDest'){
              print(cbind(id=1:2, DisaggEvent_local[intersections_ij[1],c('polygon', 'source_type', 'source_date', 'building_type', impact_type, 'notes', paste0(impact_type,'_term'))]))
            } else {
              print(cbind(id=1:2, DisaggEvent_local[intersections_ij[2],c('polygon', 'source_type', 'source_date', impact_type, 'notes', paste0(impact_type,'_term'))]))
            }
            id_chosen <- as.integer(readline(prompt='id selected: '))
            source_removed <- intersections_ij[-id_chosen]
            intersections <- intersections[-source_removed,-source_removed]
            notnans_local <- notnans_local[-source_removed]
            next
          }
          smaller_ij <- which.min(as.numeric(st_area(polygons_list[[DisaggEvent_local$polygon_id[intersect_ij[1]]]]),st_area(polygons_list[[DisaggEvent_local$polygon_id[intersect_ij[2]]]])))
          intersections <- intersections[-intersect_ij[smaller_ij],-intersect_ij[smaller_ij]]
          notnans_local <- notnans_local[-intersect_ij[smaller_ij]]
        }
        
        for (i in notnans_local){
          impact %<>% rbind(DisaggEvent_local[i,] %>% transmute(iso3=iso3,
                                                                    observed=!!sym(impact_type), 
                                                                    qualifier=!!sym(paste0(impact_type, '_qualifier')),
                                                                    polygon=polygon_id) %>% 
                              cbind(impact=impact_type)) #need to match polygon = 0 with national in ODD object creation   
        }
        
        #need to match iso3 here
        if (length(notnans_national)>0){
          for (iso3_filter in DisaggEvent_national[national_source_selected,]$iso3){
            impact %<>% rbind(DisaggEvent_national[national_source_selected,] %>% filter(iso3==iso3_filter) %>% 
                                                  transmute(iso3=iso3,
                                                            observed=!!sym(impact_type)-sum(impact %>% 
                                                                                              filter((polygon>0) & (impact== impact_type) & (iso3==iso3_filter)) %>% 
                                                                                              pull(observed)), 
                                                                                            qualifier=!!sym(paste0(impact_type, '_qualifier'))) %>% 
                                cbind(polygon=-1, impact=impact_type))
            if(any((impact %>% filter(impact==impact_type) %>% pull(observed)) < 0)){
              warning(paste('Sum of disaggregated impacts for', impact_type ,'is larger than aggregated'))
            }
          }
        }
     }
  }
  return(list(impact=impact, polygons_list=polygons_list))
}

updateODDdisagg <- function(dir, ODDy, event_date, disagg_file='IIDIPUS_Input/EQ_Subnational.xlsx'){
  #Update ODD object using subnational data:
  #   - Edit 'impact' slot to include subnational data from EQ_subnational.xlsx
  #   - Adds columns for each impact (e.g. polyMort) to ODD@data which identifies the polygon to which each pixel belongs for that impact
  DisaggData <- read.xlsx(paste0(dir, disagg_file), colNames = TRUE , na.strings = c("","NA"))
  DisaggData$source_date <- openxlsx::convertToDate(DisaggData$source_date)
  DisaggData$sdate <- openxlsx::convertToDate(DisaggData$sdate)
  DisaggData$fdate <- openxlsx::convertToDate(DisaggData$fdate)
  DisaggData$mortality <- as.integer(DisaggData$mortality)
  DisaggData$displacement <- as.integer(DisaggData$displacement)
  DisaggData$buildDam <- as.integer(DisaggData$buildDam)
  DisaggData$buildDest <- as.integer(DisaggData$buildDest)
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  DisaggData_match <- DisaggData %>% group_by(event_name, sdate) %>% filter(
    length(intersect(iso3,unique(ODDy@data$ISO3C)))>0,
    sdate > (event_date - 3) &  fdate < (event_date + 3) #NEED TO FIX HAZDATES
  ) %>% group_split()
  
  id_chosen <- 1
  if(length(DisaggData_match)>1){
    cat(paste('Please select between the following',length(DisaggData_match),'events by typing the id of the chosen source and pressing return.\n'))
    print('Desired Event:')
    print(paste('Event Date:', event_date, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
    for(i in 1:length(DisaggData_match)){
      print(paste0('Option ', i,':'))
      print(DisaggData_match[[i]][1, c('event_name', 'iso3', 'sdate', 'fdate', 'notes')])
    }
    id_chosen <- as.integer(readline(prompt='id selected: '))
  }
  DisaggEvent <- DisaggData_match[[id_chosen]]
  
  DisaggImpact <- getDisaggImpact(DisaggEvent, disaggregate=TRUE)
  ODDy@impact <- DisaggImpact$impact
  #ODDy@impact$observed <- as.integer(ODDy@impact$observed)
  ODDy <- addODDPolygons(ODDy, DisaggImpact$polygons_list)
  
  return(ODDy)
}

addODDPolygons <- function(ODDy, polygons_list){
  #Input: 
  #  - ODD Object containing a slot named 'impact' with rows corresponding to different impacts in different polygons
  # Output
  #  - ODD Object with each pixel labelled with the polygon ids to which it is assigned for each impact
  
  ODDy@data[,Model$impacts$polynames[which(Model$impacts$labels %in% ODDy@impact$impact)]] <- -1 #initially allocate to the 'remainder'/'total'
  
  for (i in 1:NROW(ODDy@impact)){
    impact_type <- ODDy@impact[i,'impact']
    polygon_id <- ODDy@impact[i,'polygon']
    poly_label <- Model$impacts$polynames[which(Model$impacts$labels==impact_type)]
    if (polygon_id > 0){
      inPolyInds <- inPoly((polygons_list[[polygon_id]] %>% as("Spatial"))@polygons[[1]], pop = ODDy@coords)
      if(length(inPolyInds)==0) warning(paste('Region not found in impacted area. Polygon id: ', polygon_id))
      if(any(ODDy@data[inPolyInds,poly_label] != -1)) warning('Trying to assign a pixel to two polygons for the same impact')
      ODDy@data[inPolyInds,poly_label] <- ODDy@impact[i,'polygon']
    } else if (polygon_id == 0){
      inPolyInds <- which(ODDy@data$ISO3C==ODDy@impact[i,'iso3'])
      if(any(ODDy@data[inPolyInds,poly_label] != -1)) warning('Trying to assign a pixel to two polygons for the same impact')
      ODDy@data[inPolyInds,poly_label] <- ODDy@impact[i,'polygon']
    }
  }
  
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

updateAllODDdisagg <- function(dir){
  #updates all the current ODD objects:
  
  #folderin<-paste0(dir, "IIDIPUS_Input/ODDobjects/")
  folderin<-"/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/"
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
  for (ODD_file in ufiles[4]){
    event_date <- as.Date(substr(ODD_file, 3, 10), "%Y%m%d") #seems to be some issues with ODD@hazdate
    #ODDy <- readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/", ODD_file))
    ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/', ODD_file))
    ODDy <- updateODDdisagg(dir, ODDy, event_date)
  }
}

getDisaggData <- function(dir, haz="EQ",extractedData=T){
  #works through EQ_Subnational.xlsx and, for each event, either updates the current ODD object or creates a new ODD object
  
  #existingODDloc <- paste0(dir, "IIDIPUS_Input/ODDobjects/")
  existingODDloc <-"/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/"
  existingODDfiles <- na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T)) 
  
  disagg_file='IIDIPUS_Input/EQ_Subnational.xlsx'
  DisaggData <- read.xlsx(paste0(dir, disagg_file), colNames = TRUE , na.strings = c("","NA"))
  DisaggData$source_date <- openxlsx::convertToDate(DisaggData$source_date)
  DisaggData$sdate <- openxlsx::convertToDate(DisaggData$sdate)
  DisaggData$fdate <- openxlsx::convertToDate(DisaggData$fdate)
  DisaggData$mortality <- as.integer(DisaggData$mortality)
  DisaggData$displacement <- as.integer(DisaggData$displacement)
  DisaggData$buildDam <- as.integer(DisaggData$buildDam)
  DisaggData$buildDest <- as.integer(DisaggData$buildDest)
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  DisaggDataByEvent <- DisaggData %>% group_by(event_name, sdate) %>% group_split()
  
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (i in 1:length(DisaggDataByEvent)){
    
    # Subset displacement and disaster database objects
    miniDam<-DisaggDataByEvent[[i]]
    
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
      
      ODDy <- updateODDdisagg(dir, ODDy, miniDam$sdate[1])
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
    
    ODDy <- updateODDdisagg(dir, ODDy, miniDam$sdate[1])
    
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
disagg_file='IIDIPUS_Input/EQ_Subnational.xlsx'
DisaggData <- read.xlsx(paste0(dir, disagg_file), colNames = TRUE , na.strings = c("","NA"))
DisaggData$source_date <- openxlsx::convertToDate(DisaggData$source_date)
DisaggData$sdate <- openxlsx::convertToDate(DisaggData$sdate)
DisaggData$fdate <- openxlsx::convertToDate(DisaggData$fdate)
DisaggData$mortality <- as.integer(DisaggData$mortality)
DisaggData$displacement <- as.integer(DisaggData$displacement)
DisaggData$buildDam <- as.integer(DisaggData$buildDam)
DisaggData$buildDest <- as.integer(DisaggData$buildDest)

# Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
DisaggDataByEvent <- DisaggData %>% group_by(event_name, sdate) %>% group_split()

for (i in 1:length(DisaggDataByEvent)){
  
  miniDam<-DisaggDataByEvent[[i]]
  print(miniDam[, c('sdate', 'iso3', 'Region', 'Subregion', 'source', 'mortality', 'displacement', 'buildDam', 'buildDest')])
  miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3), sdate=miniDam$sdate[1], 
                                  fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
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
  
  ODDy <- updateODDdisagg(dir, ODDy, miniDam$sdate[1])
  
  DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)
}

# -------------------------------------------------------------

