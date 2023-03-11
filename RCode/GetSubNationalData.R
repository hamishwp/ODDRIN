################
### ODDpolys ###
################



library(openxlsx)
library(tidyxl)
library(rgeos)
library(sfheaders)

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

readSubNatData <- function(subnat_file){
  # Clean data from xlsx file by converting data to the appropriate formats/data types
  SubNatData <-  read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  x <- xlsx_cells(paste0(dir, 'IIDIPUS_Input/', subnat_file))
  formats <- xlsx_formats(paste0(dir, 'IIDIPUS_Input/', subnat_file))

  red_cells <- x %>% filter(local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == "FFFF0000")) %>% dplyr::select(row, col)
  red_cells$row <- red_cells$row - 1
  red_cells <- red_cells[-which(red_cells$col > NCOL(SubNatData)), ]
  SubNatData[as.matrix(red_cells)] <- NA
  
  SubNatData <- SubNatData[rowSums(is.na(SubNatData)) != ncol(SubNatData), ]
  
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
  if (country=='TOTAL'){
    return(list(polygon_name = polygon_name, sf_polygon = NULL))
  }
  
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
    polygon <- tryCatch(getbb(paste(GADM_array, collapse=','), format_out='sf_polygon'), error=function(e) NULL) 
    if (!is.null(polygon)){
      if(!is.null(polygon$multipolygon)) polygon <- polygon$multipolygon
      if (NROW(polygon)> 1){
        warning(paste('Multiple Polygons for', polygon_name))
        polygon <- polygon[1,]
      }
      polygon <- polygon %>% as("Spatial")
      print('Polygon found using getbb() but not getGADM()') 
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
      polygons_list[[i]] <- list(polygon_name=ifelse(iso3='TOT', 'TOTAL', countrycode(iso3, origin='iso3c', destination='country.name')), sf_polygon=NULL) #sf_polygon not required for countries
      i = i + 1
    }
  } else { #work with the subnational data
    SubNatEvent$polygon_name <- paste0(ifelse(is.na(SubNatEvent$Subregion),'',paste0(SubNatEvent$Subregion,', ')), 
                                             ifelse(is.na(SubNatEvent$Region),'',paste0(SubNatEvent$Region,', ')), 
                                             ifelse(SubNatEvent$country=='TOTAL','TOTAL',SubNatEvent$country))
    
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
  
  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)
  
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
    print(paste('Event Date:', event_sdate, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
    for(i in 1:length(SubNatData_match)){
      print(paste0('Option ', i,':'))
      print(SubNatData_match[[i]][1, c('event_name', 'iso3', 'sdate', 'fdate', 'notes')])
    }
    id_chosen <- as.integer(readline(prompt='id selected: '))
  }
  SubNatEvent <- SubNatData_match[[id_chosen]]
  
  SubNatImpact <- getSubNatImpact(SubNatEvent, subnational=TRUE) #populate the 'impact' slot in ODDy using the national and subnational data
  
  #remove overlapping polygons using SubNatImpact$impact and SubNatImpact$polygons_list
  #find the polygons and complements we want to use for each impact type and add to ODDy@impact
  
  ODDy <- addODDPolygons(ODDy, SubNatImpact$polygons_list) #populate the 'polygons' slot in ODDy using the regions for which we have impact data
  ODDy <- addODDImpact(ODDy, SubNatImpact$impact)
  
  return(ODDy)
}

pixelsInPoly <- function(poly, ODDy){
  poly_sp <- SpatialPolygons(poly)
  proj4string(poly_sp)<- ODDy@proj4string
  poly_sf <- st_as_sf(poly_sp)
  
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
    # pixel <- SpatialPolygons(list(Polygons(list(Polygon(rbind(c(coords[j,1]-lon_cellsize/2, coords[j,2]-lat_cellsize/2),
    #                        c(coords[j,1]-lon_cellsize/2, coords[j,2]+lat_cellsize/2), 
    #                        c(coords[j,1]+lon_cellsize/2, coords[j,2]+lat_cellsize/2), 
    #                        c(coords[j,1]+lon_cellsize/2, coords[j,2]-lat_cellsize/2), 
    #                        c(coords[j,1]-lon_cellsize/2, coords[j,2]-lat_cellsize/2)
    #                        ))), 1)), proj4string=ODDy@proj4string)
    # pixel_sf <- st_as_sf(pixel)
    
    pixel_sf <- sfheaders::sf_polygon(data.frame(longitude = c(coords[j,1]-lon_cellsize/2, 
                                                   coords[j,1]-lon_cellsize/2,
                                                   coords[j,1]+lon_cellsize/2, 
                                                   coords[j,1]+lon_cellsize/2, 
                                                   coords[j,1]-lon_cellsize/2), 
                                     latitude = c(coords[j,2]-lat_cellsize/2, 
                                                  coords[j,2]+lat_cellsize/2,
                                                  coords[j,2]+lat_cellsize/2, 
                                                  coords[j,2]-lat_cellsize/2,
                                                  coords[j,2]-lat_cellsize/2)), x='longitude',y='latitude')
    sf::st_crs(pixel_sf) = ODDy@proj4string
    
    if(st_intersects(pixel_sf, poly_sf, sparse=F)[1,1]){
      insidepoly[j] <- T
    }
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
        inPolyInds <- which(!is.na(ODDy$ISO3C))
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      } else { #if national then add pixels with appropriate ISO3C value to the polygon
        iso3 <- countrycode(polygons_list[[i]]$polygon_name, origin='country.name', destination='iso3c')
        inPolyInds <- which(ODDy@data$ISO3C == iso3)
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      }
    } else {
      if(!is.null(polygons_list[[i]]$sf_polygon)){
        #inPolyInds <- intersect(which(inPoly((polygons_list[[i]]$sf_polygon)@polygons[[1]], pop = ODDy@coords)$indies), which(ODDy@data$ISO3C==iso3))
        inPolyInds <- intersect(pixelsInPoly((polygons_list[[i]]$sf_polygon)@polygons, ODDy), which(!is.na(ODDy$ISO3C)))
        if(length(inPolyInds)==0) warning(paste('Region not found in impacted area. Polygon name: ', polygons_list[[i]]$polygon_name))
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      }
    }
  }
  
  polygons_indexes <- remove_overlaps(polygons_indexes, polygons_list, ODDy@coords)
  ODDy@polygons <- polygons_indexes
  
  return(ODDy)
  
}

remove_overlaps <- function(polygons_indexes, polygons_list, coords){
  #find which polygons have intersecting indexes but not intersecting sf polygons
  #this may be the case for pixels lying on the border which have been allocated to two polygons
  #in this case, allocate the pixel to the polygon in which its center lies
  if (length(polygons_indexes) == 1) return(polygons_indexes)
  
  for(i in 2:length(polygons_indexes)){
    for(j in 1:(i-1)){
      if (length(polygons_indexes[[i]]$indexes) == 0 | length(polygons_indexes[[j]]$indexes) == 0) next
      #WHAT TO DO IF POLYGON OR ISO3
      if (!grepl(',', polygons_list[[i]]$polygon_name, fixed = TRUE) | !grepl(',', polygons_list[[j]]$polygon_name, fixed = TRUE)){
        next
      } 
      polygons_intersect <- gIntersects(polygons_list[[i]]$sf_polygon, polygons_list[[j]]$sf_polygon)
      polygons_touch <- gTouches(polygons_list[[i]]$sf_polygon, polygons_list[[j]]$sf_polygon)
      indexes_intersect <- intersect(polygons_indexes[[i]]$indexes, polygons_indexes[[j]]$indexes)
      
      if((!polygons_intersect | polygons_touch) & (length(indexes_intersect) > 0) ){
        print('intersecting polys')
        allocated_poly_i <- c()
        allocated_poly_j <- c()
        for (index in indexes_intersect){
          spatial_points <- SpatialPoints(coords, proj4string = CRS(proj4string(ODDy)))
          if (gDistance(spatial_points[index], polygons_list[[i]]$sf_polygon) < 
              gDistance(spatial_points[index], polygons_list[[j]]$sf_polygon)){
            allocated_poly_i %<>% append(index)
          } else {allocated_poly_j %<>% append(index)}
        }
        polygons_indexes[[i]]$indexes <- polygons_indexes[[i]]$indexes[! polygons_indexes[[i]]$indexes %in% allocated_poly_j]
        polygons_indexes[[j]]$indexes <- polygons_indexes[[j]]$indexes[! polygons_indexes[[j]]$indexes %in% allocated_poly_i]
      }
    }
  }

  return(polygons_indexes)
  
}

addODDImpact <- function(ODDy, impact){
  impact_updated <- impact[0,]
  impact_by_type <- impact %>% group_by(impact) %>% group_split()
  for (k in 1:length(impact_by_type)){
    impact_selected <- impact_by_type[[k]]
    if (NROW(impact_selected) == 1){
      impact_updated %<>% rbind(impact_selected)
      next
    } 
    intersections <- matrix(FALSE, nrow=length(ODDy@polygons), ncol=length(ODDy@polygons))
    for(i in impact_selected$polygon){
      for(j in impact_selected$polygon[which(impact_selected$polygon<i)]){
        intersections[i,j] <- length(intersect(ODDy@polygons[[i]]$indexes, ODDy@polygons[[j]]$indexes)) > 0
      }
    }
    while(!all(!intersections)){ #while there are still overlapping polygons
      intersect_ij <- which(intersections, arr.ind=TRUE)[1,]
      i <- intersect_ij[1]; j <- intersect_ij[2]
      poly_i_indexes <- ODDy@polygons[[i]]$indexes
      poly_j_indexes <- ODDy@polygons[[j]]$indexes
      intersecting_indexes <- intersect(poly_i_indexes, poly_j_indexes)
      smaller_ij <- c(i,j)[which.min(c(length(poly_i_indexes),length(poly_j_indexes)))] #find smaller region
      larger_ij <- c(i,j)[which(!c(i,j) %in% smaller_ij)] #find larger region
      if (length(intersecting_indexes)!= length(ODDy@polygons[[smaller_ij]]$indexes)){
        stop('One polygon does not completely consume the other')
      }
      complement_name <- paste0(ODDy@polygons[[larger_ij]]$name, ' excluding (', ODDy@polygons[[smaller_ij]]$name,')')
      complement_indexes <- ODDy@polygons[[larger_ij]]$indexes[which(! ODDy@polygons[[larger_ij]]$indexes %in% intersecting_indexes)]
      complement_value <- impact_selected$observed[which(impact_selected$polygon==larger_ij)] - impact_selected$observed[which(impact_selected$polygon==smaller_ij)]
      #find if polygon exists
      if (complement_name %in% sapply(ODDy@polygons, function(x) x$name)){
        new_polygon_id <- which(complement_name == sapply(ODDy@polygons, function(x) x$name))
      } else {
        new_polygon_id <- length(ODDy@polygons) + 1
        ODDy@polygons[[new_polygon_id]] <- list(name=complement_name, indexes=complement_indexes)
      }
      impact_selected$observed[which(impact_selected$polygon==larger_ij)] <- complement_value
      impact_selected$polygon[which(impact_selected$polygon==larger_ij)] <- new_polygon_id
      
      intersections <- matrix(FALSE, nrow=length(ODDy@polygons), ncol=length(ODDy@polygons))
      for(i in impact_selected$polygon){
        for(j in impact_selected$polygon[which(impact_selected$polygon<i)]){
          intersections[i,j] <- length(intersect(ODDy@polygons[[i]]$indexes, ODDy@polygons[[j]]$indexes)) > 0
        }
      }
    }
    indexes_impact <- list()
    for (i in NROW(impact_selected)){
      indexes_impact[[i]] <- ODDy@polygons[[impact_selected$polygon[i]]]$indexes
    }
    if (length(Reduce(intersect, indexes_impact)) > 0){
      stop('Overlapping polygons')
    }
    impact_updated %<>% rbind(impact_selected)
  }
  ODDy@impact <- impact_updated
  
  #remove obsolete polygons
  polygons_obsolete <- which(!1:length(ODDy@polygons) %in% ODDy@impact$polygon)
  id_matches <- data.frame(polygon= (1:length(ODDy@polygons))[! 1:length(ODDy@polygons) %in% polygons_obsolete],
                           polygon_new= 1:(length(ODDy@polygons)-length(polygons_obsolete)))
  
  ODDy@impact %<>% merge(id_matches)
  ODDy@impact$polygon <- NULL
  colnames(ODDy@impact)[which(colnames(ODDy@impact)=='polygon_new')] <- 'polygon'
  ODDy@polygons <- ODDy@polygons[which(! 1:length(ODDy@polygons) %in% polygons_obsolete)]
    
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
  existingODDloc <-"/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_BuildingDat/ODDobjects"
  existingODDfiles <- na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T)) 
  existingODDloc <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NoBuildingDat/ODDobjects'
  existingODDfiles <- append(existingODDfiles, na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T)) )
  
  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% group_split()
  
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  
  which(sapply(SubNatDataByEvent, function(x) return('TWN' %in% x$iso3)))
  
  for (i in 1:length(SubNatDataByEvent)){
    
    # Subset displacement and disaster database objects
    miniDam<-SubNatDataByEvent[[i]]
    print(miniDam)
    miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3), sdate=miniDam$sdate[1], 
                          fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
    
    #CHECK IF EXISTING ODD OBJECT AND UPDATE INSTEAD:
    filename_options <- paste0( miniDam$hazard[1], rep(gsub("-", "", seq(miniDam$sdate[1]-2, miniDam$sdate[1]+2, by='days')), each=length(miniDamSimplified$iso3)), miniDamSimplified$iso3)
    
    if( any(filename_options %in% substr(existingODDfiles, 1, 13))){
      matched_index <- match(filename_options, substr(existingODDfiles, 1, 13))
      file_match <- existingODDfiles[matched_index[which(!is.na(matched_index ))]]
      ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_BuildingDat/ODDobjects/', file_match))
      print('Updating the event')
      event_date <- as.Date(substr(file_match, 3, 10), "%Y%m%d") 
      print(paste('Event Date:', event_date, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
      print('Using the data:')
      print(head(miniDam))
      
      #ODDy <- updateODDSubNat(dir, ODDy, miniDam$sdate[1])
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
    bbox<-countriesbbox(unique(miniDamSimplified$iso3))
    bbox[1] <- max(bbox[1]-5, -180)
    bbox[2] <- max(bbox[2]-5, -90)
    bbox[3] <- min(bbox[3]+5, 180)
    bbox[4] <- min(bbox[4]+5, 90)
    # lhazSDF <- GetDisaster(miniDamSimplified, bbox = c(160, -11.8, 163, -9.7))
    lhazSDF<-GetDisaster(miniDamSimplified,bbox=bbox) #, bbox=c(140,-9,147.5,-3)) #, EQparams=list(I0=4.5,minmag=4.5))
    
    if(is.null(lhazSDF)) {
      print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
                   " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
      next
    }
    
    # Create the ODD object:
    ODDy<-new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified)
    if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
    
    ODDy_with_impact <- updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1])
    ODDy <- ODDy_with_impact
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDamSimplified$iso3)[1],
                  "_",ODDy@eventid)
    
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_Input_NoBuildingDat/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    HAZARDpath<-paste0(dir,"IIDIPUS_Input_NoBuildingDat/HAZARDobjects/",namer)
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
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",i, " ",unique(miniDamSimplified$iso3)[1]," ", unique(miniDamSimplified$sdate)[1])) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input_NoBuildingDat/BDobjects/",namer)
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

# Re-do subnational impact data without double counting
redoSubNat <- function(dir, haz="EQ", extractedData=T, subnat_file= 'EQ_SubNational.xlsx'){
  # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if
  # no corresponding existing ODD object can be found, creates a new ODD object.

  #existingODDloc <- paste0(dir, "IIDIPUS_Input/ODDobjects/")
  existingODDloc <-"/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All/ODDobjects"
  existingODDfiles <- na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T))

  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)

  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% group_split()

  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  DispData <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/DispData_EQ_V2.Rdata')

  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (i in 1:10){

    # Subset displacement and disaster database objects
    miniDam<-SubNatDataByEvent[[i]]
    print(miniDam)
    miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3), sdate=miniDam$sdate[1],
                                    fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])

    #CHECK IF EXISTING ODD OBJECT AND UPDATE INSTEAD:
    filename_options <- paste0( miniDam$hazard[1], rep(gsub("-", "", seq(miniDam$sdate[1]-2, miniDam$sdate[1]+2, by='days')), each=length(miniDamSimplified$iso3)), miniDamSimplified$iso3)

    if( any(filename_options %in% substr(existingODDfiles, 1, 13))){
      matched_index <- match(filename_options, substr(existingODDfiles, 1, 13))
      file_match <- existingODDfiles[matched_index[which(!is.na(matched_index ))]]
      ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All/ODDobjects/', file_match))
      print('Checking hazards for the event')
      event_date <- as.Date(substr(file_match, 3, 10), "%Y%m%d")
      print(paste('Event Date:', event_date, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
      print('Using the data:')
      print(head(miniDam))
    } else {
      print(paste('No ODD object match found for ', miniDamSimplified))
      stop()
    }
    
    
    ODDy_with_impact <- updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1])
    
    if(is.null(ODDy_with_impact)) stop()
    ODDy <- ODDy_with_impact
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDamSimplified$iso3)[which(unique(miniDamSimplified$iso3) !='TOT')][1],
                  "_",i)
    
    iso3_ODDy <- unique(ODDy$ISO3C)
    iso3_ODDy <- iso3_ODDy[which(!is.na(iso3_ODDy))]
    
    iso3_impact <- unique(ODDy@impact$iso3)
    
    if(any(!iso3_ODDy %in% iso3_impact)){
      stop(paste('No impact data found for an affected country'))
    }
    
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_Input_All_ODD_Edited/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    maxdate<-miniDamSimplified$sdate[1]-5
    if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
    
    DispData_event <- filter(DispData, (iso3 %in% iso3_ODDy) &  (as.Date(sdate) >= maxdate) & (as.Date(sdate) <= mindate))
    if (NROW(DispData_event)==0){
      stop('Event not found in disp data')
    }
    lhazSDF <- GetUSGS_id(DispData_event$USGSid)
    cellsize <- 0.008333333333333
    if (lhazSDF@bbox[1,1] < (ODDy@bbox[1,1]-cellsize) | lhazSDF@bbox[1,2] > (ODDy@bbox[1,2]+cellsize) | lhazSDF@bbox[2,1] < (ODDy@bbox[2,1]-cellsize) | lhazSDF@bbox[2,2] > (ODDy@bbox[2,2]+cellsize)){
      stop('BBOX not contained')
    }
    
    # Save some RAM
    rm(ODDy,ODDy_with_impact)
  }

  return(path)

}

library(randomcoloR)
par(mfrow=c(2,2))

impact_by_type <- ODDy@impact %>% group_by(impact) %>% group_split()
for (i in 1:length(impact_by_type)){
  impact_filtered <- impact_by_type[[i]]
  ncols <- NROW(impact_filtered)
  palette <- distinctColorPalette(ncols)
  plot(ODDy@coords[ODDy@polygons[[impact_filtered$polygon[1]]]$indexes,], col='black', 
       xlim=ODDy@bbox[1,], ylim=ODDy@bbox[2,], main=impact_filtered$impact[1])
  if (ncols>1){
    for (ij in 2:ncols){
      points(ODDy@coords[ODDy@polygons[[impact_filtered$polygon[ij]]]$indexes,], col=palette[ij])
    }
  }
}


#download Philippines ODD objects from scratch
# getSubNatDataFiltered <- function(dir, haz="EQ",subnat_file= 'EQ_SubNational.xlsx', iso_filter='IRN'){
#   # Works through EQ_Subnational.xlsx and creates a new ODD object for each event in The Philippines
#   
#   #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
#   SubNatData <- readSubNatData(subnat_file)
#   
#   # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
#   SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% filter(iso_filter %in% iso3) %>% group_split()
#   
#   # Extract all building damage points
#   Damage<-ExtractBDfiles(dir = dir,haz = haz)
#   
#   # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
#   path<-data.frame()
#   for (i in 1:length(SubNatDataByEvent)){
#     
#     # Subset displacement and disaster database objects
#     miniDam<-SubNatDataByEvent[[i]]
#     print(miniDam)
#     miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3)[!unique(miniDam$iso3) %in% 'TOTAL'], sdate=miniDam$sdate[1], 
#                                     fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
#     
#     
#     # Set some confining dates for the rest of the data to be assigned to this event
#     maxdate<-miniDamSimplified$sdate[1]-5
#     if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
#     
#     # Match displacement and hazard data and extract hazard maps
#     # HazSDF includes SpatialPixelDataFrame object of hazmean & hazsd per date 
#     # (list of each, bbox-cropped to remove M < minmag)
#     #lhazSDF<-tryCatch(GetDisaster(miniDam,miniDACS),error=function(e) NULL)
#     #lhazSDF<-tryCatch(GetDisaster(miniDamSimplified),error=function(e) NULL)
#     lhazSDF <- GetDisaster(miniDamSimplified)
#     #lhazSDF <- GetDisaster(miniDamSimplified, bbox = c(114, -10, 120, -6.5))
#     
#     if(is.null(lhazSDF)) {
#       print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
#                    " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
#       next
#     }
#     
#     # Create the ODD object:
#     #ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified),error=function(e) NULL)
#     ODDy <- new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified)
#     if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
#     
#     ODDy_with_impact <- updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1])
#     ODDy <- ODDy_with_impact
#     # Create a unique hazard event name
#     namer<-paste0(ODDy@hazard,
#                   str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
#                   unique(miniDamSimplified$iso3)[1],
#                   "_",ODDy@eventid)
#     
#     # Save out objects to save on RAM
#     ODDpath<-paste0(dir,"IIDIPUS_Input_NoBuildingDat/ODDobjects/",namer)
#     saveRDS(ODDy,ODDpath)
#     
#     HAZARDpath<-paste0(dir,"IIDIPUS_Input_NoBuildingDat/HAZARDobjects/",namer)
#     saveRDS(lhazSDF,HAZARDpath)
#     rm(lhazSDF)
#     
#     #ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
#     
#     # Building damage subset
#     miniDam<-Damage%>%filter(iso3%in%unique(miniDamSimplified$iso3) & 
#                                sdate<mindate & sdate>maxdate)
#     
#     # Get building damage data and filter to matched hazard events
#     BDpath=NA_character_
#     if(nrow(miniDam)>0) {
#       # Make building damage object BD
#       BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
#       if(is.null(BDy)) {print(paste0("BD FAIL: ",i, " ",unique(miniDamSimplified$iso3)[1]," ", unique(miniDamSimplified$sdate)[1])) ;next}
#       BDpath <-paste0(dir,"IIDIPUS_Input_NoBuildingDat/BDobjects/",namer)
#       # Save it out!
#       saveRDS(BDy, BDpath)
#     }
#     
#     # Path to file and ID
#     path%<>%rbind(data.frame(ODDpath=ODDpath,
#                              BDpath=BDpath,
#                              eventid=ODDy@eventid))
#     # Save some RAM
#     rm(ODDy,BDy,miniDam)
#   }
#   
#   return(path)
#   
# }

# # check all hazards have been retrieved for the events
# checkAllHazards <- function(dir, haz="EQ",extractedData=T, subnat_file= 'EQ_SubNational.xlsx'){
#   # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if 
#   # no corresponding existing ODD object can be found, creates a new ODD object.
#   
#   #existingODDloc <- paste0(dir, "IIDIPUS_Input/ODDobjects/")
#   existingODDloc <-"/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All/ODDobjects"
#   existingODDfiles <- na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T)) 
#   
#   #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
#   SubNatData <- readSubNatData(subnat_file)
#   
#   # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
#   SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% group_split()
#   
#   # Extract all building damage points
#   Damage<-ExtractBDfiles(dir = dir,haz = haz)
#   
#   # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
#   path<-data.frame()
#   for (i in 165:length(SubNatDataByEvent)){
#     
#     # Subset displacement and disaster database objects
#     miniDam<-SubNatDataByEvent[[i]]
#     print(miniDam)
#     miniDamSimplified <- data.frame(iso3=unique(miniDam$iso3), sdate=miniDam$sdate[1], 
#                                     fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
#     
#     #CHECK IF EXISTING ODD OBJECT AND UPDATE INSTEAD:
#     filename_options <- paste0( miniDam$hazard[1], rep(gsub("-", "", seq(miniDam$sdate[1]-2, miniDam$sdate[1]+2, by='days')), each=length(miniDamSimplified$iso3)), miniDamSimplified$iso3)
#     
#     if( any(filename_options %in% substr(existingODDfiles, 1, 13))){
#       matched_index <- match(filename_options, substr(existingODDfiles, 1, 13))
#       file_match <- existingODDfiles[matched_index[which(!is.na(matched_index ))]]
#       ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All/ODDobjects/', file_match))
#       print('Checking hazards for the event the event')
#       event_date <- as.Date(substr(file_match, 3, 10), "%Y%m%d") 
#       print(paste('Event Date:', event_date, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
#       print('Using the data:')
#       print(head(miniDam))
#     } else {
#       print(paste('No ODD object match found for ', miniDamSimplified))
#       stop()
#     }
#     
#     # Set some confining dates for the rest of the data to be assigned to this event
#     maxdate<-miniDamSimplified$sdate[1]-5
#     if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
#     
#     if(!is.na(miniDam$GetDisaster_Args[1])){
#       next
#       lhazSDF<- eval(parse(text=paste0('GetDisaster(miniDamSimplified, ', miniDam$GetDisaster_Args[1], ')')))
#     } else {
#       bbox<-countriesbbox(unique(miniDamSimplified$iso3))
#       bbox[1] <- max(bbox[1]-5, -180)
#       bbox[2] <- max(bbox[2]-5, -90)
#       bbox[3] <- min(bbox[3]+5, 180)
#       bbox[4] <- min(bbox[4]+5, 90)
#       lhazSDF <- GetDisaster(miniDamSimplified, bbox = bbox)
#     }
#   
#     if(is.null(lhazSDF)) {
#       print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
#                    " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
#       stop()
#     }
#     
#     nhaz_ODD <- length(grep("hazMean",names(ODDy),value = T))
#     if(length(lhazSDF) != nhaz_ODD + 1){
#       print('Hazard number mismatch for event')
#       stop()
#     }
#     
#     # Save some RAM
#     rm(ODDy,miniDam)
#   }
#   
#   return(path)
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# cr <- colorRamp(c("green", "black"))
# plot(ODDy@coords[which(!is.na(ODDy$Population)),], col=rgb(cr(ODDy$Population[which(!is.na(ODDy$Population))]/max(ODDy$Population, na.rm=T)), max=255))
# plot(ODDy@coords[which(!is.na(ODDy$Population)),], col=rgb(cr(ODDy$Vs30[which(!is.na(ODDy$Population))]/max(ODDy$Vs30, na.rm=T)), max=255))
# plot(ODDy@coords[which(!is.na(ODDy$Population)),], col=rgb(cr(ODDy$EQFreq[which(!is.na(ODDy$Population))]/10), max=255))
# plotODDy(ODDy)
# plot(ODDy@coords[which(!is.na(ODDy$Population)),], col=rgb(cr(ODDy$GNIc[which(!is.na(ODDy$Population))]/max(ODDy$GNIc, na.rm=T)), max=255))
# 
# points(ODDy@coords[which(is.na(ODDy@data$EQFreq)),], col='red')
# points(ODDy@coords[which(!is.na(ODDy$Population)),], col='blue')
# 
# 
# 
# library(randomcoloR)
# ncols <- length(ODDy@polygons)
# palette <- distinctColorPalette(ncols)
# plot(ODDy@coords[ODDy@polygons[[45]]$indexes,], col='black')
# 
# for (ij in c(9, 16, 30, 43, 66, 56)){
#   Sys.sleep(2)
#   points(ODDy@coords[ODDy@polygons[[ij]]$indexes,], col=palette[ij])
# }
# plot(ODDy@coords[which(!is.na(ODDy$ISO3C)),])


# -----------------------------------------------------------------------------------------

# # # Loop through and check EQ Disagg Data for errors/issues
# # 
# subnat_file='IIDIPUS_Input/EQ_SubNational.xlsx'
# #SubNatData <- read.xlsx(paste0(dir, subnat_file), colNames = TRUE , na.strings = c("","NA"))
# SubNatData <- readSubNatData(subnat_file)
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
# #SubNatData <- read.xlsx(subnat_file, colNames = TRUE , na.strings = c("","NA"))
# SubNatData <- readSubNatData(subnat_file)
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

