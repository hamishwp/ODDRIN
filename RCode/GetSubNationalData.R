################
### ODDpolys ###
################



library(openxlsx)
library(tidyxl)
library(rgeos)
library(sfheaders)
library(ecochange)

# cleanSubNatData <- function(SubNatData){
#   # Clean data from xlsx file by converting data to the appropriate formats/data types
#   SubNatData$source_date <- openxlsx::convertToDate(SubNatData$source_date)
#   SubNatData$sdate <- openxlsx::convertToDate(SubNatData$sdate)
#   SubNatData$fdate <- openxlsx::convertToDate(SubNatData$fdate)
#   SubNatData$mortality <- as.integer(SubNatData$mortality)
#   SubNatData$displacement <- as.integer(SubNatData$displacement)
#   SubNatData$buildDam <- as.integer(SubNatData$buildDam)
#   SubNatData$buildDest <- as.integer(SubNatData$buildDest)
#   return(SubNatData)
# }

folder_write <- 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/'

readSubNatData <- function(subnat_file){
  # Clean data from xlsx file by converting data to the appropriate formats/data types
  
  SubNatData <-  read.xlsx('/home/manderso/Downloads/EQ_SubNational.xlsx', colNames = TRUE , na.strings = c("","NA"))
  x <-xlsx_cells('/home/manderso/Downloads/EQ_SubNational.xlsx')
  formats <- xlsx_formats('/home/manderso/Downloads/EQ_SubNational.xlsx')
  
  # SubNatData <-  read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  # x <-xlsx_cells(paste0(dir, 'IIDIPUS_Input/', subnat_file))
  # formats <- xlsx_formats(paste0(dir, 'IIDIPUS_Input/', subnat_file))

  red_cells <- x %>% filter(local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == "FFFF0000")) %>% dplyr::select(row, col)
  red_cells$row <- red_cells$row - 1
  red_cells <- red_cells[-which(red_cells$col > NCOL(SubNatData)), ]
  SubNatData[as.matrix(red_cells)] <- NA
  
  pink_cells <- x %>% filter(local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == "FFFF00FF")) %>% dplyr::select(row, col)
  SubNatData$buildDamInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='buildDam_exlusion_reason'))]-1)
  SubNatData$buildDestInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='buildDest_exlusion_reason'))]-1)
  SubNatData$buildDamDestInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='buildDamDest_exlusion_reason'))]-1)
  SubNatData$mortalityInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='mortality_exlusion_reason'))]-1)
  SubNatData$displacementInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='displacement_exlusion_reason'))]-1)
  
  SubNatData <- SubNatData[!is.na(SubNatData$iso3), ]
  
  SubNatData$source_date <- openxlsx::convertToDate(SubNatData$source_date)
  SubNatData$sdate <- openxlsx::convertToDate(SubNatData$sdate)
  SubNatData$fdate <- openxlsx::convertToDate(SubNatData$fdate)
  SubNatData$GetDisaster_sdate <- openxlsx::convertToDate(SubNatData$GetDisaster_sdate)
  SubNatData$GetDisaster_fdate <- openxlsx::convertToDate(SubNatData$GetDisaster_fdate)
  SubNatData$mortality <- as.integer(SubNatData$mortality)
  SubNatData$displacement <- as.integer(SubNatData$displacement)
  SubNatData$buildDam <- as.integer(SubNatData$buildDam)
  SubNatData$buildDest <- as.integer(SubNatData$buildDest)
  SubNatData$buildDamDest <- as.integer(SubNatData$buildDamDest)
  SubNatData$EventID <- suppressWarnings(as.integer(SubNatData$EventID)) #will set non-integers to NA but suppress the warning for this
  SubNatData$Region <- ifelse(trimws(SubNatData$Region) == "", NA, SubNatData$Region)
  SubNatData$Subregion <- ifelse(trimws(SubNatData$Subregion) == "", NA, SubNatData$Subregion)
  
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
  if (is.na(GADM_array[1]) & is.na(GADM_array[2])){ #national data, can use iso3 code already in ODD object
    return(list(polygon_name = polygon_name, sf_polygon = NULL))
    #GADM_level <- 0
    #GADM_array <- GADM_array[3]
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
                       observed = numeric(), qualifier = character(), build_type = character())
  
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
        if (impact_type == 'buildDam' || impact_type == 'buildDest' || impact_type == 'buildDamDest'){
          if((length(unique(pull(SubNatEvent_by_polygon[[i]][,impact_type])))==1) & (length(unique(pull(SubNatEvent_by_polygon[[i]][,'building_type'])))==1)){ #all matching, so just take first
            sources_selected %<>% append(SubNatEvent_by_polygon[[i]]$source_id[1])
            next
          }
          #if different building types then take all
          if (length(unique(SubNatEvent_by_polygon[[i]] %>% pull(impact_type))) <= length(unique(SubNatEvent_by_polygon[[i]]$building_type))){
            source_ids <- SubNatEvent_by_polygon[[i]]$source_id
            if (any(duplicated(SubNatEvent_by_polygon[[i]][,c(impact_type,'building_type')]))){
              source_ids <- source_ids[-which(duplicated(SubNatEvent_by_polygon[[i]][,c(impact_type,'building_type')]))]
            }
            sources_selected %<>% append(source_ids)
            next
          } 
        }
        if(length(unique(pull(SubNatEvent_by_polygon[[i]][,impact_type])))==1){ #all matching, so just take first
          sources_selected %<>% append(SubNatEvent_by_polygon[[i]]$source_id[1])
          next
        }
        #conflicting sources:
        cat(paste('Please select between the following sources, by typing the id of the chosen source and pressing return.\n',
                  'If you would like to select more than one source, please separate each id with a comma and we will use the mean \n')) 
        if(impact_type == 'buildDam' || impact_type == 'buildDest' || impact_type == 'buildDamDest'){
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
                                                                         inferred=!!sym(paste0(impact_type, 'Inferred')),
                                                                         build_type=ifelse(rep(impact_type, length(sources_selected)) %in% c('buildDam', 'buildDest', 'buildDamDest'), building_type, NA),
                                                                         polygon=polygon_id))
    }
  }
  return(list(impact=impact, polygons_list=polygons_list))
}

updateODDSubNat <- function(dir, ODDy, event_sdate, event_fdate, event_id, subnat_file='EQ_SubNational.xlsx'){
  # Update an existing ODD object (ODDy) using subnational data:
  #   - Edit the 'impact' slot to include subnational data from the provided spreadsheet (default is EQ_subnational.xlsx)
  #   - Edit the 'polygons' slot to identify the pixels belonging to each polygon for which we have impact data
  
  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)
  SubNatData <- SubNatData[which(!is.na(SubNatData$event_name)), ]
  
  SubNatData_match <- SubNatData %>% group_by(event_name, sdate) %>% filter(EventID==event_id) %>% group_split()
  # Use event date and iso3 to identify the event in subnat data that correspond to the ODD object
  # (note that all rows for the same event in SubNatData should have the same name and sdate!)
  # SubNatData_match <- SubNatData %>% group_by(event_name, sdate) %>% filter(
  #   length(intersect(trimws(unlist(strsplit(source_listed_iso3, ";"))),unique(ODDy@data$ISO3C)))>0,
  #   sdate > (event_sdate - 1) &  fdate < (event_fdate + 1) #LOOSEEND: HAZDATES FOR EXISTING ODD OBJECTS ARE DODGEY
  # ) %>% group_split()
  
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
  
  iso3_unique <- unique(ODDy$ISO3C)[which(!is.na(unique(ODDy$ISO3C)))]
  if (length(iso3_unique)==1){
    #replace 'TOTAL' with name of the country
    SubNatEvent$iso3[which(SubNatEvent$iso3=='TOT')] = iso3_unique
    SubNatEvent$country[which(SubNatEvent$country=='TOTAL')] = countrycode(iso3_unique, origin='iso3c', destination='country.name')
  }

  SubNatImpact <- getSubNatImpact(SubNatEvent, subnational=TRUE) #populate the 'impact' slot in ODDy using the national and subnational data
  
  #remove overlapping polygons using SubNatImpact$impact and SubNatImpact$polygons_list
  #find the polygons and complements we want to use for each impact type and add to ODDy@impact
  
  ODDy <- addODDPolygons(ODDy, SubNatImpact$polygons_list) #populate the 'polygons' slot in ODDy using the regions for which we have impact data
  ODDy <- addODDImpact(ODDy, SubNatImpact$impact)
  
  #additional_poly_check(ODDy, event_id, print_to_xl=T)
  ODDy <- reweight_pixels(ODDy)
  
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
  
  if (length(indies)==0){return(integer())}
  
  check_intersection <- function(j){
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
      return(T)
    } else return(F)
  }
  
  inside_poly <- unlist(mclapply(indies, check_intersection, mc.cores=2))
  
  if (length(which(inside_poly))==0){return(integer())}
  
  return(indies[which(inside_poly)])
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
        if (polygons_list[[i]]$polygon_name=='Kosovo'){ #kosovo doesn't seem to work with countrycode
          iso3='KOS'
        } else {
          iso3 <- countrycode(polygons_list[[i]]$polygon_name, origin='country.name', destination='iso3c')
        }
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

reweight_pixels <- function(ODDy){
  # If a pixel has a total weight less than 1 across all polygons (for a specific administrative level), 
  # then either part of the pixel lies in the water or in a region for which we do not have impact data. 
  # We assume the former, and scale the weights so that they sum to one. 
  # Is there a way of knowing if it is the latter, and not rescaling? 
  pixels_of_interest <- which(!is.na(ODDy$ISO3C))
  for (i in pixels_of_interest){
    weight_sum = c(0, 0, 0)
    polygon_matches = list(c(), c(), c())
    for (p in 1:length(ODDy@polygons)){
      poly <- ODDy@polygons[[p]]
      match_i <- which(poly$indexes == i)
      if (length(match_i) > 0){
        gadm_level <- str_count(poly$name, ',')
        polygon_matches[[gadm_level+1]] %<>% rbind(c(p, match_i))
        weight_sum[gadm_level+1] <- weight_sum[gadm_level+1] + poly$weights[match_i]
      }
    }
    for (g in 2:3){
      if (is.null(polygon_matches[[g]])) next
      if(round(weight_sum[g],3) > 1){
        file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
        writeLines(paste("Region:", ODDy@polygons[[polygon_matches[[g]][1,1]]]$name, "Event Date:", ODDy@hazdates[1], '. Sum of pixel weights for a polygon is larger than 1'), file_conn)
        close(file_conn)
      } 
      if (round(weight_sum[g],6) < 1 & weight_sum[g] > 0){
        for (pix_matched_i in 1:NROW(polygon_matches[[g]])){
          ODDy@polygons[[polygon_matches[[g]][pix_matched_i,1]]]$weights[polygon_matches[[g]][pix_matched_i,2]] <- ODDy@polygons[[polygon_matches[[g]][pix_matched_i,1]]]$weights[polygon_matches[[g]][pix_matched_i,2]] / weight_sum[g]
        }
      }
    }
  }
  return(ODDy)
}

remove_overlaps <- function(polygons_indexes, polygons_list, coords){
  #find which polygons have intersecting indexes but not intersecting sf polygons
  #this may be the case for pixels lying on the border which have been allocated to two polygons
  #in this case, allocate the pixel to the polygon in which its center lies
  for (i in 1:length(polygons_indexes)){
    polygons_indexes[[i]]$weights <- rep(1, length(polygons_indexes[[i]]$indexes))
  }
  
  if (length(polygons_indexes) == 1) return(polygons_indexes)
  
  lon_cellsize <- ODDy@grid@cellsize[1]
  lat_cellsize <- ODDy@grid@cellsize[2]
  
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
        for (index in indexes_intersect){
          #to allocate a weight to the pixel based on how much of the polygon it contains
          pixel_sf <- sfheaders::sf_polygon(data.frame(longitude = c(coords[index,1]-lon_cellsize/2, 
                                                                     coords[index,1]-lon_cellsize/2,
                                                                     coords[index,1]+lon_cellsize/2, 
                                                                     coords[index,1]+lon_cellsize/2, 
                                                                     coords[index,1]-lon_cellsize/2), 
                                                       latitude = c(coords[index,2]-lat_cellsize/2, 
                                                                    coords[index,2]+lat_cellsize/2,
                                                                    coords[index,2]+lat_cellsize/2, 
                                                                    coords[index,2]-lat_cellsize/2,
                                                                    coords[index,2]-lat_cellsize/2)), x='longitude',y='latitude')
          sf::st_crs(pixel_sf) = ODDy@proj4string
          pixel_sp <- as(pixel_sf, 'Spatial')
          intersection_i <- gIntersection(pixel_sp,polygons_list[[i]]$sf_polygon)
          intersection_j <- gIntersection(pixel_sp,polygons_list[[j]]$sf_polygon)
          polygons_indexes[[i]]$weights[which(polygons_indexes[[i]]$indexes==index)] <- ifelse(is.null(intersection_i), 0, gArea(intersection_i)) / (lat_cellsize * lon_cellsize)
          polygons_indexes[[j]]$weights[which(polygons_indexes[[j]]$indexes==index)] <- ifelse(is.null(intersection_j), 0, gArea(intersection_j)) / (lat_cellsize * lon_cellsize)
        }
        
        #to simply allocate to the closest polygon:
        
        # allocated_poly_i <- c()
        # allocated_poly_j <- c()
        # for (index in indexes_intersect){
        #   spatial_points <- SpatialPoints(coords, proj4string = CRS(proj4string(ODDy)))
        #   if (gDistance(spatial_points[index], polygons_list[[i]]$sf_polygon) < 
        #       gDistance(spatial_points[index], polygons_list[[j]]$sf_polygon)){
        #     allocated_poly_i %<>% append(index)
        #   } else {allocated_poly_j %<>% append(index)}
        # }
        # polygons_indexes[[i]]$indexes <- polygons_indexes[[i]]$indexes[! polygons_indexes[[i]]$indexes %in% allocated_poly_j]
        # polygons_indexes[[j]]$indexes <- polygons_indexes[[j]]$indexes[! polygons_indexes[[j]]$indexes %in% allocated_poly_i]
      
        }
    }
  }
  
  return(polygons_indexes)
  
}

remove_double_counting <- function(impact, ODDy){
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
  return(impact_updated)
}

addODDImpact <- function(ODDy, impact){
  #impact %<>% remove_double_counting(ODDy)
  ODDy@impact <- impact
  
  #if polygon is 'TOTAL' but only one country is exposed, can rename the 'TOTAL' polygon with the country
  unique_iso3 <- unique(ODDy$ISO3C)[!is.na(unique(ODDy$ISO3C))]
  if (length(unique_iso3)==1){
    i_TOTAL <- which(sapply(ODDy@polygons, function(poly){ifelse(poly$name=='TOTAL', T, F)}))
    i_iso3 <- which(sapply(ODDy@polygons, function(poly){ifelse(grepl(",",poly$name), F, T)}))
    i_iso3 <- setdiff(i_iso3, i_TOTAL)
    ODDy@impact$iso3[which(ODDy@impact$polygon==i_TOTAL)] <- unique_iso3
    ODDy@impact$polygon[which(ODDy@impact$polygon==i_TOTAL)] <- i_iso3
  }
  
  #remove obsolete polygons
  polys_no_indexes <-  which(unlist(lapply(ODDy@polygons, function(x){length(x$indexes)==0})))
  impact_rows_remove <- which(ODDy@impact$polygon %in% polys_no_indexes)
  
  if (any(ODDy@impact$observed[impact_rows_remove] > 0)){
    missing_polys <- which(ODDy@impact$observed[impact_rows_remove] > 0)
    missing_poly_names <- unique(unlist(lapply(ODDy@polygons[ODDy@impact$polygon[impact_rows_remove[missing_polys]]], function(x){x$name})))
    file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
    writeLines(paste("Country:", unique_iso3[1], "Event Date:", ODDy@hazdates[1], ', region', missing_poly_names, 'with impact greater than 0 not included in modelled region.'), file_conn)
    close(file_conn)
  }
  
  if (length(impact_rows_remove) > 0){
    ODDy@impact <- ODDy@impact[-impact_rows_remove,]
  }
  
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
    miniDamSimplified <- data.frame(iso3=unique(miniDam$source_listed_iso3), sdate=miniDam$sdate[1], 
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
    
    ODDy_with_impact <- updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1], i)
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

#Re-do subnational impact data without double counting
redoSubNat <- function(dir, haz="EQ", extractedData=T, subnat_file= 'EQ_SubNational.xlsx'){
  # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if
  # no corresponding existing ODD object can be found, creates a new ODD object.

  #existingODDloc <- paste0(dir, "IIDIPUS_Input/ODDobjects/")
  existingODDloc <-"/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All_Edited/ODDobjects/"
  existingODDfiles <- na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T))
  
  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)

  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(sdate, event_name) %>% group_split() 

  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)

  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (i in c(1:5, 7:94, 96:162)){

    # Subset displacement and disaster database objects
    miniDam<-SubNatDataByEvent[[i]]
    miniDamSimplified <- data.frame(iso3=unique(trimws(unlist(strsplit(miniDam$source_listed_iso3, ";")))), 
                                    sdate=miniDam$sdate[1], fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
    
    #CHECK IF EXISTING ODD OBJECT AND UPDATE INSTEAD:
    filename_options <- paste0( miniDam$hazard[1], rep(gsub("-", "", seq(miniDam$sdate[1]-1, miniDam$fdate[1]+1, by='days')), each=length(miniDamSimplified$iso3)), miniDamSimplified$iso3)
    
    if( any(filename_options %in% substr(existingODDfiles, 1, 13))){
      matched_index <- match(filename_options, substr(existingODDfiles, 1, 13))
      file_match <- existingODDfiles[matched_index[which(!is.na(matched_index ))]]
      ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All_Edited/ODDobjects/', file_match))
      event_date <- as.Date(substr(file_match, 3, 10), "%Y%m%d")
    } else {
      print(paste('No ODD object match found for ', miniDamSimplified))
      stop()
    }
    maxdate<-miniDamSimplified$sdate[1]-5
    if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
    
    miniDam<-Damage%>%filter(iso3%in%unique(miniDamSimplified$iso3) & 
                               sdate<mindate & sdate>maxdate)
    
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDamSimplified$iso3)[which(unique(miniDamSimplified$iso3) !='TOT')][1],
                  "_",i)
    
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",i, " ",unique(miniDamSimplified$iso3)[1]," ", unique(miniDamSimplified$sdate)[1])) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input_All_Edited/BDobjects/",namer)
      # Save it out!
      saveRDS(BDy, BDpath)
    }
  }

  return(path)

}

# 
# # Re-do subnational impact data without double counting
# redoSubNat <- function(dir, haz="EQ", extractedData=T, subnat_file= 'EQ_SubNational.xlsx'){
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
#   DispData <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/DispData_EQ_V2.Rdata')
# 
#   # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
#   path<-data.frame()
#   for (i in 1:10){
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
#       print('Checking hazards for the event')
#       event_date <- as.Date(substr(file_match, 3, 10), "%Y%m%d")
#       print(paste('Event Date:', event_date, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
#       print('Using the data:')
#       print(head(miniDam))
#     } else {
#       print(paste('No ODD object match found for ', miniDamSimplified))
#       stop()
#     }
#     
#     
#     ODDy_with_impact <- updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1])
#     
#     if(is.null(ODDy_with_impact)) stop()
#     ODDy <- ODDy_with_impact
#     
#     # Create a unique hazard event name
#     namer<-paste0(ODDy@hazard,
#                   str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
#                   unique(miniDamSimplified$iso3)[which(unique(miniDamSimplified$iso3) !='TOT')][1],
#                   "_",i)
#     
#     iso3_ODDy <- unique(ODDy$ISO3C)
#     iso3_ODDy <- iso3_ODDy[which(!is.na(iso3_ODDy))]
#     
#     iso3_impact <- unique(ODDy@impact$iso3)
#     
#     if(any(!iso3_ODDy %in% iso3_impact)){
#       stop(paste('No impact data found for an affected country'))
#     }
#     
#     # Save out objects to save on RAM
#     ODDpath<-paste0(dir,"IIDIPUS_Input_All_ODD_Edited/ODDobjects/",namer)
#     saveRDS(ODDy,ODDpath)
#     
#     maxdate<-miniDamSimplified$sdate[1]-5
#     if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
#     
#     DispData_event <- filter(DispData, (iso3 %in% iso3_ODDy) &  (as.Date(sdate) >= maxdate) & (as.Date(sdate) <= mindate))
#     if (NROW(DispData_event)==0){
#       stop('Event not found in disp data')
#     }
#     lhazSDF <- GetUSGS_id(DispData_event$USGSid)
#     cellsize <- 0.008333333333333
#     if (lhazSDF@bbox[1,1] < (ODDy@bbox[1,1]-cellsize) | lhazSDF@bbox[1,2] > (ODDy@bbox[1,2]+cellsize) | lhazSDF@bbox[2,1] < (ODDy@bbox[2,1]-cellsize) | lhazSDF@bbox[2,2] > (ODDy@bbox[2,2]+cellsize)){
#       stop('BBOX not contained')
#     }
#     
#     # Save some RAM
#     rm(ODDy,ODDy_with_impact)
#   }
# 
#   return(path)
# 
# }
# 
# # Rename sub-national impact data
# redoSubNat <- function(dir, haz="EQ", extractedData=T, subnat_file= 'EQ_SubNational.xlsx'){
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
#   
#   for (file in existingODDfiles[49]){
#     ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All/ODDobjects/', file))
#     
#     i <-  which(sapply(SubNatDataByEvent, function(x) return(any((ODDy$ISO3C %>% unique()) %in% x$iso3) & (ODDy@hazdates[1] >= x$sdate[1]-1) & (ODDy@hazdates[1] <= x$fdate[1]+1))))
#     if (length(i) != 1){
#       stop()
#     }
#    filename_new <-  gsub("\\_.*",paste0('_',i),file)
#    saveRDS(ODDy, paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Renamed/ODDobjects/', filename_new)) 
#    if(file.exists(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All/BDobjects/', file))){
#      BDy <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All/BDobjects/', file))
#      saveRDS(BDy, paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Renamed/BDobjects/', filename_new)) 
#    }
#   }
# }


# Re-do subnational impact data with double counting
# Check impact seems to be MAR
# redoSubNat <- function(dir, haz="EQ", extractedData=T, subnat_file= 'EQ_SubNational.xlsx'){
#   # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if
#   # no corresponding existing ODD object can be found, creates a new ODD object.
#   
#   #existingODDloc <- paste0(dir, "IIDIPUS_Input/ODDobjects/")
#   existingODDloc <-"/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NMAR/ODDobjects"
#   existingODDfiles <- na.omit(list.files(path=existingODDloc,pattern=Model$haz,recursive = T, ignore.case = T))
#   
#   #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
#   SubNatData <- readSubNatData(subnat_file)
#   
#   # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
#   SubNatDataByEvent <- SubNatData %>% group_by(sdate, event_name) %>% group_split()
#   
#   # Extract all building damage points
#   #Damage<-ExtractBDfiles(dir = dir,haz = haz)
#   
#   DispData <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/DispData_EQ_V2.Rdata')
#   
#   # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
#   path<-data.frame()
#   for (i in 1:10){
#     
#     # Subset displacement and disaster database objects
#     miniDam<-SubNatDataByEvent[[i]]
#     print(miniDam)
#     miniDamSimplified <- data.frame(iso3=unique(trimws(unlist(strsplit(miniDam$source_listed_iso3, ";")))), 
#                                     sdate=miniDam$sdate[1], fdate=miniDam$fdate[1], eventid=i, hazard=miniDam$hazard[1])
#     
#     #CHECK IF EXISTING ODD OBJECT AND UPDATE INSTEAD:
#     filename_options <- paste0( miniDam$hazard[1], rep(gsub("-", "", seq(miniDam$sdate[1]-1, miniDam$sdate[1]+1, by='days')), each=length(miniDamSimplified$iso3)), miniDamSimplified$iso3)
#     
#     if( any(filename_options %in% substr(existingODDfiles, 1, 13))){
#       matched_index <- match(filename_options, substr(existingODDfiles, 1, 13))
#       file_match <- existingODDfiles[matched_index[which(!is.na(matched_index ))]]
#       ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NMAR/ODDobjects/', file_match))
#       print('Checking hazards for the event')
#       event_date <- as.Date(substr(file_match, 3, 10), "%Y%m%d")
#       print(paste('Event Date:', event_date, 'Countries:', paste(unique(ODDy@data$ISO3C[!is.na(ODDy@data$ISO3C)], na.rm=TRUE), collapse=" ")))
#       print('Using the data:')
#       print(head(miniDam))
#     } else {
#       print(paste('No ODD object match found for ', miniDamSimplified))
#       stop()
#       
#       maxdate<-miniDamSimplified$sdate[1]-5
#       if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
#       
#       bbox<-countriesbbox(unique(miniDamSimplified$iso3))
#       bbox[1] <- max(bbox[1]-5, -180)
#       bbox[2] <- max(bbox[2]-5, -90)
#       bbox[3] <- min(bbox[3]+5, 180)
#       bbox[4] <- min(bbox[4]+5, 90)
#       # lhazSDF <- GetDisaster(miniDamSimplified, bbox = c(160, -11.8, 163, -9.7))
#       lhazSDF<-GetDisaster(miniDamSimplified,bbox=bbox, EQparams=list(I0=4, minmag=4.5)) #, bbox=c(140,-9,147.5,-3)) #, EQparams=list(I0=4.5,minmag=4.5))
#       
#       if(is.null(lhazSDF)) {
#         print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
#                      " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
#         next
#       }
#       
#       # Create the ODD object:
#       ODDy<-new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified)
#       if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
#     
#     }
#     
#     
#     ODDy_with_impact <- updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1])
#     
#     if(is.null(ODDy_with_impact)) stop()
#     ODDy <- ODDy_with_impact
#     
#     # Create a unique hazard event name
#     namer<-paste0(ODDy@hazard,
#                   str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
#                   unique(miniDamSimplified$iso3)[which(unique(miniDamSimplified$iso3) !='TOT')][1],
#                   "_",i)
#     
#     iso3_ODDy <- unique(ODDy$ISO3C)
#     iso3_ODDy <- iso3_ODDy[which(!is.na(iso3_ODDy))]
#     
#     iso3_impact <- unique(ODDy@impact$iso3)
#     
#     if(any(!iso3_ODDy %in% iso3_impact)){
#       stop(paste('No impact data found for an affected country'))
#     }
#     
#     # Save out objects to save on RAM
#     ODDpath<-paste0(dir,"IIDIPUS_Input_All_ODD_Edited/ODDobjects/",namer)
#     saveRDS(ODDy,ODDpath)
#     
#     maxdate<-miniDamSimplified$sdate[1]-5
#     if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
#     
#     DispData_event <- filter(DispData, (iso3 %in% iso3_ODDy) &  (as.Date(sdate) >= maxdate) & (as.Date(sdate) <= mindate))
#     if (NROW(DispData_event)==0){
#       stop('Event not found in disp data')
#     }
#     lhazSDF <- GetUSGS_id(DispData_event$USGSid[1])
#     cellsize <- 0.008333333333333
#     if (lhazSDF@bbox[1,1] < (ODDy@bbox[1,1]-cellsize) | lhazSDF@bbox[1,2] > (ODDy@bbox[1,2]+cellsize) | lhazSDF@bbox[2,1] < (ODDy@bbox[2,1]-cellsize) | lhazSDF@bbox[2,2] > (ODDy@bbox[2,2]+cellsize)){
#       stop('BBOX not contained')
#     }
#     
#     plot_impact(ODDy)
#     
#     # Save some RAM
#     rm(ODDy,ODDy_with_impact)
#   }
#   
#   return(path)
#   
# }

library(openxlsx)
additional_poly_check <- function(ODDy, i, print_to_xl=F){
  
  noteworthy <- F
  
  impacts_split <- ODDy@impact %>% group_by(impact) %>% group_split()
  
  if(print_to_xl){wb <- createWorkbook()}
  
  for (j in 1:length(impacts_split)){
    print(j)
    if(print_to_xl){sheet <- addWorksheet(wb, impacts_split[[j]]$impact[1])}
    
    row = 1
    
    print(impacts_split[[j]]$impact[1])
    #check that all regions are on the same GADM admin level:
    
    polygons_impact <- ODDy@polygons[impacts_split[[j]]$polygon]
    polygon_names <- sapply(polygons_impact, function(x) x$name)
    # -1 = TOTAL
    # 0 = NATIONAL
    # 1 = SUBNATIONAL
    # 2 = SUB-SUB NATIONAL
    gadm_levels <- sapply(polygon_names, function(x) ifelse(grepl(",", x), length(gregexpr(",", x)[[1]]), ifelse(x=='TOTAL', -1, 0)))
    
    if (length(unique(gadm_levels)) > 1){
      noteworthy = T
      if (print_to_xl){
        writeData(wb, sheet, 'More than one admin level', startCol = 1, startRow = row)
        row = row + 1
      }
      print('More than one admin level. Go to data and manually delete if necessary.')
    }
    
    for (gadm_level in sort(unique(gadm_levels))){
      writeData(wb, sheet, paste('GADM Level ', gadm_level), startCol = 1, startRow = row)
      row = row + 1
      if (gadm_level == -1){
        writeData(wb, sheet, 'TOTAL', startCol = 1, startRow = row)
        row = row + 1
      }
      if (gadm_level == 0) { #handle countries
        iso3_incl <- c()
        for (k in 1:NROW(impacts_split[[j]])){
          if (gadm_levels[k] == gadm_level){
            iso3_incl %<>% append(impacts_split[[j]]$iso3[k])
          }
        }
        non_na_iso3 <- unique(ODDy$ISO3C)[which(!is.na(unique(ODDy$ISO3C)))]
        iso3_missing <- non_na_iso3[which(!(non_na_iso3 %in% iso3_incl))]
        for (iso3_miss in iso3_missing){
          writeData(wb, sheet, iso3_miss, startCol = 1, startRow = row)
          row = row + 1
          noteworthy=T
        }
        next
      }
      
      #check if there are any pixels not covered by the polygons:
      pixels_unallocated <- which(!is.na(ODDy$ISO3C))
      for (k in 1:NROW(impacts_split[[j]])){
        if (gadm_levels[k] == gadm_level)
          pixels_unallocated <- setdiff(pixels_unallocated, ODDy@polygons[[impacts_split[[j]]$polygon[k]]]$indexes)
      }
      if (length(pixels_unallocated) == 0){
        next
      }
      
      for (iso3 in unique(impacts_split[[j]]$iso3)){
        if (iso3 == 'TOT'){
          next
        }
        if (gadm_level == 2){
          regions_gadm_level <- getGADM(level=gadm_level, country=iso3)
          regions_gadm_level <- sapply(seq_along(regions_gadm_level), function(i) paste0(regions_gadm_level[i], ', ', names(regions_gadm_level)[i]))
          regions_gadm_level <- setdiff(regions_gadm_level, sub(",[^,]*$", "", polygon_names))
        } else if (gadm_level == 1){
          regions_gadm_level <- getGADM(level=gadm_level, country=iso3)
          regions_gadm_level <- setdiff(regions_gadm_level, sub(",[^,]*$", "", polygon_names))
        } 
        
        for (r in regions_gadm_level){
          if (gadm_level==2){
            r <- strsplit(r,',')[[1]]
            r[2] <- substring(r[2], 2)
          }
          r_poly <- getGADM(r, level = gadm_level, country=iso3)
          overlap <- (max(r_poly@bbox[1,]) >= ODDy@bbox[1] &
                        min(r_poly@bbox[1,]) <= ODDy@bbox[3] &
                        max(r_poly@bbox[2,]) >= ODDy@bbox[2] &
                        min(r_poly@bbox[2,]) <= ODDy@bbox[4])
          if (!overlap) next
          
          if (length(pixels_unallocated)==0) next
          
          spdf_pixels_unallocated <- SpatialPointsDataFrame(coords=ODDy@coords[pixels_unallocated,, drop=F],
                                                            data=ODDy@data[pixels_unallocated,1:2, drop=F], #this is arbitrary, function just seems to require data
                                                            proj4string=r_poly@proj4string)
          
          contained <- gContains(r_poly, spdf_pixels_unallocated, byid=T)
          if (any(contained)){
            print(paste('Polygon', r, 'contains an unallocated pixel'))
            pixels_unallocated <- pixels_unallocated[-which(contained)]
            if(print_to_xl){ 
              r <- rev(r)
              writeData(wb, sheet, iso3, startCol = 1, startRow = row)
              for (ll in 1:length(r)){
                writeData(wb, sheet, r[ll], startCol = ll+1, startRow = row)
              }
              row = row + 1
              noteworthy=T
            }
          }
        }
      }
    }
  }
  if (print_to_xl & noteworthy){saveWorkbook(wb, paste0(folder_write, "Missing Regions/Event_",i,".xlsx"))}
}

createODD <- function(dir, haz="EQ", subnat_file= 'EQ_SubNational.xlsx'){
  # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if
  # no corresponding existing ODD object can be found, creates a new ODD object.
  
  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(EventID) %>% group_split()
  
  DispData <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/DispData_EQ_V2.Rdata')
  
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (i in c(69, 169, 170)){
    # if (!i %in% c(8, 13, 25, 28, 31, 39, 47, 52, 54, 59, 65, 67, 72, 74, 79, 80, 81, 85, 89, 
    #               109, 111, 116, 118, 123, 124, 125, 128, 131, 132, 133, 135, 137, 144, 145,
    #               151, 157, 164, 166)){
    #   next
    # }
    if (i %in% c(1, 126)){
      next
    }
    print(i)
    #Needs redoing: 69, 169, 170
    #Needs redoing: c(20, 52, 68, 70, 91, 92, 137, 139, 140, 144, 162, 170)
    #Needs redoing: c(20, 52, 67, 68, 70, 75, 76, 81, 82, 84, 91:122, 125, 130, 139, 140, 144, 162, 170)
    #failed: 67, 68, 137, 169 
    ###89 (lhazsdf not found), 78, 68, 69, 45, 35, 21, 16, 131,
    #lhazsdf_dispdat not found: 8, 9, 76, 84
    #no dispdat: 3,4, 30, 34, 45, 55, 71, 77, 82, 90, 91, 92, 129, 136, 164, 165, 166, 167, 168, 169, 170
    #ODD bbox not contained by disp data bbox: 7, 48, 52, 66, 67, 86, 104, 107, 108, 109, 117 (x2), 133, 145
    #Impact larger than 0 outside of I > I0: 92, 97, 114
    
    
    #building counts not added: 24, 31, 36, 55, 56, 57, 81, 99, 101, 114
    #missing bing building count: peru 2013-02-09, peru 2016-04-16, 2018-09-07
    #check weights:
    #   - Philippines 2013-10-15 Iloilo, Jordan
    #BD failed: 2019-08-08
    
    # Subset displacement and disaster database objects
    SubNatDataByEvent[[i]] <- SubNatDataByEvent[[i]][rowSums(is.na(SubNatDataByEvent[[i]][,-1])) != ncol(SubNatDataByEvent[[i]])-1, ]
    miniDam<-SubNatDataByEvent[[i]]
    
    miniDamSimplified <- data.frame(iso3=unique(trimws(unlist(strsplit(miniDam$source_listed_iso3, ";")))), 
                                    sdate=data.table::fifelse(is.na(miniDam$GetDisaster_sdate[1]), miniDam$sdate[1]-5, miniDam$GetDisaster_sdate[1]), #miniDam$sdate[which(!is.na(miniDam$sdate))[1]], 
                                    fdate=data.table::fifelse(is.na(miniDam$GetDisaster_fdate[1]), miniDam$fdate[1]+21, miniDam$GetDisaster_fdate[1]), #miniDam$fdate[which(!is.na(miniDam$sdate))[1]], 
                                    eventid=i, hazard=miniDam$hazard[which(!is.na(miniDam$sdate))[1]])
    
    sdate <- miniDamSimplified$sdate[1]
    maxdate<-miniDamSimplified$sdate[1]-5
    if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
    
    if (!is.na(SubNatDataByEvent[[i]]$GetDisasterArgs_bbox[1])){
      bbox = eval(parse(text = SubNatDataByEvent[[i]]$GetDisasterArgs_bbox[1]))
    } else {
      bbox<-countriesbbox(unique(miniDamSimplified$iso3))
      bbox[1] <- ifelse( (bbox[1]-5)  < (-180), 180 - (-180-(bbox[1]-5)), bbox[1]-5)
      bbox[2] <- max(bbox[2]-5, -90)
      bbox[3] <- ifelse( (bbox[3]+5)  > 180, -180 + (bbox[3]+5-180), bbox[3] + 5)
      bbox[4] <- min(bbox[4]+5, 90)
    }
    
    if (!is.na(SubNatDataByEvent[[i]]$GetDisasterArgs_EQ_Params[1])){
      EQparams <- eval(parse(text=SubNatDataByEvent[[i]]$GetDisasterArgs_EQ_Params[1]))
    } else {
      EQparams=list(I0=4.3, minmag=5)
    }
    
    
    # lhazSDF <- GetDisaster(miniDamSimplified, bbox = c(160, -11.8, 163, -9.7))
    lhazSDF<-tryCatch(GetDisaster(miniDamSimplified,bbox=bbox, EQparams = EQparams),error=function(e) NULL)
    if(is.null(lhazSDF)) {
      file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', lhazSDF not found.'), file_conn)
      close(file_conn) 
      next
    }
    
    
    # Create the ODD object:
    ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified, agg_level=1),error=function(e) NULL)
    if(is.null(ODDy)) {
      file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', ODD object not created.'), file_conn)
      close(file_conn) 
      next
    }
    
    #Fetch building count data:
    #ODDy_build <- tryCatch(AddBuildingCounts(ODDy, i, paste0(folder_write, 'Building_count_notes')), error=function(e) NULL)
    ODDy_build <- tryCatch(getBingBuildingsGlobal(ODDy, i, paste0(folder_write, 'Building_count_notes')), error=function(e) NULL)
    if(is.null(ODDy_build)) {
      file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', Building counts not added.'), file_conn)
      close(file_conn)
    } else {
      ODDy <- ODDy_build
      rm(ODDy_build)
    }
    
    iso3_ODDy <- unique(ODDy$ISO3C)
    
    # CHECK TO MAKE SURE OLD DISP DATA IS CONTAINED IN THE EVENT
    DispData_event <- filter(DispData, (iso3 %in% iso3_ODDy) &  (as.Date(sdate) >= maxdate) & (as.Date(sdate) <= mindate))
    if (NROW(DispData_event)==0){
      file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', no DispData for event.'), file_conn)
      close(file_conn)
    } else {
      for (j in 1:length(DispData_event$USGSid)){
        lhazSDF_DispData <- tryCatch(GetUSGS_id(DispData_event$USGSid[j]),error=function(e) NULL)
        if(is.null(lhazSDF_DispData)) {
          file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
          writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', DispData lhazSDF not found.'), file_conn)
          close(file_conn)
        } else {
          cellsize <- 0.008333333333333
          if (lhazSDF_DispData@bbox[1,1] < (lhazSDF$hazard_info$bbox[1]-cellsize) | lhazSDF_DispData@bbox[1,2] > (lhazSDF$hazard_info$bbox[3]+cellsize) | lhazSDF_DispData@bbox[2,1] < (lhazSDF$hazard_info$bbox[2]-cellsize) | lhazSDF_DispData@bbox[2,2] > (lhazSDF$hazard_info$bbox[4]+cellsize)){
            file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
            writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', ODD bbox does not contain by DispData bbox.'), file_conn)
            close(file_conn)
          }
        }
      }
      rm(lhazSDF_DispData)
    }
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDamSimplified$iso3)[which(unique(miniDamSimplified$iso3) !='TOT')][1],
                  "_",i)
    HAZARDpath<-paste0(dir,folder_write, "HAZARDobjects/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    
    ODDy_with_impact <- tryCatch(updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1], i),error=function(e) NULL)
    
    if(is.null(ODDy_with_impact)) {
      file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', impact not added.'), file_conn)
      close(file_conn)
      next
    }
    
    ODDy <- ODDy_with_impact
    rm(ODDy_with_impact)
    
    
    # Save out objects to save on RAM
    ODDpath<-paste0(dir, folder_write, "ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    additional_poly_check(ODDy, i, print_to_xl=T)
    
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDamSimplified$iso3) & 
                               sdate<mindate & sdate>maxdate)
    
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {
        file_conn <- file(paste0(folder_write, 'ODD_creation_notes'), open = "a")
        writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', BD object creation failed.'), file_conn)
        close(file_conn)
      }
      BDpath <-paste0(dir, folder_write, "BDobjects/",namer)
      # Save it out!
      saveRDS(BDy, BDpath)
      rm(BDy)
    }
    
    
    # Save some RAM
    rm(ODDy)
  }

  return(path)
  
}


#creates hazard objects with additional info e.g. depth, magnitude, max_mmi, time of day
createHAZARD <- function(dir, haz="EQ", subnat_file= 'EQ_SubNational.xlsx'){
  # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if
  # no corresponding existing ODD object can be found, creates a new ODD object.
  
  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(EventID) %>% group_split()
  
  DispData <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/DispData_EQ_V2.Rdata')
  
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (i in 138:170){
    # if (!i %in% c(8, 13, 25, 28, 31, 39, 47, 52, 54, 59, 65, 67, 72, 74, 79, 80, 81, 85, 89, 
    #               109, 111, 116, 118, 123, 124, 125, 128, 131, 132, 133, 135, 137, 144, 145,
    #               151, 157, 164, 166)){
    #   next
    # }
    if (i %in% c(1, 126)){
      next
    }
    print(i)
    #Needs redoing: 69, 169, 170
    #Needs redoing: c(20, 52, 68, 70, 91, 92, 137, 139, 140, 144, 162, 170)
    #Needs redoing: c(20, 52, 67, 68, 70, 75, 76, 81, 82, 84, 91:122, 125, 130, 139, 140, 144, 162, 170)
    #failed: 67, 68, 137, 169 
    ###89 (lhazsdf not found), 78, 68, 69, 45, 35, 21, 16, 131,
    #lhazsdf_dispdat not found: 8, 9, 76, 84
    #no dispdat: 3,4, 30, 34, 45, 55, 71, 77, 82, 90, 91, 92, 129, 136, 164, 165, 166, 167, 168, 169, 170
    #ODD bbox not contained by disp data bbox: 7, 48, 52, 66, 67, 86, 104, 107, 108, 109, 117 (x2), 133, 145
    #Impact larger than 0 outside of I > I0: 92, 97, 114
    
    
    #building counts not added: 24, 31, 36, 55, 56, 57, 81, 99, 101, 114
    #missing bing building count: peru 2013-02-09, peru 2016-04-16, 2018-09-07
    #check weights:
    #   - Philippines 2013-10-15 Iloilo, Jordan
    #BD failed: 2019-08-08
    
    # Subset displacement and disaster database objects
    SubNatDataByEvent[[i]] <- SubNatDataByEvent[[i]][rowSums(is.na(SubNatDataByEvent[[i]][,-1])) != ncol(SubNatDataByEvent[[i]])-1, ]
    miniDam<-SubNatDataByEvent[[i]]
    
    miniDamSimplified <- data.frame(iso3=unique(trimws(unlist(strsplit(miniDam$source_listed_iso3, ";")))), 
                                    sdate=data.table::fifelse(is.na(miniDam$GetDisaster_sdate[1]), miniDam$sdate[1]-5, miniDam$GetDisaster_sdate[1]), #miniDam$sdate[which(!is.na(miniDam$sdate))[1]], 
                                    fdate=data.table::fifelse(is.na(miniDam$GetDisaster_fdate[1]), miniDam$fdate[1]+21, miniDam$GetDisaster_fdate[1]), #miniDam$fdate[which(!is.na(miniDam$sdate))[1]], 
                                    eventid=i, hazard=miniDam$hazard[which(!is.na(miniDam$sdate))[1]])
    
    sdate <- miniDamSimplified$sdate[1]
    maxdate<-miniDamSimplified$sdate[1]-5
    if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
    
    if (!is.na(SubNatDataByEvent[[i]]$GetDisasterArgs_bbox[1])){
      bbox = eval(parse(text = SubNatDataByEvent[[i]]$GetDisasterArgs_bbox[1]))
    } else {
      bbox<-countriesbbox(unique(miniDamSimplified$iso3))
      bbox[1] <- ifelse( (bbox[1]-5)  < (-180), 180 - (-180-(bbox[1]-5)), bbox[1]-5)
      bbox[2] <- max(bbox[2]-5, -90)
      bbox[3] <- ifelse( (bbox[3]+5)  > 180, -180 + (bbox[3]+5-180), bbox[3] + 5)
      bbox[4] <- min(bbox[4]+5, 90)
    }
    
    if (!is.na(SubNatDataByEvent[[i]]$GetDisasterArgs_EQ_Params[1])){
      EQparams <- eval(parse(text=SubNatDataByEvent[[i]]$GetDisasterArgs_EQ_Params[1]))
    } else {
      EQparams=list(I0=4.3, minmag=5)
    }
    
    
    # lhazSDF <- GetDisaster(miniDamSimplified, bbox = c(160, -11.8, 163, -9.7))
    lhazSDF<-tryCatch(GetDisaster(miniDamSimplified,bbox=bbox, EQparams = EQparams),error=function(e) NULL)
    if(is.null(lhazSDF)) {
      stop(i)
    }
    print(lhazSDF$hazard_info$max_mmi)
    lhazSDF$hazard_info$first_event <- check_preceding_hazards(lhazSDF)
    lhazSDF$hazard_info$first_event
    
    # Create a unique hazard event name
    namer<-paste0('EQ',
                  str_remove_all(as.character.Date(min(lhazSDF$hazard_info$eventdates)),"-"),
                  unique(miniDamSimplified$iso3)[which(unique(miniDamSimplified$iso3) !='TOT')][1],
                  "_",i)
    HAZARDpath<-paste0(dir,folder_write, "HAZARDobjects_PAGERfull/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
  }
  
  return(path)
  
}



editHAZARDS <- function(dir, haz="EQ", subnat_file= 'EQ_SubNational.xlsx'){
  
  folderin <- paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/HAZARDobjects_additionalInfo3/')
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
  for (namer in ufiles){
    hazobj <- readRDS(paste0(folderin, namer))
    #hazobj$hazard_info$first_event <- check_preceding_hazards(hazobj)
    saveRDS(lhazSDF,paste0(folderin, namer))
  }
  return(path)
  
}

createBDonly <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
  BD_folderin<-paste0(dir, folder_in, '/BDobjects/')
  BD_folderout1<-paste0(dir, folder_in, '/BDobjects_newVuln/')
  BD_folderout2<-paste0(dir, folder_in, '/BDobjects_newVuln_MARbySource/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  
  options(timeout=100)
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  for (jj in 1:length(ufiles)){
    #issues with 69
    file <- ufiles[jj]
    ODDy <- readRDS(paste0(ODD_folderin, file))

    #BDy <- readRDS(paste0(BD_folderin, file))
    
    miniDam<-Damage%>%filter(iso3%in%unique(ODDy$ISO3C) & 
                               sdate>(min(ODDy@hazdates)-1) & sdate<(max(ODDy@hazdates)+1))
    
    if (NROW(miniDam)==0) next
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {
        stop('Failed to create BD object')
      }
      
      #BDpath1 <-paste0(BD_folderout1,file)
      #saveRDS(BDy, BDpath1)
      
      file_conn <- file(paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/BD_creation_notes'), open = "a")
      writeLines(paste0("Event: ", file), file_conn)
      close(file_conn) 
      
      #BD_new <- tryCatch(BD_increase_coverage_bing(BDy, ODDy),error=function(e) NULL)
      BD_new <- select_MAR_polys(BDy, ODDy)
      if (is.null(BD_new)){
        file_conn <- file(paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/BD_creation_notes'), open = "a")
        writeLines('Error in creation', file_conn)
        close(file_conn) 
      }
      BDpath2 <-paste0(BD_folderout2,file)
      # Save it out!
      saveRDS(BD_new, BDpath2)
      
      rm(BDy)
      rm(BD_new)
    }
  }
}

fix_Zero_cIndies <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects_cIndiesFixed/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  for (file in ufiles){
    ODDy <- readRDS(paste0(ODD_folderin, file))
    if(any(ODDy@cIndies$value==0)){
      print(paste(ODDy@hazdates, unique(ODDy$ISO3C)))
      ODDy@cIndies$value[which(ODDy@cIndies$value==0)] = 0.0001
    }
    saveRDS(ODDy, paste0(ODD_folderout, file))
  }
  
  BD_folderin<-paste0(dir, folder_in, '/BDobjects/')
  BD_folderout<-paste0(dir, folder_in, '/BDobjects_cIndiesFixed//')
  ufiles<-list.files(path=BD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  for (file in ufiles){
    BDy <- readRDS(paste0(BD_folderin, file))
    if(any(BDy@cIndies$value==0)){
      print(paste(BDy@hazdates, unique(BDy$ISO3C)))
      BDy@cIndies$value[which(BDy@cIndies$value==0)] = 0.0001
    }
    saveRDS(BDy, paste0(BD_folderout, file))
  }
}



ODD <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_June20/ODDobjects/EQ20180908CHN_96')
ODD <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_June20/ODDobjects/EQ20171112IRN_70')


removeWeights <- function(ODD){
  pixels_of_interest <- which(!is.na(ODD$ISO3C))
  for (i in pixels_of_interest){
    weight_sum = c(0, 0, 0)
    polygon_matches = list(c(), c(), c())
    for (p in 1:length(ODD@polygons)){
      poly <- ODD@polygons[[p]]
      match_i <- which(poly$indexes == i)
      if (length(match_i) > 0){
        gadm_level <- str_count(poly$name, ',')
        polygon_matches[[gadm_level+1]] %<>% rbind(c(p, match_i))
        weight_sum[gadm_level+1] <- weight_sum[gadm_level+1] + poly$weights[match_i]
      }
    }
    for (g in 2:3){
      if (is.null(polygon_matches[[g]])) next
      if (NROW(polygon_matches[[g]])==1) next
      if (round(weight_sum[g],3)>1) stop()
      weights <- rep(0, NROW(polygon_matches[[g]]))
      for (j in 1:NROW(polygon_matches[[g]])){
        weights[j] <- ODD@polygons[[polygon_matches[[g]][j,1]]]$weights[polygon_matches[[g]][j,2]]
      }
      j_max_weight <- which.max(weights)
      for (j in 1:NROW(polygon_matches[[g]])){
        if (j == j_max_weight){
          ODD@polygons[[polygon_matches[[g]][j,1]]]$weights[polygon_matches[[g]][j,2]] = 1
        } else {
          ODD@polygons[[polygon_matches[[g]][j,1]]]$indexes <- ODD@polygons[[polygon_matches[[g]][j,1]]]$indexes[-polygon_matches[[g]][j,2]]
          ODD@polygons[[polygon_matches[[g]][j,1]]]$weights <- ODD@polygons[[polygon_matches[[g]][j,1]]]$weights[-polygon_matches[[g]][j,2]]
        }
      }
    }
  }
  return(ODD)
}

increaseAggregation <- function(ODD){
  ODD$PDens <- exp(round(log(ODD$PDens)))
  ODD$Vs30 <- round(ODD$Vs30, -2)
  ODD$EQFreq <- exp(round(log(ODD$EQFreq+0.1)/0.2, 0)*0.2)-0.1
  hrange<-grep("hazMean",names(ODD),value = T)
  ODD@data[,hrange] <- log(round(1.3^(1.5*ODD@data[,hrange])),base=1.3)/1.5  #round(ODD@data[,hrange]/0.2)*0.2 #
  
  polyMatch <- list()
  polyMatch[1:NROW(ODD@data)] <- ''
  for (i in 1:length(ODD@polygons)){
    polyMatch[ODD@polygons[[i]]$indexes] <- paste0(polyMatch[ODD@polygons[[i]]$indexes], i,',')
  }
  ODD@data$polyMatch <- unlist(polyMatch)
  grouped_by_covar <- ODD@data %>% group_by(across(all_of(c(hrange, 'ISO3C',  'PDens', 'Vs30', 'EQFreq', 'AveSchYrs', 'LifeExp', 'GNIc', 'SHDI', 'polyMatch'))))
  
  if (!is.null(ODD@data$nBuildings)){
    summarised <- grouped_by_covar %>% summarize(Population=sum(Population, na.rm=T), nBuildings=sum(nBuildings, na.rm=T))
  } else {
    summarised <- grouped_by_covar %>% summarize(Population=sum(Population, na.rm=T))
  }
  
  ODDagg <- ODD
  ODDagg@data <- as.data.frame(summarised)
  for (i in 1:length(ODDagg@polygons)){
    ODDagg@polygons[[i]]$indexes <- c()
  }
  
  for (i in 1:NROW(ODDagg@data)){
    match_polys <- as.numeric(unlist(str_extract_all(ODDagg@data$polyMatch[i], "\\d+")))
    for (j in match_polys){
      ODDagg@polygons[[j]]$indexes <- c(ODDagg@polygons[[j]]$indexes, i)
    }
  }
  for (i in 1:length(ODDagg@polygons)){
    ODDagg@polygons[[i]]$weights <- rep(1, length(ODDagg@polygons[[i]]$indexes))
  }
  return(ODDagg)
  
  #notnans <- which(!apply(is.na(ODD@data[,c('Vs30', 'EQFreq', 'AveSchYrs', 'LifeExp', 'GNIc')]), 1, any) | !apply(is.na(ODD@data[,hrange]), 1, all))
  #covar_df <- ODD@data[,c(hrange, 'Vs30', 'EQFreq', 'AveSchYrs', 'LifeExp', 'GNIc')]
  #unique_covars <- unique(covar_df)
  #new_to_old_indexes <- list()
  #rows_remaining <- 1:NROW(covar_df)
  # group_by(covar_df, all_of(hrange))
  # for (i in 1:NROW(unique_covars)){
  #   print(i)
  #   index_match <- rows_remaining
  #   for (j in 1:NCOL(unique_covars)){
  #     index_match <- intersect(index_match, rows_remaining[which(covar_df[rows_remaining,j] == unique_covars[i,j] | (is.na(unique_covars[i,j]) & is.na(covar_df[rows_remaining,j])))])
  #   }
  #   new_to_old_indexes[[i]] <- index_match
  #   rows_remaining <- setdiff(rows_remaining, index_match)
  # }   
    
}

increaseAggregation_all <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_folderout<-paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_Aug31_Agg', '/ODDobjects/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  for (file in ufiles[1:length(ufiles)]){
    event_id <- as.numeric(strsplit(file, "_")[[1]][2])
    ODDy <- readRDS(paste0(ODD_folderin, file))
    saveRDS(ODDy, paste0(dir, folder_in, '/ODDobjects/', file))
    ODDyAgg <- tryCatch(increaseAggregation(removeWeights(ODDy)),error=function(e) NULL)
    if(is.null(ODDyAgg)){
      print(event_id)
      next
    }
    saveRDS(ODDyAgg, paste0(ODD_folderout, file))
  }
  
}

moveTestData <- function(folder_in='IIDIPUS_Input'){
  ODD_folderall<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_foldertest<-paste0(dir, folder_in, '/ODDobjects/Test/')
  ufiles<-list.files(path=ODD_folderall,pattern=Model$haz,recursive = T,ignore.case = T)
  i <- 0
  for (file in ufiles){
    i <- i + 1
    if (i %%3 != 0){next}
    file.copy(from = paste0(ODD_folderall, file),
              to = paste0(ODD_foldertest, file))
    file.remove(from = paste0(ODD_folderall, file))
  }
  BD_folderall<-paste0(dir, folder_in, '/BDobjects/')
  BD_foldertest<-paste0(dir, folder_in, '/BDobjects/Test/')
  ufiles<-list.files(path=BD_folderall,pattern=Model$haz,recursive = T,ignore.case = T)
  i <- 0
  for (file in ufiles){
    i <- i + 1
    if (i %%3 != 0){next}
    file.copy(from = paste0(BD_folderall, file),
              to = paste0(BD_foldertest, file))
    file.remove(from = paste0(BD_folderall, file))
  }
  
}


# compare impact for aggregated and non-aggregated approaches
ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/EQ20151207TJK_35')
ODDyAgg <- ODDy %>% removeWeights() %>% increaseAggregation()
for (i in 1:length(ODDyAgg@polygons)){
  ODDyAgg@polygons[[i]]$weights <- rep(1, length(ODDyAgg@polygons[[i]]$indexes))
}


Disps <- DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, Method = AlgoParams)
DispsAgg <- DispX(ODD = ODDyAgg,Omega = Omega %>% addTransfParams(),center = Model$center, Method = AlgoParams)

samples <- data.frame(impact= Disps[[1]]$impact, observed= Disps[[1]]$observed)
samples_agg <- data.frame(impact= DispsAgg[[1]]$impact, observed= DispsAgg[[1]]$observed)
for (i in 1:length(Disps)){
  samples %<>% add_column(Disps[[i]]$sampled)
  samples_agg %<>% add_column(DispsAgg[[i]]$sampled)
}

plot_df <- data.frame(
  sampled_med = apply(samples[,3:NCOL(samples)], 1, median),
  sampled_min = apply(samples[,3:NCOL(samples)], 1, min),
  sampled_max = apply(samples[,3:NCOL(samples)], 1, max),
  agg_sampled_med = apply(samples_agg[,3:NCOL(samples_agg)], 1, median),
  agg_sampled_min = apply(samples_agg[,3:NCOL(samples_agg)], 1, min),
  agg_sampled_max = apply(samples_agg[,3:NCOL(samples_agg)], 1, max)
)
plot_df$poly <- 1:NROW(plot_df)

ggplot(plot_df) + 
  geom_point(aes(x=poly, y=log(sampled_med+0.1))) + 
  geom_errorbar(aes(x=poly, ymin=log(sampled_min+0.1), ymax=log(sampled_max+0.1), width=0.2)) + 
  geom_point(aes(x=poly+0.2, y=log(agg_sampled_med+0.1)), col='red') +
  geom_errorbar(aes(x=poly+0.2, ymin=log(agg_sampled_min+0.1), ymax=log(agg_sampled_max+0.1), width=0.2), col='red')



remove_partially_missing_buildings <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July21_Agg'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects_InferredFlag/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  SubNatData <- readSubNatData(subnat_file)
  for (file in ufiles){
    print(file)
    ODDy <- readRDS(paste0(ODD_folderin, file))
    impact <- ODDy@impact
    SubNatData_filt <- SubNatData %>% filter(EventID==as.numeric(strsplit(file, "_")[[1]][2]))
    SubNatData_filt$poly_name <- SubNatData_filt$country
    if(length(unique(impact$iso3))==1){
      SubNatData_filt$poly_name <- ifelse(SubNatData_filt$poly_name=='TOTAL',  ifelse(unique(impact$iso3)=='TOT', 'TOTAL', countrycode(unique(impact$iso3), origin='iso3c', destination='country.name')),SubNatData_filt$poly_name)
    }
    SubNatData_filt$poly_name <- ifelse(is.na(SubNatData_filt$Region), SubNatData_filt$poly_name, paste0(SubNatData_filt$Region, ', ', SubNatData_filt$poly_name))
    SubNatData_filt$poly_name <- ifelse(is.na(SubNatData_filt$Subregion), SubNatData_filt$poly_name, paste0(SubNatData_filt$Subregion, ', ', SubNatData_filt$poly_name))
    impact$inferred <- NA
    for (i in 1:NROW(impact)){
      matched_row <- SubNatData_filt %>% filter(poly_name==ODDy@polygons[[impact$polygon[i]]]$name & SubNatData_filt[[impact$impact[i]]]==impact$observed[i])
      if(NROW(matched_row) > 1){
        if(length(unique(matched_row[[paste0(impact$impact[i], 'Inferred')]]))==1){
          matched_row <- matched_row[1,]
        } else {
          stop('More than one matching row')
        }
      }
      impact$inferred[i] <- matched_row[[paste0(impact$impact[i], 'Inferred')]]
    }
    ODDy@impact <- impact
    
    saveRDS(ODDy, paste0(ODD_folderout, file))
  }
}

fix_Zero_cIndies <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July21_Agg'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects_withDest/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  for (j in 168:length(ufiles)){
    file <- ufiles[j]
    ODDy <- readRDS(paste0(ODD_folderin, file))
    event_id <- as.numeric(strsplit(file, "_")[[1]][2])
    
    if (event_id %in% c(22,40, 45, 46, 47, 53, 57, 60, 75, 83, 84, 87, 88, 92, 99, 104, 107, 119, 121, 134, 140, 169)){
      stop()
    }
    
    buildDam_impact = ODDy@impact %>% filter(impact == 'buildDam' | impact=='buildDest') %>% group_by(iso3, sdate, build_type, polygon) %>% 
      summarise(observed = sum(observed), impact='buildDam', qualifier = ifelse(length(unique(qualifier)) > 1, 'mixed', first(qualifier)), 
                inferred=ifelse(any(inferred), T, F))
    
    if (length(buildDam_impact$polygon) != length(unique(buildDam_impact$polygon))){
      stop()
    }
    ODDy@impact <- rbind(ODDy@impact %>% filter(impact != 'buildDam' & impact != 'buildDest'), buildDam_impact)
    ODDy@impact$impact[which(ODDy@impact$impact=='buildDamDest')] = 'buildDam'
  
    saveRDS(ODDy, paste0(ODD_folderout, file))
  }
}

removeDestroyed <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  for (j in 1:length(ufiles)){
    file <- ufiles[j]
    ODDy <- readRDS(paste0(ODD_folderin, file))
    event_id <- as.numeric(strsplit(file, "_")[[1]][2])
    
    if (any(ODDy@impact$impact %in% c('buildDest', 'buildDamDest'))){
       stop()
    } else {
      next
    }
    if (event_id %in% c(46, 53, 75, 83, 88, 99, 104, 119, 140, 169)){
      stop()
    }
    
    if (event_id %in% c(45, 57,  60, 87, 92, 107, 134)){
      ODDy@impact %<>% filter(impact != 'buildDam' & impact != 'buildDest')
      saveRDS(ODDy, paste0(ODD_folderout, file))
      next
    }
    
    if (event_id %in% c(22, 40, 47)){
      ODDy@impact %<>% filter(impact != 'buildDest')
    }
    
    buildDam_impact = ODDy@impact %>% filter(impact == 'buildDam' | impact=='buildDest') %>% group_by(iso3, sdate, build_type, polygon) %>% 
      summarise(observed = sum(observed), impact='buildDam', qualifier = ifelse(length(unique(qualifier)) > 1, 'mixed', first(qualifier)), 
                inferred=ifelse(any(inferred), T, F))
    
    if (length(buildDam_impact$polygon) != length(unique(buildDam_impact$polygon))){
      stop()
    }
    ODDy@impact <- rbind(ODDy@impact %>% filter(impact != 'buildDam' & impact != 'buildDest'), buildDam_impact)
    ODDy@impact$impact[which(ODDy@impact$impact=='buildDamDest')] = 'buildDam'
    
    saveRDS(ODDy, paste0(ODD_folderout, file))
  }
}



updateVuln <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects_OldVuln/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  ufiles_haz<-na.omit(list.files(path=paste0(folder_in, '/HAZARDobjects_additionalInfo2/'),pattern=Model$haz,recursive = T,ignore.case = T))

  for (j in 1:length(ufiles)){
    
    event_id <- as.numeric(strsplit(ufiles[j], "_")[[1]][2])
    
    if (!event_id %in% c(30,37,52,79,137,144)){
      next
    }
    ODDy <- readRDS(paste0(ODD_folderin, ufiles[j]))
    ODDy$AveSchYrs <- NULL
    ODDy$LifeExp <- NULL
    ODDy$EQFreq <- NULL
    ODDy$Vs30 <- NULL
    ODDy$GNIc <- NULL
    ODDy %<>% AddVuln()
    event_id <- as.numeric(strsplit(ufiles[j], "_")[[1]][2])
    
    haz_match <- grep(paste0("_", event_id, "\\b"),  ufiles_haz, value = TRUE)
    HAZARDobj <- readRDS(paste0(folder_in, '/HAZARDobjects_additionalInfo2/', haz_match))
    
    ODDy@hazinfo <- HAZARDobj$hazard_info
    saveRDS(ODDy, paste0(ODD_folderout, ufiles[j]))
  }
}

updateHaz2 <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects_OldFirstHaz/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  ufiles_haz<-na.omit(list.files(path=paste0(folder_in, '/HAZARDobjects_additionalInfo3/'),pattern=Model$haz,recursive = T,ignore.case = T))
  
  for (j in 1:length(ufiles)){
    
    ODDy <- readRDS(paste0(ODD_folderin, ufiles[j]))
    event_id <- as.numeric(strsplit(ufiles[j], "_")[[1]][2])
    
    haz_match <- grep(paste0("_", event_id, "\\b"),  ufiles_haz, value = TRUE)
    HAZARDobj <- readRDS(paste0(folder_in, '/HAZARDobjects_additionalInfo3/', haz_match))
    
    #HAZARDobj$hazard_info$first_event = check_preceding_hazards(HAZARDobj)
    
    ODDy@hazinfo <- HAZARDobj$hazard_info
    saveRDS(ODDy, paste0(ODD_folderout, ufiles[j]))
    # saveRDS(HAZARDobj, paste0(folder_in, '/HAZARDobjects_additionalInfo2/', haz_match))
  }
}



# folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_June20'
# ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
# ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
# i <- i + 1
# print(i)
# ODDy <- readRDS(paste0(ODD_folderin, ufiles[i]))
# plotODDy(ODDy, var='EQFreq')



removeDestroyed <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July21_Agg/'){
  BD_folderin<-paste0(dir, folder_in, 'BDobjects_withDest/')
  BD_folderout<-paste0(dir, folder_in, 'BDobjects/')
  ufiles<-list.files(path=BD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  magnitudes <- c()
  for (file in ufiles){
    BDy <- readRDS(paste0(BD_folderin, file))
    print(unique(BDy$grading))
    BDy$grading[which(BDy$grading %in% c('moderate', 'severe','destroyed'))] = 'Damaged'
    print(unique(BDy$grading))
    saveRDS(BDy, paste0(BD_folderout, file))
  }

}





















SubNatData <- readSubNatData(subnat_file)

# Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
SubNatDataByEvent <- SubNatData %>% group_by(EventID) %>% group_split()


while (i < 169){
  
  folderin<- paste0(dir, 'IIDIPUS_Input_I0_equals_4/ODDobjects/') 
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
  
  file_in <- ufiles[which(sub(".*_", "", ufiles) == i)]
  
  ODDy <- readRDS(paste0(dir,"IIDIPUS_Input_I0_equals_4/ODDobjects/",file_in))
  HAZy <- readRDS(paste0(dir,"IIDIPUS_Input_I0_equals_4/HAZARDobjects/",file_in))
  
  SubNatDataByEvent[[i]] <- SubNatDataByEvent[[i]][rowSums(is.na(SubNatDataByEvent[[i]][,-1])) != ncol(SubNatDataByEvent[[i]])-1, ]
  miniDam<-SubNatDataByEvent[[i]]
  
  miniDamSimplified <- data.frame(iso3=unique(trimws(unlist(strsplit(miniDam$source_listed_iso3, ";")))), 
                                  sdate=data.table::fifelse(is.na(miniDam$GetDisaster_sdate[1]), miniDam$sdate[1]-5, miniDam$GetDisaster_sdate[1]), #miniDam$sdate[which(!is.na(miniDam$sdate))[1]], 
                                  fdate=data.table::fifelse(is.na(miniDam$GetDisaster_fdate[1]), miniDam$fdate[1]+21, miniDam$GetDisaster_fdate[1]), #miniDam$fdate[which(!is.na(miniDam$sdate))[1]], 
                                  eventid=i, hazard=miniDam$hazard[which(!is.na(miniDam$sdate))[1]])
  
  sdate <- miniDamSimplified$sdate[1]
  maxdate<-miniDamSimplified$sdate[1]-5
  if(is.na(miniDamSimplified$fdate[1])) mindate<-miniDamSimplified$sdate[1]+3 else mindate<-miniDamSimplified$fdate[1]+3
  
  if (!is.na(SubNatDataByEvent[[i]]$GetDisasterArgs_bbox[1])){
    bbox = eval(parse(text = SubNatDataByEvent[[i]]$GetDisasterArgs_bbox[1]))
  } else {
    bbox<-countriesbbox(unique(miniDamSimplified$iso3))
    bbox[1] <- ifelse( (bbox[1]-5)  < (-180), 180 - (-180-(bbox[1]-5)), bbox[1]-5)
    bbox[2] <- max(bbox[2]-5, -90)
    bbox[3] <- ifelse( (bbox[3]+5)  > 180, -180 + (bbox[3]+5-180), bbox[3] + 5)
    bbox[4] <- min(bbox[4]+5, 90)
  }
  
  if (!is.na(SubNatDataByEvent[[i]]$GetDisasterArgs_EQ_Params[1])){
    EQparams <- eval(parse(text=SubNatDataByEvent[[i]]$GetDisasterArgs_EQ_Params[1]))
  } else {
    EQparams=list(I0=4, minmag=5)
  }
  
  lhazSDF<-tryCatch(GetDisaster(miniDamSimplified,bbox=bbox, EQparams = EQparams),error=function(e) NULL)
  if(is.null(lhazSDF)) {
    file_conn <- file('ODD_creation_notes/Updated_bbox', open = "a")
    writeLines(paste0("lhazSDF not found: ", i), file_conn)
    close(file_conn) 
    #print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3), " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
    i <- i + 1
    next
  }
  
  length(lhazSDF$hazard_info$eventdates)
  
  miniDam$sdate[1]
  plot(0,0, xlim=bbox[c(1,3)],  ylim=bbox[c(2,4)])
  
  points(lhazSDF[[2]], col='red')
  lhazSDF[[1]]$eventdates[1]
  max(lhazSDF[[2]]$mean, na.rm=T)
  
  points(lhazSDF[[3]], col='blue')
  lhazSDF[[1]]$eventdates[2]
  max(lhazSDF[[3]]$mean, na.rm=T)
  
  points(lhazSDF[[4]], col='green')
  lhazSDF[[1]]$eventdates[3]
  max(lhazSDF[[4]]$mean, na.rm=T)
  
  points(lhazSDF[[5]], col='pink')
  lhazSDF[[1]]$eventdates[4]
  max(lhazSDF[[5]]$mean, na.rm=T)
  
  points(lhazSDF[[6]], col='yellow')
  max(lhazSDF[[6]]$mean, na.rm=T)
  lhazSDF[[1]]$eventdates[5]
  
  points(lhazSDF[[7]], col='orange')
  max(lhazSDF[[7]]$mean, na.rm=T)
  lhazSDF[[1]]$eventdates[6]
  
  points(lhazSDF[[8]], col='purple')
  max(lhazSDF[[8]]$mean, na.rm=T)
  lhazSDF[[1]]$eventdates[7]
  
  points(lhazSDF[[9]], col='grey')
  max(lhazSDF[[9]]$mean, na.rm=T)
  lhazSDF[[1]]$eventdates[8]
  
  points(lhazSDF[[10]], col='salmon')
  max(lhazSDF[[10]]$mean, na.rm=T)
  lhazSDF[[1]]$eventdates[9]
  
  points(lhazSDF[[11]], col='brown')
  max(lhazSDF[[11]]$mean, na.rm=T)
  lhazSDF[[1]]$eventdates[10]
  
  points(lhazSDF[[12]], col='maroon')
  max(lhazSDF[[12]]$mean, na.rm=T)
  lhazSDF[[1]]$eventdates[11]
  
  
  if (length(grep("hazMean",names(ODDy),value = T)) < length(lhazSDF) - 1){
    upd_bbox <- c(ODDy@bbox) + c(-0.5,-0.5, 0.5, 0.5)
    lhazSDF<-tryCatch(GetDisaster(miniDamSimplified,bbox=upd_bbox, EQparams = EQparams),error=function(e) NULL)
    
    if (length(grep("hazMean",names(ODDy),value = T)) == length(lhazSDF) - 1){
      file_conn <- file('ODD_creation_notes/Updated_bbox', open = "a")
      writeLines(paste0("Index: ", i, ', c(',paste(round(upd_bbox,2), collapse=','), ')'), file_conn)
      close(file_conn) 
    } else {
      file_conn <- file('ODD_creation_notes/Updated_bbox', open = "a")
      writeLines(paste0("Need to redo hazards for ODD index: ", i), file_conn)
      close(file_conn) 
      i <- i + 1
      next
    }
  }
  if (length(grep("hazMean",names(ODDy),value = T)) > length(lhazSDF) - 1){
    file_conn <- file('ODD_creation_notes/Updated_bbox', open = "a")
    writeLines(paste0("Less hazards found than on original ODD object ", i), file_conn)
    close(file_conn) 
    i <- i + 1
    next
  }
  
  if (!all(lhazSDF$hazard_info$bbox == HAZy$hazard_info$bbox)){
    file_conn <- file('ODD_creation_notes/Updated_bbox', open = "a")
    writeLines(paste0("Same number of hazards but different hazards?: ", i), file_conn)
    close(file_conn) 
  }
  i <- i + 1
}

miniDamSimplified
lhazSDF<-tryCatch(GetDisaster(miniDamSimplified,bbox=bbox, EQparams = EQparams),error=function(e) NULL)
length(lhazSDF$hazard_info$eventdates)

miniDamSimplified2 <- miniDamSimplified
miniDamSimplified2$fdate <- miniDamSimplified$fdate -2
lhazSDF2<-tryCatch(GetDisaster(miniDamSimplified2,bbox=bbox, EQparams = EQparams),error=function(e) NULL)
length(lhazSDF2$hazard_info$eventdates)

plot(lhazSDF[[2]])
points(lhazSDF[[3]])
points(lhazSDF[[4]], col='red')
points(lhazSDF2[[2]], col='green')
points(lhazSDF2[[3]], col='yellow')
  
  miniDam$sdate[1]
  plot(0,0, xlim=bbox[c(1,3)],  ylim=bbox[c(2,4)])
  points(ODDy@coords)
  points(lhazSDF[[2]], col='red')
  lhazSDF[[1]]$eventdates[1]
  points(lhazSDF[[3]], col='blue')
  lhazSDF[[1]]$eventdates[2]
  points(lhazSDF[[4]], col='green')
  lhazSDF[[1]]$eventdates[3]
  points(lhazSDF[[5]], col='pink')
  lhazSDF[[1]]$eventdates[4]
  points(lhazSDF[[6]], col='yellow')
  lhazSDF[[1]]$eventdates[5]
  points(lhazSDF[[7]], col='orange')
  lhazSDF[[1]]$eventdates[7]
  points(lhazSDF[[8]], col='purple')
  lhazSDF[[1]]$eventdates[8]
  
  namer<-paste0(ODDy@hazard,
                str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                unique(miniDamSimplified$iso3)[which(unique(miniDamSimplified$iso3) !='TOT')][1],
                "_",i)
  
  
  ODDpath<-paste0(dir,"IIDIPUS_Input_I0_equals_4/ODDobjects/",namer)
  saveRDS(ODDy,ODDpath)
  
  # Building damage subset
  miniDam<-Damage%>%filter(iso3%in%unique(miniDamSimplified$iso3) & 
                             sdate<mindate & sdate>maxdate)
  
  HAZARDpath<-paste0(dir,"IIDIPUS_Input_I0_equals_4/HAZARDobjects/",namer)
  saveRDS(lhazSDF,HAZARDpath)
  
  BDpath=NA_character_
  if(nrow(miniDam)>0) {
    # Make building damage object BD
    BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
    if(is.null(BDy)) {
      file_conn6 <- file('ODD_creation_notes/BD_write_failed', open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", sdate), file_conn6)
      close(file_conn6) 
    }
    BDpath <-paste0(dir,"IIDIPUS_Input_I0_equals_4/BDobjects/",namer)
    # Save it out!
    saveRDS(BDy, BDpath)
  }
}

AddBingBuildingCounts(ODDy, isos_bingbuildings, plot_only=T)
plotODDy(ODDy, var='nBuildings')
























rbPal <- colorRampPalette(c('red','blue'))
cols <- rbPal(10)[as.numeric(cut(ODDy@coords[,1]> 123.6 & ODDy@coords[,1]< 124.2 & ODDy@coords[,2]< 13 & ODDy@coords[,2] >  10 & ODDy$Population > 0)],breaks = 100))]

plot(ODDy@coords[which(ODDy@coords[,1]> 123.6 & ODDy@coords[,1]< 124.2 & ODDy@coords[,2]< 13 & ODDy@coords[,2] >  10 & ODDy$Population > 0),], col=cols)

library(randomcoloR)
library(scales)
plot_impact <- function(ODDy){
  impact_by_type <- ODDy@impact %>% group_by(impact) %>% group_split()
  if (length(impact_by_type) > 1){
    par(mfrow=c(2,2))
  }
  for (i in 1:length(impact_by_type)){
    impact_filtered <- impact_by_type[[i]]
    ncols <- NROW(impact_filtered)
    palette <- distinctColorPalette(ncols)
    plot(ODDy@coords[which(is.na(ODDy$ISO3C)),], col='black', xlim=ODDy@bbox[1,], ylim=ODDy@bbox[2,], main=impact_filtered$impact[1], pch=20)
    subnat <- F
    for (k in 1:length(ODDy@polygons)){
      if (grepl(',', ODDy@polygons[[k]]$name, fixed = TRUE)) subnat <- T
    }
    for (ij in 1:ncols){
      if (!grepl(',', ODDy@polygons[[impact_filtered$polygon[ij]]]$name, fixed = TRUE) & subnat==T) next
      points(ODDy@coords[ODDy@polygons[[impact_filtered$polygon[ij]]]$indexes,], col=alpha(palette[ij], 1), pch=20, cex=0.5)
    }
  }
  par(mfrow=c(1,1))
}


plot(ODDy@coords[which(!is.na(ODDy$Population)),], col='black')
points(ODDy@coords[ODDy@polygons[[9]]$indexes,, drop=F], col='green')
points(ODDy@coords[ODDy@polygons[[16]]$indexes,, drop=F], col='blue')
points(ODDy@coords[ODDy@polygons[[17]]$indexes,, drop=F], col='pink')
points(ODDy@coords[ODDy@polygons[[19]]$indexes,, drop=F], col='orange')
points(ODDy@coords[ODDy@polygons[[20]]$indexes,, drop=F], col='yellow')
points(ODDy@coords[ODDy@polygons[[22]]$indexes,, drop=F], col='purple')
points(ODDy@coords[ODDy@polygons[[24]]$indexes,, drop=F], col='grey')
points(ODDy@coords[ODDy@polygons[[8]]$indexes,, drop=F], col='cyan')
points(ODDy@coords[ODDy@polygons[[9]]$indexes,, drop=F], col='green', pch=4)
points(ODDy@coords[ODDy@polygons[[10]]$indexes,, drop=F], col='blue', pch=4)
points(ODDy@coords[ODDy@polygons[[11]]$indexes,, drop=F], col='pink', pch=4)
points(ODDy@coords[ODDy@polygons[[12]]$indexes,, drop=F], col='orange', pch=4)
points(ODDy@coords[ODDy@polygons[[13]]$indexes,, drop=F], col='yellow', pch=4)
points(ODDy@coords[ODDy@polygons[[14]]$indexes,, drop=F], col='purple', pch=4)
points(ODDy@coords[ODDy@polygons[[15]]$indexes,, drop=F], col='grey', pch=4)
points(ODDy@coords[ODDy@polygons[[16]]$indexes,, drop=F], col='cyan', pch=4)

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

