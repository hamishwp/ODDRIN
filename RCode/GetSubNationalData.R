############################
### GetSubNationalData.R ###
############################

# Create all HAZARD, ODD and BD objects for events in EQ_SubNational.xlsx

# Functions:
#    -readSubNatData(): reads and cleans xlsx file with subnational data
#    -getPolyData(): gets SpatialPolygonsDataFrame from polygon name
#    -pixelsinpoly(): Determine which ODD object pixels are inside a SpatialPolygonsDataFrame
#    -reweight_pixels(): For pixels that have been allocated to more than one polygon, reweight based on the proportion of the pixel lying in the polygon
#    -pixel_weight_border_correction(): Increases weight of pixels on the boundary (e.g it's been given weight of 0.5 but should be 1 as half the pixel is in the sea)
#    -addODDPolygons(): Populates the polygons slot in the ODD object
#    -addODDImpact(): APopulates the impact slot of an ODD object


# Where would we like to log error messages?
#folder_write <- 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_Sept26/'

readSubNatData <- function(subnat_file){
  # Clean data from xlsx file by converting data to the appropriate formats/data types
  
  SubNatData <-  read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  x <-xlsx_cells(paste0(dir, 'IIDIPUS_Input/', subnat_file))
  formats <- xlsx_formats(paste0(dir, 'IIDIPUS_Input/', subnat_file))

  #Set the values in all red cells to NA:
  red_cells <- x %>% filter(local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == "FFFF0000")) %>% dplyr::select(row, col)
  red_cells$row <- red_cells$row - 1
  red_cells <- red_cells[-which(red_cells$col > NCOL(SubNatData)), ]
  SubNatData[as.matrix(red_cells)] <- NA
  
  #Mark which observations are 'inferred' (e.g. assumed 0s due to non-reporting):
  pink_cells <- x %>% filter(local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == "FFFF00FF")) %>% dplyr::select(row, col)
  SubNatData$buildDamInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='buildDam_exlusion_reason'))]-1)
  SubNatData$buildDestInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='buildDest_exlusion_reason'))]-1)
  SubNatData$buildDamDestInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='buildDamDest_exlusion_reason'))]-1)
  SubNatData$mortalityInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='mortality_exlusion_reason'))]-1)
  SubNatData$displacementInferred <- 1:NROW(SubNatData) %in% (pink_cells$row[which(pink_cells$col==which(names(SubNatData)=='displacement_exlusion_reason'))]-1)
  
  SubNatData <- SubNatData[!is.na(SubNatData$iso3), ]
  
  #ensure data is in the correct format:
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
  
  SubNatData %<>% combineDamagedDestroyedBuildings()
  return(SubNatData)
}

combineDamagedDestroyedBuildings <- function(SubNatData){
  # We treat 'Damaged Buildings' as Damaged + Destroyed buildings, but in the spreadsheet they are separated into
  # columns for Damaged only, Destroyed only or (Damaged + Destroyed) depending on what data is available
  for (i in 1:NROW(SubNatData)){
     if (!is.na(SubNatData$buildDamDest[i])){
       #if buildDamDest is available use that
       SubNatData$buildDam[i] = SubNatData$buildDamDest[i]
       SubNatData$buildDamInferred[i] = SubNatData$buildDamDestInferred[i]
       SubNatData$buildDam_qualifier[i] = SubNatData$buildDamDest_qualifier[i]
     } else if (!is.na(SubNatData$buildDam[i])){
         #otherwise sum over buildDam and buildDest
         SubNatData$buildDam[i]=sum(SubNatData$buildDam[i], SubNatData$buildDest[i], na.rm=T)
         SubNatData$buildDamInferred[i] = any(SubNatData$buildDamInferred[i], ifelse(is.na(SubNatData$buildDestInferred[i]), F, SubNatData$buildDestInferred[i]))
         SubNatData$buildDam_qualifier[i] = SubNatData$buildDam_qualifier[i]
      }
  }
  return(SubNatData)
}

getPolyData <- function(polygon_name, subregion, region, country, iso3){
  # Uses getGADM() or getbb() to retrieve the polygon of a region
  # Details:
  #    - Polygon_name is just used to label the polygon at the end (and isn't used to actually find the region)
  #    - Any missing values for subregion or region should be set to NA (e.g. when looking at national data)
  #    - Returns a list containing:
  #         - $polygon_name set to polygon_name
  #         - $sf_polygon set to a polygon (that works with the sf package) corresponding to this region, or NULL if polygon not found
  #    - Attempts to use ecochange::getGADM() function but if this fails will attempt the osmdata::getbb() function
  
  if (country=='TOTAL'){
    return(list(polygon_name = polygon_name, sf_polygon = NULL))
  }
  
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
    GADM_nowhitespace <- gsub(" ", "",GADM_array, fixed = TRUE)
    polygon <- tryCatch(getGADM(unit.nm=GADM_nowhitespace, level=GADM_level, country=GADM_iso3),error=function(e) NULL) 
  }
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

pixelsInPoly <- function(poly, ODDy){
  #returns the indexes of the ODDy pixels are inside the polygon 'poly'
    
  poly_sp <- SpatialPolygons(poly)
  proj4string(poly_sp)<- crs(ODDy)
  poly_sf <- st_as_sf(poly_sp)
  
  indexes_rast <- ODDy[['Population']]
  values(indexes_rast) <- 1:ncell(ODDy)
  if (is.null(intersect(extent(ext(indexes_rast)[1], ext(indexes_rast)[2], ext(indexes_rast)[3], ext(indexes_rast)[4]), extent(poly_sf)))){return(c())}
  indexes_rast <- crop(indexes_rast, extent(poly_sf), snap='out')
  pixels_as_poly <- as.polygons(rast(indexes_rast),aggregate=F)
  intersects <- as.matrix(st_intersects(st_make_valid(st_as_sf(pixels_as_poly)), st_make_valid(poly_sf)))
  return(values(indexes_rast)[which(intersects)])
  
  # insidepoly<-rep(FALSE,ncell(ODDy))
  # lon_cellsize <- res(ODDy)[1]
  # lat_cellsize <- res(ODDy)[2]
  # 
  # dummy_rast <- ODDy[['Population']]
  # values(dummy_rast) <- F
  # values(dummy_rast)[insidepoly] <- T
  # 
  # subst(rast(ODDy), 1, 1, others=NA)
  # ODDy[['Population']] %>% terra::as.polygons() %>% st_as_sf() %>% filter(value==1)
  # 
  # # Get rid of values outside the bounding box first
  # coords <- xyFromCell(ODDy, 1:ncell(ODDy))
  # indies<- which(coords[,1] >= (poly_sp@bbox[1,1]-lon_cellsize/2) &
  #   coords[,1]<= (poly_sp@bbox[1,2]+ lon_cellsize/2) &
  #   coords[,2]>= (poly_sp@bbox[2,1]-lat_cellsize/2) &
  #   coords[,2]<= (poly_sp@bbox[2,2]+lon_cellsize/2))
  # 
  # if (length(indies)==0){return(integer())}
  # 
  # #sf_use_s2(FALSE) #Error being raised with st_intersects when this is set to TRUE
  # 
  # sf_polygon
  # 
  # 
  # check_intersection <- function(j){
  #   pixel_sf <- sfheaders::sf_polygon(data.frame(longitude = c(coords[j,1]-lon_cellsize/2, 
  #                                                              coords[j,1]-lon_cellsize/2,
  #                                                              coords[j,1]+lon_cellsize/2, 
  #                                                              coords[j,1]+lon_cellsize/2, 
  #                                                              coords[j,1]-lon_cellsize/2), 
  #                                                latitude = c(coords[j,2]-lat_cellsize/2, 
  #                                                             coords[j,2]+lat_cellsize/2,
  #                                                             coords[j,2]+lat_cellsize/2, 
  #                                                             coords[j,2]-lat_cellsize/2,
  #                                                             coords[j,2]-lat_cellsize/2)), x='longitude',y='latitude')
  #   sf::st_crs(pixel_sf) = proj4string(ODDy)
  #   
  #   if(st_intersects(pixel_sf, poly_sf, sparse=F)[1,1]){
  #     return(T)
  #   } else return(F)
  # }
  # 
  # inside_poly <- unlist(mclapply(indies, check_intersection, mc.cores=2))
  # 
  # sf_use_s2(TRUE) #Change back
  # 
  # if (length(which(inside_poly))==0){return(integer())}
  # 
  # return(indies[which(inside_poly)])
}

reweight_pixels <- function(polygons_indexes, polygons_list, ODDy){
  #find which polygons have intersecting indexes but not intersecting sf polygons
  #this may be the case for pixels lying on the border which have been allocated to two polygons
  #in this case, allocate the pixel to the polygon in which its center lies
  
  coords <- xyFromCell(ODDy, 1:(NROW(ODDy)*NCOL(ODDy)))
  for (i in 1:length(polygons_indexes)){
    polygons_indexes[[i]]$weights <- rep(1, length(polygons_indexes[[i]]$indexes))
  }
  
  if (length(polygons_indexes) == 1) return(polygons_indexes)
  
  lon_cellsize <- res(ODDy)[1]
  lat_cellsize <- res(ODDy)[2]
  
  indexes_rast <- ODDy[['Population']]
  values(indexes_rast) <- 1:ncell(ODDy)
  pixels_as_poly <- as.polygons(indexes_rast,aggregate=F)
  
  for(i in 2:length(polygons_indexes)){
    for(j in 1:(i-1)){
      if (length(polygons_indexes[[i]]$indexes) == 0 | length(polygons_indexes[[j]]$indexes) == 0) next
      
      #If either is a country or total then skip
      if (!grepl(',', polygons_list[[i]]$polygon_name, fixed = TRUE) | !grepl(',', polygons_list[[j]]$polygon_name, fixed = TRUE)){
        next
      } 
      #If at different admin levels (e.g. one is admin level 2, one is admin level 1) then skip
      if (str_count(polygons_list[[i]]$polygon_name, ',') != str_count(polygons_list[[j]]$polygon_name, ',')){
        next
      } 
      
      
      #polygons_intersect <- gIntersects(polygons_list[[i]]$sf_polygon, polygons_list[[j]]$sf_polygon)
      #polygons_touch <- gTouches(polygons_list[[i]]$sf_polygon, polygons_list[[j]]$sf_polygon)
      indexes_intersect <- intersect(polygons_indexes[[i]]$indexes, polygons_indexes[[j]]$indexes)
      
      #if((!polygons_intersect | polygons_touch) & (length(indexes_intersect) > 0) ){
      if (length(indexes_intersect) > 0){
        print(paste('intersecting polys',polygons_indexes[[i]]$name,polygons_indexes[[j]]$name))
        for (index in indexes_intersect){
          #to allocate a weight to the pixel based on how much of the polygon it contains
          # pixel_sf <- sfheaders::sf_polygon(data.frame(longitude = c(coords[index,1]-lon_cellsize/2, 
          #                                                            coords[index,1]-lon_cellsize/2,
          #                                                            coords[index,1]+lon_cellsize/2, 
          #                                                            coords[index,1]+lon_cellsize/2, 
          #                                                            coords[index,1]-lon_cellsize/2), 
          #                                              latitude = c(coords[index,2]-lat_cellsize/2, 
          #                                                           coords[index,2]+lat_cellsize/2,
          #                                                           coords[index,2]+lat_cellsize/2, 
          #                                                           coords[index,2]-lat_cellsize/2,
          #                                                           coords[index,2]-lat_cellsize/2)), x='longitude',y='latitude')

          pixel <- pixels_as_poly[which(values(pixels_as_poly)==index),]
          
          #expanse(intersect(pixels_as_poly[index,], vect(polygons_list[[i]]$sf_polygon)))

          intersection_i <- expanse(intersect(pixel,vect(polygons_list[[i]]$sf_polygon)), unit='km')
          intersection_j <- expanse(intersect(pixel,vect(polygons_list[[j]]$sf_polygon)), unit='km')
          if (length(intersection_i)==0) intersection_i=0
          if (length(intersection_j)==0) intersection_j=0
          #combined_area = intersection_i + intersection_j
          polygons_indexes[[i]]$weights[which(polygons_indexes[[i]]$indexes==index)] <- intersection_i / expanse(pixel, unit='km')
          polygons_indexes[[j]]$weights[which(polygons_indexes[[j]]$indexes==index)] <- intersection_j / expanse(pixel, unit='km')
          #polygons_indexes[[i]]$weights[which(polygons_indexes[[i]]$indexes==index)] <- ifelse(is.null(intersection_i), 0, area(intersection_i)) / area(pixel_sp)
          #polygons_indexes[[j]]$weights[which(polygons_indexes[[j]]$indexes==index)] <- ifelse(is.null(intersection_j), 0, area(intersection_j)) / area(pixel_sp)
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

pixel_weight_border_correction <- function(ODDy, folder_write){
  # If a pixel has a total weight less than 1 across all polygons (for a specific administrative level), 
  # then either part of the pixel lies in the water or in a region for which we do not have impact data. 
  # We assume the former, and scale the weights so that they sum to one. 
  # Is there a way of knowing if it is the latter, and not rescaling? 
  pixels_of_interest <- which(!is.na(values(ODDy[['ISO3C']])))
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
      if(round(weight_sum[g],3) > 1.1){
        file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
        writeLines(paste("Region:", ODDy@polygons[[polygon_matches[[g]][1,1]]]$name, "Event Date:", ODDy@hazdates[1], '. Sum of pixel weights for a polygon is larger than 1'), file_conn)
        print(paste("Region:", paste(unlist(lapply(ODDy@polygons[polygon_matches[[g]][,1]], function(x) x$name)), collapse=' and '), "Event Date:", ODDy@hazdates[1],'. Sum of pixel weights for a polygon is', round(weight_sum[g],3), ' (larger than 1.1)'), file_conn)
        
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

addODDPolygons <- function(ODDy, polygons_list){
  # Populates the polygons slot in ODDy with a list, where each list element is another list containing:
  #   - polygon_name: the same as the polygon name chosen in getSubNatImpact
  #   - indexes: the indexes of the pixels in ODDy which belong to the polygon
  
  polygons_indexes <- list()
  
  for (i in 1:length(polygons_list)){
    #print(i)
    if(!grepl(',', polygons_list[[i]]$polygon_name, fixed = TRUE)){ #check if subnational by searching for comma in polygon name
      if (polygons_list[[i]]$polygon_name=='TOTAL'){ #if 'TOTAL' then add all pixels to the polygon
        inPolyInds <- which(!is.na(values(ODDy[['ISO3C']])))
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      } else { #if national then add pixels with appropriate ISO3C value to the polygon
        if (polygons_list[[i]]$polygon_name=='Kosovo'){ #kosovo doesn't seem to work with countrycode
          iso3='KOS'
        } else {
          iso3 <- countrycode(polygons_list[[i]]$polygon_name, origin='country.name', destination='iso3c')
        }
        inPolyInds <- which(values(ODDy[['ISO3C']]) == which(levels(ODDy[['ISO3C']])[[1]]$VALUE==iso3))
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      }
    } else {
      if(!is.null(polygons_list[[i]]$sf_polygon)){
        #inPolyInds <- intersect(which(inPoly((polygons_list[[i]]$sf_polygon)@polygons[[1]], pop = ODDy@coords)$indies), which(ODDy@data$ISO3C==iso3))
        inPolyInds <- intersect(pixelsInPoly((polygons_list[[i]]$sf_polygon)@polygons, ODDy), which(!is.na(values(ODDy[['ISO3C']]))))
        if(length(inPolyInds)==0){
          warning(paste('Region not found in impacted area. Polygon name: ', polygons_list[[i]]$polygon_name))
          file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
          #writeLines(paste("Region:", ODDy@polygons[[polygon_matches[[g]][1,1]]]$name, "Event Date:", ODDy@hazdates[1], '. Sum of pixel weights for a polygon is larger than 1'), file_conn)
          print(paste('Region not found in impacted area. Polygon name: ', polygons_list[[i]]$polygon_name), file_conn)
          close(file_conn)
        }
        polygons_indexes[[i]] <- list(name=polygons_list[[i]]$polygon_name, indexes = inPolyInds)
      }
    }
  }
  
  polygons_indexes <- reweight_pixels(polygons_indexes, polygons_list, ODDy)
  ODDy@polygons <- polygons_indexes
  
  return(ODDy)
  
}

addODDImpact <- function(ODDy, impact){
  #impact %<>% remove_double_counting(ODDy)
  ODDy@impact <- impact
  
  #if polygon is 'TOTAL' but only one country is exposed, can rename the 'TOTAL' polygon with the country
  unique_iso3 <- unique(ODDy$ISO3C)$ISO3C[!is.na(unique(ODDy$ISO3C)$ISO3C)]
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
    file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
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


updateODDSubNat <- function(dir, ODDy, event_sdate, event_fdate, event_id, folder_write, subnat_file='EQ_SubNational.xlsx'){
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
    print(paste('Event Date:', event_sdate, 'Countries:', paste(unique(ODDy$ISO3C)$ISO3C[!is.na(unique(ODDy$ISO3C)$ISO3C)], collapse=" ")))
    for(i in 1:length(SubNatData_match)){
      print(paste0('Option ', i,':'))
      print(SubNatData_match[[i]][1, c('event_name', 'iso3', 'sdate', 'fdate', 'notes')])
    }
    id_chosen <- as.integer(readline(prompt='id selected: '))
  }
  SubNatEvent <- SubNatData_match[[id_chosen]]
  
  iso3_unique <- unique(ODDy$ISO3C)$ISO3C[which(!is.na(unique(ODDy$ISO3C)$ISO3C))]
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
  
  ODDy <- pixel_weight_border_correction(ODDy, folder_write)
  
  #additional_poly_check(ODDy, event_id, print_to_xl=T)
  
  return(ODDy)
}


additional_poly_check <- function(ODDy, i, folder_write, print_to_xl=F){
  # Perform some checks on the polygons:
  #   - For each impact type, are all the observations at the same GADM level?
  #   - Are there any exposed pixels that haven't been allocated to a polygon? If so, what GADM region do they belong to.
  
  noteworthy <- F
  impacts_split <- ODDy@impact %>% group_by(impact) %>% group_split()
  
  if(print_to_xl){wb <- createWorkbook()}
  
  haz_max <- apply(as.data.frame(ODDy[[grep('hazMean', names(ODDy))]], na.rm=F), 1, function(row){ if(all(is.na(row))) return(0) else return(max(row, na.rm=T))})
  exposed_cells <-which(haz_max > 4.5)
  
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

        non_na_iso3 <-unique(ODDy$ISO3C[exposed_cells])$ISO3C[which(!is.na(unique(ODDy$ISO3C[exposed_cells])$ISO3C))]
        iso3_missing <- non_na_iso3[which(!(non_na_iso3 %in% iso3_incl))]
        for (iso3_miss in iso3_missing){
          max_exposed_int = max(haz_max[which(values(ODDy[['ISO3C']]==iso3_miss))], na.rm=T)
          writeData(wb, sheet, iso3_miss, startCol = 1, startRow = row)
          writeData(wb, sheet, max_exposed_int, startCol = 2, startRow = row)
          writeData(wb, sheet, length(which(values(ODDy[['ISO3C']]==iso3_miss)[exposed_cells])), startCol = 3, startRow = row)
          row = row + 1
          noteworthy=T
        }
        next
      }
      
      #check if there are any pixels not covered by the polygons:
      pixels_unallocated <- intersect(which(values(!is.na(ODDy[['ISO3C']]))), exposed_cells)
      for (k in 1:NROW(impacts_split[[j]])){
        if (gadm_levels[k] == gadm_level){
          pixels_in_poly <- ODDy@polygons[[impacts_split[[j]]$polygon[k]]]$indexes
          pixels_unallocated <- setdiff(pixels_unallocated, pixels_in_poly)
        }
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
          overlap <- (max(r_poly@bbox[1,]) >= ext(ODDy)[1] &
                        min(r_poly@bbox[1,]) <= ext(ODDy)[2] &
                        max(r_poly@bbox[2,]) >= ext(ODDy)[3] &
                        min(r_poly@bbox[2,]) <= ext(ODDy)[4])
          if (!overlap) next
          
          if (length(pixels_unallocated)==0) next
          
          sf_pixels_unallocated <- st_as_sf(data.frame(xyFromCell(ODDy, pixels_unallocated)), coords = c("x", "y"), crs = crs(ODDy))
          # spdf_pixels_unallocated <- SpatialPointsDataFrame(coords=ODDy@coords[pixels_unallocated,, drop=F],
          #                                                   data=ODDy@data[pixels_unallocated,1:2, drop=F], #this is arbitrary, function just seems to require data
          #                                                   proj4string=r_poly@proj4string)
          
          contained <- c(st_intersects(sf_pixels_unallocated, st_make_valid(st_as_sf(r_poly)), sparse=F)) #gContains(r_poly, spdf_pixels_unallocated, byid=T)
          if (any(contained)){
            print(paste('Polygon', r, 'contains an unallocated pixel'))
            max_exposed_int <- max(haz_max[pixels_unallocated[which(contained)]], na.rm=T)
            n_pixels <- length(which(contained))
            pixels_unallocated <- pixels_unallocated[-which(contained)]
            if(print_to_xl){ 
              r <- rev(r)
              writeData(wb, sheet, iso3, startCol = 1, startRow = row)
              for (ll in 1:length(r)){
                writeData(wb, sheet, r[ll], startCol = ll+1, startRow = row)
              }
              writeData(wb, sheet, max_exposed_int, startCol = ll+2, startRow = row)
              writeData(wb, sheet, n_pixels, startCol = ll+3, startRow = row)
              row = row + 1
              noteworthy=T
            }
          }
        }
      }
    }
  }
  
  if (print_to_xl & noteworthy) {
    file_path <- paste0(dir, folder_write, "MissingRegions2/Event_", i, ".xlsx")
    
    # Check if the file exists and delete it if so
    if (file.exists(file_path)) {
      file.remove(file_path)
    }
    
    # Save the new workbook
    saveWorkbook(wb, file_path)
  }
}

GetDataAll <- function(dir, haz="EQ", subnat_file= 'EQ_SubNational.xlsx', folder_write='IIDIPUS_Input_Alternatives/Aug24/'){
  # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if
  # no corresponding existing ODD object can be found, creates a new ODD object.
  
  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(EventID) %>% group_split()
  
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  options(timeout = 500)
  for (i in c(68)){#31,73,75,83,91,109,121,122)){#1:169){#c(89, 119, 122, 127, 133, 139, 150,151,152,164,165,166,167, 168,169)){#c(7,8,9,11,12,13,14,48,49,67,68,73,74,75,85,88,92,93,98,99,100,104,114, 128:163)){
    print(i)
    if (i==126) next
    # Subset displacement and disaster database objects to not all NA
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
      #stop()
      file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', lhazSDF not found.'), file_conn)
      close(file_conn) 
      next
    }
    
    # Create the ODD object:
    ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified, agg_level=1),error=function(e) NULL)
    if(is.null(ODDy)) {
      #stop()
      file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', ODD object not created.'), file_conn)
      close(file_conn) 
      next
    }
    
    
    #Fetch building count data:
    #ODDy_build <- tryCatch(AddBuildingCounts(ODDy, i, paste0(dir, folder_write, 'Building_count_notes')), error=function(e) NULL)
    ODDy_build <- tryCatch(getBingBuildingsGlobal(ODDy, i, paste0(dir, folder_write, 'Building_count_notes')), error=function(e) NULL)
    if(is.null(ODDy_build)) {
      #stop()
      file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', Building counts not added.'), file_conn)
      close(file_conn)
    } else {
      ODDy <- ODDy_build
      rm(ODDy_build)
    }
    
    iso3_ODDy <- unique(ODDy$ISO3C)$ISO3C
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDamSimplified$iso3)[which(unique(miniDamSimplified$iso3) !='TOT')][1],
                  "_",i)
    HAZARDpath<-paste0(dir,folder_write, "HAZARDobjects/",namer)
    saveHAZ(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    
    ODDy_with_impact <- tryCatch(updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1], i, folder_write),error=function(e) NULL)
    
    if(is.null(ODDy_with_impact)) {
      #stop()
      file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', impact not added.'), file_conn)
      close(file_conn)
      next
    }
    
    ODDy <- ODDy_with_impact
      
    rm(ODDy_with_impact)
    
    
    # Save out objects to save on RAM
    ODDpath<-paste0(dir, folder_write, "ODDobjects/",namer)
    saveODD(ODDy,ODDpath)
    
    additional_poly_check(ODDy, i, folder_write, print_to_xl=T)
    
    next
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDamSimplified$iso3) & 
                               sdate<mindate & sdate>maxdate)
    
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {
        file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
        writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', BD object creation failed.'), file_conn)
        close(file_conn)
      }
      BDpath <-paste0(dir, folder_write, "BDobjects/",namer)
      # Save it out!
      saveBD(BDy, BDpath)
      rm(BDy)
    }
    
    
    # Save some RAM
    rm(ODDy)
  }
  return(path)
}


UpdateBDData <- function(dir, haz="EQ", subnat_file= 'EQ_SubNational.xlsx', folder_write='IIDIPUS_Input_Alternatives/Dec24/'){
  # Works through EQ_Subnational.xlsx and, for each event, either updates the existing ODD object or, if
  # no corresponding existing ODD object can be found, creates a new ODD object.
  
  #SubNatData <- read.xlsx(paste0(dir, 'IIDIPUS_Input/', subnat_file), colNames = TRUE , na.strings = c("","NA"))
  SubNatData <- readSubNatData(subnat_file)
  
  # Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
  SubNatDataByEvent <- SubNatData %>% group_by(EventID) %>% group_split()
  
  # Extract all building damage points
  #Damage<-ExtractBDfiles(dir = dir,haz = haz)
  
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  options(timeout = 500)
  
  ufiles<-na.omit(list.files(path='IIDIPUS_Input_Alternatives/Aug24/ODDobjects/',pattern=Model$haz,recursive = T,ignore.case = T))
  
  for (i in 1:170){#31,73,75,83,91,109,121,122)){#1:169){#c(89, 119, 122, 127, 133, 139, 150,151,152,164,165,166,167, 168,169)){#c(7,8,9,11,12,13,14,48,49,67,68,73,74,75,85,88,92,93,98,99,100,104,114, 128:163)){
    print(i)
    if (i==126) next
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
    
    namer <- grep(paste0("_", i, "$"), ufiles, value = TRUE)
    #Fetch building count data:
    #ODDy_build <- tryCatch(AddBuildingCounts(ODDy, i, paste0(dir, folder_write, 'Building_count_notes')), error=function(e) NULL)
    ODDy <- readODD(paste0('IIDIPUS_Input_Alternatives/Aug24/ODDobjects/', namer))
    #ODDy@impact <- NULL
    #ODDy@polygons <- NULL
    
    ODDy_with_impact <- tryCatch(updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1], i, folder_write),error=function(e) NULL)
    
    if(is.null(ODDy_with_impact)) {
      #stop()
      file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
      writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', impact not added.'), file_conn)
      close(file_conn)
      next
    }
    
    ODDy <- ODDy_with_impact
    
    rm(ODDy_with_impact)
    
    
    # Save out objects to save on RAM
    ODDpath<-paste0(dir, folder_write, "ODDobjects/",namer)
    saveODD(ODDy,ODDpath)
    
    additional_poly_check(ODDy, i, folder_write, print_to_xl=T)
    
    next
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDamSimplified$iso3) & 
                               sdate<mindate & sdate>maxdate)
    
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {
        file_conn <- file(paste0(dir, folder_write, 'ODD_creation_notes'), open = "a")
        writeLines(paste("Index:", i, "Event Name:", SubNatDataByEvent[[i]]$event_name[1], "Event Date:", SubNatDataByEvent[[i]]$sdate[1], ', BD object creation failed.'), file_conn)
        close(file_conn)
      }
      BDpath <-paste0(dir, folder_write, "BDobjects/",namer)
      # Save it out!
      saveBD(BDy, BDpath)
      rm(BDy)
    }
    
    
    # Save some RAM
    rm(ODDy)
  }
  return(path)
}

# additional_poly_check_all <- function(input_folder='IIDIPUS_Input_Alternatives/Aug24/ODDobjects/', folder_write='IIDIPUS_Input_Alternatives/Aug24/'){
#   ufiles<-list.files(path=paste0(dir, input_folder),pattern=Model$haz,recursive = T,ignore.case = T)
#   for (file in ufiles){
#      print(file)
#      i <- as.numeric(regmatches(file, gregexpr("[0-9]+", file))[[1]][2])
#      if (!(i %in% c(1,16,22,68, 70, 89, 127, 138, 149))){
#      #if (!(i %in% c(16))){
#        next
#      }
#      ODDy <- readODD(paste0(dir,input_folder, file))
#      #ODDy_with_impact <- tryCatch(updateODDSubNat(dir, ODDy, miniDam$sdate[1], miniDam$fdate[1], i,folder_write),error=function(e) NULL)
#      #if (is.null(ODDy_with_impact)) stop()
#      
#      #ODDy <- ODDy_with_impact
#      
#      additional_poly_check(ODDy, i, folder_write, print_to_xl=T)
#      
#      #saveODD(ODDy, paste0(dir, input_folder, file, '_MAR'))
#   }
# }



# fill_missing_EQFreq <- function(ODD){
#   ODD$EQFreq <- focal(ODD$EQFreq, fun="modal", na.policy="only")
#   missing <- which(is.na(values(ODD$EQFreq)) & !is.na(values(ODD$Population)))
#   if (length(missing) > 0){
#     ODD$indexes <- 1:ncell(ODD)
#     dists <- terra::distance(crds(ODD, na.rm=F)[missing,, drop=F], crds(ODD, na.rm=F)[!is.na(values(ODD$EQFreq)),])
#     ODD$EQFreq[missing] = ODD$EQFreq[ODD$indexes[!is.na(values(ODD$EQFreq))][apply(dists,1, which.min),1]]
#     ODD$indexes <- NULL
#   }
#   return(ODD)
# }
# 
# input_folder <- 'IIDIPUS_Input_Alternatives/Aug24Agg/ODDobjects/'
# output_folder <- 'IIDIPUS_Input_Alternatives/Aug24Agg/ODDobjects_EQFreqFull/'
# fill_missing_EQFreq_all <- function(input_folder, output_folder){
#   ufiles<-list.files(path=paste0(dir, input_folder),pattern=Model$haz,recursive = T,ignore.case = T)
#   for (file in ufiles){
#     ODDy <- readODD(paste0(dir,input_folder, file))
#     ODDy %<>% fill_missing_EQFreq
#     saveODD(ODDy, paste0(dir,output_folder, file))
#   }
# }




moveTestData <- function(folder_in='IIDIPUS_Input_Alternatives/Nov24Agg'){
    ODD_folderall<-paste0(dir, folder_in, '/ODDobjects/')
    ODD_foldertest<-paste0(dir, folder_in, '/ODDobjects/Test/')
    ufiles<-list.files(path=ODD_folderall,pattern=Model$haz,recursive = T,ignore.case = T)
    #set.seed(1)
    #ufiles <- ufiles[sample(1:length(ufiles), length(ufiles), replace=F)]
    total_mortalities <- c()
    #sort by mortality
    for (file in ufiles){
      ODD <- readODD(paste0(ODD_folderall, file))
      total_mortalities <- c(total_mortalities, max(values(ODD[[grep('hazMean', names(ODD))]])[values(!is.na(ODD$Population) & ODD$Population > 0),], na.rm=t))
      # polygon_names <- unlist(lapply(ODD@polygons[ODD@impact$polygon], function(x) x$name))
      # if (any(tolower(polygon_names[which(ODD@impact$impact=='mortality')]) %in% c('tot', 'total'))){
      #   nonmatch <- which(!tolower(polygon_names[which(ODD@impact$impact=='mortality')]) %in% c('tot', 'total'))
      #   if (length(nonmatch)>0){
      #     ODD@impact <- ODD@impact[-which(ODD@impact$impact=='mortality')[nonmatch],] # in the case of total and subnational data, remove the subnational
      #   }
      # }
      # total_mortality <- sum(ODD@impact$observed[which(ODD@impact$impact=='mortality')])
      # total_mortalities <- c(total_mortalities, total_mortality)
    }
    ufiles <- ufiles[order(total_mortalities, decreasing=T)] # sort by date

    for (i in 1:length(ufiles)){
      file <- ufiles[i]
      if (i %%3 != 1){next}
      options(warn = 2)
      file.copy(from = paste0(ODD_folderall, file),
                to = paste0(ODD_foldertest, file))
      file.remove(from = paste0(ODD_folderall, file))
      options(warn = 1)
    }
  # BD_folderall<-paste0(dir, folder_in, '/BDobjects/')
  # BD_foldertest<-paste0(dir, folder_in, '/BDobjects/Test/')
  # ufiles<-list.files(path=BD_folderall,pattern=Model$haz,recursive = T,ignore.case = T)
  # i <- 0
  # for (file in ufiles){
  #   i <- i + 1
  #   if (i %%3 != 0){next}
  #   file.copy(from = paste0(BD_folderall, file),
  #             to = paste0(BD_foldertest, file))
  #   file.remove(from = paste0(BD_folderall, file))
  # }
}


# --------------------------------------------------------------------------------
#-------- For updating ODD objects rather than building them from scratch: -------
# --------------------------------------------------------------------------------


updateAllODDSubNat <- function(dir, folder_write='IIDIPUS_Input/', subnat_file='EQ_SubNational.xlsx'){
  # Takes all the ODD objects currently in IIDIPUS_Input/ODDobjects and updates them using the data from the subnational data spreadsheet
  
  folderin<- paste0(dir, 'IIDIPUS_Input/ODDobjects/') #"/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/"
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
  for (ODD_file in ufiles){
    event_date <- as.Date(substr(ODD_file, 3, 10), "%Y%m%d") #seems to be some issues with ODD@hazdate
    ODDy <- readRDS(paste0(folderin, ODD_file))
    ODDy <- updateODDSubNat(dir, ODDy, event_date, folder_write, subnat_file)
    #saveRDS(ODDy, paste0(folderin, ufiles[i])) 
  }
}


