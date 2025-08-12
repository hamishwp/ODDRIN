#install.packages(c('data.table', 'jsonlite', 'geojsonio'))
# library(data.table)
# library(jsonlite)
# library(geojsonio)

merge_rastered_spdf <- function(raster1_spdf, raster2_spdf, added_var_name){
  data <- raster1_spdf@data
  data$Longitude <-  round(raster1_spdf@coords[,1], 8)
  data$Latitude <- round(raster1_spdf@coords[,2], 8)
  data$id <- 1:NROW(data)
  data %<>% merge(data.frame(Longitude=round(raster2_spdf@coords[,1], 8), Latitude=round(raster2_spdf@coords[,2], 8), 
                             added_var = pull(raster2_spdf@data)), 
                  by=c('Latitude', 'Longitude'), all.x = TRUE)
  colnames(data)[which(colnames(data)=='added_var')] <- added_var_name
  data <- data[order(data$id),]
  data <- data[-which(colnames(data) %in% c('id', 'Latitude', 'Longitude'))]
  return(data)
}


# wkt polygons for events:
# i = 44: POLYGON((30.5 -2, 30.5 0.5, 32.5 0.5, 32.5 -2, 30.5 -2))
# i = 82: POLYGON((34.8 -17.5, 34.8 -16.2, 36.1 -16.2, 36.1 -17.5, 34.8 -17.5))
# i = 146: POLYGON((5.9 36.1,5.9 36.9,6.8 36.9,6.8 36.1, 5.9 36.1))
AddOpenBuildingCounts <- function(ODD, isos_openbuildings, i, national_coverage=T){
  lon_min <- ODD@hazinfo$bbox[1]; lon_max <- ODD@hazinfo$bbox[3]
  lat_min <- ODD@hazinfo$bbox[2]; lat_max <- ODD@hazinfo$bbox[4]
  
  indies_open_buildings <- which(levels(ODD$ISO3C)[[1]]$VALUE[values(ODD$ISO3C)] %in% isos_openbuildings) # need to check
  
  iso3_unique <- levels(ODD$ISO3C)[[1]]$VALUE
  
  if (national_coverage){
    open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/open_buildings_v2_points_your_own_wkt_polygon_',iso3_unique[1],'.csv')
    if(!file.exists(open_buildings_file)){
      stop('Download Open Buildings data from: https://colab.research.google.com/github/google-research/google-research/blob/master/building_detection/open_buildings_download_region_polygons.ipynb')
    }
    #When downloading, insert the polygon bounding the country of interest into your_own_wkt_polygon field
    # alternatively, can select the region_border_source and country, although this doesn't seem to work well when working with small islands in Philippines (and potentially other small regions)
    # select 'points' in the data_type field
  } else {
    open_buildings_files <- paste0(dir, 'Demography_Data/Buildings/open_buildings_v2_points_your_own_wkt_polygon_i',i,'.csv')
    if(!any(file.exists(open_buildings_files))){
      stop('Download Open Buildings data from: https://colab.research.google.com/github/google-research/google-research/blob/master/building_detection/open_buildings_download_region_polygons.ipynb')
    }
    open_buildings_file <- open_buildings_files[which(file.exists(open_buildings_files))]
  }

  # doesn't work once files exceed a certain size:
  # building_locs <- read.csv.sql(open_buildings_file,
  #                               paste0("select longitude, latitude from file where latitude > ", lat_min, ' AND longitude > ', lon_min, 
  #                                      ' AND latitude < ',lat_max, ' AND longitude < ', lon_max))
  
  building_locs <- data.frame(latitude=double(), longitude=double())
  i <- 1
  nrows_file <- as.integer(strsplit(system(paste0('wc -l ', open_buildings_file), intern=T), ' ')[[1]][1])
  nrow_tmp <- 50000000
  nchunks <- ceiling(nrows_file/nrow_tmp)
  for (j in 1:nchunks){
    tmp <- fread(open_buildings_file,skip=i, nrows=nrow_tmp, select=c(1,2), col.names=c('latitude', 'longitude'))
    building_locs %<>% rbind(tmp %>% filter(latitude > lat_min, longitude > lon_min, latitude < lat_max, longitude < lon_max)) 
    i <- i + nrow_tmp
  }
  
  if (length(iso3_unique) > 1 & national_coverage==T){
    for (i in 2:length(iso3_unique)){
      open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/open_buildings_v2_points_your_own_wkt_polygon_',iso3_unique[i],'.csv')
      if(!file.exists(open_buildings_file)){
        stop('Download Open Buildings data from: https://colab.research.google.com/github/google-research/google-research/blob/master/building_detection/open_buildings_download_region_polygons.ipynb')
        #follow download instructions above
      }
      i <- 1
      nrows_file <- as.integer(strsplit(system(paste0('wc -l ', open_buildings_file), intern=T), ' ')[[1]][1])
      nrow_tmp <- 50000000
      nchunks <- ceiling(nrows_file/nrow_tmp)
      for (j in 1:nchunks){
        tmp <- fread(open_buildings_file,skip=i, nrows=nrow_tmp, select=c(1,2), col.names=c('latitude', 'longitude'))
        building_locs %<>% rbind(tmp %>% filter(latitude > lat_min, longitude > lon_min, latitude < lat_max, longitude < lon_max))
        i <- i + nrow_tmp
      }
    }
  }
  
  building_locs <- building_locs[,c('longitude', 'latitude')]
  
  rastered_buildings <- rasterize(building_locs, raster(ODD), 1, fun='count')
  
  ODD[['OpenBuildings']] = rast(rastered_buildings)
  # rastered_buildings_spdf <- as(rastered_buildings, "SpatialPixelsDataFrame")
  # 
  # ODD@data <- merge_rastered_spdf(ODD, rastered_buildings_spdf, 'nBuildings')
  # if (national_coverage){
  #   ODD$nBuildings[indies_open_buildings[which(is.na(ODD$nBuildings[indies_open_buildings]))]] <- 0
  #   ODD$nBuildings[-indies_open_buildings] <- NA
  # } else {
  #   ODD$nBuildings[which(is.na(ODD$nBuildings))] <- 0
  # }
  # 
  # sedacs2020 <- GetPopulationBbox(dir, ODD@bbox, yr=2020)
  # population2020 <- merge(ODD@coords, cbind(sedacs2020@coords, population2020=sedacs2020$Population), by=c('Longitude', 'Latitude'), all.x=T, sort=F)
  # nonzero_pop <- which(population2020$population2020 > 0)
  # ODD$nBuildings[nonzero_pop] <- round(ODD$nBuildings[nonzero_pop] * (ODD$Population[nonzero_pop] / population2020$population2020[nonzero_pop]))
  # 
  return(ODD)
}

AddBuildingCounts <- function(ODD, i, file_write='IIDIPUS_Input/Building_count_notes'){
  isos_openbuildings <- c('IDN', 'PHL') 
  events_openbuildings <- c(34, 44, 82, 131, 146)
  #have checked for coverage
  #events_bingbuildings <- c(5, 39, 40, 53, 95, 142) #not covered: 15 #partially covered: 42, 50, 78
  
  iso3_unique <- unique(ODD$ISO3C)[!is.na(unique(ODD$ISO3C))]
  
  if (i %in% events_openbuildings){
    ODD %<>% AddOpenBuildingCounts(isos_openbuildings, i, national_coverage=F)
  } else if (any(iso3_unique %in% isos_openbuildings)){
    ODD %<>% AddOpenBuildingCounts(isos_openbuildings, i, national_coverage=T)
  } else if (i %in% events_bingbuildings){#(any(iso3_unique %in% isos_bingbuildings)){
    ODD %<>% AddBingBuildingCounts()
    # file_conn8 <- file('ODD_creation_notes/DoBuildingCountsManually', open = "a")
    # writeLines(paste("Bing Build Count Missing", paste(iso3_unique, sep='_'), "Event Date:", ODD@hazdates[1]), file_conn8)
    # close(file_conn8)
    #stop('Double check ODD object to ensure that Bing Building Footprints provides full coverage')
  } else {
    ODD %<>% getBingBuildingsGlobal()
  }
  if (is.null(ODD$nBuildings)){
    file_conn <- file(file_write, open = "a")
    writeLines(paste("Event", i, ": No building counts"), file_conn)
    close(file_conn) 
    return(ODD)
  } 
  missing_building_counts <- which(!is.na(values(ODD$ISO3C)) & is.na(values(ODD$nBuildings)))
  if (length(missing_building_counts)>0){
    file_conn <- file(file_write, open = "a")
    writeLines(paste("Build Count Missing for pixels in countries", paste(unique(ODD$ISO3C[missing_building_counts]), sep='_'), "Event Date:", ODD@hazdates[1]), file_conn)
    close(file_conn) 
  }
  return(ODD)
}

ReplaceBuildingCounts <- function(){
  folderin<-"/home/manderso/Documents/GitHub/ODDRIN/PHL_IIDIPUS_Input/ODDobjects/"
  ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  for (i in 1:length(ufiles)){
    print(i)
    ODDy<-readRDS(paste0(folderin,ufiles[i]))
    ODDy@data <- ODDy@data[,-which(colnames(ODDy@data)=='nBuildings')]
    ODDy %<>% AddBuildingCounts()
    saveRDS(ODDy, paste0(folderin,ufiles[i]))
    rm(ODDy)
  }
}

resave_geojsonl <- function(iso3, USA_state=NULL){
  #reads in geojsonl file for a country and saves a simplified RDS
  if (iso3 != 'USA'){
    open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/',iso3,'.geojsonl')
    if (!file.exists(open_buildings_file)) stop('Download Bing Building Footprint data from: https://www.microsoft.com/en-us/maps/building-footprints')
    tmp_df <- new.env()
    stream_in(file(open_buildings_file), handler = function(df){
      idx <- as.character(length(tmp_df) + 1)
      tmp_df[[idx]] <- sapply(df$coordinates, function(poly_coords) poly_coords[1,1,])
      #just take the first coordinate of the polygon 
      # might not be quite so accurate but shouldn't make too much difference and does speed up the process
    })
    building_footprints <- t(do.call(cbind, as.list(tmp_df)))
    colnames(building_footprints) <- c('longitude', 'latitude')
    saveRDS(building_footprints, paste0(dir, 'Demography_Data/Buildings/',iso3,'.RDS'))
  } else {
    open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/',USA_state,'.geojson')
    spdf <- geojson_read(open_buildings_file,  what = "sp")
    building_footprints <- t(sapply(spdf@polygons, function(poly) poly@labpt))
    building_footprints %<>% data.frame() 
    colnames(building_footprints) <- c('longitude', 'latitude')
    saveRDS(building_footprints, paste0(dir, 'Demography_Data/Buildings/',USA_state,'.RDS'))
  }
}

AddBingBuildingCounts <- function(ODD, plot_only = F){
  isos_bingbuildings <- c('COL', 'ECU', 'USA', 'PER')
  
  lon_min <- ODD@bbox[1,1]; lon_max <- ODD@bbox[1,2]
  lat_min <- ODD@bbox[2,1]; lat_max <- ODD@bbox[2,2]
  
  iso3_unique <- unique(ODD$ISO3C)[!is.na(unique(ODD$ISO3C))]
  
  if (any(!iso3_unique %in% isos_bingbuildings) ){
    file_conn <- file('IIDIPUS_Input_June20/ODD_creation_notes', open = "a")
    writeLines(paste("Bing buildings missing for a country. Event countries:", paste(unique(ODD$ISO3C), sep='_'), "Event Date:", ODD@hazdates[1]), file_conn)
    close(file_conn) 
  }
  
  indies_bing_buildings <- which(ODD$ISO3C %in% isos_bingbuildings)
  
  open_buildings_files <- c()
  for (iso3 in intersect(iso3_unique, isos_bingbuildings)){
    if (iso3 !='USA'){
      open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/',countrycode(iso3, origin='iso3c', destination='country.name'),'.RDS')
      if (!file.exists(open_buildings_file)){
        resave_geojsonl(iso3)
      }
    } else {
      if (lon_min > -114.06 & lon_max < -109.04 & lat_min > 37 & lat_max < 42){
        open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/Utah.RDS')
      } else {
        stop('So far only have building count data matching one event in Utah')
      }
    }
    open_buildings_files <- append(open_buildings_files, open_buildings_file)
  }
  
  building_locs <- data.frame(latitude=double(), longitude=double())
  
  for(fff in open_buildings_files){
    building_footprints <- data.frame(readRDS(fff))
  
    nrow_tmp <- 500000
    nrow_bf <- NROW(building_footprints)
    nchunks <- ceiling(nrow_bf/nrow_tmp)
    for (j in 1:nchunks){
      #tmp <- t(sapply(building_footprints$coordinates[((j-1)*nrow_tmp+1):min(j*nrow_tmp,nrow_bf)], function(poly_coords) poly_coords[1,1,]))
      #colnames(tmp) <- c('longitude', 'latitude')
      tmp = building_footprints[((j-1)*nrow_tmp+1):min(j*nrow_tmp,nrow_bf),]
      building_locs %<>% rbind(tmp %>% filter(latitude > lat_min, longitude > lon_min, latitude < lat_max, longitude < lon_max)) 
    }
  }
  
  if (plot_only){ 
    #plot convex hull of building_locs and compare to region of ODD object, to see if there is reasonable coverage
    plot(ODD@coords[which(ODD$hazMean1>4 & !is.na(ODD$ISO3C)),])
    hull <- chull(building_locs)
    lines(rbind(building_locs[hull,], building_locs[hull[1],]), col='red')
    return(NULL)
  }
  
  rastered_buildings <- rasterize(building_locs, raster(ODD), 1, fun='count')
  rastered_buildings_spdf <- as(rastered_buildings, "SpatialPixelsDataFrame")
  
  ODD@data <- merge_rastered_spdf(ODD, rastered_buildings_spdf, 'nBuildings')

  ODD$nBuildings[indies_bing_buildings[which(is.na(ODD$nBuildings[indies_bing_buildings]))]] <- 0
  ODD$nBuildings[-indies_bing_buildings] <- NA
  
  sedacs2020 <- GetPopulationBbox(dir, ODD@bbox, yr=2020)
  population2020 <- merge(ODD@coords, cbind(sedacs2020@coords, population2020=sedacs2020$Population), by=c('Longitude', 'Latitude'), all.x=T, sort=F)
  nonzero_pop <- which(population2020$population2020 > 0)
  ODD$nBuildings[nonzero_pop] <- round(ODD$nBuildings[nonzero_pop] * (ODD$Population[nonzero_pop] / population2020$population2020[nonzero_pop]))
  
  return(ODD)
}


# filename <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NoBuildingDat/ODDobjects/EQ20150425NPL_14'
# ODDy <- readRDS(filename)
# ODDy %<>% AddBuildingCounts()
# # #
# cr <- colorRamp(c("green", "black"))
# plot(ODDy@coords[which(!is.na(ODDy$Population)),], col=rgb(cr(ODDy$nBuildings[which(!is.na(ODDy$Population))]/max(ODDy$nBuildings, na.rm=T)), max=255))
# points(ODDy@coords[which(ODDy$nBuildings==0),], col='red')
# #
# saveRDS(ODDy, filename)














AddFBBuildingEstimates <- function(ODD, isos_osm){
  indies_osm <- which(ODD$ISO3C %in% isos_osm)

  bbox <- ODD@bbox 
  #resize bbox to avoid loading unnecessary FB data
  bbox[1,1] <- min(ODD@coords[indies_osm,1]); bbox[1,2] <- max(ODD@coords[indies_osm,1])
  bbox[2,1] <- min(ODD@coords[indies_osm,2]); bbox[2,2] <- max(ODD@coords[indies_osm,2])
  
  FBpop<-readFBpop(bbox) # Get Population estimates from Meta Data for Good
  rastered_FBpop_builtup <- rasterize(FBpop, raster(ODD), 'Population', fun='count')
  data_with_builtup <- merge_rastered_spdf(SEDAC, as(rastered_FBpop_builtup, "SpatialPixelsDataFrame"), 'nBuiltup')
  ODD@data$nBuiltup <- NA
  ODD@data[indies_osm] <- ODD@data_with_builtup[indies_osm,]
  return(ODD)
}


sampleBuildingsFromBuiltup <- function(nBuiltup, iso3c, reg_coefs, Np){
  
  n_pixels <- NROW(nBuildings_sample)
  reg_coefs <- merge(iso3c, reg_coefs, all.x <- T, sort = F)
  
  nBuildings_sample <- matrix(rlnorm(Np * n_pixels, rep(reg_coefs[,2] + reg_coefs[,3] * log(nBuiltup), Np), rep(reg_coefs[,4], Np) ), nrow=n_pixels)
  
  return(nBuildings_sample)
}

getBuildSamplingCoefs <- function(filename){
  nBuildSamplingCoefs <- readRDS(paste0(dir, filename))
  return(nSamplingModelCoefs)
}

addRegCoefs <- function(iso3, bbox_list, filename){
  nBuildSamplingCoefs <- readRDS(paste0(dir, filename))
  
  #read in OSM data for the polygon
  #lognormal regression against nBuiltup
  #save output to file 
  
  train_dat <- data.frame(iso3 = character(), nBuiltup = integer(), nBuildings_OSM = integer(), bbox_i = integer())
  for(i in 1:length(bbox_list)){
    SEDACS_grid <- GetPopulationBbox(dir, bbox_list[[i]])
    #make sure all pixels are fully inside bbox
    SEDACS_grid <- SEDACS_grid[-which(SEDACS_grid@coords[,1] - SEDACS_grid@grid@cellsize/2 < bbox_list[[i]][1,1] |
                       SEDACS_grid@coords[,1] + SEDACS_grid@grid@cellsize/2 > bbox_list[[i]][1,2] |                                               
                       SEDACS_grid@coords[,2] - SEDACS_grid@grid@cellsize/2 < bbox_list[[i]][2,1] |     
                       SEDACS_grid@coords[,2] + SEDACS_grid@grid@cellsize/2 > bbox_list[[i]][2,2]),] 
    
    FBpop <-readFBpop( bbox_list[[i]])
    rastered_FBpop_builtup <- rasterize(FBpop, raster(SEDACS_grid), 'Population', fun='count') 
    FB_SEDACS_merged <- merge_rastered_spdf(SEDACS_grid, as(rastered_FBpop_builtup, "SpatialPixelsDataFrame"), 'nBuiltup')
    
    OSM_buildings<-GetOSMbuildingsBbox(bbox_list[[i]], timeout=60) 
    rastered_OSM_count <- rasterize(cbind(OSM_buildings$Longitude, OSM_buildings$Latitude), raster(SEDACS_grid), fun='count') # Get OSM building counts
    OSM_SEDACS_merged <- merge_rastered_spdf(SEDACS_grid, as(rastered_OSM_count, "SpatialPixelsDataFrame"), 'nBuildings_OSM')
    
    #SEDACS_grid %<>% AddBuildingCounts()
    train_dat %<>% rbind(data.frame(iso3=SEDACS_grid$ISO3C, nBuiltup=FB_SEDACS_merged$nBuiltup, nBuildings_OSM=OSM_SEDACS_merged$nBuildings_OSM, bbox_i = i))
    
  }
  
  lnorm_fit <- lm(log(nBuildings_OSM) ~ log(nBuiltup), data=train_dat)
  
  iso3_coefs <- c(iso3, lnorm_fit$coefficients[1], lnorm_fit$coefficients[2], sigma(lnorm_fit))
  
  match_index <- which(iso3==nBuildSamplingCoefs$iso3c)
  if (length(match_index) > 0){
    print(paste('Updating Building Count Sampling Coefficients for ', iso3))
    nBuildSamplingCoefs[match_index,] <- iso3_coefs
  } else {
    nBuildSamplingCoefs %<>% rbind(iso3_coefs)
  }
  saveRDS(nBuildSamplingCoefs, paste0(dir, filename))
  
  return(nBuildSamplingCoefs)
}


# ----------------- Bing building footprints global --------------------------


# ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_June24_AggFactor5/ODDobjects/EQ20131015PHL_16')

getMercantileTile <- function(latitude, longitude, zoom){
  sinLatitude = sin(latitude * pi/180)
  pixelX = ((longitude + 180) / 360) * 256 * (2^zoom)
  pixelY = (0.5 - log((1 + sinLatitude) / (1 - sinLatitude)) / (4 * pi)) * 256 * (2^zoom)
  tileX <- floor(pixelX/256)
  tileY <- floor(pixelY/256)
  return(c(tileX, tileY, zoom))
}

getMercantileTiles <- function(west, south, east, north, zoom){
  minTile <- getMercantileTile(north, west, zoom)
  maxTile <- getMercantileTile(south, east, zoom)
  tiles_contained <- expand.grid(x=minTile[1]:maxTile[1], y=minTile[2]:maxTile[2])
  tiles_contained <- add_column(tiles_contained, z=zoom)
  return(tiles_contained)
}

getQuadKey <- function(Tile) {
  tileX <- as.integer(Tile[1])
  tileY <- as.integer(Tile[2])
  zoom <- as.integer(Tile[3])
  binX <- rev(as.integer(intToBits(tileX))[1:zoom])
  binY <- rev(as.integer(intToBits(tileY))[1:zoom])
  
  # Alternate the letters from each string
  #combined[seq(1, 2*zoom, by = 2)] <- substr(binY, 1, zoom)
  #combined[seq(2, 2*zoom, by = 2)] <- substr(binX, 1, zoom)
  
  result <- character(zoom)
  for(i in 1:length(result)){
    result[i] <- 2*binY[i] + binX[i]
  }
  
  return(as.integer(paste(result, collapse='')))
}

getBingBuildingsFromTiles <- function(bbox, event_id=NA, file_write='IIDIPUS_Input/Building_count_notes'){
  tiles <- getMercantileTiles(bbox[1], bbox[2], bbox[3], bbox[4], 9)
  quad_keys <- list()
  for (i in 1:NROW(tiles)){
    quad_keys[[i]] <- getQuadKey(tiles[i,])
  }
  
  cat(sprintf("The input area spans %d tiles: %s\n", length(quad_keys), paste(quad_keys, collapse = ", ")))
  
  df <- read.csv("https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv")
  
  # Download the GeoJSON files for each tile that intersects the input geometry
  build_coords <- array(dim=c(0,2))
  
  missing_quadkeys <- unlist(quad_keys)[!unlist(quad_keys) %in% df$QuadKey]
  # if (length(missing_quadkeys)> (0.2 * length(quad_keys))){
  #   file_conn <- file(file_write, open = "a")
  #   writeLines(paste("Event:", event_id, ", Missing", length(missing_quadkeys)/length(quad_keys)*100, "percent of quad keys, not adding building data."), file_conn)
  #   close(file_conn) 
  #   return(ODD)
  # } else if (length(missing_quadkeys)> 0){
  #   file_conn <- file(file_write, open = "a")
  #   writeLines(paste("Event:", event_id, ", Missing", length(missing_quadkeys)/length(quad_keys)*100, "percent of quad keys, but still adding building data."), file_conn)
  #   close(file_conn) 
  # }
  if (length(missing_quadkeys)> 0){
    if (is.null(file_write)){
      print(paste("Event:", event_id, ", Missing", length(missing_quadkeys)/length(quad_keys)*100, "percent of quad keys, not adding building data."))
    } else {
      file_conn <- file(file_write, open = "a")
      writeLines(paste("Event:", event_id, ", Missing", length(missing_quadkeys)/length(quad_keys)*100, "percent of quad keys, not adding building data."), file_conn)
      close(file_conn)
    }
    return(NULL)
  }
  
  missing_quadkeys_flag <- F
  
  for (quad_key in quad_keys) {
    rows <- df[df$QuadKey == quad_key, ]
    if (nrow(rows) >= 1) {
      for (url in rows$Url){
        tmp <- tempfile()
        #download.file(url, destfile =tmp,quiet = FALSE, mode = "wb")
        #out <- lapply(readLines(tmp), fromJSON)
        tryCatch({
          download.file(url, destfile = tmp, quiet = FALSE, mode = "wb")
          out <- lapply(readLines(tmp), fromJSON)
        }, error = function(e) {
          warning("An error occurred: ", e$message)
          out <- NULL  # or some fallback value
        })
        
        #convert building polygons to points by just taking the first in the coordinate 
        build_coords %<>% rbind(t(sapply(out, function(build) return(build$geometry$coordinates[1,1,]))))
      }
    } else {
      missing_quadkeys_flag <- T
      # file_conn <- file('IIDIPUS_Input_NonFinal/IIDIPUS_Input_June24_AggFactor5/BingGlobalNotes', open = "a")
      # writeLines(paste("               Event Date:", ODD@hazdates[1], ". No building data for quadkey", quad_key), file_conn)
      # close(file_conn) 
      # print(paste("QuadKey not found in dataset:", quad_key))
      #stop(paste("QuadKey not found in dataset:", quad_key))
    }
  }
  inside_bbox <- which(build_coords[,1] > bbox[1] & build_coords[,1] < bbox[3] & build_coords[,2] > bbox[2] & build_coords[,2] < bbox[4])
  building_locs <- build_coords[inside_bbox,]
  
  if (!missing_quadkeys_flag){
    if (is.null(file_write)){
      print(paste("Event:", event_id, ", Complete Building Count from Global Bing Building Footprints"))
    }
    else {
      file_conn <- file(file_write, open = "a")
      writeLines(paste("Event:", event_id, ", Complete Building Count from Global Bing Building Footprints"), file_conn)
      close(file_conn)
    }
  }
  
  colnames(building_locs) <- c('Longitude', 'Latitude')
  
  return(building_locs)
}

getBingBuildingsGlobal <- function(ODD, event_id, file_write='IIDIPUS_Input/Building_count_notes', aggregate=T){
  # Note that if you choose to proceed with some missing quadkeys, the values in these pixels will be 0 rather than NA (making them indistinguishable from true 0s)
  # This could be addressed by determining the polygons of the missing quadkeys and assigning the intersecting pixels to NA, but haven't had time to implement this yet. 
  
  zoom <- 9
  
  building_locs =  getBingBuildingsFromTiles(c(ext(ODD)[1], ext(ODD)[3], ext(ODD)[2], ext(ODD)[4]), event_id=NA, file_write='IIDIPUS_Input/Building_count_notes')
  
  if(is.null(building_locs)) {return(ODD)}
  
  if (!aggregate){return(building_locs)}
                                           
  ODD[['nBuildings']] <- rasterize(building_locs, ODD, 1, fun='count', background=0)
  #rastered_buildings_spdf <- as(rastered_buildings, "SpatialPixelsDataFrame")
  
  #ODD@data <- merge_rastered_spdf(ODD, rastered_buildings_spdf, 'nBuildings')
  
  #sedacs2020 <- GetPopulationBbox(dir, ODD@bbox, yr=2020)
  #population2020 <- merge(ODD@coords, cbind(sedacs2020@coords, population2020=sedacs2020$Population), by=c('Longitude', 'Latitude'), all.x=T, sort=F)
  #nonzero_pop <- which(population2020$population2020 > 0)
  #ODD$nBuildings[nonzero_pop] <- round(ODD$nBuildings[nonzero_pop] * (ODD$Population[nonzero_pop] / population2020$population2020[nonzero_pop]))
  
  return(ODD)
}


#------------------------------------------------------------------------------------------
#----------------------COMPARE BUILDING COUNTS FOR PHL EVENT-------------------------------
#------------------------------------------------------------------------------------------

compare_building_counts <- function(){
  # ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Nov24Agg/ODDobjects/Train/EQ20191029PHL_125')
  # ODDy %<>% AddBuildingCounts(125)
  # ODDy %<>% GetOSMbuildingsODD()
  # saveODD(ODDy, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/CaseStudies/EQ20191029PHL_125_3BuildingCounts')
  
  ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/CaseStudies/EQ20191029PHL_125_3BuildingCounts')
  
  plot_df <- as.data.frame(ODDy[[c('nBuildings', 'OpenBuildings', 'OSM')]])
  plot_df$OpenBuildings[which(is.na(plot_df$OpenBuildings))] = 0
  plot_df$OSM[which(is.na(plot_df$OSM))] = 0
  ggplot(as.data.frame(ODDy[[c('nBuildings', 'OpenBuildings', 'OSM')]])) + 
    geom_point(aes(x = nBuildings, y = OSM, color = "OSM"), cex = 1) + 
    geom_point(aes(x = nBuildings, y = OpenBuildings, color = "Microsoft Open Buildings"), cex = 1) + 
    geom_abline(slope = 1, intercept = 0, color="#D62728") + 
    scale_color_manual(values = c("OSM" = "#440154", "Microsoft Open Buildings" = "#1FA187")) + 
    xlab('Building Count: Bing') + 
    ylab('Building Count: OSM / Microsoft') + 
    theme_minimal() + 
      labs(color = "") +
    theme(
      axis.title.y = element_text(family = "Times New Roman", size = 12),
      axis.text.x = element_text(family = "Times New Roman", size = 12),
      axis.text.y = element_text(family = "Times New Roman", size = 12),
      axis.title.x = element_text(family = "Times New Roman", size = 12),
      plot.title = element_text(family = "Times New Roman", size = 14),
      #panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = unit(c(0, 20, 0, 15), "pt"),
      legend.text = element_text(family = "Times New Roman", size = 12),
      legend.title = element_blank()  # Remove legend title
    ) 
  
  #4 x 6, BuildingCountComparison
  
  plot(values(ODDy$nBuildings), values(ODDy$OpenBuildings), pch=19, cex=0.3, xlab='Building Count: Bing', ylab='Building Count: OSM/Microsoft')
  points(values(ODDy$nBuildings), values(ODDy$OSM), pch=19, cex=0.3, col='red')
  abline(0,1)
}



#-----------------------

# folderin<- paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_June24_AggFactor5/ODDobjects/') #"/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/"
# ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
# for (ODD_file in ufiles[5:length(ufiles)]){
#   print(ODD_file)
#   ODDy <- readRDS(paste0(folderin, ODD_file))
#   if (is.null(ODDy$nBuildings)){
#     ODDy <- getBingBuildingsGlobal(ODDy)
#     saveRDS(ODDy, paste0(folderin, ODD_file))
#   }
#   #saveRDS(ODDy, paste0(folderin, ufiles[i])) 
# }










# plot(train_dat$nBuiltup, train_dat$nBuildings_OSM, col=train_dat$bbox_i)
# plot(log(train_dat$nBuiltup), log(train_dat$nBuildings_OSM), col=train_dat$bbox_i)
# 
# #INDIA: 
# bbox_list <- list(rbind(c(88.69,88.86), c(22.14,22.29)), 
#                   rbind(c(76.99, 77.098), c(28.412,28.511)),
#                   rbind(c(75.34, 75.462), c(11.947,12.054)),
#                   rbind(c(72.997, 73.0335), c(20.232,20.307)))
#

# INDIA (by road completeness):
# bbox_list <- list(rbind(c(73.07,73.20), c(19.16,19.275)), 
#                   rbind(c(73.724, 73.865), c(19.95,20.036)),
#                   rbind(c(76.143, 76.30), c(10.492,10.622)))


# #PHILIPPINES:
# bbox_list <- list(rbind(c(120.414, 120.658), c(15.534, 15.659)),
#                   rbind(c(120.877, 120.984), c(14.28, 14.36)))

  #BANGLADESH:

# loon_GADM <- getGADM('Loon', level=2, country='PHL')
# building_locs_spdf <- SpatialPointsDataFrame(building_locs, data.frame(id=1:NROW(building_locs)))
# proj4string(building_locs_spdf) <- proj4string(loon_GADM)
# loon_buildings <- building_locs_spdf[loon_GADM,]
# 
# plot(loon_buildings)
# points(ODDy@coords[which(ODDy@data$Population > 0),], ylim=c(9.7,9.9), xlim=c(123.7,123.9), col='green', pch=16)
# points(ODDy@coords[intersect(which(inPoly(loon_GADM@polygons[[1]], pop = ODDy@coords)$indies),which(ODDy@data$Population > 0)),], ylim=c(9.7,9.9), xlim=c(123.7,123.9), col='red', pch=16)
# points(ODDy@coords[intersect(which(inPoly(loon_GADM@polygons[[1]], pop = ODDy@coords)$indies), which(ODDy@data$ISO3C=='PHL')),], col='blue', pch=16)
# points(ODDy@coords[intersect(which(inPoly(loon_GADM@polygons[[1]], pop = ODDy@coords)$indies), which(is.na(ODDy@data$ISO3C))),], col='yellow')


#create a smaller bbox to increase speed:
# max.lat <- 10.8
# max.lon <- 122.7
# min.lon <- 122.5
# min.lat <- 10.6
# 
# bbox <- rbind(c(min.lon, max.lon), c(min.lat,max.lat))
# 
# SEDAC <- GetPopulationBbox(dir, bbox=bbox, yr='2020') #Get SEDAC Population
# 
# OSM_buildings<-GetOSMbuildingsBbox(bbox, timeout=60) #haven't tested this yet
# rastered_OSM_count <- rasterize(cbind(buildings$Longitude, buildings$Latitude), raster(SEDAC), fun='count') # Get OSM building counts
# data_merged <- merge_rastered_spdf(SEDAC, as(rastered_OSM_count, "SpatialPixelsDataFrame"), 'nBuildings_OSM')
# 
# FBpop<-readFBpop(bbox) # Get Population estimates from Meta Data for Good
# rastered_FBpop <- rasterize(FBpop, raster(SEDAC), 'Population', fun='sum') 
# rastered_FBpop <- dropLayer(rastered_FBpop, 1)
# data_merged <- merge_rastered_spdf(SEDAC, as(rastered_FBpop, "SpatialPixelsDataFrame"), 'FBpop2')
# 
# openBuildingsDat <- read.csv('/home/manderso/Downloads/data')
# rastered_openBuildings_count <- rasterize(cbind(openBuildingsDat$longitude, openBuildingsDat$latitude), raster(SEDAC), fun='count')
# data_merged <- merge_rastered_spdf(SEDAC, as(rastered_openBuildings_count, "SpatialPixelsDataFrame"), 'nBuildings_OB')
# ggplot(data_merged %>% filter(nBuildings_OB < 100), aes(x=as.factor(nBuildings_OB), y=Population)) + geom_violin()
# 
# #on fb pop granularity
# xmn = min(FBpop@coords[,1]); xmx = max(FBpop@coords[,1])
# ymn = min(FBpop@coords[,2]); ymx = max(FBpop@coords[,2])
# FBpop_raster <- rasterize(FBpop, raster(nrows=(ymx-ymn)/(1/3600), ncols=(xmx-xmn)/(1/3600), xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx), fun='sum')
# #FBpop_raster <- dropLayer(FBpop_raster, 1)
# rastered_openBuildings_count <- rasterize(cbind(openBuildingsDat$longitude, openBuildingsDat$latitude), FBpop_raster, fun='count')
# #data_merged <- merge_rastered_spdf(FBpop, as(rastered_openBuildings_count, "SpatialPixelsDataFrame"), 'nBuildings_OB')
# raster_stack <- stack(list(FBpop=FBpop_raster, nBuildings_OB=rastered_openBuildings_count))
# plot_df <- data.frame(
#   Population = raster_stack[[2]]@data@values,
#   NBuildings = raster_stack[[3]]@data@values
# )
# ggplot(plot_df, aes(x=as.factor(NBuildings), y=Population)) + geom_violin()
# plot(raster_stack[[3]]@data@values, raster_stack[[2]]@data@values, xlab='# of Buildings', ylab='Population')
# 
# 
# #FBpop_raster2 <- rasterize(as(FBpop_raster, "SpatialPixelsDataFrame"), raster(SEDAC), fun='sum')
# #FBpop_raster2 <- dropLayer(FBpop_raster2, 1)
# #data_merged <- merge_rastered_spdf(SEDAC, as(FBpop_raster2, "SpatialPixelsDataFrame"), 'FBpop')
# rastered_openBuildings2 <- rasterize(cbind(openBuildingsDat$longitude, openBuildingsDat$latitude), raster(SEDAC), fun='count')
# data_merged <- merge_rastered_spdf(SEDAC, as(rastered_openBuildings2, "SpatialPixelsDataFrame"), 'nBuildings_OB')
# 
# #Philippines Poisson Comparison
# xgr <- seq(0,10,1)
# y_pois <- dpois(xgr, 4.23)
# plot(xgr, y_pois, xlab='Household Size', ylab='Probability')
# points(c(1,2.5,4.5), c(.0918, .15025, .1841), col='red')
# 
# pop_vs_building_count_prob <- matrix(0, 100, 30)
# intervals_90 <- matrix(0,100, 2)
# for (i in 1:NROW(pop_vs_building_count_prob)){
#   prob_sum <- 0
#   for (j in 1:min(i, 30)){
#     pop_vs_building_count_prob[i,j] <- dpois(i, 4.23 * j)
#   }
#   pop_vs_building_count_prob[i,] <- pop_vs_building_count_prob[i,] / sum(pop_vs_building_count_prob[i,] )
#   intervals_90[i, 1] <- min(which(cumsum(pop_vs_building_count_prob[i,])>0.05))
#   intervals_90[i, 2] <- min(which(cumsum(pop_vs_building_count_prob[i,])>0.95))-1
# }
# intervals_90[1,] <- c(1,1) #adjust first row
# 
# 
# hist(rastered_openBuildings_count@data@values, breaks=100)
# 
# xgr <- seq(1,25,1)
# y_pois <- dgeom(xgr-1, 1/3)
# foo = hist(rastered_openBuildings_count@data@values, freq=F, xlab='Number of Houses per 30 x 30 arcseconds (Philippines)', ylab='Probability', breaks=seq(0.5,30.5,1))
# axis(side=1, at=foo$mids, labels=seq(1,30,1))
# points(xgr, y_pois, col='red')
# 
# xgr <- seq(1,25,1)
# y_pois <- dnbinom(xgr-1, 2, 1/3)
# plot(xgr, y_pois)
# 
# 
# rastered_FBpop_builtup <- rasterize(FBpop, raster(SEDAC), 'Population', fun='count') 
# data_merged <- merge_rastered_spdf(SEDAC, as(rastered_FBpop_builtup, "SpatialPixelsDataFrame"), 'nBuiltup')
# 
# probs_nbinom <- matrix(0, nrow=810, ncol=4000)
# for (s in 1:NROW(probs_nbinom)){
#   print(s)
#   for (j in s:NCOL(probs_nbinom)){
#     probs_nbinom[s,j] <- pnbinom(j-s, s, 1/3)
#   }
# }
# 
# notnans <- which(!is.na(data_merged$nBuiltup))
# data_merged$building_lb <- 0
# data_merged$building_lb[notnans] <- apply(probs_nbinom[data_merged$nBuiltup[notnans],], 1, function(x){min(which(x>0.05))} )
# data_merged$building_ub <- 0
# data_merged$building_ub[notnans] <- apply(probs_nbinom[data_merged$nBuiltup[notnans],], 1, function(x){max(which(x<0.95))} )
# data_merged$building_med <- 0
# data_merged$building_med[notnans] <- apply(probs_nbinom[data_merged$nBuiltup[notnans],], 1, function(x){max(which(x<0.5))} )
# 
# rastered_openBuildings <- rasterize(cbind(openBuildingsDat$longitude, openBuildingsDat$latitude), raster(SEDAC), fun='count')
# data_merged2 <- merge_rastered_spdf(SEDAC, as(rastered_openBuildings, "SpatialPixelsDataFrame"), 'nBuildings_OB')
# 
# data_merged$nBuildings_OB <- data_merged2$nBuildings_OB
# 
# data_merged_sort <- data_merged[order(data_merged[,'nBuildings_OB']),]
# ggplot(data_merged_sort[notnans,], aes(x=1:length(notnans))) +
#   geom_line(aes(y=nBuildings_OB)) + 
#   geom_line(aes(y=building_med), col='red') + 
#   geom_ribbon(aes(ymin = building_lb, ymax = building_ub), alpha = 0.1) + 
#   xlab('Index') + ylab('Number of Buildings')
# 
# npoints_plot <- 1:50
# plot(npoints_plot, data_merged2$nBuildings_OB[notnans[npoints_plot]])
# points(npoints_plot, data_merged$building_med[notnans[npoints_plot]], col='red')
# 
# fit1 <- lm(nBuildings_OB ~ nBuiltup, data=data_merged)
# hist(fit1$residuals, breaks=50)
# data_merged$building_med[notnans] <- fit1$coefficients[1] + fit1$coefficients[2] * data_merged$nBuiltup[notnans]
# data_merged$building_lb[notnans] <- qnorm(0.05, data_merged$building_med[notnans], sd= sigma(fit1))
# data_merged$building_ub[notnans] <- qnorm(0.95, data_merged$building_med[notnans], sd= sigma(fit1))
# 
# fit2 <- lm(log(nBuildings_OB) ~ log(nBuiltup), data=data_merged)
# data_merged$building_med[notnans] <- qlnorm(0.5, fit2$coefficients[1] + fit2$coefficients[2] * log(data_merged$nBuiltup[notnans]), sd= sigma(fit2))
# data_merged$building_lb[notnans] <- qlnorm(0.05, fit2$coefficients[1] + fit2$coefficients[2] * log(data_merged$nBuiltup[notnans]), sd= sigma(fit2))
# data_merged$building_ub[notnans] <- qlnorm(0.95, fit2$coefficients[1] + fit2$coefficients[2] * log(data_merged$nBuiltup[notnans]), sd= sigma(fit2))
# 
# plot(log(data_merged$nBuiltup), log(data_merged$nBuildings_OB), xlab='# Pixels Builtup', ylab='# Buildings')
# abline(fit2$coefficients[1], fit2$coefficients[2])
# 
# # TEST ON A DIFFERENT AREA: 
# max.lat.test <- 11.3
# max.lon.test <- 123.2
# min.lon.test <- 123.0
# min.lat.test <- 11.1
# bbox_test <- rbind(c(min.lon.test, max.lon.test), c(min.lat.test,max.lat.test))
# SEDAC_test <- GetPopulationBbox(dir, bbox=bbox_test, yr='2020') #Get SEDAC Population
# FBpop_test<-readFBpop(bbox_test)
# 
# rastered_FBpop_builtup_test <- rasterize(FBpop_test, raster(SEDAC_test), 'Population', fun='count') 
# data_merged_test <- merge_rastered_spdf(SEDAC_test, as(rastered_FBpop_builtup_test, "SpatialPixelsDataFrame"), 'nBuiltup')
# 
# rastered_openBuildings_test <- rasterize(cbind(openBuildingsDat$longitude, openBuildingsDat$latitude), raster(SEDAC_test), fun='count')
# data_merged2_test <- merge_rastered_spdf(SEDAC_test, as(rastered_openBuildings_test, "SpatialPixelsDataFrame"), 'nBuildings_OB')
# 
# data_merged_test$nBuildings_OB <- data_merged2_test$nBuildings_OB
# plot(log(data_merged_test$nBuiltup), log(data_merged_test$nBuildings_OB), xlab='# Pixels Builtup', ylab='# Buildings')
# abline(fit2$coefficients[1], fit2$coefficients[2])
# 
# notnans_test <- which(!is.na(data_merged_test$nBuildings_OB))
# data_merged_test$building_med[notnans_test] <- qlnorm(0.5, fit2$coefficients[1] + fit2$coefficients[2] * log(data_merged_test$nBuiltup[notnans_test]), sd= sigma(fit2))
# data_merged_test$building_lb[notnans_test] <- qlnorm(0.05, fit2$coefficients[1] + fit2$coefficients[2] * log(data_merged_test$nBuiltup[notnans_test]), sd= sigma(fit2))
# data_merged_test$building_ub[notnans_test] <- qlnorm(0.95, fit2$coefficients[1] + fit2$coefficients[2] * log(data_merged_test$nBuiltup[notnans_test]), sd= sigma(fit2))
# 
# data_merged_sort_test <- data_merged_test[order(data_merged_test[,'nBuildings_OB']),]
# data_merged_sort_test[-notnans_test, c('building_med', 'building_lb', 'building_ub')] <- 0
# ggplot(data_merged_sort_test[notnans_test,], aes(x=1:length(notnans_test))) +
#   geom_line(aes(y=nBuildings_OB)) + 
#   geom_line(aes(y=building_med), col='red') + 
#   geom_ribbon(aes(ymin = building_lb, ymax = building_ub), alpha = 0.1) + 
#   xlab('Index') + ylab('Number of Buildings')
# 
# folderin <- '/home/manderso/Documents/GitHub/IIDIPUS_InputRealUnedited/BDobjects/'
# ufiles<-na.omit(list.files(path=folderin,pattern='EQ',recursive = T,ignore.case = T)) #looseend
# classifications <- c('destroyed', 'notaffected', 'possible')
# sum = 0
# for (filer in ufiles){
#   BDy<-readRDS(paste0(folderin,filer))
#   sum = sum + NROW(BDy@data)
# }
# sum


