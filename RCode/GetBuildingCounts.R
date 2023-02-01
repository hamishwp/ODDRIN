library(sqldf)


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
  data <- data[-which(colnames(data)=='id')]
  return(data)
}

AddOpenBuildingCounts <- function(ODD){
  lon_min <- ODD@bbox[1,1]; lon_max <- ODD@bbox[1,2]
  lat_min <- ODD@bbox[2,1]; lat_max <- ODD@bbox[2,2]
  iso3_unique <- unique(ODD$ISO3C)[!is.na(unique(ODD$ISO3C))]
  open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/open_buildings_v2_points_your_own_wkt_polygon_',iso3_unique[1],'.csv')
  if(!file.exists(open_buildings_file)){
    stop('Download Open Buildings data from: https://colab.research.google.com/github/google-research/google-research/blob/master/building_detection/open_buildings_download_region_polygons.ipynb')
          #When downloading, insert the polygon bounding the country of interest into your_own_wkt_polygon field
          # alternatively, can select the region_border_source and country, although this doesn't seem to work well when working with small islands in Philippines
          # select 'points' in the data_type field
  }
  building_locs <- read.csv.sql(open_buildings_file,
                                paste0("select longitude, latitude from file where latitude > ", lat_min, ' AND longitude > ', lon_min, 
                                       ' AND latitude < ',lat_max, ' AND longitude < ', lon_max))
  
  if (length(iso3_unique) > 1){
    for (i in 2:length(iso3_unique)){
      open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/open_buildings_v2_points_your_own_wkt_polygon_',iso3_unique[i],'.csv')
      if(!file.exists(open_buildings_file)){
        stop('Download Open Buildings data from: https://colab.research.google.com/github/google-research/google-research/blob/master/building_detection/open_buildings_download_region_polygons.ipynb')
        #follow download instructions above
      }
      building_locs %<>% rbind(read.csv.sql(open_buildings_file,
                                    paste0("select longitude, latitude from file where latitude > ", lat_min, ' AND longitude > ', lon_min, 
                                           ' AND latitude < ',lat_max, ' AND longitude < ', lon_max)))
    }
  }
  
  rastered_buildings <- rasterize(building_locs, raster(ODD), 1, fun='count')
  rastered_buildings_spdf <- as(rastered_buildings, "SpatialPixelsDataFrame")
  
  ODD@data <- merge_rastered_spdf(ODD, rastered_buildings_spdf, 'nBuildings')
  return(ODD)
}

AddBuildingCounts <- function(ODD){
  iso3_unique <- unique(ODD$ISO3C)[!is.na(unique(ODD$ISO3C))]
  if (all(iso3_unique %in% c('BGD', 'COD', 'IDN', 'NPL', 'PHL'))){
    ODD %<>% AddOpenBuildingCounts()
  } else if (all(iso3_unique %in% c('CHN', 'IND'))){
    stop('Functionality not yet added for building counts in countries not covered by Open Buildings')
    ODD %<>% AddFBBuildingEstimates()
  } else {
    return(ODD)
  }
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

AddFBBuildingEstimates <- function(ODD){
  bbox <- ODD@bbox
  FBpop<-readFBpop(bbox) # Get Population estimates from Meta Data for Good
  rastered_FBpop_builtup <- rasterize(FBpop, raster(ODD), 'Population', fun='count')
  ODD@data <- merge_rastered_spdf(SEDAC, as(rastered_FBpop_builtup, "SpatialPixelsDataFrame"), 'nBuiltup')
  return(ODD)
}

sampleBuildingsFromBuiltup <- function(nBuiltup, reg_int, reg_coef, reg_sd, Np, notnans){
  nBuildings_sample <- array(0, dim=c(length(nBuiltup),Np))
  notnans <- notnans[!(notnans %in% which(nBuiltup==0))] #could there also be any NA values? 
  for (ij in notnans){
    nBuildings_sample[ij,] <- round(rlnorm(Np, reg_int+reg_coef*log(nBuiltup[ij]), reg_sd))
  }
  return(nBuildings_sample)
}

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


