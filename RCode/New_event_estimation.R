
source('RCode/GetSubNationalData.R')
source('RCode/AggregateODD.R')

#Test:
#USGSid = 'us6000m0xl'
iso3= 'MAR' #'JPN'
date = '2023-09-08' #'2024-01-01'
#lhazSDF_JPN20240101 <- extractHAZARD(USGSid)
lhazSDF_JPN20240101 <- extractHAZARD(iso3=iso3, sdate=as.Date(date))
plotHAZARD(lhazSDF_JPN20240101)

ODDy_JPN20240101 <- createODD(lhazSDF_JPN20240101)
gadm_regions <- find_subnat_regions(ODDy_JPN20240101, -1)
polygons_list <- create_polygons_list(gadm_regions)

ODDy_JPN20240101 <- addODDPolygons(ODDy_JPN20240101, polygons_list)
ODD_JPN20240101 <- addPolygonsToImpact(ODDy_JPN20240101)
plotODDy(ODDy_JPN20240101)
plotODDy_GADM(ODDy_JPN20240101)

AggLevel <- 3
ODDyAgg <- aggregateODDbyX(ODDy_JPN20240101, AggLevel)
namer<-paste0(ODDy@hazard, str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),iso3, "_",i, '_AggLevel', AggLevel)
folder_write <- 'IIDIPUS_Input_NewEvents/'
saveRDS(ODDyAgg, paste0(dir, folder_write, "ODDAggobjects/",namer))

AggLevel <- 5
ODDyAgg <- aggregateODDbyX(ODDy_JPN20240101, AggLevel)
namer<-paste0(ODDy@hazard, str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),iso3, "_",i, '_AggLevel', AggLevel)
folder_write <- 'IIDIPUS_Input_NewEvents/'
saveRDS(ODDyAgg, paste0(dir, folder_write, "ODDAggobjects/",namer))

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-10_192933_alpha0.9_M60_Npart1000RealAgg5')
ODDyWithImpact <- sampleImpactODD(ODDyAgg, AlgoResults)

plotODDy(ODDyWithImpact, var='mortality.mean') # doesn't show up on map very well. Try to improve plot
plotODDy(ODDyWithImpact, var='displacement.mean')


ODDyWithAggImpact <- aggregateSampledImpact(ODDyWithImpact)
saveRDS(ODDyWithAggImpact, paste0(dir, folder_write, "ODDAggobjects/",namer, '_WithImpact100'))
#ODDyWithAggImpact <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20240101JPN_WithImpact_1000')


extractHAZARD <- function(USGSid = NULL, iso3 = NULL, bbox = NULL, sdate = NULL, fdate=NULL){
  #Need to specify either:
  #      - GDACS_id
  #      - bbox, sdate, fdate
  if (is.null(USGSid) & any(is.null(iso3), (is.null(sdate)))) stop('Specify either USGSid or iso3 and sdate')
  
  haz="EQ"; EQparams=list(I0=4.3, minmag=5)
  if (is.null(USGSid)){
    if(is.null(fdate)) fdate<-sdate+3
    if (is.null(bbox)) bbox<-countriesbbox(unique(iso3))
    miniDamSimplified <- data.frame(iso3=iso3, sdate=sdate, fdate=fdate, eventid=NA, hazard=haz)
    bbox<-countriesbbox(unique(miniDamSimplified$iso3))
    lhazSDF<-tryCatch(GetDisaster(miniDamSimplified,bbox=bbox, EQparams = EQparams),error=function(e) NULL)
    if(is.null(lhazSDF)){
      stop('Failed to download hazard using provided bounding box and dates')
    }
  } else {
    lhazSDF <- tryCatch(GetUSGS(USGSid,bbox=NULL,sdate=NULL,fdate=NULL,titlz="tmp",I0=EQparams$I0,minmag=EQparams$I0),error=function(e) NULL)
    if(is.null(lhazSDF)){
      stop('Failed to download hazard using provided USGSid')
    }
  }
  return(lhazSDF)
}

plotHAZARD <- function(lhazSDF){
  #This is just to check that the hazards are in relatively the same location 
  # and haven't picked up hazards from other parts of the country
  bbox <- lhazSDF[[1]]$bbox
  hazSDF_raster <- lapply(lhazSDF_JPN20240101[2:length(lhazSDF_JPN20240101)], raster)
  
  plot(NA, NA, xlim=bbox[c(1,3)], ylim=bbox[c(2,4)])
  for (i in 1:length(hazSDF_raster)){
    contour(hazSDF_raster[[i]], add=TRUE, levels=seq(4,9,0.5))
  }
  # plot(raster_layer, main="Contour Plot")
  # contour(raster_layer, add=TRUE)
  # contour(raster(lhaz_list[[8]]), add=TRUE)
  # lhaz_list <- lapply(lhazSDF_JPN20240101[2:length(lhazSDF_JPN20240101)], raster)
  # 
  # lhaz_list <- lhazSDF_JPN20240101[2:length(lhazSDF_JPN20240101)]
  # common_grid <- raster(lhaz_list[[1]])
  # raster_list <- lapply(lhaz_list, function(spdf) {
  #   r <- raster(spdf)
  #   resample(r, common_grid, method = "bilinear")
  # })
  # raster_stack <- stack(raster_list)
  # max_raster <- calc(raster_stack, fun = max, na.rm = TRUE)
  # plot(max_raster)

}

createODD <- function(lhazSDF, eventid = -1){
  #Create the ODD object:
  #Set event_d = -1 for new events to be predicted
  
  miniDamSimplified = data.frame(iso3=NA, sdate=min(lhazSDF[[1]]$eventdates), fdate=max(lhazSDF[[1]]$eventdates), 
             eventid=eventid, hazard=lhazSDF[[1]]$hazard)
  
  ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified, agg_level=1),error=function(e) NULL)
  if(is.null(ODDy)) {
    print('Failed to create ODD object')
  }
  
  #Fetch building count data:
  #ODDy_build <- tryCatch(AddBuildingCounts(ODDy, i, paste0(folder_write, 'Building_count_notes')), error=function(e) NULL)
  ODDy_build <- tryCatch(getBingBuildingsGlobal(ODDy, eventid, file_write=NULL), error=function(e) NULL)
  if(is.null(ODDy_build)) {
    print('No/insufficient building data available for the event from Bing Building Footprints')
  } else {
    ODDy <- ODDy_build
    rm(ODDy_build)
  }
  
  # name event using country code exposed to the highest hazard intensity
  iso3_main <- ODDy$ISO3C[which(ODDy@data[,grep('hazMean', names(ODDy@data))]==max(ODDy@data[,grep('hazMean', names(ODDy@data))], na.rm=T), arr.ind=T)[1]]
  
  namer<-paste0(ODDy@hazard,
                str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                iso3_main,
                "_",i)
  
  folder_write <- 'IIDIPUS_Input_NewEvents/'
  
  HAZARDpath<-paste0(dir,folder_write, "HAZARDobjects/",namer)
  saveRDS(lhazSDF,HAZARDpath)
  
  ODDpath<-paste0(dir, folder_write, "ODDobjects/",namer)
  saveRDS(ODDy,ODDpath)
  
  return(ODDy)
}

find_subnat_regions <- function(ODDy, eventid, gadm_levels=c(1,2), print_to_xl=F){
  
  iso3_unique <- unique(ODDy$ISO3C)[which(!is.na(unique(ODDy$ISO3C)))]
  gadm_regions_1 <- data.frame('Country'=character(), 'Region'=character(), 'Subregion'=character())
  gadm_regions_2 <- data.frame('Country'=character(), 'Region'=character(), 'Subregion'=character())
  for (gadm_level in gadm_levels){
    pixels_unallocated <- which(!is.na(ODDy$ISO3C))
    for (iso3 in iso3_unique){
      regions_gadm_level <- tryCatch(getGADM(level=gadm_level, country=iso3),error=function(e) NULL)
      if(is_null(regions_gadm_level)){
        print(paste('No gadm data found for', iso3))
        next
      }
      if (gadm_level == 2){
        regions_gadm_level <- sapply(seq_along(regions_gadm_level), function(i) paste0(regions_gadm_level[i], ', ', names(regions_gadm_level)[i]))
      }
      for (r in regions_gadm_level){
        if (gadm_level==2){
          r <- strsplit(r,',')[[1]]
          r[2] <- substring(r[2], 2)
        }
        r_poly <- tryCatch(getGADM(r, level = gadm_level, country=iso3),error=function(e) NULL)
        if(is_null(r_poly)){
          print(paste('No polygon data found for', r))
          next
        }
        overlap <- (max(r_poly@bbox[1,]) >= ODDy@bbox[1] &
                      min(r_poly@bbox[1,]) <= ODDy@bbox[3] &
                      max(r_poly@bbox[2,]) >= ODDy@bbox[2] &
                      min(r_poly@bbox[2,]) <= ODDy@bbox[4])
        if (!overlap) next
        
        if (length(pixels_unallocated)==0) next
        
        spdf_pixels_unallocated <- SpatialPointsDataFrame(coords=ODDy@coords[pixels_unallocated,, drop=F],
                                                          data=ODDy@data[pixels_unallocated,1:2, drop=F], #this is arbitrary, function just seems to require data
                                                          proj4string=r_poly@proj4string)
        spdf_pixels_unallocated$id = 1:length(pixels_unallocated)
        contained <- tryCatch(terra::intersect(r_poly, spdf_pixels_unallocated)$id, error=function(e) c())
        if (length(contained)>0){
          #print(paste('Polygon', r, 'contains an unallocated pixel'))
          pixels_unallocated <- pixels_unallocated[-contained]
          r <- rev(r)
          if (gadm_level==1){
            gadm_regions_1 <- add_row(gadm_regions_1, Country=countrycode(iso3, origin='iso3c', destination='country.name'),
                                      Region=r,
                                      Subregion=NA)
          } 
          else {gadm_regions_2 <- add_row(gadm_regions_2, Country=countrycode(iso3, origin='iso3c', destination='country.name'),
                                          Region=r[1],
                                          Subregion=r[2])}
        }
      }
    }
  }
  return(rbind(gadm_regions_1, gadm_regions_2))
}

create_polygons_list <- function(gadm_regions){
  polygons_list <- list()
  nans_polygon_id <- c()
  
  gadm_regions$polygon_name <- paste0(ifelse(is.na(gadm_regions$Subregion),'',paste0(gadm_regions$Subregion,', ')), 
                                      ifelse(is.na(gadm_regions$Region),'',paste0(gadm_regions$Region,', ')), 
                                      ifelse(gadm_regions$Country=='TOTAL','TOTAL',gadm_regions$Country))
  
  # group by polygon name and assign each an id
  gadm_regions$polygon_id <- 1:NROW(gadm_regions)
  
  for (i in gadm_regions$polygon_id){
    #collect polygon information for each polygon
    polygon_name <- gadm_regions$polygon_name[i]
    polygons_list[[i]] <- getPolyData(polygon_name, gadm_regions$Subregion[i], gadm_regions$Region[i], gadm_regions$Country[i], countrycode(gadm_regions$Country[i], origin='country.name', destination='iso3c'))
    
    if (is.null(polygons_list[[i]]$sf_polygon) & (!is.na(gadm_regions$Subregion[i]) | !is.na(gadm_regions$Region[i]))) nans_polygon_id %<>% append(i) #remove empty non-national polygons
    # should probably plot polygons here to make sure they all seem reasonable!
  }
  
  
  return(polygons_list)
  
}

addPolygonsToImpact <- function(ODDy){
  ODDy@impact <- data.frame()
  impact_types <- c('mortality', 'displacement')
  if (!is.null(ODDy$nBuildings)){impacts <- c(impacts, 'buildDam')}
  for (i in 1:length(ODDy@polygons)){
    ODDy@impact %<>% rbind(data.frame(iso3=ifelse(length(grep(',',ODDy@polygons[[i]]$name))>0, 
                                                                     countrycode::countrycode(sourcevar = trimws(substring(polygons_list[[j]]$polygon_name, max(unlist(gregexpr(",", polygons_list[[j]]$polygon_name))) + 1)),
                                                                                              origin = "country.name",
                                                                                              destination = "iso3c"),
                                                                     countrycode::countrycode(sourcevar = ODDy@polygons[[i]]$name,
                                                                                              origin = "country.name",
                                                                                              destination = "iso3c") ),
                                                         sdate=ODDy@hazinfo$sdate,
                                                         impact=impact_types,
                                                         observed=NA,
                                                         qualifier=NA,
                                                         build_type=NA,
                                                         polygon=i,
                                                         inferred=F))
  }
  
  iso3_unique <- unique(ODDy$ISO3C)[which(!is.na(unique(ODDy$ISO3C)))]
  if (length(iso3_unique) == 1){
    tot_poly <- length(ODDy@polygons) + 1
    ODDy@impact %<>% rbind(data.frame(iso3=iso3_unique,
                              sdate=ODDy@hazinfo$sdate,
                              impact=impact_types,
                              observed=NA,
                              qualifier=NA,
                              build_type=NA,
                              polygon=tot_poly,
                              inferred=F))
    ODDy@polygons[[tot_poly]] = list(name=countrycode::countrycode(sourcevar = iso3_unique,origin = "iso3c",destination = "country.name"), 
                                                 indexes=1:NROW(ODDy@data), weights=rep(1, NROW(ODDy@data)))
  } else {
    for (j in 1:length(iso3_unique)){
      added_poly <-  length(ODDy@polygons) + 1
      ODDy@impact %<>% rbind(data.frame(iso3=iso3_unique[j],
                                sdate=ODDy@hazinfo$sdate,
                                impact=impact_types,
                                observed=NA,
                                qualifier=NA,
                                build_type=NA,
                                polygon=added_poly,
                                inferred=F))
      ODDy@polygons[[added_poly]] = list(name=countrycode::countrycode(sourcevar = iso3_unique[j],origin = "iso3c",destination = "country.name"),
                                       indexes=which(ODDy$ISO3C==iso3_unique[j]), weights=rep(1, sum(ODDy$ISO3C==iso3_unique[j],na.rm=T)))
    }
    tot_poly <- length(ODDy@polygons) + 1
    ODDy@impact %<>% rbind(data.frame(iso3='TOT',
                              sdate=ODDy@hazinfo$sdate,
                              impact=impact_types,
                              observed=NA,
                              qualifier=NA,
                              build_type=NA,
                              polygon=tot_poly,
                              inferred=F))
    ODDy@polygons[[tot_poly]] = list(name='TOTAL',indexes=1:NROW(ODDy@data), weights=rep(1, NROW(ODDy@data)))
  }
  return(ODDy)
}

#AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-10_192933_alpha0.9_M60_Npart1000RealAgg5')

sampleImpactODD <- function(ODDy, AlgoResults){
  AlgoResults %<>% addAlgoParams()
  AlgoParamsPostPredict <- list(Np = 1, m_CRPS=1, cores=1)
  N_samples <- 100
  params_sampled <- rep(NA, N_samples)
  sampled <- array(NA, dim=c(NROW(ODDy@data),3, N_samples))
  for (i in 1:N_samples){
    print(paste('Sample', i))
    param_sample <- sample(1:AlgoResults$Npart, 1, replace=T, prob=AlgoResults$W[,AlgoResults$s_finish])
    params_sampled[i] <- param_sample
    proposed = AlgoResults$Omega_sample_phys[param_sample,,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton) %>% addTransfParams()
    #proposed$eps$hazard_cor <- 0.999
    #proposed$eps$local <-
    impact <- DispX(ODDy, proposed, Model$center, Method = AlgoParamsPostPredict, output='SampledFull')
    sampled[,,i] <- impact
  }
  
  for (i in 1:ifelse(is.null(ODDy$nBuildings), 2,3)){
    # i = 1: Displacement; i = 2: Mortality; i = 3: BuildDam. 
    impact_types = c('displacement', 'mortality', 'buildDam')
    impact_sampled <- sampled[,i,]
    colnames(impact_sampled) = paste0(impact_types[i],'.s', 1:N_samples)
    impact_sampled %<>% as.data.frame()
    impact_sampled[,paste0(impact_types[i], '.mean')] = apply(sampled[,i,], 1, mean )
    impact_sampled[,paste0(impact_types[i], '.median')] = apply(sampled[,i,], 1, median )
    impact_sampled[,paste0(impact_types[i], '.q05')] = apply(sampled[,i,], 1, quantile,0.05 )
    impact_sampled[,paste0(impact_types[i], '.q95')] = apply(sampled[,i,], 1, quantile, 0.95 )
    ODDy@data %<>% cbind(impact_sampled)
  } 
  return(ODDy)
}

aggregateSampledImpact <- function(ODDyWithImpact){
  
  n_samples <- length(grep(paste0('mortality', '.s'),colnames(ODDyWithImpact@data)))
  sampled_df <- data.frame(polygon=integer(), impact=character())
  for (i in 1:n_samples){ sampled_df[,paste0('sampled.', i)] = integer()}
  
  # i = 1: Displacement; i = 2: Mortality; i = 3: BuildDam
  impact_types = c('displacement', 'mortality', 'buildDam')
  for (i in 1:ifelse(is.null(ODDyWithImpact$nBuildings), 2,3)){
    impact_sampled <-ODDyWithImpact@data[,grep(paste0(impact_types[i], '.s'),colnames(ODDyWithImpact@data))]
    for (j in 1:length(ODDyWithImpact@polygons)){
        sampled_df[nrow(sampled_df)+1,] = c(j, impact_types[i], colSums(impact_sampled[ODDyWithImpact@polygons[[j]]$indexes,]))
    }
  }
  ODDyWithImpact@impact %<>% merge(sampled_df, by=c('polygon', 'impact'))
  return(ODDyWithImpact)
}

# plotImpactODD <- function(ODDy){
#   plot_table = data.frame(polygon=integer(), impact=character(), sampled=character())
#   # Find the polygon matching each region
#   
#   for (i in 1:ifelse(is.null(ODDy$nBuildings), 2,3)){
#     # i = 1: Displacement; i = 2: Mortality; i = 3: BuildDam
#     impact_types = c('displacement', 'mortality', 'buildDam')
#     impact_sampled <- ODDy@data[,grep(paste0(impact_types[i], '.s'),colnames(ODDy@data))]
#     for (i in 1:length(ODDy@polygons)){
#     }
#     
#   }
#     #region_name <- paste0(gadm_iso$NAME_2[region_i], ', ', gadm_iso$NAME_1[region_i])
#     match <- which(unlist(lapply(ODDyAgg@polygons, function(x) ifelse(length(grep(gadm_names[i], x$name)>0), T, F))))
#     if (length(match)==0) next
#     #impact_region <- colSums(impact_median[ODDyAgg@polygons[[match]]$indexes,1,])
#     plot_table %<>% add_row(id=as.character(i), var=sum(ODDyAgg@data[ODDyAgg@polygons[[match]]$indexes, var]), 
#                             region_name=gadm_names[i])
#   }
# }

plot_correlated_impacts <- function(ODDyWithAggImpact, region1, region2, impact_type='mortality', baseR=F){
  match1 <- which(unlist(lapply(ODDyWithAggImpact@polygons, function(x) ifelse(length(grep(region1, x$name)>0), T, F))))
  match2 <- which(unlist(lapply(ODDyWithAggImpact@polygons, function(x) ifelse(length(grep(region2, x$name)>0), T, F))))
  impact_row1 <- which(ODDyWithAggImpact@impact$polygon == match1 & ODDyWithAggImpact@impact$impact == impact_type)
  impact_row2 <- which(ODDyWithAggImpact@impact$polygon == match2 & ODDyWithAggImpact@impact$impact == impact_type)
  combined_samples <- t(ODDyWithAggImpact@impact[c(impact_row1, impact_row2),grep('sampled.', colnames(ODDyWithAggImpact@impact))])
  colnames(combined_samples) <- c(region1, region2)
  combined_samples %<>% as.data.frame()
  if(baseR){
    return(plot(as.numeric(combined_samples[,1]), as.numeric(combined_samples[,2]), pch=19, xlab=region1, ylab=region2))
  }
  return(ggplot(combined_samples, aes(x=as.numeric(.data[[region1]]), y=as.numeric(.data[[region2]]))) + geom_point()) #+ 
      #scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breakings=c(0,1, 10,100,1000)) + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) + xlab(region1) + ylab(region2))
}

# regions_max_impact <- c('Jōetsu, Niigata', 'Wajima, Ishikawa', 'Suzu, Ishikawa', 'Shika, Ishikawa', 'Noto, Ishikawa', 
#                         'Nanao, Ishikawa', 'Nakanoto, Ishikawa', 'Kanazawa, Ishikawa', 'Hakui, Ishikawa', 'Anamizu, Ishikawa')
# real_observed <- c(12,180, 104, 75, 110, 94, 17, 23, 23, 87)
# par(mfrow=c(6,6))
# for (i in 2:8){
#   for (j in 1:(i-1)){
#     plot_correlated_impacts(ODDyWithAggImpact, regions_max_impact[i], regions_max_impact[j], baseR=T)
#     points(x=real_observed[i], y=real_observed[j], col='red', pch=19)
#   }
# }

# plot_correlated_impacts(ODDyWithAggImpact, 'Jōetsu, Niigata', 'Wajima, Ishikawa') + geom_point(x=0, y=130, col='red')
# plot_correlated_impacts(ODDyWithAggImpact, 'Suzu, Ishikawa', 'Wajima, Ishikawa') + geom_point(x=114, y=130, col='red')
# plot_correlated_impacts(ODDyWithAggImpact, 'Shika, Ishikawa', 'Wajima, Ishikawa') + geom_point(x=2, y=130, col='red')
# plot_correlated_impacts(ODDyWithAggImpact, 'Shika, Ishikawa', 'Suzu, Ishikawa') + geom_point(x=2, y=110, col='red')


# plot_correlated_impacts(ODDyWithAggImpact, 'Morocco', 'Algeria') + geom_point(x=1684 , y=15, col='red')
# 
# 
# 
# plot_correlated_impacts(ODDyWithAggImpact, 'AlHaouz, Marrakech-Tensift-AlHaouz', 'Taroudannt, Souss-Massa-Draâ') + geom_point(x=log(1684+10) , y=log(980+10), col='red')
# ggplot(combined_samples, aes(x=as.numeric(.data[[region1]]), y=as.numeric(.data[[region2]]))) + geom_point() + 
#   xlab(region1) + ylab(region2) + scale_x_log10() + scale_y_log10() + xlab('Al Haouz (Mortality)') + ylab('Taroudant (Mortality)') + geom_point(x=log10(1684), y=log10(980), col='red', pch=4, size=5)

plotODDyAgg <- function(ODDyAgg, ODDy=NULL, zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.7,map="terrain", gadm_level=2){

  bbox <- ODDyAgg@bbox
  iso3_unique <- unique(ODDyAgg$ISO3C)
  iso3_unique <- iso3_unique[!is.na(iso3_unique)]
  gadm_iso <- as(geodata::gadm(country=iso3_unique[1], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial')
  if (length(iso3_unique) > 1){
    for (i in 2:length(iso3_unique)){
      gadm_iso %<>% rbind(as(geodata::gadm(country=iso3_unique[i], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial'))
    }
  }
  #gadm_iso <- gSimplify(gadm_iso, 0.01)
  gadm_iso <- intersect(gadm_iso, bbox)
  gadm_map <- fortify(gadm_iso)
  
  ii <- 1
  if (is.null(gadm_iso$NAME_2)){
    gadm_names <- paste0(gadm_iso$NAME_1)
  } else {
    #gadm_names <- paste0(gadm_iso$NAME_2, gadm_iso$NAME_1)
    gadm_names <- paste0(gsub(" ", "",gadm_iso$NAME_2, fixed=T), ', ', gsub(" ", "",gadm_iso$NAME_1, fixed=T))
  }
  
  plot_table = data.frame(id=character(), var=numeric(), region_name=character())
  # Find the polygon matching each region
  for (i in 1:length(gadm_names)){
    #region_name <- paste0(gadm_iso$NAME_2[region_i], ', ', gadm_iso$NAME_1[region_i])
    match <- which(unlist(lapply(ODDyAgg@polygons, function(x) ifelse(length(grep(gadm_names[i], x$name)>0), T, F))))
    if (length(match)==0) next
    #impact_region <- colSums(impact_median[ODDyAgg@polygons[[match]]$indexes,1,])
    plot_table %<>% add_row(id=as.character(i), var=sum(ODDyAgg@data[ODDyAgg@polygons[[match]]$indexes, var]), 
                            region_name=gadm_names[i])
  }
  
  #add the modeled displacement and mortality
  gadm_map[[var]] <- plyr::join(gadm_map, plot_table, by="id", type='left')$var
  #gadm_map$Mortality <- plyr::join(gadm_map, plot_table, by="id", type='left')$mort
  gadm_map[is.na(gadm_map[,var]), var] = 0
  #gadm_map$Displacement[is.na(gadm_map$Displacement)] = 0
  
  options(scipen=10000)
  gg <- ggplot()
  gg <- gg + geom_map(map=gadm_map, data=gadm_map, aes(map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
  gg <- gg + coord_map() + geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group, fill=!!ensym(var)), color='black', lwd=0.2) + 
    scale_fill_gradient(low = "white",high = "blue",na.value = "transparent", name=var)
  p <- gg + xlab("Longitude") + ylab("Latitude")# + theme(legend.position = "none")  + scale_fill_discrete(name = "New Legend Title")
  
  if (!is.null(ODDy)){
    hazard<-rep(NA_real_,length(ODDy@data$hazMean1))
    for (variable in names(ODDy)[grepl("Mean",names(ODDy))]){
      tmp<-ODDy[variable]
      tmp$hazard<-hazard
      hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
    }
    ODDy@data$hazard<-hazard
    brks<-seq(9,ceiling(2*max(hazard,na.rm = T)),by=1)/2
    
    p<-p+geom_contour(data = as.data.frame(ODDy),
                      mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
                      alpha=1.0,breaks = brks, lwd=0.8) +
      scale_colour_gradient(low = "#FFFF99",high = "red",na.value = "transparent") + 
      labs(colour = "Hazard Intensity")
    
  }
  return(p)
}

#This plot displays much better but code needs a bit of tidying
plotODDy_GADM <- function(ODDy,zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.7,map="terrain"){
  
  gadm_level=2
  bbox <- ODDy@bbox
  iso3_unique <- unique(ODDy$ISO3C)
  iso3_unique <- iso3_unique[!is.na(iso3_unique)]
  gadm_iso <- as(geodata::gadm(country=iso3_unique[1], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial')
  if (length(iso3_unique) > 1){
    for (i in 2:length(iso3_unique)){
      gadm_iso %<>% rbind(as(geodata::gadm(country=iso3_unique[i], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial'))
    }
  }
  #gadm_iso <- gSimplify(gadm_iso, 0.01)
  gadm_iso <- intersect(gadm_iso, bbox)
  gadm_map <- fortify(gadm_iso)
  
  gg <- ggplot()
  gg <- gg + geom_map(map=gadm_map, data=gadm_map, aes(map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
  gg <- gg + coord_map() + geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='white', color='black')
  p <- gg + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
  
  hazard<-rep(NA_real_,length(ODDy@data$hazMean1))
  for (variable in names(ODDy)[grepl("Mean",names(ODDy))]){
    tmp<-ODDy[variable]
    tmp$hazard<-hazard
    hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
  }
  ODDy@data$hazard<-hazard
  brks<-seq(9,ceiling(2*max(hazard,na.rm = T)),by=1)/2
  
  if(var!="hazard")  {
    ODDy@data[is.na(ODDy@data$ISO3C),var]<-NA
    
    # p<-p+geom_contour_filled(data = as.data.frame(ODDy),
    #                          mapping = aes(Longitude,Latitude,z=ODDy@data[[var]], fill=..level..),
    #                          breaks=c(0, 1,5,10,50,100,500,1000, 2000, 5000, 50000), alpha=alpha) + 
    #   #scale_fill_viridis_c() + 
    #   scale_fill_brewer(type = 'div', palette = 'Blues') +
    #   labs(fill = 'Displacement')
    fill_colours <- c('white', 'turquoise', 'deepskyblue', 'purple', 'purple4')
    # fill_colours <- c("white", "#Dcf9bc", "#94D840FF", "#74D055FF",  "#56C667FF", # "#FDE725FF", 
    #                   "#3CBC75FF", "#29AF7FFF", "#20A386FF", "#1F968BFF","#238A8DFF", "#287D8EFF", 
    #                   "#32648EFF", "#39558CFF", "#3F4788FF", "#453781FF", "#482677FF", "#481568FF", 
    #                   "#440154FF")
    # fill_values_quantiles <- seq(0, 1, length.out = length(fill_colours))
    ## use this for a vector of your quantile breaks for the labels (!)
    
    if (is.null(breakings)) {
      breakings_max = max(ODDy@data[[var]], na.rm=T) 
      log_breakings_max = round(log(breakings_max, 10))
      breakings = c(0, 10^seq(1, log_breakings_max, 1))
    }
    options(scipen=10000)
    p <- p+geom_raster(data=as.data.frame(ODDy) %>% filter(!is.na(!!ensym(var))),aes(Longitude,Latitude,fill=!!ensym(var)), #as.numeric(!!ensym(var))),
                       alpha=0.9,interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
      #scale_fill_gradientn( colours = fill_colours, breaks = fill_values_quantiles, labels = round(quants, 3)) + 
      #scale_fill_viridis_c(trans = "pseudo_log", breaks=c(0,10,100,1000, 10000), na.value = "transparent") +
      #scale_fill_gradientn(colors=c('azure','turquoise', 'deepskyblue', 'violet', 'purple', 'purple4'), breaks=breakings, na.value = "transparent", trans = scales::log_trans(base = 10)) +
      scale_fill_gradientn(colors=fill_colours, breaks=breakings, na.value = "transparent", trans = 'pseudo_log') +
      #scale_fill_gradientn(colors=c('azure','paleturquoise', 'deepskyblue', 'midnightblue'), trans = "pseudo_log",breaks=c(0,10,100,1000, 10000), limits=c(0.01,max(ODDy@data[[var]], na.rm=T)),na.value = "transparent") +
      #scale_fill_gradient2(low='white',mid = "blue", high = "midnightblue",trans = "pseudo_log", midpoint=3.9, breaks=c(0,10,100,1000, 10000),na.value = "transparent")+
      labs(fill = var)+xlab("Longitude") + ylab("Latitude") + geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='transparent', color='grey', lwd=0.1)
    
    # geom_contour_filled(data = as.data.frame(ODDy),
    #                     mapping = aes(x=Longitude,y=Latitude,z=ODDy@data[[var]], fill=..level..),
    #                     breaks=c(0, 1,5,10,100, 50000), alpha=alpha) + 
    # #scale_fill_viridis_c() + 
    # scale_fill_brewer( type = 'div', palette = 'Blues') +
    # labs(fill = 'Displacement')
    
    p<-p+geom_contour(data = as.data.frame(ODDy),
                      mapping = aes(Longitude,Latitude,z=hazard,colour=after_stat(level)),
                      alpha=1,breaks = brks, lwd=0.7) +
      scale_colour_gradient(low = "#FFFF99",high = "red",na.value = "transparent") + 
      labs(colour = "Hazard Intensity")
    
    # p+geom_contour_filled(data = as.data.frame(ODDy),
    #                       mapping = aes(Longitude,Latitude,z=1-ODDy@data$tmp),
    #                       fill="green",alpha=alpha)+ 
    #   labs(fill = "Hazard>5")
    
    return(p)
  }
  
  ODDy@data$hazard[ODDy@data$hazard==0]<-NA
  p <- p+geom_contour_filled(data = as.data.frame(ODDy),
                             mapping = aes(Longitude,Latitude,z=hazard),
                             alpha=1,breaks = brks, lwd=0.5) +
    scale_colour_gradient(low = "#FFFF99",high = "red",na.value = "transparent") + 
    labs(fill = "Hazard Intensity")
  
  # p<-p+geom_contour(data = as.data.frame(ODDy),
  #                   mapping = aes(Longitude,Latitude,z=hazard),
  #                   alpha=0.8,breaks = c(6.0),colour="red")
  
}





  
  