# Example:
# ODDyAgg <- prepareODD(iso3='MAR', sdate='2023-09-08')
# PosteriorImpactPred(ODDyAgg)

ODDyAgg <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908MAR_40_AggLevel5_WithImpact100')

# ODDyAgg <- prepareODD(iso3='AFG', sdate='2023-10-07')
# PosteriorImpactPred(ODDyAgg)
# names(ODDyWithAggImpact@data)[which(names(ODDyWithAggImpact@data)=='mortality.mean')] = 'Mean Mortality'
# plotODDy_GADM(ODDyWithAggImpact, var='Mean Mortality')
# #MoroccoMort.pdf, 8 x 12


prepareODD <- function(USGSid = NULL, iso3 = NULL, bbox = NULL, sdate = NULL, fdate=NULL, folder_write='IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/'){
  
  if (!is.null(sdate)){ sdate %<>% as.Date()}
  if (!is.null(fdate)){ fdate %<>% as.Date()}
  
  #create hazard object:
  lhazSDF <- extractHAZARD(USGSid, iso3, bbox, sdate, fdate)
  # check that hazard seems ok:
  plotHAZARD(lhazSDF) 
  
  if (is.null(iso3)){
    iso3 = na.omit(names(sort(table(coords2country(lhazSDF[[2]]@coords)), decreasing=T)))[1]
  }
  namer<-paste0(lhazSDF$hazard_info$hazard, str_remove_all(as.character.Date(min(lhazSDF$hazard_info$eventdates)),"-"),iso3, "_-1")
  saveRDS(lhazSDF, paste0(dir, folder_write, 'HAZARDobjects/', namer))
  
  #create ODD object:
  ODDy <- createODD(lhazSDF)
  saveRDS(ODDy, paste0(dir, folder_write, 'ODDobjects/', namer))
  
  #add polygons for each administrative region
  gadm_regions <- find_subnat_regions(ODDy, -1)
  polygons_list <- create_polygons_list(gadm_regions)
  saveRDS(list(gadm_regions=gadm_regions, polygons_list=polygons_list), paste0(dir, folder_write, 'tmp/', namer, 'polygons')) # save intermediary objects just in case something crashes
  ODDy <- addODDPolygons(ODDy, polygons_list)
  saveRDS(ODDy, paste0(dir, folder_write, 'ODDobjects/', namer)) # again, save frequently just in case something crashes
  ODDy <- addPolygonsToImpact(ODDy)
  saveRDS(ODDy, paste0(dir, folder_write, 'ODDobjects/', namer))
  
  #plot to make sure all is ok:
  plotODDy(ODDy)
  plotODDy_GADM(ODDy)
  
  #aggregate ODD object:
  AggLevel <- 5
  ODDyAgg <- aggregateODDbyX(ODDy, AggLevel)
  
  namer<-paste0(namer, '_AggLevel', AggLevel)
  saveRDS(ODDyAgg, paste0(dir, folder_write, "ODDAggobjects/",namer))
}

PosteriorImpactPred <- function(ODDy, AlgoResultsFilename='HPC/abcsmc_2024-07-10_192933_alpha0.9_M60_Npart1000RealAgg5'){
  #aggregation of ODD object should match that of the ODD objects used to fit AlgoResults
  
  AlgoResults <- paste0(dir, 'IIDIPUS_Results/', AlgoResultsFilename)
  ODDyWithImpact <- sampleImpactODD(ODDy, AlgoResults)
  
  plotODDy(ODDyWithImpact, var='mortality.mean') # doesn't show up on map very well. Try to improve plot
  plotODDy(ODDyWithImpact, var='displacement.mean')
  
  ODDyWithAggImpact <- aggregateSampledImpact(ODDyWithImpact)
  saveRDS(ODDyWithAggImpact, paste0(dir, folder_write, "ODDAggobjects/",namer, '_WithImpact100'))
}

extractHAZARD <- function(USGSid = NULL, iso3 = NULL, bbox = NULL, sdate = NULL, fdate=NULL){
  #Need to specify either:
  #      - GDACS_id
  #      - iso3, sdate, fdate (can further specify bbox to constrain)
  if (is.null(USGSid) & any(is.null(iso3), (is.null(sdate)))) stop('Specify either USGSid or iso3 and sdate')
  
  haz="EQ"; EQparams=list(I0=4.3, minmag=5)
  if (is.null(USGSid)){
    if(is.null(fdate)) fdate<-sdate+3
    if(is.null(bbox)) bbox<-countriesbbox(unique(iso3))
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
  hazSDF_raster <- lapply(lhazSDF[2:length(lhazSDF)], raster)
  
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
  
  return(ODDy)
}

find_subnat_regions <- function(ODDy, eventid, gadm_levels=c(1,2), print_to_xl=F){
  # Find the subnational regions at the provided administrative levels (gadm level 1 or 2)
  # that are exposed to the hazard
  # Find the pixels in the ODD object that are in each region. 
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
  # Create the list for ODD@polygons that, for each list element, contains the polygon name
  # and pixels allocated to the polygon
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
  # Create the dataframe ODD@impact that contains a row for each impact type + polygon combination
  # with the observed impact left blank
  ODDy@impact <- data.frame()
  impact_types <- c('mortality', 'displacement')
  if (!is.null(ODDy$nBuildings)){impact_types <- c(impact_types, 'buildDam')}
  for (i in 1:length(ODDy@polygons)){
    ODDy@impact %<>% rbind(data.frame(iso3=ifelse(length(grep(',',ODDy@polygons[[i]]$name))>0, 
                                                                     countrycode::countrycode(sourcevar = trimws(substring(ODDy@polygons[[i]]$name, max(unlist(gregexpr(",", ODDy@polygons[[i]]$name))) + 1)),
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
      #Add rows for national impact:
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
    #Add rows for total impact:
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

sampleImpactODD <- function(ODDy, AlgoResults, N_samples = 100){
  #Sample from the posterior predictive distribution for the event's impact
  AlgoResults %<>% addAlgoParams()
  AlgoParamsPostPredict <- list(Np = 1, m_CRPS=1, cores=1)
  params_sampled <- rep(NA, N_samples)
  sampled <- array(NA, dim=c(NROW(ODDy@data),3, N_samples))
  for (i in 1:N_samples){
    print(paste('Sample', i))
    param_sample <- sample(1:AlgoResults$Npart, 1, replace=T, prob=AlgoResults$W[,AlgoResults$s_finish])
    params_sampled[i] <- param_sample
    proposed = AlgoResults$Omega_sample_phys[param_sample,,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton) %>% addTransfParams()
    impact <- DispX(ODDy, proposed, Model$center, Method = AlgoParamsPostPredict, output='SampledFull')
    sampled[,,i] <- impact
  }
  
  for (i in 1:ifelse(is.null(ODDy$nBuildings), 2,3)){
    # i = 1: Displacement; i = 2: Mortality; i = 3: BuildDam. 
    # Add a column to ODDy@data for each posterior predictive impact sample, the mean, median, and the 5th and 95th percentiles
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
  #Aggregate the sampled impact to add posterior predictive samples to ODDy@impact:
  
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#==============================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#====== Functions for adding observed and comparing to predictions ============
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#==============================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addPolygonNames <- function(ODDy){
  #Adds polygon names to ODDy@impact
  polygon_names <- unlist(lapply(ODDy@polygons, function(x) x$name))
  polygon_names_id <- data.frame(polygon_name=polygon_names, id=1:length(polygon_names))
  ODDy@impact <- merge(ODDy@impact, polygon_names_id, by.x='polygon', by.y='id')
  return(ODDy)
}

addObs <- function(ODDyWithAggImpact, observed_data, impact_type = 'mortality', add_zeros=T){
  #Adds observations to ODDy@impact (with observations provided in the data frame observed_data)
  ODDyWithAggImpact %<>% addPolygonNames()
  impact_with_obs <- ODDyWithAggImpact@impact %>% filter(impact==impact_type)
  impact_with_obs$observed <- NULL
  impact_with_obs <- merge(impact_with_obs, observed_data, by.x=c('polygon_name', 'impact'), by.y=c('name', 'impact'), all.x=T)
  if(add_zeros) {impact_with_obs$observed[which(is.na(impact_with_obs$observed))] = 0}
  return(impact_with_obs)
}

plot_AggImpacts <- function(impact_with_obs, regions=NULL){
  #For each admin region, compares a histogram of the posterior predictive samples to the real observation
  sampled_cols <- grep('sampled.',colnames(impact_with_obs))
  if (is.null(regions)){
    regions_i = order(apply(impact_with_obs[, sampled_cols],1,function(x) mean(as.numeric(x))), decreasing=T)
    # Filter to GADM Level 2:
    regions_i <- regions_i[unlist(lapply(gregexpr(",", impact_with_obs$polygon_name[regions_i]), function(x) length(x) > 1))]
    regions_i <- regions_i[1:12]
  } else {
    stop('add some way to select named regions')
  }
  
  plots_list <- list()
  j <- 1
  for(i in regions_i){
    plots_list[[j]] <- ggplot(data.frame(sampled=as.numeric(impact_with_obs[i,sampled_cols])), aes(x = sampled))  +
      geom_histogram(aes(y= ..density..),color = "black", alpha = 0.7) +
      ggtitle(ODDyWithAggImpact@polygons[[impact_with_obs[i, 'polygon']]]$name) + scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), 
                                                                                                 breaks=c(0,10,100, 1000,10000), limits=c(-1,10000))
      #labs(x='', y='') + 
      #scale_fill_discrete("Legend", values = dd.col) + 
      #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    plots_list[[j]] = plots_list[[j]] + geom_vline(xintercept = impact_with_obs[i,'observed'], col='red')
    #scale_x_log10()+
    j <- j+1
  }
  do.call(grid.arrange,plots_list)
}

plot_correlated_impacts <- function(impact_with_obs, i, j, impact_type='mortality', baseR=F){
  # ODDyWithAggImpact is output of aggregateSampledImpact(sampleImpactODD(ODDy, AlgoResults))
  
  combined_samples <- t(impact_with_obs[c(i, j),grep('sampled.', colnames(impact_with_obs))])
  colnames(combined_samples) <- c(impact_with_obs$polygon_name[i], impact_with_obs$polygon_name[j])
  combined_samples %<>% as.data.frame()
  combined_samples <- rbind(combined_samples, observed= c(impact_with_obs[i, 'observed'],  impact_with_obs[j, 'observed']))
  if(baseR){
    cols <- rep('black', NROW(combined_samples)-1)
    cols <- c(cols, 'red')
    return(plot(as.numeric(combined_samples[,1]), as.numeric(combined_samples[,2]), pch=19, xlab=impact_with_obs$polygon_name[i], ylab=impact_with_obs$polygon_name[j], col=cols))
  }
  return(ggplot(combined_samples[1:(NROW(combined_samples)-1),], aes(x=as.numeric(.data[[region1]]), y=as.numeric(.data[[region2]]))) + geom_point()) #+ 
}

#Example of adding observation data and plotting comparison to posterior predictive.
#Unfortunately, can't automate the data retrieval as tends to come from different sources.
addObsData_MAR20230908 <- function(ODDyWithAggImpact){
  # ODDyWithAggImpact is output of aggregateSampledImpact(sampleImpactODD(ODDy, AlgoResults))
  
  observed_data_MAR20230908 <- data.frame(name=character(), impact=character(), observed=integer())
  observed_data_MAR20230908 %<>% add_row(name='Marrakech, Marrakech-Tensift-AlHaouz, Morocco', impact='mortality', observed = 15)
  observed_data_MAR20230908 %<>% add_row(name='AlHaouz, Marrakech-Tensift-AlHaouz, Morocco', impact='mortality', observed = 1684)
  observed_data_MAR20230908 %<>% add_row(name='Taroudannt, Souss-Massa-Draâ, Morocco', impact='mortality', observed = 980)
  observed_data_MAR20230908 %<>% add_row(name='Chichaoua, Marrakech-Tensift-AlHaouz, Morocco', impact='mortality', observed = 202)
  observed_data_MAR20230908 %<>% add_row(name='Azilal, Tadla-Azilal, Morocco', impact='mortality', observed = 11)
  observed_data_MAR20230908 %<>% add_row(name='Ouarzazate, Souss-Massa-Draâ, Morocco', impact='mortality', observed = 41)
  observed_data_MAR20230908 %<>% add_row(name='Morocco', impact='mortality', observed = 2946)
  observed_data_MAR20230908 %<>% add_row(name='TOTAL', impact='mortality', observed = 2946)
  
  impact_with_obs <- addObs(ODDyWithAggImpact, observed_data_MAR20230908, 'mortality', add_zeros=T )
  plot_AggImpacts(impact_with_obs)

  #quants <- apply(impact_with_obs[,c(grep('observed', names(impact_with_obs)),grep('sampled.', names(impact_with_obs)))], 1, sample_quant)
  #cbind(impact_with_obs$polygon_name, impact_with_obs$observed, quants)
  # grep(',',impact_with_obs$polygon_name[2])
  # sum(sapply(gregexpr(",", impact_with_obs$polygon_name), function(x) sum(x > 0))==2)
  # 
  # plot(as.numeric(ODDyWithAggImpact@impact[75,grep('sampled.', names(ODDyWithAggImpact@impact))]), as.numeric(ODDyWithAggImpact@impact[76,grep('sampled.', names(ODDyWithAggImpact@impact))]))
  # 
  # plot(log(as.numeric(ODDyWithAggImpact@impact[75,grep('sampled.', names(ODDyWithAggImpact@impact))])+10), log(as.numeric(ODDyWithAggImpact@impact[76,grep('sampled.', names(ODDyWithAggImpact@impact))])+10))
  # points(log(146000+10),log(2960+10), col='red', pch=19)
  
  # which(impact_with_obs$polygon==43 & impact_with_obs$impact=='mortality')
  # which(impact_with_obs$polygon==34 & impact_with_obs$impact=='mortality')
  # plot(as.numeric(impact_with_obs[42,grep('sampled.', names(impact_with_obs))]), as.numeric(impact_with_obs[30,grep('sampled.', names(impact_with_obs))]))
  # points(impact_with_obs[42,grep('observed', names(impact_with_obs))], impact_with_obs[30,grep('observed', names(impact_with_obs))], col='red')
  # par(mfrow=c(5,5))
  # for (i in 2:7){
  #   for (j in 1:(i-1)){
  #     plot_correlated_impacts(impact_with_obs, i, j, baseR=T)
  #     #points(x=ODDyWithAggImpact@impact$observed[i], y=ODDyWithAggImpact@impact$observed[j], col='red', pch=19)
  #   }
  # }
  # par(mfrow=c(1,1))
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



  
  