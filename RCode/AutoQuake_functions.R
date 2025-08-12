#
# -------------- AutoQuake_functions.R ----------------------------
#
# Generates ODD objects for new events and issues impact predictions
# Main Functions:
# -- prepareODD() takes USGS id OR country and date and creates an ODD object for the event
# -- posteriorImpactPred() takes the resulting ODD object and posterior parameter samples to issue predictions

# Examples:
# ODDyAgg <- prepareODD(iso3='MAR', sdate='2023-09-08') #readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908MAR_-1_AggLevel5')
# PosteriorImpactPred(ODDyAgg)
#
# ODDyAgg <- prepareODD(iso3='AFG', sdate='2023-10-07')
# PosteriorImpactPred(ODDyAgg)
# 
# ODDyAgg <- prepareODD(iso3='PHL', sdate='2013-10-15')
# PosteriorImpactPred(ODDyAgg)

# ODDyAgg <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908MAR_-1_AggLevel5')
# PosteriorImpactPred(ODDyAgg, 
#                     AlgoResultsFilename='HPC/mcmc_2025-02-27_162932.543295_MCMC_BDScore.05nocorr_LR40_M100_Npart1000NovAgg5_RandomFieldThree',
#                     folder_write='IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/')

input<-list(iso3='MMR',
            sdate = '2025-03-27',
            fdate = '2025-03-29',
            USGSid = 'us7000pn9s')

input<-list(
  sdate=as.Date("2019-12-13"), # "YYYY-MM-DD"
  fdate=as.Date("2019-12-17"), # "YYYY-MM-DD"
  iso3="PHL", # Country code in ISO-3C form
  datadir=dir, # Location of the main folder to access the data 
  plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

# Or extract the data purely based on the USGS id number
input<-list(USGSid="usp000huvq",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

prepareODD <- function(dir, input, getGADMregions=T, folder_write='IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/'){
  
  #if directory to save ODD object doesn't exist, then create one:
  write_loc <- paste0(dir, folder_write)
  if (!dir.exists(write_loc)) {
    dir.create(write_loc, recursive = TRUE)
  }
  if (!dir.exists(paste0(write_loc,'HAZARDobjects'))){
    dir.create(paste0(write_loc,'HAZARDobjects'))
  }
  if (!dir.exists(paste0(write_loc,'ODDobjects'))){
    dir.create(paste0(write_loc,'ODDobjects'))
  } 
  if (!dir.exists(paste0(write_loc,'tmp_polygons'))){
    dir.create(paste0(write_loc,'tmp_polygons'))
  }
  if (!dir.exists(paste0(write_loc,'ODDobjectsWithImpact'))){
    dir.create(paste0(write_loc,'ODDobjectsWithImpact'))
  }
  if (!dir.exists(paste0(write_loc,'ODDobjectsWithImpact/sampled_full'))){
    dir.create(paste0(write_loc,'ODDobjectsWithImpact/sampled_full'))
  }
  
  #getGADM
  if (!is.null(input$sdate)){ input$sdate %<>% as.Date()}
  if (!is.null(input$fdate)){ input$fdate %<>% as.Date()}
  
  #create hazard object:
  print('Downloading shakemap(s) from USGS')
  lhazSDF <- extractHAZARD(input$USGSid, input$iso3, input$bbox, input$sdate, input$fdate)
  
  # check that hazard seems ok:
  # plotHAZARD(lhazSDF) 
  
  if (is.null(input$iso3)){
    #find most common country to name the ODD object
    coords = xyFromCell(lhazSDF[[2]], 1:ncell(lhazSDF[[2]]))
    input$iso3 = na.omit(names(sort(table(coords2country( coords)), decreasing=T)))[1]
  } 
  namer<-paste0(lhazSDF$hazard_info$hazard, str_remove_all(as.character.Date(min(lhazSDF$hazard_info$eventdates)),"-"),input$iso3, "_-1")
  saveHAZ(lhazSDF, paste0(dir, folder_write, 'HAZARDobjects/', namer))
  
  #create ODD object:
  print('Creating full ODD object using hazard')
  ODDy <- createODD(lhazSDF, agg_level = 5)
  saveODD(ODDy, paste0(dir, folder_write, 'ODDobjects/', namer))
  
  if (getGADMregions){
    #add polygons for each administrative region
    print('Dividing spatial grid-cells into administrative regions')
    gadm_regions <- find_subnat_regions(ODDy, -1)
    polygons_list <- create_polygons_list(gadm_regions)
    saveRDS(list(gadm_regions=gadm_regions, polygons_list=polygons_list), paste0(dir, folder_write, 'tmp_polygons/', namer, 'polygons')) # save intermediary objects just in case something crashes
    ODDy <- addODDPolygons(ODDy, polygons_list)
    saveODD(ODDy, paste0(dir, folder_write, 'ODDobjects/', namer)) # again, save frequently just in case something crashes
    ODDy <- addPolygonsToImpact(ODDy)
    saveODD(ODDy, paste0(dir, folder_write, 'ODDobjects/', namer))
  } else {
    ODDy@impact= data.frame(
      iso3='TOT',
      sdate=ODDy@hazdates[1],
      impact=c('mortality', 'displacement', 'buildDam'),
      observed=NA,
      qualifier=NA,
      inferred=F,
      MAR=F,
      build_type=NA,
      polygon=1
    )
    ODDy@polygons = list(list(name='TOTAL', indexes=1:ncell(ODDy), weights=rep(1, ncell(ODDy))))
    saveODD(ODDy, paste0(dir, folder_write, 'ODDobjects/', namer))
  }
  
  #plot to make sure all is ok:
  #plotODDy(ODDy)
  #plotODDy_GADM(ODDy)
  
  #aggregate ODD object:
  # AggLevel <- 5
  # ODDyAgg <- aggregateODDbyX(ODDy, AggLevel)
  # namer<-paste0(namer, '_AggLevel', AggLevel)
  #saveODD(ODDyAgg, paste0(dir, folder_write, "ODDAggobjects/",namer))
  
  return(list(ODDy=ODDy, namer=namer))
}

# ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908MAR_-1_AggLevel5')
# AlgoResultsFilename = 'HPC/mcmc_2025-02-27_162932.543295_MCMC_BDScore.05nocorr_LR40_M100_Npart1000NovAgg5_RandomFieldThree'
# ODDyWithImpact <- sampleImpactODD(ODDy, readRDS(paste0(dir, 'IIDIPUS_Results/', AlgoResultsFilename)))
# saveODD(ODDyWithImpact, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908MAR_-1_AggLevel5_WithImpact100_Mar25')

PosteriorImpactPred <- function(ODDy, 
                                AlgoResultsFilename='IIDIPUS_Results/mcmc_2025-04-10',#'HPC/abcsmc_2024-08-20_051627_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2',
                                folder_write='IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/', 
                                namer = namer,
                                N_samples = 100){
  #AlgoResultsFilename should direct to a file within .../IIDIPUS_Results/ that contains an AlgoResults object
  #The AlgoResults object must have been be fitted using the same model as the one currently used to issue predictions 
  #The aggregation of ODD objects used here should also match that of the ODD objects used to fit AlgoResults
  
  AlgoResults <- readRDS(paste0(dir, AlgoResultsFilename))
  #AlgoResults$Omega_sample_phys <- AlgoResults$Omega_sample_phys[, c(1:6, 10:22),]
  
  ODDy_with_sampledfull <- sampleImpactODD(ODDy, AlgoResults, N_samples = N_samples)
  ODDyWithImpact = ODDy_with_sampledfull$ODDy
  sampled_full = ODDy_with_sampledfull$sampled_full
  saveODD(ODDyWithImpact, paste0(dir,folder_write, 'ODDobjectsWithImpact/',namer, '_WithImpactSummaries'))
  saveRDS(sampled_full, paste0(dir,folder_write,'ODDobjectsWithImpact/sampled_full/', namer, '_fullSampledImpact'))
  
  ODDyWithAggImpact <- aggregateSampledImpact(ODDyWithImpact, sampled_full)
  saveODD(ODDyWithAggImpact, paste0(dir, folder_write, 'ODDobjectsWithImpact/',namer, '_WithImpact200'))
  
  return(list(sampled_full=sampled_full, ODDyWithImpact=ODDyWithAggImpact))
  
}

MakeODDPredPlots <- function(ODDy_with_sampled, input = input,  folder_write=folder_write, namer = namer){
  
  ODDy = ODDy_with_sampled$ODDyWithImpact
  sampled_full = ODDy_with_sampled$sampled_full
  
  #plotODDy_GADM(ODDyWithImpact, var='mortality.mean') # doesn't show up on map very well. Try to improve plot
  #plotODDy(ODDyWithImpact, var='displacement.mean')
  
  # polygons_list = readRDS(paste0(dir, folder_write, 'tmp/', namer, 'polygons'))$polygons_list
  # 
  # p_regional_impact = plot_GADM_impact_polygons(ODDy, polygons_list, gadm_level = 2, impact_type = 'mortality', summary = 'mean')
  # print(p_regional_impact)
  
  p_panel = plot_predictive_panel(ODDy, sampled_full)
  print(p_panel)
  
  # Define output folder
  output_dir <- paste0(dir, folder_write, input$plotdir)
  
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save plots
  #ggsave(filename = file.path(output_dir, namer, "_regional_impact_plot.png"), plot = p_regional_impact, width = 8, height = 6)
  ggsave(filename = paste0(output_dir, namer, "_predictive_panel_plot.png"), plot = p_panel, width = 20, height = 16)
  return(p_panel)
}

extractHAZARD <- function(USGSid = NULL, iso3 = NULL, bbox = NULL, sdate = NULL, fdate=NULL){
  #Need to specify either:
  #      - GDACS_id
  #      - iso3, sdate (can further specify bbox to constrain)
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

createODD <- function(lhazSDF, eventid = -1, agg_level = 1){
  #Create the ODD object:
  #Set event_d = -1 for new events to be predicted
  
  miniDamSimplified = data.frame(iso3=NA, sdate=min(lhazSDF[[1]]$eventdates), fdate=max(lhazSDF[[1]]$eventdates), 
             eventid=eventid, hazard=lhazSDF[[1]]$hazard)
  
  ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDamSimplified, agg_level=agg_level),error=function(e) NULL)
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
  sf_use_s2(FALSE)
  iso3_unique <- unique(ODDy$ISO3C)$ISO3C[!is.na(unique(ODDy$ISO3C)$ISO3C)]
  gadm_regions_1 <- data.frame('Country'=character(), 'Region'=character(), 'Subregion'=character())
  gadm_regions_2 <- data.frame('Country'=character(), 'Region'=character(), 'Subregion'=character())
  for (gadm_level in gadm_levels){
    pixels_unallocated <- which(values(!is.na(ODDy$ISO3C)))
    for (iso3 in iso3_unique){
      print(paste(iso3, gadm_level))
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
        overlap <- (max(r_poly@bbox[1,]) >= ODDy@hazinfo$bbox[1] &
                      min(r_poly@bbox[1,]) <= ODDy@hazinfo$bbox[3] &
                      max(r_poly@bbox[2,]) >= ODDy@hazinfo$bbox[2] &
                      min(r_poly@bbox[2,]) <= ODDy@hazinfo$bbox[4])
        if (!overlap) next
        
        if (length(pixels_unallocated)==0) next
        
        crds <- crds(ODDy, na.rm=F)
        
        # spdf_pixels_unallocated <- SpatialPointsDataFrame(coords=crds[pixels_unallocated,, drop=F],
        #                                                   data=ODDy@data[pixels_unallocated,1:2, drop=F], #this is arbitrary, function just seems to require data
        #                                                   proj4string=r_poly@proj4string)
        # spdf_pixels_unallocated$id = 1:length(pixels_unallocated)
        # contained <- tryCatch(terra::intersect(r_poly, spdf_pixels_unallocated)$id, error=function(e) c())
        
        sf_pixels_unallocated <- st_as_sf(data.frame(xyFromCell(ODDy, pixels_unallocated)), coords = c("x", "y"), crs = crs(ODDy))
        contained <- c(st_intersects(sf_pixels_unallocated, st_make_valid(st_as_sf(r_poly)), sparse=F)) 
        
        if (sum(contained)>0){
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
  sf_use_s2(TRUE)
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
  if ('nBuildings' %in% names(ODDy)){impact_types <- c(impact_types, 'buildDam')}
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
  
  iso3_unique <- unique(ODDy$ISO3C)$ISO3C[!is.na(unique(ODDy$ISO3C)$ISO3C)]
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
                                                 indexes=1:ncell(ODDy), weights=rep(1, ncell(ODDy)))
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
                                       indexes=which(values(ODDy$ISO3C==iso3_unique[j])), weights=rep(1, sum(values(ODDy$ISO3C==iso3_unique[j]),na.rm=T)))
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
    ODDy@polygons[[tot_poly]] = list(name='TOTAL',indexes=1:ncell(ODDy), weights=rep(1, ncell(ODDy)))
  }
  return(ODDy)
}

sampleImpactODD <- function(ODDy, AlgoResults, N_samples = 100){
  #Sample from the posterior predictive distribution for the event's impact
  AlgoParamsPostPredict <- list(Np = 1, m_CRPS=1, cores=1, NestedCores=1)
  params_sampled <- rep(NA, N_samples)
  sampled <- array(NA, dim=c(ncell(ODDy),3, N_samples))
  for (i in 1:N_samples){
    print(paste('Sample', i))
    if (length(dim(AlgoResults$Omega_sample_phys))==2){ # MCMC
      param_sample <- sample(1:min(1000, which(is.na(AlgoResults$Omega_sample_phys[1,]))[1]-1),1) #sample within the last 1000 iterations
      proposed = AlgoResults$Omega_sample_phys[,which(is.na(AlgoResults$Omega_sample_phys[1,]))[1]-param_sample] %>% relist(skeleton=Model$skeleton) %>% addTransfParams() #MCMC
      # while(proposed$Lambda2$nu > 12.5){
      #   print('Reject')
      #   param_sample <- sample(1:1000,1) #sample within the last 1000 iterations
      #   proposed = AlgoResults$Omega_sample_phys[,which(is.na(AlgoResults$Omega_sample_phys[1,]))[1]-param_sample] %>% relist(skeleton=Model$skeleton) %>% addTransfParams() #MCMC
      # }
    } else { #ABC-SMC
      AlgoResults %<>% addAlgoParams()
      param_sample <- sample(1:AlgoResults$Npart, 1, replace=T, prob=AlgoResults$W[,AlgoResults$s_finish])
      params_sampled[i] <- param_sample
      proposed = AlgoResults$Omega_sample_phys[param_sample,,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton) %>% addTransfParams()
    }
    #param_sample <- sample(1:2000,1) #MCMC
    #param_sample <- sample(1:AlgoResults$Npart, 1, replace=T, prob=AlgoResults$W[,AlgoResults$s_finish])
    #params_sampled[i] <- param_sample
    #proposed = AlgoResults$Omega_sample_phys[,which(is.na(AlgoResults$Omega_sample_phys[1,]))[1]-param_sample] %>% relist(skeleton=Model$skeleton) %>% addTransfParams() #MCMC
    #proposed = AlgoResults$Omega_sample_phys[param_sample,,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton) %>% addTransfParams()
    impact <- DispX(ODDy, proposed, Model$center, Method = AlgoParamsPostPredict, output='SampledFull')
    sampled[,,i] <- impact
  }
  
  impact_sampled_summaries = data.frame(row.names = 1:nrow(sampled))
  for (i in 1:ifelse('nBuildings' %in% names(ODDy), 3,2)){
    # i = 1: Displacement; i = 2: Mortality; i = 3: BuildDam. 
    # Add a column to ODDy@data for each posterior predictive impact sample, the mean, median, and the 5th and 95th percentiles
    impact_types = c('displacement', 'mortality', 'buildDam')
    #impact_sampled_full <- sampled[,i,]
    #colnames(impact_sampled) = paste0(impact_types[i],'.s', 1:N_samples)
    #impact_sampled %<>% as.data.frame()
    impact_sampled_summaries[,paste0(impact_types[i], '.mean')] = apply(sampled[,i,], 1, mean )
    impact_sampled_summaries[,paste0(impact_types[i], '.median')] = apply(sampled[,i,], 1, median )
    impact_sampled_summaries[,paste0(impact_types[i], '.q05')] = apply(sampled[,i,], 1, quantile,0.05 )
    impact_sampled_summaries[,paste0(impact_types[i], '.q95')] = apply(sampled[,i,], 1, quantile, 0.95 )
    impact_sampled_summaries[,paste0(impact_types[i], '.q10')] = apply(sampled[,i,], 1, quantile,0.1 )
    impact_sampled_summaries[,paste0(impact_types[i], '.q90')] = apply(sampled[,i,], 1, quantile, 0.9 )
    ODDy[[names(impact_sampled_summaries)]] = impact_sampled_summaries
  } 
  return(list(ODDy=ODDy, sampled_full=sampled))
}

aggregateSampledImpact <- function(ODDyWithImpact, sampled_full){
  #Aggregate the sampled impact to add posterior predictive samples to ODDy@impact:
  
  n_samples <- dim(sampled_full)[3] #length(grep(paste0('mortality', '.s'),names(ODDyWithImpact)))
  sampled_df <- data.frame(polygon=integer(), impact=character(), polygon_name = character())
  for (i in 1:n_samples){ sampled_df[,paste0('sampled.', i)] = integer()}
  
  # i = 1: Displacement; i = 2: Mortality; i = 3: BuildDam
  impact_types = c('displacement', 'mortality', 'buildDam')
  if(nrow(ODDyWithImpact@impact)==0){
    ODDyWithImpact@impact= data.frame(
      iso3='TOT',
      sdate=ODDyWithImpact@hazdates[1],
      impact=c('mortality', 'displacement', 'buildDam'),
      observed=NA,
      qualifier=NA,
      inferred=F,
      MAR=F,
      build_type=NA,
      polygon=1
    )
    ODDyWithImpact@polygons = list(list(name='TOTAL', indexes=1:ncell(ODDyWithImpact), weights=rep(1, ncell(ODDyWithImpact))))
  }
  for (i in 1:ifelse('nBuildings' %in% names(ODDyWithImpact), 3,2)){
    impact_sampled <- sampled_full[,i,]#as.data.frame(ODDyWithImpact[[names(ODDyWithImpact)[grep(paste0(impact_types[i], '.s'),names(ODDyWithImpact))]]], na.rm=F)
    for (j in 1:length(ODDyWithImpact@polygons)){
      sampled_df[nrow(sampled_df)+1,] = c(j, impact_types[i], ODDyWithImpact@polygons[[j]]$name, colSums(impact_sampled[ODDyWithImpact@polygons[[j]]$indexes,, drop=F]))
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

MAR_plots <- function(ODDyWithAggImpact){
  p1 <- plotODDy_GADM(ODDyWithAggImpact, var='mortality.mean', haz_legend=F, var_legend=T) + labs(fill = "Mortality Mean")
  p2 <- plotODDy_GADM(ODDyWithAggImpact, var='mortality.q05', haz_legend=F, var_legend=T) + labs(fill = "Mortality 5%")
  p3 <- plotODDy_GADM(ODDyWithAggImpact, var='mortality.q95', haz_legend=F, var_legend=T ) + labs(fill = "Mortality 95%")
  p4 <- plotODDy_GADM(ODDyWithAggImpact, var='displacement.mean', haz_legend=F, var_legend=T, log_legend=T) + labs(fill = "Displacement Mean")
  p5 <- plotODDy_GADM(ODDyWithAggImpact, var='displacement.q05', haz_legend=F, var_legend=T) + labs(fill = "Displacement 5%")
  p6 <- plotODDy_GADM(ODDyWithAggImpact, var='displacement.q95', haz_legend=F, var_legend=T, log_legend=T) + labs(fill = "Displacement 95%")
  
  plot_with_haz_legend <- plotODDy_GADM(ODDyWithAggImpact, var='mortality.mean', haz_legend=T, var_legend=F) +
                            theme(legend.position="bottom", legend.key.size = unit(15, 'pt'),
                                  legend.key.width= unit(30, 'pt'),
                                  legend.text = element_text(margin = margin(t=2, r = 10)))
  haz_legend <- get_plot_component(plot_with_haz_legend, 'guide-box', return_all=T)[[3]]
  
  plots_all <- plot_grid( p1, p4, nrow=1,  align='v', labels = c('(a)', '(b)'),
                          label_fontface = "plain",
                          label_fontfamily = "Liberation Serif")
  plots_with_legend_mean <- plot_grid(plots_all, haz_legend, nrow=2, rel_heights=c(1,0.1))
  plots_with_legend_mean
  
  #MAR_Mean.pdf, 5 x 12
  
  plots_all <- plot_grid( p2, p3,p5, p6, nrow=2,  align='v', labels = c('(a)', '(b)', '(c)', '(d)'),
                          label_fontface = "plain",
                          label_fontfamily = "Liberation Serif")
  plots_with_legend <- plot_grid(plots_all, haz_legend, nrow=2, rel_heights=c(1,0.1))
  plots_with_legend
  #MAR_quantiles.pdf, 10 x 12
 
  #plotODDy_GADM(ODDyWithAggImpact, var='buildDam.mean', haz_legend=T, var_legend=T) + labs(fill = "Building Damage Mean")
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
  observed_data_MAR20230908 %<>% add_row(name='TOTAL', impact='displacement', observed = 500000)
  
  impact_with_obs <- addObs(ODDyWithAggImpact, observed_data_MAR20230908, 'mortality', add_zeros=T )
  impact_with_obs_disp <- addObs(ODDyWithAggImpact, observed_data_MAR20230908, 'displacement', add_zeros=F) %>% filter(!is.na(observed))
  
  plot_AggImpacts(impact_with_obs)
  
  sampled_cols <- grep('sampled.', names(impact_with_obs ))
  impact_with_obs$sampled.median = apply(impact_with_obs[,sampled_cols ], 1, function(x) median(as.numeric(x)))
  impact_with_obs$sampled.q05 = apply(impact_with_obs[,sampled_cols ], 1, function(x) quantile(as.numeric(x), 0.025))
  impact_with_obs$sampled.q95 = apply(impact_with_obs[,sampled_cols ], 1, function(x) quantile(as.numeric(x), 0.975))
  impact_with_obs[impact_with_obs$sampled.q95>0, c('polygon_name', 'observed', 'sampled.q05', 'sampled.median', 'sampled.q95')]
  
  #plot_correlated_impacts(rbind(impact_with_obs, impact_with_obs_disp), 35, 33, baseR=T)
  #plot_correlated_impacts(impact_with_obs, 4, 33, baseR=T)
  # par(mfrow=c(5,5))
  # for (i in 2:7){
  #   for (j in 1:(i-1)){
  #     plot_correlated_impacts(impact_with_obs, i, j, baseR=T)
  #     #points(x=ODDyWithAggImpact@impact$observed[i], y=ODDyWithAggImpact@impact$observed[j], col='red', pch=19)
  #   }
  # }
  # par(mfrow=c(1,1))
}

plotWithoutObs = function(ODDyWithAggImpact){
  #ODDyWithAggImpact <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20250328MMR_-1_WithImpact500')
  impact_type = 'mortality'
  polyName1 = 'Thailand'
  polyName2 = 'Myanmar (Burma)'
  
  impact_filt = ODDyWithAggImpact@impact[which(ODDyWithAggImpact@impact$impact==impact_type),]
  poly_row1 = which(impact_filt$polygon_name == polyName1)
  poly_row2 = which(impact_filt$polygon_name == polyName2)
  
  plot(as.numeric(impact_filt[poly_row1, grep('sampled', names(impact_filt))]), 
       as.numeric(impact_filt[poly_row2, grep('sampled', names(impact_filt))]),
       xlab=paste(polyName1, impact_type),
       ylab=paste(polyName2, impact_type), main='500 Posterior Samples')
  points(80, 3000,col='red', pch=19)
  
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#==============================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#====== Functions for plotting predictions and comparing to PAGER =============
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#==============================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ODDyAggWithImpact <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908MAR_-1_AggLevel5_WithImpact100_Mar25')

plotODDyImpactQuantile <- function(ODDy, impact_type, quantile, max_break = 500, gadm_level=1){
  
  impact_type_samples = paste0(impact_type, '.q')

  plot_df <- as.data.frame(ODDy, xy=T, na.rm=F)
  #plot_df[,grep(impact_type_samples, names(plot_df))] = t(apply(as.matrix(plot_df[, grep(impact_type_samples, names(plot_df)), drop = FALSE]), MARGIN=1, sort))
  
  var = paste0(impact_type, '.q', quantile)
  if (quantile == 50){
    var = paste0(impact_type, '.median')
  } else if (quantile == 5){
    var = paste0(impact_type, '.q05')
  }
  plot_df[which(is.na(plot_df$ISO3C)), var] <- NA
  names(plot_df)[which(names(plot_df)=='x')] = 'Longitude'
  names(plot_df)[which(names(plot_df)=='y')] = 'Latitude'
  
  bbox <- matrix(ODDy@hazinfo$bbox, nrow=2)
  #gadm_iso <- getData("GADM", country="NZL", level=2)
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
  
  p <- ggplot() + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
  
  p <- p + geom_map(map=gadm_map, data=gadm_map, aes(map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
  p <- p + coord_map() + geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='white', color='black')
  
  if (impact_type == 'mortality'){
    breaks = c(-Inf, 1, 3, 5, 10, 20, 50, 100, 500, Inf)
    if (max_break < 500 & max_break > 100){
      breaks = c(-Inf, 1, 3, 5, 10, 20, 50, 100, Inf)
    }
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
  } else {
    #breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, 10000, 50000, Inf)
    if (max_break < 50000 & max_break > 10000){
      breaks = c(-Inf, 1, 100, 500, 1000, 5000, 10000, Inf)
    }
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
  }
  
  bin_colors <- viridis::viridis(length(levels(plot_df$bin_category)))
  
  
  # Modify the plot
  map_plot <- p + 
    geom_raster(data = plot_df, aes(x = Longitude, y = Latitude, fill = bin_category), alpha = 0.75) + 
    coord_equal() +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal() + 
    theme(
      axis.title = element_text(family = "Liberation Serif", size = 12),  
      legend.text = element_text(family = "Liberation Serif", size = 11),
      legend.title = element_text(family = "Liberation Serif", size = 12)
    ) + 
    geom_contour(data = plot_df, aes(Longitude, Latitude, z = hazMean1, colour = ..level..),
                 alpha = 0.7, lwd = 0.8) +
    scale_color_gradientn(colors = c("transparent", "#fc9272", "#ef3b2c")) +
    labs(colour = "Hazard Intensity") + 
    scale_fill_manual(name = paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type))),#paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type)),' ', quantile,'% Quantile'), 
                      values = bin_colors, na.value = "transparent") # Discrete color scale
  
  return(map_plot)
  
  
  # pseudo_log_sigma = 100
  # max_var = max(plot_df[,var], na.rm=T)
  # if (max_var > 100000){
  #   breakings = c(0, 1000, 10000, 100000, 1000000)
  # } else if (max_var > 10000){
  #   breakings = c(0, 1000, 10000)
  # } else if (max_var > 1000){
  #   breakings = c(0, 500, 1000)
  # } else if (max_var > 100){
  #   breakings = c(0, 20,50, 90, floor(max_var/10)*10)
  #   pseudo_log_sigma=10
  # } else if (max_var > 50){
  #   breakings = c(0, 10, 50)
  # } else if (max_var > 20){
  #   breakings = c(0, 2, 5, 10, 20, floor(max_var/10)*10)
  #   pseudo_log_sigma = 1
  # } else if (max_var > 10){
  #   breakings = c(0, 5, 10)
  # } else if (max_var > 5){
  #   breakings = c(0, 5, 10)
  # } else if (max_var == 0){
  #   breakings = c(0,1)
  # }
  # 
  # p <- p + 
  #   geom_raster(dat=plot_df, aes(x=Longitude, y=Latitude, fill=!!sym(var)), alpha=0.75) + 
  #   coord_equal() +
  #   scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  #   theme_minimal() + 
  #   theme(
  #     axis.title = element_text(family = "Liberation Serif", size=12),  
  #     legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
  #     legend.title = element_text(family = "Liberation Serif", size=12)
  #   ) + 
  #   geom_contour(dat=plot_df, aes(Longitude,Latitude,z=hazMean1,colour=..level..),
  #                alpha=0.7, lwd=0.8) +
  #   scale_color_gradientn(colors = c("transparent","#fc9272", "#ef3b2c"))  + labs(colour = "Hazard Intensity      ") + 
  #   scale_fill_viridis( trans = scales::pseudo_log_trans(base = 10, sigma = pseudo_log_sigma), breaks=breakings, labels = function(x) scales::comma(x))
  # p
  # if (all(plot_df[, var]==0, na.rm=T)){
  #         p <- p + scale_fill_viridis_c( limits = c(0, 1),breaks=c(0,1), labels=c(0,1))
  # } else {
  #         p <- p + scale_fill_viridis_c()
  # }
  # p
}

#plot_df$Population[which(is.na(plot_df$Population))] = 0
# breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, 10000, 50000, Inf)
# plot_df$bin_category <- cut(
#   plot_df[['Population']], 
#   breaks = breaks, 
#   labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
#   right = FALSE
# )
# bin_colors <- viridis::viridis(length(levels(plot_df$bin_category)))
# p + 
#   geom_raster(data = plot_df, aes(x = Longitude, y = Latitude, fill = Population), alpha = 0.75) + 
#   coord_equal() +
#   scale_x_continuous(expand = c(0,0)) + 
#   scale_y_continuous(expand = c(0,0)) +
#   theme_minimal() + 
#   theme(
#     axis.title = element_text(family = "Liberation Serif", size = 12),  
#     legend.text = element_text(family = "Liberation Serif", size = 11),
#     legend.title = element_text(family = "Liberation Serif", size = 12)
#   ) + 
#   geom_contour(data = plot_df, aes(Longitude, Latitude, z = hazMean1, colour = ..level..),
#                alpha = 0.7, lwd = 0.8) +
#   scale_color_gradientn(colors = c("transparent", "#fc9272", "#ef3b2c")) +
#   labs(colour = "Hazard Intensity") + 
#   scale_fill_viridis(trans = "log10", labels = function(x) sprintf("%.0f", x))
#   scale_fill_manual(name = paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type))),#paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type)),' ', quantile,'% Quantile'), 
#                     values = bin_colors, na.value = "transparent") # Discrete color scale # Discrete color scale

#plotODDyImpactQuantile(ODDyAggWithImpact, impact_type='mortality', quantile=5, max_break = 500)

#plot mortality histogram and compare to PAGER:
#ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908_MAR_-1AggLevel5_WithImpact100_Dec24')


ODDRIN_vs_PAGER_hists <- function(ODDy, sampled_full, impact_type = 'mortality'){
  #impact_type_samples = paste0(impact_type, '.s')
  
  plot_df <- as.data.frame(ODDy, xy=T, na.rm=F)
  #tot_impact = colSums(as.matrix(plot_df[,grep(impact_type_samples, names(plot_df))]))
  impact_i = which(c('displacement', 'mortality', 'buildDam') == impact_type)
  tot_impact = colSums(sampled_full[,impact_i,])
  folderin_haz <- paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/HAZARDobjects/')
  ufiles_haz <- na.omit(list.files(path=folderin_haz,pattern=Model$haz,recursive = T,ignore.case = T))
  

  file_match <- grep(gsub("-", "", as.character(ODDy@hazinfo$eventdates)),  ufiles_haz, value = TRUE)[1]
  if (length(file_match)==0){ stop()}
  HAZy <- readHAZ(paste0(folderin_haz, file_match ))
  which.max.mmi <- which.max(sapply(HAZy[2:length(HAZy)], function(x) max(values(x$mean), na.rm=T)))
  if (!all(HAZy$hazard_info$bbox == ODDy@hazinfo$bbox)){stop('Mismatched HAZARD and ODD objects')}

  binned_preds = data.frame(bin_lower=numeric(),
                            bin_upper=numeric(), 
                            ODDRIN_prob = numeric(),
                            PAGER_prob = numeric())
  
  bins = HAZy$hazard_info$alertfull[[which.max.mmi]]$bins
  
  if (is.null(bins)){
    bins = list(list(min=0, max=1, probability=0),
                list(min=1, max=10, probability=0),
                list(min=10, max=100, probability=0),
                list(min=100, max=1000, probability=0),
                list(min=1000, max=10000, probability=0),
                list(min=10000, max=100000, probability=0),
                list(min=100000, max=10000000, probability=0))
    prob_mult = 1
  } else {
    prob_mult <- 1/sum(unlist(lapply(bins, function(x) return(x$probability)))) # some sum to 1 and some to 100
  }
  
  for (j in 1:7){
    oddrin_preds <- tot_impact
    oddrin_prob = mean((oddrin_preds >= bins[[j]]$min) & (oddrin_preds <= bins[[j]]$max))
    
    binned_preds %<>% add_row(
      bin_lower = bins[[j]]$min,
      bin_upper = bins[[j]]$max,
      ODDRIN_prob =  mean((oddrin_preds >= bins[[j]]$min) & (oddrin_preds <= bins[[j]]$max)),
      PAGER_prob = bins[[j]]$probability * prob_mult
    )
  }
  binned_preds$bin_upper[1] = 0
  
  binned_preds$color <- cut(binned_preds$bin_lower, 
                            breaks = c(-Inf, 1, 100, 1000, Inf), 
                            labels = c("green", "yellow", "orange", "red"), 
                            right = FALSE)
  
  # Create the histogram
  hist_ODDRIN <- ggplot(binned_preds, aes(x = factor(paste(bin_lower, ' - ', bin_upper)), y = ODDRIN_prob, fill = color)) +
    geom_bar(width=1,stat = "identity", color = "black") +  # Black border for clarity
    scale_fill_identity() +  # Use pre-defined colors
    labs(x = "Fatalities", y = "Probability", title=paste('ODDRIN', impact_type)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  
  hist_PAGER <- ggplot(binned_preds, aes(x = factor(paste(bin_lower, ' - ', bin_upper)), y = PAGER_prob, fill = color)) +
    geom_bar(width=1,stat = "identity", color = "black") +  # Black border for clarity
    scale_fill_identity() +  # Use pre-defined colors
    labs(x = "Fatalities", y = "Probability", title=paste('PAGER', impact_type)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  
  return(list(hist_ODDRIN, hist_PAGER))
  
  #plot_df[,grep(impact_type_samples, names(plot_df))] = t(apply(as.matrix(plot_df[, grep(impact_type_samples, names(plot_df)), drop = FALSE]), MARGIN=1, sort))
  
}

ODDRIN_hist <- function(sampled_full, impact_type = 'displacement'){
  #impact_type_samples = paste0(impact_type, '.s')
  
  #plot_df <- as.data.frame(ODDy, xy=T, na.rm=F)
  impact_order = c('displacement', 'mortality', 'buildDam')
  i_impact = which(impact_order == impact_type)
  tot_impact = colSums(sampled_full[,i_impact,])#as.matrix(plot_df[,grep(impact_type_samples, names(plot_df))]))
  
  if (impact_type == 'displacement'){
    binned_preds = data.frame(
      bin_lower = c(0, 10, 100, 1000, 10000, 100000, 1000000),
      bin_upper = c(9, 99, 999, 9999, 99999, 999999, Inf),
      ODDRIN_prob = NA)
  } else {
    binned_preds = data.frame(
      bin_lower = c(0, 1, 10, 100, 1000, 10000, 100000),
      bin_upper = c(0, 9, 99, 999, 9999, 99999, Inf),
      ODDRIN_prob = NA)
  }
  
  for (j in 1:nrow(binned_preds)){
    oddrin_preds <- tot_impact
    binned_preds$ODDRIN_prob[j] = mean((oddrin_preds >= binned_preds$bin_lower[j]) & (oddrin_preds <= binned_preds$bin_upper[j]))
  }
  #binned_preds$bin_upper[1] = 0
  
  if (impact_type == 'mortality'){
    binned_preds$color <- cut(binned_preds$bin_lower, 
                              breaks = c(-Inf, 1, 100, 1000, Inf), 
                              labels = c("green", "yellow", "orange", "red"), 
                              right = FALSE)
  } else {
    binned_preds$color <- cut(binned_preds$bin_lower, 
                              breaks = c(-Inf, 10, 1000, 100000, Inf), 
                              labels = c("green", "yellow", "orange", "red"), 
                              right = FALSE)
  }
  
  # Create the histogram
  # hist_ODDRIN <- ggplot(binned_preds, aes(x = factor(paste(bin_lower, ' - ', bin_upper)), y = ODDRIN_prob, fill = color)) +
  #   geom_bar(width=1,stat = "identity", color = "black") +  # Black border for clarity
  #   scale_fill_identity() +  # Use pre-defined colors
  #   labs(x = paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type))), y = "Probability", title=paste('ODDRIN', impact_type)) +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5),
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank()) 
  
  binned_preds$percent <- (binned_preds$ODDRIN_prob / sum(binned_preds$ODDRIN_prob)) * 100
  
  hist_ODDRIN <- ggplot(binned_preds, aes(x = factor(paste(bin_lower, ' - ', bin_upper)), 
                                          y = ODDRIN_prob, fill = color)) +
    geom_bar(width = 1, stat = "identity", color = "black") +  # Black border for clarity
    geom_text(aes(label = paste0(round(percent, 1), "%")), 
              vjust = -0.5, size = 5) +  # Position labels above bars
    scale_fill_identity() +  # Use pre-defined colors
    labs(x = paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type))), 
         y = "Probability", 
         title = paste('ODDRIN', impact_type)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + coord_cartesian(clip="off")  
  
  return(hist_ODDRIN)
  
}


#ODDyAggWithImpact <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908MAR_-1_AggLevel5_WithImpact100_Mar25')
plot_predictive_panel <- function(ODDyAggWithImpact, sampled_full){
  
  ODDRIN_mortpred.05 = plotODDyImpactQuantile(ODDyAggWithImpact, impact='mortality', quantile=5, max_break = 150) + ggtitle("Mortality 5% Quantile")
  ODDRIN_mortpred.50 = plotODDyImpactQuantile(ODDyAggWithImpact, impact='mortality', quantile=50, max_break = 150) + ggtitle("Mortality 50% Quantile")
  ODDRIN_mortpred.95 = plotODDyImpactQuantile(ODDyAggWithImpact, impact='mortality', quantile=95, max_break = 150) + ggtitle("Mortality 95% Quantile") 
  
  ODDRIN_disppred.05 = plotODDyImpactQuantile(ODDyAggWithImpact, impact='displacement', quantile=5, max_break=10001) + ggtitle("Displacement 5% Quantile")
  ODDRIN_disppred.50 = plotODDyImpactQuantile(ODDyAggWithImpact, impact='displacement', quantile=50, max_break=10001) + ggtitle("Displacement 50% Quantile")
  ODDRIN_disppred.95 = plotODDyImpactQuantile(ODDyAggWithImpact, impact='displacement', quantile=95, max_break=10001) + ggtitle("Displacement 95% Quantile") 
  
  legendmort <- get_plot_component(ODDRIN_mortpred.95, 'guide-box', return_all=T)[[1]]
  
  #grid.arrange(grid.arrange(ODDRIN_mortpred.05 + theme(legend.position="none"), ODDRIN_mortpred.50 + theme(legend.position="none"),
  #             ODDRIN_mortpred.95 + theme(legend.position="none"), legendmort, ncol=4, widths=c(1,1,1,0.35)))
  
  legenddisp <- get_plot_component(ODDRIN_disppred.95, 'guide-box', return_all=T)[[1]]
  
  p_mort_compare = ODDRIN_vs_PAGER_hists(ODDyAggWithImpact, sampled_full, impact='mortality')
  p_disp = ODDRIN_hist(sampled_full, impact='displacement')
  p_buildDam = ODDRIN_hist(sampled_full, impact='buildDam')
  
  histograms_mort = grid.arrange(
        p_mort_compare[[1]],
        p_mort_compare[[2]],
        nrow = 2, heights=c(0.5,0.5)
      )
  histograms_dispbuild = grid.arrange(
    p_disp,
    p_buildDam,
    nrow = 2, heights=c(0.5,0.5)
  )
  
  plot_panel = grid.arrange(grid.arrange(legendmort, ODDRIN_mortpred.05 + theme(legend.position="none"), ODDRIN_mortpred.50 + theme(legend.position="none"), 
                            ODDRIN_mortpred.95 + theme(legend.position="none"), histograms_mort,
                            legenddisp, ODDRIN_disppred.05 + theme(legend.position="none"), ODDRIN_disppred.50 + theme(legend.position="none"), 
                            ODDRIN_disppred.95 + theme(legend.position="none"), histograms_dispbuild,
                            ncol=5,nrow=2, widths=c(0.3,1,1,1,0.7)))
  
  return(plot_panel)
  

  # Arrange plots with relative heights
  #grid.arrange(p_mort_compare[[1]], p_mort_compare[[2]], nrow = 2)
  
  # grid.arrange(p_mort_compare[[1]], p_mort_compare[[2]], p_disp, ncol=3)
  # 
  # row1 <- grid.arrange(
  #   ODDRIN_pred.05 + theme(legend.position="none"),
  #   ODDRIN_pred.50 + theme(legend.position="none"),
  #   ODDRIN_pred.95 + theme(legend.position="none"),
  #   legend,
  #   ncol = 4, widths = c(1,1,1,0.4)
  # )
  # row2 <- grid.arrange(
  #   p_mort_compare[[1]],
  #   p_mort_compare[[2]],
  #   ncol = 2
  # )
  # grid.arrange(row1, row2, nrow = 2, heights=c(0.7, 0.4))

  # histogram_col = grid.arrange(
  #     p_mort_compare[[1]],
  #     p_mort_compare[[2]],
  #     nrow = 2, heights=c(0.5,0.5)
  #   )
  # grid.arrange(
  #     legend,
  #     ODDRIN_pred.05 + theme(legend.position="none"),
  #     ODDRIN_pred.50 + theme(legend.position="none"),
  #     ODDRIN_pred.95 + theme(legend.position="none"),
  #     histogram_col,
  #     ncol = 5, widths = c(0.35,1,1,1,0.5)
  #   )
  
  # xx <- p_mort_compare[[1]]
  # plot_grid(
  #   legend,
  #   ODDRIN_pred.05 + theme(legend.position="none"),
  #   ODDRIN_pred.50 + theme(legend.position="none"),
  #   ODDRIN_pred.95 + theme(legend.position="none"),
  #   p_mort_compare[[1]], ncol=5, rel_heights=c(1,1,1,1,0.2))
  
}

# plot_GADM_impact <- function(ODDy, gadm_level = 2, impact_type = 'mortality'){
#   #GADM_level 0 = national
#   #GADM_level 1 = admin level 1
#   #GADM_level 2 = admin level 2
#   
#   ODDy$polygon_id = 0
#   for (i in 1:length(ODDy@polygons)){
#     if (str_count(ODDy@polygons[[i]]$name, ',') != gadm_level | ODDy@polygons[[i]]$name == 'TOTAL') next
#     ODDy$polygon_id[ODDy@polygons[[i]]$indexes] = i
#   }
#   
#   impact_filt = ODDy@impact %>% filter(impact == impact_type)
#   impact_filt$median = apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) median(as.numeric(x)))
#   
#   poly_ids = as.data.frame(ODDy$polygon_id)
#   poly_ids$order = 1:nrow(poly_ids)
#   poly_ids_impact = merge(poly_ids,impact_filt[,c('polygon', 'median')], by.x='polygon_id', by.y='polygon', all.x=T)
#   poly_ids_impact = poly_ids_impact[order(poly_ids_impact$order),]
#   ODDy[[paste0(impact_type, '_median')]] = poly_ids_impact['median']
#   plot(ODDy[[paste0(impact_type, '_median')]])
#   
# }

plot_GADM_impact_polygons <- function(ODDy, polygons_list, gadm_level = 2, impact_type = 'mortality', summary = 'mean'){
  #GADM_level 0 = national
  #GADM_level 1 = admin level 1
  #GADM_level 2 = admin level 2
  #summary = mean, median, q5, q10, q20, q80, q90, q95
  
  ODDy$polygon_id = 0
  for (i in 1:length(ODDy@polygons)){
    if (str_count(ODDy@polygons[[i]]$name, ',') != gadm_level | ODDy@polygons[[i]]$name == 'TOTAL') next
    ODDy$polygon_id[ODDy@polygons[[i]]$indexes] = i
  }
  
  impact_filt = ODDy@impact[ODDy@impact$impact == impact_type,]
  impact_filt$median = apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) median(as.numeric(x)))
  impact_filt$mean = apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) mean(as.numeric(x)))
  impact_filt[paste0('q', c(05, 10, 20, 80, 90, 95))] = t(apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) quantile(as.numeric(x), c(.05, 0.10, 0.20, 0.80, 0.90, 0.95))))
  
  # poly_ids = as.data.frame(ODDy$polygon_id)
  # poly_ids$order = 1:nrow(poly_ids)
  # poly_ids_impact = merge(poly_ids,impact_filt[,c('polygon', 'median')], by.x='polygon_id', by.y='polygon', all.x=T)
  # poly_ids_impact = poly_ids_impact[order(poly_ids_impact$order),]
  # ODDy[[paste0(impact_type, '_median')]] = poly_ids_impact['median']
  # plot(ODDy[[paste0(impact_type, '_median')]])
  
  # Convert SpatialPolygonsDataFrame to a tidy format for ggplot2
  polygon_list_plot <- list()
  for (i in 1:length(polygons_list)){
    if (is.null(polygons_list[[i]]$sf_polygon$NAME_2)) next #polygons_list$polygons_list[[i]]$sf_polygon$NAME_2 = ''
    polygon_list_plot[[i]] = st_as_sf(polygons_list[[i]]$sf_polygon)
    polygon_list_plot[[i]]$impact = pull(impact_filt[summary])[which(impact_filt$polygon_name == polygons_list[[i]]$polygon_name)]
  }
  combined_polygons <- do.call(rbind, polygon_list_plot)
  
  ggplot(data = combined_polygons) +
    geom_sf(aes(fill = impact), color = "black", size = 0.2) +
    scale_fill_gradientn(colors = c("#FFF5E1", "#FDBB84", "#D73027"), trans = "sqrt")+  # Color scale
    theme_minimal() +
    labs(title = "Impact by Region", fill = paste(impact_type, summary)) + 
    xlim(ext(ODDy)[c(1,2)]) + ylim(ext(ODDy)[c(3,4)])
  
}
# library(shiny)
# library(ggplot2)
# library(dplyr)
# library(viridis)
# library(geodata)
# library(scales)
# 
# # Define UI
# ui <- fluidPage(
#   titlePanel("Interactive Impact Plot"),
#   
#   sidebarLayout(
#     sidebarPanel(
#       selectInput("impact_type", "Select Impact Type:",
#                   choices = c("mortality", "displacement"),
#                   selected = "mortality"),
#       
#       sliderInput("quantile", "Select Quantile:",
#                   min = 1, max = 100, value = 95, step = 1)
#     ),
#     
#     mainPanel(
#       plotOutput("impactPlot")
#     )
#   )
# )
# 
# # Define Server
# server <- function(input, output) {
#   output$impactPlot <- renderPlot({
#     # Reactive variables
#     impact_type <- input$impact_type
#     impact_type_samples <- paste0(impact_type, ".s")
#     quantile <- input$quantile
#     
#     # Data processing
#     
#     #ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908_MAR_-1AggLevel5_WithImpact100_Dec24')
#     
#     plot_df <- as.data.frame(ODDy, xy = TRUE, na.rm = FALSE)
#     plot_df[, grep(impact_type_samples, names(plot_df))] <- t(apply(
#       as.matrix(plot_df[, grep(impact_type_samples, names(plot_df)), drop = FALSE]), 
#       MARGIN = 1, sort
#     ))
#     
#     var <- paste0(impact_type, ".s", quantile)
#     plot_df[which(is.na(plot_df$ISO3C)), var] <- NA
#     names(plot_df)[which(names(plot_df) == "x")] <- "Longitude"
#     names(plot_df)[which(names(plot_df) == "y")] <- "Latitude"
#     
#     # Bounding Box
#     bbox <- matrix(ODDy@hazinfo$bbox, nrow = 2)
#     iso3_unique <- unique(ODDy$ISO3C)
#     iso3_unique <- iso3_unique[!is.na(iso3_unique)]
#     gadm_iso <- gadm_iso <- as(geodata::gadm(country = iso3_unique[1], level = 2, path = tempdir()), "Spatial")#as(geodata::gadm(country = iso3_unique[1], level = 2), "Spatial")
#     
#     my_path <- tempdir()  # Use a temporary directory or specify your own
#     
#     gadm_iso <- as(geodata::gadm(country = iso3_unique[1], level = 2, path = my_path), "Spatial")
#     if (length(iso3_unique) > 1) {
#       for (i in 2:length(iso3_unique)) {
#         gadm_iso <- rbind(gadm_iso, as(geodata::gadm(country = iso3_unique[i], level = 2, path = my_path), "Spatial"))
#       }
#     }
#     
#     gadm_iso <- intersect(gadm_iso, bbox)
#     gadm_map <- fortify(gadm_iso)
#     
#     # Define breakpoints
#     max_var <- max(plot_df[, var], na.rm = TRUE)
#     pseudo_log_sigma <- 100
#     breakings <- if (max_var > 100000) {
#       c(0, 1000, 10000, 100000, 1000000)
#     } else if (max_var > 10000) {
#       c(0, 1000, 10000)
#     } else if (max_var > 1000) {
#       c(0, 500, 1000)
#     } else if (max_var > 100) {
#       c(0, 50, 100)
#     } else if (max_var > 50) {
#       c(0, 10, 50)
#     } else if (max_var > 20) {
#       pseudo_log_sigma <- 1
#       c(0, 2, 5, 10, 20, floor(max_var / 10) * 10)
#     } else if (max_var > 10) {
#       c(0, 5, 10)
#     } else if (max_var > 5) {
#       c(0, 5, 10)
#     } else {
#       c(0, 1)
#     }
#     
#     # Generate the plot
#     ggplot() +
#       geom_map(map = gadm_map, data = gadm_map, aes(map_id = id, group = id)) +
#       geom_polygon(data = gadm_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
#       geom_raster(data = plot_df, aes(x = Longitude, y = Latitude, fill = !!sym(var)), alpha = 0.75) +
#       geom_contour(data = plot_df, aes(Longitude, Latitude, z = hazMean1, colour = ..level..), alpha = 0.7, lwd = 0.8) +
#       scale_fill_viridis(trans = scales::pseudo_log_trans(base = 10, sigma = pseudo_log_sigma), breaks = breakings, labels = function(x) comma(x)) +
#       scale_color_gradientn(colors = c("transparent", "#fc9272", "#ef3b2c")) +
#       labs(fill = paste0("Impact: ", impact_type), colour = "Hazard Intensity") +
#       coord_equal() +
#       scale_x_continuous(expand = c(0, 0)) +
#       scale_y_continuous(expand = c(0, 0)) +
#       theme_minimal() +
#       theme(
#         axis.title = element_text(family = "Liberation Serif", size = 12),
#         legend.text = element_text(family = "Liberation Serif", size = 11),
#         legend.title = element_text(family = "Liberation Serif", size = 12)
#       )
#   })
# }
# 
# # Run the Shiny App
# shinyApp(ui = ui, server = server)
# 
#   
#   