

setClass("ODDSim", contains="ODD")
setClass("BDSim", contains="BD")
# Get ODDSim objects to inherit the methods and attributes of ODD objects

# Make a new method for the initialisation of the ODDSim object. Much is the same as for an ODD object,
# but some steps involve simulation rather than retrieving the necessary data. 
setMethod(f="initialize", signature="ODDSim",
          definition=function(.Object,lhazSDF=NULL,DamageData=NULL, PopSim=NULL, GDPSim = NULL, EQFreqSim=NULL, 
                              Vs30Sim=NULL, dir="./",Model=list(
            WID_perc<-   c("p0p10", # Bottom 10% share of Income Distribution
                           "p10p20", # Income share held by 10th - 20th percentiles
                           "p20p30", # Income share held by 20th - 30th percentiles
                           "p30p40", # Income share held by 30th - 40th percentiles
                           "p40p50", # Income share held by 40th - 50th percentiles
                           "p50p60", # Income share held by 50th - 60th percentiles
                           "p60p70", # Income share held by 60th - 70th percentiles
                           "p70p80", # Income share held by 70th - 80th percentiles
                           "p80p90", # Income share held by 80th - 90th percentiles
                           "p90p100"), # top 10% share of Income Distribution
            impacts <- list(labels = c('mortality', 'displacement', 'buildDam', 'buildDest'), 
                            qualifiers = c('qualifierMort', 'qualifierDisp', 'qualifierBuildDam', 'qualifierBuildDest'),
                            sampled = c('mort_sampled', 'disp_sampled', 'buildDam_sampled', 'buildDest_sampled'))
            )) {
            
            if(is.null(lhazSDF)) return(.Object)
            if(!class(lhazSDF[[length(lhazSDF)]])[1]=="HAZARD") return(.Object)
            
            if(lhazSDF$hazard_info$hazard=="EQ"){ Model$INFORM_vars%<>%c("HA.NAT.EQ")
            } else if(lhazSDF$hazard_info$hazard=="TC"){Model$INFORM_vars%<>%c("HA.NAT.TC")
            } else if(lhazSDF$hazard_info$hazard=="FL"){Model$INFORM_vars%<>%c("HA.NAT.FL")
            } else stop("Not currently prepared for hazards other than EQ, TC or FL")
            
            
            .Object@dir<-dir
            .Object@hazard<-lhazSDF$hazard_info$hazard
            
            
            
            if(length(unique(DamageData$eventid))==1) .Object@eventid<-unique(DamageData$eventid)
        
            # This bounding box is taken as the minimum region that encompasses all hazard events in HAZARD object:
            bbox<-lhazSDF$hazard_info$bbox
            dater<-min(lhazSDF$hazard_info$sdate)
            .Object@hazdates<-as.Date(lhazSDF$hazard_info$eventdates)
            .Object@hazinfo <-lhazSDF$hazard_info
            
            year<-AsYear(dater)
            
            obj<-PopSim
            .Object@data <- obj@data
            .Object@coords.nrs <-obj@coords.nrs
            .Object@grid <-obj@grid
            .Object@grid.index <-obj@grid.index
            .Object@coords <-obj@coords
            .Object@bbox <-obj@bbox 
            .Object@proj4string <- PopSim@proj4string # crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            # Including minshake polygon per hazard event using getcontour from adehabitatMA package
            # LOOSEEND
            .Object%<>%AddHazSDF(lhazSDF)
            
            # Extract empty indices to save time
            inds<-!is.na(.Object$Population)
            
            .Object$GDP <- GDPSim
            .Object@data$ISO3C<-NA_character_
            .Object@data$ISO3C[inds]<-coords2country(obj@coords)[inds]
            iso3c<-unique(.Object@data$ISO3C) ; iso3c<-iso3c[!is.na(iso3c)]
            
            
            if (runif(1) < 1){ # with probability 0.025, allocate all pixels to polygons and have a lot of subnational data
              subnat_complete <- T
              grouping_mechanism <- runif(1)
              if (grouping_mechanism < 0.65 | max(.Object@data[,grep('hazMean', colnames(.Object@data))], na.rm=T) < 7){
                n_col_groups <- round(rpois(1,0.3))+1
                n_row_groups <- 1
              } else if (grouping_mechanism < 0.825){
                n_col_groups <- round(runif(1,4,20))
                n_row_groups <- 1
                if (n_col_groups==4){ n_col_groups=2; n_row_groups=2}
                if (n_col_groups==6){ n_col_groups=3; n_row_groups=2}
                if (n_col_groups==8){ n_col_groups=4; n_row_groups=2}
                if (n_col_groups==9){ n_col_groups=3; n_row_groups=3}
                if (n_col_groups==10){ n_col_groups=5; n_row_groups=2}
                if (n_col_groups==12){ n_col_groups=4; n_row_groups=3}
                if (n_col_groups==14){ n_col_groups=7; n_row_groups=2}
                if (n_col_groups==15){ n_col_groups=5; n_row_groups=3}
                if (n_col_groups==18){ n_col_groups=6; n_row_groups=3}
                if (n_col_groups==20){ n_col_groups=5; n_row_groups=4}
              } else {
                n_col_groups <- round(runif(1,3,15))#round(rweibull(1,0.3,0.15)+1) # round(runif(1,2, 15))
                n_row_groups <- round(runif(1,3,15)) #n_col_groups#round(rweibull(1,0.5,1.2)+1) # round(runif(1,3, 15))
                while ((n_col_groups*n_row_groups) > 170){
                  n_col_groups <- round(runif(1,3,15)) #round(rweibull(1,0.3,0.15)+1) # round(runif(1,2, 15))
                  n_row_groups <- round(runif(1,3,15)) #n_col_groups#round(rweibull(1,0.5,1.2)+1) #round(runif(1,3, 15))
                }
              }
              n_polygons <- n_col_groups * n_row_groups
              
            } else { # with remaining probability, allocate less pixels and have minimal subnational data
              subnat_complete <- F
              n_polygons <- ifelse(runif(1)>0.5, 1, round(rexp(1,0.3))+1)
            }
            
            polygons_list <- list()
            if (subnat_complete){ 
              ncol = sqrt(NROW(.Object))
              indexes_matrix <- t(matrix(1:NROW(.Object), ncol=ncol))
              indexes_matrix <- indexes_matrix[NROW(indexes_matrix):1,]
              nrow = ncol
              
              if (n_polygons > 5){
                #take the centre portion to get rid of some of the low intensity areas
                indexes_matrix <- indexes_matrix #[(nrow/4):(3*nrow/4),(ncol/4):(3*ncol/4)]
              }
              ncol <- NCOL(indexes_matrix); nrow=NROW(indexes_matrix)
              
              break.pts_col <- sort(sample(1:ncol, n_col_groups - 1, replace = FALSE))
              break.len_col <- diff(c(0, break.pts_col, ncol))
              groups_col <- rep(1:n_col_groups, times = break.len_col)
              
              break.pts_row <- sort(sample(1:nrow, n_row_groups - 1, replace = FALSE))
              break.len_row <- diff(c(0, break.pts_row, nrow))
              groups_row <- rep(1:n_row_groups, times = break.len_row)
              for (i in 1:n_polygons){
                col_group <- (i-1) %% n_col_groups + 1
                col_group <- ifelse(col_group==0, n_col_groups, col_group)
                row_group <- (i-1) %/% n_col_groups +1
                indexes <- indexes_matrix[groups_row==row_group,groups_col==col_group]
                polygons_list[[i]] <- list(name=paste('Polygon',i), 
                                           indexes = c(indexes),
                                           weights = rep(1, length(indexes)))
              }
              
            } else {
              polygons_list[[1]] <- list(name='Polygon 1', indexes = 1:length(.Object@data$ISO3C), weights=rep(1, length(.Object@data$ISO3C)))
              if (n_polygons > 1){
                pixels_allocated <- c()
                for (i in 2:n_polygons){ 
                  nPixels_in_poly <- c()
                  while(length(nPixels_in_poly)==0){
                    nPixels_in_poly <-  runif(1, 1, length(.Object@data$ISO3C))
                    Pixels_in_poly <- sample(1:length(.Object@data$ISO3C),nPixels_in_poly, replace=F)
                    Pixels_in_poly <- Pixels_in_poly[!(which(Pixels_in_poly) %in% pixels_allocated)] # avoid double counting
                  }
                  pixels_allocated <- c(pixels_allocated, Pixels_in_poly)
                  polygons_list[[i]] <- list(name=paste('Polygon',i), 
                                             indexes = Pixels_in_poly,
                                             weights = rep(1, length(Pixels_in_poly)))
                }
              }
            }
            
            
            .Object@data$nBuildings <- round(runif(1,0.2,1) * PopSim$Population + rnorm(length(PopSim$Population),0,20))
            .Object@data$nBuildings[.Object@data$nBuildings < 0] <- 0
            
            if(length(iso3c)==0){return(NULL)}
            
            if(.Object@hazard%in%c("EQ","TC")){
              .Object@impact <- data.frame(iso3=character(), polygon=numeric(), impact=character(),
                                           observed=numeric(), qualifier = character(), inferred=logical(), sdate=as.Date(character()))
              for (j in 1:length(Model$impacts$labels[1:3])){ #remove buildDest and buildDamDest
                n_poly_observed <- n_polygons
                if (subnat_complete & (Model$impacts$labels[j] == 'displacement')){ # less likely to have full coverage
                  if (runif(1) < 0.15){
                    n_poly_observed = 0
                    next
                  } else if (n_poly_observed < 10){
                    for (i in 1:n_poly_observed){
                      .Object@impact %<>% add_row(
                        iso3=sample(na.omit(.Object@data$ISO3C),1),
                        polygon=i,
                        impact=Model$impacts$labels[j],
                        observed= 0,
                        qualifier=ifelse(rbinom(1,1,0.8), 'Total', NA_character_), 
                        inferred=F,
                        sdate=dater)
                    }
                    next
                  } else if (n_poly_observed > 10){
                    # add a single total polygon
                    .Object@impact %<>% add_row(
                      iso3='TOT',
                      polygon=n_polygons + 1,
                      impact=Model$impacts$labels[j],
                      observed= 0,
                      qualifier=ifelse(rbinom(1,1,0.8), 'Total', NA_character_), 
                      inferred=F,
                      sdate=dater)
                    polygons_list[[n_polygons+1]] <- list(name='Total', 
                                               indexes = 1:length(.Object$ISO3C),
                                               weights = rep(1, length(.Object$ISO3C)))
                    next
                  }
                }
                if (Model$impacts$labels[j] == 'buildDam'){
                  if (runif(1)<0.2){ # no building damage data
                    next
                  }
                }
                for (i in 1:n_poly_observed){
                  
                  .Object@impact %<>% add_row(
                    iso3=sample(na.omit(.Object@data$ISO3C),1),
                    polygon=i,
                    impact=Model$impacts$labels[j],
                    observed= 0,
                    qualifier=ifelse(rbinom(1,1,0.8), 'Total', NA_character_), 
                    inferred=F,
                    sdate=dater)
                }
              }
              #remove about third as currently simulated data is very populated compared to real data
              # if (!subnat_complete){
              #   .Object@impact <- .Object@impact[sample(1:NROW(.Object@impact),round(NROW(.Object@impact)*2/3), replace=F), ]
              # }
              
              # .Object@IDPs<-DamageData[,c("sdate","gmax","qualifierDisp")]%>%
              #   transmute(date=sdate,IDPs=gmax,qualifier=qualifierDisp)
            } else {
              stop('Not setup for other hazards')
              # THIS IS READY FOR THE IPM APPROACH FOR MID/LONG DURATION HAZARDS
              # .Object@IDPs<-DamageData%>%group_by(sdate)%>%summarise(IDPs=max(IDPs),.groups = 'drop_last')%>%
              #   transmute(date=sdate,IDPs=IDPs)
              # Note that I leave .Object@gmax intentionally blank
            }
            
            # Skip the interpolation of population and GDP values for the simulations
            
            # # Simulate INFORM (Joint Research Center - JRC) data:
            # # Sample the indicators from a uniform distribution on [0,10]
            # # Set some to 0 each with a probability of 0.1 (to mimic presence of zero values in real data).
            # INFORM <- data.frame(iso3='ABC',
            #                      value=runif(length(Model$INFORM_vars),0,10),
            #                      variable=Model$INFORM_vars)
            # INFORM_missing <- rbernoulli(nrow(INFORM), 0.1)
            # INFORM$value[INFORM_missing]<-0
            
            #Simulate World Income Database (WID) data
            #Simulate a cubic distribution with max given by max_inc (WID data appears close to cubic)
            WID <- data.frame(variable=numeric(), iso3=character(), value=numeric())
            for (iso3 in unique(na.omit(.Object@data$ISO3C))){
              max_inc <- runif(1,0.4,0.8)
              perc <- seq(0.1,1,0.1)
              cubic <- perc^3+runif(1,0.3,1)*perc^2
              WID %<>% rbind(data.frame(
                percentile=Model$WID_perc,
                iso3=iso3,
                value=cubic*max_inc/max(cubic)
              ))
              mins<-WID%>%group_by(iso3)%>%summarise(mins=min(value),.groups = 'drop_last')
              for(iso3c in mins$iso3) WID$value[WID$iso3==iso3c & WID$variable!="p0p10"]%<>%subtract(WID$value[WID$iso3==iso3c & WID$variable!="p90p100"])
              for(iso3c in mins$iso3) WID$value[WID$iso3==iso3c] = (WID$value[WID$iso3==iso3c]/sum(WID$value[WID$iso3==iso3c]))
            }
            
            # Bind it all together!
            .Object@cIndies<-WID
            # .Object@fIndies<-Model$fIndies
            # 
            # linp<-rep(list(0.),length(unique(.Object@cIndies$iso3)))
            # names(linp)<-unique(.Object@cIndies$iso3)
            # .Object@modifier<-linp
            # 
            
            
            
            .Object@polygons <- polygons_list
            
            #Add linear predictor data
            #need to tidy this up, not very reflective of real data to have this much variation within such a small region!
            ulist <- unique(GDPSim)
            AveSchYrs_vals <- runif(length(ulist), 3, 18)
            LifeExp_vals <- runif(length(ulist), 30, 85)
            GNIc_vals <- exp(runif(length(ulist), 6, 12))
            # Vs30_vals <- runif(length(ulist), 98, 2197)
            # EQFreq_vals <- runif(length(ulist), 1, 10)
            SHDI_vals <- runif(length(ulist), 0.1, 0.9)
            for(i in 1:length(ulist)){
              .Object@data$AveSchYrs[GDPSim %in% ulist[i]] <- AveSchYrs_vals[1] + rbeta(1, 1,2) * (AveSchYrs_vals[i] - AveSchYrs_vals[1]) #cluster the values towards the first
              .Object@data$LifeExp[GDPSim %in% ulist[i]] <- LifeExp_vals[1] + rbeta(1, 1,2) * (LifeExp_vals[i] - LifeExp_vals[1])
              .Object@data$GNIc[GDPSim %in% ulist[i]] <- GNIc_vals[1] + rbeta(1, 1,2) * (GNIc_vals[i] - GNIc_vals[1])
              # .Object@data$Vs30[GDPSim %in% ulist[i]] <- Vs30_vals[1] + rbeta(1, 1,5) * (Vs30_vals[i] - Vs30_vals[1])
              # .Object@data$EQFreq[GDPSim %in% ulist[i]] <- EQFreq_vals[1] + rbeta(1, 1,5) * (EQFreq_vals[i] - EQFreq_vals[1])
              .Object@data$SHDI[GDPSim %in% ulist[i]] <- SHDI_vals[1] + rbeta(1, 1,2) * (SHDI_vals[i] - SHDI_vals[1])
            }
            .Object@data$PDens <- .Object@data$Population
            .Object@data$EQFreq <- EQFreqSim$EQFreq
            .Object@data$Vs30 <- Vs30Sim$Vs30
            
            print("Checking ODDSim values")
            checkODD(.Object)
            
            
            return(.Object)
          }
)

setMethod(f="initialize", signature="BDSim",
          definition=function(.Object,ODD=NULL) {
            
            if(is.null(ODD)) return(.Object)            
            
            .Object@hazard<-ODD@hazard
            #.Object@cIndies<-ODD@cIndies
            .Object@I0<-ODD@I0
            .Object@hazdates<-ODD@hazdates
            .Object@eventid<-ODD@eventid
            #.Object@fIndies<-ODD@fIndies
            
            .Object@buildingsfile<-paste0("./IIDIPUS_Input/OSM_Buildings_Objects/",unique(ODD@eventid)[1])
            
            Damage = data.frame(Latitude=double(),Longitude=double(), grading=character(), Confidence=character())
            cells_with_sat <- sample(1:NROW(ODD@coords), runif(1,3,20))
            
            for (ij in cells_with_sat){
              lonmin <- ODD@coords[ij,'Longitude']-ODD@grid@cellsize['Longitude']/2
              latmin <- ODD@coords[ij,'Latitude']-ODD@grid@cellsize['Latitude']/2
              lonmax <- ODD@coords[ij,'Longitude']+ODD@grid@cellsize['Longitude']/2
              latmax <- ODD@coords[ij,'Latitude']+ODD@grid@cellsize['Latitude']/2
              if(ODD$nBuildings[ij] == 0){
                next
              }
              for (k in 1:min(20, ODD$nBuildings[ij])){ #avoid having more than 20 buildings per grid in simulated data
                                                        #to keep the computational requirements low
                Damage %<>% add_row(
                  Longitude = runif(1,lonmin, lonmax),
                  Latitude = runif(1, latmin, latmax),
                  grading = NA, 
                  Confidence = NA
                )
              }
            }
            
            print("Forming SpatialPointsDataFrame from building damage data")
            Damage<-SpatialPointsDataFrame(coords = Damage[,c("Longitude","Latitude")],
                                           data = Damage[,c("grading","Confidence")],
                                           proj4string = ODD@proj4string)
            
            .Object@data <- Damage@data
            .Object@coords.nrs <-Damage@coords.nrs
            .Object@coords <-Damage@coords
            .Object@bbox <-Damage@bbox
            .Object@proj4string <-Damage@proj4string #crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            rm(Damage)
            
            print("Interpolating population density, hazard & GDP-PPP data")
            ODD$spatial_pixel <- 1:NROW(ODD)
            .Object%<>%BDinterpODD(ODD=ODD)
            
            print("Filter spatial data per country")
            # We could just copy Damage$iso3 directly, but I don't believe in anyone or anything...
            if(NROW(.Object@data)>0) .Object@data$ISO3C<-coords2country(.Object@coords)
            
            print("Accessing OSM to sample building height & area")
            # ExtractOSMbuildVol(.Object,ODD)
            
            # linp<-rep(list(0.),length(unique(ODD@cIndies$iso3)))
            # names(linp)<-unique(ODD@cIndies$iso3)
            # .Object@modifier<-linp
            
            print("Checking BD values")
            # checkBD(.Object)
            
            return(.Object)
          }
)

simulateEvent <- function(r, I0 = 4.5){ 
  # Input: 
  # - An empty raster field r (over the grid of interest)
  # Output:
  # - A spatial pixel data frame hazsdf containing the hazard intensity in each grid cell
  # Details:
  # - The earthquake intensity follows a gaussian kernel centered at (0,0)
  # - The maximum magnitude varies randomly in [5, 9.3] and the spread in [15, 25]
  # - The earthquake standard deviation in each cell is random uniform in [0.8,1.1]
  
  min_mmi =4.45
  max_mmi = max(4.7, rnorm(1,6.5, 1.1))
   # rnorm(1, 1, 0.28)
  sigma = 0.2 * (max_mmi-4)^3 + runif(1,-0.5,0.5) + 1#0.7*(max_mmi - 4.5)^2.2 + (max_mmi - 6) + 2 + runif(1,-0.1,0.1) #runif(1, 3,8.5)
  r <- setValues(r, spatialEco::gaussian.kernel(sigma=sigma, s=r@nrows))
  #r <- exp(r-min_mmi)
  r <- r * ((max_mmi-min_mmi)/r@data@max) + min_mmi
 
  cells_above_I0 <- xyFromCell(r, which(values(r>4.5)))
  crop_extent <- c(min(cells_above_I0[,1]), max(cells_above_I0[,1]), min(cells_above_I0[,2]), max(cells_above_I0[,2]))
  r <- crop(r, crop_extent, snap='out')
  #values(r)[values(r) < I0] = NA
  #r <- trim(r)
  #r <- (r - r@data@min) / (r@data@max-r@data@min) * (max_mmi-min_mmi) + min_mmi
  sd <- setValues(r, runif(r@ncols*r@nrows, 0.8, 1.1))
  #r <- 4.51 + exp(r)
  names(r) <- 'mmi_mean'
  r$mmi_sd <- sd
  hazsdf <- as(r, 'SpatialPixelsDataFrame')
  #hazsdf<-hazsdf[hazsdf$mmi_mean>I0,]
  colnames(hazsdf@coords)<-rownames(hazsdf@bbox)<-c("Longitude","Latitude")

  return(hazsdf)
}

simulatePopulation <- function(r){
  # Input: 
  # - An empty raster field r (over the grid of interest)
  # Output: 
  # - PopSim: a spatial pixels data frame containing the population per grid cell according to a simulation
  # Details:
  # - The population is simulated according to a Gaussian process with a spherical semi-variogram,
  #   and is then scaled onto [0, maxPopDens] using the cumulative weibull distribution.
  
  maxPopDens = runif(1, 5000, 50000)
  Field = as.data.frame(rasterToPoints(r))
  names(Field)=c('Longitude','Latitude')
  Pop_modelling=gstat(formula=Population~1, 
                      locations=~Latitude+Longitude,
                      dummy=T,    
                      beta=30,  
                      model=vgm(psill=0.3,
                                range=0.5,
                                nugget=0.01,
                                model='Sph'), 
                      nmax=40)
  
  PopSim = predict(Pop_modelling, newdata=Field, nsim=1)
  minval = min(PopSim[,3])
  maxval = max(PopSim[,3])
  PopSim[,3] = pweibull( (PopSim[,3]-minval)/(maxval-minval),5,2) /pweibull(1,5,2)* maxPopDens
  PopSim %<>% acast(Longitude~Latitude, value.var='sim1') %>% convMat2SPDF(name='Population')
  return(PopSim)
}

simulateVs30 <- function(r){
  # - works same as for population, but is a bit smoother

  Vs30_range = runif(2, 180, 880)
  while(abs(diff(Vs30_range))<30){
    Vs30_range = runif(2, 180, 880)
  }
  Field = as.data.frame(rasterToPoints(r))
  names(Field)=c('Longitude','Latitude')
  Vs30_modelling=gstat(formula=Population~1, 
                      locations=~Latitude+Longitude,
                      dummy=T,    
                      beta=30,  
                      model=vgm(psill=0.3,
                                range=2,
                                nugget=0.005,
                                model='Sph'), 
                      nmax=40)
  Vs30Sim = predict(Vs30_modelling, newdata=Field, nsim=1)
  minval = min(Vs30Sim[,3])
  maxval = max(Vs30Sim[,3])
  Vs30Sim[,3] = (Vs30Sim[,3]-minval)/(maxval-minval) * (abs(diff(Vs30_range))) + min(Vs30_range)
  Vs30Sim %<>% acast(Longitude~Latitude, value.var='sim1') %>% convMat2SPDF(name='Vs30')
  return(Vs30Sim)
}

simulateEQFreq <- function(r){
  # - works same as for Vs30, but even smoother
  EQFreq_min = rnorm(1, 200, 50)
  EQFreq_max = EQFreq_min + rbeta(1, 0.5, 5) * 350
 
  Field = as.data.frame(rasterToPoints(r))
  names(Field)=c('Longitude','Latitude')
  EQFreq_modelling=gstat(formula=Population~1, 
                       locations=~Latitude+Longitude,
                       dummy=T,    
                       beta=30,  
                       model=vgm(psill=0.00001, model= "Sph", range= 5,  kappa = 5), 
                       nmax=40)
  EQFreqSim = predict(EQFreq_modelling, newdata=Field, nsim=1)
  minval = min(EQFreqSim[,3])
  maxval = max(EQFreqSim[,3])
  EQFreqSim[,3] = (EQFreqSim[,3]-minval)/(maxval-minval) * (EQFreq_max-EQFreq_min) + EQFreq_min
  EQFreqSim %<>% acast(Longitude~Latitude, value.var='sim1') %>% convMat2SPDF(name='EQFreq')
  return(EQFreqSim)
}

simulateGDP <- function(r){
  # Input: 
  # - An empty raster field r (over the grid of interest)
  # Output: 
  # - A vector GDP containing the simulated GDP for each cell in r. 
  # Details:
  # - The GDP is simulated by generating a number of points and assigning each a GDP. Each cell in r is then assigned 
  #   the GDP of the closest point, therefore partitioning the data. This is an attempt to mimic the administrative 
  #   boundaries present in the actual GDP data.
  
  nRegions = rpois(1, 5) + 1
  regionPoints = cbind(runif(nRegions,r@extent@xmin, r@extent@xmax),runif(nRegions,r@extent@ymin, r@extent@ymax))
  GDPperRegion = round(sort(runif(nRegions, 100, 100000)))
  r_mat <- rasterToPoints(r)
  GDP = array(0, NROW(r_mat))
  for (cell in 1:NROW(r_mat)){
    GDP[cell] = GDPperRegion[which.min(distmat(r_mat[cell,], regionPoints))]
  }
  
  # ggplot() + geom_point(aes(x=r_mat[,1],y=r_mat[,2],col=as.factor(GDP)))
  
  return(GDP)
  
}

simulateODDSim <- function(miniDam, Model, I0=4.5){ 
  # Simulates hazard, population and GDP data over 120 x 120 grid cells (over the region [-0.5,0.5] x [-0.5,0.5])
  # INPUT: 
  #  - I0: minimum intensity resulting in damage 
  #  - miniDam: contains hazard info such as sdate, fdate, eventid, etc.
  # OUTPUT: 
  #  - ODDSim: an object of class ODDSim (the simulated equivalent of an ODD object)
  
  xmn <- round(runif(1, -180, 179.5)/0.05)*0.05
  ymn <- round(runif(1, -60, 79.5)/0.05)*0.05
  xsize <- 25/12
  ysize <- 25/12
  r <- raster(ncol=50, nrow=50, xmn=xmn, xmx=xmn+xsize, ymn=ymn,ymx=ymn+ysize, crs="+proj=longlat +datum=WGS84") #each cell is 30 arcseconds x 30 arcseconds
  
  while(mean(is.na(coords2country(as(r, 'SpatialPixelsDataFrame')@coords)))>0.5){ #repeat until over a country
    xmn <- round(runif(1, -180, 179.5)/0.05)*0.05
    ymn <- round(runif(1, -60, 79.5)/0.05)*0.05
    r <- raster(ncol=50, nrow=50, xmn=xmn, xmx=xmn+xsize, ymn=ymn,ymx=ymn+ysize, crs="+proj=longlat +datum=WGS84") #each cell is 30 arcseconds x 30 arcseconds
  }
  
  lenny = ifelse(runif(1)<0.5, rgeom(1, 0.6) + 1, 1) #generate number of events according to a geometric distribution
  bbox = c(xmn,ymn,xmn+xsize,ymn+ysize) 
  lhazdat<-list(hazard_info=list(bbox=bbox,sdate=min(miniDam$sdate),fdate=max(miniDam$fdate),
                                 NumEvents=lenny,hazard="EQ", I0=I0, eventdates=c(), depths=c(), 
                                magnitudes=c(), max_mmi=c(), eventtimes=c(), usgs_ids=c(), alertfull=list()))
  
  for(i in 1:lenny){
    hazsdf <- simulateEvent(r, I0)
    while(NROW(hazsdf)<=9 | length(which(hazsdf@data$mmi_mean > I0))<=1 | is.null(hazsdf@data$mmi_mean)){
      hazsdf <- simulateEvent(r, I0)
    }
  
    
    haz <- new("HAZARD",
                  obj=hazsdf,
                  hazard="EQ",
                  dater= min(miniDam$sdate),
                  I0=I0,
                  alertlevel=ifelse(max(hazsdf$mmi_mean, na.rm=T)>7.5, 'red', ifelse(max(hazsdf$mmi_mean, na.rm=T)>6, 'orange', 'green')), #assign pretty arbitrarily, don't think this is used in the model, 
                  alertscore=0,
                  depth= runif(1,0,50), #doesn't matter as not currently included in model
                  magnitude=runif(1,4,10),  
                  max_mmi=max(hazsdf@data$mmi_mean, na.rm=T),
                  USGS_id=paste0('USGS', '_',miniDam$eventid[1], '_', round(runif(1,1000,9000))), #doesn't matter as not included in model
                  eventtime=format(as.POSIXct(sample(as.numeric(as.POSIXct("00:00:00", format = "%H:%M:%S")):as.numeric(as.POSIXct("23:59:59", format = "%H:%M:%S")), 1), origin = "1970-01-01", tz = "UTC"), "%H:%M:%S"),
                  alertfull=list()
    )
    lhazdat <- append(lhazdat, haz)
    
    lhazdat$hazard_info$alertlevel%<>%c(haz@alertlevel)
    lhazdat$hazard_info$eventdates%<>%c(as.character(haz@eventdate))
    lhazdat$hazard_info$depths%<>%c(haz@depth)
    lhazdat$hazard_info$magnitudes%<>%c(haz@magnitude)
    lhazdat$hazard_info$max_mmi%<>%c(haz@max_mmi)
    lhazdat$hazard_info$eventtimes%<>%c(haz@eventtime)
    lhazdat$hazard_info$usgs_ids%<>%c(haz@USGS_id)
    lhazdat$hazard_info$alertfull[[length(lhazdat$hazard_info$alertfull)+1]] <- haz@alertfull
  }
  #make the first event in the set the 'first event' with probability 0.8
  #this somewhat mimics the real dataset in which there may be preceding events not included in the ODD object
  first_event <- order(paste(lhazdat$hazard_info$eventdates, lhazdat$hazard_info$eventtimes))[1]
  lhazdat$hazard_info$first_event <- rep(F, length(lhazdat$hazard_info$eventdates))
  lhazdat$hazard_info$first_event[first_event] = as.logical(rbinom(1, 1, 0.8)) 

  
  if (length(lhazdat)== 1){
    print('Simulation failed: affected region under simulated event is too small')
    return(NULL)
  }
  
  max_extent <- extent(lhazdat[[2]])
  if (lenny > 1){
    for (i in 2:lenny){
      max_extent[c(1,3)] <- pmin(max_extent[c(1,3)], extent(lhazdat[[i]])[c(1,3)])
      max_extent[c(2,4)] <- pmax(max_extent[c(2,4)], extent(lhazdat[[i]])[c(2,4)])
    }
  }
 
  PopSim <- simulatePopulation(crop(r, max_extent))
  GDPSim <- simulateGDP(crop(r, max_extent))
  EQFreqSim <- simulateEQFreq(crop(r, max_extent))
  Vs30Sim <- simulateVs30(crop(r, max_extent))
  
  ODDSim <- tryCatch(new('ODDSim', lhazSDF=lhazdat,DamageData=miniDam, PopSim=PopSim, GDPSim=GDPSim, EQFreqSim=EQFreqSim, # no longer use GDP but use the spatial groupings of GDPSim for grouping vulnerability covariates
                         Vs30Sim=Vs30Sim, Model=Model), error=function(e) NULL)
  return(ODDSim)
}

simulateDataSet <- function(nEvents, Omega, Model, dir, outliers = FALSE, I0=4.5, folder_write='IIDIPUS_SimInput/'){
  # Input:
  # - nEvents: The number of ODDSim objects to generate
  # - Omega: The model parameterisation
  # - outliers: If true, includes an outlier at a high intensity with higher impact than simulated by the model (to check robustness of algorithm)
  # Output: 
  # - center: the center values of PDens and GDP under the simulated data
  # Details:
  # - Saves each ODDSim object to IIDIPUS_SimInput/ODDobjects/
  
  haz='EQ'
  
  # Rather than simulating eventids and sdates, just take real events but remove all data except 
  # sdate, fdate, eventid, and displacement/mortality qualifiers
  #DamageData = GetDisplacements(haz, saved=F, GIDD=F, EMDAT=T)
  #DamageData %<>% head(n=nEvents) %>% mutate(iso3='ABC', gmax=NA, mortality=NA) 
  
  DamageData <- data.frame(iso3='ABC', hazard='EQ', eventid=1:nEvents, 
                           sdate = sort(sample(seq(as.Date('2012/01/01'), as.Date('2024/01/01'), by="day"), nEvents)))
  DamageData$fdate <- DamageData$sdate + runif(nEvents, 0,14)
  
  set.seed(1)
  seeds_events <- sample(1:1000, nEvents, replace=F)
  for(i in 1:nEvents){
    set.seed(seeds_events[i]) # helps with debugging if issue with a particular event
    ev = unique(DamageData$eventid)[i]
    miniDam<-DamageData%>%filter(eventid==ev)
    ODDSim <- simulateODDSim(miniDam, Model, I0=I0)
    if (is.null(ODDSim) | is.null(ODDSim$hazMean1)){
      next
    }
    #BDSim <- new('BDSim', ODDSim)
    namer<-paste0(ODDSim@hazard,
                  str_remove_all(as.character.Date(min(ODDSim@hazdates)),"-"),
                  unique(miniDam$iso3)[1],
                  "_",ODDSim@eventid)
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,folder_write, "ODDobjects/",namer)
    saveRDS(ODDSim,ODDpath)
    
    #BDpath<-paste0(dir,folder_write,"BDobjects/",namer) 
    #saveRDS(BDSim,BDpath)
    print(paste0('Saved simulated hazard to ', namer))
  }
  
  #calculate Model$center based on the simulated events
  Model$center <- ExtractCentering(dir, input_folder=folder_write, saver=F)
  
  ODDpaths <-na.omit(list.files(path=paste0(folder_write,"ODDobjects/")))
  BDpaths <-na.omit(list.files(path=paste0(folder_write, "BDobjects/")))
  k <- 10
  intensities <- c() #store eq intensities
  #now loop through each event and simulate the displacement, mortality, and building destruction using DispX()
  for(i in 1:length(ODDpaths)){
    ODDSim <- readRDS(paste0(folder_write,"ODDobjects/",ODDpaths[i]))
    intensities <- append(intensities, max(ODDSim@data$hazMean1, na.rm=TRUE))
    #simulate displacement, mortality and building destruction using DispX
    ODDSim %<>% DispX(Omega %>% addTransfParams(), Model$center, output='ODDwithSampled',
                      Method=list(Np=1,cores=1, NestedCores=1))
    #take these simulations as the actual values
    ODDSim@impact <- data.frame(impact=character(), iso3=character(), qualifier=character(), value=numeric())
    ODDSim@impact <- ODDSim@predictDisp %>% mutate(observed=sampled, qualifier = ifelse(is.na(sampled), NA, 'total'))
    ODDSim@impact %<>% dplyr::select(-sampled)
    
    #overwrite ODDSim with the updated
    saveRDS(ODDSim, paste0(folder_write, "ODDobjects/",ODDpaths[i]))
    # if (i < 100){
    #   saveRDS(ODDSim, paste0("IIDIPUS_Input/ODDobjects/Train/",ODDpaths[i]))
    # } else {
    #   saveRDS(ODDSim, paste0("IIDIPUS_Input/ODDobjects/Test/",ODDpaths[i]))
    # }
    
  }
  # for (i in 1:length(BDpaths)){
  #   BDSim <- readRDS(paste0(folder_write,"BDobjects/", BDpaths[i]))
  #   BDSim %<>% BDX(Omega %>% addTransfParams(), Model, Method=list(Np=1,cores=1), output='sim')
  #   #take these simulations as the actual values
  #   BDSim@data$grading <- BDSim@data$ClassPred
  #   #switch those affected to 'Damaged' with probability 0.5
  #   affected <- which(BDSim@data$grading != 'notaffected')
  #   for (j in affected){
  #     BDSim@data$grading[j] = ifelse(runif(1)>0.5, BDSim@data$grading[j], 'Damaged')
  #   }
  #   if(NROW(BDSim@data)>100){
  #     BDSim@data <- BDSim@data[1:100,]
  #   }
  #   saveRDS(BDSim,paste0(folder_write, "BDobjects/",BDpaths[i]))
  # }
  
  if (outliers){
    #double the impact of the highest intensity earthquake
    i <- which.max(intensities)
    ODDSim <- readRDS(paste0(folder_write, "ODDobjects/",ODDpaths[i]))
    ODDSim@gmax %<>% mutate(
      gmax = 2 * gmax,
      mortality = 2 * mortality,
      buildDestroyed = 2 * buildDestroyed
    )
    saveRDS(ODDSim, paste0(folder_write, "ODDobjects/",ODDpaths[i]))
  }
  return(Model$center)
}

perturb_impacts <- function(d=1900, AlgoParams){
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  impacts_list <- list()
  for (i in 1:length(ufiles)){
    ufile <- ufiles[i]
    impacts_list[[i]] <- readRDS(paste0(folderin,ufile))@impact
  }
  LL_sum = d + 1
  d_recorded <- c()
  while (LL_sum > d){
    LL_sum <- 0
    for(j in 1:length(impacts_list)){
      print(j)
      impacts_list[[j]]$sampled <- NA
      for (k in 1:NROW(impacts_list[[j]])){
        cap <- rexp(1, 0.2)
        added_LL <- cap + 1
        while (added_LL > cap){
          impacts_list[[j]]$sampled[k] <- round(runif(1, 0, impacts_list[[j]]$observed[k] * 100+100))
          added_LL <- CalcPolyDist(impacts_list[[j]][k,], impact_weights=AlgoParams$impact_weights, kernel=AlgoParams$kernel, cap=AlgoParams$cap)
        }
        LL_sum <- LL_sum + added_LL
      }
    }
    d_recorded <- append(d_recorded, LL_sum)
  }
  
  for (i in 1:length(ufiles)){
    ufile <- ufiles[i]
    ODDy <- readRDS(paste0(folderin,ufile))
    ODDy@impact$observed <- impacts_list[[i]]$sampled
    saveRDS(ODDy, paste0("IIDIPUS_InputPerturbed/ODDobjects/",ufile))
  }
  
}


simulateRealData <- function(input_folder, Omega, Model, dir, outliers = FALSE, I0=4.5){
  #simulate impact data for real events
  input_folder <- "IIDIPUS_Input_RealAgg5/ODDobjects/"
  ODDpaths <-na.omit(list.files(path=input_folder, recursive=T))
  k <- 10
  #now loop through each event and simulate the displacement, mortality, and building damage using DispX()
  for(i in 1:length(ODDpaths)){
    ODD <- readRDS(paste0(input_folder,ODDpaths[i]))

    #simulate displacement, mortality and building destruction using DispX
    impact_sampled = DispX(ODD, Omega %>% addTransfParams(), Model$center, output='SampledAgg',
                      Method=list(Np=1,cores=1))[[1]]
    impact_sampled$observed = impact_sampled$sampled
    impact_sampled$sampled = NULL
    ODD@impact = impact_sampled

    saveRDS(ODD, paste0("IIDIPUS_SimInput_RealAgg5/ODDobjects/",ODDpaths[i]))
    
  }
}


#plot the S-curve for a given parameterisation (and optionally compare to a second)
plot_S_curves <- function(Omega, Omega_curr=NULL){
  Intensity <- seq(4,12,0.01)
  I0 = 4.5
  D_MortDisp <- D_MortDisp_calc(h_0(Intensity, I0, Omega), Omega %>% addTransfParams())
  D_Mort <- D_MortDisp[1,]
  D_Disp <- D_MortDisp[2,]
  D_Dam <- D_Dam_calc(h_0(Intensity, I0, Omega), Omega %>% addTransfParams())
  
  plot(Intensity, D_Mort, col='red', type='l', ylim=c(0,1)); lines(Intensity, D_Disp, col='orange');
  lines(Intensity, D_Dam, col='blue', type='l', ylim=c(0,1)); #lines(Intensity, D_Dam, col='dark green', type='l');
  
  if(!is.null(Omega_curr)){
    
    D_MortDisp_curr <- D_MortDisp_calc(h_0(Intensity, I0, Omega_curr), Omega_curr %>% addTransfParams())
    D_Mort_curr <- D_MortDisp_curr[1,]
    D_Disp_curr <- D_MortDisp_curr[2,]
    D_DestDam_curr <- D_DestDam_calc(h_0(Intensity, I0, Omega_curr), Omega_curr %>% addTransfParams())
    D_Dest_curr <- D_DestDam_curr[1,]
    D_Dam_curr <- D_DestDam_curr[2,]
    
    lines(Intensity, D_Mort_curr, col='pink', lty=2); lines(Intensity, D_Disp_curr, col='yellow', lty=2); 
    lines(Intensity, D_Dest_curr, col='cyan', lty=2, lwd=2); lines(Intensity, D_Dam_curr, col='green', lty=2, lwd=2);
  }
}


# # -------------------------------------------------------------------------
# 
# #Plot the simulated data
# #Names: ODDSim.png, Sim DispMortBD.png  
# #Size: 1500 x 700
# 


plotODDy_sim <- function(ODDy, zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.7,
                         log_plot=F, haz_legend=F, var_legend=T, discrete=F){
  #Plots background as GADM regions rather than terrain:
  
  bbox <- ODDy@bbox
  
  gg <- ggplot()
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
    
    p <- p+geom_raster(data = as.data.frame(ODDy),
                       mapping = aes(Longitude,Latitude,fill=ODDy@data[[var]]),alpha=0.8) + 
      labs(fill = ifelse(GetVarName(var)=='NULL', var, GetVarName(var))) + 
      scale_x_continuous(expand = c(0, 0)) +#f7fcb9
      scale_y_continuous(expand = c(0, 0)) 
    #scale_fill_distiller(palette = "Spectral", direction=-1)
    # +
    #scale_fill_gradientn(colors = c("#e0ecf4", "#9ebcda", "#8856a7", "#810f7c", "#4d004b"))
    #e0ecf4
    
    
    if (log_plot){
      p <- p + scale_fill_viridis( trans = "log10", labels = function(x) sprintf("%.0f", x))
    } else if (discrete){
      p <- p + scale_fill_viridis_d()
    } else if (var %in% c('SHDI', 'EQFreq')){
      #p <- p + scale_fill_viridis(option='mako')
      if (var_legend){
        p <- p + scale_fill_gradientn(colors=viridis_pal()(60)[1:6*10], 
                                      values=c(0, 0.1, 0.4, 0.6, 0.9, 1))
      } else {
        p <- p + scale_fill_gradientn(colors=viridis_pal()(60)[1:6*10], 
                                      values=c(0, 0.1, 0.4, 0.6, 0.9, 1), guide='none')
      }
      
    } else {
      if (var_legend){
        p <- p + scale_fill_gradientn(colors=viridis_pal()(10),
                                      values=c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.85,  1))
      } else {
        p <- p + scale_fill_gradientn(colors=viridis_pal()(10),
                                      values=c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.85,  1), guide='none')#(colors = c("#f7fcf5","#c7e9c0", "#addd8e", "#78c679", "#238443", "#006837","#004529")
      }
    } 
    
    p<-p+geom_contour(data = as.data.frame(ODDy),
                      mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
                      alpha=0.7,breaks = brks, lwd=0.8)
    if (haz_legend){
      p <- p + scale_color_gradientn(colors = c("transparent","#fc9272", "#ef3b2c")) + labs(colour = "Hazard Intensity      ")
    } else {
      p <- p + scale_color_gradientn(colors = c("transparent","#fc9272", "#ef3b2c"), guide='none')
    }
    
    #scale_color_gradientn(colors = c("","#fc9272", "#cb181d","#a50f15"))
    #scale_color_distiller(palette='RdBu')
    #scale_colour_gradient(low = "transparent",high = "#ef3b2c",na.value = "transparent") + 
    
    
    # p+geom_contour_filled(data = as.data.frame(ODDy),
    #                       mapping = aes(Longitude,Latitude,z=1-ODDy@data$tmp),
    #                       fill="green",alpha=alpha)+ 
    #   labs(fill = "Hazard>5")
    
    return(p)
  }
  
  ODDy@data$hazard[ODDy@data$hazard==0]<-NA
  
  p<-p+geom_contour_filled(data = as.data.frame(ODDy),
                           mapping = aes(Longitude,Latitude,z=hazard),
                           alpha=alpha,breaks = brks) +
    # scale_fill_discrete(low = "transparent",high = "red",na.value = "transparent") + 
    labs(fill = "Hazard Intensity")
  p<-p+geom_contour(data = as.data.frame(ODDy),
                    mapping = aes(Longitude,Latitude,z=hazard),
                    alpha=0.8,breaks = c(6.0),colour="red")
  
}

plot_sim_event <- function(ODDy){
  #ODDSim <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_SimInput2/ODDobjects/Train/EQ20120115ABC_1')
  # ODDSim <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_SimInput/ODDobjects/EQ20130409ABC_-3')
  # ODDSim <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_SimInput/ODDobjects/EQ20130416ABC_-4')
  # 
  plot_with_haz_legend <- plotODDy_sim(ODDy, var='Population', haz_legend=T, var_legend=F) + coord_fixed()+
    theme(legend.position="bottom", legend.key.size = unit(15, 'pt'),
          legend.key.width= unit(30, 'pt'),
          legend.text = element_text(margin = margin(t=2, r = 10)))
  haz_legend <- get_plot_component(plot_with_haz_legend, 'guide-box', return_all=T)[[3]]
  
  plots_all <- plot_grid( plotODDy_sim(ODDy, var='Population') + coord_fixed(),
                          plotODDy_sim(ODDy, var='nBuildings') + coord_fixed(), 
                          plotODDy_sim(ODDy, var='SHDI') + coord_fixed(), 
                          plotODDy_sim(ODDy, var='Vs30') + coord_fixed(), 
                          plotODDy_sim(ODDy, var='BuildDam') + coord_fixed(),
                          plotODDy_sim(ODDy, var='BuildDamAgg', discrete=T) + coord_fixed(), 
                          nrow=2, rel_heights=c(1,1), align='v', labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'),
                          label_fontface = "plain")
  
  plots_with_legend <- plot_grid(plots_all, haz_legend, nrow=2, rel_heights=c(1,0.1))
  plots_with_legend
}

# ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_SimInput3/ODDobjects/Train/EQ20200502ABC_118')
# ODDy$BuildDamAgg = NA
# ODDy$BuildDamAgg[ODDy@polygons[[1]]$indexes] <- ODDy@impact$observed[1]
# ODDy$BuildDamAgg[ODDy@polygons[[2]]$indexes] <- ODDy@impact$observed[4]
# ODDy$BuildDamAgg <- as.factor(ODDy$BuildDamAgg)
# plot_sim_event(ODDy)
#SimEvent.pdf , 6 x 12 inches


# grid.arrange(plotODDy(ODDSim, var='Population') + coord_fixed(),
#               plotODDy(ODDSim, var='nBuildings') + coord_fixed(),
#               plotODDy(ODDSim, var='GNIc') + coord_fixed(), nrow=1) # 10 x 4, Input
# 
# grid.arrange(plotODDy(ODDSim, var='Mort') + coord_fixed(),
#              plotODDy(ODDSim, var='BuildDam') + coord_fixed(), nrow=1)
# 
# plotODDy(ODDSim, var='Population') + coord_fixed() # 3.5 x 3.5, Population
# plotODDy(ODDSim, var='nBuildings') + coord_fixed() # 3.5 x 3.5, BuildCount
# plotODDy(ODDSim, var='GNIc') + coord_fixed()  # 3.5 x 3.5, GNIc
# 
# plotODDy(ODDSim, var='Disp') + coord_fixed() #5.5 x 5.5, Displacement
# plotODDy(ODDSim, var='BuildDest') + coord_fixed() #5.5 x 5.5, BuildDest
# 
# # \begin{figure}
# # \centering
# # \makebox[\textwidth][c]{\begin{subfigure}{1.4\textwidth}
# #   \includegraphics[width=\textwidth]{Images/ODDSim.pdf}
# #   \caption{From left to right, the simulated population, building count, and GNI per capita for a simulated event.}
# #   \label{fig:ODDSim}
# #   \end{subfigure}}
# # \makebox[\textwidth][c]{\begin{subfigure}{1.4\textwidth}
# #   \includegraphics[width=\textwidth]{Images/Sim DispMortBD.png}
# #   \caption{Simulated Damage Data}
# #   \label{fig:DamSim}
# #   \end{subfigure}}
# # \caption{Simulated Data}
# # \end{figure}
# 
# plotODDy<-function(ODDy,zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.5,map="terrain"){
#   
#   if(is.null(bbox)) bbox<-ODDy@bbox
#   
#   mad_map <- get_stamenmap(bbox,source = "stamen",maptype = map,zoom=zoomy)
#   p<- ggplot() #plggmap(mad_map) + xlab("Longitude") + ylab("Latitude")
#   
#   hazard<-rep(NA_real_,length(ODDy@data$hazMean1))
#   for (variable in names(ODDy)[grepl("Mean",names(ODDy))]){
#     tmp<-ODDy[variable]
#     tmp$hazard<-hazard
#     hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
#   }
#   ODDy@data$hazard<-hazard
#   brks<-seq(9,ceiling(2*max(hazard,na.rm = T)),by=1)/2
#   
#   if (var=="GNIc"){
#     ODDy@data[is.na(ODDy@data$ISO3C),var]<-NA
#     dat <- as.data.frame(ODDy)
#     dat <- dat[which(!is.na(dat$GNIc)),]
#     p<-p+geom_tile(data = dat,
#                    mapping = aes(Longitude,Latitude,fill=as.factor(round(dat[[var]]))),alpha=alpha)+ 
#       labs(fill = 'GNI') +
#       scale_fill_brewer(palette = "RdYlGn")
#     #p<-p+geom_contour(data = as.data.frame(ODDy),
#     #                   mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
#     #                   alpha=1.0,breaks = brks) +
#     #   scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent", guide='none') + 
#     #   labs(colour = "Hazard Intensity") 
#     
#     return(p)
#   } else if(var!="hazard")  {
#     ODDy@data[is.na(ODDy@data$ISO3C),var]<-NA
#     
#     p<-p+geom_contour_filled(data = as.data.frame(ODDy),
#                              mapping = aes(Longitude,Latitude,z=as.numeric(ODDy@data[[var]])),alpha=alpha)+ 
#       scale_fill_brewer(palette = "RdYlGn", direction=-1) + 
#       labs(fill = GetVarName(var))
#     if (var %in% c('Disp', 'Mort', 'BuildDam', 'BuildDest')){
#       p<-p+geom_contour(data = as.data.frame(ODDy),
#                         mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
#                         alpha=1.0,breaks = brks) +
#         scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
#         labs(colour = "Hazard Intensity") + guides(fill = guide_legend(order = 1))
#     }
#     
#     
#     # p+geom_contour_filled(data = as.data.frame(ODDy),
#     #                       mapping = aes(Longitude,Latitude,z=1-ODDy@data$tmp),
#     #                       fill="green",alpha=alpha)+ 
#     #   labs(fill = "Hazard>5")
#     
#     return(p)
#   }
#   
#   ODDy@data$hazard[ODDy@data$hazard==0]<-NA
#   
#   p<-p+geom_contour_filled(data = as.data.frame(ODDy),
#                            mapping = aes(Longitude,Latitude,z=hazard),
#                            alpha=alpha,breaks = brks) +
#     scale_fill_brewer(palette = "RdYlGn", direction=-1) + 
#     labs(fill = "Hazard Intensity")
#   
#   # scale_fill_discrete(low = "transparent",high = "red",na.value = "transparent") + 
#   
#   
#   return(p)
#   
# }
# 
# grid.arrange(plotODDy(ODDSim, var='Population') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
#              plotODDy(ODDSim, var='GNI')+ xlim(-0.25,0.25) + ylim(-0.25,0.25), 
#              plotODDy(ODDSim, var='nBuildings')+ xlim(-0.25,0.25) + ylim(-0.25,0.25), nrow=1)
# 
# grid.arrange(plotODDy(ODDSim, var='Disp') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
#              plotODDy(ODDSim, var='Mort') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
#              plotODDy(ODDSim, var='nBD') + xlim(-0.25,0.25) + ylim(-0.25,0.25), nrow=1)



