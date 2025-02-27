#######################################################################
####################### ODDobj CLASS DEFINITION #######################
########### (Oxford-University Disaster Displacement Object) ##########
#######################################################################
# FIELDS:
#   - Helix id number
#   - Data - population, GDP, country ISO3C code, mean hazard intensity, stand.dev hazard intensity
#   - Hazard type
#   - Hazard mapping source
#   - Dates corresponding to hazard intensity (initialised with NA)
#   - Start date (could be evacuation initialisations)
#   - Per Iso3 extract and store indicators
#   - Displacement total data
# METHODS:
#   - Sample hazard intensity (mean,sd) truncated normal distribution
#   - plotGDP, plotPop, plotHaz
#   - PosteriorTotal, PosteriorSample, PosteriorTime
#   - ODDfill: AddGDP, FindISOs, ISO indicators, hazard(haztype)
#   - NOTE: hazard(haztype) does kriging for irregular grid or non-gridded data 
#           and cubic spline interpolation otherwise
#   - AffectedPop(magnitude)
#   - SampleBuildHeight(object@Population,object@buildheightpars)
#   - SampleBuildArea(object@Population,object@buildareapars)
#   - Extract hazard intensity values (I0)
#######################################################################

# source('RCode/Functions.R')
# source('RCode/Model.R')
# source('RCode/GetPopDemo.R')
# source('RCode/GetSocioEconomic.R')
# source('RCode/AddVulnerability.R')
# library(devtools)
# library(parallel)
# library(doParallel)
# library(foreach)
# library(parallelsugar)

checkODD<-function(object) {
  
  if('Population' %in% names(object)) if (any(values(object[['Population']])<0,na.rm = T)) return(F) 
  if('GDP' %in% names(object)) if (any(values(object[['GDP']]) <0,na.rm = T)) return(F) 
  if(any(is.na(object@cIndies$value))) 
    print(paste0("WARNING: missing country indicator elements for ",object@cIndies$iso3[is.na(object@cIndies$value)]))
  
  TRUE
}

Genx0y0<-function(ODDobj){
  
  xo<-ODDobj@coords[1:ODDobj@grid@cells.dim[1],1]
  yo<-ODDobj@coords[1:ODDobj@grid@cells.dim[2]*ODDobj@grid@cells.dim[1]-ODDobj@grid@cells.dim[1]+1,2]
  
  return(list(xo=xo,yo=yo))
}

setClass("ODD", 
         slots = c(dir="character",
                   hazard="character",
                   cIndies="data.frame",
                   #fIndies="list",
                   IDPs="data.frame", # includes measurement dates
                   impact="data.frame",
                   alerts="data.frame",
                   I0="numeric",
                   hazdates="Date",
                   eventid="numeric",
                   predictDisp="data.frame",
                   #modifier="list",
                   polygons="list", 
                   hazinfo="list"),
         contains = "SpatRaster")

ExtractI0poly<-function(HazDat,ODD){

  # Extract contours
  pcontour<-adehabitatMA::getcontour(raster::subset(HazDat,mean>=ODD@I0))
  pcontour<-adehabitatMA::getcontour(HazDat)
  
  
  conts<-data.frame()
  id<-1
  # For each contour, extract points within only if it has a large enough area inside
  for(k in 1:length(pcontour@polygons)) {
    if(!pcontour@polygons[[k]]@Polygons[[1]]@area<1e-3) {
      conts%<>%rbind(data.frame(id=rep(id,length(pcontour@polygons[[k]]@Polygons[[1]]@coords[,1])),
                                Longitude=pcontour@polygons[[k]]@Polygons[[1]]@coords[,1],
                                Latitude=pcontour@polygons[[k]]@Polygons[[1]]@coords[,2]))
      id<-id+1
    }
    # Check for holes in the closed contours: WE DON'T LIKE DONUTS
    if(pcontour@polygons[[k]]@Polygons[[1]]@hole) 
      print(paste0("WARNING: hole in polygon of area ",
                   pcontour@polygons[[k]]@Polygons[[1]]@area,
                   " for ", ODD@hazard," event in countries: ",unique(ODD@data$ISO3C)))
  }
  
  return(conts)
}

interp_overlay <- function(layer, ODD){
  layer_df <- layer$z
  rownames(layer_df) <- layer$x
  colnames(layer_df) <- layer$y

  layer_df <- as.data.frame(as.table(layer_df))
  colnames(layer_df) <- c('Longitude', 'Latitude', 'var')
  layer_df$Longitude <- as.numeric(as.character(layer_df$Longitude))
  layer_df$Latitude <- as.numeric(as.character(layer_df$Latitude))
  
  coords_df <- as.data.frame(ODD@coords)
  coords_df$order <- 1:NROW(coords_df)
  coords_df %<>% merge(layer_df, by=c('Latitude', 'Longitude'))
  coords_df %<>% arrange(order)
  return(coords_df$var)
}

# Add hazard data to the ODD object by interpolating onto the grid using bilinear interpolation
setGeneric("AddHazSDF", function(ODD,lhazSDF) 
  standardGeneric("AddHazSDF") )
setMethod("AddHazSDF", "ODD", function(ODD,lhazSDF){
  
  ODD@I0<-lhazSDF$hazard_info$I0
  lenny<-length(lhazSDF) ; start<-2
  alertscores<-alertlevels<-c() ; dates<-rep(lhazSDF$hazard_info$sdate,lenny-start+1)
  
  polysave<-array(F,dim=c(nrow(ODD),(lenny-start+1)))
  
  for (i in start:lenny){
    print(i-start+1)
    HazDat<-lhazSDF[[i]]
    # Extract detail of specific hazard
    dates[i-start+1]<-HazDat@eventdate
    alertlevels%<>%c(HazDat@alertlevel)
    alertscores%<>%c(HazDat@alertscore)
    
    if(lhazSDF$hazard_info$hazard=="TC"){
      stop('Not yet set up for TC hazards')
    } else {
      haz_interpolated <- resample(HazDat,ODD, method = "bilinear")
      ODD[[paste0("hazMean", i-start+1)]] <- haz_interpolated[['mean']]
      ODD[[paste0("hazSD", i-start+1)]] <- haz_interpolated[['sd']]
    }
  }
  
  ODD@alerts<-data.frame(alertscores=alertscores,alertlevels=alertlevels)
  
  return(ODD)
  
})

readODD <- function(path){
  .Object = new('ODD')
  ODD_list = readRDS(path)
  slotnames <- slotNames(.Object)
  pointer_slot <- ifelse('ptr' %in% slotnames, 'ptr', ifelse('pnt' %in% slotnames, 'pnt', 'pntr'))
  for (slot in slotnames[slotnames!=pointer_slot]){
    slot(.Object, slot) = ODD_list[[slot]]
  }
  slot(.Object, pointer_slot) <- slot(unwrap(ODD_list$spatrast), pointer_slot)
  #.Object@ptr = unwrap(ODD_list$spatrast)@ptr
  names(.Object)[which(names(.Object)=='VALUE')] = 'ISO3C' #sometimes error reading in factor layer
  return(.Object)
}

saveODD <- function(ODD, path){
  ODD_list = list()
  slotnames = slotNames(ODD)
  pointer_slot <- ifelse('ptr' %in% slotnames, 'ptr', ifelse('pnt' %in% slotnames, 'pnt', 'pntr'))
  for (slot in slotnames[slotnames!=pointer_slot]){
    ODD_list[[slot]] = slot(ODD, slot)
  }
  ODD_list$spatrast <- wrap(ODD)
  saveRDS(ODD_list, path)
}

# saveODD <- function(ODD, path){
#   ODD_no_ptr <- ODD
#   ODD_no_ptr$ptr <- NULL
#   saveRDS(list(ODD_no_ptr, wrap(ODD)), path)
# }
# 
# readODD <- function(path){
#   ODD_read <- readRDS(path)
#   ODDy_read[[1]]@ptr = unwrap(ODDy_read[[2]])@ptr
#   return(ODDy_read[[1]])
# }



setMethod(f="initialize", signature="ODD",
          # definition=function(.Object,bbox,lhazSDF,dater=NULL,dir=directory,
          definition=function(.Object,lhazSDF=NULL,DamageData=NULL,dir="./",Model=list(
            WID_perc=   c("p0p10", # Bottom 10% share of Income Distribution
                          "p10p20", # Income share held by 10th - 20th percentiles
                          "p20p30", # Income share held by 20th - 30th percentiles
                          "p30p40", # Income share held by 30th - 40th percentiles
                          "p40p50", # Income share held by 40th - 50th percentiles
                          "p50p60", # Income share held by 50th - 60th percentiles
                          "p60p70", # Income share held by 60th - 70th percentiles
                          "p70p80", # Income share held by 70th - 80th percentiles
                          "p80p90", # Income share held by 80th - 90th percentiles
                          "p90p100") # top 10% share of Income Distribution
          ), agg_level=1) {
            
            if(is.null(lhazSDF)) return(.Object)
            if(!class(lhazSDF[[length(lhazSDF)]])[1]=="HAZARD") return(.Object)
            
            .Object@dir<-dir
            .Object@hazard<-lhazSDF$hazard_info$hazard
            
            if(length(unique(DamageData$eventid))==1) .Object@eventid<-unique(DamageData$eventid)
            if(.Object@hazard%in%c("EQ","TC")){
              if(!is.null(DamageData$gmax)){
                #If using subnational data, @impact is overwritten separately later in GetSubNationalData.R
                .Object@impact<-DamageData%>%group_by(iso3)%>%
                  summarise(gmax=max(gmax),qualifier=qualifierDisp[which.max(gmax)]) 
                .Object@IDPs<-DamageData[,c("sdate","gmax","qualifierDisp")]%>%
                  transmute(date=sdate,IDPs=gmax,qualifier=qualifierDisp)
              }
            } else {
              # THIS IS READY FOR THE IPM APPROACH FOR MID/LONG DURATION HAZARDS
              .Object@IDPs<-DamageData%>%group_by(sdate)%>%summarise(IDPs=max(IDPs),.groups = 'drop_last')%>%
                transmute(date=sdate,IDPs=IDPs)
              # Note that I leave .Object@gmax intentionally blank
            }
            # This bounding box is taken as the minimum region that encompasses all hazard events in HAZARD object:
            bbox<-lhazSDF$hazard_info$bbox
            dater<-min(lhazSDF$hazard_info$sdate)
            .Object@hazdates<-lhazSDF$hazard_info$eventdates
            .Object@hazinfo <-lhazSDF$hazard_info
            
            year<-AsYear(dater)
            
            print("Fetching population data")
            #obj <-GetPopulationBbox(.Object@dir,bbox=bbox)
            obj <- getWorldPop_ODD(.Object@dir, year, bbox, agg_level)
            .Object@ptr <- obj@ptr
            # .Object@file <- obj@file
            # .Object@data <- obj@data
            # .Object@legend <- obj@legend
            # .Object@title <- obj@title
            # .Object@extent <- obj@extent
            # .Object@rotated <- obj@rotated
            # .Object@rotation <- obj@rotation
            # .Object@ncols <- obj@ncols
            # .Object@nrows <- obj@nrows
            # .Object@crs <- obj@crs
            # .Object@srs <- obj@srs
            # .Object@history <- obj@history
            # .Object@z <- obj@z
            # .Object@proj4string <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            #levels(pop[['ISO3C']])[[1]]$VALUE[values(pop[['ISO3C']])]
            print("Adding hazard events")
            # Including minshake polygon per hazard event using getcontour from adehabitatMA package
            .Object%<>%AddHazSDF(lhazSDF)
            
            # Extract empty indices to save time
            inds<-which(!is.na(values(.Object[['Population']])))
            .Object[['ISO3C']][-inds] <- NA
            #levels(pop[['ISO3C']])[[1]]$VALUE[values(pop[['ISO3C']])]
            
            print("Filter spatial data per country")
            #.Object@data$ISO3C<-NA_character_
            #.Object@data$ISO3C[inds]<-coords2country(.Object@coords[inds,])
            
            
            #print("Interpolate population values")
            
            # Note there are as many values returned as iso3c codes (returns as data.frame with columns 'iso3' and 'factor')
            # Popfactors<-InterpPopWB(iso3c,dater, normdate=as.Date(paste(year, "2015-01-01"))
            # 
            # for (iso in iso3c){
            #   indie<-.Object@data$ISO3C==iso & !is.na(.Object@data$ISO3C)
            #   .Object@data$Population[indie]%<>%
            #     multiply_by(Popfactors$factor[Popfactors$iso3==iso])
            # }
            
            # World Income Database (WID) data:
            if(year==AsYear(Sys.Date())) year<-AsYear(Sys.Date())-1
            print("Extract country indicators - WID:")
            
            WID<-GetWID_perc(Model$WID_perc,unique(.Object[['ISO3C']])$ISO3C,year)
            
            
            #stop("Add the full variables to the cIndies data.frame")
            # Bind it all together!
            .Object@cIndies<-WID
            .Object@cIndies$value[which(.Object@cIndies$value==0)] = 0.0001 # 0 values cause issues when applying log transform to GNIc
            
            # Here we add the vulnerabilities used in the linear predictor
            #.Object@ISO3C <- levels(.Object[['ISO3C']])[[1]][[1]]$VALUE
            
            .Object%<>%AddVuln()
            
            
            #print("Fetching GDP-PPP data")
            #.Object%<>%AddGDP(inds)
            
            print("Checking ODD values")
            checkODD(.Object)
            
            return(.Object)
          }
)

# setReplaceMethod( "$", signature("ODD"), function(ODD, name, value) {
#   if( name=="sides" ){
#     ODD@sides <- value
#   }
#   x
# })
# setReplaceMethod( "[", signature("ODD"), function(ODD, name, value) {
#   if( name=="sides" ){
#     ODD@data[] <- value
#   }
#   x
# })

# setGeneric("ExtractCIndy", function(ODD,iso,var)
#   standardGeneric("ExtractCIndy") )
# setMethod("ExtractCIndy", "ODD", function(ODD,iso = NULL,var=NULL){
#   cIndies<-ODD@cIndies
#   if(!is.null(iso)) cIndies%<>%filter(iso3%in%iso)
#   if(!is.null(var)) cIndies%<>%filter(variable%in%var)
#   cIndies
# })

ExtractCIndy<- function(ODD,iso = NULL,var=NULL){
  cIndies<-ODD@cIndies
  if(!is.null(iso)) cIndies%<>%filter(iso3%in%iso)
  if(!is.null(var)) cIndies%<>%filter(percentile%in%var)
  cIndies
}

FormParams<-function(ODD,listy){
  
  
  return(c(listy,list(I0=ODD@I0)))
  #return(c(listy,list(I0=ODD@I0,fIndies=ODD@fIndies)))
  
  
  
  # listy%<>%c(list(I0=ODD@I0,fIndies=ODD@fIndies))
  # Ivars<-unique(ODD@cIndies$variable)
  # Params<-listy
  # tParams<-list()
  # for (iso3c in unique(ODD@cIndies$iso3)){
  #   # Extract the income distribution stochastic diffusion enhancement variable
  #   tParams$Ik<-ExtractCIndy(ODD,iso = iso3c,var = "Ik")$value
  #   # Extract all the constant country specific variables
  #   tParams$var<-ExtractCIndy(ODD,iso = iso3c,
  #                             var = Ivars[!(Ivars=="Ik" | endsWith(Ivars,"p100"))])%>%
  #     dplyr::select(-iso3)
  #   # tParams$var%<>%rbind(data.frame(value=NA,variable="dollar"))
  #   Params[[iso3c]]<-tParams
  # }
  # # names(Params)[(length(listy)+1):length(Params)]<-unique(ODD@cIndies$iso3)
  # return(Params)
}

setGeneric("DispX", function(ODD,Omega,center, Method, output='SampledAgg', event_i=NA)
  standardGeneric("DispX") )
# Code that calculates/predicts the total human displacement 
setMethod("DispX", "ODD", function(ODD,Omega,center,
                                   Method=list(Np=20,cores=8,cap=-300, 
                                               impact_weights=list(displacement=1,mortality=7,buildDam=0.6), 
                                               kernel='energy_score'), output='SampledAgg', event_i = NA
){
  # Samples from the model to simulate the impact of an event.
  # Output type options:
  #     - SampledFull: Simulates and returns the impact in every pixel
  #     - SampledAgg: Simulates the impact in every pixel, but returns a data frame only containing
  #                   simulations aggregated at the same level as the true observations
  #     - SampledTotal: Simulates the impact in every pixel, but returns a data frame containing
  #                   only the total for each impact type (aggregated over all pixels)
  #     - ODDwithSampled: Returns the ODD object with columns Mort, Disp, and BuildDam for the simulated impact
  
  #elapsed_time <- c()
  #start_time <- Sys.time()
  
  # Extract 0D parameters & speed up loop
  Params<-FormParams(ODD,list(Np=Method$Np,center=center))
  Params$I0 <- Model$I0 #some objects have different I0 but don't want this to affect the model
  
  # Income distribution percentiles & extract income percentile  
  SincN<-paste0('p',seq(10,80,10), 'p', seq(20,90,10))
  Sinc<-ExtractCIndy(ODD,var = SincN)
  
  #------------------------- GAUSSIAN PROCESS OVER ERRORS:
  
  # local_errors_all <- matrix(NA, nrow=ncell(ODD), ncol=Method$Np)
  # for (sim in 1:Method$Np){
  #   grid<- list( x= seq( 0,5,length.out=nrow(ODD)), y= seq(0,5,length.out=ncol(ODD))) 
  #   obj<-matern.image.cov( grid=grid, theta=.5, setup=TRUE, aRange=1, smoothness=1)
  #   local_errors <- sim.rf( obj)
  #   local_errors_all[,sim] <- c(t(local_errors))
  # }
  # #image.plot(grid$x, grid$y, local_errors)
  # local_errors2 <-  sim.rf( obj)
  # local_errors3 <-  local_errors + 0.2 * local_errors2
  # cor(c(local_errors3), c(local_errors2))
  # image.plot(grid$x, grid$y, local_errors + local_errors2)
  # 
  #-------------------------- ATTEMPT WITH GSTAT:
  
  #r <- rast(ODD)
  if(is.na(event_i)){
    grid <- as.data.frame(xyFromCell(ODD, 1:ncell(ODD)))  # Extract grid coordinates
    names(grid) <- c("x", "y")  # Name the columns
    vgm_model <- vgm(psill = 1, model = "Mat", range = 0.5, kappa = 1)
    gstat_mod = gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 3)
    
    # set.seed(1)
    #RF_local <- as.matrix(predict(gstat_mod, newdata = grid, nsim = Method$Np)[, 3:(2+Method$Np)])
    
    RF_mort <- Omega$eps$local * as.matrix(predict(gstat_mod, newdata = grid, nsim = Method$Np)[, 3:(2+Method$Np)])
    RF_disp <- Omega$eps$hazard_cor * RF_mort + sqrt(1-Omega$eps$hazard_cor^2) * Omega$eps$local * Omega$eps$hazard_disp / Omega$eps$hazard_mort * as.matrix(predict(gstat_mod, newdata = grid, nsim = Method$Np)[, 3:(2+Method$Np)])
    RF_bd <- Omega$eps$hazard_cor * RF_mort + sqrt(1-Omega$eps$hazard_cor^2) * Omega$eps$local * Omega$eps$hazard_bd / Omega$eps$hazard_mort * as.matrix(predict(gstat_mod, newdata = grid, nsim = Method$Np)[, 3:(2+Method$Np)])
    
  } else {
    #RF_local = Omega$u_local[[event_i]]
    
    RF_mort <- Omega$eps$local * Omega$u_local[[event_i]][,,1]
    RF_disp <- Omega$eps$hazard_cor * RF_mort + sqrt(1-Omega$eps$hazard_cor^2) * Omega$eps$local * Omega$eps$hazard_disp / Omega$eps$hazard_mort * Omega$u_local[[event_i]][,,2]
    RF_bd <- Omega$eps$hazard_cor * RF_mort + sqrt(1-Omega$eps$hazard_cor^2) * Omega$eps$local * Omega$eps$hazard_bd / Omega$eps$hazard_mort * Omega$u_local[[event_i]][,,3]
  }
  # RF_mort <- Omega$eps$local * RF_local
  # RF_disp <- Omega$eps$local * Omega$eps$hazard_disp / Omega$eps$hazard_mort * RF_local
  # RF_bd <- Omega$eps$local * Omega$eps$hazard_bd / Omega$eps$hazard_mort * RF_local
  
  # ODD$mort = RF_local[,1]
  # plot(ODD$mort)
  # RF_disp <- predict(gstat_mod, newdata = grid, nsim = Method$Np)[, 3:(2+Method$Np)]
  # RF_bd <- predict(gstat_mod, newdata = grid, nsim = Method$Np)[, 3:(2+Method$Np)]
  # 
  # errors <- aperm(abind(
  #   Omega$eps$hazard_mort * RF_mort, 
  #   Omega$eps$hazard_disp * (Omega$eps$hazard_cor * RF_mort + RF_disp),
  #   Omega$eps$hazard_bd * (Omega$eps$hazard_cor * RF_mort + RF_bd),
  #   along=3
  # ), c(1,3,2))
  
  #------------------------------------------------------
  
  ODD_df <- as.data.frame(ODD, na.rm=F)
  
  #ODD_df$ISO3C <- levels(ODD[['ISO3C']])[[1]]$VALUE[ODD_df$ISO3C]
 
  # Speed-up calculation (through accurate cpu-work distribution) to only values that are not NA
  nans_haz <- apply(ODD_df[,grep('hazMean', names(ODD_df)), drop=F], 1, function(x) all(is.na(x)))
  
  notnans<-which(!(is.na(ODD_df$Population) | is.na(ODD_df$ISO3C) | is.na(ODD_df$SHDI) | nans_haz))

  # Calculate non-local linear predictor values
  LP<-GetLP(ODD_df,Omega,Params,Sinc,notnans, split_GNI=F)
  LP_buildings <- GetLP(ODD_df,Omega,Params,Sinc,notnans, split_GNI=F)
  
  #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, GetLP = finish_time-start_time); start_time <- Sys.time()

  BD_data_present <- !is.null(ODD_df$nBuildings)
  hrange<-grep("hazMean",names(ODD_df),value = T)
  hrange_order <- order(paste(ODD@hazinfo$eventdates, ODD@hazinfo$eventtimes))
  
  event_lp <- getLP_event(ODD@hazinfo, Omega, Params)
  
  cov_mort_disp = Omega$eps$hazard_cor * Omega$eps_adj$hazard_mort * Omega$eps_adj$hazard_disp
  cov_mort_bd = Omega$eps$hazard_cor * Omega$eps_adj$hazard_mort * Omega$eps_adj$hazard_bd
  cov_disp_bd = Omega$eps$hazard_cor * Omega$eps_adj$hazard_disp * Omega$eps_adj$hazard_bd
  covar_matrix = cbind(c(Omega$eps_adj$hazard_mort^2, cov_mort_disp, cov_mort_bd), c(0, Omega$eps_adj$hazard_disp^2, cov_disp_bd), c(0, 0, Omega$eps_adj$hazard_bd^2))
  covar_matrix[upper.tri(covar_matrix)] = covar_matrix[lower.tri(covar_matrix)]
  
  
  covar_matrix_local = covar_matrix * Omega$eps_adj$local
  # covar_matrix_local[1,2] = covar_matrix_local[2,1]  = 1 * sqrt(covar_matrix_local[1,1] * covar_matrix_local[2,2])
  # covar_matrix_local[1,3] = covar_matrix_local[3,1]  =  1 * sqrt(covar_matrix_local[1,1] * covar_matrix_local[3,3])
  # covar_matrix_local[2,3] = covar_matrix_local[3,2]  =  1 * sqrt(covar_matrix_local[2,2] * covar_matrix_local[3,3])
  # 
  # eps_local_ij <- array(0, dim=c(length(hrange), 3, Method$Np))
  # for (i in 1:Method$Np){
  #   eps_local_ij[,,i] <- rmvnorm(length(hrange), rep(0,3), sigma=covar_matrix_local)
  # }
  
  #eps_event <- array(0, dim=c(3, Method$Np))
  if (is.na(event_i)){
    eps_event <- t(rmvnorm(Method$Np, rep(0, 3), sigma=covar_matrix))
  } else {
    eps_event <-  chol(covar_matrix) %*% t(Omega$u[event_i,,])
  }
  #eps_event[,] <- 0
  
  #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); start_time <- Sys.time()

  #Function to predict damage per gridpoint
  CalcDam<-function(ij){
    
    # elapsed_time <- c()
    # start_time <- Sys.time()

    #finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochasticregen = finish_time-start_time); start_time <- Sys.time()
    
    #eps_local_long <- rmvnorm(length(hrange)*Method$Np, rep(0,3), sigma=covar_matrix_local)
    #eps_local_ij <- aperm(array(eps_local_long, dim=c(length(hrange), Method$Np, 3)), c(1,3,2))
    #eps_local_ij[1,,] <- cbind(local_errors_all[ij,],local_errors_all[ij,],local_errors_all[ij,])
    #eps_local_ij <- array(0, dim=c(length(hrange),3,Method$Np))
    
    #eps_local_ij <- RF_local[ij,]
    eps_local_ij <- rbind(RF_mort[ij,], RF_disp[ij,], RF_bd[ij,])
    
    
    ##SLOWER:
    # for (i in 1:Method$Np){
    #   eps_local_ij[,,i] <- rmvnorm(length(hrange), rep(0,3), sigma=covar_matrix_local)
    # }
    # notnans_ij <- which(notnans==ij)
    #eps_local_ij <- adrop(eps_local_transf[,,,notnans_ij, drop=F], drop=4)
    
    locallinp<- LP[ij] # LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]]*LP$Plinp[ij]*LP$linp[[iso3c]] 
    
    # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochasticreshape = finish_time-start_time); start_time <- Sys.time()
    
    # Sample population per income distribution (Assumes 8 percentiles):
    # Population is split evenly between the income quantiles, with remainders randomly allocated between
    lPopS <- rep(floor(ODD_df$Population[ij]), Method$Np) #+ rbinom(Method$Np,1, ODD_df$Population[ij] %% 1) 
    #lPopS <- SplitSamplePop(Pop=ODD@data$Population[ij],Method$Np) #matrix(round(ODD@data$Population[ij]/length(locallinp)), nrow=length(locallinp), ncol = Method$Np)
    #lPopS <- matrix(round(ODD@data$Population[ij]/ 8), nrow=8, ncol=Method$Np)
    
    # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); start_time <- Sys.time()
    
    
    lPopDisp <- rep(0, Method$Np)
    lPopMort <- rep(0, Method$Np)
    tPop <-array(0,c(3, Method$Np)) #row 1 = tDisp, #row 2 = tMort, #row 3 = tRem
    tPop[3,]=lPopS
    
    # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); 
    
    for(h_i in hrange_order){
      
      #eps_local_ij[h_i,,] <- errors[ij,,] #array(NA, dim=c(1, 3, Method$Np))
      #eps_local_ij[h_i,,] <- as.matrix(rbind(RF_local[ij,], RF_local[ij,], RF_local[ij,]))
      
      start_time <- Sys.time()
      h <- hrange[h_i]

      if(is.na(ODD_df[ij,h])) next
      
      nonzero_pop <- which(lPopS != 0)
      if (length(nonzero_pop)==0) next
      
      
      # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); start_time <- Sys.time()
      
      ## Sample hazard Intensity 
      ## Doesn't work very well as don't have the covariance of the errors, so must instead assume they are independent
      ## So instead just use the mean value
      # I_ij<-rnorm(n = Method$Np,
      #             mean = ODD@data[ij,paste0("hazMean",h)],
      #             sd = ODD@data[ij,paste0("hazSD",h)]/10)
      I_ij<-ODD_df[ij,h]
      Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=1),Omega) + locallinp + event_lp[h_i], error=function(e) NA) #+ rep(eps_local[h_i,], each=8), error=function(e) NA)
      
      # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, damcalc = finish_time-start_time); start_time <- Sys.time()
      
      #D_MortDisp <- D_MortDisp_calc(Damage, Omega, eps_event[1:2, nonzero_pop[,2], drop=F] + adrop(eps_local_ij[h_i,1:2,nonzero_pop[,2], drop=F], drop = 1)) #First row of D_MortDisp is D_Mort, second row is D_Disp
      #D_MortDisp <- D_MortDisp_calc(Damage, Omega, sweep(eps_event[1:2, nonzero_pop[,2], drop=F], 2, eps_local_ij[nonzero_pop[,2], drop=F], "+")) #First row of D_MortDisp is D_Mort, second row is D_Disp
      D_MortDisp <- D_MortDisp_calc(rep(Damage, length(nonzero_pop)), Omega, eps_event[1:2, nonzero_pop, drop=F] + eps_local_ij[1:2, nonzero_pop, drop=F]) #First row of D_MortDisp is D_Mort, second row is D_Disp
      
      
      # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, mortdisp = finish_time-start_time); start_time <- Sys.time()
      
      D_Rem <- pmax(0, 1 - D_MortDisp[2,] - D_MortDisp[1,]) #probability of neither death nor displacement. Use pmax to avoid errors caused by numerical accuracy.

      # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, d_rem = finish_time-start_time); start_time <- Sys.time()
      
      Dam <- Fbdam(lPopS[nonzero_pop], D_MortDisp[2,], D_MortDisp[1,], D_Rem)
      
      # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, fbdam = finish_time-start_time); start_time <- Sys.time()
      
      lPopS[nonzero_pop] <- Dam[3,]
      lPopDisp[nonzero_pop] <- lPopDisp[nonzero_pop] + Dam[1,]
      lPopMort[nonzero_pop] <- lPopMort[nonzero_pop] + Dam[2,]
      tPop[3,] <- lPopS
      
      # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, sum_pops = finish_time-start_time); 
      
      # # This is a bit clearer than the above but slower:
      # Accumulate the number of people displaced/deceased, but don't accumulate the remaining population
      # tPop[3,ind]<-0
      # for (s in 1:length(SincN)){ #Separate into income distributions (as each have 10% of population, order doesn't matter)
      #   if(all(lPopS[s,]==0)) next
      #   # Predict damage at coordinate {i,j} (vector with MC particles)
      #   Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=sum(ind)),Omega) + locallinp[s] + eps_event[h_i,ind], error=function(e) NA)
      #   if(any(is.na(Damage))) print(ij)
      # 
      #   #LOOSEEND: Include [ind] here
      #   D_MortDisp <- D_MortDisp_calc(Damage, Omega) #First row of D_MortDisp is D_Mort, second row is D_Disp
      #   D_Rem <- pmax(0, 1 - D_MortDisp[2,] - D_MortDisp[1,]) #probability of neither death nor displacement. Use pmax to avoid errors caused by numerical accuracy.
      # 
      #   tPop[,ind]<-tPop[,ind] + Fbdam(lPopS[s,ind],D_MortDisp[2,], D_MortDisp[1,], D_Rem)
      # }
    }
    
    # start_time <- Sys.time()
    
    tPop[1,] <- lPopDisp
    tPop[2,] <- lPopMort
    
    #ensure the total displaced, deceased or remaining does not exceed total population
    tPop[tPop>ODD_df$Population[ij]] <- floor(ODD_df$Population[ij])
    
    # if (length(elapsed_time)==9) {return(elapsed_time)
    # } else {(return(rep(0,9)))}
    
    #if no building destruction data:
    if(!BD_data_present) return(list(samples=rbind(tPop[1:2,, drop=FALSE], rep(NA, Method$Np))))#, 

    locallinp_buildings <- LP_buildings[ij]
    
    nUnaff = rep(ODD_df$nBuildings[ij], Method$Np)
    nDam = rep(0, Method$Np)
    
    #finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); start_time <- Sys.time()
    
    for (h_i in hrange_order){
      h <- hrange[h_i]
      if(is.na(ODD_df[ij,h])) next
      if(all(nUnaff==0)) break #if no remaining buildings, skip modelling

      I_ij<-ODD_df[ij,h]
      Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Method$Np),Omega) + locallinp_buildings + event_lp[h_i], error=function(e) NA) #+ eps_local[h_i,], error=function(e) NA) #calculate unscaled damage (excluding GDP)
 
      D_Dam <- D_Dam_calc(Damage, Omega, eps_event[3,] + eps_local_ij[3,])
      
      # Accumulate the number of buildings damaged/destroyed, but not the number of buildings remaining
      nDam_new <- rbinom(Method$Np, nUnaff, D_Dam)
      nUnaff <- nUnaff - nDam_new
      nDam <- nDam + nDam_new
    
    }
    
    #finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time);
    
    #return(elapsed_time)
    
    
    return(list(samples = rbind(tPop[1:2,,drop=FALSE], nDam[1:Method$Np])))
  }
  
  Dam<-array(0,c(nrow(ODD_df),Method$Np,3)) # Dam[,,1] = Displacement, Dam[,,2] = Mortality, Dam[,,3] = Buildings Damaged
  Dam_means<-array(0,c(nrow(ODD_df),3))
  
  if(Method$NestedCores>1) { 
    CalcDam_out <- mclapply(X = notnans,FUN = CalcDam,mc.cores = Method$NestedCores)
  } else  {CalcDam_out <- lapply(X = notnans,FUN = CalcDam)}
  
  Dam[notnans,,]<-aperm(simplify2array(lapply(CalcDam_out, function(x) x$samples)), perm=c(3,2,1))
  
  if (output=='SampledFull'){
    return(Dam)
  } else if (output == 'SampledTotal'){
    SampledTot <- colSums(Dam)
    df_SampledTot <- list()
    impact_types <- unique(ODD@impact$impact)
    
    for (impact_type in impact_types){
      polygon_names <- unlist(lapply(ODD@polygons[ODD@impact$polygon], function(x) x$name))
      if (any(tolower(polygon_names[which(ODD@impact$impact==impact_type)]) %in% c('tot', 'total'))){
        nonmatch <- which(!tolower(polygon_names[which(ODD@impact$impact==impact_type)]) %in% c('tot', 'total'))
        if (length(nonmatch)>0){
          ODD@impact <- ODD@impact[-which(ODD@impact$impact==impact_type)[nonmatch],] # in the case of total and subnational data, remove the subnational
        }
      }
    }
    
    #Many ODD objects don't contain 'total' impact values (to avoid double counting), so we need to obtain these. 
    #We combine the polygons with observations so long as the proportion of overlapping pixels is less than 10%
    #and the total coverage of the exposed area is greater than 90%. Then sum the observations across these polygons.
    observed_total=rep(NA, length(impact_types))
    exposed_haz <- which(apply(ODD_df[,grep('hazMean', names(ODD_df)), drop=F], 1, function(row) any(!is.na(row))) & !is.na(ODD_df$ISO3C) & (ODD_df$Population>0))
    get_overlap_coverage <- function(impact_type){
      indexes_list <- lapply(ODD@polygons[ODD@impact$polygon[which(ODD@impact$impact==impact_type)]], function(x) x$indexes)
      universal_set <- Reduce(union, indexes_list)
      overlap <- Reduce(intersect, indexes_list)
      overlap <- intersect(overlap, exposed_haz)
      prop_overlap <- ifelse(length(indexes_list)>1,length(overlap)/length(universal_set),0)
      prop_coverage <- length(intersect(unique(universal_set), exposed_haz))/length(exposed_haz)#sum(!is.na(ODD_df$ISO3C))
      return(c(prop_overlap, prop_coverage))
    }
    
    overlap_coverage <- sapply( impact_types,get_overlap_coverage)
      
    for (i in 1:length(impact_types)){
      observed_total[i] = sum(ODD@impact$observed[which(ODD@impact$impact==impact_types[i])])
    }
  
    
    for (i in 1:NROW(SampledTot)){
      df_SampledTot[[i]] <- data.frame(event_id = ODD@impact$event_id[1],
                                       polygon=0,
                                       iso3 = paste0(unique(ODD@impact$iso3), collapse=' '),
                                       impact=impact_types, 
                                       sampled=SampledTot[i,match(impact_types, c('displacement', 'mortality','buildDam'))], 
                                       observed=observed_total,
                                       qualifier='total',
                                       overlap = overlap_coverage[1,],
                                       coverage = overlap_coverage[2,],
                                       inferred=F)
      df_SampledTot[[i]] %<>% filter(overlap < 0.1 & coverage > 0.9)
    }

    
    return(df_SampledTot)
  }
  
  
  #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, dam_sample = finish_time-start_time); start_time <- Sys.time()
  
  funcy<-function(i, impact_weights=Method$impact_weights, kernel=Method$kernel, cap=Method$cap, polygons_indexes=ODD@polygons) {
    
    # funcy() aggregates the simulate impact (for each impact type) across each polygon, when
    # there is a matching observation for that impact type and polygon
    
    tmp<-data.frame(iso3=ODD_df$ISO3C, displacement=Dam[,i,1], mortality=Dam[,i,2], buildDam=Dam[,i,3])#, 
                                    #mort_mean=Dam_means[,1], disp_mean=Dam_means[,2], buildDam_mean=Dam_means[,3])
    impact_sampled<-data.frame(polygon = numeric(), impact = character(), sampled = numeric(), mean=numeric())
    
    for (polygon_id in unique(ODD@impact$polygon)){
      polygon_impacts <- ODD@impact$impact[which(ODD@impact$polygon==polygon_id)]
      for (impact in polygon_impacts){
        if (impact == 'mortality'){
          impact_sampled <- rbind(impact_sampled, data.frame(polygon=polygon_id,
                                                 impact=impact,
                                                 sampled=floor(sum(tmp[polygons_indexes[[polygon_id]]$indexes,impact]  * polygons_indexes[[polygon_id]]$weights, na.rm=T))))
        } else if (impact=='displacement'){
          impact_sampled <- rbind(impact_sampled, data.frame(polygon=polygon_id,
                                                 impact=impact,
                                                 sampled=floor(sum(tmp[polygons_indexes[[polygon_id]]$indexes,impact]  * polygons_indexes[[polygon_id]]$weights, na.rm=T))))
        } else {
          impact_sampled <- rbind(impact_sampled, data.frame(polygon=polygon_id,
                                                 impact=impact,
                                                 sampled=floor(sum(tmp[polygons_indexes[[polygon_id]]$indexes,impact]  * polygons_indexes[[polygon_id]]$weights, na.rm=T))))
        }
      }
    }
    impact_obs_sampled <- arrange(merge(impact_sampled, ODD@impact, by=c("polygon", "impact")),desc(observed)) 
    return(impact_obs_sampled)
  }
  
  if (output == 'ODDwithSampled'){ #usually used for generating simulated data
    ODD[['Disp']]<- Dam[,1,1]  
    ODD[['Mort']]<- Dam[,1,2]
    ODD[['BuildDam']] <- Dam[,1,3]  
    ODD@predictDisp<-funcy(1) 
    return(ODD)
  }
  
  for (i in 1:length(ODD@polygons)){ #weights can be used to downweight pixels that are only partly in a polygon
    if(is.null(ODD@polygons[[i]]$weights)){
      ODD@polygons[[i]]$weights <- rep(1, length=length(ODD@polygons[[i]]$indexes))
    }
  }
  
  ## IF TIMING: --------------------------------
  # dummy <- lapply(1:Method$Np, funcy, LLout=F)
  # finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, agg_dam = finish_time-start_time); start_time <- Sys.time()
  # return(elapsed_time)
  # --------------------------------------------
  
  if(is.null(ODD_df$nBuildings)){
    ODD@impact <- ODD@impact[!ODD@impact$impact %in% c('buildDam', 'buildDest', 'buildDamDest'),]
  }
  
  if(output == 'SampledAgg'){
    return(lapply(1:Method$Np, funcy))
  }
  
  stop('Output type selected as "output" not recognised, please select either SampledAgg, SampledTotal,
       SampledFull, or ODDwithSampled')
  
  return(ODD)
  
})


GroupODDyBoundaries<-function(ODDy,boundaries){
  
  for (i in 1:length(boundaries$geometry)){

    insiders<-rep(F,length(ODDy$Population))
    
    nc <- as(boundaries$geometry[i], "Spatial")
    for(j in 1:length(nc@polygons[[1]]@Polygons)) {
      
      coordz<-nc@polygons[[1]]@Polygons[[j]]@coords
      
      insiders[!is.na(ODDy$Population)]<-insiders[!is.na(ODDy$Population)] |
      point.in.polygon(ODDy@coords[!is.na(ODDy$Population),1],
                       ODDy@coords[!is.na(ODDy$Population),2],
                       coordz[,1],
                       coordz[,2])>0

    }
    
    # Allocate damaged buildings to each region boundary
    ODDy$Region[insiders]<-boundaries$name[i]
    
  }
  
  ODDy$Region[is.na(ODDy$Population)]<-NA
  
  return(ODDy)
  
}

ODDyAggPerRegion<-function(ODDy,featurez="Disp",name="RegionAggDisp"){
  
  ODDy@data$featurez<-ODDy@data[,featurez]
  # Need to create a DF with aggregates per region, with 'Region' as a column name
  aggy<-ODDy@data%>%group_by(Region)%>%summarise(Aggy=sum(featurez,na.rm = T),.groups="drop")
  ODDy@data$featurez<-NULL
  # Then need to join the DF with ODDy$data to be able to plot it
  output<-plyr::join(ODDy@data,aggy,by="Region")
  output[name]<-output$Aggy; output$Aggy<-NULL
  
  ODDy@data<-output
  return(ODDy)
  
}

ExtractDispData_ODD<-function(haz){
  
  # Load ODD files
  folderin<-paste0(dir,"IIDIPUS_Results/ODDobjects")
  ufiles<-list.files(path=folderin,pattern=haz,recursive = T,ignore.case = T)
  
  DispData<-data.frame()
  for(i in 1:length(ufiles)){
    # Extract the ODD object
    ODDy<-readRDS(paste0(dir,"IIDIPUS_Results/ODDobjects/",ufiles[i]))
    # iso3, start date and event_id values for each hazard
    DispData%<>%rbind(data.frame(iso3=ODDy@gmax$iso3[1],
                                 sdate=ODDy@hazdates[1],
                                 event_id=ODDy@eventid))
  }
  
  return(DispData)
  
}

GetVarName<-function(varname){
  
  lookup<-list("Population"="CIESIN Population Data",
               "Disp" = "Displaced Population",
               "Mort" = "Mortality",
               "GDP" = "GDP [USD-2015]",
               "FBPop" = "HRSL Population",
               "RWI" = "FB Relative Wealth Index",
               "DamBDosm" = "Damaged Buildings (OSM)",
               "DamBDFB" = "Damaged Buildings (HRSL)")
  
  #Max's lookup for easier plotting
  lookup<-list("Population"="Population  ",
               "Disp" = "Displacement ",
               "Mort" = "Mortality   ",
               "nBuildings" = "# Buildings ",
               "nBD" = "# Build Dest",
               "GDP" = "GDP         ",
               "FBPop" = "HRSL Population",
               "RWI" = "FB Relative Wealth Index",
               "DamBDosm" = "Damaged Buildings (OSM)",
               "DamBDFB" = "Damaged Buildings (HRSL)",
               "MedianDisp" = "Median", 
               "RangeDisp" = "Range")
  
  as.character(lookup[varname])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#==============================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#====================Functions for plotting ODD objects =======================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#==============================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Plot the ODD object with the GADM administrative boundaries as the background:
# plotODDy_GADM <- function(ODDy,zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.7){
#   gadm_level=2
#   if (is.null(bbox)) bbox <- ODDy@bbox
#   iso3_unique <- unique(ODDy$ISO3C)
#   iso3_unique <- iso3_unique[!is.na(iso3_unique)]
#   gadm_iso <- as(geodata::gadm(country=iso3_unique[1], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial')
#   if (length(iso3_unique) > 1){
#     for (i in 2:length(iso3_unique)){
#       gadm_iso %<>% rbind(as(geodata::gadm(country=iso3_unique[i], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial'))
#     }
#   }
#   #gadm_iso <- gSimplify(gadm_iso, 0.01)
#   gadm_iso <- intersect(gadm_iso, bbox)
#   gadm_map <- fortify(gadm_iso)
#   
#   gg <- ggplot()
#   gg <- gg + geom_map(map=gadm_map, data=gadm_map, aes(map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
#   gg <- gg + coord_map() + geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='white', color='black')
#   p <- gg + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
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
#   if(var!="hazard")  {
#     ODDy@data[is.na(ODDy@data$ISO3C),var]<-NA
#     
#     fill_colours <- c('white', 'turquoise', 'deepskyblue', 'purple', 'purple4')
#     
#     if (is.null(breakings)) {
#       breakings_max = max(ODDy@data[[var]], na.rm=T) 
#       log_breakings_max = round(log(breakings_max, 10))
#       breakings = c(0, 10^seq(1, log_breakings_max, 1))
#     }
#     options(scipen=10000)
#     p <- p+geom_raster(data=as.data.frame(ODDy) %>% filter(!is.na(!!ensym(var))),aes(Longitude,Latitude,fill=!!ensym(var)), #as.numeric(!!ensym(var))),
#                        alpha=0.9,interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
#       scale_fill_gradientn(colors=fill_colours, breaks=breakings, na.value = "transparent", trans = 'pseudo_log') +
#       labs(fill = var)+xlab("Longitude") + ylab("Latitude") + geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='transparent', color='grey', lwd=0.1)
#     
#     p<-p+geom_contour(data = as.data.frame(ODDy),
#                       mapping = aes(Longitude,Latitude,z=hazard,colour=after_stat(level)),
#                       alpha=1,breaks = brks, lwd=0.7) +
#       scale_colour_gradient(low = "#FFFF99",high = "red",na.value = "transparent") + 
#       labs(colour = "Hazard Intensity")
#     
#     
#     return(p)
#   }
#   
#   ODDy@data$hazard[ODDy@data$hazard==0]<-NA
#   p <- p+geom_contour_filled(data = as.data.frame(ODDy),
#                              mapping = aes(Longitude,Latitude,z=hazard),
#                              alpha=1,breaks = brks, lwd=0.5) +
#     scale_colour_gradient(low = "#FFFF99",high = "red",na.value = "transparent") + 
#     labs(fill = "Hazard Intensity")
#   
# }


# plotODDy_GADM <- function(ODDy,zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.7,map="terrain"){
#   
#   if(is.null(breakings) & (var=="Population" | var=="Disp" | var=='Population2')) breakings<-c(0,1,5,10,50,100,500,1000, 2000, 5000, 50000)
#   
#   
#   bbox <- ODDy@bbox
#   gadm_iso <- getData("GADM", country="NZL", level=2)
#   gadm_iso <- gadm_iso#gSimplify(gadm_iso, 0.01)
#   gadm_iso <- intersect(gadm_iso, bbox)
#   gadm_map <- fortify(gadm_iso)
# 
#   gg <- ggplot()
#   gg <- gg + geom_map(map=gadm_map, data=gadm_map, aes(x=long, y=lat, map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
#   gg <- gg + coord_map() + geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='white', color='black')
#   p <- gg + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
# 
#   p
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
#   if(var!="hazard")  {
#     ODDy@data[is.na(ODDy@data$ISO3C),var]<-NA
#     
#     p<-p+geom_contour_filled(data = as.data.frame(ODDy),
#                              mapping = aes(Longitude,Latitude,z=ODDy@data[[var]]),
#                              breaks=breakings,alpha=alpha)+ 
#       labs(fill = GetVarName(var))
#     p<-p+geom_contour(data = as.data.frame(ODDy),
#                       mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
#                       alpha=1.0,breaks = brks, lwd=1.2) +
#       scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
#       labs(colour = "Hazard Intensity")
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
#     # scale_fill_discrete(low = "transparent",high = "red",na.value = "transparent") + 
#     labs(fill = "Hazard Intensity")
#   p<-p+geom_contour(data = as.data.frame(ODDy),
#                     mapping = aes(Longitude,Latitude,z=hazard),
#                     alpha=0.8,breaks = c(6.0),colour="red")
#   
#   
# }

plotODDy <-function(ODDy,zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.5,map="terrain", api_key_loc='/home/manderso/Documents/Miscellaneous/API_Keys/stadiamap_key_mal'){
  
  if(is.null(breakings) & (var=="Population" | var=="Disp" | var=='Population2')) breakings<-c(0,1,5,10,50,100,500,1000, 2000, 5000, 50000)
  
  if(is.null(bbox)) bbox<-ODDy@bbox
  
  if(!file.exists(api_key_loc)){
    warning('You need an API key from StadiaMaps to get terrain background. 
            Get one here: https://client.stadiamaps.com/signup/ and save your key (as a string) as an RDS in a local dir. 
            Pass the location of that directory using the api_key_loc argument')
    return()
  }
  stadiamap_api_key <- readRDS(api_key_loc)
  register_stadiamaps(key = stadiamap_api_key)
  
  #mad_map <- get_stamenmap(bbox,source = "terrain",maptype = map,zoom=zoomy)
  mad_map <- get_stadiamap(bbox = bbox, zoom = zoomy, maptype = "stamen_terrain_background")
  p<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")
  
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
    
    p<-p+geom_contour_filled(data = as.data.frame(ODDy),
                             mapping = aes(Longitude,Latitude,z=ODDy@data[[var]]),
                             breaks=breakings,alpha=alpha)+ 
      labs(fill = GetVarName(var))
    p<-p+geom_contour(data = as.data.frame(ODDy),
                      mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
                      alpha=1.0,breaks = brks, size=1) +
      scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
      labs(colour = "Hazard Intensity")
    
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
  
  
  return(p)
  
}

# plotODDy_GADM <- function(ODDy, zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.7, gadm_level=2){
#   #Plots background as GADM regions rather than terrain:
#   
#   bbox <- ODDy@bbox
#   #gadm_iso <- getData("GADM", country="NZL", level=2)
#   iso3_unique <- unique(ODDy$ISO3C)
#   iso3_unique <- iso3_unique[!is.na(iso3_unique)]
#   gadm_iso <- as(geodata::gadm(country=iso3_unique[1], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial')
#   if (length(iso3_unique) > 1){
#     for (i in 2:length(iso3_unique)){
#       gadm_iso %<>% rbind(as(geodata::gadm(country=iso3_unique[i], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial'))
#     }
#   }
#   #gadm_iso <- gSimplify(gadm_iso, 0.01)
#   gadm_iso <- intersect(gadm_iso, bbox)
#   gadm_map <- fortify(gadm_iso)
#   
#   gg <- ggplot()  p <- gg + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
#   
#   gg <- gg + geom_map(map=gadm_map, data=gadm_map, aes(x=long, y=lat, map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
#   gg <- gg + coord_map() + geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='white', color='black')
#   p <- gg + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
#   
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
#   if(var!="hazard")  {
#     ODDy@data[is.na(ODDy@data$ISO3C),var]<-NA
#     
#     p <- p+geom_tile(data = as.data.frame(ODDy),
#                      mapping = aes(Longitude,Latitude,fill=ODDy@data[[var]]+0.1),alpha=0.8, width=ODDy@grid@cellsize[1]*5,
#                      height=ODDy@grid@cellsize[2]*5) + scale_fill_viridis( trans = "log10", labels = function(x) sprintf("%.0f", x)) + labs(fill = ifelse(GetVarName(var)=='NULL', var, GetVarName(var)))
#     
#     p<-p+geom_contour(data = as.data.frame(ODDy),
#                       mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
#                       alpha=1.0,breaks = brks, lwd=0.8) +
#       scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
#       labs(colour = "Hazard Intensity")
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
#     # scale_fill_discrete(low = "transparent",high = "red",na.value = "transparent") + 
#     labs(fill = "Hazard Intensity")
#   p<-p+geom_contour(data = as.data.frame(ODDy),
#                     mapping = aes(Longitude,Latitude,z=hazard),
#                     alpha=0.8,breaks = c(6.0),colour="red")
#   
#   
# }

plotODDy_GADM <- function(ODDy, var, gadm_level=2, haz_legend=F, var_legend=T, var_discrete=F, log_legend=F){
  
  plot_df <- as.data.frame(ODDy, xy=T, na.rm=F)
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
  
  # if (log_legend){
  #   p <- p + geom_raster(dat=plot_df, aes(x=Longitude, y=Latitude, fill=!!sym(var)+0.1), alpha=0.75) # quick fix but must be a better way to do log scale with 0s
  # } else {
  #   p <- p + geom_raster(dat=plot_df, aes(x=Longitude, y=Latitude, fill=!!sym(var)), alpha=0.75)
  # }
  p <- p + geom_raster(dat=plot_df, aes(x=Longitude, y=Latitude, fill=!!sym(var)), alpha=0.75)
  p <- p + coord_equal() +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    theme_minimal() + 
    theme(
      axis.title = element_text(family = "Liberation Serif", size=12),  
      legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
      legend.title = element_text(family = "Liberation Serif", size=12)
    ) + 
    geom_contour(dat=plot_df, aes(Longitude,Latitude,z=hazMean1,colour=..level..),
                 alpha=0.7, lwd=0.8) +
    scale_color_gradientn(colors = c("transparent","#fc9272", "#ef3b2c"))
  if (haz_legend){
    p <- p + scale_color_gradientn(colors = c("transparent","#fc9272", "#ef3b2c")) + labs(colour = "Hazard Intensity      ")
  } else {
    p <- p + scale_color_gradientn(colors = c("transparent","#fc9272", "#ef3b2c"), guide='none')
  }
  if (var_discrete){
    p <- p + scale_fill_viridis_d()
  } else {
    if (var_legend){
      if (log_legend){
        if (max(plot_df[,var], na.rm=T)> 100000){
          breakings = c(0, 1000, 10000, 100000, 1000000)
        } else {
          breakings = c(0, 1000, 10000, 100000)
        }
        p <- p + scale_fill_viridis( trans = scales::pseudo_log_trans(base = 10, sigma = 100), 
                                     breaks=breakings, labels = function(x) scales::comma(x))
      }
      else {
        if (all(plot_df[, var]==0, na.rm=T)){
          p <- p + scale_fill_viridis_c( limits = c(0, 1),breaks=c(0,1), labels=c(0,1))
        } else {
          p <- p + scale_fill_viridis_c()
        }
        
      }
    } else {
      p <- p + scale_fill_viridis_c(guide='none')
    }
  }
  return(p)
}

plotODDPolygons <-function(ODDy,admin_level=2, bbox=NULL){

  #if(is.null(breakings) & (var=="Population" | var=="Disp" | var=='Population2')) breakings<-c(0,1,5,10,50,100,500,1000, 2000, 5000, 50000)
  
  if(is.null(bbox)) bbox<-as.numeric(c(ext(ODDy)[1], ext(ODDy)[3], ext(ODDy)[2], ext(ODDy)[4]))

  poly_match <- rep(NA, ncell(ODDy))
  match_polygon <- function(index){
    matches <- unlist(lapply(ODDy@polygons, function(x){
      if(str_count(x$name, ',') != admin_level) return(c())
      poly_match <- which(x$indexes==index)
      if (length(poly_match)==0) return(c())
      if (x$weights[poly_match]>0.5){
        return(x$name)
      } else {
        return(c())
      }
    } ))
    return(ifelse(length(matches)>0,matches,NA))
  }
  
  for (i in 1:length(poly_match)){
    if (i %% 100==0){print(paste('Identifying poly of pixel', i))}
    poly_match[i] <- match_polygon(i)
  }
  
  ODDy[['poly_match']] = poly_match
  p <- plot(ODDy[['poly_match']])
  
  return(p)
}

plotODDyAgg <- function(ODDyAgg, ODDy=NULL, zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.7,map="terrain", gadm_level=2){
  # Plots aggregated impact/exposure across GADM regions rather than by pixel
  # Will sum over the variable (so not appropriate for covariates e.g. SHDI or EQFreq)
  
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


# MakeODDPlots<-function(ODDy, input=NULL,zoomer=7, bbox=NULL){
#   
#     ODDy@data$Disp[ODDy@data$Disp<1]<-0
#   
#     mad_map <- get_stamenmap(ODDy@bbox,source = "stamen",maptype = "terrain",zoom=zoomer)
#     q<-ggmap(mad_map)
#     p1<-q+ geom_raster(data=as.data.frame(ODDy),aes(Longitude,Latitude,fill=Disp),
#                        alpha=0.5,interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
#       scale_fill_gradient2(low = "blue",mid="blue",high = "red",trans = "log",
#                            # breaks=c(0,1,10,100),
#                            na.value = "transparent")+#,limits=c(0,100)) +
#       # labs(fill = "Displaced Pop. / Total Pop.")+xlab("Longitude") + ylab("Latitude");p
#       labs(fill = "Displaced Population")+xlab("Longitude") + ylab("Latitude");p1
# 
#     p2<-q+ geom_raster(data = as.data.frame(ODDy),
#                                mapping = aes(Longitude,Latitude,fill=hazMean1),
#                                alpha=0.5,  inherit.aes = FALSE) + coord_cartesian() +
#       scale_fill_gradient2(low = "darkorchid2",mid="olivedrab4",high = "yellow",
#                            midpoint = 6,breaks=c(4.5,5,5.5,6,6.5,7,7.5,8),
#                            na.value = "transparent")+#,limits=c(0,100)) +
#       labs(fill = "Hazard Intensity")+xlab("Longitude") + ylab("Latitude");p2    
#     # p2<-q+ geom_contour_filled(data = as.data.frame(ODDy),
#     #                   mapping = aes(Longitude,Latitude,z=hazMean1,fill=..level..),
#     #                   alpha=0.5, inherit.aes = FALSE,bins = 20,) + coord_cartesian() +
#     #   scale_fill_discrete(low = "transparent",high = "purple",na.value = "transparent") + 
#     #   labs(fill = "Hazard Intensity")+xlab("Longitude") + ylab("Latitude");p2
#     
#     plotty<-gridExtra::grid.arrange(p1,p2,nrow=1); plotty
#     
#     if(is.null(input)) return(plotty)
#         
#     ggsave(paste0(input$datadir,input$plotdir,ODDy@hazard,input$sdate,input$iso3,".png"), plotty)
#     
#     return(plotty)
#   
# }

convODDy2raster<-function(ODDy,filename=NULL){
  
  ODDy$ISO3C<-NULL
  class(ODDy)<-"SpatialPixelsDataFrame"
  ODDy%<>%brick()
  
  if(!is.null(filename)) {
    x <- terra::rast(raster(ODDy))
    terra::writeRaster(x,filename=filename,overwrite=TRUE)
  }
  return(ODDy)
  
}

OutputODDyFile<-function(ODDy,filer=NULL){
  
  # if(sum(grepl("FBPop",colnames(ODDy@data),fixed = T))>0) {
    # ODDy$Population<-NULL
  #   
  # }
  ODDy$ISO3C<-NULL
  
  ODDy@modifier$HTI<- -0.16
  ODDy<-DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)
  ODDy@predictDisp
  ODDy$LowCIDisp<-ODDy$Disp
  ODDy@modifier$HTI<- 0.4
  ODDy<-DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)
  ODDy@predictDisp
  ODDy$HiCIDisp<-ODDy$Disp
  ODDy@modifier$HTI<- 0
  ODDy<-DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)
  ODDy@predictDisp
  ODDy$ExpDisp<-ODDy$Disp
  ODDy$Disp<-NULL
  
  ODDy<-convODDy2raster(ODDy,filename=filer)
  
  return(ODDy)
  
}




# setMethod("DispX", "ODD", function(ODD,Omega,center,saver=F,
#                                    Method=list(Np=20,cores=4,
#                                                clNeeds=c("SplitSamplePop","ExtractCIndy","fDamUnscaled","BinR","rbinom"),
#                                                packages=c("dplyr","magrittr"))
# ){
#   
#   # Extract 0D parameters
#   Params<-FormParams(ODD,list(Np=Method$Np,center=center))
#   # Income distribution percentiles
#   SincN<-seq(10,90,by = 10)
#   
#   Disp<-c()
#   notnans<-(1:nrow(ODD))[!(is.na(ODD$Population) | is.na(ODD$ISO3C) | is.na(ODD$GDP))]
#   
#   CalcDisp<-function(ij,ODD,Omega,center,saver,Method,Params,SincN){
#     iso3c<-ODD$ISO3C[ij]
#     tParams<-Params[[iso3c]]
#     # Sample population per income distribution (Assumes 9 percentiles)
#     lPopS<-SplitSamplePop(Pop=ODD$Population[ij],Method$Np) 
#     tDisp<-array(0,Method$Np)
#     # for(h in 1:sum(grepl("hazMean",names(ODD)))){
#     for(h in c(1)){
#       if(is.na(ODD@data[ij,paste0("hazMean",h)])) next
#       # Resample population based on who is already displaced
#       ind<-(colSums(lPopS)-tDisp)>0
#       if(h>1) {
#         if(ncol(lPopS[,ind])==0) break
#         if(sum(ind>1)) sumz<-colSums(lPopS[,ind])
#         else sumz<-sum(lPopS[,ind])
#         lPopS[,!ind]<-0
#         lPopS[,ind]<-SplitSamplePop(Pop=(sumz-tDisp[ind]))
#       }
#       # Sample hazard Intensity 
#       # the uncertainty is too high... so I scale it to get some interpretable results (I know, I'm not really a statistician, I don't even have a degree, I was cleaning when they confused me for the interviewee. I didn't have the heart to say anything. You don't hate me as much as I do)
#       # I_ij<-rnorm(n = Method$Np,
#       #             mean = ODD@data[ij,paste0("hazMean",h)],
#       #             sd = ODD@data[ij,paste0("hazSD",h)]/10) 
#       I_ij<-ODD@data[ij,paste0("hazMean",h)]
#       
#       # Separate into income distributions
#       for (s in 1:length(SincN)){
#         if(all(lPopS[s,]==0)) next
#         # Extract income percentile  
#         Sinc<-ExtractCIndy(ODD,iso = iso3c,var = paste0("p",SincN[s],"p100"))
#         # Multiply-in the GDP factor to Sinc
#         dollar<-pull(dplyr::select(Sinc,-iso3),value); dollar%<>%multiply_by(ODD$GDP[ij])
#         tParams$var$value[tParams$var$variable=="dollar"]<-dollar
#         # Predict damage at coordinate {i,j} (vector with MC particles)
#         Damage<-fDamUnscaled(I_ij,tParams,Omega)
#         if(any(is.na(Damage))) {print("WARNING: NA produced by damage calculation")}
#         # Scaled damage
#         Dprime<-BinR(Omega$Lambda$nu*Damage + Omega$Lambda$omega,Omega$zeta)
#         # Local displacement additions
#         tDisp[ind]%<>%add(mapply(function(size,p) 
#           rbinom(n = 1,size,p),lPopS[s,ind],Dprime[ind]))
#       }  
#     }
#     # Add to total displaced for this country
#     return(tDisp)
#   }
#   
#   cl<-makeCluster(Method$cores,outfile="./Rcode/OUTPUT.txt")
#   registerDoParallel(cl)
#   clusterExport(cl, Method$clNeeds)
#   
#   Disp <- foreach(ij%in%notnans, .packages = Method$packages, .combine = c) %dopar%{
#     to.Disp <- tryCatch({CalcDisp(ij,ODD,Omega,center,saver,Method,Params,SincN)}, error=function(e) NA)
#   }
#   
#   return(Disp)
#   
#   # Find the maximum number of displaced persons per country
#   funcy<-function(i) data.frame(iso3=ODD$ISO3C,IDPs=Disp[i,])%>%
#     group_by(iso3)%>%summarise(gmax=sum(IDPs,na.rm = T),.groups = 'drop_last')%>%
#     filter(iso3==ODD@gmax$iso3)%>%pull(gmax)
#   predictDisp<-vapply(1:Method$Np,funcy,numeric(1))
#   
#   # If the IDMC estimate is foreseen to be a lower or upper bound, or a generally poor estimate
#   predictDisp%<>%qualifierDisp(qualifier = ODD@gmax$qualifier,mu = Omega$mu)
#   
#   # Ensure that values are scaled
#   tpredictDisp<-apply(cbind(predictDisp,rep(sum(ODD$Population,na.rm = T),Method$Np)),1,min)
#   if(!saver) return(tpredictDisp) 
#   
#   # Scale the spatial values
#   Disp%<>%multiply_by(tpredictDisp/predictDisp)%>%round()
#   # Find the best fit solution
#   MLE<-which.min((rowSums(Disp)-ODD@gmax$gmax)^2)
#   # Save into ODD object
#   ODD$Displacements<-Disp[MLE,]
#   # I know this is kind of repeating code, but I want the other function as fast as possible
#   ODD@predictDisp<-data.frame(iso3=ODD$ISO3C,IDPs=Disp[MLE,])%>%
#     group_by(iso3)%>%summarise(gmax=sum(IDPs,na.rm = T),.groups = 'drop_last')
#   
# })

# Convert from ISO3C to country ("GBR"="United Kingdom")
# transIso(ODD)
# do.call(selectMethod("transIso","ODD"),list(ODD=ODD))
setGeneric("transIso", function(ODD) 
  standardGeneric("transIso") )
setMethod("transIso", "ODD", function(ODD)
  return(countrycode::countrycode(sourcevar = as.character(unique(ODD@cIndies$iso3)),
                                  origin = "iso3c",
                                  destination = "country.name")))
