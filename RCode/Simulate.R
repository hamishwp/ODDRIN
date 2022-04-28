library(gstat)

setClass("ODDSim", contains="ODD")
setClass("BDSim", contains="BD")
# Get ODDSim objects to inherit the methods and attributes of ODD objects

# Make a new method for the initialisation of the ODDSim object. Much is the same as for an ODD object,
# but some steps involve simulation rather than retrieving the necessary data. 
setMethod(f="initialize", signature="ODDSim",
          definition=function(.Object,lhazSDF=NULL,DamageData=NULL, PopSim=NULL, GDPSim = NULL, dir="./",Model=list(
            INFORM_vars=c("CC.INS.GOV.GE", # Government Effectiveness
                          "VU.SEV.AD", # Economic Dependency Vulnerability
                          "CC.INS.DRR", # Disaster Risk Reduction
                          "VU.SEV.PD", # Multi-dimensional Poverty
                          "CC.INF.PHY" # Physical Infrastructure
            ), 
            fIndies=list(CC.INS.GOV.GE=returnX, # Government Effectiveness
                         VU.SEV.AD=returnX, # Economic Dependency Vulnerability
                         CC.INS.DRR=returnX, # Disaster Risk Reduction
                         VU.SEV.PD=returnX, # Multi-dimensional Poverty
                         CC.INF.PHY=returnX, # Physical Infrastructure
                         dollar=returnX, # IncomeDistribution*GDP
                         Pdens=returnX), # IncomeDistribution*GDP
            WID_perc=   c("p10p100", # top 90% share of Income Distribution
                          "p20p100", # top 80% share of Income Distribution
                          "p30p100", # top 70% share of Income Distribution
                          "p40p100", # top 60% share of Income Distribution
                          "p50p100", # top 50% share of Income Distribution
                          "p60p100", # top 40% share of Income Distribution
                          "p70p100", # top 30% share of Income Distribution
                          "p80p100", # top 20% share of Income Distribution
                          "p90p100" # top 10% share of Income Distribution
            ))) {
            
            if(is.null(lhazSDF)) return(.Object)
            if(!class(lhazSDF[[length(lhazSDF)]])[1]=="HAZARD") return(.Object)
            
            if(lhazSDF$hazard_info$hazard=="EQ") Model$INFORM_vars%<>%c("HA.NAT.EQ")
            else if(lhazSDF$hazard_info$hazard=="TC") Model$INFORM_vars%<>%c("HA.NAT.TC")
            else if(lhazSDF$hazard_info$hazard=="FL") Model$INFORM_vars%<>%c("HA.NAT.FL")
            else stop("Not currently prepared for hazards other than EQ, TC or FL")
            
            .Object@dir<-dir
            .Object@hazard<-lhazSDF$hazard_info$hazard
            
            if(length(unique(DamageData$eventid))==1) .Object@eventid<-unique(DamageData$eventid)
            if(.Object@hazard%in%c("EQ","TC")){
              .Object@gmax<-DamageData%>%group_by(iso3)%>%
                summarise(gmax=max(gmax),
                          qualifier=ifelse(all(is.na(gmax)), NA_character_, qualifierDisp[which.max(gmax)]),
                          mortality=max(mortality),
                          qualifierMort=ifelse(all(is.na(mortality), NA_character_, qualifierMort[which.max(mortality)])),
                          buildDestroyed=max(buildDestroyed), 
                          qualifierBD = ifelse(all(is.na(buildDestroyed), NA_character_, qualifierMort[which.max(mortality)])))
              .Object@IDPs<-DamageData[,c("sdate","gmax","qualifierDisp")]%>%
                transmute(date=sdate,IDPs=gmax,qualifier=qualifierDisp)
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
            
            year<-AsYear(dater)
            
            obj<-PopSim
            .Object@data <- obj@data
            .Object@coords.nrs <-obj@coords.nrs
            .Object@grid <-obj@grid
            .Object@grid.index <-obj@grid.index
            .Object@coords <-obj@coords
            .Object@bbox <-obj@bbox
            .Object@proj4string <-crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            
            .Object@data$nBuildings <- round(runif(1,0.2,1) * PopSim$Population + rnorm(length(PopSim$Population),0,20))
            .Object@data$nBuildings[.Object@data$nBuildings < 0] <- 0
            
            # Including minshake polygon per hazard event using getcontour from adehabitatMA package
            # LOOSEEND
            .Object%<>%AddHazSDF(lhazSDF)
            
            # Extract empty indices to save time
            inds<-!is.na(.Object$Population)
            
            .Object$GDP <- GDPSim
            
            .Object@data$ISO3C<-NA_character_
            .Object@data$ISO3C[inds]<-'ABC'
            iso3c<-unique(.Object@data$ISO3C) ; iso3c<-iso3c[!is.na(iso3c)]
            
            # Skip the interpolation of population and GDP values for the simulations
            
            # Simulate INFORM (Joint Research Center - JRC) data:
            # Sample the indicators from a uniform distribution on [0,10]
            # Set some to 0 each with a probability of 0.1 (to mimic presence of zero values in real data).
            INFORM <- data.frame(iso3='ABC', 
                                 value=runif(length(Model$INFORM_vars),0,10),
                                 variable=Model$INFORM_vars)
            INFORM_missing <- rbernoulli(nrow(INFORM), 0.1)
            INFORM$value[INFORM_missing]<-0
            
            #Simulate World Income Database (WID) data
            #Simulate a cubic distribution with max given by max_inc (WID data appears close to cubic)
            max_inc <- runif(1,0.4,0.8)
            perc <- seq(0.1,0.9,0.1)
            cubic <- perc^3+runif(1,0.3,1)*perc^2
            WID <- data.frame(
              variable=Model$WID_perc,
              iso3='ABC',
              value=cubic*max_inc/max(cubic)
            )
            
            # Bind it all together!
            .Object@cIndies<-rbind(INFORM,WID)
            .Object@fIndies<-Model$fIndies
            
            linp<-rep(list(1.),length(unique(.Object@cIndies$iso3)))
            names(linp)<-unique(.Object@cIndies$iso3)
            .Object@modifier<-linp
            
            print("Checking ODDSim values")
            checkODD(.Object)
            
            return(.Object)
          }
)

setMethod(f="initialize", signature="BDSim",
          definition=function(.Object,ODD=NULL) {
            
            if(is.null(ODD)) return(.Object)            
            
            .Object@hazard<-ODD@hazard
            .Object@cIndies<-ODD@cIndies[ODD@cIndies$iso3%in%'ABC',]
            .Object@I0<-ODD@I0
            .Object@hazdates<-ODD@hazdates
            .Object@eventid<-ODD@eventid
            .Object@fIndies<-ODD@fIndies
            
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
              for (k in 1:ODD$nBuildings[ij]){
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
                                           proj4string = crs("+proj=longlat +datum=WGS84 +ellps=WGS84"))
            
            .Object@data <- Damage@data
            .Object@coords.nrs <-Damage@coords.nrs
            .Object@coords <-Damage@coords
            .Object@bbox <-Damage@bbox
            .Object@proj4string <-crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            rm(Damage)
            
            print("Interpolating population density, hazard & GDP-PPP data")
            .Object%<>%BDinterpODD(ODD=ODD)
            
            print("Filter spatial data per country")
            # We could just copy Damage$iso3 directly, but I don't believe in anyone or anything...
            if(NROW(.Object@data)>0) .Object@data$ISO3C<-'ABC'
            
            print("Accessing OSM to sample building height & area")
            # ExtractOSMbuildVol(.Object,ODD)
            
            linp<-rep(list(1.),length(unique(ODD@cIndies$iso3)))
            names(linp)<-unique(ODD@cIndies$iso3)
            .Object@modifier<-linp
            
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
  
  maxMag = runif(1, 4.5, 10)
  sigma = runif(1, 20, 38)
  r <- setValues(r, spatialEco::gaussian.kernel(sigma=sigma, n=r@nrows)) 
  r <- r * (maxMag/r@data@max)
  sd <- setValues(r, runif(r@ncols*r@nrows, 0.8,1.1))
  names(r) <- 'mmi_mean'
  r$mmi_sd <- sd
  hazsdf <- as(r, 'SpatialPixelsDataFrame')
  hazsdf<-hazsdf[hazsdf$mmi_mean>I0,]
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
  
  maxPopDens = runif(1, 100, 1000)
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
  PopSim[,3] = pweibull((PopSim[,3]-minval)/(maxval-minval),5,0.8) * maxPopDens
  PopSim %<>% acast(Longitude~Latitude, value.var='sim1') %>% convMat2SPDF(name='Population')
  return(PopSim)
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

simulateODDSim <- function(miniDam, I0=4.5){ 
  # Simulates hazard, population and GDP data over 120 x 120 grid cells (over the region [-0.5,0.5] x [-0.5,0.5])
  # INPUT: 
  #  - I0: minimum intensity resulting in damage 
  #  - miniDam: contains hazard info such as sdate, fdate, eventid, etc.
  # OUTPUT: 
  #  - ODDSim: an object of class ODDSim (the simulated equivalent of an ODD object)
  
  r <- raster(ncol=100, nrow=100, xmn=-0.125, xmx=0.125, ymn=-0.125,ymx=0.125, crs="+proj=longlat +datum=WGS84") #each cell is 30 arcseconds x 30 arcseconds
  
  lenny = rgeom(1, 0.8) + 1 #generate number of events according to a geometric distribution
  bbox = c(-0.125,-0.125,0.125,0.125) 
  lhazdat<-list(hazard_info=list(bbox=bbox,sdate=min(miniDam$sdate),fdate=max(miniDam$fdate),
                                 NumEvents=lenny,hazard="EQ", I0=I0, eventdates=rep(miniDam$sdate, lenny)))
  for(i in 1:lenny){
    hazsdf <- simulateEvent(r, I0)
    if(is.null(hazsdf@data$mmi_mean)){
      next
    }
    lhazdat <- append(lhazdat, new("HAZARD",
        obj=hazsdf,
        hazard="EQ",
        dater=min(miniDam$sdate),
        I0=I0,
        alertlevel=ifelse(max(hazsdf$mmi_mean)>7.5, 'red', ifelse(max(hazsdf$mmi_mean)>6, 'orange', 'green')), #assign pretty arbitrarily, don't think this is used in the model
        alertscore=0))
  }
  if (length(lhazdat)== 1){
    print('Simulation failed: affected region under simulated event is too small')
    return(NULL)
  }
  PopSim <- simulatePopulation(r)
  GDPSim <- simulateGDP(r)
  
  ODDSim <- new('ODDSim', lhazSDF=lhazdat,DamageData=miniDam, PopSim=PopSim, GDPSim=GDPSim)
  
  return(ODDSim)
}


simulateDataSet <- function(nEvents, Omega, Model, dir, I0=4.5, cap=-300){
  # Input:
  # - nEvents: The number of ODDSim objects to generate
  # - Omega: The model parameterisation
  # Output: None
  # Details:
  # - Saves each ODDSim object to IIDIPUS_SimInput/ODDobjects/
  set.seed(round(runif(1,0,100000)))
  
  haz='EQ'
  
  # Rather than simulating eventids and sdates, just take real events but remove all data except 
  # sdate, fdate, eventid, and displacement/mortality qualifiers
  DamageData = GetDisplacements(haz, saved=F, GIDD=F, EMDAT=T)
  DamageData %<>% head(n=nEvents) %>% mutate(iso3='ABC', gmax=NA, mortality=NA, buildDestroyed=NA, qualifierBD=NA) 
  
  ODDpaths = c()
  BDpaths = c()
  for(ev in unique(DamageData$eventid)){
    miniDam<-DamageData%>%filter(eventid==ev)
    ODDSim <- simulateODDSim(miniDam, I0=I0)
    if (is.null(ODDSim)){
      next
    }
    BDSim <- new('BDSim', ODDSim)
    namer<-paste0(ODDSim@hazard,
                  str_remove_all(as.character.Date(min(ODDSim@hazdates)),"-"),
                  unique(miniDam$iso3)[1],
                  "_",ODDSim@eventid)
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_SimInput/ODDobjects/",namer)
    saveRDS(ODDSim,ODDpath)
    
    BDpath<-paste0(dir,"IIDIPUS_SimInput/BDobjects/",namer) 
    saveRDS(BDSim,BDpath)
    print(paste0('Saved simulated hazard to ', namer))
  }
  
  #calculate Model$center based on the simulated events
  Model$center <- ExtractCentering(dir, input_folder='IIDIPUS_SimInput/', saver=F)
  
  ODDpaths <-na.omit(list.files(path="IIDIPUS_SimInput/ODDobjects/"))
  BDpaths <-na.omit(list.files(path="IIDIPUS_SimInput/BDobjects/"))
  k <- 10
  #now loop through each event and simulate the displacement, mortality, and building destruction using DispX()
  for(i in 1:length(ODDpaths)){
    ODDSim <- readRDS(paste0("IIDIPUS_SimInput/ODDobjects/",ODDpaths[i]))
    #simulate displacement, mortality and building destruction using DispX
    ODDSim %<>% DispX(Omega, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
    
    #take these simulations as the actual values
    ODDSim@gmax <- ODDSim@predictDisp %>% transmute(iso3='ABC',
                                                    gmax = round(rlnormTrunc(1,log(disp_predictor+k), sdlog=0.1, min=k)) - k,
                                                    qualifier = ifelse(is.na(gmax), NA, 'total'),
                                                    mortality = round(rlnormTrunc(1,log(mort_predictor+k), sdlog=0.03, min=k)) - k,
                                                    qualifierMort = ifelse(is.na(mort_predictor), NA, 'total'), 
                                                    buildDestroyed = round(rlnormTrunc(1,log(nBD_predictor+k), sdlog=0.1, min=k)) - k,
                                                    qualifierBD = ifelse(is.na(buildDestroyed), NA, 'total'))
                                        
    #overwrite ODDSim with the updated
    saveRDS(ODDSim, paste0("IIDIPUS_SimInput/ODDobjects/",ODDpaths[i]))
  }
  for (i in 1:length(BDpaths)){
    BDSim <- readRDS(paste0("IIDIPUS_SimInput/BDobjects/", BDpaths[i]))
    BDSim %<>% BDX(Omega, Model, LL=FALSE, Method=list(Np=1,cores=1))
    #take these simulations as the actual values
    BDSim@data$grading <- BDSim@data$ClassPred
    #switch those affected to 'Damaged' with probability 0.5
    affected <- which(BDSim@data$grading != 'notaffected')
    for (j in affected){
      BDSim@data$grading[j] = ifelse(runif(1)>0.5, BDSim@data$grading[j], 'Damaged')
    }
    saveRDS(BDSim,paste0("IIDIPUS_SimInput/BDobjects/",BDpaths[i]))
  }
  return(Model$center)
}

