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
#   - Displacement date data
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

source('RCode/Functions.R')
source('RCode/Model.R')
source('RCode/GetPopDemo.R')
source('RCode/GetSocioEconomic.R')
source('RCode/AddVulnerability.R')
library(devtools)
library(parallel)
library(doParallel)
library(foreach)
library(parallelsugar)

checkODD<-function(object) {
  
  if(!is.null(object$Population) & any(object$Population<0,na.rm = T)) return(F) 
  if(!is.null(object$GDP) & any(object$GDP<0,na.rm = T)) return(F) 
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
                   polygons="list"),
         contains = "SpatialPixelsDataFrame")

ExtractI0poly<-function(hsdf,ODD){

  # Extract contours
  pcontour<-adehabitatMA::getcontour(raster::subset(hsdf,mean>=ODD@I0))
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

# Add GDP data to the ODD object by interpolating onto the grid using cubic splines
setGeneric("AddHazSDF", function(ODD,lhazSDF) 
  standardGeneric("AddHazSDF") )
setMethod("AddHazSDF", "ODD", function(ODD,lhazSDF){
  
  ODD@I0<-lhazSDF$hazard_info$I0
  # interpolate data onto the grid
  coords<-Genx0y0(ODD)
  lenny<-length(lhazSDF) ; start<-2
  alertscores<-alertlevels<-c() ; dates<-rep(lhazSDF$hazard_info$sdate,lenny-start+1)
  
  polysave<-array(F,dim=c(nrow(ODD),(lenny-start+1)))
  
  for (i in start:lenny){
    print(i-start+1)
    hsdf<-lhazSDF[[i]]
    # Extract detail of specific hazard
    dates[i-start+1]<-hsdf@eventdate
    alertlevels%<>%c(hsdf@alertlevel)
    alertscores%<>%c(hsdf@alertscore)
    
    if(lhazSDF$hazard_info$hazard=="TC"){
      
      layer<-with(as.data.frame(hsdf),akima::interp(x=Longitude,y=Latitude,z=mean,
                                     xo=coords$xo,yo=coords$yo,
                                     linear=F,extrap = F))

      
      layer<-c(layer$z)
      layer[layer<ODD@I0]<-NA
      
      ODD@data[paste0("hazMean",i-start+1)]<-layer
      
    } else {
      
      # extract polycontour of I<I0
      pcontour<-ExtractI0poly(hsdf=hsdf,ODD=ODD)
      # Find all ODD coordinates not inside polycontour
      insidepoly<-rep(F,nrow(ODD))
      if (length(unique(pcontour$id)) > 0){
        for(p in 1:length(unique(pcontour$id))){
          tcont<-filter(pcontour,id==p)
          insidepoly<-insidepoly | sp::point.in.polygon(ODD@coords[,1],
                                                        ODD@coords[,2],
                                                        tcont$Longitude,
                                                        tcont$Latitude)>0
        }
      }
      rm(tcont)
      hsdf%<>%as.data.frame
      # Interpolate BOTH MEAN & SD onto the ODD grid
      print("mean")
      layer<-with(hsdf,akima::interp(x=Longitude,y=Latitude,z=mean,
                                     xo=coords$xo,yo=coords$yo,
                                     linear=F,extrap = F))
      layer<-c(layer$z)
      layer[!insidepoly]<-NA
      if(all(is.na(layer))) next
      
      ODD@data[paste0("hazMean",i-start+1)]<-layer
      
      print("sd")
      layer<-with(hsdf,akima::interp(x=Longitude,y=Latitude,z=sd,
                                     xo=coords$xo,yo=coords$yo,
                                     linear=F,extrap = F))
      layer<-c(layer$z)
      layer[!insidepoly]<-NA
      ODD@data[paste0("hazSD",i-start+1)]<-layer
      
      polysave[,i-start+1]<-insidepoly
    }
  }
  
  if(lhazSDF$hazard_info$hazard=="TC") { 
    ind<-unname(apply(ODD@data,2, function(x) sum(!is.na(x))))
    ODD@data<-ODD@data[,ind>0]
  }else ODD$Population[rowSums(polysave)==0]<-NA
  
  # ODD@hazdates<-dates
  ODD@alerts<-data.frame(alertscores=alertscores,alertlevels=alertlevels)
  
  return(ODD)
  
})

# Initialisation of the ODD object
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
            )) {
            
            if(is.null(lhazSDF)) return(.Object)
            if(!class(lhazSDF[[length(lhazSDF)]])[1]=="HAZARD") return(.Object)
            
            
            
            
 
            
            #stop("Also remove INFORM crap from Model.R")
            
            

            
            
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

            year<-AsYear(dater)
            
            print("Fetching population data")
            obj <-GetPopulationBbox(.Object@dir,bbox=bbox)
            .Object@data <- obj@data
            .Object@coords.nrs <-obj@coords.nrs
            .Object@grid <-obj@grid
            .Object@grid.index <-obj@grid.index
            .Object@coords <-obj@coords
            .Object@bbox <-obj@bbox
            .Object@proj4string <-crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            print("Adding hazard events")
            # Including minshake polygon per hazard event using getcontour from adehabitatMA package
            .Object%<>%AddHazSDF(lhazSDF)
            
            # Extract empty indices to save time
            inds<-which(!is.na(.Object$Population))
            .Object@data$ISO3C[-inds] <- NA
            
            print("Filter spatial data per country")
            #.Object@data$ISO3C<-NA_character_
            #.Object@data$ISO3C[inds]<-coords2country(.Object@coords[inds,])
            iso3c<-unique(.Object@data$ISO3C) ; iso3c<-iso3c[!is.na(iso3c)]
            

            print("Interpolate population values")
            # Note there are as many values returned as iso3c codes (returns as data.frame with columns 'iso3' and 'factor')
            Popfactors<-InterpPopWB(iso3c,dater)
            
            for (iso in iso3c){
              indie<-.Object@data$ISO3C==iso & !is.na(.Object@data$ISO3C)
              .Object@data$Population[indie]%<>%
                multiply_by(Popfactors$factor[Popfactors$iso3==iso])
            }
            
            # World Income Database (WID) data:
            if(year==AsYear(Sys.Date())) year<-AsYear(Sys.Date())-1
            print("Extract country indicators - WID:")
            
            WID<-GetWID_perc(Model$WID_perc,iso3c,year)
            
            
            #stop("Add the full variables to the cIndies data.frame")
            # Bind it all together!
            .Object@cIndies<-WID
            
            # Here we add the vulnerabilities used in the linear predictor
            .Object%<>%AddVuln()
            
            
            #print("Fetching GDP-PPP data")
            #.Object%<>%AddGDP(inds)
            
            
            # # print("Fetching building count data")
            # .Object%<>%AddBuildingCounts()
            
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

setGeneric("DispX", function(ODD,Omega,center, BD_params, LL, sim=F, Method)
  standardGeneric("DispX") )
# Code that calculates/predicts the total human displacement 
setMethod("DispX", "ODD", function(ODD,Omega,center, BD_params, LL=F, sim=F, 
                                   Method=list(Np=20,cores=8,cap=-300, 
                                               kernel_sd=list(displacement=0.15,mortality=0.03,buildDam=0.15,buildDest=0.1), kernel='lognormal')
){
  # ... Function description ...
  # LL: Returns 'likelihood' if true or data simulated from model if false
  # Sim: Set to true when generating data for a simulated ODD object. 
  
  # Extract 0D parameters & speed up loop
  Params<-FormParams(ODD,list(Np=Method$Np,center=center))
  # Income distribution percentiles & extract income percentile  
  SincN<-paste0('p',seq(10,80,10), 'p', seq(20,90,10))
  Sinc<-ExtractCIndy(ODD,var = SincN)
  # Speed-up calculation (through accurate cpu-work distribution) to only values that are not NA
  if(!LL) {notnans<-which(!(is.na(ODD$Population) | is.na(ODD$ISO3C) | is.na(ODD$ExpSchYrs)))
  } else notnans<-which(!(is.na(ODD$Population) | is.na(ODD$ISO3C) ))#!ODD$ISO3C%in%ODD@impact$iso3))
  
  nBuildings_sample <- array(0, dim=c(NROW(ODD@data), Method$Np))
  if (!is.null(ODD@data$nBuildings)){
    indies_openbuildings <- which(!is.na(ODD@data$nBuildings))
    nBuildings_sample[indies_openbuildings,] <- matrix(ODD@data$nBuildings[indies_openbuildings] , length(ODD@data$nBuildings[indies_openbuildings]), Method$Np)
    BD_data_present <- T
  } else if (!is.null(ODD@data$nBuiltup)){
    indies_osm <- which(ODD@data$nBuiltup > 0)
    nBuildings_sample[indies_osm,] <- sampleBuildingsFromBuiltup(ODD$nBuiltup[indies_osm], ODD$ISO3C[indies_osm], reg_coefs, Method$Np)
    BD_data_present <- T
  } else {
    BD_data_present <- F
  }
  
  # Calculate non-local linear predictor values
  LP<-GetLP(ODD,Omega,Params,Sinc,notnans)
  LP_buildings <- GetLP(ODD,Omega,Params,Sinc,notnans, split_GNI=F) #could be sped up by combining
  
  # Speed things up a little
  hrange<-grep("hazMean",names(ODD),value = T)
  
  # Function to predict displacement per gridpoint
  CalcDam<-function(ij){
    # Calculate local linear predictor (NOTE: is a vector due to income distribution)
    locallinp<- LP[ij,] # LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]]*LP$Plinp[ij]*LP$linp[[iso3c]] 
    #locallinp<-rep(1,10) #reduces parameter space and removes demographic covariates
    locallinp_buildings <- LP_buildings[ij]
    # Sample population per income distribution (Assumes 9 percentiles):
    lPopS <- SplitSamplePop(Pop=ODD@data$Population[ij],Method$Np) 
    tPop <-array(0,c(3, Method$Np)) #row 1 = tDisp, #row 2 = tMort, #row 3 = tRem
    tPop[3,]=colSums(lPopS)
    for(h in hrange){
      # for(h in c(1)){
      if(is.na(ODD@data[ij,h])) next
      # Resample population based on who is remaining
      ind<-tPop[3,]>0
      if(h!=hrange[1]) {
        if(sum(ind)==0) break #if no remaining population, skip modelling
        if(length(lPopS[,ind])==0) break #if no remaining population, skip modelling
        #if(sum(ind)>1) sumz<-colSums(lPopS[,ind])
        #else sumz<-sum(lPopS[,ind])
        #lPopS[,!ind]<-0
        lPopS[,ind]<-SplitSamplePop(Pop=tPop[3,ind])
      }
      # Sample hazard Intensity 
      # the uncertainty is too high... so I scale it to get some interpretable results (I know, I'm not really a statistician, I don't even have a degree, I was actually just having a look around the department when they confused me for the interviewee. I didn't have the heart to say anything. You don't hate me as much as I do)
      # I_ij<-rnorm(n = Method$Np,
      #             mean = ODD@data[ij,paste0("hazMean",h)],
      #             sd = ODD@data[ij,paste0("hazSD",h)]/10)
      I_ij<-ODD@data[ij,h]
      
      # Separate into income distributions (as each have 10% of population, order doesn't matter)
      for (s in 1:length(SincN)){
        if(all(lPopS[s,]==0)) next
        # Predict damage at coordinate {i,j} (vector with MC particles)
        Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Params$Np),Omega)*locallinp[s], error=function(e) NA)
        if(any(is.na(Damage))) print(ij)
        
        #LOOSEEND: Include [ind] here 
        D_MortDisp <- D_MortDisp_calc(Damage, Omega) #First row of D_MortDisp is D_Mort, second row is D_Disp
        D_Rem <- pmax(0, 1 - D_MortDisp[2,] - D_MortDisp[1,]) #probability of neither death nor displacement. Use pmax to avoid errors caused by numerical accuracy.
        
        # Accumulate the number of people displaced/deceased, but don't accumulate the remaining population
        tPop[3,ind]<-0
        tPop[,ind]<-tPop[,ind] + Fbdam(lPopS[s,ind],D_MortDisp[2,ind], D_MortDisp[1,ind], D_Rem[ind]) 
      } 
    }
    #ensure the total displaced, deceased or remaining does not exceed total population
    tPop[tPop>ODD@data$Population[ij]]<-floor(ODD@data$Population[ij])
    
    #if no building destruction data:
    if(!BD_data_present) return(rbind(tPop[1:2,, drop=FALSE], rep(NA, Method$Np))) #return total displacement and mortality, set number of buildings destroyed to NA
    
    #otherwise, sample from the model for the number of buildings destroyed:
    #we take locallinp[5] which corresponds to locallinp for the median GDP
    
    nBuild <-array(0,c(3, Method$Np)) #row 1 = nDam, #row 2 = nDest, #row 3 = nRem
    nBuild[3,]= nBuildings_sample[ij,]
    first_haz = T
    for (h in hrange){
      if(is.na(ODD@data[ij,h])) next
      if(h!=hrange[1]) {
        first_haz = F
        if(all(nBuild[3,]==0)) break #if no remaining buildings, skip modelling LOOSEEND
      }

      I_ij<-ODD@data[ij,h]
      Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Params$Np),Omega)*locallinp_buildings, error=function(e) NA) #calculate unscaled damage (excluding GDP)
      
      D_DestDam <- D_DestDam_calc(Damage, Omega, first_haz, Model$DestDam_modifiers) #First row of D_DestDam is D_Dest, second row is D_Dam
      D_Rem <- pmax(0, 1 - D_DestDam[1,] - D_DestDam[2,]) #probability of neither damaged nor destroyed. Use pmax to avoid errors caused by numerical accuracy.
      
      if (h != hrange[1]){
        n_DamtoDest <- fBD(nBuild[2,], D_DestDam[1,] ^ (Model$DestDam_modifiers[3]/Model$DestDam_modifiers[2]))
        nBuild[1,] <- nBuild[1, ] + n_DamtoDest
        nBuild[2,] <- nBuild[2, ] - n_DamtoDest
      }
      
      # Accumulate the number of buildings damaged/destroyed, but not the number of buildings remaining
      nRem <- nBuild[3,]
      nBuild[3,] <- 0 
      nBuild <- nBuild + Fbdam(nRem, D_DestDam[2,], D_DestDam[1,], D_Rem) 
    }
    return(rbind(tPop[1:2,,drop=FALSE], nBuild[1:2,,drop=FALSE]))
  }

  Dam<-array(0,c(nrow(ODD),Method$Np,4)) # Dam[,,1] = Displacement, Dam[,,2] = Mortality, Dam[,,3] = Buildings Damaged, Dam[,,4] = Buildings Destroyed
  
  #Method$cores is equal to AlgoParams$NestedCores (changed in Model file)
  if(Method$cores>1) { Dam[notnans,,]<-aperm(simplify2array(mclapply(X = notnans,FUN = CalcDam,mc.cores = Method$cores)), perm=c(3,2,1))
  } else  Dam[notnans,,]<- aperm(simplify2array(lapply(X = notnans,FUN = CalcDam)), perm=c(3,2,1))

  # return(Disp)

  # If the IDMC estimate is foreseen to be a lower or upper bound, or a generally poor estimate
  # for(c in ODD@gmax$iso3){
  #   ind<-ODD@data$ISO3C==c  & !is.na(ODD@data$ISO3C)
  #   Disp[ind,]%<>%qualifierDisp(qualifier = ODD@gmax$qualifier[ODD@gmax$iso3==c],mu = Omega$mu)
  # }
  
  
  funcy<-function(i,LLout=T, kernel_sd=Method$kernel_sd, kernel=Method$kernel, cap=Method$cap,polygons_indexes=ODD@polygons) {
    
    tmp<-data.frame(iso3=ODD$ISO3C, displacement=Dam[,i,1], mortality=Dam[,i,2], buildDam=Dam[,i,3], buildDest=Dam[,i,4])
    impact_sampled<-data.frame(polygon = numeric(), impact = character(), sampled = numeric())
    
    for (polygon_id in unique(ODD@impact$polygon)){
      polygon_impacts <- ODD@impact$impact[which(ODD@impact$polygon==polygon_id)]
      for (impact in polygon_impacts){
        impact_sampled %<>% rbind(data.frame(polygon=polygon_id, 
                                             impact=impact,
                                             sampled= floor(sum(tmp[polygons_indexes[[polygon_id]]$indexes,impact], na.rm=T))))
      }
    }
    impact_obs_sampled <- merge(impact_sampled, ODD@impact, by=c("polygon", "impact")) %>% arrange(desc(observed)) 
    if(LLout) return(LL_IDP(impact_obs_sampled, kernel_sd,  kernel, cap))
    return(impact_obs_sampled)
    
  }
  
  if (sim == T){ 
    ODD@data$Disp<-Dam[,1,1]  #should this be named data$DispPred or something instead?
    ODD@data$Mort<-Dam[,1,2]
    ODD@data$BuildDam<-Dam[,1,3]
    ODD@data$BuildDest<-Dam[,1,4]
    ODD@predictDisp<-funcy(1,LLout=F) 
    return(ODD)
  }
  
  if(LL == F){
    return(lapply(1:Method$Np, funcy, LLout=F))
  }
  
  outer<-vapply(1:Method$Np,funcy,numeric(1), kernel_sd=Method$kernel_sd)
  outer[outer < (-745)]<- -745
  
  
  # Find the best fit solution, not sure what this is for ...
  # if(NROW(outer)>1) {
  #   if(LL)  return(log(rowMeans(exp(outer),na.rm=T)))
  #   MLE<-which.max(log(colSums(exp(outer),na.rm=T)))
  # }  else {
  #   return(outer) #SMC-CHANGE
  #   if(LL)  return(log(mean(exp(outer),na.rm=T)))
  #   MLE<-which.max(log(exp(outer)))
  # }
  
  if(LL) return(outer)
  
  if(Method$Np == 1){
    MLE=1
  } else {
    MLE <- which.max(outer)
  }
  # Save into ODD object
  # ODD@data$Disp<-Disp[,MLE]*sum(ODD@gmax$gmax)/mean(sum(Disp[,MLE])) %>% round()
  ODD@data$Disp<-Dam[,MLE,1]  #should this be named data$DispPred or something instead?
  ODD@data$Mort<-Dam[,MLE,2]
  ODD@data$BuildDam<-Dam[,MLE,3]
  ODD@data$BuildDest<-Dam[,MLE,4]
  # I know this is kind of repeating code, but I want the other function as fast as possible
  ODD@predictDisp<-funcy(MLE,LLout=F) 
  
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
               "DamBDFB" = "Damaged Buildings (HRSL)")
  
  as.character(lookup[varname])
  
}

plotODDy<-function(ODDy,zoomy=7,var="Population",breakings=NULL,bbox=NULL,alpha=0.5,map="terrain"){
  
  if(is.null(breakings) & (var=="Population" | var=="Disp")) breakings<-c(0,1,5,10,50,100,500,1000)
  
  if(is.null(bbox)) bbox<-ODDy@bbox
  
  mad_map <- get_stamenmap(bbox,source = "stamen",maptype = map,zoom=zoomy)
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
                        alpha=1.0,breaks = brks) +
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

MakeODDPlots<-function(ODDy, input=NULL,zoomer=7, bbox=NULL){
  
    ODDy@data$Disp[ODDy@data$Disp<1]<-0
  
    mad_map <- get_stamenmap(ODDy@bbox,source = "stamen",maptype = "terrain",zoom=zoomer)
    q<-ggmap(mad_map)
    p1<-q+ geom_raster(data=as.data.frame(ODDy),aes(Longitude,Latitude,fill=Disp),
                       alpha=0.5,interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
      scale_fill_gradient2(low = "blue",mid="blue",high = "red",trans = "log",
                           # breaks=c(0,1,10,100),
                           na.value = "transparent")+#,limits=c(0,100)) +
      # labs(fill = "Displaced Pop. / Total Pop.")+xlab("Longitude") + ylab("Latitude");p
      labs(fill = "Displaced Population")+xlab("Longitude") + ylab("Latitude");p1

    p2<-q+ geom_raster(data = as.data.frame(ODDy),
                               mapping = aes(Longitude,Latitude,fill=hazMean1),
                               alpha=0.5,  inherit.aes = FALSE) + coord_cartesian() +
      scale_fill_gradient2(low = "darkorchid2",mid="olivedrab4",high = "yellow",
                           midpoint = 6,breaks=c(4.5,5,5.5,6,6.5,7,7.5,8),
                           na.value = "transparent")+#,limits=c(0,100)) +
      labs(fill = "Hazard Intensity")+xlab("Longitude") + ylab("Latitude");p2    
    # p2<-q+ geom_contour_filled(data = as.data.frame(ODDy),
    #                   mapping = aes(Longitude,Latitude,z=hazMean1,fill=..level..),
    #                   alpha=0.5, inherit.aes = FALSE,bins = 20,) + coord_cartesian() +
    #   scale_fill_discrete(low = "transparent",high = "purple",na.value = "transparent") + 
    #   labs(fill = "Hazard Intensity")+xlab("Longitude") + ylab("Latitude");p2
    
    plotty<-gridExtra::grid.arrange(p1,p2,nrow=1); plotty
    
    if(is.null(input)) return(plotty)
        
    ggsave(paste0(input$datadir,input$plotdir,ODDy@hazard,input$sdate,input$iso3,".png"), plotty)
    
    return(plotty)
  
}

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
