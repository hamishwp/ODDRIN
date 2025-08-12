#######################################################################
####################### BDobj CLASS DEFINITION #######################
########### (Oxford-University Disaster Displacement Object) ##########
#######################################################################
# FIELDS:
#   - Helix id number
#   - Data - population, GDP, country ISO3C code, mean hazard intensity, stand.dev hazard intensity
#   - Hazard type
#   - Dates corresponding to hazard intensity (initialised with NA)
#   - Start date (could be evacuation initialisations)
#   - Per Iso3 extract and store indicators
#   - Building height GAMMA pars: gradient, intercept and sigma values glm(formula = y ~ 1, family = Gamma)
#   - Building area GAMMA pars: gradient, intercept and sigma values (>0,>0,>0)
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
# source('RCode/GetPopDemo.R')
# source('RCode/GetSocioEconomic.R')
# source('RCode/GetINFORM.R')

checkBD<-function(object) {
  stop("modify checks function")
  if(!str_length(object@iso3)==3&!is.na(object@iso3)) stop("ISO3 string length is not equal to 3 or NA")
  if(!is_empty(object@cIndies)|ncol(object@cIndies)!=ifelse(is.na(object@iso3),0,(length(object@iso3)+1))) stop("Country indicator data.frame column length not equal to number of countries (iso3's) or empty")
  
  TRUE
}

Genx0y0<-function(ODDobj){
  xo<-seq(ODDobj@bbox["Longitude","min"],
          ODDobj@bbox["Longitude","max"],
          length.out=ODDobj@grid@cells.dim[1])
  yo<-seq(ODDobj@bbox["Latitude","min"],
          ODDobj@bbox["Latitude","max"],
          length.out=ODDobj@grid@cells.dim[2])
  return(list(xo=xo,yo=yo))
}

setClass("BD", 
         slots = c(dir="character",
                   hazard="character",
                   buildingsfile='character',
                   cIndies="data.frame",
                   #fIndies="list",
                   alerts="data.frame",
                   I0="numeric",
                   hazdates="Date",
                   eventid="numeric",
                   #modifier="list",
                   hazinfo="list"),
         contains = "SpatVector")

setMethod(f="initialize", signature="BD",
          definition=function(.Object,Damage=NULL,ODD=NULL) {
            
            if(is.null(Damage)|is.null(ODD)) return(.Object)            
            
            .Object@hazard<-ODD@hazard
            #.Object@cIndies<-ODD@cIndies[ODD@cIndies$iso3%in%unique(Damage$iso3),]
            .Object@I0<-ODD@I0
            .Object@hazdates<-ODD@hazdates
            .Object@eventid<-ODD@eventid
            .Object@hazinfo<-ODD@hazinfo
            #.Object@fIndies<-ODD@fIndies
            
            .Object@buildingsfile<-paste0("./IIDIPUS_Input/OSM_Buildings_Objects/",unique(Damage$event)[1])
            
            print("Forming SpatVect object from building damage data")
            data_var <- c('grading')
            if (!is.null(Damage$Confidence)) data_var %<>% c('Confidence')
            if (!is.null(Damage$source)) data_var %<>% c('source')
            # Damage<-SpatialPointsDataFrame(coords = Damage[,c("Longitude","Latitude")],
            #                             data = Damage[,data_var],
            #                             proj4string = ODD@proj4string) #crs("+proj=longlat +datum=WGS84 +ellps=WGS84"))
            
            BD_points <- vect(x = as.matrix(Damage[,c("Longitude","Latitude")]), type='points', 
                 atts = Damage[,data_var], crs=crs(ODD))
            
            .Object@ptr <- BD_points@ptr
            # .Object@data <- Damage@data
            # .Object@coords.nrs <-Damage@coords.nrs
            # .Object@coords <-Damage@coords
            # .Object@bbox <-Damage@bbox
            # .Object@proj4string <- ODD@proj4string #crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            rm(Damage)
            
            print("Interpolating population density, hazard & GDP-PPP data")
            ODD[['spatial_pixel']] <- 1:ncell(ODD)
            .Object%<>%BDinterpODD(ODD=ODD)
            
            print("Filter spatial data per country")
            # We could just copy Damage$iso3 directly, but I don't believe in anyone or anything...
            .Object$ISO3C<-coords2country(crds(.Object))
            
            .Object@cIndies <- ODD@cIndies
            
            print("Accessing OSM to sample building height & area")
            # ExtractOSMbuildVol(.Object,ODD)
            
            #linp<-rep(list(1.),length(unique(ODD@cIndies$iso3)))
            #names(linp)<-unique(ODD@cIndies$iso3)
            #.Object@modifier<-linp
            
            print("Checking BD values")
            #checkBD(.Object)
            
            return(.Object)
          }
)

saveBD <- function(BD, path){
  BD_list = list()
  slotnames = slotNames(BD)
  pointer_slot <- ifelse('ptr' %in% slotnames, 'ptr', ifelse('pnt' %in% slotnames, 'pnt', 'pntr'))
  for (slot in slotnames[slotnames!=pointer_slot]){
    BD_list[[slot]] = slot(BD, slot)
  }
  BD_list$spatvector <- wrap(BD)
  saveRDS(BD_list, path)
}

readBD <- function(path){
  .Object = new('BD')
  BD_list = readRDS(path)
  slotnames <- slotNames(.Object)
  pointer_slot <- ifelse('ptr' %in% slotnames, 'ptr', ifelse('pnt' %in% slotnames, 'pnt', 'pntr'))
  for (slot in slotnames[slotnames!=pointer_slot]){
    slot(.Object, slot) = BD_list[[slot]]
  }
  slot(.Object, pointer_slot) = slot(terra::unwrap(BD_list$spatvector), pointer_slot) 
  names(.Object)[which(names(.Object)=='VALUE')] = 'ISO3C'
  return(.Object)
}

ODDI0poly<-function(ODD,var){
  
  isos<-unique(ODD@data$ISO3C)
  hazard<-ODD@hazard
  I0<-ODD@I0
  ODD<-ODD[var]
  names(ODD)<-"mean"
  
  pcontour<-suppressWarnings(tryCatch(adehabitatMA::getcontour(raster::subset(ODD,mean>=I0)),
                                      error=function(e) NULL))
  if(is.null(pcontour)) stop("Error in extracting contour by adehabitatMA::getcontour")
  
  conts<-data.frame()
  id<-1
  for(k in 1:length(pcontour@polygons)) {
    if(!pcontour@polygons[[k]]@Polygons[[1]]@area<1e-4) {
      conts%<>%rbind(data.frame(id=rep(id,length(pcontour@polygons[[k]]@Polygons[[1]]@coords[,1])),
                                Longitude=pcontour@polygons[[k]]@Polygons[[1]]@coords[,1],
                                Latitude=pcontour@polygons[[k]]@Polygons[[1]]@coords[,2]))
      id<-id+1
    }
    if(pcontour@polygons[[k]]@Polygons[[1]]@hole) 
      print(paste0("WARNING: hole in polygon of area ",
                   pcontour@polygons[[k]]@Polygons[[1]]@area,
                   " for ", hazard," event in countries: ",isos))
  }
  
  return(conts)
}

# This entire method disgusts me, I'm sorry if you are reading this.
# I'm running out of time, and apparently the capacity to think.
setGeneric("BDinterpODD", function(BD,ODD) 
  standardGeneric("BDinterpODD") )
setMethod("BDinterpODD", "BD", function(BD,ODD){
  
  if(sum(grepl("hazMean",names(ODD)))>0 & sum(grepl("hazSD",names(ODD)))>0 &
     sum(grepl("hazMean",names(ODD))) != sum(grepl("hazSD",names(ODD)))) 
    stop("Number of hazard mean values is not equal to number of hazard SD values")
  
  if(ODD@hazard=="TC"){
    stop('Not prepared for TC hazard type yet')
    for (var in names(ODD)){
      if(grepl("ISO3C",var)) next
      tmp<-ODD[var]%>%raster%>%raster::extract(BD@coords)%>%data.frame; colnames(tmp)<-var
      BD@data%<>%cbind(tmp)
    }
    ind<-unname(apply(BD@data,2, function(x) sum(!is.na(x))))
    BD@data<-BD@data[,ind>0]
    BDy@data<-BDy@data[(apply(BDy@data[,grep("hazMean",names(BDy),value = T)],1,function(x) !all(is.na(x)))),]
  } else { 
    
    polysave<-array(F,dim=c(nrow(BD),sum(grepl("hazMean",names(ODD)))))
    
    for (var in names(ODD)){
      if(grepl("ISO3C",var)) next
      #if(!grepl("haz",var)) {
      tmp<-ODD[[var]]%>%extract(crds(BD))%>%data.frame; colnames(tmp)<-var
      BD[[var]] <- tmp
      #}
    }
  # unsure what the rest of this was for?:
      # if(grepl("hazSD",var)) {
      #   # We deal with these later (just incase the ordering of ODD columns has been muddled about)
      #   next
      # }
      # # Find polygons of I>=I0
      # pcontour<-ODDI0poly(ODD,var)
      # 
      # insidepoly<-rep(F,nrow(BD))
      # if (length(unique(pcontour$id)) > 0){
      #   for(p in 1:length(unique(pcontour$id))){
      #     tcont<-filter(pcontour,id==p)
      #     insidepoly<-insidepoly | sp::point.in.polygon(BD@coords[,1, drop=F],
      #                                                   BD@coords[,2, drop=F],
      #                                                   tcont$Longitude,
      #                                                   tcont$Latitude)>0
      #   }
      # }
      # 
      # if(sum(insidepoly)==0) next
      # # Add the hazard to the building damage object
      # subBD<-ODD[var]%>%raster%>%raster::extract(BD@coords[insidepoly,, drop=F])
      # tmp<-rep(NA_real_,nrow(BD)); tmp[insidepoly]<-subBD; tmp%<>%data.frame; colnames(tmp)<-var
      # BD@data%<>%cbind(tmp)
      # polysave[,extractnumbers(var)]<-insidepoly
      
    # }
    
    # Add the standard deviations of the hazard intensity
  #   for (var in grep("hazSD",names(ODD),value = T)){
  #     if(sum(polysave[,extractnumbers(var)])==0) next
  #     subBD<-ODD[var]%>%raster%>%raster::extract(BD@coords[polysave[,extractnumbers(var)],, drop=F])
  #     tmp<-rep(NA_real_,nrow(BD)); tmp[polysave[,extractnumbers(var)]]<-subBD; tmp%<>%data.frame; colnames(tmp)<-var
  #     BD@data%<>%cbind(tmp)
  #   }
  #   
  #   if(length(sum(rowSums(polysave)==0))>0){
  #     print(paste0(sum(rowSums(polysave)==0)," building damage values were not exposed to a hazard... "))#,unique(miniDam$event)))
  #     BD%<>%raster::subset(rowSums(polysave)!=0)
  #   }
  # }
  } 
  return(BD)
  
})

Likefunc<-function(val,buildingV,coefs){
  xx<-val-buildingV
  return(1/(abs(coefs[5])*(xx)^8 + 
              abs(coefs[4])*(xx)^6 +   
              abs(coefs[3])*(xx)^4 + 
              abs(coefs[2])*(xx)^2 + 
              abs(coefs[1])))
}

OSMfunctionFIT<-function(buildings){
  
  findcoefsA<-function(coefs) {  
    Likelihood<-0
    indies<-(1:nrow(buildings))[!is.na(buildings$Pdensity)]
    sizer<-nrow(buildings[indies,])
    testers<-sample(indies,min(100,sizer))
    for (ij in testers){
      
      # Unweighted sample from all buildings
      miniBD<-sample(indies,size = min(300,sizer),replace = F)
      # Sample the building surface area from these buildings based weighted by the population density
      proby<-Likefunc(buildings$Pdensity[ij],buildings$Pdensity[miniBD],coefs)
      area<-buildings$area[sample(miniBD,1,prob=proby)]
      # Sample the building height from the buildings weighted by the building surface area
      # BD$Hb[ij]<-buildings$building.levels[sample(Ibuildlevel,1,replace = T,prob=1/((buildings$area[Ibuildlevel]-BD$Ab[ij])^3 + 1.2))]
      
      Likelihood%<>%add((area-buildings$area[ij])^2/sizer)
      
    }
    return(Likelihood)
  }
  
  findcoefsH<-function(coefs) {  
    Likelihood<-0
    indies<-(1:nrow(buildings))[!is.na(buildings$area) & !is.na(buildings$building.levels)]
    sizer<-nrow(buildings[indies,])
    testers<-sample(indies,min(100,sizer))
    for (ij in testers){
      # Unweighted sample from all buildings
      miniBD<-sample(indies,size = min(300,sizer),replace = F)
      # Sample the building surface area from these buildings based weighted by the population density
      proby<-Likefunc(buildings$area[ij],buildings$area[miniBD],coefs)
      hlevels<-buildings$building.levels[sample(miniBD,1,prob=proby)]
      
      Likelihood%<>%add((hlevels-buildings$building.levels[ij])^2/sizer)
      
    }
    return(Likelihood)
  }
  
  icoef<-c(1.1,1,0.1,0.01,0.001)
  outputA<-optim(icoef,findcoefsA)
  # outputH<-optim(icoef,findcoefsH,control = list(maxit = 10000))
  
  # return(list(coefsA=outputA$par,coefsH=outputH$par))
  return(outputA$par)
  
}

setGeneric("SampleBuildings", function(BD,buildings,samplecoefs) 
  standardGeneric("SampleBuildings") )
setMethod("SampleBuildings", "BD", function(BD,buildings=NULL,samplecoefs=F){
  
  BD@data<-BD@data[!is.na(BD@data$GDP),]
  BD@data$Population[is.na(BD@data$Population)]<-0
  
  if(is.null(buildings)) buildings<-readRDS(BD@buildingsfile)
  
  if(samplecoefs | is.null(BD@coefs)){
    coefs<-abs(OSMfunctionFIT(buildings))
    BD@coefs<-coefs
    
    return(BD)
  }
  
  # Do we not have enough OSM buildings to sample without replacement?
  repl<-min(300,nrow(buildings))==nrow(buildings) | nrow(BD)>nrow(buildings)
  
  BD$Ab<-rep(NA,nrow(BD))
  Ibuildlevel<-which(!is.na(buildings$building.levels))
  Ipop<-which(!is.na(buildings$Pdensity))
  lPop<-nrow(buildings[Ipop,]); llevel<-nrow(buildings[Ibuildlevel,])
  for (ij in 1:nrow(BD)){
    # Unweighted sample from all buildings
    miniBD<-sample(Ipop,size = min(300,lPop),replace = repl)
    # Sample the building surface area from these buildings based weighted by the population density
    proby<-Likefunc(buildings$Pdensity[miniBD],BD@data$Population[ij],BD@coefs)
    BD$Ab[ij]<-buildings$area[sample(miniBD,1,replace = F,prob=proby)]
    
    # miniBD<-sample(Ibuildlevel,size = min(300,llevel),replace = repl)
    # Sample the building height from the buildings weighted by the building surface area
    # proby<-Likefunc(buildings$area[miniBD],BD$Ab[ij],coefs$coefsH)
    # BD$Hb[ij]<-buildings$building.levels[sample(Ibuildlevel,1,prob=proby)]
    
  }
  
  return(BD)
  
})

# Note that bbox is most likely larger for ODD than BD object, hence bbox in function input
setGeneric("ExtractOSMbuildVol", function(BD,ODD) 
  standardGeneric("ExtractOSMbuildVol") )
setMethod("ExtractOSMbuildVol", "BD", function(BD,ODD){
  
  minnum<-max(nrow(BD)*0.05,30)
  buildings<-GetOSMbuildings(BD=BD,bbox=ODD@bbox,minnum=minnum)
  
  tmp<-ODD["Population"]%>%raster%>%raster::extract(buildings@coords)%>%data.frame; colnames(tmp)<-"Pdensity"
  buildings@data%<>%cbind(tmp) 
  buildings$Pdensity[is.na(buildings$Pdensity)]<-0
  
  print(paste0("Saving the buildings samples out to the file ",BD@buildingsfile))
  saveRDS(buildings,BD@buildingsfile)
  # BD%<>%SampleBuildings(buildings)
  
})

# Code that calculates/predicts the total human displacement 
setGeneric("BDX", function(BD,Omega,Model,Method, output='LL')
  standardGeneric("BDX") )
setMethod("BDX", "BD", function(BD,Omega,Model,Method=list(Np=20,cores=8), output='LL'){
  # Only calculate buildings with all key parameters
  
  #LL_cap <- -40
  hrange<-grep("hazMean",names(BD),value = T)
  
  elapsed_time <- c()
  start_time <- Sys.time()
  
  BD$order <- 1:ncell(BD)
  BD_df <- as.data.frame(BD)
  
  if(output != 'LL') {notnans<-which(!(is.na(BD_df$Population) | is.na(BD_df$ISO3C) | apply(BD_df[,hrange, drop=F], 1, function(row) all(is.na(row)))))
  } else notnans<-which(!(is.na(BD_df$Population) | is.na(BD_df$ISO3C) | is.na(BD_df$grading) | apply(BD_df[,hrange, drop=F], 1, function(row) all(is.na(row)))))
  
  #LOOSEEND: This would be faster if moved outside of BDX and have the BD objects resaved without possibly the damaged buildings
  #also can replace 'Minor', 'Moderate' and 'Severe' with just 'Damaged' outside and resave BD objects
  possiblyDamaged_ij <- which(BD_df$grading=='possible') 
  notnans <- notnans[!notnans %in% possiblyDamaged_ij]
  
  BD<-BD[notnans,]
  BD_df<-BD_df[notnans,] ; notnans<-1:nrow(BD_df)
  
  if(nrow(BD_df) ==0){
    if(output=='LL'){return(0)}
    else if (output=='sim'){
      return(BD)
    } else {
      return(array(0, dim=c(4, Method$Np)))
    }
  }
  # Get parameters for model
  
  # finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, GetNans = finish_time-start_time); start_time <- Sys.time()
  
  Params<-FormParams(BD,list(Np=Method$Np,center=Model$center))
  Params$I0 <- Model$I0 #some objects have different I0 but don't want this to affect the model
  
  # Income distribution percentiles & extract income percentile  
  #SincN<-paste0('p',seq(10,80,10), 'p', seq(20,90,10))
  #Sinc<-ExtractCIndy(BD,var = SincN)
  # Load buildings file
  # buildings<-readRDS(BD@buildingsfile)
  # Sample income distribution by area*building height?
  # BD%<>%SampleBuildings(buildings,F)
  # Calculate non-local linear predictor values
  LP<-GetLP(BD_df,Omega,Params,Sinc=NULL, notnans, split_GNI=F)
  
  event_lp <- getLP_event(BD@hazinfo, Omega, Params)
  
  # finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, GetLP = finish_time-start_time); start_time <- Sys.time()
  
  eps_event <- array(0, dim=c(length(hrange), Method$Np))
  for (h in 1:length(hrange)){
    eps_event[h,] <- stochastic(Method$Np,Omega$eps_adj$hazard_bd)
  }
  
  n_groups <- length(unique(BD_df$spatial_pixel))
  BD_df$spatial_pixel_id <- as.numeric(factor(BD_df$spatial_pixel, levels = unique(BD_df$spatial_pixel)))
  eps_local <- array(0, dim=c(n_groups, Method$Np))
  for (g in 1:n_groups){
    eps_local[g,] <- stochastic(Method$Np, Omega$eps_adj$local)
  }
  
  BD_df$LP <- LP
  BD_pixellated <- BD_df %>%
    group_by(spatial_pixel) %>%                    
    dplyr::summarise(
      tot_obs = n(),                     
      obs_prop_bd = 1 - mean(grading == 'notaffected'),
      LP = unique(LP),
      spatial_pixel_id = unique(spatial_pixel_id),
      across(starts_with("hazMean"), unique)
    )
  

  CalcBD_pixel <- function(ij){
    pUnaff <- rep(0, Method$Np)
    pUnaff_sampled <- rep(1, Method$Np)
    for (h_i in 1:length(hrange)){
      h <- hrange[h_i]
      #if(length(BD@data[ij,h])==0) next
      if(is.na(BD_pixellated[ij,h])) next
  
      I_ij <- BD_pixellated[ij,h]
      Damage <-fDamUnscaled_BD(I_ij,list(I0=Params$I0, Np=Method$Np),Omega) + as.numeric(BD_pixellated[ij, 'LP']) + eps_event[h_i,] + eps_local[as.numeric(BD_pixellated[ij,'spatial_pixel_id']), ] + event_lp[h_i]
      D_Dam <- D_Dam_calc(Damage, Omega)
      pUnaff_sampled <- pUnaff_sampled * (1-D_Dam)
    }
    return(1-pUnaff_sampled)
  }
  
  pDam_sampled <- matrix(unlist(mclapply(X = 1:nrow(BD_pixellated),FUN = CalcBD_pixel,mc.cores = Method$cores)), ncol=nrow(BD_pixellated), byrow=F)
  BD_pixellated[,paste0('pDamSamp.', 1:Method$Np)] = t(pDam_sampled)
  
  if (output=='results_analysis'){
    return(BD_pixellated)
  } 
  
  else if (output=='sim'){
    BD_df_withGrading <- merge(BD_pixellated[,c('spatial_pixel', 'pDamSamp.1')], BD_df, by='spatial_pixel')
    BD_df_withGrading$grading <- ifelse(rbern(NROW(BD_df_withGrading), BD_df_withGrading$pDamSamp.1) == 1, 'Damaged', 'notaffected')
    BD$grading <- BD_df_withGrading$grading
    return(BD)
  }
  # finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, GetEpsEvent = finish_time-start_time); start_time <- Sys.time()
  
  # to instead work building by building individually:
  
  # for each building in list,
  CalcBD<-function(ij){
    # Calculate local linear predictor (NOTE: is a scalar - we randomly sample one value)
    #locallinp<-sample(LP[ij,], Method$Np, replace=T)
    
    locallinp <- LP[ij]
    #tryCatch(sample(LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]],size=Method$Np, replace=TRUE)*
    #LP$Plinp[ij]*LP$linp[[iso3c]],         error=function(e) NA) #LOOSEEND: Assumes that a house is equally likely to be from each income bracket. 
    #locallinp<- LP$dGDP$linp[5]* LP$Plinp[ij]*LP$linp[[iso3c]]
    # if(is.na(locallinp)) stop(ij)
    # locallinp<-1.
    bDamage<-0
    bDamage<-rep(0, Method$Np) #0 = notaffected, 1 = damaged
    ind <- rep(T, Method$Np)
    p_notaff <- rep(1, Method$Np)
    for(h_i in 1:length(hrange)){
      h <- hrange[h_i]
      #if(length(BD@data[ij,h])==0) next
      if(is.na(BD_df[ij,h])) next
      
      # calculate the sampled hazard intensity I_ij
      # I_ij<-rnorm(n = Method$Np,
      #             mean = BD@data[ij,h],
      #             sd = BD@data[ij,paste0("hazSD",h)]/10)
      
      I_ij <- BD_df[ij,h]
      Damage <-fDamUnscaled_BD(I_ij,list(I0=Params$I0, Np=sum(ind)),Omega) + locallinp + eps_event[h_i,ind] + eps_local[BD_df$spatial_pixel_id[ij], ind] + event_lp[h_i]
      D_Dam <- D_Dam_calc(Damage, Omega)
      
      if (output=='LL' | output == 'results_analysis'){
        p_notaff <- p_notaff * (1-D_Dam)
      } else {
        bDamage[ind] <- rbinom(sum(ind), rep(1, sum(ind)), D_Dam) #vapply(1:sum(ind), sampleDamDest, numeric(1))
        ind <- bDamage != 1
        if(all(!ind)) break
        #bDamage[ind] <- apply(D_DestDamUnaf, 2, function(p){3- which(rmultinom(1, 1, p)==1)}) #slower
        #bDamage[ind] <- apply(D_DestDamUnaf, 2, sample, x=2:0, size=1, replace=F) #even slower: note that positional matching of arguments to sample
        #is overridden by matching names
        # if (h_i == 3){
        #   finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, pmaxbDam = finish_time-start_time); start_time <- Sys.time()
        # }
      }
      
    }
    
    if(output=='LL'){
      if (BD_df$grading[ij] == 'notaffected'){
        log_L <- log(p_notaff)
        return(ifelse(is.finite(log_L), log_L, LL_cap))
      } else {
        log_L <- log(1-p_notaff)
        return(ifelse(is.finite(log_L),log_L, LL_cap))
      }
    } else if (output=='results_analysis'){
      #return Long, Lat, Grading, Max Intensity, Simulated Probability of Damage
      coords <- crds(BD, na.rm=F)
      return(c(coords[ij,1], coords[ij,2], BD_df$grading[ij], max(BD_df[ij,hrange], na.rm=T), p_notaff))
    }
    #if(LL) return(bPred!=BD@data$grading[ij]) #return 1 if incorrectly classified
    
    bPred <- ifelse(bDamage == 1, 'Damaged','notaffected')
    
    if(output != 'sim'){
      if(BD_df$grading[ij] == 'notaffected'){ # S: refers to simulated value, O: refers to observed value
        return(ifelse(bPred == 'Damaged', 'S:Dam,O:NAff','S:Naff,O:NAff'))
      } else {
        return(ifelse(bPred == 'Damaged', 'S:Dam,O:Dam', 'S:Naff,O:Dam'))
      }
    }
    if(output == 'sim') return(bPred)
  }
  
  #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, DefCalcBD = finish_time-start_time); start_time <- Sys.time()
  #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, findPossiblyDam = finish_time-start_time); start_time <- Sys.time()
  
  if(output == 'LL') { # return the log likelihood
    if(Method$cores>1) {return(colSums(t(matrix(unlist(mclapply(X = notnans,FUN = CalcBD,mc.cores = Method$cores)),ncol=length(notnans)))))
    } else return(colSums(t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans)))))
  } 
  
  if(output == 'results_analysis') { #return the proportion of buildings correctly classified
    if(Method$cores>1) {return(t(matrix(unlist(mclapply(X = notnans,FUN = CalcBD,mc.cores = Method$cores)),ncol=length(notnans))))
    } else return(t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans))))
  }
  
  # classified<-t(matrix(unlist(mclapply(X = notnans,FUN = predBD,mc.cores = Method$cores)),ncol=length(notnans)))
  classified<-t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans)))
  
  #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, calcBD = finish_time-start_time); start_time <- Sys.time()
  
  if(output != 'sim'){ 
    #find values in contingency table between simulated and observed
    N11 <- colSums(classified=='S:Naff,O:NAff') #simulated: notaffected, observed: notaffected
    N12 <- colSums(classified=='S:Dam,O:NAff') #simulated: damaged, observed: notaffected
    N21 <- colSums(classified=='S:Naff,O:Dam')
    N22 <- colSums(classified=='S:Dam,O:Dam')
    return(rbind(N11, N12, N21, N22))
  }
  
  
  #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, classifications = finish_time-start_time); start_time <- Sys.time()
  
  #return(elapsed_time)
  
  
  # Therefore, sim==F and LL==F
  # Save into the file
  BD$ClassPred<-classified
  
  return(BD)
  
})

# # Code that calculates/predicts the total human displacement 
# setGeneric("BDX", function(BD,Omega,Model,Method, output='LL')
#   standardGeneric("BDX") )
# setMethod("BDX", "BD", function(BD,Omega,Model,Method=list(Np=20,cores=8), output='LL'){
#   # Only calculate buildings with all key parameters
#   
#   #LL_cap <- -40
#   hrange<-grep("hazMean",names(BD),value = T)
#   
#   elapsed_time <- c()
#   start_time <- Sys.time()
#   
#   
#   BD_df <- as.data.frame(BD)
#   
#   if(output != 'LL') {notnans<-which(!(is.na(BD$Population) | is.na(BD$ISO3C) | apply(BD_df[,hrange], 1, function(row) all(is.na(row)))))
#   } else notnans<-which(!(is.na(BD@data$Population) | is.na(BD@data$ISO3C) | is.na(BD@data$grading) | apply(BD_df[,hrange], 1, function(row) all(is.na(row)))))
#   
#   #LOOSEEND: This would be faster if moved outside of BDX and have the BD objects resaved without possibly the damaged buildings
#   #also can replace 'Minor', 'Moderate' and 'Severe' with just 'Damaged' outside and resave BD objects
#   possiblyDamaged_ij <- which(BD@data$grading=='possible') 
#   notnans <- notnans[!notnans %in% possiblyDamaged_ij]
#   
#   BD<-BD[notnans,] ; notnans<-1:nrow(BD)
#   
#   if(nrow(BD) ==0){
#     if(output=='LL'){return(0)}
#     else if (output=='sim'){
#       return(BD)
#     } else {
#       return(array(0, dim=c(4, Method$Np)))
#     }
#   }
#   # Get parameters for model
#   
#   # finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, GetNans = finish_time-start_time); start_time <- Sys.time()
#   
#   Params<-FormParams(BD,list(Np=Method$Np,center=Model$center))
#   Params$I0 <- Model$I0 #some objects have different I0 but don't want this to affect the model
#   
#   # Income distribution percentiles & extract income percentile  
#   #SincN<-paste0('p',seq(10,80,10), 'p', seq(20,90,10))
#   #Sinc<-ExtractCIndy(BD,var = SincN)
#   # Load buildings file
#   # buildings<-readRDS(BD@buildingsfile)
#   # Sample income distribution by area*building height?
#   # BD%<>%SampleBuildings(buildings,F)
#   # Calculate non-local linear predictor values
#   LP<-GetLP(BD,Omega,Params,Sinc=NULL, notnans, split_GNI=F)
#   
#   # finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, GetLP = finish_time-start_time); start_time <- Sys.time()
#   
#   eps_event <- array(0, dim=c(length(hrange), Method$Np))
#   for (h in 1:length(hrange)){
#     eps_event[h,] <- stochastic(Method$Np,Omega$eps_adj$hazard_bd)
#   }
#   
#   n_groups <- length(unique(BD$spatial_pixel))
#   BD$spatial_pixel_id <- as.numeric(factor(BD$spatial_pixel, levels = unique(BD$spatial_pixel)))
#   eps_local <- array(0, dim=c(n_groups, Method$Np))
#   for (g in 1:n_groups){
#     eps_local[g,] <- stochastic(Method$Np, Omega$eps_adj$local)
#   }
#   
#   
#   # finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, GetEpsEvent = finish_time-start_time); start_time <- Sys.time()
#   
#   # for each building in list,
#   CalcBD<-function(ij){
#     # Calculate local linear predictor (NOTE: is a scalar - we randomly sample one value)
#     #locallinp<-sample(LP[ij,], Method$Np, replace=T)
#     
#     locallinp <- LP[ij]
#     #tryCatch(sample(LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]],size=Method$Np, replace=TRUE)*
#     #LP$Plinp[ij]*LP$linp[[iso3c]],         error=function(e) NA) #LOOSEEND: Assumes that a house is equally likely to be from each income bracket. 
#     #locallinp<- LP$dGDP$linp[5]* LP$Plinp[ij]*LP$linp[[iso3c]]
#       # if(is.na(locallinp)) stop(ij)
#     # locallinp<-1.
#     bDamage<-0
#     bDamage<-rep(0, Method$Np) #0 = notaffected, 1 = damaged
#     ind <- rep(T, Method$Np)
#     p_notaff <- rep(1, Method$Np)
#     for(h_i in 1:length(hrange)){
#       h <- hrange[h_i]
#       #if(length(BD@data[ij,h])==0) next
#       if(is.na(BD@data[ij,h])) next
#       
#       # calculate the sampled hazard intensity I_ij
#       # I_ij<-rnorm(n = Method$Np,
#       #             mean = BD@data[ij,h],
#       #             sd = BD@data[ij,paste0("hazSD",h)]/10)
#       
#       I_ij <- BD@data[ij,h]
#       Damage <-fDamUnscaled_BD(I_ij,list(I0=Params$I0, Np=sum(ind)),Omega) + locallinp + eps_event[h_i,ind] + eps_local[BD@data$spatial_pixel_id[ij], ind]
#       D_Dam <- D_Dam_calc(Damage, Omega)
#       
#       if (output=='LL' | output == 'results_analysis'){
#         p_notaff <- p_notaff * (1-D_Dam)
#       } else {
#         bDamage[ind] <- rbinom(sum(ind), rep(1, sum(ind)), D_Dam) #vapply(1:sum(ind), sampleDamDest, numeric(1))
#         ind <- bDamage != 1
#         if(all(!ind)) break
#         #bDamage[ind] <- apply(D_DestDamUnaf, 2, function(p){3- which(rmultinom(1, 1, p)==1)}) #slower
#         #bDamage[ind] <- apply(D_DestDamUnaf, 2, sample, x=2:0, size=1, replace=F) #even slower: note that positional matching of arguments to sample
#         #is overridden by matching names
#         # if (h_i == 3){
#         #   finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, pmaxbDam = finish_time-start_time); start_time <- Sys.time()
#         # }
#       }
#       
#     }
#     
#     if(output=='LL'){
#       if (BD@data$grading[ij] == 'notaffected'){
#         log_L <- log(p_notaff)
#         return(ifelse(is.finite(log_L), log_L, LL_cap))
#       } else {
#         log_L <- log(1-p_notaff)
#         return(ifelse(is.finite(log_L),log_L, LL_cap))
#       }
#     } else if (output=='results_analysis'){
#       #return Long, Lat, Grading, Max Intensity, Simulated Probability of Damage
#       return(c(BD@coords[ij,1], BD@coords[ij,2], BD@data$grading[ij], max(BD@data[ij,hrange], na.rm=T), p_notaff))
#     }
#     #if(LL) return(bPred!=BD@data$grading[ij]) #return 1 if incorrectly classified
#     
#     bPred <- ifelse(bDamage == 1, 'Damaged','notaffected')
#     
#     if(output != 'sim'){
#       if(BD@data$grading[ij] == 'notaffected'){ # S: refers to simulated value, O: refers to observed value
#         return(ifelse(bPred == 'Damaged', 'S:Dam,O:NAff','S:Naff,O:NAff'))
#       } else {
#         return(ifelse(bPred == 'Damaged', 'S:Dam,O:Dam', 'S:Naff,O:Dam'))
#       }
#     }
#     if(output == 'sim') return(bPred)
#   }
#   
#   #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, DefCalcBD = finish_time-start_time); start_time <- Sys.time()
#   #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, findPossiblyDam = finish_time-start_time); start_time <- Sys.time()
#   
#   if(output == 'LL') { # return the log likelihood
#     if(Method$cores>1) {return(colSums(t(matrix(unlist(mclapply(X = notnans,FUN = CalcBD,mc.cores = Method$cores)),ncol=length(notnans)))))
#     } else return(colSums(t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans)))))
#   } 
#   
#   if(output == 'results_analysis') { #return the proportion of buildings correctly classified
#     if(Method$cores>1) {return(t(matrix(unlist(mclapply(X = notnans,FUN = CalcBD,mc.cores = Method$cores)),ncol=length(notnans))))
#     } else return(t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans))))
#   }
#   
#   # classified<-t(matrix(unlist(mclapply(X = notnans,FUN = predBD,mc.cores = Method$cores)),ncol=length(notnans)))
#   classified<-t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans)))
#   
#   #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, calcBD = finish_time-start_time); start_time <- Sys.time()
#   
#   if(output != 'sim'){ 
#     #find values in contingency table between simulated and observed
#     N11 <- colSums(classified=='S:Naff,O:NAff') #simulated: notaffected, observed: notaffected
#     N12 <- colSums(classified=='S:Dam,O:NAff') #simulated: damaged, observed: notaffected
#     N21 <- colSums(classified=='S:Naff,O:Dam')
#     N22 <- colSums(classified=='S:Dam,O:Dam')
#     return(rbind(N11, N12, N21, N22))
#   }
#   
#   
#   #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, classifications = finish_time-start_time); start_time <- Sys.time()
#   
#   #return(elapsed_time)
#   
#   
#   # Therefore, sim==F and LL==F
#   # Save into the file
#   BD$ClassPred<-classified
#   
#   return(BD)
#   
# })
