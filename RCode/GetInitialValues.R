
ExtractCentering<-function(dir, haz="EQ",saver=T, input_folder='IIDIPUS_Input/'){
  
  if(saver & file.exists(paste0(dir, input_folder, "centerings"))) 
        return(readRDS(paste0(dir, input_folder, "centerings")))
  
  # path<-paste0("/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/")
  # ufiles<-list.files(path=path,pattern=haz,recursive = T,ignore.case = T)
  # ufiles<-ufiles[grepl(ufiles,pattern = haz)]
  # PDens<-c()
  # Vs30 <- c()
  # for(fff in ufiles){
  #   ODDy<-readRDS(paste0(path,fff))
  #   PDens<-append(PDens, ODDy@data$Population[!is.na(ODDy@data$Population)])
  #   Vs30<-append(Vs30, ODDy@data$Vs30[!is.na(ODDy@data$Vs30)])
  # }
  # mean(log(PDens+0.1)); sd(log(PDens+0.1))
  # mean(log(Vs30)); sd(log(Vs30))
  
  PDens_mean <- 2.197932 #mean(log(PDens+0.1))
  PDens_sd <- 2.349322 #sd(log(PDens+0.1))
  
  Vs30_mean <- 6.23811 #old: cellStats(stiff,'mean')
  Vs30_sd <- 0.4035713 #old: cellStats(stiff,'sd')
  
  # # Read in Stiff data and calculate the mean and standard deviation, again using all regions in the dataset:
  # if(!file.exists(paste0(dir,"Hazard_Data/global_vs30_tif/global_vs30.tif"))) stop("Please download the VS30 dataset (geotiff and auxiliary files) here https://earthquake.usgs.gov/data/vs30/.")
  #stiff<-raster(paste0(dir,"Hazard_Data/global_vs30_tif/global_vs30.tif"))
  #stiffAgg <- aggregate(stiff, 5) #Need to aggregate a bit otherwise quantile crashes R
  #raster::quantile(stiffAgg, probs=c(0.01,0.99))
  
  # # Read in Global Data Lab data and calculate the mean and standard deviation of each variable
  # # Note that we calculate the mean and sd using all regions (not just those in the training set)
  GDLdata <- readGlobalDataLab()
  
  AveSchYrs_mean <- 7.014779 #mean(GDLdata$AveSchYrs, na.rm=T); 
  AveSchYrs_sd <- 3.407637 #sd(GDLdata$AveSchYrs);
  LifeExp_mean <- 68.30403 #mean(GDLdata$LifeExp); 
  LifeExp_sd <- 9.607136 #sd(GDLdata$LifeExp);
  GNIc_mean <- 8.913694 #mean(log(GDLdata$GNIc)); 
  GNIc_sd <- 1.18729 #sd(log(GDLdata$GNIc));
  SHDI_mean <- 0.6453444#mean(GDLdata$SHDI, na.rm=T);
  SHDI_sd <- 0.1712748#sd(GDLdata$SHDI, na.rm=T);
  
  
 
  
  # # Read in GDPA data and calculate the mean and standard deviation, again using all regions in the dataset:
  # 
  # if(!file.exists(paste0(dir,"Hazard_Data/gdpga/gdpga.asc"))) stop("Please download the hazard frequency data here: https://sedac.ciesin.columbia.edu/data/set/ndh-earthquake-distribution-peak-ground-acceleration/data-download.  Note you'll have to create a SEDAC account.")
  # pga<-SortDemoData(paste0(dir,"Hazard_Data/gdpga/gdpga.asc"))
  # pga%<>%convMat2SPDF(name="PGA")
  
  #pga <- raster(paste0(dir,"Hazard_Data/gdpga/pga_475y.tif"))
  #pga_vals <- values(pga)
  
  # EQfreqs <- c()
  # Mags <- c()
  # for (file in ufiles){
  #   ODDy <- readRDS(paste0(folderin, file))
  #   EQfreqs <- c(EQfreqs, log(median(ODDy$EQFreq)+0.1))
  #   #Mags <- c(Mags, max(ODDy@hazinfo$magnitudes))
  # }
  # mean(EQfreqs)
  # sd(EQfreqs)
  #pga <- rast(paste0(dir,"Hazard_Data/GEM-GSHM_PGA-475y-rock_v2023/v2023_1_pga_475_rock_3min.tif"))
  #pga_vals <- values(pga)
  #pga_vals_nonzero <- pga_vals[pga_vals>0]
  EQFreq_mean <- -3.123102 #mean(log(pga_vals_nonzero+0.001), na.rm=T) # mean(EQfreqs) #2.568421 #mean(log(pga_vals+1)[which(pga_vals > 0.01)]) 5.411167 # mean(pga$PGA, na.rm=T) 
  EQFreq_sd <- 1.223276 #sd(log(pga_vals_nonzero+0.001), na.rm=T) # sd(EQfreqs) #1.860136 #sd(log(pga_vals+1)[which(pga_vals > 0.01)]) 2.918439 # sd(pga$PGA, na.rm=T)
  
  #mean and sd of magnitudes in data (max magnitude taken from each event)
  Mag_mean <- 6.187425
  Mag_sd <- 0.7543592
  
  FirstHaz_mean <- 0.5
  FirstHaz_sd <- 0.5
  
  Night_mean <- 1/3
  Night_sd <- sqrt(2/9)
  
  FirstHaz.Night_mean <- 1/6
  FirstHaz.Night_sd <- sqrt(434/3125)
  
    
  center<-list(PDens=list(mean=PDens_mean, sd=PDens_sd), #LOOSEEND: CAME UP WITH THIS SD. CHECK
               AveSchYrs=list(mean=AveSchYrs_mean, sd=AveSchYrs_sd),
               LifeExp=list(mean=LifeExp_mean, sd=LifeExp_sd),
               GNIc=list(mean=GNIc_mean, sd=GNIc_sd),
               Vs30=list(mean=Vs30_mean, sd=Vs30_sd),
               EQFreq=list(mean=EQFreq_mean, sd=EQFreq_sd),
               Mag=list(mean=Mag_mean, sd=Mag_sd), 
               SHDI=list(mean=SHDI_mean, sd=SHDI_sd),
               FirstHaz = list(mean=FirstHaz_mean, sd=FirstHaz_sd), 
               Night = list(mean = Night_mean, sd=Night_sd),
               FirstHaz.Night = list(mean=FirstHaz.Night_mean, sd=FirstHaz.Night_sd))
  
  # center<-list(Gov=98.7,Vuln=51,CC=44,MPI=53.7,Pinf=103,Pexp=112,Sinc=0.2152956,Ik=0.4,A=3.6,H=1.65)
  print(unlist(center))
  
  saveRDS(center,paste0(dir,input_folder, "centerings"))
  
  return(center)
}

GetInitVals<-function(ODDpath, Model, AlgoParams, usePastPost=T, pastPostFile='mcmc_2025-04-10'){
  # Produces initial proposal and perturbation covariance for the MCMC chains
  # if usePastPost = T, then uses the posterior from another chain to estimate the proposal covariance
  # if usePastPost = F, then uses prior samples to estimate the proposal covariance
  
  if (usePastPost){
    if(length(pastPostFile)==0){
      stop('Please provide AlgoResults file from which to obtain the posterior')
    }
    post_samples <- t(readRDS(paste0(dir, 'IIDIPUS_Results/', pastPostFile))$Omega_sample)
    s_finish = which(is.na(post_samples)[,1])[1] - 1
    post_samples = post_samples[max(1, s_finish-1000):s_finish,]
    propCOV = cov(post_samples)
    init_value = post_samples[nrow(post_samples),]
    iVals = list(x0 = init_value, COV=propCOV/10)
    
  } else {
    n_samples <- 500
    hp_samples <- array(NA, dim=c(n_samples,length(Model$par_lb)))
    for (i in 1:n_samples){
      hp_samples[i, ] = HLPrior_sample(Model, AlgoParams)
    }
    iVals = list(x0=hp_samples[1,],COV=cov(hp_samples)/10)
  }
  return(iVals)
}

HLPrior_sample <- function(Model, AlgoParams){
  #Draws a sample from the prior that accepts the Higher Level Prior
  #Returns sample on the transformed space (not the physical space)
  n_x <- length(unlist(Model$links))
  HP <- AlgoParams$ABC + 1
  sample <- rep(0, n_x)
  while (HP > AlgoParams$ABC){
    s <- 1
    for (i in 1:length(Model$Priors)){
      for (j in 1:length(Model$Priors[[i]])){
        prior_dist <- Model$Priors[[i]][[j]]
        sample[s] <- do.call(match.fun(paste0('r', prior_dist$dist)), c(list(n=1), prior_dist[2:length(prior_dist)]))
        s <- s + 1
      }
    }
    HP <- Model$HighLevelPriors(relist(sample,skeleton=Model$skeleton) %>% addTransfParams(),Model) #check higher level prior of the proposal
  }
  return(sample %>% relist(Model$skeleton)%>% Physical2Proposed(Model) %>% unlist())
}

# 
# params_store <- array(0, dim=c(1000,19))
# sums <- c()
# for (i in 1:1000){
#   print(i)
#   params_store[i, ] <- HLPrior_sample(Model, AlgoParams)
#   sums <- c(sums, sum(abs(params_store[i,13:19])))
# }
# saveRDS(params_store, paste0(dir, 'IIDIPUS_Input/HLPriorSamples2'))

# plot(params_store[,1], params_store[,2])
# omega_samples <- array(0, dim=c(500, 14))
# for (i in 1:500){
#   print(i)
#   #Draws a sample from the prior that accepts the Higher Level Prior
#   #Returns sample on the transformed space (not the physical space)
#   n_x <- length(unlist(Model$links))
#   HP <- AlgoParams$ABC + 1
#   sample <- rep(0, n_x)
#   while (HP > AlgoParams$ABC){
#     sample <- runif(n_x, Model$par_lb, Model$par_ub) #generate proposal on the physical space
#     sample[c(2,4,6,8)] <- (exp(sample[c(2,4,6,8)]*sample[13])-1)/6
#     HP <- Model$HighLevelPriors(relist(sample,skeleton=Model$skeleton) %>% addTransfParams(),Model) #check higher level prior of the proposal
#   }
#   omega_samples[i,] <- sample
# }
# 
# pairs(omega_samples)
# hist(omega_samples[,2])
# plot(log(6*omega_samples2[,2]+1)/omega_samples2[,13])
# 
# plot(omega_samples[,13], log(6*omega_samples[,8]+1)/omega_samples[,13])
# 
# plot(seq(4.5,10, 0.01), pnorm(exp(omega_samples[1,13]*(seq(4.5,10, 0.01) - 4.5))-1,exp(omega_samples[1,13]*(omega_samples[1,3]-4.5))-1,omega_samples[1,4]), type='l')
# for(i in 59:59){
#   lines(seq(4.5,10, 0.01), pnorm(exp(omega_samples[i,13]*(seq(4.5,10, 0.01) - 4.5))-1,exp(omega_samples[i,13]*(omega_samples[i,3]-4.5))-1,omega_samples[i,4]), type='l')
# }
# plot(omega_samples[,2])
# Model$par_ub




















# For hazards where only one event is available
# SingleEventInits<-function(ODDy,BDy,Omega,haz="EQ", Model, AlgoParams){
#   
#   dollar<-Omega$beta$dollar
#   Omega$beta[1:8]<-0 ; Omega$beta$dollar<-dollar
#   Omega$eps$xi<-0
#   
#   Fopty<-function(vals){
#     
#     vals%<>%exp()
#     vals[8]<--vals[8]
#     
#     print(vals)
#     
#     Omega$Lambda[1:3]<-vals[1:3]
#     Omega$zeta[1:2]<-vals[4:5]
#     Omega$theta[1]<-vals[6]
#     Omega$eps[1]<-vals[7]
#     Omega$beta$dollar<-vals[8]
#     
#     HP<-0
#     # Apply higher order priors
#     if(!is.null(Model$HighLevelPriors)){
#       HP<-Model$HighLevelPriors(Omega,Model)
#       print(paste0("Higher Level Priors = ",HP))
#     }
#     
#     LL<-HP
#     ### ODD Object! ###
#     tLL<-tryCatch(DispX(ODD = ODDy,Omega = Omega,center = Model$center,LL = T,Method = AlgoParams),
#                   error=function(e) NA)
#     if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));return(-Inf)}
#     if(expLL) LL<-LL+max(log(mean(exp(tLL),na.rm=T)),AlgoParams$cap,na.rm = T)
#     else LL<-LL+max(mean(tLL,na.rm=T),AlgoParams$cap,na.rm = T)
#     LL<-LL
#     print(paste0("LL Disp = ",LL-HP))
#     sLL<-LL
#     
#     ### BD Object! ###
#     tLL<-tryCatch(BDX(BD = BDy,Omega = Omega,center = Model$center,Method=AlgoParams),
#                   error=function(e) NA)
#     # tLL<-BDX(BD = BDy,Omega = Omega,center = Model$center,Method=AlgoParams)
#     # if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));next}
#     if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));return(-Inf)}
#     # print(paste0(ufiles[i]," BD LL = ",mean(tLL,na.rm = T)))
#     if(expLL) LL<-LL+max(log(mean(exp(tLL),na.rm=T)),AlgoParams$cap,na.rm = T)
#     else LL<-LL+max(mean(tLL,na.rm=T),AlgoParams$cap,na.rm = T)
#     
#     print(paste0("LL BD = ",LL-sLL))
#     
#     # Add priors
#     if(!is.null(Model$priors)){
#       LL<-LL+sum(Priors(Omega,Model$priors),na.rm=T)
#     }
#     LL<-LL+HP
#     
#     print(LL)
#     print("...")
#     return(LL)
#     
#   } 
#   
#   ivals<-log(unname(unlist(Omega[c("Lambda","zeta","theta","eps")])))
#   ivals[8]<-log(-Omega$beta$dollar)
#   
#   output<-optim(par=ivals,
#                 fn = Fopty,control = list(maxit = 150,fnscale=-1,reltol=1.5e-4),hessian = T)  
#   
# }



# ExtractCentering<-function(dir,Model,haz="EQ",saver=T){
#   
#   if(saver) return(readRDS(paste0(dir,"IIDIPUS_Input/centerings")))
#   
#   # Indicators<-INFORMmean(Model$INFORM_vars); Indicators%<>%arrange(variable)
#   # print("Re-extracting the INFORM indicator values for 2019 and replacing interpolated values")
#   
#   lennie<-length(Model$INFORM_vars)
#   
#   ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/ODDobjects_count"),pattern=haz,recursive = T,ignore.case = T)
#   # ufiles<-str_split_fixed(ufiles,"_",2) ; ufiles<-as.integer(ufiles[,2])
#   DF<-rep(0,length(Model$INFORM_vars)+1) # (including Ik too!)
#   GDP<-c(dollar=0)
#   Pdens<-c(Pdens=0)
#   nGDP<-nDF<-nPdens<-0
#   for(i in 1:length(ufiles)){
#     print(paste0(i," :   ",(str_split_fixed(ufiles[i],"_",2))[1]))
#     ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects_count/",ufiles[i]))
#     ODDy@gmax%<>%as.data.frame.list()
#     saveRDS(ODDy,paste0(dir,"IIDIPUS_Input/ODDobjects_count/",ufiles[i]))
#     
#     ind<-ODDy@cIndies$iso3%in%unique(ODDy@data$ISO3C[!is.na(ODDy@data$Population)])
#     isos<-unique(ODDy@cIndies$iso3[ind])
#     
#     cIndies<-filter(ODDy@cIndies,variable%in%
#                       c(paste0("p",(1:9)*10,"p100"),
#                         "Ik")) 
#     # tmp<-filter(Indicators,iso3%in%isos)
#     # if(nrow(tmp)!=length(isos)*lennie) print(paste0("WARNING: missing INFORM variables for ", ufiles[i]))
#     
#     if(nrow(cIndies)>0){
#       checker<-as.data.frame(table(cIndies$variable))
#       checker%<>%filter(Freq<max(Freq))%>%pull(Var1)%>%as.character()
#     } else {checker<-Model$WID_perc}
#     
#     if(length(checker)>0){
#       print(paste0("WARNING: missing ",checker," variables for ", ufiles[i]))
#       # cisos<-isos[!isos%in%cIndies$iso3[cIndies$variable==checker]]
#       # WID<-GetWID_perc(checker,cisos,AsYear(ODDy@hazdates[1]))
#       # if(nrow(WID)==0) {print(paste0("WARNING: skipping earthquake ", ufiles[i])) ; next}
#       # cIndies%<>%rbind(WID)
#     }
#     
#     # ODDy@cIndies<-rbind(cIndies,filter(Indicators,iso3%in%isos))%>%
#     #               arrange(variable)
#     
#     # saveRDS(ODDy,paste0(dir,"IIDIPUS_Input/ODDobjects/",ufiles[i]))
#     
#     aind<-ODDy@cIndies$iso3%in%unique(ODDy@data$ISO3C[!is.na(ODDy@data$Population)])
#     ind<-aind & !ODDy@cIndies$variable%in%paste0("p",(1:9)*10,"p100")
#     
#     DF<-DF+colSums(matrix(ODDy@cIndies$value[ind],
#                           nrow = length(isos),
#                           dimnames = list(isos,unique(ODDy@cIndies$variable[ind]))))
#     nDF<-nDF+length(isos)
#     
#     Pdens<-Pdens+sum(log(ODDy@data$Population[ODDy@data$Population>0]),na.rm=T)
#     nPdens<-nPdens+length(ODDy@data$Population[ODDy@data$Population>0 & !is.na(ODDy@data$Population)])
#     
#     GDP<-GDP+sum(log(ODDy@data$GDP[ODDy@data$GDP>0]),na.rm=T)
#     nGDP<-nGDP+length(ODDy@data$GDP[ODDy@data$GDP>0 & !is.na(ODDy@data$GDP)])
#     # tmp<-unique(ODDy@data$GDP);tmp<-tmp[!is.na(tmp)]
#     # tmp<-log(tmp%o%ODDy@cIndies$value[aind & ODDy@cIndies$variable%in%paste0("p",(1:9)*10,"p100")])
#     # GDP<-GDP+sum(log(tmp))
#     # nGDP<-nGDP+length(tmp)
#     
#     rm(ODDy)
#     
#   }
#   
#   # Do the same for building damage objects BD
#   ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/BDobjects"),pattern=haz,recursive = T,ignore.case = T)
#   A<-c(A=0)
#   H<-c(H=0)
#   nA<-nH<-0
#   for(i in 1:length(ufiles)){
#     print(paste0(i," :   ",(str_split_fixed(ufiles[i],"_",2))[1]))
#     BDy<-readRDS(paste0(dir,"IIDIPUS_Input/BDobjects/",ufiles[i]))
#     # ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/",ufiles[i]))
#     
#     # BDy@cIndies<-ODDy@cIndies
#     BDy@fIndies<-Model$fIndies
#     saveRDS(BDy,paste0(dir,"IIDIPUS_Input/BDobjects/",ufiles[i]))
#     
#     H<-H+sum(log(BDy@data$Population[BDy@data$Population>0]),na.rm=T)
#     nH<-nH+length(BDy@data$Population[BDy@data$Population>0 & !is.na(BDy@data$Population)])
#     
#     buildings<-readRDS(BDy@buildingsfile)
#     A<-A+sum(log(buildings$area[buildings$area>0]),na.rm=T)
#     nA<-nA+length(buildings$area[!is.na(buildings$area)])
#     rm(buildings)
#     
#   }
#   
#   center<-as.list(c(DF[Model$INFORM_vars]/nDF,GDP/nGDP,log(301),DF["Ik"]/nDF,A/nA,H/nH))
#   # center<-list(Gov=98.7,Vuln=51,CC=44,MPI=53.7,Pinf=103,Pexp=112,Sinc=0.2152956,Ik=0.4,A=3.6,H=1.65)
#   print(unlist(center))
#   
#   saveRDS(center,paste0(dir,"IIDIPUS_Input/centerings"))
#   
#   return(center)
# }


# path<-paste0("/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NoBuildingDat/ODDobjects/")
# ufiles<-list.files(path=path,pattern=haz,recursive = T,ignore.case = T)
# ufiles<-ufiles[grepl(ufiles,pattern = haz)]
# 
# for(fff in ufiles){
#   ODDy<-readRDS(paste0(path,fff))
# 
#   if(any(is.na(ODDy$GNIc)&!is.na(ODDy$ISO3C))) print(fff)
# }


# rbPal <- colorRampPalette(c('red','blue'))
# plot(ODDy@coords[ODDy@data$Population>0,], col=rbPal(20)[as.numeric(cut(ODDy$EQFreq[ODDy@data$Population>0],breaks = 20))])

# for(fff in ufiles){
#   ODDy<-readRDS(paste0(path,fff))
#   missing_GDL_observed_Pop <- which(!is.na(ODDy@data$Population) & is.na(ODDy@data$AveSchYrs))
#   missing_GDL <- which(is.na(ODDy@data$AveSchYrs))
#   #assign these with the GDL of the closest pixel if within a distance of 2 arcminutes
#   for (i in missing_GDL_observed_Pop){
#     closest <- which.min(rowSums(sweep(ODDy@coords[-missing_GDL,], 2, ODDy@coords[i,], FUN='-')^2 ))
#     dist <- sqrt(sum((ODDy@coords[i,] - ODDy@coords[-missing_GDL,][closest,])^2))
#     if (dist < 0.05){
#       ODDy@data$AveSchYrs[i] <- ODDy@data$AveSchYrs[-missing_GDL][closest]
#       ODDy@data$LifeExp[i] <- ODDy@data$LifeExp[-missing_GDL][closest]
#       ODDy@data$GNIc[i] <- ODDy@data$GNIc[-missing_GDL][closest]
#     } else {
#       stop(paste('GDL data not found for coordinate',ODDy@coords[i,1], ODDy@coords[i,2]))
#     }
#   }
#   
#   saveRDS(ODDy, paste0(path,fff))
# }


