
ExtractCentering<-function(dir,haz="EQ",saver=T){
  
  if(saver & file.exists(paste0(dir,"IIDIPUS_Input/centerings_V2"))) 
        return(readRDS(paste0(dir,"IIDIPUS_Input/centerings_V2")))
  
  path<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  ufiles<-list.files(path=path,pattern=haz,recursive = T,ignore.case = T)
  ufiles<-ufiles[grepl(ufiles,pattern = haz)]
  GDP<-nGDP<-0
  for(fff in ufiles){
    ODDy<-readRDS(paste0(path,fff))

    GDP<-GDP+sum(log(ODDy@data$GDP[ODDy@data$GDP>0]),na.rm=T)
    nGDP<-nGDP+length(ODDy@data$GDP[ODDy@data$GDP>0 & !is.na(ODDy@data$GDP)])
    
  }
  rm(ODDy)
  
  center<-list(Pdens=log(301),dollar=GDP/nGDP)
  # center<-list(Gov=98.7,Vuln=51,CC=44,MPI=53.7,Pinf=103,Pexp=112,Sinc=0.2152956,Ik=0.4,A=3.6,H=1.65)
  print(unlist(center))
  
  saveRDS(center,paste0(dir,"IIDIPUS_Input/centerings_V2"))
  
  return(center)
}

GetInitVals<-function(ODDpath, Model, AlgoParams, optimiser=F){
  Omega<-list(
    Lambda=list(kappa=0.03786569,nu=0.19768294,omega=0.09954531),
    zeta=list(k=2.02265867,lambda=5.40390699), # zeta=list(k=2.5,lambda=1.6),
    # beta=list(CC.INS.GOV.GE=0,VU.SEV.AD=0,CC.INS.DRR=0,VU.SEV.PD=0,CC.INF.PHY=0,HA.NAT.EQ=0),
    Pdens=list(M=0.05,k=1.),
    dollar=list(M=0.05,k=1.),
    theta=list(e=0.67431138), #list(e=0.25),
    # rho=list(A=0,H=0),
    eps=list(eps=0.12709657)#,xi=3.52269924)
    # mu=list(muplus=1,muminus=1,sigplus=0.001,sigminus=0.001)
  )

  lenny<-length(unlist(Omega))
  # Note that this is in proposal space, not physical space, hence the logarithms
  propCOV<-diag(rep(lenny/60,lenny))
  # propCOV<-diag(abs(c(log(0.17/0.12),log(0.04/0.02),log(0.09/0.055),
  #                     log(2/1.5),log(5.5/4.5),
  #                     log(0.67/0.5),log(0.13/0.09),
  #                     log(3.5/4.5)
  #                     )))*0.5
  
  iVals<-list(x0=Omega,COV=propCOV)
  return(iVals)
  
}






























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