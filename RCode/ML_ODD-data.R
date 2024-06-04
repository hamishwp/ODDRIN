# Extract Environment Variables
# dir<-directory<-"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/";setwd(directory); packred<-T
# Download and install the necessary packages:
# source('RCode/GetODDPackages.R')
library(dplyr)
library(magrittr)
library(tidyverse)
library(boot)
library(MASS)
library(pscl)
library(multiColl)
library(FactoMineR)
library(factoextra)
library(parallel)
library(doParallel)
library(caret)
library(sp)
library(tensorflow)
# tf$config$list_physical_devices("GPU")
library(keras)

source('RCode/ODDobj.R')
source('RCode/Functions.R')
# install.packages(c("ggcorrplot","vip","pdp","ggcorrplot"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE GLM MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
GetModel<-function(GLMer,weights,mvm=NULL){
  # Setup the GLM here:
  if(!is.null(mvm) & GLMer%in%c("LM","lognorm")){
    # Change between the normal (linear) model and the log-normal
    if(GLMer=="LM") {
      lognorm=F 
      fncy<-function(x) log(x+10)
    } else {
      lognorm=T
      fncy<-function(x) x
    }
    
    modely<-function(equation,datar,responsers,modout=F){
      # If we are using a lognormal model
      if(lognorm) datar[,responsers]<-log(datar[,responsers]+10)
      # if we just want the model output:
      if(modout) return(lm(formula = as.formula(equation),data = datar,weights = weights))
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      out<-do.call(rbind, lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-lm(formula = as.formula(equation),data = train,weights = weights[!1:nrow(datar)%in%indies[[i]]])
        as.data.frame(t(colMeans((abs(fncy(pmax(predict(modeler,test),-9))-fncy(test[,responsers])))*weights[indies[[i]]]/(sum(weights[indies[[i]]])*nrow(test)),na.rm = T)))
      }))
      c(colMeans(out),apply(out,2,sd))
    }
    
  } else if(GLMer=="LM"){
    modely<-function(equation,datar,responsers="Y",modout=F) {
      # Train the model on all the data in order to check the BIC
      modeler<-glm(formula = as.formula(equation),data = datar,family = gaussian(link = "identity"),weights = weights)
      # if we just want the model output:
      if(modout) return(modeler)
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      StandErr<-unlist(lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-glm(formula = as.formula(equation),data = train,family = gaussian(link = "identity"),weights = weights[!1:nrow(datar)%in%indies[[i]]])
        return(sum(abs(log(pmax(predict(modeler,test),-9)+10)-log(test$Y+10))*weights[indies[[i]]])/(sum(weights[indies[[i]]])*nrow(test)))
      }))
      
      data.frame(StandErr=mean(StandErr),
                 StandErrSD=sd(StandErr),
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="pois"){
    modely<-function(equation,datar,responsers="Y",modout=F) {
      # Train the model on all the data in order to check the BIC
      modeler<-glm(formula = as.formula(equation),data = datar,family = poisson(),weights = weights)
      # if we just want the model output:
      if(modout) return(modeler)
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      StandErr<-unlist(lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-glm(formula = as.formula(equation),data = train,family =poisson(),weights = weights[!1:nrow(datar)%in%indies[[i]]])
        return(sum(abs(log(predict(modeler,test)+10)-log(test$Y+10))*weights[indies[[i]]])/(sum(weights[indies[[i]]])*nrow(test)))
      }))
      
      data.frame(StandErr=mean(StandErr),
                 StandErrSD=sd(StandErr),
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
      
    }
    
  } else if(GLMer=="lognorm"){
    modely<-function(equation,datar,responsers="Y",modout=F) {
      # Convert to log scale
      datar$Y<-log(datar$Y+10)
      # Train the model on all the data in order to check the BIC
      modeler<-glm(formula = as.formula(equation),data = datar,family = gaussian(link = "identity"),weights = weights)
      # if we just want the model output:
      if(modout) return(modeler)
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      StandErr<-unlist(lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-glm(formula = as.formula(equation),data = train,family = gaussian(link = "identity"),weights = weights[!1:nrow(datar)%in%indies[[i]]])
        return(sum(abs(predict(modeler,test)-test$Y)*weights[indies[[i]]])/(sum(weights[indies[[i]]])*nrow(test)))
      }))
      
      data.frame(StandErr=mean(StandErr),
                 StandErrSD=sd(StandErr),
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
      
    }
    
  } else if(GLMer=="HurdlePois"){
    modely<-function(equation,datar,responsers="Y",modout=F) {
      # Train the model on all the data in order to check the BIC
      modeler<-pscl::hurdle(formula = as.formula(equation),data = datar,dist="poisson",weights = weights)
      # if we just want the model output:
      if(modout) return(modeler)
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      StandErr<-unlist(lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-pscl::hurdle(formula = as.formula(equation),data = train,dist="poisson",weights = weights[!1:nrow(datar)%in%indies[[i]]])
        return(sum(abs(log(predict(modeler,test)+10)-log(test$Y+10))*weights[indies[[i]]])/(sum(weights[indies[[i]]])*nrow(test)))
      }))
      
      data.frame(StandErr=mean(StandErr),
                 StandErrSD=sd(StandErr),
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="HurdleNegBin"){
    modely<-function(equation,datar,responsers="Y",modout=F) {
      # Train the model on all the data in order to check the BIC
      modeler<-pscl::hurdle(formula = as.formula(equation),data = datar,dist="negbin",weights = weights)
      # if we just want the model output:
      if(modout) return(modeler)
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      StandErr<-unlist(lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-pscl::hurdle(formula = as.formula(equation),data = train,dist="negbin",weights = weights[!1:nrow(datar)%in%indies[[i]]])
        return(sum(abs(log(predict(modeler,test)+10)-log(test$Y+10))*weights[indies[[i]]])/(sum(weights[indies[[i]]])*nrow(test)))
      }))
      
      data.frame(StandErr=mean(StandErr),
                 StandErrSD=sd(StandErr),
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="ZInegbin"){
    modely<-function(equation,datar,responsers="Y",modout=F) {
      # Train the model on all the data in order to check the BIC
      modeler<-pscl::zeroinfl(formula = as.formula(equation),data = datar,dist="negbin",weights = weights)
      # if we just want the model output:
      if(modout) return(modeler)
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      StandErr<-unlist(lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-pscl::zeroinfl(formula = as.formula(equation),data = train,dist="negbin",weights = weights[!1:nrow(datar)%in%indies[[i]]])
        return(sum(abs(log(predict(modeler,test)+10)-log(test$Y+10))*weights[indies[[i]]])/(sum(weights[indies[[i]]])*nrow(test)))
      }))
      
      data.frame(StandErr=mean(StandErr),
                 StandErrSD=sd(StandErr),
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="ZIpois"){
    modely<-function(equation,datar,responsers="Y",modout=F) {
      # Train the model on all the data in order to check the BIC
      modeler<-pscl::zeroinfl(formula = as.formula(equation),data = datar,dist="poisson",weights = weights)
      # if we just want the model output:
      if(modout) return(modeler)
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      StandErr<-unlist(lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-pscl::zeroinfl(formula = as.formula(equation),data = train,dist="poisson",weights = weights[!1:nrow(datar)%in%indies[[i]]])
        return(sum(abs(log(predict(modeler,test)+10)-log(test$Y+10))*weights[indies[[i]]])/(sum(weights[indies[[i]]])*nrow(test)))
      }))
      
      data.frame(StandErr=mean(StandErr),
                 StandErrSD=sd(StandErr),
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else stop("Regression model not recognised")
  
  return(modely)
  
}

# Used to correlate the vulnerability variables with the modifier term
LMFeatureSelection<-function(output,Nb=30,GLMer="LM",mvm=NULL,intercept=F,fn="+",nlim=3,weights=NULL,ncores=4){
  # Remove all the NAs
  output%<>%na.omit()
  
  if(is.null(weights)) weights<-rep(1,nrow(output))
  modely<-GetModel(GLMer,weights,mvm)
  
  # Use 80% of the observations as training and 20% for testing
  Ns <- floor(0.80*nrow(output))
  # For multivariate response models
  if(!is.null(mvm)){
    # The parts of the LM that never change
    if(!intercept) {eqn_base<-c(paste0("cbind(",paste0(mvm,collapse = ","),") ~ "),"+0") 
    } else eqn_base<-c(paste0("cbind(",paste0(mvm,collapse = ","),") ~ "),"")
    # Store the names of the MV responses
    responsers<-mvm
  } else {
    # The parts of the LM that never change
    if(!intercept) eqn_base<-c("Y ~ ","+0") else eqn_base<-c("Y ~ ","")
    responsers<-"Y"
  }
  # List of variables in the output data.frame
  vars<-colnames(output)[!colnames(output)%in%responsers]
  # Store the LM outputs
  prediction<-data.frame()
  # How many variables to use at a time?
  eqeqeq<-c()
  for(n in 1:nlim){
    eqn<-combn(vars,n)
    if(n==1) {eqn%<>%as.character()
    }else eqn<-apply(eqn,2,function(x) pracma::strcat(x,collapse = fn))
    eqeqeq%<>%c(eqeqeq,unname(sapply(eqn,function(eee) paste0(eqn_base[1],eee,eqn_base[2]))))
  }
  
  eqeqeq<-unique(eqeqeq)
  
  prediction<-mclapply(eqeqeq,mc.cores = ncores,function(eee){
    
    predictors<-data.frame()
    for (i in 1:Nb){
      # Run the model
      ressies<-tryCatch(modely(eee,output,responsers),
                        error = function(e) rep(NA,ifelse(is.null(mvm),4,length(responsers)*2L)))
      
      predictors%<>%rbind(ressies)
    }
    
    return(c(t(colMeans(predictors))))
    
  })
  # Fudging between univariate and multivariate response models
  if(is.null(mvm)) {
    prediction<-as.data.frame(t(matrix(unlist(prediction),nrow = 4)))
    colnames(prediction)<-c("StandErr","StandErrSD","AIC","BIC")
  } else {
    prediction<-as.data.frame(t(matrix(unlist(prediction),nrow = length(responsers)*2L)))
    colnames(prediction)<-c(responsers,paste0(responsers,"SD"))#[c(sapply(1:4,function(n) c(n,n+length(responsers))))]
  }
  
  return(cbind(data.frame(eqn=eqeqeq),prediction))
  
}

# Extract the data from the aggregated impact object and spatially aggregate to 0D
convODD<-function(out,ODDy,Event){
  # Make MaxHaz function over all EQ fore and aftershocks:
  if(length(names(ODDy)[grepl("hazMean",names(ODDy))])==1){
    ODDy@data$hazMax<-ODDy@data$hazMean1
  } else {
    hazard<-rep(NA_real_,length(ODDy@data$hazMean1))
    for (variable in names(ODDy)[grepl("Mean",names(ODDy))]){
      tmp<-ODDy[variable]
      tmp$hazard<-hazard
      hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
    }
    ODDy@data$hazMax<-hazard
  }
  # And for the S.D. of the hazard intensity
  if(length(names(ODDy)[grepl("hazSD",names(ODDy))])==1){
    ODDy@data$hazSD<-ODDy@data$hazSD1
  } else {
    for (variable in names(ODDy)[grepl("hazSD",names(ODDy))]){
      tmp<-ODDy[variable]
      tmp$hazard<-hazard
      hazard<-apply(tmp@data,1,function(x) mean(x,na.rm=T))
    }
    ODDy@data$hazSD<-hazard
    rm(hazard)
  }
  
  for (iso in unique(ODDy@impact$iso3)){
    
    if(!iso%in%ODDy@cIndies$iso3) next
    
    inISO<-!is.na(ODDy@data$ISO3C) & ODDy@data$ISO3C==iso
    
    tmp<-ODDy@cIndies%>%filter(iso3==iso)%>%dplyr::select(c(1,2))
    WID<-tmp$value%>%t()%>%as.data.frame(); colnames(WID)<-tmp$percentile; rm(tmp)
    
    for (ppp in unique(ODDy@impact$polygon[ODDy@impact$iso3==iso])){
      
      # extract which grid points lie within the given polygon
      inP<-rep(F,nrow(ODDy@data))
      inP[ODDy@polygons[[ppp]]$indexes] <-T
      indies<-inISO & inP
      
      # Extract impact data per polygon
      impacts<-ODDy@impact%>%filter(polygon==ppp)%>%summarise(mortality=ifelse(length(observed[impact=="mortality"]==0),sum(observed[impact=="mortality"]),NA),
                                                              displacement=ifelse(length(observed[impact=="displacement"]==0),sum(observed[impact=="displacement"]),NA),
                                                              buildDam=ifelse(length(observed[impact=="buildDam"]==0),sum(observed[impact=="buildDam"]),NA),
                                                              buildDest=ifelse(length(observed[impact=="buildDest"]==0),sum(observed[impact=="buildDest"]),NA))
      
      # Weighted (by population) vulnerability values
      vulny<-ODDy@data[inP & ODDy@data$hazMax>5,]%>%dplyr::select(c(ExpSchYrs,LifeExp,GNIc,Vs30,EQFreq,hazSD,Population))%>%
        summarise(ExpSchYrs=weighted.mean(ExpSchYrs,Population),
                  LifeExp=weighted.mean(LifeExp,Population),
                  GNIc=weighted.mean(GNIc,Population),
                  Vs30=weighted.mean(Vs30,Population),
                  EQFreq=weighted.mean(EQFreq,Population),
                  hazSD=weighted.mean(hazSD,Population))
      
      out%<>%rbind(
        cbind(impacts,data.frame(
          iso3=iso,
          date=ODDy@hazdates[1],
          Event=Event,
          Exp5=sum(ODDy@data$Population[ODDy@data$hazMax>=5 & indies],na.rm = T),
          Exp5.5=sum(ODDy@data$Population[ODDy@data$hazMax>=5.5 & indies],na.rm = T),
          Exp6=sum(ODDy@data$Population[ODDy@data$hazMax>=6 & indies],na.rm = T),
          Exp6.5=sum(ODDy@data$Population[ODDy@data$hazMax>=6.5 & indies],na.rm = T),
          Exp7=sum(ODDy@data$Population[ODDy@data$hazMax>=7 & indies],na.rm = T),
          Exp7.5=sum(ODDy@data$Population[ODDy@data$hazMax>=7.5 & indies],na.rm = T),
          Exp8=sum(ODDy@data$Population[ODDy@data$hazMax>=8 & indies],na.rm = T),
          Exp8.5=sum(ODDy@data$Population[ODDy@data$hazMax>=8.5 & indies],na.rm = T),
          Exp9=sum(ODDy@data$Population[ODDy@data$hazMax>=9 & indies],na.rm = T),
          maxHaz=max(ODDy@data$hazMax[indies],na.rm = T)
        ),WID,vulny)
      )
    }
  }
  
  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRACT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
folder<-"./IIDIPUS_Input/IIDIPUS_Input_All_2023May19/ODDobjects/"
filez<-list.files(folder)

out<-data.frame()
for(fff in filez){
  # Exceptions... sigh
  if(fff=="EQ20170529IDN_136") next
  # Read in the hazard & impact object
  ODDy<-readRDS(paste0(folder,fff))
  # Convert it into the necessary form
  out%<>%convODD(ODDy,fff)
  
  print(paste0("Finished EQ: ",fff))
}

out$time<-as.numeric(out$date-min(out$date))

saveRDS(out,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/InputData_ODD.RData")

out<-readRDS("./IIDIPUS_Results/SpatialPolygons_ML-GLM/InputData_ODD.RData")

# avHS<-read.csv("~/Downloads/GDL-Population-(2021)-data.csv")%>%
#   filter(Level=="National")%>%dplyr::select(ISO_Code,HH.size)
# 
# colnames(avHS)[1]<-"iso3"
# 
# out%<>%left_join(avHS)
# 
# out$buildDisp<-rowSums(cbind(out$buildDam,out$buildDest),na.rm = T)*out$HH.size
# 
# out$buildDisp[is.na(out$buildDam) & is.na(out$buildDest)]<-NA
# 
# 
# nnn<-sum(!is.na(out$buildDisp) & !is.na(out$displacement))
# rsq<-summary(lm(log(displacement+10) ~ log(buildDisp+10) + 0,out))$adj.r.squared
# 
# p<-out%>%ggplot(aes(displacement, buildDisp))+geom_point()+
#   scale_x_log10(limits=c(100,1.3e6))+scale_y_log10(limits=c(100,1.3e6)) + 
#   geom_abline(slope = 1,intercept = 0)+
#   xlab("Displacement")+ylab("Building Damage x Av. Household Size")+
#   annotate(geom="text", x=8e2, y=5e5, size=6,
#            label=paste0("Adj-R-sq = ",signif(rsq,2)));p
#   
# ggsave("DispHH-Size.eps",p,path="./Plots/IIDIPUS_Results/",width=6,height=5,device = grDevices::cairo_ps)  

scaleIMPs<-function(out){
  # Per impact, do the analysis (for all impacts, even without building data)
  outred<-dplyr::select(out,-c("iso3","date")); rm(out)
  outred[,-(1:5)]<-scale(outred[,-(1:5)])
  # This reveals 95% of variance is contained within 4 dimensions of the PCA
  ExpPCA<-outred%>%dplyr::select(colnames(outred)[str_starts(colnames(outred),"Exp") & colnames(outred)!="ExpSchYrs"])%>%FactoMineR::PCA(ncp = 5,graph = F)
  iExp<-sum(ExpPCA$eig[,3]<95); iExp<-2
  tmpExp<-as.data.frame(ExpPCA$ind$coord[,1:iExp]); colnames(tmpExp)<-paste0("ExpDim",1:iExp)
  outred%<>%cbind(tmpExp)%>%dplyr::select(-colnames(outred)[str_starts(colnames(outred),"Exp") & colnames(outred)!="ExpSchYrs"])
  rm(ExpPCA,iExp,tmpExp)
  # This reveals 95% of variance is contained within 2 dimensions of the PCA
  WIDPCA<-outred%>%dplyr::select(colnames(outred)[str_starts(colnames(outred),"p")])%>%FactoMineR::PCA(ncp = 5,graph = F)
  iWID<-sum(WIDPCA$eig[,3]<95)
  tmpWID<-as.data.frame(WIDPCA$ind$coord[,1:iWID]); colnames(tmpWID)<-paste0("WIDDim",1:iWID)
  outred%<>%cbind(tmpWID)%>%dplyr::select(-colnames(outred)[str_starts(colnames(outred),"p")])
  rm(WIDPCA,iWID,tmpWID)
  # Remove NAs from the object
  outred<-outred[apply(outred[,-(1:4)],1,function(x) !any(is.na(x))),]
  # Check the Variance Inflation Factor
  print(paste0("Checking the VIF of the remaining variables: percentage that didn't pass = ",sum(multiColl::VIF(outred[,-(1:5)])>6)/ncol(outred[,-(1:5)])))
  
  return(outred)
}

outred<-scaleIMPs(out); rm(out)
covariates<-colnames(outred)[-c(1:5)]
allimps<-c("mortality","displacement","buildDam","buildDest")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNIVARIATE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

ExtractGLMresults<-function(algo,impact,othimps=NULL){
  
  print(paste0("Working on ",algo," for impact=",impact))  
  
  if(is.null(othimps)){
    outFrame<-dplyr::select(outred,-allimps[impact!=allimps])  
  } else {
    outFrame<-dplyr::select(outred,-allimps[!allimps%in%c(impact,othimps)])
    ind<-which(colnames(outFrame)==impact)
    outFrame<-outFrame[,c(ind,(1:ncol(outFrame))[-ind])]
    
    for(i in 1:length(othimps)) outFrame[,othimps[i]]<-log(outFrame[,othimps[i]]+10)
  }
  
  names(outFrame)[1]<-"Y"
  outFrame%<>%na.omit()
  # Make weights from the different events to make sure that no single event dominates the model parameterisation
  weights<-outFrame%>%group_by(Event)%>%summarise(www=1/length(time))%>%merge(outFrame)%>%pull(www)
  # Remove the variable Event after weighting is calculated
  outFrame%<>%dplyr::select(-Event)
  
  # Function from the file CorrelateModifier.R:
  out<-tryCatch(LMFeatureSelection(outFrame,Nb=1,intercept=T,fn="+",nlim=(ncol(outFrame)-1),
                                   GLMer = algo, weights = weights, ncores=60),
                error=function(e) NA)
  # out$model<-algo
  saveRDS(out,paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_",algo,"_",impact,ifelse(is.null(othimps),"",paste0("_with-",othimps)),".RData"))
  
  return(out)
  
}

minimods<-c("LM","pois","lognorm","HurdlePois","HurdleNegBin","ZIpois","ZInegbin")

# Run it!
ODD_ML<-lapply(allimps, function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractGLMresults(algo,impact),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
  })})

# Analyse the results:
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/"); filez<-filez[!grepl(filez,pattern = "_with")]
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"GLM_",simplify = T)[,2]; namerz<-namerz[!namerz%in%c("Mortality","BuildDam","BuildDest","Displacement","")]

predictions<-data.frame()
for(i in 1:length(filez)) {
  if(filez[i]=="InputDataGLM.RData" | grepl(filez[i],pattern = "MVGLM")) next
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/",filez[i]))
  if(all(is.na(tmp))) {print(paste0("Failed for ",filez[i]));next}
  if("model"%in%colnames(tmp)) tmp%<>%dplyr::select(-"model")
  predictions%<>%rbind(cbind(tmp,data.frame(model=namerz[i])))
}

tmp<-str_split(predictions$model,"_",simplify = T)
predictions$impact<-tmp[,2]
predictions$algo<-tmp[,1]
predictions$model<-NULL
table(predictions$impact)
table(predictions$algo)

predictions%>%arrange(StandErr)%>%group_by(impact,algo)%>%slice(1)%>%View()
# predictions%>%arrange(BIC)%>%group_by(impact)%>%slice(1:5)
# 
# predictions%>%group_by(impact,algo)%>%
#   summarise(StandErr=min(StandErr),BIC=BIC[which.min(StandErr)],StandErrSD=StandErrSD[which.min(StandErr)])%>%
#   ggplot(aes(StandErr,BIC))+geom_point(aes(colour=algo,shape=impact),size=3)+
#   geom_errorbar(aes(xmin=StandErr-StandErrSD, 
#                                       xmax=StandErr+StandErrSD,colour=algo))
#   scale_y_log10()+scale_x_log10()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%% CONDITIONAL UNIVARIATE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

minimods<-c("lognorm")

# Run it!
ODD_ML<-lapply(allimps[2:4], function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractGLMresults(algo,impact,othimps = "mortality"),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
  })})

ODD_ML<-lapply(allimps[c(1,3,4)], function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractGLMresults(algo,impact,othimps = "displacement"),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
  })})

ODD_ML<-lapply(allimps[c(1,2,4)], function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractGLMresults(algo,impact,othimps = "buildDam"),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
  })})

ODD_ML<-lapply(allimps[c(1,2,3)], function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractGLMresults(algo,impact,othimps = "buildDest"),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
  })})

ODD_ML<-lapply(allimps[3:4], function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractGLMresults(algo,impact,othimps = c("mortality","displacement")),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
  })})

ODD_ML<-lapply(allimps[1:2], function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractGLMresults(algo,impact,othimps = c("buildDam","buildDest")),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
  })})

# Find the best equation, then compare the error with the best GLM without the impact terms
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/"); filez<-filez[grepl(filez,pattern = "_with")]
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"GLM_",simplify = T)[,2]; namerz<-namerz[!namerz%in%c("Mortality","BuildDam","BuildDest","Displacement","")]

predictionsCUV<-data.frame()
for(i in 1:length(filez)) {
  if(filez[i]=="InputDataGLM.RData" | grepl(filez[i],pattern = "MVGLM")) next
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/",filez[i]))
  if(all(is.na(tmp))) {print(paste0("Failed for ",filez[i]));next}
  if("model"%in%colnames(tmp)) tmp%<>%dplyr::select(-"model")
  predictionsCUV%<>%rbind(cbind(tmp,data.frame(model=namerz[i])))
}

predictionsCUV$otherimp<-str_split(predictionsCUV$model,"-",simplify = T)[,2]
predictionsCUV$impact<-str_split(predictionsCUV$model,"_",simplify = T)[,2]
predictionsCUV$impineq<-unlist(lapply(1:nrow(predictionsCUV),function(i) grepl(predictionsCUV$otherimp[i],predictionsCUV$eqn[i])))

predictionsCUV%<>%filter(otherimp%in%allimps)

predictionsCUV%<>%na.omit()

condout<-predictionsCUV%>%group_by(impact,otherimp)%>%
  summarise(bestCond=min(StandErr[impineq],na.rm = T),
            bestCondSD=StandErrSD[which.min(StandErr[impineq])],
            bestNoCond=min(StandErr[!impineq],na.rm = T),
            bestNoCondSD=StandErrSD[which.min(StandErr[!impineq])])
condout$improvedPred<-condout$bestCond<condout$bestNoCond

condout%<>%group_by(impact,otherimp)%>%mutate(nn=nrow(na.omit(outred[,c(unique(impact),unique(otherimp))])))

condout%<>%mutate(minSD=bestCondSD^2/unique(nn),
                                      thisSD=bestNoCondSD^2/unique(nn))%>%
  mutate(df=(nn-1)*(thisSD+minSD)^2/(thisSD^2+minSD^2))%>%dplyr::select(-c(minSD,thisSD))

condout%<>%mutate(tval=sqrt(nn)*(bestCond-bestNoCond)/
                                        sqrt(bestCondSD^2+bestNoCondSD^2))
condout%<>%mutate(pval=dt(tval,df))%>%dplyr::select(-c(df,nn,tval,bestCondSD,bestNoCondSD))

condout$improvedPred<-condout$improvedPred & condout$pval<0.05

condtab<-xtable::xtable(condout,digits = 3,
               "Increase in model performance by including other impact types as covariates.")
names(condtab)<-c("Response Impact","Covariate Impact","Best MADL w Covariate","Best MADL w/o Covariate","Improved? (Stat. Sig. Only)","P-value [95%]")

print(condtab,include.rownames=F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

costie<-function(data, lev = NULL, model = NULL) {
  out<-mean(abs(data$obs-data$pred)*data$weights/sum(data$weights),na.rm = T)
  names(out)<-"RelativeAbs"
  return(out)
}
  
train_control <- caret::trainControl(method="repeatedcv", number=10, repeats=3,
                                     search = "random",
                                     summaryFunction=costie)

allimps<-c("mortality","displacement","buildDam","buildDest")

parallelML<-function(algo,impact,modRet=F,othimps=NULL) {
  
  print(paste0("Working on ",algo," for impact=",impact))  
  
  if(is.null(othimps)){
    outFrame<-dplyr::select(outred,-allimps[impact!=allimps])  
  } else {
    outFrame<-dplyr::select(outred,-allimps[!allimps%in%c(impact,othimps)])
    ind<-which(colnames(outFrame)==impact)
    outFrame<-outFrame[,c(ind,(1:ncol(outFrame))[-ind])]
  }
  
  names(outFrame)[1]<-"Y"
  outFrame$Y<-log(outFrame$Y+10)
  outFrame%<>%na.omit()
  # Make weights from the different events to make sure that no single event dominates the model parameterisation
  weights<-outFrame%>%group_by(Event)%>%summarise(www=1/length(time))%>%merge(outFrame)%>%pull(www)
  # Remove the variable Event after weighting is calculated
  outFrame%<>%dplyr::select(-Event)
  
  modeler<-caret::train(Y~., data = outFrame, method = algo, metric="RelativeAbs",
                        tuneLength = 12, trControl = train_control,linout = TRUE,
                        weights = weights, preProcess = c("center","scale"))
  
  print(modeler$results)
  
  if(modRet) return(modeler)
  
  out<-cbind(dplyr::select(filter(modeler$results,RelativeAbs==min(RelativeAbs)),RelativeAbs,RelativeAbsSD),
          t(as.data.frame((t(as.data.frame(varImp(modeler, useModel=F, nonpara=F, scale=FALSE)$importance))[1,]))))
  rownames(out)<-NULL
  
  saveRDS(out,paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/ML_",algo,"_",impact,ifelse(is.null(othimps),"",paste0("_with-",othimps)),".RData"))
  
  return(out)
  
}

# tabmod<-getModelInfo()
# carmods<-unlist(sapply(tabmod,function(x) x$type%in%"Regression"))
# carmods<-data.frame(algorithm=names(carmods),regression=unname(carmods))
# carmods%<>%filter(regression)%>%pull(algorithm)
# # Check that we have all that we need to run each model
# checkerz<-unlist(lapply(carmods,function(stst) ifelse(is.null(tryCatch(checkInstall(getModelInfo(stst)$library),error=function(e) NA)),T,F)))
# carmods<-carmods[checkerz]; rm(checkerz)

minimods<-c("nnet","brnn","rf","glmnet","lm","lmStepAIC","svmLinear","svmRadial","svmPoly")
            # "blasso","bridge","glm.nb","lasso","nnls","pcr")
# ,"gbm"

# Parallelise
ncores<-60
cl <- makePSOCKcluster(ncores)  # Create computing clusters
registerDoParallel(cl)
getDoParWorkers()
# Run it!
ODD_ML<-lapply(allimps, function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(parallelML(algo,impact),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
    })})

# Remember to close the computing cluster
stopCluster(cl)
registerDoSEQ()

# CONDITIONAL MODELS

cl <- makePSOCKcluster(ncores)  # Create computing clusters
registerDoParallel(cl)
getDoParWorkers()

ODD_ML_Cond<-lapply(allimps, function(impact) {
  lapply(allimps[allimps!=impact],function(oth){
    lapply(minimods,function(algo) {
      out<-tryCatch(parallelML(algo,impact,othimp = oth),error=function(e) NA)
      if(any(is.na(out))) return(NULL)
    })
  })
})

# Remember to close the computing cluster
stopCluster(cl)
registerDoSEQ()

# Let's have a look! :)
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/"); filez<-filez[!grepl(filez,pattern = "_with")]
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"ML_",simplify = T)[,2]
impact<-str_split(namerz,"_",simplify = T)[,2]
namerz<-str_split(namerz,"_",simplify = T)[,1]

predictionsML<-data.frame()
for(i in 1:length(filez)) {
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/",filez[i]))
  predictionsML%<>%rbind(cbind(data.frame(model=namerz[i],impact=impact[i]),tmp))
}
predictionsML[,covariates]<-100*predictionsML[,covariates]/rowSums(predictionsML[,covariates])

table(predictionsML$model)
table(predictionsML$impact)

View(predictionsML)

# Find the best equation, then compare the error with the best GLM without the impact terms
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/"); filez<-filez[grepl(filez,pattern = "_with")]
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"ML_",simplify = T)[,2]; namerz<-namerz[!namerz%in%c("Mortality","BuildDam","BuildDest","Displacement","")]

predictionsCUVML<-data.frame()
for(i in 1:length(filez)) {
  if(filez[i]=="InputDataGLM.RData" | grepl(filez[i],pattern = "MVGLM")) next
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/",filez[i]))
  if(all(is.na(tmp))) {print(paste0("Failed for ",filez[i]));next}
  if("model"%in%colnames(tmp)) tmp%<>%dplyr::select(-"model")
  predictionsCUVML%<>%rbind(cbind(tmp[,1:2],data.frame(model=namerz[i])))
}

predictionsCUVML$otherimp<-str_split(predictionsCUVML$model,"-",simplify = T)[,2]
predictionsCUVML$impact<-str_split(predictionsCUVML$model,"_",simplify = T)[,2]
predictionsCUVML$model<-str_split(predictionsCUVML$model,"_",simplify = T)[,1]

predictionsCUVML%<>%filter(otherimp%in%allimps)

predictionsCUVML%<>%na.omit()%>%dplyr::select(impact,otherimp,model,everything())
colnames(predictionsCUVML)[4:5]<-c("bestCond","bestCondSD")

tmp<-predictionsML%>%dplyr::select(model,impact,RelativeAbs,RelativeAbsSD)
colnames(tmp)[3:4]<-c("bestNoCond","bestNoCondSD")

predictionsCUVML%<>%full_join(tmp)#,by=c("model","impact"))

condout<-predictionsCUV%>%group_by(impact,otherimp)%>%
  summarise(bestCond=min(StandErr[impineq],na.rm = T),
            bestCondSD=StandErrSD[which.min(StandErr[impineq])],
            bestNoCond=min(StandErr[!impineq],na.rm = T),
            bestNoCondSD=StandErrSD[which.min(StandErr[!impineq])])

condout%<>%rbind(dplyr::select(predictionsCUVML,-model))

condout%<>%group_by(impact,otherimp)%>%
  summarise(bestCond=min(bestCond,na.rm = T),
            bestCondSD=bestCondSD[which.min(bestCond)],
            bestNoCond=min(bestNoCond,na.rm = T),
            bestNoCondSD=bestNoCondSD[which.min(bestNoCond)])

condout$improvedPred<-condout$bestCond<condout$bestNoCond

condout%<>%group_by(impact,otherimp)%>%mutate(nn=nrow(na.omit(outred[,c(unique(impact),unique(otherimp))])))

condout%<>%mutate(minSD=bestCondSD^2/unique(nn),
                  thisSD=bestNoCondSD^2/unique(nn))%>%
  mutate(df=(nn-1)*(thisSD+minSD)^2/(thisSD^2+minSD^2))%>%dplyr::select(-c(minSD,thisSD))

condout%<>%mutate(tval=sqrt(nn)*(bestCond-bestNoCond)/
                    sqrt(bestCondSD^2+bestNoCondSD^2))
condout%<>%mutate(pval=dt(tval,df))%>%dplyr::select(-c(df,nn,tval,bestCondSD,bestNoCondSD))

condout$improvedPred<-condout$improvedPred & condout$pval<0.05

condtab<-xtable::xtable(condout,
                        "Increase in model performance by including other impact types as covariates.")
names(condtab)<-c("Response Impact","Covariate Impact","Best MADL w Covariate","Best MADL w/o Covariate","Improved? (Stat. Sig. Only)","P-value [95%]")

print(condtab,include.rownames=F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%% COMPARE UNIVARIATE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

fuller<-rbind(transmute(predictions,RelativeAbsDiff=StandErr,RelativeAbsDiffSD=StandErrSD,
                        impact=impact,algo=algo, equation=eqn),
              transmute(predictionsML,RelativeAbsDiff=RelativeAbs,RelativeAbsDiffSD=RelativeAbsSD,
                        impact=impact,algo=model,equation="Y ~ maxHaz + ExpSchYrs + LifeExp + GNIc + Vs30 + EQFreq + hazSD + time + ExpDim1 + ExpDim2 + WIDDim1 + WIDDim2"))
fuller%<>%na.omit()

# For the top-five models in terms of performance across all impacts
fuller%>%group_by(impact,algo)%>%filter(RelativeAbsDiff==min(RelativeAbsDiff))%>%
  ungroup()%>%group_by(algo)%>%summarise(Cost=prod(RelativeAbsDiff))%>%
  ungroup()%>%mutate(Cost=Cost/min(Cost))%>%arrange(Cost)



View(fuller)     

predictions%<>%group_by(impact)%>%mutate(BestDiff=StandErr-min(StandErr,na.rm = T))

p<-predictions%>%filter(StandErr<0.6)%>%ggplot(aes(StandErr,group=algo)) +
  geom_density(aes(fill=algo),alpha=0.4)+scale_x_log10()+
  xlab("Relative Absolute Log-Difference Errors - Log Scale ")+ylab("Density")+labs(fill="Model")+
  facet_wrap(.~impact,scales = "free");p
ggsave("GLM_Errors.eps",p,path="./Plots/IIDIPUS_Results/GLM-ML_Work/",width=10,height=5,device = grDevices::cairo_ps)  


# colnames(allUV)<-c("Impact","Model","Avg. MADL Value", "S.D. MADL Value")

# For each of the top models, calculate the feature importance using vip package
# See here for more info: https://cran.r-project.org/web/packages/vip/vignettes/vip-introduction.pdf
# Then compare the values between the top 10-models (over all model types) for each impact

# FOR ALL MODELS NOT STAT. SIG. DIFF FROM BEST MODEL, PER GLM, MEASURE THE MODEL-AGNOSTIC VIP

topOeach<-fuller%>%group_by(impact,algo)%>%arrange(RelativeAbsDiff)%>%dplyr::select(3,4,1,2,5)
topOeach%<>%group_by(impact)%>%mutate(nn=sum(!is.na(outred[,unique(impact)])))

topOeach%<>%group_by(impact)%>%mutate(minSD=RelativeAbsDiffSD[which.min(RelativeAbsDiff)]^2/unique(nn),
                                      thisSD=RelativeAbsDiffSD^2/unique(nn))%>%
  mutate(df=(nn-1)*(thisSD+minSD)^2/(thisSD^2+minSD^2))

topOeach%<>%group_by(impact)%>%mutate(tval=sqrt(nn)*(RelativeAbsDiff-min(RelativeAbsDiff))/
                                        sqrt(RelativeAbsDiffSD^2+RelativeAbsDiffSD[which.min(RelativeAbsDiffSD)]^2))
topOeach%<>%mutate(pval=dt(tval,df))

topOeach%<>%dplyr::select(-c(nn,minSD,thisSD,df,tval))

topOeach%<>%filter(pval>0.05)

table(topOeach$impact[topOeach$pval>0.05])
table(topOeach$algo[topOeach$pval>0.05])

GLMmods<-predictions%>%filter(algo%in%unique(topOeach$algo))%>%group_by(impact,algo)%>%
  filter(eqn%in%topOeach$equation[topOeach$algo==unique(algo) & topOeach$impact==unique(impact)])

colnames(topOeach)[5]<-"eqn"

GLMmods%<>%merge(dplyr::select(topOeach,impact,algo,eqn,pval))

GLMvarImp<-function(GLMmods,tester=NULL,predictionsML=NULL){
  # line by line model development and variable information extraction
  resultsUV<-do.call(rbind,lapply(1:nrow(GLMmods),function(i){
    inpy<-as.data.frame(GLMmods[i,])
    # Sort out the data
    outFrame<-dplyr::select(outred,-allimps[inpy$impact!=allimps])
    names(outFrame)[1]<-"Y"
    outFrame%<>%na.omit()
    # Make weights from the different events to make sure that no single event dominates the model parameterisation
    weights<-outFrame%>%group_by(Event)%>%summarise(www=1/length(time))%>%merge(outFrame)%>%pull(www)
    # Remove the variable Event after weighting is calculated
    outFrame%<>%dplyr::select(-Event)
    # Find the GLM model
    modeler<-GetModel(inpy$algo,weights)
    # Extract the trained model
    outmod<-modeler(inpy$eqn,outFrame,"Y",modout=T)
    # In case we don't want feature importance but the performance on a separate validation dataset
    if(!is.null(tester)){
      return(data.frame(pred=exp(predict(outmod,tester))+10,
                        true=tester[[inpy$impact]],
                        impact=inpy$impact,model=inpy$algo))
    }
    # Which variables to include
    vars<-names(outmod$model)[-1];vars<-vars[-length(vars)]
    # Feature importance (model-agnostic) calculation
    vippy<-as.data.frame(vip::vi(outmod,scale=T,ice=T,
                                 feature_names=vars))%>%dplyr::select(Variable,Importance)
    # Add missing columns
    missies<-colnames(outFrame)[-1][!colnames(outFrame)[-1]%in%vippy$Variable]
    if(length(missies)>0) vippy%<>%rbind(data.frame(Variable=missies,Importance=0))
    # Reformulate to column form 
    rownames(vippy)<-vippy$Variable; vippy%<>%dplyr::select(Importance)%>%t()%>%as.data.frame()
    # colnames in alphabetical order and add the model information to the data frame
    vippy%<>%dplyr::select(sort(colnames(vippy)))%>%cbind(inpy)%>%
      dplyr::select(impact,algo,StandErr,StandErrSD,everything())
    # output it all!
    return(vippy)
  }))
  
  if(!is.null(tester)) return(resultsUV)
  # Reorder the dataframe
  resultsUV%<>%dplyr::select(algo,impact,StandErr,StandErrSD,everything())%>%dplyr::select(-c(AIC,BIC,pval,eqn))
  colnames(resultsUV)[1:4]<-colnames(predictionsML)[1:4]
  rownames(resultsUV)<-NULL
  # Scale to percentage the feature importance
  resultsUV[,covariates]<-100*resultsUV[,covariates]/rowSums(resultsUV[,covariates])
  # Potentially add the ML results as well
  if(!is.null(predictionsML)) resultsUV%<>%rbind(predictionsML)
  # Check for statistical significance.
  # Sample size
  resultsUV%<>%group_by(impact)%>%mutate(nn=sum(!is.na(outred[,unique(impact)])))
  # degrees of freedom
  resultsUV%<>%group_by(impact)%>%mutate(minSD=RelativeAbsSD[which.min(RelativeAbs)]^2/unique(nn),
                                         thisSD=RelativeAbsSD^2/unique(nn))%>%
    mutate(df=(nn-1)*(thisSD+minSD)^2/(thisSD^2+minSD^2))
  # t-value from t-distribution
  resultsUV%<>%group_by(impact)%>%mutate(tval=sqrt(nn)*(RelativeAbs-min(RelativeAbs))/
                                           sqrt(RelativeAbsSD^2+RelativeAbsSD[which.min(RelativeAbsSD)]^2))
  # p-value associated to the t-values
  resultsUV%<>%mutate(pval=dt(tval,df))
  
  return(resultsUV)
}

resultsUV<-GLMvarImp(GLMmods)

resultsUV%<>%filter(pval>0.05)

# Order the variables by their importance, weighted by each models error value
varimport<-resultsUV%>%group_by(impact)%>%
  reframe(across(all_of(c("RelativeAbs","RelativeAbsSD",covariates)), ~ weighted.mean(.x,pval)))%>%
  dplyr::select(-c(RelativeAbs,RelativeAbsSD))

# To account for the variance in the estimates
wavg<-(varimport[varimport$impact=="mortality",-1] +
  varimport[varimport$impact=="displacement",-1]/3.22 + 
  varimport[varimport$impact=="buildDam",-1]/3.06 +
  varimport[varimport$impact=="buildDest",-1]/4)/(1+1/3.22+1/3.06+1/4)

varimport%<>%rbind(cbind(data.frame(impact="average"),wavg))
                         
varimport[,-1]<-100*varimport[,-1]/rowSums(as.matrix(varimport[,-1]))

ordy<-colnames(varimport[,-1])[rev(order(as.numeric(varimport[varimport$impact=="average",-1])))]

varimport%<>%reshape2::melt("impact")

varimport$variable%<>%factor(levels=ordy)

pal <- c(
  "mortality" = scales::hue_pal()(4)[1],
  "displacement" = scales::hue_pal()(4)[2], 
  "buildDam" = scales::hue_pal()(4)[3], 
  "buildDest" = scales::hue_pal()(4)[4],
  "average" = "black"
)

saveRDS(varimport,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/varimport.RData")

p<-varimport%>%
  ggplot(aes(variable,value,group=impact))+
  geom_point(aes(colour=impact,shape=impact),size=3)+
  geom_line(aes(colour=impact),linewidth=0.7, alpha=0.5,linetype="dotdash")+
  scale_shape_manual(values=c(15:18,1),breaks=c(allimps,"average"))+
  scale_colour_manual(values = pal,limits = names(pal))+
  xlab("Model Covariate") + ylab("Feature Importance (%)")+
  labs(colour="Impact Type",shape="Impact Type")+
  theme(axis.text.x = element_text(angle = 45, hjust=1));p

ggsave("UV_FeatImp.eps",p,path="./Plots/IIDIPUS_Results/",width=9,height=4.,device = grDevices::cairo_ps)  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MULTIVARIATE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

ExtractMVGLMresults<-function(algo,impact){
  
  print(paste0("Working on ",algo," for impacts=",paste0(impact,collapse = " & ")))
  # Get rid of NAs
  outFrame<-outred%>%na.omit()%>%dplyr::select(-allimps[!allimps%in%impact])
  # Make weights from the different events to make sure that no single event dominates the model parameterisation
  weights<-outFrame%>%group_by(Event)%>%summarise(www=1/length(time))%>%merge(outFrame)%>%pull(www)
  # Remove the variable Event after weighting is calculated
  outFrame%<>%dplyr::select(-Event)
  # Function from the file CorrelateModifier.R:
  out<-LMFeatureSelection(outFrame,Nb=1,intercept=T,fn="+",nlim=(ncol(outFrame)-length(impact)),
                          GLMer = algo, weights = weights, mvm = impact, ncores=60)
  # out<-tryCatch(LMFeatureSelection(outFrame,Nb=1,intercept=T,fn="+",nlim=(ncol(outFrame)-1),
  #                                  GLMer = algo, weights = weights, mvm = impact, ncores=60),
  #               error=function(e) NA)
  out$model<-algo
  
  saveRDS(out,paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/MVGLM_",algo,"_",paste0(impact,collapse = ""),".RData"))
  
  return(out)
  
}

minimods<-c("LM","lognorm")

impcombs<-list(c("mortality","displacement"),
               c("buildDam","buildDest"),
               c("mortality","buildDest"),
               c("mortality","buildDam"),
               c("displacement","buildDest"),
               c("displacement","buildDam"),
               c("mortality","displacement","buildDam","buildDest"))

# Run it!
ODD_ML<-lapply(1:length(impcombs), function(jj) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractMVGLMresults(algo,unlist(impcombs[[jj]])),error=function(e) NA)
    if(any(is.na(out))) {print(paste0("Fail for ",algo,unlist(impcombs[[jj]]))); return(NULL)}
  })})

# Let's have a look! :)
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/")
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"GLM_",simplify = T)[,2]; namerz<-namerz[namerz!=""]

predictionsMV<-data.frame()
for(i in 1:length(filez)) {
  if(filez[i]=="InputData_BD.RData") next
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/",filez[i]))
  if("model"%in%colnames(tmp)) tmp%<>%dplyr::select(-"model")
  missies<-allimps[!allimps%in%colnames(tmp)]
  if(length(missies)>0) tmp[,c(missies,paste0(allimps[!allimps%in%colnames(tmp)],"SD"))]<-NA
  tmp<-tmp[apply(tmp,1,function(x)!all(is.na(x))),]
  predictionsMV%<>%rbind(cbind(tmp,data.frame(model=namerz[i])))
}
predictionsMV$Cost<-apply(predictionsMV[,allimps],1,prod,na.rm=T)

tmp<-str_split(predictionsMV$model,"_",simplify = T)
predictionsMV$impact<-tmp[,2]
predictionsMV$algo<-tmp[,1]
predictionsMV$model<-NULL

predictionsMV%>%arrange(Cost)%>%dplyr::select(-c(allimps,paste0(allimps,"SD")))%>%
  group_by(algo,impact)%>%slice(1:5)%>%View()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%% COMPARE MULTIVARIATE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

predictionsMV%>%arrange(Cost)%>%group_by(impact,algo)%>%slice(1)

MVGLMmods<-predictionsMV%>%arrange(Cost)%>%group_by(impact,algo)%>%slice(1)
MVvars<-str_split(str_split(str_split((str_split(unique(MVGLMmods$eqn),"~",simplify = T)[,1]),"\\(",simplify = T)[,2],"\\)",simplify = T)[,1],",")

MVGLMmods%<>%filter(algo=="lognorm")

resultsMV<-lapply(1:nrow(MVGLMmods),function(i){
  inpy<-as.data.frame(MVGLMmods[i,])
  if(inpy$algo=="LM") return(NA)
  # Remove NAs
  outFrame<-outred%>%na.omit()%>%dplyr::select(-allimps[!allimps%in%MVvars[[i]]])
  # Make sure log-normal model is taken into account
  if(inpy$algo=="lognorm") outFrame[,MVvars[[i]]]<-log(outFrame[,MVvars[[i]]]+10)
  # Remove the variable Event after weighting is calculated
  outFrame%<>%dplyr::select(-Event)
  # Have the full equation
  eqn<-paste0("cbind(",paste0(MVvars[[i]],collapse=","),") ~ ",paste0(covariates,collapse = " + "))
  # Extract the trained model
  outmod<-lm(formula = as.formula(eqn),data = outFrame)
  # Feature importance (model-agnostic) calculation
  vippy<-as.data.frame(vip::vi(outmod,method="firm",scale=T,ice=T,
                               feature_names=covariates))
  # Add missing columns
  missies<-colnames(outFrame)[!colnames(outFrame)%in%c(vippy$Variable,MVvars[[i]])]
  vippy%<>%rbind(data.frame(Variable=missies,Importance=0))
  # Reformulate to column form 
  rownames(vippy)<-vippy$Variable; vippy%<>%dplyr::select(Importance)%>%t()%>%as.data.frame()
  # colnames in alphabetical order and add the model information to the data frame
  vippy%<>%dplyr::select(sort(colnames(vippy)))%>%cbind(dplyr::select(inpy,Cost,impact,algo))
  # output it all!
  return(list(vip=vippy,covvy=vcov(outmod)))
})

# Displacement and Mortality
covvy<-resultsMV[[1]]$covvy
ind<-which(grepl("Intercept",colnames(covvy)))
covvy<-covvy[-ind,]; covvy<-covvy[,-ind]

covvy<-covvy[(1:(ncol(covvy)/2)),
      ((ncol(covvy)/2)+1):ncol(covvy)]

covvy<-covvy/max(abs(covvy))

ggcorrplot::ggcorrplot(covvy, 
                       type = "lower",show.diag = T,
                       lab = TRUE)  

# Building damage & destruction
covvy<-resultsMV[[2]]$covvy
ind<-which(grepl("Intercept",colnames(covvy)))
covvy<-covvy[-ind,]; covvy<-covvy[,-ind]

covvy<-covvy[(1:(ncol(covvy)/2)),
             ((ncol(covvy)/2)+1):ncol(covvy)]

covvy<-covvy/max(abs(covvy))

ggcorrplot::ggcorrplot(covvy, 
                       type = "lower",show.diag = T,
                       lab = TRUE)  

outer<-as.data.frame(do.call(rbind,lapply(1:(length(MVvars)-1),function(i){
  eqn<-paste0("cbind(",paste0(MVvars[[i]],collapse=","),") ~ ",paste0(covariates,collapse = " + "))
  covvy<-vcov(lm(formula = as.formula(eqn),data = outred))
  
  max(covvy)
  
  # covvy<-resultsMV[[i]]$covvy
  ind<-which(grepl("Intercept",colnames(covvy)))
  covvy<-covvy[-ind,]; covvy<-covvy[,-ind]
  
  covvy<-covvy[(1:(ncol(covvy)/2)),
               ((ncol(covvy)/2)+1):ncol(covvy)]
  
  diag(covvy/max(abs(covvy)))

})))

colnames(outer)<-covariates
outer<-cbind(data.frame(impacts=unlist(lapply(1:(length(MVvars)-1),function(i) paste0(MVvars[[i]],collapse="-")))),outer)

outer%<>%reshape2::melt("impacts")

outer$variable%<>%factor(levels=ordy)

pal <- c(
  "mortality" = scales::hue_pal()(4)[1],
  "displacement" = scales::hue_pal()(4)[2], 
  "buildDam" = scales::hue_pal()(4)[3], 
  "buildDest" = scales::hue_pal()(4)[4]
)

p<-outer%>%
  ggplot(aes(variable,value,group=impacts))+
  geom_point(aes(colour=impacts,shape=impacts),size=3)+
  geom_line(aes(colour=impacts),linewidth=0.7, alpha=0.5,linetype="dotdash")+
  scale_shape_manual(values=c(15:18,1),breaks=c(allimps,"average"))+
  scale_colour_manual(values = pal,limits = names(pal))+
  xlab("Model Covariate") + ylab("Inter-impact Parameter Variance [Normalised]")+
  labs(colour="Impact Interactions",shape="Impact Interactions")+
  theme(axis.text.x = element_text(angle = 45, hjust=1));p

ggsave("MV_FeatCor.eps",p,path="./Plots/IIDIPUS_Results/",width=9,height=4.,device = grDevices::cairo_ps)  



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%% CONVOLUTIONAL NEURAL NETWORKS %%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%% EXTRACT DATA %%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
folder<-"./IIDIPUS_Input/IIDIPUS_Input_NMAR/ODDobjects/"
filez<-list.files(folder)
filez<-filez[filez!="EQ20170529IDN_136"]

mnlo<-mnla<-10000; mxlo<-mxla<-0
mnlof<-mnlaf<-c(); mxlof<-mxlaf<-c()
for(fff in filez){
  
  ODDy<-readRDS(paste0(folder,fff))
  
  mnlo<-min(mnlo,ODDy@grid@cells.dim[1]);mnla<-min(mnla,ODDy@grid@cells.dim[2])
  mxlo<-max(mxlo,ODDy@grid@cells.dim[1]);mxla<-max(mxla,ODDy@grid@cells.dim[2])
  
  mnlof<-c(mnlof,ODDy@grid@cells.dim[1]);mnlaf<-c(mnlaf,ODDy@grid@cells.dim[2])
  mxlof<-c(mxlof,ODDy@grid@cells.dim[1]);mxlaf<-c(mxlaf,ODDy@grid@cells.dim[2])
  
}

reddie<-2#min(mnlo,mnla)
CNNdim<-3
padding<-CNNdim-1 # this is the size of the filter -1

finDim<-c((mxlo+reddie-mxlo%%reddie)/reddie,
          (mxla+reddie-mxla%%reddie)/reddie)



nevs<-length(filez)

padArray<-function(arraz,padding=2){
  out<-array(NA_real_,c(dim(arraz)[1:2]+padding*2L,dim(arraz)[3]))
  out[(padding+1):(dim(arraz)[1]+padding),(padding+1):(dim(arraz)[2]+padding),]<-arraz
  return(out)
}

resizeArray<-function(arraz,reddie=4,fn=mean){
  
  # indexes in the input array to apply the function over
  indie<-seq.int(1,max(dim(arraz)[1:2]),reddie)
  
  if(length(dim(arraz))==3) {
    
    out<-array(NA_real_,c(dim(arraz)[1:2]/reddie,dim(arraz)[3]))
    
    for(k in 1:dim(arraz)[3]){
      for(i in 1:nrow(out)){
        for(j in 1:ncol(out)){
          # Highlight the square
          out[i,j,k]<-fn(arraz[indie[i]:(indie[i]+reddie-1),indie[j]:(indie[j]+reddie-1),k],na.rm=T)
        }
      }    
    }
    
  } else {
    
    out<-array(NA_real_,c(dim(arraz)[1:2]/reddie))
    
    for(i in 1:nrow(out)){
      for(j in 1:ncol(out)){
        # Highlight the square
        out[i,j]<-fn(arraz[indie[i]:(indie[i]+reddie-1),indie[j]:(indie[j]+reddie-1)],na.rm=T)
      }
    }
  }
  
  return(out)
  
}


# outer<-array(NA,dim = c(nevs,finDim+2L*padding,8))
outer<-c()
impies<-data.frame()
for(i in seq_along(filez)){
  
  fff<-filez[i]
  # Read in the ODD object
  ODDy<-readRDS(paste0(folder,fff))
  # Get rid of unnecessary others
  ODDy@data%<>%dplyr::select_if(!names(.) %in% c("Longitude","Latitude","nBuildings","ISO3C","nBuiltup"))
  # Make sure impact data is on log-scale
  ODDy@impact$observed<-log(ODDy@impact$observed+10)
  # Put it into the impacts file
  impies%<>%rbind(ODDy@impact)
  # Make MaxHaz function over all EQ fore and aftershocks:
  if(length(names(ODDy)[grepl("hazMean",names(ODDy))])==1){
    ODDy@data$hazMax<-ODDy@data$hazMean1
  } else {
    hazard<-rep(NA_real_,length(ODDy@data$hazMean1))
    for (variable in names(ODDy)[grepl("Mean",names(ODDy))]){
      tmp<-ODDy[variable]
      tmp$hazard<-hazard
      hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
    }
    ODDy@data$hazMax<-hazard-ODDy@I0
  }
  # And for the S.D. of the hazard intensity
  if(length(names(ODDy)[grepl("hazSD",names(ODDy))])==1){
    hazard<-ODDy@data$hazSD1
  } else {
    for (variable in names(ODDy)[grepl("hazSD",names(ODDy))]){
      tmp<-ODDy[variable]
      tmp$hazard<-hazard
      hazard<-apply(tmp@data,1,function(x) mean(x,na.rm=T))
    }
  }
  # Get rid of the other hazard info and leave only hazMax
  ODDy@data<-ODDy@data[,!grepl("hazSD",names(ODDy)) & 
                         !grepl("hazMean",names(ODDy))]
  # Now we add the hazSD column back on
  ODDy@data$hazSD<-hazard
  rm(hazard)
  # Reorder things for later
  datar<-ODDy@data%>%dplyr::select(Population,hazMax,everything())
  # Extract impact data per polygon
  for(im in 1:nrow(ODDy@impact)){
    impacts<-ODDy@impact[im,]
    # Remove all elements that lie outside of the area
    datar[!1:nrow(datar)%in%ODDy@polygons[[impacts$polygon]]$indexes,]<-NA
    # Remove all nan values
    removers<-is.na(datar$hazMax)     |
      is.na(datar$Population) |
      datar$hazMax<ODDy@I0    |
      datar$Population<=0
    # Need to differentiate for the array resizing and padding done later
    datar$Population[removers]<-0
    datar[removers,-1]<-NA
    # Convert from ODD object to SPDF to Array
    tmp<-SpatialPixelsDataFrame(ODDy@coords,datar,grid = ODDy@grid,proj4string = ODDy@proj4string)%>%
      convSPDF2Array()
    # What are the missing rows to make it divisible by the reduction factor?
    addie<-reddie-ODDy@grid@cells.dim%%reddie
    
    # Add these values as padding, with zeros
    dimmie<-unname(c(addie+ODDy@grid@cells.dim,ncol(datar)-1))
    # Handle population data first
    pop<-array(0,dim = dimmie[1:2])
    pop[1:nrow(tmp),1:ncol(tmp)]<-tmp[,,1]
    # Resize it using the sum of the population
    pop%<>%resizeArray(reddie,fn=sum)
    # Now handle the remaining columns, by taking the average
    oth<-array(NA_real_,dim = c(dimmie[1:2],dimmie[3]))
    oth[1:nrow(tmp),1:ncol(tmp),]<-tmp[,,-1]; rm(tmp)
    # Resize it using the sum of the population
    oth%<>%resizeArray(reddie,fn=mean)
    # Combine both arrays into one
    out<-array(NA_real_,dim = c(finDim,dim(oth)[3]+1))
    out[1:nrow(pop),1:ncol(pop),1]<-pop; out[1:nrow(pop),1:ncol(pop),-1]<-oth
    # Cleaning!
    rm(oth,pop)
    # Now pad me out!
    out%<>%padArray(padding=padding)
    # # Replace all NAs in hazMax with 0
    # tmp<-out[,,2]; tmp[is.na(tmp)]<-0
    # out[,,2]<-tmp; rm(tmp)
    # Convert all the NAs in the remaining columns with the average values, to not confuse the CNN
    out[is.na(out)]<-0
    # The mother array
    # outer[i,,,]<-out
    # outer%<>%c(list(out=out, impact=impacts))
    outer%<>%c(list(out))
    
  }
  saveRDS(list(outer=outer,impies=impies),paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/Input-CNN-Data_ODD_2.RData"))
  print(paste0("Finished EQ: ",fff))
}

innies<-c()
for(i in 1:nrow(impies)) if(!all(is.na(outer[[i]])) & !is.na(impies$observed[i])) innies%<>%c(T) else innies%<>%c(F)

if(sum(innies)!=length(innies)) stop("Errors in the data wrangling, remove some entries from outer")

inpy<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/Input-CNN-Data_ODD_2.RData"))
outer<-array(NA,dim = c(length(inpy$outer),dim(inpy$outer[[1]])))
for(i in 1:length(inpy$outer)) outer[i,,,]<-inpy$outer[[i]]

saveRDS(list(outer=outer,impies=inpy$impies),paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/Input-CNN-Data_ODD_full.RData"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%% MODEL CNNS! #%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

Hyperparams<-list(cnnfilters=1,
                  kerneldim=c(3,3),
                  poolsize=2,
                  denselayers=3,
                  epocher=100,
                  SCV=10,
                  nvul=1:2)

# Let's do this!
oddCNN<-function(outer,impies,cnnfilters,poolsize,denselayers,droppie=0.2){
  
  Hyperparams$cnnfilters<-cnnfilters
  Hyperparams$poolsize<-poolsize
  Hyperparams$denselayers<-denselayers
  
  performance<-data.frame()
  cnamers<-c("Loss")
  for(j in 1:3){
    indies <- caret::createFolds(1:dim(outer)[1], k = Hyperparams$SCV, list = T, returnTrain = FALSE)
    #@@@@@@@@@@@@@@@@@@@@@ STRATIFIED CROSS-VALIDATION @@@@@@@@@@@@@@@@@@@@@#
    for(cv in 1:Hyperparams$SCV){
      # Split it up!
      xtrain <- outer[!1:dim(outer)[1]%in%indies[[cv]],,,]
      xtest  <- outer[indies[[cv]],,,]
      ytrain <- impies$observed[!1:dim(outer)[1]%in%indies[[cv]]]
      ytest  <- impies$observed[indies[[cv]]]
      ################# CNN SECTION #################
      cnn_model <- keras_model_sequential() %>%
        # Data augmentation
        layer_random_flip() %>%
        layer_random_rotation(0.2)%>%
        layer_conv_2d(filters = Hyperparams$cnnfilters, 
                      kernel_size = Hyperparams$kerneldim,
                      # activation = actie, 
                      input_shape = c(dim(outer)[2:3],length(Hyperparams$nvul))) %>%
        layer_max_pooling_2d(pool_size = c(Hyperparams$poolsize, Hyperparams$poolsize)) %>%
        layer_flatten() %>%
        layer_dropout(droppie)%>%
        # layer_dense(units = Hyperparams$denselayers) %>%
        layer_dense(units = 1)
      # summary(cnn_model)
      # Compile it
      cnn_model %>% compile(
        loss = 'mean_absolute_error', # for some reason this was preferred over binary
        optimizer = optimizer_adam() # could use optimizer_adadelta() or optimizer_adam() or optimizer_sgd()
      )
      # Fit the model!
      cnn_history <- cnn_model %>% fit(
        xtrain, ytrain,
        # Batch size taken from DeFine
        batch_size = length(ytrain),
        # epochs taken from DeFine
        epochs = Hyperparams$epocher,
        validation_split = 0.0,
        verbose=0,
        callbacks = list(callback_early_stopping(monitor = "loss", patience = 20, restore_best_weights = TRUE))
      )
      
      tmp<-cnn_model%>%evaluate(xtest,ytest)
      
      # Bind to the output file
      performance%<>%rbind(cbind(as.data.frame(as.list(t(tmp)),col.names=cnamers),
                                 data.frame(filters=Hyperparams$cnnfilters,
                                            denselayers=Hyperparams$denselayers,
                                            poolsize=Hyperparams$poolsize,
                                            CVfold=cv,j=j)))
      # }
    }
  }
  # colnames(performance)<-cnamers
  out<-data.frame(avLoss=mean(performance$Loss),
                  sdLoss=sd(performance$Loss),
                  filters=Hyperparams$cnnfilters,
                  denselayers=Hyperparams$denselayers,
                  poolsize=Hyperparams$poolsize,
                  dropout=droppie)
  
}

inpy<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/Input-CNN-Data_ODD_full.RData"))

outer<-inpy$outer[inpy$impies$impact=="mortality",,,Hyperparams$nvul]
impies<-inpy$impies[inpy$impies$impact=="mortality",]

maxies<-apply(outer,4,max)
for(i in seq_along(maxies)) outer[,,,i]<-outer[,,,i]/maxies[i]

# impies$observed<-(impies$observed-min(impies$observed))/(max(impies$observed)-min(impies$observed))
impies$observed<-impies$observed-min(impies$observed)

performance<-data.frame()

# for(ac in c("relu","sigmoid")){
  for(fff in c(1)){
    for(ps in c(4,6,8)){
      for(dl in c(1)){  
        for(dp in c(0.2,0.4,0.6)){  
          # out<-tryCatch(oddCNN(outer,fff,ps,dl,dp),error=function(e) NA)
          out<-oddCNN(outer,impies,fff,ps,dl,dp)
          # if(any(is.na(out))) next
          performance%<>%rbind(out)
          saveRDS(performance,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/CNN_performance_3.RData")
        }
      }
    }
  }
# }

morties<-data.frame()
# Final formulations
for(dp in c(0.2,0.4,0.6)){  
  out<-oddCNN(outer,impies,1,5,1,0.5)
  morties%<>%rbind(out)
}

outer<-inpy$outer[inpy$impies$impact=="displacement",,,Hyperparams$nvul]
impies<-inpy$impies[inpy$impies$impact=="displacement",]
maxies<-apply(outer,4,max)
for(i in seq_along(maxies)) outer[,,,i]<-outer[,,,i]/maxies[i]
impies$observed<-impies$observed-min(impies$observed)

dispies<-data.frame()
# Final formulations
for(ps in c(4,6,8)){
  for(dp in c(0.2,0.4,0.6)){  
    out<-oddCNN(outer,impies,1,ps,1,dp)
    dispies%<>%rbind(out)
  }
}

outer<-inpy$outer[inpy$impies$impact=="buildDam",,,Hyperparams$nvul]
impies<-inpy$impies[inpy$impies$impact=="buildDam",]
maxies<-apply(outer,4,max)
for(i in seq_along(maxies)) outer[,,,i]<-outer[,,,i]/maxies[i]
impies$observed<-impies$observed-min(impies$observed)

bdamies<-data.frame()
# Final formulations
for(ps in c(4,6,8)){
  for(dp in c(0.2,0.4,0.6)){  
    out<-oddCNN(outer,impies,1,ps,1,dp)
    bdamies%<>%rbind(out)
  }
}

outer<-inpy$outer[inpy$impies$impact=="buildDest",,,Hyperparams$nvul]
impies<-inpy$impies[inpy$impies$impact=="buildDest",]
maxies<-apply(outer,4,max)
for(i in seq_along(maxies)) outer[,,,i]<-outer[,,,i]/maxies[i]
impies$observed<-impies$observed-min(impies$observed)

bdesties<-data.frame()
# Final formulations
for(ps in c(4,6,8)){
  for(dp in c(0.2,0.4,0.6)){  
    out<-oddCNN(outer,impies,1,ps,1,dp)
    bdesties%<>%rbind(out)
  }
}

morties$impact<-"mortality"
dispies$impact<-"displacement"
bdamies$impact<-"buildDam"
bdesties$impact<-"buildDest"

outout<-rbind(morties,dispies,bdamies,bdesties)

write_csv(outout,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/CNN_allimpacts.csv")

outout<-read_csv("./IIDIPUS_Results/SpatialPolygons_ML-GLM/CNN_allimpacts.csv")

allUV<-fuller%>%arrange(RelativeAbsDiff)%>%group_by(impact,algo)%>%slice(1)%>%
  dplyr::select(impact,algo,RelativeAbsDiff,RelativeAbsDiffSD)%>%filter(algo!="lm")

colnames(allUV)
outout$algo<-"CNN"
outout%<>%arrange(avLoss)%>%group_by(impact)%>%slice(1)%>%
  dplyr::select(impact,algo,avLoss,sdLoss)

colnames(outout)[3:4]<-c("RelativeAbsDiff","RelativeAbsDiffSD")

allUV%<>%rbind(outout)

pal <- c(
  "buildDam" = scales::hue_pal()(4)[3], 
  "buildDest" = scales::hue_pal()(4)[4],
  "displacement" = scales::hue_pal()(4)[2], 
  "mortality" = scales::hue_pal()(4)[1]
)

p<-allUV%>%ggplot(aes(algo,RelativeAbsDiff,group=impact))+
  # geom_point(aes(colour=impact,shape=impact),size=4) +
  geom_pointrange(aes(ymin=RelativeAbsDiff-RelativeAbsDiffSD,
                      ymax=RelativeAbsDiff+RelativeAbsDiffSD,colour=impact,
                      shape=impact),
                  size=1)+theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_y_log10()+scale_shape_manual(values=15:18,breaks=allimps)+
  xlab("Model") + ylab("Mean Absolute Deviation of Logs (MADL)")+
  scale_colour_manual(values = pal,limits = names(pal))+
  labs(colour="Impact Type",shape="Impact Type");p
ggsave("UV_MADL_CNN.eps",p,path="./Plots/IIDIPUS_Results/",width=10,height=4.,device = grDevices::cairo_ps)  


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%% MODEL TURKEY! #%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

folder<-"./IIDIPUS_Input/IIDIPUS_Input_All_2023May19/ODDobjects/"
# Read in the Turkey earthquake
ODDy<-readRDS(paste0(folder,"EQ20230206TUR_169"))
# Convert it into the necessary form
outT<-convODD(data.frame(),ODDy,"EQ20230206TUR_169")%>%filter(outT%in%c("LBN","TUR","SYR"))
# Read in the rest of the events
out<-readRDS("./IIDIPUS_Results/SpatialPolygons_ML-GLM/InputData_ODD.RData")
# Add the time component before binding
outT$time<-as.numeric(outT$date-min(out$date))
# Bind them
out%<>%rbind(outT)
# Scale and adjust
outred<-scaleIMPs(out)

tursyr<-outred%>%filter(Event=="EQ20230206TUR_169")
outred%<>%filter(Event!="EQ20230206TUR_169")

ncores<-60
cl <- makePSOCKcluster(ncores)  # Create computing clusters
registerDoParallel(cl)
getDoParWorkers()
# Run it!
ODD_RF<-do.call(rbind,lapply(allimps, function(impact) {
    out<-parallelML("rf",impact,modRet = T)
    return(data.frame(pred=(exp(unname(predict(out,tursyr)))-10),
                      true=tursyr[[impact]],
                      impact=impact,
                      ISO3C=outT$iso3,
                      model="rf"))
  }))

ODD_Rad<-do.call(rbind,lapply(allimps, function(impact) {
  out<-parallelML("svmRadial",impact,modRet = T)
  return(data.frame(pred=(exp(unname(predict(out,tursyr)))-10),
                    true=tursyr[[impact]],
                    impact=impact,
                    ISO3C=outT$iso3,
                    model="svmRadial"))
}))
# Remember to close the computing cluster
stopCluster(cl)
registerDoSEQ()

ODD_ML<-rbind(ODD_RF,ODD_Rad)
ODD_ML%<>%mutate(MADL=abs(log(pred+10)-log(true+10)))

View(ODD_ML)

saveRDS(ODD_ML,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/TUR-SYR_preds.RData")

TURpreds<-do.call(rbind,lapply(1:nrow(tursyr),function(i){
  out<-GLMvarImp(GLMmods,tursyr[i,])
  return(cbind(out,data.frame(ISO3C=outT$iso3[i])))
}))

TURpreds%<>%group_by(impact,ISO3C)%>%
  summarise(pred=exp(mean(log(pred+10))-10),
            true=mean(true),
            MADL=abs(log(pred+10)-log(true+10)),
            model="GLM")

finfin<-rbind(TURpreds,ODD_ML)

finfin%>%group_by(model)%>%
  summarise(avMADL=mean(MADL,na.rm=T))

finfin%>%group_by(model)%>%filter(ISO3C=="TUR")%>%View()

ExpHaz<-rbind(out,outT)%>%dplyr::select(grep("Exp",colnames(out),value = T),Event,iso3)%>%
  dplyr::select(-LifeExp,-ExpSchYrs)%>%reshape2::melt()
ExpHaz$variable%<>%extractnumbers()%>%as.character()

p<-ExpHaz%>%ggplot()+geom_boxplot(aes(variable,value,fill=variable))+scale_y_log10()+
  geom_point(data = filter(ExpHaz,iso3=="TUR" & Event=="EQ20230206TUR_169"),aes(variable,value),colour="red")+
  scale_fill_ordinal()+labs(fill="Intensity")+
  xlab("Min. Exposed Earthquake Intensity [MMI]")+ylab("Exposed Population");p
ggsave("TUR_HazExp.eps",p,path="./Plots/IIDIPUS_Results/",width=10,height=5,device = grDevices::cairo_ps)  

outred%>%ggplot()+geom_histogram(aes(ExpDim2))+scale_x_log10()
# 

ncores<-60
cl <- makePSOCKcluster(ncores)  # Create computing clusters
registerDoParallel(cl)
getDoParWorkers()
# Run it!
ODD_all<-do.call(rbind,lapply(allimps, function(impact) {
  out<-parallelML("rf",impact,modRet = T)
  return(data.frame(pred=(exp(unname(predict(out,outred)))-10),
                    true=outred[[impact]],
                    impact=impact,
                    Event=outred$Event,
                    ISO3C=outT$iso3,
                    model="rf"))
}))
# Remember to close the computing cluster
stopCluster(cl)
registerDoSEQ()

ODD_all%<>%mutate(MADL=abs(log(pred+10)-log(true+10)))

p<-ODD_all%>%ggplot()+geom_boxplot(aes(impact,MADL,fill=impact))+
  geom_point(data = filter(ODD_ML,model=="rf" & ISO3C=="TUR"),aes(impact,MADL),colour="red")+
  xlab("Impact Type")+ylab("MADL Prediction Error")+labs(fill="Impact Type")+
  scale_fill_manual(values = pal,limits = names(pal))+
  ggtitle("Random Forest (Top Model) Predictions")+
  theme(plot.title = element_text(hjust = 0.5));p
ggsave("agg_RF_Predictions.eps",p,path="./Plots/IIDIPUS_Results/",width=8,height=5,device = grDevices::cairo_ps)  


MMIlevels<-5000:9500/1000

# Proof that Turkiye EQ event is an outlier
expPop_MMI<-function(ODDy,Event){
  # Make MaxHaz function over all EQ fore and aftershocks:
  if(length(names(ODDy)[grepl("hazMean",names(ODDy))])==1){
    ODDy@data$hazMax<-ODDy@data$hazMean1
  } else {
    hazard<-rep(NA_real_,length(ODDy@data$hazMean1))
    for (variable in names(ODDy)[grepl("Mean",names(ODDy))]){
      tmp<-ODDy[variable]
      tmp$hazard<-hazard
      hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
    }
    ODDy@data$hazMax<-hazard
  }
  
  # Calculate the population exposed per MMI level
  return(do.call(rbind,parallel::mclapply(MMIlevels,function(mmi){
    return(data.frame(Event=Event, MMI=mmi,
                      PopExp=sum(ODDy@data$Population[ODDy@data$hazMax>=mmi],
                                 na.rm = T)))
  },mc.cores = ncores)))

}

# Run it across all events
expop<-data.frame()
for(fff in filez){
  # Exceptions... sigh
  if(fff=="EQ20170529IDN_136") next
  # Read in the hazard & impact object
  ODDy<-readRDS(paste0(folder,fff))
  # Convert it into the necessary form
  expop%<>%rbind(expPop_MMI(ODDy,fff))
  
  print(paste0("Finished EQ: ",fff))
}
ODDy<-readRDS(paste0(folder,"EQ20230206TUR_169"))
expop%<>%rbind(expPop_MMI(ODDy,"EQ20230206TUR_169"))

p <- expop %>% ggplot() + geom_line(aes(MMI, PopExp, colour=Event),alpha=0.5,size=0.5) + 
  theme(legend.position = "none") + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits=c(1e3,1e8)) + annotation_logticks(sides="l") +
  xlab("USGS Shakemap Intensity [MMI]") + ylab("Exposed Population [Cumulative]") +
  geom_point(data = expop%>%filter(Event=="EQ20230206TUR_169"),mapping = aes(MMI,PopExp), colour="black", size=2)

ggsave("TUR_outlier_ExpPop-MMI.eps",p,path="./Plots/IIDIPUS_Results/",width=8,height=5,device = grDevices::cairo_ps)  


costie<-function(data, lev = NULL, model = NULL) {
  out<-mean(abs(data$obs-data$pred)*data$weights/sum(data$weights),na.rm = T)
  names(out)<-"RelativeAbs"
  return(out)
}

train_control <- caret::trainControl(method="repeatedcv", number=10, repeats=3,
                                     search = "random",
                                     summaryFunction=costie)
# MADL example with random forest algorithm, by plotting the observed vs predicted impacts
rfres<-do.call(rbind,lapply(allimps,function(impact){
  
  outFrame<-dplyr::select(outred,-allimps[impact!=allimps])
  names(outFrame)[1]<-"Y"
  outFrame$Y<-log(outFrame$Y+10)
  outFrame%<>%na.omit()
  # Make weights from the different events to make sure that no single event dominates the model parameterisation
  weights<-outFrame%>%group_by(Event)%>%summarise(www=1/length(time))%>%merge(outFrame)%>%pull(www)
  # Remove the variable Event after weighting is calculated
  outFrame%<>%dplyr::select(-Event)
  # Train the model!
  modeler<-caret::train(Y~., data = outFrame, method = "rf", metric="RelativeAbs",
                        tuneLength = 12, trControl = train_control,linout = TRUE,
                        weights = weights, preProcess = c("center","scale"))
  # Get the predictions
  ybar=unname(predict(modeler,outFrame))
  
  data.frame(y=exp(outFrame$Y)-10, 
             ybar=exp(ybar)-10,
             MADL=abs(outFrame$Y-ybar),
             impact=impact)
}))

# Due to messing around with log(x+10) we need to round the observed values
rfres$y%<>%round()

# Function to print out the numbers how we want them
scinote <- function(x) {
  # Get the exponent of the number
  exponent <- floor(log10(abs(x)))
  # Round the number to 2 significant figures
  rounded_value <- round(x / 10^exponent, 2)
  # Format the number in scientific notation with 1 significant figures
  formatted_value <- sprintf("%.1f", rounded_value)
  # Construct the string in scientific notation
  scientific_notation <- paste0(formatted_value, "e", exponent)
  
  return(scientific_notation)
}

# Create the text to go onto the plot, showing the adj-R^2 and the L1 + L2 norms
texty<-do.call(rbind,lapply(allimps,function(imp){
  # Filter the specificimpact type
  tmp<-rfres%>%filter(impact==imp)
  # Extract the metrics of interest
  adjR2<-summary(lm(y ~ ybar , data=tmp))$adj.r.squared
  L1<-mean(abs(tmp$y-tmp$ybar))
  L2<-mean((tmp$y-tmp$ybar)^2)
  # Find a nice way to calculate the y-value of the text location
  yv<-quantile(tmp$ybar,probs=0.99); yv<-c(yv,yv*(0.7-0.05*log10(yv[1])))
  y_diff <- (log10(yv[2]) - log10(yv[1]))  # Difference between the first two entries in log space
  yv <- c(yv, 10^(log10(yv[2]) + y_diff)) 
  
  data.frame(impact=imp,
    text=c(
      paste0("adj-R^2 = ",signif(adjR2,2)),
      paste0("L1 / N = ",scinote(L1)),
      paste0("L2 / N = ",scinote(L2))
    ),
    variable=c("adj-R^2","L1","L2"),
    x_pos=1, y_pos=yv)
}))

rfres%<>%mutate(impact=factor(impact,levels(),ordered=T))
# Plot it!
p<-rfres%>%ggplot()+geom_point(aes(y,ybar,size=MADL,colour=impact))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits=c(1,NA))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits=c(1,NA))+
  geom_abline(slope = 1,intercept = 0,colour="black")+
  scale_size_continuous(breaks = c(0.05,0.1,0.5,1,5)) +
  annotation_logticks() +
  labs(colour="Impact")+xlab("Observed Impact")+ylab("Predicted Impact")+
  ggtitle("Top Model: Random Forest")+theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~impact,nrow=2,scales = "free")+
  geom_text(data = texty,
            aes(label = text, x = 1, y = y_pos),
            hjust = 0, vjust = 0, size = 3);p

ggsave("RF_Perf_wMADL.eps",p,path="./Plots/IIDIPUS_Results/",width=10,height=8,device = grDevices::cairo_ps)  

# Melt the dataframe to plot
# resultsUV%<>%reshape2::melt(id.vars=1:4)
# # Make a weighting variable for the alpha plotting
# resultsUV$normError<-1/resultsUV$StandErr
# resultsUV%<>%group_by(impact)%>%mutate(normError=normError/max(normError))%>%ungroup()
# # Plot
# resultsUV%>%ggplot(aes(as.factor(variable),value,group=algo))+geom_point(aes(colour=algo,alpha=normError))+
#   facet_wrap(.~impact,nrow = 4)

#  HYPOTHESIS TESTING ON BEST PERFORMING MODEL
# resultsUV%<>%group_by(impact)%>%mutate(nn=sum(!is.na(outred[,unique(impact)])))
# resultsUV%<>%group_by(impact)%>%mutate(withBest=dnorm(RelativeAbs-min(RelativeAbs),0,
#    sqrt(RelativeAbsSD^2+RelativeAbsSD[which.min(RelativeAbs)]^2)/sqrt(nn))<0.05)%>%
#   dplyr::select(-nn)
# 
# resultsUV%<>%group_by(impact)%>%mutate(withBest=dnorm(RelativeAbs-min(RelativeAbs),0,
#                                                       sqrt(RelativeAbsSD^2+RelativeAbsSD[which.min(RelativeAbs)]^2)/sqrt(nn))<0.05)%>%
#   dplyr::select(-nn)

# tabUV<-lapply(unique(resultsUV$impact), function(imp){
#   xtable::xtable(filter(resultsUV,impact==imp & withBest),"Feature importance ranking measure (FIRM), scaled as a percentage, for the \'best-performing\' models. The definition of \'best performing\' here is defined as any model that had a cost that was not statistically significantly different from the overall best performing model for each specific impact. When applying the FIRM, Individual Conditional Expectation (ICE) curves are used.")})
# 
# for(i in 1:length(tabUV)) print(tabUV[[i]],include.rownames=F)
# 
# # The models that made it as one of the best performing, per impact:
# resultsUV[resultsUV$withBest,]%>%group_by(impact)%>%reframe(models=unique(model))%>%View()
# resultsUV[resultsUV$withBest,]%>%group_by(impact)%>%reframe(models=length(model))
# nrow(resultsUV)
# 
# # Order the variables by their importance, weighted by each models error value
# varimport<-do.call(rbind,lapply(allimps,function(imps) {tmp<-filter(resultsUV,impact==imps);apply(tmp[tmp$withBest,-c(1:4,ncol(tmp),(ncol(tmp)-1))],2,function(x) weighted.mean(x,1/tmp$RelativeAbs[tmp$withBest]))}))
# varimport<-100*varimport/rowSums(varimport)
# varimport%<>%rbind(as.numeric(colMeans(varimport)))
# rownames(varimport)<-c(allimps,"average")
# varimport<-100*varimport/rowSums(varimport)
# varimport%<>%t()
# varimport%<>%as.data.frame(row.names = rownames(varimport))%>%mutate(Covariate=rownames(varimport))
# 
# varimport%<>%reshape2::melt(id.vars=6)
# colnames(varimport)[2]<-"impact"
# 
# pal <- c(
#   "mortality" = "red",
#   "displacement" = "blue", 
#   "buildDam" = "forestgreen", 
#   "buildDest" = "purple",
#   "average" = "black"
# )
# varimport%>%ggplot(aes(Covariate,value,group=impact))+geom_point(aes(colour=impact,shape=impact),size=2)+
#   geom_line(aes(colour=impact),alpha=0.25)+scale_colour_manual(values = pal,limits = names(pal))
# 
# 
# 
# 
# 
# 
# 
# colnames(resultsUV)[2:4]<-c("RelativeAbs","RelativeAbsSD","maxHaz")
# 
# predictionsML%<>%rbind(resultsUV)
# 
# varimp<-predictionsML%>%arrange(RelativeAbs)%>%group_by(model)%>%
#   summarise(Cost=prod(RelativeAbs),
#             maxHaz=mean(maxHaz),
#             ExpSchYrs=mean(ExpSchYrs),
#             LifeExp=mean(LifeExp),
#             GNIc=mean(GNIc),
#             Vs30=mean(Vs30),
#             EQFreq=mean(EQFreq),
#             hazSD=mean(hazSD),
#             time=mean(time),
#             ExpDim1=mean(ExpDim1),
#             ExpDim2=mean(ExpDim2),
#             WIDDim1=mean(WIDDim1),
#             WIDDim2=mean(WIDDim2))
# 
# varimp[3:ncol(varimp)]<-100*varimp[3:ncol(varimp)]/
#   rowSums(varimp[3:ncol(varimp)])
# # 
# varimp%>%filter(model%in%c("rf","svmRadial","svmPoly","lognorm","nnet"))%>%xtable::xtable()%>%print(include.rownames=F)

