# Extract Environment Variables
# dir<-directory<-"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/";setwd(directory); packred<-T
# Download and install the necessary packages:
# source('RCode/GetODDPackages.R')
source('RCode/ODDobj.R')
library(dplyr)
library(magrittr)
library(tidyverse)
library(boot)
library(MASS)
library(pscl)
library(FactoMineR)
library(factoextra)
library(parallel)
library(doParallel)

# Used to correlate the vulnerability variables with the modifier term
LMFeatureSelection<-function(output,Nb=30,GLMer="LM",mvm=NULL,intercept=F,fn="+",nlim=3,weights=NULL,ncores=4){
  # Nb - Number of times LM is run with different samples of training vs test data
  # intercept - do we include an intercept in the equation?
  
  # Setup the GLM here:
  if(!is.null(mvm) & GLMer%in%c("LM","lognorm")){
    # Change between the normal (linear) model and the log-normal
    if(GLMer=="LM") {
      lognorm=F 
      fncy<-function(x) x
    } else {
      lognorm=T
      fncy<-exp
    }
  
    modely<-function(equation,datar,responsers){
      # If we are using a lognormal model
      if(lognorm) datar[,responsers]<-log(datar[,responsers]+10)
      # Remove all the NAs
      datar%<>%na.omit()
      # Set up the cross validation folds
      indies <- createFolds(1:nrow(datar), k = 10, list = T, returnTrain = FALSE)
      # Train a model for each fold
      as.data.frame(t(colMeans(do.call(rbind, lapply(1:10,function(i){
        test<-datar[indies[[i]],]
        train<-datar[!1:nrow(datar)%in%indies[[i]],]
        modeler<-lm(formula = as.formula(equation),data = train)
        as.data.frame(t(colMeans(abs(fncy(predict(modeler,test))-fncy(test[,responsers]))/(fncy(test[,responsers])+1),na.rm = T)))
      })))))
    }
    
  } else if(GLMer=="LM"){
    modely<-function(equation,datar,responsers="Y") {
      datar%<>%na.omit()
      modeler<-glm(formula = as.formula(equation),data = datar,family = gaussian(link = "identity"))
      StandErr<-cv.glm(data = datar, glmfit = modeler, K = 10, cost=function(y,yhat) mean(abs(y-yhat)/(y+1),na.rm = T))$delta[1]
      
      data.frame(StandErr=StandErr,
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="pois"){
    modely<-function(equation,datar,responsers="Y") {
      datar%<>%na.omit()
      modeler<-glm(formula = as.formula(equation),data = datar,family = poisson())
      StandErr<-cv.glm(data = datar, glmfit = modeler, K = 10, cost=function(y,yhat) mean(abs(y-yhat)/(y+1),na.rm = T))$delta[1]
      
      data.frame(StandErr=StandErr,
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="lognorm"){
    modely<-function(equation,datar,responsers="Y") {
      datar%<>%na.omit()
      datar$Y<-log(datar$Y+10)
      modeler<-glm(formula = as.formula(equation),data = datar,family = gaussian(link = "identity"))
      StandErr<-cv.glm(data = datar, glmfit = modeler, K = 10, cost=function(y,yhat) mean(abs(exp(y)-exp(yhat))/(y+1),na.rm = T))$delta[1]
      
      data.frame(StandErr=StandErr,
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="HurdlePois"){
    modely<-function(equation,datar,responsers="Y") {
      datar%<>%na.omit()
      modeler<-pscl::hurdle(formula = as.formula(equation),data = datar,dist="poisson")
      StandErr<-cv.glm(data = datar, glmfit = modeler, K = 10, cost=function(y,yhat) mean(abs(y-yhat)/(y+1),na.rm = T))$delta[1]
      
      data.frame(StandErr=StandErr,
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="HurdleNegBin"){
    modely<-function(equation,datar,responsers="Y") {
      datar%<>%na.omit()
      modeler<-pscl::hurdle(formula = as.formula(equation),data = datar,dist="negbin")
      StandErr<-cv.glm(data = datar, glmfit = modeler, K = 10, cost=function(y,yhat) mean(abs(y-yhat)/(y+1),na.rm = T))$delta[1]
      
      data.frame(StandErr=StandErr,
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="ZInegbin"){
    modely<-function(equation,datar,responsers="Y") {
      datar%<>%na.omit()
      modeler<-pscl::zeroinfl(formula = as.formula(equation),data = datar,dist="negbin")
      StandErr<-cv.glm(data = datar, glmfit = modeler, K = 10, cost=function(y,yhat) mean(abs(y-yhat)/(y+1),na.rm = T))$delta[1]
      
      data.frame(StandErr=StandErr,
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else if(GLMer=="ZIpois"){
    modely<-function(equation,datar,responsers="Y") {
      datar%<>%na.omit()
      modeler<-pscl::zeroinfl(formula = as.formula(equation),data = datar,dist="poisson")
      StandErr<-cv.glm(data = datar, glmfit = modeler, K = 10, cost=function(y,yhat) mean(abs(y-yhat)/(y+1),na.rm = T))$delta[1]
      
      data.frame(StandErr=StandErr,
                 AIC=AIC(modeler),
                 BIC=BIC(modeler))
    }
    
  } else stop("Regression model not recognised")
  
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
      ressies<-tryCatch(modely(eee,output,responsers),error = function(e) NA)
      if(all(is.na(ressies))) next
      
      predictors%<>%rbind(ressies)
    }
    
    return(c(t(colMeans(predictors))))
    
  })
  # Fudging between univariate and multivariate response models
  if(is.null(mvm)) {
    prediction<-as.data.frame(t(matrix(unlist(prediction),nrow = 3)))
    colnames(prediction)<-c("StandErr","AIC","BIC")
  } else {
    prediction<-as.data.frame(t(matrix(unlist(prediction),nrow = length(responsers))))
    colnames(prediction)<-responsers
  }
  
  return(cbind(data.frame(eqn=eqeqeq),prediction))
  
}

# folder<-"./IIDIPUS_Input/IIDIPUS_Input_NMAR/ODDobjects/"
# filez<-list.files(folder)
# 
# out<-data.frame()
# for(fff in filez){
# 
#   if(fff=="EQ20170529IDN_136") next
# 
#   ODDy<-readRDS(paste0(folder,fff))
# 
#   # Make MaxHaz function over all EQ fore and aftershocks:
#   if(length(names(ODDy)[grepl("hazMean",names(ODDy))])==1){
#     ODDy@data$hazMax<-ODDy@data$hazMean1
#   } else {
#     hazard<-rep(NA_real_,length(ODDy@data$hazMean1))
#     for (variable in names(ODDy)[grepl("Mean",names(ODDy))]){
#       tmp<-ODDy[variable]
#       tmp$hazard<-hazard
#       hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
#     }
#     ODDy@data$hazMax<-hazard
#   }
#   # And for the S.D. of the hazard intensity
#   if(length(names(ODDy)[grepl("hazSD",names(ODDy))])==1){
#     ODDy@data$hazSD<-ODDy@data$hazSD1
#   } else {
#     for (variable in names(ODDy)[grepl("hazSD",names(ODDy))]){
#       tmp<-ODDy[variable]
#       tmp$hazard<-hazard
#       hazard<-apply(tmp@data,1,function(x) mean(x,na.rm=T))
#     }
#     ODDy@data$hazSD<-hazard
#     rm(hazard)
#   }
# 
#   for (iso in unique(ODDy@impact$iso3)){
# 
#     if(!iso%in%ODDy@cIndies$iso3) next
# 
#     inISO<-!is.na(ODDy@data$ISO3C) & ODDy@data$ISO3C==iso
# 
#     tmp<-ODDy@cIndies%>%filter(iso3==iso)%>%dplyr::select(c(1,2))
#     WID<-tmp$value%>%t()%>%as.data.frame(); colnames(WID)<-tmp$percentile; rm(tmp)
# 
#     for (ppp in unique(ODDy@impact$polygon[ODDy@impact$iso3==iso])){
# 
#       # extract which grid points lie within the given polygon
#       inP<-rep(F,nrow(ODDy@data))
#       inP[ODDy@polygons[[ppp]]$indexes] <-T
#       indies<-inISO & inP
# 
#       # Extract impact data per polygon
#       impacts<-ODDy@impact%>%filter(polygon==ppp)%>%summarise(mortality=ifelse(length(observed[impact=="mortality"]==0),sum(observed[impact=="mortality"]),NA),
#                                                      displacement=ifelse(length(observed[impact=="displacement"]==0),sum(observed[impact=="displacement"]),NA),
#                                                      buildDam=ifelse(length(observed[impact=="buildDam"]==0),sum(observed[impact=="buildDam"]),NA),
#                                                      buildDest=ifelse(length(observed[impact=="buildDest"]==0),sum(observed[impact=="buildDest"]),NA))
# 
#       # Weighted (by population) vulnerability values
#       vulny<-ODDy@data[inP & ODDy@data$hazMax>5,]%>%dplyr::select(c(ExpSchYrs,LifeExp,GNIc,Vs30,EQFreq,hazSD,Population))%>%
#         summarise(ExpSchYrs=weighted.mean(ExpSchYrs,Population),
#                   LifeExp=weighted.mean(LifeExp,Population),
#                   GNIc=weighted.mean(GNIc,Population),
#                   Vs30=weighted.mean(Vs30,Population),
#                   EQFreq=weighted.mean(EQFreq,Population),
#                   hazSD=weighted.mean(hazSD,Population))
# 
#       out%<>%rbind(
#         cbind(impacts,data.frame(
#           iso3=iso,
#           date=ODDy@hazdates[1],
#           Exp5=sum(ODDy@data$Population[ODDy@data$hazMax>=5 & indies],na.rm = T),
#           Exp5.5=sum(ODDy@data$Population[ODDy@data$hazMax>=5.5 & indies],na.rm = T),
#           Exp6=sum(ODDy@data$Population[ODDy@data$hazMax>=6 & indies],na.rm = T),
#           Exp6.5=sum(ODDy@data$Population[ODDy@data$hazMax>=6.5 & indies],na.rm = T),
#           Exp7=sum(ODDy@data$Population[ODDy@data$hazMax>=7 & indies],na.rm = T),
#           Exp7.5=sum(ODDy@data$Population[ODDy@data$hazMax>=7.5 & indies],na.rm = T),
#           Exp8=sum(ODDy@data$Population[ODDy@data$hazMax>=8 & indies],na.rm = T),
#           Exp8.5=sum(ODDy@data$Population[ODDy@data$hazMax>=8.5 & indies],na.rm = T),
#           Exp9=sum(ODDy@data$Population[ODDy@data$hazMax>=9 & indies],na.rm = T),
#           maxHaz=max(ODDy@data$hazMax,na.rm = T)
#         ),WID,vulny)
#       )
#     }
#   }
#   print(paste0("Finished EQ: ",fff))
# }
# 
# out$time<-as.numeric(out$date-min(out$date))
# 
# saveRDS(out,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/InputDataGLM.RData")

out<-readRDS("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/InputDataGLM.RData")

# Per impact, do the analysis (for all impacts, even without building data)
outred<-dplyr::select(out,-c("iso3","date")); rm(out)
outred[,-(1:4)]<-scale(outred[,-(1:4)])

# This reveals 95% of variance is contained within 4 dimensions of the PCA
ExpPCA<-outred%>%dplyr::select(colnames(outred)[str_starts(colnames(outred),"Exp") & colnames(outred)!="ExpSchYrs"])%>%PCA(ncp = 5,graph = F)
iExp<-sum(ExpPCA$eig[,3]<95); iExp<-2
tmpExp<-as.data.frame(ExpPCA$ind$coord[,1:iExp]); colnames(tmpExp)<-paste0("ExpDim",1:iExp)
outred%<>%cbind(tmpExp)%>%dplyr::select(-colnames(outred)[str_starts(colnames(outred),"Exp") & colnames(outred)!="ExpSchYrs"])
rm(ExpPCA,iExp,tmpExp)
# This reveals 95% of variance is contained within 2 dimensions of the PCA
WIDPCA<-outred%>%dplyr::select(colnames(outred)[str_starts(colnames(outred),"p")])%>%PCA(ncp = 5,graph = F)
iWID<-sum(WIDPCA$eig[,3]<95)
tmpWID<-as.data.frame(WIDPCA$ind$coord[,1:iWID]); colnames(tmpWID)<-paste0("WIDDim",1:iWID)
outred%<>%cbind(tmpWID)%>%dplyr::select(-colnames(outred)[str_starts(colnames(outred),"p")])
rm(WIDPCA,iWID,tmpWID)

# Remove NAs from the object
outred<-outred[apply(outred[,-(1:4)],1,function(x) !any(is.na(x))),]
# Check the Variance Inflation Factor
print(paste0("Checking the VIF of the remaining variables: percentage that didn't pass = ",sum(multiColl::VIF(outred[,-(1:4)])>6)/ncol(outred[,-(1:4)])))
colnames(outred)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MORTALITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
outFrame<-dplyr::select(outred,-c("displacement","buildDam","buildDest"))
names(outFrame)[1]<-"Y"

# Function from the file CorrelateModifier.R:
fuller<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "LM", ncores=40),
                      error=function(e) NA)
fuller$model<-"LM"

saveRDS(fuller,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Mortality_LM.RData")

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "pois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Mortality_pois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="pois")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "lognorm", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Mortality_lognorm.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="lognorm")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "HurdlePois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Mortality_HurdlePois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdlePois")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "HurdleNegBin", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Mortality_HurdleNegBin.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdleNegBin")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "ZInegbin", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Mortality_ZInegbin.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZInegbin")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "ZIpois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Mortality_ZIpois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZIpois")))

saveRDS(fuller,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Mortality.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISPLACEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
outFrame<-dplyr::select(outred,-c("mortality","buildDam","buildDest"))
names(outFrame)[1]<-"Y"

# Function from the file CorrelateModifier.R:
fuller<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                    GLMer = "LM", ncores=40),
                 error=function(e) NA)
fuller$model<-"LM"
saveRDS(fuller,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Displacement_LM.RData")

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "pois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Displacement_pois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="pois")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "lognorm", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Displacement_lognorm.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="lognorm")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "HurdlePois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Displacement_HurdlePois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdlePois")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "HurdleNegBin", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Displacement_HurdleNegBin.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdleNegBin")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "ZInegbin", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Displacement_ZInegbin.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZInegbin")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "ZIpois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Displacement_ZIpois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZIpois")))

saveRDS(fuller,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_Displacement.RData")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUILDING DAMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
outFrame<-dplyr::select(outred,-c("mortality","displacement","buildDest"))
names(outFrame)[1]<-"Y"

# Function from the file CorrelateModifier.R:
fuller<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                    GLMer = "LM", ncores=40),
                 error=function(e) NA)
fuller$model<-"LM"
saveRDS(fuller,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDam_LM.RData")

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "pois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDam_pois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="pois")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "lognorm", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDam_lognorm.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="lognorm")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "HurdlePois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDam_HurdlePois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdlePois")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "HurdleNegBin", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDam_HurdleNegBin.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdleNegBin")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=3,intercept=T,fn="+",nlim=12,
                                         GLMer = "ZInegbin", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDam_ZInegbin.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZInegbin")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "ZIpois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDam_ZIpois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZIpois")))

saveRDS(fuller,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDam.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUILDING DESTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
outFrame<-dplyr::select(outred,-c("mortality","buildDam","displacement"))
names(outFrame)[1]<-"Y"

# Function from the file CorrelateModifier.R:
fuller<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                    GLMer = "LM", ncores=40),
                 error=function(e) NA)
fuller$model<-"LM"
saveRDS(fuller,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDest_LM.RData")

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "pois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDest_pois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="pois")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "lognorm", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDest_lognorm.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="lognorm")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "HurdlePois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDest_HurdlePois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdlePois")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "HurdleNegBin", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDest_HurdleNegBin.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdleNegBin")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "ZInegbin", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDest_ZInegbin.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZInegbin")))

predictions<-tryCatch(LMFeatureSelection(outFrame,Nb=15,intercept=T,fn="+",nlim=12,
                                         GLMer = "ZIpois", ncores=40),
                      error=function(e) NA)
saveRDS(predictions,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDest_ZIpois.RData")

if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZIpois")))

saveRDS(fuller,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_BuildDest.RData")



# Analyse the results:
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/"); filez<-filez[!filez%in%c("GLM_Mortality.RData","GLM_BuildDam.RData","GLM_BuildDest.RData","GLM_Displacement.RData","InputDataGLM.RData")]
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"GLM_",simplify = T)[,2]; namerz<-namerz[!namerz%in%c("Mortality","BuildDam","BuildDest","Displacement","")]

predictions<-data.frame()
for(i in 1:length(filez)) {
  if(filez[i]=="InputDataGLM.RData") next
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/",filez[i]))
  if("model"%in%colnames(tmp)) tmp%<>%dplyr::select(-"model")
  predictions%<>%rbind(cbind(tmp,data.frame(model=namerz[i])))
}

tmp<-str_split(predictions$model,"_",simplify = T)
predictions$impact<-tmp[,1]
predictions$algo<-tmp[,2]
predictions$model<-NULL
table(predictions$impact)
table(predictions$algo)

predictions%>%arrange(StandErr)%>%group_by(impact)%>%slice(1:5)
predictions%>%arrange(BIC)%>%group_by(impact)%>%slice(1:5)

predictions%>%group_by(impact,algo)%>%
  summarise(minnie=min(StandErr),BIC=BIC[which.min(StandErr)])%>%
  ggplot(aes(minnie,BIC))+geom_point(aes(colour=algo,shape=impact),size=3)+
  scale_y_log10()+scale_x_log10()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MULTIVARIATE MODELS @%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
allimps<-c("mortality","displacement","buildDam","buildDest")

imps<-c("mortality","displacement")
predictionsMV<-tryCatch(LMFeatureSelection(dplyr::select(outred,-allimps[!allimps%in%imps]),
                                           Nb=15,intercept=T,fn="+",nlim=12,ncores=40,
                                           mvm = imps,GLMer = "lognorm"),
                      error=function(e) NA)
saveRDS(predictionsMV,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/MVGLM_MortDisp_lognorm.RData")

imps<-c("buildDam","buildDest")
predictionsMV<-tryCatch(LMFeatureSelection(dplyr::select(outred,-allimps[!allimps%in%imps]),
                                           Nb=15,intercept=T,fn="+",nlim=12,ncores=40,
                                           mvm = imps,GLMer = "lognorm"),
                        error=function(e) NA)
saveRDS(predictionsMV,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/MVGLM_buildDambuildDest_lognorm.RData")

# ALL!
predictionsMV<-tryCatch(LMFeatureSelection(outred,
                                           Nb=15,intercept=T,fn="+",nlim=12,ncores=40,
                                           mvm = allimps,GLMer = "lognorm"),
                        error=function(e) NA)
saveRDS(predictionsMV,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/MVGLM_MortDispbuildDambuildDest_lognorm.RData")

# Linear models

imps<-c("mortality","displacement")
predictionsMV<-tryCatch(LMFeatureSelection(dplyr::select(outred,-allimps[!allimps%in%imps]),
                                           Nb=15,intercept=T,fn="+",nlim=12,ncores=40,
                                           mvm = imps,GLMer = "LM"),
                        error=function(e) NA)
saveRDS(predictionsMV,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/MVGLM_MortDisp_LM.RData")

imps<-c("buildDam","buildDest")
predictionsMV<-tryCatch(LMFeatureSelection(dplyr::select(outred,-allimps[!allimps%in%imps]),
                                           Nb=15,intercept=T,fn="+",nlim=12,ncores=40,
                                           mvm = imps,GLMer = "LM"),
                        error=function(e) NA)
saveRDS(predictionsMV,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/MVGLM_buildDambuildDest_LM.RData")

# ALL!
predictionsMV<-tryCatch(LMFeatureSelection(outred,
                                           Nb=15,intercept=T,fn="+",nlim=12,ncores=40,
                                           mvm = allimps,GLMer = "LM"),
                        error=function(e) NA)
saveRDS(predictionsMV,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/MVGLM_MortDispbuildDambuildDest_LM.RData")

# Let's have a look! :)
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/")
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"GLM_",simplify = T)[,2]

predictionsMV<-data.frame()
for(i in 1:length(filez)) {
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/",filez[i]))
  if("model"%in%colnames(tmp)) tmp%<>%dplyr::select(-"model")
  tmp[,allimps[!allimps%in%colnames(tmp)]]<-NA
  predictionsMV%<>%rbind(cbind(tmp,data.frame(model=namerz[i])))
}
predictionsMV$Cost<-apply(predictionsMV[,2:5],1,prod,na.rm=T)

tmp<-str_split(predictionsMV$model,"_",simplify = T)
predictionsMV$impact<-tmp[,1]
predictionsMV$algo<-tmp[,2]
predictionsMV$model<-NULL

predictionsMV%>%arrange(Cost)%>%dplyr::select(-allimps)%>%
  group_by(algo,impact)%>%slice(1:5)%>%View()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

costie<-function(data, lev = NULL, model = NULL) {
  out<-mean(abs(data$obs-data$pred)/(data$obs+1),na.rm = T)
  names(out)<-"RelativeAbs"
  return(out)
}
  
train_control <- caret::trainControl(method="repeatedcv", number=10, repeats=5,
                                     search = "random",
                                     summaryFunction=costie)


allimps<-c("mortality","displacement","buildDam","buildDest")

parallelML<-function(algo,impact) {
  
  print(paste0("Working on ",algo," for impact=",impact))  
  
  outFrame<-dplyr::select(outred,-allimps[impact!=allimps])
  names(outFrame)[1]<-"Y"
  outFrame$Y<-log(outFrame$Y+10)
  outFrame%<>%na.omit()
   
  modeler<-caret::train(Y~., data = outFrame, method = algo, metric="RelativeAbs",
                        tuneLength = 12, trControl = train_control,linout = TRUE,
                        preProcess = c("center","scale"))
  
  print(modeler$results)
  
  out<-cbind(dplyr::select(filter(modeler$results,RelativeAbs==min(RelativeAbs)),RelativeAbs),
          t(as.data.frame((t(as.data.frame(varImp(modeler, scale=FALSE)$importance))[1,]))))
  rownames(out)<-NULL
    
  saveRDS(out,paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/ML_",algo,"_",impact,".RData"))
  
  return(out)
  
}

tabmod<-getModelInfo()
carmods<-unlist(sapply(tabmod,function(x) x$type%in%"Regression"))
carmods<-data.frame(algorithm=names(carmods),regression=unname(carmods))
carmods%<>%filter(regression)%>%pull(algorithm)
# Check that we have all that we need to run each model
checkerz<-unlist(lapply(carmods,function(stst) ifelse(is.null(tryCatch(checkInstall(getModelInfo(stst)$library),error=function(e) NA)),T,F)))
carmods<-carmods[checkerz]; rm(checkerz)

minimods<-c("nnet","brnn","rf","glmnet","lm","lmStepAIC","rvmLinear","rvmRadial")
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

# Let's have a look! :)
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/")
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"ML_",simplify = T)[,2]
impact<-str_split(namerz,"_",simplify = T)[,2]
namerz<-str_split(namerz,"_",simplify = T)[,1]

predictionsML<-data.frame()
for(i in 1:length(filez)) {
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/",filez[i]))
  predictionsML%<>%rbind(cbind(data.frame(model=namerz[i],impact=impact[i]),tmp))
}
predictionsML

# Building damage assessment - classification with spatial element using kriging

# CNNs and maybe some others (RBF-NN, ResNet?)

# Then do multivariate model for displacement and mortality for CNNs

# Then do multivariate model for displacement and mortality for CNNs






