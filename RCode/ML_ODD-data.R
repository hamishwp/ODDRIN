# Extract Environment Variables
# dir<-directory<-"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/";setwd(directory); packred<-T
# Download and install the necessary packages:
# source('RCode/GetODDPackages.R')
# source('RCode/ODDobj.R')
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
library(caret)

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
    
    modely<-function(equation,datar,responsers){
      # If we are using a lognormal model
      if(lognorm) datar[,responsers]<-log(datar[,responsers]+10)
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
    modely<-function(equation,datar,responsers="Y") {
      # Train the model on all the data in order to check the BIC
      modeler<-glm(formula = as.formula(equation),data = datar,family = gaussian(link = "identity"),weights = weights)
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
    modely<-function(equation,datar,responsers="Y") {
      # Train the model on all the data in order to check the BIC
      modeler<-glm(formula = as.formula(equation),data = datar,family = poisson(),weights = weights)
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
    modely<-function(equation,datar,responsers="Y") {
      # Convert to log scale
      datar$Y<-log(datar$Y+10)
      # Train the model on all the data in order to check the BIC
      modeler<-glm(formula = as.formula(equation),data = datar,family = gaussian(link = "identity"),weights = weights)
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
    modely<-function(equation,datar,responsers="Y") {
      # Train the model on all the data in order to check the BIC
      modeler<-pscl::hurdle(formula = as.formula(equation),data = datar,dist="poisson",weights = weights)
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
    modely<-function(equation,datar,responsers="Y") {
      # Train the model on all the data in order to check the BIC
      modeler<-pscl::hurdle(formula = as.formula(equation),data = datar,dist="negbin",weights = weights)
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
    modely<-function(equation,datar,responsers="Y") {
      # Train the model on all the data in order to check the BIC
      modeler<-pscl::zeroinfl(formula = as.formula(equation),data = datar,dist="negbin",weights = weights)
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
    modely<-function(equation,datar,responsers="Y") {
      # Train the model on all the data in order to check the BIC
      modeler<-pscl::zeroinfl(formula = as.formula(equation),data = datar,dist="poisson",weights = weights)
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
      ressies<-tryCatch(modely(eee,output,responsers),error = function(e) NA)
      if(all(is.na(ressies))) next
      
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRACT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
folder<-"./IIDIPUS_Input/IIDIPUS_Input_NMAR/ODDobjects/"
filez<-list.files(folder)

out<-data.frame()
for(fff in filez){

  if(fff=="EQ20170529IDN_136") next

  ODDy<-readRDS(paste0(folder,fff))

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
          Event=fff,
          Exp5=sum(ODDy@data$Population[ODDy@data$hazMax>=5 & indies],na.rm = T),
          Exp5.5=sum(ODDy@data$Population[ODDy@data$hazMax>=5.5 & indies],na.rm = T),
          Exp6=sum(ODDy@data$Population[ODDy@data$hazMax>=6 & indies],na.rm = T),
          Exp6.5=sum(ODDy@data$Population[ODDy@data$hazMax>=6.5 & indies],na.rm = T),
          Exp7=sum(ODDy@data$Population[ODDy@data$hazMax>=7 & indies],na.rm = T),
          Exp7.5=sum(ODDy@data$Population[ODDy@data$hazMax>=7.5 & indies],na.rm = T),
          Exp8=sum(ODDy@data$Population[ODDy@data$hazMax>=8 & indies],na.rm = T),
          Exp8.5=sum(ODDy@data$Population[ODDy@data$hazMax>=8.5 & indies],na.rm = T),
          Exp9=sum(ODDy@data$Population[ODDy@data$hazMax>=9 & indies],na.rm = T),
          maxHaz=max(ODDy@data$hazMax,na.rm = T)
        ),WID,vulny)
      )
    }
  }
  print(paste0("Finished EQ: ",fff))
}

out$time<-as.numeric(out$date-min(out$date))

saveRDS(out,"./IIDIPUS_Results/SpatialPolygons_ML-GLM/InputData_ODD.RData")

out<-readRDS("./IIDIPUS_Results/SpatialPolygons_ML-GLM/InputData_ODD.RData")

# Per impact, do the analysis (for all impacts, even without building data)
outred<-dplyr::select(out,-c("iso3","date")); rm(out)
outred[,-(1:5)]<-scale(outred[,-(1:5)])

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
print(paste0("Checking the VIF of the remaining variables: percentage that didn't pass = ",sum(multiColl::VIF(outred[,-(1:5)])>6)/ncol(outred[,-(1:5)])))
colnames(outred)

allimps<-c("mortality","displacement","buildDam","buildDest")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNIVARIATE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

ExtractGLMresults<-function(algo,impact){
  
  print(paste0("Working on ",algo," for impact=",impact))  
  
  outFrame<-dplyr::select(outred,-allimps[impact!=allimps])
  names(outFrame)[1]<-"Y"
  outFrame%<>%na.omit()
  # Make weights from the different events to make sure that no single event dominates the model parameterisation
  weights<-outFrame%>%group_by(Event)%>%summarise(www=1/length(time))%>%merge(outFrame)%>%pull(www)
  # Remove the variable Event after weighting is calculated
  outFrame%<>%dplyr::select(-Event)
  
  # Function from the file CorrelateModifier.R:
  out<-tryCatch(LMFeatureSelection(outFrame,Nb=1,intercept=T,fn="+",nlim=12,
                                   GLMer = algo, weights = weights, ncores=60),
                error=function(e) NA)
  # out$model<-algo
  saveRDS(out,paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/GLM_",algo,"_",impact,".RData"))
  
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
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/GLM_Models/"); filez<-filez[!filez%in%c("GLM_Mortality.RData","GLM_BuildDam.RData","GLM_BuildDest.RData","GLM_Displacement.RData","InputDataGLM.RData")]
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"GLM_",simplify = T)[,2]; namerz<-namerz[!namerz%in%c("Mortality","BuildDam","BuildDest","Displacement","")]

predictions<-data.frame()
for(i in 1:length(filez)) {
  if(filez[i]=="InputDataGLM.RData" | grepl(filez[i],pattern = "MVGLM")) next
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

predictions%>%arrange(StandErr)%>%group_by(impact)%>%slice(1:5)%>%View()
# predictions%>%arrange(BIC)%>%group_by(impact)%>%slice(1:5)
# 
# predictions%>%filter(algo!="lognorm")%>%group_by(impact,algo)%>%
#   summarise(minnie=min(StandErr),BIC=BIC[which.min(StandErr)])%>%
#   ggplot(aes(minnie,BIC))+geom_point(aes(colour=algo,shape=impact),size=3)+
#   scale_y_log10()+scale_x_log10()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MULTIVARIATE MODELS @%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

ExtractMVGLMresults<-function(algo,impact){
  
  print(paste0("Working on ",algo," for impacts=",impact))  
  # Make weights from the different events to make sure that no single event dominates the model parameterisation
  weights<-outred%>%group_by(Event)%>%summarise(www=1/length(time))%>%merge(outred)%>%pull(www)
  # Remove the variable Event after weighting is calculated
  outFrame%<>%dplyr::select(-Event)
  # Function from the file CorrelateModifier.R:
  out<-tryCatch(LMFeatureSelection(outred,Nb=1,intercept=T,fn="+",nlim=12,
                                   GLMer = algo, weights = weights, mvm = imps, ncores=60),
                error=function(e) NA)
  out$model<-algo
  
  saveRDS(out,paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/MVGLM_",algo,"_",impact,".RData"))
  
  return(out)
  
}

minimods<-c("LM","lognorm")

# Run it!
ODD_ML<-lapply(allimps, function(impact) {
  lapply(minimods,function(algo) {
    out<-tryCatch(ExtractGLMresults(algo,impact),error=function(e) NA)
    if(any(is.na(out))) return(NULL)
  })})

# Let's have a look! :)
filez<-list.files("./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/")
namerz<-str_split(str_split(filez,".RData",simplify = T)[,1],"GLM_",simplify = T)[,2]; namerz<-namerz[namerz!=""]

predictionsMV<-data.frame()
for(i in 1:length(filez)) {
  if(filez[i]=="InputData_BD.RData") next
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/MV_GLM_Models/",filez[i]))
  if("model"%in%colnames(tmp)) tmp%<>%dplyr::select(-"model")
  tmp[,allimps[!allimps%in%colnames(tmp)]]<-NA
  tmp<-tmp[apply(tmp,1,function(x)!all(is.na(x))),]
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

# costie<-function(data, lev = NULL, model = NULL) {
costie<-function(data) {
  out<-mean(abs(data$obs-data$pred)*data$weights/sum(data$weights),na.rm = T)
  names(out)<-"RelativeAbs"
  return(out)
}
  
train_control <- caret::trainControl(method="repeatedcv", number=10, repeats=3,
                                     search = "random",
                                     summaryFunction=costie)

allimps<-c("mortality","displacement","buildDam","buildDest")

parallelML<-function(algo,impact) {
  
  print(paste0("Working on ",algo," for impact=",impact))  
  
  outFrame<-dplyr::select(outred,-allimps[impact!=allimps])
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
  
  out<-cbind(dplyr::select(filter(modeler$results,RelativeAbs==min(RelativeAbs)),RelativeAbs),
          t(as.data.frame((t(as.data.frame(varImp(modeler, scale=FALSE)$importance))[1,]))))
  rownames(out)<-NULL
    
  saveRDS(out,paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/ML_",algo,"_",impact,".RData"))
  
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
View(predictionsML)

# Building damage assessment - classification with spatial element using kriging

# CNNs and maybe some others (RBF-NN, ResNet?)

# Then do multivariate model for displacement and mortality for CNNs

# Then do multivariate model for displacement and mortality for CNNs






