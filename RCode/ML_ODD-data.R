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

# install.packages("ggcorrplot")

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

ExtractGLMresults<-function(algo,impact,othimps=NULL){
  
  print(paste0("Working on ",algo," for impact=",impact))  
  
  if(is.null(othimps)){
    outFrame<-dplyr::select(outred,-allimps[impact!=allimps])  
  } else {
    outFrame<-dplyr::select(outred,-allimps[!allimps%in%c(impact,othimps)])
    ind<-which(colnames(outFrame)==impact)
    outFrame<-outFrame[,c(ind,(1:ncol(outFrame))[-ind])]
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

parallelML<-function(algo,impact,modRet=F) {
  
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
  
  if(modRet) return(modeler)
  
  print(modeler$results)
  
  out<-cbind(dplyr::select(filter(modeler$results,RelativeAbs==min(RelativeAbs)),RelativeAbs,RelativeAbsSD),
          t(as.data.frame((t(as.data.frame(varImp(modeler, useModel=F, nonpara=F, scale=FALSE)$importance))[1,]))))
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
  if(filez[i]=="ML_brnn_mortality.RData") next
  tmp<-readRDS(paste0("./IIDIPUS_Results/SpatialPolygons_ML-GLM/NoSpace_ML_models/",filez[i]))
  predictionsML%<>%rbind(cbind(data.frame(model=namerz[i],impact=impact[i]),tmp))
}

table(predictionsML$model)
table(predictionsML$impact)

View(predictionsML)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%% COMPARE UNIVARIATE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

fuller<-rbind(transmute(predictions,RelativeAbsDiff=StandErr,RelativeAbsDiffSD=StandErrSD,
                        impact=impact,algo=algo, equation=eqn),
              transmute(predictionsML,RelativeAbsDiff=RelativeAbs,RelativeAbsDiffSD=RelativeAbsSD,
                        impact=impact,algo=model,equation="Y ~ maxHaz + ExpSchYrs + LifeExp + GNIc + Vs30 + EQFreq + hazSD + time + ExpDim1 + ExpDim2 + WIDDim1 + WIDDim2"))
View(fuller)     

predictions%<>%group_by(impact)%>%mutate(BestDiff=StandErr-min(StandErr,na.rm = T))

p<-predictions%>%filter(StandErr<0.6)%>%ggplot(aes(StandErr,group=algo)) +
  geom_density(aes(fill=algo),alpha=0.4)+scale_x_log10()+
  xlab("Relative Absolute Log-Difference Errors - Log Scale ")+ylab("Density")+labs(fill="Model")+
  facet_wrap(.~impact,scales = "free");p
ggsave("GLM_Errors.eps",p,path="./Plots/IIDIPUS_Results/GLM-ML_Work/",width=10,height=5,device = grDevices::cairo_ps)  


fuller%<>%na.omit()

allUV<-fuller%>%arrange(RelativeAbsDiff)%>%group_by(impact,algo)%>%slice(1)%>%
  dplyr::select(impact,algo,RelativeAbsDiff,RelativeAbsDiffSD)%>%filter(algo!="lm")

# colnames(allUV)<-c("Impact","Model","Avg. MADL Value", "S.D. MADL Value")
pal <- c(
  "mortality" = scales::hue_pal()(4)[1],
  "displacement" = scales::hue_pal()(4)[2], 
  "buildDam" = scales::hue_pal()(4)[3], 
  "buildDest" = scales::hue_pal()(4)[4],
  "average" = "black"
)

p<-allUV%>%ggplot(aes(algo,RelativeAbsDiff,group=impact))+
  # geom_point(aes(colour=impact,shape=impact),size=4) +
  geom_pointrange(aes(ymin=RelativeAbsDiff-RelativeAbsDiffSD,
                    ymax=RelativeAbsDiff+RelativeAbsDiffSD,colour=impact,
                    shape=impact),
                size=1)+theme(axis.text.x = element_text(angle = 90))+
  scale_y_log10()+scale_shape_manual(values=15:18,breaks=allimps)+
  xlab("Model") + ylab("Mean Absolute Deviation of Logs (MADL)")+
  scale_colour_manual(values = pal,limits = names(pal))+
  labs(colour="Impact Type",shape="Impact Type");p
ggsave("UV_MADL.eps",p,path="./Plots/IIDIPUS_Results/",width=10,height=4.,device = grDevices::cairo_ps)  

varimport<-predictionsML%>%group_by(impact)%>%summarise_at(colnames(predictionsML)[-c(1:4)],mean)
varimport%<>%rbind(as.data.frame(cbind(impact="average",as.data.frame(t(colMeans(varimport[,-1]))))))
varimport[,-1]<-100*varimport[,-1]/rowSums(as.matrix(varimport[,-1]))
varimport<-varimport[,rev(c(1,order(as.numeric(varimport[varimport$impact=="average",-1]))+1))]

p<-varimport%>%reshape2::melt("impact")%>%
  ggplot(aes(variable,value,group=impact))+
  geom_point(aes(colour=impact,shape=impact),size=3)+
  geom_line(aes(colour=impact),alpha=0.5,linetype="dotdash")+
  scale_shape_manual(values=15:19,breaks=allimps)+
  scale_colour_manual(values = pal,limits = names(pal))+
  xlab("Model Covariate") + ylab("Feature Importance [%]")+
  labs(colour="Impact Type",shape="Impact Type")+
  theme(axis.text.x = element_text(angle = 90));p

ggsave("UV_FeatImp.eps",p,path="./Plots/IIDIPUS_Results/",width=10,height=4.,device = grDevices::cairo_ps)  
ggsave("UV_FeatImp.png",p,path="./Plots/IIDIPUS_Results/",width=10,height=4.)  

# tabUV<-lapply(unique(allUV$Impact), function(imp){
#   print(xtable::xtable(filter(allUV,Impact==imp),
#                        paste0("MADL distances for the univariate response, multiple regression models, for ",imp," observational data only")), row.names = FALSE)
# })


#  HYPOTHESIS TESTING ON BEST PERFORMING MODEL
# fuller%<>%group_by(impact)%>%mutate(nn=sum(!is.na(outred[,unique(impact)])))
# fuller%<>%group_by(impact)%>%mutate(BestDiff=dnorm(RelativeAbsDiff-min(RelativeAbsDiff),0,
#                                                       sqrt(RelativeAbsDiffSD^2+RelativeAbsDiffSD[which.min(RelativeAbsDiff)]^2)/sqrt(nn))<0.05)%>%
#   dplyr::select(-nn)
  
# fuller%<>%group_by(impact)%>%mutate(minSD=RelativeAbsDiffSD[which.min(RelativeAbsDiff)]^2/unique(nn),
#                                     thisSD=RelativeAbsDiffSD^2/nn)%>%
#   mutate(df=(thisSD+minSD)^2/( (thisSD^2+minSD^2) / (nn-1) ))
  
# fuller%<>%group_by(impact)%>%mutate(withBest=dt(RelativeAbsDiff-min(RelativeAbsDiff),df))

# topOeach<-fuller%>%group_by(impact,algo)%>%arrange(RelativeAbsDiff)%>%slice(1)%>%dplyr::select(3,4,1,2)
# 
# topOeach%>%ggplot(aes(algo,RelativeAbsDiff,group=impact)) +
#   geom_point(aes(colour=impact,shape=impact),size=4)+
#   geom_errorbar(aes(ymin=RelativeAbsDiff-RelativeAbsDiffSD, 
#                     ymax=RelativeAbsDiff+RelativeAbsDiffSD,colour=impact),
#                 width=.4,alpha=0.7)+
#   scale_y_log10()+scale_shape_manual(values=15:18,breaks=allimps)

# For each of the top models, calculate the feature importance using vip package
# See here for more info: https://cran.r-project.org/web/packages/vip/vignettes/vip-introduction.pdf
# Then compare the values between the top 10-models (over all model types) for each impact

# FOR ALL MODELS NOT STAT. SIG. DIFF FROM BEST MODEL, PER GLM, MEASURE THE MODEL-AGNOSTIC VIP

GLMmods<-predictions%>%arrange(StandErr)%>%group_by(impact,algo)%>%slice(1)%>%
  dplyr::select(eqn,StandErr,StandErrSD,impact,algo)

GLMmods%<>%filter(algo=="lognorm")

resultsUV<-do.call(rbind,lapply(1:nrow(GLMmods),function(i){
  inpy<-as.data.frame(GLMmods[i,])
  inpy$impact
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
  vippy%<>%dplyr::select(sort(colnames(vippy)))%>%cbind(inpy)%>%dplyr::select(-eqn)%>%
    dplyr::select(impact,algo,StandErr,StandErrSD,everything())
  # output it all!
  return(vippy)
}))

# Order the variables by their importance, weighted by each models error value
varimport<--apply(resultsUV[,-(1:4)],2,function(x) weighted.mean(x,1/resultsUV$StandErr))
# Reorder the dataframe
resultsUV%<>%dplyr::select(algo,impact,StandErr,StandErrSD,names(sort(varimport)))
colnames(resultsUV)<-colnames(predictionsML)
rownames(resultsUV)<-NULL

resultsUV%<>%rbind(predictionsML)

tabUV<-lapply(unique(resultsUV$impact), function(imp){
  xtable::xtable(filter(resultsUV,impact==imp),"Feature importance ranking measure (FIRM), scaled as a percentage, for the best performing model-formulation for each GLM or ML model. When applying the FIRM, Individual Conditional Expectation (ICE) curves are used.")})

for(i in 1:length(tabUV)) print(tabUV[[i]],include.rownames=F)
# Melt the dataframe to plot
# resultsUV%<>%reshape2::melt(id.vars=1:4)
# # Make a weighting variable for the alpha plotting
# resultsUV$normError<-1/resultsUV$StandErr
# resultsUV%<>%group_by(impact)%>%mutate(normError=normError/max(normError))%>%ungroup()
# # Plot
# resultsUV%>%ggplot(aes(as.factor(variable),value,group=algo))+geom_point(aes(colour=algo,alpha=normError))+
#   facet_wrap(.~impact,nrow = 4)

#  HYPOTHESIS TESTING ON BEST PERFORMING MODEL
resultsUV%<>%group_by(impact)%>%mutate(nn=sum(!is.na(outred[,unique(impact)])))
resultsUV%<>%group_by(impact)%>%mutate(withBest=dnorm(RelativeAbs-min(RelativeAbs),0,
   sqrt(RelativeAbsSD^2+RelativeAbsSD[which.min(RelativeAbs)]^2)/sqrt(nn))<0.05)%>%
  dplyr::select(-nn)

resultsUV%<>%group_by(impact)%>%mutate(withBest=dnorm(RelativeAbs-min(RelativeAbs),0,
                                                      sqrt(RelativeAbsSD^2+RelativeAbsSD[which.min(RelativeAbs)]^2)/sqrt(nn))<0.05)%>%
  dplyr::select(-nn)

tabUV<-lapply(unique(resultsUV$impact), function(imp){
  xtable::xtable(filter(resultsUV,impact==imp & withBest),"Feature importance ranking measure (FIRM), scaled as a percentage, for the \'best-performing\' models. The definition of \'best performing\' here is defined as any model that had a cost that was not statistically significantly different from the overall best performing model for each specific impact. When applying the FIRM, Individual Conditional Expectation (ICE) curves are used.")})

for(i in 1:length(tabUV)) print(tabUV[[i]],include.rownames=F)

# The models that made it as one of the best performing, per impact:
resultsUV[resultsUV$withBest,]%>%group_by(impact)%>%reframe(models=unique(model))%>%View()
resultsUV[resultsUV$withBest,]%>%group_by(impact)%>%reframe(models=length(model))
nrow(resultsUV)

# Order the variables by their importance, weighted by each models error value
varimport<-do.call(rbind,lapply(allimps,function(imps) {tmp<-filter(resultsUV,impact==imps);apply(tmp[tmp$withBest,-c(1:4,ncol(tmp),(ncol(tmp)-1))],2,function(x) weighted.mean(x,1/tmp$RelativeAbs[tmp$withBest]))}))
varimport<-100*varimport/rowSums(varimport)
varimport%<>%rbind(as.numeric(colMeans(varimport)))
rownames(varimport)<-c(allimps,"average")
varimport<-100*varimport/rowSums(varimport)
varimport%<>%t()
varimport%<>%as.data.frame(row.names = rownames(varimport))%>%mutate(Covariate=rownames(varimport))

varimport%<>%reshape2::melt(id.vars=6)
colnames(varimport)[2]<-"impact"

pal <- c(
  "mortality" = "red",
  "displacement" = "blue", 
  "buildDam" = "forestgreen", 
  "buildDest" = "purple",
  "average" = "black"
)
varimport%>%ggplot(aes(Covariate,value,group=impact))+geom_point(aes(colour=impact,shape=impact),size=2)+
  geom_line(aes(colour=impact),alpha=0.25)+scale_colour_manual(values = pal,limits = names(pal))







colnames(resultsUV)[2:4]<-c("model","RelativeAbs","RelativeAbsSD")

predictionsML%<>%rbind(resultsUV)

varimp<-predictionsML%>%arrange(RelativeAbs)%>%group_by(model)%>%
  summarise(Cost=prod(RelativeAbs),
            maxHaz=mean(maxHaz),
            ExpSchYrs=mean(ExpSchYrs),
            LifeExp=mean(LifeExp),
            GNIc=mean(GNIc),
            Vs30=mean(Vs30),
            EQFreq=mean(EQFreq),
            hazSD=mean(hazSD),
            time=mean(time),
            ExpDim1=mean(ExpDim1),
            ExpDim2=mean(ExpDim2),
            WIDDim1=mean(WIDDim1),
            WIDDim2=mean(WIDDim2))

varimp[3:ncol(varimp)]<-100*varimp[3:ncol(varimp)]/
rowSums(varimp[3:ncol(varimp)])
# 
varimp%>%filter(model%in%c("rf","svmRadial","svmPoly","lognorm","nnet"))%>%xtable::xtable()%>%print(row.names = FALSE)







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
  out<-LMFeatureSelection(outFrame,Nb=1,intercept=T,fn="+",nlim=(ncol(outFrame)-1),
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
MVvars<-str_split(str_split(str_split((str_split(MVGLMmods$eqn,"~",simplify = T)[,1]),"\\(",simplify = T)[,2],"\\)",simplify = T)[,1],",")

MVGLMmods%<>%filter(algo=="lognorm")

resultsMV<-lapply(1:nrow(MVGLMmods),function(i){
  inpy<-as.data.frame(MVGLMmods[i,])
  if(inpy$algo=="LM") return(NA)
  # Remove NAs
  outFrame<-outred%>%na.omit()%>%dplyr::select(-allimps[allimps%in%inpy$impact])
  # Make sure log-normal model is taken into account
  if(inpy$algo=="lognorm") outFrame[,MVvars[[i]]]<-log(outFrame[,MVvars[[i]]]+10)
  # Remove the variable Event after weighting is calculated
  outFrame%<>%dplyr::select(-Event)
  # Extract the trained model
  outmod<-lm(formula = as.formula(inpy$eqn),data = outFrame)
  # Which variables to include
  vars<-names(outmod$model)[-1];vars<-vars[-length(vars)]
  # Feature importance (model-agnostic) calculation
  vippy<-as.data.frame(vip::vi(outmod,method="firm",scale=T,ice=T,
                               feature_names=vars))
  # Add missing columns
  missies<-colnames(outFrame)[!colnames(outFrame)%in%c(vippy$Variable,inpy$impact)]
  vippy%<>%rbind(data.frame(Variable=missies,Importance=0))
  # Reformulate to column form 
  rownames(vippy)<-vippy$Variable; vippy%<>%dplyr::select(Importance)%>%t()%>%as.data.frame()
  # colnames in alphabetical order and add the model information to the data frame
  vippy%<>%dplyr::select(sort(colnames(vippy)))%>%cbind(dplyr::select(inpy,Cost,impact,algo))
  # output it all!
  return(list(vip=vippy,covvy=vcov(outmod)))
})

covvy<-resultsMV[[1]]$covvy;covvy<-covvy-min(covvy);covvy<-covvy/max(covvy);covvy<-2*covvy - 1

ggcorrplot::ggcorrplot(covvy, type = "lower",
                       lab = TRUE)
ncl<-(ncol(covvy)/2L)
ccrr<-data.frame()
for(i in 1:ncl) {
  ccrr%<>%rbind(data.frame(corii=covvy[i,i+ncl],name=str_split(colnames(covvy)[i],":",simplify = T)[,2]))
}

# Compare StandErr for lognormal against LM to show why you will only use lognormal afterwards
# Mention which covariates were in the highest-performing models
# Run the model with all covariates for the lognormal model
# Look at which covariates are stat. sig. for both models
# Do the same for the highest performing models
# For which covariates do you print out the covariance?
#     those that are stat.sig in both?
#     those that are present in the highest performing MV models?
#     those that are present in the highest performing UV models?

# I think most likely the last one. Make sure to have done the analysis on the other 









# TO DO:
# 1) Table of top-5 models, per impact
#    FOR BD AND ODD MODELS
# 2) Plot comparison of StandErr between GLM and ML models, with error bar
# 4) Application of GLM and ML models on the MV impact data to make the comparison
# 5) 


























# Building damage assessment - classification with spatial element using kriging

# CNNs and maybe some others (RBF-NN, ResNet?)

# Then do multivariate model for displacement and mortality for CNNs

# Then do multivariate model for displacement and mortality for CNNs






