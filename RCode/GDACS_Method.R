# Extract Environment Variables
dir<-directory<-"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/";setwd(directory); packred<-T
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
source('RCode/ODDobj.R')
library(boot)
library(MASS)
install.packages("pscl")
library(pscl)

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
          Exp5=sum(ODDy@data$Population[ODDy@data$hazMax>=5 & indies],na.rm = T),
          Exp5.5=sum(ODDy@data$Population[ODDy@data$hazMax>=5.5 & indies],na.rm = T),
          Exp6=sum(ODDy@data$Population[ODDy@data$hazMax>=6 & indies],na.rm = T),
          Exp6.5=sum(ODDy@data$Population[ODDy@data$hazMax>=6.5 & indies],na.rm = T),
          Exp7=sum(ODDy@data$Population[ODDy@data$hazMax>=7 & indies],na.rm = T),
          Exp7.5=sum(ODDy@data$Population[ODDy@data$hazMax>=7.5 & indies],na.rm = T),
          Exp8=sum(ODDy@data$Population[ODDy@data$hazMax>=8 & indies],na.rm = T),
          Exp8.5=sum(ODDy@data$Population[ODDy@data$hazMax>=8.5 & indies],na.rm = T),
          Exp9=sum(ODDy@data$Population[ODDy@data$hazMax>=9 & indies],na.rm = T)
        ),WID,vulny)
      ) 
    }
  }
  print(paste0("Finished EQ: ",fff))
}

out$time<-as.numeric(out$date-min(out$date))

saveRDS(out,"./IIDIPUS_Results/GLM_Models/InputDataGLM.RData")

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


#$$$$$$$$$$$$$$$$$$$ TO DO $$$$$$$$$$$$$$$$$$$$$$#
# for each impact

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MORTALITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
outMort<-dplyr::select(outred,-c("displacement","buildDam","buildDest"))
names(outMort)[1]<-"Y"

# Function from the file CorrelateModifier.R:
fuller<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "LM", ncores=40),
                      error=function(e) NA)
fuller$model<-"LM"

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "pois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="pois")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "lognorm", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="lognorm")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "HurdlePois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdlePois")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "HurdleNegBin", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdleNegBin")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "ZInegbin", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZInegbin")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "ZIpois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZIpois")))

saveRDS(predictions,"./IIDIPUS_Results/GLM_Models/GLM_Mortality.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISPLACEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
outMort<-dplyr::select(outred,-c("mortality","buildDam","buildDest"))
names(outMort)[1]<-"Y"

# Function from the file CorrelateModifier.R:
fuller<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                    GLMer = "LM", ncores=40),
                 error=function(e) NA)
fuller$model<-"LM"

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "pois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="pois")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "lognorm", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="lognorm")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "HurdlePois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdlePois")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "HurdleNegBin", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdleNegBin")))

# predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
#                                          GLMer = "ZInegbin", ncores=40),
#                       error=function(e) NA)
# if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZInegbin")))
# 
# predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
#                                          GLMer = "ZIpois", ncores=40),
#                       error=function(e) NA)
# if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZIpois")))

saveRDS(predictions,"./IIDIPUS_Results/GLM_Models/GLM_Displacement.RData")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUILDING DAMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
outMort<-dplyr::select(outred,-c("mortality","displacement","buildDest"))
names(outMort)[1]<-"Y"

# Function from the file CorrelateModifier.R:
fuller<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                    GLMer = "LM", ncores=40),
                 error=function(e) NA)
fuller$model<-"LM"

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "pois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="pois")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "lognorm", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="lognorm")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "HurdlePois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdlePois")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "HurdleNegBin", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdleNegBin")))

# predictions<-tryCatch(LMFeatureSelection(outMort,Nb=3,intercept=T,fn="+",nlim=6,
#                                          GLMer = "ZInegbin", ncores=40),
#                       error=function(e) NA)
# if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZInegbin")))
# 
# predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
#                                          GLMer = "ZIpois", ncores=40),
#                       error=function(e) NA)
# if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZIpois")))

saveRDS(predictions,"./IIDIPUS_Results/GLM_Models/GLM_BuildDam.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUILDING DESTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
outMort<-dplyr::select(outred,-c("mortality","buildDam","displacement"))
names(outMort)[1]<-"Y"

# Function from the file CorrelateModifier.R:
fuller<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                    GLMer = "LM", ncores=40),
                 error=function(e) NA)
fuller$model<-"LM"

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "pois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="pois")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "lognorm", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="lognorm")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "HurdlePois", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdlePois")))

predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
                                         GLMer = "HurdleNegBin", ncores=40),
                      error=function(e) NA)
if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="HurdleNegBin")))

# predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
#                                          GLMer = "ZInegbin", ncores=40),
#                       error=function(e) NA)
# if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZInegbin")))
# 
# predictions<-tryCatch(LMFeatureSelection(outMort,Nb=15,intercept=T,fn="+",nlim=6,
#                                          GLMer = "ZIpois", ncores=40),
#                       error=function(e) NA)
# if(!all(is.na(predictions))) fuller%<>%rbind(cbind(predictions,data.frame(model="ZIpois")))

saveRDS(predictions,"./IIDIPUS_Results/GLM_Models/GLM_BuildDest.RData")


print("MAKE SURE TO REPLACE THE LAST ROUND OF model=ZIpois IN GLM_MORTALITY")


# Building damage assessment - classification with and without spatial element (GLM, SVM, kriging)

# Graph NN for building damage data?

# CNNs and maybe some others (RBF-NN, ResNet?)

# Then do multivariate model for displacement and mortality for both CNNs and GLMs








