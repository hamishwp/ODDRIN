folder<-"./IIDIPUS_Input/BDobjects/"
filez<-list.files(folder)
RWIfolder<-"./Demography_Data/RWI/"
haz<-"EQ"
if(haz=="EQ") funcyfun<-exp else funcyfun<-returnX

BDs<-data.frame()
for(fff in filez){
 
  BDy<-readRDS(paste0(folder,fff))
  if(length(grep(names(BDy@data),pattern = "hazMean",value = T))==1) hazMax<-BDy$hazMean1
  else hazMax<-apply(BDy@data[,grep(names(BDy@data),pattern = "hazMean",value = T)],1,max,na.rm=T)
  
  iso<-names(which.max(table(BDy$ISO3C)))
  
  fRWI<-grep(list.files(RWIfolder),pattern = str_to_lower(iso),value = T)
  if(length(fRWI)==0) next
 
  print(fff)
   
  RWI<-read.csv(paste0(RWIfolder,fRWI))
  names(RWI)<-c("Latitude","Longitude","Intensity","error")

  # OPTION TO KRIG VALUES INSTEAD:  
  insideBBOX<-RWI$Longitude>1.1*BDy@bbox[1] &
    RWI$Latitude>1.1*BDy@bbox[2] &
    RWI$Longitude<1.1*BDy@bbox[3] &
    RWI$Latitude<1.1*BDy@bbox[4]
  
  if(sum(insideBBOX)<1000) insideBBOX<-T
  
  print(paste0("Number of RWI values used = ",sum(insideBBOX)," / ",nrow(RWI)))
  
  RWI<-RWI[insideBBOX,]
  
  interpRWI<-apply(as.data.frame(BDy)[,c("Longitude","Latitude")],1,function(x){
    disty<-geosphere::distHaversine(x,RWI[,c("Longitude","Latitude")])
    return(list(Exp=weighted.mean(RWI$Intensity,1./(disty^2),na.rm=TRUE),
                Error=weighted.mean(RWI$error,1./(disty^2),na.rm=TRUE)))
  })
  
  tmpRWI<-data.frame(Exp=rep(NA,length(interpRWI)),Error=rep(NA,length(interpRWI)))
  for(i in 1:length(interpRWI)){
    tmpRWI$Exp[i]<-unlist(interpRWI[[i]][1])
    tmpRWI$Error[i]<-unlist(interpRWI[[i]][2])
  }
  
  # interpRWI<-KrigMeUp(RWI[insideBBOX,],as.data.frame(BDy)[,c("Longitude","Latitude")])
  # Place hard bounding limits on the values interpolated from the GPR
  # interpRWI[interpRWI>max(RWI,na.rm = T)]<-max(RWI,na.rm = T)
  # interpRWI[interpRWI<min(RWI,na.rm = T)]<-min(RWI,na.rm = T)
  
  BDs%<>%rbind(
    data.frame(
      Event=fff,
      iso3=iso,
      date=BDy@hazdates[1],
      Grading=BDy@data$grading,
      PopDens=BDy@data$Population,
      HazardIntensity=funcyfun(hazMax-BDy@I0),
      GDP=BDy@data$GDP,
      RWI=tmpRWI$Exp,
      RWIerror=tmpRWI$Error
    )
  )
   
}

tmp<-dplyr::select(BDs,c(Grading,PopDens,HazardIntensity,GDP,RWI))  
tmp$PopDens<-log(tmp$PopDens)
tmp$GDP<-log(tmp$GDP)
tmp%<>%filter(!is.na(tmp$Grading))

n<-15

damage<-sampleBDdamage(tmp$Grading,n = n)
damage<-logit(damage)

donnees<-data.frame()
for (j in 1:n){
  donnees<-tmp%>%mutate(damage=damage[j,])%>%rbind(donnees)
}
donnees%<>%dplyr::select(-c(Grading))

# monmod<-monreg::monreg(donnees, tmp$damage)

# library(tfprobability)

poly<-donnees%>%transmute(Intensity=damage,Longitude=HazardIntensity,Latitude=RWI)
KrigMeUp()





# install.packages("neuralnet")
library(neuralnet)
library(plyr) 
set.seed(450)
cv.error <- NULL
k <- 10
maxs <- apply(donnees, 2, max) 
mins <- apply(donnees, 2, min)
scaled <- as.data.frame(scale(donnees, center = mins, scale = maxs - mins))
n <- names(scaled[,])
f <- as.formula(paste("damage ~", paste(n[!n %in% "damage"], collapse = " + ")))
pbar <- create_progress_bar('text')
pbar$init(k)
for(i in 1:k){
  index <- sample(1:nrow(donnees),round(0.1*nrow(donnees)))
  train.cv <- scaled[index,]
  test.cv <- scaled[-index,]
  nn <- neuralnet(f,data=train.cv,hidden=c(5,2),linear.output=T)   
  pr.nn <- compute(nn,test.cv)
  pr.nn <- pr.nn$net.result*(max(donnees$damage)-min(donnees$damage))+min(donnees$damage)   
  test.cv.r <- (test.cv$damage)*(max(donnees$damage)-min(donnees$damage))+min(donnees$damage)   
  cv.error[i] <- sum((test.cv.r - pr.nn)^2)/nrow(test.cv)    
  pbar$step()
}

nn<-neuralnet("damage ~ HazardIntensity + RWI",data=scaled,hidden=c(5,2),linear.output=T)  






BDs<-data.frame()
for(fff in filez){
  
  BDy<-readRDS(paste0(folder,fff))
  if(length(grep(names(BDy@data),pattern = "hazMean",value = T))==1) hazMax<-BDy$hazMean1
  else hazMax<-apply(BDy@data[,grep(names(BDy@data),pattern = "hazMean",value = T)],1,max,na.rm=T)
  
  BDs%<>%rbind(
    data.frame(
      Event=fff,
      date=BDy@hazdates[1],
      Grading=BDy@data$grading,
      PopDens=BDy@data$Population,
      HazardIntensity=funcyfun(hazMax-BDy@I0),
      GDP=BDy@data$GDP
    )
  )
  
}



























