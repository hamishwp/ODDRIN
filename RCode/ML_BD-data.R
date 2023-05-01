# install.packages("caret",repos = "http://cran.r-project.org", 
#                  dependencies = c("Depends", "Imports","Suggests"))
# devtools::install_github("souravc83/fastAdaboost")
# devtools::install_github("davpinto/fastknn")
# install.packages("ROCit")

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
library(fastknn)
library(caTools)
library(gstat)
library(sf)
library(sp)
source('RCode/BDobj.R')

# folder<-"./IIDIPUS_Input/IIDIPUS_Input_NMAR/BDobjects/"
# filez<-list.files(folder)
# haz<-"EQ"
# if(haz=="EQ") funcyfun<-exp else funcyfun<-returnX
# 
# BDs<-data.frame()
# for(fff in filez){
# 
#   # if(fff%in%c("EQ20180225PNG_68","EQ20190925IDN_95","EQ20201030TUR_88")) next
# 
#   BDy<-readRDS(paste0(folder,fff))
#   if(length(BDy@data[,-grep(names(BDy@data),pattern = "haz",value = F)])<9) next
# 
#   if(length(grep(names(BDy@data),pattern = "hazMean",value = T))==1) {
#     hazMax<-BDy$hazMean1
#     hazSD<-BDy$hazSD1
#   } else {
#     hazMax<-apply(BDy@data[,grep(names(BDy@data),pattern = "hazMean",value = T)],1,max,na.rm=T)
#     hazSD<-apply(BDy@data[,grep(names(BDy@data),pattern = "hazSD",value = T)],1,median,na.rm=T)
#   }
#   BDy@data%<>%dplyr::select(-c(grep(names(BDy@data),pattern = "hazMean",value = T),
#                                grep(names(BDy@data),pattern = "hazSD",value = T),
#                                grep(names(BDy@data),pattern = "itude",value = T),
#                                grep(names(BDy@data),pattern = "nBuildings",value = T)))%>%
#     mutate(Event=fff,date=BDy@hazdates[1],hazMax=funcyfun(hazMax-BDy@I0),hazSD=hazSD)
# 
#   BDy@data$Longitude<-BDy@coords[,1]
#   BDy@data$Latitude<-BDy@coords[,2]
# 
#   BDs%<>%rbind(BDy@data)
# 
# }
# rm(BDy)
# 
# BDs$time<-as.numeric(BDs$date-min(BDs$date))
# 
# BDs%<>%dplyr::select(-c("date","Confidence","ISO3C"))
# BDs$Population<-log(BDs$Population+1)
# BDs$GNIc<-log(BDs$GNIc)
# BDs%<>%filter(!as.character(grading)%in%c("possible") & apply(BDs,1,function(x) !any(is.na(x))))
# BDs$Damage<-abs(1-as.integer(BDs$grading=="notaffected"))
# # BDs%<>%dplyr::select(-"grading")
# 
# BDs%>%ggplot(aes(log(hazMax),group=as.factor(grading)))+geom_density(aes(colour=as.factor(grading),fill=as.factor(grading)),alpha=0.1)
# 
# tmp<-BDs%>%group_by(Event)%>%summarise(weighting=1/length(GNIc))
# BDs%<>%merge(tmp,by="Event")
# BDs%>%group_by(grading)%>%summarise(wmean=weighted.mean(hazMax,weighting),wcov=modi::weighted.var(hazMax,weighting))
# 
# # model<-e1071::svm(BDs[1:10000,-c(1,5)], BDs$Damage[1:10000], kernel = "polynomial", cost = 5, scale = FALSE)
# BDs$www<-BDs$weighting
# BDs$www[BDs$Damage==1]<-BDs$www[BDs$Damage==1]/sum(BDs$Damage==1)
# BDs$www[BDs$Damage==0]<-BDs$www[BDs$Damage==0]/sum(BDs$Damage==0)
# BDs$www<-BDs$www/sum(BDs$www)
# 
# BDs$Damage%<>%as.factor()
# levels(BDs$Damage)<-c("Unaffected","Damaged")
# BDs$hazMax%<>%unname();BDs$hazSD%<>%unname()
# 
# saveRDS(BDs,"./IIDIPUS_Results/SpatialPoints_ML-GLM/InputData_BD.RData")

BDs<-readRDS("./IIDIPUS_Results/SpatialPoints_ML-GLM/InputData_BD.RData")

BDcoords<-BDs%>%dplyr::select(Longitude,Latitude)
BDs%<>%dplyr::select(-Longitude,-Latitude)

tmp<-apply(BDs[,3:11],2,function(x) length(unique(x)))<30
BDs%<>%dplyr::select(-names(tmp[tmp]))
BDs$hazMax%<>%log()

covariates<-BDs%>%dplyr::select(colnames(BDs),-c(Event,grading,www,weighting,Damage))%>%colnames()

train_control <- caret::trainControl(method="repeatedcv", number=8, repeats=2,
                                     search = "random",classProbs=T,
                                     summaryFunction=twoClassSummary)

# parallelML<-function(algo) {
#   # Remove pesky columns
#   datar<-BDs%>%dplyr::select(-c("Event","grading","weighting","www"))
#   # Run the model!
#   modeler<-caret::train(Damage~., data = datar, method = algo, metric="ROC",
#                         tuneLength = 12, trControl = train_control,
#                         preProcess = c("center","scale"))
#   
#   return(cbind(modeler$results[-1],
#                t(as.data.frame((t(as.data.frame(varImp(modeler, scale=FALSE)$importance))[1,])))))
#   
# }

parallelML_balanced<-function(algo,splitties=NULL,ncores=4,retmod=F) {
  
  # How many damaged buildings are there?
  numun<-round(table(BDs$Damage)["Damaged"]*1.5)
  # By default, add all events that have less than 500 points
  permys<-BDs%>%filter(Event%in%names(table(BDs$Event)[table(BDs$Event)<500]))
  # Adding all damaged buildings, too
  permys<-BDs%>%filter(Damage=="Damaged" & 
                       !Event%in%names(table(BDs$Event)[table(BDs$Event)<500]))%>%
                rbind(permys)
  # Lets reduce the bias towards predicting well only the unaffected buildings
  indies<-which(BDs$Damage=="Unaffected" & !BDs$Event%in%names(table(BDs$Event)[table(BDs$Event)<500]))
  # How many times to split the dataset and model
  if(is.null(splitties)) splitties<-floor(length(indies)/numun)
  print(paste0("Working with the ",algo," model with data split into ",splitties," groups")) 
  # Now split the remaining into groups of indices
  indies <- createFolds(indies, k = splitties, list = T, returnTrain = FALSE)
  # Parallelise
  cl <- makePSOCKcluster(ncores)  # Create computing clusters
  registerDoParallel(cl)
  getDoParWorkers()
  # CV-split and model the damaged buildings
  out<-lapply(seq.int(1,length(indies),by=3),function(i){
    datar<-rbind(BDs[indies[[i]],],permys)
    datar%<>%dplyr::select(-c("Event","grading","weighting","www"))
    # Run the model!
    modeler<-caret::train(Damage~., data = datar, method = algo, metric="ROC",
                          tuneLength = 12, trControl = train_control,
                          preProcess = c("center","scale"))
    # Let me know!
    print(paste0(signif(100*i/splitties,2),"% done"))
    
    if(retmod) return(as.data.frame(vip::vi(modeler,
                                            feature_names=covariates))%>%dplyr::select(Variable,Importance))

    
    return(cbind(filter(modeler$results[-1],ROC==max(ROC)),
                 t(as.data.frame((t(as.data.frame(varImp(modeler, scale=FALSE)$importance))[1,])))))
  })
  
  # Remember to close the computing cluster
  stopCluster(cl)
  registerDoSEQ()
  # Save out, then get out!
  saveRDS(out,paste0("./IIDIPUS_Results/SpatialPoints_ML-GLM/NoSpace_ML_models/ML-",algo,"_2.RData"))
  
  return(out)
}

# tabmod<-getModelInfo()
# carmods<-unlist(sapply(tabmod,function(x) x$type%in%"Classification"))
# carmods<-data.frame(algorithm=names(carmods),classification=unname(carmods))
# carmods%<>%filter(classification)%>%pull(algorithm)
# # Check that we have all that we need to run each model
# checkerz<-unlist(lapply(carmods,function(stst) ifelse(is.null(tryCatch(checkInstall(getModelInfo(stst)$library),error=function(e) NA)),T,F)))
# carmods<-carmods[checkerz]; rm(checkerz)

minimods<-c("svmLinear","svmRadial","svmPoly","naive_bayes","rf","glmnet","ada","adaboost")

ncores<-60

# Run ALL THE MODELLLLLSSS
ML_BDs<-lapply(minimods,function(stst) tryCatch(parallelML_balanced(stst,ncores = ncores),error=function(e) NA))

# Parallelise
# cl <- makePSOCKcluster(60)  # Create computing clusters
# registerDoParallel(cl)
# getDoParWorkers()
# # Run ALL THE MODELLLLLSSS
# ML_BDs<-lapply(carmods,function(stst) tryCatch(parallelML_balanced(stst),error=function(e) NA))
# # Remember to close the computing cluster
# stopCluster(cl)
# registerDoSEQ()

# plot(varImp(modeler, scale=FALSE))
# tmp<-confusionMatrix(predict(modeler, BDs), BDs$Damage)

datar<-BDs%>%dplyr::select(-c("Event","grading","weighting","www"))
# datar<-datar[sample(1:nrow(datar),50000,F),]
ind<-colnames(datar)=="Damage"
# Let's run it!
set.seed(123)
yhat <- fastknnCV(data.matrix(datar[,!ind]),as.factor(datar[,ind]),
                  k = 1:15, method = "vote", folds = 8, eval.metric ="auc")

knn.out <- fastknn::fastknn(data.matrix(datar[,!ind]),
                            as.factor(datar[,ind]),
                            data.matrix(datar[,!ind]), k = yhat$best_k) 

# split <- sample(1:nrow(datar), size = 0.8*nrow(datar),replace = F)
# 
# featies<-knnDecision(data.matrix(datar[split,!ind]), as.factor(datar[split,ind]),
#                      data.matrix(datar[-split,!ind]), as.factor(datar[-split,ind]),
#                      k = yhat$best_k)

ROC<-max(yhat$cv_table$mean)
ROCSD<-sd(yhat$cv_table[which.max(yhat$cv_table$mean),1:8])

knn.out$pred<-factor(apply(knn.out$prob,1,which.max)-1,levels = c("Unaffected"=0,"Damaged"=1))

score<-knn.out$prob[,1]
class<-as.numeric(knn.out$class)-1

perfy<-ROCit::measureit(score = score, class = class,
          measure = c("SENS", "SPEC"))

ROCit::rocit(score = score, negref = "-",
      class = class, 
      method = "emp")

head(perfy)

namerz<-c("ROC","Sens","Spec","ROCSD","SensSD","SpecSD",colnames(datar)); namerz<-namerz[1:(length(namerz)-1)]
out<-as.data.frame(t(data.frame(namerz))); colnames(out)<-out[1,]
out[1,]<-NA
out$ROC<-ROC
out$ROCSD<-ROCSD

saveRDS(out,"./IIDIPUS_Results/SpatialPoints_ML-GLM/NoSpace_ML_models/ML-knn_2.RData")
# 
loccy<-"./IIDIPUS_Results/SpatialPoints_ML-GLM/NoSpace_ML_models/"
filez<-list.files(loccy);filez<-filez[grepl("_2",filez)];filez<-filez[!grepl("GPR",filez)]

out<-do.call(rbind,lapply(seq_along(filez),
                           function(i){
                             print(filez[i])
                             if(grepl("knn",filez[i])) return(cbind(data.frame(algo=filez[i]),dplyr::select(readRDS(paste0(loccy,filez[i])),namerz)))
                             cbind(data.frame(algo=filez[i]),
                                   t(colMeans(dplyr::select(do.call(rbind,readRDS(paste0(loccy,filez[i]))),namerz))))
                           }))

row.names(out)<-NULL

for(i in 8:ncol(out)) out[,i]%<>%as.numeric()
out[,-c(1:7)]<-100*out[,-c(1:7)]/rowSums(out[,-c(1:7)],na.rm = T)

# out%>%ggplot(aes(algo,ROC))+geom_point()


# tmp<-do.call(rbind,lapply(seq_along(out), function(i) out[[i]]$vippy))
# outer<-data.frame(variable=colnames(tmp),meany=colMeans(tmp),vary=apply(tmp,2,sd))%>%
#   arrange(desc(meany))







# outer<-do.call(rbind,parallelML_balanced("ada",ncores = 60,retmod = T))%>%
#   group_by(Variable)%>%
#   reframe(ImportanceSD=sd(Importance),
#           Importance=mean(Importance))%>%
#   dplyr::select(Variable,Importance,ImportanceSD)

outer<-out[1,-(1:7)]%>%reshape2::melt()

colnames(outer)<-c("Variable","Importance")

outer$Importance<-100*outer$Importance/sum(outer$Importance)

outer$Variable%<>%factor(levels=outer$Variable[rev(order(outer$Importance))])

namerz<-1:nrow(outer); names(namerz)<-outer$Variable

p<-outer%>%
  ggplot(aes(Variable,Importance))+
  geom_point(size=3)+
  # scale_shape_manual(values=15:19,breaks=allimps)+
  # scale_colour_manual(values = pal,limits = names(pal))+
  xlab("Model Covariate") + ylab("Feature Importance [%]")+
  ylim(c(0,25))+
  # scale_x_discrete(labels=namerz)+
  # labs(colour="Impact Type",shape="Impact Type")+
  theme(axis.text.x = element_text(angle = 90));p
ggsave("BD_Feat.eps",p,path="./Plots/IIDIPUS_Results/",width=10,height=4.,device = grDevices::cairo_ps)

outer%>%dplyr::select(meany)%>%t()%>%as.data.frame()%>%
  ggplot(aes())

p<-BDs%>%dplyr::select(Population,Vs30,EQFreq,hazMax,hazSD,Damage)%>%
  # mutate(hazMax=exp(hazMax))%>%
  reshape2::melt("Damage")%>%
  group_by(variable)%>%reframe(value=value/max(value),Damage=Damage)%>%
  ggplot(aes(variable,value))+geom_violin(aes(fill=Damage),scale = "width");p
ggsave("BD_Violin.eps",p,path="./Plots/IIDIPUS_Results/",width=12,height=4.,device = grDevices::cairo_ps)
# 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%% GAUSSIAN PROCESS REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

BDs%<>%dplyr::select(-c("grading","weighting","www"))

BDs%<>%cbind(BDcoords)

bbs<-do.call(rbind,lapply(unique(BDs$Event), function(ev){
  BDs%>%filter(Event==ev)%>%
    reframe(bbox1=min(Longitude,na.rm = T),
            bbox2=min(Latitude,na.rm = T),
            bbox3=max(Longitude,na.rm = T),
            bbox4=max(Latitude,na.rm = T))
}))

maxdist<-0.25 #max(0.5*sqrt((bbs[,3]-bbs[,1])*(bbs[,4]-bbs[,2])))/5

# How many damaged buildings are there?
numun<-round(table(BDs$Damage)["Damaged"]*1.5)
# By default, add all events that have less than 500 points
permys<-BDs%>%filter(Event%in%names(table(BDs$Event)[table(BDs$Event)<500]))
# Adding all damaged buildings, too
permys<-BDs%>%filter(Damage=="Damaged" & 
                       !Event%in%names(table(BDs$Event)[table(BDs$Event)<500]))%>%
  rbind(permys)
# Lets reduce the bias towards predicting well only the unaffected buildings
indies<-which(BDs$Damage=="Unaffected" & !BDs$Event%in%names(table(BDs$Event)[table(BDs$Event)<500]))
# How many times to split the dataset and model
splitties<-8
if(is.null(splitties)) splitties<-floor(length(indies)/numun)
print(paste0("Working with the GPR model with data split into ",splitties," groups")) 
# Now split the remaining into groups of indices
indies <- createFolds(indies, k = splitties, list = T, returnTrain = FALSE)
# Initial variogram
vario <- gstat::variogram(Damage~1, BDs, cutoff=maxdist)
fit <- gstat::fit.variogram(vario, model=gstat::vgm(0.05, "Gau", 0.18, nugget = min(vario$gamma)))
fit[1,2]<-0.004; fit[2,2:3]<-c(0.045,0.08); #plot(vario,fit)
# Formulate the equations
eqns<-c("Damage ~ 1")
for(n in 1:length(covariates)){
  eqn<-combn(covariates,n)
  if(n==1) {eqn%<>%as.character()
  }else eqn<-apply(eqn,2,function(x) pracma::strcat(x,collapse = "+"))
  eqns%<>%c(eqns,unname(sapply(eqn,function(eee) paste0("Damage ~ ",eee,""))))
}
eqns<-unique(eqns)
# Ensure all unique values in BD object
BDs[,covariates]<-BDs[,covariates]+rnorm(nrow(BDs),0,1e-5)
# Parallelise
ncores<-50
cl <- makePSOCKcluster(ncores)  # Create computing clusters
registerDoParallel(cl)
getDoParWorkers()
# CV-split and model the damaged buildings
out<-lapply(eqns, function(eee) {print(eee);do.call(rbind,lapply(1:length(indies),function(i){
  print(i)
  test<-BDs[indies[[i]],]
  train<-BDs[!(1:nrow(BDs))%in%indies[[i]],]
  # datar<-BDs
  test%<>%dplyr::select(-c("Event")); train%<>%dplyr::select(-c("Event"))
  test$Damage<-as.numeric(test$Damage)-1; train$Damage<-as.numeric(train$Damage)-1
  coordinates(test) <- ~ Longitude + Latitude; coordinates(train) <- ~ Longitude + Latitude
  
  # split <- sample(1:nrow(datar), size = 0.8*nrow(datar),replace = F)
  # train <- datar[split,]
  # test <- datar[-split,]
  
  # Run the model!
  # vario <- gstat::variogram(Damage~Population + Vs30 + EQFreq + hazMax + hazSD, train, cutoff=maxdist)
  
  parts <- split(x = 1:nrow(test), f = 1:ncores)
  
  clusterExport(cl = cl, varlist = c("test", "train", "fit", "parts", "eee"), envir = environment())
  clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
  
  out<-do.call(rbind,parLapply(cl = cl, X = 1:ncores, 
                         fun = function(x) cbind(as.data.frame(krige(formula = as.formula(eee),nmax=300,
                                                    locations = train, newdata = test[parts[[x]],],
                                                    model = fit),col.names=c("prediction"))$var1.pred,as.data.frame(test[parts[[x]],]))))
  
  return(ROCit::rocit(score = out[,1], class = out$Damage)$AUC)

}))})
# Remember to close the computing cluster
stopCluster(cl)
registerDoSEQ()

outer<-do.call(rbind,lapply(1:length(out), function(i) data.frame(ROC=mean(out[[i]]),ROCSD=sd(out[[i]]),equation=eqns[i])))

outer$nn<-8

outer%<>%mutate(minSD=ROCSD[which.min(ROC)]^2/unique(nn),
                                       thisSD=ROCSD^2/unique(nn))%>%
  mutate(df=(nn-1)*(thisSD+minSD)^2/(thisSD^2+minSD^2))

outer%<>%mutate(tval=sqrt(nn)*(ROC-min(ROC))/
                                         sqrt(ROCSD^2+ROCSD[which.min(ROC)]^2))
outer%<>%mutate(pval=dt(tval,df))
# 



# Save out, then get out!
saveRDS(outer,paste0("./IIDIPUS_Results/SpatialPoints_ML-GLM/NoSpace_ML_models/ML-GPR_2.RData"))

GPR<-outer%>%arrange(desc(ROC))%>%slice(1)%>%dplyr::select(ROC,ROCSD)%>%
  cbind(data.frame(Sens=NA,
                   Spec=NA,
                   SensSD=NA,
                   SpecSD=NA))
# 

out%<>%rbind(cbind(dplyr::select(GPR,namerz),data.frame(algo="ML-GPR_2.RData")))

out$algo<-str_split(str_split(out$algo,"_2.RData",simplify = T)[,1],"-",simplify = T)[,2]
rownames(out)<-NULL
out[,1:2]

ordz<-out%>%arrange(desc(ROC))%>%pull(algo)

out$algo%<>%factor(levels=ordz)

p<-out%>%ggplot(aes(algo,ROC)) +
  geom_point(aes(colour=algo),size=3)+ylim(c(0.7,1))+
  geom_errorbar(aes(ymin=ROC-ROCSD, ymax=ROC+ROCSD,colour=algo),width=.4)+
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Model")+ylab("Area Under ROC Curve (AUC)")+
  labs(colour="Model");p

ggsave("BD_Perf.eps",p,path="./Plots/IIDIPUS_Results/",width=8,height=5.,device = grDevices::cairo_ps)












# vario<-geoR::variog(as.array(BDs[,c("Longitude","Latitude")]),option = "cloud",max.dist = maxdist)
# 
# tmpK<-krige.cv(Damage~Population + Vs30 + EQFreq + hazMax + hazSD, 
#                test, fit, maxdist=maxdist, nfold=8)
# 
# 
# G <- auto_basis(manifold = plane(), # 2D plane
#                 data = train, # meuse data
#                 nres = 2, # number of resolutions
#                 type = "Gaussian", # type of basis function
#                 regular = 1)
# 
# BAUs <- auto_BAUs(manifold = plane(),
#                   type = "grid",
#                   data = train,
#                   nonconvex_hull = FALSE)
# 
# BAUs@data%<>%cbind(dplyr::select(train@data,Population , Vs30 , EQFreq , hazMax , hazSD))
# 
# outFRK<-FRK(f = Damage~Population + Vs30 + EQFreq + hazMax + hazSD, # Formula to FRK
#          list(train), # All datasets are supplied in list
#          nres = 2, # Low-rank model to reduce run-time
#          response = "gaussian", # data model
#          link = "identity", # link function
#          normalise_wts = T,
#          nonconvex_hull = FALSE)
# 
# predy <- predict(outFRK, newdata=test)
# 
# ROCit::rocit(score = predy@data$mu, class = as.numeric(predy@data$Damage))
# 
# outFRK_full<-FRK(f = Damage~Population + Vs30 + EQFreq + hazMax + hazSD, # Formula to FRK
#             list(train), # All datasets are supplied in list
#             nres = 2, # Low-rank model to reduce run-time
#             response = "gaussian", # data model
#             link = "identity", # link function
#             normalise_wts = T,
#             nonconvex_hull = FALSE)
# 
# predy_full <- predict(outFRK_full, newdata=test)
# 
# ROCit::rocit(score = predy_full@data$mu, class = as.numeric(predy_full@data$Damage))


#
# 
# 





# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# source('RCode/GetODDPackages.R')
# library(tidyverse)
# library(reticulate)
# library(tensorflow)
# library(keras)
# use_implementation("tensorflow")
# 
# # install.packages("tfprobability")
# library(tfprobability)
# 
# BDs<-readRDS("./IIDIPUS_Input/BD_nonparametric.Rdata")
# tmp<-dplyr::select(BDs,c(Grading,PopDens,HazardIntensity,GDP,RWI))  
# tmp$PopDens<-log(tmp$PopDens+1)
# tmp$GDP<-log(tmp$GDP)
# tmp%<>%filter(as.character(Grading)%in%c("destroyed","moderate","severe","notaffected"))
# tmp$Damage<-as.integer(tmp$Grading=="notaffected")
# tmp%<>%filter(!(is.na(tmp$Grading)|is.na(tmp$Damage)|is.na(tmp$GDP)|is.na(tmp$PopDens)))
# 
# model<-glm(formula = "Damage~HazardIntensity*PopDens*GDP+0",tmp,family = "binomial")
# # model<-glm("Damage~1",tmp,family = "binomial")
# # model_stepAIC<-MASS::stepAIC(model,direction = "forward")
# summary(model)
# data.frame(pred=model$fitted.values,Damage=tmp$Damage)%>%ggplot(aes(pred,group=Damage))+geom_density(aes(fill=tmp$Damage))
# 
# n<-15
# 
# damage<-sampleBDdamage(tmp$Grading,n = n)
# damage<-logit(damage)
# 
# donnees<-data.frame()
# for (j in 1:n){
#   donnees<-tmp%>%mutate(damage=damage[j,])%>%rbind(donnees)
# }
# donnees%<>%dplyr::select(-c(Grading))
# 
# bt <- import("builtins")
# RBFKernelFn <- reticulate::PyClass(
#   "KernelFn",
#   inherit = tensorflow::tf$keras$layers$Layer,
#   list(
#     `__init__` = function(self, ...) {
#       kwargs <- list(...)
#       super()$`__init__`(kwargs)
#       dtype <- kwargs[["dtype"]]
#       self$`_amplitude` = self$add_variable(initializer = initializer_zeros(),
#                                             dtype = dtype,
#                                             name = 'amplitude')
#       self$`_length_scale` = self$add_variable(initializer = initializer_zeros(),
#                                                dtype = dtype,
#                                                name = 'length_scale')
#       NULL
#     },
#     
#     call = function(self, x, ...) {
#       x
#     },
#     
#     kernel = bt$property(
#       reticulate::py_func(
#         function(self)
#           tfp$math$psd_kernels$ExponentiatedQuadratic(
#             amplitude = tf$nn$softplus(array(0.1) * self$`_amplitude`),
#             length_scale = tf$nn$softplus(array(2) * self$`_length_scale`)
#           )
#       )
#     )
#   )
# )
# 
# num_inducing_points <- 50
# k_set_floatx("float64")
# 
# maxs <- apply(donnees, 2, max) 
# mins <- apply(donnees, 2, min)
# scaled <- as.data.frame(scale(donnees, center = mins, scale = maxs - mins))
# # n <- names(scaled[,])
# # f <- as.formula(paste("damage ~", paste(n[!n %in% "damage"], collapse = " + ")))
# split <- sample(1:nrow(scaled), size = 0.8*nrow(scaled),replace = F)
# train <- scaled[split,]
# test <- scaled[-split,]
# 
# sample_dist <- tfd_uniform(low = 1, high = nrow(train) + 1)
# sample_ids <- sample_dist %>%
#   tfd_sample(num_inducing_points) %>%
#   tf$cast(tf$int32) %>%
#   as.numeric()
# sampled_points <- train[sample_ids, 1:4]
# 
# model <- keras_model_sequential() %>%
#   layer_dense(units = 2,
#               input_shape = 2,
#               use_bias = FALSE) %>%
#   layer_variational_gaussian_process(
#     num_inducing_points = num_inducing_points,
#     kernel_provider = RBFKernelFn(),
#     event_shape = 1,
#     inducing_index_points_initializer = initializer_constant(as.matrix(sampled_points)),
#     unconstrained_observation_noise_variance_initializer =
#       initializer_constant(array(0.1))
#   )
# 
# 
# 
# 
# 
# 
# 
# 
