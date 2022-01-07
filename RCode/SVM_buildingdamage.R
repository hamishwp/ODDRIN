source('RCode/GetODDPackages.R')

folder<-"./IIDIPUS_Input/BDobjects/"
filez<-list.files(folder)
haz<-"EQ"
if(haz=="EQ") funcyfun<-exp else funcyfun<-returnX

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
rm(BDy)

tmp<-dplyr::select(BDs,c(Grading,PopDens,HazardIntensity,GDP))  
tmp$PopDens<-log(tmp$PopDens+1)
tmp$GDP<-log(tmp$GDP)
tmp%<>%filter(as.character(Grading)%in%c("destroyed","moderate","severe","notaffected"))
tmp$Damage<-abs(1-as.integer(tmp$Grading=="notaffected"))
tmp%<>%filter(!(is.na(tmp$Grading)|is.na(tmp$Damage)|is.na(tmp$GDP)|is.na(tmp$PopDens)))
tmp%<>%dplyr::select(-"Grading")

BDs%>%ggplot(aes(HazardIntensity,RWI,group=Damage))+geom_point(aes(colour=Damage))


library(caret)
# model<-e1071::svm(tmp[1:10000,-c(1,5)], tmp$Damage[1:10000], kernel = "polynomial", cost = 5, scale = FALSE)
train_control <- trainControl(method="repeatedcv", number=5, repeats=5)
svm1 <- caret::train(Damage ~., data = tmp[1:5000,-1], method = "svmLinear", trControl = train_control)#,  preProcess = c("center","scale"),)

pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 8)

# tmp<-tmp[1:10000,]

set.seed(2016)

# USING DEFAULT RBF-KERNEL
gttune<-Rgtsvm::tune.svm(tmp[,1:3],tmp$Damage,
                         gamma=c(2.25,2.5,2.75),
                         cost=c(30,35,40,45),
                         tunecontrol=tune.control(sampling="cross",cross=8,rough.cross = 3),
                         scale=F)
# mdl <- Rgtsvm::svm(tmp[,1:3],tmp$Damage, type = "C-classification", kernel = "radial", degree=3, cost = 5, gamma = 1, probability = TRUE)

saveRDS(gttune,"./IIDIPUS_Results/SVM/SVM_tune_gamma2_2.25_2.5_2.75_cost_30_35_40_45.Rdata")

isplit<-sample(1:nrow(tmp),size = round(0.8*nrow(tmp)))
mdl <- Rgtsvm::svm(tmp[isplit,1:3],tmp$Damage[isplit], kernel = "radial", 
                   cost = gttune$best.parameters$cost, gamma = gttune$best.parameters$gamma, 
                   tunecontrol=tune.control(sampling="cross",cross=8,rough.cross = 3),
                   probability = TRUE)

saveRDS(mdl,"./IIDIPUS_Results/SVM/SVM_tuned_model.Rdata")

pred <- predict(mdl, tmp[-isplit,1:3], decision.values = TRUE, probability = TRUE)
table(pred==tmp$Damage[-isplit])
pred<-as.integer(as.character(pred))





tmp$fold <- caret::createFolds(1:nrow(tmp), k = 8, list = FALSE)
### PARAMETER LIST ###
cost <- c(5,10,20)
gamma <- c(1,2,3)
parms <- expand.grid(cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:max(tmp$fold), .combine = rbind, .inorder = FALSE) %dopar% {
    deve <- tmp[tmp$fold != j, ]
    test <- tmp[tmp$fold == j, ]
    mdl <- Rgtsvm::svm(Damage~., data = deve, type = "C-classification", kernel = "radial", cost = c, gamma = g, probability = TRUE)
    pred <- predict(mdl, test, decision.values = TRUE, probability = TRUE)
    data.frame(y = test$Damage, prob = attributes(pred)$probabilities[, 2])
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- pROC::roc(as.factor(out$y), out$prob) 
  data.frame(parms[i, ], roc = roc$auc[1])
}

mdl_radial<-mdl

result2 <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:max(tmp$fold), .combine = rbind, .inorder = FALSE) %dopar% {
    deve <- tmp[tmp$fold != j, ]
    test <- tmp[tmp$fold == j, ]
    mdl <- Rgtsvm::svm(Damage~., data = deve, type = "C-classification", kernel = "polynomial", degree=3, cost = c, gamma = g, probability = TRUE)
    pred <- predict(mdl, test, decision.values = TRUE, probability = TRUE)
    data.frame(y = test$Damage, prob = attributes(pred)$probabilities[, 2])
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- pROC::roc(as.factor(out$y), out$prob) 
  data.frame(parms[i, ], roc = roc$auc[1])
}
