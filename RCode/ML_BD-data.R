# install.packages("caret",repos = "http://cran.r-project.org", 
#                  dependencies = c("Depends", "Imports","Suggests"))
# devtools::install_github("souravc83/fastAdaboost")
# devtools::install_github("davpinto/fastknn")

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

tmp<-apply(BDs[,3:11],2,function(x) length(unique(x)))<30
BDs%<>%dplyr::select(-names(tmp[tmp]))
BDs$hazMax%<>%log()

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

parallelML_balanced<-function(algo,splitties=NULL,ncores=4) {
  
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

minimods<-c("svmLinear","svmRadial","svmPoly","naive_bayes","rf","glmnet","AdaBoost")

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

namerz<-c("ROC","Sens","Spec","ROCSD","SensSD","SpecSD",colnames(datar)); namerz<-namerz[1:(length(namerz)-1)]
out<-as.data.frame(t(data.frame(namerz))); colnames(out)<-out[1,]
out[1,]<-NA
out$ROC<-ROC
out$ROCSD<-ROCSD

saveRDS(out,"./IIDIPUS_Results/SpatialPoints_ML-GLM/NoSpace_ML_models/ML-knn_2.RData")
# 
# loccy<-"./IIDIPUS_Results/SpatialPoints_ML-GLM/NoSpace_ML_models/"
# filez<-list.files(loccy)
# 
# out<-readRDS(paste0(loccy,filez[i]))
# 
# out<-do.call(rbind, lapply(seq_along(filez),
#                            function(i){ cbind(data.frame(algo=filez[i]),
#                                               dplyr::select(readRDS(paste0(loccy,filez[i])),namerz))}))
# 
# tmp<-do.call(rbind,lapply(seq_along(out), function(i) out[[i]]$vippy))
# outer<-data.frame(variable=colnames(tmp),meany=colMeans(tmp),vary=apply(tmp,2,sd))%>%
#   arrange(desc(meany))
# outer$meany<-100*outer$meany/sum(outer$meany)
# 
# outer$variable%<>%factor(levels=outer$variable[rev(order(outer$meany))])
# 
# namerz<-1:nrow(outer); names(namerz)<-outer$variable
# 
# p<-outer%>%
#   ggplot(aes(variable,meany))+
#   geom_point(size=3)+
#   geom_line(alpha=0.5,linetype="dotdash")+
#   # scale_shape_manual(values=15:19,breaks=allimps)+
#   # scale_colour_manual(values = pal,limits = names(pal))+
#   xlab("Model Covariate") + ylab("Feature Importance [%]")+
#   # scale_x_discrete(labels=namerz)+
#   # labs(colour="Impact Type",shape="Impact Type")+
#   theme(axis.text.x = element_text(angle = 90));p
# ggsave("BD_Feat.eps",p,path="./Plots/IIDIPUS_Results/",width=10,height=4.,device = grDevices::cairo_ps)  
# 
# outer%>%dplyr::select(meany)%>%t()%>%as.data.frame()%>%
#   ggplot(aes())
# 
# out$algo<-str_split(str_split(out$algo,".RData",simplify = T)[,1],"-",simplify = T)[,2]
# rownames(out)<-NULL
# out[,1:2]
# 
# out%>%ggplot(aes(algo,ROC)) + 
#   geom_point(aes(colour=algo),size=3)+
#   geom_errorbar(aes(ymin=ROC-ROCSD, ymax=ROC+ROCSD,colour=algo),width=.4)
# # 
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
