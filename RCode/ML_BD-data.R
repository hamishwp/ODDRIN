# install.packages(c("bit64", "snow", "SparseM", "Matrix"))
# devtools::install_github("Danko-lab/Rgtsvm/Rgtsvm", args="--configure-args='--with-cuda-home=YOUR_CUDA_PATH --with-boost-home=YOU_BOOST_PATH'")
install.packages("modi")
# source('RCode/GetODDPackages.R')

folder<-"./IIDIPUS_Input/IIDIPUS_Input_NMAR/BDobjects/"
filez<-list.files(folder)
haz<-"EQ"
if(haz=="EQ") funcyfun<-exp else funcyfun<-returnX

BDs<-data.frame()
for(fff in filez){
  
  # if(fff%in%c("EQ20180225PNG_68","EQ20190925IDN_95","EQ20201030TUR_88")) next
  
  BDy<-readRDS(paste0(folder,fff))
  if(length(BDy@data[,-grep(names(BDy@data),pattern = "haz",value = F)])<9) next
  
  if(length(grep(names(BDy@data),pattern = "hazMean",value = T))==1) {
    hazMax<-BDy$hazMean1
    hazSD<-BDy$hazSD1
  } else {
    hazMax<-apply(BDy@data[,grep(names(BDy@data),pattern = "hazMean",value = T)],1,max,na.rm=T)
    hazSD<-apply(BDy@data[,grep(names(BDy@data),pattern = "hazSD",value = T)],1,median,na.rm=T)
  }
  BDy@data%<>%dplyr::select(-c(grep(names(BDy@data),pattern = "hazMean",value = T),
                               grep(names(BDy@data),pattern = "hazSD",value = T),
                               grep(names(BDy@data),pattern = "itude",value = T),
                               grep(names(BDy@data),pattern = "nBuildings",value = T)))%>%
    mutate(Event=fff,date=BDy@hazdates[1],hazMax=funcyfun(hazMax-BDy@I0),hazSD=hazSD)
  
  BDs%<>%rbind(BDy@data)
  
}
rm(BDy)

BDs$time<-as.numeric(BDs$date-min(BDs$date))

BDs%<>%dplyr::select(-c("date","Confidence","ISO3C"))  
BDs$Population<-log(BDs$Population+1)
BDs$GNIc<-log(BDs$GNIc)
BDs%<>%filter(!as.character(grading)%in%c("possible") & apply(BDs,1,function(x) !any(is.na(x))))
BDs$Damage<-abs(1-as.integer(BDs$grading=="notaffected"))
# BDs%<>%dplyr::select(-"grading")

BDs%>%ggplot(aes(log(hazMax),group=as.factor(grading)))+geom_density(aes(colour=as.factor(grading),fill=as.factor(grading)),alpha=0.1)

tmp<-BDs%>%group_by(Event)%>%summarise(weighting=1/length(GNIc))
BDs%<>%merge(tmp,by="Event")
BDs%>%group_by(grading)%>%summarise(wmean=weighted.mean(hazMax,weighting),wcov=modi::weighted.var(hazMax,weighting))

library(caret)
# model<-e1071::svm(BDs[1:10000,-c(1,5)], BDs$Damage[1:10000], kernel = "polynomial", cost = 5, scale = FALSE)
BDs$www<-BDs$weighting
BDs$www[BDs$Damage==1]<-BDs$www[BDs$Damage==1]/sum(BDs$Damage==1)
BDs$www[BDs$Damage==0]<-BDs$www[BDs$Damage==0]/sum(BDs$Damage==0)
BDs$www<-BDs$www/sum(BDs$www)

BDs$Damage%<>%as.factor()
levels(BDs$Damage)<-c("Unaffected","Damaged")
BDs$hazMax%<>%unname();BDs$hazSD%<>%unname()

train_control <- trainControl(method="repeatedcv", number=10, repeats=5,classProbs=T)#,summaryFunction=twoClassSummary)

parallelML<-function(algo) {
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
  # Now split the remaining into groups of indices
  indies <- createFolds(indies, k = floor(length(indies)/numun), list = T, returnTrain = FALSE)
  # CV-split and model the damaged buildings
  as.data.frame(t(colMeans(do.call(rbind, lapply(1:length(indies),function(i){
    datar<-rbind(BDs[indies[[i]],],permys)
    datar%<>%dplyr::select(-c("Event","grading","weighting","www"))
    modeler<-caret::train(Damage~., data = datar, method = algo, metric="ROC", 
                          trControl = train_control,  preProcess = c("center","scale"))
    
    
    # What is it you want to extract? ROC? Kappa?
    
  })))))
  
  
}

cl <- makePSOCKcluster(12)  # Create 8 clusters
registerDoParallel(cl)
getDoParWorkers()

svm1<-parallelML("svmLinear")

stopCluster(cl)
registerDoSEQ()

plot(varImp(svm1, scale=FALSE))
confusionMatrix(predict(svm1, BDs), BDs$Damage)

saveRDS(svm1,"./IIDIPUS_Results/ML_Models/BD_svm_initial.RData")

tabmod<-getModelInfo()
carmods<-unlist(sapply(tabmod,function(x) x$type%in%"Classification"))
carmods<-data.frame(algorithm=names(carmods),classification=unname(carmods))
carmods%<>%filter(classification)
#








# BDs<-BDs[1:10000,]

set.seed(2016)

# USING DEFAULT RBF-KERNEL
gttune<-BDs%>%dplyr::select(-c("Event","grading","weighting"))%>%
  Rgtsvm::tune.svm(BDs$Damage,
                   gamma=c(2.25,2.5,2.75),
                   cost=c(30,35,40,45),
                   tunecontrol=tune.control(sampling="cross",cross=8,rough.cross = 3),
                   scale=F)
# mdl <- Rgtsvm::svm(BDs[,1:3],BDs$Damage, type = "C-classification", kernel = "radial", degree=3, cost = 5, gamma = 1, probability = TRUE)

saveRDS(gttune,"./IIDIPUS_Results/SVM/SVM_tune_gamma2_2.25_2.5_2.75_cost_30_35_40_45.Rdata")

isplit<-sample(1:nrow(BDs),size = round(0.8*nrow(BDs)))
mdl <- Rgtsvm::svm(BDs[isplit,1:3],BDs$Damage[isplit], kernel = "radial", 
                   cost = gttune$best.parameters$cost, gamma = gttune$best.parameters$gamma, 
                   tunecontrol=tune.control(sampling="cross",cross=8,rough.cross = 3),
                   probability = TRUE)

saveRDS(mdl,"./IIDIPUS_Results/SVM/SVM_tuned_model.Rdata")

pred <- predict(mdl, BDs[-isplit,1:3], decision.values = TRUE, probability = TRUE)
table(pred==BDs$Damage[-isplit])
pred<-as.integer(as.character(pred))





BDs$fold <- caret::createFolds(1:nrow(BDs), k = 8, list = FALSE)
### PARAMETER LIST ###
cost <- c(5,10,20)
gamma <- c(1,2,3)
parms <- expand.grid(cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:max(BDs$fold), .combine = rbind, .inorder = FALSE) %dopar% {
    deve <- BDs[BDs$fold != j, ]
    test <- BDs[BDs$fold == j, ]
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
  out <- foreach(j = 1:max(BDs$fold), .combine = rbind, .inorder = FALSE) %dopar% {
    deve <- BDs[BDs$fold != j, ]
    test <- BDs[BDs$fold == j, ]
    mdl <- Rgtsvm::svm(Damage~., data = deve, type = "C-classification", kernel = "polynomial", degree=3, cost = c, gamma = g, probability = TRUE)
    pred <- predict(mdl, test, decision.values = TRUE, probability = TRUE)
    data.frame(y = test$Damage, prob = attributes(pred)$probabilities[, 2])
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- pROC::roc(as.factor(out$y), out$prob) 
  data.frame(parms[i, ], roc = roc$auc[1])
}































































source('RCode/GetODDPackages.R')
library(tidyverse)
library(reticulate)
library(tensorflow)
library(keras)
use_implementation("tensorflow")

# install.packages("tfprobability")
library(tfprobability)

BDs<-readRDS("./IIDIPUS_Input/BD_nonparametric.Rdata")
tmp<-dplyr::select(BDs,c(Grading,PopDens,HazardIntensity,GDP,RWI))  
tmp$PopDens<-log(tmp$PopDens+1)
tmp$GDP<-log(tmp$GDP)
tmp%<>%filter(as.character(Grading)%in%c("destroyed","moderate","severe","notaffected"))
tmp$Damage<-as.integer(tmp$Grading=="notaffected")
tmp%<>%filter(!(is.na(tmp$Grading)|is.na(tmp$Damage)|is.na(tmp$GDP)|is.na(tmp$PopDens)))

model<-glm(formula = "Damage~HazardIntensity*PopDens*GDP+0",tmp,family = "binomial")
# model<-glm("Damage~1",tmp,family = "binomial")
# model_stepAIC<-MASS::stepAIC(model,direction = "forward")
summary(model)
data.frame(pred=model$fitted.values,Damage=tmp$Damage)%>%ggplot(aes(pred,group=Damage))+geom_density(aes(fill=tmp$Damage))

n<-15

damage<-sampleBDdamage(tmp$Grading,n = n)
damage<-logit(damage)

donnees<-data.frame()
for (j in 1:n){
  donnees<-tmp%>%mutate(damage=damage[j,])%>%rbind(donnees)
}
donnees%<>%dplyr::select(-c(Grading))

bt <- import("builtins")
RBFKernelFn <- reticulate::PyClass(
  "KernelFn",
  inherit = tensorflow::tf$keras$layers$Layer,
  list(
    `__init__` = function(self, ...) {
      kwargs <- list(...)
      super()$`__init__`(kwargs)
      dtype <- kwargs[["dtype"]]
      self$`_amplitude` = self$add_variable(initializer = initializer_zeros(),
                                            dtype = dtype,
                                            name = 'amplitude')
      self$`_length_scale` = self$add_variable(initializer = initializer_zeros(),
                                               dtype = dtype,
                                               name = 'length_scale')
      NULL
    },
    
    call = function(self, x, ...) {
      x
    },
    
    kernel = bt$property(
      reticulate::py_func(
        function(self)
          tfp$math$psd_kernels$ExponentiatedQuadratic(
            amplitude = tf$nn$softplus(array(0.1) * self$`_amplitude`),
            length_scale = tf$nn$softplus(array(2) * self$`_length_scale`)
          )
      )
    )
  )
)

num_inducing_points <- 50
k_set_floatx("float64")

maxs <- apply(donnees, 2, max) 
mins <- apply(donnees, 2, min)
scaled <- as.data.frame(scale(donnees, center = mins, scale = maxs - mins))
# n <- names(scaled[,])
# f <- as.formula(paste("damage ~", paste(n[!n %in% "damage"], collapse = " + ")))
split <- sample(1:nrow(scaled), size = 0.8*nrow(scaled),replace = F)
train <- scaled[split,]
test <- scaled[-split,]

sample_dist <- tfd_uniform(low = 1, high = nrow(train) + 1)
sample_ids <- sample_dist %>%
  tfd_sample(num_inducing_points) %>%
  tf$cast(tf$int32) %>%
  as.numeric()
sampled_points <- train[sample_ids, 1:4]

model <- keras_model_sequential() %>%
  layer_dense(units = 2,
              input_shape = 2,
              use_bias = FALSE) %>%
  layer_variational_gaussian_process(
    num_inducing_points = num_inducing_points,
    kernel_provider = RBFKernelFn(),
    event_shape = 1,
    inducing_index_points_initializer = initializer_constant(as.matrix(sampled_points)),
    unconstrained_observation_noise_variance_initializer =
      initializer_constant(array(0.1))
  )








