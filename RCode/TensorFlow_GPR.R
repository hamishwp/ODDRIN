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








