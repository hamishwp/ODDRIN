dir<-directory<- "/home/manderso/Documents/GitHub/ODDRIN/" #"./" #"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/"
# Set the working directory from your environment variables
setwd(directory)

# Extract Environment Variables
source('RCode/GetEnv.R')
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
# Sourcing the data:
source('RCode/GetData.R')
# Extract model functions and priors
source('RCode/Model.R')
# Extract the model parameterisation algorithm, default = Adaptive MCMC
source('RCode/Method.R')
#Extract the functions for generating the simulations
source('RCode/Simulate.R')

AlgoParams$Np = 2
AlgoParams$m_CRPS = 2
AlgoParams$smc_Npart = 10
AlgoParams$cores = 4
AlgoParams$n_nodes = 1
out <- delmoral_parallel(AlgoParams, Model, unfinished=F, oldtag='')

#Parameterise the model and simulate the data:

Omega <- list(Lambda1 = list(nu=9.6, kappa=1.9049508),
              Lambda2 = list(nu=9.8626452, kappa=1.9049508),
              Lambda3 = list(nu=9.3, kappa=0.9),
              Lambda4 = list(nu=9.9, kappa=1.9),
              theta= list(theta1=0.6),
              eps = list(local=2.2292053, hazard_mort=0.8383464, hazard_disp=0.9, hazard_bd=0.9, hazard_cor=0.55),
              vuln_coeff = list(PDens=0, SHDI=-0.5, GNIc=-0.1, Vs30=0.1, EQFreq=-0.1, FirstHaz=0.05, Night=0.05, FirstHaz.Night=0.1),
              check = list(check=0.5))

Model$HighLevelPriors(Omega %>% addTransfParams(),Model)

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects_SubNat/Train/EQ20191118PHL_130')

xx <- DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, Method = AlgoParams, output='SampledAgg')

for(var in c('PDens', 'EQFreq', 'GNIc', 'Vs30', 'SHDI', 'hazMean1')){
  print(range(ODDyOld[[var]], na.rm=T))
  print(range(ODDyAgg[[var]], na.rm=T))
}

 
 Omega %<>% addTransfParams()
 Intensity <- seq(4.5, 10, 0.01)
 Dfun<-function(I_ij, Omega) h_0(I = I_ij,I0 = 4.3,Omega=Omega)
 Damage <- Dfun(Intensity, Omega)
 D_MortDisp <-  D_MortDisp_calc(Damage, Omega)
 D_MortDisp <-  D_MortDisp_calc(Damage, Omega)
 D_BD <- D_Dam_calc(Damage, Omega)
 # D_BD <-  plnorm(Damage, Omega$Lambda3$nu, Omega$Lambda3$omega)
 #plot(Intensity, D_MortDisp[1,], col='red', type='l', ylab='Proportion', lwd=3);
 plot(Intensity, D_MortDisp[2,], col='blue', type='l', ylab='Proportion', lwd=3);
 #lines(Intensity, log(D_BD), col='green', type='l', ylab='Proportion', lwd=3);
 for (i in 1:NROW(params_store)){
   Omega <- params_store[i,] %>% relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model)
   Omega %<>% addTransfParams()
   Intensity <- seq(4.5, 10, 0.01)
   Dfun<-function(I_ij, Omega) h_0(I = I_ij,I0 = 4.3,Omega=Omega)
   Damage <- Dfun(Intensity, Omega)
   D_MortDisp <-  D_MortDisp_calc(Damage, Omega)
   D_MortDisp <-  D_MortDisp_calc(Damage, Omega)
   D_BD <- D_Dam_calc(Damage, Omega)
   # D_BD <-  plnorm(Damage, Omega$Lambda3$nu, Omega$Lambda3$omega)
   #lines(Intensity, D_MortDisp[1,], col='red', type='l', ylab='Proportion', lwd=3);
   lines(Intensity, D_MortDisp[2,], col='blue', type='l', ylab='Proportion', lwd=3);
   #lines(Intensity, log(D_BD), col='green', type='l', ylab='Proportion', lwd=3);
}
#  
 
 
 AlgoParams$cores <- 1
 AlgoParams$NestedCores <- 6
 AlgoParams$Np <- 2
 AlgoParams$m_CRPS <- 5
 start_time <- Sys.time()
 impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams)
 finish_time <- Sys.time()
 finish_time-start_time
 
 dist_samp <- CalcDist(impact_sample, AlgoParams)
 
# 
# AlgoParams$cores <- 1
# AlgoParams$NestedCores <- 6
# AlgoParams$Np <- 2
# AlgoParams$m_CRPS <- 15
# start_time <- Sys.time()
# impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams)
# finish_time <- Sys.time()
# finish_time-start_time
# 
# d_prop <- CalcDist(impact_sample, AlgoParams)
# 
# AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-08-24_074550')
# AlgoResults$d_full = array(NA,  dim=c(AlgoParams$smc_Npart, NCOL(d_prop), AlgoParams$Np * AlgoParams$m_CRPS, AlgoParams$smc_steps))
# AlgoResults$sampled_full = array(NA,  dim=c(AlgoParams$smc_Npart, NROW(impact_sample$poly[[1]]), AlgoParams$Np * AlgoParams$m_CRPS, AlgoParams$smc_steps))
# 
# start_time <- Sys.time()
# SamplePointImpact(dir,Model, Omega %>% addTransfParams(),AlgoParams)
# finish_time <- Sys.time()
# finish_time - start_time
# 
# impact_filt <- 'mortality'
# plot_df <- data.frame(
#   obs = impact_sample$poly[[1]] %>% filter(impact==impact_filt) %>% pull(observed),
#   sample1 = impact_sample$poly[[1]] %>% filter(impact==impact_filt) %>% pull(sampled),
#   sample2 = impact_sample$poly[[2]] %>% filter(impact==impact_filt) %>% pull(sampled))
# 
# ggplot(plot_df) + geom_errorbar(aes(x=log(obs+0.1), ymin=log(sample1+0.1), ymax=log(sample2+0.1), width=0.2)) + 
#   geom_abline(slope=1, intercept=0) + ggtitle(impact_filt)
# 
# plot(log(impact_sample_filt$observed+0.1), log(impact_sample_filt$sampled+0.1))
# 
# sum(elapsed_time)
# 
# 
# ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_June24_AggFactor5/ODDobjects/EQ20131015PHL_16')
# DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)
# 
# ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_June20/ODDobjects/EQ20151026AFG_170')
# DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)
# 
# SamplePointImpact(dir,Model,proposed = Omega %>% addTransfParams(),AlgoParams,expLL=T)
#   
# # start_time <- Sys.time()
# impact_sample <- SampleImpact(dir = dir,Model = Model,
#                           proposed = Omega %>% addTransfParams(),
#                           AlgoParams = AlgoParams)
# 
# 
# 
# 
# d_i <- CalcDist(impact_sample, AlgoParams)
# # end_time <- Sys.time()
# # print(paste('Time:', end_time-start_time))
# 
# Model$center <- simulateDataSet(5, Omega %>% addTransfParams(), Model=Model, dir = dir, outliers = FALSE)
# 
# #After generating the simulated data, need to move 'from 'centerings' 
# #and 'ODDobjects' from 'IIDIPUS_SimInput' to 'IIDIPUS_Input'
# 
# main_simulated <- function(){
#   #ODDpaths<-ExtractData(Model$haz,dir,extractedData)
#   # Extract initial parameterisation estimate
#   
#   # Parameterise... Here we go!
#   
#   iVals<-GetInitVals(dir,Model,AlgoParams)
# 
#   output <- delmoral_parallel(dir=dir,
#                       Model=Model,
#                       iVals=iVals,
#                       AlgoParams=AlgoParams, 
#                       unfinished=F)
#   
#   
#   Omega_MAP = output$PhysicalValues
#   
#   # Save the output
#   saveRDS(object = output,
#           file = paste0(dir,"IIDIPUS_Results/output_",
#                         DateTimeString(),"10simulated.Rdata"))
#   
#   # Calculate the single linear predictor term that would make each event prediction perfectly accurate
#   modifiers<-SingleEventsModifierCalc(dir,Model,output$PhysicalValues,AlgoParams)
#   # Correlate the modifiers to systemic vulnerability, note that systemic vulnerability variables must be defined in the model
#   vulnerability<-CorrelateModifier(modifiers,Model)
#   
#   return(list(DataDirectory=ODDpaths,
#               Parameterisation=output,
#               fullvulnerability=vulnerability))
# }
