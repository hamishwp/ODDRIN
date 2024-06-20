#dir<-directory<- "/home/manderso/Documents/GitHub/ODDRIN/"
dir<-directory<- "/home/kebl6973/ODDRIN/" #"./" #"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/"
# Set the working directory from your environment variables
setwd(directory)

packred<-T
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

AlgoParams$Np = 1
AlgoParams$m_CRPS = 60
AlgoParams$smc_Npart = 1
AlgoParams$cores = 1
AlgoParams$n_nodes = 1

#Parameterise the model and simulate the data:

Omega <- Omega_true <- list(Lambda1 = list(nu=9.05, kappa=1.2),
                            Lambda2 = list(nu=9.7, kappa=1.25), #list(nu=10.65, kappa=1.5), #
                            Lambda3 = list(nu=8.9, kappa=1.1),
                            Lambda4 = list(nu=9.9, kappa=1.6),
                            theta= list(theta1=0.6),
                            eps=list(local=0.6, hazard_mort=0.3, hazard_disp=0.5, hazard_bd=0.4, hazard_cor=0.55),
                            #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                            vuln_coeff = list(PDens=0, SHDI=-0.18, GNIc=-0.05, Vs30=0.1, EQFreq=-0.25, FirstHaz=0.05, Night=0, FirstHaz.Night=0.1),
                            check = list(check=0.5))

Model$HighLevelPriors(Omega %>% addTransfParams(),Model)

#out <- delmoral_parallel(AlgoParams, Model, unfinished=F, oldtag='')

#ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Train/EQ20220701IRN_166')
#DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F, Method = AlgoParams)
 
AlgoParams$cores <- 4
AlgoParams$NestedCores <- 1
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 60
AlgoParams$smc_Npart <- 25
AlgoParams$n_nodes <- 39
AlgoParams$smc_steps <- 100
AlgoParams$rel_weightings <- c(0,1)
AlgoParams$kernel_sd <- list(displacement = 1, mortality = 7, buildDam=0.6,
                             buildDest = 0, buildDamDest = 0)

start_time <- Sys.time()

tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_M', AlgoParams$m_CRPS, '_Npart', AlgoParams$smc_Npart, '_150events_simulatedfull_15by15')
out <- delmoral_parallel(AlgoParams, Model, unfinished=F, tag_notes=tag_notes)
#impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams)
finish_time <- Sys.time()
print(finish_time-start_time)
 
#CalcDist(impact_sample, AlgoParams)

