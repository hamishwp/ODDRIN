dir<-directory<- "/home/manderso/Documents/GitHub/ODDRIN/"
#dir<-directory<- "/home/kebl6973/ODDRIN/" #"./" #"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/"
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
AlgoParams$m_CRPS = 1
AlgoParams$smc_Npart = 1
AlgoParams$cores = 1
AlgoParams$n_nodes = 1

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

#out <- delmoral_parallel(AlgoParams, Model, unfinished=F, oldtag='')

#ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Train/EQ20220701IRN_166')
#DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F, Method = AlgoParams)
 
AlgoParams$cores <- 4
AlgoParams$NestedCores <- 1
AlgoParams$Np <- 2
AlgoParams$m_CRPS <- 5
AlgoParams$smc_Npart <- 500
AlgoParams$n_nodes <- 25
AlgoParams$smc_steps <- 200
start_time <- Sys.time()
out <- delmoral_parallel(AlgoParams, Model, unfinished=T, oldtag='HPC/abcsmc_2023-11-01_133620')
#impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams)
finish_time <- Sys.time()
print(finish_time-start_time)
 
#CalcDist(impact_sample, AlgoParams)

