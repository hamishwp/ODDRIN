dir<-directory<-"/home/ubuntu/ODDRIN/"
# Set the working directory from your environment variables
setwd(directory)
packred<-F

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

Omega <- list(Lambda1 = list(nu=0,omega=0.5),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=0.4,omega=0.6),
              zeta = list(k=2.978697, lambda=1.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.2359788),
              eps = list(eps=0.01304351))

AlgoParams$kernel <- 'lognormal' # 'loglaplace'
AlgoParams$AllParallel <- TRUE
AlgoParams$cores <- 8
AlgoParams$NestedCores <- 1
AlgoParams$Np <- 10
AlgoParams$ABC <- 1
AlgoParams$epsilon_max <- AlgoParams$epsilon_min
AlgoParams$smc_Npart <- 100

abcSmc(AlgoParams, Model, unfinished=T, oldtag='2022-07-15_104248')
