#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%% IIDIPUS - Integrated Internal DIsplaced PopUlation Sampler %%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%% Coded by Hamish Patten and Max Anderson Loake %%%%%%%%%%#
#%%%%%%%%%% Collaboration between the University of Oxford %%%%%%%%%%#
#%%%%%%%%%% and the Internal Displacement Monitoring Centre %%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%% Started January 2020 %%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%% Funded by the EPSRC - Impact Acceleration Account %%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#change the working directory to the location of cloned ODDRIN repository:
dir<-directory<- "/home/manderso/Documents/GitHub/ODDRIN/" #"./" #"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/"
# Set the working directory from your environment variables
setwd(directory)

packred<-T
loadRmpi<-F

# Extract Environment Variables
# source('RCode/GetEnv.R')
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
# Sourcing the data:
source('RCode/GetData.R')
# Extract model functions and priors
source('RCode/Model.R')
# Extract the model parameterisation algorithm, default = Adaptive MCMC
source('RCode/Method.R')

# Methodology parameters required
AlgoParams<-list(#Data:
  input_folder = 'IIDIPUS_Input_Alternatives/Apr25Agg/', # Directory of input data, allows you to vary between simulated input, aggregated input, etc.
  
  #Compute Parameters:
  Np=1, # Number of Monte Carlo particles
  cores=20, # Number of parallelised threads per event
  NestedCores=6, # How many cores are to be used inside the ODD displacement calculations?
  n_nodes=1, #parallelise between nodes on HPC
  AllParallel=T, # Do you want fully-nested (data) parallelisation?
  
  # Loss function:
  ABC=0, # Approximate Bayesian Computation rejection
  kernel='energy_score', #Distance kernel between simulated and observed data
  learning_rate = 40, #power of likelihood 
  impact_weights=list(displacement=1,mortality=7,buildDam=0.8), #Relative weights of the different impact types when obtaining the energy score
  m_CRPS = 100, # number of draws to estimate CRPS/Energy Score for each particle. Number of samples from model therefore becomes Np * m_CRPS
  BD_weight = 0, #scales from 0 to 1: 0 means point data has no weighting, 1 means point data is included to same weighting (by inverse MAD) as mortality
  log_offset = 10, #where log_offset = k, we compare log(x+k) with log(y+k) for simulated data x with observed data y
  W_rankhist = 0.05, #weighting of Anderson Darling statistic over Band Depth rank histograms vs energy score
  
  # SMC hyperparameters: 
  tol0=12000, # should be larger than largest expected distance
  smc_steps = 200, #Number of steps in the ABC-SMC algorithm
  smc_Npart = 1000, #Number of particles in the ABC-SMC algorithm
  smc_alpha = 0.9, #Proportion of particles that we aim to keep 'alive' between steps
  propCOV_multiplier = 0.2, #mutiplied by `optimal' covariance obtained using Filippi et al. 
  
  # Adaptive MCMC hyperparameters:
  N_steps = 10000, # number of MCMC steps
  rho = 0.95, # correlation between subsequent event-level error terms
  rho_local = 0.95, # correlation between subsequent local error terms
  lambda_rate = 0.9, #parameter of adaptive scaling and shaping algorithm
  alpha_star = 0.234, #parameter of adaptive scaling and shaping algorithm
  v_0 = 15 #parameter of adaptive scaling and shaping algorithm
)

IIDIPUSModelTraining<-function(extractedData=T){
 
  # Extract Displacement & Hazard data and match & reduce
  ODDpaths<-ExtractData(Model$haz,dir,extractedData)
  # Extract initial parameterisation estimate
  iVals<-GetInitVals(dir,Model,AlgoParams,usePastPost=T)
  
  tag_notes = paste0('_Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot')
  # Parameterise... Here we go!
  output <- AMCMC(dir=dir,
                  AlgoParams=AlgoParams,
                  Model=Model,
                  iVals=iVals,
                  unfinished = F,
                  oldtag=NULL, 
                  tag_notes=NULL)
  
  #Alternatively, using ABCSMC:
  # output <- ABCSMC(dir=dir,
  #                 AlgoParams=AlgoParams,
  #                 Model=Model)
  
  return(list(DataDirectory=ODDpaths,
              Parameterisation=output))
  
}

ODDRIN<-IIDIPUSModelTraining()
