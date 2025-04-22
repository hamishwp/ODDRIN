
#change the working directory to the location of cloned ODDRIN repository:
dir<-directory<- "/home/manderso/Documents/ODDRIN-master/" #"./" #"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/"
# Set the working directory from your environment variables
setwd(directory)

packred<-T
loadRmpi<-F

# Extract Environment Variables
#source('RCode/GetEnv.R')
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
# Extract model functions and priors
source('RCode/Model.R')
# Extract the model parameterisation algorithm, default = Adaptive MCMC
source('RCode/Method.R')

# This file holds all the required data - hazard, exposure, vulnerability, as well as event information and observed displacement estimates in the 'ODD' class format
ODDy<-readODD(paste0(dir,"IIDIPUS_Input/ODDobjects/EQ20210814HTI_164_example"))

# Plot the ODD object:
plot(ODDy) 

# Retrieve a past posterior fit, trained on earthquakes at a global level
AlgoResults <- readRDS(paste0(dir,"IIDIPUS_Results/mcmc_2025-04-10"))

# Sample a posterior parameter set:
post_samples <- t(AlgoResults$Omega_sample_phys)
s_finish = which(is.na(post_samples)[,1])[1] - 1 #find stopping point
post_samples = post_samples[max(1, s_finish-1000):s_finish,] #take the last 1000 samples as posterior samples
Omega = post_samples[sample(1:nrow(post_samples),1),] %>% relist(skeleton=Model$skeleton) %>% addTransfParams()

# Test to see if the displacement prediction calculations are working
ODDy%<>%DispX(Omega = Omega,center = Model$center, Method = AlgoParams, output='ODDyWithSampled')

plot(ODDy$SampledMortality)
plot(ODDy$SampledDisplacement)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Now check ability to download new ODD objects (and run Autoquake.R)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

packred<-F
loadRmpi<-F

source('RCode/GetODDPackages.R')
source('RCode/AutoQuake.R')

input<-list(
  sdate=as.Date("2019-12-13"), # "YYYY-MM-DD"
  fdate=as.Date("2019-12-17"), # "YYYY-MM-DD"
  iso3="PHL", # Country code in ISO-3C form
  datadir=dir, # Location of the main folder to access the data 
  plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

ODDy <- AutoQuake(input)




