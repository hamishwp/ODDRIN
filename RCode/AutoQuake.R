#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%% AutoQuake - automatically generate object from earthquake data %%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%% Developed by Hamish Patten %%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%% Collaboration between the University of Oxford %%%%%%%%%%%%#
#%%%%%%%%%%%% and the Internal Displacement Monitoring Centre %%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%% Funded by the EPSRC - Impact Acceleration Account %%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%% Started March 2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#change the working directory to the location of cloned ODDRIN repository:
dir<-directory<- "/home/manderso/Documents/GitHub/ODDRIN/" #"./" #"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/"
# Set the working directory from your environment variables
setwd(directory)

packred<-F
loadRmpi<-F

# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
# Sourcing the data:
source('RCode/GetData.R')
# Extract model functions and priors
source('RCode/Model.R')
# Extract the model parameterisation algorithm, default = Adaptive MCMC
source('RCode/Method.R')
# Load functions for setting up event objects and sampling from model
source('RCode/AutoQuake_functions.R')


#%%%%%%%%%%%%% User defined input - bare minimum required %%%%%%%%%%%%%%%#
input<-list(
  sdate=as.Date("2019-12-13"), # "YYYY-MM-DD"
  fdate=as.Date("2019-12-17"), # "YYYY-MM-DD"
  iso3="PHL", # Country code in ISO-3C form
  datadir=dir, # Location of the main folder to access the data 
  plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

# Or extract the data purely based on the USGS id number
input<-list(USGSid="usp000huvq",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

#%%%%%%%%%%%%% Variables and functions from IIDIPUS files %%%%%%%%%%%%%%%#
input%<>%append(list(Model=Model, # Model parameters - equations, parameters, ... (Model.R)
                     PosteriorFileLoc='IIDIPUS_Results/mcmc_2025-04-10', # Location of fitted model (AlgoResults)
                     Method=AlgoParams)) # Number of CPUs, particles in SMC, etc. (Method.R)

#%%%%%%% Function to extract data, predict displacement & plot %%%%%%%%%%#
AutoQuake<-function(input,predImpact=T, folder_write='IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/'){
  
  # Collect hazard, population and vulnerability data:
  options(timeout = 500)
  ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
  ODDy = ODDy_with_namer$ODDy
  namer = ODDy_with_namer$namer
  
  if(!predImpact) return(ODDy)
  
  ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                            AlgoResultsFilename=input$PosteriorFileLoc,
                            folder_write=folder_write,
                            namer=namer,
                            N_samples = 100)
  
  #sampled_full = ODDy_with_sampled$sampled_full
  #ODDy = ODDy_with_sampled$ODDyWithImpact
  #namer = ODDy_with_sampled$namer
  
  # Make plots and store them in a specific folder (ODDobj.R)
  
  # Makes plots and saves them in folder_write/input$plotdir
  pout<-MakeODDPredPlots(ODDy_with_sampled, 
                         input = input,
                         folder_write=folder_write, 
                         namer = namer)
  
  pout
  
  return(ODDy)
}

# RUN IT! (Cross your fingers that the EQ exist in USGS)
ODDy<-AutoQuake(input)


