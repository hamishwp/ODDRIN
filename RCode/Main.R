#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%% IIDIPUS - Integrated Internal DIsplaced PopUlation Sampler %%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%% Coded by Hamish Patten %%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%% Collaboration between the University of Oxford %%%%%%%%%%#
#%%%%%%%%%% and the Internal Displacement Monitoring Centre %%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%% Started January 2020 %%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%% Funded by the EPSRC - Impact Acceleration Account %%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

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

IIDIPUSModelTraining<-function(extractedData=T){
 
  # Extract Displacement & Hazard data and match & reduce
  ODDpaths<-ExtractData(Model$haz,dir,extractedData)
  # Extract initial parameterisation estimate
  iVals<-GetInitVals(dir,Model,AlgoParams)
  # Parameterise... Here we go!
  output <- Algorithm(dir=dir,
                      Model=Model,
                      iVals=iVals,
                      AlgoParams=AlgoParams)
  # Save the output
  saveRDS(object = output,
          file = paste0(dir,"IIDIPUS_Results/output_",
                        DateTimeString(),".Rdata"))
  
  # Calculate the single linear predictor term that would make each event prediction perfectly accurate
  modifiers<-SingleEventsModifierCal(dir,Model,output$PhysicalValues,AlgoParams)
  # Correlate the modifiers to systemic vulnerability, note that systemic vulnerability variables must be defined in the model
  vulnerability<-CorrelateModifier(modifiers,Model)
  
  return(list(DataDirectory=ODDpaths,
              Parameterisation=output,
              fullvulnerability=vulnerability))
   
}

ODDRIN<-IIDIPUSModelTraining()
