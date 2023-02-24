
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


#Parameterise the model and simulate the data:
Omega <- list(Lambda1 = list(nu=8,omega=4.9),
              Lambda2 = list(nu= 9.8, omega=5.7),
              Lambda3 = list(nu=7.5,omega=4.1),
              Lambda4 = list(nu=9.4, omega=5.2),
              theta = list(e=0.2359788),
              eps = list(eps=0.01304351),
              vuln_coeff = list(itc=1, PDens=0.01, ExpSchYrs=0,LifeExp=-0.02, GNIc=-0.045, Vs30=-0.001, EQFreq=-0.01)) 
              #intercept term itc is redundant

Model$HighLevelPriors(Omega %>% addTransfParams(), Model)

Model$center <- simulateDataSet(20, Omega %>% addTransfParams(), Model=Model, dir = dir, outliers = FALSE)

#After generating the simulated data, need to move 'from 'centerings' 
#and 'ODDobjects' from 'IIDIPUS_SimInput' to 'IIDIPUS_Input'

main_simulated <- function(){
  #ODDpaths<-ExtractData(Model$haz,dir,extractedData)
  # Extract initial parameterisation estimate
  
  # Parameterise... Here we go!
  
  iVals<-GetInitVals(dir,Model,AlgoParams)

  output <- delmoral_parallel(dir=dir,
                      Model=Model,
                      iVals=iVals,
                      AlgoParams=AlgoParams, 
                      unfinished=F)
  
  
  Omega_MAP = output$PhysicalValues
  
  # Save the output
  saveRDS(object = output,
          file = paste0(dir,"IIDIPUS_Results/output_",
                        DateTimeString(),"10simulated.Rdata"))
  
  # Calculate the single linear predictor term that would make each event prediction perfectly accurate
  modifiers<-SingleEventsModifierCalc(dir,Model,output$PhysicalValues,AlgoParams)
  # Correlate the modifiers to systemic vulnerability, note that systemic vulnerability variables must be defined in the model
  vulnerability<-CorrelateModifier(modifiers,Model)
  
  return(list(DataDirectory=ODDpaths,
              Parameterisation=output,
              fullvulnerability=vulnerability))
}


# -------------------------------------------------------------------------

#Plot the simulated data
#Names: ODDSim.png, Sim DispMortBD.png  
#Size: 1500 x 700
grid.arrange(plotODDy(ODDSim, var='Population') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
             plotODDy(ODDSim, var='GDP')+ xlim(-0.25,0.25) + ylim(-0.25,0.25), 
             plotODDy(ODDSim, var='nBuildings')+ xlim(-0.25,0.25) + ylim(-0.25,0.25), nrow=1)
grid.arrange(plotODDy(ODDSim, var='Disp') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
             plotODDy(ODDSim, var='Mort') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
             plotODDy(ODDSim, var='nBD') + xlim(-0.25,0.25) + ylim(-0.25,0.25), nrow=1)
