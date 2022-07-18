
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
Omega <- list(Lambda1 = list(nu=0,omega=0.5),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=0.4,omega=0.6),
              zeta = list(k=2.978697, lambda=1.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.2359788),
              eps = list(eps=0.01304351))

Model$center <- simulateDataSet(100, Omega, Model=Model, dir = dir, outliers = TRUE)
#After generating the simulated data, need to move 'from 'centerings' 
#and 'ODDobjects' from 'IIDIPUS_SimInput' to 'IIDIPUS_Input'



main_simulated <- function(){
  #ODDpaths<-ExtractData(Model$haz,dir,extractedData)
  # Extract initial parameterisation estimate
  
  # Parameterise... Here we go!
  
  AlgoParams$AllParallel <- TRUE
  AlgoParams$cores <- 8
  AlgoParams$Np <- 10
  AlgoParams$ABC <- 0.1
  AlgoParams$itermax <- 10000
  AlgoParams$epsilon_max <- AlgoParams$epsilon_min
  
  iVals<-GetInitVals(dir,Model,AlgoParams)
  AlgoParams$ABC <- 1
  output <- AMCMC(dir=dir,
                      Model=Model,
                      iVals=iVals,
                      AlgoParams=AlgoParams, 
                      unfinished=F)
  
  # output <- AMCMC(dir=dir,
  #                 Model=Model,
  #                 iVals=iVals,
  #                 AlgoParams=AlgoParams, 
  #                 unfinished=T, 
  #                 tag='2022-06-22_141021')
  
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
