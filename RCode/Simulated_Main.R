
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
Omega <- list(Lambda1 = list(nu=0.9,omega=0.1),
              Lambda2 = list(nu= 0.3, omega=0.7),
              Lambda3 = list(nu=0.8,omega=0.1),
              zeta = list(k=2.978697, lambda=1.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.2359788),
              eps = list(eps=0.01304351))

simulateDataSet(10, Omega, Model=Model, dir = dir)
#After generating the simulated data, need to move 'from 'centerings' 
#and 'ODDobjects' from 'IIDIPUS_SimInput' to 'IIDIPUS_Input'


main_simulated <- function(){
  #ODDpaths<-ExtractData(Model$haz,dir,extractedData)
  # Extract initial parameterisation estimate
  
  iVals<-GetInitVals(dir,Model,AlgoParams)
  # Parameterise... Here we go!
  
  AlgoParams$AllParallel <- TRUE
  AlgoParams$ABC <- 2
  AlgoParams$itermax <- 80
  
  output <- Algorithm(dir=dir,
                      Model=Model,
                      iVals=iVals,
                      AlgoParams=AlgoParams)
  
  Omega_MAP = output$PhysicalValues
  
  # Save the output
  saveRDS(object = output,
          file = paste0(dir,"IIDIPUS_Results/output_",
                        DateTimeString(),"simulated.Rdata"))
  
  # Calculate the single linear predictor term that would make each event prediction perfectly accurate
  modifiers<-SingleEventsModifierCalc(dir,Model,output$PhysicalValues,AlgoParams)
  # Correlate the modifiers to systemic vulnerability, note that systemic vulnerability variables must be defined in the model
  vulnerability<-CorrelateModifier(modifiers,Model)
  
  return(list(DataDirectory=ODDpaths,
              Parameterisation=output,
              fullvulnerability=vulnerability))
}

# -----------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------- MISCELLANEOUS PLOTS -------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

#Plot the simulated data
#Names: ODDSim.png, Sim DispMortBD.png  
#Size: 1500 x 700
grid.arrange(plotODDy(ODDSim, var='Population') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
             plotODDy(ODDSim, var='GDP')+ xlim(-0.25,0.25) + ylim(-0.25,0.25), 
             plotODDy(ODDSim, var='nBuildings')+ xlim(-0.25,0.25) + ylim(-0.25,0.25), nrow=1)
grid.arrange(plotODDy(ODDSim, var='Disp') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
             plotODDy(ODDSim, var='Mort') + xlim(-0.25,0.25) + ylim(-0.25,0.25), 
             plotODDy(ODDSim, var='nBD') + xlim(-0.25,0.25) + ylim(-0.25,0.25), nrow=1)


output <- readRDS('/home/manderso/Documents/GitHub/ODDRINfork/IIDIPUS_Results/output_2022-04-15_220417')

Intensity <- seq(0,10,0.1)
Dfun<-function(I_ij, theta) h_0(I = I_ij,I0 = 4.5,theta = Omega$theta) 

# Plot S-curves for the actual and MAP parameterisation
D_MortDisp <- BinR( Omega$Lambda1$nu * Dfun(Intensity, theta=Omega$theta) + Omega$Lambda1$omega, Omega$zeta)
D_MortDisp_sample <- BinR( Omega_MAP$Lambda1$nu * Dfun(Intensity, theta=Omega_MAP$theta) + Omega_MAP$Lambda1$omega, Omega_MAP$zeta)
D_Mort <- BinR(Omega$Lambda2$nu * Dfun(Intensity, theta=Omega$theta) + Omega$Lambda2$omega , Omega$zeta) * D_MortDisp
D_Mort_sample <- BinR(Omega_MAP$Lambda2$nu * Dfun(Intensity, theta=Omega_MAP$theta) + Omega_MAP$Lambda2$omega , Omega_MAP$zeta) * D_MortDisp_sample
D_Disp <- D_MortDisp - D_Mort
D_Disp_sample <- D_MortDisp_sample - D_Mort_sample
D_BD <- BinR(Omega$Lambda3$nu * Dfun(Intensity, theta=Omega$theta) + Omega$Lambda3$omega, Omega$zeta)
D_BD_sample <- BinR(Omega_MAP$Lambda3$nu * Dfun(Intensity, theta=Omega_MAP$theta) + Omega_MAP$Lambda3$omega, Omega_MAP$zeta)
plot(Intensity, D_Mort, col='red', type='l', ylim=c(0,1)); lines(Intensity, D_Mort_sample, col='red', lty=2)
lines(Intensity, D_Disp, col='blue'); lines(Intensity, D_Disp_sample, col='blue', lty=2)
lines(Intensity, D_BD, col='pink', type='l'); lines(Intensity, D_BD_sample, col='pink', lty=2, lwd=2)

logTarget(dir = dir,Model = Model, proposed = Omega_MAP,AlgoParams = AlgoParams)
logTarget(dir = dir,Model = Model, proposed = Omega,AlgoParams = AlgoParams)
