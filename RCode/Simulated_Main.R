  
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
Omega <- list(Lambda1 = list(mu=9, sigma=1.1),
              Lambda2 = list(mu=10.8, sigma=0.88),
              Lambda3 = list(mu=8.5, sigma=1.1),
              Lambda4 = list(mu=9.8, sigma=0.85),
              eps = list(local=0.05, hazard=0.1),
              vuln_coeff = list(PDens=0, AveSchYrs=0,LifeExp=0, GNIc=0.03, Vs30=0, EQFreq=0),
              check = list(check=0.5)) 

Model$HighLevelPriors(Omega %>% addTransfParams(), Model)

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_June24_AggFactor5/ODDobjects/EQ20150425NPL_31')
DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)



start_time <- Sys.time()
sampleBDDist(dir,Model, Omega %>% addTransfParams(),AlgoParams)
finish_time <- Sys.time()
finish_time - start_time


sum(elapsed_time)


ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_June24_AggFactor5/ODDobjects/EQ20131015PHL_16')
DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_June20/ODDobjects/EQ20151026AFG_170')
DispX(ODD = ODDy,Omega = Omega %>% addTransfParams(),center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams)

sampleBDDist(dir,Model,proposed = Omega %>% addTransfParams(),AlgoParams,expLL=T)
  
# start_time <- Sys.time()
dist_sample <- sampleDist(dir = dir,Model = Model,
                          proposed = Omega %>% addTransfParams(),
                          AlgoParams = AlgoParams)




d_i <- logTarget2(dist_sample, AlgoParams)
# end_time <- Sys.time()
# print(paste('Time:', end_time-start_time))

Model$center <- simulateDataSet(5, Omega %>% addTransfParams(), Model=Model, dir = dir, outliers = FALSE)

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