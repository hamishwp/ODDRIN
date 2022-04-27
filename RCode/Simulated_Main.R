
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
Omega <- list(Lambda1 = list(nu=1,omega=0.1),
              Lambda2 = list(nu= 0.15, omega=0.75),
              Lambda3 = list(nu=0.7,omega=0.05),
              zeta = list(k=2.978697, lambda=1.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.2359788),
              eps = list(eps=0.01304351))

simulateDataSet(150, Omega, Model=Model, dir = dir)
#After generating the simulated data, need to move 'from 'centerings' 
#and 'ODDobjects' from 'IIDIPUS_SimInput' to 'IIDIPUS_Input'


main_simulated <- function(){
  #ODDpaths<-ExtractData(Model$haz,dir,extractedData)
  # Extract initial parameterisation estimate
  
  iVals<-GetInitVals(dir,Model,AlgoParams)
  # Parameterise... Here we go!
  
  AlgoParams$AllParallel <- TRUE
  AlgoParams$cores <- 6
  AlgoParams$Np <- 5
  AlgoParams$ABC <- 1
  AlgoParams$itermax <- 2500
  
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



output <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/output_2022-04-26_234345')


Omega_MAP <- output[which.max(output[,1]),2:ncol(output)] %>% 
  relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model)


# Plot S-curves for the actual and MAP parameterisation

Intensity <- seq(0,10,0.1)
Dfun<-function(I_ij, theta) h_0(I = I_ij,I0 = 4.5,theta = Omega$theta) 

D_extent <- BinR(Dfun(Intensity, theta=Omega$theta) , Omega$zeta)
D_extent_sample <- BinR(Dfun(Intensity, theta=Omega_MAP$theta) , Omega_MAP$zeta)
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
lines(Intensity, D_extent, col='green', type='l'); lines(Intensity, D_extent_sample, col='green', lty=2, lwd=2)

#check that the likelihood from the true parameters is not greater than the MAP estimate
logTarget(dir = dir,Model = Model,proposed = Omega_MAP, AlgoParams = AlgoParams)
logTarget(dir = dir,Model = Model,proposed = Omega, AlgoParams = AlgoParams)


#Magnitude vs contribution to likelihood
folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
LLs <- c()
#LLs2 <- c()
#LLs3 <- c()
mags <- c()
for(i in 1:length(ufiles)){
  # Extract the ODD object
  ODDy<-readRDS(paste0(folderin,ufiles[i]))
  # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
  ODDy@fIndies<-Model$fIndies
  ODDy@gmax%<>%as.data.frame.list()
  # Apply DispX
  LLs <- append(LLs, DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams))
  #LLs2 <- append(LLs2, DispX(ODD = ODDy,Omega = Omega_MAP,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams))
  #LLs3 <- append(LLs3, DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams))
  hrange<-grep("hazMean",names(ODDy),value = T)
  mags <- append(mags, max(ODDy[hrange]@data, na.rm=TRUE))
}
plot(1:142, LLs)
#points(1:142, LLs2,col='blue')
#points(1:142, LLs3, col='red')

#Misc plots to compare predictions from the true and MAP parameters
# folderin<-paste0(dir,"IIDIPUS_SimInput/ODDobjects/")
# ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
# x_actual = c()
# x_pred_TRUE = c()
# x_pred_MAP = c()
# y_actual = c()
# y_pred_TRUE = c()
# y_pred_MAP = c()
# z_actual = c()
# z_pred_TRUE = c()
# z_pred_MAP = c()
# 
# for(i in 1:length(ufiles)){
#   ODDSim <- readRDS(paste0(folderin,ufiles[i]))
#   #simulate displacement, mortality and building destruction using DispX
#   ODDSim_Omega <- ODDSim %>% DispX(Omega, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
#   ODDSim_MAP <- ODDSim %>% DispX(Omega_MAP, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
#   x_actual = append(x_actual, ODDSim@gmax$gmax)
#   x_pred_TRUE = append(x_pred_TRUE, ODDSim_Omega@predictDisp$disp_predictor)
#   x_pred_MAP = append(x_pred_MAP, ODDSim_MAP@predictDisp$disp_predictor)
#   y_actual = append(y_actual, ODDSim@gmax$mortality)
#   y_pred_TRUE = append(y_pred_TRUE, ODDSim_Omega@predictDisp$mort_predictor)
#   y_pred_MAP = append(y_pred_MAP, ODDSim_MAP@predictDisp$mort_predictor)
#   z_actual = append(z_actual, ODDSim@gmax$buildDestroyed)
#   z_pred_TRUE = append(z_pred_TRUE, ODDSim_Omega@predictDisp$nBD_predictor)
#   z_pred_MAP = append(z_pred_MAP, ODDSim_MAP@predictDisp$nBD_predictor)
# }
# 
# par(mfrow=c(1,1))
# plot(log(x_actual), log(x_actual))
# points(x_actual%>% log(), x_pred_TRUE%>% log(), col='blue')
# points(x_actual %>% log(), x_pred_MAP%>% log(), col='red')
# 
# par(mfrow=c(1,1))
# plot(log(y_actual), log(y_actual), ylim=c(2,10))
# points(y_actual%>% log(), y_pred_TRUE%>% log(), col='blue')
# points(y_actual %>% log(), y_pred_MAP%>% log(), col='red')
# 
# 
# par(mfrow=c(1,1))
# plot(log(z_actual), log(z_actual), ylim=c(6.5,7.5))
# points(z_actual%>% log(), z_pred_TRUE%>% log(), col='blue')
# points(z_actual %>% log(), z_pred_MAP%>% log(), col='red')
# 
# 
# i = 10
# samples = c()
# folderin<-paste0(dir,"IIDIPUS_SimInput/ODDobjects/")
# ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
# for (k in 1:50){
#   ODDSim <- readRDS(paste0(folderin, ufiles[i]))
#   ODDSim_Omega <- ODDSim %>% DispX(Omega, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
#   samples <- append(sampels, ODDSim_Omega@predictDisp$disp_predictor)
# }
# plot(samples, ylim=c(1800,2000))
# abline(h=ODDSim@gmax$gmax)
# 
# 
# par(mfrow=c(2,2))
# for(i in 1:4){
#   plot(output$OptimisationOut[,i+1], type='l')
#   abline(h=unlist(Omega)[i], col='red')

# }

BDSim <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_SimInput/BDobjects/EQ20130409ABC_-3')
BDSim@data
output <- readRDS('/home/manderso/Documents/GitHub/ODDRINfork/IIDIPUS_Results/output_2022-04-25_135511')

LLs1 <- c()
LLs2 <- c()
folderin<-paste0(dir,"IIDIPUS_Input/BDobjects/")
ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
for(i in 1:length(ufiles)){
  # Extract the BD object
  BDy<-readRDS(paste0(folderin,ufiles[i]))
  # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
  BDy@fIndies<-Model$fIndies
  # Apply BDX
  tLL1<-mean(tryCatch(BDX(BD = BDy,Omega = Omega,Model = Model,Method=AlgoParams, LL=T),
                error=function(e) NA))
  tLL2<-mean(tryCatch(BDX(BD = BDy,Omega = Omega_MAP,Model = Model,Method=AlgoParams, LL=T),
                  error=function(e) NA))
  LLs1 <- append(LLs1, tLL1)
  LLs2 <- append(LLs2, tLL2)
  # If all is good, add the LL to the total LL
  # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
}

plot(1:142, LLs1)
points(1:142, LLs2, col='red')