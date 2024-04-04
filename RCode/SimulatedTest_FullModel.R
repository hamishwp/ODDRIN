#------------------------------------------------------------------------------------------------
#------------------------------------ SETUP -----------------------------------------------------
#-------------------(copied from GetEnv.R and Main.R for ease)-----------------------------------
#------------------------------------------------------------------------------------------------

# Where is the main folder with all the code and data
dir<-directory<- "/home/manderso/Documents/GitHub/ODDRIN/" #"./" #"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/"
# Set the working directory from your environment variables
setwd(directory)
# Directory of the Data for Good data, e.g. Disaster Mapping, 4G connectivity, etc
FBdirectory<-'/home/patten/Documents/IDMC/Facebook_Data/'
# Do you want only the reduced packages or all? Choose via packred
packred<-F

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

#Run SMC-ABC algorithm with Simulated Data generated using the full model

#Choose true Omega for Simulated Data
Omega <- list(Lambda1 = list(nu=8.5, kappa=1.9049508),
                   Lambda2 = list(nu=9.8626452, kappa=1.8049508),
                   Lambda3 = list(nu=9.3, kappa=1.2),
                   Lambda4 = list(nu=9.9, kappa=1.6),
                   theta= list(theta1=0.6),
                   eps = list(local=2.2292053, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                   vuln_coeff = list(PDens=0.05, SHDI=-0.5, GNIc=-0.1, Vs30=0.1, EQFreq=-0.1, FirstHaz=0.05, Night=0.05, FirstHaz.Night=0.1),
                   check = list(check=0.5))
  
set.seed(1)
simulateDataSet(150, Omega, Model, dir)



#------------------------------------------------------------------------------------------------
#----------------------------------- PLOT SIMULATED DATA ----------------------------------------
#--------------------------------(and compare to real data)--------------------------------------
#------------------------------------------------------------------------------------------------

# Collect mortality, building damage, and displacement data for simulated data:

ODDsim_paths <-na.omit(list.files(path="IIDIPUS_SimInput/ODDobjects/"))
df_SimImpact <- data.frame(observed=numeric(),
                           impact=character(),
                           polygon=integer(),
                           event=integer(),
                           I_max=numeric())
for(i in 1:length(ODDsim_paths)){
  ODDSim <- readRDS(paste0("IIDIPUS_SimInput/ODDobjects/",ODDsim_paths[i]))
  if (length(ODDSim@impact$impact=='mortality')>0){
    for (j in 1:NROW(ODDSim@impact)){
      df_SimImpact %<>% add_row(observed=ODDSim@impact$observed[j], impact=ODDSim@impact$impact[j], polygon=ODDSim@impact$polygon[j],
                                event=i, I_max=max(ODDSim$hazMean1[ODDSim@polygons[[ODDSim@impact$polygon[j]]]$indexes],  na.rm=T))
    }
  }
}
ggplot(df_SimImpact %>% filter(impact=='mortality'), aes(x=I_max, y=observed)) + geom_point()

# Collect mortality, building damage, and displacement data for real data:

ODDpath <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/'
ODDpaths <-na.omit(list.files(path=ODDpath))
df_Impact <- data.frame(observed=numeric(), impact=character(),
                        polygon=integer(), event=integer(), I_max=numeric())
for(i in 1:length(ODDpaths)){
  ODD <- readRDS(paste0(ODDpath,ODDpaths[i]))
  if (length(ODD@impact$impact=='mortality')>0){
    for (j in 1:NROW(ODD@impact)){
      df_Impact %<>% add_row(observed=ODD@impact$observed[j], impact=ODD@impact$impact[j], polygon=ODD@impact$polygon[j],
                                event=i, I_max=max(ODD$hazMean1[ODD@polygons[[ODD@impact$polygon[j]]]$indexes],  na.rm=T))
    }
  }
}

ggplot(df_Impact %>% filter(impact=='mortality'), aes(x=I_max, y=observed)) + geom_point()

#Plots:
# plot of mortality (simulated) vs displacement (simulated) for the same  polygon, compare to true
impact_type='displacement'
ggplot() + 
  geom_histogram(data=df_Impact %>% filter(impact==impact_type), aes(x=log(observed),y=after_stat(count)), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  geom_histogram(data=df_SimImpact %>% filter(impact==impact_type), aes(x=log(observed),y=after_stat(count)), alpha=0.3, col='yellow', lwd=0.2, fill='yellow')

#compare the number of observations per event:
ggplot() +
  geom_histogram(data=df_Impact %>% filter(impact==impact_type) %>% group_by(event) %>% tally(), aes(x=n, y=after_stat(count)), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  geom_histogram(data=df_SimImpact %>% filter(impact==impact_type) %>% group_by(event) %>% tally(), aes(x=n, y=after_stat(count)), alpha=0.3, col='yellow', lwd=0.2, fill='yellow')

#compare correlation between impact types for simulated and true data:
impact_type1 = 'mortality'
impact_type2 = 'displacement'
cor_true <- ggplot(data=merge(df_Impact%>% filter(impact==impact_type1), df_Impact %>% filter(impact==impact_type2), by=c('event', 'polygon')),
       aes(x=observed.x, y=observed.y)) + xlab(impact_type1) + ylab(impact_type2) +
  geom_point() + scale_x_log10() + scale_y_log10() + ggtitle('Simulated Impact Data, with each point representing an observation for one event and region')

cor_sim <- ggplot(data=merge(df_SimImpact%>% filter(impact==impact_type1), df_SimImpact %>% filter(impact==impact_type2), by=c('event', 'polygon')),
       aes(x=observed.x, y=observed.y)) + xlab(impact_type1) + ylab(impact_type2) +
  geom_point() + scale_x_log10() + scale_y_log10() + ggtitle('Observed Impact Data, with each point representing an observation for one event and region')

grid.arrange(cor_true, cor_sim, ncol=1)

# histogram of intensities (simulated) vs intensities (true)
impact_type='displacement'
ggplot() + 
  geom_histogram(data=df_Impact %>% filter(impact==impact_type), aes(x=I_max,y=after_stat(count)), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  geom_histogram(data=df_SimImpact %>% filter(impact==impact_type), aes(x=I_max,y=after_stat(count)), alpha=0.3, col='yellow', lwd=0.2, fill='yellow')



