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
#source('RCode/GetData.R')
# Extract model functions and priors
source('RCode/Model.R')
# Extract the model parameterisation algorithm, default = Adaptive MCMC
source('RCode/Method.R')
#Extract the functions for generating the simulations
source('RCode/Simulate.R')

#Run SMC-ABC algorithm with Simulated Data generated using the full model

#Choose true Omega for Simulated Data

Omega <- Omega_true <-  list(Lambda1 = list(nu=7.5, kappa=0.6),
                                    Lambda2 = list(nu=11.4, kappa=0.775), #list(nu=10.65, kappa=1.5), #
                                    Lambda3 = list(nu=7.9, kappa=0.68),
                                    #Lambda4 = list(nu=9.9, kappa=1.6),
                                    #theta= list(theta1=0.6),
                                    eps=list(local=1.2, hazard_mort=0.4, hazard_disp=0.65, hazard_bd=0.5, hazard_cor=0.1),
                                    #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                                    vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15))

Model$HighLevelPriors(Omega %>% addTransfParams(), Model)

# Omega <- Omega_true <-  list(Lambda1 = list(nu=7.5, kappa=0.6),
#                                     Lambda2 = list(nu=10.4, kappa=0.775), #list(nu=10.65, kappa=1.5), #
#                                     Lambda3 = list(nu=7.9, kappa=0.68),
#                                     #Lambda4 = list(nu=9.9, kappa=1.6),
#                                     #theta= list(theta1=0.6),
#                                     eps=list(local=0.3, hazard_mort=0.55, hazard_disp=0.65, hazard_bd=0.5, hazard_cor=0.1),
#                                     #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
#                                     vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15))
#   
  
  
# Omega = list(Lambda1 = list(nu=7.5, kappa=0.6),
#              Lambda2 = list(nu=10.7, kappa=0.775), #list(nu=10.65, kappa=1.5), #
#              Lambda3 = list(nu=7.9, kappa=0.68),
#              #Lambda4 = list(nu=9.9, kappa=1.6),
#              #theta= list(theta1=0.6),
#              eps=list(local=0.4, hazard_mort=0.45, hazard_disp=0.65, hazard_bd=0.5, hazard_cor=0.65),
#              #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
#              vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15))
                            # check = list(check=0.5))


HLPrior_samples <- readRDS(paste0(dir, 'IIDIPUS_Input/HLPriorSamples_MCMCOut'))
propCOV <- cov(HLPrior_samples)/5
init_val_phys <- Proposed2Physical(HLPrior_samples[1,] %>% relist(skeleton=Model$skeleton) %>% unlist(), Model)
xx <- SampleImpact(dir, Model, init_val_phys %>% addTransfParams(), AlgoParams)


Model$HighLevelPriors(Omega %>% addTransfParams(), Model)
plot_S_curves(Omega_true)

Model$HighLevelPriors(Omega_ln %>% addTransfParams(), Model)

# plot(seq(5,10,0.05),log(pnorm(seq(5,10,0.05)-4.5, 4,0.5)), type='l')
# points(seq(5,10,0.05),log(pnorm(exp(0.3*seq(5,10,0.05)), 12.1, 1.3)), col='blue', type='l')
# points(seq(5,10,0.05), log(pnorm(exp(0.5*seq(5,10,0.05)), 65.1, 11)), col='red', type='l')
# points(seq(5,10,0.05), log(pnorm(exp(1*seq(5,10,0.05)), 3565.1, 850)), col='pink', type='l')
# points(seq(5,10,0.05),log(pnorm(exp(0.01*seq(5,10,0.05)), 1.087, 0.005)), col='green', type='l')
# plot(seq(5,10,0.05), exp(0.9*seq(5,10,0.05)))

Model$HighLevelPriors(Omega %>% addTransfParams(), Model)

#MCMC:
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput/' #'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg5/'
AlgoParams$N_steps <- 1000
AlgoParams$learning_rate <- 1000
AlgoParams$m_CRPS <- 60
AlgoParams$N_steps <- 1000
AlgoParams$learning_rate <- 40
AlgoParams$Np <- 1
AlgoParams$rel_weightings <- c(1,0)
AlgoParams$rho <- 0.95
#AlgoResults <- readRDS(paste0(dir, 'IIDIPUS_Results/abcsmc_2024-07-29_140215_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2'))
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-29_140215_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2')
init_val_phys <- AlgoResults$Omega_sample_phys[401,,130] %>% relist(skeleton=Model$skeleton)
propCOV <- cov(AlgoResults$Omega_sample[,,130])/4


tag_notes <- paste0('MCMC_RealAgg5_Trial')
AlgoResults <- correlated_MCMC(AlgoParams, Model, unfinished = F, propCOV = propCOV, init_val_phys = init_val_phys, tag_notes=tag_notes)

tag_notes <- paste0('MCMC_RealAgg5_Trial')
AlgoResults <- correlated_MCMC(AlgoParams, Model, unfinished = T, propCOV = propCOV, init_val_phys = init_val_phys, tag_notes=tag_notes,
                               oldtag='mcmc_2024-08-05_155604.382051_MCMC_RealAgg5_Trial')


# Generate samples from higher level prior to initialise MCMC:
n_samples <- 500
samples <- array(NA, dim=c(n_samples, length(unlist(Omega))))
for (i in 1:n_samples){
  if (i %% 10 == 0){
    print(paste('Sample:',i))
  }
  samples[i,] <- HLPrior_sample(Model, AlgoParams)
}
saveRDS(samples, paste0(dir, 'IIDIPUS_Results/HLPriorSamples'))

diag(cov(samples))




set.seed(1)
simulateDataSet(167, Omega, Model, dir, folder_write='IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9/')


AlgoParams$smc_steps <- 2
AlgoParams$smc_Npart <- 250
AlgoParams$m_CRPS <- 60
AlgoParams$Np <- 1
AlgoParams$smc_alpha <- 0.9
AlgoParams$rel_weightings <- c(1,0)


AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-24_121215_alpha0.95_M60_Npart990RealAgg3')
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg3/'
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg3/'
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, M=200, output='SampledTotal')

AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput/'

tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_RealAgg5_asprev')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F, tag_notes=tag_notes)


Model$Priors <- list( #All uniform so currently not included in the acceptance probability. 
  Lambda1=list(nu=list(dist='unif', min=6.5, max=10.5), 
               kappa=list(dist='unif', min=0.25, max=2) #, alpha=list(dist='unif', min=-0.1, max=0.5)
  ), 
  Lambda2=list(nu=list(dist='unif', min=9, max=12.5), 
               kappa=list(dist='unif', min=0.25, max=2)),
  Lambda3=list(nu=list(dist='unif', min=6.5, max=10), 
               kappa=list(dist='unif', min=0.25, max=2)),
  Lambda4=list(nu=list(dist='unif', min=8, max=12.5), 
               kappa=list(dist='unif', min=0.25, max=2.5)),
  theta=list(theta1=list(dist='unif', min=0, max=1)),
  eps=list(local=list(dist='unif', min=0.01, max=2),
           hazard_mort=list(dist='unif', min=0, max=1.5),
           hazard_disp=list(dist='unif', min=0, max=1.5),
           hazard_bd=list(dist='unif', min=0, max=1.5),
           hazard_cor=list(dist='unif', min=0, max=1)),
  vuln_coeff=list(PDens=list(dist='laplace', location=0, scale=0.25),
                  SHDI=list(dist='laplace', location=0, scale=0.25),
                  GNIc=list(dist='laplace', location=0, scale=0.25),
                  Vs30=list(dist='laplace', location=0, scale=0.25),
                  EQFreq=list(dist='laplace', location=0, scale=0.25),
                  #Mag=list(dist='laplace', location=0, scale=0.25),
                  FirstHaz=list(dist='laplace', location=0, scale=0.25),
                  Night=list(dist='laplace', location=0, scale=0.25),
                  FirstHaz.Night=list(dist='laplace', location=0, scale=0.25)),
  check=list(check=list(dist='unif', min=0, max=1))
)

Model$links <-Model$unlinks <- Model$acceptTrans <- Model$skeleton

#Currently, all parameters use a lower and upper bound
for (i in 1:length(Model$links)){
  if (is.list(Model$links[[i]])){
    for (j in 1:length(Model$links[[i]])){
      Model$links[[i]][[j]] <- 'ab_bounded'
      Model$unlinks[[i]][[j]] <- 'ab_bounded_inv'
      Model$acceptTrans[[i]][[j]] <- 'ab_bounded_acc'
    }
  } else {
    Model$links[[i]] <- 'ab_bounded'
    Model$unlinks[[i]] <- 'ab_bounded_inv'
    Model$acceptTrans[[i]] <- 'ab_bounded_acc'
  }
}

#Set lower and upper bounds for the parameters
Model$par_lb <- c()
Model$par_ub <- c()

for (i in 1:length(Model$Priors)){
  if (is.list(Model$Priors[[i]])){
    for (j in 1:length(Model$Priors[[i]])){
      if(Model$Priors[[i]][[j]]$dist == 'unif'){
        Model$par_lb = c(Model$par_lb, Model$Priors[[i]][[j]]$min)
        Model$par_ub = c(Model$par_ub, Model$Priors[[i]][[j]]$max)
      } else if (Model$Priors[[i]][[j]]$dist == 'norm'){
        Model$par_lb = c(Model$par_lb, Model$Priors[[i]][[j]]$mean - 6 * Model$Priors[[i]][[j]]$sd)
        Model$par_ub = c(Model$par_ub, Model$Priors[[i]][[j]]$mean + 6 * Model$Priors[[i]][[j]]$sd)
      } else if (Model$Priors[[i]][[j]]$dist == 'laplace'){
        Model$par_lb = c(Model$par_lb, Model$Priors[[i]][[j]]$location - 15 * Model$Priors[[i]][[j]]$scale)
        Model$par_ub = c(Model$par_ub, Model$Priors[[i]][[j]]$location + 15 * Model$Priors[[i]][[j]]$scale)
      } else {
        stop('Please update Method.R to adjust acceptance probability to account for other priors before continuing.')
      }
    }
  } else {
    if(Model$Priors[[i]]$dist == 'unif'){
      Model$par_lb = c(Model$par_lb, Model$Priors[[i]]$min)
      Model$par_ub = c(Model$par_ub, Model$Priors[[i]]$max)
    } else if (Model$Priors[[i]]$dist == 'norm'){
      Model$par_lb = c(Model$par_lb, Model$Priors[[i]]$mean - 6 * Model$Priors[[i]]$sd)
      Model$par_ub = c(Model$par_ub, Model$Priors[[i]]$mean + 6 * Model$Priors[[i]]$sd)
    } else if (Model$Priors[[i]]$dist == 'laplace'){
      Model$par_lb = c(Model$par_lb, Model$Priors[[i]]$location - 15 * Model$Priors[[i]]$scale)
      Model$par_ub = c(Model$par_ub, Model$Priors[[i]]$location + 15 * Model$Priors[[i]]$scale)
    } else {
      stop('Please update Method.R to adjust acceptance probability to account for other priors before continuing.')
    }
  }
}

addTransfParams <- function(Omega, I0=Model$I0){
  #Omega$theta$theta1 <- 0.6
  Omega$Lambda1$loc <- h_0(Omega$Lambda1$nu, I0, Omega)
  Omega$Lambda2$loc <- h_0(Omega$Lambda2$nu, I0, Omega)
  Omega$Lambda3$loc <- h_0(Omega$Lambda3$nu, I0, Omega)
  Omega$Lambda4$loc <- h_0(Omega$Lambda4$nu, I0, Omega)
  h_10_minus_h_4.5 = h_0(10, I0, Omega) - h_0(4.5, I0, Omega)
  Omega$Lambda1$scale <- h_10_minus_h_4.5 / (6 * Omega$Lambda1$kappa)
  Omega$Lambda2$scale <- h_10_minus_h_4.5 / (6 * Omega$Lambda2$kappa)
  Omega$Lambda3$scale <- h_10_minus_h_4.5 / (6 * Omega$Lambda3$kappa)
  Omega$Lambda4$scale <- h_10_minus_h_4.5 / (6 * Omega$Lambda4$kappa)
  Omega$vuln_coeff_adj <- lapply(Omega$vuln_coeff, function(x) x * Omega$Lambda2$scale)
  Omega$eps_adj <- lapply(Omega$eps, function(x) x * Omega$Lambda2$scale)
  Omega$eps_adj$local <- Omega$eps$local / Omega$eps$hazard_mort
  return(Omega)
}

tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_RealAgg5_localerrupd')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F, tag_notes=tag_notes)



tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulatedfull_energyscore_150events_Npart50_M60')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F, tag_notes=tag_notes)



#Speed Test:

AlgoParams$cores <- 1
AlgoParams$NestedCores <- 4
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 5
AlgoParams$smc_Npart <- 50
AlgoParams$n_nodes <- 1
AlgoParams$smc_steps <- 100
AlgoParams$rel_weightings <- c(1,1)
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/Nov24Agg/'

tag_notes <- paste0('alpha', AlgoParams$smc_alpha, 'test_ucorr')
AlgoResults <- delmoral_parallel_corr(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)


start_time <- Sys.time()
impact_sample <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams)
CalcDist(impact_sample, AlgoParams)
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)

start_time <- Sys.time()
AlgoParams$Np <- 60
#ODD <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Train/EQ20180114PER_78')
impactsamp <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams)
AlgoParams$Np <- 1
end_time <- Sys.time()
execution_time <- end_time - start_time
execution_time

HLPrior_sample <- saveRDS(readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-12-10_060355_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2_further')$Omega_sample[,,170], '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/HLPriorSamples_SMCOut')


#------------------------------------------------------------------------------------------------
#----------------------------------- PLOT SIMULATED DATA ----------------------------------------
#--------------------------------(and compare to real data)--------------------------------------
#------------------------------------------------------------------------------------------------

moveTestData <- function(folder_in='IIDIPUS_Input_Alternatives/IIDIPUS_SimInput'){
  ODD_folderall<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_foldertest<-paste0(dir, folder_in, '/ODDobjects/Test/')
  ufiles<-list.files(path=ODD_folderall,pattern=Model$haz,recursive = T,ignore.case = T)
  i <- 1
  
  # total_mortalities <- c()
  # for (file in ufiles){
  #   ODD <- readODD(paste0(ODD_folderall, file))
  #   total_mortalities <- c(total_mortalities, max(values(ODD[[grep('hazMean', names(ODD))]])[values(!is.na(ODD$Population) & ODD$Population > 0),], na.rm=t))
  # }
  # ufiles <- ufiles[order(total_mortalities, decreasing=T)] # sort by date
  # 
  # for (i in 1:length(ufiles)){
  #   file <- ufiles[i]
  #   if (i %%3 != 1){next}
  #   options(warn = 2)
  #   file.copy(from = paste0(ODD_folderall, file),
  #             to = paste0(ODD_foldertest, file))
  #   file.remove(from = paste0(ODD_folderall, file))
  #   options(warn = 1)
  # }
  
  for (file in ufiles){
    i <- i + 1
    if (i %%3 != 2){next}
    options(warn = 2)
    file.copy(from = paste0(ODD_folderall, file),
              to = paste0(ODD_foldertest, file))
    file.remove(from = paste0(ODD_folderall, file))
    options(warn = 1)
  }
  # BD_folderall<-paste0(dir, folder_in, '/BDobjects/')
  # BD_foldertest<-paste0(dir, folder_in, '/BDobjects/Test/')
  # ufiles<-list.files(path=BD_folderall,pattern=Model$haz,recursive = T,ignore.case = T)
  # i <- 0
  # for (file in ufiles){
  #   i <- i + 1
  #   if (i %%3 != 0){next}
  #   file.copy(from = paste0(BD_folderall, file),
  #             to = paste0(BD_foldertest, file))
  #   file.remove(from = paste0(BD_folderall, file))
  # }
  
}

moveTestData('IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9')

# Collect mortality, building damage, and displacement data for simulated data:
ODDsim_paths <-na.omit(list.files(path="IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9/ODDobjects/", recursive=T))
df_SimImpact <- data.frame(observed=numeric(),
                           impact=character(),
                           polygon=integer(),
                           exposure=numeric(),
                           event=integer(),
                           I_max=numeric())
nHazSim <- c()
maxIntSim <- c()
for(i in 1:length(ODDsim_paths)){
  ODDSim <- readODD(paste0("IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9/ODDobjects/",ODDsim_paths[i]))
  if (length(ODDSim@impact$impact)>0){
    nHazSim <- c(nHazSim, length(grep('hazMean', names(ODDSim))))
    maxIntSim <- c(maxIntSim, max(values(ODDSim[[grep('hazMean', names(ODDSim))]]),  na.rm=T))
    # if (length(grep('hazMean', names(ODDSim@data))) ==1 & (length(grep('-4', ODDsim_paths[i]))==0)& (length(grep('-5', ODDsim_paths[i]))==0)){
    #   stop()
    # }
    if(!is.finite(maxIntSim[length(maxIntSim)])) stop()
    for (j in 1:NROW(ODDSim@impact)){
      df_SimImpact %<>% add_row(observed=ODDSim@impact$observed[j], impact=ODDSim@impact$impact[j], polygon=ODDSim@impact$polygon[j],
                                exposure=ifelse(impact=='buildDam', sum(ODDSim[['nBuildings']][ODDSim@polygons[[ODDSim@impact$polygon[j]]]$indexes]), sum(ODDSim[['Population']][ODDSim@polygons[[ODDSim@impact$polygon[j]]]$indexes])),
                                event=i, I_max=max(values(ODDSim[[grep('hazMean', names(ODDSim))]])[ODDSim@polygons[[ODDSim@impact$polygon[j]]]$indexes,],  na.rm=T))
    }
  }
}
ggplot(df_SimImpact %>% filter(impact=='mortality'), aes(x=I_max, y=observed)) + geom_point()

# Collect mortality, building damage, and displacement data for real data:
#ODDpath <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/'

ODDpath <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg5/ODDobjects/'
ODDpaths <-na.omit(list.files(path=ODDpath, recursive=T))
df_Impact <- data.frame(observed=numeric(), impact=character(),
                        polygon=integer(), exposure=numeric(), event=integer(), I_max=numeric())

# EQFreq_min <- c()
# EQFreq_max <- c()

nHazReal <- c()
maxIntReal <- c()

map_data_df <- data.frame(lon=numeric(), lat=numeric(), max_mmi=numeric())

for(i in 1:length(ODDpaths)){
  ODD <- readRDS(paste0(ODDpath,ODDpaths[i]))
  loc_mmi <- c(ODD@coords[which(ODD@data[,grep('hazMean',names(ODD@data))]==max(ODD@data[,grep('hazMean', names(ODD@data))], na.rm=T), arr.ind=T)[1],],max(ODD@data[,grep('hazMean', names(ODD@data))], na.rm=T))
  names(loc_mmi) <- c('lon', 'lat', 'max_mmi')
  map_data_df[nrow(map_data_df)+1,] = loc_mmi
  #iso3 <- c(iso3, unique(ODD$ISO3C))
  if (length(ODD@impact$impact)>0){
    nHazReal <- c(nHazReal, length(grep('hazMean', names(ODD@data))))
    maxIntReal <- c(maxIntReal, max(ODD@data[, grep('hazMean', colnames(ODD@data))],  na.rm=T))
    for (j in 1:NROW(ODD@impact)){
      df_Impact %<>% add_row(observed=ODD@impact$observed[j], impact=ODD@impact$impact[j], polygon=ODD@impact$polygon[j],
                             exposure=ifelse(impact=='buildDam', sum(ODD@data[ODD@polygons[[ODD@impact$polygon[j]]]$indexes,'nBuildings']), sum(ODD@data[ODD@polygons[[ODD@impact$polygon[j]]]$indexes,'Population'])),
                              event=i, I_max=max(ODD@data[ODD@polygons[[ODD@impact$polygon[j]]]$indexes, grep('hazMean', colnames(ODD@data))],  na.rm=T))
      # EQFreq_min <- c(EQFreq_min, min(ODD$EQFreq[ODD@polygons[[ODD@impact$polygon[j]]]$indexes], na.rm=T))
      # EQFreq_max <- c(EQFreq_max, max(ODD$EQFreq[ODD@polygons[[ODD@impact$polygon[j]]]$indexes], na.rm=T))
    }
  }
}

ggplot() +
  borders("world", colour = "gray80", fill = "gray80") +
  geom_point(data = map_data_df, aes(x = lon, y = lat, color = max_mmi), 
             size = 2.5, alpha = 0.8) +
  scale_color_gradientn(colors = c("darkgreen", 'chartreuse3', "gold", 'orange', "red")) +
  labs(x = "Longitude", y = "Latitude", col = "Max MMI") +
  theme_minimal() +
  theme(axis.title = element_text(family = "Liberation Serif", size=12),  
  legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
  legend.title = element_text(family = "Liberation Serif", size=12))
#event_map.pdf, 3.8 x 8

ggplot(df_Impact %>% filter(impact=='mortality'), aes(x=I_max, y=observed)) + geom_point()

#Plots:
# plot of mortality (simulated) vs mortality (observed)
library(cowplot)
library(scales)
plot_true_vs_simulated_obsvals <- function(impact_type){
  p <- ggplot()  + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 3, 10, 30, 100, 300, 1000), labels = label_comma(), expand = expand_scale(mult = c(0,0.1)))  + 
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000), labels=function(x) {ifelse(x==0, "0", ifelse(x==10, "10",parse(text=gsub("[+]", "", gsub("1e", "10^", scientific_format()(x))))))}, expand = expand_scale(mult = c(0,0.1))) +
    geom_histogram(data=df_SimImpact %>% filter(impact==impact_type), aes(x=observed,y=after_stat(count)), alpha=0.4, col="blue", lwd=0.2, fill="blue") +
    geom_histogram(data=df_Impact %>% filter(impact==impact_type), aes(x=observed,y=after_stat(count)), alpha=0.4, col='red', lwd=0.2, fill='red') +
    theme_bw() + ylab('Count') + xlab('Oberved') +
    theme(axis.title = element_text(family = "Liberation Serif", size=12),  
          legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
          legend.title = element_text(family = "Liberation Serif", size=12),
          panel.grid.minor = element_blank())
  return(p)
}
plot_true_vs_simulated_obscount <- function(impact_type){
  df_sim_tally <- df_SimImpact %>% group_by(event) %>% summarise(n = sum(impact == impact_type))
  df_real_tally <- df_Impact %>% group_by(event) %>% summarise(n = sum(impact == impact_type)) 
  if (impact_type=='mortality' | impact_type=='buildDam'){ 
    max_val <- max(max(df_sim_tally$n), max(df_real_tally$n))
    #breaks <- seq(2.5, max_val + 5 + (5-max_val %%5), 5)
    
    p <- ggplot() +
      scale_x_continuous(expand = expand_scale(add = c(2.5,7.5))) + 
      scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 3, 10, 30, 100, 300), labels = label_comma(), expand = expand_scale(mult = c(0,0.1)))  + 
      geom_histogram(data=df_sim_tally, aes(x=n, y=after_stat(count),fill='Simulated Data', col='Simulated Data'), binwidth=4, center=3,alpha=0.3, col='blue', lwd=0.2) + 
      geom_histogram(data=df_real_tally, aes(x=n, y=after_stat(count),fill='Real Data', col='Real Data'), binwidth=4, center=3, alpha=0.3, col='red', lwd=0.2) +
      ylab('Count') + xlab('Number of Observations per event') + theme_bw() +
      scale_fill_manual(values = c("Real Data" = "red", "Simulated Data" = "blue"))+ labs(fill = "") +
      guides(fill = guide_legend(override.aes = list(color = NULL))) +
      theme(axis.title = element_text(family = "Liberation Serif", size=12),  
            legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
            legend.title = element_text(family = "Liberation Serif", size=12),
            panel.grid.minor = element_blank())
    p
  } else {
    p <- ggplot() +
      scale_x_continuous(expand = expand_scale(mult = c(0,0.1))) + 
      scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 3, 10, 30, 100, 300), labels = label_comma(), expand = expand_scale(mult = c(0,0.1)))  + 
      geom_histogram(data=df_sim_tally, aes(x=n, y=after_stat(count),fill='Simulated Data', col='Simulated Data'), binwidth=1, center=0,alpha=0.3, col='blue', lwd=0.2) + 
      geom_histogram(data=df_real_tally, aes(x=n, y=after_stat(count),fill='Real Data', col='Real Data'), binwidth=1, center=0, alpha=0.3, col='red', lwd=0.2) +
      ylab('Count') + xlab('Number of Observations per event') + theme_bw() +
      scale_fill_manual(values = c("Real Data" = "red", "Simulated Data" = "blue"))+ labs(fill = "") +
      guides(fill = guide_legend(override.aes = list(color = NULL))) +
      theme(axis.title = element_text(family = "Liberation Serif", size=12),  
            legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
            legend.title = element_text(family = "Liberation Serif", size=12),
            panel.grid.minor = element_blank())
  }
  p
  return(p)
}
p_mort_obsvals <- plot_true_vs_simulated_obsvals('mortality') + xlab('Observed Mortality') 
p_disp_obsvals <- plot_true_vs_simulated_obsvals('displacement') + xlab('Observed Displacement')
p_bd_obsvals <- plot_true_vs_simulated_obsvals('buildDam') + xlab('Observed Building Damage')
p_mort_obscount <- plot_true_vs_simulated_obscount('mortality') + xlab('Number of Mortality Observations in Event') + scale_x_continuous(breaks=c(1,40,80, 120, 160), labels=c(1,40,80,120,160),expand = expand_scale(add=c(1, 0), mult = c(0,0.1)))
p_disp_obscount <- plot_true_vs_simulated_obscount('displacement') + xlab('Number of Displacement Observations in Event') + scale_x_continuous(breaks=0:7, labels=0:7, expand = expand_scale(mult = c(0,0.1)))
p_bd_obscount <- plot_true_vs_simulated_obscount('buildDam') + xlab('Number of Building Damage Observations in Event')
legend <- get_plot_component(p_mort_obscount +
                               theme(legend.position="bottom"), 'guide-box', return_all=T)[[3]]

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p_1 <- plot_grid( plot_grid( p_mort_obsvals, p_disp_obsvals, p_bd_obsvals, align = 'vh', nrow = 1),legend, ncol = 1, rel_heights=c(1,0.1))

p_2 <- plot_grid( plot_grid( p_mort_obscount + theme(legend.position="none"),
                             p_disp_obscount + theme(legend.position="none"), 
                             p_bd_obscount + theme(legend.position="none"), align = 'vh', nrow = 1), legend,ncol = 1, rel_heights=c(1,0.1))

plot_grid(p_1, p_2, ncol=1)

plot_grid(legend,plot_grid(p_mort_obsvals, p_disp_obsvals, p_bd_obsvals,
          p_mort_obscount + theme(legend.position="none"),
          p_disp_obscount + theme(legend.position="none"), 
          p_bd_obscount + theme(legend.position="none"), align = 'vh'), ncol = 1, rel_heights=c(0.075,1))

#SimVsObs.pdf, 7 x 12 inches


files <-  paste0(dir, "IIDIPUS_Input_Alternatives/IIDIPUS_SimInput/ODDobjects/")

ufiles<-na.omit(list.files(path=files,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
for (file in ufiles){
  print(file)
  ODD <- readRDS(paste0(files,file))
  print(NROW(ODD@impact))
}

#compare the number of hazards per event:
hist(nHazSim)
hist(nHazReal)
#compare the number of hazards per event:
hist(maxIntSim)
hist(maxIntReal)

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

merged_df <- merge(df_Impact%>% filter(impact==impact_type1), df_Impact %>% filter(impact==impact_type1), by=c('event'))
merged_df %<>% filter(polygon.x < polygon.y)
cor_true <- ggplot(merged_df, aes(x=observed.x, y=observed.y)) + xlab(impact_type1) + ylab(impact_type1) +
  geom_point() + scale_x_log10() + scale_y_log10() + ggtitle('Observed Impact Data, with each point representing an observation for one event and region')


merged_df_sim <- merge(df_SimImpact%>% filter(impact==impact_type1), df_SimImpact %>% filter(impact==impact_type1), by=c('event'))
merged_df_sim %<>% filter(polygon.x < polygon.y)
cor_sim <- ggplot(merged_df_sim, aes(x=observed.x, y=observed.y)) + xlab(impact_type1) + ylab(impact_type1) +
  geom_point() + scale_x_log10() + scale_y_log10() + ggtitle('Simulated Impact Data, with each point representing an observation for one event and region')

grid.arrange(cor_true, cor_sim, ncol=1)

# histogram of intensities (simulated) vs intensities (true)
impact_type='mortality'
ggplot() + 
  geom_histogram(data=df_Impact %>% filter(impact==impact_type), aes(x=I_max,y=after_stat(count)), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  geom_histogram(data=df_SimImpact %>% filter(impact==impact_type), aes(x=I_max,y=after_stat(count)), alpha=0.3, col='yellow', lwd=0.2, fill='yellow')  #+
  #scale_y_continuous(breaks = seq(0, 90, by = 10))

xx <- df_Impact %>% group_by(event) %>%
  summarise(unique_prop = n_distinct(polygon) / n())
yy <- df_SimImpact %>% group_by(event) %>%
  summarise(unique_prop = n_distinct(polygon) / n())

plot((df_Impact %>% group_by(event) %>% summarise(n_obs = n()))$n_obs)
points((df_SimImpact %>% group_by(event) %>% summarise(n_obs = n()))$n_obs, col='red')

ggplot() + 
  geom_histogram(data=df_Impact %>% filter(impact==impact_type), aes(x=exposure,y=after_stat(count)), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  geom_histogram(data=df_SimImpact %>% filter(impact==impact_type), aes(x=exposure,y=after_stat(count)), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') + 
  scale_x_log10() #+
#scale_y_continuous(breaks = seq(0, 90, by = 10))


#------------------------------------------------------------------------------------------------
#------------------------------------ CHECK RANKS/PAIRPLOTS -------------------------------------
#------------------------------------------------------------------------------------------------

#source('RCode/Model.R')
#source('RCode/Model Variations/ODDobjJune.R')

Omega=list(Lambda1 = list(nu=7.5, kappa=0.6),
           Lambda2 = list(nu=10.4, kappa=0.775), #list(nu=10.65, kappa=1.5), #
           Lambda3 = list(nu=7.9, kappa=0.68),
           #Lambda4 = list(nu=9.9, kappa=1.6),
           #theta= list(theta1=0.6),
           eps=list(local=0.4, hazard_mort=0.45, hazard_disp=0.65, hazard_bd=0.5, hazard_cor=0.65),
           #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
           vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15))

flattenImpactSample2 <- function(sampled_out){
  df <- data.frame(iso3 = sampled_out[[1]]$iso3,
                   polygon = sampled_out[[1]]$polygon,
                   impact = sampled_out[[1]]$impact,
                   observed = sampled_out[[1]]$observed)#,
  #n_pixels = impact_sample$poly[[1]]$n_pixels,
  #max_intensity = impact_sample$poly[[1]]$max_intensity)
  for (j in 1:length(sampled_out)){
    df %<>% cbind(sampled_out[[j]]$sampled)
  }
  colnames(df)[grep('sampled',names(df))] <- paste0('sampled.', 1:length(sampled_out))
  return(df)
}

i = 1
ufiles = na.omit(list.files(path=paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput6/ODDobjects/'),pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
ufiles = grep('^Train/' , ufiles, value = TRUE)
x <- file.info(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput6/ODDobjects/',ufiles))
ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
event_ids_all <- as.numeric(sub(".*_(\\d+)$", "\\1", ufiles))


par(mfrow=c(5,5), mar=c(0.2, 0.2, 0.2, 0.2))
for (i in event_ids_all[1:27]){
  print(i)
  file = ufiles[which(event_ids_all==i)]
  ODD <- readODD(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_SimInput7/ODDobjects/', file))
  
  Omega$eps$hazard_cor = 0.1
  sampled_out_single <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                         replace(which(names(AlgoParams)==c('Np')), 1), 
                       output='SampledAgg')
  
  sampled_out1 <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                                replace(which(names(AlgoParams)==c('Np')), 100), 
                              output='SampledAgg')
  
  Omega$eps$hazard_cor =0.9
  sampled_out <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                         replace(which(names(AlgoParams)==c('Np')), 100), 
                       output='SampledAgg')
  
  impact_mort1 = which(ODD@impact$impact=='buildDam')[1]
  impact_disp1 = which(ODD@impact$impact=='displacement')[1]
  if (is.na(impact_mort1) | is.na(impact_disp1)) next
  
  impact_sampled_flattened1 = flattenImpactSample2(sampled_out1)
  impact_sampled_flattened = flattenImpactSample2(sampled_out)
  
  combined_df <- cbind(
    c(log(sampled_out_single[[1]]$sampled[impact_mort1] + 10),
      log(as.numeric(impact_sampled_flattened[impact_mort1, 5:104]) + 10)),
    c(log(sampled_out_single[[1]]$sampled[impact_disp1] + 10),
      log(as.numeric(impact_sampled_flattened[impact_disp1, 5:104]) + 10))
  )
  print(get_banddepth_rank_single(t(combined_df)))
  
  # Compute MMDs
  mmds <- c()
  for (j in 1:101) {
    mmds <- c(mmds, mmds_sample(combined_df[j, ], t(combined_df[-j, ])))
  }
  
  # Normalize MMDs to a color scale
  colors <- colorRampPalette(c("blue", "red"))(101)
  color_indices <- as.numeric(cut(mmds, breaks = 101))
  point_colors <- colors[color_indices]
  
  # Plot base
  plot(
    log(as.numeric(impact_sampled_flattened[impact_mort1, 5:104]) + 10),
    log(as.numeric(impact_sampled_flattened[impact_disp1, 5:104]) + 10),
    xlab = '', ylab = '', main = file,
    col = point_colors[2:length(point_colors)], pch = 16
  )
  
  points(
    log(as.numeric(impact_sampled_flattened1[impact_mort1, 5:104]) + 10),
    log(as.numeric(impact_sampled_flattened1[impact_disp1, 5:104]) + 10),
    xlab = '', ylab = '', main = file,
    col = 'green', pch = 16
  )
  
  # Highlight first point
  points(
    log(sampled_out_single[[1]]$sampled[impact_mort1] + 10),
    log(sampled_out_single[[1]]$sampled[impact_disp1] + 10),
    pch = 1, cex = 2, col=point_colors[1]
  )
}



Omega <- Omega_true <-  list(Lambda1 = list(nu=7.5, kappa=0.6),
                             Lambda2 = list(nu=11.4, kappa=0.775), #list(nu=10.65, kappa=1.5), #
                             Lambda3 = list(nu=7.9, kappa=0.68),
                             #Lambda4 = list(nu=9.9, kappa=1.6),
                             #theta= list(theta1=0.6),
                             eps=list(local=1.2, hazard_mort=0.4, hazard_disp=0.65, hazard_bd=0.5, hazard_cor=0.1),
                             #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                             vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15))

testing <- function(){
  ufiles = na.omit(list.files(path=paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput6/ODDobjects/'),pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  ufiles = grep('^Train/' , ufiles, value = TRUE)
  x <- file.info(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput6/ODDobjects/',ufiles))
  ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
  event_ids_all <- as.numeric(sub(".*_(\\d+)$", "\\1", ufiles))
  for (file in ufiles){
    ODD <- readODD(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput6/ODDobjects/', file))
    Omega$eps$hazard_cor = -0.4
    sampled_out_single <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                                  replace(which(names(AlgoParams)==c('Np')), 2), 
                                output='SampledAgg')
    
    ODD@impact$id <- 1:nrow(ODD@impact)
    
    #flats = flattenImpactSample2(sampled_out_single)
    #plot(log(t(flats[1,5:104])+10), log(t(flats[2,5:104])+10))
    
    merged_df <- merge(
      ODD@impact,
      sampled_out_single[[1]],
      by = c('polygon', 'impact', 'iso3', 'observed', 'inferred', 'qualifier', 'sdate')
    )
  
    merged_df <- merged_df[order(merged_df$id), ]
    
    ODD@impact$observed = merged_df$sampled
    
    #ODD@impact$observed = sampled_out_single[[1]]$sampled
    saveODD(ODD, paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput7/ODDobjects/', file))
  }
  
  Omega$eps$hazard_cor = -0.4
  AlgoParams$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput7/'
  AlgoParams$Np = 1
  AlgoParams$m_CRPS = 100
  AlgoParams$cores =2
  AlgoParams$NestedCores =1
  impact_sample = SampleImpact(dir, Model, Omega %>% addTransfParams, AlgoParams, dat='Train', output='SampledAgg')
  CalcDist(impact_sample, AlgoParams)
  
  flattened = flattenImpactSample(impact_sample)
  for (i in event_ids_all[1:27]){
    combined_df = flattened %>% filter(event_id == i)
    
    impact_mort1 = which(combined_df$impact=='buildDam')[1]
    impact_disp1 = which(combined_df$impact=='displacement')[1]
    if (is.na(impact_mort1) | is.na(impact_disp1)) next
    
    combined_df = combined_df[c(impact_mort1, impact_disp1), 7:NCOL(combined_df)]
    print(get_banddepth_rank_single(t(combined_df)))
    
   plot(log(t(combined_df[1,])+10), log(t(combined_df[2,])+10))
   points(log(combined_df[1,1]+10), log(combined_df[2,1]+10),col='red', pch=19)
  }
  
  
  
  impact_sample = SampleImpact(dir, Model, Omega %>% addTransfParams, AlgoParams, dat='Train', output='SampledAgg')
  CalcDist(impact_sample, AlgoParams)
  
  impact_sample = SampleImpact(dir, Model, Omega %>% addTransfParams, AlgoParams, dat='Train', output='SampledAgg')
  CalcDist(impact_sample, AlgoParams)
  
  
  
}



#### MOST RIGOROUS AS OF 21 JULY #######

Omega <- Omega_true <-  list(Lambda1 = list(nu=7.5, kappa=0.6),
                             Lambda2 = list(nu=11.4, kappa=0.775), #list(nu=10.65, kappa=1.5), #
                             Lambda3 = list(nu=7.9, kappa=0.68),
                             #Lambda4 = list(nu=9.9, kappa=1.6),
                             #theta= list(theta1=0.6),
                             eps=list(local=1.2, hazard_mort=0.4, hazard_disp=0.65, hazard_bd=0.5, hazard_cor=0.1),
                             #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                             vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15))

source(paste0(dir, 'RCode/Unfinished/Model_allRanks.R'))
testing2 <- function(){
  #Omega$eps$local = 0.4
  hazard_cors <- c(0.05, 0.5, 0.95)
  n_dim <- 9   # dimensions of CalcDist output
  n_rep <- 20
  n_hz <- length(hazard_cors)
  
  # Initialize 4D array: [dim, repeat, test, gen]
  results_array <- array(NA, dim = c(n_dim, n_rep, n_hz, n_hz),
                         dimnames = list(NULL,
                                         paste0("rep", 1:n_rep),
                                         paste0("test_", hazard_cors),
                                         paste0("gen_", hazard_cors)))
  
  # Load and order training files
  ufiles <- na.omit(list.files(path = paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9/ODDobjects/'),
                               pattern = Model$haz, recursive = TRUE, ignore.case = TRUE))
  ufiles <- grep('^Train/', ufiles, value = TRUE)
  x <- file.info(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9/ODDobjects/', ufiles))
  ufiles <- na.omit(ufiles[match(length(ufiles):1, rank(x$size))])
  event_ids_all <- as.numeric(sub(".*_(\\d+)$", "\\1", ufiles))
  
  
  for (gen_idx in seq_along(hazard_cors)) {
    gen_hz <- hazard_cors[gen_idx]
    cat("Generating data for hazard_cor =", gen_hz, "\n")
    
    for (rep in 1:n_rep){
      
      Omega$eps$hazard_cor <- gen_hz
      
      for (file in ufiles) {
        ODD <- readODD(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9/ODDobjects/', file))
        #ODD$Population = ifelse(values(ODD$Population) >7500 & values(ODD$hazMean1)> 8, 15000, 0)
        #ODD$Population = values(ODD$Population) * exp(values(ODD$Population)/max(values(ODD$Population)))/exp(1)
        
        # Omega$eps$hazard_cor = 0.1
        AlgoParams$cores <- 2
        AlgoParams$NestedCores <- 1
        suppressMessages(suppressWarnings(invisible(capture.output({result = DispX(ODD,
                                                                                   Omega %>% addTransfParams(),
                                                                                   Model$center,
                                                                                   AlgoParams %>% replace(which(names(AlgoParams) == 'm_CRPS'), 1) %>%
                                                                                     replace(which(names(AlgoParams) == 'Np'), 2),
                                                                                   output = 'SampledAgg')}))))
        
        sampled_out_single = result
        
        # Omega$eps$hazard_cor = 0.9
        # suppressMessages(suppressWarnings(invisible(capture.output({result = DispX(ODD,
        #                                                                            Omega %>% addTransfParams(),
        #                                                                            Model$center,
        #                                                                            AlgoParams %>% replace(which(names(AlgoParams) == 'm_CRPS'), 1) %>%
        #                                                                              replace(which(names(AlgoParams) == 'Np'), 100),
        #                                                                            output = 'SampledAgg')}))))
        # 
        # sampled_out_single2 = result
        
        
        
        ODD@impact$id <- 1:nrow(ODD@impact)
        
        # flattened1 = flattenImpactSample2(sampled_out_single)
        # flattened2 = flattenImpactSample2(sampled_out_single2)
        # pairs(log(rbind(t(flattened1[,5:104]), t(flattened2[,5:104]))+10), col=c(rep('red',100), rep('blue',100)))
        # 
        
        merged_df <- merge(
          ODD@impact,
          sampled_out_single[[1]],
          by = c('polygon', 'impact', 'iso3', 'observed', 'inferred', 'qualifier', 'sdate')
        )
        
        merged_df <- merged_df[order(merged_df$id), ]
        
        ODD@impact$observed = merged_df$sampled
        
        saveODD(ODD, paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput10/ODDobjects/', file))
      }
      
      for (test_idx in seq_along(hazard_cors)) {
        test_hz <- hazard_cors[test_idx]
        cat("Testing with hazard_cor =", test_hz, "on gen =", gen_hz, "\n")
        
          Omega$eps$hazard_cor <- test_hz
          AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput10/'
          AlgoParams$Np <- 1
          AlgoParams$m_CRPS <- 100
          AlgoParams$cores <- 2
          AlgoParams$NestedCores <- 1
          suppressMessages(suppressWarnings(invisible(capture.output({impact_sample = SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams, dat = 'Train', output = 'SampledAgg')}))))
          dist_output <- CalcDist(impact_sample, AlgoParams)
          
          print(dist_output)
          
          # Store in 4D array
          results_array[, rep, test_idx, gen_idx] <- as.numeric(dist_output)
          saveRDS(results_array, file = paste0("CalcDist_4D_matrix", Sys.Date(), ".rds"))
        }
    }
  }
}


# results_array <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/CalcDist_4D_matrix2025-07-23.rds')
par(mfrow=c(2,2), mar = c(4, 4, 2, 1), family="Times")
method_colors <- c('#440154', '#FDE725FF', '#1FA187')
method_labels <- c('Method 1', 'Method 2', 'Method 3')

titles <- c('Energy Score', '(a) Multivariate Rank', '(b) MST Rank','MVR MortDisp', '(c) BD Rank', 'BD MortDisp', '(d) Average Rank', 'Ave Rank MortDisp', 'MST rank mortdisp')

for (i in c(2,3,5, 7)) {
  box_data <- list()
  box_names <- list()
  for (j in 1:3) {
    for (k in 1:3) {
      box_data[[length(box_data)+1]] = results_array[i, ,j,k] / 0.05
      box_names[[length(box_names) + 1]] <- substitute(rho[Sim] == x, list(x = hazard_cors[k]))
    }
  }
  
  group_labels <- lapply(hazard_cors, function(x) substitute(rho[Obs] == y, list(y = x)))
  
  ylims <- range(results_array[i,,,]/0.05, na.rm = TRUE) 
  
  # Add 10% vertical space on top
  y_margin <- 0.1 * diff(ylims)
  ylims[2] <- ylims[2] + y_margin
  
  bp <- boxplot(box_data, 
                col = rep(method_colors, times = 3),
                xaxt = "n",
                main = "", #titles[i],
                ylim = ylims,
                ylab = "",
                cex.axis = 1)
  mtext(titles[i], side = 3, line = 0.5, adj = 0, cex = 0.9, font = 1)
  
  # Tick labels for each box
  tick_labels <- rep(hazard_cors, times = 3)
  axis(1, at = 1:9, labels = tick_labels, las = 1, cex.axis = 1)
  
  # Add x-axis label
  mtext(expression(rho[Sim]), side = 1, line = 2.5, cex = 1)
  
  # Add rho[Obs] inside plot above groups
  group_positions <- c(mean(1:3), mean(4:6), mean(7:9))
  usr <- par("usr")
  y_pos <- usr[4] - 0.07 * diff(usr[3:4])  # inside plot, near top
  
  text(x = group_positions, y = y_pos, labels = do.call(expression, group_labels), cex = 1.25)
  
  # Optional: Add group separators
  abline(v = c(3.5, 6.5), lty = 2, col = "grey")
}
#8.5 x 6, PreRankComparison.pdf


testing3 <- function(){
  hazard_cors <- c(0.1, 0.5, 0.9)
  n_dim <- 8   # dimensions of CalcDist output
  n_rep <- 6
  n_hz <- length(hazard_cors)
  
  # Initialize 4D array: [dim, repeat, test, gen]
  results_array <- array(NA, dim = c(n_dim, n_rep, n_hz),
                         dimnames = list(NULL,
                                         paste0("rep", 1:n_rep),
                                         paste0("test_", hazard_cors)))
    
  for (rep in 1:n_rep){
    for (test_idx in seq_along(hazard_cors)) {
      test_hz <- hazard_cors[test_idx]
      cat("Testing with hazard_cor =", test_hz, "on gen =", 0.1, "\n")
      
      Omega$eps$hazard_cor <- test_hz
      AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput8/'
      AlgoParams$Np <- 1
      AlgoParams$m_CRPS <- 100
      AlgoParams$cores <- 1
      AlgoParams$NestedCores <- 1
      suppressMessages(suppressWarnings(invisible(capture.output({ impact_sample = SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams, dat = 'Train', output = 'SampledAgg')}))))
      dist_output <- CalcDist(impact_sample, AlgoParams)
      
      print(dist_output)
      
      # Store in 4D array
      results_array[, rep, test_idx] <- as.numeric(dist_output)
      saveRDS(results_array, file = paste0("CalcDist_4D_matrix_fixedGenTrain", Sys.Date(), ".rds"))
    }
  }
}

#results_array = readRDS('CalcDist_4D_matrix.rds')
titles = c('ES', 'MVR Rank', 'MMDS', 'MVR Rank MortDisp', 'BD Rank', 'BD MortDisp', 'MMDs Rank', 'MMDS Ranks MortDisp')

par(mfrow=c(4,3))
for (j in c(2,5,6,7)){
  for (gen_idx in 1:3){
    box_data <- as.list(rep(NA, length(hazard_cors)))
    for (test_idx in seq_along(hazard_cors)) {
      box_data[[test_idx]] <- results_array[j, , test_idx, gen_idx]
    }
    
    boxplot(box_data, names = round(hazard_cors, 2), main = paste0('Hazard gen: ', hazard_cors[gen_idx]),
            xlab = "Hazard Test", ylab = titles[j])
    
    #plot(rep(hazard_cors,each=20), c(results_array[j,, ,gen_idx]), main=hazard_cors[gen_idx])
    #abline(v=hazard_cors[gen_idx])
  }
}

par(mfrow=c(3,2))
for (gen_idx in 1:6){
  plot(1:n_rep, c(results_array[2,, 1,gen_idx]), ylim=range(results_array[2,,,gen_idx]), col='yellow')
  #points(1:n_rep, c(results_array[2,,2,gen_idx]), col='orange')
  points(1:n_rep, c(results_array[2,,3,gen_idx]), col='red')
  points(1:n_rep, c(results_array[2,,4,gen_idx]), col='darkred')
  points(1:n_rep, c(results_array[2,,5,gen_idx]), col='purple')
  points(1:n_rep, c(results_array[2,,6,gen_idx]), col='blue')
  abline(v=hazard_cors[gen_idx])
}

par(mfrow=c(4,3))
for(i in 1:10){
  for (gen_idx in 1:1){
    plot(rep(hazard_cors,each=10), c(results_array[i,, ]), main=hazard_cors[gen_idx])
    abline(v=hazard_cors[gen_idx])
  }
}

titles = c('ES', 'MVR Rank', 'MMDS', 'MVR Rank MortDisp', 'BD Rank', 'BD MortDisp', 'MMDs Rank', 'MMDS Ranks MortDisp')

par(mfrow = c(3, 4))
for (i in 1:10) {
  box_data <- as.list(rep(NA, length(hazard_cors)))
  for (gen_idx in seq_along(hazard_cors)) {
    box_data[[gen_idx]] <- results_array[i, , gen_idx]
  }
  
  boxplot(box_data, names = round(hazard_cors, 2), main = titles[i],
          xlab = "hazard correlation", ylab = "Value")
}


#------------------------------------------------------------------------------------------------
#------------------------------------ CHECK HIGH LEVEL PRIORS -----------------------------------
#------------------------------------------------------------------------------------------------

n_samples <- 300
samples <- array(NA, dim=c(n_samples, length(unlist(Omega))))
for (i in 1:n_samples){
  if (i %% 10 == 0){
    print(paste('Sample:',i))
  }
  samples[i,] <- Proposed2Physical(HLPrior_sample(Model, AlgoParams), Model) %>% unlist()
}

par(mfrow=c(3, 2))
pairings <- rbind(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12))
for (p in 1:NROW(pairings)){
  plot(samples[,pairings[p,1]], samples[,pairings[p,2]], 
       xlab=names(unlist(Model$skeleton))[pairings[p,1]],
       ylab=names(unlist(Model$skeleton))[pairings[p,2]])
  if (!is.null(Omega)){
    points(unlist(Omega)[pairings[p,1]], unlist(Omega)[pairings[p,2]], col='red', pch=4, cex=2, lwd=4)
  }
}
pairings <- rbind(c(13,14), c(15,16), c(17,18), c(19,20), c(21,22))
for (p in 1:NROW(pairings)){
  plot(samples[,pairings[p,1]], samples[,pairings[p,2]], 
       xlab=names(unlist(Model$skeleton))[pairings[p,1]],
       ylab=names(unlist(Model$skeleton))[pairings[p,2]])
  if (!is.null(Omega)){
    points(unlist(Omega)[pairings[p,1]], unlist(Omega)[pairings[p,2]], col='red', pch=4, cex=2, lwd=4)
  }
}
par(mfrow=c(1, 1))

Omega_test <- Omega
Omega_test$vuln_coeff$SHDI <- -0.15
Omega_test$Lambda1$kappa <- 1.9
Omega_test$eps$local <- 1.2
Model$HighLevelPriors(Omega_test %>% addTransfParams(), Model)


#------------------------------------------------------------------------------------------------
#------------------------------------ RESULTS ANALYSIS ------------------------------------------
#------------------------------------------------------------------------------------------------

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-05_202654_alpha0.9_simulated1D_v3_crps_noplus10_200events_M10')
# Notes: Take everything with a big grain of salt as only 50 particles
#           - essstore freezes from steps 48 - 65ish, something odd going on
#           - acceptance rate stays within 0.1 - 0.45, doesn't seem to be dropping
#           - simulated data for some events is unexpectedly 0 due to a bug!

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-08_185840_alpha0.9_simulatedfull_energyscore_150events_Npart50_M60')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-04-09_135858_alpha0.9_M60_Npart1000_150events_simulatedfull')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-11_114910_alpha0.9_M60_Npart1000_150events_simulatedfull')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-15_100950_alpha0.9_M60_Npart1000_150events_simulatedfull')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-16_172900_alpha0.9_M60_Npart1000_150events_simulatedfull_smallereps')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-19_093811_alpha0.9_M60_Npart1000_150events_simulatedfull_smallereps')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-22_214949_alpha0.9_M60_Npart1000_150events_simulatedfull_smallereps')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-23_130232_alpha0.9_M60_Npart1000_150events_simulatedfull_smallereps')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-23_130232_alpha0.9_M60_Npart1000_150events_simulatedfull_smallereps')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-25_093630_alpha0.9_M60_Npart1000_150events_simulatedfull_smallereps')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-28_175817_alpha0.9_M60_Npart1000_150events_simulatedfull_smallereps')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-28_175817_alpha0.9_M60_Npart1000_150events_simulatedfull_smallereps')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-04-30_174728_alpha0.9_M60_Npart1000_150events_simulatedfull_15by15')
#     - have now reduced to a 15x 15 grid

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-05-03_172241_alpha0.9_M60_Npart1000_150events_simulatedfull_15by15')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-05-05_155223_alpha0.9_M60_Npart1000_150events_simulatedfull_15by15')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-05-16_222830_alpha0.9_M60_Npart1000_150events_simulatedfull_15by15')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-05-17_140645_alpha0.9_M60_Npart1000_150events_simulatedfull_15by15')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-05-31_140728_alpha0.9_M60_Npart1000_150events_simulatedfull_15by15_wRankScores')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-05-31_162704_alpha0.9_M60_Npart1000_150events_simulatedfull_15by15_wRankScores')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-06-05_120111_alpha0.9_M60_Npart990_150events_simulatedfull_15by15_wRankScores')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-06-12_163924_alpha0.9_M60_Npart990_150events_simulatedfull_15by15_wBothRankScores')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-12_185056_alpha0.9_M60_Npart990RealAgg3')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-12_163924_alpha0.9_M60_Npart990_150events_simulatedfull_15by15_wBothRankScores')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-06-13_124500_alpha0.9_M60_Npart990RealAgg3')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-15_194119_alpha0.95_M60_Npart990_150events_simulatedfull_15by15_wBothRankScores_alphaincreased')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-15_194119_alpha0.95_M60_Npart990_150events_simulatedfull_15by15_wBothRankScores_alphaincreased')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-16_083342_alpha0.9_M60_Npart990RealAgg3')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-16_083342_alpha0.9_M60_Npart990RealAgg3')

plot_correlated_posteriors(AlgoResults, Omega=Omega, pairings=rbind(c(1,2), c(3,4), c(5,6),c(9,10), c(11,12), c(13,14)))
plot_correlated_posteriors(AlgoResults, Omega=Omega, pairings=rbind(c(15,16), c(17,18), c(19,20), c(21,22), c(7,8), c(22,23)))

plot_corr_posterior_vs_d(AlgoResults, Omega=Omega, pairing=c(21,22))


plot_corr_transf_posterior_vs_d(AlgoResults, Omega=Omega, pairing=c(41,44))

plot_posteriors(AlgoResults, 10:15)

plot(AlgoResults$Omega_sample_phys[,2,30], AlgoResults$Omega_sample_phys[,3,30])

unique_part <- c()
for (s in 1:AlgoResults$s_finish){
  unique_part <- c(unique_part, length(unique(AlgoResults$Omega_sample_phys[,1,s])))
}
plot(unique_part)

plot(AlgoResults$Omega_sample_phys[,3,195], AlgoResults$d_full[,1,3,195])

plot_acc_prob(AlgoResults); abline(h=0.05)
plot_d_vs_step(AlgoResults, 5); abline(h=2.25, col='red')
plot(AlgoResults$essstore)

impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams)
CalcDist(impact_sample, AlgoParams)

# compare different particles: 
cor_seq <- c(0.1, 0.25, 0.4, 0.55, 0.7, 0.95)
n_repeats <- 4
results <- array(0, dim=c(n_repeats, length(cor_seq), 2))
for (i in 1:n_repeats){
  # Omega_high_cor <- AlgoResults$Omega_sample_phys[which.max(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  # Omega_low_cor <- AlgoResults$Omega_sample_phys[which.min(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  # Omega_true_cor <- AlgoResults$Omega_sample_phys[which.min(abs(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]-0.55)),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  # 
  for (j in 1:length(cor_seq)){
    om_cor <- cor_seq[j]
    Omega <- Omega_true
    Omega$eps$hazard_cor <- om_cor
    impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams)
    dist <- CalcDist(impact_sample, AlgoParams)
    results[i,j,1] <- dist[4]
    results[i,j,2] <- dist[5]
  }
}

cor_seq <- c(0.1, 0.55,0.95)
n_repeats <- 20
results2 <- array(0, dim=c(n_repeats, length(cor_seq), 2))
for (i in 1:n_repeats){
  # Omega_high_cor <- AlgoResults$Omega_sample_phys[which.max(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  # Omega_low_cor <- AlgoResults$Omega_sample_phys[which.min(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  # Omega_true_cor <- AlgoResults$Omega_sample_phys[which.min(abs(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]-0.55)),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  # 
  for (j in 1:length(cor_seq)){
    om_cor <- cor_seq[j]
    Omega <- Omega_true
    Omega$eps$hazard_cor <- om_cor
    impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams)
    dist <- CalcDist(impact_sample, AlgoParams)
    print(dist)
    results2[i,j,1] <- dist[4]
    results2[i,j,2] <- dist[5]
  }
}
plot(rep(cor_seq, each=20), results2[,,1], xlab='Correlation', ylab='VS Distance')
#========================================

folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")

ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
ufiles <- grep('^Train/' , ufiles, value = TRUE)
AlgoParams$Np <- AlgoParams$Np * AlgoParams$m_CRPS
Omega_highcor <- Omega_true
Omega_highcor$eps$hazard_cor <- 0.99
for (i in 1:length(ufiles)){
  ODDy<-readRDS(paste0(folderin,ufiles[i]))
  tLL_truecor <- DispX(ODD = ODDy,Omega = Omega_true %>% addTransfParams(),center = Model$center, Method = AlgoParams, output='SampledAgg')
  tLL_highcor <- DispX(ODD = ODDy,Omega = Omega_highcor %>% addTransfParams(),center = Model$center, Method = AlgoParams, output='SampledAgg')
  obs <- tLL_truecor[[1]]$observed
  sampled_truecor <- matrix(unlist(lapply(tLL_truecor, function(x){x$sampled})), ncol=NROW(ODDy@impact), byrow=T)
  sampled_highcor <- matrix(unlist(lapply(tLL_highcor, function(x){x$sampled})), ncol=NROW(ODDy@impact), byrow=T)
  plot(sampled_truecor[,c(1,2)], xlab='Sampled Displacement', ylab='Sampled Mortality')
  points(sampled_highcor[,c(1,2)], col='red')
  points(obs[1], obs[2], col='blue', pch=19)
  tLL_truecor[[1]]
  w_vs <- matrix(unlist(AlgoParams$impact_weights[tLL_truecor[[1]]$impact]) %*% t(unlist(AlgoParams$impact_weights[tLL_truecor[[1]]$impact])), ncol=NROW(tLL_truecor[[1]]))
  vs_true <- vs_sample(log(obs+10), log(t(sampled_truecor)+10), w_vs=w_vs)
  vs_highcor <- vs_sample(log(obs+10), log(t(sampled_highcor)+10), w_vs=w_vs)
  
  vs_sample(log(obs+10),t(log(sampled_truecor+10)))
  vs_sample(log(obs+10), t(log(sampled_highcor+10)))
  
  vs_sample(obs,t(sampled_truecor))
  vs_sample(obs, t(sampled_highcor))
  
  print(paste(vs_true, vs_highcor))
  vs_true_unlog <- vs_sample(obs, t(sampled_truecor))
  vs_highcor_unlog <- vs_sample(obs, t(sampled_highcor))
  print(paste(vs_true_unlog, vs_highcor_unlog))
}

#compare conditional predictions: with vs without correlation term
AlgoParams$Np <- 1
preds_omega_corpost <- matrix(NA, nrow=0, ncol=2)
preds_omega_corprior <- matrix(NA, nrow=0, ncol=2)
for (i in 1:500){
  print(i)
  Omega_i <- sample(1:990, 1)
  Omega <- AlgoResults$Omega_sample_phys[Omega_i,,220] %>% relist(skeleton=Model$skeleton)
  
  #ODDyAgg <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Train/EQ20191215PHL_135')
  #event <- 'EQ20230206TUR_169' #'EQ20191215PHL_135'
  #ODDyAgg <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects_RealAgg3/Train/', event))
  ODDyAgg <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Train/EQ20130409ABC_-3')
  AggImpact <- DispX(ODDyAgg, Omega %>% addTransfParams(), Model$center, AlgoParams, output='SampledAgg')
  impAgg_sampled <- AggImpact[[1]]$sampled[c(1,2)]
  preds_omega_corpost <- rbind(preds_omega_corpost, impAgg_sampled)
  
  Omega$eps$hazard_cor <- runif(1,0,1)
  AggImpact <- DispX(ODDyAgg, Omega %>% addTransfParams(), Model$center, AlgoParams, output='SampledAgg')
  impAgg_sampled <- AggImpact[[1]]$sampled[c(1,2)]
  preds_omega_corprior <- rbind(preds_omega_corprior, impAgg_sampled)
}



df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=T, particle_best = F)

fitpred_mort <- plot_df_postpredictive(df_postpredictive_sampled_best, 'mortality')
fitpred_mort
fitpred_disp <- plot_df_postpredictive(df_postpredictive_sampled_best, 'displacement')
fitpred_disp
fitpred_bd <-plot_df_postpredictive(df_postpredictive_sampled_best, 'buildDam')
fitpred_bd

df_postpredictive_true <- create_df_postpredictive(AlgoResults, single_particle=T, Omega = Omega_true)

truepred_mort <- plot_df_postpredictive(df_postpredictive_true, 'mortality')
truepred_disp <- plot_df_postpredictive(df_postpredictive_true, 'displacement')
truepred_bd <- plot_df_postpredictive(df_postpredictive_true, 'buildDam')

grid.arrange(fitpred_mort, truepred_mort, nrow=2)

Omega_adj <- Omega
Omega_adj$vuln_coeff$PDens <- 0.8
df_postpredictive_adj <- create_df_postpredictive(AlgoResults, single_particle=T, Omega = Omega_adj)

adjpred_mort <- plot_df_postpredictive(df_postpredictive_adj, 'mortality')
adjpred_disp <- plot_df_postpredictive(df_postpredictive_adj, 'displacement')
adjpred_bd <- plot_df_postpredictive(df_postpredictive_adj, 'buildDam')

grid.arrange(truepred_mort, adjpred_mort, nrow=2)

#compare true and best distance:
impact_sample_trueOmega <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(),AlgoParams)
CalcDist(impact_sample_trueOmega, AlgoParams)
min(AlgoResults$d, na.rm=T)

eucldist2truepart <- function(part,Omega_true){
  sqrt(sum((part[c(1:6, 9:23)]-Omega_true[c(1:6, 9:23)])^2))
}

eucl_dists <- array(NA, dim=c(AlgoResults$s_finish, length(AlgoResults$Omega_sample_phys[,1,20])))
for (s in 1:AlgoResults$s_finish){
  eucl_dists[s,] <-apply(AlgoResults$Omega_sample_phys[,,s], 1, eucldist2truepart, unlist(Omega_true))
}
plot(rep(1:AlgoResults$s_finish, 1000), eucl_dists)

plot(AlgoResults$Omega_sample_phys[,15,150], AlgoResults$Omega_sample_phys[,9,150])

plot(apply(AlgoResults$Omega_sample_phys[,,AlgoResults$s_finish], 1, eucldist2truepart, unlist(Omega_true)))
points(apply(AlgoResults$Omega_sample_phys[,,1], 1, eucldist2truepart, unlist(Omega_true)), col='blue')
eucl_dist <- apply(AlgoResults$Omega_sample_phys[,,AlgoResults$s_finish], 1, eucldist2truepart, unlist(Omega_true))
plot(eucl_dist,AlgoResults$d[,1,AlgoResults$s_finish])
eucl_dist10 <- apply(AlgoResults$Omega_sample_phys[,,20], 1, eucldist2truepart, unlist(Omega_true))
points(eucl_dist10,AlgoResults$d[,1,20], col='blue')

#investigate broad spread of correlation covariate:
Omega_high_cor <- AlgoResults$Omega_sample_phys[which.max(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
Omega_low_cor <- AlgoResults$Omega_sample_phys[which.min(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
df_postpredictive_sampled_highcor <- create_df_postpredictive(AlgoResults, Omega=Omega_high_cor, single_particle=T)
df_postpredictive_sampled_lowcor <- create_df_postpredictive(AlgoResults, Omega=Omega_high_cor, single_particle=T)


Omegatrue_high_cor <- Omega_true
Omegatrue_high_cor$eps$hazard_cor <- 0.95
Omegatrue_low_cor <- Omega_true
Omegatrue_low_cor$eps$hazard_cor <- 0.05
df_postpredictive_true <- create_df_postpredictive(AlgoResults, single_particle=T, Omega = Omega_true, M=60)
df_postpredictive_truehighcor <- create_df_postpredictive(AlgoResults, single_particle=T, Omega = Omegatrue_high_cor, M=60)
df_postpredictive_truelowcor <- create_df_postpredictive(AlgoResults, single_particle=T, Omega = Omegatrue_low_cor, M=60)

i_1 <- 1; i_2 <- 4
plot(as.numeric(df_postpredictive_true[i_1,6:55]), as.numeric(df_postpredictive_true[ i_2,6:55]))
points(as.numeric(df_postpredictive_sampled_highcor[i_1,6:55]), as.numeric(df_postpredictive_sampled_highcor[ i_2,6:55]), col='green')
points(as.numeric(df_postpredictive_sampled_lowcor[i_1,6:55]), as.numeric(df_postpredictive_sampled_lowcor[ i_2,6:55]), col='blue')
points(df_postpredictive_true[i_1,5],df_postpredictive_true[ i_2,5], col='red' ,pch=12, cex=2)

i_1 <- 1; i_2 <- 3
plot(as.numeric(df_postpredictive_true[i_1,6:55]), as.numeric(df_postpredictive_true[ i_2,6:55]))
points(as.numeric(df_postpredictive_truehighcor[i_1,6:55]), as.numeric(df_postpredictive_truehighcor[ i_2,6:55]), col='green')
points(as.numeric(df_postpredictive_truelowcor[i_1,6:55]), as.numeric(df_postpredictive_truelowcor[ i_2,6:55]), col='blue')
points(df_postpredictive_true[i_1,5],df_postpredictive_true[ i_2,5], col='red' ,pch=12, cex=2)

points(as.numeric(df_postpredictive_truelowcor[i_1,6:155]), as.numeric(df_postpredictive_truelowcor[ i_2,6:155]), col='blue')
points(as.numeric(df_postpredictive_truehighcor[i_1,6:155]), as.numeric(df_postpredictive_truehighcor[ i_2,6:155]), col='green')
points(df_postpredictive_true[i_1,5],df_postpredictive_true[ i_2,5], col='red' ,pch=12, cex=2)
#particle_min.d <- which(min(AlgoResults$d, na.rm=T)==AlgoResults$d, arr.ind=T)

ii <- 2:3
plot(t(as.matrix(sweep(log(df_postpredictive_true[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_true$impact[ii]]), "*"))))
points(t(as.matrix(sweep(log(df_postpredictive_truelowcor[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_truelowcor$impact[ii]]), "*"))), col='blue')
points(t(as.matrix(sweep(log(df_postpredictive_truehighcor[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_truehighcor$impact[ii]]), "*"))), col='green')
points(t(as.numeric(log(df_postpredictive_true$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_true$impact[ii]]))), col='red', cex=2, pch=12)


i_1 <- 1153; i_2 <- 1154
ii <- c(i_1,i_2)

mean_mort <- c()
n_impacts <- c()
es_true <- c()
es_low <- c()
es_high <- c()
for (event_id in unique(df_postpredictive_true$event_id)){
  ii <- which(df_postpredictive_true$event_id==event_id)
  
  es_true <- c(es_true, es_sample(as.numeric(log(df_postpredictive_true$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_true$impact[ii]])), 
            as.matrix(sweep(log(df_postpredictive_true[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_true$impact[ii]]), "*"))))
  
  es_low <- c(es_low, es_sample(as.numeric(log(df_postpredictive_truelowcor$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_truelowcor$impact[ii]])), 
            as.matrix(sweep(log(df_postpredictive_truelowcor[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_truelowcor$impact[ii]]), "*"))))
  es_high <- c(es_high, es_sample(as.numeric(log(df_postpredictive_truehighcor$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_truehighcor$impact[ii]])), 
                                as.matrix(sweep(log(df_postpredictive_truehighcor[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_truehighcor$impact[ii]]), "*"))))
  
  n_impacts <- c(n_impacts, length(ii))
  mean_mort <- c(mean_mort, mean(df_postpredictive_true[ii[which(df_postpredictive_true$impact[ii]=='mortality')], 5]))
  
}

df_postpredictive_true[ii,]

i_1 <-241
i_2 <-248
plot(as.numeric(df_postpredictive_true[i_1,6:55]), as.numeric(df_postpredictive_true[ i_2,6:55]), 
     xlim=range(c(as.numeric(df_postpredictive_true[i_1,6:55]),as.numeric(df_postpredictive_truehighcor[i_1,6:55]),as.numeric(df_postpredictive_truelowcor[i_1,6:55]))),
     ylim=range(c(as.numeric(df_postpredictive_true[ i_2,6:55]),as.numeric(df_postpredictive_truehighcor[ i_2,6:55]),as.numeric(df_postpredictive_truelowcor[ i_2,6:55]))))
points(as.numeric(df_postpredictive_truehighcor[i_1,6:55]), as.numeric(df_postpredictive_truehighcor[ i_2,6:55]), col='green')
points(as.numeric(df_postpredictive_truelowcor[i_1,6:55]), as.numeric(df_postpredictive_truelowcor[ i_2,6:55]), col='blue')
points(df_postpredictive_true[i_1,5],df_postpredictive_true[ i_2,5], col='red' ,pch=12, cex=2)

es_sample(as.numeric(log(df_postpredictive_sampled_highcor$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_sampled_highcor$impact[ii]])), 
          as.matrix(sweep(log(df_postpredictive_sampled_highcor[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_sampled_highcor$impact[ii]]), "*")))

pred_true <- as.matrix(sweep(log(df_postpredictive_true[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_true$impact[ii]]), "*"))
pred_highcor <- as.matrix(sweep(log(df_postpredictive_sampled_highcor[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_sampled_highcor$impact[ii]]), "*"))
obs <- as.numeric(log(df_postpredictive_sampled_highcor$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_sampled_highcor$impact[ii]]))


flattenImpactSample <- function(impact_sample){
  df <- data.frame(event_id = impact_sample$poly[[1]]$event_id,
                   iso3 = impact_sample$poly[[1]]$iso3,
                   polygon = impact_sample$poly[[1]]$polygon,
                   impact = impact_sample$poly[[1]]$impact,
                   observed = impact_sample$poly[[1]]$observed)
  for (j in 1:length(impact_sample$poly)){
    df %<>% cbind(impact_sample$poly[[j]]$sampled)
  }
  colnames(df)[6:NCOL(df)] <- paste0('sampled.', 1:length(impact_sample$poly))
  return(df)
}


for (i in 1:5){
  Omegatrue_high_cor <- Omega_true
  Omegatrue_high_cor$eps$hazard_cor <- 0.95
  Omegatrue_low_cor <- Omega_true
  Omegatrue_low_cor$eps$hazard_cor <- 0.05
  impact_sample_high_cor <- SampleImpact(dir, Model, Omegatrue_high_cor %>% addTransfParams(), AlgoParams, dat='all')
  dist_high <- CalcDist(impact_sample_high_cor, AlgoParams)
  print(paste('high',dist_high))
  df_high3 <- flattenImpactSample(impact_sample_high_cor)
  
  impact_sample_low_cor <- SampleImpact(dir, Model, Omegatrue_low_cor %>% addTransfParams(), AlgoParams, dat='all')
  dist_low <- CalcDist(impact_sample_low_cor, AlgoParams)
  print(paste('low',dist_low))
  df_low3 <- flattenImpactSample(impact_sample_low_cor)
  
  impact_sample_true_cor <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams, dat='all')
  dist_true <- CalcDist(impact_sample_true_cor, AlgoParams)
  print(paste('true', dist_true))
  df_true3 <- flattenImpactSample(impact_sample_true_cor)
  
  get_ks_stat_rank_histogram(df)
  print(paste('true',dist_true))
}

for (i in 1:5){
  #Omega_high_cor <- AlgoResults$Omega_sample_phys[which.max(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  #Omega_low_cor <- AlgoResults$Omega_sample_phys[which.min(AlgoResults$Omega_sample_phys[,14,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  Omega_high_cor <- Omega_true
  Omega_high_cor$eps$hazard_cor <- 0.95
  Omega_low_cor <- Omega_true
  Omega_low_cor$eps$hazard_cor <- 0.05
  df_postpredictive_sampled_highcor <- create_df_postpredictive(AlgoResults, Omega=Omega_high_cor, single_particle=T)
  df_postpredictive_sampled_lowcor <- create_df_postpredictive(AlgoResults, Omega=Omega_low_cor, single_particle=T)
  df_postpredictive_true <- create_df_postpredictive(AlgoResults, Omega = Omega_true, single_particle=T)
  
  rank_high <- get_ks_stat_rank_histogram(df_postpredictive_sampled_highcor)
  rank_low <- get_ks_stat_rank_histogram(df_postpredictive_sampled_lowcor)
  rank_true <- get_ks_stat_rank_histogram(df_postpredictive_true)
  hist(rank_true)
  
  print(get_ks_stat_rank_histogram(df_postpredictive_sampled_highcor))
  print(get_ks_stat_rank_histogram(df_postpredictive_sampled_lowcor))
  print(get_ks_stat_rank_histogram(df_postpredictive_true))
}

get_ks_stat_rank_histogram <- function(df){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    mat <- cbind(df[groupings[[i]],5], df[groupings[[i]],6:55])
    z_j <- c()
    for (j in 1:NCOL(mat)){
      n_greater <- 0
      for (k in 1:NCOL(mat[,-j])){
        if(all(mat[,j] <= mat[,-j][,k])){
          n_greater <- n_greater + 1
        }
      }
      z_j <- c(z_j, n_greater)
    }
    z_all <- c(z_all, rank(z_j,  ties.method ='random')[1])
  }
  return(z_all)
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

get_band_depth_rank <- function(df){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    print(i)
    mat <- cbind(df[groupings[[i]],5], df[groupings[[i]],6:55])
    m <- NCOL(mat)
    d <- NROW(mat)
    pre_ranks <- c()
    for (j in 1:m){
      sum <- 0
      for (k in 1:d){
        rank_k <- sum(mat[k,j] <= mat[k,])
        sum <- sum + rank_k * (m - rank_k) + (rank_k-1) * sum(mat[d,j]==mat[k,])
      }
      pre_ranks <- c(pre_ranks, sum /d)
    }
    z_all <- c(z_all, rank(pre_ranks,  ties.method ='random')[1])
  }
  random_runif <-((z_all-runif(length(z_all),0,1))/51)
  hist(random_runif)
  AndersonDarlingTest(random_runif, null='punif')
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

saveRDS(list(
  df_high1 = df_high1,
  df_high2 = df_high2,
  df_high3 = df_high3,
  df_low1 = df_low1,
  df_low2 = df_low2,
  df_low3 = df_low3,
  df_true1 = df_true1,
  df_true2 = df_true2,
  df_true3 = df_true3
), 'ImpactSample_varyingcorrelation')



get_rank_histogram_pairwise <- function(df){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    print(i)
    pairs <- combn(NROW(groupings[[i]]), 2)
    for (combn in 1:NCOL(pairs)){
      mat <- cbind(df[groupings[[i]][pairs[,combn]],5], df[groupings[[i]][pairs[,combn]],6:55])
      z_j <- c()
      for (j in 1:NCOL(mat)){
        z_j <- c(z_j, sum(colSums(mat[,j]<= mat[,-j])==NROW(mat)))
      }
      z_all <- c(z_all, sum(z_j[1] <= z_j[-1]))
      #z_all <- c(z_all, rank(z_j,  ties.method ='random')[1])
    }
  }
  return(z_all)
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

prop_zeros <- c()
for (i in 1:length(groupings)){
  prop_zeros <- c(prop_zeros, mean(df[groupings[[i]],5]<0.5))
}
plot(prop_zeros, z_all)
hist(z_all[prop_zeros<0.5]/51)
AndersonDarlingTest((z_all[prop_zeros<0.5]-runif(length(z_all)))/51, null='punif')$statistic


hist(z_j_true)
abline(v=z_j_true[1], col='red')

hist(z_j_highcor)
abline(v=z_j_highcor[1], col='red')

es_sample(as.numeric(log(df_postpredictive_truelowcor$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_truelowcor$impact[ii]])), 
          as.matrix(sweep(log(df_postpredictive_truelowcor[ii,6:155]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_truelowcor$impact[ii]]), "*")))

es_sample(as.numeric(log(df_postpredictive_truehighcor$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_truehighcor$impact[ii]])), 
          as.matrix(sweep(log(df_postpredictive_truehighcor[ii,6:155]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_truehighcor$impact[ii]]), "*")))


crps_sample(as.numeric(log(df_postpredictive_true$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_true$impact[ii]]))[2], as.matrix(sweep(log(df_postpredictive_true[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_true$impact[ii]]), "*"))[2,])
crps_sample(as.numeric(log(df_postpredictive_sampled_highcor$observed[ii]+AlgoParams$log_offset)*unlist(AlgoParams$impact_weights[df_postpredictive_sampled_highcor$impact[ii]]))[2], as.matrix(sweep(log(df_postpredictive_sampled_highcor[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_sampled_highcor$impact[ii]]), "*"))[2,])

ii <- 801
plot(as.matrix(sweep(log(df_postpredictive_sampled_highcor[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_sampled_highcor$impact[ii]]), "*"))[1,])
points(as.matrix(sweep(log(df_postpredictive_true[ii,6:55]+AlgoParams$log_offset), 1, as.numeric(AlgoParams$impact_weights[df_postpredictive_true$impact[ii]]), "*"))[1,], col='blue')
#Yes, samples are correct:
# ii <- which(impact_sample$poly[[1]]$impact=='displacement' & impact_sample$poly[[1]]$observed==0)[1:3]
# plot(1:length(ii), impact_sample$poly[[1]]$observed[ii], col='red', ylim=c(0, 100000))
# for (i in 1:50){
#   points(1:length(ii), impact_sample$poly[[i]]$sampled[ii], col='blue')
# }
# points(1:length(ii), impact_sample$poly[[1]]$observed[ii], col='red', pch=19)

# ii <- which(df_poly_train$impact=='mortality' & df_poly_train$observed>1)[1:100]
# plot(1:length(ii), df_poly_train$observed[ii], col='red')
# for (i in 1:50){
#   points(1:length(ii), df_poly_train[[paste0('sampled.', i)]][ii], col='blue')
# }
# points(1:length(ii), df_poly_train$observed[ii], col='red', pch=19)

#What do the higher level priors look like?

HLPrior_sample(Model, AlgoParams)


# plot(seq(4.5,10,0.01), log(pnorm(seq(4.5,10,0.01), 10,1.05)))
# lines(seq(4.5,10,0.01), -20+2*seq(4.5,10,0.01))
























