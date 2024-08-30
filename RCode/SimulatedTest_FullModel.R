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
packred<-T


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

Omega <- Omega_true <- list(Lambda1 = list(nu=8.75, kappa=0.6),
                            Lambda2 = list(nu=11.7, kappa=0.75), #list(nu=10.65, kappa=1.5), #
                            Lambda3 = list(nu=8.7, kappa=0.7),
                            Lambda4 = list(nu=9.9, kappa=1.6),
                            theta= list(theta1=0.6),
                            eps=list(local=0.8, hazard_mort=0.45, hazard_disp=0.6, hazard_bd=0.5, hazard_cor=0.55),
                            #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                            vuln_coeff = list(PDens=0, SHDI=-0.18, GNIc=-0.05, Vs30=0.1, EQFreq=-0.12, FirstHaz=0.05, Night=0, FirstHaz.Night=0.1),
                            check = list(check=0.5))

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
n_samples <- 300
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
simulateDataSet(170, Omega, Model, dir, folder_write='IIDIPUS_Input_Alternatives/IIDIPUS_SimInput/')


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
AlgoParams$m_CRPS <- 60
AlgoParams$smc_Npart <- 500
AlgoParams$n_nodes <- 1
AlgoParams$smc_steps <- 100
AlgoParams$rel_weightings <- c(1,1)
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput/'

tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '500parttest_0.1propcov')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)


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


#------------------------------------------------------------------------------------------------
#----------------------------------- PLOT SIMULATED DATA ----------------------------------------
#--------------------------------(and compare to real data)--------------------------------------
#------------------------------------------------------------------------------------------------

moveTestData <- function(folder_in='IIDIPUS_Input_Alternatives/IIDIPUS_SimInput'){
  ODD_folderall<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_foldertest<-paste0(dir, folder_in, '/Test/')
  ufiles<-list.files(path=ODD_folderall,pattern=Model$haz,recursive = T,ignore.case = T)
  i <- 0
  for (file in ufiles){
    i <- i + 1
    if (i %%3 != 0){next}
    file.copy(from = paste0(ODD_folderall, file),
              to = paste0(ODD_foldertest, file))
    file.remove(from = paste0(ODD_folderall, file))
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
moveTestData('IIDIPUS_Input_Alternatives/IIDIPUS_SimInput')

# Collect mortality, building damage, and displacement data for simulated data:
ODDsim_paths <-na.omit(list.files(path="IIDIPUS_Input_Alternatives/IIDIPUS_SimInput/ODDobjects/"))
df_SimImpact <- data.frame(observed=numeric(),
                           impact=character(),
                           polygon=integer(),
                           exposure=numeric(),
                           event=integer(),
                           I_max=numeric())
nHazSim <- c()
maxIntSim <- c()
for(i in 1:length(ODDsim_paths)){
  ODDSim <- readRDS(paste0("IIDIPUS_Input_Alternatives/IIDIPUS_SimInput/ODDobjects/",ODDsim_paths[i]))
  if (length(ODDSim@impact$impact)>0){
    nHazSim <- c(nHazSim, length(grep('hazMean', names(ODDSim@data))))
    maxIntSim <- c(maxIntSim, max(ODDSim@data[, grep('hazMean', colnames(ODDSim@data))],  na.rm=T))
    # if (length(grep('hazMean', names(ODDSim@data))) ==1 & (length(grep('-4', ODDsim_paths[i]))==0)& (length(grep('-5', ODDsim_paths[i]))==0)){
    #   stop()
    # }
    for (j in 1:NROW(ODDSim@impact)){
      df_SimImpact %<>% add_row(observed=ODDSim@impact$observed[j], impact=ODDSim@impact$impact[j], polygon=ODDSim@impact$polygon[j],
                                exposure=ifelse(impact=='buildDam', sum(ODDSim@data[ODDSim@polygons[[ODDSim@impact$polygon[j]]]$indexes,'nBuildings']), sum(ODDSim@data[ODDSim@polygons[[ODDSim@impact$polygon[j]]]$indexes,'Population'])),
                                event=i, I_max=max(ODDSim@data[ODDSim@polygons[[ODDSim@impact$polygon[j]]]$indexes, grep('hazMean', colnames(ODDSim@data))],  na.rm=T))
    }
  }
}
ggplot(df_SimImpact %>% filter(impact=='mortality'), aes(x=I_max, y=observed)) + geom_point()

# Collect mortality, building damage, and displacement data for real data:
#ODDpath <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/'

ODDpath <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg5/ODDobjects/Train/'
ODDpaths <-na.omit(list.files(path=ODDpath))
df_Impact <- data.frame(observed=numeric(), impact=character(),
                        polygon=integer(), exposure=numeric(), event=integer(), I_max=numeric())

nHazReal <- c()
maxIntReal <- c()
for(i in 1:length(ODDpaths)){
  ODD <- readRDS(paste0(ODDpath,ODDpaths[i]))
  if (length(ODD@impact$impact)>0){
    nHazReal <- c(nHazReal, length(grep('hazMean', names(ODD@data))))
    maxIntReal <- c(maxIntReal, max(ODD@data[, grep('hazMean', colnames(ODD@data))],  na.rm=T))
    for (j in 1:NROW(ODD@impact)){
      df_Impact %<>% add_row(observed=ODD@impact$observed[j], impact=ODD@impact$impact[j], polygon=ODD@impact$polygon[j],
                             exposure=ifelse(impact=='buildDam', sum(ODD@data[ODD@polygons[[ODD@impact$polygon[j]]]$indexes,'nBuildings']), sum(ODD@data[ODD@polygons[[ODD@impact$polygon[j]]]$indexes,'Population'])),
                              event=i, I_max=max(ODD@data[ODD@polygons[[ODD@impact$polygon[j]]]$indexes, grep('hazMean', colnames(ODD@data))],  na.rm=T))
    }
  }
}

ggplot(df_Impact %>% filter(impact=='mortality'), aes(x=I_max, y=observed)) + geom_point()

#Plots:
# plot of mortality (simulated) vs mortality (observed)
library(cowplot)
library(scales)
plot_true_vs_simulated_obsvals <- function(impact_type){
  p <- ggplot()  + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000), labels = label_comma())  + 
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000), labels = label_comma()) +
    geom_histogram(data=df_SimImpact %>% filter(impact==impact_type), aes(x=observed,y=after_stat(count)), alpha=0.4, col="blue", lwd=0.2, fill="blue") +
    geom_histogram(data=df_Impact %>% filter(impact==impact_type), aes(x=observed,y=after_stat(count)), alpha=0.4, col='red', lwd=0.2, fill='red') +
    theme_bw() + ylab('Count') + xlab('Oberved')
  return(p)
}
plot_true_vs_simulated_obscount <- function(impact_type){
  p <- ggplot() +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100), labels = label_comma())  + 
    geom_histogram(data=df_SimImpact %>% filter(impact==impact_type) %>% group_by(event) %>% tally(), aes(x=n, y=after_stat(count), fill='Simulated Data', col='Simulated Data'), alpha=0.3, col='blue', lwd=0.2) + 
    geom_histogram(data=df_Impact %>% filter(impact==impact_type) %>% group_by(event) %>% tally(), aes(x=n, y=after_stat(count), fill='Real Data', col='Real Data'), alpha=0.3, col='red', lwd=0.2) +
    ylab('Count') + xlab('Number of Observations per event') + theme_bw() +
    scale_fill_manual(values = c("Real Data" = "red", "Simulated Data" = "blue"))+ labs(fill = "") +
    guides(fill = guide_legend(override.aes = list(color = NULL)))
  return(p)
}

p_mort_obsvals <- plot_true_vs_simulated_obsvals('mortality') + xlab('Observed Mortality')
p_disp_obsvals <- plot_true_vs_simulated_obsvals('displacement') + xlab('Observed Displacement')
p_bd_obsvals <- plot_true_vs_simulated_obsvals('buildDam') + xlab('Observed Building Damage')
p_mort_obscount <- plot_true_vs_simulated_obscount('mortality') + xlab('Number of Mortality Observations per Event') 
p_disp_obscount <- plot_true_vs_simulated_obscount('displacement') + xlab('Number of Displacement Observations per Event')
p_bd_obscount <- plot_true_vs_simulated_obscount('buildDam') + xlab('Number of Building Damage Observations per Event')
legend <- get_legend(p_mort_obscount + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
plot_grid(legend, plot_grid( p_mort_obsvals, p_disp_obsvals, p_bd_obsvals, 
                          p_mort_obscount + theme(legend.position="none"), p_disp_obscount + theme(legend.position="none"), 
                          p_bd_obscount + theme(legend.position="none"),
                          align = 'vh', hjust = -1, nrow = 2
), ncol = 1, rel_heights=c(0.1,1))

grid.arrange(p_mort_obsvals, p_disp_obsvals, p_bd_obsvals, 
             p_mort_obscount, p_disp_obscount, p_bd_obscount, ncol=3, nrow=2)

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


