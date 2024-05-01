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
Omega <- Omega_true <- list(Lambda1 = list(nu=8.85, kappa=1.6),
                   Lambda2 = list(nu=9.7, kappa=1.8), #list(nu=10.65, kappa=1.5), #
                   Lambda3 = list(nu=8.9, kappa=1.1),
                   Lambda4 = list(nu=9.9, kappa=1.6),
                   theta= list(theta1=0.6),
                   eps=list(local=0.6, hazard_mort=0.3, hazard_disp=0.5, hazard_bd=0.4, hazard_cor=0.55),
                   #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                   vuln_coeff = list(PDens=0, SHDI=-0.18, GNIc=-0.05, Vs30=0.1, EQFreq=-0.25, FirstHaz=0.05, Night=0, FirstHaz.Night=0.1),
                   check = list(check=0.5))

plot_S_curves(Omega_true)

Model$HighLevelPriors(Omega %>% addTransfParams(), Model)

set.seed(1)
simulateDataSet(150, Omega, Model, dir)

#copy into IIDIPUS_Input and create Training + Testing folders



AlgoParams$smc_steps <- 100
AlgoParams$smc_Npart <- 50
AlgoParams$m_CRPS <- 60
AlgoParams$Np <- 1
AlgoParams$smc_alpha <- 0.9
AlgoParams$rel_weightings <- c(0,1)
AlgoParams$kernel_sd <- list(displacement = 1, mortality = 16, buildDam=1.2,
                             buildDest = 0.9, buildDamDest = 1)

tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulatedfull_energyscore_150events_Npart50_M60')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F, tag_notes=tag_notes)

#Speed Test:

start_time <- Sys.time()
impact_sample <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams)
CalcDist(impact_sample, AlgoParams)
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)



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

#compare the number of observations per event:


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
  geom_histogram(data=df_SimImpact %>% filter(impact==impact_type), aes(x=I_max,y=after_stat(count)), alpha=0.3, col='yellow', lwd=0.2, fill='yellow')  +
  scale_y_continuous(breaks = seq(0, 90, by = 10), limits = c(0, 90))

# POINT DATA COMPARISON:


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

plot_correlated_posteriors(AlgoResults, Omega=Omega)
plot_correlated_posteriors(AlgoResults, Omega=Omega, pairings=rbind(c(13,14), c(15,16), c(17,18), c(19,20), c(21,22)))

plot_corr_posterior_vs_d(AlgoResults, Omega=Omega, pairing=c(21,22))

plot_corr_transf_posterior_vs_d(AlgoResults, Omega=Omega, pairing=c(21,22))

plot(AlgoResults$Omega_sample_phys[,2,30], AlgoResults$Omega_sample_phys[,3,30])

plot_acc_prob(AlgoResults); abline(h=0.05)
plot_d_vs_step(AlgoResults, 12)
plot(AlgoResults$essstore)

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
  sqrt(sum((part[c(10:14)]-Omega_true[c(10:14)])^2))
}
eucl_dist <- apply(AlgoResults$Omega_sample_phys[,,AlgoResults$s_finish], 1, eucldist2truepart, unlist(Omega_true))
plot(eucl_dist,AlgoResults$d[,1,AlgoResults$s_finish], ylim=c(6.5, 20))
eucl_dist10 <- apply(AlgoResults$Omega_sample_phys[,,20], 1, eucldist2truepart, unlist(Omega_true))
points(eucl_dist10,AlgoResults$d[,1,20], col='blue')
#particle_min.d <- which(min(AlgoResults$d, na.rm=T)==AlgoResults$d, arr.ind=T)


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


