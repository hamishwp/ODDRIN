
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-24_121215_alpha0.95_M60_Npart990RealAgg3')

AlgoResults$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg3/'

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-10_192933_alpha0.9_M60_Npart1000RealAgg5')
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-29_140215_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2')
AlgoResults$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg5/'

AlgoResultsSMC <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-08-18_204305_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2')

AlgoResultsMCMC <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-08-19_075743_MCMC_RealAgg5_Trial_LR40_Rho0.975propCOVdiv12')

plot(AlgoResultsMCMC$loss, type='l')

plot(AlgoResultsSMC$Omega_sample_phys[,16,1], AlgoResultsSMC$Omega_sample_phys[,17,1], col='blue', xlab='SHDI.coef', ylab='GNIc.coef')
points(AlgoResultsSMC$Omega_sample_phys[,16,163], AlgoResultsSMC$Omega_sample_phys[,17,163], pch=19)
points(AlgoResultsMCMC$Omega_sample_phys[16,], AlgoResultsMCMC$Omega_sample_phys[17,], col='red', pch=19)


#compare MCMC and ABC-SMC:
AlgoResults_smc <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-10_192933_alpha0.9_M60_Npart1000RealAgg5')
AlgoResults_mcmc <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-08-11_123917_MCMC_RealAgg5_Trial_LR30_Rho0.975propCOVdiv12backup')

plot(AlgoResults_smc$Omega_sample_phys[,16,1], AlgoResults_smc$Omega_sample_phys[,17,1], col='blue',xlab='SHDI', ylab='GNIc')
points(AlgoResults_smc$Omega_sample_phys[,16,120], AlgoResults_smc$Omega_sample_phys[,17,120], col='red')
points(AlgoResults_mcmc$Omega_sample_phys[16,], AlgoResults_mcmc$Omega_sample_phys[17,], col='green')


#which(AlgoResults$d==min(AlgoResults$d), arr.ind=T)
Omega_best <- AlgoResults$Omega_sample_phys[192,,61] %>% relist(skeleton=Model$skeleton)
Omega_best <- AlgoResults$Omega_sample_phys[466,,71] %>% relist(skeleton=Model$skeleton)
Omega_best <- AlgoResults$Omega_sample_phys[466,,71] %>% relist(skeleton=Model$skeleton)
Omega_best <- AlgoResults$Omega_sample_phys[626,,80] %>% relist(skeleton=Model$skeleton) # min mst
Omega_best <- AlgoResults$Omega_sample_phys[489,,80] %>% relist(skeleton=Model$skeleton)

Omega_best <- AlgoResults$Omega_sample_phys[401,,130] %>% relist(skeleton=Model$skeleton)

loss_store <- c()
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 60
for (i in 1:10){
  imp_best <- SampleImpact(dir, Model, Omega_best %>% addTransfParams(), AlgoParams, dat='Train')
  es_sample <- CalcDist(imp_best, AlgoParams)[1]
  loss_store <- c(loss_store, es_sample)
}
var(exp(-8*loss_store[1:5])/exp(-8*loss_store[6:10])) # learning rate of 8 seems pretty reasonable, but only estimating variance using a sample of 5!

  
Omega_best <- AlgoResults$Omega_sample_phys[984,,120] %>% relist(skeleton=Model$skeleton)
Omega_best <- AlgoResults$Omega_sample_phys[401,,130] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder='IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg5/'
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
imp_best <- SampleImpact(dir, Model, Omega_best %>% addTransfParams(), AlgoParams, dat='All')

imp_best_plot <- create_df_postpredictive_from_impact_sample(imp_best)
plot_df_postpredictive(imp_best_plot, 'mortality')

#generic plots
plot_posteriors(AlgoResults, 15:22)
plot_d_vs_step(AlgoResults, ymax=6)
plot_acc_prob(AlgoResults); abline(h=0.05)

plot_correlated_posteriors(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(9,10), c(11,12), c(13,14)))

#plot posterior predictive distributions

#df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, M=100, output='SampledTotal')
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=T, 
                                                           Omega = AlgoResults$Omega_sample_phys[401,,130] %>% relist(skeleton=Model$skeleton),
                                                           M=100, output='SampledTotal')

#saveRDS(df_postpredictive_sampled_best, paste0(dir, 'IIDIPUS_Results/', 'sampledpost', gsub(gsub(gsub(Sys.time(), pattern = " ", replacement = "_"), pattern = ":", replacement = ""), pattern = "\\..*", replacement = "")))

plot_df_postpredictive(df_postpredictive_sampled_best, 'mortality')
plot_df_postpredictive(df_postpredictive_sampled_best, 'displacement')
plot_df_postpredictive(df_postpredictive_sampled_best, 'buildDam')

plot_cor_posts_poster(AlgoResults)

plot_df_postpredictive_PAGER_coloured(df_postpredictive_sampled_best, 'mortality') #8 x 6.pdf


# obs_sampled <- as.numeric(df_postpredictive_sampled_best[106,c(grep('observed',names(df_postpredictive_sampled_best)),grep('sampled',names(df_postpredictive_sampled_best)))])
# plot(obs_sampled, col=c('red',rep('blue', length(obs_sampled)-1)), pch=19)

