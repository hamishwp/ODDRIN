
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-24_121215_alpha0.95_M60_Npart990RealAgg3')

AlgoResults$input_folder = 'IIDIPUS_Input_RealAgg3/'

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-10_192933_alpha0.9_M60_Npart1000RealAgg5')
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-29_140215_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2')
AlgoResults$input_folder = 'IIDIPUS_Input_RealAgg5/'


Omega_best <- AlgoResults$Omega_sample_phys[192,,61] %>% relist(skeleton=Model$skeleton)
Omega_best <- AlgoResults$Omega_sample_phys[466,,71] %>% relist(skeleton=Model$skeleton)
Omega_best <- AlgoResults$Omega_sample_phys[466,,71] %>% relist(skeleton=Model$skeleton)
Omega_best <- AlgoResults$Omega_sample_phys[626,,80] %>% relist(skeleton=Model$skeleton) # min mst
Omega_best <- AlgoResults$Omega_sample_phys[489,,80] %>% relist(skeleton=Model$skeleton)
  
Omega_best <- AlgoResults$Omega_sample_phys[984,,120] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder='IIDIPUS_Input_RealAgg5/'
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
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, M=100, output='SampledTotal')
#saveRDS(df_postpredictive_sampled_best, paste0(dir, 'IIDIPUS_Results/', 'sampledpost', gsub(gsub(gsub(Sys.time(), pattern = " ", replacement = "_"), pattern = ":", replacement = ""), pattern = "\\..*", replacement = "")))

plot_df_postpredictive(df_postpredictive_sampled_best, 'mortality')
plot_df_postpredictive(df_postpredictive_sampled_best, 'displacement')
plot_df_postpredictive(df_postpredictive_sampled_best, 'buildDam')

plot_cor_posts_poster(AlgoResults)

plot_df_postpredictive_PAGER_coloured(df_postpredictive_sampled_best, 'mortality') #8 x 6.pdf


# obs_sampled <- as.numeric(df_postpredictive_sampled_best[106,c(grep('observed',names(df_postpredictive_sampled_best)),grep('sampled',names(df_postpredictive_sampled_best)))])
# plot(obs_sampled, col=c('red',rep('blue', length(obs_sampled)-1)), pch=19)