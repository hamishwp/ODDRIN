
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-19_101900_alpha0.9_M60_Npart990RealAgg3')

AlgoResults$input_folder = 'IIDIPUS_Input_RealAgg3/'

#generic plots
plot_posteriors(AlgoResults, 15:20)
plot_d_vs_step(AlgoResults, ymax=7)
plot_acc_prob(AlgoResults); abline(h=0.05)

#plot posterior predictive distributions
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, M=100, output='SampledTotal')
#saveRDS(df_postpredictive_sampled_best, paste0(dir, 'IIDIPUS_Results/', 'sampledbest', gsub(gsub(gsub(Sys.time(), pattern = " ", replacement = "_"), pattern = ":", replacement = ""), pattern = "\\..*", replacement = "")))

plot_df_postpredictive(df_postpredictive_sampled_best, 'mortality')
plot_df_postpredictive(df_postpredictive_sampled_best, 'displacement')
plot_df_postpredictive(df_postpredictive_sampled_best, 'buildDam')

plot_cor_posts_poster(AlgoResults)

plot_df_postpredictive_PAGER_coloured(df_postpredictive_sampled_best, 'mortality')

# obs_sampled <- as.numeric(df_postpredictive_sampled_best[106,c(grep('observed',names(df_postpredictive_sampled_best)),grep('sampled',names(df_postpredictive_sampled_best)))])
# plot(obs_sampled, col=c('red',rep('blue', length(obs_sampled)-1)), pch=19)