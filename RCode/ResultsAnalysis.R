
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#--------------------------------------------------------------------------------------------------
#------------------------------ Simulated Data Results Analysis -----------------------------------
#--------------------------------------------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source('RCode/ResultsAnalysisFunctions.R')
library(patchwork)

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2025-10-18_084656.908241_SimInput13_NormalCDF_ESplus0.1AVE_RF_range1')
#plot(AlgoResults$tolerancestore[2:length(AlgoResults$tolerancestore)]/AlgoResults$tolerancestore[1:(length(AlgoResults$tolerancestore)-1)], ylim=c(0.99,1))

AlgoParams$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput13/'


plot_correlated_posteriors_ggplot(AlgoResults,  pairings = rbind(c(1,2), c(3,4), c(5,6), c(8,7), c(9,10), c(9,11), c(15,12), c(13, 14), c(16,17), c(18,19)), 
                                  Omega = list(Lambda1 = list(nu=7.5, kappa=0.6),
                                               Lambda2 = list(nu=12.3, kappa=1.4),#Lambda2 = list(nu=10.6, kappa=0.95),
                                               Lambda3 = list(nu=7.9, kappa=0.68),
                                               eps=list(local=0.3, hazard_mort=0.4, hazard_disp=0.55, hazard_bd=0.5, hazard_cor=0.6),
                                               vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15)), subfig_title_adj=-0.8,
                                  s_finish=350)
#SimPosteriors2.pdf, 5 x 12 inches

AlgoResults_real <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1')
plot(AlgoResults$tolerancestore[2:length(AlgoResults$tolerancestore)]/AlgoResults$tolerancestore[1:(length(AlgoResults$tolerancestore)-1)], ylim=c(0.99, 1))
points(AlgoResults_real$tolerancestore[2:length(AlgoResults_real$tolerancestore)]/AlgoResults_real$tolerancestore[1:(length(AlgoResults_real$tolerancestore)-1)], ylim=c(0.99, 1), col='red')

plot(AlgoResults$accrate_store, ylim=c(0, 0.2))
points(AlgoResults_real$accrate_store, col='red')

plot(AlgoResults$accrate_store, ylim=c(0, 0.2))
plot(AlgoResults$tolerancestore, ylim=c(0, 5))


AlgoParams$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput13/'
AlgoResults$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput13/'

df_postpredictive_sampled <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                      M=200, output='SampledTotal', s_finish = 350)

df_postpredictive_sampled_true <- create_df_postpredictive(AlgoResults, single_particle=T, 
                                                           list(Lambda1 = list(nu=7.5, kappa=0.6),
                                                                Lambda2 = list(nu=12.3, kappa=1.4),#Lambda2 = list(nu=10.6, kappa=0.95),
                                                                Lambda3 = list(nu=7.9, kappa=0.68),
                                                                eps=list(local=0.3, hazard_mort=0.4, hazard_disp=0.55, hazard_bd=0.5, hazard_cor=0.6),
                                                                vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15)),
                                                           M=200, output='SampledTotal')

saveRDS(df_postpredictive_sampled_true, 'IIDIPUS_Results/SimPostPredictiveOct22_2025')
saveRDS(df_postpredictive_sampled, 'IIDIPUS_Results/SimTruePredictiveOct22_2025')

#df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SimPostPredictiveAug2025')
#df_postpredictive_sampled_true <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/SimTruePredictiveAug2025')

plot_df_postpredictive(df_postpredictive_sampled_true, 'mortality')

df_plot <- bind_rows(
  df_postpredictive_sampled_true %>%
    filter(impact == "mortality") %>%
    transmute(val = log10(observed+10), source = "postpredictive"),
  sampled_best_realdat %>%
    filter(impact == "mortality") %>%
    transmute(val = log10(observed+10), source = "realdat")
)

ggplot(df_plot, aes(val, fill = source)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(x = "log10(observed)", y = "count", title = "Overlayed histograms") +
  theme_minimal()

plot_df_postpredictive_compare(df_postpredictive_sampled_true %>% filter(train_flag=='TEST'),
                               df_postpredictive_sampled %>% filter(train_flag=='TEST'), 'mortality')

#PostPredictiveCompare.pdf, 8 x 5

df_postpredictive_sampled_best %<>% filter(observed>0)
df_postpredictive_sampled_true %<>% filter(observed>0)
plot_df_postpredictive_compare(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),
                               df_postpredictive_sampled_true %>% filter(train_flag=='TEST'), 'mortality')
#PostPredSim.pdf, 5 x 9 inches

# mat <- df_postpredictive_sampled_true %>%
#   filter(impact == "mortality") %>%
#   dplyr::select(5:105) %>%
#   as.matrix()
# 
# # rank of the FIRST row within each column
# r_first <- apply(mat, 1, function(x) rank(x, ties.method = "random")[1])



#--------------------------------------------------------------------------------------------------
#-------------------------------------- Simulated Trace Plots -------------------------------------
#--------------------------------------------------------------------------------------------------

source('RCode/ResultsAnalysisFunctions.R')
library(patchwork)

mcmc_sim1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-10-27_155104.365377__Sim13_NormalCDF_ESplus.1AVErank_RFwithTot_range1_Chain1_backup')
mcmc_sim2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-10-27_155036.357656__Sim13_NormalCDF_ESplus.1AVErank_RFwithTot_range1_Chain2_backup')
mcmc_sim3 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-10-27_155247.811179__Sim13_NormalCDF_ESplus.1AVErank_RFwithTot_range1_Chain3')

mcmc_trace(mcmc_sim1, mcmc_sim2, mcmc_sim3, var_plot=1:19, xlim=c(1,25000), Omega=Omega_true)
#MCMC_Sim_Trace.pdf, 9 x 9


library(coda)

# matrices: iterations x variables for each chain
mat1 <- t(mcmc_sim1$Omega_sample_phys[, seq(20001, 25000, 10), drop = FALSE])
mat2 <- t(mcmc_sim2$Omega_sample_phys[, seq(20001, 25000, 10), drop = FALSE])
mat3 <- t(mcmc_sim3$Omega_sample_phys[, seq(20001, 25000, 10), drop = FALSE])

gd <- gelman.diag(mcmc.list(mcmc(mat1), mcmc(mat2), mcmc(mat3)),
                  autoburnin = FALSE, multivariate = FALSE)
rhat_vec <- gd$psrf[, "Point est."]                # one per parameter
rhat_vec

plot_correlated_posteriors_ggplot(AlgoResults, include_priors=T, Omega = list(Lambda1 = list(nu=7.5, kappa=0.6),
                                                                                    Lambda2 = list(nu=12.3, kappa=1.4),#Lambda2 = list(nu=10.6, kappa=0.95),
                                                                                    Lambda3 = list(nu=7.9, kappa=0.68),
                                                                                    eps=list(local=0.3, hazard_mort=0.4, hazard_disp=0.55, hazard_bd=0.5, hazard_cor=0.6),
                                                                                    vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15)),
                                  pairings = rbind(c(1,2), c(3,4), c(5,6), c(8,7), c(9,10), c(9,11), c(15,12), c(13, 14), c(16,17), c(18,19)), s_finish=350, 
                                  AlgoResultsMCMC=list(mcmc_sim1, mcmc_sim2, mcmc_sim3), subfig_title_adj=-0.3, mcmc_warmup = 20000, plots_ncol=4)

# MCMC_SMC_Compare.pdf, 9 x 7


#---------------------Compare with and without rank histograms-----------


ESonly1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-10-23_142412.52814__Sim13_NormalCDF_ESonly_RFwithTot_range1_Chain1')
ESonly2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-10-23_142548.329482__Sim13_NormalCDF_ESonly_RFwithTot_range1_Chain2')

mcmc_trace(ESonly1, ESonly2, var_plot=1:19, xlim=c(1, 15000))


xlim = 20001:25000
xlim_esonly = 8001:12000
hist(c(ESonly1$Omega_sample_phys[11,xlim_esonly], ESonly2$Omega_sample_phys[11,xlim_esonly]))
hist(c(mcmc_sim1$Omega_sample_phys[11, xlim], mcmc_sim2$Omega_sample_phys[11, xlim], mcmc_sim3$Omega_sample_phys[11, xlim]))



plot_correlated_posteriors_ggplot(AlgoResults, include_priors=T, 
                                  pairings = rbind(c(1,2), c(3,4), c(5,6), c(8,7), c(9,10), c(9,11), c(15,12), c(13, 14), c(16,17), c(18,19)), 
                                  Omega = list(Lambda1 = list(nu=7.5, kappa=0.6),
                                               Lambda2 = list(nu=12.3, kappa=1.4),#Lambda2 = list(nu=10.6, kappa=0.95),
                                               Lambda3 = list(nu=7.9, kappa=0.68),
                                               eps=list(local=0.3, hazard_mort=0.4, hazard_disp=0.55, hazard_bd=0.5, hazard_cor=0.6),
                                               vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15)), 
                                  s_finish=350, 
                                  AlgoResultsMCMC=list(ESonly1, ESonly2), subfig_title_adj=-0.3, mcmc_warmup = 10000)

plot_compare_two_mcmc_posteriors_ggplot(AlgoResults,
                                        include_priors = TRUE,
                                        Omega = list(Lambda1 = list(nu=7.5, kappa=0.6),
                                                     Lambda2 = list(nu=12.3, kappa=1.4),#Lambda2 = list(nu=10.6, kappa=0.95),
                                                     Lambda3 = list(nu=7.9, kappa=0.68),
                                                     eps=list(local=0.3, hazard_mort=0.4, hazard_disp=0.55, hazard_bd=0.5, hazard_cor=0.6),
                                                     vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15)),
                                        pairings = rbind(c(1,2), c(3,4), c(5,6), c(8,7), c(9,10), c(9,11), c(15,12), c(13, 14), c(16,17), c(18,19)),
                                        s_finish = NULL,
                                        AlgoResultsMCMC_ES = list(ESonly1, ESonly2),   # Energy Score Only
                                        AlgoResultsMCMC_RH = list(mcmc_sim1, mcmc_sim2, mcmc_sim3),   # with Rank Histogram
                                        mcmc_warmup = 20000,
                                        mcmc_warmup_es = 15000,
                                        subfig_title_adj = -0.4,
                                        plots_ncol = 4)
# MCMC_RankHist_Compare.pdf, 9 x 7

p <- make_param_hist_overlay_linear(
  ESonly1, ESonly2, mcmc_sim1, mcmc_sim2, mcmc_sim3,
  i = 11, xlim_esonly = 15001:20000, xlim_mcmc = 20001:25000,
  nbins = 20
)
print(p)

post_samples = 15001:20000
low_rho_i = which.min(ESonly1$Omega_sample_phys[11,post_samples])
Omega_low_rho = ESonly1$Omega_sample_phys[,post_samples[low_rho_i]] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder = "IIDIPUS_Input_Alternatives/IIDIPUS_SimInput13/"
AlgoParams$cores = 1
AlgoParams$NestedCores = 1
sampled_impact = SampleImpact(dir, Model, Omega_low_rho %>% addTransfParams(), AlgoParams)
ranks = ReturnRankHist(sampled_impact$poly, AlgoParams)

ranks1 = plot_ranks(ranks)

post_samples = 18001:20000
high_rho_i = which.max(ESonly1$Omega_sample_phys[11,post_samples])
Omega_high_rho = ESonly1$Omega_sample_phys[,post_samples[high_rho_i]] %>% relist(skeleton=Model$skeleton)
sampled_impact2 = SampleImpact(dir, Model, Omega_high_rho %>% addTransfParams(), AlgoParams)
ranks_high = ReturnRankHist(sampled_impact2$poly, AlgoParams)

ranks2 = plot_ranks(ranks_high)

plot_grid(p, plot_grid(ranks1, ranks2, nrow=1,
                       labels=c('(b)', '(c)'),
                       label_fontfamily = "Times New Roman",
                       label_fontface = "plain",
                       label_x = -0.02, 
                       label_y = 1), nrow=2, rel_heights = c(0.6, 0.4),
          labels=c('(a)', ''),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain",
          label_x = 0, 
          label_y = 0.96)

# RankHistResults.pdf, 8.5 x 7.5

#------------------------------------------------------------------------


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#--------------------------------------------------------------------------------------------------
#------------------------------ Real Data Results Analysis ----------------------------------------
#--------------------------------------------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('RCode/ResultsAnalysisFunctions.R')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1')
AlgoResults$s_finish = 250

#df <- as.data.frame(AlgoResults$Omega_sample_phys[, 12:19,AlgoResults$s_finish])
#ggpairs(df)

#points(AlgoResults$tolerancestore[2:length(AlgoResults$tolerancestore)]/AlgoResults$tolerancestore[1:(length(AlgoResults$tolerancestore)-1)], ylim=c(0.99,1), col='red')

#check variability in loss at posterior:
AlgoParams$cores = 1
AlgoParams$NestedCores = 1
nreps = 10
dists = array(NA, dim=nreps)
AlgoParams$input_folder = 'IIDIPUS_Input_Alternatives/Sept25/'
for (i in 1:nreps){
  impact_sample = SampleImpact(dir, Model, proposed = AlgoResults$Omega_sample_phys[2,,250] %>% 
                                                          relist(skeleton=Model$skeleton) %>%
                                                          addTransfParams(), AlgoParams)
  dists[i] = sum(CalcDist(impact_sample, AlgoParams)[1:2])
}
# source('RCode/Model Variations/ODDobj_NoTot.R')
# 
# source('RCode/Model Variations/ODDobjJune.R')
# source('RCode/Model Variations/Model_Logistic2.R')
# AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2025-06-12_115400.392519__Apr25Agg_LogLinMort_ESplus0.05BD_RFwithTot_range.5_kappa.5_withHLPriors_flexHLP_backup')
# #AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-12-10_060355_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2')
# 
# 
# source('RCode/ODDobj.R')
# source('RCode/Model Variations/Model_MVRrank.R')
# AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-06-17_151849.678089__Apr25Agg_LogNormal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_withHLPriors_flexHLP')
# #AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-12-10_060355_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2')
# 
# source('RCode/ODDobj.R')
# source('RCode/Model.R')
# AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-08_110537.972824__Apr25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_backup')

plot_correlated_posteriors_ggplot(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(8, 7), c(8,9), c(10,11)), plots_ncol=3, s_finish=250)
#RealPosteriors.pdf, 9 x 5.5 inches




library(grid)
plot_vuln_posteriors(AlgoResults)
#VulnPosteriors.pdf, 9 x 4.5 inches

#AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-13_234229.627104__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain2_backup')
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/Sept25/'
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=200, output='SampledTotal')
saveRDS(df_postpredictive_sampled_best, 'IIDIPUS_Results/fit_MVR_Oct18_Total')
#saveRDS(df_postpredictive_sampled_best, 'SampledTotal8June')

df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=200, output='SampledAgg')
saveRDS(df_postpredictive_sampled_best, 'IIDIPUS_Results/fit_MVR_Oct18_Agg')

#df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/fit_MVR_June18')

quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='mortality'& df_postpredictive_sampled_best$train_flag=='TRAIN'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/(length(grep('sampled.', names(df_postpredictive_sampled_best)))+1)
hist(quants)
mean(quants>0.025 & quants < 0.975)

quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='mortality'& df_postpredictive_sampled_best$train_flag=='TEST'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/(length(grep('sampled.', names(df_postpredictive_sampled_best)))+1)
hist(quants)
mean(quants>0.025 & quants < 0.975)
quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='displacement'& df_postpredictive_sampled_best$train_flag=='TEST'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/(length(grep('sampled.', names(df_postpredictive_sampled_best)))+1)
hist(quants)
mean(quants>0.025 & quants < 0.975)
quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='buildDam'& df_postpredictive_sampled_best$train_flag=='TEST'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/(length(grep('sampled.', names(df_postpredictive_sampled_best)))+1)
hist(quants)
mean(quants>0.025 & quants < 0.975)

library(grid)
mort_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'mortality')  + guides(color="none") 
disp_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'displacement')  + guides(color="none")
buildDam_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'buildDam')  + guides(color="none")# + 
  #scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma, limits=c(50,100000))


p1 <- arrangeGrob(mort_test, top = textGrob("(a)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p2 <- arrangeGrob(disp_test, top = textGrob("(b)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p3 <- arrangeGrob(buildDam_test, top = textGrob("(c)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))

plot_grid(mort_test, disp_test, buildDam_test, nrow=2, labels=c('(a)', '(b)', '(c)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain",
          label_x = -0.025, 
          label_y = 0.99)

#PostPredReal.pdf, 8 x 8

grid.arrange(p1, p2, p3, ncol=3)
#PostPredReal.pdf, 4 x 12.1 inches



#--------------------------------------------------------------------------------------------------
#------------------------------ Residuals by Region -----------------------------------------------
#--------------------------------------------------------------------------------------------------

df_postpredictive_sampled_best <- readRDS('IIDIPUS_Results/fit_MVR_Oct17')

df_post_with_covar = add_covar( df_postpredictive_sampled_best, covars=c('hazMean', 'EQFreq', 'Vs30', 'SHDI', 'PDens', 'region', 'shakeSD'), dat='all')

df_post_with_covar %<>% filter(train_flag=='TEST')
impact_type ='mortality'

resample_obs_quant = function(df_poly){
  df_poly$obs_quant = apply(df_poly[, c('observed', grep('sampled.', names(df_poly), value=T, fixed=T))],1, function(x) sample_quant(x))
  df_poly$obs_quant = (df_poly$obs_quant - runif(length(df_poly$obs_quant), 0, 1))/length(c('observed', grep('sampled.', names(df_poly), value=T, fixed=T)))
  return(df_poly)
}

df_post_with_covar %<>% resample_obs_quant()
ggplot(df_post_with_covar %>% filter(impact==impact_type), aes(y = obs_quant)) +
  geom_point(aes(x = shakeSD,            color = "Maximum uncertainty")) +
  geom_point(aes(x = shakeSD_epicentre,  color = "Uncertainty at epicentre")) +
  scale_color_manual(
    name = "",
    values = c("Maximum uncertainty" = "#440154", "Uncertainty at epicentre" = "#1FA187")
  ) +
  ggtitle('') +
  ylab('Observation Quantile') + xlab('ShakeMap Uncertainty') +
  theme_minimal(base_family = "Times New Roman", base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.text      = ggplot2::element_text(family = "Liberation Serif", size = 13)
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 0, 10, 10)) +
  scale_y_continuous(expand = c(0, 0))

#ShakeSDvsQuantile.pdf, 7.5 x 4


#df_post_with_covar$sampled_median = apply(df_post_with_covar[, grep('sampled', names(df_post_with_covar))], 1, median)

#df_post_with_covar$log_dif = log(df_post_with_covar$observed+10) - log(df_post_with_covar$sampled_median+10)

p1 = ggplot(df_post_with_covar %>% filter(impact==impact_type), aes(x=region, y=obs_quant)) +
  #geom_point(aes(color=as.factor(train_flag))) + 
  geom_point(col='#440154') +
  #geom_point() + 
  scale_color_viridis() +
  ggtitle('') +
  ylab('Observation Quantile') + xlab('') +
  theme_minimal(base_family = "Times New Roman", base_size = 14) + 
  theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1), 
        legend.position='none',
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))+ #geom_text(hjust=0, vjust=0.5, size=3, alpha=0.5) + 
  coord_cartesian(clip="off")+
  theme(plot.margin = margin(10, 50, 10, 10)) + scale_y_continuous(expand=c(0,0))
# Region vs Residual, 8.5 x 5

p2 <- ggplot() +
  borders("world", colour = "gray80", fill = "gray80") +
  geom_point(
    data = df_post_with_covar,
    aes(x = lon, y = lat, color = obs_quant),
    size = 2.5, alpha = 0.8
  ) +
  scale_color_viridis(name = "Observation\nQuantile") +
  #scale_color_gradientn(
    #colors = c("darkgreen", "chartreuse3", "gold", "orange", "red"),
    #colors = c("purple", "white", "green"),
    #name = "Observation\nQuantile"
  #) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_family = "Liberation Serif", base_size = 14) +
  theme(
    #axis.title  = element_text(size = 12.5),
    #legend.title = element_text(size = 12.5),
    #legend.text  = element_text(size = 12.5),
    
    # place legend INSIDE the plotting area
    legend.position = c(0.025, 0.1),          # (x, y) in NPC; tweak as needed
    legend.justification = c(0, 0),           # anchor bottom-left of the legend box
    
    # add border and background to the legend box
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      colour = "black",
      linewidth = 0.5
    ),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(4, 6, 4, 6),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  guides(color = guide_colorbar(
    title.position = "top",
    barheight = unit(50, "pt"),
    barwidth  = unit(8, "pt"),
    ticks = TRUE
  )) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

plot_grid(p1 +
            theme(
              plot.margin = margin(t = -10, b = -20, l = 7, r = 7) # 20pt top and bottom margin
            )
          ,p2,  nrow=2, rel_heights = c(1.1, 0.8),
          labels=c('(a)', '(b)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain",
          label_x = c(-0.01, -0.01), 
          label_y = c(1,1))

#Residual_Map.pdf, 7.5 x 9.5 portrait

#RSS 2:

plot_grid(p1 +
            theme(
              plot.margin = margin(t = -10, b = -20, l = 7, r = 7) # 20pt top and bottom margin
            )
          ,p2,  nrow=1, rel_widths = c(1.1, 0.8),
          labels=c('(a)', '(b)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain",
          label_x = c(-0.01, -0.01), 
          label_y = c(1,1))

#Residual_Map_Flat.pdf, 18 x 5


#----------- Residuals by date, RSS: ------------

subnatData = readSubNatData('EQ_SubNational.xlsx')
subnatData_dates = subnatData[, c('EventID', 'sdate')] %>% unique()

df_post_with_covar %<>% merge(subnatData_dates, by.x='event_id', by.y='EventID')

ggplot(df_post_with_covar %>% filter(impact==impact_type), aes(x=sdate, y=obs_quant)) +
  #geom_point(aes(color=as.factor(train_flag))) + 
  geom_point(col='#440154') +
  #geom_point() + 
  scale_color_viridis() +
  scale_x_date(
    limits = as.Date(c("2011-01-01", "2023-12-31")),
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  ggtitle('') +
  ylab('Observation Quantile') + xlab('Event Date') +
  theme_minimal(base_family = "Times New Roman", base_size = 14) + 
  theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1), 
        legend.position='none',
        panel.border = element_rect(color = "black", fill = NA, size = 0.4))+ #geom_text(hjust=0, vjust=0.5, size=3, alpha=0.5) + 
  coord_cartesian(clip="off")+
  theme(plot.margin = margin(10, 50, 10, 10)) + scale_y_continuous(expand=c(0,0))

#Residual_Time.pdf, 9 x 5

#------------------------------------------------



#--------------------------------------------------------------------------------------------------
#------------------------------ ABC-MCMC Trace (and comparison to ABC-SMC)-------------------------
#--------------------------------------------------------------------------------------------------

SR1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-03_092825.199826__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain1_backup')
#SR2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-16_213626.504342__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain2')
#SR3 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-22_094221.160343__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain3')
SR3 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-06_112257.598195__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain3')
SR4 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-03_101635.044837__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain4')
#/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-22_094221.160343__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain3
#/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-23_121212.470563__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain4

#SR1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2025-07-31_094827.645236__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain2_backup')
plot_mcmc_compare(SR1, SR3, SR4, xlim=c(1, 15000))

plot_mcmc_compare = function(SR1, SR2, SR3, xlim=c(1, 10000)){
  
  # Define subfigure labels (a), (b), (c), ...
  subfig_labels <- letters[1:19]
  plot_list <- list()
  var_plot <- 1:19
  var_plot[12] = 15
  var_plot[13:15] = 12:14
  # Loop through the variables and generate each plot
  for (i in seq_along(var_plot)){#seq_along(vuln_var)) {
    #var <- vuln_var[i]
    var <- var_plot[i]
    
    post_samples <- data.frame(
      x=(xlim[1]:xlim[2])-xlim[1],
      SR1_param <- SR1$Omega_sample_phys[var,xlim[1]:xlim[2]],
      SR2_param <-SR2$Omega_sample_phys[var,xlim[1]:xlim[2]],
      SR3_param <-SR3$Omega_sample_phys[var,xlim[1]:xlim[2]])
    
    post_samples_long <- post_samples %>%
      pivot_longer(cols = -x, names_to = "Parameter", values_to = "Value")
    
    if (length(grep('vuln', names(unlist(Omega))[var]))>0){
      ylim=c(-1,1)
    } else {
      ylim=c(Model$par_lb[var], Model$par_ub[var])
    }
    # Plot the lines
    plot <- ggplot(data = post_samples_long, aes(x = x, y = Value, color = Parameter)) +
      labs(y = get_greek_titles(names(unlist(Model$skeleton))[var]),
           x= 'Iteration') +
      scale_y_continuous(limits = ylim, expand = c(0.01, 0.01)) +
      geom_line(alpha=0.8) +
      theme_minimal() +
      scale_color_viridis(discrete = TRUE) +
      theme(
        axis.title.y = element_text(family = "Times New Roman", size = 12),  # Change y-axis title font
        axis.text.x = element_text(family = "Times New Roman", size = 12),   # Change x-axis text font
        axis.text.y = element_text(family = "Times New Roman", size = 12),  # Remove y-axis text (numbers)
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.title.x = element_text(family = "Times New Roman", size = 12),  # Change x-axis title font
        plot.title = element_text(family = "Times New Roman", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = unit(c(0, 20, 0, 15), "pt"),
        legend.position = "none"
      )
    
    # Add subfigure label to each plot
    plot_with_label <- arrangeGrob(plot, #padding = unit(c(0, 20), "pt"),
                                   top = textGrob(paste0("(", subfig_labels[i], ")"),
                                                  x = unit(0.0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
                                                  gp = gpar(fontsize = 14, fontfamily = "Times New Roman"))
    )
    
    # Store the labeled plot in the list
    plot_list[[i]] <- plot_with_label
  }
  
  # Arrange all 8 plots in a 2x4 grid
  grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)
  #Trace_MCMC.pdf, 12 x 9
}

AlgoResultsABC <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2025-08-03_093559.602748__July25Agg_NormalCDF_ESplus0.05MVR_RF_kappa0.5_range0.5')

SR1$HLP_vals = rep(Inf, length(SR1$loss))
SR3$HLP_vals = rep(Inf, length(SR3$loss))
SR4$HLP_vals = rep(Inf, length(SR4$loss))
for (i in 1:15000){
  if (!is.na(SR1$Omega_sample_phys[1,i])){
    SR1$HLP_vals[i] <- Model$HighLevelPriors(SR1$Omega_sample_phys[,i] %>% 
                                               relist(skeleton=Model$skeleton) %>% 
                                               addTransfParams(),
                                             Model)
  }
  if (!is.na(SR3$Omega_sample_phys[1,i])){
    SR3$HLP_vals[i] <- Model$HighLevelPriors(SR3$Omega_sample_phys[,i] %>% 
                                                   relist(skeleton=Model$skeleton) %>% 
                                                   addTransfParams(),
                                                 Model)
  }
  if (!is.na(SR4$Omega_sample_phys[1,i])){
    SR4$HLP_vals[i] <- Model$HighLevelPriors(SR4$Omega_sample_phys[,i] %>% 
                                                   relist(skeleton=Model$skeleton) %>% 
                                                   addTransfParams(),
                                                 Model)
  }
}

plot_correlated_posteriors_ggplot(AlgoResultsABC, include_priors=T, Omega=NULL,
                                  pairings = rbind(c(1,2), c(3,4), c(5,6), c(8,7), c(8,9),  c(10,11), c(15,12), c(13, 14), c(16,17), c(18,19)), s_finish=NULL, 
                                  AlgoResultsMCMC=list(SR1, SR3, SR4), subfig_title_adj=-0.3, mcmc_warmup = 13000)
#MCMCvsSMC.pdf, 5 x 13
# 7 x 13 for presentation plot

#------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------- Euclidean Distance Plots -------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------

source('RCode/ResultsAnalysisFunctions.R')
EucDist1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-10-07_164759.048116__Sept25_NormalCDF_EuclDist_RFwithTot_range1_Chain1_Np10_backup')
EucDist2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-10-07_164114.99222__Sept25_NormalCDF_EuclDist_RFwithTot_range1_Chain2_Np10_backup')
EucDist3 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-10-08_100845.107303__Sept25_NormalCDF_EuclDist_RFwithTot_range1_Chain3_Np10_backup')

EucDist_combined = EucDist1
EucDist_combined$Omega_sample_phys = EucDist_combined$Omega_sample_phys[,1:15000]
EucDist_combined$Omega_sample_phys[,1:5000] = EucDist2$Omega_sample_phys[,5001:10000]
EucDist_combined$Omega_sample_phys[,5001:10000] = EucDist3$Omega_sample_phys[,5001:10000]
EucDist_combined$Omega_sample_phys[,1:15000] = EucDist_combined$Omega_sample_phys[,sample(1:1500, 1500, replace=F)]
EucDist_combined$Omega_sample_phys[,15000] = NA

df_postpredictive_euc <- create_df_postpredictive_MCMC(EucDist_combined, single_particle=F, Omega=NULL, particle_best=F,
                                                                M=100, output='SampledTotal')

saveRDS(df_postpredictive_euc, 'IIDIPUS_Results/df_postpredictive_euc_Oct25')

p_mort = plot_df_postpredictive(df_postpredictive_euc %>% filter(train_flag=='TEST'), impact_type='mortality')
p_disp = plot_df_postpredictive(df_postpredictive_euc %>% filter(train_flag=='TEST'), impact_type='displacement')

p_mort + theme(legend.position = "none")
p_disp + theme(legend.position = "none")

plot_grid(p_mort + theme(legend.position = "none"), p_disp + theme(legend.position = "none"), 
          labels=c('(a)', '(b)'), label_fontfamily = "Times New Roman", label_fontface = "plain", 
          label_size = 13, label_x = -0.015, label_y = 0.98)

# EucDistPred.pdf, 8 x 4.1

mcmc_trace(EucDist1, EucDist2, EucDist3, var_plot=c(1:19), xlim=c(1,10000))


quants <- apply(df_postpredictive_euc[which(df_postpredictive_euc$impact=='mortality'& df_postpredictive_euc$train_flag=='TEST'),c(grep('observed', names(df_postpredictive_euc)), grep('sampled.', names(df_postpredictive_euc)))], 1, function(x){sample_quant(x)})/(length(grep('sampled.', names(df_postpredictive_euc)))+1)
hist(quants)
mean(quants>0.025 & quants < 0.975)
quants <- apply(df_postpredictive_euc[which(df_postpredictive_euc$impact=='displacement'& df_postpredictive_euc$train_flag=='TEST'),c(grep('observed', names(df_postpredictive_euc)), grep('sampled.', names(df_postpredictive_euc)))], 1, function(x){sample_quant(x)})/(length(grep('sampled.', names(df_postpredictive_euc)))+1)
hist(quants)
mean(quants>0.025 & quants < 0.975)


#---------------------------------------------------------------------------------------------------

traces = plot_mcmc_compare(EucDist1, EucDist2, SR1, SR4, xlimEuc = c(5000, 10000), xlimSR = c(10000, 15000))


plot_mcmc_compare = function(EucDist1, EucDist2, SR1, SR2, xlimEuc=c(8000, 10000), xlimSR=c(13000, 15000)){
  
  # Define subfigure labels (a), (b), (c), ...
  subfig_labels <- letters[1:8]
  plot_list <- list()
  var_plot <- c(8,9,10, 7)
  # Loop through the variables and generate each plot
  for (i in seq_along(var_plot)){#seq_along(vuln_var)) {
    #var <- vuln_var[i]
    var <- var_plot[i]
    
    post_samples <- data.frame(
      x=(xlimSR[1]:xlimSR[2])-xlimSR[1],
      Euc1_param <- EucDist1$Omega_sample_phys[var,xlimEuc[1]:xlimEuc[2]],
      Euc2_param <- EucDist2$Omega_sample_phys[var,xlimEuc[1]:xlimEuc[2]],
      SR1_param <- SR1$Omega_sample_phys[var,xlimSR[1]:xlimSR[2]],
      SR2_param <-SR2$Omega_sample_phys[var,xlimSR[1]:xlimSR[2]])
    
    post_samples_long <- post_samples %>%
      pivot_longer(cols = -x, names_to = "Parameter", values_to = "Value")
    
    # Plot the lines
    plot <- ggplot(data = post_samples_long, aes(x = x, y = Value, color = Parameter)) +
      labs(y = get_greek_titles(names(unlist(Model$skeleton))[var]),
           x= 'Post warmup iteration') +
      scale_y_continuous(limits = c(Model$par_lb[var], Model$par_ub[var]), expand = c(0.01, 0.01)) +
      geom_line(alpha=0.8) +
      theme_minimal() +
      scale_color_viridis(discrete = TRUE) +
      theme(
        axis.title.y = element_text(family = "Times New Roman", size = 12),  # Change y-axis title font
        axis.text.x = element_text(family = "Times New Roman", size = 12),   # Change x-axis text font
        axis.text.y = element_text(family = "Times New Roman", size = 12),  # Remove y-axis text (numbers)
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.title.x = element_text(family = "Times New Roman", size = 12),  # Change x-axis title font
        plot.title = element_text(family = "Times New Roman", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = unit(c(0, 20, 0, 15), "pt"),
        legend.position = "none"
      )
    
    # Add subfigure label to each plot
    plot_with_label <- arrangeGrob(plot, #padding = unit(c(0, 20), "pt"),
                                   top = textGrob(paste0("(", subfig_labels[i], ")"),
                                                  x = unit(0.0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
                                                  gp = gpar(fontsize = 14, fontfamily = "Times New Roman"))
    )
    
    # Store the labeled plot in the list
    plot_list[[i]] <- plot_with_label
  }
  
  # Arrange all 8 plots in a 2x4 grid
  grid.arrange(grobs = plot_list, ncol = 2, nrow = 2)
  #Trace_EuclDist.pdf, 12 x 9 (or 13 x 3 for single row)
}


df_postpredictive_sampled_best <- create_df_postpredictive_MCMC(EucDist1, single_particle=F, Omega=NULL, particle_best=F,
                                                                M=100, output='SampledTotal')

#saveRDS(df_postpredictive_sampled_best, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/EuclDistAug7')
df_postpredictive_sampled_bestEuc <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/EuclDistAug7')

df_postpredictive_sampled_bestSR = df_postpredictive_sampled
df_postpredictive_sampled_bestSR <- readRDS('fit_MVR_August6')

mort_test <- plot_df_postpredictive_compare(df_postpredictive_sampled_bestEuc %>% filter(train_flag=='TEST'),
                                            df_postpredictive_sampled_bestSR %>% filter(train_flag=='TEST'), 'mortality')  + guides(color="none") 
disp_test <- plot_df_postpredictive_compare(df_postpredictive_sampled_bestEuc %>% filter(train_flag=='TEST'),
                                            df_postpredictive_sampled_bestSR %>% filter(train_flag=='TEST'), 'displacement')+ guides(color="none") 
buildDam_test <- plot_df_postpredictive_compare(df_postpredictive_sampled_bestEuc %>% filter(train_flag=='TEST'),
                                                df_postpredictive_sampled_bestSR %>% filter(train_flag=='TEST'), 'buildDam')+ guides(color="none") + ylab('Sampled building damage') + xlab('Observed building damage')

p1 <- arrangeGrob(mort_test + theme(plot.margin = unit(c(0, 0.1, 1, 1), "lines")), top = textGrob("(a)", x = unit(0.025, "npc"), y = unit(0, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p2 <- arrangeGrob(disp_test + theme(plot.margin = unit(c(0, 0.1, 1, 1), "lines")), top = textGrob("(b)", x = unit(0.025, "npc"), y = unit(0, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p3 <- arrangeGrob(buildDam_test + theme(plot.margin = unit(c(0, 0.1, 1, 1), "lines")), top = textGrob("(c)", x = unit(0.025, "npc"), y = unit(0, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))

grid.arrange(p1, p2, p3, ncol=3)
#EucDist_PostPred.pdf, 17 x 5.5






















################################################################################################################################################








AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-05-20_053958.556514__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain')
AlgoResults_Logistic <- AlgoResults

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-04-27_014135.935892__Apr25Agg_Lognormal_ESplus0.05BD_RFwithTot_backup')
AlgoResults_Logistic <- AlgoResults

get_tangent_7.5 <- function(mu, sigma){
  x0 <- 7.5
  z <- (x0 - mu) / sigma
  
  g <- log(pnorm(z))                # log CDF value at 7.5
  gp <- dnorm(z) / (sigma * pnorm(z))  # derivative at 7.5
  
  a <- g - gp * (x0 - 7)  # intercept at (x - 7)
  b <- gp                # slope
  return(c(a,b))
}


D_MortDisp_calc <- function(Damage, Omega, stoch_event=rbind(0,0)){
  D_Mort_and_Disp <- pnorm(Damage + stoch_event[2,], mean = Omega$Lambda1$loc, sd = Omega$Lambda1$scale)
  D_Mort <- pmin(exp(Omega$Lambda2$nu + Omega$Lambda2$kappa * (Damage + stoch_event[1,] - 7)), 1-1e-6)
  D_Disp <- D_Mort_and_Disp - D_Mort
  D_Disp <- ifelse(D_Disp<0, 10^(-60), D_Disp)
  return(rbind(D_Mort, D_Disp))
}

fit_log_linear_mort_to_pnorm <- function(omega_pnorm, damage_range = seq(6, 9, 0.1)) {
  # True probabilities from pnorm model
  D_Mort_and_Disp <- pnorm(damage_range, mean = omega_pnorm[1], sd = omega_pnorm[2])
  D_Mort <- pnorm(damage_range, mean = omega_pnorm[3], sd = omega_pnorm[4])
  D_Disp <- pmax(D_Mort_and_Disp - D_Mort, 1e-12)
  D_True <- rbind(D_Mort, D_Disp)
  
  # Loss function: linear-on-log mortality
  loss_fn <- function(par) {
    a <- par[1]  # intercept for log mortality
    b <- par[2]  # slope for log mortality
    nu_disp <- par[3]
    kappa_disp <- par[4]
    
    # Log-linear mortality
    P_mort <- pmin(exp(a + b * (damage_range - 7)), 1-1e-6)
    
    # Logistic model for displacement (residual probability from non-mortality)
    P_dispmort = pnorm(damage_range, nu_disp, kappa_disp)
    # 1 / (1 + exp(-(nu_disp + kappa_disp * (damage_range - 7))))
    P_disp = pmax(P_dispmort-P_mort, 1e-12)
    
    D_hat <- rbind(P_mort, P_disp)
    loss <- sum((log(D_hat) - log(D_True))^2)
    return(loss)
  }
  
  # Initial guess
  init <- c(-8, 1.5, 7.5, 2)
  
  # Optimization
  fit <- optim(init, loss_fn, method = "BFGS")
  
  # Return parameters
  list(
    mortality = list(a = fit$par[1], b = fit$par[2]),
    displacement = list(nu = fit$par[3], kappa = fit$par[4]),
    convergence = fit$convergence,
    loss = fit$value
  )
}
  
for (i in 1:NCOL(AlgoResults_Logistic$Omega_sample_phys)){
  print(i)
  fitted_vals = fit_log_linear_mort_to_pnorm(AlgoResults$Omega_sample_phys[,i])
  
    # Get tangents and update logistic Omega to match slope & intercept at x = 7.5
  tang_params_disp <- get_tangent_7.5(AlgoResults$Omega_sample_phys[1,i], AlgoResults$Omega_sample_phys[2,i])
  AlgoResults_Logistic$Omega_sample_phys[1,i] <- fitted_vals$displacement$nu #-4 #tang_params_disp[1]
  AlgoResults_Logistic$Omega_sample_phys[2,i] <- fitted_vals$displacement$kappa #4.5 #tang_params_disp[2]
  
  tang_params_mort <- get_tangent_7.5(AlgoResults$Omega_sample_phys[3,i], AlgoResults$Omega_sample_phys[4,i])
  AlgoResults_Logistic$Omega_sample_phys[3,i] <- fitted_vals$mortality$a#tang_params_mort[1]
  AlgoResults_Logistic$Omega_sample_phys[4,i] <- fitted_vals$mortality$b# tang_params_mort[2]
  
  #tang_params_bd <- get_tangent_7.5(AlgoResults$Omega_sample_phys[5,i], AlgoResults$Omega_sample_phys[6,i])
  #AlgoResults_Logistic$Omega_sample_phys[5,i] <- tang_params_bd[1]
  #AlgoResults_Logistic$Omega_sample_phys[6,i] <- tang_params_bd[2]
    
    # Plot for visual comparison
  # curve_x <- seq(4.5, 10, 0.1)
  # plot(curve_x, log(pnorm(curve_x, AlgoResults$Omega_sample_phys[3,i], AlgoResults$Omega_sample_phys[4,i])),
  #      type = "l", col = "black", ylab = "log prob", main = paste("Sample", i), ylim=c(-20,0))
  # lines(curve_x, log(pnorm(curve_x, AlgoResults$Omega_sample_phys[1,i], AlgoResults$Omega_sample_phys[2,i])),
  #       type = "l", col = "black")
  # # Compute Mort/Disp from logistic softmax
  # MortDisp <- D_MortDisp_calc(curve_x, AlgoResults_Logistic$Omega_sample_phys[,i] %>% 
  #                               relist(skeleton = Model$skeleton) %>%
  #                               addTransfParams())
  # lines(curve_x, log(MortDisp[1,]), col = "red")
  # lines(curve_x, log(MortDisp[2,]), col = "blue")
}


hlps <- c()
for (i in 1:NCOL(AlgoResults_Logistic$Omega_sample_phys)){
  hlps <- c(hlps, Model$HighLevelPriors(AlgoResults_Logistic$Omega_sample_phys[,i] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), Model))
}
plot(hlps)

saveRDS(AlgoResults_Logistic, paste0(dir, 'IIDIPUS_Results/Logistic_Artificial'))

Logistic_Artificial = readRDS(paste0(dir, 'IIDIPUS_Results/Logistic_Artificial'))
Logistic_Artificial$Omega_sample_phys[8,]= Logistic_Artificial$Omega_sample_phys[8,] / 3 + 0.5 
saveRDS(Logistic_Artificial,  'IIDIPUS_Results/Logistic_Artificial2')


AlgoResults_Logistic = readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-04-24_181652.03706__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup')
hlps <- c()
for (i in 1:NCOL(AlgoResults_Logistic$Omega_sample_phys)){
  hlps <- c(hlps, Model$HighLevelPriors(AlgoResults_Logistic$Omega_sample_phys[,i] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), Model))
}


moveSubFiles <- function(filepath='IIDIPUS_Input_Alternatives/Sept25/ODDobjects'){
  all_files <- list.files(path = filepath, recursive = TRUE, full.names = TRUE)
  file_paths <- all_files[file.info(all_files)$isdir == FALSE]
  
  # Move each file to the current directory
  file.rename(from = file_paths,
              to = file.path(filepath, basename(file_paths)))
}


moveTestData <- function(folder_in='IIDIPUS_Input_Alternatives/Sept25'){
  
    if (!dir.exists(paste0(dir, folder_in, '/ODDobjects/Test/'))) {
      dir.create(paste0(dir, folder_in, '/ODDobjects/Test/'))
    } 
    if (!dir.exists(paste0(dir, folder_in, '/ODDobjects/Train/'))){
      dir.create(paste0(dir, folder_in, '/ODDobjects/Train/'))
    }
    
    ODD_folderall<-paste0(dir, folder_in, '/ODDobjects/')
    ODD_foldertest<-paste0(dir, folder_in, '/ODDobjects/Test/')
    ODD_foldertrain<-paste0(dir, folder_in, '/ODDobjects/Train/')
    
    ufiles<-list.files(path=ODD_folderall,pattern=Model$haz,recursive = T,ignore.case = T)
    #set.seed(1)
    #ufiles <- ufiles[sample(1:length(ufiles), length(ufiles), replace=F)]
    total_mortalities <- c()
    #sort by mortality
    for (file in ufiles){
      ODD <- readODD(paste0(ODD_folderall, file))
      #total_mortalities <- c(total_mortalities, max(values(ODD[[grep('hazMean', names(ODD))]])[values(!is.na(ODD$Population))], na.rm=T))
      #[values(!is.na(ODD$Population) & ODD$Population > 100),], na.rm=T))
      #total_mortalities <- c(total_mortalities, max(values(ODD[[grep('hazMean', names(ODD))]])[values(!is.na(ODD$Population))], na.rm=T))#[values(!is.na(ODD$Population) & ODD$Population > 100),], na.rm=T))
      polygon_names <- unlist(lapply(ODD@polygons[ODD@impact$polygon], function(x) x$name))
      if (any(tolower(polygon_names[which(ODD@impact$impact=='mortality')]) %in% c('tot', 'total'))){
        nonmatch <- which(!tolower(polygon_names[which(ODD@impact$impact=='mortality')]) %in% c('tot', 'total'))
        if (length(nonmatch)>0){
          ODD@impact <- ODD@impact[-which(ODD@impact$impact=='mortality')[nonmatch],] # in the case of total and subnational data, remove the subnational
        }
      }
      total_mortality <- sum(ODD@impact$observed[which(ODD@impact$impact=='mortality')])
      total_mortalities <- c(total_mortalities, total_mortality)
    }
    ufiles <- ufiles[order(total_mortalities, decreasing=T)] # sort by date

    for (i in 1:length(ufiles)){
      file <- ufiles[i]
      #if (i < 2/3 * length(ufiles)){ #next
      if (i %%3 != 1){
        options(warn = 2)
        file.copy(from = paste0(ODD_folderall, file),
                  to = paste0(ODD_foldertrain, file))
        file.remove(from = paste0(ODD_folderall, file))
        options(warn = 1)
        next
      }
      options(warn = 2)
      file.copy(from = paste0(ODD_folderall, file),
                to = paste0(ODD_foldertest, file))
      file.remove(from = paste0(ODD_folderall, file))
      options(warn = 1)
    }
}

df_postpredictive_sampled_best = readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/fit_MVR_August6')
df_postpredictive_sampled_best %<>% add_row( event_id = 171, 
                                            iso3='MAR',
                                            polygon=0,
                                            impact='mortality', 
                                            observed = 2933,
                                            obs_id = 1000,
                                            crps=100,
                                            train_flag='TRAIN', 
                                            train_flag2='TRAIN')
df_postpredictive_sampled_best[nrow(df_postpredictive_sampled_best), grep('sampled', names(df_postpredictive_sampled_best))] = rnorm(100, 1000, 100)
df_postpredictive_sampled_best$train_flag2 <- 'TRAIN'
ufiles<-list.files(path=paste0(dir, 'IIDIPUS_Input_Alternatives/Sept25/ODDobjects/Test'),pattern=Model$haz,recursive = T,ignore.case = T)
test_id = as.numeric(sub(".*_", "", ufiles))
for (i in 1:length(df_postpredictive_sampled_best$train_flag2)){
  if (df_postpredictive_sampled_best$event_id[i] %in% test_id){
    df_postpredictive_sampled_best$train_flag2[i] <- 'TEST'
  }
}

df_postpredictive_sampled_best$train_flag=df_postpredictive_sampled_best$train_flag2
mort_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'mortality' )  + guides(color="none") 
mort_all <- plot_df_postpredictive(df_postpredictive_sampled_best,'mortality' )  + guides(color="none") 
mort_all
grid.arrange(mort_test, mort_all, ncol=2)
disp_test <- plot_df_postpredictive(df_postpredictive_sampled_best,'displacement')  + guides(color="none")
buildDam_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'buildDam')  + guides(color="none") #+ 
  #scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma, limits=c(50,100000))


p1 <- arrangeGrob(mort_test, top = textGrob("(a)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p2 <- arrangeGrob(disp_test, top = textGrob("(b)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p3 <- arrangeGrob(buildDam_test, top = textGrob("(c)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))

grid.arrange(p1,p2,p3, nrow=1)

mort_test
mort_all
disp_test
buildDam_test








plot_df_postpredictive(df_postpredictive_sampled_best,'mortality', label_points=T)


plot(seq(4.5, 10,0.01), log(pnorm(seq(4.5, 10,0.01), AlgoResults$Omega_sample_phys[1,3,120],  AlgoResults$Omega_sample_phys[1,4,120])), ylim=c(-20,-1))
for (i in 1:100){
  lines(seq(4.5, 10,0.01), log(pnorm(seq(4.5, 10,0.01), AlgoResults$Omega_sample_phys[i,3,120],  AlgoResults$Omega_sample_phys[i,4,120])))
  lines(seq(4.5, 10,0.01), 2.6*seq(4.5, 10,0.01)-27.5, col='red')
}

mort_df = df_postpredictive_sampled_best %>% filter(impact=='mortality')
mort_df$sampled.median= apply(mort_df[,grep('sampled.', names(mort_df))], 1, median)
mort_df$obs_quant = apply(cbind(mort_df$observed, mort_df[,grep('sampled.', names(mort_df))]), 1, function(x) rank(x, ties.method='random')[1])
hist(mort_df$obs_quant[mort_df$sampled.median>100])

plot(log(mort_df$observed+10), log(mort_df$sampled.median+10), col=ifelse(mort_df$train_flag=='TRAIN', 'blue', 'red'))
abline(0, 1)


plot(log(mort_df$sampled.median+10), mort_df$obs_quant, col=ifelse(mort_df$train_flag=='TRAIN', 'blue', 'red'))


ODDy <-readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Apr25Agg/ODDobjects/Test/EQ20220622AFG_165')
ODDy2 <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Apr25Agg/ODDobjects/Train/EQ20210213JPN_155')
plot(values(ODD1$hazMean1), log(values(ODD1$Population)), ylim=c(2, 15))
points(values(ODD2$hazMean1), values(log(ODD2$Population)), col='red')



AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2025-06-11_114651.469474__June25Agg_LogLinMort_ESplus0.05BD_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain3')

source('RCode/Model Variations/ODDobj_NoTot.R')
source('RCode/Model Variations/Model_Logistic2.R')

AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/Apr25Agg/'
df_postpredictive_sampled_best <- create_df_postpredictive_MCMC(AlgoResults, single_particle=F, 
                                                           M=50, output='SampledTotal')

mort_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag == 'TEST'),'mortality' )  + guides(color="none") 
mort_test

disp_test <- plot_df_postpredictive(df_postpredictive_sampled_best%>% filter(train_flag == 'TEST'),'displacement', label_points=T)  + guides(color="none")
disp_test

bd_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'buildDam')  + guides(color="none")
bd_test


#--------------------------------------------------------------------------------------------------
#------------------------------ Model Comparison --------------------------------------------------
#--------------------------------------------------------------------------------------------------
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

#Normal CDF (misnamed)
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-06-17_151849.678089__Apr25Agg_LogNormal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_withHLPriors_flexHLP')

AlgoResultsLogLin <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2025-06-12_115400.392519__Apr25Agg_LogLinMort_ESplus0.05BD_RFwithTot_range.5_kappa.5_withHLPriors_flexHLP_backup')
plot(AlgoResultsLogLin$loss)

source('RCode/Model Variations/Model_Logistic2.R')
source('RCode/ODDobj.R')

AlgoParams$input_folder = 'IIDIPUS_Input_Alternatives/Apr25Agg/'
xx1 = SampleImpact(dir, Model, proposed = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
dist1 = CalcDist(xx1, AlgoParams)
print(dist1)

xx2 = SampleImpact(dir, Model, proposed = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
dist2 = CalcDist(xx2, AlgoParams)
print(dist2)

xx3 = SampleImpact(dir, Model, proposed = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
dist3 = CalcDist(xx3, AlgoParams)
print(dist3)

xx4 = SampleImpact(dir, Model, proposed = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
dist4 = CalcDist(xx4, AlgoParams)
print(dist4)

yy1 = SampleImpact(dir, Model, proposed = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
disty1 = CalcDist(yy1, AlgoParams)
print(disty1)

yy2 = SampleImpact(dir, Model, proposed = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
disty2 = CalcDist(yy2, AlgoParams)
print(disty2)

Omega_alt = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton)
Omega_alt$Lambda2$nu = -12.3
Omega_alt$Lambda2$kappa = 3.45
plot(seq(4.5, 12, 0.01), log( 1/(1+exp(-(-12.02085 + 2.447839 * (seq(4.5, 12, 0.01)- 7))))))
lines(seq(4.5, 12, 0.01), log( 1/(1+exp(-(-12.3 + 3.45 * (seq(4.5, 12, 0.01)- 7))))))

zz1 = SampleImpact(dir, Model, proposed = Omega_alt %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
distz1 = CalcDist(zz1, AlgoParams)
print(distz1)


Omega_alt2 = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton)
#Omega_alt2$Lambda2$nu = -12.3
#Omega_alt2$Lambda2$kappa = 3.45
plot(seq(4.5, 12, 0.01), log( 1/(1+exp(-(-12.02085 + 2.447839 * (seq(4.5, 12, 0.01)- 7))))))
lines(seq(4.5, 12, 0.01), log( 1/(1+exp(-(Omega_alt2$Lambda2$nu + Omega_alt2$Lambda2$kappa * (seq(4.5, 12, 0.01)- 7))))))
lines(seq(4.5, 12, 0.01), log( 1/(1+exp(-(Omega_alt2$Lambda2$nu + Omega_alt2$Lambda2$kappa * (seq(4.5, 12, 0.01)+1- 7))))), col='red')
lines(seq(4.5, 12, 0.01), log( 1/(1+exp(-(Omega_alt2$Lambda2$nu + Omega_alt2$Lambda2$kappa * (seq(4.5, 12, 0.01)-1- 7))))), col='red')
a = 0.001
b = 6.93
Damage_tranf = exp(a*seq(4.5, 12, 0.01)+b)
Damage_transf_perturbed = Damage_tranf + 1
Damage_transf_perturbed2 = Damage_tranf - 1
Damage_untransf = (log(Damage_transf_perturbed)-b)/a
Damage_untransf2 = (log(Damage_transf_perturbed2)-b)/a
lines(seq(4.5, 12, 0.01), log( 1/(1+exp(-(Omega_alt2$Lambda2$nu + Omega_alt2$Lambda2$kappa * (Damage_untransf- 7))))), col='blue')
lines(seq(4.5, 12, 0.01), log( 1/(1+exp(-(Omega_alt2$Lambda2$nu + Omega_alt2$Lambda2$kappa * (Damage_untransf2 - 7))))), col='blue')


aa1 = SampleImpact(dir, Model, proposed = Omega_alt2 %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
dista1 = CalcDist(aa1, AlgoParams)
print(dista1)



mort_test <- plot_df_postpredictive(flattenImpactSample2(zz1$poly) %>% add_column(train_flag='TRAIN'),'mortality')  + guides(color="none")
mort_test2 <- plot_df_postpredictive(flattenImpactSample2(xx1$poly) %>% add_column(train_flag='TRAIN'),'mortality')  + guides(color="none")
mort_test3 <- plot_df_postpredictive(flattenImpactSample2(aa1$poly) %>% add_column(train_flag='TRAIN'),'mortality')  + guides(color="none") 
grid.arrange(mort_test2, mort_test3, ncol=2)

zz2 = SampleImpact(dir, Model, proposed = Omega_alt %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
distz2 = CalcDist(zz2, AlgoParams)
print(distz2)

xx3 = SampleImpact(dir, Model, proposed = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
dist3 = CalcDist(xx3, AlgoParams)
print(dist3)

xx4 = SampleImpact(dir, Model, proposed = AlgoResultsLogLin$Omega_sample_phys[,940] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), AlgoParams, dat = 'Train', output='SampledAgg')
dist4 = CalcDist(xx4, AlgoParams)
print(dist4)







#AlgoResults <- readRDS('IIDIPUS_Results/mcmc_2025-05-05_225121.802448__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5')
#AlgoResults2 <- readRDS('IIDIPUS_Results/HPC/mcmc_2025-05-15_211354.026585__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain2')
#AlgoResults3 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2025-06-12_115400.392519__Apr25Agg_LogLinMort_ESplus0.05BD_RFwithTot_range.5_kappa.5_withHLPriors_flexHLP_backup')

# plot(AlgoResults2$loss, ylim=c(3, 5))
# lines(AlgoResults$loss)
# 
# param=13
# plot(AlgoResults2$Omega_sample_phys[param,], ylim=c(Model$par_lb[param],Model$par_ub[param] ))
# lines(AlgoResults$Omega_sample_phys[param,], col='red')
# 
# 
# plot(0, 0, xlim=c(4.5, 12), ylim=c(-15, 0))
# for (i in 1:1300){
#   lines(seq(4.5, 12, 0.01), log(pnorm(seq(4.5, 12, 0.01), AlgoResults$Omega_sample_phys[3,i], AlgoResults$Omega_sample_phys[4,i])))
#   lines(seq(4.5, 12, 0.01), log(pnorm(seq(4.5, 12, 0.01), AlgoResults2$Omega_sample_phys[3,i+3000], AlgoResults2$Omega_sample_phys[4,i+3000])), col='blue')
#   lines(seq(4.5, 12, 0.01), log(1/(1+exp(-(AlgoResults3$Omega_sample_phys[3,i]+AlgoResults3$Omega_sample_phys[4,i]*(seq(4.5, 12, 0.01)-7))))), col='red')
#   # if ( AlgoResults$Omega_sample_phys[3,i] > 12.5){
#   #   lines(seq(4.5, 12, 0.01), log(pnorm(seq(4.5, 12, 0.01), AlgoResults$Omega_sample_phys[3,i], AlgoResults$Omega_sample_phys[4,i])), col='blue')
#   # } else {
#   #   lines(seq(4.5, 12, 0.01), log(pnorm(seq(4.5, 12, 0.01), AlgoResults$Omega_sample_phys[3,i], AlgoResults$Omega_sample_phys[4,i])))
#   # }
# }
# lines(seq(4.5, 12,0.01), -10 + 3 * (seq(4.5, 12,0.01)-7))
# abline(v=8.7)

#--------------------------------------------------------------------------------------------------
#------------------------------ Euclidean Distance Comparison -------------------------------------
#--------------------------------------------------------------------------------------------------

# MCMC Compare Euclidean distance to Scoring Rule Posterior:
EucDist1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-05-24_022553.767998__Apr25Agg_NormalCDF_EuclDist_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain1')
EucDist2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-05-24_023053.892028__Apr25Agg_NormalCDF_EuclDist_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain2')
SR1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-05-20_053958.556514__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain')
SR2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-05-20_054516.483032__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain2')

plot_mcmc_compare = function(EucDist1, EucDist2, SR1, SR2, xlim=c(5000, 10000)){
  
  # Define subfigure labels (a), (b), (c), ...
  subfig_labels <- letters[1:8]
  plot_list <- list()
  var_plot <- c(7, 9, 10)
  # Loop through the variables and generate each plot
  for (i in seq_along(var_plot)){#seq_along(vuln_var)) {
    #var <- vuln_var[i]
    var <- var_plot[i]
    
    post_samples <- data.frame(
      x=(xlim[1]:xlim[2])-xlim[1],
      Euc1_param <- EucDist1$Omega_sample_phys[var,xlim[1]:xlim[2]],
      Euc2_param <- EucDist2$Omega_sample_phys[var,xlim[1]:xlim[2]],
      SR1_param <- SR1$Omega_sample_phys[var,xlim[1]:xlim[2]],
      SR2_param <-SR2$Omega_sample_phys[var,xlim[1]:xlim[2]])
    
    post_samples_long <- post_samples %>%
      pivot_longer(cols = -x, names_to = "Parameter", values_to = "Value")
    
    # Plot the lines
    plot <- ggplot(data = post_samples_long, aes(x = x, y = Value, color = Parameter)) +
      labs(y = get_greek_titles(names(unlist(Model$skeleton))[var]),
           x= 'Post warmup iteration') +
      scale_y_continuous(limits = c(Model$par_lb[var], Model$par_ub[var]), expand = c(0.01, 0.01)) +
      geom_line(alpha=0.8) +
      theme_minimal() +
      scale_color_viridis(discrete = TRUE) +
      theme(
        axis.title.y = element_text(family = "Times New Roman", size = 12),  # Change y-axis title font
        axis.text.x = element_text(family = "Times New Roman", size = 12),   # Change x-axis text font
        axis.text.y = element_text(family = "Times New Roman", size = 12),  # Remove y-axis text (numbers)
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.title.x = element_text(family = "Times New Roman", size = 12),  # Change x-axis title font
        plot.title = element_text(family = "Times New Roman", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = unit(c(0, 20, 0, 15), "pt"),
        legend.position = "none"
      )
    
    # Add subfigure label to each plot
    plot_with_label <- arrangeGrob(plot, #padding = unit(c(0, 20), "pt"),
                                   top = textGrob(paste0("(", subfig_labels[i], ")"),
                                                  x = unit(0.0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
                                                  gp = gpar(fontsize = 14, fontfamily = "Times New Roman"))
    )
    
    # Store the labeled plot in the list
    plot_list[[i]] <- plot_with_label
  }
  
  # Arrange all 8 plots in a 2x4 grid
  grid.arrange(grobs = plot_list, ncol = 3, nrow = 1)
  #Trace_EuclDist.pdf, 10 x 3
}

#Postpredictive comparison
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/Apr25Agg/'
source('RCode/Model Variations/ODDobj_NoTot.R')

df_postpredictive_sampled_best <- create_df_postpredictive_MCMC(EucDist1, single_particle=F, Omega=NULL, particle_best=F,
                                                                M=100, output='SampledTotal')
saveRDS(df_postpredictive_sampled_best, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/EuclDistJune2025')
df_postpredictive_sampled_bestEuc <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/EuclDistJune2025')

df_postpredictive_sampled_bestSR <- create_df_postpredictive_MCMC(SR1, single_particle=F, Omega=NULL, particle_best=F,
                                                                M=100, output='SampledTotal')
saveRDS(df_postpredictive_sampled_bestSR, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SRJune2025')
df_postpredictive_sampled_bestSR <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SRJune2025')

mort_test <- plot_df_postpredictive_compare(df_postpredictive_sampled_bestEuc %>% filter(train_flag=='TEST'),
                                            df_postpredictive_sampled_bestSR %>% filter(train_flag=='TEST'), 'mortality')  + guides(color="none") 
disp_test <- plot_df_postpredictive_compare(df_postpredictive_sampled_bestEuc %>% filter(train_flag=='TEST'),
                                            df_postpredictive_sampled_bestSR %>% filter(train_flag=='TEST'), 'displacement')+ guides(color="none") 
# buildDam_test <- plot_df_postpredictive_compare(df_postpredictive_sampled_bestEuc %>% filter(train_flag=='TEST'),
#                                                 df_postpredictive_sampled_bestSR %>% filter(train_flag=='TEST'), 'buildDam')

p1 <- arrangeGrob(mort_test, top = textGrob("(a)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p2 <- arrangeGrob(disp_test, top = textGrob("(b)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
#p3 <- arrangeGrob(buildDam_test, top = textGrob("(c)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))

grid.arrange(p1, p2, ncol=2)
#EucDist_PostPred.pdf, 12 x 6

#---------------------------------------------------------------------------------------------
#------------  Compare Fit with and without rank histogram tests -----------------------------
#---------------------------------------------------------------------------------------------

source('RCode/Model Variations/ODDobj_NoTot.R')

AlgoResults_bd <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-05-20_053958.556514__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain')
Omega_bd <- AlgoResults_bd$Omega_sample_phys[,5001] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Apr25Agg/"
AlgoResults_bd$input_folder <- "IIDIPUS_Input_Alternatives/Apr25Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
set.seed(123)
impact_sample_bd <- SampleImpact(dir, Model, Omega_bd %>% addTransfParams(), AlgoParams, dat='Train')
#impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train', output='SampledTotal')
#saveRDS(impact_sample_rf, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
df_bd <- flattenImpactSample(impact_sample_bd)
ranks_bdpth <- get_banddepth_rank(df_bd, log=T)
hist(ranks_bdpth, main='Banddepth Rank Scores')



AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-05-24_023914.521443__Apr25Agg_NormalCDF_ESOnly_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain1')
Omega <- AlgoResults$Omega_sample_phys[,4001] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Apr25Agg/"
AlgoResults$input_folder <- "IIDIPUS_Input_Alternatives/Apr25Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams, dat='Train')
#impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train', output='SampledTotal')
#saveRDS(impact_sample_rf, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
df_low <- flattenImpactSample(impact_sample)
ranks_bdpth_og <- get_banddepth_rank(df_low, log=T)
hist(ranks_bdpth_og, main='Banddepth Rank Scores')


par(family = "Times")
par(mfrow=c(1,2))
par(mar = c(4.5, 4, 0, 2)) 
hist(ranks_bdpth, col = "#c3b1c7",  # light purple fill
     border = "#440154",
     xlab='BD Ranks (loss function with rank testing)', freq=T, ylab='Frequency', main='',cex.main = 1.2, font.main = 1,  ylim=c(0,32))
mtext("(a)", side = 3,  line = -1, adj = -0.2,  font = 1, cex = 1.2)

hist(ranks_bdpth_og, col = "#c3b1c7",  # light purple fill
     border = "#440154", xlab='BD Ranks (loss function without rank testing)', freq=T, ylab='Frequency', main='',font.main = 1, ylim=c(0,32))
mtext("(b)", side = 3, line = -1, adj = -0.2, font = 1, cex = 1.2)


