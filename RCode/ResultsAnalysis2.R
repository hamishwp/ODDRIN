
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#--------------------------------------------------------------------------------------------------
#------------------------------ Simulated Data Results Analysis -----------------------------------
#--------------------------------------------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(patchwork)

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2025-08-06_112005.101225_SimInput9_NormalCDF_ESplus0.05MVR_RF_kappa0.5_range0.5')
#plot(AlgoResults$tolerancestore[2:length(AlgoResults$tolerancestore)]/AlgoResults$tolerancestore[1:(length(AlgoResults$tolerancestore)-1)], ylim=c(0.99,1))
mean(AlgoResults$d[,1,1])
mean(AlgoResults$d[,1,250])

AlgoResults$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9/'
AlgoParams$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput9/'

plot_correlated_posteriors_ggplot(AlgoResults,  pairings = rbind(c(1,2), c(3,4), c(5,6), c(8,9), c(10,11), c(8,7), c(15,12), c(13, 14), c(16,17), c(18,19)), 
                                  Omega = list(Lambda1 = list(nu=7.5, kappa=0.6),
                                               Lambda2 = list(nu=11.4, kappa=0.775), #list(nu=10.65, kappa=1.5), #
                                               Lambda3 = list(nu=7.9, kappa=0.68),
                                               #Lambda4 = list(nu=9.9, kappa=1.6),
                                               #theta= list(theta1=0.6),
                                               eps=list(local=1.2, hazard_mort=0.4, hazard_disp=0.65, hazard_bd=0.5, hazard_cor=0.1),
                                               #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                                               vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15)), subfig_title_adj=-0.8)
#SimPosteriors2.pdf, 5 x 12 inches



df_postpredictive_sampled_true <- create_df_postpredictive(AlgoResults, single_particle=T, 
                                                           Omega = list(Lambda1 = list(nu=7.5, kappa=0.6),
                                                                        Lambda2 = list(nu=11.4, kappa=0.775), #list(nu=10.65, kappa=1.5), #
                                                                        Lambda3 = list(nu=7.9, kappa=0.68),
                                                                        #Lambda4 = list(nu=9.9, kappa=1.6),
                                                                        #theta= list(theta1=0.6),
                                                                        eps=list(local=1.2, hazard_mort=0.4, hazard_disp=0.65, hazard_bd=0.5, hazard_cor=0.1),
                                                                        #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                                                                        vuln_coeff = list(PDens=0, SHDI=-0.15, GNIc=-0.03, Vs30=-0.03, EQFreq=-0.05, FirstHaz=0.05, Night=0, FirstHaz.Night=0.15)),
                                                           M=100, output='SampledTotal')

df_postpredictive_sampled <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=100, output='SampledTotal')

saveRDS(df_postpredictive_sampled_true, 'SimPostPredictiveAug2025')
saveRDS(df_postpredictive_sampled, 'SimTruePredictiveAug2025')

df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/SimPostPredictiveAug2025')
df_postpredictive_sampled_true <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/SimTruePredictiveAug2025')

plot_df_postpredictive_compare(df_postpredictive_sampled_true %>% filter(train_flag=='TEST'),
                               df_postpredictive_sampled %>% filter(train_flag=='TEST'), 'mortality')

df_postpredictive_sampled_best %<>% filter(observed>0)
df_postpredictive_sampled_true %<>% filter(observed>0)
plot_df_postpredictive_compare(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),
                               df_postpredictive_sampled_true %>% filter(train_flag=='TEST'), 'mortality')
#PostPredSim.pdf, 5 x 9 inches

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#--------------------------------------------------------------------------------------------------
#------------------------------ Real Data Results Analysis ----------------------------------------
#--------------------------------------------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2025-08-03_093559.602748__July25Agg_NormalCDF_ESplus0.05MVR_RF_kappa0.5_range0.5')
#points(AlgoResults$tolerancestore[2:length(AlgoResults$tolerancestore)]/AlgoResults$tolerancestore[1:(length(AlgoResults$tolerancestore)-1)], ylim=c(0.99,1), col='red')


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

plot_correlated_posteriors_ggplot(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(8,9), c(10,11), c(8, 7)))
#RealPosteriors.pdf, 3 x 15 inches

library(grid)
plot_vuln_posteriors(AlgoResults)
#VulnPosteriors.pdf, 6.5 x 12 inches

#AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-13_234229.627104__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain2_backup')
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/July25Agg/'
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=100, output='SampledTotal')
saveRDS(df_postpredictive_sampled_best, 'fit_MVR_August6')
#saveRDS(df_postpredictive_sampled_best, 'SampledTotal8June')

df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=100, output='SampledAgg')
#saveRDS(df_postpredictive_sampled_best, 'SampledAgg8June')

df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/fit_MVR_June18')

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

grid.arrange(p1, p2, p3, ncol=3)
#PostPredReal.pdf, 4 x 12.1 inches



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

#---------------------------------------------------------------------------------------------------------------------------------------------------

EucDist1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-07_112730.188674__July25Agg_Normal_EuclDist_RFwithTot_range.5_kappa.5_Chain1')
EucDist2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-29_132814.361033__July25Agg_Normal_EuclDist_RFwithTot_range.5_kappa.5_Chain2')
SR1 <- SR1
SR2 <- SR3

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


moveSubFiles <- function(filepath='IIDIPUS_Input_Alternatives/July25Agg/ODDobjects'){
  all_files <- list.files(path = filepath, recursive = TRUE, full.names = TRUE)
  file_paths <- all_files[file.info(all_files)$isdir == FALSE]
  
  # Move each file to the current directory
  file.rename(from = file_paths,
              to = file.path(filepath, basename(file_paths)))
}


moveTestData <- function(folder_in='IIDIPUS_Input_Alternatives/July25Agg'){
  
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
      total_mortalities <- c(total_mortalities, max(values(ODD[[grep('hazMean', names(ODD))]])))#[values(!is.na(ODD$Population) & ODD$Population > 100),], na.rm=T))
      # polygon_names <- unlist(lapply(ODD@polygons[ODD@impact$polygon], function(x) x$name))
      # if (any(tolower(polygon_names[which(ODD@impact$impact=='mortality')]) %in% c('tot', 'total'))){
      #   nonmatch <- which(!tolower(polygon_names[which(ODD@impact$impact=='mortality')]) %in% c('tot', 'total'))
      #   if (length(nonmatch)>0){
      #     ODD@impact <- ODD@impact[-which(ODD@impact$impact=='mortality')[nonmatch],] # in the case of total and subnational data, remove the subnational
      #   }
      # }
      # total_mortality <- sum(ODD@impact$observed[which(ODD@impact$impact=='mortality')])
      #total_mortalities <- c(total_mortalities, total_mortality)
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

df_postpredictive_sampled_best$train_flag2 <- 'TRAIN'
ufiles<-list.files(path=paste0(dir, 'IIDIPUS_Input_Alternatives/July25Agg/ODDobjects/Test'),pattern=Model$haz,recursive = T,ignore.case = T)
test_id = as.numeric(sub(".*_", "", ufiles))
for (i in 1:length(df_postpredictive_sampled_best$train_flag2)){
  if (df_postpredictive_sampled_best$event_id[i] %in% test_id){
    df_postpredictive_sampled_best$train_flag2[i] <- 'TEST'
  }
}

df_postpredictive_sampled_best$train_flag=df_postpredictive_sampled_best$train_flag2
mort_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'mortality' )  + guides(color="none") 
mort_all <- plot_df_postpredictive(df_postpredictive_sampled_best,'mortality' )  + guides(color="none") 
disp_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'displacement')  + guides(color="none")
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


