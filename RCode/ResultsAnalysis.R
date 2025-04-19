
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-24_121215_alpha0.95_M60_Npart990RealAgg3')

AlgoResults$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg3/'

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-10_192933_alpha0.9_M60_Npart1000RealAgg5')
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-07-29_140215_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2')
AlgoResults$input_folder = 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg5/'

AlgoResultsSMC <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-08-18_204305_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2')

AlgoResultsMCMC <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-08-19_075743_MCMC_RealAgg5_Trial_LR40_Rho0.975propCOVdiv12')


AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-09-07_081752_alphaAdaptive_M100_Npart1000_Sim50by50_propCOVmult0.2')
AlgoParams$cores <- 1
AlgoParams$NestedCores <- 4
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
AlgoParams$smc_Npart <- 50
AlgoParams$rel_weightings <- c(1,1)
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput2/'
dists <- array(NA, dim=c(100, 2))
for (i in 1:100){
  print(i)
  proposed<- Omega_true 
  proposed$eps$hazard_cor = i/100
  impact_sample <- SampleImpact(dir, Model, proposed %>% addTransfParams(), AlgoParams)
  dists[i,] = CalcDist(impact_sample, AlgoParams)[1:2]
}

#Simulated correlated vuln
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-09-22_151452_alphaAdaptive_M100_Npart1000_Sim50by50_propCOVmult0.2')
plot_correlated_posteriors(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(10,11), c(12,13), c(11, 14), c(18,15), c(16,17), c(19,20), c(21,22), c(9,23)), 
                           Omega=list(Lambda1 = list(nu=8.75, kappa=0.6),
                                      Lambda2 = list(nu=11.7, kappa=0.75), #list(nu=10.65, kappa=1.5), 
                                      Lambda3 = list(nu=8.55, kappa=0.8),
                                      Lambda4 = list(nu=9.9, kappa=1.6),
                                      theta= list(theta1=0.6),
                                      eps=list(local=0.8, hazard_mort=0.48, hazard_disp=0.6, hazard_bd=0.5, hazard_cor=0.55),
                                      #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                                      vuln_coeff = list(PDens=0, SHDI=-0.3, GNIc=-0.05, Vs30=0.1, EQFreq=-0.12, FirstHaz=0.05, Night=0, FirstHaz.Night=0.1),
                                      check = list(check=0.5)))
plot_d_vs_step(AlgoResults, ymax=6)


#AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-09-15_204947_alphaAdaptive_M100_Npart1000_Sim50by50_propCOVmult0.2')
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-09-18_080708_alphaAdaptive_M100_Npart1000_Sim50by50_propCOVmult0.2')
plot_correlated_posteriors(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(10,11), c(12,13), c(11, 14), c(18,15), c(16,17), c(19,20), c(21,22), c(9,23)), 
                           Omega=list(Lambda1 = list(nu=8.75, kappa=0.6),
                                      Lambda2 = list(nu=11.7, kappa=0.75), #list(nu=10.65, kappa=1.5), #
                                      Lambda3 = list(nu=8.7, kappa=0.7),
                                      Lambda4 = list(nu=9.9, kappa=1.6),
                                      theta= list(theta1=0.6),
                                      eps=list(local=0.8, hazard_mort=0.45, hazard_disp=0.6, hazard_bd=0.5, hazard_cor=0.55),
                                      #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                                      vuln_coeff = list(PDens=0, SHDI=-0.18, GNIc=-0.05, Vs30=0.1, EQFreq=-0.12, FirstHaz=0.05, Night=0, FirstHaz.Night=0.1),
                                      check = list(check=0.5)))
#SimPosteriors.pdf, 7 x 9 inches

plot_d_vs_step(AlgoResults, ymax=6)
points(AlgoResults$tolerancestore, pch=20, col='red')
#DvsStep.pdf, 6 x 9 inches

df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=100, output='SampledAgg')


df_postpredictive_sampled_true <- create_df_postpredictive(AlgoResults, single_particle=T, 
                                                           Omega =list(Lambda1 = list(nu=8.75, kappa=0.6),
                                                                       Lambda2 = list(nu=11.7, kappa=0.75), #list(nu=10.65, kappa=1.5), #
                                                                       Lambda3 = list(nu=8.7, kappa=0.7),
                                                                       Lambda4 = list(nu=9.9, kappa=1.6),
                                                                       theta= list(theta1=0.6),
                                                                       eps=list(local=0.8, hazard_mort=0.45, hazard_disp=0.6, hazard_bd=0.5, hazard_cor=0.55),
                                                                       #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                                                                       vuln_coeff = list(PDens=0, SHDI=-0.18, GNIc=-0.05, Vs30=0.1, EQFreq=-0.12, FirstHaz=0.05, Night=0, FirstHaz.Night=0.1),
                                                                       check = list(check=0.5)),
                                                           M=100, output='SampledAgg')

saveRDS(df_postpredictive_sampled_best, 'SimPostPredictive23rdSept')
saveRDS(df_postpredictive_sampled_true, 'SimTruePredictive23rdSept')

df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SimPostPredictive16thSept')
df_postpredictive_sampled_true <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SimTruePredictive16thSept')

plot_df_postpredictive_compare(df_postpredictive_sampled_true %>% filter(train_flag=='TEST'),
                               df_postpredictive_sampled_best %>% filter(train_flag=='TEST'), 'mortality')
#PostPredSim.pdf, 5 x 9 inches

#----------------------------------------------------------------------------------------------------
#------------------------------ OLD REAL DATA RESULTS -----------------------------------------------
#----------------------------------------------------------------------------------------------------

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-08-20_051627_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2')
plot_correlated_posteriors(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(10,11), c(12,13), c(11, 14), c(9,23)))#, c(18,15), c(16,17), c(19,20), c(21,22), c(9,23)))
#RealPosteriors.pdf, 5 x 10 inches


#plot_vuln(AlgoResults)#, c(18,15), c(16,17), c(19,20), c(21,22), c(9,23)))
#RealPosteriors.pdf, 5 x 10 inches
library(grid)
Model$skeleton <- list(Lambda1=list(nu=NA,kappa=NA), Lambda2=list(nu=NA,kappa=NA),
                       Lambda3=list(nu=NA,kappa=NA), Lambda4=list(nu=NA,kappa=NA),
                       eps=list(local=NA,hazard_mort=NA, hazard_disp=NA, hazard_bd=NA, hazard_cor=NA),
                       theta=list(theta=NA),
                       vuln_coeff=list(PDens=NA, SHDI=NA, GNIc=NA, Vs30=NA, EQFreq=NA, #Mag=NA,
                       FirstHaz=NA, Night=NA, FirstHaz.Night=NA), check=list(check=NA)
)
plot_vuln_posteriors(AlgoResults)
#VulnPosteriors.pdf, 6.5 x 12 inches


AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg5/'
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=50, output='SampledTotal')

df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledBest18thSept')
df_postpredictive_sampled_best1 = df_postpredictive_sampled_best
df_postpredictive_sampled_best1$observed[which(df_postpredictive_sampled_best1$observed == 17899)] = 8957
plot_df_postpredictive_PAGER_coloured(df_postpredictive_sampled_best1, 'mortality')

quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='mortality'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/101
hist(quants)
mean(quants>0.05 & quants < 0.95)
quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='displacement'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/101
hist(quants)
mean(quants>0.05 & quants < 0.95)
quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='buildDam'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/101
hist(quants)
mean(quants>0.05 & quants < 0.95)

hist(quants/101)
df_postpredictive_sampled_best[1,]

library(grid)
mort_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'mortality')  + guides(color="none") 
disp_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'displacement')  + guides(color="none")
buildDam_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'buildDam')  + guides(color="none") + 
                      scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma, limits=c(50,100000))
  

p1 <- arrangeGrob(mort_test, top = textGrob("(a)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p2 <- arrangeGrob(disp_test, top = textGrob("(b)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p3 <- arrangeGrob(buildDam_test, top = textGrob("(c)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))

grid.arrange(p1, p2, p3, ncol=3)
#PostPredReal.pdf, 4 x 12.1 inches

AlgoResults %<>% addAlgoParams(s_finish=NULL)
for (i in c(16,)){
  p_H_given_y = sum(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish]>0)/length(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish])
  p_H = sum(AlgoResults$Omega_sample_phys[,i,1]>0)/length(AlgoResults$Omega_sample_phys[,i,1])
  print(paste('Bayes factor for', names(unlist(Omega))[i], ':', round(p_H_given_y * (1-p_H)/((1-p_H_given_y)*p_H), 3)))
  print(paste('Post. prob of greater than 0:', sum(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish]>0)/length(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish])))
}


#----------------------------------------------------------------------------------------------------
#------------------------------ UPDATED REAL DATA RESULTS -------------------------------------------
#----------------------------------------------------------------------------------------------------

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-12-10_060355_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2')
plot_correlated_posteriors(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(8,11)))
#RealPosteriors.pdf, 5 x 10 inches

library(grid)
plot_vuln_posteriors(AlgoResults)
#VulnPosteriors.pdf, 6.5 x 12 inches

AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/Nov24Agg/'
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=50, output='SampledTotal')
df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledBest10thDec')

#df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledBest18thSept')

quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='mortality' & df_postpredictive_sampled_best$train_flag=='TEST'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/51
hist(quants)
mean(quants>0.05 & quants < 0.95)
quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='displacement'& df_postpredictive_sampled_best$train_flag=='TEST'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/51
hist(quants)
mean(quants>0.05 & quants < 0.95)
quants <- apply(df_postpredictive_sampled_best[which(df_postpredictive_sampled_best$impact=='buildDam'& df_postpredictive_sampled_best$train_flag=='TRAIN'),c(grep('observed', names(df_postpredictive_sampled_best)), grep('sampled.', names(df_postpredictive_sampled_best)))], 1, function(x){sample_quant(x)})/51
hist(quants)
mean(quants>0.05 & quants < 0.95)

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


#MCMC:

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-11_141748_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1')
AlgoResults$Omega_sample_phys <- AlgoResults$Omega_sample_phys[c(1:6, 10:22),]

#plot_correlated_posteriors(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(10,11), c(12,13), c(11, 14), c(9,23)))#, c(18,15), c(16,17), c(19,20), c(21,22), c(9,23)))
#RealPosteriors.pdf, 5 x 10 inches


#plot_vuln(AlgoResults)#, c(18,15), c(16,17), c(19,20), c(21,22), c(9,23)))
#RealPosteriors.pdf, 5 x 10 inches
#library(grid)
#plot_vuln_posteriors(AlgoResults)
#VulnPosteriors.pdf, 6.5 x 12 inches
AlgoResults$Sigma_store <- AlgoResults$Sigma_store[c(1:6, 10:22),c(1:6, 10:22),]

init_cov <- cov(readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/HLPriorSamples'))
init_cov2 <- cov(readRDS(paste0(dir, 'IIDIPUS_Input/HLPriorSamples2')))
init_cov <- init_cov[c(1:6, 10:22),c(1:6, 10:22)]
param <- 12
plot(AlgoResults$Sigma_store[param,param,], ylim=range(c(init_cov[param,param], init_cov2[param, param], AlgoResults$Sigma_store[param,param,]), na.rm=T))
points(1, init_cov[param,param], col='red')
abline(h=init_cov2[param, param]/5, col='blue')

AlgoParams$input_folder <- AlgoResults$input_folder
df_postpredictive_sampled_best <- create_df_postpredictive_MCMC(AlgoResults, single_particle=F, Omega=NULL, particle_best=F,
                                                           M=20, output='SampledTotal')

df_postpredictive_sampled_best_old <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledBest18thSept')

#hist(quants/101)
#df_postpredictive_sampled_best[1,]

library(grid)

df_postpredictive_sampled_best$train_flag2 <- 'TRAIN'
ufiles<-list.files(path=paste0(dir, 'IIDIPUS_Input_Alternatives/Nov24Agg/ODDobjects/Test'),pattern=Model$haz,recursive = T,ignore.case = T)
test_id = as.numeric(sub(".*_", "", ufiles))
for (i in 1:length(df_postpredictive_sampled_best$train_flag2)){
  if (df_postpredictive_sampled_best$event_id[i] %in% test_id){
    df_postpredictive_sampled_best$train_flag2[i] <- 'TEST'
  }
}

df_postpredictive_sampled_best$train_flag=df_postpredictive_sampled_best$train_flag2
mort_test <- plot_df_postpredictive(df_postpredictive_sampled_best%>% filter(train_flag=='TEST'),'mortality' )  + guides(color="none") 
disp_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'displacement')  + guides(color="none")
buildDam_test <- plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'buildDam')  + guides(color="none") + 
  scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma, limits=c(50,100000))


p1 <- arrangeGrob(mort_test, top = textGrob("(a)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p2 <- arrangeGrob(disp_test, top = textGrob("(b)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))
p3 <- arrangeGrob(buildDam_test, top = textGrob("(c)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 13)))

mort_test
disp_test
buildDam_test


#-----
source('RCode/ODDobj_noRF.R')
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-12-10_060355_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2')
AlgoResults %<>% addAlgoParams()
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/Mar25Agg/'

Omega_trial <- list(Lambda1 = list(nu=8.75, kappa=0.6),
                            Lambda2 = list(nu=12, kappa=1), #list(nu=10.65, kappa=1.5), #
                            Lambda3 = list(nu=9.55, kappa=0.68),
                            #Lambda4 = list(nu=9.9, kappa=1.6),
                            #theta= list(theta1=0.6),
                            eps=list(local=0.25, hazard_mort=1, hazard_disp=0.6, hazard_bd=0.5, hazard_cor=0.55),
                            #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                            vuln_coeff = list(PDens=-0.1239836, SHDI=-0.6717716, GNIc=0.1994072, Vs30=0.473589, EQFreq=-0.5352366, FirstHaz=0.08874202, Night=0.03315047, FirstHaz.Night=-0.1316748))

df_postpredictive_sampled_total <- create_df_postpredictive(AlgoResults, single_particle=T, 
                                                           M=100, output='SampledTotal', 
                                                           Omega = Omega_trial %>% addTransfParams())

plot_df_postpredictive(df_postpredictive_sampled_total %>% filter(train_flag=='TRAIN'),'mortality')  + guides(color="none") 

#----------------------------------------------------------------------------------------------------
#------------------------------ Band Depth Score Results  -------------------------------------------
#----------------------------------------------------------------------------------------------------

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-27_162932.543295_MCMC_BDScore.05nocorr_LR40_M100_Npart1000NovAgg5_RandomFieldThree')
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/Nov24Agg/'
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/Mar25Agg/'
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/Mar25Agg/'
df_postpredictive_sampled <- create_df_postpredictive_MCMC(AlgoResults, single_particle=F, 
                                                           M=100, output='SampledTotal')

#saveRDS(df_postpredictive_sampled, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledAggBandDepth24thMarch')
#saveRDS(df_postpredictive_sampled, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledTotalBandDepth25thMarch')


library(plotly)
df_postpredictive_sampled <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledTotalBandDepth25thMarch')
plot_df_postpredictive_PAGER_coloured(df_postpredictive_sampled, 'mortality', interactive=T) 



plot(density(rnorm(10000)), freq=F)
lines(density(rt(10000,3)*0.7), freq=F, col='red')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-27_162932.543295_MCMC_BDScore.05nocorr_LR40_M100_Npart1000NovAgg5_RandomFieldThree')
AlgoParams$input_folder = 'IIDIPUS_Input_Alternatives/Nov24Agg/'
proposed = AlgoResults$Omega_sample_phys[,1] %>% relist(skeleton=Model$skeleton) %>% addTransfParams()
xx <- SampleImpact(dir, Model, proposed, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 2) %>% 
               replace(which(names(AlgoParams)==c('Np')), 1),)

n_post_samples = 10
ODD <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDAggobjects/EQ20230908MAR_-1_AggLevel5')
for (i in 1:n_post_samples){
  
  #plot: 
  #grid.arrange(plotODDy_GADM(ODD, 'ISO3C', gadm_level=2, haz_legend=T, var_legend=T, var_discrete=T),plotODDy_GADM(ODD, 'Population', gadm_level=2, haz_legend=T, var_legend=T, var_discrete=F, log_legend=T), ncol=2)
  Omega = AlgoResults$Omega_sample_phys[,1] %>% relist(skeleton=Model$skeleton)
  sampled_out <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                         replace(which(names(AlgoParams)==c('Np')), 1), 
                       output='Sampled')
}


#-------------------------------------
#------------No RF Total--------------
#-------------------------------------

AlgoResults_noTotErr <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-04-13_120740.076549_MCMC_BDScore.05nocorr_M100_Npart1000NovAgg5_RandomFieldThree_rfNoTot')
AlgoResults_noTotErr$input_folder <- 'IIDIPUS_Input_Alternatives/Mar25Agg/'

plot(AlgoResults$loss, type='l')
lines(AlgoResults_noTotErr$loss, col='blue')


mcmc_trace(AlgoResults_noTotErr, AlgoResults,  var_plot=1:19, xlim=c(1, 2000))



AlgoResults_LogLinear <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-04-12_195835_MCMC_BDScore.05nocorr_M100_Npart1000MarAgg5_RandomFieldThree_LogLinear')
plot(AlgoResults$loss, type='l')
lines(AlgoResults_LogLinear$loss, type='l', col='red')







#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

mcmc_trace <- function(..., var_plot = 7:10, xlim = c(4000, 8000)) {
  
  chains <- list(...)
  num_chains <- length(chains)
  
  if (num_chains < 1) {
    stop("At least one chain must be provided.")
  }
  
  # Define subfigure labels (a), (b), (c), ...
  subfig_labels <- letters[1:20]
  plot_list <- list()
  
  for (i in seq_along(var_plot)) {
    var <- var_plot[i]
    
    # Create base data frame with x values
    post_samples <- data.frame(x = (xlim[1]:xlim[2]) - xlim[1])
    
    # Loop through each chain and add its samples to the data frame
    for (j in seq_along(chains)) {
      chain_name <- paste0("Chain", j)
      post_samples[[chain_name]] <- chains[[j]]$Omega_sample_phys[var, xlim[1]:xlim[2]]
    }
    
    # Convert to long format
    post_samples_long <- post_samples %>%
      pivot_longer(cols = -x, names_to = "Parameter", values_to = "Value")
    
    
    if (length(grep('vuln_coeff.', names(unlist(Model$skeleton))[var])>0)){
      y_lb = -1
      y_ub = 1
    } else {
      y_lb = Model$par_lb[var]
      y_ub = Model$par_ub[var]
    }
    
    # Plot the lines
    plot <- ggplot(data = post_samples_long, aes(x = x, y = Value, color = Parameter)) +
      labs(y = get_greek_titles(names(unlist(Model$skeleton))[var]),
           x = 'Post warmup iteration') +
      scale_y_continuous(limits = c(y_lb, y_ub), expand = c(0.01, 0.01)) +
      geom_line(alpha = 0.8) +
      theme_minimal() +
      scale_color_viridis(discrete = TRUE) +
      theme(
        axis.title.y = element_text(family = "Times New Roman", size = 12),
        axis.text.x = element_text(family = "Times New Roman", size = 12),
        axis.text.y = element_text(family = "Times New Roman", size = 12),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(family = "Times New Roman", size = 12),
        plot.title = element_text(family = "Times New Roman", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = unit(c(0, 20, 0, 15), "pt"),
        legend.position = "none"
      )
    
    # Add subfigure label
    plot_with_label <- arrangeGrob(
      plot,
      top = textGrob(paste0("(", subfig_labels[i], ")"),
                     x = unit(0.0, "npc"), y = unit(1, "npc"),
                     just = c("left", "top"),
                     gp = gpar(fontsize = 14, fontfamily = "Times New Roman"))
    )
    
    plot_list[[i]] <- plot_with_label
  }
  
  # Arrange plots in a 1-row grid
  grid.arrange(grobs = plot_list, ncol = 4)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------


#MCMC Single:

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-10-21_120105_MCMC_RealAgg5_LR40_Rho0.9_adaptive')
plot(AlgoResults$loss[100:which(!is.finite(AlgoResults$loss))[1]], type='l')

#SMC:
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-10-01_223810_alphaAdaptivePropCOV0.2_M100_Npart1000SimInput')

#MCMC:
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-07_103308_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP')

AlgoResults2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2024-11-06_143442_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV')

AlgoResults2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-10-12_222810_MCMC_SimCorrVuln_Trial_LR40_Rho0.9_adaptive')
plot(AlgoResults$loss, type='l', ylim=c(2.2, 5))
points(AlgoResults$loss, type='l', col='red')

par(mfrow=c(2,2))
params <- c(9,23)
plot(AlgoResults2$Omega_sample_phys[params[1],], AlgoResults2$Omega_sample_phys[params[2],], type='l')
points(AlgoResults2$Omega_sample_phys[params[1],1], AlgoResults2$Omega_sample_phys[params[2],1], col='blue', pch=19)

params <- c(9,23)
plot(AlgoResults$Omega_sample_phys[params[1],2000:5000], AlgoResults$Omega_sample_phys[params[2],2000:5000], type='l')
plot(AlgoResults$Omega_sample_phys[params[2],], type='l', xlim=c(4000,5000))
plot(AlgoResults$Sigma_store[params[2],params[2],]* AlgoResults$lambda_store, type='l', ylim=c(0,0.5))










plot(AlgoResults2$lambda_store * AlgoResults2$Sigma_store[19,19,], type='l', xlim=c(0,1000))

xlim=range(AlgoResults2$Omega_sample_phys[params[1],], AlgoResults$Omega_sample_phys[params[1],],na.rm=T)
ylim=range(AlgoResults2$Omega_sample_phys[params[2],], AlgoResults$Omega_sample_phys[params[2],], na.rm=T)
plot(AlgoResults2$Omega_sample_phys[params[1],10000:12000], AlgoResults2$Omega_sample_phys[params[2],10000:12000], xlim=xlim, ylim=ylim,type='l')
points(AlgoResults$Omega_sample_phys[params[1],1:6000], AlgoResults$Omega_sample_phys[params[2],1:6000], type='l', col='red')

xlim=range(AlgoResults$Omega_sample_phys[params[1],], AlgoResults$Omega_sample_phys[params[1],],na.rm=T)
ylim=range(AlgoResults2$Omega_sample_phys[params[2],], AlgoResults$Omega_sample_phys[params[2],], na.rm=T)
plot(AlgoResults2$Omega_sample_phys[params[1],10000:12000], AlgoResults2$Omega_sample_phys[params[2],10000:12000], xlim=xlim, ylim=ylim,type='l')
points(AlgoResults$Omega_sample_phys[params[1],1:6000], AlgoResults$Omega_sample_phys[params[2],1:6000], type='l', col='red')
xlim=range(AlgoResults2$Omega_sample[params[1],], AlgoResults$Omega_sample[params[1],],na.rm=T)
ylim=range(AlgoResults2$Omega_sample[params[2],], AlgoResults$Omega_sample[params[2],], na.rm=T)
plot(AlgoResults2$Omega_sample[params[1],], AlgoResults2$Omega_sample[params[2],], xlim=xlim, ylim=ylim,type='l')
points(AlgoResults$Omega_sample[params[1],1:6000], AlgoResults$Omega_sample[params[2],1:6000], type='l', col='red')
plot(AlgoResults2$lambda_store * AlgoResults2$Sigma_store[params[1],params[1],],type='l', ylim=c(0,0.5))
lines(AlgoResults$lambda_store * AlgoResults$Sigma_store[params[1],params[1],],type='l', col='red')


params <- c(22,23)
xlim=range(AlgoResults2$Omega_sample_phys[params[1],], AlgoResults$Omega_sample_phys[params[1],],na.rm=T)
ylim=range(AlgoResults2$Omega_sample_phys[params[2],], AlgoResults$Omega_sample_phys[params[2],], na.rm=T)
plot(AlgoResults2$Omega_sample_phys[params[1],10000:14000], AlgoResults2$Omega_sample_phys[params[2],10000:14000], xlim=xlim, ylim=ylim,type='l')
points(unlist(Omega)[params[1]], unlist(Omega)[params[2]], col='red', pch=19)

plot(AlgoResults$lambda_store)
points(AlgoResults2$lambda_store, col='red')

par(mfrow=c(2,1))
plot(AlgoResults$lambda_store * AlgoResults$Sigma_store[1,1,])
points(AlgoResults2$lambda_store * AlgoResults2$Sigma_store[1,1,], col='red')
plot(AlgoResults$lambda_store * AlgoResults$Sigma_store[7,7,])
points(AlgoResults2$lambda_store * AlgoResults2$Sigma_store[7,7,], col='red')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-09-24_121551_MCMC_SimCorrVuln_Trial_LR40_Rho0.9_adaptive_backup')
Prior <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HLPriorSamples')
prior_phys <- Prior
for (i in 1:NROW(Prior)){
  prior_phys[i, ] <- Proposed2Physical(Prior[i,] %>% relist(skeleton=Model$skeleton) %>% unlist(), Model) %>% unlist()
}


params <- c(15,16)
plot(prior_phys[,params[1]], prior_phys[,params[2]])
points(AlgoResults$Omega_sample_phys[params[1],1], AlgoResults$Omega_sample_phys[params[2],1], col='red', pch=19)
lines(AlgoResults$Omega_sample_phys[params[1],], AlgoResults$Omega_sample_phys[params[2],], col='blue')


#plot_df_postpredictive(df_postpredictive_sampled_true %>% filter(train_flag=='TEST'), 'mortality')
#plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'), 'mortality')


plot(AlgoResultsMCMC$loss, type='l')

plot(AlgoResultsSMC$Omega_sample_phys[,16,1], AlgoResultsSMC$Omega_sample_phys[,17,1], col='blue', xlab='SHDI.coef', ylab='GNIc.coef')
points(AlgoResultsSMC$Omega_sample_phys[,16,163], AlgoResultsSMC$Omega_sample_phys[,17,163], pch=19)
points(AlgoResultsMCMC$Omega_sample_phys[16,], AlgoResultsMCMC$Omega_sample_phys[17,], col='red', pch=19)


#compare MCMC and ABC-SMC:
AlgoResults_smc <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-11-14_224431_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2')
AlgoResults_mcmc <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-16_172107_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1_backup')
stop <- which(is.na(AlgoResults_mcmc$Omega_sample_phys[1,]))[1]
params <- c(12,15)
plot(AlgoResults_smc$Omega_sample_phys[,params[1],1], AlgoResults_smc$Omega_sample_phys[,params[2],1], col='blue',xlab=names(unlist(Omega))[params[1]], ylab=names(unlist(Omega))[params[2]])
points(AlgoResults_smc$Omega_sample_phys[,params[1],60], AlgoResults_smc$Omega_sample_phys[,params[2],60], col='red')
points(AlgoResults_mcmc$Omega_sample_phys[params[1],(stop-3000):stop], AlgoResults_mcmc$Omega_sample_phys[params[2],(stop-3000):stop], col='green')


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

#MCMC single:
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-10_084817_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP')
params <- c(9,23)
stop1 <- which(is.na(AlgoResults1$Omega_sample_phys[params[1],]))[1]-1
plot(AlgoResults1$Omega_sample_phys[params[1],(stop1-3000):stop1], AlgoResults1$Omega_sample_phys[params[2],(stop1-3000):stop1],type='l',
     xlab ='Dummy 3', ylab='Dummy 4') #names(unlist(Omega))[params[1]], ylab = names(unlist(Omega))[params[2]])
abline(h=0)
abline(v=0)
plot(AlgoResults2$Omega_sample_phys[params[2],], type='l')

#MCMC compare:
AlgoResults1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-20_163934_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1')
AlgoResults2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-20_163933_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2_backup')

param = 3
plot(AlgoResults1$lambda_store* AlgoResults1$Sigma_store[param, param,], xlim=c(0, 9000))
lines(AlgoResults2$lambda_store* AlgoResults2$Sigma_store[param, param, ], col='red')

plot(AlgoResults1$loss, type='l', ylim=c(3.5, 5), xlim=c(0,8000))
points(AlgoResults2$loss, type='l', col='red')

plot(AlgoResults2$Omega_sample_phys[5,1000:2000], AlgoResults2$loss[1000:2000])

par(mfrow=c(2,2))
params <- c(11,12)
xlim=range(AlgoResults1$Omega_sample_phys[params[1],], AlgoResults2$Omega_sample_phys[params[1],],na.rm=T)
ylim=range(AlgoResults1$Omega_sample_phys[params[2],], AlgoResults2$Omega_sample_phys[params[2],], na.rm=T)
stop1 <- which(is.na(AlgoResults1$Omega_sample_phys[params[1],]))[1]-1
stop2 <- which(is.na(AlgoResults2$Omega_sample_phys[params[1],]))[1]-1
plot(AlgoResults1$Omega_sample_phys[params[1],(stop1-4000):stop1], AlgoResults1$Omega_sample_phys[params[2],(stop1-4000):stop1], xlim=xlim, ylim=ylim,type='l',
     xlab = names(unlist(Omega))[params[1]], ylab = names(unlist(Omega))[params[2]] )
#points(AlgoResults1$Omega_sample_phys[params[1],1], AlgoResults1$Omega_sample_phys[params[2],1], xlim=xlim, ylim=ylim,pch=19, col='green')
points(AlgoResults2$Omega_sample_phys[params[1],(stop2-4000):stop2], AlgoResults2$Omega_sample_phys[params[2],(stop2-4000):stop2], type='l', col='blue')
#points(AlgoResults2$Omega_sample_phys[params[1],1], AlgoResults2$Omega_sample_phys[params[2],1], col='green', pch=19)
points(AlgoResults1$Omega_sample_phys[params[1],stop1], AlgoResults1$Omega_sample_phys[params[2],stop1], xlim=xlim, ylim=ylim,pch=19, col='red')
points(AlgoResults2$Omega_sample_phys[params[1],stop2], AlgoResults2$Omega_sample_phys[params[2],stop2], xlim=xlim, ylim=ylim,pch=19, col='red')

params[1] <- 4
plot(AlgoResults1$Omega_sample_phys[params[1],], type='l')
lines(AlgoResults2$Omega_sample_phys[params[1],], type='l', col='blue')


#compare to previous abc smc results:
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_start_step_2024-11-26_030403_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2')
par(mfrow=c(2,1))
params[1] = 2
hist(AlgoResults$Omega_sample_phys[,params[1],115], type='l', freq=F)
lines(density(AlgoResults$Omega_sample_phys[,params[1],1]), type='l', freq=F, col='blue')
lines(density(AlgoResults1$Omega_sample_phys[params[1],4000:6000]), type='l', col='red')
lines(density(AlgoResults2$Omega_sample_phys[params[1],4000:6000]), type='l', col='red')
plot(AlgoResults$Omega_sample_phys[,params[1],115], AlgoResults$d[,1,115])




xlim=range(AlgoResults1$Omega_sample[params[1],], AlgoResults2$Omega_sample[params[1],],na.rm=T)
ylim=range(AlgoResults1$Omega_sample[params[2],], AlgoResults2$Omega_sample[params[2],], na.rm=T)
plot(AlgoResults1$Omega_sample[params[1],3500:4500], AlgoResults1$Omega_sample[params[2],3500:4500], xlim=xlim, ylim=ylim,type='l')
points(AlgoResults2$Omega_sample[params[1],], AlgoResults2$Omega_sample[params[2],], type='l', col='red')
points(AlgoResults1$Omega_sample[params[1],stop], AlgoResults1$Omega_sample[params[2],stop], xlim=xlim, ylim=ylim,pch=19, col='blue')

HLPrior_sample <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/HLPriorSamples')
param <- 16
plot(seq(1,3600, length.out=500), HLPrior_sample[,param], ylab='SHDI Coefficient (untransformed)')
lines(AlgoResults1$Omega_sample[param,], col='red', type='l')


plot(AlgoResults1$loss, type='l', ylim=c(3.8, 4.5))
lines(AlgoResults2$loss, type='l', col='red')
AlgoResults1$loss 
plot(AlgoResults1$lambda_store * AlgoResults1$Sigma_store[23,23,], type='l', ylim=c(0,3.5))
points(AlgoResults2$lambda_store * AlgoResults2$Sigma_store[23,23,], col='red', type='l')

mean(AlgoResults1$accprob_store[1:500], na.rm=T)

par(mfrow=c(2,1))
plot(log(AlgoResults1$lambda_store), xlim=c(0,2500), ylim=c(-2.5,0))
#plot(log(AlgoResults1$lambda_store), xlim=c(0,2000))
abline(h=-log(2.5))
#abline(h=-log(3))
abline(h=-2*log(2.5))
abline(v=125)
lines(log(AlgoResults2$lambda_store), col='red')
accprob_mean <- cumsum(AlgoResults1$accprob_store)/seq_along(AlgoResults1$accprob_store)
plot(rollmean(AlgoResults1$accprob_store, 21, fill=NA, align='center'), xlim=c(0, 6000), type='l', ylab='Rolling Mean of Acceptance Probability')
abline(h=0.234)

plot(AlgoResults1$lambda_store * AlgoResults1$Sigma_store[23,23,])
#plot(log(AlgoResults1$lambda_store), xlim=c(0,2000))
lines(AlgoResults2$lambda_store* AlgoResults2$Sigma_store[23,23,], col='red')

plot(log(AlgoResults1$lambda_store))

HLP_vals <- c()
for (i in 1:8000){
  HLP_vals <- c(HLP_vals,Model$HighLevelPriors(AlgoResultsMCMC$Omega_sample_phys[,i] %>% 
                                      relist(skeleton=Model$skeleton) %>% 
                                      addTransfParams(),
                                    Model))
}
plot(HLP_vals, type='l')

#MCMC Compare Euclidean distance to Scoring Rule Posterior
EucDist1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-25_151709_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat_EuclDist_backup')
EucDist2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2024-12-01_010730.812193_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat_EuclDist2')
SR1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180047_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2')
SR2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180050_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1_backup')

plot_mcmc_compare = function(EucDist1, EucDist2, SR1, SR2, xlim=c(4000, 8000)){
  
  # Define subfigure labels (a), (b), (c), ...
  subfig_labels <- letters[1:8]
  plot_list <- list()
  var_plot <- 7:10
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
  grid.arrange(grobs = plot_list, ncol = 4, nrow = 1)
  #Trace_EuclDist.pdf, 12 x 9 (or 13 x 3 for single row)
}

#Postpredictive comparison
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/Nov24Agg/'

df_postpredictive_sampled_best <- create_df_postpredictive_MCMC(EucDist1, single_particle=F, Omega=NULL, particle_best=F,
                                                                M=50, output='SampledTotal')
saveRDS(df_postpredictive_sampled_best, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/EuclDist25thNov')
df_postpredictive_sampled_bestEuc <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/EuclDist25thNov')

df_postpredictive_sampled_bestSR <- create_df_postpredictive_MCMC(SR1, single_particle=F, Omega=NULL, particle_best=F,
                                                                M=50, output='SampledTotal')
saveRDS(df_postpredictive_sampled_bestSR, '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SR25thNov')
df_postpredictive_sampled_bestSR <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SR25thNov')

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


# df_postpredictive_sampled_true <- create_df_postpredictive(AlgoResults, single_particle=T, 
#                                                            Omega =list(Lambda1 = list(nu=8.75, kappa=0.6),
#                                                                        Lambda2 = list(nu=11.7, kappa=0.75), #list(nu=10.65, kappa=1.5), #
#                                                                        Lambda3 = list(nu=8.7, kappa=0.7),
#                                                                        Lambda4 = list(nu=9.9, kappa=1.6),
#                                                                        theta= list(theta1=0.6),
#                                                                        eps=list(local=0.8, hazard_mort=0.45, hazard_disp=0.6, hazard_bd=0.5, hazard_cor=0.55),
#                                                                        #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
#                                                                        vuln_coeff = list(PDens=0, SHDI=-0.18, GNIc=-0.05, Vs30=0.1, EQFreq=-0.12, FirstHaz=0.05, Night=0, FirstHaz.Night=0.1),
#                                                                        check = list(check=0.5)),
#                                                            M=100, output='SampledAgg')
# 
# saveRDS(df_postpredictive_sampled_best, 'SimPostPredictive23rdSept')
# saveRDS(df_postpredictive_sampled_true, 'SimTruePredictive23rdSept')

# df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SimPostPredictive16thSept')
# df_postpredictive_sampled_true <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/SimTruePredictive16thSept')

plot_df_postpredictive_compare(df_postpredictive_sampled_bestEuc %>% filter(train_flag=='TEST'),
                               df_postpredictive_sampled_bestSR %>% filter(train_flag=='TEST'), 'buildDam')

# param = 5
# for (param in 1:length(unlist(Omega))){
#   plot(EucDist1$Omega_sample_phys[param,], type='l', lwd=1,col='#440154', xlim=c(5000,10000)) # Dark Purple
#   lines(EucDist2$Omega_sample_phys[param,], type='l', col='#31688E') # Blue
#   lines(SR1$Omega_sample_phys[param,], type='l', col='#35B779') # Greenish Teal
#   lines(SR2$Omega_sample_phys[param,], type='l', col='#FDE725') # Yellow
# }

#Compare ABC-SMC and MCMC

SR1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-12-08_171036_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2')
SR2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180050_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1_backup')
SR3 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-12-08_180223_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat3_backup')



plot_mcmc_compare = function(SR1, SR2, SR3, xlim=c(1, 9000)){
  
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
      SR2_param <-SR2$Omega_sample_phys[var,xlim[1]:xlim[2]])#,
      #SR3_param <-SR3$Omega_sample_phys[var,xlim[1]:xlim[2]])
    
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

AlgoResultsABC <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-12-10_060355_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2_further')
param=10
hist(AlgoResultsABC$Omega_sample_phys[,param,163], freq=F)
lines(density(SR1$Omega_sample_phys[param,], na.rm=T))
lines(density(SR2$Omega_sample_phys[param,], na.rm=T), col='blue')

AlgoResultsMCMC <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180050_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1_backup')
#AlgoResultsMCMC <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180047_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2')
AlgoResultsMCMC$HLP_vals <- rep(1, length(AlgoResultsMCMC$loss))
for (i in 1:8000){
  AlgoResultsMCMC$HLP_vals[i] = Model$HighLevelPriors(AlgoResultsMCMC$Omega_sample_phys[,i] %>% 
                                                 relist(skeleton=Model$skeleton) %>% 
                                                 addTransfParams(),
                                               Model)
}

darker_colors <- c("yellow", "red", "blue")  # Darker shades of blue, green, purple

# Ensure HLP_vals is mapped to these colors
point_colors <- darker_colors[as.factor(AlgoResultsMCMC$HLP_vals)]
plot(
  AlgoResultsMCMC$Omega_sample_phys[5,],
  AlgoResultsMCMC$Omega_sample_phys[6,],
  col = point_colors)

plot_correlated_posteriors(AlgoResultsABC, include_priors=T, Omega=NULL,
                                pairings=rbind(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,10), c(15,12), c(13,14), c(16, 17), c(18,19)), s_finish=NULL, 
                                AlgoResultsMCMC, subfig_title_adj=-0.3)
#MCMCvsSMC.pdf, 9 x 10.5 (7 x 12 for wider)
# 7 x 13 for presentation plot





#Compare LogNormal and Normal Models

Normal1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180047_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2')
Normal2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180050_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1_backup')
LogNormal <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2024-12-07_203532.216188_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1_LogNormal_backup')

plot(LogNormal$loss, ylim=c(3.5, 5), type='l')
lines(Normal1$loss, ylim=c(3.5, 5), type='l', col='red')
lines(Normal2$loss, ylim=c(3.5, 5), type='l', col='green')

#Distance measure exploration:
plot_distance_measure_exp <- function(){
  
  impact_type='mortality'
  df_poly_jitter <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledBest18thSept') %>% filter(train_flag=='TEST')
  #plot_df_postpredictive(df_postpredictive_sampled_best %>% filter(train_flag=='TEST'),'mortality')  + guides(color="none") 
  
  jitter_val <- function(x){
    return_arr <- c()
    for (x_i in x){
      if(x_i==0){return_arr <- c(return_arr, runif(1,0.2, 0.3))}
      else {return_arr <- c(return_arr, runif(1, x_i-0.5, x+0.5))}
    }
    return(return_arr)
  }
  #df_poly_jitter[which(df_poly_jitter==0, arr.ind=T)] = 0.1
  df_poly_jitter$observed <- sapply(df_poly_jitter$observed, jitter_val)
  df_poly_jitter[,grep('sampled', names(df_poly_jitter))] <- apply(df_poly_jitter[,grep('sampled', names(df_poly_jitter))],1:2,jitter_val)
  
  
  df_poly_jitter[,grep('sampled', names(df_poly_jitter))] <- t(apply(df_poly_jitter[,grep('sampled', names(df_poly_jitter))],1,sort))
  M <- length(grep('sampled', names(df_poly_jitter)))
  M_lower <- round(quantile(1:M, 0.05))
  M_lower <- ifelse(M_lower==0, 1, M_lower)
  M_median <- round(quantile(1:M, 0.5))
  M_upper <- round(quantile(1:M, 0.95))
  df_poly_jitter$medians_sampled <- apply(df_poly_jitter[,grep("sampled",names(df_poly_jitter),value = T)], 1, median)
  
  
  
  df_poly_filt <- df_poly_jitter %>% filter(impact==impact_type)
  mean(df_poly_filt$observed <= df_poly_filt[,paste0('sampled.', M_lower)]) + mean(df_poly_filt$observed >= df_poly_filt[,paste0('sampled.', M_upper)])
  
  # 
  SS_res_raw = sum((df_poly_filt$observed-df_poly_filt$medians_sampled)^2)
  R_squared_raw = 1-SS_res_raw/SS_tot_raw
  print(paste('R squared raw', R_squared_raw))
  
  k <- 10
  SS_res_log = sum((log(df_poly_filt$observed+k)-log(df_poly_filt$medians_sampled+k))^2)
  SS_tot_log = sum((log(df_poly_filt$observed+k)-mean(log(df_poly_filt$observed+k)))^2)
  R_squared_log = 1-SS_res_log/SS_tot_log
  print(paste('R squared log', R_squared_log))
  
  
  ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=get(paste0('sampled.', M_median)), 
                                                                     ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
    #ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=means_sampled, ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
    geom_point(aes(col='black')) + 
    scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    scale_y_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    #geom_pointrange(aes(col=train_flag)) + 
    geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
    ylab(paste('Sampled', ifelse(impact_type=='buildDam', 'building damage', impact_type))) + xlab(paste('Observed', ifelse(impact_type=='buildDam', 'building damage', impact_type))) + scale_color_manual(values = c('red', 'blue')) + 
    theme_bw() +
    theme(axis.title = element_text(family = "Liberation Serif", size=12),  
          legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
          legend.title = element_text(family = "Liberation Serif", size=12)) + guides(color="none") 
  
  g1 <- ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=get(paste0('sampled.', M_median)))) + 
                                                                     #ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
    #ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=means_sampled, ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
    geom_point() + 
    scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    scale_y_continuous(trans='log10', limits=c(0.25, 100000),breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    #geom_pointrange(aes(col=train_flag)) + 
    geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
    ylab(paste('Sampled', ifelse(impact_type=='buildDam', 'building damage', impact_type))) + xlab(paste('Observed', ifelse(impact_type=='buildDam', 'building damage', impact_type))) +#+ scale_color_manual(values = c('red', 'blue')) + 
    theme_bw() +
    theme(axis.title = element_text(family = "Liberation Serif", size=12),  
          legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
          legend.title = element_text(family = "Liberation Serif", size=12)) + guides(color="none") +
    geom_point(data=data.frame(x=17899, y=2058), aes(x=x, y=y), col='red') +
    geom_point(data=data.frame(x=17899, y=rlnorm(15, log(2058), 0.5)), aes(x=x, y=y), col='red')
  
  g2 <- ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=get(paste0('sampled.', M_median)))) + 
    #ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
    #ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=means_sampled, ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
    geom_point() + 
    scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    scale_y_continuous(trans='log10', limits=c(0.25, 100000),breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    #geom_pointrange(aes(col=train_flag)) + 
    geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
    ylab(paste('Sampled', ifelse(impact_type=='buildDam', 'building damage', impact_type))) + xlab(paste('Observed', ifelse(impact_type=='buildDam', 'building damage', impact_type))) +#+ scale_color_manual(values = c('red', 'blue')) + 
    theme_bw() +
    theme(axis.title = element_text(family = "Liberation Serif", size=12),  
          legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
          legend.title = element_text(family = "Liberation Serif", size=12)) + guides(color="none") +
    geom_point(data=data.frame(x=17899, y=2058), aes(x=x, y=y), col='red') +
    geom_point(data=data.frame(x=17899, y=rlnorm(15, log(2058), 2)), aes(x=x, y=y), col='red')
  
  grid.arrange(g1, g2, ncol=2)
  #5 x 11, Sim_vs_obs_calibration
  
  
}

#Scoring rule with and without variogram score
SR1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180047_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2')
SR2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180050_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat1_backup')
SR_vs <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-01-10_181501_MCMC_VariogramScore_M100_Npart1000NovAgg5_propCOVmult0.2')

plot_mcmc_compare = function(SR1, SR2, SR_vs, xlim=c(1, 10000)){
  
  # Define subfigure labels (a), (b), (c), ...
  subfig_labels <- letters[1:16]
  plot_list <- list()
  var_plot <- 12:19
  # Loop through the variables and generate each plot
  for (i in seq_along(var_plot)){#seq_along(vuln_var)) {
    #var <- vuln_var[i]
    var <- var_plot[i]
    
    post_samples <- data.frame(
      x=(xlim[1]:xlim[2])-xlim[1],
      SR1_param <- SR1$Omega_sample_phys[var,xlim[1]:xlim[2]],
      SR2_param <-SR2$Omega_sample_phys[var,xlim[1]:xlim[2]],
      SRvs_param <-SR_vs$Omega_sample_phys[var,xlim[1]:xlim[2]])
    
    post_samples_long <- post_samples %>%
      pivot_longer(cols = -x, names_to = "Parameter", values_to = "Value")
    
    # Plot the lines
    plot <- ggplot(data = post_samples_long, aes(x = x, y = Value, color = Parameter)) +
      labs(y = get_greek_titles(names(unlist(Model$skeleton))[var]),
           x= 'Post warmup iteration') +
      scale_y_continuous(limits = c(-1,1), expand = c(0.01, 0.01)) +
      #scale_y_continuous(limits = c(Model$par_lb[var], Model$par_ub[var]), expand = c(0.01, 0.01)) +
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
  grid.arrange(grobs = plot_list, ncol = 4, nrow = 2)
  #Trace_EuclDist.pdf, 12 x 9 (or 13 x 3 for single row)
}

# AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-01-25_125459_MCMC_VariogramScore_M100_Npart1000NovAgg5_RandomField_backup')
# HLPrior_samples <- t(AlgoResults$Omega_sample[,300:900])
# cov(HLPrior_samples)
# Proposed2Physical(HLPrior_samples[300,] %>% relist(skeleton=Model$skeleton) %>% unlist(), Model)
# saveRDS(HLPrior_samples, paste0(dir, 'IIDIPUS_Input/HLPriorSamples_MCMCOut'))

#------------------ Prediction updating Myanmar using Microsoft building damage -----------------------------

ODDyWithImpact <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/EQ20250328MMR_-1_WithBuildings_WithImpactSummaries')
sampled_full <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp/EQ20250328MMR_-1_WithBuildings_fullSampledImpact')

mandalay_dam = st_read("/home/manderso/Downloads/skysat_20250329_080228_predictions_merged.gpkg")
bbox_obsdam = project(ext(mandalay_dam), from=crs(mandalay_dam), to="EPSG:4326")

overlapping_cells = cells(ODDyWithImpact, bbox_obsdam)
plot(t(sampled_full[overlapping_cells[c(1,2)],3,]))

ODD_cropped <- crop(ODDyWithImpact, bbox_obsdam)

project(mandalay_dam, from=crs(mandalay_dam), to="EPSG:4326")
mandalay_dam_v <- project(vect(mandalay_dam),crs(ODD_cropped))
mandalay_dam_v$cond_1 <- as.numeric(mandalay_dam_v$damage_pct_0m > 0 & mandalay_dam_v$unknown_pct < .50)
mandalay_dam_v$cond_2 <- as.numeric(mandalay_dam_v$unknown_pct > .50)

dam_agg = data.frame(
  cell = numeric(),
  n_buildings = numeric(),
  damage_count = numeric(),
  unknown_count = numeric()
)

for (cell_num in overlapping_cells){

  xy <- xyFromCell(ODDyWithImpact, cell=cell_num)
  
  # Cell resolution
  res_x <- res(ODDyWithImpact)[1]
  res_y <- res(ODDyWithImpact)[2]
  
  # Calculate extent: xmin, xmax, ymin, ymax
  cell_ext <- ext(
    xy[1] - res_x / 2,
    xy[1] + res_x / 2,
    xy[2] - res_y / 2,
    xy[2] + res_y / 2
  )
  
  cropped_builds = crop(mandalay_dam_v, cell_ext)
  dam_agg %<>% add_row(
    cell = cell_num,
    n_buildings = length(cropped_builds$cond_1),
    damage_count = sum(cropped_builds$cond_1, na.rm=T),
    unknown_count = sum(cropped_builds$cond_2, na.rm=T))
  
}
cells_plot = c(15,6)
plot(t(sampled_full[overlapping_cells[cells_plot],3,]/values(ODDyWithImpact$nBuildings)[overlapping_cells[cells_plot]]),
     xlab='Prop Dam Cell 1',
     ylab='Prop Dam Cell 2')
true_prop = dam_agg[cells_plot,'damage_count']/dam_agg[cells_plot,'n_buildings']
points(ifelse(is.na(true_prop[1]), 0, true_prop[1]), ifelse(is.na(true_prop[2]), 0, true_prop[2]), col='red', pch=19)

plot(colSums(sampled_full[overlapping_cells, 3, ])/sum(values(ODDyWithImpact$nBuildings)[overlapping_cells], na.rm=T), 
     colSums(sampled_full[,3,]),
     xlab='Prop Dam Mandalay', ylab='Total Buildings Damaged')
abline(v=sum(dam_agg$damage_count)/sum(dam_agg$n_buildings), col='red')

plot(colSums(sampled_full[overlapping_cells, 3, ])/sum(values(ODDyWithImpact$nBuildings)[overlapping_cells], na.rm=T), 
     colSums(sampled_full[,1,]),
     xlab='Prop Dam Mandalay', ylab='Total Population Displacement')
abline(v=sum(dam_agg$damage_count)/sum(dam_agg$n_buildings), col='red')

ODDyWithImpact <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/EQ20250328MMR_-1_WithBuildings_WithImpactSummaries')
sampled_full <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp/EQ20250328MMR_-1_WithBuildings_fullSampledImpact')


#-----------------------------------

naypyidaw_dam = st_read("/home/manderso/Downloads/naypyidaw_myanmar_results_03313025_skysat.gpkg")
bbox_loc2 = project(ext(naypyidaw_dam), from=crs(naypyidaw_dam), to="EPSG:4326")

overlapping_cells2 = cells(ODDyWithImpact, bbox_loc2)
plot(t(sampled_full[overlapping_cells2[c(1,2)],3,]))

ODD_cropped2 <- crop(ODDyWithImpact, bbox_loc2)

project(naypyidaw_dam, from=crs(naypyidaw_dam), to="EPSG:4326")
naypyidaw_dam_v <- project(vect(naypyidaw_dam),crs(ODD_cropped))
naypyidaw_dam_v$cond_1 <- as.numeric(naypyidaw_dam_v$damage_pct_0m > 0 & naypyidaw_dam_v$unknown_pct < .50)
naypyidaw_dam_v$cond_2 <- as.numeric(naypyidaw_dam_v$unknown_pct > .50)

dam_agg2 = data.frame(
  cell = numeric(),
  n_buildings = numeric(),
  damage_count = numeric(),
  unknown_count = numeric()
)

for (cell_num in overlapping_cells2){
  
  xy <- xyFromCell(ODDyWithImpact, cell=cell_num)
  
  # Cell resolution
  res_x <- res(ODDyWithImpact)[1]
  res_y <- res(ODDyWithImpact)[2]
  
  # Calculate extent: xmin, xmax, ymin, ymax
  cell_ext <- ext(
    xy[1] - res_x / 2,
    xy[1] + res_x / 2,
    xy[2] - res_y / 2,
    xy[2] + res_y / 2
  )
  
  cropped_builds = crop(naypyidaw_dam_v, cell_ext)
  dam_agg2 %<>% add_row(
    cell = cell_num,
    n_buildings = length(cropped_builds$cond_1),
    damage_count = sum(cropped_builds$cond_1, na.rm=T),
    unknown_count = sum(cropped_builds$cond_2, na.rm=T))
  
}

plot(colSums(sampled_full[overlapping_cells, 3,])/sum(values(ODDyWithImpact$nBuildings)[overlapping_cells]), 
     colSums(sampled_full[overlapping_cells2, 3,])/sum(values(ODDyWithImpact$nBuildings)[overlapping_cells2]), 
     xlab='Mandalay Dam Prop', ylab='Naypyidaw Dam Prop')
points(sum(dam_agg$damage_count)/sum(dam_agg$n_buildings), sum(dam_agg2$damage_count)/sum(dam_agg2$n_buildings),
     col='red', pch=19)


