
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
                                      Lambda2 = list(nu=11.7, kappa=0.75), #list(nu=10.65, kappa=1.5), #
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

# Real data results:

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-08-20_051627_alpha0.9_M60_Npart1000RealAgg5_propCOVmult0.2')
plot_correlated_posteriors(AlgoResults, pairings = rbind(c(1,2), c(3,4), c(5,6), c(10,11), c(12,13), c(11, 14), c(18,15), c(16,17), c(19,20), c(21,22), c(9,23)))
#RealPosteriors.pdf, 7 x 9 inches

AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_RealAgg5/'
df_postpredictive_sampled_best <- create_df_postpredictive(AlgoResults, single_particle=F, 
                                                           M=100, output='SampledTotal')

df_postpredictive_sampled_best <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/sampledBest18thSept')

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
  

p1 <- arrangeGrob(mort_test, top = textGrob("(a)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 14)))
p2 <- arrangeGrob(disp_test, top = textGrob("(b)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 14)))
p3 <- arrangeGrob(buildDam_test, top = textGrob("(c)", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 14)))

grid.arrange(p1, p2, p3, ncol=3)
#PostPredReal.pdf, 4 x 12.1 inches

AlgoResults %<>% addAlgoParams(s_finish=NULL)
for (i in c(16,)){
  p_H_given_y = sum(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish]>0)/length(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish])
  p_H = sum(AlgoResults$Omega_sample_phys[,i,1]>0)/length(AlgoResults$Omega_sample_phys[,i,1])
  print(paste('Bayes factor for', names(unlist(Omega))[i], ':', round(p_H_given_y * (1-p_H)/((1-p_H_given_y)*p_H), 3)))
  print(paste('Post. prob of greater than 0:', sum(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish]>0)/length(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish])))
}



#SMC:
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-10-01_223810_alphaAdaptivePropCOV0.2_M100_Npart1000SimInput')

#MCMC:
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-09-30_223103_MCMC_SimCorrVuln_Trial_LR40_Rho0.9_adaptive')

AlgoResults2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-10-08_213516_MCMC_SimCorrVuln_Trial_LR40_Rho0.9_adaptive')

AlgoResults2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-10-12_222810_MCMC_SimCorrVuln_Trial_LR40_Rho0.9_adaptive')
plot(AlgoResults2$loss, type='l', ylim=c(2.2, 5))
points(AlgoResults$loss, type='l', col='red')

par(mfrow=c(2,2))
params <- c(7,8)
xlim=range(AlgoResults2$Omega_sample_phys[params[1],], AlgoResults$Omega_sample_phys[params[1],],na.rm=T)
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

