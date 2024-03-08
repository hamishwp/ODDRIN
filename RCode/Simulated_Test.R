
#Run SMC-ABC algorithm with Simulated Data with features that mimic the real data

#Choose true Omega for Simulated Data
Omega_true <- list(Lambda1 = list(nu=9.6, kappa=1.9049508),
                   Lambda2 = list(nu=9.8626452, kappa=1.9049508),
                   Lambda3 = list(nu=9.3, kappa=0.6),
                   Lambda4 = list(nu=9.9, kappa=1.6),
                   theta= list(theta1=0.6),
                   eps = list(local=2.2292053, hazard_mort=0.8383464, hazard_disp=0.9, hazard_bd=0.9, hazard_cor=0.55),
                   vuln_coeff = list(PDens=0.05, SHDI=-0.5, GNIc=-0.1, Vs30=0.1, EQFreq=-0.1, FirstHaz=0.05, Night=0.05, FirstHaz.Night=0.1),
                   check = list(check=0.5))
  
#Simulate observed data:
n_events <-200
set.seed(3)
means <- exp(seq(-4,2.3, length.out=n_events))
data_y <- rbind(data.frame(polygon=1, impact='mortality', sampled=NA,
                      iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                      build_type=NA, inferred=F, event_id = 1:length(means),
                      observed = rlnorm(length(means), 10*Omega_true$vuln_coeff$Vs30* (means-mean(means)) + mean(means), 5*Omega_true$Lambda3$kappa)))
                      #observed = abs(rnorm(1, 100 * Omega_true$vuln_coeff$Vs30 * (means[i]-mean(means))+1500, Omega_true$Lambda3$kappa*500))) #round(rt(1,Omega_true$Lambda1$nu,5)*Omega_true$Lambda1$kappa))

#Remove Higher Level Priors
Model$HighLevelPriors <- function(Omega, Model, modifier=NULL){
  return(0)
}

#Replace SampleImpact with appropriate function for this model
prior_tightening = 1 # improve the prior so that the algorithm reaches stopping point / particle degeneracy faster (doesn't actually change prior but the meaning of the parameter)
                       # Closer to 0 = tighter prior
                       # Closer to 1 = closer to the original (relatively uninformative) prior
SampleImpact <- function(dir, Model, proposed, AlgoParams){
  impact_sample <- list()
  for (i in 1:(AlgoParams$Np*AlgoParams$m_CRPS)){
    impact_sample[[i]] <- data_y
    impact_sample[[i]]$sampled <-  rlnorm(length(means),  10* (prior_tightening*(proposed$vuln_coeff$Vs30-Omega_true$vuln_coeff$Vs30)+Omega_true$vuln_coeff$Vs30)* (means-mean(means)) + mean(means), 
                                          5*(prior_tightening*(proposed$Lambda3$kappa-Omega_true$Lambda3$kappa)+Omega_true$Lambda3$kappa))
  }
  return(list(poly=impact_sample))
}

#Plot data:
impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 100))

plot_impact_sample <- function(impact_sample, impact_filter='mortality'){
  sampled_df <- impact_sample$poly[[1]]$sampled
  for (i in 2:length(impact_sample$poly)){
    sampled_df <- cbind(sampled_df, impact_sample$poly[[i]]$sampled)
  }
  
  plot_data_jitter <- sampled_df
  jitter_val <- function(x){
    return_arr <- c()
    for (x_i in x){
      if(x_i<0.7){return_arr <- c(return_arr, runif(1,0.2, 0.3))}
      else {return_arr <- c(return_arr, runif(1, x_i-0.5, x+0.5))}
    }
    return(return_arr)
  }
  #df_poly_jitter[which(df_poly_jitter==0, arr.ind=T)] = 0.1
  plot_data_jitter <- apply(plot_data_jitter,1:2,jitter_val)
  
  plot_data <- data.frame(observed = sapply(impact_sample$poly[[1]]$observed, jitter_val),
                          sampled_median = apply(plot_data_jitter, 1, median), 
                          sampled_mean = apply(plot_data_jitter, 1, mean), 
                          sampled_quant2.5 = apply(plot_data_jitter, 1, quantile, 0.025),
                          sampled_quant97.5 = apply(plot_data_jitter, 1, quantile, 0.975),
                          impact_type = impact_sample$poly[[1]]$impact)
  
  p <- ggplot(plot_data %>% filter(impact_type==impact_filter), mapping=aes(x=observed, y=sampled_median, ymin=sampled_quant2.5, ymax=sampled_quant97.5)) + 
    geom_errorbar() + geom_point(col='red') + 
    scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    scale_y_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    #geom_pointrange(aes(col=train_flag)) + 
    geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1)
  return(p)
}
plot_impact_sample(impact_sample)

#Fit the model using ABC-SMC Algorithm:
AlgoParams$smc_steps <- 100
AlgoParams$smc_Npart <- 500
AlgoParams$m_CRPS <- 10
AlgoParams$smc_alpha <- 0.9
AlgoParams$rel_weightings <- c(0,1)
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulated1D_v3_crps_noplus10_200events_M10')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F, tag_notes=tag_notes)
AlgoParams$m_CRPS <- 5
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulated1D_v3_crps_noplus10_200events_M5')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F, tag_notes=tag_notes)

AlgoParams$m_CRPS <- 20
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulated1D_v3_crps_withplus10_200events_M20')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F, tag_notes=tag_notes)
AlgoParams$m_CRPS <- 25
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulated1D_v3_crps_noplus10_200events_M25')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F, tag_notes=tag_notes)

AlgoParams$rel_weightings <- c(1,10)
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulated1D_v3_WM1WV10_500events')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)
AlgoParams$rel_weightings <- c(1,100)
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulated1D_v3_WM1WV100')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)

#tag key: 
#1D = 1 dimensional (i.e. single impact type of mortality)
#v2 = more closely resembling real data, more zeros , v3 = even further resembling real data, larger errors
#WM = Weighting median component, WV = Weighting variance component

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_start_step_2024-02-27_121727_alpha0.9_simulated1D_v3_crps_noplus10_200events_M25')

library(rstanarm)
glm_out <- log(data_y$observed)-mean(means)
glm_in <- means-mean(means)
fit_stan <- stan_glm(glm_out ~ glm_in - 1, 
                     prior=NULL)


fit_samples <- as.matrix(fit_stan)
plot(AlgoResults$Omega_sample_phys[,18,1], AlgoResults$Omega_sample_phys[,6,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),18,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),6,s_finish], col='green')
points(fit_samples[1:500,1]/10, fit_samples[1:500,2]/5, col='blue')
points(0.1, 0.6, col='red', pch=19)
plot(density(AlgoResults$Omega_sample_phys[,18,s_finish]), xlab='Mean Parameter', col='blue')
lines(density(fit_samples[1:500,1]/10), col='green')
abline(v=0.1,col='red')

Omega <- AlgoResults$Omega_sample_phys[1,,140] %>% relist(skeleton=Model$skeleton)
impact_sample <- SampleImpact(dir, Model, Omega, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 100))
p1 <- plot_impact_sample(impact_sample)
impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 100))
p2 <- plot_impact_sample(impact_sample)

grid.arrange(p1, p2, nrow = 2)

plot(density(AlgoResults$Omega_sample_phys[,6,s_finish]))
lines(density(fit_samples[1:500,2]/5), col='blue')
lines(density(AlgoResults2$Omega_sample_phys[,6,120]), col='blue')
lines(density(AlgoResults3$Omega_sample_phys[,6,150]), col='green')

AlgoResultsM5 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-04_200737_alpha0.9_simulated1D_v3_crps_noplus10_200events_M5')
AlgoResultsM10 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-04_191934_alpha0.9_simulated1D_v3_crps_noplus10_200events_M10')
AlgoResultsM15 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-04_095738_alpha0.9_simulated1D_v3_crps_noplus10_200events_M15')
AlgoResultsM20 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-04_110930_alpha0.9_simulated1D_v3_crps_noplus10_200events_M20')
AlgoResultsM25 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-04_122602_alpha0.9_simulated1D_v3_crps_noplus10_200events_M25')

plot(density(AlgoResultsM15$Omega_sample_phys[,6,50]),col= 'red')
lines(density(AlgoResultsM20$Omega_sample_phys[,6,50]), col='purple')
lines(density(AlgoResultsM25$Omega_sample_phys[,6,50]), col='blue')
lines(density(AlgoResultsM10$Omega_sample_phys[,6,50]), col='orange')
lines(density(AlgoResultsM5$Omega_sample_phys[,6,50]), col='yellow')
abline(v=0.6)

# -------------------------------------------------------------------------------------------------

#Compare 3 different relative weightings and the resulting posterior:
AlgoResults_WM1_WV1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-02-19_180046_alpha0.95_simulated1D_v3_WM1WV1_500events')
AlgoResults_WM1_WV10 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-02-19_191506_alpha0.95_simulated1D_v3_WM1WV10_500events')
plot(AlgoResults_WM1_WV1$Omega_sample_phys[,18,1], AlgoResults_WM1_WV1$Omega_sample_phys[,6,1], col='darkgrey', xlab='Location parameter', ylab='Variance Parameter')
points(AlgoResults_WM1_WV10$Omega_sample_phys[,18,100], AlgoResults_WM1_WV10$Omega_sample_phys[,6,100], col='orange')
points(AlgoResults_WM1_WV1$Omega_sample_phys[,18,100], AlgoResults_WM1_WV1$Omega_sample_phys[,6,100], col='red')
points(fit_samples[1:300,1]/10, fit_samples[1:300,2]/5, col='blue')
points(0.1, 0.6, col='green', pch=19)
  
cr <- colorRamp(c("red", "blue"))
cols <- rgb(cr(AlgoResults_WM1_WV1$d_full[,1,1,60] / max(AlgoResults_WM1_WV1$d_full[,1,1,60])), max=255)
plot(AlgoResults_WM1_WV1$Omega_sample_phys[,6,60], AlgoResults_WM1_WV1$Omega_sample_phys[,18,60], col=cols, pch=19)

par_i <- 18
plot(density(AlgoResults_WM1_WV1$Omega_sample_phys[,par_i,100]), col='red', ylim=c(0,80), xlab='Posterior for variance parameter', main='')
lines(density(AlgoResults_WM1_WV10$Omega_sample_phys[,par_i,140]), col='orange') # run for longer to see if this becomes like the others
lines(density(fit_samples[1:300,1]/10), col='blue')
lines(density(AlgoResults_WM1_WV1$Omega_sample_phys[,par_i,1]), col='grey')
abline(v=0.1)

plot(AlgoResults$Omega_sample_phys[,6,100], AlgoResults$d_full[,1,4,100])
Omega_fake <- AlgoResults$Omega_sample_phys[which.min(AlgoResults$d_full[,1,4,100]), ,100] %>% relist(skeleton=Model$skeleton)

d_true <- c()
d_fake <- c()
for (i in 1:100){
  M <- 5
  impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M))
  dist <- CalcDist(impact_sample, AlgoParams  %>% replace(which(names(AlgoParams)==c('m_CRPS')), M))
  d_true <- c(d_true, dist[1,4])
  impact_sample <- SampleImpact(dir, Model, Omega_fake, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M))
  dist <- CalcDist(impact_sample, AlgoParams  %>% replace(which(names(AlgoParams)==c('m_CRPS')), M))
  d_fake <- c(d_fake, dist[1,4])
}

unif_true <- c()
for(i in 1:5000){
  unif_true <- c(unif_true, AndersonDarlingTest(runif(500), null='punif')$statistic)
}

plot(d_true - d_fake)
plot(density(d_true/16), col='red')
lines(density(d_fake/16), col='blue')
lines(density(unif_true), col='grey')




# -------------------------------------------------------------------------------------------------

# As size of mean component increases (e.g. by more noisy data), we need weighting of variance term to decrease. 
# e.g : 
AlgoResults_WM1_WV1_v2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-02-14_113905_alpha0.9_simulated1D_v2_WM1WV1')
AlgoResults_WM1_WV1_v3 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-02-14_161454_alpha0.9_simulated1D_v3_WM1WV1')

# v3 with much noiser data than v2
#Admittedly though, the data is different so it's not completely clear that this isn't just due to variation in the data
plot(AlgoResults_WM1_WV1_v3$Omega_sample_phys[,18,50], AlgoResults_WM1_WV1_v3$Omega_sample_phys[,6,50], col='blue', ylim=c(0.38, 0.68))
points(AlgoResults_WM1_WV1_v2$Omega_sample_phys[,18,50], AlgoResults_WM1_WV1_v2$Omega_sample_phys[,6,50])
points(0.1, 0.6, col='red', pch=19)

# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- MULTIVARIATE MODEL ----------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
n_events <- 200
set.seed(1)

# so I can avoid changing the actual priors:
mean_coeffs <- c(10, 8, 13)
var_coeffs <- c(5, 1.2, 1.6)

means_mort <- exp(seq(-4,2.3, length.out=n_events))
mort_impacts <- data.frame(polygon=j, impact='mortality', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = i,
                           observed = rlnorm(length(means_mort), mean_coeffs[1]*Omega_true$vuln_coeff$Vs30* (means_mort-mean(means_mort)) + mean(means_mort), var_coeffs[1]*Omega_true$Lambda3$kappa))

means_disp <- exp(seq(-4,2.6, length.out=n_events))
disp_impacts <- data.frame(polygon=j, impact='displacement', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = i,
                           observed = rlnorm(length(means_disp), mean_coeffs[2]*Omega_true$vuln_coeff$SHDI* (means_disp-mean(means_disp)) + mean(means_disp), var_coeffs[2]*Omega_true$Lambda2$kappa))

means_bd <- exp(seq(-4,2.8, length.out=n_events))
bd_impacts <- data.frame(polygon=j, impact='buildDam', sampled=NA,
                         iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                         build_type=NA, inferred=F, event_id = i,
                         observed = rlnorm(length(means_bd), mean_coeffs[3]*Omega_true$vuln_coeff$PDens* (means_bd-mean(means_bd)) + mean(means_bd), var_coeffs[3]*Omega_true$Lambda4$kappa))

data_y <- rbind(mort_impacts, disp_impacts, bd_impacts)
plot(data_y$observed[which(data_y$impact=='buildDam')])

prior_tightening <- 1
SampleImpact <- function(dir, Model, proposed, AlgoParams){
  impact_sample <- list()
  for (i in 1:(AlgoParams$Np*AlgoParams$m_CRPS)){
    impact_sample[[i]] <- data_y
    #impact_sample[[i]]$sampled <-  abs(rnorm(length(means), 100*proposed$vuln_coeff$Vs30 * (means-mean(means))+1500, proposed$Lambda3$kappa*500)) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='mortality')] <-  rlnorm(length(means_mort),  mean_coeffs[1]*(prior_tightening*(proposed$vuln_coeff$Vs30-Omega_true$vuln_coeff$Vs30)+Omega_true$vuln_coeff$Vs30) * (means_mort-mean(means_mort)) + mean(means_mort), var_coeffs[1]*proposed$Lambda3$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='displacement')] <-  rlnorm(length(means_disp),  mean_coeffs[2]*(prior_tightening*(proposed$vuln_coeff$SHDI-Omega_true$vuln_coeff$SHDI)+Omega_true$vuln_coeff$SHDI) * (means_disp-mean(means_disp)) + mean(means_disp), var_coeffs[2]*proposed$Lambda2$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='buildDam')] <-  rlnorm(length(means_bd), mean_coeffs[3]*(prior_tightening*(proposed$vuln_coeff$PDens-Omega_true$vuln_coeff$PDens)+Omega_true$vuln_coeff$PDens) * (means_bd-mean(means_bd)) + mean(means_bd), var_coeffs[3]*proposed$Lambda4$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
  }
  return(list(poly=impact_sample))
}

impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 100))
plot_impact_sample(impact_sample, impact_filter='buildDam')

AlgoParams$smc_steps <- 120
AlgoParams$smc_Npart <- 500
AlgoParams$m_CRPS <- 5
AlgoParams$smc_alpha <- 0.95
AlgoParams$rel_weightings <- c(1,1)
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulated3D_v3_WM1WV1')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)


# -----------------------------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- MULTIVARIATE MODEL WITH CORRELATION -----------------------------------
# -------------------------------------------------------------------------------------------------------------
n_events <- 200
set.seed(1)

# so I can avoid changing the actual priors:
mean_coeffs <- c(10, 8, 13)
var_coeffs <- var_coeffs <- c(2, 0.5, 0.5) #c(5, 1.2, 1.6)
means_mort <- exp(seq(-4,2.3, length.out=n_events))
means_disp <- exp(seq(-4,2.6, length.out=n_events))
means_bd <- exp(seq(-4,2.8, length.out=n_events))

cov_mort_disp = Omega_true$eps$hazard_cor * var_coeffs[1]*Omega_true$Lambda3$kappa * var_coeffs[2]*Omega_true$Lambda2$kappa
cov_mort_bd = Omega_true$eps$hazard_cor * var_coeffs[1]*Omega_true$Lambda3$kappa * var_coeffs[3]*Omega_true$Lambda4$kappa
cov_disp_bd = Omega_true$eps$hazard_cor * var_coeffs[2]*Omega_true$Lambda2$kappa * var_coeffs[3]*Omega_true$Lambda4$kappa
covar_matrix = cbind(c((var_coeffs[1]*Omega_true$Lambda3$kappa)^2, cov_mort_disp, cov_mort_bd), c(0, (var_coeffs[2]*Omega_true$Lambda2$kappa)^2, cov_disp_bd), c(0, 0, (var_coeffs[3]*Omega_true$Lambda4$kappa)^2))
covar_matrix[upper.tri(covar_matrix)] = covar_matrix[lower.tri(covar_matrix)]

covar_matrix <- diag(c(var_coeffs[1]*Omega_true$Lambda3$kappa, var_coeffs[2]*Omega_true$Lambda2$kappa, var_coeffs[3]*Omega_true$Lambda4$kappa))

observed_multivar <- exp(rmvnorm(length(means_mort), rep(0, 3), covar_matrix) + cbind(
                          mean_coeffs[1]*Omega_true$vuln_coeff$Vs30* (means_mort-mean(means_mort)) + mean(means_mort), 
                          mean_coeffs[2]*Omega_true$vuln_coeff$SHDI* (means_disp-mean(means_disp)) + mean(means_disp),
                          mean_coeffs[3]*Omega_true$vuln_coeff$PDens* (means_bd-mean(means_bd)) + mean(means_bd)))

mort_impacts <- data.frame(polygon=j, impact='mortality', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = i,
                           observed = observed_multivar[,1])


disp_impacts <- data.frame(polygon=j, impact='displacement', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = i,
                           observed = observed_multivar[,2])

bd_impacts <- data.frame(polygon=j, impact='buildDam', sampled=NA,
                         iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                         build_type=NA, inferred=F, event_id = i,
                         observed = observed_multivar[,3])

data_y <- rbind(mort_impacts, disp_impacts, bd_impacts)
plot(data_y$observed[which(data_y$impact=='buildDam')])

prior_tightening <- 0.1
SampleImpact <- function(dir, Model, proposed, AlgoParams){
  impact_sample <- list()
  for (i in 1:(AlgoParams$Np*AlgoParams$m_CRPS)){
    impact_sample[[i]] <- data_y
    #impact_sample[[i]]$sampled <-  abs(rnorm(length(means), 100*proposed$vuln_coeff$Vs30 * (means-mean(means))+1500, proposed$Lambda3$kappa*500)) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
   
    #cov_mort_disp = proposed$eps$hazard_cor * var_coeffs[1]*proposed$Lambda3$kappa * var_coeffs[2]*proposed$Lambda2$kappa
    #cov_mort_bd = proposed$eps$hazard_cor * var_coeffs[1]*proposed$Lambda3$kappa * var_coeffs[3]*proposed$Lambda4$kappa
    #cov_disp_bd = proposed$eps$hazard_cor * var_coeffs[2]*proposed$Lambda2$kappa * var_coeffs[3]*proposed$Lambda4$kappa
    #covar_matrix = cbind(c((var_coeffs[1]*proposed$Lambda3$kappa)^2, cov_mort_disp, cov_mort_bd), c(0, (var_coeffs[2]*proposed$Lambda2$kappa)^2, cov_disp_bd), c(0, 0, (var_coeffs[3]*proposed$Lambda4$kappa)^2))
    #covar_matrix[upper.tri(covar_matrix)] = covar_matrix[lower.tri(covar_matrix)]
    
    covar_matrix <- diag(c(var_coeffs[1]*proposed$Lambda3$kappa, var_coeffs[2]*proposed$Lambda2$kappa, var_coeffs[3]*proposed$Lambda4$kappa))

    means_mort_upd <- mean_coeffs[1]*(prior_tightening*(proposed$vuln_coeff$Vs30-Omega_true$vuln_coeff$Vs30)+Omega_true$vuln_coeff$Vs30) * (means_mort-mean(means_mort)) + mean(means_mort)
    means_disp_upd <- mean_coeffs[2]*(prior_tightening*(proposed$vuln_coeff$SHDI-Omega_true$vuln_coeff$SHDI)+Omega_true$vuln_coeff$SHDI) * (means_disp-mean(means_disp)) + mean(means_disp)
    means_bd_upd <- mean_coeffs[3]*(prior_tightening*(proposed$vuln_coeff$PDens-Omega_true$vuln_coeff$PDens)+Omega_true$vuln_coeff$PDens) * (means_bd-mean(means_bd)) + mean(means_bd)
    
    sampled_multivar <- exp(rmvnorm(length(means_mort), rep(0, 3), covar_matrix) + cbind(means_mort_upd, means_disp_upd, means_bd_upd))
    
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='mortality')] <-  sampled_multivar[,1] #rlnorm(length(means_mort),  mean_coeffs[1]*(prior_tightening*(proposed$vuln_coeff$Vs30-Omega_true$vuln_coeff$Vs30)+Omega_true$vuln_coeff$Vs30) * (means_mort-mean(means_mort)) + mean(means_mort), var_coeffs[1]*proposed$Lambda3$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='displacement')] <-  sampled_multivar[,2] #rlnorm(length(means_disp),  mean_coeffs[2]*(prior_tightening*(proposed$vuln_coeff$SHDI-Omega_true$vuln_coeff$SHDI)+Omega_true$vuln_coeff$SHDI) * (means_disp-mean(means_disp)) + mean(means_disp), var_coeffs[2]*proposed$Lambda2$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='buildDam')] <- sampled_multivar[,3] #rlnorm(length(means_bd), mean_coeffs[3]*(prior_tightening*(proposed$vuln_coeff$PDens-Omega_true$vuln_coeff$PDens)+Omega_true$vuln_coeff$PDens) * (means_bd-mean(means_bd)) + mean(means_bd), var_coeffs[3]*proposed$Lambda4$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
  }
  return(list(poly=impact_sample))
}

impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 100))
plot_impact_sample(impact_sample, impact_filter='buildDam')

AlgoParams$smc_steps <- 160
AlgoParams$smc_Npart <- 500
AlgoParams$m_CRPS <- 20
AlgoParams$smc_alpha <- 0.9
AlgoParams$rel_weightings <- c(0,1)
AlgoParams$kernel_sd$mortality <- 1
AlgoParams$kernel_sd$buildDam <- 1
AlgoParams$kernel_sd$buildDest <- 1
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_simulated3D_withcorr_noplus10_smallvar_CRPS_M20_equalimpactweightings_priortightening')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)

plot(AlgoResults$Omega_sample_phys[,18,1], AlgoResults$Omega_sample_phys[,6,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),18,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),6,s_finish], col='green')
points(0.1, 0.6, col='red', pch=19)

hist(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),6,s_finish])

plot(AlgoResults$Omega_sample_phys[,16,1], AlgoResults$Omega_sample_phys[,4,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),16,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),4,s_finish], col='green')
points(-0.5, 1.9, col='red', pch=19)

plot(AlgoResults$Omega_sample_phys[,15,1], AlgoResults$Omega_sample_phys[,8,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),15,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),8,s_finish], col='green')
points(0.05, 1.6, col='red', pch=19)

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-06_120756_alpha0.9_simulated3D_withcorr_noplus10_smallvar_CRPS_M20')

# par_i <- 18
# s_finish <- 160
# plot(AlgoResults$Omega_sample_phys[, par_i, 135], AlgoResults$d[,1,135])
# 
# library(rstanarm)
# glm_out <- log(data_y$observed[which(data_y$impact=='mortality')])-mean(means_mort)
# glm_in <- means_mort-mean(means_mort)
# fit_stan <- stan_glm(glm_out ~ glm_in - 1, 
#                      prior=NULL)
# 
# fit_samples <- as.matrix(fit_stan)
# 
# 
# abline(v=0.6, col='red')

# -----------------------------------------------------------------------------------------------------------------------------------------



AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-02-18_092654_alpha0.95_simulated3D_v3_WM1WV1')
s_finish <- 120

#mortality:
plot(AlgoResults$Omega_sample_phys[,18,1], AlgoResults$Omega_sample_phys[,6,1], xlab='Location Param', ylab='Variance Param')
points(AlgoResults$Omega_sample_phys[,18,s_finish], AlgoResults$Omega_sample_phys[,6,s_finish], col='blue')
points(Omega_true$vuln_coeff$Vs30, Omega_true$Lambda3$kappa, col='red', pch=19)

#displacement:
plot(AlgoResults$Omega_sample_phys[,16,1], AlgoResults$Omega_sample_phys[,4,1], xlab='Location Param', ylab='Variance Param')
points(AlgoResults$Omega_sample_phys[,16,s_finish], AlgoResults$Omega_sample_phys[,4,s_finish], col='blue')
points(Omega_true$vuln_coeff$SHDI, Omega_true$Lambda2$kappa, col='red', pch=19)

#building damage:
plot(AlgoResults$Omega_sample_phys[,15,1], AlgoResults$Omega_sample_phys[,8,1], xlab='Location Param', ylab='Variance Param')
points(AlgoResults$Omega_sample_phys[,15,s_finish], AlgoResults$Omega_sample_phys[,8,s_finish], col='blue')
points(Omega_true$vuln_coeff$PDens, Omega_true$Lambda4$kappa, col='red', pch=19)

#contributions to distance function:
plot(AlgoResults$Omega_sample_phys[,8,s_finish], AlgoResults$d_full[,1,6,s_finish] * AlgoParams$rel_weightings[1])
points(AlgoResults$Omega_sample_phys[,8,s_finish], AlgoResults$d_full[,1,3,s_finish] * AlgoParams$rel_weightings[2], col='red')
points(AlgoResults$Omega_sample_phys[,8,s_finish], AlgoResults$d_full[,1,6,s_finish] * AlgoParams$rel_weightings[1] + 
                                                   AlgoResults$d_full[,1,3,s_finish] * AlgoParams$rel_weightings[2], col='blue', pch=19)

plot(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d_full[,1,4,s_finish] * AlgoParams$rel_weightings[1])
points(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d_full[,1,1,s_finish] * AlgoParams$rel_weightings[2], col='red')
points(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d_full[,1,1,s_finish] * AlgoParams$rel_weightings[1] + 
         AlgoResults$d_full[,1,4,s_finish] * AlgoParams$rel_weightings[2], col='blue', pch=19)

#check correlation between them? 





AlgoResults$propCOV[6,6,2]
AlgoResults$propCOV[18,18,2]

AlgoResults$propCOV[6,6,2]
AlgoResults$propCOV[18,18,2]
AlgoResults$propCOV[6,18,2]
s_finish <-50
plot(AlgoResults$Omega_sample_phys[,18,1], AlgoResults$Omega_sample_phys[,6,1], col='blue')
points(AlgoResults$Omega_sample_phys[,18,s_finish], AlgoResults$Omega_sample_phys[,6,s_finish])
points(Omega_true$vuln_coeff$Vs30, Omega_true$Lambda3$kappa, col='red',pch=19)


s_finish <-45
plot(AlgoResults$Omega_sample_phys[,16,1], AlgoResults$Omega_sample_phys[,4,1], col='blue')
points(AlgoResults$Omega_sample_phys[,16,s_finish], AlgoResults$Omega_sample_phys[,4,s_finish])
points(Omega_true$vuln_coeff$SHDI, Omega_true$Lambda2$kappa, col='red',pch=19)

plot(AlgoResults$Omega_sample_phys[,18,s_finish], AlgoResults$d[,1,s_finish])
plot(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d[,1,s_finish])
plot(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d_full[,1,4,s_finish], col='blue')
plot(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d_full[,1,1,s_finish]* AlgoResults$rel_weights[1], col='red')
plot(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d_full[,1,4,s_finish], col='blue')
plot(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d_full[,1,4,s_finish] * 41, col='blue')

s_plot <- 50
plot(AlgoResults$Omega_sample_phys[,18,s_plot], AlgoResults$d_full[,1,1,s_plot] * AlgoResults$rel_weights[1]/16, col='red', ylim=c(0,8))
points(AlgoResults$Omega_sample_phys[,18,s_plot], AlgoResults$d_full[,1,4,s_plot] * AlgoResults$rel_weights[4]/16, col='blue')
points(AlgoResults$Omega_sample_phys[,18,s_plot], AlgoResults$d_full[,1,1,s_plot] * AlgoResults$rel_weights[1]/16 + 
         AlgoResults$d_full[,1,4,s_plot] * AlgoResults$rel_weights[4]/16, col='black')


y_obs <- 100
y_samp <- 1:1000
dist <- c()
for(x in y_samp){
  dist <- c(dist, abs((log(x+10)-log(y_obs+10))))
}
plot(y_samp, dist)

#points(AlgoResults$Omega_sample_phys[,18,20], AlgoResults$d[,1,20], col='black')


plot(AlgoResults$Omega_sample_phys[,6,s_finish], AlgoResults$d_full[,1,1,s_finish]* AlgoResults$rel_weights[1])
plot(AlgoResults$Omega_sample_phys[,18,s_finish], AlgoResults$d_full[,1,1,s_finish]* AlgoResults$rel_weights[1])
plot(AlgoResults$Omega_sample_phys[,4,s_finish], AlgoResults$d_full[,1,5,s_finish]* AlgoResults$rel_weights[5])
points(AlgoResults$Omega_sample_phys[,4,s_finish], AlgoResults$d_full[,1,2,s_finish]* AlgoResults$rel_weights[2], col='red')

plot(AlgoResults$d_full[,1,2,s_finish]* AlgoResults$rel_weights[2], col='red')


plot(AlgoResults$Omega_sample_phys[,6,1])
points(AlgoResults$Omega_sample_phys[,6,s_finish], col='red')
abline(h=Omega_true$Lambda3$kappa)

#Investigate the effect of different PropCOVs
unique_part <- c()
for (i in 1:120){
  unique_part <- c(unique_part, length(unique(AlgoResults$Omega_sample_phys[AlgoResults$W[,i]>0,1,i])))
}
plot(unique_part)
plot(apply(AlgoResults$Omega_sample_phys[,1,], 2, function(x) length(unique(x)))) # propCOV = propCOV/12
points(apply(AlgoResults2$Omega_sample_phys[,1,], 2, function(x) length(unique(x))), col='red') # propCOV = propCOV
points(apply(AlgoResults3$Omega_sample_phys[,1,], 2, function(x) length(unique(x))), col='red') # propCOV = propCOV /36

plot(density(AlgoResults2$Omega_sample_phys[,1,1]))
lines(density(AlgoResults2$Omega_sample_phys[,1,30]), col='red')

plot(AlgoResults$Omega_sample_phys[,6,20], AlgoResults$d_full[,1,1,20])
plot(AlgoResults$Omega_sample_phys[,6,1], AlgoResults$d_full[,1,4,1]/16)
plot(AlgoResults$Omega_sample_phys[,18,20], AlgoResults$d_full[,1,4,20]/16)
plot(AlgoResults$Omega_sample_phys[,18,20], AlgoResults$d_full[,1,1,20])
plot(AlgoResults$Omega_sample_phys[,18,20],AlgoResults$Omega_sample_phys[,6,20])
plot(AlgoResults$Omega_sample_phys[,18,20],AlgoResults$d[,1,20])


plot(AlgoResults$Omega_sample_phys[,1,30],AlgoResults$Omega_sample_phys[,2,30])
abline(v=0.9, col='red')


diffs_min1 <- c()
diffs_min2 <- c()
diffs_min3 <- c()
data_y$observed <- rlnorm(length(means), means, 0.7)
for (i in 1:1000){
  y_samp <- rlnorm(length(means), means, 0.9)
  y_samp2 <-  rlnorm(length(means), means, 0.7)
  y_samp3 <-  rlnorm(length(means), means, 0.5)
  for (m in 2:M){
    y_samp <- cbind(y_samp, rlnorm(length(means), means, 0.9))
    y_samp2 <- cbind(y_samp2, rlnorm(length(means), means, 0.7))
    y_samp3 <- cbind(y_samp3, rlnorm(length(means), means, 0.5))
  }
  diffs <- abs(sweep(y_samp, 1, data_y$observed))
  diffs_min1 <- c(diffs_min1, sum(apply(diffs, 1, min)))
  diffs <- abs(sweep(y_samp2, 1, data_y$observed))
  diffs_min2 <- c(diffs_min2, sum(apply(diffs, 1, min)))
  diffs <- abs(sweep(y_samp3, 1, data_y$observed))
  diffs_min3 <- c(diffs_min3, sum(apply(diffs, 1, min)))
}
plot(diffs_min1, ylim=c(min(diffs_min1, diffs_min2, diffs_min3), max(diffs_min1, diffs_min2, diffs_min3)))
points(diffs_min2, col='red')
points(diffs_min3, col='blue')
# ----------------------------------------- Justifying choice of relative weights ----------------------------------------------


plot(AlgoResults$Omega_sample_phys[,6,30], AlgoResults$d_full[,1,4,30], xlab='Variance parameter', ylab='Distance Component 2 (Quantile component)')
abline(v=0.9, col='red')

samps <- array(rnorm(100, Omega_true$Lambda1$nu, Omega_true$Lambda1$kappa), dim=c(5, 20))
plot(rep(1:NCOL(samps), each=NROW(samps)), samps)
points(1:NCOL(samps), apply(samps, 2, median), col='red', pch=19)
samps2 <- array(rnorm(100, 9.75, 1.6), dim=c(5, 20))
points(rep(1:NCOL(samps2), each=NROW(samps2)), samps2, col='blue')
points(1:NCOL(samps2), apply(samps2, 2, median), col='blue', pch=19)

#take a model run (not quite the final model nor algorithm but should be a good example still)
AlgoResults2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-01-17_214246')
AlgoResults2$sampled_full <- NULL

hist(AlgoResults$Omega_sample_phys[,2,30])
plot(AlgoResults$Omega_sample_phys[,1,30], AlgoResults$Omega_sample_phys[,2,30])

test_stats <- c()
for (i in 1:10000){
  test.stat <- AndersonDarlingTest(runif(10, 0,1), null='punif')$statistic
  test_stats <- c(test_stats, test.stat)
}
hist(test_stats)
mean(test_stats)

store <- c()
for (x in 1:500){
  obs <- rlnorm(1, 4, 0.75)
  samp <- apply(matrix(rlnorm(1500, 4, 0.75), nrow=300),1, median)
  store <- c(store, mean(abs(log(samp+10)-log(obs+10))))
}

P_unif = ecdf(test_stats)
P_unif_test <- function(x){
  p <- 1-P_unif(x)
  if (p==1) return (0.999999999)
  if (p==0) return (0.000000001)
  return(p)
}


hist(AlgoResults2$d_full[,1,4,40]/16, breaks=30)

mean(AlgoResults$d_full[,1,4,1])
hist(AlgoResults$d_full[,1,1,1], breaks=30)
mad(AlgoResults$d_full[order(AlgoResults$d_full[,1,1,1])[1:50],1,1,1]/15)

median_dist <- c()
var_dist <- c()
for (i in 1:100){
  xx <- SampleImpact(dir, Model, Omega_true, AlgoParams)
  dist_true <- CalcDist(xx, AlgoParams)
  median_dist <- c(median_dist, dist_true[1,1])
  var_dist <- c(var_dist, dist_true[1,4])
}
hist(median_dist)
hist(var_dist)
mad(median_dist)
mad(var_dist)

#simulate data like real data:
N_events <- 100
y_true <- rep(0, N_events)
means <- revd(100, 4.5, 0.85)
y_true <- rlnorm(100, means, 0.75)
#Try 1:1?

# -------------------------------------------------------------------------------------------------------------------------------

AlgoResults$rel_weights

mad(AlgoResults$d_full[,1,1,20])/mad(AlgoResults$d_full[,1,4,20])

mad(AlgoResults$d_full[order(AlgoResults$d_full[,1,1,1])[1:15],1,1,1])/mad(AlgoResults$d_full[order(AlgoResults$d_full[,1,4,1])[1:15],1,4,1])


ll_true <- c()
ll_fitted <- c()
for (j in 1:1000){
  xx <- SampleImpact(dir, Model, Omega_true, AlgoParams)
  ll_true <- c(ll_true, rowSums(CalcDist(xx, AlgoParams)))
  xx <- SampleImpact(dir, Model, AlgoResults$Omega_sample_phys[2,,20] %>% relist(skeleton=Model$skeleton), AlgoParams)
  ll_fitted <- c(ll_fitted, rowSums(CalcDist(xx, AlgoParams)))
}
plot(ll_true, col='blue')
points(ll_fitted, col='red')


median_dist <- c()
var_dist <- c()
for (j in 1:100){
  xx <- SampleImpact(dir, Model, Omega_true, AlgoParams)
  dist_true <- CalcDist(xx, AlgoParams)
  median_dist <- c(median_dist, dist_true[1,1])
  var_dist <- c(var_dist, dist_true[1,4])
}
mad(median_dist) / mad(var_dist)

plot(ll_true, col='blue')
points(ll_fitted, col='red')
  