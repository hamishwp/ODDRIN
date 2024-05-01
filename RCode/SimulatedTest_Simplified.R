
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
Omega_true <- list(Lambda1 = list(nu=1, kappa=0.6),
                   Lambda2 = list(nu=1.2, kappa=0.9),
                   Lambda3 = list(nu=-0.5, kappa=0.2),
                   Lambda4 = list(nu=9.9, kappa=1.6),
                   theta= list(theta1=0.6),
                   eps = list(local=2.2292053, hazard_mort=0.8383464, hazard_disp=0.9, hazard_bd=0.9, hazard_cor=0.55),
                   vuln_coeff = list(PDens=0.05, SHDI=-0.5, GNIc=-0.1, Vs30=0.1, EQFreq=-0.1, FirstHaz=0.05, Night=0.05, FirstHaz.Night=0.1),
                   check = list(check=0.5))

Model$Priors <- list( #All uniform so currently not included in the acceptance probability. 
  Lambda1=list(mu=list(dist='laplace', location=1, scale=0.1), 
               sigma=list(dist='unif', min=0, max=1)), 
  Lambda2=list(mu=list(dist='laplace', location=1.2, scale=0.1), 
               sigma=list(dist='unif', min=0, max=1.5)),
  Lambda3=list(mu=list(dist='laplace', location=-0.5, scale=0.1), 
               sigma=list(dist='unif', min=0, max=1)),
  Lambda4=list(mu=list(dist='unif', min=8, max=12.5), 
               sigma=list(dist='unif', min=0.25, max=2.5)),
  theta=list(theta1=list(dist='unif', min=0, max=1)),
  eps=list(local=list(dist='unif', min=0, max=1.5),
           hazard_mort=list(dist='unif', min=0, max=1.5),
           hazard_disp=list(dist='unif', min=0, max=1.5),
           hazard_bd=list(dist='unif', min=0, max=1.5),
           hazard_cor=list(dist='unif', min=0, max=1)),
  vuln_coeff=list(PDens=list(dist='laplace', location=0, scale=0.25),
                  EQFreq=list(dist='laplace', location=0, scale=0.25),
                  SHDI=list(dist='laplace', location=0, scale=0.25),
                  GNIc=list(dist='laplace', location=0, scale=0.25),
                  Vs30=list(dist='laplace', location=0, scale=0.25),
                  Mag=list(dist='laplace', location=0, scale=0.25),
                  FirstHaz=list(dist='laplace', location=0, scale=0.25),
                  Night=list(dist='laplace', location=0, scale=0.25),
                  FirstHaz.Night=list(dist='laplace', location=0, scale=0.25)),
  check=list(check=list(dist='unif', min=0, max=1))
)

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

Model$HighLevelPriors <- function(Omega, Model, modifier=NULL){
  return(0)
}


n_events <- 200
set.seed(2)

# so I can avoid changing the actual priors:
means_mort <- exp(seq(-4,2.3, length.out=n_events))
means_disp <- exp(seq(2,2.6, length.out=n_events))
means_bd <- exp(seq(-3,4, length.out=n_events))

cov_mort_disp = Omega_true$eps$hazard_cor * Omega_true$Lambda1$kappa * Omega_true$Lambda2$kappa
cov_mort_bd = Omega_true$eps$hazard_cor * Omega_true$Lambda1$kappa * Omega_true$Lambda3$kappa
cov_disp_bd = Omega_true$eps$hazard_cor * Omega_true$Lambda2$kappa * Omega_true$Lambda3$kappa
covar_matrix = cbind(c((Omega_true$Lambda1$kappa)^2, cov_mort_disp, cov_mort_bd), c(0, (Omega_true$Lambda2$kappa)^2, cov_disp_bd), c(0, 0, (Omega_true$Lambda3$kappa)^2))
covar_matrix[upper.tri(covar_matrix)] = covar_matrix[lower.tri(covar_matrix)]

#covar_matrix <- diag(c(Omega_true$Lambda1$kappa^2, Omega_true$Lambda2$kappa^2, Omega_true$Lambda3$kappa^2))

observed_multivar <- exp(rmvnorm(length(means_mort), rep(0, 3), covar_matrix) + cbind(
                          Omega_true$Lambda1$nu* (means_mort-mean(means_mort)) + mean(means_mort), 
                          Omega_true$Lambda2$nu* (means_disp-mean(means_disp)) + mean(means_disp),
                          Omega_true$Lambda3$nu* (means_bd-mean(means_bd)) + mean(means_bd)))

mort_impacts <- data.frame(polygon=1, impact='mortality', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = 1:n_events,
                           observed = observed_multivar[,1])

disp_impacts <- data.frame(polygon=1, impact='displacement', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = 1:n_events,
                           observed = observed_multivar[,2])

bd_impacts <- data.frame(polygon=1, impact='buildDam', sampled=NA,
                         iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                         build_type=NA, inferred=F, event_id = 1:n_events,
                         observed = observed_multivar[,3])

data_y <- rbind(mort_impacts, disp_impacts, bd_impacts)
plot(data_y$observed[which(data_y$impact=='buildDam')])

#prior_tightening <- 0.1
SampleImpact <- function(dir, Model, proposed, AlgoParams){
  impact_sample <- list()
  for (i in 1:(AlgoParams$Np*AlgoParams$m_CRPS)){
    impact_sample[[i]] <- data_y
    #impact_sample[[i]]$sampled <-  abs(rnorm(length(means), 100*proposed$vuln_coeff$Vs30 * (means-mean(means))+1500, proposed$Lambda3$kappa*500)) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
   
    cov_mort_disp = proposed$eps$hazard_cor * proposed$Lambda1$kappa * proposed$Lambda2$kappa
    cov_mort_bd = proposed$eps$hazard_cor * proposed$Lambda1$kappa * proposed$Lambda3$kappa
    cov_disp_bd = proposed$eps$hazard_cor * proposed$Lambda2$kappa * proposed$Lambda3$kappa
    covar_matrix = cbind(c((proposed$Lambda1$kappa)^2, cov_mort_disp, cov_mort_bd), c(0, (proposed$Lambda2$kappa)^2, cov_disp_bd), c(0, 0, (proposed$Lambda3$kappa)^2))
    covar_matrix[upper.tri(covar_matrix)] = covar_matrix[lower.tri(covar_matrix)]
    
    #covar_matrix <- diag(c(proposed$Lambda1$kappa^2, proposed$Lambda2$kappa^2, proposed$Lambda3$kappa^2))

    means_mort_upd <- proposed$Lambda1$nu * (means_mort-mean(means_mort)) + mean(means_mort)
    means_disp_upd <- proposed$Lambda2$nu * (means_disp-mean(means_disp)) + mean(means_disp)
    means_bd_upd <- proposed$Lambda3$nu * (means_bd-mean(means_bd)) + mean(means_bd)
    
    sampled_multivar <- exp(rmvnorm(length(means_mort), rep(0, 3), covar_matrix) + cbind(means_mort_upd, means_disp_upd, means_bd_upd))
    
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='mortality')] <-  sampled_multivar[,1] #rlnorm(length(means_mort),  mean_coeffs[1]*(prior_tightening*(proposed$vuln_coeff$Vs30-Omega_true$vuln_coeff$Vs30)+Omega_true$vuln_coeff$Vs30) * (means_mort-mean(means_mort)) + mean(means_mort), var_coeffs[1]*proposed$Lambda3$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='displacement')] <-  sampled_multivar[,2] #rlnorm(length(means_disp),  mean_coeffs[2]*(prior_tightening*(proposed$vuln_coeff$SHDI-Omega_true$vuln_coeff$SHDI)+Omega_true$vuln_coeff$SHDI) * (means_disp-mean(means_disp)) + mean(means_disp), var_coeffs[2]*proposed$Lambda2$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$impact=='buildDam')] <- sampled_multivar[,3] #rlnorm(length(means_bd), mean_coeffs[3]*(prior_tightening*(proposed$vuln_coeff$PDens-Omega_true$vuln_coeff$PDens)+Omega_true$vuln_coeff$PDens) * (means_bd-mean(means_bd)) + mean(means_bd), var_coeffs[3]*proposed$Lambda4$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
  }
  return(list(poly=impact_sample))
}

impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 100))
#plot_impact_sample(impact_sample, impact_filter='buildDam')

AlgoParams$smc_steps <- 105
AlgoParams$smc_Npart <- 2500
AlgoParams$m_CRPS <- 60
AlgoParams$smc_alpha <- 0.9
AlgoParams$rel_weightings <- c(0,1)
AlgoParams$kernel_sd$mortality <- 5
AlgoParams$kernel_sd$buildDam <- 0.1
AlgoParams$log_offset <- 0
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_Npart', AlgoParams$smc_Npart,'_M', AlgoParams$m_CRPS,'_simulated3D_withcorr_noplus10_scroremves_weightings5M1D0.1BD')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)

AlgoParams$kernel_sd$mortality <- 1
AlgoParams$kernel_sd$buildDam <- 1
AlgoParams$log_offset <- 10
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_Npart', AlgoParams$smc_Npart,'_M', AlgoParams$m_CRPS,'_simulated3D_withcorr_withplus10_scroremves')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)

AlgoParams$kernel_sd$mortality <- 1
AlgoParams$kernel_sd$buildDam <- 1
AlgoParams$log_offset <- 0
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_Npart', AlgoParams$smc_Npart,'_M', AlgoParams$m_CRPS,'_simulated3D_withcorr_noplus10_scroremves')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)


plot(AlgoResults$Omega_sample_phys[,1,1], AlgoResults$Omega_sample_phys[,2,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),1,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),2,s_finish], col='green')
points(1, 0.6, col='red', pch=19)

plot(AlgoResults$Omega_sample_phys[,3,1], AlgoResults$Omega_sample_phys[,4,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),3,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),4,s_finish], col='green')
points(1.2, 0.9, col='red', pch=19)

plot(AlgoResults$Omega_sample_phys[,5,1], AlgoResults$Omega_sample_phys[,6,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),5,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),6,s_finish], col='green')
points(-0.5, 0.2, col='red', pch=19)

#Some of the uniform distributions aren't going great in tails still e.g. 
plot(density(AlgoResults$Omega_sample_phys[,9,1], bw=0.03))
lines(density(AlgoResults$Omega_sample_phys[,9,86], bw=0.03), col='red')

par_i <- 18
plot(density(AlgoResults$Omega_sample_phys[,par_i,1], bw=0.03))
lines(density(AlgoResults$Omega_sample_phys[,par_i,86], bw=0.03), col='red')
lines(density(AlgoResults2$Omega_sample_phys[,par_i,63], bw=0.03), col='green')
lines(density(rLaplace(2000, 0,0.25)), col='blue')
plot(AlgoResults$propCOV[18,18,]);points(AlgoResults2$propCOV[18,18,], col='red')

plot(density(AlgoResults$Omega_sample_phys[,14,1], bw=0.03))
lines(density(AlgoResults$Omega_sample_phys[,14,90], bw=0.03), col='red')

plot(density(AlgoResults$Omega_sample_phys[,2,1], bw=0.03))
lines(density(AlgoResults$Omega_sample_phys[,2,70], bw=0.03), col='red')

propCOR <- array(NA, dim=c(24, 24, s_finish))
for (i in 2:s_finish){
  propCOR[,,i] <- cov2cor(AlgoResults$propCOV[,,i])
}
plot(propCOR[9,3,])

plot(AlgoResults$Omega_sample[,3,180])
points(AlgoResults$Omega_sample[,3,210], col='green')

plot(AlgoResults$Omega_sample[,9,200])
points(AlgoResults$Omega_sample[,9,210], col='green')

d1 <- runif(1000)
d2 <- runif(1000)
cor(d1,d2)



#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#----------------------------------COMPARE DIFFERENT NPART---------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------


AR2500 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-30_093349_alpha0.9_simulated3D_withcorr_noplus10_newparams_CRPS_M60_scroremves_Npart2500')
AR2500_smallcov <- AlgoResults
AR10000 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-29_120705_alpha0.9_simulated3D_withcorr_noplus10_M60_scroremves_Npart10000')

plot(AR2500$tolerancestore,ylim=c(0.5,2.1))
points(AR10000$tolerancestore, col='red')
#almost identical!

plot(AR2500$accrate_store)
points(AR10000$accrate_store, col='red')

#why is propcov[,,1] for 10000 particles so different? 
#pre_resample_steps <- which(AR2500$essstore == 2500) -1 
AR2500$s_finish <- which(AR2500$accrate_store<0.05)[1] #pre_resample_steps[which(pre_resample_steps > which(AR2500$accrate_store<0.05)[1])][1]

#pre_resample_steps <- which(AR10000$essstore == 10000) -1 
AR10000$s_finish <- which(AR10000$accrate_store<0.05)[1] #pre_resample_steps[which(pre_resample_steps > which(AR10000$accrate_store<0.05)[1])][1]
AR2500_smallcov$s_finish <- which(AR2500_smallcov$accrate_store<0.05)[1]

ggplot() + 
  geom_histogram(data=data.frame(v=AR2500$Omega_sample_phys[,15,AR2500$s_finish], w=AR2500$W[,AR2500$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR10000$Omega_sample_phys[,15,AR10000$s_finish], w=AR10000$W[,AR10000$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  geom_histogram(data=data.frame(v=AR2500_smallcov$Omega_sample_phys[,15,AR2500_smallcov$s_finish], w=AR2500_smallcov$W[,AR2500_smallcov$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='red', lwd=0.2, fill='red') #+
#geom_vline(xintercept=1, col='red')

ggplot() + 
  geom_histogram(data=data.frame(v=AR2500$Omega_sample_phys[,2,AR2500$s_finish], w=AR2500$W[,AR2500$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR10000$Omega_sample_phys[,2,AR10000$s_finish], w=AR10000$W[,AR10000$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  geom_vline(xintercept=0.6, col='red')

ggplot() + 
  geom_histogram(data=data.frame(v=AR2500$Omega_sample_phys[,6,AR2500$s_finish], w=AR2500$W[,AR2500$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR10000$Omega_sample_phys[,6,AR10000$s_finish], w=AR10000$W[,AR10000$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  geom_vline(xintercept=0.6, col='red')

ggplot() + 
  geom_histogram(data=data.frame(v=AR2500$Omega_sample_phys[,16,AR2500$s_finish], w=AR2500$W[,AR2500$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR10000$Omega_sample_phys[,16,AR10000$s_finish], w=AR10000$W[,AR10000$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  stat_function(fun = function(x) dLaplace(x, 0, 0.25), col='red')

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#-------------------------------COMPARE DIFFERENT NEVENTS ---------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

AR200Ev <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-30_093349_alpha0.9_simulated3D_withcorr_noplus10_newparams_CRPS_M60_scroremves_Npart2500')
AR600Ev <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-31_224922_alpha0.9_simulated3D_withcorr_noplus10_newparams_CRPS_M60_scroremves_Npart2500_600events')

plot(AR200Ev$tolerancestore,ylim=c(0.5,2.1))
points(AR600Ev$tolerancestore, col='red')
#almost identical!

plot(AR200Ev$accrate_store)
points(AR600Ev$accrate_store, col='red')

#why is propcov[,,1] for 10000 particles so different? 
#pre_resample_steps <- which(AR2500$essstore == 2500) -1 
AR200Ev$s_finish <- which(AR200Ev$accrate_store<0.05)[1] #pre_resample_steps[which(pre_resample_steps > which(AR2500$accrate_store<0.05)[1])][1]

#pre_resample_steps <- which(AR10000$essstore == 10000) -1 
AR600Ev$s_finish <- which(AR600Ev$accrate_store<0.05)[1] #pre_resample_steps[which(pre_resample_steps > which(AR10000$accrate_store<0.05)[1])][1]

ggplot() + 
  geom_histogram(data=data.frame(v=AR200Ev$Omega_sample_phys[,5,AR200Ev$s_finish], w=AR200Ev$W[,AR200Ev$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR600Ev$Omega_sample_phys[,5,AR600Ev$s_finish], w=AR600Ev$W[,AR600Ev$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow')

ggplot() + 
  geom_histogram(data=data.frame(v=AR200Ev$Omega_sample_phys[,16,AR200Ev$s_finish], w=AR200Ev$W[,AR200Ev$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR600Ev$Omega_sample_phys[,16,AR600Ev$s_finish], w=AR600Ev$W[,AR600Ev$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  stat_function(fun = function(x) dLaplace(x, 0, 0.25), col='red')


library(rstanarm)
glm_out <- log(data_y$observed[which(data_y$impact=='mortality')])-mean(means_mort)
glm_in <- means_mort-mean(means_mort)
fit_stan <- stan_glm(glm_out ~ glm_in - 1,
                     prior=NULL)
fit_samples <- as.matrix(fit_stan)

glm_out <- log(data_y$observed[which(data_y$impact=='mortality')[1:200*3]])-mean(means_mort)
glm_in <- means_mort[1:200*3]-mean(means_mort)
fit_stan200_ev <- stan_glm(glm_out ~ glm_in - 1,
                           prior=NULL)
fit_samples_200ev <- as.matrix(fit_stan200_ev)
hist(fit_samples_200ev[,1], freq=F)
lines(density(fit_samples[,1])) # a bit different but not thaaaat different

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#-------------------------------COMPARE log(x+10) vs log(x) -------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

ARplus10 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-01_235846_alpha0.9_Npart2500_M60_simulated3D_withcorr_withplus10_scroremves')
AR <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-02_040704_alpha0.9_Npart2500_M60_simulated3D_withcorr_noplus10_scroremves')

plot(ARplus10$tolerancestore,ylim=c(0.5,2.1))
points(AR$tolerancestore, col='red')
#almost identical!

plot(ARplus10$accrate_store)
points(AR$accrate_store, col='red')

#why is propcov[,,1] for 10000 particles so different? 
#pre_resample_steps <- which(AR2500$essstore == 2500) -1 
ARplus10$s_finish <- which(ARplus10$accrate_store<0.05)[1] #pre_resample_steps[which(pre_resample_steps > which(AR2500$accrate_store<0.05)[1])][1]

#pre_resample_steps <- which(AR10000$essstore == 10000) -1 
AR$s_finish <- which(AR$accrate_store<0.05)[1] #pre_resample_steps[which(pre_resample_steps > which(AR10000$accrate_store<0.05)[1])][1]

ggplot() + 
  geom_histogram(data=data.frame(v=ARplus10$Omega_sample_phys[,5,ARplus10$s_finish], w=ARplus10$W[,ARplus10$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR$Omega_sample_phys[,5,AR$s_finish], w=AR$W[,AR$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow')

ggplot() + 
  geom_histogram(data=data.frame(v=ARplus10$Omega_sample_phys[,16,ARplus10$s_finish], w=ARplus10$W[,ARplus10$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR$Omega_sample_phys[,16,AR$s_finish], w=AR$W[,AR$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  stat_function(fun = function(x) dLaplace(x, 0, 0.25), col='red')

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#-------------------------------COMPARE RELATIVE WEIGHTINGS -------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

AReven <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-02_040704_alpha0.9_Npart2500_M60_simulated3D_withcorr_noplus10_scroremves')
ARweighted <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-01_193532_alpha0.9_Npart2500_M60_simulated3D_withcorr_noplus10_scroremves_weightings5M1D0.1BD')

plot(AReven$tolerancestore,ylim=c(0.5,2.1))
points(ARweighted$tolerancestore, col='red')
#almost identical!

plot(AReven$accrate_store)
points(ARweighted$accrate_store, col='red')

#why is propcov[,,1] for 10000 particles so different? 
#pre_resample_steps <- which(AR2500$essstore == 2500) -1 
AReven$s_finish <- which(AReven$accrate_store<0.05)[1] #pre_resample_steps[which(pre_resample_steps > which(AR2500$accrate_store<0.05)[1])][1]

#pre_resample_steps <- which(AR10000$essstore == 10000) -1 
ARweighted$s_finish <- which(ARweighted$accrate_store<0.05)[1] #pre_resample_steps[which(pre_resample_steps > which(AR10000$accrate_store<0.05)[1])][1]

ggplot() + 
  geom_histogram(data=data.frame(v=AReven$Omega_sample_phys[,5,AReven$s_finish], w=AReven$W[,AReven$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=ARweighted$Omega_sample_phys[,5,ARweighted$s_finish], w=ARweighted$W[,ARweighted$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow')

ggplot() + 
  geom_histogram(data=data.frame(v=AReven$Omega_sample_phys[,15,AReven$s_finish], w=AReven$W[,AReven$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=ARweighted$Omega_sample_phys[,15,ARweighted$s_finish], w=ARweighted$W[,ARweighted$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  stat_function(fun = function(x) dLaplace(x, 0, 0.25), col='red')




# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- SIMPLE MODEL BUT HIGH DIMENSIONALITY ----------------------------------
# -------------------------------------------------------------------------------------------------------------
Omega <- Omega_true <- list(Lambda1 = list(nu=9, kappa=1.4),
                            Lambda2 = list(nu=10.5, kappa=1.7), #list(nu=10.65, kappa=1.5), #
                            Lambda3 = list(nu=9, kappa=1.2),
                            Lambda4 = list(nu=9.9, kappa=1.6),
                            theta= list(theta1=0.6),
                            eps=list(local=0.6, hazard_mort=0.3, hazard_disp=0.5, hazard_bd=0.4, hazard_cor=0.55),
                            #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                            vuln_coeff = list(PDens=0.05, SHDI=-0.1, GNIc=-0.05, Vs30=0.1, EQFreq=-0.15, FirstHaz=0.05, Night=0.05, FirstHaz.Night=0.1),
                            check = list(check=0.5))

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
  eps=list(local=list(dist='unif', min=0, max=1.5),
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

Model$HighLevelPriors <- function(Omega, Model, modifier=NULL){
  return(0)
}


n_events <- 200
set.seed(1)

cov_mort_disp = Omega_true$eps$hazard_cor * Omega_true$eps$hazard_mort * Omega_true$eps$hazard_disp
cov_mort_bd = Omega_true$eps$hazard_cor * Omega_true$eps$hazard_mort * Omega_true$eps$hazard_bd
cov_disp_bd = Omega_true$eps$hazard_cor * Omega_true$eps$hazard_disp * Omega_true$eps$hazard_bd
covar_matrix = cbind(c((Omega_true$eps$hazard_mort)^2, cov_mort_disp, cov_mort_bd), c(0, (Omega_true$eps$hazard_disp)^2, cov_disp_bd), c(0, 0, (Omega_true$eps$hazard_bd)^2))
covar_matrix[upper.tri(covar_matrix)] = covar_matrix[lower.tri(covar_matrix)]

eps_ev <- rmvnorm(n_events, rep(0,3), covar_matrix) * 3
vuln <- rmvnorm(n_events, rep(0, 8), diag(100, nrow=8))
vuln_adj <- vuln %*% unlist(Omega_true$vuln_coeff)

mort_impacts <- data.frame(polygon=1, impact='mortality', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = 1:n_events,
                           observed = rlnorm(n_events, Omega_true$Lambda1$nu + vuln_adj + eps_ev[,1], Omega_true$Lambda1$kappa))

disp_impacts <- data.frame(polygon=1, impact='displacement', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = 1:n_events,
                           observed = rlnorm(n_events, Omega_true$Lambda2$nu + vuln_adj + eps_ev[,2], Omega_true$Lambda2$kappa))

bd_impacts <- data.frame(polygon=1, impact='buildDam', sampled=NA,
                         iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                         build_type=NA, inferred=F, event_id = 1:n_events,
                         observed = rlnorm(n_events, Omega_true$Lambda3$nu + vuln_adj + eps_ev[,3], Omega_true$Lambda3$kappa))

data_y <- rbind(mort_impacts, disp_impacts, bd_impacts)
plot(data_y$observed[which(data_y$impact=='buildDam')])

#prior_tightening <- 0.1
SampleImpact <- function(dir, Model, proposed, AlgoParams){
  impact_sample <- list()
  for (i in 1:(AlgoParams$Np*AlgoParams$m_CRPS)){
    impact_sample[[i]] <- data_y
    #impact_sample[[i]]$sampled <-  abs(rnorm(length(means), 100*proposed$vuln_coeff$Vs30 * (means-mean(means))+1500, proposed$Lambda3$kappa*500)) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    cov_mort_disp = proposed$eps$hazard_cor * proposed$eps$hazard_mort * proposed$eps$hazard_disp
    cov_mort_bd = proposed$eps$hazard_cor * proposed$eps$hazard_mort * proposed$eps$hazard_bd
    cov_disp_bd = proposed$eps$hazard_cor * proposed$eps$hazard_disp * proposed$eps$hazard_bd
    covar_matrix = cbind(c((proposed$eps$hazard_mort)^2, cov_mort_disp, cov_mort_bd), c(0, (proposed$eps$hazard_disp)^2, cov_disp_bd), c(0, 0, (proposed$eps$hazard_bd)^2))
    covar_matrix[upper.tri(covar_matrix)] = covar_matrix[lower.tri(covar_matrix)]
    
    eps_ev <- rmvnorm(n_events, rep(0,3), covar_matrix)
    
    vuln_adj <- vuln %*% unlist(proposed$vuln_coeff)

    impact_sample[[i]]$sampled[1:n_events] <- rlnorm(n_events, proposed$Lambda1$nu + vuln_adj, proposed$Lambda1$kappa) + eps_ev[,1] #rlnorm(length(means_mort),  mean_coeffs[1]*(prior_tightening*(proposed$vuln_coeff$Vs30-Omega_true$vuln_coeff$Vs30)+Omega_true$vuln_coeff$Vs30) * (means_mort-mean(means_mort)) + mean(means_mort), var_coeffs[1]*proposed$Lambda3$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[(n_events+1): (2*n_events)] <-  rlnorm(n_events, proposed$Lambda2$nu  + vuln_adj, proposed$Lambda2$kappa) + eps_ev[,2] #rlnorm(length(means_disp),  mean_coeffs[2]*(prior_tightening*(proposed$vuln_coeff$SHDI-Omega_true$vuln_coeff$SHDI)+Omega_true$vuln_coeff$SHDI) * (means_disp-mean(means_disp)) + mean(means_disp), var_coeffs[2]*proposed$Lambda2$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[(2*n_events+1):(3*n_events)] <- rlnorm(n_events, proposed$Lambda3$nu + vuln_adj, proposed$Lambda3$kappa) + eps_ev[,3] #rlnorm(length(means_bd), mean_coeffs[3]*(prior_tightening*(proposed$vuln_coeff$PDens-Omega_true$vuln_coeff$PDens)+Omega_true$vuln_coeff$PDens) * (means_bd-mean(means_bd)) + mean(means_bd), var_coeffs[3]*proposed$Lambda4$kappa) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    impact_sample[[i]]$sampled[which(impact_sample[[i]]$sampled < 0)] <- 0
  }
  return(list(poly=impact_sample))
}

impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 60))
CalcDist(impact_sample, AlgoParams)
#plot_impact_sample(impact_sample, impact_filter='buildDam')

#impact_sample <- SampleImpact(dir, Model, AlgoResults$Omega_sample_phys[1,,70] %>% relist(skeleton=Model$skeleton), AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 60))
#plot_impact_sample(impact_sample, impact_filter='buildDam')
#CalcDist(impact_sample, AlgoParams)


AlgoParams$smc_steps <- 250
AlgoParams$smc_Npart <- 1000
AlgoParams$m_CRPS <- 60
AlgoParams$smc_alpha <- 0.9
AlgoParams$rel_weightings <- c(0,1)
AlgoParams$kernel_sd$mortality <- 5
AlgoParams$kernel_sd$buildDam <- 0.1
AlgoParams$log_offset <- 10
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_Npart', AlgoParams$smc_Npart,'_M', AlgoParams$m_CRPS,'_', 'SimpleHighDimensionalTestingLargerPropCOV')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes, oldtag='')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-30_174026_alpha0.9_Npart1000_M60_SimpleHighDimensionalTestingLargerPropCOV')

plot_correlated_posteriors(AlgoResults, Omega=Omega)
plot_correlated_posteriors(AlgoResults, Omega=Omega, pairings=rbind(c(13,14), c(15,16), c(17,18), c(19,20), c(21,22)))

plot_corr_posterior_vs_d(AlgoResults, Omega=Omega, pairing=c(19,20))

plot_corr_transf_posterior_vs_d(AlgoResults, Omega=Omega, pairing=c(5,6))

plot(AlgoResults$Omega_sample_phys[,2,30], AlgoResults$Omega_sample_phys[,3,30])

plot_acc_prob(AlgoResults); abline(h=0.05)
plot_d_vs_step(AlgoResults, 30); abline(h=5.1, col='red')
plot(AlgoResults$essstore)


# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- COMPARE DIFFERENT PROPCOV MULT ----------------------------------------
# -------------------------------------------------------------------------------------------------------------

AlgoResults_smallMult <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-29_220522_alpha0.9_Npart1000_M60_SimpleHighDimensionalTesting')

AlgoResults_largeMult <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-04-30_174026_alpha0.9_Npart1000_M60_SimpleHighDimensionalTestingLargerPropCOV')

par(mfrow=c(2,1))
plot_d_vs_step(AlgoResults_smallMult, 15); abline(h=5.1, col='red')
plot_d_vs_step(AlgoResults_largeMult, 15); abline(h=5.1, col='red')

plot_acc_prob(AlgoResults_smallMult); abline(h=0.05)
plot_acc_prob(AlgoResults_largeMult); abline(h=0.05)
par(mfrow=c(1,1))

n_unique_large <- c()
n_unique_small <- c()
for (i in 1:200){
  n_unique_large <- c(n_unique_large, length(unique(AlgoResults_largeMult$Omega_sample_phys[,1,i])))
  n_unique_small <- c(n_unique_small, length(unique(AlgoResults_smallMult$Omega_sample_phys[,1,i])))
}
plot(n_unique_small)
points(n_unique_large, col='blue')

#smaller seems to do better - larger doesn't seem to be any faster but less particles from pretty much the outset

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------


#Still drifting more from priors than would be desirable
#OPTIONS:
#  - PropCOV is too small / large to properly explore tails
#  - Not enough particles to properly explore tails
#  - Not enough particles so picking up wrong correlation between parameters, unsure how this then translates
#  - As this post becomes more clustered, this is what you'd expect
#  - Undesired correlation could be getting produced by duplicated particles?

plot(AlgoResults$Omega_sample_phys[,24,1])
points(AlgoResults$Omega_sample_phys[,24,63], col='red')

plot(AlgoResults$Omega_sample_phys[,7,1], freq=F)
points(AlgoResults$Omega_sample_phys[,7,90], col='red')

#--------------------------------------------
#------------TEST IN STAN--------------------
#--------------------------------------------

library(rstan)

fit1 <- stan(
  file = "RCode/simulated_test_stan.stan",  # Stan program
  data = list(N_events = n_events,
              observed_mort=data_y$observed[which(data_y$impact=='mortality')],
              observed_disp=data_y$observed[which(data_y$impact=='displacement')],
              observed_bd=data_y$observed[which(data_y$impact=='buildDam')],
              means_mort=means_mort,
              means_disp=means_disp,
              means_bd=means_bd
              ),    # named list of data
  chains = 3,             # number of Markov chains
  warmup = 4000,          # number of warmup iterations per chain
  iter = 6000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
)

plot(fit1, pars=c('mu_mort', 'sigma_mort', 'mu_disp', 'sigma_disp', 'mu_bd', 'sigma_bd'))
traceplot(fit1, pars=c('mu_mort', 'sigma_mort', 'mu_disp', 'sigma_disp', 'mu_bd', 'sigma_bd'))

summary(fit1)

fit1[['mu_disp']
#-------------------------------------------
#-------------------------------------------
#-------------------------------------------




AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-06_120756_alpha0.9_simulated3D_0corr_noplus10_midvar_reparam_CRPS_M20')

par_i <- 4
s_finish <- 75
plot(AlgoResults$Omega_sample_phys[, par_i, s_finish], AlgoResults$d[,1,s_finish])

library(rstanarm)
glm_out <- log(data_y$observed[which(data_y$impact=='mortality')])-mean(means_mort)
glm_in <- means_mort-mean(means_mort)
fit_stan <- stan_glm(glm_out ~ glm_in - 1,
                     prior=NULL)

fit_samples <- as.matrix(fit_stan)
hist(fit_samples[,1], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,1,s_finish]))
hist(fit_samples[,2], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,2,s_finish]))

glm_out <- log(data_y$observed[which(data_y$impact=='displacement')])-mean(means_disp)
glm_in <- means_disp-mean(means_disp)
fit_stan <- stan_glm(glm_out ~ glm_in - 1,
                     prior=NULL)

fit_samples <- as.matrix(fit_stan)
hist(fit_samples[,1], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,3,s_finish]))
hist(fit_samples[,2], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,4,s_finish]))

abline(v=0.6, col='red')

fit1 <- stan(
  file = "simulated_test_stan.stan",  # Stan program
  data = schools_data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = 0             # no progress shown
)

# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# --------------- SPEEDY MODEL TO CHECK WHY NUISANCE PARAMETERS ARE NOT STAYING AT PRIOR  ---------------------
# -------------------------------------------------------------------------------------------------------------

Model$HighLevelPriors <- function(Omega, Model, modifier=NULL){
  return(0)
}
Omega_true <- list(Lambda1 = list(nu=1, kappa=0.6),
                   Lambda2 = list(nu=1.2, kappa=0.9),
                   Lambda3 = list(nu=-0.5, kappa=0.2),
                   Lambda4 = list(nu=9.9, kappa=1.6),
                   theta= list(theta1=0.6),
                   eps = list(local=2.2292053, hazard_mort=0.8383464, hazard_disp=0.9, hazard_bd=0.9, hazard_cor=0.55),
                   vuln_coeff = list(PDens=0.05, SHDI=-0.5, GNIc=-0.1, Vs30=0.1, EQFreq=-0.1, FirstHaz=0.05, Night=0.05, FirstHaz.Night=0.1),
                   check = list(check=0.5))

Model$Priors <- list( #All uniform so currently not included in the acceptance probability. 
  Lambda1=list(mu=list(dist='laplace', location=1, scale=0.1), 
               sigma=list(dist='unif', min=0, max=1)), 
  Lambda2=list(mu=list(dist='laplace', location=1.2, scale=0.1), 
               sigma=list(dist='unif', min=0, max=1.5)),
  Lambda3=list(mu=list(dist='laplace', location=-0.5, scale=0.1), 
               sigma=list(dist='unif', min=0, max=1)),
  Lambda4=list(mu=list(dist='unif', min=8, max=12.5), 
               sigma=list(dist='unif', min=0.25, max=2.5)),
  theta=list(theta1=list(dist='unif', min=0, max=1)),
  eps=list(local=list(dist='unif', min=0, max=1.5),
           hazard_mort=list(dist='unif', min=0, max=1.5),
           hazard_disp=list(dist='unif', min=0, max=1.5),
           hazard_bd=list(dist='unif', min=0, max=1.5),
           hazard_cor=list(dist='unif', min=0, max=1)),
  vuln_coeff=list(PDens=list(dist='laplace', location=0, scale=0.25),
                  EQFreq=list(dist='laplace',  location=0, scale=0.25),
                  SHDI=list(dist='laplace', location=0, scale=0.25),
                  GNIc=list(dist='laplace', location=0, scale=0.25),
                  Vs30=list(dist='laplace', location=0, scale=0.25),
                  Mag=list(dist='laplace', location=0, scale=0.25),
                  FirstHaz=list(dist='laplace', location=0, scale=0.25),
                  Night=list(dist='laplace', location=0, scale=0.25),
                  FirstHaz.Night=list(dist='laplace', location=0, scale=0.25)),
  check=list(check=list(dist='unif', min=0, max=1))
)

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

n_events <- 100
set.seed(1)

data_y <- data.frame(polygon=j, impact='mortality', sampled=NA,
                           iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                           build_type=NA, inferred=F, event_id = i,
                           observed = rlnorm(n_events, Omega_true$Lambda1$nu, Omega_true$Lambda1$kappa))

#plot(mort_impacts$observed)

#prior_tightening <- 0.1
SampleImpact <- function(dir, Model, proposed, AlgoParams){
  impact_sample <- list()
  for (i in 1:(AlgoParams$Np*AlgoParams$m_CRPS)){
    impact_sample[[i]] <- data_y
    #impact_sample[[i]]$sampled <-  abs(rnorm(length(means), 100*proposed$vuln_coeff$Vs30 * (means-mean(means))+1500, proposed$Lambda3$kappa*500)) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    
    impact_sample[[i]]$sampled <-  rlnorm(n_events, proposed$Lambda1$nu, proposed$Lambda1$kappa)
  }
  return(list(poly=impact_sample))
}

impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 100))
#plot_impact_sample(impact_sample, impact_filter='mortality')

AlgoParams$smc_steps <- 50
AlgoParams$smc_Npart <- 10000
AlgoParams$m_CRPS <- 60
AlgoParams$smc_alpha <- 0.9
AlgoParams$rel_weightings <- c(0,1)
AlgoParams$kernel_sd$mortality <- 1
AlgoParams$kernel_sd$buildDam <- 1
AlgoParams$kernel_sd$buildDest <- 1
tag_notes <- paste0('Npart', AlgoParams$smc_Npart,'_alpha', AlgoParams$smc_alpha, '_Np', AlgoParams$Np, '_M',AlgoParams$m_CRPS, '_AllParamsSpeedy_PropCOVmult0.1')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = F,tag_notes=tag_notes)

AlgoResults1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_110124_alpha0.9_simulatedspeed_M60_Np500_alpha0.9_PropCOVmult1')
AlgoResults2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_111427_alpha0.9_simulatedspeed_M60_Np500_alpha0.9_PropCOVmult0.1')
AlgoResults3 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_113832_alpha0.9_simulatedspeed_M60_Np500_alpha0.9_PropCOVmult0.01')

plot(density(AlgoResults$Omega_sample[,6,1], bw=0.05))
lines(density(AlgoResults$Omega_sample[,6,50], bw=0.05), col='blue')

plot(density(AlgoResults1$Omega_sample[,1,1], bw=0.01))
lines(density(AlgoResults1$Omega_sample[,1,20], bw=0.01), col='green')
lines(density(AlgoResults2$Omega_sample[,1,20], bw=0.01), col='blue')
lines(density(AlgoResults3$Omega_sample[,1,20], bw=0.01), col='red')
lines(density(AlgoResults2$Omega_sample[,1,50], bw=0.01), col='cyan')

plot(density(AlgoResults1$Omega_sample[,2,1], bw=0.1), xlim=c(-2,2))
lines(density(AlgoResults1$Omega_sample[,2,20], bw=0.1), col='green')
lines(density(AlgoResults2$Omega_sample[,2,20], bw=0.1), col='blue')
lines(density(AlgoResults3$Omega_sample[,2,20], bw=0.01), col='red')
lines(density(AlgoResults2$Omega_sample[,2,50], bw=0.1), col='cyan')

plot(density(AlgoResults1$Omega_sample[,10,1], bw=0.1))
lines(density(AlgoResults1$Omega_sample[,10,20], bw=0.1), col='green')
lines(density(AlgoResults2$Omega_sample[,10,40], bw=0.1), col='blue')
lines(density(AlgoResults3$Omega_sample[,10,50], bw=0.1), col='red')

plot(density(AlgoResults1$Omega_sample[,16,1], bw=0.01))
lines(density(AlgoResults1$Omega_sample[,16,20], bw=0.01), col='green')
lines(density(AlgoResults2$Omega_sample[,16,40], bw=0.01), col='blue')
lines(density(AlgoResults2$Omega_sample[,16,50], bw=0.01), col='cyan')
lines(density(AlgoResults3$Omega_sample[,16,50], bw=0.01), col='red')


hist(AlgoResults1$Omega_sample_phys[,10,1])

library(rstanarm)
fit_stan <- stan_glm(log(data_y$observed) ~ 1,
                     prior=NULL)

fit_samples <- as.matrix(fit_stan)
hist(fit_samples[,1], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,1,30]))

hist(fit_samples[,2], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,2,30]))

par_i <- 1
plot(density(AlgoResults$Omega_sample_phys[,par_i,1], bw=0.02))
plot(density(rLaplace(4312, -0.5, 0.1), bw=0.02))
lines(density(AlgoResults$Omega_sample_phys[,par_i,30], bw=0.02), col='red')
lines(density(AlgoResults2$Omega_sample_phys[,par_i,45], bw=0.02), col='green')
lines(density(AlgoResults3$Omega_sample_phys[,par_i,45], bw=0.02), col='blue')

plot(AlgoResults$propCOV[17,17,])
points(AlgoResults2$propCOV[17,17,], col='green')



n_unique <- c()
for (i in 1:50){
  n_unique <- c(n_unique, length(unique(AlgoResults$Omega_sample_phys[,16,i])))
}
plot(n_unique)

#even with Np = 2000 and alpha = 0.95 the tails of those that aren't uniform are too small. 
#uniform parameters seem to be doing ok
#don't think it's an issue with the number of particles as doesn't improve with larger Np. 
#don't think it's an issue with the prior for nu being tighter than kappa, as a wider prior on nu doesn't help
#don't think it's an issue with propCOV too large, as allowing it to go lower doesn't help
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_start_step_2024-03-20_124532_alpha0.95_simulatedspeed_M60_Np500_alpha0.95')
plot(density(AlgoResults$Omega_sample_phys[,17,1], bw=0.02))
lines(density(AlgoResults$Omega_sample_phys[,17,40], bw=0.02), col='red')

hist(AlgoResults$Omega_sample_phys[,1,25], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,1,50]), col='blue')


plot(AlgoResults$Omega_sample_phys[,1,1], AlgoResults$Omega_sample_phys[,2,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),1,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),2,s_finish], col='green')
points(1, 0.6, col='red', pch=19)

plot(AlgoResults$Omega_sample_phys[,3,1], AlgoResults$Omega_sample_phys[,4,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),3,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),4,s_finish], col='green')
points(1.2, 0.9, col='red', pch=19)

plot(AlgoResults$Omega_sample_phys[,5,1], AlgoResults$Omega_sample_phys[,6,1], xlab='Mean Parameter', ylab='Variance Parameter')
points(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),5,s_finish], AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s_finish]>0),6,s_finish], col='green')
points(-0.5, 0.2, col='red', pch=19)

hist(AlgoResults$Omega_sample_phys[,7,1], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,7,90]), col='red')

plot(AlgoResults$Omega_sample_phys[,7,1], freq=F)
points(AlgoResults$Omega_sample_phys[,7,90], col='red')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-06_120756_alpha0.9_simulated3D_0corr_noplus10_midvar_reparam_CRPS_M20')

par_i <- 4
s_finish <- 75
plot(AlgoResults$Omega_sample_phys[, par_i, s_finish], AlgoResults$d[,1,s_finish])

library(rstanarm)
glm_out <- log(data_y$observed[which(data_y$impact=='mortality')])-mean(means_mort)
glm_in <- means_mort-mean(means_mort)
fit_stan <- stan_glm(glm_out ~ glm_in - 1,
                     prior=NULL)

fit_samples <- as.matrix(fit_stan)
hist(fit_samples[,1], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,1,s_finish]))
hist(fit_samples[,2], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,2,s_finish]))

glm_out <- log(data_y$observed[which(data_y$impact=='displacement')])-mean(means_disp)
glm_in <- means_disp-mean(means_disp)
fit_stan <- stan_glm(glm_out ~ glm_in - 1,
                     prior=NULL)

fit_samples <- as.matrix(fit_stan)
hist(fit_samples[,1], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,3,s_finish]))
hist(fit_samples[,2], freq=F)
lines(density(AlgoResults$Omega_sample_phys[,4,s_finish]))



abline(v=0.6, col='red')


#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------- REDUCED DIMENSIONALITY -----------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

Omega_true <- list(Lambda1 = list(nu=1, kappa=0.6),
                   Lambda2 = list(nu=1.2, kappa=0.9),
                   Lambda3 = list(nu=-0.5, kappa=0.2),
                   Lambda4 = list(nu=9.9, kappa=1.6),
                   theta= list(theta1=0.6),
                   eps = list(local=2.2292053, hazard_mort=0.8383464, hazard_disp=0.9, hazard_bd=0.9, hazard_cor=0.55),
                   vuln_coeff = list(PDens=0.05, SHDI=-0.5, GNIc=-0.1, Vs30=0.1, EQFreq=-0.1, FirstHaz=0.05, Night=0.05, FirstHaz.Night=0.1),
                   check = list(check=0.5))

Model$skeleton <- list(
  Lambda1=list(nu=NA,kappa=NA#,alpha=NA
  ), 
  Lambda2=list(nu=NA,kappa=NA)
)


Model$HighLevelPriors <- function(Omega, Model, modifier=NULL){
  return(0)
}
Omega_true <- list(Lambda1 = list(nu=1, kappa=0.6),
                   Lambda2 = list(nu=1.2, kappa=0.9))

Model$Priors <- list( #All uniform so currently not included in the acceptance probability. 
  Lambda1=list(mu=list(dist='laplace', location=1, scale=0.1), 
               sigma=list(dist='unif', min=0, max=1)), 
  Lambda2=list(mu=list(dist='laplace', location=1.2, scale=0.1), 
               sigma=list(dist='unif', min=0, max=1.5))
)

#Set lower and upper bounds for the parameters
Model$par_lb <- c()
Model$par_ub <- c()

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

addTransfParams <- function(Omega, I0=Model$I0){
  return(Omega)
}

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

n_events <- 100
set.seed(1)

data_y <- data.frame(polygon=j, impact='mortality', sampled=NA,
                     iso3='ABC', sdate=as.Date('16-03-1999'), qualifier=NA,
                     build_type=NA, inferred=F, event_id = i,
                     observed = rlnorm(n_events, Omega_true$Lambda1$nu, Omega_true$Lambda1$kappa))

#plot(mort_impacts$observed)

#prior_tightening <- 0.1
SampleImpact <- function(dir, Model, proposed, AlgoParams){
  impact_sample <- list()
  for (i in 1:(AlgoParams$Np*AlgoParams$m_CRPS)){
    impact_sample[[i]] <- data_y
    #impact_sample[[i]]$sampled <-  abs(rnorm(length(means), 100*proposed$vuln_coeff$Vs30 * (means-mean(means))+1500, proposed$Lambda3$kappa*500)) #rt(NROW(data_y), proposed$Lambda1$nu, 3) * proposed$Lambda1$kappa #round(rt(NROW(data_y), proposed$Lambda1$nu, 5) * proposed$Lambda1$kappa)
    
    impact_sample[[i]]$sampled <-  rlnorm(n_events, proposed$Lambda1$nu, proposed$Lambda1$kappa)
  }
  return(list(poly=impact_sample))
}

impact_sample <- SampleImpact(dir, Model, Omega_true, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 100))
#plot_impact_sample(impact_sample, impact_filter='mortality')

AlgoParams$smc_steps <- 100
AlgoParams$smc_Npart <- 10000
AlgoParams$m_CRPS <- 60
AlgoParams$smc_alpha <- 0.9
AlgoParams$rel_weightings <- c(0,1)
AlgoParams$kernel_sd$mortality <- 1
AlgoParams$kernel_sd$buildDam <- 1
AlgoParams$kernel_sd$buildDest <- 1
tag_notes <- paste0('alpha', AlgoParams$smc_alpha, '_Np', AlgoParams$smc_Npart, '_M',AlgoParams$m_CRPS, '_4params_PropCOVmult0.1')
AlgoResults <- delmoral_parallel(AlgoParams, Model, unfinished = T,tag_notes=tag_notes, oldtag='abcsmc_start_step_2024-03-26_160110_alpha0.9_4params_M60_Np1000_alpha0.9_PropCOVmult0.1')

s <- 50
hist(AlgoResults$Omega_sample_phys[,3,s], breaks=50, freq=F)
dens_xgr <- seq(min(AlgoResults$Omega_sample_phys[,3,s]), max(AlgoResults$Omega_sample_phys[,3,s]), 0.01)
lines(dens_xgr, dLaplace(dens_xgr, 1.2, 0.1))

100-AlgoResults$tolerancestore[2:50]/AlgoResults$tolerancestore[1:49]*100

ggplot() + 
  geom_histogram(data=data.frame(v=AlgoResults$Omega_sample_phys[,3,s], w=AlgoResults$W[,s]), aes(x=v,y=..density.., weight = w)) + 
  stat_function(fun = function(x) dLaplace(x, 1.2, 0.1))

ggplot() + 
  geom_histogram(data=data.frame(v=AlgoResults$Omega_sample_phys[,4,s], w=AlgoResults$W[,s]), aes(x=v,y=..density.., weight = w)) + 
  scale_x_continuous(limits = c(0, 1.5)) +
  stat_function(fun = function(x) dunif(x, 0, 1.5))

#-----------------------------------------------------------------------------------------
#-------------------------- Compare different Np and alpha: ------------------------------
#-----------------------------------------------------------------------------------------

AR_Np1000_alpha0.9 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_190945_alpha0.9_4params_M60_Np1000_alpha0.9_PropCOVmult0.1')
AR_Np1000_alpha0.99 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_191925_alpha0.99_4params_M60_Np1000_alpha0.99_PropCOVmult0.1')
AR_Np10000_alpha0.9 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-28_122520_alpha0.9_Np10000_M60_4params_PropCOVmult0.1')

calc_s_finish <- function(AlgoResults){
  tol_percentage_decrease <- 100-AlgoResults$tolerancestore[2:length(AlgoResults$tolerancestore)]/AlgoResults$tolerancestore[1:(length(AlgoResults$tolerancestore)-1)]*100
  AlgoResults$s_finish <- which(tol_percentage_decrease<0.1)[1]
  return(AlgoResults)
}
AR_Np1000_alpha0.9 %<>% calc_s_finish()
AR_Np1000_alpha0.99 %<>% calc_s_finish()
AR_Np10000_alpha0.9 %<>% calc_s_finish()

ggplot() + 
  geom_histogram(data=data.frame(v=AR_Np1000_alpha0.9$Omega_sample_phys[,1,AR_Np1000_alpha0.9$s_finish], w=AR_Np1000_alpha0.9$W[,AR_Np1000_alpha0.9$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,1,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR_Np10000_alpha0.9$Omega_sample_phys[,1,AR_Np10000_alpha0.9$s_finish], w=AR_Np10000_alpha0.9$W[,AR_Np10000_alpha0.9$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  geom_vline(xintercept=1, col='red')

ggplot() + 
  geom_histogram(data=data.frame(v=AR_Np1000_alpha0.9$Omega_sample_phys[,2,AR_Np1000_alpha0.9$s_finish], w=AR_Np1000_alpha0.9$W[,AR_Np1000_alpha0.9$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2, fill='blue') +
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,2,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') +
  geom_histogram(data=data.frame(v=AR_Np10000_alpha0.9$Omega_sample_phys[,2,AR_Np10000_alpha0.9$s_finish], w=AR_Np10000_alpha0.9$W[,AR_Np10000_alpha0.9$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='yellow', lwd=0.2, fill='yellow') +
  geom_vline(xintercept=0.6, col='red')

ggplot() + 
  #geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,3,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='green', lwd=0.2, fill='green') + 
  geom_histogram(data=data.frame(v=AR_Np1000_alpha0.9$Omega_sample_phys[,3,AR_Np1000_alpha0.9$s_finish], w=AR_Np1000_alpha0.9$W[,AR_Np1000_alpha0.9$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.05, col='blue', lwd=0.5, fill='blue') + 
  geom_histogram(data=data.frame(v=AR_Np10000_alpha0.9$Omega_sample_phys[,3,AR_Np10000_alpha0.9$s_finish], w=AR_Np10000_alpha0.9$W[,AR_Np10000_alpha0.9$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.05, col='yellow', lwd=0.5, fill='yellow') +
  xlim(qLaplace(0.001, 1.2,0.1), qLaplace(0.999, 1.2,0.1)) +
  stat_function(fun = function(x) dLaplace(x, 1.2, 0.1), col='red')

ggplot() + 
  geom_histogram(data=data.frame(v=AR_Np1000_alpha0.99$Omega_sample_phys[,4,AR_Np1000_alpha0.99$s_finish], w=AR_Np1000_alpha0.99$W[,AR_Np1000_alpha0.99$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.1, col='green', lwd=0.2, fill='green') + 
  geom_histogram(data=data.frame(v=AR_Np1000_alpha0.9$Omega_sample_phys[,4,AR_Np1000_alpha0.9$s_finish], w=AR_Np1000_alpha0.9$W[,AR_Np1000_alpha0.9$s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.1, col='blue', lwd=0.2, fill='blue') + 
  xlim(0,1.5) +
  stat_function(fun = function(x) dunif(x, 0, 1.5), col='red')

ggplot() + 
  geom_histogram(data=data.frame(v=AR_Np1000_alpha0.9$Omega_sample_phys[,4,s], w=AR_Np1000_alpha0.9$W[,s]), aes(x=v,y=..density.., weight = w)) + 
  scale_x_continuous(limits = c(0, 1.5)) +
  stat_function(fun = function(x) dunif(x, 0, 1.5))

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-------------------------- Compare different stopping rules: ----------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

AR <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_start_step_2024-03-28_122520_alpha0.9_Np10000_M60_4params_PropCOVmult0.1')
#AR <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_190945_alpha0.9_4params_M60_Np1000_alpha0.9_PropCOVmult0.1')
#AR <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_191925_alpha0.99_4params_M60_Np1000_alpha0.99_PropCOVmult0.1')

s_finish1 <- 45
s_finish2 <- 75
s_finish3 <- 50

ggplot() + 
  geom_histogram(data=data.frame(v=AR$Omega_sample_phys[,1,s_finish1], w=AR$W[,s_finish1]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2) +
  geom_histogram(data=data.frame(v=AR$Omega_sample_phys[,1,s_finish2], w=AR$W[,s_finish2]), aes(x=v,y=..density.., weight = w), alpha=0.3, fill='green', col='green', lwd=0.2) +
  #geom_histogram(data=data.frame(v=AR_Np10000_alpha0.9$Omega_sample_phys[,1,s_finish3], w=AR_Np10000_alpha0.9$W[,s_finish3]), aes(x=v,y=..density.., weight = w), alpha=0.3, fill='red', col='red', lwd=0.2) +
  geom_vline(xintercept=1, col='black')

ggplot() + 
  #geom_histogram(data=data.frame(v=AR_Np10000_alpha0.9$Omega_sample_phys[,2,s_finish1], w=AR_Np10000_alpha0.9$W[,s_finish1]), aes(x=v,y=..density.., weight = w), alpha=0.3, col='blue', lwd=0.2) +
  geom_histogram(data=data.frame(v=AR$Omega_sample_phys[,2,s_finish2], w=AR$W[,s_finish2]), aes(x=v,y=..density.., weight = w), alpha=0.3, fill='green', col='green', lwd=0.2) +
  geom_histogram(data=data.frame(v=AR$Omega_sample_phys[,2,s_finish3], w=AR$W[,s_finish3]), aes(x=v,y=..density.., weight = w), alpha=0.3, fill='red', col='red', lwd=0.2) +
  geom_vline(xintercept=0.6, col='black')


plot_dens_s <- function(s, col='red'){
  ps <- ggplot() + 
          geom_histogram(data=data.frame(v=AR$Omega_sample_phys[,3,s], w=AR$W[,s]), aes(x=v,y=..density.., weight = w), alpha=0.3, fill=col, col=col, lwd=0.2) +
          xlim(qLaplace(0.001, 1.2,0.1), qLaplace(0.999, 1.2,0.1)) +
          stat_function(fun = function(x) dLaplace(x, 1.2, 0.1), col='black')
  return(ps)
}

s_breaks <- c(30,45,60,75)
p_s1 <- plot_dens_s(s_breaks[1], col='yellow')
p_s2 <- plot_dens_s(s_breaks[2], col='green')
p_s3 <- plot_dens_s(s_breaks[3], col='red')
p_s4 <- plot_dens_s(s_breaks[4], col='blue')

s_finish <- 75
unique_part <- c()
for (i in 1:s_finish){
  unique_part <- c(unique_part, length(unique(AR$Omega_sample[AR$W[,i]>0,1,i])))
}
  
plot_df <- data.frame(s = 1:s_finish,
                      acc_rate = AR$accrate_store[1:s_finish], 
                      unique_alive_part = unique_part, 
                      tolerance_store_perc_change = 1:s_finish )

accrate_plot <- ggplot(plot_df, aes(x=s, y=acc_rate)) + geom_point() + geom_vline(xintercept=s_breaks[1], col='yellow') + 
  geom_vline(xintercept=s_breaks[2], col='green') + geom_vline(xintercept=s_breaks[3], col='red') + geom_vline(xintercept=s_breaks[4], col='blue') +
  geom_hline(yintercept=0.05, col='red')

unique_part_plot <- ggplot(plot_df, aes(x=s, y=unique_alive_part)) + geom_point() + geom_vline(xintercept=s_breaks[1], col='yellow') + 
  geom_vline(xintercept=s_breaks[2], col='green') + geom_vline(xintercept=s_breaks[3], col='red') + geom_vline(xintercept=s_breaks[4], col='blue')

grid.arrange(accrate_plot, unique_part_plot, p_s1, p_s2, p_s3, p_s4, nrow = 3)

#ggplot() + 
#  geom_histogram(data=data.frame(v=rLaplace(3500,1.2,0.1)), aes(x=v,y=..density..), alpha=0.3, lwd=0.2) +
#  xlim(qLaplace(0.001, 1.2,0.1), qLaplace(0.999, 1.2,0.1)) +
#  stat_function(fun = function(x) dLaplace(x, 1.2, 0.1), col='black')

#p_s3
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-------------------------- Compare different PropCOV_mult: ----------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

AR.1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_160110_alpha0.9_4params_M60_Np10000_alpha0.9_PropCOVmult0.1')
AR.01 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_175418_alpha0.9_4params_M60_Np10000_alpha0.9_PropCOVmult0.01')
AR1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2024-03-26_171206_alpha0.9_4params_M60_Np10000_alpha0.9_PropCOVmult1')

plot(AR.1$tolerancestore, ylim=c(3, 5.5))
points(AR.01$tolerancestore, col='blue')
points(AR1$tolerancestore, col='red')
#has almost no effect on the descent in tolerances ??

plot(AR.1$accrate)
points(AR.01$accrate, col='blue')
points(AR1$accrate, col='red')
abline(h=0.05)
#doesn't make a huge difference to acceptance rate. 

s_finish <- 50
unique_part.1 <- c()
unique_part.01 <- c()
unique_part1 <- c()
for (i in 1:s_finish){
  unique_part.1 <- c(unique_part.1, length(unique(AR.1$Omega_sample[AR.1$W[,i]>0,1,i])))
  unique_part.01 <- c(unique_part.01, length(unique(AR.01$Omega_sample[AR.01$W[,i]>0,1,i])))
  unique_part1 <- c(unique_part1, length(unique(AR1$Omega_sample[AR1$W[,i]>0,1,i])))
}
plot(unique_part.1, ylab='Unique Alive Part')
points(unique_part.01, col='blue')
points(unique_part1, col='red')

ggplot() + 
  geom_histogram(data=data.frame(v=AR.1$Omega_sample_phys[,1,s_finish], w=AR.1$W[,s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.2, col='blue', fill='blue', lwd=0.2) +
  #geom_histogram(data=data.frame(v=AR1$Omega_sample_phys[,1,s_finish], w=AR1$W[,s_finish]), aes(x=v,y=..density.., weight = w), alpha=0.2, fill='red', col='red', lwd=0.2) +
  geom_histogram(data=data.frame(v=AR.01$Omega_sample_phys[,1,s_finish], w=AR.01$W[,s_finish3]), aes(x=v,y=..density.., weight = w), alpha=0.3, fill='yellow', col='yellow', lwd=0.2) +
  geom_vline(xintercept=1, col='black')












#-----------------------------------------------------------------------------------------



plot(density(AlgoResults$Omega_sample[,4,1]))
lines(density(AlgoResults$Omega_sample[,4,50]), col='blue')

# 10000 particles, alpha 0.9: everything seems to work
# 1000 particles, alpha 0.9: propcov_mult 0.01 seems to be pretty haphazard, whereas propcov_mult 0.1 seems to work better
#

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
  