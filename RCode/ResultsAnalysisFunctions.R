
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-07-30_185759'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-07-22_193115'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-08-04_134932'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-08-11_164714'
# results_file <-'/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-08-16_160839'
# results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-08-28_192239'
# results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-08-29_170442'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-09-03_171706'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-09-13_182006'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2023-11-02_170925'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2023-11-04_225255'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2023-12-08_091753'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2023-12-11_020113'
#results_file <-'/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-01-06_090554'
# results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-01-08_102907_3'
# results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_start_step_2024-01-17_214246'
# 
# AlgoResults <- readRDS(results_file)
# AlgoResults$sampled_full <- NULL

compare_sampled_full_observed <- function(){
  Omega_min.d <- AlgoResults$Omega_sample_phys[1,,1] %>% relist(skeleton=Model$skeleton)
  impact_sample <- SampleImpact(dir, Model, Omega_min.d  %>% addTransfParams(), AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), 1))
  #observed <- impact_sample$poly[[1]]$observed
  #impact_type <- impact_sample$poly[[1]]$impact
  
  p <- 346
  sampled <- AlgoResults$sampled_full[p,,,AlgoResults$s_finish]
  plot_df <- data.frame(sampled)
  names(plot_df) <- paste0('sampled', 1:NCOL(plot_df))
  plot_df$impact_type=impact_sample$poly[[1]]$impact
  plot_df$observed=impact_sample$poly[[1]]$observed
  plot_df$event_id=impact_sample$poly[[1]]$event_id
  
  ggplot(plot_df %>% filter(impact_type=='buildDam' & event_id==135), aes(x=log(observed+1), y=log(sampled3+1))) + 
    geom_point() + geom_abline(slope=1, intercept=0)
  
}

check_changes_1_event <- function(){
  ODD <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/EQ20191016PHL_124')
  Omega <- list(Lambda1 = list(nu=9.6, kappa=1.9049508),
       Lambda2 = list(nu=9.8626452, kappa=1.9049508),
       Lambda3 = list(nu=7.8, kappa=3),
       Lambda4 = list(nu=9.9, kappa=1),
       theta= list(theta1=0.6),
       eps = list(local=0.001, hazard_mort=0.8383464, hazard_disp=0.9, hazard_bd=0.001, hazard_cor=0.55),
       vuln_coeff = list(PDens=0, SHDI=0, GNIc=0, Vs30=0, EQFreq=0, FirstHaz=0, Night=0, FirstHaz.Night=0),
       check = list(check=0.5))
  impact <- DispX(ODD, Omega %>% addTransfParams(), Model$center, Model$BD_params, AlgoParams %>% replace(which(names(AlgoParams)==c('Np')), 5), output='SampledAgg')
  plot_df <- impact[[1]][, c('observed', 'impact')]
  for (i in 1:length(impact)){
    plot_df[,paste0('sampled', i)] <- impact[[i]]$sampled
  }
  plot_df$sampled_min <- apply(plot_df[,grep("sampled", names(plot_df))], 1, min)
  plot_df$sampled_max <- apply(plot_df[,grep("sampled", names(plot_df))], 1, max)
  plot_df$sampled_median <- apply(plot_df[,grep("sampled", names(plot_df))], 1, median)
  ggplot(plot_df %>% filter(impact=='buildDam'), aes(x=log(observed+1), y=log(sampled_median+1), ymin=log(sampled_min+1), ymax=log(sampled_max+1))) + 
    geom_point() + geom_abline(slope=1, intercept=0) #+ geom_errorbar()
  
  ints <- seq(6,8, 0.05)
  plot(ints, D_Dam_calc(h_0(ints,Model$I0,Omega %>% addTransfParams()), Omega %>% addTransfParams()))
  
  bd_impacts <- data.frame(buildings=numeric(), intensity=numeric(), observed=numeric())
  for (i in 1:NROW(ODD@impact)){
    if(ODD@impact$impact[i] != 'buildDam') next
    buildings <- sum(ODD$nBuildings[ODD@polygons[[ODD@impact$polygon[i]]]$indexes])
    intensity <- max(ODD$hazMean1[ODD@polygons[[ODD@impact$polygon[i]]]$indexes], na.rm=T)
    observed <- ODD@impact$observed[i]
    bd_impacts %<>% add_row(buildings=buildings, intensity=intensity, observed=observed)
  }
  plot(bd_impacts$intensity, jitter(bd_impacts$observed/bd_impacts$buildings))
  which(bd_impacts$intensity > 7 & bd_impacts$observed==0)
}

# plot(ints, h_0(ints,Model$I0,Omega %>% addTransfParams()))
# h_0<-function(I,I0, Omega){
#   I0 <- 7
#   ind<-I>I0
#   h<-rep(0,length(I))
#   #h[ind]<- I[ind]-I0
#   h[ind] <- exp(Omega$theta$theta1*I[ind])
#   return(h)
# }
# points(ints, h_0(ints,Model$I0,Omega %>% addTransfParams()), col='red')

library(zoo)

addAlgoParams <- function(AlgoResults, s_finish=NULL){
  if(is.null(s_finish)){
    AlgoResults$s_finish <- which(is.na(AlgoResults$Omega_sample_phys[1,1,]))[1]-1
    if(is.na(AlgoResults$s_finish)){
      acc = c()
      for (s in 2:length(AlgoResults$accrate_store)){
        W_nonzero <- which(AlgoResults$W[,s-1] != 0)
        acc <- c(acc, sum(AlgoResults$Omega_sample_phys[W_nonzero,1,s] != AlgoResults$Omega_sample_phys[W_nonzero,1,s-1])/length(W_nonzero))
      }
      AlgoResults$s_finish <- which(rollapply(acc<0.025, 5, all))[1]
    }
    if (is.na(AlgoResults$s_finish)){
      AlgoResults$s_finish <- length(AlgoResults$Omega_sample_phys[1,1,])
    }
  } else {
    AlgoResults$s_finish <- s_finish
  }
  #AlgoResults$s_finish <- 220
  AlgoResults$n_x <- length(unlist(Model$skeleton))
  AlgoResults$Npart <- NROW(AlgoResults$W)
  return(AlgoResults)
}

plot_density_vs_step = function(AlgoResults, Omega=NULL, s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  par(mfrow=c(5,4))
  plot_titles <- names(unlist(Model$skeleton))
  plot_titles[1:8] <- paste0(ifelse(1:8 %% 2 == 0, 'sigma_', 'mu_'), rep(1:4,each=2))
  for (i in 1:length(unlist(Model$skeleton))){
    ymin= min(AlgoResults$Omega_sample_phys[,i,1:AlgoResults$s_finish])
    ymax= max(AlgoResults$Omega_sample_phys[,i,1:AlgoResults$s_finish])
    if(!is.null(Omega)){
      ymin <- min(ymin, unlist(Omega)[i])
      ymax <- max(ymax, unlist(Omega)[i])
    }
    plot(rep(1:AlgoResults$s_finish,each=AlgoResults$Npart), AlgoResults$Omega_sample_phys[,i,1:AlgoResults$s_finish], ylim=c(ymin, ymax), main=plot_titles[i], xlab='Step', ylab='')
    if (!is.null(Omega)){
      abline(h=unlist(Omega)[i], col='red')
    }
  }
  par(mfrow=c(1,1)) # 800 x 800
}

plot_d_vs_step = function(AlgoResults, ymax=NULL, s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  
  ymin=min(AlgoResults$d, na.rm=T)
  ymax=ifelse(is.null(ymax), max(AlgoResults$d[which(is.finite(AlgoResults$d))], na.rm=T), ymax)
  plot(rep(1, AlgoResults$Npart), apply(adrop(AlgoResults$d[,,1, drop=F], drop=3), 1, min, na.rm=T), xlim=c(1, AlgoResults$s_finish), xlab='Algorithm step', ylab='Mean energy score for each particle', ylim=c(ymin, ymax), cex=0.5)
  for (s in 2:AlgoResults$s_finish){
    nonzero_weights <- which(AlgoResults$W[,s] != 0)
    points(rep(s, length(nonzero_weights)), apply(adrop(AlgoResults$d[nonzero_weights,,s, drop=F], drop=3), 1, min, na.rm=T), cex=0.5)
  }
}

plot_mean_d_vs_step = function(AlgoResults, ymax=NULL, s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  
  ymin=min(AlgoResults$d, na.rm=T)
  ymax=ifelse(is.null(ymax), max(AlgoResults$d[which(is.finite(AlgoResults$d))], na.rm=T), ymax)
  plot(1, mean(apply(adrop(AlgoResults$d[,,1, drop=F], drop=3), 1, min, na.rm=T)), xlim=c(1, AlgoResults$s_finish), xlab='Algorithm Step', ylab='Mean Energy Score', ylim=c(ymin, 35))
  for (s in 2:AlgoResults$s_finish){
    nonzero_weights <- which(AlgoResults$W[,s] != 0)
    points(s,mean(apply(adrop(AlgoResults$d[nonzero_weights,,s, drop=F], drop=3), 1, min, na.rm=T)))
  }
}


plot_acc_prob = function(AlgoResults, s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  
  acc = c()
  for (s in 2:AlgoResults$s_finish){
    W_nonzero <- which(AlgoResults$W[,s-1] != 0)
    acc <- c(acc, sum(AlgoResults$Omega_sample_phys[W_nonzero,1,s] != AlgoResults$Omega_sample_phys[W_nonzero,1,s-1])/length(W_nonzero))
  }
  plot(acc, ylim=c(0,0.5))
}

get_greek_titles = function(var_name){
  if (var_name == 'Lambda1.nu') return(expression(mu['Disp']))
  if (var_name == 'Lambda1.kappa') return(expression(kappa['Disp']))
  if (var_name == 'Lambda2.nu') return(expression(mu['Mort']))
  if (var_name == 'Lambda2.kappa') return(expression(kappa['Mort']))
  if (var_name == 'Lambda3.nu') return(expression(mu['BuildDam']))
  if (var_name == 'Lambda3.kappa') return(expression(kappa['BuildDam']))
  if (var_name == 'eps.hazard_mort') return(expression(sigma['Mort']))
  if (var_name == 'eps.hazard_disp') return(expression(sigma['Disp']))
  if (var_name == 'eps.hazard_bd') return(expression(sigma['BuildDam']))
  if (var_name == 'eps.hazard_cor') return(expression(rho))
  if (var_name == 'eps.local') return(expression(sigma['Local'['Mort']]))
  if (var_name == 'vuln_coeff.PDens') return(expression(paste(beta[2], '  (PDens)')))
  if (var_name == 'vuln_coeff.Vs30') return(expression(paste(beta[1], '  (Vs30)')))
  if (var_name == 'vuln_coeff.SHDI') return(expression(paste(beta[3], '  (SHDI)')))
  if (var_name == 'vuln_coeff.EQFreq') return(expression(paste(beta[5], '  (EQFreq)')))
  if (var_name == 'vuln_coeff.GNIc') return(expression(paste(beta[4], '  (GNIc)')))
  if (var_name == 'vuln_coeff.FirstHaz') return(expression(paste(beta[6], '  (FirstHaz)')))
  if (var_name == 'vuln_coeff.Night') return(expression(paste(beta[7], '  (Night)')))
  if (var_name == 'vuln_coeff.FirstHaz.Night') return(expression(paste(beta[8], '  (FirstHaz x Night)')))
  if (var_name == 'theta.theta1') return('Dummy 1')
  if (var_name == 'check.check') return('Dummy 2')
  else return(var_name)
    #expression(paste("Diameter of aperture ( ", mu, " )")),xlim=c(xmin, xmax),
    #
}

plot_correlated_posteriors = function(AlgoResults, include_priors=T, Omega=NULL,
                                      pairings=rbind(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)), s_finish=NULL,
                                      AlgoResultsMCMC=NULL, subfig_title_adj=-0.4){
  
  AlgoResults %<>% addAlgoParams(s_finish)
  post_samples <- AlgoResults$Omega_sample_phys[which(AlgoResults$W[,AlgoResults$s_finish]>0),,AlgoResults$s_finish]
  if (include_priors) prior_samples <- AlgoResults$Omega_sample_phys[,,1]
  

  if(NROW(pairings)>7){
    par(mfrow=c(3,4), mai = c(0.6, 0.7, 0.2, 0.1), family='Liberation Serif', cex.lab=1.25)
  } else {
    par(mfrow=c(2,4), mai = c(0.6, 0.7, 0.2, 0.1), family='Liberation Serif', cex.lab =1.25)
  }
  for (p in 1:NROW(pairings)){
    xmin= min(post_samples[,pairings[p,1]]); xmax= max(post_samples[,pairings[p,1]])
    ymin= min(post_samples[,pairings[p,2]]); ymax= max(post_samples[,pairings[p,2]])
    if(include_priors){
      xmin <- min(xmin, prior_samples[,pairings[p,1]]); xmax <- max(xmax, prior_samples[,pairings[p,1]])
      ymin <- min(ymin, prior_samples[,pairings[p,2]]); ymax <- max(ymax, prior_samples[,pairings[p,2]])
    }
    
    
    if (include_priors){
      plot(prior_samples[,pairings[p,1]], prior_samples[,pairings[p,2]], col='blue', 
           xlab=get_greek_titles(names(unlist(Model$skeleton))[pairings[p,1]]), xlim=c(xmin, xmax),
           ylab=get_greek_titles(names(unlist(Model$skeleton))[pairings[p,2]]), ylim=c(ymin, ymax), cex=0.75)
      if (!is.null(AlgoResultsMCMC) & p < 7){
        MCMC_s_finish <- which(is.infinite(AlgoResultsMCMC$loss))[1]-1
        post_s <- ceil(MCMC_s_finish/2):MCMC_s_finish
        post_s <- post_s[which(AlgoResultsMCMC$HLP_vals[post_s]==0 & !is.na(AlgoResultsMCMC$HLP_vals[post_s]))]
        MCMC_samples <- post_s[round(seq(1, length(post_s), length.out=NROW(AlgoResults$Omega_sample_phys[,1,1])))] #round(seq(MCMC_s_finish/2, MCMC_s_finish, length.out=NROW(AlgoResults$Omega_sample_phys[,1,1])))
        points(AlgoResultsMCMC$Omega_sample_phys[pairings[p,1],MCMC_samples], AlgoResultsMCMC$Omega_sample_phys[pairings[p,2],MCMC_samples], col='red')
      }
      points(post_samples[,pairings[p,1]], post_samples[,pairings[p,2]], cex=0.75,col='black')
      if (!is.null(AlgoResultsMCMC) & p >= 7){
        MCMC_s_finish <- which(is.infinite(AlgoResultsMCMC$loss))[1]-1
        post_s <- ceil(MCMC_s_finish/2):MCMC_s_finish
        post_s <- post_s[which(AlgoResultsMCMC$HLP_vals[post_s]==0 & !is.na(AlgoResultsMCMC$HLP_vals[post_s]))]
        MCMC_samples <- post_s[round(seq(1, length(post_s), length.out=NROW(AlgoResults$Omega_sample_phys[,1,1])))] #round(seq(MCMC_s_finish/2, MCMC_s_finish, length.out=NROW(AlgoResults$Omega_sample_phys[,1,1])))
        points(AlgoResultsMCMC$Omega_sample_phys[pairings[p,1],MCMC_samples], AlgoResultsMCMC$Omega_sample_phys[pairings[p,2],MCMC_samples], col=alpha('red', 0.3))
      }
    } else {
      plot(post_samples[,pairings[p,1]], post_samples[,pairings[p,2]], 
           xlab=names(unlist(Model$skeleton))[pairings[p,1]], xlim=c(xmin, xmax),
           ylab=names(unlist(Model$skeleton))[pairings[p,2]], ylim=c(ymin, ymax))
    }
  
    if (!is.null(Omega)){
      points(unlist(Omega)[pairings[p,1]], unlist(Omega)[pairings[p,2]], col='red', pch=4, cex=1.5, lwd=2.5)
    }
    subfig_label <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l')[p]
    mtext(paste0('(', subfig_label, ')'), side = 3, line = 0, adj = subfig_title_adj, cex = 1, font = 1)
  }
  par(mfrow=c(1,1), mai=c(1,1,1,1))
}


plot_vuln_posteriors = function(AlgoResults, include_priors=T, Omega=NULL,
                                      pairings=rbind(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)), s_finish=NULL){
  
  AlgoResults %<>% addAlgoParams(s_finish)
  post_samples <- AlgoResults$Omega_sample_phys[,,AlgoResults$s_finish]
  if (include_priors) prior_samples <- AlgoResults$Omega_sample_phys[,,1]
  
  vuln_var <- grep('vuln', names(unlist(Model$skeleton)))

  vuln_var <- vuln_var[c(4,1,2,3,5,6,7,8)]
  plot_list <- list()
  
  # Define subfigure labels (a), (b), (c), ...
  subfig_labels <- letters[1:8]
  
  # Loop through the variables and generate each plot
  for (i in seq_along(vuln_var)) {
    var <- vuln_var[i]
    
    
    p_H_given_y = sum(post_samples[,var]>0)/length(post_samples[,var])
    p_H = sum(prior_samples[,var]>0)/length(prior_samples[,var])
    print(paste('Post Weight Positive', names(unlist(Model$skeleton))[var], ':', round(p_H_given_y, 3)))
    print(paste('Bayes factor for', names(unlist(Model$skeleton))[var], ':', round(p_H_given_y * (1-p_H)/((1-p_H_given_y)*p_H), 3)))
    
    post_median = median(post_samples[,var])
    
    # Create each plot and store it in the plot_list
    plot <- ggplot() +
      # Add the histogram from the posterior samples
      geom_histogram(aes(x = post_samples[, var], y = ..density..), 
                     bins = 30, fill = "lightgrey", color='black', lwd=0.5) +
      geom_density(aes(x = prior_samples[, var]), color = "blue", size = 0.7) +
      geom_vline(xintercept = 0, color = "blue", linetype = "dashed", size = 0.7) +
      geom_vline(xintercept = post_median, color = "black", linetype = "dashed", size = 0.7) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = get_greek_titles(names(unlist(Model$skeleton))[var])) +
      xlim(c(-1, 1)) +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),  # Change y-axis title font
        axis.text.x = element_text(family = "Times New Roman", size = 12),   # Change x-axis text font
        axis.text.y = element_blank(),  # Remove y-axis text (numbers)
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.title.x = element_text(family = "Times New Roman", size = 12),  # Change x-axis title font
        plot.title = element_text(family = "Times New Roman", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = unit(c(0, 20, 0, 15), "pt")
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
  

}


plot_corr_posterior_vs_d = function(AlgoResults, include_priors=F, Omega=NULL,
                                      pairing=c(1,2), s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  post_samples <- AlgoResults$Omega_sample_phys[,,AlgoResults$s_finish]
  if (include_priors) prior_samples <- AlgoResults$Omega_sample_phys[,,1]
  
  xmin= min(post_samples[,pairing[1]]); xmax= max(post_samples[,pairing[1]])
  ymin= min(post_samples[,pairing[2]]); ymax= max(post_samples[,pairing[2]])
  if(include_priors){
    xmin <- min(xmin, prior_samples[,pairing[p,1]]); xmax <- max(xmax, prior_samples[,pairing[p,1]])
    ymin <- min(ymin, prior_samples[,pairing[p,2]]); ymax <- max(ymax, prior_samples[,pairing[p,2]])
  }
  plot_df <- data.frame(x=post_samples[,pairing[1]],
                        y=post_samples[,pairing[2]],
                        d=AlgoResults$d[,1,AlgoResults$s_finish])
  
  ggplot(plot_df, aes(x=x, y=y, col=d)) + geom_point() + 
    geom_point(aes(x=unlist(Omega)[pairing[1]], y=unlist(Omega)[pairing[2]]), col='green', shape=4, size=3) + 
    xlab(names(unlist(Omega))[pairing[1]]) + ylab(names(unlist(Omega))[pairing[2]]) +
    scale_color_viridis_c(option = "magma", direction=-1)
  
  # cr <- colorRamp(c("red", "blue"))
  # col = rgb(cr(AlgoResults$d[,1,AlgoResults$s_finish] / max(AlgoResults$d[,1,AlgoResults$s_finish])), max=255)
  # 
  # plot(post_samples[,pairing[1]], post_samples[,pairing[2]], 
  #        xlab=names(unlist(Model$skeleton))[pairing[1]], xlim=c(xmin, xmax),
  #        ylab=names(unlist(Model$skeleton))[pairing[2]], ylim=c(ymin, ymax), col=col)
  # 
  # if (include_priors) points(prior_samples[,pairing[1]], prior_samples[,pairing[2]], col='blue')
  # if (!is.null(Omega)){
  #   points(unlist(Omega)[pairing[1]], unlist(Omega)[pairing[2]], col='red', pch=4, cex=2, lwd=4)
  # }
}

plot_corr_transf_posterior_vs_d = function(AlgoResults, include_priors=F, Omega=NULL,
                                    pairing=c(1,2), s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  relist_transf <- function(Omega_phys){
    unlist(addTransfParams(relist(Omega_phys, skeleton=Model$skeleton)))
  }
  post_transf <- apply(AlgoResults$Omega_sample_phys[,,AlgoResults$s_finish], 1,relist_transf)
  prior_samples <- apply(AlgoResults$Omega_sample_phys[,,1], 1,relist_transf)
  
  xmin= min(post_transf[pairing[1],]); xmax= max(post_transf[pairing[1],])
  ymin= min(post_transf[pairing[2],]); ymax= max(post_transf[pairing[2],])
  if(include_priors){
    xmin <- min(xmin, prior_samples[pairing[1],]); xmax <- max(xmax, prior_samples[pairing[1],])
    ymin <- min(ymin, prior_samples[pairing[2],]); ymax <- max(ymax, prior_samples[pairing[2],])
  }
  plot_df <- data.frame(x=post_transf[pairing[1],],
                        y=post_transf[pairing[2],],
                        d=AlgoResults$d[,1,AlgoResults$s_finish])
  
  ggplot(plot_df, aes(x=x, y=y, col=d)) + geom_point() + 
    geom_point(aes(x=prior_samples[pairing[1],], y=prior_samples[pairing[2],]), col='grey', alpha=0.5) + 
    geom_point(aes(x=unlist(addTransfParams(Omega))[pairing[1]], y=unlist(addTransfParams(Omega))[pairing[2]]), col='green', shape=4, size=3) + 
    xlab(names(unlist(addTransfParams(Omega)))[pairing[1]]) + ylab(names(unlist(addTransfParams(Omega)))[pairing[2]]) +
    scale_color_viridis_c(option = "magma", direction=-1)+ theme_bw()
  
  # cr <- colorRamp(c("red", "blue"))
  # col = rgb(cr(AlgoResults$d[,1,AlgoResults$s_finish] / max(AlgoResults$d[,1,AlgoResults$s_finish])), max=255)
  # 
  # plot(post_samples[,pairing[1]], post_samples[,pairing[2]], 
  #        xlab=names(unlist(Model$skeleton))[pairing[1]], xlim=c(xmin, xmax),
  #        ylab=names(unlist(Model$skeleton))[pairing[2]], ylim=c(ymin, ymax), col=col)
  # 
  # if (include_priors) points(prior_samples[,pairing[1]], prior_samples[,pairing[2]], col='blue')
  # if (!is.null(Omega)){
  #   points(unlist(Omega)[pairing[1]], unlist(Omega)[pairing[2]], col='red', pch=4, cex=2, lwd=4)
  # }
}

crps <- function(sample, obs){
  sample <- sort(sample)
  m <- length(sample)
  crps <- 0
  for (i in 1:m){
    crps <- crps + (sample[i]-obs)*(m*as.numeric(obs<sample[i]) -i + 0.5)
  }
  crps <- (crps*2)/(m^2)
  return(crps)
}

sample_post_predictive <- function(AlgoResults, M, s, dat='Train', single_particle=F, Omega=NULL, particle_i = NULL, 
                                   return_type='Poly', output='SampledAgg'){
  AlgoParams$input_folder <- AlgoResults$input_folder
  if (is.null(Omega)){
    sampled_part <- sample(1:AlgoResults$Npart, M, prob=AlgoResults$W[, s], replace=T)
    impact_sample <- SampleImpact(dir = dir,Model = Model,
                                  proposed = AlgoResults$Omega_sample_phys[sampled_part[1],,s] %>% relist(Model$skeleton) %>% addTransfParams(), 
                                  AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                  dat=dat, output=output)
    poly_sampled <- impact_sample$poly[[1]][,c('event_id', 'iso3', 'polygon', 'impact', 'observed', 'sampled')] # 'sdate', 'inferred',
    point_sampled <- impact_sample$point
    
    if (M>1){
      for (m in 2:M){
        print(m)
        impact_sample <- SampleImpact(dir = dir,Model = Model,
                                      proposed = AlgoResults$Omega_sample_phys[sampled_part[m],,s] %>% relist(Model$skeleton) %>% addTransfParams(), 
                                      AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                      dat=dat, output=output)
        poly_sampled <- cbind(poly_sampled, impact_sample$poly[[1]]$sampled)
        point_sampled <- cbind(point_sampled, impact_sample$point)
      }
    }
  } else {
    if(is.null(Omega)){
      sampled_part <- ifelse(is.null(particle_i), sample(1:AlgoResults$Npart, 1, prob=AlgoResults$W[, s],replace=T), particle_i)
      Omega <- AlgoResults$Omega_sample_phys[sampled_part[1],,s] %>% relist(Model$skeleton) 
    } 
    impact_sample <- SampleImpact(dir = dir,Model = Model,
                                  proposed = Omega %>% addTransfParams(), 
                                  AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                  dat=dat, output=output)
    poly_sampled <- impact_sample$poly[[1]][,c('event_id', 'iso3', 'polygon', 'impact', 'observed', 'sampled')] #'sdate', 'inferred'
    point_sampled <- impact_sample$point
    
    if (M>1){
      for (m in 2:M){
        poly_sampled <- cbind(poly_sampled, impact_sample$poly[[m]]$sampled)
      }
    }
  }
  
  df_poly <- poly_sampled
  names(df_poly)[grep("sampled", names(df_poly))] = paste0('sampled.', 1:M)
  df_poly$obs_id[order(df_poly$observed)] <- 1:NROW(df_poly)
  df_poly$crps <- NA
  
  for (i in 1:NROW(df_poly)){
    if (is.na(df_poly[i,'sampled.1'])) next 
    df_poly$crps[i] <- crps(log(as.numeric(df_poly[i, grep("sampled", names(df_poly))])+10), log(df_poly$observed[i]+10)) * 
      unlist(AlgoParams$impact_weights)[df_poly$impact[i]]
  }
  
  if (tolower(dat)=='train'){
    df_poly$train_flag = 'TRAIN'
  } else if (tolower(dat)=='test'){
    df_poly$train_flag = 'TEST'
  } else {
    folderin<-paste0(dir, AlgoResults$input_folder, "ODDobjects/")
    ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
    ufiles_train <- grep('^Train/' , ufiles, value = TRUE)
    i_train <- as.numeric(gsub(".*_(\\d+)", "\\1", ufiles_train))
    df_poly$train_flag <- ifelse(df_poly$event_id %in% i_train, 'TRAIN', 'TEST')
  }
  
  if(tolower(return_type)=='poly'){
    return(df_poly)
  } else {
    df_point = point_sampled
    colnames(df_point) = 'count'
    df_point <- cbind(df_point, 0)
    
    df_point[,2] = ifelse(rownames(df_point) %in% c('N12', 'N21', 'N23', 'N32'), 0.5 * df_point[,1], ifelse(rownames(df_point) %in% c('N13', 'N31'), df_point[,1], 0)) * 0.1
    return(list(poly=df_poly, point=df_point))
  }
}

plot_posteriors <- function(AlgoResults, pars_i=NULL, Omega=NULL, s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  n_params <- NCOL(AlgoResults$Omega_sample_phys[,,1])
  #ncol = 4
  #nrow = n_params %/% ncol+1
  plots_list <- list()
  plot_df <- data.frame(AlgoResults$Omega_sample_phys[,, AlgoResults$s_finish])
  priorplot_df <- data.frame(AlgoResults$Omega_sample_phys[,, 1])
  colnames(plot_df) <- names(priorplot_df) <- names(Omega_true %>% unlist())
  j <- 1
  if (is.null(pars_i)) pars_i <- 1:n_params
  for(i in pars_i){
    var <- colnames(plot_df)[i]
    plots_list[[j]] <- ggplot(plot_df, aes(x = !!ensym(var)))  +
      geom_histogram(aes(y= ..density..),color = "black", alpha = 0.7) +
      ggtitle(var)+
      #labs(x='', y='') + 
      #scale_fill_discrete("Legend", values = dd.col) + 
      #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      geom_density(data=priorplot_df, aes(x=!!ensym(var), y=..density..), col='blue')
    if (!is.null(Omega)){plots_list[[j]] = plots_list[[j]] + geom_vline(xintercept = unlist(Omega)[var], col='red')}
      #scale_x_log10()+
    j <- j+1
  }
  do.call(grid.arrange,plots_list)
}

sample_post_predictive_BDam <- function(AlgoResults, M, s, dat='Train', single_particle=F, particle_i = NULL, Omega=NULL){
  
  if (!single_particle){
    sampled_part <- sample(1:AlgoResults$Npart, M, prob=AlgoResults$W[, s], replace=T)
    proposed <- AlgoResults$Omega_sample_phys[sampled_part[1],,s] %>% relist(Model$skeleton) %>% addTransfParams()
    AlgoParams %>% replace(which(names(AlgoParams)==c('Np')), 1)
    
    folderin<-paste0(dir, AlgoResults$input_folder, "BDobjects/")
    ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
    if (tolower(dat)=='train'){
      ufiles <- grep('^Train/' , ufiles, value = TRUE)
    } else if (tolower(dat)=='test'){
      ufiles <- grep('^Test/' , ufiles, value = TRUE)
    }
    x <- file.info(paste0(folderin,ufiles))
    ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
    
    plot_dat <- data.frame(event_name=character(), 
                           obs_p_dam = numeric(), 
                           mean_p_mod = numeric(), 
                           mean_p_modDam = numeric(),
                           mean_p_modUnaff = numeric())
    
    for (filer in ufiles){
      BDy<-readRDS(paste0(folderin,filer))
      
      if(is.null(BDy)){return(rep(0, AlgoParams$Np))}
      if(nrow(BDy@data)==0){return(rep(0, AlgoParams$Np))}
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      #BDy@fIndies<-Model$fIndies
      # Apply BDX
      AlgoParams$Np <- 1
      BD_with_p <- BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, output='results_analysis')
      
      BD_with_p %<>% as.data.frame()
      names(BD_with_p) <- c('Longitude', 'Latitude', 'grading', 'Intensity','pDam')
      for (var_nam in c('Longitude', 'Latitude', 'Intensity', 'pDam')){
        BD_with_p[,var_nam] <- as.numeric(BD_with_p[,var_nam])
      } 
      BD_with_p$pDam <- 1-BD_with_p$pDam
      BD_with_p$order <- 1:NROW(BD_with_p)
      
      BD_with_p$roundedLongitude = round(BD_with_p$Longitude, 1)
      BD_with_p$roundedLatitude = round(BD_with_p$Latitude, 1)
      BD_with_p$roundedIntensity = round(BD_with_p$Intensity, 1)
      
      BD_with_p <- arrange(BD_with_p, Latitude)
      breaks_condition <- c(TRUE, diff(sort(BD_with_p$Latitude)) > 0.01)
      BD_with_p$groupLat <- cumsum(breaks_condition)
      
      BD_with_p <- arrange(BD_with_p, Longitude)
      breaks_condition <- c(TRUE, diff(sort(BD_with_p$Longitude)) > 0.01)
      BD_with_p$groupLon <- cumsum(breaks_condition)
      
      #BD_with_p <- arrange(BD_with_p, Intensity)
      #breaks_condition <- c(TRUE, diff(sort(BD_with_p$Intensity)) > 0.1)
      BD_with_p$groupp <- round(BD_with_p$Intensity / 0.1) * 0.1
      
      BD_with_p_grouped <- BD_with_p %>% group_by(groupLon, groupLat) %>% 
        summarise(obs_p_dam = mean(grading!='notaffected'), 
                  mean_p_mod = mean(pDam), 
                  mean_p_modDam = mean(pDam[which(grading!='notaffected')]),
                  mean_p_modUnaff = mean(pDam[which(grading=='notaffected')]))
      
      BD_with_group <- BD_with_p %>% group_by(groupLon, groupLat) %>% mutate(group = cur_group_id())
      BD_with_group <- arrange(BD_with_group, order)
      
      
      plot_dat <- plot_dat %>% add_row(BD_with_p_grouped[,-c(1,2,3)] %>% add_column(event_name=filer))
      
      
      # BS = mean(ifelse(BD_with_p$grading!='notaffected', (BD_with_p$pDam-1)^2, BD_with_p$pDam^2))
      # print(filer)
      # print(BS)
      # print(mean(BD_with_p$pDam[BD_with_p$grading=='notaffected']))
      # print(mean(BD_with_p$pDam[BD_with_p$grading!='notaffected']))
      # print(mean(BD_with_p$grading!='notaffected'))
      BD_with_p$grading[which(BD_with_p$grading != 'notaffected' & BD_with_p$grading != 'possible')] = 'Damaged'
      ggplot(BD_with_p %>% arrange(desc(grading)), aes(x=Intensity, y=pDam, col=grading)) + geom_jitter(width=0.05, height=0.05) + ylab('Modelled Probability of Damage')
      # #ggplot(BD_with_p %>% filter(grading!='notaffected'), aes(x=Intensity, y=pDam, col=grading)) + geom_jitter(width=0.1, height=0.05)
      # 
      # plot_bd_df_int <- data.frame(intensity=numeric(), bs_prop_bd=numeric(), tot_obs=integer())
      # increment <- 0.2
      # for(int in seq(5.4,7.1, increment)){
      #   plot_bd_df_int %<>% add_row(intensity=int + increment/2, 
      #                               bs_prop_bd=1-length(which(BD_with_p$grading[which(BD_with_p$Intensity>int & BD_with_p$Intensity<(int+increment))]=='notaffected'))/length(which(BD_with_p$Intensity>int & BD_with_p$Intensity<(int+increment))),
      #                               tot_obs=length(which(BD_with_p$Intensity>int & BD_with_p$Intensity<(int+increment))))
      # }
      
      
    }
  }
}

compare_CRPS_breakdown_at_different_s <- function(AlgoResults, M, s_finish=NULL){
  M <- 1
  AlgoResults %<>% addAlgoParams(s_finish)
  df_poly1 <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, single_particle=T, particle_i = 1)
  df_poly2 <- sample_post_predictive(AlgoResults, M, 20, single_particle=T, particle_i = 1)
  
  
  ix <- which(df_poly1$impact=='mortality')
  ymin=min(df_poly1$crps[ix], df_poly2$crps[ix])
  ymax=max(df_poly1$crps[ix], df_poly2$crps[ix])
  plot(df_poly1$crps[ix], xlim=c(310, 350))
  points(df_poly2$crps[ix], col='red')
  
}

plot_predictive <- function(AlgoResults, dat='Train', s_finish=NULL){
  
  M <- 10
  AlgoResults %<>% addAlgoParams(s_finish)
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat=dat, single_particle=T, particle=particle_min.d[1])
  
  impact_type <- 'displacement'
  k <- 10
  ix <- which(df_poly$impact==impact_type)
  inferred_cols <- ifelse(df_poly$inferred[ix], 'red', 'black')
  plot(log(df_poly$observed[ix]+k), log(df_poly$sampled.1[ix]+k), 
       xlab=paste0('log(observed+',k,')'), ylab=paste0('log(sampled+',k,')'), col=inferred_cols)
  abline(0,1)
  
  #large_discrepancies
  df_poly[which(abs(log(df_poly$sampled.1+1)-log(df_poly$observed+1)) >7),] %>% filter(impact==impact_type)
  
  #plot covariates against discrepancy between sampled and observed
  df_poly$sampled_median <- apply(df_poly[grep("sampled", names(df_poly))], 1, median)
  df_poly <- add_covar(df_poly, covar='EQFreq', dat=dat)
  plot(df_poly$EQFreq, df_poly$sampled_median)
  
  abline(0,1)
}

create_df_postpredictive <- function(AlgoResults, single_particle = F, Omega = NULL, particle_best=T, M=50,
                                     s_finish=NULL, output='SampledAgg'){
  # If single_particle == T:
  #       If Omega = NULL, finds the parameters that produced the smallest distance and samples using only these parameters
  #       If Omega is given, samples uses Omega
  # If single_particle == F, samples from the full posterior
  
  start_time <- Sys.time()
  if (length(dim(AlgoResults$Omega_sample_phys))==3){ #SMC
    AlgoResults %<>% addAlgoParams(s_finish)
  } else { #MCMC
    AlgoResults$s_finish = NA
  }
  
  if (single_particle){
    if(is.null(Omega)){
      if (particle_best){ # best particle at final step
        particle_i <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
      } else { # worst particle at final step (get a sense for how well we're doing at s_finish in the worst case)
        particle_i <- which(AlgoResults$d[,,AlgoResults$s_finish] == max(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
      }
      Omega <- AlgoResults$Omega_sample_phys[particle_i[1],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
    }
  } 
  df_poly_train <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat='Train', single_particle=single_particle, Omega = Omega, output=output) #T, particle=particle_min.d[1])
  df_poly_test <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat='Test', single_particle=single_particle, Omega = Omega, output=output) #T, particle=particle_min.d[1])
  #df_poly_jitter <- SampleImpact(dir, Model, Omega_min.d %>% addTransfParams(), AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M) %>% replace(which(names(AlgoParams)==c('Np')), 1))
    
  finish_time <- Sys.time()
  finish_time-start_time
  
  if (M > 1){
    df_poly <- rbind(df_poly_train, df_poly_test)
    if (length(which(is.na(df_poly$sampled.1))) > 0){
      df_poly <- df_poly[-which(is.na(df_poly$sampled.1)),]
    }
    return(df_poly)
    # jitter_val <- function(x){
    #   return_arr <- c()
    #   for (x_i in x){
    #     if(x_i==0){return_arr <- c(return_arr, runif(1,0.2, 0.3))}
    #     else {return_arr <- c(return_arr, runif(1, x_i-0.5, x+0.5))}
    #   }
    #   return(return_arr)
    # }
    # #df_poly_jitter[which(df_poly_jitter==0, arr.ind=T)] = 0.1
    # df_poly_jitter$observed <- sapply(df_poly_jitter$observed, jitter_val)
    # df_poly_jitter[,grep('sampled', names(df_poly_jitter))] <- apply(df_poly_jitter[,grep('sampled', names(df_poly_jitter))],1:2,jitter_val)
    #df_poly_jitter[,grep('sampled', names(df_poly_jitter))] <- t(apply(df_poly_jitter[,grep('sampled', names(df_poly_jitter))],1,sort))
    
    #return(df_poly_jitter)
    
    M_lower <- round(quantile(1:M, 0.05))
    M_lower <- ifelse(M_lower==0, 1, M_lower)
    M_median <- round(quantile(1:M, 0.5))
    M_upper <- round(quantile(1:M, 0.95))
    
    impact_type='displacement'
    ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=get(paste0('sampled.', M_median)), ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
      geom_errorbar() + geom_point(aes(col=train_flag)) + 
      scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
      scale_y_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
      #geom_pointrange(aes(col=train_flag)) + 
      geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
      ylab(paste('Sampled', impact_type)) + xlab(paste('Observed', impact_type)) + scale_color_manual(values = c('red', 'blue'))
    
    
    data_filt <- df_poly_test %>% filter(impact==impact_type)
    data_filt$sampled_median <- apply(data_filt[,grep('sampled', names(data_filt))],1,median)
    cor(data_filt$observed,data_filt$sampled_median)
    k <- 10
    cor(log(data_filt$observed+k),log(data_filt$sampled_median+k))
    
    SS_res_raw = sum((data_filt$observed-data_filt$sampled_median)^2)
    SS_tot_raw = sum((data_filt$observed-mean(data_filt$observed))^2)
    R_squared_raw = 1-SS_res_raw/SS_tot_raw
    R_squared_raw
    
    SS_res_log = sum((log(data_filt$observed+k)-log(data_filt$sampled_median+k))^2)
    SS_tot_log = sum((log(data_filt$observed+k)-mean(log(data_filt$observed+k)))^2)
    R_squared_log = 1-SS_res_log/SS_tot_log
    R_squared_log
    
    ggplot(df_poly_jitter %>% filter(impact=='buildDam'), mapping=aes(x=log(observed_jitter+1), y=log(get(paste0('sampled.', M_median))+1), ymin=log(get(paste0('sampled.', M_lower))+1), ymax=log(get(paste0('sampled.', M_upper))+1))) + 
      geom_errorbar(position = position_dodge(width = 0.1)) +
      geom_point(aes(col=inferred), position = position_dodge(width = 0.1)) + geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
      ylab('log(Sampled Impact+1)') + xlab('log(Observed Impact+1)')+ scale_color_manual(values = c('red', 'blue'))
    
    ggplot(df_poly_jitter %>% filter(impact=='displacement'), mapping=aes(x=log(observed_jitter+1), y=log(get(paste0('sampled.', M_median))+1), ymin=log(get(paste0('sampled.', M_lower))+1), ymax=log(get(paste0('sampled.', M_upper))+1))) + 
      geom_errorbar(position = position_dodge(width = 0.1)) +
      geom_point(aes(col=train_flag), position = position_dodge(width = 0.1)) + geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
      ylab('log(Sampled Impact+1)') + xlab('log(Observed Impact+1)')+ scale_color_manual(values = c('red', 'blue'))
  
    quants <- (apply(df_poly_jitter[,c(grep('observed', names(df_poly_jitter)),grep('sampled', names(df_poly_jitter)))], 1, sample_quant)-runif(NROW(df_poly_jitter),0,1))/(length(grep('sampled', names(df_poly_jitter)))+1)
    hist(quants[impact_type=='mortality'])
  }
  
  impact_type <- 'mortality'
  k <- 1
  ix_train <- which(df_poly_train$impact==impact_type)
  df_poly <- rbind(df_poly_train, df_poly_test)
  ix <- which(df_poly$impact==impact_type)
  cols = c(rep('black', length(ix_train)), rep('red', length(ix)-length(ix_train)))
  plot(log(df_poly$observed[ix]+k), log(df_poly$sampled.1[ix]+k), xlab=paste0('log(observed+',k,')'), ylab=paste0('log(sampled+',k,')'), col=cols)
  abline(0,1)
  
  #large_discrepancies
  df_poly[which(abs(log(df_poly$sampled.1+1)-log(df_poly$observed+1)) >7),] %>% filter(impact==impact_type)
}

create_df_postpredictive_MCMC <- function(AlgoResults, single_particle = F, Omega = NULL, particle_best=T, 
                                          M=50, output='SampledAgg'){
  # If single_particle == T:
  #       If Omega = NULL, finds the parameters that produced the smallest distance and samples using only these parameters
  #       If Omega is given, samples uses Omega
  # If single_particle == F, samples from the full posterior
  
  if (single_particle){
    if(is.null(Omega)){
      if (particle_best){ # best particle at final step
        particle_i <- which(AlgoResults$loss == min(AlgoResults$loss))
      } else { # worst particle at final step (get a sense for how well we're doing at s_finish in the worst case)
        stop()
        #particle_i <- which(AlgoResults$d[,,AlgoResults$s_finish] == max(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
      }
      Omega <- AlgoResults$Omega_sample_phys[,particle_i[1]] %>% relist(skeleton=Model$skeleton)
    }
  } 
  #get AlgoResults into ABC-SMC form
  AlgoResults$s_finish <- 1
  AlgoResults$W <- matrix(1/1000, nrow=1000, ncol=2)
  Omega_sample_phys_new <- array(0, dim=c(1000, nrow(AlgoResults$Omega_sample_phys), 2))
  stop_point <- which(is.na(AlgoResults$Omega_sample_phys[1,]))[1] -1
  Omega_sample_phys_new[,,1] <- t(AlgoResults$Omega_sample_phys[,(stop_point-999):stop_point])
  AlgoResults$Omega_sample_phys <- Omega_sample_phys_new
  AlgoResults$Npart <- 1000
  
  
  df_poly_train <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat='Train', single_particle=single_particle, Omega = Omega, output=output) #T, particle=particle_min.d[1])
  df_poly_test <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat='Test', single_particle=single_particle, Omega = Omega, output=output) #T, particle=particle_min.d[1])
  #df_poly_jitter <- SampleImpact(dir, Model, Omega_min.d %>% addTransfParams(), AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M) %>% replace(which(names(AlgoParams)==c('Np')), 1))
  
  # finish_time <- Sys.time()
  # finish_time-start_time
  
  if (M > 1){
    df_poly <- rbind(df_poly_train, df_poly_test)
    if (length(which(is.na(df_poly$sampled.1))) > 0){
      df_poly <- df_poly[-which(is.na(df_poly$sampled.1)),]
    }
    return(df_poly)
    # jitter_val <- function(x){
    #   return_arr <- c()
    #   for (x_i in x){
    #     if(x_i==0){return_arr <- c(return_arr, runif(1,0.2, 0.3))}
    #     else {return_arr <- c(return_arr, runif(1, x_i-0.5, x+0.5))}
    #   }
    #   return(return_arr)
    # }
    # #df_poly_jitter[which(df_poly_jitter==0, arr.ind=T)] = 0.1
    # df_poly_jitter$observed <- sapply(df_poly_jitter$observed, jitter_val)
    # df_poly_jitter[,grep('sampled', names(df_poly_jitter))] <- apply(df_poly_jitter[,grep('sampled', names(df_poly_jitter))],1:2,jitter_val)
    #df_poly_jitter[,grep('sampled', names(df_poly_jitter))] <- t(apply(df_poly_jitter[,grep('sampled', names(df_poly_jitter))],1,sort))
    
    #return(df_poly_jitter)
    
    M_lower <- round(quantile(1:M, 0.05))
    M_lower <- ifelse(M_lower==0, 1, M_lower)
    M_median <- round(quantile(1:M, 0.5))
    M_upper <- round(quantile(1:M, 0.95))
    
    impact_type='displacement'
    ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=get(paste0('sampled.', M_median)), ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
      geom_errorbar() + geom_point(aes(col=train_flag)) + 
      scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
      scale_y_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
      #geom_pointrange(aes(col=train_flag)) + 
      geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
      ylab(paste('Sampled', impact_type)) + xlab(paste('Observed', impact_type)) + scale_color_manual(values = c('red', 'blue'))
    
    
    data_filt <- df_poly_test %>% filter(impact==impact_type)
    data_filt$sampled_median <- apply(data_filt[,grep('sampled', names(data_filt))],1,median)
    cor(data_filt$observed,data_filt$sampled_median)
    k <- 10
    cor(log(data_filt$observed+k),log(data_filt$sampled_median+k))
    
    SS_res_raw = sum((data_filt$observed-data_filt$sampled_median)^2)
    SS_tot_raw = sum((data_filt$observed-mean(data_filt$observed))^2)
    R_squared_raw = 1-SS_res_raw/SS_tot_raw
    R_squared_raw
    
    SS_res_log = sum((log(data_filt$observed+k)-log(data_filt$sampled_median+k))^2)
    SS_tot_log = sum((log(data_filt$observed+k)-mean(log(data_filt$observed+k)))^2)
    R_squared_log = 1-SS_res_log/SS_tot_log
    R_squared_log
    
    ggplot(df_poly_jitter %>% filter(impact=='buildDam'), mapping=aes(x=log(observed_jitter+1), y=log(get(paste0('sampled.', M_median))+1), ymin=log(get(paste0('sampled.', M_lower))+1), ymax=log(get(paste0('sampled.', M_upper))+1))) + 
      geom_errorbar(position = position_dodge(width = 0.1)) +
      geom_point(aes(col=inferred), position = position_dodge(width = 0.1)) + geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
      ylab('log(Sampled Impact+1)') + xlab('log(Observed Impact+1)')+ scale_color_manual(values = c('red', 'blue'))
    
    ggplot(df_poly_jitter %>% filter(impact=='displacement'), mapping=aes(x=log(observed_jitter+1), y=log(get(paste0('sampled.', M_median))+1), ymin=log(get(paste0('sampled.', M_lower))+1), ymax=log(get(paste0('sampled.', M_upper))+1))) + 
      geom_errorbar(position = position_dodge(width = 0.1)) +
      geom_point(aes(col=train_flag), position = position_dodge(width = 0.1)) + geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
      ylab('log(Sampled Impact+1)') + xlab('log(Observed Impact+1)')+ scale_color_manual(values = c('red', 'blue'))
    
    quants <- (apply(df_poly_jitter[,c(grep('observed', names(df_poly_jitter)),grep('sampled', names(df_poly_jitter)))], 1, sample_quant)-runif(NROW(df_poly_jitter),0,1))/(length(grep('sampled', names(df_poly_jitter)))+1)
    hist(quants[impact_type=='mortality'])
  }
  
  impact_type <- 'mortality'
  k <- 1
  ix_train <- which(df_poly_train$impact==impact_type)
  df_poly <- rbind(df_poly_train, df_poly_test)
  ix <- which(df_poly$impact==impact_type)
  cols = c(rep('black', length(ix_train)), rep('red', length(ix)-length(ix_train)))
  plot(log(df_poly$observed[ix]+k), log(df_poly$sampled.1[ix]+k), xlab=paste0('log(observed+',k,')'), ylab=paste0('log(sampled+',k,')'), col=cols)
  abline(0,1)
  
  #large_discrepancies
  df_poly[which(abs(log(df_poly$sampled.1+1)-log(df_poly$observed+1)) >7),] %>% filter(impact==impact_type)
}

flattenImpactSample <- function(impact_sample){
  df <- data.frame(event_id = impact_sample$poly[[1]]$event_id,
                   iso3 = impact_sample$poly[[1]]$iso3,
                   polygon = impact_sample$poly[[1]]$polygon,
                   impact = impact_sample$poly[[1]]$impact,
                   observed = impact_sample$poly[[1]]$observed, 
                   train_flag = impact_sample$poly[[1]]$train_flag)
  for (j in 1:length(impact_sample$poly)){
    df %<>% cbind(impact_sample$poly[[j]]$sampled)
  }
  colnames(df)[7:NCOL(df)] <- paste0('sampled.', 1:length(impact_sample$poly))
  return(df)
}


create_df_postpredictive_from_impact_sample = function(impact_sample){
  
  impact_sample_flattened <- flattenImpactSample(impact_sample)
  
  df_poly_jitter <- impact_sample_flattened
  if (any(is.na(df_poly_jitter$sampled.1))){
    df_poly_jitter <- df_poly_jitter[-which(is.na(df_poly_jitter$sampled.1)),]
  }
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
  #df_poly_jitter[,grep('sampled', names(df_poly_jitter))] <- t(apply(df_poly_jitter[,grep('sampled', names(df_poly_jitter))],1,sort))
  return(df_poly_jitter)
}

plot_df_postpredictive <- function(df_poly_jitter, impact_type){
  
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
  SS_tot_raw = sum((df_poly_filt$observed-mean(df_poly_filt$observed))^2)
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
    geom_errorbar() + geom_point(aes(col=train_flag)) + 
    scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    scale_y_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    #geom_pointrange(aes(col=train_flag)) + 
    geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
    ylab(paste('Sampled', ifelse(impact_type=='buildDam', 'building damage', impact_type))) + xlab(paste('Observed', ifelse(impact_type=='buildDam', 'building damage', impact_type))) + scale_color_manual(values = c('red', 'blue')) + 
    theme_bw() +
    theme(axis.title = element_text(family = "Liberation Serif", size=12),  
                legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
                legend.title = element_text(family = "Liberation Serif", size=12))
}

plot_df_postpredictive_compare <- function(df_poly1, df_poly2, impact_type){
  df_poly1[,grep('sampled', names(df_poly1))] <- t(apply(df_poly1[,grep('sampled', names(df_poly1))],1,sort))
  df_poly2[,grep('sampled', names(df_poly2))] <- t(apply(df_poly2[,grep('sampled', names(df_poly2))],1,sort))
  
  M <- length(grep('sampled', names(df_poly1)))
  M_lower <- round(quantile(1:M, 0.05))
  M_lower <- ifelse(M_lower==0, 1, M_lower)
  M_median <- round(quantile(1:M, 0.5))
  M_upper <- round(quantile(1:M, 0.95))
  df_poly1$means_sampled <- apply(df_poly1[,grep("sampled",names(df_poly1),value = T)], 1, median)
  df_poly2$means_sampled <- apply(df_poly2[,grep("sampled",names(df_poly2),value = T)], 1, median)
  
  
  df_poly1_jitter = df_poly1 %>% filter(impact==impact_type)
  df_poly2_jitter = df_poly2 %>% filter(impact==impact_type)
  
  #jitter
  # for (repeat_jitter in 1:5){
  #   for (i in 1:NROW(df_poly1_jitter)){
  #     for (j in 1:NROW(df_poly1_jitter)){
  #       if (i == j) next
  #       # if ((df_poly1_jitter$observed[i] <= 2) | (df_poly2_jitter$observed[j] <= 2)){
  #       #   if(abs(df_poly1_jitter$observed[i] - df_poly2_jitter$observed[j])<2){
  #       #     stop()
  #       #     k <- c(i,j)[which.max(c(i,j))][1]
  #       #     print(paste(df_poly2_jitter$observed[k],  df_poly1_jitter$observed[k] + df_poly2_jitter$observed[k]/5))
  #       #     df_poly1_jitter$observed[k] = df_poly1_jitter$observed[k] + df_poly2_jitter$observed[k]/5
  #       #     df_poly2_jitter$observed[k] = df_poly2_jitter$observed[k] + df_poly2_jitter$observed[k]/5
  #       #   }
  #       #   next
  #       # }
  #       if(abs(df_poly1_jitter$observed[i]/20 - df_poly2_jitter$observed[j]/20)<2){
  #         k <- c(i,j)[which.max(c(i,j))][1]
  #         #print(paste(df_poly2_jitter$observed[k],  df_poly1_jitter$observed[k] + df_poly2_jitter$observed[k]/5))
  #         df_poly1_jitter$observed[k] = df_poly1_jitter$observed[k] + df_poly2_jitter$observed[k]/10
  #         df_poly2_jitter$observed[k] = df_poly2_jitter$observed[k] + df_poly2_jitter$observed[k]/10
  #       }
  #     }
  #   }
  # }
  # df_poly1_jitter %<>% filter(observed > 1)
  # df_poly2_jitter %<>% filter(observed > 1)
  
  #df_poly1_jitter <- df_poly1_jitter[order(df_poly1_jitter$observed),][1:(NROW(df_poly1_jitter)/2)*2,]
  #df_poly2_jitter <- df_poly2_jitter[order(df_poly2_jitter$observed),][1:(NROW(df_poly2_jitter)/2)*2,]

  ggplot() +
    #ggplot(df_poly_jitter %>% filter(impact==impact_type), mapping=aes(x=observed, y=means_sampled, ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
    #geom_errorbar(aes(x=observed, y=get(paste0('sampled.', M_median)), ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + geom_point(aes(x=observed, y=get(paste0('sampled.', M_median)),col=train_flag)) + 
    geom_errorbar(mapping=aes(x=df_poly2_jitter$observed-df_poly2_jitter$observed/30, y=df_poly2_jitter[,paste0('sampled.', M_median)], ymin=df_poly2_jitter[,paste0('sampled.', M_lower)], ymax=df_poly2_jitter[,paste0('sampled.', M_upper)], col='Posterior Predictive'), lwd=1) +
    geom_errorbar(mapping=aes(x=df_poly1_jitter$observed, y=df_poly1_jitter[,paste0('sampled.', M_median)], ymin=df_poly1_jitter[,paste0('sampled.', M_lower)], ymax=df_poly1_jitter[,paste0('sampled.', M_upper)], col='True Predictive'), lwd=0.75) +
    geom_point(aes(x=df_poly2_jitter$observed-df_poly2_jitter$observed/30, y=df_poly2_jitter[,paste0('sampled.',M_median)],col='Posterior Predictive'), cex=2) + 
    geom_point(aes(x=df_poly1_jitter$observed, y=df_poly1_jitter[,paste0('sampled.',M_median)],col='True Predictive'), cex=2) + 
    scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    scale_y_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    #geom_pointrange(aes(col=train_flag)) + 
    geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) + 
    ylab(paste('Sampled', impact_type)) + xlab(paste('Observed', impact_type)) + scale_color_manual(values = c('black', 'red')) + theme_bw()  + labs(color='Median and 90% Credible Interval') +
    theme(axis.title = element_text(family = "Liberation Serif", size=12),  
          legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
          legend.title = element_text(family = "Liberation Serif", size=12))
  
}

plot_df_postpredictive_PAGER_coloured <- function(df_poly_jitter, impact_type){
  
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
  df_poly_jitter$means_sampled <- apply(df_poly_jitter[,grep("sampled",names(df_poly_jitter),value = T)], 1, mean)
  
  folderin_haz <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/HAZARDobjects_wMaxMMIDiff/'
  ufiles_haz <- na.omit(list.files(path=folderin_haz,pattern=Model$haz,recursive = T,ignore.case = T))
  df_poly_jitter$alertlevel <- ''
  for(i in 1:NROW(df_poly_jitter)){
    file_match <- grep(paste0("_", df_poly_jitter$event_id[i], "\\b"),  ufiles_haz, value = TRUE)
    HAZy <- readRDS(paste0(folderin_haz, file_match ))
    which.max.mmi <- which.max(sapply(HAZy[2:length(HAZy)], function(x) max(x$mean)))
    df_poly_jitter$alertlevel[i] <- HAZy$hazard_info$alertlevel[which.max.mmi]
  }
  
  df_poly_jitter$alertlevel <- factor(df_poly_jitter$alertlevel, levels=c('green', 'yellow', 'orange', 'red'))
  ggplot(df_poly_jitter %>% filter(impact==impact_type & train_flag=='TEST' & alertlevel !='null'), mapping=aes(x=observed, y=get(paste0('sampled.', M_median)), ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
  #ggplot(df_poly_jitter %>% filter(impact==impact_type & train_flag=='TEST' & alertlevel !='null'), mapping=aes(x=observed, y=means_sampled, ymin=get(paste0('sampled.', M_lower)), ymax=get(paste0('sampled.', M_upper)))) + 
    geom_errorbar() + 
    geom_abline(slope=1, intercept=0) + theme(aspect.ratio=1) +
    geom_point(col='red', size=2.5) + #geom_point(aes(col=alertlevel), size=2.5) + 
    geom_point(shape = 1, size = 3,colour = "black") + 
    #scale_x_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) log(x+10), labels = scales::trans_format("log10")), labels = scales::comma) + 
    scale_x_continuous(trans=scales::pseudo_log_trans(sigma=5,base = 10), breaks = c(0,20,100,1000, 10000, 100000), labels= scales::comma_format(),  minor_breaks =NULL) + 
    scale_y_continuous(trans=scales::pseudo_log_trans(sigma=5,base = 10), breaks = c(0,20,100,1000, 10000, 100000), labels= scales::comma_format(),  minor_breaks =NULL) + 
    #scale_y_continuous(trans='log10', breaks = scales::trans_breaks("log10", function(x) 10^x, labels = scales::trans_format("log10")), labels = scales::comma) + 
    #geom_pointrange(aes(col=train_flag)) + 
    ylab(paste('ODDRIN Posterior Predictive', impact_type)) + xlab(paste('Observed', impact_type)) + 
    #scale_color_manual(values = c('green'='green', 'yellow'='yellow', 'orange'='orange', 'red'='red'), labels = c('green'='0', 'yellow'='1 - 99', 'orange'='100 - 999', 'red'='1000+'),name = "PAGER Predicted Mortality") + 
    theme_bw() + theme(plot.background = element_rect(fill = rgb(0.95,0.95,0.95)))
}

gdacs_df <- rbind(
  c(31, 'red'),
  c(70, 'red'),
  c(7, 'green'),
  c(170, 'red'),
  c(54, 'green'),
  c(39, 'orange'),
  c(140, 'orange'),
  c(92, 'green'),
  c(155, 'green'),
  c(164, 'red'),
  c(22, 'green'),
  c(36, 'red'),
  c(16, 'red'),
  c(109, 'orange'),
  c(94, 'green'),
  c(167, 'orange'),
  c(112, 'green'),
  c(149, 'orange'),
  c(124, 'orange'),
  c(118, 'green'),
  c(131, 'green'),
  c(51, 'orange'),
  c(63, 'green'),
  c(25, 'green'),
  c(106, 'orange'),
  c(48, 'green'),
  c(152, 'orange'),
  c(136, 'orange'),
  c(143, 'orange'),
  c(100, 'green'),
  c(79, 'green'),
  c(88, 'green'),
  c(57, 'green'),
  c(121, 'orange'),
  c(115, 'green'),
  c(128, 'null'),
  c(137, 'orange'),
  c(28, 'green'),
  c(4, 'green'),
  c(82, 'orange'),
  c(97, 'green'),
  c(161, 'green'),
  c(13, 'green'),
  c(13, 'green'),
  c(42, 'green'),
  c(158, 'green'),
  c(146, 'green'),
  c(85, 'green'),
  c(76, 'green'),
  c(60, 'green'),
  c(45, 'green'),
  c(103, 'green'),
  c(10, 'green'))
gdacs_df = data.frame(gdacs_df)
colnames(gdacs_df) <- c('event_id', 'gdacs')

gdacs_df <- rbind(
  c(169, 'red', 'TUR', '2023-02-06'),
  c(7, 'green', 'IRN', '2013-04-16'),
  c(67, 'orange', 'MEX', '2017-09-07'),
  c(54, 'green', 'CHL', '2016-12-25'),
  c(81, 'red', 'PNG', '2018-02-26'),
  c(14, 'red', 'PAK', '2013-09-24'),
  c(145, 'green', 'MEX', '2020-06-23'),
  c(107, 'orange', 'PNG', '2019-05-07'),
  c(150, 'orange', 'TUR', '2020-10-30'),
  c(133, 'orange', 'ALB', '2019-11-26'),
  c(16, 'red', 'PHL', '2013-10-15'),
  c(64, 'orange', 'GTM', '2017-06-14'),
  c(78, 'green', 'PER', '2018-01-14'),
  c(38, 'green', 'JPN', '2016-04-14'),
  c(15, 'green', 'PER', '2013-09-25'),
  c(40, 'orange', 'ECU', '2016-05-18'),
  c(23, 'null', 'GTM', '2014-07-07'),
  c(26, 'green', 'IRN', '2014-08-18'),
  c(130, 'orange', 'PHL', '2019-11-18'),
  c(118, 'green', 'IDN', '2019-07-14'),
  c(3, 'null', 'UZB', '2011-07-20'),
  c(63, 'green', 'GRC', '2017-06-12'),
  c(62, 'orange', 'IDN', '2017-05-29'),
  c(165, 'red', 'AFG', '2022-06-22'),
  c(56, 'orange', 'PHL', '2017-02-10'),
  c(37, 'null', 'TWN', '2016-02-06'),
  c(154, 'orange', 'PHL', '2021-02-07'),
  c(47, 'green', 'JPN', '2016-10-21'),
  c(12, 'green', 'CHN', '2013-07-22'),
  c(100, 'green', 'CHN', '2018-12-16'),
  c(71, 'green', 'CRI', '2017-11-12'),
  c(79, 'green', 'TWN', '2018-02-04'),
  c(55, 'green', 'BGD', '2017-01-03'),
  c(147, 'green', 'USA', '2020-08-09'),
  c(77, 'green', 'MMR', '2018-01-12'),
  c(57, 'green', 'ZMB', '2017-02-24'),
  c(58, 'orange', 'PHL', '2017-04-08'),
  c(119, 'green', 'TUR', '2019-08-08'),
  c(115, 'green', 'IRN', '2019-07-08'),
  c(18, 'green', 'GRC', '2014-01-26'),
  c(117, 'green', 'PHL', '2019-07-13'),
  c(4, 'green', 'BGR', '2012-05-22'),
  c(132, 'green', 'CHN', '2019-11-25'),
  c(82, 'orange', 'MOZ', '2018-03-08'),
  c(72, 'green', 'KOR', '2017-11-15'),
  c(146, 'green', 'DZA', '2020-08-07'),
  c(110, 'green', 'ALB', '2019-06-01'),
  c(60, 'green', 'CHN', '2017-05-11'),
  c(86, 'green', 'SLV', '2018-05-06'),
  c(103, 'green', 'AZE', '2019-02-05'),
  c(102, 'green', 'ITA', '2018-12-26'),
  c(30, 'green', 'CYP', '2015-04-15'),
  c(84, 'green', 'IDN', '2018-04-18')
)
gdacs_df = data.frame(gdacs_df)
colnames(gdacs_df) <- c('event_id', 'gdacs', 'iso3c', 'date')

pager_compare <- function(df_poly_jitter, impact_type){
  
  df_poly_jitter <- df_postpredictive_sampled_best %>% filter(impact=='mortality' & train_flag=='TEST')
  folderin_haz <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Aug24/HAZARDobjects/'
  ufiles_haz <- na.omit(list.files(path=folderin_haz,pattern=Model$haz,recursive = T,ignore.case = T))
  
  df_poly_jitter <- add_landslide_flag(df_poly_jitter)
  df_poly_jitter$date <- as.Date('2010-01-01')
  for (i in 1:NROW(df_poly_jitter)){
    file_match <- grep(paste0("_", df_poly_jitter$event_id[i], "\\b"),  ufiles_haz, value = TRUE)
    HAZy <- readRDS(paste0(folderin_haz, file_match ))
    which.max.mmi <- which.max(sapply(HAZy[2:length(HAZy)], function(x) max(values(unwrap(x$spatrast)$mean), na.rm=T)))
    df_poly_jitter$date[i] <- HAZy$hazard_info$eventdates[which.max.mmi]
  }
  
  df_poly_jitter$Landslide.Fatalities <- ifelse(is.na(df_poly_jitter$Landslide.Fatalities), 0, df_poly_jitter$Landslide.Fatalities)
  df_poly_jitter$observed_LS_adj <- df_poly_jitter$observed - df_poly_jitter$Landslide.Fatalities
  df_poly_jitter$alertlevel <- ''
  
  df_poly_jitter$obs_prob_pager <- NA
  df_poly_jitter$obs_prob_oddrin <- NA
  df_poly_jitter$oddrin_rps <- 0
  df_poly_jitter$pager_rps <- 0
  df_poly_jitter$oddrin_rps <- 0
  df_poly_jitter$pager_brier <- 0
  df_poly_jitter$oddrin_brier <- 0
  for(i in 1:NROW(df_poly_jitter)){
    file_match <- grep(paste0("_", df_poly_jitter$event_id[i], "\\b"),  ufiles_haz, value = TRUE)
    HAZy <- readRDS(paste0(folderin_haz, file_match ))
    which.max.mmi <- which.max(sapply(HAZy[2:length(HAZy)], function(x) max(values(unwrap(x$spatrast)$mean), na.rm=T)))
    df_poly_jitter$alertlevel[i] <- HAZy$hazard_info$alertlevel[which.max.mmi]
    fatality_alert <- HAZy$hazard_info$alertfull[[which.max.mmi]]$alert_level
    if (is.null(fatality_alert)){
      df_poly_jitter$alertlevel[i] <- 'null'
      next
    }
    df_poly_jitter$alertlevel[i] <- fatality_alert
    
    if (length(HAZy$hazard_info$alertfull[[which.max.mmi]]$bins) == 0) next
    prob_mult <- 1/sum(unlist(lapply(HAZy$hazard_info$alertfull[[which.max.mmi]]$bins, function(x) return(x$probability)))) # some sum to 1 and some to 100
    cum_prob_oddrin <- 0
    cum_prob_pager <- 0
    cum_obs_oddrin <- 0
    cum_obs_pager <- 0
    if (length(HAZy$hazard_info$alertfull[[which.max.mmi]]$bins) != 7) warning(paste('Check alert bins for event id:',df_poly_jitter$event_id[i]))
    for (j in 1:7){
      oddrin_preds <- df_poly_jitter[i,grep('sampled.', names(df_poly_jitter))]
      oddrin_prob = mean((oddrin_preds >= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$min) & (oddrin_preds <= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$max))
      df_poly_jitter$oddrin_brier[i] = df_poly_jitter$oddrin_brier[i] + (oddrin_prob - ifelse(df_poly_jitter$observed[i] >= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$min & df_poly_jitter$observed[i] <= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$max, 1, 0))^2
      df_poly_jitter$pager_brier[i] = df_poly_jitter$pager_brier[i] + (HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$probability * prob_mult - ifelse(df_poly_jitter$observed_LS_adj[i] >= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$min & df_poly_jitter$observed_LS_adj[i] <= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$max, 1, 0))^2
      #df_poly_jitter$oddrin_rps[i] = df_poly_jitter$oddrin_rps[i] + (oddrin_prob - ifelse(df_poly_jitter$observed[i] < HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$max, 1, 0))^2
      #df_poly_jitter$pager_rps[i] = df_poly_jitter$pager_rps[i] + (HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$probability * prob_mult - ifelse(df_poly_jitter$observed_LS_adj[i]  <= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$max, 1, 0))^2
      cum_prob_oddrin <- cum_prob_oddrin + oddrin_prob
      cum_prob_pager <- cum_prob_pager + HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$probability * prob_mult
      cum_obs_oddrin <- ifelse(df_poly_jitter$observed[i] >= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$min & df_poly_jitter$observed[i] <= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$max, 1, cum_obs_oddrin)
      cum_obs_pager <- ifelse(df_poly_jitter$observed[i] >= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$min & df_poly_jitter$observed[i] <= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$max, 1, cum_obs_pager)
      df_poly_jitter$oddrin_rps[i] = df_poly_jitter$oddrin_rps[i] + (cum_prob_oddrin - cum_obs_oddrin)^2
      df_poly_jitter$pager_rps[i] = df_poly_jitter$pager_rps[i] + (cum_prob_pager - cum_obs_pager)^2
      
      if ((df_poly_jitter$observed[i] >= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$min) & (df_poly_jitter$observed[i] <= HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$max)){
        
        df_poly_jitter$obs_prob_pager[i] = -log2(HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$probability * prob_mult)#HAZy$hazard_info$alertfull[[which.max.mmi]]$bins[[j]]$probability * prob_mult
        if(!is.finite(df_poly_jitter$obs_prob_pager[i])){ df_poly_jitter$obs_prob_pager[i] = -log2(0.01)}
        df_poly_jitter$obs_prob_oddrin[i] = -log2(oddrin_prob)
      }
    }
  }
  mean(df_poly_jitter$oddrin_rps)
  mean(df_poly_jitter$pager_rps)
  
  plot(df_poly_jitter$oddrin_rps, ylim=c(0,3))
  points(df_poly_jitter$pager_rps, col='red')
  
  plot(df_poly_jitter$oddrin_rps)
  points(df_poly_jitter$pager_rps, col='red')
  
  mean(df_poly_jitter$oddrin_rps)
  mean(df_poly_jitter$pager_rps)
  sum(df_poly_jitter$oddrin_rps[order( df_poly_jitter$oddrin_rps - df_poly_jitter$pager_rps, decreasing=T)][1:51])
  sum(df_poly_jitter$pager_rps[order(df_poly_jitter$oddrin_rps - df_poly_jitter$pager_rps , decreasing=T)][1:51])
  
  mean(df_poly_jitter$oddrin_brier)
  mean(df_poly_jitter$pager_brier)
  
  mean(df_poly_jitter$oddrin_rps)
  mean(df_poly_jitter$pager_rps)
  sort(df_poly_jitter$oddrin_rps-df_poly_jitter$pager_rps)
  
  df_poly_jitter %<>% filter(impact==impact_type)
  
  
  df_poly_jitter$median_sampled <- apply(df_poly_jitter[,grep("sampled",names(df_poly_jitter),value = T)], 1, median)
  df_poly_jitter$oddrin_alert <- ifelse(df_poly_jitter$median_sampled>=1000, 'red', 
                                        ifelse(df_poly_jitter$median_sampled>=100, 'orange', 
                                               ifelse(df_poly_jitter$median_sampled>=1, 'yellow', 'green')))
  df_poly_jitter$real_alert <- ifelse(df_poly_jitter$observed>=1000, 'red', 
                                        ifelse(df_poly_jitter$observed>=100, 'orange', 
                                               ifelse(df_poly_jitter$observed>=1, 'yellow', 'green')))
  
  sum(df_poly_jitter$alertlevel!= 'null')
  mean(df_poly_jitter$oddrin_alert[df_poly_jitter$alertlevel!= 'null']==df_poly_jitter$real_alert[df_poly_jitter$alertlevel!= 'null'])
  mean(df_poly_jitter$alertlevel[df_poly_jitter$alertlevel!= 'null']==df_poly_jitter$real_alert[df_poly_jitter$alertlevel!= 'null'])
  
  df_poly_jitter[df_poly_jitter$alertlevel!= 'null',c('event_id','real_alert','observed', 'alertlevel','oddrin_alert', 'landslide_flag', 'Landslide.Fatalities', 'oddrin_rps', 'pager_rps')]
  
  df_poly_jitter <- merge(df_poly_jitter, gdacs_df, by='event_id')
  df_poly_jitter$real_alert_gdacs <- ifelse(df_poly_jitter$observed>=100, 'red', 
                                             ifelse(df_poly_jitter$observed>=10, 'orange', 'green'))
  df_poly_jitter$alert_oddrin_gdacs <- ifelse(df_poly_jitter$median_sampled>=100, 'red', 
                                            ifelse(df_poly_jitter$median_sampled>=10, 'orange', 'green'))
  
  mean(df_poly_jitter$gdacs[df_poly_jitter$gdacs!= 'null']== df_poly_jitter$real_alert_gdacs[df_poly_jitter$gdacs!= 'null'])
  mean(df_poly_jitter$alert_oddrin_gdacs[df_poly_jitter$gdacs!= 'null']== df_poly_jitter$real_alert_gdacs[df_poly_jitter$gdacs!= 'null'])
  
  plot(log(df_poly_jitter$observed+10), log(df_poly_jitter$means_sampled+10), col=df_poly_jitter$oddrin_alert)
  plot(log(df_poly_jitter$observed+10), log(df_poly_jitter$means_sampled+10), col=ifelse(df_poly_jitter$alertlevel=='null', 'black', df_poly_jitter$alertlevel))
  plot(log(df_poly_jitter$observed+10), log(df_poly_jitter$means_sampled+10), col=df_poly_jitter$real_alert)
  
}

check_quants <- function(AlgoResults, s_finish=NULL){
  #check quants over multiple simulations to check for consistency
  sim_n <- 5
  M <- 20
  AlgoResults %<>% addAlgoParams(s_finish)
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  Omega_min.d <- AlgoResults$Omega_sample_phys[particle_min.d[1],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  quants_store <- NULL
  for (i in 1:sim_n){
    print(i)
    impact_sample <- SampleImpact(dir, Model, Omega_min.d  %>% addTransfParams(), AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M) %>% replace(which(names(AlgoParams)==c('Np')), 1))
    observed <- impact_sample$poly[[1]]$observed
    impact_type <- impact_sample$poly[[1]]$impact
    samples_allocated <- 1:length(impact_sample$poly)
    samples_combined <- sapply(impact_sample$poly[samples_allocated], function(x){x$sampled}) #doesn't work if samples_allocated is length 1
    medians <- apply(samples_combined, 1, median)
    quants <- (apply(cbind(observed,samples_combined), 1, sample_quant)-runif(length(observed),0,1))/(NCOL(samples_combined)+1)
    quants_store <- cbind(quants_store, quants)
  }
  par(mfrow=c(round(sqrt(sim_n)), ceiling(sim_n / round(sqrt(sim_n)))))
  for (i in 1:sim_n){
    hist(quants_store[impact_type=='mortality',i], breaks=10, main='', xlab='')
  }
  par(mfrow=c(1,1))
  
  observed <- impact_sample$poly[[1]]$observed
  impact_type <- impact_sample$poly[[1]]$impact
  samples_allocated <- 1:length(impact_sample$poly)
  samples_combined <- sapply(impact_sample$poly[samples_allocated], function(x){x$sampled}) #doesn't work if samples_allocated is length 1
  medians <- apply(samples_combined, 1, median)
  sum((log(medians[which(impact_type=='mortality')]+10)-log(observed[which(impact_type=='mortality')]+10))^2 * unlist(AlgoParams$impact_weights['mortality']))
  sum((log(medians[which(impact_type=='displacement')]+10)-log(observed[which(impact_type=='displacement')]+10))^2 * unlist(AlgoParams$impact_weights['displacement']))
  sum((log(medians[which(impact_type=='buildDam')]+10)-log(observed[which(impact_type=='buildDam')]+10))^2 * unlist(AlgoParams$impact_weights['buildDam']))
  
  hist((apply(cbind(observed,samples_combined), 1, sample_quant))[impact_type=='mortality' & medians!=0])
  
  quants <- (apply(cbind(observed,samples_combined), 1, sample_quant)-runif(length(observed),0,1))/(NCOL(samples_combined)+1)
  hist(quants[impact_type=='mortality'])
  hist(quants[impact_type=='mortality' & medians != 0])
  hist(quants[impact_type=='displacement'])
  hist(quants[impact_type=='buildDam' & !impact_sample$poly[[1]]$inferred])
  AndersonDarlingTest(quants[impact_type=='mortality'], 'punif')$statistic * unlist(AlgoParams$impact_weights['mortality'])
  AndersonDarlingTest(quants[impact_type=='displacement'], null='punif')$statistic * unlist(AlgoParams$impact_weights['displacement'])
  AndersonDarlingTest(quants[impact_type=='buildDam'], null='punif')$statistic * unlist(AlgoParams$impact_weights['buildDam'])

  finish_time <- Sys.time()
  finish_time-start_time
  
  #compare to different parameters
  Omega_min.d <- AlgoResults$Omega_sample_phys[particle_min.d[1],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  Omega_min.d$eps$hazard_mort <- 0.8
  impact_sample2 <- SampleImpact(dir, Model, Omega_min.d  %>% addTransfParams(), AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M) %>% replace(which(names(AlgoParams)==c('Np')), 1))
  
  observed2 <- impact_sample2$poly[[1]]$observed
  impact_type2 <- impact_sample2$poly[[1]]$impact
  samples_allocated2 <- 1:length(impact_sample2$poly)
  samples_combined2 <- sapply(impact_sample2$poly[samples_allocated2], function(x){x$sampled}) #doesn't work if samples_allocated is length 1
  medians2 <- apply(samples_combined2, 1, median)
  sum((log(medians2[which(impact_type2=='mortality')]+10)-log(observed2[which(impact_type2=='mortality')]+10))^2 * unlist(AlgoParams$impact_weights['mortality']))

  quants2 <- (apply(cbind(observed2,samples_combined2), 1, sample_quant)-runif(length(observed2),0,1))/(NCOL(samples_combined2)+1)
  hist(quants[impact_type=='mortality'])

  quants <- (apply(cbind(observed,samples_combined), 1, sample_quant)-runif(length(observed),0,1))/(NCOL(samples_combined)+1)
  quants2 <- (apply(cbind(observed2,samples_combined2), 1, sample_quant)-runif(length(observed2),0,1))/(NCOL(samples_combined2)+1)
  hist(quants2)
  cvm.test(quants[impact_type=='mortality'], 'punif')$statistic
  cvm.test(quants2[impact_type2=='mortality'], 'punif')$statistic
  
  10*AndersonDarlingTest(quants[impact_type=='displacement'], null='punif')$statistic * unlist(AlgoParams$impact_weights['displacement'])
  10*AndersonDarlingTest(quants[impact_type=='buildDam'], null='punif')$statistic * unlist(AlgoParams$impact_weights['buildDam'])
  
}

plot_cor_posts_poster <- function(AlgoResults, pars=c(16,17)){
  AlgoResults %<>% addAlgoParams()
  p1 <- ggplot() + 
        geom_point(aes(x=AlgoResults$Omega_sample_phys[1:450,pars[1],1], y=AlgoResults$Omega_sample_phys[1:450,pars[2],1]), col='blue') +
        geom_point(aes(x=AlgoResults$Omega_sample_phys[,pars[1],AlgoResults$s_finish], y=AlgoResults$Omega_sample_phys[,pars[2],AlgoResults$s_finish])) +
        theme_bw() + theme(plot.background = element_rect(fill = rgb(0.95,0.95,0.95))) + 
        xlab(names(unlist(Omega))[pars[1]]) + ylab(names(unlist(Omega))[pars[2]])
  return(p1)
}

check_outliers <- function(s_finish=NULL){
  start_time <- Sys.time()
  AlgoResults %<>% addAlgoParams(s_finish)
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  Omega_min.d <- AlgoResults$Omega_sample_phys[particle_min.d[1],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  impact_sample <- SampleImpact(dir, Model, Omega_min.d  %>% addTransfParams(), AlgoParams)
  
  impact_sample$poly[[1]] %>% filter(impact=='buildDam' & log(sampled+1)>(log(observed+1)+3))
}


add_covar <- function(df_poly, covars='EQFreq', dat='all'){

  folderin<-paste0(dir,AlgoResults$input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  for (covar in covars){
    df_poly[[covar]] <- NA
  }
  
  for(i in 1:NROW(df_poly)){
    file_match <- grep(paste0("_", df_poly$event_id[i], "\\b"),  ufiles, value = TRUE)
    if (length(file_match) > 1) stop('Multiple ODD files for event with the same ID.')
    ODDy <- readRDS(paste0(folderin,file_match))
    poly_indexes <- ODDy@polygons[[df_poly$polygon[i]]]$indexes
    if ('hazMean' %in% covars){
      poly_pop_restricted <- poly_indexes[which(ODDy@data$Population[poly_indexes] > 1000)]
      df_poly[['hazMean']][i] <- max(ODDy@data[poly_pop_restricted,grep("hazMean",names(ODDy),value = T)], na.rm=T)
    } 
    #df_poly[[covar]][i] <- median(ODDy@data[[covar]][poly_indexes])
    #find mode:
    for (covar in covars[which(covars != 'hazMean')]){
      tbl <- table(ODDy@data[[covar]][poly_indexes])
      if (length(tbl)==0){df_poly[[covar]][i] <- NA; next }
      df_poly[[covar]][i] <- as.numeric(names(tbl)[which.max(tbl)])
    }
  }
  return(df_poly)
}

add_covar2 <- function(df_poly, covars=c('NightFlag', ''), dat='all'){
  
  folderin<-paste0(dir,AlgoResults$input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  impact_sample_flattened <- flattenImpactSample(imp_best)
  columns_obs_sampled <- c(grep('observed', colnames(impact_sample_flattened)), grep('sampled.', colnames(impact_sample_flattened)))
  impact_sample_flattened$obs_quantile <- apply(impact_sample_flattened, 1, function(x) sample_quant(x[columns_obs_sampled]))
  impact_sample_flattened$samp_median <- apply(impact_sample_flattened, 1, function(x) median(as.numeric(x[grep('sampled.', colnames(impact_sample_flattened))])))
  #hist(impact_sample_flattened$obs_quantile[impact_sample_flattened$samp_median>0], breaks=c(0,10,20,30, 40, 50, 60, 70, 80, 90, 101), xlab='Observation Quantile', main='')
  #hist(impact_sample_flattened$obs_quantile[impact_sample_flattened$samp_median], breaks=c(0,10,20,30, 40, 50, 60, 70, 80, 90, 101), xlab='Observation Quantile', main='')
  df_plot <- impact_sample_flattened
  hist(df_plot$obs_quantile[impact_sample_flattened$samp_median>0 &  df_plot$impact=='mortality'], breaks=c(0,10,20,30, 40, 50, 60, 70, 80, 90, 101), xlab='Observation Quantile', main='')
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  covar <- 'SHDI'
  df_plot[[covar]] <- NA
  df_plot[['night_main_hazard']] <- NA
  
  for(i in 1:NROW(df_plot)){
    file_match <- grep(paste0("_", df_plot$event_id[i], "\\b"),  ufiles, value = TRUE)
    if (length(file_match) > 1) stop('Multiple ODD files for event with the same ID.')
    ODDy <- readRDS(paste0(folderin,file_match))
    poly_indexes <- ODDy@polygons[[df_plot$polygon[i]]]$indexes
    # if ('hazMean' %in% covars){
    #   poly_pop_restricted <- poly_indexes[which(ODDy@data$Population[poly_indexes] > 1000)]
    #   df_poly[['hazMean']][i] <- max(ODDy@data[poly_pop_restricted,grep("hazMean",names(ODDy),value = T)], na.rm=T)
    # } 
    #df_poly[[covar]][i] <- median(ODDy@data[[covar]][poly_indexes])
    #find mode:
    max_haz <- which(as.matrix(ODDy@data[,grep("hazMean",names(ODDy@data))] == max(ODDy@data[,grep("hazMean",names(ODDy@data))], na.rm=T)), arr.ind=T)[1,2]
    hour <- as.numeric(substr(ODDy@hazinfo$eventtimes, 1, 2))
    night_flag <- ifelse(hour>=22 | hour < 6, 1, 0)
    df_plot[i, 'night_main_hazard'] <- night_flag[max_haz]
    exposed_restrict <- which(ODDy@data[poly_indexes,paste0('hazMean', max_haz)] > 6 & ODDy@data$Population[poly_indexes] > 100)
    print(length(exposed_restrict))
    if (length(exposed_restrict)==0){
      exposed_restrict = poly_indexes
    }
    df_plot[[covar]][i] <- mean(ODDy@data[[covar]][poly_indexes[exposed_restrict]], na.rm=T)
  }
  ggplot(df_plot %>% filter(samp_median > 0 & night_main_hazard==1), aes(x=SHDI, y=obs_quantile)) + geom_point()
  plot(df_plot[[covar]], df_plot$obs_quantile, col=ifelse(df_plot$night_main_hazard, 'red', 'blue'), xlab='SHDI', ylab='Observation Quantile')
  
  #return(df_poly)
}



# ------------- Regionalise countries according to PAGER proposed regionalisation --------------

region_list <- list(
  list(name = "Australia, USA, and Canada", countries = c("Australia", "Canada", "United States (without California)", "Mexico", "Saint Pierre and Miquelon", "United States Minor Outlying Islands")),
  list(name = "New Zealand and California", countries = c("New Zealand", "California (USA)")),
  list(name = "Central America", countries = c("Costa Rica", "Panama")),
  list(name = "South Central America", countries = c("Dominican Republic", "Jamaica", "Guadeloupe", "El Salvador")),
  list(name = "Caribbean and Central America", countries = c("Guatemala", "Belize", "Honduras", "Nicaragua", "Haiti", "Puerto Rico", "Cayman Islands", "Turks and Caicos Islands", "Anguilla", "Montserrat", "Cuba", "Bahamas", "Saint Kitts and Nevis", "Saint Lucia", "Antigua and Barbuda", "Trinidad and Tobago", "Aruba", "Netherlands Antilles", "Dominica", "Grenada", "Saint Vincent and the Grenadines", "Martinique", "British Virgin Islands", "U.S. Virgin Islands", "Barbados", "Saint Barthelemy", "Saint Martin (France)")),
  list(name = "Western South America", countries = c("Colombia", "Ecuador", "Peru", "Chile", "Argentina", "Uruguay", "Brazil", "Paraguay")),
  list(name = "Eastern South America", countries = c("Venezuela", "Bolivia", "Brazil", "Uruguay", "Guyana", "Suriname", "Paraguay", "French Guiana")),
  list(name = "North Africa", countries = c("Algeria", "Egypt", "Tunisia", "Western Sahara")),
  list(name = "South-central Africa", countries = c("Botswana", "Namibia", "South Africa", "Swaziland", "Zimbabwe", "Morocco", "Sudan", "Chad", "Central African Republic", "Cameroon", "Congo", "DRP Congo", "Gabon", "Equatorial Guinea", "Sao Tome and Principe", "Angola", "Mauritania", "Senegal", "Gambia", "Guinea-Bissau", "Sierra Leone", "Liberia", "Cote d'Ivoire", "Ghana", "Togo", "Benin", "Niger", "Nigeria", "Mali", "Burkina Faso", "Guinea", "Yemen", "Eritrea", "Djibouti", "Ethiopia", "Somalia", "Kenya", "Uganda", "Rwanda", "Burundi", "United Republic of Tanzania", "Malawi", "Madagascar", "Mozambique", "Zambia", "Lesotho")),
  list(name = "Italy", countries = c("Italy", "Holy See", "Malta", "San Marino")),
  list(name = "Northern Europe", countries = c("Norway", "Sweden", "Finland", "Denmark", "Germany", "Belgium", "France", "Austria", "Switzerland", "Aland Islands", "Monaco", "Poland", "Bouvet Island", "United Kingdom", "Ireland", "Guernsey", "Isle of Man", "Jersey", "Falkland Islands (Malvinas)", "Saint Helena", "South Georgia and the South Sandwich Islands", "Iceland", "Faroe Islands", "Greenland", "Svalbard and Jan Mayen", "Liechtenstein", "Luxembourg", "Netherlands", "Greece", "Spain", "Portugal", "Gibraltar", "Cape Verde", "Andorra")),
  list(name = "Eastern Europe", countries = c("Czech Republic", "Slovenia", "Slovakia", "Hungary", "Bosnia and Herzegovina", "Croatia", "Serbia", "Montenegro", "Romania", "Albania", "Former Yugoslav Republic of Macedonia", "Bulgaria", "Republic of Moldova")),
  list(name = "Baltic States and Russia", countries = c("Estonia", "Latvia", "Lithuania", "Belarus", "Ukraine", "Russian Federation", "Georgia", "Armenia", "Azerbaijan")),
  list(name = "Central Asia", countries = c("Kazakhstan", "Uzbekistan", "Turkmenistan", "Kyrgyzstan", "Tajikistan")),
  list(name = "Arabian Peninsula", countries = c("Turkey", "Oman", "United Arab Emirates", "Qatar", "Saudi Arabia", "Bahrain", "Kuwait", "Lebanon", "Jordan", "Palestinian Territory", "Syrian Arab Republic", "Israel", "Cyprus", "Libyan Arab Jamahiriya")),
  list(name = "Iran & Iraq", countries = c("Iran", "Iraq", "Afghanistan", "Pakistan")),
  list(name = "Chinese Peninsula", countries = c("Brunei Darussalam", "China", "North Korea", "South Korea", "Macao", "Mongolia")),
  list(name = "Philippines and Malaysian Peninsula", countries = c("Singapore", "Thailand", "Hong Kong", "Malaysia", "Philippines")),
  list(name = "Indian Peninsula", countries = c("India", "Sri Lanka", "Bangladesh", "Nepal", "Bhutan", "Myanmar")),
  list(name = "Indonesian Peninsula", countries = c("Cambodia", "Laos", "Timor-Leste", "Vietnam", "Papua New Guinea", "American Samoa", "Samoa", "Tokelau", "Tuvalu", "Fiji", "Tonga", "Vanuatu", "Wallis and Futuna", "Niue", "Nauru", "New Caledonia", "Solomon Islands", "Palau", "Guam", "Northern Mariana Islands", "Marshall Islands", "Federated States of Micronesia", "Kiribati", "Cook Islands", "French Polynesia", "Norfolk Island", "Pitcairn", "British Indian Ocean Territory", "Christmas Island", "Cocos (Keeling) Islands", "French Southern Territories", "Heard Island and McDonald Islands", "Maldives", "Comoros", "Mauritius", "Mayotte", "Reunion", "Seychelles", "Bermuda", "Antarctica", "Indonesia")),
  list(name = "Japan & South Korea", countries = c("Japan", "Taiwan"))
)

#simplify region list given that we just know country, and so that it matches output of countrycode()
# region_list[[1]]$countries[3] <- 'United States'
# region_list[[2]]$countries <- region_list[[2]]$countries[1]
# region_list[[9]]$countries[which(region_list[[9]]$countries=='United Republic of Tanzania')] = 'Tanzania'
# region_list[[9]]$countries[which(region_list[[9]]$countries=='DRP Congo')] = 'Congo - Kinshasa'
# region_list[[12]]$countries[which(region_list[[12]]$countries=='Former Yugoslav Republic of Macedonia')] = 'North Macedonia'
# region_list[[15]]$countries[11] <- 'Syria'
# region_list[[19]]$countries[which(region_list[[19]]$countries=='Myanmar')] = 'Myanmar (Burma)'


iso3_to_region <- function(iso3, region_list) {
  #having some issues with getting some country name matches with countrycode() but this seems to work ok for now:
  if (iso3 %in% c('TTO', 'BIH')){
    country <- countrycode(iso3, origin='iso3c', destination='iso.name.en')
  } else {
    country <- ifelse(iso3=='KOS', 'Serbia', countrycode(iso3, origin='iso3c', destination='country.name'))
  }
  for (i in seq_along(region_list)) {
    sublist <- region_list[[i]]$countries
    if (country %in% sublist) {
      return(region_list[[i]]$name)
    }
  }
  return(NA)  # Country not found in the list
}

add_regions <- function(df_poly){
  df_poly$region <- NA
  for (i in 1:NROW(df_poly)){
    iso3=df_poly$iso3[i]
    if (df_poly$iso3[i] == 'TOT'){
      df_poly_event <- filter(df_poly, event_id==df_poly$event_id[i] & iso3 != 'TOT')
      if (NROW(df_poly_event)==0){
        if (df_poly$event_id[i]==81){ iso3='PNG'}
        if (df_poly$event_id[i]==90){ iso3='CRI'}
        if (iso3=='TOT') {stop(df_poly$event_id[i])}
      } else {
        iso3=df_poly_event$iso3[which.max(df_poly_event$observed)]
      }
    }
    df_poly$region[i] <- iso3_to_region(iso3, region_list)
  }
  return(df_poly)
}

add_landslide_flag <- function(df_poly){
  EQIL <- read.csv('/home/manderso/Downloads/EQIL_Database_2022.csv', stringsAsFactors=FALSE, fileEncoding="latin1")[, c('Event', 'Country', 'Year', 'Month', 'Day', 'Mw', 'Fault.Type', 'Depth..km.', 'Ms', 'Landslide.Fatalities')]
  df_poly$landslide_flag <- F
  df_poly$Landslide.Fatalities <- NA
  df_poly$Mw <- NA
  df_poly$Fault.Type <- NA
  df_poly$Depth <- NA
  df_poly$Ms <- NA
  
  EQIL$date <- as.Date(paste(EQIL$Year, EQIL$Month, EQIL$Day, sep = "-"))
  EQIL %<>% filter(date > '2007-01-01')
  for (i in unique(df_poly$event_id)){
    EQIL_match <- c()
    corr_rows <- which(df_poly$event_id==i)
    date_match <- which(EQIL$date > df_poly$date[corr_rows[1]] - 1 & EQIL$date < df_poly$date[corr_rows[1]] + 1)
    for(j in date_match){
      if(EQIL$Country[j] %in% countrycode(strsplit(df_poly$iso3[corr_rows], ' ')[[1]], origin='iso3c', destination='country.name')){
        EQIL_match <- c(EQIL_match, j)
      }
    }
    if (length(EQIL_match)>1){
      stop()
    }
    if (length(EQIL_match)==1){
      df_poly$landslide_flag[corr_rows] <- T
      df_poly$Landslide.Fatalities[corr_rows] <- EQIL$Landslide.Fatalities[EQIL_match]
      df_poly$Mw[corr_rows] <- EQIL$Mw[EQIL_match]
      df_poly$Fault.Type[corr_rows] <- EQIL$Fault.Type[EQIL_match]
      df_poly$Depth[corr_rows] <- EQIL$Depth..km.[EQIL_match]
      df_poly$Ms[corr_rows] <- EQIL$Ms[EQIL_match]
    }
  }
  return(df_poly)
}

add_hazard_info <- function(df_poly){
  folderin <- paste0(dir,"IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/HAZARDobjects_wMaxMMIDiff/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
  
  df_poly$depth <- NA
  df_poly$max_mmi <- NA
  df_poly$magnitude <- NA
  df_poly$time <- NA
  df_poly$main_shock_first <- NA
  
  for(i in 1:NROW(df_poly)){
    file_match <- grep(paste0("_", df_poly$event_id[i], "\\b"),  ufiles, value = TRUE)
    HAZARDobj <- readRDS(paste0(folderin, file_match))
    max_mmi_i <- which.max(HAZARDobj$hazard_info$max_mmi)
    df_poly$max_mmi[i] <- HAZARDobj$hazard_info$max_mmi[max_mmi_i]
    df_poly$depth[i] <- HAZARDobj$hazard_info$depth[max_mmi_i]
    df_poly$magnitude[i] <- max(HAZARDobj$hazard_info$magnitude)
    df_poly$time[i] <- HAZARDobj$hazard_info$eventtime[max_mmi_i]
    datetimes <- paste0(HAZARDobj$hazard_info$eventdates, HAZARDobj$hazard_info$eventtimes)
    df_poly$main_shock_first[i] <- (max_mmi_i == ifelse(length(which(HAZARDobj$hazard_info$first_event))>0,which(HAZARDobj$hazard_info$first_event), 0 ))
  }
  
  return(df_poly)
}

addSHDI <- function(df_poly){ 
  GDLdata <- readGlobalDataLab()
  df_poly$shdi <- NA
  df_poly$healthindex <- NA
  df_poly$ExpSchYrs <- NA
  df_poly$sgdi <- NA
  for (i in 1:NROW(df_poly)){
    GDLdata_i <- which(GDLdata$AveSchYrs==df_poly$AveSchYrs[i] & GDLdata$LifeExp==df_poly$LifeExp[i] & GDLdata$GNIc==df_poly$GNIc[i])
    if (length(GDLdata_i)> 1){
      stop(paste('Duplicate Matches:', i))
    } 
    if (length(GDLdata_i)==0){
      print(paste('No Match', df_poly$iso3[i]))
      next
    }
    df_poly$shdi[i] <- GDLdata[GDLdata_i, 'shdi']
    df_poly$healthindex[i] <- GDLdata[GDLdata_i, 'healthindex']
    df_poly$ExpSchYrs[i] <- GDLdata[GDLdata_i, 'ExpSchYrs']
    df_poly$sgdi[i] <- GDLdata[GDLdata_i, 'sgdi']
  }
  
}



# -----------------------------------------------------------------------------------------------

manual_modify_params = function(AlgoResults, s_finish=NULL){
  M=3
  AlgoResults %<>% addAlgoParams(s_finish)
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  particle_min.d <- 40
  omega1 <- AlgoResults$Omega_sample_phys[particle_min.d, , AlgoResults$s_finish] %>% relist(Model$skeleton)
  omega1$eps$local = 0
  omega1$eps$hazard = 0
  omega2 <- omega1
  omega2$vuln_coeff$AveSchYrs = -0.1
  omega2$vuln_coeff$GNIc <- 0
  omega2$Lambda2$mu <- 11.68
  #omega2$vuln_coeff$GNIc = -0.05
  Model$HighLevelPriors(omega1 %>% addTransfParams(), Model)
  Model$HighLevelPriors(omega2 %>% addTransfParams(), Model)
  
  df_all1 <- sample_post_predictive(AlgoResults, 3, AlgoResults$s_finish, dat='train', single_particle=T, Omega = omega1)
  df_all2 <- sample_post_predictive(AlgoResults, 3, AlgoResults$s_finish, dat='train', single_particle=T, Omega = omega2)
  
  df_poly1 <- df_all1$poly
  df_poly2 <- df_all2$poly
  sum(df_poly1$crps); sum(df_poly2$crps)
  
  df_point1 <- df_all1$point
  df_point2 <- df_all2$point
  sum(df_point1[,2]); sum(df_point2[,2])
  
  impact_sample1 <- SampleImpact(dir = dir,Model = Model,
                                proposed = omega1 %>% addTransfParams(), 
                                AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 3) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                dat='Train')
  
  impact_sample2 <- SampleImpact(dir = dir,Model = Model,
                                proposed = omega2 %>% addTransfParams(), 
                                AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 3) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                dat='Train')
  
  sum(df_poly1$crps)
  sum(df_poly2$crps)
  
  dist1 <- CalcDist(impact_sample1, AlgoParams)
  dist2 <- CalcDist(impact_sample2, AlgoParams)
  hist(AlgoResults$d[,1,AlgoResults$s_finish])
  Model$HighLevelPriors(omega2 %>% addTransfParams(), Model, AlgoParams)
  
  ix <- which(impact_sample2$poly[[1]]$impact=='mortality')
  plot(log(impact_sample1$poly[[1]]$observed[ix] + k), log(impact_sample1$poly[[1]]$sampled[ix]+k) - (log(impact_sample1$poly[[1]]$observed[ix] + k)))
  points(log(impact_sample1$poly[[1]]$observed[ix] + k), log(impact_sample2$poly[[1]]$sampled[ix]+k) - (log(impact_sample2$poly[[1]]$observed[ix] + k)), col='red')
  
  
}



plot_covar_vs_error = function(AlgoResults, covar='EQFreq', dat='all', s_finish=NULL){
  #plot covariates against discrepancy between sampled and observed
  
  M <- 10
  AlgoResults %<>% addAlgoParams(s_finish)
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  proposed <- relist(AlgoResults$Omega_sample_phys[particle_min.d[1],,AlgoResults$s_finish], skeleton=Model$skeleton)
  #proposed$vuln_coeff$Vs30=0
  # proposed$vuln_coeff$Night=0
  # proposed$vuln_coeff$FirstHaz.Night=0
  # proposed$check$check = 0.5
  #df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat=dat, single_particle=T, Omega = proposed)
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat=dat, single_particle=T, particle_i = particle_min.d[1])
  
  #df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat=dat, single_particle=T, Omega = proposed)
  df_poly$sampled_median <- apply(df_poly[grep("sampled", names(df_poly))], 1, median)
  df_poly$sampled_mean <- apply(df_poly[grep("sampled", names(df_poly))], 1, mean)
  
  
  impact_type <- 'buildDam'
  k <- 10
  
  df_poly %>% filter(impact=='mortality' & log(df_poly$sampled_median+10) < 3 & log(df_poly$observed+10)>5)
  df_poly %>% filter(impact=='mortality' & abs(log(df_poly$sampled_median+10) - log(df_poly$observed+10))>2.8)
  df_poly %>% filter(impact=='mortality' & log(df_poly$sampled_median+10) - log(df_poly$observed+10)>1.5)
  
  ggplot(df_poly %>% filter(impact==impact_type & inferred==F), aes(x= log(observed+k), y=log(sampled_median+k))) +
    geom_point(aes(color=train_flag)) +theme_minimal() + ggtitle(impact_type) +
    scale_color_manual(values = c("TRAIN" = "black", "TEST" = "red")) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + coord_fixed(ratio = 1)
  
  
  df_poly %>% filter(impact=='mortality' & train_flag=='TEST' & log(observed+k)>5 & log(sampled_median+k) < 3)
  
  df_poly_filt <- df_poly %>% filter(impact==impact_type & inferred==F)
  # df_poly_filt2 <- df_poly2 %>% filter(impact==impact_type  & inferred==F)
  # abline(0,1)
  # plot(log(df_poly_filt$observed+10), log(df_poly_filt$sampled_median+10))
  # points(log(df_poly_filt2$observed+10), log(df_poly_filt2$sampled_median+10), col='red')
  # 
  # plot(log(df_poly_filt2$sampled_median+10)-log(df_poly_filt2$observed+10), log(df_poly_filt$sampled_median+10)-log(df_poly_filt$observed+10) )
  # abline(0,1)
  # 
  # df_poly_filt2[ log(df_poly_filt$sampled_median+10)-log(df_poly_filt$observed+10)>2.1,]
  # df_poly_filt[ log(df_poly_filt$sampled_median+10)-log(df_poly_filt$observed+10)>2.1,]

  cor(df_poly_filt$sampled_median, df_poly_filt$observed)
  cor(log(df_poly_filt$sampled_median+k), log(df_poly_filt$observed+k))

  SS_res_raw = sum((df_poly_filt$observed-df_poly_filt$sampled_median)^2)
  SS_tot_raw = sum((df_poly_filt$observed-mean(df_poly_filt$observed))^2)
  R_squared_raw = 1-SS_res_raw/SS_tot_raw
  R_squared_raw

  SS_res_log = sum((log(df_poly_filt$observed+k)-log(df_poly_filt$sampled_median+k))^2)
  SS_tot_log = sum((log(df_poly_filt$observed+k)-mean(log(df_poly_filt$observed+k)))^2)
  R_squared_log = 1-SS_res_log/SS_tot_log
  R_squared_log
  
  
  df_poly %>% filter(impact==impact_type & abs((log(sampled_median+k)-log(observed+k)))>4)
  
  df_poly %>% filter(impact==impact_type)
  
  df_poly %<>% add_hazard_info()
  df_poly %<>% add_covar(covars=c('hazMean', 'EQFreq', 'GNIc', 'Vs30', 'AveSchYrs', 'LifeExp', 'SHDI'), dat=dat)
  
  
  plot(df_poly$SHDI, log(df_poly$GNIc))
  
  ggplot(df_poly %>% filter(impact==impact_type & inferred==F & train_flag=='TRAIN'), aes(x= as.POSIXct(time, format = "%H:%M:%S"), y=log(sampled_median+k)-log(observed+k))) +
    geom_point(aes(color=main_shock_first)) +theme_minimal() + ggtitle(impact_type) #+
    #scale_color_manual(values = c("TRAIN" = "black", "TEST" = "red"))
  
  ggplot(df_poly %>% filter(impact==impact_type & inferred==F), aes(x= Vs30, y=log(sampled_median+k)-log(observed+k))) +
    geom_point(aes(color=train_flag)) +theme_minimal() + ggtitle(impact_type) +
    scale_color_manual(values = c("TRAIN" = "black", "TEST" = "red"))
  
 df_poly %>% filter(impact==impact_type & inferred==F & main_shock_first & log(sampled_median+k)-log(observed+k) > 1)
         
  time2 <- as.numeric(format(as.POSIXct(df_poly$time, format = "%H:%M:%S"), "%H")) + as.numeric(format(as.POSIXct(df_poly$time, format = "%H:%M:%S"), "%M") )/60
  
  df_poly$sin_time = cos(2*pi * (time2 - 2) / 24)
  df_poly$night_flag = ifelse(as.numeric(format(as.POSIXct(df_poly$time, format = "%H:%M:%S"), "%H")) < 6 | as.numeric(format(as.POSIXct(df_poly$time, format = "%H:%M:%S"), "%H")) > 20, T, F)
  
  ggplot(df_poly %>% filter(impact==impact_type & inferred==F & train_flag=='TRAIN'), aes(x= -sin_time, y=log(sampled_median+k)-log(observed+k))) +
    geom_point(aes(color=main_shock_first)) +theme_minimal() + ggtitle(impact_type) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red"))
  
  ggplot(df_poly %>% filter(impact==impact_type & inferred==F), aes(x= night_flag, y=log(sampled_median+k)-log(observed+k))) +
    geom_point(aes(color=main_shock_first)) +theme_minimal() + ggtitle(impact_type) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red"))
  
  ggplot(df_poly %>% filter(impact==impact_type & inferred==F), aes(x= main_shock_first, y=log(sampled_median+k)-log(observed+k))) +
    geom_point(aes(color=night_flag)) +theme_minimal() + ggtitle(impact_type) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red"))
  
  ggplot(df_poly %>% filter(impact==impact_type & inferred==F), aes(x= Vs30, y=log(sampled_median+k)-log(observed+k))) +
    geom_point(aes(color=train_flag)) +theme_minimal() + ggtitle(impact_type) +
    scale_color_manual(values = c("TRAIN" = "black", "TEST" = "red"))
  
  
  
  df_poly %>% filter(impact==impact_type & inferred==F & night_flag & (log(sampled_mean+k)-log(observed+k)>2))
  
  ggplot(df_poly %>% filter(impact==impact_type), aes(x= log(observed+k), y=log(sampled_mean+k))) + 
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color=AveSchYrs))  +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  #df_poly %>% filter(magnitude < 6.5 & (log(sampled_median+k) - log(observed+k) < -2.5) & impact==impact_type)
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = log(EQFreq+0.1), y = log(sampled_mean + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = log(EQFreq))) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = log(GNIc), y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = train_flag)) + scale_color_manual(values = c("TRAIN" = "black", "TEST" = "red"))
    #scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = AveSchYrs, y = log(sampled_mean + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = train_flag)) + scale_color_manual(values = c("TRAIN" = "black", "TEST" = "red"))
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type & train_flag=='TRAIN'), aes(x = LifeExp, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = LifeExp)) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = Vs30, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = Vs30)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn"))) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = LifeExp, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point(aes(color = train_flag, shape=train_flag))
  
  
  
  ggplot(df_poly %>% filter(impact == impact_type & train_flag=='TRAIN'), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = AveSchYrs)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn"))) + theme_minimal()
  
  
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = log(depth), y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(aes(color = magnitude)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn"))) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(aes(color = hazMean)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn"))) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(aes(color = log(depth))) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = magnitude, y = log(sampled_median + k) - log(observed + k))) + geom_point()
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = max_mmi, y = log(sampled_median + k) - log(observed + k))) + geom_point()
  
  df_poly %<>% add_covar(covar='hazMean', dat=dat)
  df_poly %<>% add_covar(covar='AveSchYrs', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) +
    geom_point(aes(color = AveSchYrs)) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k)))
  
  df_poly %<>% add_covar(covar='GNIc', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = log(GNIc), y = log(sampled_median + k) - log(observed + k))) + 
    geom_point()
  
  df_poly %<>% add_covar(covar='LifeExp', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = LifeExp, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point()
  
  df_poly %<>% add_covar(covar='AveSchYrs', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = AveSchYrs, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point()
  
  df_poly %<>% add_covar(covar='Vs30', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = Vs30, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point(aes(color=as.factor(event_id))) + theme(legend.position='none')
  
  df_poly %<>% add_covar(covar='EQFreq', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = EQFreq, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point(aes(color=as.factor(event_id))) + theme(legend.position='none')
  
  #df_poly[ix[which(7.7 < df_poly$hazMean[ix] & df_poly$hazMean[ix] < 8.3)],]
  #which(log(df_poly$sampled_median[ix]+k) - log(df_poly$observed[ix]+k) > 2.5)
  
  #how reliable is landslide data?
  df_poly_landslide = add_landslide_flag(df_poly)
  ggplot(df_poly_landslide %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point(aes(color=landslide_flag)) + xlab('Maximum intensity to which Population > 1000 is exposed')
  
  # ggplot(df_poly_landslide %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed -ifelse(is.na(Landslide.Fatalities), 0,Landslide.Fatalities) + k))) + 
  #   geom_point(aes(color=landslide_flag))
  
  
  
  df_poly %<>% add_regions()
  ggplot(df_poly %>% filter(impact==impact_type), aes(x=region, y=log(sampled_median+k) - log(observed+k))) +
    geom_point(aes(color=as.factor(event_id))) + ggtitle(impact_type) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='none')
  
  ggplot(df_poly %>% filter(impact==impact_type), aes(x=region, y=log(sampled_median+k) - log(observed+k))) +
    geom_point(aes(color=Vs30)) + ggtitle(impact_type) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  plot(df_poly$EQFreq[ix], log(df_poly$sampled_median[ix]+k)-log(df_poly$observed[ix]+k), xlab=covar, ylab='log(Sampled+10) - log(Observed+10)')
  
  ggplot(df_poly, aes(x=GNIc, y=log(observed+k)-log(sampled_median+k))) + geom_point() 
  #ggplot(df_poly, aes(x=countrycode(iso3, origin='iso3c', destination='continent'), y=log(observed+k)-log(sampled_median+k))) + geom_point()
  
  
  #count number of points that lie within credible intervals
  df_poly_filt <- df_poly2 %>% filter(impact=='mortality')
  count_min <- 0
  count_max <- 0
  total <- 0
  for (i in 1:NROW(df_poly_filt)){
    if (all(c(df_poly_filt$observed[i], df_poly_filt[i, grep("sampled", names(df_poly_filt))])==0)){
      next
    }
    if(all(df_poly_filt$observed[i] > df_poly_filt[i, grep("sampled", names(df_poly_filt))])){
      count_max <- count_max + 1
    }
    if(all(df_poly_filt$observed[i] < df_poly_filt[i, grep("sampled", names(df_poly_filt))])){
      count_min <- count_min + 1
    }
    total <- total + 1
  }
  (count_max+count_min)/total # vs 2/(M+1)
  count_max /total
  count_min / total
  2/11
  
  
  
  abline(0,0)
}

# par(mfrow=c(3,2))
# k <- 1
# event_dat <- as.numeric(df_poly[1000,grep("sampled", names(df_poly))])
# hist(log(event_dat+k), freq=F, breaks=min(length(unique(event_dat)), 20))
# lines(seq(0, max(log(event_dat+k)), 0.1), dnorm(seq(0, max(log(event_dat+k)), 0.1),mean(log(event_dat+k)), sd(log(event_dat+k))))

plot_impact_curves = function(AlgoResults, s_finish=NULL){
  #plot covariates against discrepancy between sampled and observed
  
  M <- 5
  AlgoResults %<>% addAlgoParams(s_finish)
  
  n_post_samples <- 500
  particle_sample <- sample(1:AlgoResults$Npart, n_post_samples, prob=AlgoResults$W[,AlgoResults$s_finish], replace=T)
  I <- seq(4.5, 5.5, 0.01)
  Omega <- relist(AlgoResults$Omega_sample_phys[particle_sample[1],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()
  #fDamUnscaled(I,list(I0=Model$I0, Np=length(I)),relist(AlgoResults$Omega_sample_phys[particle_sample[1],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()) 
  Damage <- h_0(I,Model$I0, Omega)
  plot(I, D_MortDisp_calc(Damage, Omega)[2,], type='l', xlim=c(4.5, 5.5), ylim=c(0,0.3))
  lines(I, D_MortDisp_calc(Damage, Omega)[1,], col='red')
  for (i in 2:n_post_samples){
    Omega <- relist(AlgoResults$Omega_sample_phys[particle_sample[i],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()
    #fDamUnscaled(I,list(I0=Model$I0, Np=length(I)),relist(AlgoResults$Omega_sample_phys[particle_sample[1],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()) 
    Damage <- h_0(I,Model$I0, Omega)
    lines(I, D_MortDisp_calc(Damage, Omega)[2,])
    lines(I, D_MortDisp_calc(Damage, Omega)[1,], col='red')
    lines(I, D_DestDam_calc(Damage, Omega)[2,], col='blue')
    lines(I, D_DestDam_calc(Damage, Omega)[1,], col='green')
  }
}

plot_vuln = function(AlgoResults, dat='all', s_finish=NULL){
  
  folderin<-paste0(dir,AlgoResults$input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  AlgoResults %<>% addAlgoParams(s_finish)
  n_post_samples <- 5
  df_vuln <- data.frame(
    post_sample = integer(),
    iso3 = character(),
    vuln_overall = numeric(),
    date = as.Date(character()),  # Convert the date column to Date type
    PDens = numeric(),
    AveSchYrs = numeric(),
    LifeExp = numeric(),
    GNIc = numeric(),
    Vs30 = numeric(),
    EQFreq = numeric(),
    mag = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion of character to factors
  )
  for (i in 1:n_post_samples){
    particle_sample <- sample(1:AlgoResults$Npart, 1, prob=AlgoResults$W[,AlgoResults$s_finish], replace=T)
    Omega <- relist(AlgoResults$Omega_sample_phys[particle_sample[1],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()
    EQfreqs <- c()
    for (file in ufiles){
      ODDy <- readRDS(paste0(folderin, file))
      pop_restricted = which(ODDy@data$Population > 1000)
      max_haz <- max(ODDy@data[pop_restricted, grep("hazMean",names(ODDy),value = T)], na.rm=T)
      max_haz_i <- which(ODDy@data[pop_restricted, grep("hazMean",names(ODDy),value = T), drop=F] == max_haz, arr.ind=T)[1]
      EQfreqs <- c(EQfreqs, log(ODDy@data[pop_restricted[max_haz_i], 'EQFreq']+0.1))
      vuln_overall <- GetLP_single(Omega, Model$center, vuln_terms=list(PDens=ODDy@data[pop_restricted[max_haz_i], 'PDens'], 
                                                        AveSchYrs=ODDy@data[pop_restricted[max_haz_i], 'AveSchYrs'],
                                                        LifeExp=ODDy@data[pop_restricted[max_haz_i], 'LifeExp'],
                                                        GNIc=ODDy@data[pop_restricted[max_haz_i], 'GNIc'],
                                                        Vs30=ODDy@data[pop_restricted[max_haz_i], 'Vs30'],
                                                        EQFreq=ODDy@data[pop_restricted[max_haz_i], 'EQFreq'],
                                                        Mag=6.187425))
      df_vuln %<>% add_row(post_sample=i, 
                           iso3=ODDy@data[pop_restricted[max_haz_i], 'ISO3C'], 
                           vuln_overall=vuln_overall, 
                           date=ODDy@hazdates[1],
                           PDens=Omega$vuln_coeff$PDens * ((log(ODDy@data[pop_restricted[max_haz_i], 'PDens']+1) - Model$center$PDens$mean)/Model$center$PDens$sd),
                           AveSchYrs=Omega$vuln_coeff$AveSchYrs * ((ODDy@data[pop_restricted[max_haz_i], 'AveSchYrs'] - Model$center$AveSchYrs$mean)/Model$center$AveSchYrs$sd),
                           LifeExp=Omega$vuln_coeff$LifeExp * ((ODDy@data[pop_restricted[max_haz_i], 'LifeExp'] - Model$center$LifeExp$mean)/Model$center$LifeExp$sd),
                           GNIc=Omega$vuln_coeff$GNIc * ((log(ODDy@data[pop_restricted[max_haz_i], 'GNIc']) - Model$center$GNIc$mean)/Model$center$GNIc$sd),
                           Vs30=Omega$vuln_coeff$Vs30 * ((ODDy@data[pop_restricted[max_haz_i], 'Vs30'] - Model$center$Vs30$mean)/Model$center$Vs30$sd),
                           EQFreq=Omega$vuln_coeff$EQFreq * ((log(ODDy@data[pop_restricted[max_haz_i], 'EQFreq']+0.1) - Model$center$EQFreq$mean)/Model$center$EQFreq$sd), 
                           mag=Omega$vuln_coeff$Mag * ((6.187425 - Model$center$Mag$mean)/Model$center$Mag$sd))
    }
  }
  plot(df_vuln$date, df_vuln$EQFreq, col='blue',  ylim=c(-0.3,0.3))
  points(df_vuln$date, df_vuln$GNIc, col='green')
  points(df_vuln$date, df_vuln$AveSchYrs)
  points(df_vuln$date, df_vuln$PDens, col='red')
  points(df_vuln$date, df_vuln$Vs30, col='yellow')
  points(df_vuln$date, df_vuln$LifeExp, col='orange')
  
}

plot_vuln_corr = function(AlgoResults, dat='all'){
  
  folderin<-paste0(dir,AlgoResults$input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  df_vuln <- data.frame(
    iso3 = character(),
    date = as.Date(character()),  # Convert the date column to Date type
    PDens = numeric(),
    AveSchYrs = numeric(),
    LifeExp = numeric(),
    GNIc = numeric(),
    Vs30 = numeric(),
    EQFreq = numeric(),
    mag = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion of character to factors
  )
  
  for (file in ufiles){
    ODDy <- readRDS(paste0(folderin, file))
    pop_restricted = which(ODDy@data$Population > 1000)
    max_haz <- max(ODDy@data[pop_restricted, grep("hazMean",names(ODDy),value = T)], na.rm=T)
    max_haz_i <- which(ODDy@data[pop_restricted, grep("hazMean",names(ODDy),value = T), drop=F] == max_haz, arr.ind=T)[1]
    df_vuln %<>% add_row(iso3=ODDy@data[pop_restricted[max_haz_i], 'ISO3C'], 
                         date=ODDy@hazdates[1],
                         PDens=((log(ODDy@data[pop_restricted[max_haz_i], 'PDens']+1) - Model$center$PDens$mean)/Model$center$PDens$sd),
                         AveSchYrs=((ODDy@data[pop_restricted[max_haz_i], 'AveSchYrs'] - Model$center$AveSchYrs$mean)/Model$center$AveSchYrs$sd),
                         LifeExp=((ODDy@data[pop_restricted[max_haz_i], 'LifeExp'] - Model$center$LifeExp$mean)/Model$center$LifeExp$sd),
                         GNIc=((log(ODDy@data[pop_restricted[max_haz_i], 'GNIc']) - Model$center$GNIc$mean)/Model$center$GNIc$sd),
                         Vs30=((ODDy@data[pop_restricted[max_haz_i], 'Vs30'] - Model$center$Vs30$mean)/Model$center$Vs30$sd),
                         EQFreq=((log(ODDy@data[pop_restricted[max_haz_i], 'EQFreq']+0.1) - Model$center$EQFreq$mean)/Model$center$EQFreq$sd))
  }
  plot(df_vuln[,c('PDens', 'AveSchYrs', 'LifeExp', 'GNIc', 'Vs30', 'EQFreq')], xlim=c(-4,4), ylim=c(-4,4))
}

plot_fitted_vuln_coefficients = function(AlgoResults, s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  var_posts <- data.frame(AlgoResults$Omega_sample_phys[,15:22,AlgoResults$s_finish])
  var_prior <- data.frame(AlgoResults$Omega_sample_phys[,15:22,1])
  names(var_posts) <- names(var_prior) <- sub("^[^.]*\\.\\s*", "", names(unlist(Omega))[15:22])
  par(mfrow=c(2,4))
  for (var in c('PDens', 'SHDI', 'GNIc', 'Vs30', 'EQFreq', 'FirstHaz', 'Night', 'FirstHaz.Night')){
    hist(var_posts[[var]], main=var, xlab='', xlim=c(-1,1), freq=F)
    lines(density(var_prior[[var]]), col='blue')
    #lines(seq(-1,1,0.01), dLaplace(seq(-1,1,0.01), 0,0.25))
    abline(v=0, col='red')
  }
  
  for (i in 15:22){
    p_H_given_y = sum(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish]>0)/length(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish])
    p_H = sum(AlgoResults$Omega_sample_phys[,i,1]>0)/length(AlgoResults$Omega_sample_phys[,i,1])
    print(paste('Bayes factor for', names(unlist(Omega))[i], ':', round(p_H_given_y * (1-p_H)/((1-p_H_given_y)*p_H), 3)))
    print(paste('Post. prob of greater than 0:', sum(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish]>0)/length(AlgoResults$Omega_sample_phys[,i,AlgoResults$s_finish])))
  }
  
  p <- list()
  i <- 1
  for (var in c('PDens', 'SHDI', 'GNIc', 'Vs30', 'EQFreq', 'FirstHaz', 'Night', 'FirstHaz.Night')){
    p[[i]] <- ggplot(var_posts, aes(x = !!ensym(var))) +
      geom_histogram(color = "black", alpha = 0.7) +
      ggtitle(paste0("Coefficient for ", var ))+
      labs(x='', y='') + 
      #scale_fill_discrete("Legend", values = dd.col) + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      geom_vline(xintercept = 0, col='red') + 
      #scale_x_log10()+
      guides(fill = FALSE) 
    i <- i + 1
  }
  do.call(grid.arrange, p)
  
  ggplot(var_posts, aes(x=SHDI, y=GNIc)) + geom_point() + geom_abline(slope=-1, intercept=0, col='red')
  ggplot(var_posts, aes(x=FirstHaz, y=FirstHaz.Night)) + geom_point() + geom_abline(slope=-1, intercept=0, col='red')
  
}

MAD_investigation = function(AlgoResults){
  for (i in 1:6){
    d <- AlgoResults$d_full[,i,1]
    d_mad <- sum(abs(d-mean(d)))/length(d) 
    if (i > 3){d_mad <- d_mad / 10}
    print(d_mad)
  }
}

plot_post_predictive = function(AlgoResults, M, s_finish=NULL){
  #mapply(crps, df_plot[grep("^sampled\\.", names(df_sampled))], df_plot$observed)
  AlgoResults %<>% addAlgoParams(s_finish)
  
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat='all')
  
  df_poly$q10 <- NA
  df_poly$q90 <- NA
  for (i in 1:NROW(df_poly)){
    df_poly$q10[i] <- as.numeric(quantile(as.numeric(df_poly[i, grep("sampled", names(df_poly))]), 0.1))
    df_poly$q90[i] <- as.numeric(quantile(as.numeric(df_poly[i, grep("sampled", names(df_poly))]), 0.9))
  }
  
  ggplot(df_poly %>% filter(impact=='displacement' & train_flag=='TEST'), aes(x=log(observed+10), ymin=log(q10+10), ymax=log(q90+10), col=train_flag)) + 
    geom_errorbar() + geom_abline(intercept=0,slope=1) +ylab('Sampled range across 10 simulations') +
    scale_color_manual(values = c("TRAIN" = "black", "TEST" = "red"))
  
  
  df_long_sampled <-   df_poly %>% pivot_longer(paste0('sampled.', 1:M))
  df_long_sampled$log_observed <- log(df_long_sampled$observed+10)
  df_long_sampled$log_sampled <- log(df_long_sampled$value+10)
  
  
  # Create the ggplot
  
  id_plot <- sample((df_long_sampled %>% filter(impact=='mortality' & train_flag=='TEST'))$obs_id, 50, replace=F)
  ggplot(df_long_sampled %>% filter(impact=='mortality' & train_flag=='TEST' & obs_id %in% id_plot), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
  
  ggplot(df_long_sampled %>% filter(impact=='mortality', observed != 0, obs_id >1000), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
  
  geom_density(alpha = 0.5) + geom_vline(aes(xintercept = observed), col='red') + theme_minimal()
  
  ggplot(df_sampled, aes(x=)) + + geom_density()
  
}

compare_mort_disp <- function(){
  
  df_poly$sampled_median <- apply(df_poly[grep("sampled", names(df_poly))], 1, median)
  df_poly$sampled_mean <- apply(df_poly[grep("sampled", names(df_poly))], 1, mean)
  
  merged_df <- merge(df_poly %>% filter(impact=='mortality') %>% transmute(obs_mort=observed,
                                                                        samp_mort=sampled_median,
                                                                        event_id=event_id,
                                                                        polygon=polygon),
                     df_poly %>% filter(impact=='displacement') %>% transmute(obs_disp=observed,
                                                                            sam_disp=sampled_median,
                                                                            event_id=event_id,
                                                                            polygon=polygon),
                     by=c("event_id","polygon"))
  
  ggplot(merged_df, aes(x=log(samp_mort+k)-log(obs_mort+k), y=log(sam_disp+k)-log(obs_disp+k))) + geom_point()
  
  
  
}


model_deepdive = function(AlgoResults, M, s_finish=NULL){
  #mapply(crps, df_plot[grep("^sampled\\.", names(df_sampled))], df_plot$observed)
  AlgoResults %<>% addAlgoParams(s_finish)
  
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish)
  
  df_long_sampled <-   df_poly %>% pivot_longer(paste0('sampled.', 1:M))
  df_long_sampled$log_observed <- log(df_long_sampled$observed+10)
  df_long_sampled$log_sampled <- log(df_long_sampled$value+10)
  
  # Create the ggplot
  ggplot(df_long_sampled %>% filter(impact=='mortality', obs_id %in% 900:1000), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
  
  ggplot(df_long_sampled %>% filter(impact=='mortality', observed != 0, obs_id >1000), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
  
  geom_density(alpha = 0.5) + geom_vline(aes(xintercept = observed), col='red') + theme_minimal()
  
  ggplot(df_sampled, aes(x=)) + + geom_density()
  
}


# plot_correlated_posteriors(AlgoResults, include_priors=T)
# plot_density_vs_step(AlgoResults, Omega)


plot_satellite_data = function(AlgoResults, dat='all', s_finish=NULL){

  M <- 5
  AlgoResults %<>% addAlgoParams(s_finish)
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat=dat, single_particle=T, particle_i = particle_min.d[1], return_type='all')
  df_point <- df_poly$point
  
  NROW(df_point)
  colnames(df_point)[2] <- 'distance'
  df_point_class <- rownames(df_point)
  df_point %<>% as.data.frame()
  rownames(df_point) <- NULL
  df_point$classification <- df_point_class
  df_point$impact_true <- ifelse(df_point$classification %in% c('N11', 'N12', 'N13'), 'Unaffected', ifelse(df_point$classification %in% c('N21', 'N22', 'N23'), 'Damaged', 'Destroyed'))
  df_point$impact_sampled <- ifelse(df_point$classification %in% c('N11', 'N21', 'N31'), 'Unaffected', ifelse(df_point$classification %in% c('N12', 'N22', 'N32'), 'Damaged', 'Destroyed'))
  df_point$event_i <- (1:NROW(df_point)- 1) %/% 9 + 1
  
  
  p1 <- ggplot(df_point %>% filter(impact_true=='Unaffected'), aes(fill= impact_sampled, alpha = (impact_sampled==impact_true),y=ifelse(count==0, 0, log(count)), x=event_i)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_alpha_manual(values=c(0.2, 1)) +
    xlim(0,27) + ylim(0,30) + ylab('log(count)') + 
    scale_fill_manual(values = c("Unaffected" = "forestgreen", "Damaged" = "orange", "Destroyed" = "red"))
  
  p2 <- ggplot(df_point %>% filter(impact_true=='Damaged'), aes(fill= impact_sampled, alpha = (impact_sampled==impact_true),y=ifelse(count==0, 0, log(count)), x=event_i)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_alpha_manual(values=c(0.2, 1)) +
    xlim(0,27) + ylim(0,30) + ylab('log(count)') + 
    scale_fill_manual(values = c("Unaffected" = "forestgreen", "Damaged" = "orange", "Destroyed" = "red"))
  
  p3 <- ggplot(df_point %>% filter(impact_true=='Destroyed'), aes(fill= impact_sampled, alpha = (impact_sampled==impact_true),y=ifelse(count==0, 0, log(count)), x=event_i)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_alpha_manual(values=c(0.2, 1)) +
    xlim(0,27) + ylim(0,30) + ylab('log(count)') + 
    scale_fill_manual(values = c("Unaffected" = "forestgreen", "Damaged" = "orange", "Destroyed" = "red"))
  
  ggarrange(p1, p2, p3, 
            labels = c("Unaffected Buildings", "Damaged Buildings", "Destroyed Buildings"),
            ncol = 1, nrow = 3)
  
  ggplot(df_point, aes(fill= paste('Obs:', impact_true, ', Sampled:', impact_sampled), y=ifelse(count==0, 0, log(count)), x=event_i)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_alpha_manual(values=c(0.2, 1)) +
    scale_fill_manual(values = c("Obs: Unaffected , Sampled: Unaffected" = "green", 
                                 "Obs: Unaffected , Sampled: Damaged" = "forestgreen",
                                 "Obs: Unaffected , Sampled: Destroyed" = "darkgreen",
                                 "Obs: Damaged , Sampled: Unaffected" = "cyan",
                                 "Obs: Damaged , Sampled: Damaged" = "blue",
                                 "Obs: Damaged , Sampled: Destroyed" = "darkblue",
                                 "Obs: Destroyed , Sampled: Unaffected" = "yellow",
                                 "Obs: Destroyed , Sampled: Damaged" = "orange",
                                 "Obs: Destroyed , Sampled: Destroyed" = "red"))
  
  
}

compare2events = function(){
  ODDobj1 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Train/EQ20160824ITA_43')
  ODDobj2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Test/EQ20161026ITA_48')
  df_vuln <- data.frame(
    population = integer(),
    haz_max = numeric(),
    event= integer()
  )
  df_vuln %<>% rbind(data.frame(population = ODDobj1$Population,
                        haz_max = apply(ODDobj1@data[,grep("hazMean",names(ODDobj1),value = T)], 1, max, na.rm=T)+0.01,
                        event = 1))
  df_vuln %<>% add_row(population = ODDobj2$Population,
                        haz_max = apply(ODDobj2@data[,grep("hazMean",names(ODDobj1),value = T)], 1, max, na.rm=T),
                        event = 2)
  
  ggplot(df_vuln, aes(x=haz_max, y=log(population), col=event)) + geom_point()
  
}



checkMeans <- function(AlgoResults, s_finish=NULL){
  AlgoResults %<>% addAlgoParams(s_finish)
  M <- 30
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  proposed <- relist(AlgoResults$Omega_sample_phys[particle_min.d[1],,AlgoResults$s_finish], skeleton=Model$skeleton)
  proposed$vuln_coeff$EQFreq=0
  proposed$vuln_coeff$PDens=0
  proposed$vuln_coeff$Vs30=0
  proposed$vuln_coeff$Mag=0
  proposed$vuln_coeff$FirstHaz=0
  proposed$vuln_coeff$Night=0
  proposed$vuln_coeff$FirstHaz.Night=0
  proposed$check$check = 0.5
  impact_sample <- SampleImpact(dir=dir, 
                                Model=Model, 
                                proposed= Omega %>% addTransfParams(), 
                                AlgoParams=list(Np=1, m_CRPS = 30, cores=8, kernel='crps_with_mean', AllParallel=T, NestedCores=1), dat='all')
  
  samples <- rbind(sapply(impact_sample$poly, function(x) x$sampled[1:200]))
  
  k <- 14
  hist(samples[k,])
  abline(v=impact_sample$poly[[1]]$mean[k])
  
  impact_type = 'mortality'
  impact_sample_filt = filter(impact_sample$poly[[5]], impact==impact_type)
  plot(log(impact_sample_filt$mean+10), log(impact_sample_filt$observed+10))
  plot(log(impact_sample$poly[[1]]$sampled[which(impact_sample$poly[[1]]$impact==impact_type)]+10)-log(impact_sample_filt$mean+10), xlim=c(0,100))
  points(log(impact_sample$poly[[2]]$sampled[which(impact_sample$poly[[1]]$impact==impact_type)]+10)-log(impact_sample_filt$mean+10), col='red')
  points(log(impact_sample$poly[[3]]$sampled[which(impact_sample$poly[[1]]$impact==impact_type)]+10)-log(impact_sample_filt$mean+10), col='blue')
  points(log(impact_sample$poly[[4]]$sampled[which(impact_sample$poly[[1]]$impact==impact_type)]+10)-log(impact_sample_filt$mean+10), col='green')
  points(log(impact_sample$poly[[5]]$sampled[which(impact_sample$poly[[1]]$impact==impact_type)]+10)-log(impact_sample_filt$mean+10), col='yellow')
  
}
# 
# dist1_store = c()
# dist2_store = c()
# Np <- 10
# N <- 5
# M <- 5
# for (i in 1:1000){
#   obs <- rnorm(N, 0, 2)
#   dist_m1 <- c()
#   for (np in 1:Np){
#     samp1 <- matrix(rnorm(N*M, 0, 1), NROW=M)
#     means = apply(samp1, 1, mean)
#     sds = apply(samp1, 1, sd)
#     dist_m1 <- c(dist_m1, dnorm(obs, means, sds))
#   }
#   dist_m2 <- c()
#   for (np in 1:Np){
#     samp2 <- rnorm(N, 0, 2)
#     dist_m2 <- c(dist_m2, sum(abs(samp2-obs)))
#   }
#   dist1_store <- c(dist1_store,min(dist_m1))
#   dist2_store <- c(dist2_store,min(dist_m2))
# }
# sum(dist1_store < dist2_store)
# 
# x <- AlgoResults$Omega_sample_phys[,6,1]
# wt <- AlgoResults$W[,1]
# xm <- weighted.mean(x, wt)
# sum(wt * (x - xm)^2)
# plot(x, wt)
# 
# for (i in c(5,6,7,8,16,20)){
#   plot(density(AlgoResults$Omega_sample_phys[,i,1]))
#   lines(density(AlgoResults$Omega_sample_phys[,i,45]), col='red')
# }


compare2Pager <- function(AlgoResults, s_finish=NULL){
  M <- 10
  start_time <- Sys.time()
  AlgoResults %<>% addAlgoParams(s_finish)
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  Omega_min.d <- AlgoResults$Omega_sample_phys[particle_min.d[1],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
  impact_sample <- SampleImpact(dir, Model, Omega_min.d  %>% addTransfParams(), 
                                AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M) 
                                           %>% replace(which(names(AlgoParams)==c('Np')), 1) 
                                           %>% replace(which(names(AlgoParams)==c('SubNat')), F), 
                                output='SampledAgg')
  
  folderin<-paste0(dir,AlgoResults$input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  plot_df <- impact_sample$poly[[1]][, c('observed', 'impact', 'alertlevel_fatalities', 'alertnull_fatalities', 'event_id')]
  for (i in 1:length(impact_sample$poly)){
    plot_df[,paste0('sampled', i)] <- impact_sample$poly[[i]]$sampled
  }
  
  plot_df$max_mmi <- NA
  
  folderin_haz <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/HAZARDobjects_wMaxMMIDiff/'
  ufiles_haz <- na.omit(list.files(path=folderin_haz,pattern=Model$haz,recursive = T,ignore.case = T))
  for(i in 1:NROW(plot_df)){
    file_match <- grep(paste0("_", plot_df$event_id[i], "\\b"),  ufiles_haz, value = TRUE)
    HAZy <- readRDS(paste0(folderin_haz, file_match ))
    which.max.mmi <- which.max(sapply(HAZy[2:length(HAZy)], function(x) max(x$mean)))
    print(which.max.mmi)
    plot_df$max_mmi[i] <- HAZy$hazard_info$max_mmi[which.max.mmi]
  }
  
  plot_df$sampled_min <- apply(plot_df[,grep("sampled", names(plot_df))], 1, min)
  plot_df$sampled_max <- apply(plot_df[,grep("sampled", names(plot_df))], 1, max)
  plot_df$sampled_median <- apply(plot_df[,grep("sampled", names(plot_df))], 1, median)
  plot_df$sampled_mean <- apply(plot_df[,grep("sampled", names(plot_df))], 1, mean)
  plot_df$alertlevel_ODDRIN <- ifelse(plot_df$sampled_median > 1000, 'red', ifelse(plot_df$sampled_median>100, 'orange', ifelse(plot_df$sampled_median>1, 'yellow', 'green')))
  
  ggplot(plot_df %>% filter(impact=='mortality' & alertnull==F), aes(x=max_mmi, y= log(observed+1)-log(sampled_median+1))) + geom_point() + xlab('max_mmi_diff')
  
  ggplot(plot_df %>% filter(impact=='mortality' & alertnull==F), aes(x=observed, y=sampled_median, col=max_mmi, ymin=sampled_min, ymax=sampled_max+1)) + 
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 10^seq(0,5,1)), labels = scales::comma) + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 10^seq(0,5,1)), labels = scales::comma) +
    geom_point() + geom_abline(slope=1, intercept=0) + geom_errorbar() 
  
  ggplot(plot_df %>% filter(impact=='mortality' & alertnull==F), aes(x=observed, y=sampled_median, col=alertlevel_fatalities, ymin=sampled_min, ymax=sampled_max+1)) + 
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 10^seq(0,5,1)), labels = scales::comma) + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 10^seq(0,5,1)), labels = scales::comma) +
    geom_point() + geom_abline(slope=1, intercept=0) + geom_errorbar() + scale_color_manual(values=list('red'='red', 'orange'='orange', 'yellow'='yellow', 'green'='green', 'null'='white'))
  
  ggplot(plot_df %>% filter(impact=='mortality' & alertnull_fatalities==F), aes(x=1:112, y=log(observed+1), col=alertlevel_fatalities, ymin=log(sampled_min+1), ymax=log(sampled_max+1))) + 
    geom_point() + scale_color_manual(values=list('red'='red', 'orange'='orange', 'yellow'='yellow', 'green'='green', 'null'='white'))# + geom_errorbar() 
  
  ggplot(plot_df %>% filter(impact=='mortality' & alertnull==F), aes(x=1:101, y=log(observed+1), col=alertlevel_fatalities, ymin=log(sampled_min+1), ymax=log(sampled_max+1))) + 
    geom_point() + scale_color_manual(values=list('red'='red', 'orange'='orange', 'yellow'='yellow', 'green'='green', 'null'='white'))# + geom_errorbar() 
  
  ggplot(plot_df %>% filter(impact=='mortality' & alertnull_fatalities==F), aes(x=event_id, y=log(observed+1), col=alertlevel_ODDRIN, ymin=log(sampled_min+1), ymax=log(sampled_max+1))) + 
    geom_point(size=2.5) + scale_color_manual(values=list('red'='red', 'orange'='orange', 'yellow'='yellow', 'green'='green', 'null'='white')) + # + geom_errorbar() 
    geom_abline(intercept=log(1000+1), slope=0, col='red') + geom_abline(intercept=log(100+1), slope=0, col='orange') + geom_abline(intercept=log(1+1), slope=0, col='yellow')
  
}

#Plot impact histogram for a single event:

plot_total_impact <- function(ODD_with_impact){
  tot_impact <- colSums(ODD_with_impact@data[,grep('displacement.s', names(ODD_with_impact@data))])
  ggplot(data.frame(value = tot_impact), aes(x = value)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.5, fill = "#002147", color = "black") +
    scale_x_log10(labels = scales::comma) +   # Log10 scale with full number labels
    labs(x = "Posterior Predictive Total Displacement (Morocco Earthquake)", y = "")
}

# Check correlation structure
plot_joint <- function(event_ids=c(16, 23, 31, 67, 68, 70, 89, 94, 124, 125, 135, 139, 164, 170), AlgoParams){
  folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
  
  event_ids_all <- as.numeric(sub(".*_(\\d+)$", "\\1", ufiles))
  #Selecting events:
  # for (file in ufiles){
  #   ODD <- readODD(paste0(folderin, file))
  #   impact_filt <- ODD@impact %>% filter(impact=='mortality')
  #   print(paste('Event:', file, '.Non-zero mortality observations:', sum(impact_filt$observed != 0)))
  # }
  AlgoResults_vs <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-01-10_181501_MCMC_VariogramScore_M100_Npart1000NovAgg5_propCOVmult0.2')
  Omega_i <- which(AlgoResults_vs$Omega_sample_phys[7,] <0.5 & AlgoResults_vs$Omega_sample_phys[8,] >0.5)
  Omega <- AlgoResults_vs$Omega_sample_phys[,Omega_i[Omega_i > 2000][1]] %>% relist(skeleton=Model$skeleton)
  
  AlgoResults_es <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180047_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2')
  Omega_i <- which(AlgoResults_es$Omega_sample_phys[7,] <0.5 & AlgoResults_es$Omega_sample_phys[8,] >0.8)
  Omega <- AlgoResults_es$Omega_sample_phys[,Omega_i[Omega_i > 5000][1]] %>% relist(skeleton=Model$skeleton)
  
  Omega_i <- which(AlgoResults_es$Omega_sample_phys[11,] >0.6)
  Omega <- AlgoResults_es$Omega_sample_phys[,Omega_i[Omega_i > 5000][1]] %>% relist(skeleton=Model$skeleton)
  
  i = 70
  for (i in 70){
    ODD <- readODD(paste0(folderin, ufiles[which(event_ids_all==i)]))
    
    #plot: 
    #grid.arrange(plotODDy_GADM(ODD, 'ISO3C', gadm_level=2, haz_legend=T, var_legend=T, var_discrete=T),plotODDy_GADM(ODD, 'Population', gadm_level=2, haz_legend=T, var_legend=T, var_discrete=F, log_legend=T), ncol=2)
    
    ODD@impact %<>% rbind(data.frame(iso3=c('IRN', 'IRQ'),
                                     sdate= ODD@impact$sdate[1],
                                     impact='mortality',
                                     observed=c(620, 10),
                                     qualifier=NA,
                                     inferred=F,
                                     build_type=NA, 
                                     polygon=c(99, 100)))
    
    sampled_out <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                           replace(which(names(AlgoParams)==c('Np')), 50), 
                         output='SampledAgg')
    
    #ordered_obs <- which(sampled_out[[1]]$impact=='mortality')[order(sampled_out[[1]]$observed[which(sampled_out[[1]]$impact=='mortality')], decreasing=T)]
    ordered_obs <- c(1,2,3, 8)
    par(mfrow=c(2,4))
    polygon_names <- unlist(lapply(ODD@polygons[sampled_out[[1]]$polygon], function(x) x$name))
    polygon_names[which(polygon_names=='n.a. (03), Kermanshah, Iran')] = "Salas Babajani, Kermanshah, Iran" 
    for (obs1 in 1:4){
      for (obs2 in obs1:4){
        if (obs1==obs2 ) next
        obs_plot <- ordered_obs[c(obs1,obs2)]
        obs_sampled <- sampled_out[[1]]$observed[obs_plot]
        for (j in 1:length(sampled_out)){
          obs_sampled %<>% rbind(sampled_out[[j]]$sampled[obs_plot])
        }
        plot(obs_sampled[2:nrow(obs_sampled),], xlim=range(obs_sampled[,1]), ylim=range(obs_sampled[,2]),
             xlab = paste(polygon_names[obs_plot[1]], sampled_out[[j]]$impact[obs_plot[1]]), 
             ylab = paste(polygon_names[obs_plot[2]], sampled_out[[j]]$impact[obs_plot[2]]))
        points(t(obs_sampled[1,]), col='red', pch=3)
      }
    }
    more_obs_plot <- rbind(c(4,7), c(5,6))
    for (obs_plot_i in 1:nrow(more_obs_plot)){
      obs_plot <- more_obs_plot[obs_plot_i,]
      obs_sampled <- sampled_out[[1]]$observed[obs_plot]
      for (j in 1:length(sampled_out)){
        obs_sampled %<>% rbind(sampled_out[[j]]$sampled[obs_plot])
      }
      plot(obs_sampled[2:nrow(obs_sampled),], xlim=range(obs_sampled[,1]), ylim=range(obs_sampled[,2]),
           xlab = paste(polygon_names[obs_plot[1]], sampled_out[[j]]$impact[obs_plot[1]]), 
           ylab = paste(polygon_names[obs_plot[2]], sampled_out[[j]]$impact[obs_plot[2]]))
      points(t(obs_sampled[1,]), col='red', pch=3)
    }
  }
  
  
  for (i in event_ids){
    ODD <- readODD(paste0(folderin, ufiles[which(event_ids_all==i)]))
    
    sampled_out <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                                                    replace(which(names(AlgoParams)==c('Np')), 50), 
          output='SampledAgg')
    
    #ordered_obs <- which(sampled_out[[1]]$impact=='mortality')[order(sampled_out[[1]]$observed[which(sampled_out[[1]]$impact=='mortality')], decreasing=T)]
    ordered_obs <- c(1,2,3,4)
    par(mfrow=c(4,4))
    polygon_names <- unlist(lapply(ODD@polygons[sampled_out[[1]]$polygon], function(x) x$name))
    if (i == 70) polygon_names[4] = "Salas Babajani, Kermanshah, Iran" 
    for (obs1 in 1:4){
      for (obs2 in obs1:4){
        if (obs1==obs2 ) next
        obs_plot <- ordered_obs[c(obs1,obs2)]
        obs_sampled <- sampled_out[[1]]$observed[obs_plot]
        for (j in 1:length(sampled_out)){
          obs_sampled %<>% rbind(sampled_out[[j]]$sampled[obs_plot])
        }
        plot(obs_sampled[2:nrow(obs_sampled),], xlim=range(obs_sampled[,1]), ylim=range(obs_sampled[,2]),
             xlab = paste(polygon_names[obs_plot[1]], sampled_out[[j]]$impact[obs_plot[1]]), 
             ylab = paste(polygon_names[obs_plot[2]], sampled_out[[j]]$impact[obs_plot[2]]))
        points(t(obs_sampled[1,]), col='red', pch=3)
      }
    }
  }
}
