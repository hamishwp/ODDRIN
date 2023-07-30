
results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-07-22_193115'

AlgoResults <- readRDS(results_file)

addAlgoParams <- function(AlgoResults){
  AlgoResults$s_finish <- which(is.na(AlgoResults$Omega_sample_phys[1,1,]))[1]-1
  AlgoResults$n_x <- length(unlist(Model$skeleton))
  AlgoResults$Npart <- NROW(AlgoResults$W)
  return(AlgoResults)
}

plot_density_vs_step = function(AlgoResults, Omega=NULL){
  AlgoResults %<>% addAlgoParams()
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

plot_correlated_posteriors = function(AlgoResults, include_priors=T, Omega=NULL,
                                      pairings=rbind(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12))){
  AlgoResults %<>% addAlgoParams()
  post_samples <- AlgoResults$Omega_sample_phys[,,AlgoResults$s_finish]
  if (include_priors) prior_samples <- AlgoResults$Omega_sample_phys[,,1]
    
  par(mfrow=c(2,3))
  for (p in 1:NROW(pairings)){
    xmin= min(post_samples[,pairings[p,1]]); xmax= max(post_samples[,pairings[p,1]])
    ymin= min(post_samples[,pairings[p,2]]); ymax= max(post_samples[,pairings[p,2]])
    if(include_priors){
      xmin <- min(xmin, prior_samples[,pairings[p,1]]); xmax <- max(xmax, prior_samples[,pairings[p,1]])
      ymin <- min(ymin, prior_samples[,pairings[p,2]]); ymax <- max(ymax, prior_samples[,pairings[p,2]])
    }
    plot(post_samples[,pairings[p,1]], post_samples[,pairings[p,2]], 
         xlab=names(unlist(Model$skeleton))[pairings[p,1]], xlim=c(xmin, xmax),
         ylab=names(unlist(Model$skeleton))[pairings[p,2]], ylim=c(ymin, ymax))
    if (include_priors) points(prior_samples[,pairings[p,1]], prior_samples[,pairings[p,2]], col='blue', alpha=0.5)
    if (!is.null(Omega)){
      points(unlist(Omega)[pairings[p,1]], unlist(Omega)[pairings[p,2]], col='red', pch=4, cex=2, lwd=4)
    }
  }
  
}

sample_post_predictive <- function(AlgoResults, M, s, single_particle=F){
  if (!single_particle){
    sampled_part <- sample(1:AlgoResults$Npart, M, prob=AlgoResults$W[, s])
    impact_sample <- SampleImpact(dir = dir,Model = Model,
                                  proposed = AlgoResults$Omega_sample_phys[sampled_part[1],,s] %>% relist(Model$skeleton) %>% addTransfParams(), 
                                  AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), 1))
    poly_sampled <- impact_sample$poly[[1]][,c('iso3', 'sdate', 'polygon', 'impact', 'observed', 'sampled')]
    point_sampled <- impact_sample$point
    
    for (m in 2:M){
      impact_sample <- SampleImpact(dir = dir,Model = Model,
                                    proposed = AlgoResults$Omega_sample_phys[sampled_part[m],,s] %>% relist(Model$skeleton) %>% addTransfParams(), 
                                    AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), 1))
      poly_sampled <- cbind(poly_sampled, impact_sample$poly[[1]]$sampled)
      point_sampled <- cbind(point_sampled, impact_sample$point)
    }
  } else {
    sampled_part <- sample(1:AlgoResults$Npart, 1, prob=AlgoResults$W[, s])
    impact_sample <- SampleImpact(dir = dir,Model = Model,
                                  proposed = AlgoResults$Omega_sample_phys[sampled_part[1],,s] %>% relist(Model$skeleton) %>% addTransfParams(), 
                                  AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), M))
    poly_sampled <- impact_sample$poly[[1]][,c('iso3', 'sdate', 'polygon', 'impact', 'observed', 'sampled')]
    point_sampled <- impact_sample$point
    
    for (m in 2:M){
      poly_sampled <- cbind(poly_sampled, impact_sample$poly[[m]]$sampled)
    }
  }
  
  df_poly <- poly_sampled
  names(df_poly)[grep("sampled", names(df_poly))] = paste0('sampled.', 1:M)
  df_poly$obs_id[order(df_poly$observed)] <- 1:NROW(df_poly)
  df_poly$crps <- NA
  
  for (i in 1:NROW(df_poly)){
    df_poly$crps[i] <- crps(log(as.numeric(df_poly[i, grep("sampled", names(df_poly))])+10), log(df_poly$observed[i]+10)) * 
      unlist(AlgoParams$kernel_sd)[df_poly$impact[i]]
  }
  
  return(df_poly)
}

compare_CRPS_breakdown_at_different_s <- function(AlgoResults, M){
  AlgoResults %<>% addAlgoParams()
  df_poly1 <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, single_particle=T)
  df_poly2 <- sample_post_predictive(AlgoResults, M, 20, single_particle=T)
  
  ix <- which(df_poly1$impact=='mortality')
  ymin=min(df_poly1$crps[ix], df_poly2$crps[ix])
  ymax=max(df_poly1$crps[ix], df_poly2$crps[ix])
  plot(df_poly1$crps[ix]-df_poly2$crps[ix])
  
  plot(df_poly1$crps[ix], ylim=c(ymin, ymax))
  abline(0,1)
  points(df_poly2$crps[ix], col='red')
  plot(log(df_poly1$observed[ix]+0.1), df_poly1$crps[ix])
  plot(log(df_poly1$observed[ix]+1), log(df_poly1$sampled.1[ix]+1))
  abline(0,1)
  
  impact_type <- 'buildDam'
  df_poly1[which(abs(log(df_poly1$sampled.1+1)-log(df_poly1$observed+1)) >7),] %>% filter(impact==impact_type)
  plot(log(df_poly1$observed[which(df_poly1$impact==impact_type)]+10), log(df_poly1$sampled.2[which(df_poly1$impact==impact_type)]+10), xlab='log(observed+10)', ylab='log(sampled+10)')
  abline(0,1)
  
}



plot_post_predictive = function(AlgoResults, M){
  #mapply(crps, df_plot[grep("^sampled\\.", names(df_sampled))], df_plot$observed)
  AlgoResults %<>% addAlgoParams()
  
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish)
  
  df_long_sampled <-   df_poly %>% pivot_longer(paste0('sampled.', 1:M))
  df_long_sampled$log_observed <- log(df_long_sampled$observed+0.1)
  df_long_sampled$log_sampled <- log(df_long_sampled$value+0.1)
  
  # Create the ggplot
  ggplot(df_long_sampled %>% filter(impact=='mortality', obs_id %in% 1:50), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
    
    geom_density(alpha = 0.5) + geom_vline(aes(xintercept = observed), col='red') + theme_minimal()
  
  ggplot(df_sampled, aes(x=)) + + geom_density()
  
}


plot_correlated_posteriors(AlgoResults, include_priors=T)
plot_density_vs_step(AlgoResults, Omega)
