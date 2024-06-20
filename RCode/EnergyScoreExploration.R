library(scoringRules)

# Univariate with CRPS
mean_true <- 0
sigma_true <- 1
sigma_range <- seq(0.5, 1.5, 0.1)

scoring_assess <- function(M, sigma){
  N_repeats <- 100000
  score_store <- array(NA, N_repeats)
  for (i in 1:N_repeats){
    obs <- rnorm(1, mean_true, sigma_true)
    samples <- rnorm(M, mean_true, sigma)
    score_store[i] <- crps(samples, obs)
  }
  return(c(quantile(score_store, probs=c(0.05, 0.95)),mean(score_store)))
}

M_seq <-  c(5, 15, 30, 45, 60, 75)
sigma_seq <- seq(0.5, 1.5, 0.1)
results_store <- array(NA, dim=c(length(M_seq), length(sigma_seq),3))
for (i in 1:length(M_seq)){
  for (j in 1:length(sigma_seq)){
    results_store[i, j,] <- scoring_assess(M_seq[i], sigma_seq[j]) #calculate 5th quantile and 95th quantile and mean
  }
}
colnames(results_store) <- sigma_seq
rownames(results_store) <- M_seq
#results_store %<>% as.data.frame()

results_long <- merge(merge(as.data.frame(results_store[,,1]) %>% 
                          rownames_to_column("M") %>%
                          gather(Sigma, CRPS.05, -M),
                      as.data.frame(results_store[,,2]) %>% 
                        rownames_to_column("M") %>%
                        gather(Sigma, CRPS.95, -M),
                      by=c('Sigma', 'M')),
                      as.data.frame(results_store[,,3]) %>% 
                        rownames_to_column("M") %>%
                        gather(Sigma, CRPS.mean, -M),  by=c('Sigma', 'M'))
  

cols <- c("5" = "orange", "15" = "red", "30" = "purple", "45" = "blue", "60" = "navy", "75"="black")
results_long$M <- factor(results_long$M, levels=c('5', '15', '30', '45', '60', '75'))

ggplot(results_long, aes(x=Sigma, y=CRPS.mean, group=M, color=M)) + geom_line() + ylab('Mean Continuous Ranked Probability Score')+ 
  scale_colour_manual(values = cols) + geom_point()

# Multivariate with energy score - check correlation
mean_true <- c(0,0)
sigma_true <- c(1, 1)
cor_true <- 0.5
cor_range <- seq(0, 1, 0.1)

scoring_assess <- function(M, cor){
  N_repeats <- 10000
  score_store <- array(NA, N_repeats)
  for (i in 1:N_repeats){
    obs <- rmvnorm(1, mean_true, rbind(c(sigma_true[1], cor_true), c(cor_true, sigma_true[2])))
    samples <- rmvnorm(M, mean_true, rbind(c(sigma_true[1], cor), c(cor, sigma_true[2])))
    score_store[i] <- mmds_sample(c(obs), t(samples))
  }
  return(c(quantile(score_store, probs=c(0.05, 0.95)),mean(score_store)))
}

M_seq <-  c(5, 15, 30, 45, 60, 75)
cor_seq <- seq(0, 1, 0.1)
results_store <- array(NA, dim=c(length(M_seq), length(cor_seq),3))
for (i in 1:length(M_seq)){
  for (j in 1:length(cor_seq)){
    print(paste(i, j))
    results_store[i, j,] <- scoring_assess(M_seq[i], cor_seq[j]) #calculate 5th quantile and 95th quantile and mean
  }
}
colnames(results_store) <- cor_seq
rownames(results_store) <- M_seq
#results_store %<>% as.data.frame()

results_long <- merge(merge(as.data.frame(results_store[,,1]) %>% 
                              rownames_to_column("M") %>%
                              gather(Rho, CRPS.05, -M),
                            as.data.frame(results_store[,,2]) %>% 
                              rownames_to_column("M") %>%
                              gather(Rho, CRPS.95, -M),
                            by=c('Rho', 'M')),
                      as.data.frame(results_store[,,3]) %>% 
                        rownames_to_column("M") %>%
                        gather(Rho, CRPS.mean, -M),  by=c('Rho', 'M'))


cols <- c("5" = "orange", "15" = "red", "30" = "purple", "45" = "blue", "60" = "navy", "75"="black")
results_long$M <- factor(results_long$M, levels=c('5', '15', '30', '45', '60', '75'))

ggplot(results_long, aes(x=Rho, y=CRPS.mean, group=M, color=M)) + geom_line() + ylab('Mean Continuous Ranked Probability Score')+ 
  scale_colour_manual(values = cols) + geom_point()

# Multivariate with energy score - check sd
mean_true <- c(0,0)
sigma_true <- c(1, 1)
cor_true <- 0.55
sigma_range <- seq(0.5, 1.5, 0.1)

scoring_assess <- function(M, sigma){
  N_repeats <- 50000
  score_store <- array(NA, N_repeats)
  for (i in 1:N_repeats){
    obs <- rmvnorm(1, mean_true, rbind(c(sigma_true[1], cor_true), c(cor_true, sigma_true[2])))
    samples <- rmvnorm(M, mean_true, rbind(c(sigma_true[1], cor_true), c(cor_true, sigma)))
    score_store[i] <- es_sample(c(obs), t(samples))
  }
  return(c(quantile(score_store, probs=c(0.05, 0.95)),mean(score_store)))
}

M_seq <-  c(5, 15, 30, 45, 60, 75)
sigma_seq <- seq(0.5, 1.5, 0.1)
results_store <- array(NA, dim=c(length(M_seq), length(sigma_seq),3))
for (i in 1:length(M_seq)){
  for (j in 1:length(sigma_seq)){
    print(paste(i, j))
    results_store[i, j,] <- scoring_assess(M_seq[i], sigma_seq[j]) #calculate 5th quantile and 95th quantile and mean
  }
}
colnames(results_store) <- sigma_seq
rownames(results_store) <- M_seq
#results_store %<>% as.data.frame()

results_long <- merge(merge(as.data.frame(results_store[,,1]) %>% 
                              rownames_to_column("M") %>%
                              gather(Sigma, CRPS.05, -M),
                            as.data.frame(results_store[,,2]) %>% 
                              rownames_to_column("M") %>%
                              gather(Sigma, CRPS.95, -M),
                            by=c('Sigma', 'M')),
                      as.data.frame(results_store[,,3]) %>% 
                        rownames_to_column("M") %>%
                        gather(Sigma, CRPS.mean, -M),  by=c('Sigma', 'M'))


cols <- c("5" = "orange", "15" = "red", "30" = "purple", "45" = "blue", "60" = "navy", "75"="black")
results_long$M <- factor(results_long$M, levels=c('5', '15', '30', '45', '60', '75'))

ggplot(results_long, aes(x=Sigma, y=CRPS.mean, group=M, color=M)) + geom_line() + ylab('Mean Continuous Ranked Probability Score')+ 
  scale_colour_manual(values = cols) + geom_point()


# Multivariate with MMD score with single M
mean_true <- c(0,0)
sigma_true <- c(1, 1)
cor_true <- sqrt(sigma_true[1])/2
cor_range <- seq(0, 1, 0.1)
N_repeats <- 10000

scoring_assess <- function(M, cor){
  score_store <- array(NA, N_repeats)
  obs_store <- array(NA, dim=c(2, N_repeats))
  for (i in 1:N_repeats){
    obs <- exp(rmvnorm(1, mean_true, rbind(c(sigma_true[1], cor_true), c(cor_true, sigma_true[2]))))
    obs_store[,i] <- obs
    #samples <- rmvnorm(M, mean_true, rbind(c(sigma_true[1], cor), c(cor, sigma_true[2])))
    samples <- exp(rmvnorm(M, mean_true, rbind(c(sigma_true[1] , sqrt(sigma_true[1])*cor), c(sqrt(sigma_true[1])*cor, sigma_true[2]))))
    score_store[i] <- vs_sample(c(obs), t(samples))
  }
  return(c(quantile(score_store, probs=c(0.05, 0.95)),mean(score_store)))
}
# plot(obs_store[1,1], obs_store[2,1])
# plot(samples[,1], samples[,2], col='red', pch=19)
#plot(obs_store[1,], obs_store[2,], pch=19)
#points(samples[,1], samples[,2], col='red', pch=19)

M <- 60
cor_seq <- seq(0, 1, 0.1)
results_store <- array(NA, dim=c(length(cor_seq),3))
for (j in 1:length(cor_seq)){
  print(j)
  results_store[j,] <- scoring_assess(M, cor_seq[j]) #calculate 5th quantile and 95th quantile and mean
}

plot(cor_seq, results_store[,3], type='l', xlab='Simulated Data Correlation (true=0.5)', ylab='Mean Energy Score',
     main=paste0('M=60, N=', N_repeats))

# Comparing under similar structure to real data
mean_true <- c(0,0)
sigma_larger_true <- 2
sigma_true <-1
cor_true <- 0.5
cor_range <- seq(0, 1, 0.1)
N_obs <- 100
N_repeats <- 5
dim_per_obs <- 3

scoring_assess <- function(M, cor){
  score_store <- array(NA, dim=c(N_repeats, N_obs))
  
  covar_true <- matrix(cor_true, nrow=dim_per_obs, ncol=dim_per_obs)
  diag(covar_true) <- rep(sigma_true, dim_per_obs)
  covar_true[,1] <-  covar_true[,1] * sqrt(sigma_larger_true)
  covar_true[1,] <-  covar_true[1,] * sqrt(sigma_larger_true)
  covar <- matrix(cor, nrow=dim_per_obs, ncol=dim_per_obs)
  diag(covar) <- rep(sigma_true, dim_per_obs)
  covar[,1] <-  covar[,1] * sqrt(sigma_larger_true)
  covar[1,] <-  covar[1,] * sqrt(sigma_larger_true)
  
  obs <- exp(rmvnorm(N_obs, rep(0, dim_per_obs), covar_true))
  
  for (i in 1:N_repeats){
    #samples <- rmvnorm(M, mean_true, rbind(c(sigma_true[1], cor), c(cor, sigma_true[2])))
    for (j in 1:N_obs){
      samples <- exp(rmvnorm(M, rep(0, dim_per_obs), covar))
      score_store[i,j] <- vs_sample(c(obs[j,]), t(samples))
    }
  }
  print(mean(score_store))
  return(c(quantile(score_store, probs=c(0.05, 0.95)),mean(score_store)))
}
# plot(obs_store[1,1], obs_store[2,1])
# plot(samples[,1], samples[,2], col='red', pch=19)
# plot(obs_store[1,], obs_store[2,], pch=19)
# points(samples[,1], samples[,2], col='red', pch=19)

M <- 60
cor_seq <- seq(0, 1, 0.1)
results_store <- array(NA, dim=c(length(cor_seq),3))
for (j in 1:length(cor_seq)){
  print(j)
  results_store[j,] <- scoring_assess(M, cor_seq[j]) #calculate 5th quantile and 95th quantile and mean
}

plot(cor_seq, results_store[,3], type='l', xlab='Simulated Data Correlation (true=0.5)', ylab='Mean Energy Score',
     main='M=60, 5 repeats over 100 events')

# Multivariate with MMD score with larger sigma
mean_true <- c(0,0)
sigma_true <- c(0.1, 0.1)
cor_true <- 0.5
cor_range <- seq(0, 1, 0.1)

scoring_assess <- function(M, cor){
  N_repeats <- 100000
  score_store <- array(NA, N_repeats)
  for (i in 1:N_repeats){
    obs <- rmvnorm(1, mean_true, rbind(c(sigma_true[1], cor_true), c(cor_true, sigma_true[2])))
    samples <- rmvnorm(M, mean_true, rbind(c(sigma_true[1], cor), c(cor, sigma_true[2])))
    score_store[i] <- mmds_sample(c(obs), t(samples))
  }
  return(c(quantile(score_store, probs=c(0.05, 0.95)),mean(score_store)))
}

M <- 60
cor_seq <- seq(0, 1, 0.1)
results_store <- array(NA, dim=c(length(cor_seq),3))
for (j in 1:length(cor_seq)){
  print(j)
  results_store[j,] <- scoring_assess(M, cor_seq[j]) #calculate 5th quantile and 95th quantile and mean
}

plot(cor_seq, results_store[,3], type='l')

# Dimensionality vs Energy Score
mean_true <- 0
sigma_true <- 1
cor_true <- 0.5
cor_seq <- seq(0,1,0.1)
dim_seq <- seq(1,5,1)

scoring_assess <- function(M, cor, dim){
  N_repeats <- 1000
  score_store <- array(NA, N_repeats)
  
  covar_true <- matrix(cor_true, nrow=dim, ncol=dim)
  diag(covar_true) <- sigma_true
  covar <- matrix(cor, nrow=dim, ncol=dim)
  diag(covar) <- sigma_true
  
  for (i in 1:N_repeats){
    obs <- exp(rmvnorm(1, rep(mean_true, dim), covar_true))
    samples <- exp(rmvnorm(M, rep(mean_true, dim), covar))
    score_store[i] <- vs_sample(c(obs), t(samples))
  }
  return(c(quantile(score_store, probs=c(0.05, 0.95)),mean(score_store)))
}

dim_seq <-  seq(1,5,1)
cor_seq <- seq(0,1,0.1)
results_store <- array(NA, dim=c(length(dim_seq), length(cor_seq),3))
for (i in 1:length(dim_seq)){
  for (j in 1:length(cor_seq)){
    print(paste(i, j))
    results_store[i, j,] <- scoring_assess(60, cor_seq[j], dim_seq[i]) #calculate 5th quantile and 95th quantile and mean
  }
}
colnames(results_store) <- cor_seq
rownames(results_store) <- dim_seq
#results_store %<>% as.data.frame()

results_long <- merge(merge(as.data.frame(results_store[,,1]) %>% 
                              rownames_to_column("Dim") %>%
                              gather(Cor, CRPS.05, -Dim),
                            as.data.frame(results_store[,,2]) %>% 
                              rownames_to_column("Dim") %>%
                              gather(Cor, CRPS.95, -Dim),
                            by=c('Cor', 'Dim')),
                      as.data.frame(results_store[,,3]) %>% 
                        rownames_to_column("Dim") %>%
                        gather(Cor, CRPS.mean, -Dim),  by=c('Cor', 'Dim'))


cols <- c("1" = "orange", "2" = "red", "3" = "purple", "4" = "blue", "5" = "navy")
results_long$Dim <- factor(results_long$Dim, levels=c('1', '2', '3', '4', '5'))

ggplot(results_long, aes(x=Cor, y=CRPS.mean, group=Dim, color=Dim)) + geom_line() + ylab('Mean Continuous Ranked Probability Score')+ 
  scale_colour_manual(values = cols) + geom_point()
