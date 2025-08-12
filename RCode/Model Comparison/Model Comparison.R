library(rstan)
library(arrow)
library(bayesplot)
library(tidyverse)
library(dplyr)
library(magrittr)
library(scales)
library(loo)

stan_model_normalCDF <- stan_model('/home/manderso/Documents/GitHub/ODDRIN/RCode/Model Comparison/Model_NormalCDF.stan')
stan_model_logistic <- stan_model('/home/manderso/Documents/GitHub/ODDRIN/RCode/Model Comparison/Model_Logistic.stan')

stan_model_lognormalCDF <- stan_model('/home/manderso/Documents/GitHub/ODDRIN/RCode/Model Comparison/Model_LogNormalCDF.stan')

event_summaries = read_feather('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/July25Agg/event_summaries')
event_summaries %<>% filter(impact=='mortality')
x = event_summaries[, c('Pop4', 'Pop5', 'Pop6', 'Pop7', 'Pop8', 'Pop9')]
y = event_summaries$observed

x_train = x[1:111,]
x_test = x[112:length(x),]

y_train = y[1:111]
y_test = y[112:length(y)]

stan_data <- list(
  N = nrow(x_train),
  X = x_train[,1:6],
  morts = y_train,
  phi = 100
)

fit_normalCDF <- sampling(
  stan_model_normalCDF,  
  data = stan_data,
  iter = 3000,
  chains = 3,
  warmup = 2000, 
  cores=parallel::detectCores(),
  control = list(max_treedepth=15),
  save_warmup = TRUE,
  seed = 1#,
)

fit_logistic <- sampling(
  stan_model_logistic,  
  data = stan_data,
  iter = 3000,
  chains = 3,
  warmup = 2000, 
  cores=parallel::detectCores(),
  control = list(max_treedepth=15),
  save_warmup = TRUE,
  seed = 1#,
)

fit_lognormal <- sampling(
  stan_model_lognormalCDF,  
  data = stan_data,
  iter = 3000,
  chains = 3,
  warmup = 2000, 
  cores=parallel::detectCores(),
  control = list(max_treedepth=15),
  save_warmup = TRUE, 
  seed = 1#,
)

log_lik1 <- loo::extract_log_lik(fit_normalCDF, merge_chains = FALSE)
log_lik2 <- loo::extract_log_lik(fit_logistic, merge_chains = FALSE)
log_lik3 <- loo::extract_log_lik(fit_lognormal, merge_chains = FALSE)

loo1 <- loo::loo(log_lik1)
loo2 <- loo::loo(log_lik2)
loo3 <- loo::loo(log_lik3)

loo_compare(loo1, loo2, loo3)


mcmc_trace(fit_normalCDF, pars=c('loc_param', 'sd_param', 'event_sd'))
mcmc_trace(fit_logistic, pars=c('loc_param', 'sd_param', 'event_sd'))
mcmc_trace(fit_lognormal, pars=c('loc_param', 'sd_param', 'event_sd'))

posterior = rstan::extract(fit_normalCDF)
posterior_logistic = rstan::extract(fit_logistic)
posterior_lognormal = rstan::extract(fit_lognormal)

N <- nrow(x)
n_draws <- 100
true_pred_samples <- matrix(NA, nrow = N, ncol = n_draws)
true_pred_samples_logistic <- matrix(NA, nrow = N, ncol = n_draws)
true_pred_samples_lognormal <- matrix(NA, nrow = N, ncol = n_draws)

for (i in 1:N) {
  print(i)
  for (j in 1:n_draws) {
    
    event_error_raw <- rnorm(1, 0, 1)
    
    loc <-  posterior$loc_param[j]
    spread <-  posterior$sd_param[j]
    event_sd <- posterior$event_sd[j]

    loc_logistic <-  posterior_logistic$loc_param[j]
    spread_logistic <-  posterior_logistic$sd_param[j]
    event_sd_logistic <- posterior_logistic$event_sd[j]
    
    loc_lognormal <-  posterior_lognormal$loc_param[j]
    spread_lognormal <-  posterior_lognormal$sd_param[j]
    event_sd_lognormal <- posterior_lognormal$event_sd[j]
    
    event_error = event_error_raw * event_sd
    event_error_logistic = event_error_raw * event_sd_logistic
    
    mu = 0
    mu_logistic = 0
    mu_lognormal = 0
    for (I in 4:9){
      score = I-3
      mu = mu + x[i,score] * pnorm(I + event_error, loc, spread) #n(1/(1+exp(-(loc+spread*(I_upd-7))))) 
      mu_logistic = mu_logistic + x[i,score] * (1/(1+exp(-(loc_logistic + spread_logistic *(I + event_error)))))#n(1/(1+exp(-(loc+spread*(I_upd-7))))) 
      mu_lognormal = mu_lognormal + x[i,score] * plnorm(I + event_error, log(loc_lognormal), spread_lognormal)
   }
    log_expected = log(as.numeric(mu_logistic)+10) 
    true_pred_samples[i, j] <-  rpois(1, as.numeric(mu))
    true_pred_samples_logistic[i, j] <-  rpois(1, as.numeric(mu_logistic))#rpois(1, as.numeric(mu))
    true_pred_samples_lognormal[i, j] <-  rpois(1, as.numeric(mu_lognormal))
    #ll[i, j] <- dnorm(log(y[i]+10), log(pmax(1e-6,as.numeric(mu)+10)), 0.05, log=T) #dpois(y[i], as.numeric(mu), log=T)#log(rnbinom(1, size = phi, mu = mu) + 10)
  }
}

# Compute summaries
medians <- apply(true_pred_samples, 1, median)
lower90 <- apply(true_pred_samples, 1, quantile, probs = 0.05)
upper90 <- apply(true_pred_samples, 1, quantile, probs = 0.95)

medians_logistic <- apply(true_pred_samples_logistic , 1, median)
lower90_logistic  <- apply(true_pred_samples_logistic , 1, quantile, probs = 0.05)
upper90_logistic  <- apply(true_pred_samples_logistic , 1, quantile, probs = 0.95)

medians_lognormal <- apply(true_pred_samples_lognormal, 1, median)
lower90_lognormal  <- apply(true_pred_samples_lognormal, 1, quantile, probs = 0.05)
upper90_lognormal  <- apply(true_pred_samples_lognormal, 1, quantile, probs = 0.95)

# Plot with uncertainty bands
plot_df <- data.frame(
  obs = y,
  median = medians,
  lower90 = lower90,
  upper90 = upper90,
  median_logistic = medians_logistic,
  lower90_logistic = lower90_logistic,
  upper90_logistic = upper90_logistic,
  median_lognormal = medians_lognormal,
  lower90_lognormal = lower90_lognormal,
  upper90_lognormal = upper90_lognormal
)

plot_df$train_flag = 'TEST'
plot_df$train_flag[1:nrow(x_train)] = 'TRAIN'

library(scales)
ggplot(plot_df %>% filter(train_flag=='TRAIN')) +
  geom_point(aes(x = obs, y = median), col='black') +
  #geom_point(aes(x = obs, y = median_logistic), col='blue') +
  #geom_point(aes(x = obs, y = median_lognormal), col='red') +
  # geom_errorbar(aes(ymin = lower90, ymax = upper90),
  #               width = 0, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = comma_format(),
    minor_breaks = NULL, 
    limits=c(NA, 100000)
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = comma_format(),
    minor_breaks = NULL, 
    limits=c(NA, 100000)
  ) +
  labs(
    x = "Observed",
    y = "Predicted (median)",
    title = "Posterior Predictive Summary"
  ) +
  theme_minimal()

df_train = plot_df %>% filter(train_flag=='TRAIN')
mean((log(df_train$obs+10) - log(df_train$median+10))^2)
mean((log(df_train$obs+10) - log(df_train$median_logistic+10))^2)
mean((log(df_train$obs+10) - log(df_train$median_lognormal+10))^2)



#-------------------------------------------------------
#---------- COMPARE USING OPTIMISATION -----------------
#-------------------------------------------------------

pred_normalCDF = function(x, loc_param, sigma_param){
  probs = pnorm(seq(4,9,1), loc_param, sigma_param)
  impact = rowSums( sweep(x, 2, probs, `*`))
}

fit_mse <- function(params, x_train, y_train) {
  loc_param <- params[1]
  sigma_param <- params[2]
  
  y_pred <- pred_normalCDF(x_train, loc_param, sigma_param)
  
  mse <- mean((log(y_pred+10) - log(y_train+10))^2)
  return(mse)
}

res <- optim(
  par = c(loc_param = 13, sigma_param = 1),  # starting values
  fn = fit_mse,
  x_train = x_train,
  y_train = y_train,
  method = "L-BFGS-B",
  lower = c(9, 0.01),  # lower bounds for loc and sigma
  upper = c(20, 5)      # upper bounds
)

best_loc <- res$par[1]
best_sigma <- res$par[2]
cat("Optimised loc =", best_loc, " sigma =", best_sigma, "\n")

y_pred_best = pred_normalCDF(x_train, best_loc, best_sigma)
plot_best_fit(y_pred_best, y_train)

pred_lognormalCDF = function(x, loc_param, sigma_param){
  probs = plnorm(seq(4,9,1), loc_param, sigma_param)
  impact = rowSums( sweep(x, 2, probs, `*`))
}

fit_mse_lnorm <- function(params, x_train, y_train) {
  loc_param <- params[1]
  sigma_param <- params[2]
  
  y_pred <- pred_lognormalCDF(x_train, loc_param, sigma_param)
  
  #plot(log(y_pred+10), log(y_train+10))
  
  mse <- mean((log(y_pred+10) - log(y_train+10))^2)
  return(mse)
}

res <- optim(
  par = c(loc_param = 2.7, sigma_param = 0.2),  # starting values
  fn = fit_mse_lnorm,
  x_train = x_train,
  y_train = y_train,
  method = "L-BFGS-B",
  lower = c(1, 0.01),  # lower bounds for loc and sigma
  upper = c(3, 0.3)      # upper bounds
)

best_loc <- res$par[1]
best_sigma <- res$par[2]
cat("Optimised loc =", best_loc, " sigma =", best_sigma, "\n")

y_pred_best_ln = pred_lognormalCDF(x_train, best_loc, best_sigma)
plot(log(y_pred_best+10), log(y_train+10))



pred_logistic = function(x, loc_param, sigma_param){
  probs = 1/(1+exp(-(loc_param + seq(4,9,1)*sigma_param)))
  impact = rowSums( sweep(x, 2, probs, `*`))
}

fit_mse_logistic <- function(params, x_train, y_train) {
  loc_param <- params[1]
  sigma_param <- params[2]
  
  y_pred <- pred_logistic(x_train, loc_param, sigma_param)
  
  #plot(log(y_pred+10), log(y_train+10))
  
  mse <- mean((log(y_pred+10) - log(y_train+10))^2)
  return(mse)
}

res <- optim(
  par = c(loc_param = -15, sigma_param = 2),  # starting values
  fn = fit_mse_logistic,
  x_train = x_train,
  y_train = y_train,
  method = "L-BFGS-B",
  lower = c(-30, 0.1),  # lower bounds for loc and sigma
  upper = c(-5, 3.5)      # upper bounds
)

best_loc <- res$par[1]
best_sigma <- res$par[2]
cat("Optimised loc =", best_loc, " sigma =", best_sigma, "\n")

y_pred_best_logistic = pred_logistic(x_train, best_loc, best_sigma)

ggplot() +
  geom_point(data = data.frame(y_pred=y_pred_best, y_obs=y_train), aes(y=y_pred, x=y_train), size = 2) +
  geom_point(data = data.frame(y_pred=y_pred_best_ln, y_obs=y_train), aes(y=y_pred, x=y_train), size = 2, col='red') +
  geom_point(data = data.frame(y_pred=y_pred_best_logistic, y_obs=y_train), aes(y=y_pred, x=y_train), size = 2, col='blue') +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = comma_format(),
    minor_breaks = NULL, 
    limits=c(NA, 5000)
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = comma_format(),
    minor_breaks = NULL, 
    limits=c(NA, 5000)
  ) +
  labs(
    x = "Observed",
    y = "Predicted",
    title = "Posterior Predictive Summary"
  ) +
  theme_minimal()

mean((log(y_pred_best+10) - log(y_train+10))^2)
mean((log(y_pred_best_ln+10) - log(y_train+10))^2)
mean((log(y_pred_best_logistic+10) - log(y_train+10))^2)






