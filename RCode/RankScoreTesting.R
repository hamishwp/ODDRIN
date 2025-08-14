
flattenImpactSample <- function(impact_sample){
  df <- data.frame(event_id = impact_sample$poly[[1]]$event_id,
                   iso3 = impact_sample$poly[[1]]$iso3,
                   polygon = impact_sample$poly[[1]]$polygon,
                   impact = impact_sample$poly[[1]]$impact,
                   observed = impact_sample$poly[[1]]$observed)#,
                   #n_pixels = impact_sample$poly[[1]]$n_pixels,
                   #max_intensity = impact_sample$poly[[1]]$max_intensity)
  for (j in 1:length(impact_sample$poly)){
    df %<>% cbind(impact_sample$poly[[j]]$sampled)
  }
  colnames(df)[grep('sampled',names(df))] <- paste0('sampled.', 1:length(impact_sample$poly))
  return(df)
}

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




#------------------------------------------------------------------------------------------
#--------------------------------Why is ABC failing?---------------------------------------
#------------------------------------------------------------------------------------------

n_obs_all = c(10, 50, 100, 150)
n_rep = 100
n_rep2 = 1000000
mean_df = array(Inf, dim=c(2, length(n_obs_all), n_rep))

for (i in seq_along(n_obs_all)){
  n_obs = n_obs_all[i]
  print(n_obs)
  for (r in 1:n_rep){
    obs = rnorm(n_obs, 0, 1)
    for (j in 1:n_rep2){
      y_true = rnorm(n_obs, 0, 1)
      y_missp = rnorm(n_obs, 0, 0)
      mean_df[1,i, r] = min(mean_df[1,i, r], sqrt(sum((obs-y_true)^2)))
      mean_df[2,i, r] = min(mean_df[2,i, r], sqrt(sum((obs-y_missp)^2)))
    }
  }
}
n_obs_i = 1
plot(mean_df[1,n_obs_i,], ylim=range(mean_df[1,n_obs_i,], mean_df[2,n_obs_i,]))
points(mean_df[2,n_obs_i,], col='red')

#saveRDS(mean_df, 'dists_toy_NObsIncreasing16thJuly')
mean_df = readRDS('dists_toy_NObsIncreasing16thJuly')

data_list <- list()
for (i in seq_along(n_obs_all)) {
  for (r in 1:n_rep) {
    data_list[[length(data_list) + 1]] <- data.frame(
      n_obs = n_obs_all[i],
      value = mean_df[1, i, r],
      type = "y_true"
    )
    data_list[[length(data_list) + 1]] <- data.frame(
      n_obs = n_obs_all[i],
      value = mean_df[2, i, r],
      type = "y_missp"
    )
  }
}
plot_df <- do.call(rbind, data_list)

plot_df$type <- factor(plot_df$type, levels = c("y_missp", "y_true"))

# Set numeric labels for the x-axis
x_labels <- c("y_missp" = "0.01", "y_true" = "1")

# Create individual plots for each unique n_obs
plots <- list()
n_obs_unique <- unique(plot_df$n_obs)

for (i in seq_along(n_obs_unique)) {
  sub_df <- subset(plot_df, n_obs == n_obs_unique[i])
  p <- ggplot(sub_df, aes(x = type, y = value, fill = type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = c("y_missp" = "#1FA187", "y_true" = "#440154")) +
    scale_x_discrete(labels = x_labels) +
    labs(
      x = expression(sigma[sim]),
      y = if (i == 1) expression("Min " * lambda(x, y) * " over " * 10^6 * " simulations") else NULL
    ) +
    ggtitle(bquote("("*.(letters[i])*")"~N[obs]~"="~.(n_obs_unique[i]))) +
    theme_minimal(base_family = "Times") +
    theme(
      plot.title = element_text(hjust = -0.25, face = "bold", size = 15),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size=13.5),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text( size = 13),
      axis.text.y = element_text( size = 13),
      legend.position = "none"
    )
  plots[[i]] <- p
}

# Arrange all plots in one row
gridExtra::grid.arrange(grobs = plots, nrow = 1)
#ABCbias.pdf, 10.5 x 3

#MCMC
n_obs = 150
n_rep = 100
n_rep2 = 40000
mean_df = array(Inf, dim=c(2, length(n_obs_all), n_rep))
y_true_store = array(0, dim=c(n_rep2, n_obs))
sigma = 0.5
sigma_store = c()
mean_store = c()
i = 1
for (r in 1:1){
  obs = rnorm(n_obs, 0, 1)
  u = rnorm(n_obs, 0, 1)
  y_true = u * sigma
  y_missp = rnorm(n_obs, 0, 0)
  ll_old = -200* max(abs(obs-y_true))
  for (j in 1:n_rep2){
    sigma_prop = max(0, rnorm(1, sigma, 0.05))
    u_prop = 0.999 * u + sqrt(1-0.999^2) * rnorm(n_obs, 0, 1)
    y_true_prop = u_prop * sigma_prop
    ll = -200*max(abs(obs-y_true_prop))
    if (runif(1) < exp(ll - ll_old)){
      u = u_prop
      sigma = sigma_prop
      y_true = y_true_prop
      ll_old = ll
    } 
    y_true_store[j,] = y_true
    mean_df[1,i, r] = min(mean_df[1,i, r], mean(abs(obs-y_true)))
    mean_store = c(mean_store, mean(abs(obs-y_true)))
    mean_df[2,i, r] = min(mean_df[2,i, r], mean(abs(obs-y_missp)))
    sigma_store = c(sigma_store, sigma)
    # y_missp = y_missp 
    # mean_df[1,i, r] = min(mean_df[1,i, r], mean(abs(obs-y_true)))
    # mean_df[2,i, r] = min(mean_df[2,i, r], mean(abs(obs-y_missp)))
  }
}

par(mfrow=c(2,1))
plot(sigma_store)
abline(h=1)

plot(mean_store)
abline(h=mean_df[2,i, r], col='red')
par(mfrow=c(1,1))

mean_df[1,1,1]
mean_df[2,1,1]

plot(mean_df[1,1,], ylim=range(mean_df[1,1,], mean_df[2,1,]))
points(mean_df[2,1,], col='red')

plot(y_true_store[, 2], type='l')
abline(h = obs[2], col='red')



N <- 5            # Dimensionality
n_samples <- 100000000        # Number of samples

# Function to simulate scaled chi distribution
simulate_chi_distances <- function(N, sigma, n_samples) {
  scaling <- sqrt(sigma^2 + 1)
  chi_samples <- sqrt(rchisq(n_samples, df = N))
  return(scaling * chi_samples)
}

# Simulate distances
dist_correct <- simulate_chi_distances(N, sigma = 1, n_samples)
dist_missp   <- simulate_chi_distances(N, sigma = 0.01, n_samples)

# Combine and plot
library(ggplot2)
df <- data.frame(
  distance = c(dist_correct, dist_missp),
  model = factor(rep(c("sigma = 1", "sigma = 0.01"), each = n_samples))
)

ggplot(df, aes(x = distance, fill = model)) +
  geom_density(alpha = 0.6) +
  labs(
    title = paste("Euclidean Distance Distributions (Chi, N =", N, ")"),
    x = "Distance", y = "Density", fill = "Model"
  ) +
  theme_minimal()

#------------------------------------------------------------------------------------------
#------------------- ES vs VS vs MMD vs MVR Rank Histogram --------------------------------
#------------------------------------------------------------------------------------------

get_mmds_rank_single <- function(mat){
  ranks = c()
  for (j in 1:NCOL(mat)){
    ranks = c(ranks, mmds_sample(mat[,j], mat[,-j]))
  }
  return(rank(ranks, ties.method='random')[1])
}

get_banddepth_rank_single <- function(mat, noise = NA){
  m <- NCOL(mat)
  d <- NROW(mat)
  pre_ranks <- c()
  for (j in 1:m){
    sum <- 0
    for (k in 1:d){
      rank_k <- sum(mat[k,j] <= mat[k,])
      sum <- sum + rank_k * (m - rank_k) + (rank_k-1) * sum(mat[d,j]==mat[k,])
    }
    pre_ranks <- c(pre_ranks, sum /d)
  }
  if (all(!is.na(noise))){
    pre_ranks = pre_ranks + noise / d
  }
  #pre_ranks = MBD(t(mat), manage_ties=T)
  return(rank(pre_ranks,  ties.method ='random')[1])
}

get_average_rank_single <- function(mat){
  pre_ranks <- apply(apply(mat, 1, rank), 1, mean)
  return(rank(pre_ranks, ties.method='random')[1])
}

get_multivariate_rank_single <- function(mat){
  pre_ranks <- c()
  for (j in 1:NCOL(mat)){
    pre_ranks <- c(pre_ranks, sum(colSums(mat[,j]<= mat[,-j])==NROW(mat)))
  }
  return(rank(pre_ranks, ties.method='random')[1])
}

library(vegan)
get_mst_rank_single <- function(mat){
  pre_ranks <- c()
  for (j in 1:NCOL(mat)){
    pre_ranks <- c(pre_ranks, sum(spantree(dist(t(mat[,-j])))$dist))
  }
  return(rank(pre_ranks, ties.method='random')[1])
}

standardise_and_test_ranks <- function(pre_ranks, n_sim){
  ranks_std <- (pre_ranks-runif(length(pre_ranks),0,1))/(n_sim + 1)
  test_stat <- AndersonDarlingTest(ranks_std, null='punif')$statistic
  return(test_stat)
}

n_rep = 20
gen_cor_all = c(0.05,0.5,0.95)
test_cor_all = c(0.05,0.5,0.95)
results_store = array(0, dim=c(8, length(gen_cor_all), length(test_cor_all), n_rep))
n_obs = 100
ndim = 2

for (gen_i in seq_along(gen_cor_all)){
  #loop through correlations in observation generation
  gen_cor =  gen_cor_all[gen_i]
  print(gen_cor)
  for (test_i in seq_along(test_cor_all)){
    #loop through correlations in simulation generation
    test_cor = test_cor_all[test_i]
    print(test_cor)
    
    for (rep in 1:n_rep){
      print(rep)
      
      #scoring rules:
      es_stores = c()
      vs_stores = c()
      mmds_stores = c()
      
      #set up pre_rank arrays:
      pre_ranks_mmd = c()
      pre_ranks_mst = c()
      pre_ranks_bd = c()
      pre_ranks_av = c()
      pre_ranks_mvr = c()
      
      for (obs_i in 1:n_obs){
        
        
        # tests=cbind(x1, x2)
        # testsother = cbind(x3, x4)
        # 
        # bds_single = c()
        # for (i in 1:NROW(testsother)){
        #   bds_single = c(bds_single, get_banddepth_rank_single(cbind(testsother[i,], t(tests))))
        # }
        # 
        # color_func <- colorRampPalette(c("yellow", "red"))
        # n_colors <- 100
        # colors <- color_func(n_colors)
        # bd_scaled <- as.numeric(cut(bds_single, breaks = n_colors))
        # point_colors <- colors[bd_scaled]
        # 
        # plot(tests, xlim=range(c(tests[,1], testsother[,1])), ylim=range(tests[,2], testsother[,2]))
        # points(testsother, col= point_colors)
        
        # df <- 5  # degrees of freedom for t-distribution
        # test_cor <- 0.1
        # gen_cor <- 0.9
        # 
        # # t-distribution variance:
        # t_var <- df / (df - 2)  # for df > 2
        # scale_factor <- 3*sqrt(1 / t_var)  # to rescale to variance 1
        # 
        # # First banana
        # tests1 <- rt(100, df) * scale_factor
        # tests2 <- rt(100, df) * scale_factor + 2 * test_cor * (tests1^2 - 1)
        # tests2 <- tests2 / sqrt(1 + 2 * (2 * test_cor)^2)
        # tests <- cbind(tests1, tests2)
        # 
        # # Second banana
        # tests3 <- rt(100, df) * scale_factor
        # tests4 <- rt(100, df) * scale_factor + 2 * gen_cor * (tests3^2 - 1)
        # tests4 <- tests4 / sqrt(1 + 2 * (gen_cor)^2)
        # testsother <- cbind(tests3, tests4)
        # 
        # # Plot
        # plot(tests, xlab = "X", ylab = "Y", main = "Banana Distributions")
        # points(testsother, col = 'red')
        # 
        # tests1 = rnorm(100) * sqrt(4)
        # tests2 = rnorm(100) +  2*test_cor*(tests1^2-1)
        # tests2 = tests2/sqrt(1+2*(2*test_cor)^2) * sqrt(4)
        # tests = cbind(tests1, tests2)
        # 
        # tests3 = rnorm(100) * sqrt(4)
        # tests4 = rnorm(100) + 2*(gen_cor)*(tests3^2-1)
        # tests4 = tests4/sqrt(1+2*(gen_cor)^2) * sqrt(4)
        # testsother = cbind(tests3, tests4)
        # 
        # plot(tests)
        # points(testsother, col='red')
        
        obs = t(rmvnorm(1,c(2,0.5), cbind(c(1,test_cor),c(test_cor,1))))
        obs[obs < 0] = 0
        
        tests = rmvnorm(100,c(2,0.5), cbind(c(1,gen_cor),c(gen_cor,1)))
        tests[tests < 0] = 0
        
        # obs1 <- rnorm(1) * 2  # Var = 4
        # obs2 <- (rnorm(1) + 2 * test_cor * (obs1^2 - 1)) / sqrt(1 + 2 * (2 * test_cor)^2) * 2  # Var = 4
        # obs <- c(obs1, obs2)
        # 
        # # Tests
        # tests1 <- rnorm(100) * 2
        # tests2 <- (rnorm(100) + 2 * gen_cor * (tests1^2 - 1)) / sqrt(1 + 2 * (2 * gen_cor)^2) * 2
        # tests <- cbind(tests1, tests2)
        
        # df <- 5  # degrees of freedom for t-distribution
        # t_var <- df / (df - 2)  # variance of t-distribution for df > 2
        # scale_factor <- sqrt(1 / t_var)  # to rescale to variance 1
        # 
        # obs1 <- rt(1, df) * scale_factor
        # obs2 <- rt(1, df) * scale_factor + 2 * test_cor * (obs1^2 - 1)
        # obs2 <- obs2 / sqrt(1 + 2 * (test_cor)^2)
        # obs <- c(obs1, obs2)
        # 
        # # Test distribution
        # tests1 <- rt(100, df) * scale_factor
        # tests2 <- rt(100, df) * scale_factor + 2 * gen_cor * (tests1^2 - 1)
        # tests2 <- tests2 / sqrt(1 + 2 * (gen_cor)^2)
        # tests <- cbind(tests1, tests2)
        
        # 
        # tests3 = rnorm(100)
        # tests4 = rnorm(100) + 2*(5*gen_cor)*(tests3^2-1)
        # tests4 = tests4/sqrt(1+2*(5*gen_cor)^2)
        # testsother = cbind(tests3, tests4)

        # 
        # bds_single = c()
        # for (i in 1:NROW(testsother)){
        #   bds_single = c(bds_single, get_banddepth_rank_single(cbind(testsother[i,], t(tests))))
        # }
        # 
        # color_func <- colorRampPalette(c("yellow", "red"))
        # n_colors <- 100
        # colors <- color_func(n_colors)
        # bd_scaled <- as.numeric(cut(bds_single, breaks = n_colors))
        # point_colors <- colors[bd_scaled]
        # 
        # plot(tests, xlim=range(c(tests[,1], testsother[,1])), ylim=range(tests[,2], testsother[,2]))
        # points(testsother, col= point_colors)
        
        # 
        # plot(tests)
        # points(testsother, col='red')
        
        # obs1 = rnorm(1)
        # obs2 = rnorm(1) + 2*gen_cor*(obs1^2-1)
        # obs2= obs2/sqrt(1+2*(2*gen_cor)^2)
        # obs = c(obs1, obs2)
        #tests = mvtnorm::rmvt(100, sigma = (1-test_cor)* diag(n_dim) + test_cor * ones(n_dim), df=3)
        #obs = mvtnorm::rmvt(1, sigma = (1-gen_cor)* diag(n_dim) + gen_cor * ones(n_dim), df=3)
        # N = sample(1000:1000000, 1)
        # test_error = rmvnorm(100,rep(0, n_dim), (1-test_cor)* diag(n_dim) + test_cor * ones(n_dim))
        # test_probs = t(apply(test_error, 1, function(x) pnorm(x, c(2,4))))
        # test_probs[,1] = pmax(test_probs[,1] - test_probs[,2], 0)
        # test_probs =  cbind(test_probs, 1-rowSums(test_probs))
        # tests = t(log(apply(test_probs, 1, function(x) rmultinom(1, N, x))+10)[1:2,])
        # 
        # obs_error = rmvnorm(1, rep(0, n_dim), (1-gen_cor)* diag(n_dim) + gen_cor * ones(n_dim))
        # obs_probs = pnorm(obs_error, c(2,4))
        # obs_probs[,1] = pmax(obs_probs[,1] - obs_probs[,2], 0)
        # obs_probs =  c(obs_probs, 1-sum(obs_probs))
        # obs = as.matrix(log(rmultinom(1, N, obs_probs) + 10)[1:2], nrow=1)
        
        es_stores = c(es_stores,  es_sample(c(obs), t(tests)))
        vs_stores = c(vs_stores,  vs_sample(c(obs), t(tests)))
        mmds_stores = c(mmds_stores,  mmds_sample(c(obs), t(tests)))
        dat= cbind(obs, t(tests))
        
        pre_ranks_mmd = c(pre_ranks_mmd, get_mmds_rank_single(dat))
        pre_ranks_mst = c(pre_ranks_mst, get_mst_rank_single(dat))
        pre_ranks_bd = c(pre_ranks_bd, get_banddepth_rank_single(dat))
        pre_ranks_av = c(pre_ranks_av, get_average_rank_single(dat))
        pre_ranks_mvr = c(pre_ranks_mvr, get_multivariate_rank_single(dat))
        
      }
      
      results_store[1, gen_i, test_i, rep] = mean(es_stores)
      results_store[2, gen_i, test_i, rep] = mean(vs_stores)
      results_store[3, gen_i, test_i, rep] = mean(mmds_stores)
      
      results_store[4, gen_i, test_i, rep] = standardise_and_test_ranks(pre_ranks_mmd, 100)
      results_store[5, gen_i, test_i, rep] = standardise_and_test_ranks(pre_ranks_mst, 100)
      results_store[6, gen_i, test_i, rep] = standardise_and_test_ranks(pre_ranks_bd, 100)
      results_store[7, gen_i, test_i, rep] = standardise_and_test_ranks(pre_ranks_av, 100)
      results_store[8, gen_i, test_i, rep] = standardise_and_test_ranks(pre_ranks_mvr, 100)
    }
  }
}

par(mfrow = c(2, 5), mar = c(4, 4, 2, 1))
titles = c('ES', 'VS', 'MMDS', 'MMD Rank', 'MST Rank', 'BD Rank', 'AVE Rank', 'MVR Rank')
for (i in 1:8){
  plot(rep(gen_cor_all,n_rep), c(results_store[i,,1,]), col='red', main=titles[i], ylim=range(results_store[i,,,]))
  points(rep(gen_cor_all,n_rep), results_store[i,,2,], col='orange')
  points(rep(gen_cor_all, n_rep), results_store[i,,3,], col='yellow')
}

par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
titles <- c('(b) Energy Score', '(c) Variogram Score', '(d) MMD Kernel Score', 'MMD Rank', 'MST Rank', 'BD Rank', 'AVE Rank', ' (e) MVR Rank Hist. Test')

# Colors for 3 methods
layout_matrix <- matrix(c(1, 2, 3,
                          1, 4, 5), ncol = 3, byrow = TRUE)

layout(layout_matrix, widths = c(1.5, 1))

# Adjust margins if needed:
par(mar = c(4, 4, 2, 1), family="Times")


cor_gens <- c(0.05, 0.5, 0.95)
n <- 100

tests_list <- lapply(cor_gens, function(cor) {
  
  dat = rmvnorm(n,c(2,0), cbind(c(1,cor),c(cor,1)))
  dat[dat < 0] = 0
  
  return(dat)
})

# Combine for limits
all_tests <- do.call(rbind, tests_list)
xlim <- range(all_tests[, 1]) + c(-0.5, 0.5)
ylim <- range(all_tests[, 2]) + c(-0.5, 0.5)

# Plot
plot(tests_list[[1]], col = '#440154', pch = 19,
     xlim = xlim, ylim = ylim,
     xlab = expression(x[1]), ylab = expression(x[2]), cex.lab = 1.5)

points(tests_list[[2]], col = '#FDE725FF', pch = 19)
points(tests_list[[3]], col = '#1FA187', pch = 19)

legend("topleft",
       legend = c(
         expression(rho == 0.05),
         expression(rho == 0.5),
         expression(rho == 0.95)
       ),
       col = c('#440154', '#FDE725FF', '#1FA187'),
       pch = 19,
       bty = "n",
       cex = 1.5,          
       y.intersp = 1.5)
mtext('(a) Simulated Data Under Varying Correlations', side = 3, line = 0.5, adj = 0, cex = 0.9, font = 1)

method_colors <- c('#440154', '#FDE725FF', '#1FA187')
method_labels <- c('Method 1', 'Method 2', 'Method 3')

for (i in c(1:3, 8)) {
  box_data <- list()
  box_names <- list()
  for (j in 1:3) {
    for (k in 1:3) {
      if (i <= 3) {
        box_data[[length(box_data) + 1]] <- results_store[i, j, k, ]
      } else {
        box_data[[length(box_data) + 1]] <- results_store[i, j, k, ] * 0.05
      }
      box_names[[length(box_names) + 1]] <- substitute(rho[Sim] == x, list(x = test_cor_all[k]))
    }
  }
  
  group_labels <- lapply(gen_cor_all, function(x) substitute(rho[Obs] == y, list(y = x)))
  
  ylims <- range(results_store[i,,,], na.rm = TRUE)
  if (i > 3) ylims <- ylims * 0.05
  
  # Add 10% vertical space on top
  y_margin <- 0.1 * diff(ylims)
  ylims[2] <- ylims[2] + y_margin
  
  bp <- boxplot(box_data, 
                col = rep(method_colors, times = 3),
                xaxt = "n",
                main = "", #titles[i],
                ylim = ylims,
                ylab = "",
                cex.axis = 1)
  mtext(titles[i], side = 3, line = 0.5, adj = 0, cex = 0.9, font = 1)
  
  # Tick labels for each box
  tick_labels <- rep(test_cor_all, times = 3)
  axis(1, at = 1:9, labels = tick_labels, las = 1, cex.axis = 1)
  
  # Add x-axis label
  mtext(expression(rho[Sim]), side = 1, line = 2.5, cex = 1)
  
  # Add rho[Obs] inside plot above groups
  group_positions <- c(mean(1:3), mean(4:6), mean(7:9))
  usr <- par("usr")
  y_pos <- usr[4] - 0.07 * diff(usr[3:4])  # inside plot, near top
  
  text(x = group_positions, y = y_pos, labels = do.call(expression, group_labels), cex = 1.25)
  
  # Optional: Add group separators
  abline(v = c(3.5, 6.5), lty = 2, col = "grey")
}

# for (i in c(1:3, 8)) {
#   # Collect boxplot data
#   box_data <- list()
#   box_names <- list()  # use list for expression labels
#   for (j in 1:3) {  # for each gen_cor value
#     for (k in 1:3) {  # for each method
#       if (i <= 3) {
#         box_data[[length(box_data) + 1]] <- results_store[i, j, k, ]
#       } else {
#         box_data[[length(box_data) + 1]] <- results_store[i, j, k, ] * 0.05
#       }
#       box_names[[length(box_names) + 1]] <- substitute(rho[Sim] == x, list(x = test_cor_all[k]))
#     }
#   }
#   
#   # Group labels as expressions too
#   group_labels <- lapply(gen_cor_all, function(x) substitute(rho[Obs] == y, list(y = x)))
#   
#   # Plot without x-axis labels
#   ylims <- range(results_store[i,,,], na.rm = TRUE)
#   if (i > 3) ylims <- ylims * 0.05
#   
#   bp <- boxplot(box_data, 
#                 col = rep(method_colors, times = 3),
#                 xaxt = "n",  # no x-axis labels
#                 main = titles[i],
#                 ylim = ylims,
#                 ylab = "",
#                 cex.axis = 1)
#   
#   # Add box labels (expressions)
#   axis(1, at = 1:9, labels = do.call(expression,box_names), las = 2, cex.axis = 1)
#   
#   # Add group labels (expressions)
#   group_positions <- c(mean(1:3), mean(4:6), mean(7:9))
#   axis(1, at = group_positions, labels = do.call(expression,group_labels), tick = FALSE, line = 4.25, cex.axis = 1.2)
#   
#   # Optional: Add group separators
#   abline(v = c(3.5, 6.5), lty = 2, col = "grey")
# }


#------------------------------------------------------------------------------------------
#-------------------Rank Histogram Comparison Trace Plots --------------------------------
#------------------------------------------------------------------------------------------

BD <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-09_151133.450368__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_RankBD')
MST <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-09_090443.45158__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_RankMST_backup')
AVE <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-09_091042.264797__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_RankAVE')

par(mfrow=c(2,2))
plot(BD$Omega_sample_phys[1,], type='l'); lines(MST$Omega_sample_phys[1,], type='l', col='red'); 
plot(BD$Omega_sample_phys[2,], type='l'); lines(MST$Omega_sample_phys[2,], type='l', col='red'); 
plot(BD$Omega_sample_phys[3,], type='l'); lines(MST$Omega_sample_phys[3,], type='l', col='red'); 
plot(BD$Omega_sample_phys[4,], type='l'); lines(MST$Omega_sample_phys[4,], type='l', col='red'); 

#------------------------------------------------------------------------------------------


# Create a plot that demonstrates rank histogram

# 1. Set up side-by-side plots
par(mfrow = c(1, 2))

# 2. First plot: your custom scatter
par(mar = c(0, 0, 0, 0) + 0.1)
set.seed(1)
plot(
  runif(10, 0, 1), rep(1, 10), pch = 19,
  col = c(rep('#1FA187', 9), 'tomato'),
  ylim = c(-1, 1), xlim = c(0, 1.5),
  axes = FALSE, xlab = "", ylab = ""
)
for (i in seq(0.8, -1, by = -0.2)) {
  if (i == -0.8){set.seed(14)}
  points(runif(10, 0, 1), rep(i, 10), pch = 19,
         col = c(rep('#1FA187', 9), 'tomato'))
}
ranks <- c(1, 9, 5, 3, 6, 5, 9, 10, 2, 2, 4)
ys <- seq(1, -1, length.out = 11)
text(x = 1.01, y = ys, labels = paste0("Rank = ", ranks), pos = 4)

# 3. Second plot: histogram
par(mar = c(4, 2, 2, 2))
hist(ranks, main = "Rank Histogram", breaks = seq(0.5, 10.5, by = 1), xlab = "Rank", col = 'tomato')
axis(side = 1, at = 1:10)

# 4. Draw arrow across figure AFTER both plots
# Switch to figure coordinates
par(xpd = NA)  # allow drawing outside plot region

# Get figure coordinates
fig <- par("fig")
usr <- par("usr")

# Add arrow from middle right of left plot to middle left of right plot
# For mfrow = c(1,2), each plot takes ~50% horizontally.
# So approximate:
arrows(x0 = -3, y0 = 1, x1 = -1.75, y1 = 1,
       code = 2, length = 0.1, lwd = 2)
#RankHistDemo, 13 x 4

#------------------------------------------------------------------------------------------
#----------------------- Compare random field & independent errors ------------------------
#------------------------------------------------------------------------------------------
 

i = 31
#c(16, 23, 31, 67, 68, 70, 89, 94, 124, 125, 135, 139, 164, 170)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/July25Agg/"

folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects/")
ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
event_ids_all <- as.numeric(sub(".*_(\\d+)$", "\\1", ufiles))

# for (i in 29:1){
#   print(i)
#   ODD <- readODD(paste0(folderin, ufiles[which(event_ids_all==i)]))
#   
#   obs1_i = 1
#   obs2_i = 2
#   N_postpred_samples = 50
#   
#   #Omega = AlgoResults$Omega_sample_phys[order(AlgoResults$Omega_sample_phys[,7,AlgoResults$s_finish])[111],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
#   Omega = AlgoResults$Omega_sample_phys[which.max(AlgoResults$Omega_sample_phys[,7,AlgoResults$s_finish]/AlgoResults$Omega_sample_phys[,8,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
#   #Omega = AlgoResults$Omega_sample_phys[order(AlgoResults$Omega_sample_phys[,7,AlgoResults$s_finish]/AlgoResults$Omega_sample_phys[,8,AlgoResults$s_finish])[500],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
#   Omega$eps$local
#   Omega$eps$hazard_mort
#   
#   sampled_out1 <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
#                           replace(which(names(AlgoParams)==c('Np')), N_postpred_samples), 
#                         output='SampledAgg')
#   
#   impact_sample_flattened=flattenImpactSample2(sampled_out1)
#   impact_sample_flattened$polygon_name = paste0(impact_sample_flattened$impact, impact_sample_flattened$polygon)
#   #impact_sample_flattened %<>% filter(impact=='mortality')
#   print(pairplot_regions(impact_sample_flattened[1:min(10, nrow(impact_sample_flattened)),]))#[1:10,])#[1:10,])
#   readline(prompt="Press [enter] to continue")
# }
ODD <- readODD(paste0(folderin, ufiles[which(event_ids_all==i)]))

obs1_i = 1
obs2_i = 4
N_postpred_samples = 2000

# Independent errors
source('RCode/Model Variations/ODDobj_IndErr.R')
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2025-07-31_093541.918622__July25Agg_NormalCDF_ESplus0.05MVR_RF_kappa0.5_range0.5')
AlgoResults %<>% addAlgoParams()
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/July25Agg/'

#Omega = AlgoResults$Omega_sample_phys[order(AlgoResults$Omega_sample_phys[,7,AlgoResults$s_finish])[111],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
Omega = AlgoResults$Omega_sample_phys[which.max(AlgoResults$Omega_sample_phys[,7,AlgoResults$s_finish]* AlgoResults$Omega_sample_phys[,8,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
#Omega = AlgoResults$Omega_sample_phys[order(AlgoResults$Omega_sample_phys[,7,AlgoResults$s_finish]/AlgoResults$Omega_sample_phys[,8,AlgoResults$s_finish])[500],,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
Omega$eps$local
Omega$eps$hazard_mort

sampled_out1 <- DispX(ODD, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                       replace(which(names(AlgoParams)==c('Np')), N_postpred_samples), 
                     output='SampledAgg')

impact_sample_flattened=flattenImpactSample2(sampled_out1)
impact_sample_flattened$polygon_name = paste0(impact_sample_flattened$impact, impact_sample_flattened$polygon)
#impact_sample_flattened %<>% filter(impact=='mortality')
pairplot_regions(impact_sample_flattened[1:min(10, nrow(impact_sample_flattened)),])#[1:10,])#[1:10,])

joint_pred1 = data.frame(
  sample = 1:length(sampled_out1),
  india_disp_sampled = unlist(lapply(sampled_out1, function(x) x$sampled[obs1_i])),
  china_disp_sampled = unlist(lapply(sampled_out1, function(x) x$sampled[obs2_i])),
  class = 'Independent errors'
)

#Random over errors
source('RCode/ODDobj.R')
#AlgoResults <- readRDS(paste0(dir, 'IIDIPUS_Results/mcmc_2025-04-21_005554.153622_MCMC_BDScore.05_M100AprAgg5_RandomFieldThree_DispOnly_backup'))
#AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2025-04-10')
#AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/Apr25Agg/'
Omega_rf = AlgoResults$Omega_sample_phys[which.max(AlgoResults$Omega_sample_phys[,7,AlgoResults$s_finish]* AlgoResults$Omega_sample_phys[,8,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)
#Omega_rf = AlgoResults$Omega_sample_phys[,which.max(AlgoResults$Omega_sample_phys[7,]/AlgoResults$Omega_sample_phys[8,]*AlgoResults$Omega_sample_phys[9,])] %>% relist(skeleton=Model$skeleton)
#Omega_rf = AlgoResults$Omega_sample_phys[,which.max(AlgoResults$Omega_sample_phys[7,]/AlgoResults$Omega_sample_phys[8,])] %>% relist(skeleton=Model$skeleton)
#Omega_rf = AlgoResults$Omega_sample_phys[,which(AlgoResults$Omega_sample_phys[7,]/AlgoResults$Omega_sample_phys[8,] * AlgoResults$Omega_sample_phys[9,]> 1.5 & AlgoResults$Omega_sample_phys[1,]>9.5 & AlgoResults$Omega_sample_phys[2,]<1.8)[1]] %>% relist(skeleton=Model$skeleton)
#Omega_rf = AlgoResults$Omega_sample_phys[,which(AlgoResults$Omega_sample_phys[7,]/AlgoResults$Omega_sample_phys[8,] * AlgoResults$Omega_sample_phys[9,]> 1 & AlgoResults$Omega_sample_phys[2,]<1.35)[1]] %>% relist(skeleton=Model$skeleton)
# Omega_rf$Lambda1$nu = 9
# Omega_rf$Lambda1$kappa = 0.8
# Omega_rf$eps$local = 0.7
# Omega_rf$eps$hazard_mort = 0.3
# Omega_rf$eps$hazard_disp = 0.7

#plot_joint(event_ids=70, AlgoParams, Omega=Omega_rf, N_postpred_samples = 100)

sampled_out2 <- DispX(ODD, Omega_rf %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                       replace(which(names(AlgoParams)==c('Np')), N_postpred_samples), 
                     output='SampledAgg')


joint_pred2 = data.frame(
  sample = 1:length(sampled_out2),
  india_disp_sampled = unlist(lapply(sampled_out2, function(x) x$sampled[obs1_i])),
  china_disp_sampled = unlist(lapply(sampled_out2, function(x) x$sampled[obs2_i])),
  class = 'Spatially correlated errors'
)

impact_sample_flattened2=flattenImpactSample2(sampled_out2)
impact_sample_flattened2$polygon_name = paste0(impact_sample_flattened2$impact, impact_sample_flattened2$polygon)
#impact_sample_flattened %<>% filter(impact=='mortality')
pairplot_regions(impact_sample_flattened2[1:min(10, nrow(impact_sample_flattened2)),])#[1:10,])#[1:10,])

joint_pred = rbind(joint_pred1, joint_pred2)
saveRDS(joint_pred, paste0(dir, 'IIDIPUS_Results/NPL_jointpred2000'))

#joint_pred500 = readRDS(paste0(dir, 'IIDIPUS_Results/NPL_jointpred'))
#joint_pred2000 = readRDS(paste0(dir, 'IIDIPUS_Results/NPL_jointpred2000'))

p1 <- ggplot(joint_pred, aes(x = india_disp_sampled, y = china_disp_sampled, color = class)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  geom_point(aes(x = sampled_out1[[1]]$observed[obs1_i], y = sampled_out1[[1]]$observed[obs2_i], color = "Observed"), size = 4, shape = 4) +  # Crimson red point
  xlab("Population Displacement (India)") +
  ylab("Population Displacement (China)") +
  theme(
    axis.title.y = element_text(family = "Times New Roman", size = 12),
    axis.text.x = element_text(family = "Times New Roman", size = 12),
    axis.text.y = element_text(family = "Times New Roman", size = 12),
    axis.title.x = element_text(family = "Times New Roman", size = 12),
    plot.title = element_text(family = "Times New Roman", size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0, 20, 0, 15), "pt"),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank()  # Remove legend title
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  scale_color_manual(
    values = c(
      "Independent errors" = "#1FA187",           # Viridis green
      "Spatially correlated errors" = "#440154",  # Viridis purple
      "Observed" = "#D62728"                      # Crimson red last
    ),
    breaks = c("Independent errors", "Spatially correlated errors", "Observed")  # Control order
  )

p2 <- ggplot(joint_pred, aes(x = china_disp_sampled, color = class, fill = class)) +
  geom_density(alpha = 0.2, adjust = 1) +
  scale_x_log10(labels = scales::comma_format()) +
  xlab("Predicted Population Displacement in China") +
  ylab("Density") +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "Independent errors" = "#1FA187",
      "Spatially correlated errors" = "#440154"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Independent errors" = "#1FA187",
      "Spatially correlated errors" = "#440154"
    )
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank(),
    legend.text = element_text(family = "Times New Roman", size = 12)
  ) + 
  geom_vline(xintercept = sampled_out1[[1]]$observed[obs2_i], color = "#D62728", size = 1)

filtered_df <- joint_pred %>% filter(india_disp_sampled > 10000 & india_disp_sampled < 50000)

# Plot conditional densities
p3 <- ggplot(filtered_df, aes(x = china_disp_sampled, color = class, fill = class)) +
  geom_density(alpha = 0.2, adjust = 0.8) +
  scale_x_log10(labels = scales::comma_format()) +
  xlab("Predicted Pop Disp in China (given Pop Disp in India)") +
  ylab("Density") +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "Independent errors" = "#1FA187",
      "Spatially correlated errors" = "#440154"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Independent errors" = "#1FA187",
      "Spatially correlated errors" = "#440154"
    )
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank()
  ) + 
  geom_vline(xintercept = sampled_out1[[1]]$observed[obs2_i], color = "#D62728", size = 1)

plot_grid(
  p1 + theme(
    legend.position = "bottom",
    legend.spacing = unit(2, "cm")
  ),
  plot_grid(
    p2, p3,
    nrow = 2,
    labels = c('(b)', '(c)'),
    label_fontfamily = "Times New Roman",
    label_fontface = "plain",
    align = "v"
  ),
  labels = c('(a)', '', ''),
  ncol = 2,
  nrow = 1,
  rel_widths = c(0.5, 0.5),
  align = "h",
  axis = "t",  # <-- This aligns top of p1 and p2 nicely
  label_fontfamily = "Times New Roman",
  label_fontface = "plain"
)
# NPL_Disp_Correlations.pdf,14 x 6

#-----------------------------------------------------------------------------------------
#---------------------- ESTIMATING BEST MATERN PARAMS ------------------------------------
#-----------------------------------------------------------------------------------------

library(geoR)
library(fields)
ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/July25Agg/ODDobjects/Train/EQ20140801DZA_24')

AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/July25Agg/"
folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects/")
ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
ufiles <- grep('^Train/' , ufiles, value = TRUE)

for (file in ufiles){
  ODDy = readODD(paste0(folderin, file))
  
  coords <- crds(ODDy, na.rm=F)
  values <- as.numeric(values(ODDy$hazMean1,na.rm=F))
  df <- data.frame(x = coords[,1], y = coords[,2], z = values)
  df = df[!is.na(df$z),]
  
  # Convert to spatial object
  coordinates(df) <- ~x + y
  vg_emp <- variogram(z ~ 1, data = df)
  plot(vg_emp)
  vg_fit <- fit.variogram(vg_emp, model = vgm("Mat"), fit.kappa=T)
  plot(vg_emp, vg_fit, main = "Fitted Exponential Variogram")
  vg_fit$kappa
  vg_fit$range
}

grid <- as.data.frame(xyFromCell(ODDy, 1:ncell(ODDy)))  # Extract grid coordinates
names(grid) <- c("x", "y")  # Name the columns
vgm_model <- vgm(psill = 1, model = "Mat", range = 0.5, kappa = 0.5)
gstat_mod = gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 3)

# set.seed(1)
RF_local = predict(gstat_mod, newdata = grid, nsim = 1000)

vgm_model <- vgm(psill = 1, model = "Exp", range = 0.5)
gstat_mod = gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 3)

# set.seed(1)
RF_local = predict(gstat_mod, newdata = grid, nsim = 1)[,3]
r <- rast(ODD)
r[['hi']]= RF_local
plot(r$hi)

vg_fit$kappa


vg_fit <- fit.variogram(vg_emp, model = vgm("Exp"))
plot(vg_emp, vg_fit, main = "Fitted Exponential Variogram")

set.seed(123)

fit <- spatialProcess(coords, values, Covariance = "Matern")


data_geo <- as.geodata(cbind(coords, values))


fit <- likfit(data_geo,
              ini.cov.pars = c(1, 1),   # initial: c(sill, range)
              nugget = 0.1,
              fix.nugget = FALSE,
              cov.model = "matern",
              kappa = 0.5,               # fixed smoothness
              fix.kappa = TRUE)         # only range is estimated

summary(fit)











ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/July25Agg/ODDobjects/Train/EQ20170919MEX_68')
#obs1_i = 1
#obs2_i = 4
N_postpred_samples = 100

# Independent errors
source('RCode/Model Variations/ODDobj_VarRF.R')

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2025-07-31_093541.918622__July25Agg_NormalCDF_ESplus0.05MVR_RF_kappa0.5_range0.5')
AlgoResults %<>% addAlgoParams()
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/July25Agg/'

Omega = AlgoResults$Omega_sample_phys[which.max(AlgoResults$Omega_sample_phys[,7,AlgoResults$s_finish]* AlgoResults$Omega_sample_phys[,8,AlgoResults$s_finish]),,AlgoResults$s_finish] %>% relist(skeleton=Model$skeleton)

Omega$eps$local
Omega$eps$hazard_mort
Omega$RF_range=2
Omega$RF_kappa=0.5
sampled_out1 <- DispX(ODDy, Omega %>% addTransfParams(), Model$center, AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% 
                        replace(which(names(AlgoParams)==c('Np')), N_postpred_samples), 
                      output='SampledAgg')

impact_sample_flattened=flattenImpactSample2(sampled_out1)
impact_sample_flattened$polygon_name = paste0(impact_sample_flattened$impact, impact_sample_flattened$polygon)
#impact_sample_flattened %<>% filter(impact=='mortality')
pairplot_regions(impact_sample_flattened[1:min(10, nrow(impact_sample_flattened)),])#[1:10,])#[1:10,])


#------------------------------------------------------------------------------------------
#----------------------- Compare ES fit and ES + BD fit -----------------------------------
#------------------------------------------------------------------------------------------

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-27_162932.543295_MCMC_BDScore.05nocorr_LR40_M100_Npart1000NovAgg5_RandomFieldThree')
AlgoResults$input_folder <- 'IIDIPUS_Input_Alternatives/Mar25Agg/'

ii <- which(AlgoResults$Omega_sample_phys[7,] > 0.5 & AlgoResults$Omega_sample_phys[7,] < 0.55)
Omega_rf = AlgoResults$Omega_sample_phys[,1600] %>% relist(skeleton=Model$skeleton)
Omega_rf$eps$local
Omega_rf$eps$hazard_mort
set.seed(123)
plot_joint(event_ids=70, AlgoParams, Omega=Omega, N_postpred_samples = 100)


#------------------------------------------------------------------------------------------
#----------------------- Compare MMDS sample and ES sample --------------------------------
#------------------------------------------------------------------------------------------
corr_seq = seq(-1, 1, 0.05)
n_samples = 100
es_samples = array(0, dim=c(n_samples, length(corr_seq)))
vs_samples = array(0, dim=c(n_samples, length(corr_seq)))
mmds_samples = array(0, dim=c(n_samples, length(corr_seq)))

for (i in 1:n_samples){
  print(i)
  dat = t(rmvnorm(1, c(0,0), rbind(c(1,0.5), c(0.5,1))))
  for (corr_i in 1:length(corr_seq)){
    corr = corr_seq[corr_i]
    preds = rmvnorm(100, c(0,0), rbind(c(1,corr), c(corr,1)))
    es_samples[i, corr_i] = es_sample(c(dat), t(preds))
    vs_samples[i, corr_i] = vs_sample(c(dat), t(preds))
    mmds_samples[i, corr_i] = mmds_sample(c(dat), t(preds))
  }
}
plot(rep(corr_seq, n_samples), c(es_samples)); points(corr_seq, apply(es_samples, 2, mean), col='red')
plot(rep(corr_seq, each=n_samples), c(vs_samples), ylim=c(0,2)); points(corr_seq, apply(vs_samples, 2, mean), col='red')
plot(rep(corr_seq, each=n_samples), c(mmds_samples)); points(corr_seq, apply(mmds_samples, 2, mean), col='red')

plot(corr_seq, apply(es_samples, 2, median))
plot(corr_seq, apply(mmds_samples, 2, median))

ad_samples = c()
for (i in 1:10000){
  ad_samples = c(ad_samples, AndersonDarlingTest(runif(50), null='punif')$statistic)
}
hist(ad_samples)
sd(ad_samples)

#--------------------------------------------------------------------------------------
# -----------  COMPARE FIT WITH AND WITHOUT MVR RANK HISTOGRAMS -----------------------
#--------------------------------------------------------------------------------------

getRanks <- function(impact_sample_poly, AlgoParams, corr_noise = NA, corr_noise2 = NA){
  #Calculate the mean energy score from impact_sample
  # Note that dist_poly has length 7, however, only the first element is currently used when calculating 
  # the distance. It has been kept with length greater than 1 in case we want to store up to 6 other values 
  # (e.g. the Anderson Darling Statistic of pre-ranks using the Minimum Spanning Tree), but these other
  # elements are not actually currently used when calculating the distance. 
  
  observed <- impact_sample_poly[[1]]$observed
  
  dist_poly <- array(NA, dim=c(AlgoParams$Np,8))
  impact_type <- impact_sample_poly[[1]]$impact
  impact_weightings <- unlist(AlgoParams$impact_weights[impact_type])
  event_id <- impact_sample_poly[[1]]$event_id
  grouped_events <- split(seq_along(event_id), event_id)
  
  for(n in 1:AlgoParams$Np){
    
    samples_allocated <- ((n-1)*AlgoParams$m_CRPS+1):(n*AlgoParams$m_CRPS)
    samples_combined <- sapply(impact_sample_poly[samples_allocated], function(x){x$sampled}) #doesn't work if samples_allocated is length 1
    
    dist_poly[n,1] <- 0
    
    es_store <- c()
    mmds_store <- c()
    mmds_store_nobd <- c()
    pre_ranks_bd <- c()
    pre_ranks_mmds <- c()
    pre_ranks_mvr <- c()
    pre_ranks_bd_mortdisp <- c()
    pre_ranks_mmds_mortdisp <- c()
    pre_ranks_mvr_mortdisp <- c()
    #pre_ranks_average <- c() #can also assess quantiles for uniformity based on the average pre-rank function
    #pre_ranks_mst <- c() #can also assess quantiles for uniformity based on the minimum spanning tree pre-rank function
    for (i in 1:length(grouped_events)){
      #For each event, compute the energy score of the observed data vs the 'prediction' (simulated data)
      #Each impact type is weighted differently, simply multiplying the observation and the simulations by this weight performs the weighting
      obs <- log(observed[grouped_events[[i]]]+AlgoParams$log_offset) *impact_weightings[grouped_events[[i]]]
      sims <- log(samples_combined[grouped_events[[i]],, drop=F]+AlgoParams$log_offset) * impact_weightings[grouped_events[[i]]]
      
      obs_mortdisp <- log(observed[intersect(grouped_events[[i]], which(impact_type!='buildDam'))]+AlgoParams$log_offset) *impact_weightings[intersect(grouped_events[[i]], which(impact_type!='buildDam'))]
      sims_mortdisp <- log(samples_combined[intersect(grouped_events[[i]], which(impact_type!='buildDam')),,drop=F]+AlgoParams$log_offset) * impact_weightings[intersect(grouped_events[[i]], which(impact_type!='buildDam'))]
      #obs_unweighted <- log(observed[grouped_events[[i]]]+AlgoParams$log_offset)# * impact_weightings[grouped_events[[i]]]
      #sims_unweighted <- log(samples_combined[grouped_events[[i]],,drop=F]+AlgoParams$log_offset)# * impact_weightings[grouped_events[[i]]]
      
      #nrow(sims_rh)
      #pairs(rbind(t(sims_rh), obs_rh), col=c(rep('black',ncol(sims_rh)), 'red'))
      #mort_obs1 = which(impact_type[grouped_events[[i]]] == 'mortality')[1]
      #disp_obs1 = which(impact_type[grouped_events[[i]]] == 'displacement')[1]
      # 
      # pairs(rbind(t(sims_rh), obs_rh), col=c(rep('black',ncol(sims_rh)), 'red'))
      # get_banddepth_rank_single(cbind(obs_rh[2:3], sims_rh[2:3,]))
      # 
      # mat = cbind(obs_rh[c(2,4)], sims_rh[c(2,4),])
      # vss=c()
      # for (j in 1:ncol(mat)){
      #   vss = c(vss, mmds_sample(mat[,j], mat[,-j]))
      # }
      # plot(vss)
      # vss[1]
      # xx = data.frame(x = mat[1,], y=mat[2,], col=vss)
      # ggplot(xx, aes(x=x, y=y, col=col))+geom_point()
      # if (length(obs_rh)>5) stop()
      #if(get_banddepth_rank_single(cbind(obs_rh, sims_rh)) < 2) stop()
      
      if (length(grouped_events[[i]])==1){
        #crps is the univariate case of the energy score, so use when dimension is 1. 
        es_store<- c(es_store, crps_sample(obs, sims))
        mmds_store<-c(mmds_store, 0)
        
        rank = rank(c(obs,sims), ties.method='random')[1]
        pre_ranks_bd <- c(pre_ranks_bd, rank)
        pre_ranks_mmds <- c(pre_ranks_mmds, rank)
        pre_ranks_mvr <- c(pre_ranks_mvr, rank)
        pre_ranks_bd_mortdisp <- c(pre_ranks_bd_mortdisp,rank)
        pre_ranks_mmds_mortdisp <- c(pre_ranks_mmds_mortdisp,rank)
        pre_ranks_mvr_mortdisp <- c(pre_ranks_mvr_mortdisp,rank)
        
        next
      } 
      
      # if (nrow(sims_rh) > 1){
      #   plot(t(sims_rh[1:2,]))
      #   points(t(obs_rh[1:2]), col='red', pch=19)
      # }
      
      es_store <- c(es_store, es_sample(obs, sims))
      mmds_store <- c(mmds_store, mmds_sample(obs, sims))
      
      if (length(obs_mortdisp)==1){
        obs_sims_mat = cbind(obs,sims)
        obs_sims_mat_mortdisp = c(obs_mortdisp,sims_mortdisp)
        rank = rank(obs_sims_mat_mortdisp, ties.method='random')[1]
        
        pre_ranks_bd <- c(pre_ranks_bd, get_banddepth_rank_single(obs_sims_mat))
        pre_ranks_mmds <- c(pre_ranks_mmds, get_mmds_rank_single(obs_sims_mat))
        pre_ranks_mvr <- c(pre_ranks_mvr, get_mvr_rank_single(obs_sims_mat))
        
        pre_ranks_bd_mortdisp <- c(pre_ranks_bd_mortdisp,rank)
        pre_ranks_mmds_mortdisp <- c(pre_ranks_mmds_mortdisp,rank)
        pre_ranks_mvr_mortdisp <- c(pre_ranks_mvr_mortdisp,rank)
      } else {
        obs_sims_mat = cbind(obs,sims)
        obs_sims_mat_mortdisp = cbind(obs_mortdisp,sims_mortdisp)
        
        pre_ranks_bd <- c(pre_ranks_bd, get_banddepth_rank_single(obs_sims_mat))#, noise=corr_noise[i,]))
        pre_ranks_bd_mortdisp <- c(pre_ranks_bd_mortdisp, get_banddepth_rank_single(obs_sims_mat_mortdisp))
        
        pre_ranks_mmds <- c(pre_ranks_mmds, get_mmds_rank_single(obs_sims_mat))
        pre_ranks_mmds_mortdisp <- c(pre_ranks_mmds_mortdisp, get_mmds_rank_single(obs_sims_mat_mortdisp))
        
        pre_ranks_mvr <- c(pre_ranks_mvr, get_mvr_rank_single(obs_sims_mat))#, noise=corr_noise[i,]))
        pre_ranks_mvr_mortdisp <- c(pre_ranks_mvr_mortdisp, get_mvr_rank_single(obs_sims_mat_mortdisp))
      }
      
      
    }
    
    ranks_bd_std <- (pre_ranks_bd-runif(length(pre_ranks_bd),0,1))/(AlgoParams$m_CRPS + 1)
    ranks_bd_mortdisp_std <- (pre_ranks_bd_mortdisp -runif(length(pre_ranks_bd_mortdisp),0,1))/(AlgoParams$m_CRPS + 1)
    ranks_mmds_std <- (pre_ranks_mmds-runif(length(pre_ranks_mmds),0,1))/(AlgoParams$m_CRPS + 1)
    ranks_mmds_mortdisp_std <- (pre_ranks_mmds_mortdisp-runif(length(pre_ranks_mmds_mortdisp),0,1))/(AlgoParams$m_CRPS + 1)
    ranks_mvr_std <- (pre_ranks_mvr-runif(length(pre_ranks_mvr),0,1))/(AlgoParams$m_CRPS + 1)
    ranks_mvr_mortdisp_std <- (pre_ranks_mvr_mortdisp-runif(length(pre_ranks_mvr_mortdisp),0,1))/(AlgoParams$m_CRPS + 1)
    #ranks_std <- (pre_ranks-corr_noise2)/(AlgoParams$m_CRPS + 1) #can also correlate the noise added here
    
    return(ranks_mvr_std)
  }
  
  return(dist_poly)
  
}


source('RCode/ODDobj.R')
source('RCode/Model.R')
#AlgoResults_mvr <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-06_112257.598195__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain3')
AlgoResults_mvr <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-03_101635.044837__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain4')
#AlgoResults_mvr$loss[12001]
Omega_mvr <- AlgoResults_mvr$Omega_sample_phys[,13001] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/July25Agg/"
AlgoResults_mvr$input_folder <- "IIDIPUS_Input_Alternatives/July25Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
set.seed(123)
for (i in 1:10){
  impact_sample_mvr <- SampleImpact(dir, Model, Omega_mvr %>% addTransfParams(), AlgoParams, dat='Train')
  ranks_mvr2 = getRanks(impact_sample_mvr$poly, AlgoParams)
  saveRDS(ranks_mvr2, paste0('/home/manderso/Documents/GitHub/ODDRIN/Ranks/MVRranks', i))
}
hist(ranks_mvr2)
#impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train', output='SampledTotal')
#saveRDS(impact_sample_rf, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
#CalcDist(impact_sample_mvr, AlgoParams)
#df_mvr <- flattenImpactSample(impact_sample_mvr)
#ranks_mvr <- get_multivariate_ranking(df_mvr, log=T)
#hist(ranks_mvr, main='Multivariate Rank Scores')

# ranks_mst <- get_mst_rank(df_bd, log=T)
# hist(ranks_mst, main='MST Rank Scores')


AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-08-07_173458.049397__July25Agg_Normal_ESOnly_RFwithTot_range.5_kappa.5_Chain1')
for (i in c(2001, 2251, 2501, 2751, 3001)){
  for (j in 1:5){
    Omega <- AlgoResults$Omega_sample_phys[,i] %>% relist(skeleton=Model$skeleton)
    AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/July25Agg/"
    AlgoResults$input_folder <- "IIDIPUS_Input_AlternativesJuly25Agg/"
    AlgoParams$Np <- 1
    AlgoParams$m_CRPS <- 100
    impact_sample <- SampleImpact(dir, Model, Omega %>% addTransfParams(), AlgoParams, dat='Train')
    ranks_mvr_og = getRanks(impact_sample$poly, AlgoParams)
    saveRDS(ranks_mvr_og, paste0('/home/manderso/Documents/GitHub/ODDRIN/Ranks/MVRranksOG', i, '_', j))
  }
}
hist(ranks_mvr_og)
#impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train', output='SampledTotal')
#saveRDS(impact_sample_rf, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
# df_low <- flattenImpactSample(impact_sample)
# ranks_mvr_og <- get_multivariate_ranking(df_low, log=T)
# hist(ranks_mvr_og, main='Banddepth Rank Scores')
# 
# ranks_mvr_og <- get_mst_rank(df_low, log=T)
# hist(ranks_mvr_og, main='MST Rank Scores')


ranks_mvr2 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/Ranks/MVRranks2')
par(family = "Times")
par(mfrow=c(1,2))
par(mar = c(4.5, 4, 0, 2)) 
hist(ranks_mvr2, col = "#c3b1c7",  # light purple fill
     border = "#440154",
     xlab='Ranks (loss function with rank testing)', freq=T, ylab='Frequency', main='',cex.main = 1.2, font.main = 1,  ylim=c(0,32))
mtext("(a)", side = 3,  line = -1, adj = -0.2,  font = 1, cex = 1.2)

ranks_mvr_og <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/Ranks/MVRranksOG3001_5')
hist(ranks_mvr_og, col = "#c3b1c7",  # light purple fill
     border = "#440154", xlab='Ranks (loss function without rank testing)', freq=T, ylab='Frequency', main='',font.main = 1, ylim=c(0,32))
mtext("(b)", side = 3, line = -1, adj = -0.2, font = 1, cex = 1.2)

#RankComparison.pdf, 4 x 10


# hist(ranks_mst, col = "#c3b1c7",  # light purple fill
#      border = "#440154",
#      xlab='MST Ranks (loss function with rank testing)', freq=T, ylab='Frequency', main='',cex.main = 1.2, font.main = 1,  ylim=c(0,32))
# mtext("(c)", side = 3,  line = -1, adj = -0.22,  font = 1, cex = 1.2)
# 
# hist(ranks_mst_og, col = "#c3b1c7",  # light purple fill
#      border = "#440154", xlab='MST Ranks (loss function without rank testing)', freq=T, ylab='Frequency', main='',font.main = 1, ylim=c(0,32))
# mtext("(d)", side = 3, line = -1, adj = -0.22, font = 1, cex = 1.2)

#RankComparison.pdf, 4 x 10

#OLD ------------------------------------------------------------------------------------------

# for (i in 1:5){
#   Omegatrue_high_cor <- Omega_true
#   Omegatrue_high_cor$eps$hazard_cor <- 0.95
#   Omegatrue_low_cor <- Omega_true
#   Omegatrue_low_cor$eps$hazard_cor <- 0.05
#   impact_sample_high_cor <- SampleImpact(dir, Model, Omegatrue_high_cor %>% addTransfParams(), AlgoParams, dat='all')
#   dist_high <- CalcDist(impact_sample_high_cor, AlgoParams)
#   print(paste('high',dist_high))
#   df_high3 <- flattenImpactSample(impact_sample_high_cor)
#   
#   impact_sample_low_cor <- SampleImpact(dir, Model, Omegatrue_low_cor %>% addTransfParams(), AlgoParams, dat='all')
#   dist_low <- CalcDist(impact_sample_low_cor, AlgoParams)
#   print(paste('low',dist_low))
#   df_low3 <- flattenImpactSample(impact_sample_low_cor)
#   
#   impact_sample_true_cor <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams, dat='all')
#   dist_true <- CalcDist(impact_sample_true_cor, AlgoParams)
#   print(paste('true', dist_true))
#   df_true3 <- flattenImpactSample(impact_sample_true_cor)
#   
# }

# saveRDS(list(
#   df_high1 = df_high1,
#   df_high2 = df_high2,
#   df_high3 = df_high3,
#   df_low1 = df_low1,
#   df_low2 = df_low2,
#   df_low3 = df_low3,
#   df_true1 = df_true1,
#   df_true2 = df_true2,
#   df_true3 = df_true3
# ), 'ImpactSample_varyingcorrelation')


# Omegatrue_high_locerror <- Omega_true
# Omegatrue_high_locerror$eps$local <- 0.1
# Omegatrue_low_locerror <- Omega_true
# Omegatrue_low_locerror$eps$local <- 1.5
# impact_sample_high_locerror <- SampleImpact(dir, Model, Omegatrue_high_locerror %>% addTransfParams(), AlgoParams, dat='all')
# dist_high <- CalcDist(impact_sample_high_locerror, AlgoParams)
# print(paste('high',dist_high))
# df_highloc1 <- flattenImpactSample(impact_sample_high_locerror)
# 
# impact_sample_low_locerror <- SampleImpact(dir, Model, Omegatrue_low_locerror %>% addTransfParams(), AlgoParams, dat='all')
# dist_low <- CalcDist(impact_sample_low_locerror, AlgoParams)
# print(paste('low',dist_low))
# df_lowloc1 <- flattenImpactSample(impact_sample_low_locerror)
# 
# impact_sample_true_locerror <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams, dat='all')
# dist_true <- CalcDist(impact_sample_true_locerror, AlgoParams)
# print(paste('true', dist_true))
# df_trueloc1 <- flattenImpactSample(impact_sample_true_locerror)

# NEW ------------------------------------------------------------------------------------------

# AlgoResults_es <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180047_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2')
# Omega_i <- which(AlgoResults_es$Omega_sample_phys[7,] <0.5 & AlgoResults_es$Omega_sample_phys[8,] >0.8)
# Omega_es <- AlgoResults_es$Omega_sample_phys[,Omega_i[Omega_i > 5000][1]] %>% relist(skeleton=Model$skeleton)
# AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
# AlgoParams$Np <- 1
# AlgoParams$m_CRPS <- 100
# impact_sample_es <- SampleImpact(dir, Model, Omega_es %>% addTransfParams(), AlgoParams, dat='Train')
# df_es <- flattenImpactSample(impact_sample_es)
# df_low <- df_es

# Energy score fit & independent errors, best case scenario (high local error and displacement error)
AlgoResults_es <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-11-21_180047_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP_smallerStartPropCOV_NovDat2')
Omega_i <- which(AlgoResults_es$Omega_sample_phys[7,] > 1.4 & AlgoResults_es$Omega_sample_phys[9,] > 1.3 & AlgoResults_es$Omega_sample_phys[8,] > 0.75)
Omega_es <- AlgoResults_es$Omega_sample_phys[,Omega_i[Omega_i > 5000][1]] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_es <- SampleImpact(dir, Model, Omega_es %>% addTransfParams(), AlgoParams, dat='Train')
saveRDS(impact_sample_es, paste0(dir, 'IIDIPUS_Results/impact_sample_es_Feb10_best_case_errors'))
df_es <- flattenImpactSample(impact_sample_es)
df_low <- df_es

# Variogram score fit & random field errors, with correlation 1 between random fields across impact types, best case scenario (high local error and displacement error)
source('RCode/ODDobj_rf1.R')
AlgoResults_rf1 <-readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-01-23_123740_MCMC_VariogramScore_M100_Npart1000NovAgg5_RandomField')
Omega_rf1 <- AlgoResults_rf1$Omega_sample_phys[,290] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_rf1 <- SampleImpact(dir, Model, Omega_rf1 %>% addTransfParams(), AlgoParams, dat='Train')
#saveRDS(impact_sample_rf1, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
df_rf <- flattenImpactSample(impact_sample_rf1)
df_low <- df_rf

# Variogram score fit & random field errors, with correlation <1 between random fields across impact types, best case scenario (high local error and displacement error)
source('RCode/ODDobj_rf.R')
AlgoResults_rf <-readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-01-31_172601_MCMC_VariogramScore_M100_Npart1000NovAgg5_RandomFieldThree')
Omega_rf <- AlgoResults_rf$Omega_sample_phys[,290] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train')
#saveRDS(impact_sample_rf1, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
df_rf <- flattenImpactSample(impact_sample_rf)
df_low <- df_rf

# Energy score fit only & random field errors, with correlation <1 between random fields across impact types, best case scenario (high local error and displacement error)
source('RCode/ODDobj_rf.R')
AlgoResults_rf <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-05_135459.168382_MCMC_VariogramScore_M100_Npart1000NovAgg5_RandomFieldThree')
Omega_rf <- AlgoResults_rf$Omega_sample_phys[,which(AlgoResults_rf$Omega_sample_phys[7,]>0.8)[1]] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train')
impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train', output='SampledTotal')
saveRDS(impact_sample_rf, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
df_rf <- flattenImpactSample(impact_sample_rf)
df_low <- df_rf

# Energy score fit & random field errors, with correlation <1 between random fields across impact types, not quite best case scenario 
source('RCode/ODDobj_rf.R')
AlgoResults_rf <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-05_135459.168382_MCMC_VariogramScore_M100_Npart1000NovAgg5_RandomFieldThree')
Omega_rf <- AlgoResults_rf$Omega_sample_phys[,which(is.na(AlgoResults_rf$Omega_sample_phys[7,]))-1] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train')
impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train', output='SampledTotal')
#saveRDS(impact_sample_rf, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
df_rf <- flattenImpactSample(impact_sample_rf)
df_low <- df_rf

# Energy score + 0.02 * Anderson Darling(MST Rank Histograms) with correlation <1 between random fields across impact types, better case scenario (small errors)
AlgoResults_rh <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/mcmc_2025-02-15_171140.550783_MCMC_MstScore.02_M100_Npart1000NovAgg5_RandomFieldThree')
Omega_i <- which(AlgoResults_rh$Omega_sample_phys[7,] > 0.8 & AlgoResults_rh$Omega_sample_phys[8,] > 0.5 & AlgoResults_rh$Omega_sample_phys[9,] > 0.9)
Omega_rh <- AlgoResults_rh$Omega_sample_phys[,Omega_i[Omega_i > 1000][1]] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_rh <- SampleImpact(dir, Model, Omega_rh %>% addTransfParams(), AlgoParams, dat='Train')
saveRDS(impact_sample_es, paste0(dir, 'IIDIPUS_Results/impact_sample_es_Feb10_best_case_errors'))
df_rh <- flattenImpactSample(impact_sample_rh)
df_low <- df_rh

# Energy score + 0.05 * Anderson Darling(MST Rank Histograms) with correlation <1 between random fields across impact types, better case scenario (small errors)
AlgoResults_rh <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-22_234005.938826_MCMC_MstScore.05corr_M100_Npart1000NovAgg5_RandomFieldThree')
Omega_i <- which.min(AlgoResults_rh$rh_score) #which(AlgoResults_rh$Omega_sample_phys[7,] > 0.8 & AlgoResults_rh$Omega_sample_phys[8,] > 0.5 & AlgoResults_rh$Omega_sample_phys[9,] > 0.9)
Omega_i = which(AlgoResults_rh$rh_score==unique(sort(AlgoResults_rh$rh_score, decreasing=F))[2])[1]
Omega_rh <- AlgoResults_rh$Omega_sample_phys[,Omega_i[Omega_i > 1000][1]] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_rh <- SampleImpact(dir, Model, Omega_rh %>% addTransfParams(), AlgoParams, dat='Train')
#saveRDS(impact_sample_es, paste0(dir, 'IIDIPUS_Results/impact_sample_es_Feb10_best_case_errors'))
df_rh <- flattenImpactSample(impact_sample_rh)
df_low <- df_rh

# Energy score + 0.05 * Anderson Darling(Band Depth Rank Histograms) with correlation <1 between random fields across impact types, better case scenario (small errors)
AlgoResults_rh <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-24_111536.758596_MCMC_BDScore.05nocorr_M100_Npart1000NovAgg5_RandomFieldThree')
#Omega_i <- which.min(AlgoResults_rh$rh_score) #which(AlgoResults_rh$Omega_sample_phys[7,] > 0.8 & AlgoResults_rh$Omega_sample_phys[8,] > 0.5 & AlgoResults_rh$Omega_sample_phys[9,] > 0.9)
#Omega_i = which(AlgoResults_rh$rh_score==unique(sort(AlgoResults_rh$rh_score, decreasing=F))[2])[1]
Omega_i = 1020 #which(!is.finite(AlgoResults_rh$loss))[1] - 1
Omega_rh <- AlgoResults_rh$Omega_sample_phys[,Omega_i] %>% relist(skeleton=Model$skeleton) #AlgoResults_rh$Omega_sample_phys[,Omega_i[Omega_i > 1000][1]] %>% relist(skeleton=Model$skeleton)
AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
proposed = Omega_rh %>% addTransfParams()
proposed$u = AlgoResults_rh$u[,,,ifelse(Omega_i %% 3==0, 3, Omega_i%%3)] 
proposed$u_local = AlgoResults_rh$RF_current
impact_sample_rh <- SampleImpact(dir, Model, proposed, AlgoParams, dat='Train')
CalcDist(impact_sample_rh, AlgoParams, corr_noise = pnorm(AlgoResults_rh$u_rh), corr_noise2 = pnorm(AlgoResults_rh$u_rh2))
CalcDist(impact_sample_rh, AlgoParams)

#saveRDS(impact_sample_es, paste0(dir, 'IIDIPUS_Results/impact_sample_es_Feb10_best_case_errors'))
df_rh <- flattenImpactSample(impact_sample_rh)
df_low <- df_rh




# Energy score fit & random field errors, with correlation 1 between random fields across impact types, best case scenario (high local error and displacement error)
AlgoResults_rf<- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-05_135459.168382_MCMC_VariogramScore_M100_Npart1000NovAgg5_RandomFieldThree')


AlgoResults_rf <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-01-27_115044.604171_MCMC_VariogramScore_M100_Npart1000NovAgg5_RandomFieldThree')
Omega_rf <- AlgoResults_rf$Omega_sample_phys[,200] %>% relist(skeleton=Model$skeleton)
Omega_rf$eps$local <- 0.8
Omega_rf$eps$hazard_mort <- 0.6
Omega_rf$eps$hazard_disp <- 1

AlgoResults_rf <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-02-05_135459.168382_MCMC_VariogramScore_M100_Npart1000NovAgg5_RandomFieldThree')
Omega_rf <- AlgoResults_rf$Omega_sample_phys[,which(AlgoResults_rf$loss<3.6)[1]] %>% relist(skeleton=Model$skeleton)
Omega_rf <- AlgoResults_rf$Omega_sample_phys[,which(AlgoResults_rf$Omega_sample_phys[7,]>0.8)[1]] %>% relist(skeleton=Model$skeleton)


AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train')
impact_sample_rf <- SampleImpact(dir, Model, Omega_rf %>% addTransfParams(), AlgoParams, dat='Train', output='SampledTotal')
saveRDS(impact_sample_rf, paste0(dir, 'IIDIPUS_Results/impact_sample_rf_vs_score_fitted_6Feb_TotImpact'))
df_rf <- flattenImpactSample(impact_sample_rf)
df_low <- df_rf


AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_small_local <- SampleImpact(dir, Model, Omega_small_local %>% addTransfParams(), AlgoParams, dat='Train')
df_low <- flattenImpactSample(impact_sample_small_local)
par(mfrow=c(2,2))
ranks_mvr <- get_multivariate_ranking(df_low)
hist(ranks_mvr, main='Multivariate Rank Scores')
ranks_avg <- get_average_rank(df_low)
hist(ranks_avg, main='Average Rank Scores')
ranks_mst <- get_mst_rank(df_low, log=T)
hist(ranks_mst, main='Minimum Spanning Tree Rank Scores', breaks=10)
ranks_bdpth <- get_banddepth_rank(df_low, log=T)
hist(ranks_bdpth, main='Banddepth Rank Scores')
par(mfrow=c(1,1))
ranks_uni <- get_univariate_rank(df_low)
hist(ranks_uni, main='Univariate Rank Histogram')



saveRDS(impact_sample_small_local, paste0(dir, 'IIDIPUS_Results/impact_sample_es_score_fitted_single'))

AlgoResults_vs <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-01-10_181501_MCMC_VariogramScore_M100_Npart1000NovAgg5_propCOVmult0.2')
Omega_i <- which(AlgoResults_vs$Omega_sample_phys[7,] >1.3 & AlgoResults_vs$Omega_sample_phys[8,] >0.5)
Omega_vs <- AlgoResults_vs$Omega_sample_phys[,Omega_i[Omega_i > 2000][1]] %>% relist(skeleton=Model$skeleton)

AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_vs <- SampleImpact(dir, Model, Omega_vs %>% addTransfParams(), AlgoParams, dat='Train')
saveRDS(impact_sample_vs, paste0(dir, 'IIDIPUS_Results/impact_sample_vs_score_fitted_range0.5_goodcal'))
df_vs <- flattenImpactSample(impact_sample_vs)


AlgoResults_vs <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-01-10_181501_MCMC_VariogramScore_M100_Npart1000NovAgg5_propCOVmult0.2')
Omega_i <- which(AlgoResults_vs$Omega_sample_phys[7,] >1.3 & AlgoResults_vs$Omega_sample_phys[8,] >0.5)
Omega_vs <- AlgoResults_vs$Omega_sample_phys[,Omega_i[Omega_i > 2000][1]] %>% relist(skeleton=Model$skeleton)
Omega_vs$eps$local <- 1

AlgoParams$input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
AlgoParams$Np <- 1
AlgoParams$m_CRPS <- 100
impact_sample_vs <- SampleImpact(dir, Model, Omega_vs %>% addTransfParams(), AlgoParams, dat='Train')
saveRDS(impact_sample_vs, paste0(dir, 'IIDIPUS_Results/impact_sample_vs_score_fitted_randomfield'))
df_vs <- flattenImpactSample(impact_sample_vs)

  
ImpactSample_varyingcorrelation <- readRDS(paste0(dir,'Testing_Results/ImpactSample_varyingcorrelation'))

get_univariate_rank <- function(df, log=F){
  z_all <- c()
  for (i in 1:NROW(df)){
    if (!log){
      mat <- cbind(df[i,grep('observed', names(df))], df[i,grep('sampled', names(df))])
    } else {
      mat <- cbind(log(df[i,grep('observed', names(df))]+10), log(df[i,grep('sampled', names(df))]+10))
    }
    z_j <- c()
    z_all <- c(z_all, rank(mat,  ties.method ='random')[1])
  }
  random_runif = (z_all-runif(length(z_all), 0,1))/length(mat)
  return(random_runif)
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

get_univariate_rank2 <- function(df, log=F){
  z_all <- c()
  obs_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    if (length(groupings[[i]])>1) next
    if (!log){
      mat <- cbind(df[groupings[[i]],grep('observed', names(df))], df[groupings[[i]],grep('sampled', names(df))])
    } else {
      mat <- cbind(log(df[groupings[[i]],grep('observed', names(df))]+10), log(df[groupings[[i]],grep('sampled', names(df))]+10))
    }
    for (j in 1:NROW(mat)){
      z_all <- c(z_all, rank(mat[j,],  ties.method ='random')[1])
      obs_all <- c(obs_all, mat[j,1])
    }
  }
  random_runif = (z_all-runif(length(z_all), 0,1))/length(mat)
  hist(random_runif)
}


get_multivariate_ranking <- function(df, log=F){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    if (!log){
      mat <- cbind(df[groupings[[i]],grep('observed', names(df))], df[groupings[[i]],grep('sampled', names(df))])
    } else {
      mat <- cbind(log(df[groupings[[i]],grep('observed', names(df))]+10), log(df[groupings[[i]],grep('sampled', names(df))]+10))
    }
    z_j <- c()
    for (j in 1:NCOL(mat)){
      z_j <- c(z_j, sum(colSums(mat[,j]<= mat[,-j])==NROW(mat)))
    }
    z_all <- c(z_all, rank(z_j,  ties.method ='random')[1])
  }
  random_runif = (z_all-runif(length(z_all), 0,1))/length(z_j)

  return(random_runif)
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

get_average_rank <- function(df, log=F){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    if(length(groupings[[i]])>1) next
    if (!log){
      mat <- cbind(df[groupings[[i]],grep('observed', names(df))], df[groupings[[i]],grep('sampled', names(df))])
    } else {
      mat <- cbind(log(df[groupings[[i]],grep('observed', names(df))]+10), log(df[groupings[[i]],grep('sampled', names(df))]+10))
    }
    mat = mat + runif(length(mat),-0.5,0.5)
    pre_ranks <- apply(apply(mat, 1, rank), 1, mean)
    z_all <- c(z_all, rank(pre_ranks,  ties.method ='random')[1])
  }
  random_runif <-((z_all-runif(length(z_all),0,1))/length(pre_ranks))
  return(random_runif)
  AndersonDarlingTest(random_runif, null='punif')
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

investigate_rankings <- function(){
  groupings_high = which(names(groupings) %in% names(groupings[which(random_runif > 0.9)]))
  groupings[groupings_high]
  i = groupings_high[13]
  mat <- cbind(df[groupings[[i]],5], df[groupings[[i]],7:NCOL(df)])
  mat = mat + runif(length(mat),-0.5,0.5)
  pre_ranks <- apply(apply(mat, 1, rank), 1, mean)
  rank(pre_ranks,  ties.method ='random')[1]
  plot(as.numeric(mat[1,]), as.numeric(mat[2,]))
  points(as.numeric(mat[1,1]), as.numeric(mat[2,1]), col='red',pch=19)
  groupings[groupings_high]
}

get_mst_rank <- function(df, log=F){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    print(i)
    #if(length(groupings[[i]])>1) next
    if (!log){
      mat <- cbind(df[groupings[[i]],grep('observed', names(df))], df[groupings[[i]],grep('sampled', names(df))])
    } else {
      mat <- cbind(log(df[groupings[[i]],grep('observed', names(df))]+10), log(df[groupings[[i]],grep('sampled', names(df))]+10))
    }
    w_mat <- matrix(unlist(AlgoParams$impact_weights[df$impact]) %*% t(unlist(AlgoParams$impact_weights[df$impact])), ncol=NROW(df))
    mat_weighted <- sweep(mat, 1, unlist(AlgoParams$impact_weights[df[groupings[[i]], 'impact']]), '*')
    
    pre_ranks <- c()
    for (j in 1:NCOL(mat)){
      pre_ranks <- c(pre_ranks, sum(spantree(dist(t(mat_weighted[,-j])))$dist))
    }
    z_all <- c(z_all, rank(pre_ranks,  ties.method ='random')[1])
  }
  random_runif <-((z_all-runif(length(z_all),0,1))/length(pre_ranks))
  print(AndersonDarlingTest(random_runif, null='punif')$statistic)
  return(random_runif)
  AndersonDarlingTest(random_runif, null='punif')
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

get_banddepth_rank <- function(df, log=F){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    print(i)
    #if(length(groupings[[i]])==1) next
    if (!log){
      mat <- cbind(df[groupings[[i]],grep('observed', names(df))], df[groupings[[i]],grep('sampled', names(df))])
    } else {
      mat <- cbind(log(df[groupings[[i]],grep('observed', names(df))]+10), log(df[groupings[[i]],grep('sampled', names(df))]+10))
    }
    mat %<>% as.matrix()
    m <- NCOL(mat)
    d <- NROW(mat)
    pre_ranks <- c()
    for (j in 1:m){
      sum <- 0
      for (k in 1:d){
        rank_k <- sum(mat[k,j] <= mat[k,])
        sum <- sum + rank_k * (m - rank_k) + (rank_k-1) * sum(mat[d,j]==mat[k,])
      }
      pre_ranks <- c(pre_ranks, sum /d)
    }
    z_all <- c(z_all, rank(pre_ranks,  ties.method ='random')[1])
  }
  random_runif <-((z_all-runif(length(z_all),0,1))/length(pre_ranks))
  print(AndersonDarlingTest(random_runif, null='punif')$statistic * 0.05)
  return(random_runif)
  AndersonDarlingTest(random_runif, null='punif')
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}


df_low1 <- ImpactSample_varyingcorrelation$df_low1
df_low2 <- ImpactSample_varyingcorrelation$df_low2
df_low3 <- ImpactSample_varyingcorrelation$df_low3
df_high1 <- ImpactSample_varyingcorrelation$df_high1
df_high2 <- ImpactSample_varyingcorrelation$df_high2
df_high3 <- ImpactSample_varyingcorrelation$df_high3
df_true1 <- ImpactSample_varyingcorrelation$df_true1
df_true2 <- ImpactSample_varyingcorrelation$df_true2
df_true3 <- ImpactSample_varyingcorrelation$df_true3

groupings <- split(seq_along(df_true1$event_id), df_true1$event_id)
for (i in 1:length(groupings)){
  plot(unlist(df_true2[groupings[[i]][1], 6:55]), unlist(df_true2[groupings[[i]][2], 6:55]))
  points(unlist(df_high2[groupings[[i]][1], 6:55]), unlist(df_high2[groupings[[i]][2], 6:55]), col='red')
  points(unlist(df_low2[groupings[[i]][1], 6:55]), unlist(df_low2[groupings[[i]][2], 6:55]), col='blue')
  points(unlist(df_true2[groupings[[i]][1], 5]), unlist(df_true2[groupings[[i]][2], 5]), col='green', pch=19)
}
par(mfrow=c(1,3))

for (rank_func in c('get_multivariate_ranking', 'get_average_rank', 'get_mst_rank', 'get_banddepth_rank')){
  log_flag = T
  ranks_low1 <- get(rank_func)(df_low1, log=log_flag)
  ranks_low2 <- get(rank_func)(df_low2, log=log_flag)
  ranks_low3 <- get(rank_func)(df_low3, log=log_flag)
  hist(ranks_low1); hist(ranks_low2); hist(ranks_low3);
  
  ranks_true1 <- get(rank_func)(df_true1, log=log_flag)
  ranks_true2 <- get(rank_func)(df_true2, log=log_flag)
  ranks_true3 <- get(rank_func)(df_true3, log=log_flag)
  hist(ranks_true1); hist(ranks_true2); hist(ranks_true3);
  
  ranks_high1 <- get(rank_func)(df_high1, log=log_flag)
  ranks_high2 <- get(rank_func)(df_high2, log=log_flag)
  ranks_high3 <- get(rank_func)(df_high3, log=log_flag)
  hist(ranks_high1); hist(ranks_high2); hist(ranks_high3);
  
}

#10 x 4, MST_rankhist
hist(ranks_low1, xlab='MST Rank Histogram', ylab='Density', freq=F, main='Under correlated joint distribution')
hist(ranks_true1, xlab='MST Rank Histogram', ylab='Density', freq=F, main='Correct correlation of joint distribution')
hist(ranks_high1, xlab='MST Rank Histogram', ylab='Density', freq=F, main='Over correlated joint distribution')

#### test using simulated data:
n_entries <- 900
cor_true <- 0.5
df_test <- data.frame(event_id = rep(1:(n_entries/3), each=3),
                      iso3 = 'ABC',
                      polygon = 1,
                      impact = rep(c('mortality', 'displacement', 'buildDam'), n_entries/3),
                      observed = c(t(rmvnorm(n_entries/3, c(0,0,0), cbind(c(1,cor_true,cor_true), c(cor_true,1,cor_true), c(cor_true,cor_true,1))))))
cor_sim <- 0.3
for (i in 1:60){
  df_test %<>% cbind(c(t(rmvnorm(n_entries/3, c(0,0,0), cbind(c(1,cor_sim,cor_sim), c(cor_sim,1,cor_sim), c(cor_sim,cor_sim,1))))))
}
names(df_test)[6:65] <- paste0('sampled', 1:60)

for (rank_func in c('get_multivariate_ranking', 'get_average_rank', 'get_mst_rank', 'get_banddepth_rank')){
  ranks_test <- get(rank_func)(df_test, log=F)
  hist(ranks_test)
}

#### Average Rank + MST:
n_repeats <- 3
dists <- array(0, dim=c(n_repeats, 8, 3))
Omega <- AlgoResults$Omega_sample_phys[489,,80] %>% relist(skeleton=Model$skeleton)
Omega$eps$hazard_cor <- 0.6
Omega_true <- Omega
Omega_trialled <- list(Omega_true, Omega_true, Omega_true)
Omega_trialled[[2]]$eps$hazard_cor <- 0.05
Omega_trialled[[3]]$eps$hazard_cor <- 0.95
set.seed(1)

AlgoParams$input_folder = 'IIDIPUS_SimInput_RealAgg5/'
for (i in 1:n_repeats){
  #simulateDataSet(150, Omega, Model, dir)
  simulateRealData("IIDIPUS_Input_RealAgg5/ODDobjects/", Omega_true, Model, dir)
  for (j in 1:3){
    impact_sample <- SampleImpact(dir, Model, Omega_trialled[[j]] %>% addTransfParams(), AlgoParams, dat='Train')
    dists[i,,j] <- CalcDist(impact_sample, AlgoParams)
  }
  directory_path <- paste0(dir, 'IIDIPUS_SimInput_RealAgg5')
  # # List all files in the directory (excluding directories)
  files_to_delete <- list.files(directory_path, full.names = TRUE, recursive = TRUE)
  # # Filter out directories from the list
  files_to_delete <- files_to_delete[!file.info(files_to_delete)$isdir]
  # # Remove the files
  file.remove(files_to_delete)
}
# Energy score favours undercorrelated? Not enough to really tell
# AD with average definitely punishes undercorrelated
# MST favours true then undercorrelated then punishes overcorrelated
# chi squared performs poorly, cannot identify under from correctly correlated
# cvm maybe does slightly better?
plot(rep(1:3, each=3), as.numeric(dists[,3,]), col=c(rep(c('red', 'blue', 'green'),3)), pch=19)
plot(rep(1:3, each=3), as.numeric(dists[,2,]), col=c(rep(c('red', 'blue', 'green'),3)), pch=19)


, ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l', ylim=c(0,1))
for (i in 2:5){
  lines(c(1,3), dists[i,3,c(1,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l')
}


chi_store <- c()
chi_store2 <- c()
for (j in 1:100){
  ranks_std_mst <- runif(100)
  bin_counts <- table(cut(ranks_std_mst, seq(0, 1, 0.1), include.lowest = TRUE, right = FALSE))
  chi_store <- c(chi_store, sum((bin_counts-length(ranks_std_mst)/10)^2)/(length(ranks_std_mst)/10)) #0 #ks.test(ranks_std_average, y='punif')$p.value
  bin_counts <- table(cut(ranks_std_mst, seq(0, 1, 0.02), include.lowest = TRUE, right = FALSE))
  chi_store2 <- c(chi_store2, sum((bin_counts-length(ranks_std_mst)/50)^2)/(length(ranks_std_mst)/50)) #0 #ks.test(ranks_std_average, y='punif')$p.value
}




#MST Plot:
plot(c(1,3), dists[1,3,c(1,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l', ylim=c(0,1))
for (i in 2:5){
  lines(c(1,3), dists[i,3,c(1,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l')
}
plot(c(1,2,3), dists[1,2,c(1,2,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l', ylim=c(0,1))
for (i in 2:5){
  lines(c(1,2,3), dists[i,2,c(1,2,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l')
}

plot(rep(1:3, each=5), dists[1:5,2,], ylab='p.value for Anderson Darling Test', col=rep(1:5, 3), main='Average Ranking')


dists <- array(0, dim=c(5, 8))
for (i in 1:NROW(dists)){
  impact_sample <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams)
  dists[i,] <- CalcDist(impact_sample, AlgoParams)
}
#-----------------------------------------------------------------------------------------------------------------
#--------------Explore correlation in model with different parameter configurations-------------------------------
#-----------------------------------------------------------------------------------------------------------------

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_RealAgg3/ODDobjects/Train/EQ20191029PHL_125')
ODDy@impact %<>% add_row(ODDy@impact[6,] %>% replace(which(names(ODDy@impact[6,])==c('impact')), 'mortality'))
ODDy@impact %<>% add_row(ODDy@impact[6,] %>% replace(which(names(ODDy@impact[6,])==c('impact')), 'buildDam'))

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_SimInput_RealAgg5/ODDobjects/Train/EQ20191029PHL_125')

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_RealAgg5/ODDobjects/Train/EQ20201229HRV_151')

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_RealAgg5/ODDobjects/Train/EQ20170908MEX_67')

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Test/EQ20190620ABC_5148')
ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Test/EQ20191016ABC_6061')
ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Train/EQ20140803ABC_-22')



#Omega$eps$local <- 2 
#Omega$eps$hazard_disp <- 1
#Omega$Lambda1$nu <- 4.3000000001
#Omega$Lambda1$kappa <- 0.1
#Omegatransf <- Omega %>% addTransfParams()
#Omegatransf$Lambda1$loc
#plot(seq(1,20,0.01), plnorm(seq(1,20,0.01), 2.7, 0.2))

#Omega$Lambda1$nu <- 6
#Omega$Lambda1$kappa <- 0.01
Omega$theta$theta1 <- 0.01
Omega$Lambda1$nu <- 8.3
Omega$Lambda1$kappa <- 2
Omega$eps$hazard_disp <- 1.1
Omegatransf <- Omega %>% addTransfParams()
plot(seq(4.5, 10,0.1), pnorm(h_0(seq(4.5, 10,0.1), 4.5, Omegatransf), Omegatransf$Lambda1$loc, Omegatransf$Lambda1$scale), ylim=c(0,1))

#Omega$eps$hazard_mort <- 0.001
#Omega$eps$hazard_disp <- 0.25
#Omega$eps$hazard_bd <- 0.1
#Omega$eps$hazard_cor <- 0.9
impactODD <- DispX(ODDy, Omega %>% addTransfParams(), Model$center, 
                   AlgoParams  %>% replace(which(names(AlgoParams)==c('Np')), 60), output='SampledAgg')
sampled_impact_flattened <- matrix(unlist(lapply(impactODD, function(x){x$sampled})), nrow=length(impactODD), byrow=T)

#investigate two different impact types, total aggregation
obs_of_interest <- c(1,2)
plot(log(sampled_impact_flattened[,obs_of_interest]+10))
points(log(impactODD[[1]]$observed[obs_of_interest[1]]+10), log(impactODD[[1]]$observed[obs_of_interest[2]]+10), col='red', pch=19)

#investigate two different spatial polygons, same impact type
obs_of_interest <- c(6,8)
plot(log(sampled_impact_flattened[,obs_of_interest]+10))
points(log(impactODD[[1]]$observed[obs_of_interest[1]]+10), log(impactODD[[1]]$observed[obs_of_interest[2]]+10), col='red', pch=19)


#eps_local needs to be unattached to hazard-wide error so that it can be larger and permit (lack of) correlation between regions

impact_sample <- SampleImpact(dir, Model, Omega_best %>% addTransfParams(), AlgoParams)
df_impact <- flattenImpactSample(impact_sample) #%>% filter(impact=='mortality')
mat_impact <- df_impact[,5:65]

quants <- apply(mat_impact, 1, sample_quant)
median_pred <- apply(mat_impact, 1, median)
plot(log(median_pred+10), quants)

es_store <- c()
pre_ranks_average <- c()
pre_ranks_mst <- c()
i_store <- c()
grouped_events <- split(seq_along(df_impact$event_id), df_impact$event_id)


#df_impact$group_key <- paste(df_impact$impact, df_impact$event_id, sep = "_")
#grouped_events <- split(seq_along(df_impact$group_key), df_impact$group_key)

for (i in 1:length(grouped_events)){

  obs <- log(df_impact$observed[grouped_events[[i]]]+AlgoParams$log_offset) * unlist(AlgoParams$impact_weights[df_impact$impact[grouped_events[[i]]]])
  sims <- as.matrix(log(df_impact[grouped_events[[i]], 6:65]+AlgoParams$log_offset)) * unlist(AlgoParams$impact_weights[df_impact$impact[grouped_events[[i]]]])
  if (length(obs )> 1){
    i_store <- c(i_store, i)
    es_store<- c(es_store, es_sample(obs, sims))
    pre_ranks_average <- c(pre_ranks_average, get_average_rank_single(cbind(obs, sims)))
    pre_ranks_mst <- c(pre_ranks_mst, get_mst_rank_single(cbind(obs,sims)))  
  }
  #vs_store <- c(vs_store, vs_sample(obs,sims))
  
  #mrh_store <- c(mrh_store, mrh_calc(cbind(obs, sims)))
  #crps_store <- c(crps_store, crps_sample(log(observed[i]), log(samples_combined[i,])))
  #crps_store <- c(crps_store, es_sample(c(log(observed[i]), log(observed[i+200]),log(observed[i+400])), log(samples_combined[c(i, i+200, i+400),])))
  #crps_store <- c(crps_store, crps_sample(log(observed[i]), log(samples_combined[i,])))
}
obs_of_interest <- c(1,2)
plot(sims[obs_of_interest[1],], sims[obs_of_interest[2],])
points(obs[obs_of_interest[1]], obs[obs_of_interest[2]], col='red', pch=19)
vs_sample(obs[obs_of_interest],sims[obs_of_interest,])

grouped_events[i_store[which(pre_ranks_mst<3)]]

lapply(grouped_events, function(x) df_impact$observed[x])
which(names(grouped_events)=='6')
hist((pre_ranks_average-runif(length(pre_ranks_average),0,1))/(AlgoParams$m_CRPS + 1))
hist((pre_ranks_mst-runif(length(pre_ranks_mst),0,1))/(AlgoParams$m_CRPS + 1))

  ranks_std_average <- (pre_ranks_average-runif(length(pre_ranks_average),0,1))/(AlgoParams$m_CRPS + 1)
  ranks_std_mst <- (pre_ranks_mst-runif(length(pre_ranks_mst),0,1))/(AlgoParams$m_CRPS + 1)
  dist_poly[n,1] <- mean(es_store)
  AndersonDarlingTest(ranks_std_average, null='punif')$p.value
  AndersonDarlingTest(ranks_std_mst, null='punif')$p.value

plot(quants) #not too badly calibrated, tendency to overpredict rather than underpredict


#### Investigate what are the problems with the sampled impact that's causing issues with MST histogram #######

# AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-06-24_121215_alpha0.95_M60_Npart990RealAgg3')
# 
# s_change <- 70
# plot(AlgoResults$d_full[,1,7,10], AlgoResults$d[,1,10])

#AlgoResults$d[,,s_change] <- AlgoResults$d_full[,1, 7,s_change]
impact_sample <- imp_best
impact_sample$poly <- lapply(impact_sample$poly, function(x) return(x %>% filter(train_flag=='Train')))

observed <- impact_sample$poly[[1]]$observed
dist_poly <- array(NA, dim=c(AlgoParams$Np,7))
impact_type <- impact_sample$poly[[1]]$impact
impact_weightings <- unlist(AlgoParams$impact_weights[impact_type])
event_id <- impact_sample$poly[[1]]$event_id
grouped_events <- split(seq_along(event_id), event_id)
i_low_mst <- c()
for(n in 1:AlgoParams$Np){
  samples_allocated <- ((n-1)*AlgoParams$m_CRPS+1):(n*AlgoParams$m_CRPS)
  samples_combined <- sapply(impact_sample$poly[samples_allocated], function(x){x$sampled}) #doesn't work if samples_allocated is length 1
  #medians <- apply(samples_combined, 1, mean)
  dist_poly[n,1] <- 0#mean((log(medians[which(impact_type=='mortality')]+10)-log(observed[which(impact_type=='mortality')]+10))^2) * unlist(AlgoParams$impact_weights['mortality'])
  dist_poly[n,2] <- 0#mean((log(medians[which(impact_type=='displacement')]+10)-log(observed[which(impact_type=='displacement')]+10))^2) * unlist(AlgoParams$impact_weights['displacement'])
  dist_poly[n,3] <- 0#mean((log(medians[which(impact_type=='buildDam')]+10)-log(observed[which(impact_type=='buildDam')]+10))^2) * unlist(AlgoParams$impact_weights['buildDam'])
  
  #quants <- (apply(cbind(observed,samples_combined), 1, sample_quant)-runif(length(observed),0,1))/(NCOL(samples_combined)+1)
  #AD_mort <- AndersonDarlingTest(quants[impact_type=='mortality'],null='punif')$statistic
  #dist_poly[n,4] <- AD_mort * unlist(AlgoParams$impact_weights['mortality'])
  es_store <- c()
  pre_ranks_average <- c()
  pre_ranks_mst <- c()
  #vs_store <- c()
  #mrh_store <- c()
  for (i in 1:length(grouped_events)){
    #For each event, compute the energy score of the observed data vs the 'prediction' (simulated data)
    #Each impact type is weighted differently, simply multiplying the observation and the simulations by this weight performs the weighting
    obs <- log(observed[grouped_events[[i]]]+AlgoParams$log_offset) *impact_weightings[grouped_events[[i]]]
    sims <- log(samples_combined[grouped_events[[i]],]+AlgoParams$log_offset) * impact_weightings[grouped_events[[i]]]
    #obs <- log(observed[grouped_events[[i]]]+AlgoParams$log_offset)*impact_weightings[grouped_events[[i]]]
    #sims <- log(samples_combined[grouped_events[[i]],]+AlgoParams$log_offset) * impact_weightings[grouped_events[[i]]]
    
    if (length(grouped_events[[i]])==1){
      #LOOSEEND: Double check that crps_sample is in fact the same as 
      es_store<- c(es_store, crps_sample(obs, sims))
      #mrh_store <- c(mrh_store, mrh_calc(cbind(obs, sims)))
      next
    } 
    #es_store<- c(es_store, vs_sample(obs, sims, w_vs = matrix(impact_weightings[grouped_events[[i]]] %*% t(impact_weightings[grouped_events[[i]]]), ncol=length(grouped_events[[i]]))))
    es_store<- c(es_store, es_sample(obs, sims))
    pre_ranks_average <- c(pre_ranks_average, get_average_rank_single(cbind(obs, sims)))
    pre_ranks_mst <- c(pre_ranks_mst, get_mst_rank_single(cbind(obs,sims)))
    #vs_store <- c(vs_store, vs_sample(obs,sims))
    if (pre_ranks_mst[length(pre_ranks_mst)] < 3){
      i_low_mst <- c(i_low_mst, i)
    }
    
    #mrh_store <- c(mrh_store, mrh_calc(cbind(obs, sims)))
    #crps_store <- c(crps_store, crps_sample(log(observed[i]), log(samples_combined[i,])))
    #crps_store <- c(crps_store, es_sample(c(log(observed[i]), log(observed[i+200]),log(observed[i+400])), log(samples_combined[c(i, i+200, i+400),])))
    #crps_store <- c(crps_store, crps_sample(log(observed[i]), log(samples_combined[i,])))
  }
  #logscores <- ifelse(is.finite(logscores), logscores, 600)
  ranks_std_average <- (pre_ranks_average-runif(length(pre_ranks_average),0,1))/(AlgoParams$m_CRPS + 1)
  ranks_std_mst <- (pre_ranks_mst-runif(length(pre_ranks_mst),0,1))/(AlgoParams$m_CRPS + 1)
  dist_poly[n,1] <- mean(es_store) #mean(crps_store[which(impact_type=='mortality')]) * unlist(AlgoParams$impact_weights['mortality'])
  dist_poly[n,2] <- 0.2*(1 - AndersonDarlingTest(ranks_std_average, null='punif')$p.value) #mean(vs_store) #0.5*AndersonDarlingTest(mrh_store, null='punif')$statistic #mean(crps_store[which(impact_type=='displacement')]) * unlist(AlgoParams$impact_weights['displacement'])
  dist_poly[n,3] <- 0.2*(1 - AndersonDarlingTest(ranks_std_mst, null='punif')$p.value) #mean(crps_store[which(impact_type=='buildDam')]) * unlist(AlgoParams$impact_weights['buildDam'])
  
  bin_counts <- table(cut(ranks_std_mst, seq(0, 1, 0.1), include.lowest = TRUE, right = FALSE))
  #sum((bin_counts-length(ranks_std_mst)/10)^2)/(length(ranks_std_mst)/10)
  dist_poly[n,4] <- 0#sum((bin_counts-length(ranks_std_mst)/10)^2)/(length(ranks_std_mst)/10) #0 #ks.test(ranks_std_average, y='punif')$p.value
  dist_poly[n,5] <- 0#gof.uniform(ranks_std_mst)$Usq #ks.test(ranks_std_mst, y='punif')$p.value
  dist_poly[n,6] <- 0#gof.uniform(ranks_std_mst)$Usq.pvalue #AndersonDarlingTest(ranks_std_average, null='punif')$statistic
  dist_poly[n,7] <- 0#AndersonDarlingTest(ranks_std_mst, null='punif')$statistic
  #logscores  %>% mean()
  #dist_poly[n,4] <- log(ifelse(AD_mort < 2, 2, AD_mort)+1) * unlist(AlgoParams$impact_weights['mortality'])
  #AD_disp <- AndersonDarlingTest(quants[impact_type=='displacement'], null='punif')$statistic
  #dist_poly[n,5] <- AD_disp * unlist(AlgoParams$impact_weights['displacement'])
  
  #AD_bd <- AndersonDarlingTest(quants[impact_type=='buildDam'], null='punif')$statistic
  #dist_poly[n,6] <- AD_bd * unlist(AlgoParams$impact_weights['buildDam'])
  
  #AD_mort_nonzero <- AndersonDarlingTest(quants[impact_type=='mortality' & medians != 0], null='punif')$statistic
  #dist_poly[n,7] <- AD_mort_nonzero * unlist(AlgoParams$impact_weights['mortality']) #ifelse(rbinom(1, 1, P_unif_test(AD_mort_nonzero))==1, 0, AD_mort_nonzero)
  #dist_poly[n,7] <- log(ifelse(AD_mort_nonzero < 2, 2, AD_mort_nonzero)+1) * unlist(AlgoParams$impact_weights['mortality'])
  #dist_poly[n,7] <- 0#ifelse(is.na(dist_poly[n,7]), 50 * unlist(AlgoParams$impact_weights['mortality']), dist_poly[n,7])
  #  dist_poly[n,j] <- ifelse(is.na(dist_poly[n,j]), 0, dist_poly[n,j])
  #}
}

plot_obs_vs_sims <- function(i){
  #dev.off()
  obs <- log(observed[grouped_events[[i]]]+AlgoParams$log_offset) * impact_weightings[grouped_events[[i]]]
  sims <- log(samples_combined[grouped_events[[i]],]+AlgoParams$log_offset) * impact_weightings[grouped_events[[i]]]
  polygon_id = impact_sample$poly[[1]]$polygon[grouped_events[[i]]]
  impact_types = impact_type[grouped_events[[i]]]
  obs <- obs[which(impact_types=='mortality')]
  sims <- sims[which(impact_types=='mortality'),, drop=F]
  polygon_id <- polygon_id[which(impact_types=='mortality')]
  impact_types <- impact_types[which(impact_types=='mortality')]
  if (length(obs) < 8){
    par(mfrow=c(length(obs), length(obs)))
    for (j in 1:length(obs)){
      for (k in 1:length(obs)){
        xlim=range(c(obs[j], sims[j,]))
        ylim=range(c(obs[k], sims[k,]))
        plot(sims[j,], sims[k,], xlim=xlim, ylim=ylim, xlab = paste(impact_types[j], polygon_id[j]), ylab= paste(impact_types[k],polygon_id[k]))
        points(obs[j], obs[k], col='red', pch=4)
      }
    }
  }
  if (length(obs) > 8){
    par(mfrow=c(4, 4), mai = c(0.5, 0.5, 0.5, 0.5))
    for (j in 1:4){
      for (k in 1:4){
        xlim=range(c(obs[j], sims[j,]))
        ylim=range(c(obs[k], sims[k,]))
        plot(sims[j,], sims[k,], xlim=xlim, ylim=ylim, xlab = paste(impact_types[j], polygon_id[j]), ylab= paste(impact_types[k],polygon_id[k]))
        points(obs[j], obs[k], col='red', pch=4)
      }
    }
  }
}

i_low_mst

# Weak ones: 37, 45 slightly, 90 slightly, 93, 109 well off.
# Strong ones: 112
plot_obs_vs_sims(37)




ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_RealAgg5/ODDobjects/Train/EQ20170908MEX_67')
sampled_mat <- array(0, dim=c(100, 15))
for (i in 1:100){
  print(i)
  Omega_samp <-  AlgoResults$Omega_sample_phys[sample(1:1000,1),,130] %>% relist(skeleton=Model$skeleton)
  sampled <- DispX(ODDy, Omega_samp %>% addTransfParams(),center=Model$center, 
                   AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                   output='SampledAgg')
  sampled_mat[i, ] <- sampled[[1]]$sampled
}
obs <- sampled[[1]]$observed

obs_of_interest <- c(1,2)
plot(log(sampled_mat[,c(obs_of_interest[1], obs_of_interest[2])]+10), xlim=log(range(sampled_mat[,obs_of_interest[1]], obs[obs_of_interest[1]])+10), ylim=log(range(sampled_mat[,obs_of_interest[2]], obs[obs_of_interest[2]])+10),
     xlab='MEX67Displacement', ylab='GTM67Displacement')
points(log(obs[obs_of_interest[1]]+10), log(obs[obs_of_interest[2]]+10), col='red', pch=19)




#----------------------------------------------------------------------------------------------------


flattenImpactSample <- function(impact_sample){
  df <- data.frame(event_id = impact_sample$poly[[1]]$event_id,
                   iso3 = impact_sample$poly[[1]]$iso3,
                   polygon = impact_sample$poly[[1]]$polygon,
                   impact = impact_sample$poly[[1]]$impact,
                   observed = impact_sample$poly[[1]]$observed)
  for (j in 1:length(impact_sample$poly)){
    df %<>% cbind(impact_sample$poly[[j]]$sampled)
  }
  colnames(df)[6:NCOL(df)] <- paste0('sampled.', 1:length(impact_sample$poly))
  return(df)
}

# for (i in 1:5){
#   Omegatrue_high_cor <- Omega_true
#   Omegatrue_high_cor$eps$hazard_cor <- 0.95
#   Omegatrue_low_cor <- Omega_true
#   Omegatrue_low_cor$eps$hazard_cor <- 0.05
#   impact_sample_high_cor <- SampleImpact(dir, Model, Omegatrue_high_cor %>% addTransfParams(), AlgoParams, dat='all')
#   dist_high <- CalcDist(impact_sample_high_cor, AlgoParams)
#   print(paste('high',dist_high))
#   df_high3 <- flattenImpactSample(impact_sample_high_cor)
#   
#   impact_sample_low_cor <- SampleImpact(dir, Model, Omegatrue_low_cor %>% addTransfParams(), AlgoParams, dat='all')
#   dist_low <- CalcDist(impact_sample_low_cor, AlgoParams)
#   print(paste('low',dist_low))
#   df_low3 <- flattenImpactSample(impact_sample_low_cor)
#   
#   impact_sample_true_cor <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams, dat='all')
#   dist_true <- CalcDist(impact_sample_true_cor, AlgoParams)
#   print(paste('true', dist_true))
#   df_true3 <- flattenImpactSample(impact_sample_true_cor)
#   
# }

# saveRDS(list(
#   df_high1 = df_high1,
#   df_high2 = df_high2,
#   df_high3 = df_high3,
#   df_low1 = df_low1,
#   df_low2 = df_low2,
#   df_low3 = df_low3,
#   df_true1 = df_true1,
#   df_true2 = df_true2,
#   df_true3 = df_true3
# ), 'ImpactSample_varyingcorrelation')


# Omegatrue_high_locerror <- Omega_true
# Omegatrue_high_locerror$eps$local <- 0.1
# Omegatrue_low_locerror <- Omega_true
# Omegatrue_low_locerror$eps$local <- 1.5
# impact_sample_high_locerror <- SampleImpact(dir, Model, Omegatrue_high_locerror %>% addTransfParams(), AlgoParams, dat='all')
# dist_high <- CalcDist(impact_sample_high_locerror, AlgoParams)
# print(paste('high',dist_high))
# df_highloc1 <- flattenImpactSample(impact_sample_high_locerror)
# 
# impact_sample_low_locerror <- SampleImpact(dir, Model, Omegatrue_low_locerror %>% addTransfParams(), AlgoParams, dat='all')
# dist_low <- CalcDist(impact_sample_low_locerror, AlgoParams)
# print(paste('low',dist_low))
# df_lowloc1 <- flattenImpactSample(impact_sample_low_locerror)
# 
# impact_sample_true_locerror <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams, dat='all')
# dist_true <- CalcDist(impact_sample_true_locerror, AlgoParams)
# print(paste('true', dist_true))
# df_trueloc1 <- flattenImpactSample(impact_sample_true_locerror)

ImpactSample_varyingcorrelation <- readRDS(paste0(dir,'ImpactSample_varyingcorrelation'))

get_multivariate_ranking <- function(df, log=F){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    if (!log){
      mat <- cbind(df[groupings[[i]],5], df[groupings[[i]],6:55])
    } else {
      mat <- cbind(log(df[groupings[[i]],5]+10), log(df[groupings[[i]],6:55]+10))
    }
    z_j <- c()
    for (j in 1:NCOL(mat)){
      z_j <- c(z_j, sum(colSums(mat[,j]<= mat[,-j])==NROW(mat)))
    }
    z_all <- c(z_all, rank(z_j,  ties.method ='random')[1])
  }
  random_runif = (z_all-runif(length(z_all), 0,1))/length(z_j)
  return(random_runif)
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

get_average_rank <- function(df, log=F){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    if (!log){
      mat <- cbind(df[groupings[[i]],5], df[groupings[[i]],6:55])
    } else {
      mat <- cbind(log(df[groupings[[i]],5]+10), log(df[groupings[[i]],6:55]+10))
    }
    pre_ranks <- apply(apply(mat, 1, rank), 1, mean)
    z_all <- c(z_all, rank(pre_ranks,  ties.method ='random')[1])
  }
  random_runif <-((z_all-runif(length(z_all),0,1))/length(pre_ranks))
  return(random_runif)
  AndersonDarlingTest(random_runif, null='punif')
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

get_mst_rank <- function(df, log=F){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    print(i)
    if (!log){
      mat <- cbind(df[groupings[[i]],5], df[groupings[[i]],6:55])
    } else {
      mat <- cbind(log(df[groupings[[i]],5]+10), log(df[groupings[[i]],6:55]+10))
    }
    w_mat <- matrix(unlist(AlgoParams$impact_weights[df$impact]) %*% t(unlist(AlgoParams$impact_weights[df$impact])), ncol=NROW(df))
    mat_weighted <- sweep(mat, 1, unlist(AlgoParams$impact_weights[df[groupings[[i]], 'impact']]), '*')
    
    pre_ranks <- c()
    for (j in 1:NCOL(mat)){
      pre_ranks <- c(pre_ranks, sum(spantree(dist(t(mat_weighted[,-j])))$dist))
    }
    z_all <- c(z_all, rank(pre_ranks,  ties.method ='random')[1])
  }
  random_runif <-((z_all-runif(length(z_all),0,1))/length(pre_ranks))
  return(random_runif)
  AndersonDarlingTest(random_runif, null='punif')
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}

get_banddepth_rank <- function(df, log=F){
  z_all <- c()
  groupings <- split(seq_along(df$event_id), df$event_id)
  for (i in 1:length(groupings)){
    print(i)
    if (!log){
      mat <- cbind(df[groupings[[i]],5], df[groupings[[i]],6:55])
    } else {
      mat <- cbind(log(df[groupings[[i]],5]+10), log(df[groupings[[i]],6:55]+10))
    }
    m <- NCOL(mat)
    d <- NROW(mat)
    pre_ranks <- c()
    for (j in 1:m){
      sum <- 0
      for (k in 1:d){
        rank_k <- sum(mat[k,j] <= mat[k,])
        sum <- sum + rank_k * (m - rank_k) + (rank_k-1) * sum(mat[d,j]==mat[k,])
      }
      pre_ranks <- c(pre_ranks, sum /d)
    }
    z_all <- c(z_all, rank(pre_ranks,  ties.method ='random')[1])
  }
  random_runif <-((z_all-runif(length(z_all),0,1))/length(pre_ranks))
  return(random_runif)
  AndersonDarlingTest(random_runif, null='punif')
  return(AndersonDarlingTest((z_all-runif(length(z_all), 0,1))/length(z_j), null='punif')$statistic)
}


df_low1 <- ImpactSample_varyingcorrelation$df_low1
df_low2 <- ImpactSample_varyingcorrelation$df_low2
df_low3 <- ImpactSample_varyingcorrelation$df_low3
df_high1 <- ImpactSample_varyingcorrelation$df_high1
df_high2 <- ImpactSample_varyingcorrelation$df_high2
df_high3 <- ImpactSample_varyingcorrelation$df_high3
df_true1 <- ImpactSample_varyingcorrelation$df_true1
df_true2 <- ImpactSample_varyingcorrelation$df_true2
df_true3 <- ImpactSample_varyingcorrelation$df_true3

groupings <- split(seq_along(df_true1$event_id), df_true1$event_id)
for (i in 1:length(groupings)){
  plot(unlist(df_true2[groupings[[i]][1], 6:55]), unlist(df_true2[groupings[[i]][2], 6:55]))
  points(unlist(df_high2[groupings[[i]][1], 6:55]), unlist(df_high2[groupings[[i]][2], 6:55]), col='red')
  points(unlist(df_low2[groupings[[i]][1], 6:55]), unlist(df_low2[groupings[[i]][2], 6:55]), col='blue')
  points(unlist(df_true2[groupings[[i]][1], 5]), unlist(df_true2[groupings[[i]][2], 5]), col='green', pch=19)
}
par(mfrow=c(1,3))

for (rank_func in c('get_multivariate_ranking', 'get_average_rank', 'get_mst_rank', 'get_banddepth_rank')){
  log_flag = T
  ranks_low1 <- get(rank_func)(df_low1, log=log_flag)
  ranks_low2 <- get(rank_func)(df_low2, log=log_flag)
  ranks_low3 <- get(rank_func)(df_low3, log=log_flag)
  hist(ranks_low1); hist(ranks_low2); hist(ranks_low3);
  
  ranks_true1 <- get(rank_func)(df_true1, log=log_flag)
  ranks_true2 <- get(rank_func)(df_true2, log=log_flag)
  ranks_true3 <- get(rank_func)(df_true3, log=log_flag)
  hist(ranks_true1); hist(ranks_true2); hist(ranks_true3);
  
  ranks_high1 <- get(rank_func)(df_high1, log=log_flag)
  ranks_high2 <- get(rank_func)(df_high2, log=log_flag)
  ranks_high3 <- get(rank_func)(df_high3, log=log_flag)
  hist(ranks_high1); hist(ranks_high2); hist(ranks_high3);
  
}

#### test using simulated data:
n_entries <- 900
cor_true <- 0.5
df_test <- data.frame(event_id = rep(1:(n_entries/3), each=3),
                      iso3 = 'ABC',
                      polygon = 1,
                      impact = rep(c('mortality', 'displacement', 'buildDam'), n_entries/3),
                      observed = c(t(rmvnorm(n_entries/3, c(0,0,0), cbind(c(1,cor_true,cor_true), c(cor_true,1,cor_true), c(cor_true,cor_true,1))))))
cor_sim <- 0.3
for (i in 1:60){
  df_test %<>% cbind(c(t(rmvnorm(n_entries/3, c(0,0,0), cbind(c(1,cor_sim,cor_sim), c(cor_sim,1,cor_sim), c(cor_sim,cor_sim,1))))))
}
names(df_test)[6:65] <- paste0('sampled', 1:60)

for (rank_func in c('get_multivariate_ranking', 'get_average_rank', 'get_mst_rank', 'get_banddepth_rank')){
  ranks_test <- get(rank_func)(df_test, log=F)
  hist(ranks_test)
}

#### Average Rank + MST:
n_repeats <- 8
dists <- array(0, dim=c(n_repeats, 8, 3))
Omega_trialled <- list(Omega_true, Omega_true, Omega_true)
Omega_trialled[[2]]$eps$hazard_cor <- 0.95
Omega_trialled[[3]]$eps$hazard_cor <- 0.05
set.seed(1)
for (i in 1:n_repeats){
  simulateDataSet(150, Omega, Model, dir)
  for (j in 1:3){
    impact_sample <- SampleImpact(dir, Model, Omega_trialled[[j]] %>% addTransfParams(), AlgoParams, dat='Train')
    dists[i,,j] <- CalcDist(impact_sample, AlgoParams)
  }
  directory_path <- paste0(dir, 'IIDIPUS_Input/ODDobjects')
  # List all files in the directory (excluding directories)
  files_to_delete <- list.files(directory_path, full.names = TRUE, recursive = TRUE)
  # Filter out directories from the list
  files_to_delete <- files_to_delete[!file.info(files_to_delete)$isdir]
  # Remove the files
  file.remove(files_to_delete)
}

#MST Plot:
plot(c(1,3), dists[1,3,c(1,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l', ylim=c(0,1))
for (i in 2:5){
  lines(c(1,3), dists[i,3,c(1,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l')
}
plot(c(1,2,3), dists[1,2,c(1,2,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l', ylim=c(0,1))
for (i in 2:5){
  lines(c(1,2,3), dists[i,2,c(1,2,3)], ylab='p.value for Anderson Darling Test',main='Minimum Spanning Tree Ranking', type='l')
}

plot(rep(1:3, each=5), dists[1:5,2,], ylab='p.value for Anderson Darling Test', col=rep(1:5, 3), main='Average Ranking')


dists <- array(0, dim=c(5, 8))
for (i in 1:NROW(dists)){
  impact_sample <- SampleImpact(dir, Model, Omega_true %>% addTransfParams(), AlgoParams)
  dists[i,] <- CalcDist(impact_sample, AlgoParams)
}

# Explore correlation in model with different parameter configurations:
ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_RealAgg3/ODDobjects/Train/EQ20210207PHL_154')
DispX(ODDy, Omega %>% addTransfParams(), Model$center, 
      AlgoParams  %>% replace(which(names(AlgoParams)==c('Np')), 50), output='SampledAgg')


# Which Anderson Darlings to approve?


