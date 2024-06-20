
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
    w_mat <- matrix(unlist(AlgoParams$kernel_sd[df$impact]) %*% t(unlist(AlgoParams$kernel_sd[df$impact])), ncol=NROW(df))
    mat_weighted <- sweep(mat, 1, unlist(AlgoParams$kernel_sd[df[groupings[[i]], 'impact']]), '*')
    
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
#-----------------------------------------------------------------------------------------------------------------
#--------------Explore correlation in model with different parameter configurations-------------------------------
#-----------------------------------------------------------------------------------------------------------------

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_RealAgg3/ODDobjects/Train/EQ20191029PHL_125')

ODDy@impact %<>% add_row(ODDy@impact[6,] %>% replace(which(names(ODDy@impact[6,])==c('impact')), 'mortality'))
ODDy@impact %<>% add_row(ODDy@impact[6,] %>% replace(which(names(ODDy@impact[6,])==c('impact')), 'buildDam'))


Omega$eps$local <- 30
Omega$eps$hazard_mort <- 0.15
Omega$eps$hazard_disp <- 0.25
Omega$eps$hazard_bd <- 0.2
Omega$eps$hazard_cor <- 0.9
impactODD <- DispX(ODDy, Omega %>% addTransfParams(), Model$center, 
                   AlgoParams  %>% replace(which(names(AlgoParams)==c('Np')), 60), output='SampledAgg')
sampled_impact_flattened <- matrix(unlist(lapply(impactODD, function(x){x$sampled})), ncol=NROW(ODDy@impact), byrow=T)

#investigate two different impact types, total aggregation
obs_of_interest <- c(1,2)
plot(log(sampled_impact_flattened[,obs_of_interest]+10))

#investigate two different spatial polygons, same impact type
obs_of_interest <- c(4,5)
plot(log(sampled_impact_flattened[,obs_of_interest]+10))



#eps_local needs to be unattached to hazard-wide error so that it can be larger and permit (lack of) correlation between regions


