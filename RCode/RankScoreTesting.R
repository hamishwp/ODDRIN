
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

  obs <- log(df_impact$observed[grouped_events[[i]]]+AlgoParams$log_offset) * unlist(AlgoParams$kernel_sd[df_impact$impact[grouped_events[[i]]]])
  sims <- as.matrix(log(df_impact[grouped_events[[i]], 6:65]+AlgoParams$log_offset)) * unlist(AlgoParams$kernel_sd[df_impact$impact[grouped_events[[i]]]])
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
impact_weightings <- unlist(AlgoParams$kernel_sd[impact_type])
event_id <- impact_sample$poly[[1]]$event_id
grouped_events <- split(seq_along(event_id), event_id)
i_low_mst <- c()
for(n in 1:AlgoParams$Np){
  samples_allocated <- ((n-1)*AlgoParams$m_CRPS+1):(n*AlgoParams$m_CRPS)
  samples_combined <- sapply(impact_sample$poly[samples_allocated], function(x){x$sampled}) #doesn't work if samples_allocated is length 1
  #medians <- apply(samples_combined, 1, mean)
  dist_poly[n,1] <- 0#mean((log(medians[which(impact_type=='mortality')]+10)-log(observed[which(impact_type=='mortality')]+10))^2) * unlist(AlgoParams$kernel_sd['mortality'])
  dist_poly[n,2] <- 0#mean((log(medians[which(impact_type=='displacement')]+10)-log(observed[which(impact_type=='displacement')]+10))^2) * unlist(AlgoParams$kernel_sd['displacement'])
  dist_poly[n,3] <- 0#mean((log(medians[which(impact_type=='buildDam')]+10)-log(observed[which(impact_type=='buildDam')]+10))^2) * unlist(AlgoParams$kernel_sd['buildDam'])
  
  #quants <- (apply(cbind(observed,samples_combined), 1, sample_quant)-runif(length(observed),0,1))/(NCOL(samples_combined)+1)
  #AD_mort <- AndersonDarlingTest(quants[impact_type=='mortality'],null='punif')$statistic
  #dist_poly[n,4] <- AD_mort * unlist(AlgoParams$kernel_sd['mortality'])
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
  dist_poly[n,1] <- mean(es_store) #mean(crps_store[which(impact_type=='mortality')]) * unlist(AlgoParams$kernel_sd['mortality'])
  dist_poly[n,2] <- 0.2*(1 - AndersonDarlingTest(ranks_std_average, null='punif')$p.value) #mean(vs_store) #0.5*AndersonDarlingTest(mrh_store, null='punif')$statistic #mean(crps_store[which(impact_type=='displacement')]) * unlist(AlgoParams$kernel_sd['displacement'])
  dist_poly[n,3] <- 0.2*(1 - AndersonDarlingTest(ranks_std_mst, null='punif')$p.value) #mean(crps_store[which(impact_type=='buildDam')]) * unlist(AlgoParams$kernel_sd['buildDam'])
  
  bin_counts <- table(cut(ranks_std_mst, seq(0, 1, 0.1), include.lowest = TRUE, right = FALSE))
  #sum((bin_counts-length(ranks_std_mst)/10)^2)/(length(ranks_std_mst)/10)
  dist_poly[n,4] <- 0#sum((bin_counts-length(ranks_std_mst)/10)^2)/(length(ranks_std_mst)/10) #0 #ks.test(ranks_std_average, y='punif')$p.value
  dist_poly[n,5] <- 0#gof.uniform(ranks_std_mst)$Usq #ks.test(ranks_std_mst, y='punif')$p.value
  dist_poly[n,6] <- 0#gof.uniform(ranks_std_mst)$Usq.pvalue #AndersonDarlingTest(ranks_std_average, null='punif')$statistic
  dist_poly[n,7] <- 0#AndersonDarlingTest(ranks_std_mst, null='punif')$statistic
  #logscores  %>% mean()
  #dist_poly[n,4] <- log(ifelse(AD_mort < 2, 2, AD_mort)+1) * unlist(AlgoParams$kernel_sd['mortality'])
  #AD_disp <- AndersonDarlingTest(quants[impact_type=='displacement'], null='punif')$statistic
  #dist_poly[n,5] <- AD_disp * unlist(AlgoParams$kernel_sd['displacement'])
  
  #AD_bd <- AndersonDarlingTest(quants[impact_type=='buildDam'], null='punif')$statistic
  #dist_poly[n,6] <- AD_bd * unlist(AlgoParams$kernel_sd['buildDam'])
  
  #AD_mort_nonzero <- AndersonDarlingTest(quants[impact_type=='mortality' & medians != 0], null='punif')$statistic
  #dist_poly[n,7] <- AD_mort_nonzero * unlist(AlgoParams$kernel_sd['mortality']) #ifelse(rbinom(1, 1, P_unif_test(AD_mort_nonzero))==1, 0, AD_mort_nonzero)
  #dist_poly[n,7] <- log(ifelse(AD_mort_nonzero < 2, 2, AD_mort_nonzero)+1) * unlist(AlgoParams$kernel_sd['mortality'])
  #dist_poly[n,7] <- 0#ifelse(is.na(dist_poly[n,7]), 50 * unlist(AlgoParams$kernel_sd['mortality']), dist_poly[n,7])
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
