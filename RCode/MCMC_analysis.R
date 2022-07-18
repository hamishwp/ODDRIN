

tag <- '2022-07-04_145319'
n_x <- 14
output <- readRDS(paste0(dir, 'AWS_IIDIPUS_Results/output_', tag))

endpoint <- which(is.na(output[,1]))[1] - 1 #Find the iteration that MCMC execution ended

# -----------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------- SINGLE CHAIN --------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------


#Plot likelihood over time
plot(output[1:endpoint,1], type='l', ylab='Likelihood', xlab='Iteration')

#Plot trace plots for a single chain, compare to the true values
par(mfrow=c(4,4))
for(i in 1:n_x){
  ylim=c(min(unlist(Omega)[i], output[1:endpoint,i+1]), max(unlist(Omega)[i], output[1:endpoint,i+1]))
  plot(output[1:endpoint,i+1], type='l', ylab='', ylim=ylim, main=names(unlist(Omega))[i])
  abline(h=unlist(Omega)[i], col='red')
}
par(mfrow=c(1,1))

#Compare maximum a posteri estimate with the true values
Omega_MAP <- output[which.max(output[1:endpoint,1]),2:ncol(output)] %>% 
  relist(skeleton=Model$skeleton)

# Plot S-curves for the actual and MAP parameterisation
plot_S_curves(Omega,Omega_MAP)

# Plot Posterior S-curves for the actual and MAP parameterisation
plot_posterior_S_curves <- function(Omega, output){
  Intensity <- seq(0,10,0.1)
  Dfun<-function(I_ij, theta) h_0(I = I_ij,I0 = 4.5,theta = theta)
  Damage <- Dfun(Intensity, theta=Omega$theta)
  D_extent <- BinR(Damage, Omega$zeta)
  D_MortDisp <-  D_MortDisp_calc(Damage, Omega)
  D_BD <-  plnorm(Damage, Omega$Lambda3$nu, Omega$Lambda3$omega)
  plot(Intensity, D_MortDisp[1,], col='red', type='l', ylim=c(0,1), ylab='Proportion', lwd=3); 
  
  plotted_iter <- (endpoint %/% 2):endpoint
  for (i in plotted_iter){
    Omega_i <- output[i,2:(n_x+1)] %>% relist(skeleton=Model$skeleton)
    Damage_i <- Dfun(Intensity, theta=Omega_i$theta) 
    D_extent_sample <- BinR(Damage_i, Omega_i$zeta)
    D_MortDisp_sample <-  D_MortDisp_calc(Damage_i, Omega_i)
    D_BD_sample <- plnorm( Dfun(Intensity, theta=Omega_i$theta), Omega_i$Lambda3$nu, Omega_i$Lambda3$omega)
    lines(Intensity, D_MortDisp_sample[1,], col=adjustcolor("red", alpha = 0.1)); lines(Intensity, D_MortDisp_sample[2,], col=adjustcolor("cyan", alpha = 0.1)); 
    #lines(Intensity, D_extent_sample, col=adjustcolor("green", alpha = 0.1)); lines(Intensity, D_BD_sample, col=adjustcolor("pink", alpha = 0.1)); 
  }
  lines(Intensity, D_MortDisp[1,], col='darkred', lwd=3); lines(Intensity, D_MortDisp[2,], col='blue', lwd=3); 
  #lines(Intensity, D_BD, col='hotpink', type='l', lwd=3); lines(Intensity, D_extent, col='darkgreen', type='l', lwd=3); 
  legend(x=1,y=0.7, c('D_Mort', 'D_Disp', 'D_BD', 'D_B'), col=c('red','blue','pink', 'green'), lty=1)
}

plot_posterior_S_curves(Omega, output)

#compare likelihood of the true parameters to the MAP estimate
for(i in 1:5){
  logTarget(dir = dir,Model = Model,proposed = Omega, AlgoParams = AlgoParams, epsilon=AlgoParams$epsilon_min)
  logTarget(dir = dir,Model = Model,proposed = Omega_MAP, AlgoParams = AlgoParams)
} 


# -----------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------- MULTIPLE CHAINS -----------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

tags <- c('2022-07-04_131850_laplace', '2022-07-04_145319_normal')
outputs <- list()
param_max <- unlist(Omega)
param_min <- unlist(Omega)
endpoint_max <- 1

#combine all chains into a single list
for (i in 1:length(tags)){
  output_i <- list(
    output = readRDS(paste0(dir, 'AWS_IIDIPUS_Results/output_', tags[i])),
    endpoint = which(is.na(output[,1]))[1] - 1
  )
  endpoint_max <- max(endpoint_max, output_i$endpoint)
  param_max <- apply(rbind(output_i$output[1:output_i$endpoint,2:(n_x+1)], param_max), 2, max, na.rm = TRUE)
  param_min <- apply(rbind(output_i$output[1:output_i$endpoint,2:(n_x+1)], param_min), 2, min, na.rm = TRUE)
  outputs[[i]] <- output_i
}

#Trace plots for multiple chains
par(mfrow=c(4,4))
for(i in 1:n_x){
  plot(outputs[[1]]$output[1:endpoint_max,i+1], type='l', ylab='', ylim=c(param_min[i], param_max[i]), main=names(unlist(Omega))[i])
  if (length(tags) > 1){
    for(j in 2:length(tags)){
      lines(outputs[[j]]$output[1:(outputs[[j]]$endpoint),i+1], type='l', col='blue')
    }
  }
  abline(h=unlist(Omega)[i], col='red')
}
par(mfrow=c(1,1))

# -----------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------ SMC ABC --------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

tag <- '2022-07-16_093045'
output <- readRDS(paste0(dir,"AWS_IIDIPUS_Results/abcsmc_",tag))
s_finish = max(which(colSums(is.na(output$omegastore[,1,]))==0)) #find last completed step

Omega_sample_phys <- output$omegastore[,,1:s_finish]
W <- output$W[,1:s_finish]
d <- output$d[,1:s_finish]
Npart <- NROW(W)

par(mfrow=c(4,4))
for (i in 1:n_x){
  plot(rep(1:s_finish,each=Npart),Omega_sample_phys[,i,1:s_finish], ylim=c(min(Omega_sample_phys[,i,1:s_finish], unlist(Omega)[i]), max(Omega_sample_phys[,i,1:s_finish], unlist(Omega)[i])),
       main=names(unlist(Omega))[i], xlab='step', ylab='')
  abline(h=unlist(Omega)[i], col='red')
}
par(mfrow=c(1,1))

par(mfrow=c(4,4))
for (i in 1:n_x){
  plot(rep(1, Npart), Omega_sample_phys[,i,1], ylim=c(min(Omega_sample_phys[,i,1:s_finish], unlist(Omega)[i]), max(Omega_sample_phys[,i,1:s_finish], unlist(Omega)[i])), xlim=c(1, s_finish),
       main=names(unlist(Omega))[i], xlab='step', ylab='')
  for (s in 2:s_finish){
    retain <- which(d[,s] < quantile(d[,s], 0.6))
    points(rep(s, length(retain)), Omega_sample_phys[retain,i,s], col='blue')
    points(rep(s, Npart - length(retain)), Omega_sample_phys[-(retain),i,s], col='black')
  }
  abline(h=unlist(Omega)[i], col='red')
}
par(mfrow=c(1,1))



plot(rep(1:s_finish, each=Npart), as.vector(d[,1:s_finish]), ylab='Negative Log Likelihood', xlab='step')


# Compare posterior S-curves with the S-curve under the true parameterisation
Intensity <- seq(0,10,0.1)
Dfun<-function(I_ij, theta) h_0(I = I_ij,I0 = 4.5,theta = theta)
Damage <- Dfun(Intensity, theta=Omega$theta)
D_extent <- BinR(Damage, Omega$zeta)
D_MortDisp <-  D_MortDisp_calc(Damage, Omega)
D_BD <-  plnorm(Damage, Omega$Lambda3$nu, Omega$Lambda3$omega)
plot(Intensity, D_MortDisp[1,], col='red', type='l', ylim=c(0,1), ylab='Proportion', lwd=3); 
for (i in 1:Npart){
  Omega_i <- Omega_sample_phys[i,,s_finish] %>% relist(skeleton=Model$skeleton)
  Damage_i <- Dfun(Intensity, theta=Omega_i$theta) 
  D_extent_sample <- BinR(Damage_i, Omega_i$zeta)
  D_MortDisp_sample <-  D_MortDisp_calc(Damage_i, Omega_i)
  D_BD_sample <- plnorm( Dfun(Intensity, theta=Omega_i$theta), Omega_i$Lambda3$nu, Omega_i$Lambda3$omega)
  lines(Intensity, D_MortDisp_sample[1,], col=adjustcolor("red", alpha = 0.1)); lines(Intensity, D_MortDisp_sample[2,], col=adjustcolor("cyan", alpha = 0.1)); 
}
lines(Intensity, D_MortDisp[1,], col='darkred', lwd=3); lines(Intensity, D_MortDisp[2,], col='blue', lwd=3); 
#lines(Intensity, D_BD, col='hotpink', type='l', lwd=3); lines(Intensity, D_extent, col='darkgreen', type='l', lwd=3); 
legend(x=1,y=0.7, c('D_Mort', 'D_Disp', 'D_BD', 'D_B'), col=c('red','blue','pink', 'green'), lty=1)


Omega_MAP <- Omega_sample_phys[which.min(d[,s_finish]),,s_finish] %>% relist(skeleton=Model$skeleton)
#Plot disp, mort and nBD simulations under the true and MAP parameterisations
folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
impact_MAP <- array(NA, dim=c(3, length(ufiles)))
impact_Omega <- array(NA, dim=c(3, length(ufiles)))
impact_true <- array(NA, dim=c(3, length(ufiles)))
for(i in 1:length(ufiles)){
  ODDSim <- readRDS(paste0(folderin,ufiles[i]))
  #simulate displacement, mortality and building destruction using DispX
  ODDSim_Omega <- ODDSim %>% DispX(Omega, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
  ODDSim_MAP <- ODDSim %>% DispX(Omega_MAP, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
  impact_true[1, i] = ODDSim@gmax$gmax
  impact_Omega[1, i] = ODDSim_Omega@predictDisp$disp_predictor
  impact_MAP[1, i] = ODDSim_MAP@predictDisp$disp_predictor
  impact_true[2, i] = ODDSim@gmax$mortality
  impact_Omega[2, i] = ODDSim_Omega@predictDisp$mort_predictor
  impact_MAP[2, i] = ODDSim_MAP@predictDisp$mort_predictor
  impact_true[3, i] =  ODDSim@gmax$buildDestroyed
  impact_Omega[3, i] = ODDSim_Omega@predictDisp$nBD_predictor
  impact_MAP[3, i] = ODDSim_MAP@predictDisp$nBD_predictor
}

par(mfrow=c(1,1))
plot(impact_true[2,],impact_true[2,], xlim=c(1,500000),ylim=c(1,500000))
points(impact_true[2,], impact_Omega[2,], col='blue')
points(impact_true[2,], impact_MAP[2,], col='red')

k <- 10
epsilon = 0.03
log(dlnormTrunc(impact_true[2,36]+k, log(impact_MAP[2,36]+k), sdlog=epsilon, min=k))
log(dlnormTrunc(impact_true[2,36]+k, log(impact_Omega[2,36]+k), sdlog=epsilon, min=k))


# -----------------------------------------------------------------------------------------------------------------------


LLs_disp <- array(NA, c(2, 5))
LLs_tot <- array(NA, c(2, 5))
j = 1
for (Np in c(1)){
  print(j)
  for (i in 1){
    print(i)
    AlgoParams$Np <- Np
    LLs_disp[j, i] = LL_Displacement(0, dir = dir,Model = Model,proposed = Omega, AlgoParams = AlgoParams, epsilon=AlgoParams$epsilon_min)
    LLs_tot[j,i] = LL_Buildings(LLs_disp[j,i], dir = dir,Model = Model,proposed = Omega, AlgoParams = AlgoParams)
  }
  j = j + 1
}

plot(rep(1,5), LLs_disp[1,],xlim=c(0.5,5.5), ylim=c(min(LLs_disp), max(LLs_disp)))
points(rep(2,5), LLs_disp[2,])
points(rep(3,5), LLs_disp[3,])
points(rep(4,5), LLs_disp[4,])
points(rep(5,5), LLs_disp[5,])
range(LLs_tot[5,])
# folderin<-paste0(dir,"IIDIPUS_SimInput/ODDobjects/")
# ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
# mags = c()
# mort = c()
# disp = c()
# bd = c()
# for(i in 1:length(ufiles)){
#   ODDy<-readRDS(paste0(folderin,ufiles[i]))
#   mags = append(mags, max(ODDy@data$hazMean1, na.rm=TRUE))
#   mort = append(mort, ODDy@gmax$mortality/sum(ODDy@data$Population, na.rm=TRUE))
#   disp = append(disp, ODDy@gmax$gmax/sum(ODDy@data$Population, na.rm=TRUE))
#   bd = append(bd, ODDy@gmax$buildDestroyed/sum(ODDy@data$nBD, na.rm=TRUE))
# }


#Magnitude vs contribution to likelihood
folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
LLs <- c()
LLs2 <- c()
#LLs3 <- c()
mags <- c()
for(i in 1:length(ufiles)){
  # Extract the ODD object
  print(i)
  ODDy<-readRDS(paste0(folderin,ufiles[i]))
  AlgoParams$Np <- 1
  print(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams))
  AlgoParams$Np <- 2
  print(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams))
  #LLs <- append(LLs, DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams))
  #LLs2 <- append(LLs2, DispX(ODD = ODDy,Omega = Omega_MAP,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams))
  #LLs3 <- append(LLs3, DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams))
  #hrange<-grep("hazMean",names(ODDy),value = T)
  #mags <- append(mags, max(ODDy[hrange]@data, na.rm=TRUE))
}
plot(mags, LLs, xlab='Event Magnitude', ylab='Contribution to Log Likelihood')
points(mags, LLs2,col='blue')
#points(1:142, LLs3, col='red')

#Plot Linear Predictor Terms
for(i in 1:length(ufiles)){
  ODDSim <- readRDS(paste0(folderin,ufiles[i]))
  #simulate displacement, mortality and building destruction using DispX
  Params<-FormParams(ODDSim,list(Np=5,center=Model$center))
  # Income distribution percentiles & extract income percentile  
  SincN<-seq(10,90,by = 10); Sinc<-ExtractCIndy(ODDSim,var = paste0("p",SincN,"p100"))
  notnans<-which(!(is.na(ODDSim$Population) | is.na(ODDSim$ISO3C) | is.na(ODDSim$GDP)))
  LP <- GetLP(ODDSim,Omega,Params,Sinc,notnans)
  locallinp<-LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]]*LP$Plinp[ij]*LP$linp[[iso3c]]
}

#Misc plots to compare predictions from the true and MAP parameters
folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
x_actual = c()
x_pred_TRUE = c()
x_pred_MAP = c()
y_actual = c()
y_pred_TRUE = c()
y_pred_MAP = c()
z_actual = c()
z_pred_TRUE = c()
z_pred_MAP = c()

for(i in 1:length(ufiles)){
  ODDSim <- readRDS(paste0(folderin,ufiles[i]))
  #simulate displacement, mortality and building destruction using DispX
  ODDSim_Omega <- ODDSim %>% DispX(Omega, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
  ODDSim_MAP <- ODDSim %>% DispX(Omega_MAP, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
  x_actual = append(x_actual, ODDSim@gmax$gmax)
  x_pred_TRUE = append(x_pred_TRUE, ODDSim_Omega@predictDisp$disp_predictor)
  x_pred_MAP = append(x_pred_MAP, ODDSim_MAP@predictDisp$disp_predictor)
  y_actual = append(y_actual, ODDSim@gmax$mortality)
  y_pred_TRUE = append(y_pred_TRUE, ODDSim_Omega@predictDisp$mort_predictor)
  y_pred_MAP = append(y_pred_MAP, ODDSim_MAP@predictDisp$mort_predictor)
  z_actual = append(z_actual, ODDSim@gmax$buildDestroyed)
  z_pred_TRUE = append(z_pred_TRUE, ODDSim_Omega@predictDisp$nBD_predictor)
  z_pred_MAP = append(z_pred_MAP, ODDSim_MAP@predictDisp$nBD_predictor)
}

LL_Omega <- c()
LL_MAP <- c()
for(i in 1:length(ufiles)){
  ODDSim <- readRDS(paste0(folderin,ufiles[i]))
  #simulate displacement, mortality and building destruction using DispX
  LL_Omega <- append(LL_Omega, ODDSim %>% DispX(Omega, Model$center, Model$BD_params, LL=TRUE, Method=list(Np=1,cores=1,cap=-300))
  LL_MAP <- ODDSim %>% DispX(Omega_MAP, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
}

par(mfrow=c(1,1))
plot(x_actual,x_actual, xlim=c(1,100000),ylim=c(1,100000))
points(x_actual, x_pred_TRUE, col='blue')
points(x_actual, x_pred_MAP, col='red')

par(mfrow=c(1,1))
plot(y_actual, y_actual, xlim=c(1,1000), ylim=c(1,1000))
points(y_actual, y_pred_TRUE, col='blue')
points(y_actual, y_pred_MAP, col='red')


par(mfrow=c(1,1))
plot(log(z_actual), log(z_actual))
points(z_actual%>% log(), z_pred_TRUE%>% log(), col='blue')
points(z_actual %>% log(), z_pred_MAP%>% log(), col='red')

# 
# i = 10
# samples = c()
# folderin<-paste0(dir,"IIDIPUS_SimInput/ODDobjects/")
# ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
# for (k in 1:50){
#   ODDSim <- readRDS(paste0(folderin, ufiles[i]))
#   ODDSim_Omega <- ODDSim %>% DispX(Omega, Model$center, Model$BD_params, LL=FALSE, Method=list(Np=1,cores=1,cap=-300))
#   samples <- append(sampels, ODDSim_Omega@predictDisp$disp_predictor)
# }
# plot(samples, ylim=c(1800,2000))
# abline(h=ODDSim@gmax$gmax)
# 
# 
par(mfrow=c(4,4))
for(i in 1:14){
   ylim=c(min(unlist(Omega)[i], output[1:1500,i+1]), max(unlist(Omega)[i], output[1:1500,i+1]))
   plot(output[1:1500,i+1], type='l', ylab='', ylim=ylim)
   abline(h=unlist(Omega)[i], col='red')
}

BDSim <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_SimInput/BDobjects/EQ20130409ABC_-3')
BDSim@data
output <- readRDS('/home/manderso/Documents/GitHub/ODDRINfork/IIDIPUS_Results/output_2022-04-25_135511')

LLs1 <- c()
LLs2 <- c()
folderin<-paste0(dir,"IIDIPUS_Input/BDobjects/")
ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
for(i in 1:length(ufiles)){
  # Extract the BD object
  BDy<-readRDS(paste0(folderin,ufiles[i]))
  # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
  BDy@fIndies<-Model$fIndies
  # Apply BDX
  tLL1<-mean(tryCatch(BDX(BD = BDy,Omega = Omega,Model = Model,Method=AlgoParams, LL=T),
                error=function(e) NA))
  tLL2<-mean(tryCatch(BDX(BD = BDy,Omega = Omega_MAP,Model = Model,Method=AlgoParams, LL=T),
                  error=function(e) NA))
  LLs1 <- append(LLs1, tLL1)
  LLs2 <- append(LLs2, tLL2)
  # If all is good, add the LL to the total LL
  # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
}

plot(1:142, LLs1)
points(1:142, LLs2, col='red')

output <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/output_2022-04-30_170148')

for (i in 1:NROW(output)){
  if (is.na(output[i,1])){
    output[i,] <- output[i-1,]
  }
}
output <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/output_2022-05-06_093328')

cov2cor(cov(cbind(c(1,1,2), 
                  c(1,1,3))))

