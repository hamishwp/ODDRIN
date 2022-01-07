library(pracma)
library(plyr)
library(magrittr)
# library(rstan,lib.loc="../Rstan/")
library(rstan)
rstan_options(auto_write = TRUE)

mcs<-4


my_inits <- function(chainnum){
  list(
    M_i = as.matrix(M_i),
    Delta_i = as.matrix(Delta_i),
    tauis = as.matrix(tauis_i),
    beta  =  as.vector(beta_i),
    gompertz = gompertz_i   
  )
}

standata <- list(
  Surv_T_and_delta = Surv_T_and_delta,
  N = N,
  max_nj = max(nj),
  yhomeS = as.matrix(xhomeS),
  yclinicS = as.matrix(xS) ,
  yhomeD = as.matrix(xhomeD),
  yclinicD = as.matrix(xD) ,    
  sex=sex,
  race=race,
  age= age,
  
  m_M=m_M,
  tau_M=tau_M,
  m_Delta=m_Delta,
  tau_Delta=tau_Delta,
  
  alpha_CLINIC=alpha_CLINIC,
  beta_CLINIC=beta_CLINIC,
  alpha_HOME=alpha_HOME,
  beta_HOME=beta_HOME      
)

my_compiled_stan <- stan_model(file="./RCode/GP_regression.stan")

hwpstan<- sampling(my_compiled_stan, data=standata, init=my_inits , iter=6000,
                   warmup=4500, cores=mcs, chains=mcs, sample_file="./IIDIPUS_Results/GP_regression/")