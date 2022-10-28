##########################################################################
##########################################################################
###################### Machine Learning Methodologies ####################
### This file contains the Markov Chain Monte Carlo algorithm          ###
### that is used to optimise the model parameters, given historical    ###
### events and observation data, including building damage evaluations ###
###                                                                    ###
### I chose to use the Generalised Adaptative Metropolis-Hastings      ###
### algorithm with Global Adaptative Scaling                           ###
### see DOI 10.1007/s11222-008-9110-y                                  ###
### By C. Andrieu & J. Thoms, Stat Comput (2008) 18: 343–373           ###
### 'A tutorial on adaptive MCMC'                                      ###
### or otherwise http://drops.dagstuhl.de/opus/volltexte/2010/2813     ###
### C.L. Muller, (ETH) 'Exploring the common concepts of adaptive      ###
### MCMC and covariance matrix adaptation schemes'                     ###
##########################################################################
##########################################################################
##########################################################################

# Methodology parameters required
AlgoParams<-list(Np=20, # Number of Monte Carlo particles
                 cores=8, # Number of parallelised threads per event
                 NestedCores=1, # How many cores are to be used inside the ODD displacement calculations?
                 AllParallel=F, # Do you want fully-nested (data) parallelisation?
                 itermax=2000, # How many iterations do we want?
                 ABC=-500, # Approximate Bayesian Computation rejection
                 cap=-300, # if log values are too low, then log(mean(exp(LL)))=-Inf
                 GreedyStart=0, # How sure are we of the initial covariance matrix for accepted parameters? (Larger -> more confident)
                 Pstar=0.234, # Adaptive metropolis acceptance rate
                 gamzy0=0.2, # How quickly do the rejected parameters start having an influence on the covariance? (like GreedyStart) 
                 epsilon=50, # Do we still want values at larger numbers of iterations to have an influence on the covariance?
                 minVar=1e-4, # Prevent certain parameters from being too sure of themselves
                 t_0 =200,
                 eps = 0.000000001,
                 epsilon=c(0.15,0.03,0.15,0.1),
                 kernel='lognormal', #options are lognormal or loglaplace
                 smc_steps = 200, #Number of steps in the ABC-SMC algorithm
                 smc_Npart = 100, #Number of particles in the ABC-SMC algorithm
                 smc_alpha = 0.9
                 )
		 
if(is.null(AlgoParams$AllParallel)){
  if(AlgoParams$cores>4) { AlgoParams$AllParallel<-T
  } else AlgoParams$AllParallel<-F
}
# Choose the parameterisation algorithm - the string must match the function name exactly
Algorithm<-"AMCMC" # "NelderMeadOptim", "AMCMC"

# Metropolis-Hastings proposal distribution, given old values and covariance matrix
multvarNormProp <- function(xt, propPars){
  # purpose : A multivariate Gaussian random walk proposal for Met-Hastings
  #           MCMC
  # inputs  : xt       - The value of the chain at the previous time step 
  #           propPars - The correlation structure of the proposal
  return(array(mvtnorm::rmvnorm(1, mean=xt, sigma=propPars),dimnames = list(names(xt))))
}

Proposed2Physical<-function(proposed,Model,index=NULL){
  
  Model$links%<>%unlist()
  
  if(is.null(index)) index<-1:length(Model$links)
  # Link functions to convert values into useable/physical values
  for (i in index)  {
    proposed[i] <- match.fun(Model$links[[names(proposed)[i]]])(proposed[i], Model$par_lb[i],Model$par_ub[i])
  }
  # Reshape into desired structure
  proposed%>%relist(skeleton=Model$skeleton)
  
}

Physical2Proposed<-function(proposed,Model,index=NULL){
  
  proposed%<>%unlist()
  Model$unlinks%<>%unlist()
  
  if(is.null(index)) index<-1:length(Model$unlinks)
  # Link functions to convert values into useable/physical values
  for (i in index)  {
    proposed[i] <- match.fun(Model$unlinks[[names(proposed)[i]]])(proposed[i], Model$par_lb[i],Model$par_ub[i])
  }
  # Reshape into desired structure
  proposed%>%relist(skeleton=Model$skeleton)
  
}

Array2Physical <- function(array, Model){
  return(array %>% relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model))
}

# Generate the Adaptive Metropolis Global Scaling Factor discount factor
GenerateGamzy<-function(AlgoParams){
  AlgoParams$gamzy0/(1:(AlgoParams$itermax+AlgoParams$cores))^
    (seq(from=1/(1+AlgoParams$epsilon),to=1,length.out=(AlgoParams$itermax+AlgoParams$cores)))
}

# Check that the proposed initial values, model and methodology parameters
checkLTargs<-function(Model,iVals,AlgoParams){
  if(length(unlist(Model$links))!=length(unlist(iVals$x0))) stop("Mismatching link functions and initial values of Omega space")
  if(length(unlist(iVals$x0))!=nrow(iVals$COV) | length(unlist(iVals$x0))!=ncol(iVals$COV)) stop("Mismatching initial values and initial covariance matrix")
}

modifyAcc <- function(xNew, xPrev, Model){
  xNew %<>% relist(Model$skeleton) %>% unlist() %>% Proposed2Physical(Model) %>% unlist()
  xPrev %<>% relist(Model$skeleton) %>% unlist() %>% Proposed2Physical(Model) %>% unlist()
  Model$acceptTrans%<>%unlist()
  index<-1:length(unlist(Model$unlinks))
  product <- 1
  for (i in index){
    product <- product * match.fun(Model$acceptTrans[[names(xNew)[i]]])(xNew[i], xPrev[i], Model$par_lb[i],Model$par_ub[i])
  }
  return(product)
}

AMCMC <-function(dir, Model, iVals, AlgoParams, unfinished=F, tag=''){
  
  # Set Random Number Generator (RNG) initial seed
  set.seed(round(runif(1,0,100000))) 
  # Check no mistakes have been made in the model, methodology and initial values
  checkLTargs(Model,iVals,AlgoParams)
  AlgoParams$Np <- 10
  AlgoParams$itermax <- 100000
  
  xPrev<-unlist(Omega %>% Physical2Proposed(Model=Model)) #unlist(iVals$x0)
  n_x <- length(xPrev)
  s_d = (2.38)^2/n_x
  eps = diag(AlgoParams$eps, nrow=n_x)
  
  if (!unfinished){
    xbar_tminus1 <- xPrev
    output <- matrix(NA, nrow=AlgoParams$itermax, ncol=n_x+1)
    C_0 = diag(0.0001, nrow=n_x) / exp(xPrev) #iVals$COV 
    propCOV <- diag(n_x)
    it <- 1
    epsilon <- AlgoParams$epsilon_max
  } else {
    output <- readRDS(paste0(dir, '/IIDIPUS_Results/output_', tag)) 
    xbar_tminus1 <- readRDS(paste0(dir, '/IIDIPUS_Results/xbar_tminus1_', tag)) 
    propCOV <- readRDS(paste0(dir, '/IIDIPUS_Results/covariance_', tag)) 
    it <- min(which(is.na(output[,1])))-1
    xPrev[1:n_x] <- unlist(Physical2Proposed(relist(output[it, 2:NCOL(output)], skeleton=Model$skeleton), Model))
    if((it-1) > AlgoParams$t_0){
      epsilon <- AlgoParams$epsilon_min
    } else {
      epsilon <- AlgoParams$epsilon_max - (it-1) * (AlgoParams$epsilon_max-AlgoParams$epsilon_min)/AlgoParams$t_0
    }
  }
  C_0 <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2022-09-08_134604')$propCOV / 2
  xNew <- rep(NA,n_x)
  lTargNew<-alpha<-c()
  # Create file names with unique names for storage reasons
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  
  # Find first log-target value using initial values
  lTargOld_Np <- colSums(logTarget(dir = dir,Model = Model,
                             proposed = xPrev %>% Proposed2Physical(Model),AlgoParams = AlgoParams, epsilon=epsilon))
  
  lTargOld_Max <- max(lTargOld_Np, na.rm=T)
  lTargOld <- log(mean(exp(lTargOld_Np-lTargOld_Max), na.rm=T)) + lTargOld_Max
  
  output[it,] <- c(lTargOld, unlist(xPrev  %>% Proposed2Physical(Model)))
  
  # Start the iterations!
  it <- it + 1
  while (it <= AlgoParams$itermax){
    print(it)
    t <- it - 1
  
    # Parameter proposal
    xNew <- xPrev
    if (t > AlgoParams$t_0){
      xNew <- multvarNormProp(xt=xPrev, propPars=propCOV)
      if(any(xNew < -20) || any(xNew > 20)){
        next
      }
      epsilon <- AlgoParams$epsilon_min
    } else {
      xNew <- multvarNormProp(xt=xPrev, propPars=C_0)
      epsilon <- AlgoParams$epsilon_max - t * (AlgoParams$epsilon_max-AlgoParams$epsilon_min)/AlgoParams$t_0
    }
    
    # Check proposal is within the parameter space:
    #if(any(xNew < Model$par_lb) | any(xNew > Model$par_ub) ){
    #  output[it,] <- c(lTargOld, xPrev)
    #  propCOV <- (t-1)/t * propCOV + s_d/(t+1) * (xPrev - xbar_tminus1) %*% t(xPrev - xbar_tminus1) + s_d /t * eps
    #  xbar_tminus1 <- (t * xbar_tminus1 + xPrev)/(t+1)
    #  it <- it + 1 
    #  next
    #}
    
    # Convert parameters to physical/useable values
    if(!is.null(Model$links)) xProp<-xNew%>%Proposed2Physical(Model)
    
    start.time <- Sys.time()
    
    # Calculate log-target value
    lTargNew_Np <- tryCatch(colSums(logTarget(dir = dir,Model = Model,proposed = xProp,
                                   AlgoParams = AlgoParams, epsilon= epsilon)), error=function(e) NA)
    
    lTargNew_Max <- max(lTargNew_Np, na.rm=T)
    lTargNew <- log(mean(exp(lTargNew_Np-lTargNew_Max), na.rm=T)) + lTargNew_Max
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    
    # Check if we have a NaN
    if(is.na(lTargNew)|is.infinite(lTargNew)) {
      output[it,] <- c(lTargOld, xPrev %>% Proposed2Physical(Model) %>% unlist())
      propCOV <- (t-1)/t * propCOV + s_d/(t+1) * (xPrev - xbar_tminus1) %*% t(xPrev - xbar_tminus1) + s_d /t * eps
      xbar_tminus1 <- (t * xbar_tminus1 + xPrev)/(t+1)
      it <- it + 1 
      next
    }
    
    #if proposed log likelihood is close to the old log likelihood, then resample the old log likelihood
    if (lTargNew > (lTargOld - 50) || (t <= AlgoParams$t_0 & t %% 10 == 0 )){
      lTargOld_Np <- tryCatch(colSums(logTarget(dir = dir,Model = Model,proposed = xPrev %>%Proposed2Physical(Model),
                                     AlgoParams = AlgoParams, epsilon=epsilon)), error=function(e) NA)
      lTargOld_Max <- max(lTargOld_Np, na.rm=T)
      lTargOld <- log(mean(exp(lTargOld_Np-lTargOld_Max), na.rm=T)) + lTargOld_Max
    }
    
    # Prepare for acceptance
    u <- runif(1)
    
    # Acceptance probability
    alpha <- min(exp(lTargNew - lTargOld) * modifyAcc(xNew, xPrev, Model), 1)
    
    # Metropolis Acceptance Algorithm
    if (alpha>=u) { # Accepted!
      print('ACCEPTED!')
      xPrev<-xNew
      lTargOld <- lTargNew
    } 
    
    propCOV <- (t-1)/t * propCOV + s_d/(t+1) * (xPrev - xbar_tminus1) %*% t(xPrev - xbar_tminus1) + s_d /t * eps
    xbar_tminus1 <- (t * xbar_tminus1 + xPrev)/(t+1)
    
    print(paste0(round(it*100/AlgoParams$itermax),"% done. LL = ",lTargOld))
    print(" ")
    
    output[it,] <- c(lTargOld, xPrev %>% Proposed2Physical(Model) %>% unlist())
    
    par(mfrow=c(4,4))
    # Plot S-curves for the actual and current
    plot_S_curves(Omega,  xPrev %>% Proposed2Physical(Model))
    
    #create trace plots
    for(i in 1:n_x){
      ylim=c(min(unlist(Omega)[i], output[1:it,i+1]), max(unlist(Omega)[i], output[1:it,i+1]))
      plot(output[1:it,i+1], type='l', ylab='', ylim=ylim)
      abline(h=unlist(Omega)[i], col='red')
    }
    
    plot(output[ifelse(it>500,500,1):it,1], type='l', ylab='')
    
    print(paste('Single Chain R-Hat', paste(round(apply(output[round(it/2):it,2:(n_x+1)],2, rhat, split = TRUE), digits=2), collapse=' ')))
    print(paste('ESS', paste(round(apply(output[min(1000, round(it/2)):it,2:(n_x+1)],2, ess_basic), digits=2), collapse=' ')))
    
    
    # Save log-target and parameters
    saveRDS(output,paste0(dir,"IIDIPUS_Results/output_",tag))
    # Save covariance matrix
    saveRDS(propCOV,paste0(dir,"IIDIPUS_Results/covariance_",tag))
    
    saveRDS(xbar_tminus1,paste0(dir,"IIDIPUS_Results/xbar_tminus1_",tag))
    print(cov2cor(propCOV))
    
    it <- it + 1
  }
  
  return(list(PhysicalValues=output[which.max(output[,1]),2:ncol(output)] %>% 
                relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model), # MAP value 
              OptimisationOut=output))
}

Algorithm <- match.fun('AMCMC')

abcSmc <- function(AlgoParams, Model, unfinished=F, oldtag=''){
  #Input: 
  # - AlgoParams: Parameters describing the ABC-SMC Algorithm (e.g. the ABC rejection threshold for higher level priors)
  # - Model: Describes the data simulation and calculation of the distance measure 
  # - Unfinished: If TRUE, then include oldtag - the tag (end of the filename) of an unfinished ABC-SMC run to be completed.
  # Output:
  # - A list containing the parameter values for each particle at each step of the ABC-SMC algorithm, as well as the
  #   corresponding weights, distances and final perturbation kernel. 
  # Details:
  # - Uses the algorithm described on left of Page 5 of Beaument et al., 2007: https://arxiv.org/pdf/0805.2256.pdf
  # - Uses the multivariate perturbation kernel described in Filippi et al., 2012: https://arxiv.org/pdf/1106.6280.pdf
  

  steps = AlgoParams$smc_steps
  Npart = AlgoParams$smc_Npart
  
  start_time <- Sys.time()
  n_x <- length(Model$par_lb) #n_x = number of parameters
  Omega_sample <- array(NA, dim=c(Npart, n_x, steps)) #store sampled parameters on the transformed space
  Omega_sample_phys <- array(NA, dim=c(Npart, n_x, steps)) #store sampled parameters on the untransformed space
  W <- array(NA, dim=c(Npart, steps))
  d <- array(Inf, dim=c(Npart, steps))
  
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  
  if(unfinished==F){ #Initialize and perform sampling for s=1
    W[,1] <- 1/Npart # Give each particle an equal weight
    for (n in 1:Npart){
      print(paste('Step: 1, Particle:', n))
      while(d[n,1] > 200000){ #make less arbitrary! 50% quantile?
        #draw initial particle values from the prior and ensure that they satisfy the higher level prior
        Omega_sample[n,,1] <- HLPrior_sample(Model, AlgoParams)
        Omega_sample_phys[n,,1] <-  Omega_sample[n,,1] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
        #calculate distance
        d[n,1] <- -logTarget(dir = dir,Model = Model,
                             proposed = Omega_sample_phys[n,,1] %>% relist(skeleton=Model$skeleton), 
                             AlgoParams = AlgoParams)
      } 
    }
    saveRDS(
      list(d=d, 
           omegastore=Omega_sample_phys,
           propCOV=array(0, dim=c(n_x,n_x)),
           W=W),
      paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
    )  
    s_start = 2 # continue the algorithm from s = 2
    n_start = 1
  } else { #Collect relevant information from the unfinished sample
    output_unfinished <- readRDS(paste0(dir,"AWS_IIDIPUS_Results/abcsmc_",oldtag))
    s_finish = max(which(colSums(!is.na(output_unfinished$omegastore[,1,]))>0)) #find last completed step
    n_finish = max(which(!is.na(output_unfinished$omegastore[,1,s_finish])))
    W[,1:s_finish] <- output_unfinished$W[,1:s_finish]
    d[,1:s_finish] <- output_unfinished$d[,1:s_finish]
    propCOV <- output_unfinished$propCOV
    Omega_sample_phys[,,1:s_finish] <- output_unfinished$omegastore[,,1:s_finish]
    for (n in 1:Npart){
      for (s in 1:s_finish){
        #convert particles from the physical to the transformed space
        if(!is.na(Omega_sample_phys[n,1,s]))
        Omega_sample[n,,s] <- Omega_sample_phys[n,,s] %>% relist(skeleton=Model$skeleton) %>% Physical2Proposed(Model) %>% unlist()
      }
    }
    s_start = ifelse(n_finish==Npart, s_finish+1, s_finish) #identify the appropriate step from which to continue the algorithm
    n_start = ifelse(n_finish==Npart, 1, n_finish + 1) #identify the appropriate particle from which to continue the algorithm
  }

  for (s in s_start:steps){
    
    #record and print time for each step
    end_time <- Sys.time()
    print(paste('Time:', end_time-start_time))
    start_time <- Sys.time()
    
    #update the tolerance using the 80th quantile
    tolerance <- quantile(d[,s-1], probs=0.6)
    print(paste('      Step:', s, ', New tolerance is:', tolerance))
    
    #calculate perturbation covariance using Filippi et al., 2012
    tilda_i <- which(d[,s-1] < tolerance) #identify old particles that fall within the new tolerance
    Omega_tilda <- Omega_sample[tilda_i,,s-1] 
    W_tilda <- W[tilda_i,s-1]
    W_tilda <- W_tilda/sum(W_tilda) #normalise weights
    
    propCOV <- matrix(0, n_x, n_x)
    for(i in 1:Npart){
      for(k in 1:length(tilda_i)){
        propCOV <- propCOV + W[i,s-1] * W_tilda[k] * ((Omega_tilda[k,]-Omega_sample[i,,s-1]) %*% t(Omega_tilda[k,]-Omega_sample[i,,s-1]))
      }
    }
    
    for(n in n_start:Npart){
      print(paste(' Step:', s, ', Particle:', n))
      d[n,s] <- tolerance + 1
      while (d[n,s] > tolerance){ 
        weighted_sample <- sample(1:Npart, 1, prob=W[,s-1]) #draw weighted sample from the previous population
        Omega_sample[n,,s] <- multvarNormProp(xt=Omega_sample[weighted_sample,,s-1], propPars=propCOV) #perturb the proposal
        
        Omega_sample_phys[n,,s] <-  Omega_sample[n,,s] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% 
                                          Proposed2Physical(Model) %>% unlist() #convert to physical space
        
        d[n,s] <-  -logTarget(dir = dir,Model = Model,
                              proposed = Omega_sample_phys[n,,s] %>% relist(skeleton=Model$skeleton), 
                              AlgoParams = AlgoParams) #calculate distance
      }
      W[n,s] <- 1/sum(W[,s-1]*apply(-sweep(Omega_sample[,,s-1], 2, Omega_sample[n,,s]), 1, dmvnorm, mean=rep(0,n_x), sigma=propCOV)) #update weight
      saveRDS(
        list(d=d, 
             omegastore=Omega_sample_phys,
             propCOV=propCOV,
             W=W),
        paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
      )
    }
    W[,s] <- W[,s]/sum(W[,s]) #normalise weights
    print(paste('Normalised Weights', W[,s]))
    n_start = 1
    
    par(mfrow=c(4,4))
    for (i in 1:n_x){
      plot(rep(1:s,each=Npart),Omega_sample_phys[,i,1:s], ylim=c(min(Omega_sample_phys[,i,1:s], unlist(Omega)[i]), max(Omega_sample_phys[,i,1:s], unlist(Omega)[i])),
           main=names(unlist(Omega))[i], xlab='step', ylab='')
      abline(h=unlist(Omega)[i], col='red')
    }
    par(mfrow=c(1,1))
    
    saveRDS(
      list(d=d, 
           omegastore=Omega_sample_phys,
           propCOV=propCOV,
           W=W),
      paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
    )  
  }
}

abcSmc_delmoral <- function(AlgoParams, Model, unfinished=F, oldtag=''){
  #Input: 
  # - AlgoParams: Parameters describing the ABC-SMC Algorithm (e.g. the ABC rejection threshold for higher level priors)
  # - Model: Describes the data simulation and calculation of the distance measure 
  # - Unfinished: If TRUE, then include oldtag - the tag (end of the filename) of an unfinished ABC-SMC run to be completed.
  # Output:
  # - A list containing the parameter values for each particle at each step of the ABC-SMC algorithm, as well as the
  #   corresponding weights, distances and final perturbation kernel. 
  # Details:
  # - Uses the algorithm described in Del Moral et al., 2011: https://link.springer.com/content/pdf/10.1007/s11222-011-9271-y.pdf
  # - Uses the multivariate perturbation kernel described in Filippi et al., 2012: https://arxiv.org/pdf/1106.6280.pdf
  
  steps = 400 #AlgoParams$smc_steps
  Npart = AlgoParams$smc_Npart
  M = AlgoParams$Np
  
  start_time <- Sys.time()
  n_x <- length(Model$par_lb) #n_x = number of parameters
  Omega_sample <- array(NA, dim=c(Npart, n_x, steps)) #store sampled parameters on the transformed space
  Omega_sample_phys <- array(NA, dim=c(Npart, n_x, steps)) #store sampled parameters on the untransformed space
  W <- array(NA, dim=c(Npart, steps)) #Weights
  d <- array(Inf, dim=c(Npart, M, steps)) #Distances
  d_full <- array(0, dim=c(Npart, 178, M, steps))
  npa=rep(0,Npart) #Number of alive particles (have at least one distance less than the tolerance)
  tolerance = 10000000 #Ensure all initial distances are less than the tolerance
  tolerancetarget=1 #LOOSEEND: need to add in an adaptive stopping rule
  tolerancestore=c(tolerance)
  essstore=c(Npart)
  N_T=Npart/2
  alpha= 0.9 #AlgoParams$smc_alpha
  
  tpa<-function(tolerance,x){ #calculate the total proportion of alive particles (have at least one distance less than the tolerance)
    d<-abs(x)
    for (i in 1:Npart){
      count<-0
      for (j in 1:M){
        if (d[i,j]<tolerance){
          count<-count+1
        }
      }
      npa[i]<-count
    }
    return(sum(npa>0)/Npart)
  }
  
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  
  if(unfinished==F){ #Initialize and perform sampling for s=1
    W[,1] <- 1/Npart # Give each particle an equal weight
    for (n in 1:Npart){
      print(paste('Step: 1, Particle:', n))
      #draw initial particle values from the prior and ensure that they satisfy the higher level prior
      Omega_sample[n,,1] <- HLPrior_sample(Model, AlgoParams)
      Omega_sample_phys[n,,1] <-  Omega_sample[n,,1] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
      #calculate distance
      #CHECK HLP
      dist_sample <- sampleDist(dir = dir,Model = Model,
                           proposed = Omega_sample_phys[n,,1] %>% relist(skeleton=Model$skeleton) %>% add_Loc_Params(), 
                           AlgoParams = AlgoParams)
      d[n,,1] <- logTarget2(dist_sample)
      #d_full[n,,,1] <- 
    }
    saveRDS(
      list(d=d, 
           d_full=d_full,
           omegastore=Omega_sample_phys,
           propCOV=array(0, dim=c(n_x,n_x)),
           W=W, 
           tolerancestore = tolerancestore,
           essstore = essstore),
      paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
    )  
    s_start = 2 # continue the algorithm from s = 2
    n_start = 1
  } else { #Collect relevant information from the unfinished sample
    output_unfinished <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2022-09-15_121638') #readRDS(paste0(dir,"IIDIPUS_Results/abcsmc_",oldtag))
    s_finish = max(which(colSums(!is.na(output_unfinished$omegastore[,1,]))>0)) #find last completed step
    n_finish = max(which(!is.na(output_unfinished$omegastore[,1,s_finish])))
    W[,1:s_finish] <- output_unfinished$W[,1:s_finish]
    d[,,1:s_finish] <- output_unfinished$d[,,1:s_finish]
    #d_full[,,,1:s_finish] <- output_unfinished$d_full[,,,1:s_finish]
    propCOV <- output_unfinished$propCOV
    Omega_sample_phys[,,1:s_finish] <- output_unfinished$omegastore[,,1:s_finish]
    for (n in 1:Npart){
      for (s in 1:s_finish){
        #convert particles from the physical to the transformed space
        if(!is.na(Omega_sample_phys[n,1,s]))
          Omega_sample[n,,s] <- Omega_sample_phys[n,,s] %>% relist(skeleton=Model$skeleton) %>% Physical2Proposed(Model) %>% unlist()
      }
    }
    s_start = ifelse(n_finish==Npart, s_finish+1, s_finish) #identify the appropriate step from which to continue the algorithm
    n_start = ifelse(n_finish==Npart, 1, n_finish + 1) #identify the appropriate particle from which to continue the algorithm
    tolerance = ifelse(n_finish==Npart, output_unfinished$tolerancestore[s_finish], output_unfinished$tolerancestore[s_finish-1])
  }
  
  for (s in s_start:steps){
    
    #record and print time for each step
    end_time <- Sys.time()
    print(paste('Time:', end_time-start_time))
    start_time <- Sys.time()
    
    #Find the new tolerance such that alpha proportion of the current alive particles stay alive
    toleranceold <- tolerance
    reflevel <- alpha * tpa(toleranceold, d[,,s-1])
    tolerance<-uniroot(function(tolerance) tpa(tolerance,x=d[,,s-1])-reflevel,c(0,toleranceold))$root
    print(paste('      Step:', s, ', New tolerance is:', tolerance))
    if (tolerance<tolerancetarget){
         tolerance<-tolerancetarget
    }
    
    essstore<-c(essstore, 1/sum(W[,s-1]^2)) #effective sample size
    tolerancestore<-c(tolerancestore, tolerance)
    
    #compute the associated weights
    npa_old<-rowSums(d[,,s-1]<toleranceold)
    npa<-rowSums(d[,,s-1]<tolerance)
    a<-which(npa_old>0)
    b<-which(npa_old==0)
    W[a,s]<-W[a,s-1]*npa[a]/npa_old[a]
    W[b,s]<-rep(0,length(b))
    W[,s]<-W[,s]/sum(W[,s])

    # resample if necessary
    if ((sum(W[,s]^2)*N_T)>1){
         choice<-sample(1:Npart,Npart,replace= TRUE, prob = W[,s])
         Omega_sample[,,s]<-Omega_sample[choice,,s-1] 
         Omega_sample_phys[,,s]<-Omega_sample_phys[choice,,s-1] 
         d[,,s]<-d[choice,,s-1] 
         #d_full[,,,s]<-d_full[choice,,,s-1] 
         W[,s]<-rep(1/Npart,Npart) 
    } else { #these will be perturbed later via a MCMC step. 
        Omega_sample[,,s]<-Omega_sample[,,s-1] 
        Omega_sample_phys[,,s]<-Omega_sample_phys[,,s-1] 
        d[,,s]<-d[,,s-1] 
        d_full[,,,s]<-d_full[,,,s-1] 
    }
    
    #calculate perturbation covariance based on Filippi et al., 2012
    tilda_i <- which(rowSums(d[,,s-1]<tolerance)>4) #identify old particles that fall within the new tolerance
    Omega_tilda <- Omega_sample[tilda_i,,s-1] 
    W_tilda <- W[tilda_i,s-1]
    W_tilda <- W_tilda/sum(W_tilda) #normalise weights
    
    propCOV <- matrix(0, n_x, n_x)
    for(i in 1:Npart){
      for(k in 1:length(tilda_i)){
        propCOV <- propCOV + W[i,s-1] * W_tilda[k] * ((Omega_tilda[k,]-Omega_sample[i,,s-1]) %*% t(Omega_tilda[k,]-Omega_sample[i,,s-1]))
      }
    }
    
    for(n in n_start:Npart){
      print(paste(' Step:', s, ', Particle:', n))
      if(W[n,s]>0){
        Omega_prop <- multvarNormProp(xt=Omega_sample[n,,s], propPars=propCOV) #perturb the proposal
        Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
        
        HP<- Model$HighLevelPriors(Omega_prop_phys %>% add_Loc_Params(), Model)
        if (HP> AlgoParams$ABC) next
        
        dist_sample <- sampleDist(dir = dir,Model = Model,
                                  proposed = Omega_prop_phys %>% add_Loc_Params(), 
                                  AlgoParams = AlgoParams)
        d_prop <- logTarget2(dist_sample)
        
        if(d_prop[1]==Inf){#if (d_full_prop[1]==Inf){
          d_prop <- Inf
        } 
        
        #calculate the acceptance probability:
        acc <- sum(d_prop<tolerance)/sum(d[n,,s]<tolerance) * modifyAcc(Omega_prop, Omega_sample[n,,s], Model)
        u <- runif(1)
        if(u < acc){
          Omega_sample[n,,s] <- Omega_prop
          Omega_sample_phys[n,,s] <- Omega_sample[n,,s] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
          #d_full[n,,,s] <- d_full_prop
          d[n,,s] <- d_prop
        }
        saveRDS(
          list(d=d,
               d_full=d_full,
               omegastore=Omega_sample_phys,
               propCOV=propCOV,
               W=W,
               tolerancestore=tolerancestore,
               essstore=essstore),
          paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
        )
      }
    }
    
    n_start = 1
   
    #plot the output so far:
    i_zeroweight <- which(W==0,  arr.ind=TRUE)
    plot(rep(1, Npart), apply(d[,,1], c(1), min), xlim=c(1, s), ylim=c(min(d), max(d[which(d<Inf)])))
    for (i in 2:s){
      i_s_zeroweight <- i_zeroweight[which(i_zeroweight[,2]==i),1]
      if (length(i_s_zeroweight)==0){
        points(rep(i, Npart), apply(d[,,i], c(1), min))
      }
      else{
        points(rep(i, Npart-length(i_s_zeroweight)), apply(d[-i_s_zeroweight,,i], c(1), min))
      }
    }
    
    par(mfrow=c(4,4))
    for (i in 1:n_x){
      plot(rep(1:s,each=Npart),Omega_sample_phys[,i,1:s], ylim=c(min(Omega_sample_phys[,i,1:s], unlist(Omega)[i]), max(Omega_sample_phys[,i,1:s], unlist(Omega)[i])),
           main=names(unlist(Omega))[i], xlab='step', ylab='')
      abline(h=unlist(Omega)[i], col='red')
    }
    par(mfrow=c(1,1))
    
  }
}

initialise_particles <- function(dir, Npart, n_nodes){
  particle_divisions <- split(1:Npart, sort(1:Npart%%n_nodes))
  allocated_particles <- particle_divisions[[mpi.comm.rank()]]
  n_allocated <- length(allocated_particles)
  n_x <- length(Model$par_lb)
  Omega_sample_node <- array(NA, dim=c(n_allocated, n_x))
  Omega_sample_phys_node <- array(NA, dim=c(n_allocated, n_x))
  d_node <- array(Inf, dim=c(n_allocated, AlgoParams$Np))
  for (n in 1:n_allocated){
    HLP <- 100
    while (HLP > 1){
      Omega_sample_i <- runif(n_x, min=Model$par_lb, max=Model$par_ub) %>% relist(skeleton=Model$skeleton) %>% Physical2Proposed(Model) %>% unlist()
      Omega_sample_phys_i <-  Omega_sample_i %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
      HLP = Model$HighLevelPriors(Omega_sample_phys_i %>% relist(skeleton=Model$skeleton),Model)
    }
    #calculate distance
    d_i <- -logTarget(dir = dir, Model = Model,
                      proposed = Omega_sample_phys_i %>% relist(skeleton=Model$skeleton), 
                      AlgoParams = AlgoParams)
    
    Omega_sample_node[n,] = Omega_sample_i 
    Omega_sample_phys_node[n,] = Omega_sample_phys_i
    d_node[n,] = d_i 
  }
  return(list(Omega_sample_node=Omega_sample_node,
              Omega_sample_phys_node=Omega_sample_phys_node,
              d_node=d_node))
}

perturb_particles <- function(dir, Npart, n_nodes, W_curr, Omega_curr, Omega_phys_curr, d_curr, propCOV, tolerance){
  particle_divisions <- split(1:Npart, sort(1:Npart%%n_nodes))
  allocated_particles <- particle_divisions[[mpi.comm.rank()]]
  n_allocated <- length(allocated_particles)
  n_x <- length(Model$par_lb)
  
  Omega_sample_node <- Omega_curr[allocated_particles,]
  Omega_sample_phys_node <- Omega_phys_curr[allocated_particles,]
  d_node <- d_curr[allocated_particles,]
  
  for(n in 1:n_allocated){
    if(W_curr[n]>0){
      Omega_prop <- multvarNormProp(xt=Omega_curr[n,], propPars=propCOV) #perturb the proposal
      Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
      
      d_prop <-  -logTarget(dir = dir, Model = Model,
                            proposed = Omega_prop_phys, 
                            AlgoParams = AlgoParams) #calculate distance
      acc <- sum(d_prop<tolerance)/sum(d_curr[n,]<tolerance) * modifyAcc(Omega_prop, Omega_curr[n,], Model)
      u <- runif(1)
      if(u < acc){
        Omega_sample_node[n,] <- Omega_prop
        Omega_sample_phys_node[n,] <- Omega_sample_node[n,] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
        d_node[n,] <- d_prop
      }
    }
  }
  return(list(Omega_sample_node=Omega_sample_node, 
              Omega_sample_phys_node=Omega_sample_phys_node,
              d_node=d_node))
}


delmoral_parallel <- function(AlgoParams, Model, unfinished=F, oldtag=''){
  #Input: 
  # - AlgoParams: Parameters describing the ABC-SMC Algorithm (e.g. the ABC rejection threshold for higher level priors)
  # - Model: Describes the data simulation and calculation of the distance measure 
  # - Unfinished: If TRUE, then include oldtag - the tag (end of the filename) of an unfinished ABC-SMC run to be completed.
  # Output:
  # - A list containing the parameter values for each particle at each step of the ABC-SMC algorithm, as well as the
  #   corresponding weights, distances and final perturbation kernel. 
  # Details:
  # - Uses the algorithm described in Del Moral et al., 2011: https://link.springer.com/content/pdf/10.1007/s11222-011-9271-y.pdf
  # - Uses the multivariate perturbation kernel described in Filippi et al., 2012: https://arxiv.org/pdf/1106.6280.pdf
  
  steps = AlgoParams$smc_steps
  Npart = AlgoParams$smc_Npart
  M = AlgoParams$Np
  n_nodes = AlgoParams$n_nodes
  
  start_time <- Sys.time()
  n_x <- length(Model$par_lb) #n_x = number of parameters
  Omega_sample <- array(NA, dim=c(Npart, n_x, steps)) #store sampled parameters on the transformed space
  Omega_sample_phys <- array(NA, dim=c(Npart, n_x, steps)) #store sampled parameters on the untransformed space
  W <- array(NA, dim=c(Npart, steps)) #Weights
  d <- array(Inf, dim=c(Npart, M, steps)) #Distances
  npa=rep(0,Npart) #Number of alive particles (have at least one distance less than the tolerance)
  tolerance = 10000000 #Ensure all initial distances are less than the tolerance
  tolerancetarget=1 #LOOSEEND: need to add in an adaptive stopping rule
  tolerancestore=c(tolerance)
  essstore=c(Npart)
  N_T=Npart/2
  alpha=AlgoParams$smc_alpha
  
  
  tpa<-function(tolerance,x){ #calculate the total proportion of alive particles (have at least one distance less than the tolerance)
    d<-abs(x)
    for (i in 1:Npart){
      count<-0
      for (j in 1:M){
        if (d[i,j]<tolerance){
          count<-count+1
        }
      }
      npa[i]<-count
    }
    return(sum(npa>0)/Npart)
  }
  
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  
  mpi.spawn.Rslaves(nslaves=2)
  mpi.bcast.Rfun2slave()
  mpi.bcast.Robj2slave(dir)
  mpi.bcast.Robj2slave(directory)
  mpi.bcast.Robj2slave(Model)
  mpi.bcast.Robj2slave(AlgoParams)
  mpi.remote.exec(GetODDPackages, FALSE)
  mpi.remote.exec(source('RCode/Simulate.R'))
  
  
  if(unfinished==F){ #Initialize and perform sampling for s=1
    W[,1] <- 1/Npart # Give each particle an equal weight
    
    node_return <- mpi.remote.exec(initialise_particles, dir, Npart, n_nodes)
    particle_divisions <- split(1:Npart, sort(1:Npart%%n_nodes))
    for (j in 1:length(node_return)){
      Omega_sample[particle_divisions[[j]],,1] <- node_return[[j]]$Omega_sample_node
      Omega_sample_phys[particle_divisions[[j]],,1] <- node_return[[j]]$Omega_sample_phys_node
      d[particle_divisions[[j]],,1] <- node_return[[j]]$d_node
    }
    
    s_start = 2 # continue the algorithm from s = 2
    saveRDS(
      list(d=d, 
           omegastore=Omega_sample_phys,
           propCOV=array(0, dim=c(n_x,n_x)),
           W=W, 
           tolerance = tolerance),
      paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
    )  
  } else { #Collect relevant information from the unfinished sample
    output_unfinished <- readRDS(paste0(dir,"IIDIPUS_Results/abcsmc_",oldtag))
    s_finish = max(which(colSums(!is.na(output_unfinished$omegastore[,1,]))>0)) #find last completed step
    n_finish = max(which(!is.na(output_unfinished$omegastore[,1,s_finish])))
    W[,1:s_finish] <- output_unfinished$W[,1:s_finish]
    d[,,1:s_finish] <- output_unfinished$d[,,1:s_finish]
    propCOV <- output_unfinished$propCOV
    Omega_sample_phys[,,1:s_finish] <- output_unfinished$omegastore[,,1:s_finish]
    for (n in 1:Npart){
      for (s in 1:s_finish){
        #convert particles from the physical to the transformed space
        if(!is.na(Omega_sample_phys[n,1,s]))
          Omega_sample[n,,s] <- Omega_sample_phys[n,,s] %>% relist(skeleton=Model$skeleton) %>% Physical2Proposed(Model) %>% unlist()
      }
    }
    s_start = ifelse(n_finish==Npart, s_finish+1, s_finish) #identify the appropriate step from which to continue the algorithm
    n_start = ifelse(n_finish==Npart, 1, n_finish + 1) #identify the appropriate particle from which to continue the algorithm
    tolerance = output_unfinished$tolerance
  }
  
  for (s in s_start:steps){
    
    #record and print time for each step
    end_time <- Sys.time()
    print(paste('Time:', end_time-start_time))
    start_time <- Sys.time()
    
    #Find the new tolerance such that alpha proportion of the current alive particles stay alive
    toleranceold <- tolerance
    reflevel <- alpha * tpa(toleranceold, d[,,s-1])
    tolerance<-uniroot(function(tolerance) tpa(tolerance,x=d[,,s-1])-reflevel,c(0,toleranceold))$root
    print(paste('      Step:', s, ', New tolerance is:', tolerance))
    if (tolerance<tolerancetarget){
      tolerance<-tolerancetarget
    }
    tolerancestore<-c(tolerancestore, tolerance)
    essstore<-c(essstore, 1/sum(W[,s-1]^2)) #effective sample size
    
    #compute the associated weights
    npa_old<-rowSums(d[,,s-1]<toleranceold)
    npa<-rowSums(d[,,s-1]<tolerance)
    a<-which(npa_old>0)
    b<-which(npa_old==0)
    W[a,s]<-W[a,s-1]*npa[a]/npa_old[a]
    W[b,s]<-rep(0,length(b))
    W[,s]<-W[,s]/sum(W[,s])
    
    # resample if necessary
    if ((sum(W[,s]^2)*N_T)>1){
      choice<-sample(1:Npart,Npart,replace= TRUE, prob = W[,s])
      Omega_sample[,,s]<-Omega_sample[choice,,s-1] 
      Omega_sample_phys[,,s]<-Omega_sample_phys[choice,,s-1] 
      d[,,s]<-d[choice,,s-1] 
      W[,s]<-rep(1/Npart,Npart) 
    } else { #these will be perturbed later via a MCMC step. 
      Omega_sample[,,s]<-Omega_sample[,,s-1] 
      Omega_sample_phys[,,s]<-Omega_sample_phys[,,s-1] 
      d[,,s]<-d[,,s-1] 
    }
    
    #calculate perturbation covariance based on Filippi et al., 2012
    tilda_i <- which(rowSums(d[,,s-1]<tolerance)>0) #identify old particles that fall within the new tolerance
    Omega_tilda <- Omega_sample[tilda_i,,s-1] 
    W_tilda <- W[tilda_i,s-1]
    W_tilda <- W_tilda/sum(W_tilda) #normalise weights
    
    propCOV <- matrix(0, n_x, n_x)
    for(i in 1:Npart){
      for(k in 1:length(tilda_i)){
        propCOV <- propCOV + W[i,s-1] * W_tilda[k] * ((Omega_tilda[k,]-Omega_sample[i,,s-1]) %*% t(Omega_tilda[k,]-Omega_sample[i,,s-1]))
      }
    }
    
    node_return <- mpi.remote.exec(perturb_particles, dir, Npart, n_nodes, W[,s], Omega_sample[,,s], Omega_sample_phys[,,s], d[,,s], propCOV, tolerance)
    
    particle_divisions <- split(1:Npart, sort(1:Npart%%n_nodes))
    for (j in 1:length(node_return)){
      Omega_sample[particle_divisions[[j]],,s] <- node_return[[j]]$Omega_sample_node
      Omega_sample_phys[particle_divisions[[j]],,s] <- node_return[[j]]$Omega_sample_phys_node
      d[particle_divisions[[j]],,s] <- node_return[[j]]$d_node
    }
    
    #plot the output so far:
    par(mfrow=c(4,4))
    for (i in 1:n_x){
      plot(rep(1:s,each=Npart),Omega_sample_phys[,i,1:s], ylim=c(min(Omega_sample_phys[,i,1:s], unlist(Omega)[i]), max(Omega_sample_phys[,i,1:s], unlist(Omega)[i])),
           main=names(unlist(Omega))[i], xlab='step', ylab='')
      abline(h=unlist(Omega)[i], col='red')
    }
    par(mfrow=c(1,1))
    
    saveRDS(
      list(d=d, 
           omegastore=Omega_sample_phys,
           propCOV=propCOV,
           W=W, 
           tolerance=tolerance),
      paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
    )  
  }
}


NelderMeadOptim<-function(dir,Model,iVals,AlgoParams){
  
  # We don't need an initial guess of the covariance matrix for MLE optimisation
  x0=Physical2Proposed(iVals$x0,Model)%>%unlist()
  # Only optimise over the input iVals that are not NAs
  if(is.null(AlgoParams$indices)) AlgoParams$indices<-length(x0)
  # Cost function (note that this still includes the priors and ABC rejection, so isn't purely frequentist MLE)
  Fopty<-function(vals){
    #x0[!AlgoParams$indices]<-vals
    x0 <- vals
    # Convert proposal in to physical values ODDRIN understands
    x0%<>%Proposed2Physical(Model)
    # Trust me, it's nice to know the parameters tested, this optim algorithm can struggle with stochastic target distributions
    print(unname(unlist(x0)))
    # Posterior calculation
    posterior<-logTarget(dir,Model,x0,AlgoParams,expLL = F)
    print(posterior)
    print("...")
    return(posterior)
  }
  # Optimisation algorithm - Nelder & Mead, 1965, and includes outputting the Hessian
  output<-optim(par=x0,
                fn = Fopty,control = list(maxit = AlgoParams$itermax,
                                          fnscale=-1,
                                          reltol=1.5e-3),
                hessian = T)
  
  
  x0[!AlgoParams$indices]<-output$par
  # Convert the optimised values into physical ones
  x0%<>%Proposed2Physical(Model)
  
  return(list(PhysicalValues=x0,OptimisationOut=output,optimisedParams=(1:length(unlist(x0)))[!AlgoParams$indices]))
  
}


# Using the bisection method to generate a new guess
OptimMd<-function(FFF,x,LLs){
  
  Ldiff<-xdiff<-1
  while(Ldiff>1e-4 & xdiff>1e-6){
    xp<-c(mean(x[c(1,2)]),mean(x[c(2,3)]))
    LLxp<-c(FFF(xp[1]),FFF(xp[2]))
    LLt<-c(LLs[1],LLxp[1],LLs[2],LLxp[2],LLs[3]);  xt<-c(x[1],xp[1],x[2],xp[2],x[3])
    Ldiff<-min(abs(diff(LLt))); xdiff<-min(abs(diff(xt)))
    ix<-which.max(LLt)
    if(any(is.na(x))) {print(c(x,LLs));print(c(xt,LLt))}
    x<-xt[c(ix-1,ix,ix+1)];  LLs<-LLt[c(ix-1,ix,ix+1)]
    
  }
  
  ix<-which.max(LLs)
  print(paste0("Modifier = ",x[ix],",    LLmd = ",LLs[ix]))
  
  return(x[ix])
  
}

GetModifier<-function(Fopty,init,lower,upper){
  
  skely<-rep(0,length(init))
  
  for(i in 1:length(init)){
    
    Ftmp<-function(x) {
      y<-skely; y[i]<-x
      as.numeric((Fopty(y))[i])
    }
    
    if(init[i]>0){ # Calculate values between 0 and upper
      xx<-(0:20)*upper/20
      LLMd<-vapply(xx,Ftmp,numeric(1))
      up<-upper
      while(which.max(LLMd)==length(xx)){
        up<-up*1.5
        xx<-0.6*up+0.4*up*1:10/10
        LLMd<-vapply(xx,Ftmp,numeric(1))
      }
      srtr<-sort(LLMd,decreasing = T,index.return = T); srtr$x<-srtr$x[1:3];  srtr$ix<-srtr$ix[1:3]
      srtr2<-sort(xx[srtr$ix[1:3]],index.return = T);  xxx<-srtr2$x;  LLs<-srtr$x[srtr2$ix]
      skely[i] <-OptimMd(Ftmp,xxx,LLs)
    } else { # Calculate values between lower and 0
      xx<-(0:20)*lower/20
      LLMd<-vapply(xx,Ftmp,numeric(1))
      lo<-lower
      while(which.max(LLMd)==length(xx)){
        lo<-lo*1.5
        xx<-0.6*lo+0.4*lo*1:10/10
        LLMd<-vapply(xx,Ftmp,numeric(1))
      }
      srtr<-sort(LLMd,decreasing = T,index.return = T); srtr$x<-srtr$x[1:3];  srtr$ix<-srtr$ix[1:3]
      srtr2<-sort(xx[srtr$ix[1:3]],index.return = T);  xxx<-srtr2$x;  LLs<-srtr$x[srtr2$ix]
      skely[i] <-OptimMd(Ftmp,xxx,LLs)
    }
    
  }
  
  return(as.list(skely))
  
}

SingleEventsModifierCalc<-function(dir,Model,Omega,AlgoParams){
  
  warning("You need to set the subnational linear predictor variables to their average value before applying the modifiers")
  
  modifiers<-data.frame()
  ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/ODDobjects"),pattern=Model$haz,recursive = T,ignore.case = T)
  
  if(Model$haz=="EQ"){
    lower<--0.5
    upper<-3
  } else stop("Insert lower and upper bounds for the hazard studied, only have EQ so far")
  
  for(i in 1:length(ufiles)){
    
    print(ufiles[i])
    
    ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/",ufiles[i]))
    ODDy@gmax%<>%arrange(desc(gmax))
    
    nisos<-length(unique(ODDy@gmax$iso3))
    linp<-rep(list(0.),nisos)
    names(linp)<-unique(ODDy@gmax$iso3)
    ODDy@modifier<-linp
      
    tLL<-tryCatch(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams),
                  error=function(e) NA)
    
    if(is.na(tLL)) stop(ufiles[i])
    
    init<-log(tLL@predictDisp$gmax/tLL@predictDisp$disp_predictor)*0.1
    init[init<lower]<-0.5*lower; init[init>upper]<-0.5*upper; 
    init%<>%as.list()
    names(init)<-as.character(tLL@predictDisp$iso3)
  
    Fopty<-function(Md){
      
      Md%<>%as.list()
      names(Md)<-as.character(ODDy@gmax$iso3)
      ODDy@modifier<-Md
      
      lTarg<-logTargetSingle(ODDy,Model,Omega,AlgoParams)+
        Model$HighLevelPriors(Omega,Model,modifier = list(Md))
      
      print(paste0(ufiles[i],"   :   lTarg={",strcat(as.character(lTarg),collapse = ","),"},    params= {",strcat(as.character(Md),collapse = ","),"}"))
      return(lTarg)
    }
    
    fin<-tryCatch(GetModifier(Fopty,init,lower,upper),error=function(e) NA)
    if(is.na(fin)) {
      print(paste0("Failed to calculate modifier for ",ufiles[i]))
      print("...")
      next
    }
    print(paste0("Final params= {",strcat(as.character(unlist(fin)),collapse = ","),"}"))
    print("...")
    names(fin)<-as.character(ODDy@gmax$iso3)
    ODDy@modifier<-fin
    
    ODDy<-tryCatch(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams),
                  error=function(e) NA)
    if(is.na(ODDy)) stop(ufiles[i])
    
    modifiers%<>%rbind(data.frame(iso3=ODDy@predictDisp$iso3,modifier=as.numeric(unlist(fin)),
                                  event=ufiles[i],gmax=ODDy@gmax$gmax,
                                  eventid=as.numeric(strsplit(as.character(ufiles[i]),split = "_")[[1]][2]),
                                  predictor=ODDy@predictDisp$disp_predictor))
    
    saveRDS(ODDy, paste0(dir,"IIDIPUS_Results/ODDobjects/",ufiles[i]))
    saveRDS(modifiers,paste0(dir,"IIDIPUS_Results/ODDobjects/modifiers.Rdata"))

  }

  return(modifiers)
  
}

# Used to correlate the vulnerability variables with the modifier term
LMFeatureSelection<-function(output,Nb=12,intercept=F,fn="+",nlim=3,weights=NULL){
  # Nb - Number of times LM is run with different samples of training vs test data
  # intercept - do we include an intercept in the equation?
  
  # Use 80% of the observations as training and 20% for testing
  Ns <- floor(0.80*nrow(output))
  # List of variables in the output data.frame
  vars<-colnames(dplyr::select(output,-Y))
  # The parts of the LM that never change
  if(!intercept) eqn_base<-c("Y ~ ","+0") else eqn_base<-c("Y ~ ","")
  # Store the LM outputs
  prediction<-data.frame(eqn=NULL,AIC=NULL,BIC=NULL,adjR2=NULL)
  # How many variables to use at a time?
  for(n in 1:nlim){
    print(n)
    eqn<-combn(vars,n)
    if(n==1) {eqn%<>%as.character()
    }else eqn<-apply(eqn,2,function(x) pracma::strcat(x,collapse = fn))
    
    for(eee in eqn){
      equation<-paste0(eqn_base[1],eee,eqn_base[2])
      predictors<-data.frame(AIC=NULL,BIC=NULL)
      for (i in 1:Nb){
        ind = sample(seq_len(nrow(output)),size = Ns)
        train<-output[ind,]
        test<-output[-ind,]
        if(!is.null(weights)) predy<-lm(formula = as.formula(equation),
                                        data = train, weights = weights[ind])
        else predy<-lm(formula = as.formula(equation),data = train)
        # print(paste0(ceiling(AIC(predy)),", ",ceiling(BIC(predy))))
        # predme<-predict(predy,test, interval = 'confidence')
        predictors%<>%rbind(data.frame(AIC=AIC(predy),BIC=BIC(predy),adjR2=summary(predy)$adj.r.squared))
      }
      prediction<-rbind(cbind(data.frame(eqn=rep(equation)),t(colMeans(predictors))),prediction)
      
    }
  }
  
  return(prediction)
  
}

CorrelateModifier<-function(modifiers,Model){
  
  # This function uses all the ODDobjects available to extract iso3, start date and event_id values for each hazard
  DispData<-ExtractDispData_ODD(Model$haz)
  # Now extract the vulnerability variables to be used in the correlation - these are defined in Model
  val<-vulnerabilityVars(DispData,Model)
  # Use only the variables we think are relevant
  warning("Using only Hamish's pre-selected World Bank indicators for correlation")
  indies<-read_csv("./IIDIPUS_Input/REDUCED_WB-WorldDevelopment_Indicators.csv")
  val%<>%filter(indicator%in%indies$indicator_id)
  # Merge the two data frames
  modifiers$indicator<-"modifier"
  modifiers<-val%>%group_by(eventid)%>%summarise(date=unique(date))%>%merge(modifiers,by="eventid")
  # Housekeeping
  modifiers%<>%transmute(eventid=eventid,nearval=modifier,indicator=indicator)%>%
    dplyr::select(eventid,nearval); colnames(modifiers)[3]<-"modifier"
  
  # Get data frames in the correct structure
  for(inds in unique(val$indicator[val$indicator!="modifier"])){
    tmp<-filter(val,indicator==inds)
    tn<-data.frame(value=tmp$nearval); colnames(tn)<-inds
    modifiers%<>%cbind(tn)
  }
  rm(tn,ti,tmp)
  
  # Scale and center the vulnerability variables
  modifiers[3:ncol(modifiers)]<-base::scale(modifiers[3:ncol(modifiers)])
  # Remove all NA values
  redmodies<-redmodies[!apply(redmodies,1,function(x) any(is.na(x))),]
  # Remove all vulnerability variables that correlate strongly with one another - Variance Inflation Factor
  ncor <- cor(redmodies[,4:ncol(redmodies)])
  redmodies<-redmodies[,c(1,3,which(!1:ncol(redmodies)%in%caret::findCorrelation(ncor, cutoff=0.75)))]
  redmodies<-redmodies[,c(1,3,2,4:ncol(redmodies))]
  redmodies%<>%dplyr::select(-c(eventid))
  
  # Here we go!
  vulncor<-LMFeatureSelection(final,nlim = 3,weights = weights)
  # vulncor<-LMFeatureSelection(final,nlim = 3,weights = weights,fn = ":",intercept = T)
  # vulncor<-LMFeatureSelection(final,nlim = 4,weights = weights,fn = ":",intercept = T)
  
  return(vulncor)
  
}


# HighPriorModifiers<-function(modifiers,dir,Model,Omega,AlgoParams){
#   
#   hpmods<-modifiers
#   hpmods$hpmods<-hpmods$modifier
#   
#   for(i in 1:nrow(hpmods)){
#     
#     ODDy<-readRDS(paste0(dir,"IIDIPUS_Results/ODDobjects/",hpmods$event[i]))
#     
#     n<-15
#     
#     if(hpmods$modifier[i]< 0) {xx<-seq(-0.001,-0.01,length.out = n)
#     }else if(hpmods$modifier[i]> 0) {xx<-seq(0.01,1.75,length.out = n)
#     }else next
#     
#     Fopty<-function(Md){
#       
#       ODDy@modifier[[as.character(hpmods$iso3[i])]]<-Md
#       
#       LL<-logTargetSingle(ODDy,Model,Omega,AlgoParams) +
#           Model$HighLevelPriors(Omega,Model,modifier = list(Md))
# 
#       print(paste0(hpmods$event[i],"   :   LL={",strcat(as.character(LL),collapse = ","),"},    params= {",strcat(as.character(Md),collapse = ","),"}"))
#       return(LL[which(ODDy@predictDisp$iso3==as.character(hpmods$iso3[i]))])
#       
#     }
#     
#     LLout<-vapply(xx,Fopty,numeric(1))
#       
#     hpmods$hpmods[i] <- xx[which.max(LLout)]
#     
#   }
#   
#   
#   
# }


