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
AlgoParams<-list(#Data:
                 input_folder = 'IIDIPUS_Input/', # Directory of input data, allows you to vary between simulated input, aggregated input, etc.
                 
                 #Compute Parameters:
                 Np=1, # Number of Monte Carlo particles
                 cores=2, # Number of parallelised threads per event
                 NestedCores=1, # How many cores are to be used inside the ODD displacement calculations?
                 n_nodes=1, #parallelise between nodes on HPC
                 AllParallel=T, # Do you want fully-nested (data) parallelisation?
                 
                 # Loss function:
                 ABC=0, # Approximate Bayesian Computation rejection
                 kernel='energy_score', #Distance kernel between simulated and observed data
                 learning_rate = 40, #power of likelihood 
                 impact_weights=list(displacement=1,mortality=7,buildDam=0.8), #Relative weights of the different impact types when obtaining the energy score
                 m_CRPS = 100, # number of draws to estimate CRPS/Energy Score for each particle. Number of samples from model therefore becomes Np * m_CRPS
                 BD_weight = 0, #scales from 0 to 1: 0 means point data has no weighting, 1 means point data is included to same weighting (by inverse MAD) as mortality
                 log_offset = 10, #where log_offset = k, we compare log(x+k) with log(y+k) for simulated data x with observed data y
                 W_rankhist = 0.05, #weighting of Anderson Darling statistic over Band Depth rank histograms vs energy score
                 
                 # SMC hyperparameters: 
                 tol0=12000, # should be larger than largest expected distance
                 smc_steps = 200, #Number of steps in the ABC-SMC algorithm
                 smc_Npart = 1000, #Number of particles in the ABC-SMC algorithm
                 smc_alpha = 0.9, #Proportion of particles that we aim to keep 'alive' between steps
                 propCOV_multiplier = 0.2, #mutiplied by `optimal' covariance obtained using Filippi et al. 
    
                 # Adaptive MCMC hyperparameters:
                 N_steps = 10000, # number of MCMC steps
                 rho = 0.95, # correlation between subsequent event-level error terms
                 rho_local = 0.95, # correlation between subsequent local error terms
                 lambda_rate = 0.9, #parameter of adaptive scaling and shaping algorithm
                 alpha_star = 0.234, #parameter of adaptive scaling and shaping algorithm
                 v_0 = 15 #parameter of adaptive scaling and shaping algorithm
)
		 
if(is.null(AlgoParams$AllParallel)){
  if(AlgoParams$cores>4) { AlgoParams$AllParallel<-T
  } else AlgoParams$AllParallel<-F
}

# Choose the parameterisation algorithm - the string must match the function name exactly
Algorithm<- "AMCMC" #"ABCSMC"

###################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#----------------------------- GENERIC FUNCTIONS ---------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################################

# Metropolis-Hastings proposal distribution, given old values and covariance matrix
multvarNormProp <- function(xt, propPars){
  # purpose : A multivariate Gaussian random walk proposal for Met-Hastings
  #           MCMC
  # inputs  : xt       - The value of the chain at the previous time step 
  #           propPars - The correlation structure of the proposal
  
  #ensure positive semi definite:
  eig <- eigen(propPars)
  # Floor very small or negative eigenvalues
  eig$values[eig$values < 1e-6] <- 1e-6
  # Rebuild matrix
  propPars_psd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

  return(array(mvtnorm::rmvnorm(1, mean=xt, sigma=propPars_psd),dimnames = list(names(xt))))
}

Proposed2Physical<-function(proposed,Model,index=NULL){
  #Convert parameters on transformed space to the physical space:
  
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
  #Convert parameters on physical space to the transformed space:
  
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
  #Convert an unlabelled array of parameters on the transformed space to parameters on physical space
  return(array %>% relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model))
}

modifyAcc <- function(xNew, xPrev, Model, propCOV){
  #Modify the acceptance rate to account for Jacobian of parameter transformation and priors
  
  xNew %<>% Array2Physical(Model) %>% unlist() 
  xPrev %<>% Array2Physical(Model) %>% unlist()
  Model$acceptTrans%<>%unlist()
  index<-1:length(unlist(Model$unlinks))
  
  product <- 1
  #Account for jacobian of parameter transformation:
  for (i in index){
    product <- product * match.fun(Model$acceptTrans[[names(xNew)[i]]])(xNew[i], xPrev[i], Model$par_lb[i],Model$par_ub[i])
  }
  
  #Account for priors:
  index = 1
  for (i in 1:length(Model$Priors)){
    if (is.list(Model$Priors[[i]])){
      for (j in 1:length(Model$Priors[[i]])){ # probably a more efficient way of doing this
        if(Model$Priors[[i]][[j]]$dist == 'norm'){
          product = product * (dnorm(xNew[index], Model$Priors[[i]][[j]]$mean, Model$Priors[[i]][[j]]$sd)/ dnorm(xPrev[index], Model$Priors[[i]][[j]]$mean, Model$Priors[[i]][[j]]$sd))
        }
        if(Model$Priors[[i]][[j]]$dist == 'laplace'){
          product = product * (dlaplace(xNew[index], Model$Priors[[i]][[j]]$location, Model$Priors[[i]][[j]]$scale)/ dlaplace(xPrev[index], Model$Priors[[i]][[j]]$location, Model$Priors[[i]][[j]]$scale))
        }
        if(Model$Priors[[i]][[j]]$dist == 'invgamma'){
          product = product * (dinvgamma(xNew[index], shape=Model$Priors[[i]][[j]]$shape, scale=Model$Priors[[i]][[j]]$scale)/ dinvgamma(xPrev[index], shape=Model$Priors[[i]][[j]]$shape, scale=Model$Priors[[i]][[j]]$scale))
        }
        index = index + 1
      }
    } else {
      if(Model$Priors[[i]]$dist == 'norm'){
        product = product * (dnorm(xNew[index], Model$Priors[[i]]$mean, Model$Priors[[i]]$sd)/ dnorm(xPrev[index], Model$Priors[[i]]$mean, Model$Priors[[i]]$sd))
      }
      if(Model$Priors[[i]]$dist=='laplace'){
        product = product * (dLaplace(xNew[index], Model$Priors[[i]]$location, Model$Priors[[i]]$scale)/ dLaplace(xPrev[index], Model$Priors[[i]]$location, Model$Priors[[i]]$scale))
      }
      if(Model$Priors[[i]]$dist=='invgamma'){
        product = product * (dinvgamma(xNew[index], shape=Model$Priors[[i]]$shape, scale=Model$Priors[[i]]$scale)/ dinvgamma(xPrev[index], shape=Model$Priors[[i]]$shape, scale= Model$Priors[[i]]$scale))
      }
      index = index + 1
    }
  }
  return(product)
}

###################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#----------------------- ABC-SMC PARALLELISED ACROSS NODES -----------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################################

tpa<-function(tolerance, d){ 
  #calculate the total proportion of alive particles (have at least one distance less than the tolerance)
  #d is a matrix, with element (i,j) representing the distance calculated for the jth Monte Carlo data simulation of particle i. 
  Npart <- NROW(d)
  M <- NCOL(d)
  npa <- rep(0,Npart) 
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

bcast_ODDRIN <- function(dir, Model, AlgoParams, nslaves){
  #Broadcast ODDRIN functions and objects to the nodes
  mpi.spawn.Rslaves(nslaves=nslaves)
  mpi.bcast.Rfun2slave()
  mpi.bcast.Robj2slave(dir)
  directory = dir
  mpi.bcast.Robj2slave(directory)
  mpi.bcast.Robj2slave(Model)
  mpi.bcast.Robj2slave(AlgoParams)
  mpi.bcast.Robj2slave(HLPrior_sample)
  
  #Run GetODDPackages and Simulate.R on the nodes
  mpi.remote.exec(GetODDPackages, TRUE)
  mpi.remote.exec(source('RCode/Simulate.R'))
}

initialise_particles <- function(dir, Model, AlgoParams, AlgoResults){
  # Sample initial values for each particle from the prior, obtain distances, and store.
  # This function is used when AlgoParams$n_nodes = 1. When n_nodes > 1, the initialisation is parallelised
  # across nodes using initialise_particles_Rmpi().
  iter_times <- c()
  for (n in 1:AlgoParams$smc_Npart){
    print(paste('Step: 1, Particle:', n))
    
    d_prop <- array(AlgoParams$tol0+1, dim=c(1,1))
    while(any(d_prop[,1]> AlgoParams$tol0)){
      #draw initial particle values from the prior and ensure that they satisfy the higher level prior
      AlgoResults$Omega_sample[n,,1] <- HLPrior_sample(Model, AlgoParams)
      AlgoResults$Omega_sample_phys[n,,1] <-  AlgoResults$Omega_sample[n,,1] %>% relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model) %>% unlist()
      
      start_time <- Sys.time()
      impact_sample <- SampleImpact(dir = dir,Model = Model,
                                    proposed = AlgoResults$Omega_sample_phys[n,,1] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), 
                                    AlgoParams = AlgoParams,
                                    dat=dat)
      
      d_prop <- CalcDist(impact_sample, AlgoParams)
    }
     
    AlgoResults$d[n,,1] <- rowSums(d_prop)
    min.d.part <- which.min(AlgoResults$d[n,,1])
    
    if (n == 1){
      AlgoResults$d_full = array(NA, dim=c(AlgoParams$smc_Npart, AlgoParams$Np, NCOL(d_prop), AlgoParams$smc_steps))
    }
    AlgoResults$d_full[n,,,1] <- d_prop
  
    ##currently not storing sampled_full due to the space it takes up:
    
    # if (n == 1){
    #   AlgoResults$sampled_full = array(NA, dim=c(AlgoParams$smc_Npart, NROW(impact_sample$poly[[1]]), AlgoParams$m_CRPS, AlgoParams$smc_steps))
    # }
    #AlgoResults$sampled_full[n,,,1] = matrix(unlist(lapply(impact_sample$poly[((min.d.part-1)*AlgoParams$m_CRPS+1):(min.d.part*AlgoParams$m_CRPS)], function(x) x$sampled)), nrow=NROW(impact_sample$poly[[1]])) #take sampled values for the sample with the smallest distance
    
    end_time <- Sys.time()
    
    iter_times <- append(iter_times, end_time-start_time)
  }
  AlgoResults$iter_times <- iter_times
  return(AlgoResults)
}

initialise_particles_Rmpi <- function(dir, Npart, n_nodes){
  # Sample initial values for each particle from the prior, obtain distances, and store.
  # This function is used when AlgoParams$n_nodes > 1. When n_nodes = 1, the initialisation is 
  # performed on a single node using initialise_particles().
  
  # Use mpi.comm.rank() to determine which particles the node has been allocated:
  particle_divisions <- split(1:Npart, sort(1:Npart%%n_nodes))
  allocated_particles <- particle_divisions[[mpi.comm.rank()]]
  n_allocated <- length(allocated_particles)
  
  n_x <- length(Model$par_lb) # dimension of parameter space
  
  #create empty arrays to store parameter values and distances:
  Omega_sample_node <- array(NA, dim=c(n_allocated, n_x))
  Omega_sample_phys_node <- array(NA, dim=c(n_allocated, n_x))
  d_node <- array(Inf, dim=c(n_allocated, AlgoParams$Np))
  
  iter_times <- c()
  for (n in 1:n_allocated){
    #Sample from the parameter ranges until the higher level prior is satisfied:
    d_prop <- array(AlgoParams$tol0+1, dim=c(1,1))
    while(any(d_prop[,1]> AlgoParams$tol0)){
      Omega_sample_node[n,] <- HLPrior_sample(Model, AlgoParams)
      Omega_sample_phys_node[n,] <-  Omega_sample_node[n,] %>% relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model) %>% unlist()
      start_time <- Sys.time()
      
      #calculate distance
      impact_sample <- SampleImpact(dir = dir, Model = Model,
                                    proposed = Omega_sample_phys_node[n,] %>% relist(skeleton=Model$skeleton) %>% addTransfParams(), 
                                    AlgoParams = AlgoParams)
      d_prop = CalcDist(impact_sample, AlgoParams) 
    }
    d_node[n,] = rowSums(d_prop)
    
    # min.d.part <- which.min(d_node[n,])
    # 
    if (n == 1){
       d_full = array(NA, dim=c(n_allocated, AlgoParams$Np, NCOL(d_prop)))
    #   sampled_full = array(NA, dim=c(n_allocated, NROW(impact_sample$poly[[1]]), AlgoParams$m_CRPS))
    }
    d_full[n,,] <- d_prop
    #sampled_full[n,,] = matrix(unlist(lapply(impact_sample$poly[((min.d.part-1)*AlgoParams$m_CRPS+1):(min.d.part*AlgoParams$m_CRPS)], function(x) x$sampled)), nrow=NROW(impact_sample$poly[[1]])) #take sampled values for the sample with the smallest distance
    
    end_time <- Sys.time()
    iter_times <- append(iter_times, end_time-start_time)
    
  }
  return(list(Omega_sample_node=Omega_sample_node,
              Omega_sample_phys_node=Omega_sample_phys_node,
              d_node=d_node, 
              iter_times=iter_times, 
              d_full_node=d_full,
              sampled_full_node=NULL))#sampled_full))
}


AlgoStep1 <- function(dir, AlgoParams, AlgoResults){
  #Initialise the full set of particles at the first step of the ABC-SMC algorithm, using either:
  #         - initialise_particles() when AlgoParams$n_nodes = 1
  #         - initialise_particles_Rmpi() when AlgoParams$n_nodes > 1
  
  AlgoResults$W[,1] <- 1/AlgoParams$smc_Npart # Give each particle an equal weight
  start_time <- Sys.time()
  
  if(AlgoParams$n_nodes > 1){ #parallelise across nodes
    node_return <- mpi.remote.exec(initialise_particles_Rmpi, dir, AlgoParams$smc_Npart, AlgoParams$n_nodes)
    #print(node_return)
    particle_divisions <- split(1:AlgoParams$smc_Npart, sort(1:AlgoParams$smc_Npart%%AlgoParams$n_nodes))
    AlgoResults$d_full = array(Inf, dim=c(AlgoParams$smc_Npart, length(node_return[[1]]$d_full_node[1,,1]), length(node_return[[1]]$d_full_node[1,1,]), AlgoParams$smc_steps))
    # AlgoResults$sampled_full = array(NA, dim=c(AlgoParams$smc_Npart, length(node_return[[1]]$sampled_full_node[1,,1]), length(node_return[[1]]$sampled_full_node[1,1,]), AlgoParams$smc_steps))
    
    for (j in 1:length(node_return)){
      AlgoResults$Omega_sample[particle_divisions[[j]],,1] <- node_return[[j]]$Omega_sample_node
      AlgoResults$Omega_sample_phys[particle_divisions[[j]],,1] <- node_return[[j]]$Omega_sample_phys_node
      AlgoResults$d[particle_divisions[[j]],,1] <- node_return[[j]]$d_node
      AlgoResults$d_full[particle_divisions[[j]],,,1] <- node_return[[j]]$d_full_node
      #AlgoResults$sampled_full[particle_divisions[[j]],,,1] <- node_return[[j]]$sampled_full_node
    }
    
  } else { #no parallelisation
    AlgoResults <- initialise_particles(dir, Model, AlgoParams, AlgoResults)
    print(AlgoResults$iter_times)
  }
  
  AlgoResults$tolerancestore[1] <- max(AlgoResults$d[,,1]) + 1 # set tolerance to larger than maximum distance
  AlgoResults$accrate_store[1] <- 1
  end_time <- Sys.time()
  
  ## Currently just applying an equal weight of 1 to all elements of the distance
  rel_weights <- rep(1, length(AlgoResults$d_full[1,1,,1]))
  AlgoResults$rel_weights = rel_weights

  d_full_weighted <- sweep(adrop(AlgoResults$d_full[,,,1, drop=F], drop=4), 3, rel_weights, FUN = "*")
  AlgoResults$d[,,1] <- apply(d_full_weighted, c(1,2), sum)
  
  ## Alternatively, weight distances by inverse MAD at first step
  # rel_weights <- 1/apply(AlgoResults$d_full[,1,,1], 2, stats::mad)
  # rel_weights <- ifelse(is.finite(rel_weights), rel_weights, 0)
  
  return(AlgoResults)
}

perturb_particles <- function(s, propCOV, AlgoParams, tolerance_s, W_s, Omega_sample_s, Omega_sample_phys_s,
                              d_s, d_full_s, sampled_full_s, rel_weights){
  # For each particle:
  #     - Perturb using the proposal covariance
  #     - For the proposed particle, sample from the model and calculate the distance to the observed data
  #     - Accept or reject the proposal
  # Return the set of perturbed (or non-perturbed, in the case of rejection) particles and distances
  
  for(n in 1:AlgoParams$smc_Npart){
    print(paste(' Step:', s, ', Particle:', n))
    if(W_s[n]>0){
      Omega_prop <- multvarNormProp(xt=Omega_sample_s[n,], propPars=propCOV) #perturb the proposal
      Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
      
      HP<- Model$HighLevelPriors(Omega_prop_phys %>% addTransfParams(), Model)
      if (HP> AlgoParams$ABC & Model$higherpriors) next
      
      impact_sample <- SampleImpact(dir = dir,Model = Model,
                                proposed = Omega_prop_phys %>% addTransfParams(), 
                                AlgoParams = AlgoParams)
      d_prop <- CalcDist(impact_sample, AlgoParams)
      if (AlgoParams$Np > 1){
        d_prop_tot <- rowSums(sweep(d_prop, 2, rel_weights, FUN = "*"))
      } else {
        d_prop_tot <- sum(d_prop * rel_weights)
      }
      
      print(d_prop_tot)
      
      if(d_prop[1]==Inf){#if (d_full_prop[1]==Inf){
        d_prop <- Inf
      } 
      
      d_curr <- d_s[n,which(is.finite(d_s[n,]))]
      
      #calculate the acceptance probability:
      acc <- (sum(d_prop_tot<tolerance_s)/length(d_prop_tot))/(sum(d_curr<tolerance_s)/length(d_curr)) * modifyAcc(Omega_prop, Omega_sample_s[n,], Model, propCOV)
        
      u <- runif(1)
      if(u < acc){
        min.d.part <- which.min(d_prop_tot)
        Omega_sample_s[n,] <- Omega_prop
        Omega_sample_phys_s[n,] <- Omega_sample_s[n,] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
        #sampled_full_s[n,,] <- matrix(unlist(lapply(impact_sample$poly[((min.d.part-1)*AlgoParams$m_CRPS+1):(min.d.part*AlgoParams$m_CRPS)], function(x) x$sampled)), nrow=NROW(impact_sample$poly[[1]])) #take sampled values for the sample with the smallest distance
        d_full_s[n,,] <- d_prop
        d_s[n,] <- d_prop_tot
      }
      #saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/abcsmc_",tag))
    }
  }
  return(list(Omega_sample_s = Omega_sample_s,
              Omega_sample_phys_s = Omega_sample_phys_s,
              d_full_s = d_full_s,
              sampled_full_s = sampled_full_s,
              d_s = d_s))
}


perturb_particles_Rmpi <- function(dir, Npart, n_nodes, W_curr, Omega_curr, Omega_phys_curr, d_curr, d_full_curr, sampled_full_curr, propCOV, tolerance_s, rel_weights){
  # The MPI alternative to perturb_particles()
  # For each particle:
  #     - Perturb using the proposal covariance
  #     - For the proposed particle, sample from the model and calculate the distance to the observed data
  #     - Accept or reject the proposal
  # Return the set of perturbed (or non-perturbed, in the case of rejection) particles and distances
  
  #Determine which particles have been allocated to this node:
  particle_divisions <- split(1:Npart, sort(1:Npart%%n_nodes))
  allocated_particles <- particle_divisions[[mpi.comm.rank()]]
  n_allocated <- length(allocated_particles)
  n_x <- length(Model$par_lb)
  
  W_node <- W_curr[allocated_particles]
  Omega_sample_node <- Omega_curr[allocated_particles,, drop=F]
  Omega_sample_phys_node <- Omega_phys_curr[allocated_particles,, drop=F]
  d_node <- d_curr[allocated_particles,, drop=F]
  d_full_node <- d_full_curr[allocated_particles,,, drop=F]
  #sampled_full_node <- sampled_full_curr[allocated_particles,,,drop=F]
  
  for(n in 1:n_allocated){
    if(W_node[n]>0){
      Omega_prop <- multvarNormProp(xt=Omega_sample_node[n,], propPars=propCOV) #perturb the proposal
      Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
      
      if (any(unlist(Omega_prop_phys) < Model$par_lb) | any(unlist(Omega_prop_phys) > Model$par_ub)) next
      
      HP<- Model$HighLevelPriors(Omega_prop_phys %>% addTransfParams(), Model)
      if (HP> AlgoParams$ABC & Model$higherpriors) next
      
      impact_sample <-  SampleImpact(dir = dir, Model = Model,
                            proposed = Omega_prop_phys %>% addTransfParams(), 
                            AlgoParams = AlgoParams) 
      d_prop <- CalcDist(impact_sample, AlgoParams)
      
      if(AlgoParams$Np > 1){
        d_prop_tot <- rowSums(sweep(d_prop, 2, rel_weights, FUN = "*"))
      } else {
        d_prop_tot <- sum(d_prop * rel_weights)
      }
     
      if(d_prop[1]==Inf){
        d_prop <- Inf
      } 
      d_curr <- d_node[n,which(is.finite(d_node[n,]))]
      
      acc <- (sum(d_prop_tot<tolerance_s)/length(d_prop_tot))/(sum(d_curr<tolerance_s)/length(d_curr)) * modifyAcc(Omega_prop, Omega_sample_node[n,], Model)
      u <- runif(1)
      if(u < acc){
        min.d.part <- which.min(d_prop_tot)
        Omega_sample_node[n,] <- Omega_prop
        Omega_sample_phys_node[n,] <- Omega_sample_node[n,] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
        d_node[n,] <- d_prop_tot
        #sampled_full_node[n,,] <- matrix(unlist(lapply(impact_sample$poly[((min.d.part-1)*AlgoParams$m_CRPS+1):(min.d.part*AlgoParams$m_CRPS)], function(x) x$sampled)), nrow=NROW(impact_sample$poly[[1]])) #take sampled values for the sample with the smallest distance
        d_full_node[n,,] <- d_prop
      }
    }
  }
  return(list(Omega_sample_node=Omega_sample_node, 
              Omega_sample_phys_node=Omega_sample_phys_node,
              d_node=d_node, 
              d_full_node=d_full_node,
              sampled_full_node=NULL))#sampled_full_node))
}

retrieve_UnfinishedAlgoResults <- function(dir, oldtag, Npart, AlgoResults){
  #new AlgoResults must have same M_CRPS and smc_Npart as old one
  AlgoResults_unfinished <- readRDS(paste0(dir,"IIDIPUS_Results/",oldtag))
  s_finish = max(which(colSums(is.na(AlgoResults_unfinished$Omega_sample[,1,]))==0)) #find last completed step
  s_start <- s_finish 
  
  # AlgoResults_unfinished <- readRDS(paste0(dir,"IIDIPUS_Results/",oldtag))
  # s_finish = max(which(colSums(is.na(AlgoResults_unfinished$Omega_sample[,1,]))==0))-1 #find last completed step
  # #n_finish = max(which(!is.na(AlgoResults_unfinished$Omega_sample[,1,s_finish]))) #find last completed particle
  # 
  ## carry out these steps in case s differs between the unfinished and current runs
  AlgoResults$W[,1:s_finish] <- AlgoResults_unfinished$W[,1:s_finish]
  AlgoResults$d[,,1:s_finish] <- AlgoResults_unfinished$d[,,1:s_finish]
  AlgoResults$d_full <- array(NA, dim=c(dim(AlgoResults_unfinished$d_full)[1:3], dim(AlgoResults$W)[2]))
  AlgoResults$sampled_full <- NULL #array(NA, dim=c(dim(AlgoResults_unfinished$sampled_full)[1:3], dim(AlgoResults$W)[2]))
  AlgoResults$d_full[,,,1:s_finish] <- AlgoResults_unfinished$d_full[,,,1:s_finish]
  AlgoResults$sampled_full[,,,1:s_finish] <- AlgoResults_unfinished$sampled_full[,,,1:s_finish]
  AlgoResults$Omega_sample[,,1:s_finish] <- AlgoResults_unfinished$Omega_sample[,,1:s_finish]
  AlgoResults$Omega_sample_phys[,,1:s_finish] <- AlgoResults_unfinished$Omega_sample_phys[,,1:s_finish]
  AlgoResults$tolerancestore[1:s_finish] <- AlgoResults_unfinished$tolerancestore[1:s_finish]
  AlgoResults$essstore[1:s_finish] <- AlgoResults_unfinished$essstore[1:s_finish]
  AlgoResults$propCOV[,,1:s_finish] <- AlgoResults_unfinished$propCOV[,,1:s_finish]
  AlgoResults$accrate_store[1:s_finish] <- AlgoResults_unfinished$accrate_store[1:s_finish]
  AlgoResults$rel_weights <- AlgoResults_unfinished$rel_weights
  AlgoResults$input_folder <- AlgoResults_unfinished$input_folder
  # s_start = s_finish + 1 #ifelse(n_finish==Npart, s_finish+1, s_finish) #identify the appropriate step from which to continue the algorithm
  # #n_start = ifelse(n_finish==Npart, 1, n_finish + 1) #identify the appropriate particle from which to continue the algorithm
  # 
  return(list(
    AlgoResults=AlgoResults,
    s_start=s_start
  ))
}

update_tolerance_and_weights <- function(s, alpha, AlgoResults){
  #Find the new tolerance such that alpha proportion of the current alive particles stay alive
  toleranceold <- AlgoResults$tolerancestore[s-1]
  d_old <- adrop(AlgoResults$d[,,s-1,drop=FALSE], drop = 3)
  reflevel <- alpha * tpa(toleranceold, d_old)

  tolerance<-uniroot(function(tolerance) tpa(tolerance,d=d_old)-reflevel,c(0,toleranceold*3))$root
  print(paste('      Step:', s, ', New tolerance is:', tolerance))

  AlgoResults$tolerancestore[s] <- tolerance
  AlgoResults$essstore[s-1]<- 1/sum(AlgoResults$W[,s-1]^2) #effective sample size
  
  #compute the associated weights
  npa_old<-rowSums(d_old<toleranceold)/apply(d_old, 1, function(x){sum(is.finite(x))})
  npa<-rowSums(d_old<tolerance)/apply(d_old, 1, function(x){sum(is.finite(x))})
  a<-which(npa_old>0)
  b<-which(npa_old==0)
  AlgoResults$W[a,s]<-AlgoResults$W[a,s-1]*npa[a]/npa_old[a]
  AlgoResults$W[b,s]<-rep(0,length(b))
  AlgoResults$W[,s]<-AlgoResults$W[,s]/sum(AlgoResults$W[,s])
  
  return(AlgoResults)
}

resample_particles <- function(s, N_T, Npart, AlgoResults){
  # resample if effective sample size is too low
  if ((sum(AlgoResults$W[,s]^2)*N_T)>1){
    choice<- sample(1:Npart,Npart,replace= TRUE, prob = AlgoResults$W[,s])
    AlgoResults$Omega_sample[,,s] <- AlgoResults$Omega_sample[choice,,s-1] 
    AlgoResults$Omega_sample_phys[,,s] <- AlgoResults$Omega_sample_phys[choice,,s-1] 
    AlgoResults$d[,,s] <- adrop(AlgoResults$d[choice,,s-1, drop=F], drop=3)
    AlgoResults$d_full[,,,s] <- adrop(AlgoResults$d_full[choice,,,s-1, drop=F], drop=4)
    #AlgoResults$sampled_full[,,,s] <- adrop(AlgoResults$sampled_full[choice,,,s-1, drop=F], drop=4)
    AlgoResults$W[,s] <-rep(1/Npart,Npart) 
  } else { #otherwise do not resample (the particles will be perturbed later via a MCMC step) 
    AlgoResults$Omega_sample[,,s] <- AlgoResults$Omega_sample[,,s-1] 
    AlgoResults$Omega_sample_phys[,,s] <- AlgoResults$Omega_sample_phys[,,s-1] 
    AlgoResults$d[,,s] <- adrop(AlgoResults$d[,,s-1, drop=F], drop=3)
    AlgoResults$d_full[,,,s] <- adrop(AlgoResults$d_full[,,,s-1, drop=F], drop=4)
    #AlgoResults$sampled_full[,,,s] <- adrop(AlgoResults$sampled_full[,,,s-1, drop=F], drop=4)
  }
  return(AlgoResults)
}

calc_propCOV <- function(s, n_x, Npart, AlgoResults){
  #calculate perturbation covariance based on Filippi et al., 2012
  
  tilda_i <- which(rowSums(adrop(AlgoResults$d[,,s-1, drop=F], drop=3)<AlgoResults$tolerancestore[s])>0)
  Omega_tilda <- AlgoResults$Omega_sample[tilda_i,,s-1] 
  W_tilda <- AlgoResults$W[tilda_i,s-1]
  W_tilda <- W_tilda/sum(W_tilda) #normalise weights
  
  # When Npart > 2000, sometimes calculate the covariance just using the first 2000 to speed things up a bit!
  
  #check that the indexes are right here!
  propCOV <- matrix(0, nrow=n_x, ncol=n_x)
  #tilda_i <- tilda_i[tilda_i<2000]
  #for (n in 1:min(2000, Npart)){
  for (n in 1:Npart){
    for(k in 1:length(tilda_i)){
      propCOV <- propCOV + AlgoResults$W[n,s] * W_tilda[k] * ((Omega_tilda[k,]-AlgoResults$Omega_sample[n,,s]) %*% t(Omega_tilda[k,]-AlgoResults$Omega_sample[n,,s]))
    }
  }
  
  ## Local alternative (different covariance per particle), but no longer reversible steps so doesn't seem to work:
  # propCOV <- list()
  # for (n in 1:Npart){
  #   propCOV[[n]] <- matrix(0, n_x, n_x)
  #   for(k in 1:length(tilda_i)){
  #     propCOV[[n]] <- propCOV[[n]] + W_tilda[k] * ((Omega_tilda[k,]-AlgoResults$Omega_sample[n,,s]) %*% t(Omega_tilda[k,]-AlgoResults$Omega_sample[n,,s]))
  #   }
  # }
  
  return(propCOV)
}

plot_abcsmc <- function(s, n_x, Npart, Omega_sample_phys, Omega, accrate=NULL){
  par(mfrow=c(5,4))
  for (i in 1:min(n_x, 20)){
    plot(Omega_sample_phys[,i,s], main=names(unlist(Omega))[i])
    abline(h=unlist(Omega)[i], col='red')
  }
  par(mfrow=c(1,1))
}

ABCSMC <- function(dir, AlgoParams, Model, unfinished=F, oldtag=NULL, tag_notes=NULL){
  #Input: 
  # - AlgoParams: Parameters describing the ABC-SMC Algorithm (e.g. the ABC rejection threshold for higher level priors)
  # - Model: Describes the data simulation and calculation of the distance measure 
  # - Unfinished: Continue an old MCMC chain or start a new one?
  # -                If TRUE, then include oldtag - the tag (end of the filename) of an unfinished ABC-SMC run to be completed.
  # - Tag_notes: Additional extention to the results file name (e.g. describing hyperparameters, dataset)
  # Output:
  # - A list containing the parameter values for each particle at each step of the ABC-SMC algorithm, as well as the
  #   corresponding weights, distances and final perturbation kernel
  # Details:
  # - Uses the algorithm described in Del Moral et al., 2011: https://link.springer.com/content/pdf/10.1007/s11222-011-9271-y.pdf
  # - Uses the multivariate perturbation kernel described in Filippi et al., 2012: https://arxiv.org/pdf/1106.6280.pdf
  
  n_x <- length(Model$par_lb) #n_x = number of parameters
  N_T <- AlgoParams$smc_Npart/2
  tolerancetarget=1 # At this point, no defined stopping point
  
  # Filename:
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  tag<-ifelse(is.null(tag_notes), tag, paste0(tag, '_', tag_notes))
  
  AlgoResults <- list(
    input_folder = AlgoParams$input_folder,
    Omega_sample = array(NA, dim=c(AlgoParams$smc_Npart, n_x, AlgoParams$smc_steps)), #store sampled parameters on the transformed space
    Omega_sample_phys = array(NA, dim=c(AlgoParams$smc_Npart, n_x, AlgoParams$smc_steps)), #store sampled parameters on the untransformed space
    W = array(NA, dim=c(AlgoParams$smc_Npart, AlgoParams$smc_steps)), #Weights
    d = array(Inf, dim=c(AlgoParams$smc_Npart, AlgoParams$Np, AlgoParams$smc_steps)), #Distances
    d_full = NULL, #array(Inf, dim=c(AlgoParams$smc_Npart, 3 * AlgoParams$Np, AlgoParams$smc_steps)), #Distances broken down into: poly_crps, poly_mean, point
    sampled_full = NULL,
    tolerancestore=array(NA, AlgoParams$smc_steps),
    essstore=array(NA, AlgoParams$smc_steps),
    accrate_store=array(NA, AlgoParams$smc_steps),
    propCOV=array(NA, dim=c(n_x, n_x, AlgoParams$smc_steps))
  )
  
  if (AlgoParams$n_nodes > 1) bcast_ODDRIN(dir, Model, AlgoParams, nslaves=AlgoParams$n_nodes)
  start_time <- Sys.time()
  
  if(unfinished==F){ 
    #Initialize and perform sampling for s=1
    AlgoResults <- AlgoStep1(dir, AlgoParams, AlgoResults)
    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/abcsmc_",tag))  
    s_start = 2 # continue the algorithm from s = 2
  } else { 
    #Collect relevant information from the unfinished sample
    UnfinishedAlgoResults <- retrieve_UnfinishedAlgoResults(dir, oldtag, AlgoParams$smc_Npart, AlgoResults)
    AlgoResults <- UnfinishedAlgoResults$AlgoResults
    s_start <- UnfinishedAlgoResults$s_start
    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/abcsmc_",tag)) 
  }
  
  for (s in s_start:AlgoParams$smc_steps){
    
    #record and print time for each step
    end_time <- Sys.time()
    print(paste('Time:', end_time-start_time))
    start_time <- Sys.time()
    
    if (s>2){ AlgoParams$smc_alpha <- min(0.99,max(0.9,(1 - AlgoResults$accrate_store[s-1]))) }
    AlgoResults <- update_tolerance_and_weights(s, AlgoParams$smc_alpha, AlgoResults)
    AlgoResults <- resample_particles(s, N_T, AlgoParams$smc_Npart, AlgoResults)
    propCOV <- calc_propCOV(s, n_x, AlgoParams$smc_Npart, AlgoResults) * AlgoParams$propCOV_multiplier
    
    AlgoResults_small <- list(
      tolerance_s = AlgoResults$tolerancestore[s],
      Omega_sample_s = AlgoResults$Omega_sample[,,s],
      Omega_sample_phys_s = AlgoResults$Omega_sample_phys[,,s],
      d_s = adrop(AlgoResults$d[,,s, drop=F], drop=3),
      d_full_s = adrop(AlgoResults$d_full[,,,s, drop=F], drop=4),
      sampled_full_s = NULL, #adrop(AlgoResults$sampled_full[,,,s, drop=F], drop=4),
      W_s = AlgoResults$W[,s], 
      rel_weights = AlgoResults$rel_weights
    )
    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/abcsmc_start_step_",tag))
    rm(AlgoResults) #free up some space while doing the heavy lifting
    
    if(AlgoParams$n_nodes>1){
      node_return <- mpi.remote.exec(perturb_particles_Rmpi, dir, AlgoParams$smc_Npart, AlgoParams$n_nodes, 
                                     AlgoResults_small$W_s, AlgoResults_small$Omega_sample_s, AlgoResults_small$Omega_sample_phys_s, 
                                     AlgoResults_small$d_s, AlgoResults_small$d_full_s, AlgoResults_small$sampled_full_s, propCOV, AlgoResults_small$tolerance_s, 
                                     rel_weights = AlgoResults_small$rel_weights)
      #print(node_return)
      particle_divisions <- split(1:AlgoParams$smc_Npart, sort(1:AlgoParams$smc_Npart%%AlgoParams$n_nodes))
      
      AlgoResults <- readRDS(paste0(dir,"IIDIPUS_Results/abcsmc_start_step_",tag))
      for (j in 1:length(node_return)){
        #print(node_return[[j]])
        AlgoResults$Omega_sample[particle_divisions[[j]],,s] <- node_return[[j]]$Omega_sample_node
        AlgoResults$Omega_sample_phys[particle_divisions[[j]],,s] <- node_return[[j]]$Omega_sample_phys_node
        AlgoResults$d[particle_divisions[[j]],,s] <- node_return[[j]]$d_node
        AlgoResults$d_full[particle_divisions[[j]],,,s] <- node_return[[j]]$d_full_node
        #AlgoResults$sampled_full[particle_divisions[[j]],,,s] <- node_return[[j]]$sampled_full_node
      }
    } else {
      step_s_results <- perturb_particles(s, propCOV, AlgoParams, AlgoResults_small$tolerance_s, AlgoResults_small$W_s, AlgoResults_small$Omega_sample_s, AlgoResults_small$Omega_sample_phys_s,
                        AlgoResults_small$d_s, AlgoResults_small$d_full_s, AlgoResults_small$sampled_full_s, rel_weights = AlgoResults_small$rel_weights)
      AlgoResults <- readRDS(paste0(dir,"IIDIPUS_Results/abcsmc_start_step_",tag))
      AlgoResults$Omega_sample[,,s] <- step_s_results$Omega_sample_s
      AlgoResults$Omega_sample_phys[,,s] <- step_s_results$Omega_sample_phys_s
      AlgoResults$d[,,s] <- step_s_results$d_s
      AlgoResults$d_full[,,,s] <- step_s_results$d_full_s
      #AlgoResults$sampled_full[,,,s] <- step_s_results$sampled_full_s
    }
    AlgoResults$propCOV[,,s] <- propCOV
    AlgoResults$accrate_store[s] <- mean(AlgoResults$Omega_sample_phys[which(AlgoResults$W[,s]>0),1,s]!=AlgoResults_small$Omega_sample_phys_s[which(AlgoResults$W[,s]>0),1])

    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/abcsmc_",tag))
    
    print(s)
    
    #plot_abcsmc(s, n_x, AlgoParams$smc_Npart, AlgoResults$Omega_sample_phys, Omega, accrate=AlgoResults$accrate_store)
    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/abcsmc_",tag))  
  }
  
  if (AlgoParams$n_nodes > 1) mpi.close.Rslaves()
  return(AlgoResults)
}

###################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#---------Adaptive MCMC with accelerated scaling and shaping ---------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################################


retrieve_UnfinishedAlgoResults_AMCMC <- function(dir, oldtag, AlgoResults){
  #retrieve an unfinished AMCMC chain so that it can be continued running:
  
  AlgoResults_unfinished <- readRDS(paste0(dir,"IIDIPUS_Results/",oldtag))
  s_finish = max(which(!is.na(AlgoResults_unfinished$Omega_sample[1,]))) #find last completed step
  s_start <- s_finish
  
  # # carry out these steps in case s differs between the unfinished and current runs
  AlgoResults$W[,1:s_finish] <- AlgoResults_unfinished$W[,1:s_finish]
  AlgoResults$d[,,1:s_finish] <- AlgoResults_unfinished$d[,,1:s_finish]
  
  AlgoResults$Omega_sample[,1:s_finish] <- AlgoResults_unfinished$Omega_sample[,1:s_finish]
  AlgoResults$Omega_sample_phys[,1:s_finish] <- AlgoResults_unfinished$Omega_sample_phys[,1:s_finish]
  AlgoResults$loss[1:s_finish] = AlgoResults_unfinished$loss[1:s_finish] 
  AlgoResults$u = AlgoResults_unfinished$u
  AlgoResults$u_selected[,1:s_finish] = AlgoResults_unfinished$u_selected[,1:s_finish] 
  AlgoResults$essstore[1:s_finish] <- AlgoResults_unfinished$essstore[1:s_finish]
  AlgoResults$lambda_store[1:s_finish] <- AlgoResults_unfinished$lambda_store[1:s_finish]
  AlgoResults$mu_store[,1:s_finish] <- AlgoResults_unfinished$mu_store[,1:s_finish]
  AlgoResults$Sigma_store[,,1:s_finish] <- AlgoResults_unfinished$Sigma_store[,,1:s_finish]
  AlgoResults$accprob_store[1:s_finish] <- AlgoResults_unfinished$accprob_store[1:s_finish]
  AlgoResults$es_score[1:s_finish] <- AlgoResults_unfinished$es_score[1:s_finish]
  AlgoResults$rh_score[1:s_finish] <- AlgoResults_unfinished$rh_score[1:s_finish]
  AlgoResults$RF_current <- AlgoResults_unfinished$RF_current

  AlgoResults$input_folder <- AlgoResults_unfinished$input_folder
  
  AlgoResults$lambda_start <- AlgoResults_unfinished$lambda_start
  AlgoResults$s_restart <- AlgoResults_unfinished$s_restart
  
  return(list(
    AlgoResults=AlgoResults,
    s_start=s_start
  ))
}

AMCMC <- function(dir, AlgoParams, Model, iVals, unfinished=F, oldtag=NULL, tag_notes=NULL){
  #Input: 
  # - AlgoParams: Parameters describing the AMCMC Algorithm (e.g. the ABC rejection threshold for higher level priors)
  # - Model: Describes the data simulation and calculation of the distance measure 
  # - iVals: The initial values and starting proposal covariance
  # - init_val_phys: the initial proposal (in the physical space)
  # - Unfinished: Continue an old MCMC chain or start a new one?
  # -                If TRUE, then include oldtag - the tag (end of the filename) of an unfinished ABC-SMC run to be completed.
  # - Tag_notes: Additional extention to the results file name (e.g. describing hyperparameters, dataset)
  # Output:
  # - A list containing the parameter values for each particle at each step of the ABC-SMC algorithm, as well as the
  #   corresponding weights, distances and final perturbation kernel
  # Details:
  # - Uses the adaptive algorithm described in 'Accelerating adaptation in the adaptive Metropolis-Hastings random walk algorithm', Spencer 2021
  
  n_x <- length(Model$par_lb) # number of parameters
  
  # results file name:
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "") 
  tag<-ifelse(is.null(tag_notes), tag, paste0(tag, '_', tag_notes))
  folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects/")
  
  # number of training events:
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
  
  #restrict to training data only if subdivided:
  dat='all'
  if (length(grep('^Train/' , ufiles, value = TRUE))!=0){
    dat = 'Train'
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  }
  
  n_events <- length(ufiles)
  
  # parameters of the adaptive MCMC algorithm:
  c = 2.38^2 / n_x
  A = pnorm(AlgoParams$a/2)
  delta = (1-1/n_x)*(sqrt(2*pi)*exp((A^2)/2)/(2*A)) + 1/(n_x*AlgoParams$a*(1-AlgoParams$a))
  
  iter_func = function(s){
    return(floor(s/2))
  }
  
  #set up the results file
  AlgoResults <- list(
    input_folder = AlgoParams$input_folder,
    Omega_sample = array(NA, dim=c(n_x, AlgoParams$N_steps)), #parameters on untransformed space
    Omega_sample_phys = array(NA, dim=c(n_x, AlgoParams$N_steps)), #parameters on physical space
    loss = array(Inf, dim=c( AlgoParams$N_steps)), # evaluations of the loss function at each step
    es_score = array(Inf, dim=c( AlgoParams$N_steps)), # evaluations of the energy score at each step
    rh_score = array(Inf, dim=c( AlgoParams$N_steps)), # evaluations of the rank histogram score at each step
    u = array(NA, dim=c(n_events, AlgoParams$m_CRPS*AlgoParams$Np, 3, 3)), # event-level errors at most recent 3 steps
    u_selected = array(NA, dim=c(3, AlgoParams$N_steps)), # a small sample of the event-level errors across all steps
    lambda_store=array(NA, AlgoParams$N_steps), #parameter of adaptive MCMC algorithm (see Spencer, 2021)
    RF_current = NULL, # correlated local error at current step
    mu_store=array(NA, dim=c(n_x, AlgoParams$N_steps)), #parameter of adaptive MCMC algorithm (see Spencer, 2021)
    Sigma_store=array(NA, dim=c(n_x, n_x, AlgoParams$N_steps)), #parameter of adaptive MCMC algorithm (see Spencer, 2021)
    accprob_store=array(NA, AlgoParams$N_steps), #parameter of adaptive MCMC algorithm (see Spencer, 2021)
    s_restart = round(5/((AlgoParams$a)*(1-AlgoParams$a))), #parameter of adaptive MCMC algorithm (see Spencer, 2021)
    lambda_start = 1 #parameter of adaptive MCMC algorithm (see Spencer, 2021)
  )
  
  
  #-------------------------------------------------------------
  # Store ODD object extents and resolutions to calculate the random fields
  #-------------------------------------------------------------
  x <- file.info(paste0(folderin,ufiles))
  ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
  AlgoResults$RF_current <- list()
  ODD<-readODD(paste0(folderin,ufiles[1]))
  crs_all <- crs(ODD)
  res_all <- res(ODD)
  
  ext_all <- list()
  for (filer in ufiles){
    ODD<-readODD(paste0(folderin,filer))
    event_i = which(ufiles==filer)
    
    ext_all[[event_i]] = ext(ODD)
  }
  
  if(unfinished==F){ 
    #-----------------------------------------------------------
    # STEP 1:
    #-----------------------------------------------------------
    cat('Step: 1\n')
    cat('   Sampling local and event wide error terms (this may take several minutes)\n')
    # Store ODD object extents and resolutions to calculate the random fields
    initialise_random_fields <- function(filer) {
      
      event_i <- which(ufiles == filer)
      
      ext_i <- ext_all[[event_i]]
      r <- rast(crs = crs_all, extent = ext_i, resolution = res_all)
      grid <- as.data.frame(crds(r))
      names(grid) <- c("x", "y")  # Name the columns
      
      vgm_model <- vgm(psill = 1, model = "Mat", range = 0.5, kappa = 1)
      gstat_mod <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 3)
      
      RF_local <- array(
          as.matrix(predict(gstat_mod, newdata = grid, nsim = AlgoParams$Np * AlgoParams$m_CRPS * 3)[, 3:(2 + AlgoParams$Np * AlgoParams$m_CRPS * 3)]),
          dim = c(ncell(r), AlgoParams$Np * AlgoParams$m_CRPS, 3)
      )
      
      return(list(event_i = event_i, RF_local = RF_local))
    }
    results_list <- mclapply(ufiles, initialise_random_fields, mc.cores = AlgoParams$cores, mc.preschedule=F)
    AlgoResults$RF_current <- vector("list", length(ufiles))
    for (event in results_list) {
      AlgoResults$RF_current[[event$event_i]] <- event$RF_local
    }
    
    # Sample from the model and evaluate the loss function
    AlgoResults$Omega_sample[,1] = iVals$x0
    AlgoResults$Omega_sample_phys[,1] = iVals$x0 %>% Array2Physical(Model) %>% unlist()
    AlgoResults$u[,,,1] = rnorm(length(AlgoResults$u[,,,1]))
    proposed =  AlgoResults$Omega_sample_phys[,1] %>% relist(skeleton=Model$skeleton)
    proposed$u = AlgoResults$u[,,,1] 
    proposed$u_local = AlgoResults$RF_current
    
    cat('   Sampling from model given parameterisation and errors\n')
    impact_sample = SampleImpact(dir, Model, proposed %>% addTransfParams(), AlgoParams, dat=dat)
    loss_prop_all = CalcDist(impact_sample, AlgoParams)
    AlgoResults$loss[1] = sum(loss_prop_all)
    AlgoResults$u_selected[,1] = c(AlgoResults$u[1,1,1,1],AlgoResults$u[2,2,2,1],AlgoResults$u[3,3,3,1])
    cat(paste0('   Loss: ', AlgoResults$loss[1],'\n'))
    
    #Store adaptive MCMC parameters
    AlgoResults$lambda_store[1] = AlgoResults$lambda_start
    AlgoResults$mu_store[,1] = AlgoResults$Omega_sample[,1]
    AlgoResults$Sigma_store[,,1] = iVals$COV
    AlgoResults$accprob_store[1] = AlgoParams$alpha_star
    
    #-----------------------------------------------------------
    # STEP 2:
    #-----------------------------------------------------------
    cat('Step: 2\n')
    Omega_prop <- multvarNormProp(xt=AlgoResults$Omega_sample[,1], propPars= c * AlgoResults$lambda_store[1] * AlgoResults$Sigma_store[,,1]) #perturb the proposal
    Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
    
    #----------------------------------
    
    HP<- Model$HighLevelPriors(Omega_prop_phys %>% addTransfParams(), Model)
    if (HP> AlgoParams$ABC & Model$higherpriors){
      #reject parameterisation due to higher level prior
      cat('   Rejected by higher level priors\n')
      AlgoResults$accprob_store[2] = 0
      AlgoResults$Omega_sample[,2] = AlgoResults$Omega_sample[,1]
      AlgoResults$Omega_sample_phys[,2] = AlgoResults$Omega_sample_phys[,1]
      AlgoResults$loss[2] = AlgoResults$loss[1]
      AlgoResults$u[,,,2] = AlgoResults$u[,,,1]
      AlgoResults$u_selected[,2] = AlgoResults$u_selected[,1]
    }  else {
      
      #sample from model using proposal
      proposed = Omega_prop_phys %>% addTransfParams()
      
      cat('   Sampling local and event wide error terms\n')
      
      #perturbate the local error terms:
      perturbate_random_field <- function(event_i) {
        # Create raster grid
        r <- rast(crs = crs_all, extent = ext_all[[event_i]], resolution = res_all)
        grid <- as.data.frame(crds(r))
        names(grid) <- c("x", "y")
        
        # Define variogram model
        vgm_model <- vgm(psill = 1, model = "Mat", range = 0.5, kappa = 1)
        gstat_mod <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 3)
        
        # Generate random field perturbation
        RF_local <- suppressMessages(array(
          as.matrix(predict(gstat_mod, newdata = grid, nsim = AlgoParams$Np * AlgoParams$m_CRPS * 3)[, 3:(2 + AlgoParams$Np * AlgoParams$m_CRPS * 3)]),
          dim = c(ncell(r), AlgoParams$Np * AlgoParams$m_CRPS, 3)
        ))
        
        RF_perturbed <- AlgoParams$rho_local * AlgoResults$RF_current[[event_i]] + sqrt(1 - AlgoParams$rho_local^2) * RF_local
        return(RF_perturbed)
        
      }
      RF_prop <- mclapply(seq_along(ext_all), perturbate_random_field, mc.cores = AlgoParams$cores, mc.preschedule=F)
      proposed$u_local = RF_prop
      
      #perturbate the event-wide error terms:
      epsilon <- rnorm(length(AlgoResults$u[,,,1]))
      u_prop <- AlgoParams$rho * c(AlgoResults$u[,,,1]) + sqrt(1-AlgoParams$rho^2) * epsilon 
      proposed$u = array(u_prop, dim=c(n_events, AlgoParams$m_CRPS*AlgoParams$Np, 3))
      
      #sample from the model and calculate the loss function:
      cat('   Sampling from model given parameterisation and errors\n')
      impact_sample <- SampleImpact(dir = dir,Model = Model,
                                    proposed = proposed, 
                                    AlgoParams = AlgoParams,
                                    dat=dat)
      loss_prop_all <- CalcDist(impact_sample, AlgoParams)#, corr_noise = pnorm(u_rh_prop), corr_noise2 = pnorm(u_rh2_prop))
      loss_prop <- sum(loss_prop_all)
      cat(paste0('   Loss: ', loss_prop, '\n'))
      
      #calculate the acceptance probability:
      #print(modifyAcc(Omega_prop, Omega_sample_s[n,], Model))
      min_loss <- min(AlgoParams$learning_rate * loss_prop, AlgoParams$learning_rate * AlgoResults$loss[1])
      AlgoResults$accprob_store[2] <- min(1, exp(-AlgoParams$learning_rate * loss_prop + min_loss)/exp(-AlgoParams$learning_rate * AlgoResults$loss[1] + min_loss) * modifyAcc(Omega_prop, AlgoResults$Omega_sample[,1], Model, AlgoResults$lambda_store[1] * AlgoResults$Sigma_store[,,1]))
      
      u <- runif(1)
      if(u < AlgoResults$accprob_store[2]){
        #ACCEPT
        cat('   Proposal accepted!\n')
        AlgoResults$Omega_sample[,2] = Omega_prop
        AlgoResults$Omega_sample_phys[,2] = unlist(Omega_prop_phys)
        AlgoResults$loss[2] = loss_prop
        AlgoResults$u[,,,2] = proposed$u
        AlgoResults$u_selected[,2] = c(proposed$u[1,1,1],proposed$u[2,2,2],proposed$u[3,3,3])
        AlgoResults$RF_current = RF_prop
        #AlgoResults$u_rh = u_rh_prop
        #AlgoResults$u_rh2 = u_rh2_prop
      }  else {
        cat('   Proposal rejected\n')
        AlgoResults$Omega_sample[,2] = AlgoResults$Omega_sample[,1]
        AlgoResults$Omega_sample_phys[,2] = AlgoResults$Omega_sample_phys[,1]
        AlgoResults$loss[2] = AlgoResults$loss[1]
        AlgoResults$u[,,,2] = AlgoResults$u[,,,1]
        AlgoResults$u_selected[,2] = AlgoResults$u_selected[,1]
      }
    }
    
    #Store parameters for adaptive MCMC method:
    AlgoResults$mu_store[,2] = (AlgoResults$Omega_sample[,2] +  AlgoResults$Omega_sample[,1])/2
    AlgoResults$Sigma_store[,,2] = 1 / (AlgoParams$v_0 + n_x + 3) * 
      ((AlgoResults$mu_store[,2] %*% t(AlgoResults$mu_store[,2]) + AlgoResults$mu_store[,1] %*% t(AlgoResults$mu_store[,1])) - 
         2*AlgoResults$mu_store[,2] %*% t(AlgoResults$mu_store[,2]) +
         (AlgoParams$v_0 + n_x + 1)* AlgoResults$Sigma_store[,,1] )
    AlgoResults$lambda_store[2] = max(AlgoParams$lambda_min, AlgoResults$lambda_store[1] * exp(delta / (2) * (AlgoResults$accprob_store[2]-AlgoParams$a)))
    
    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag))  
    
    s_start = 3
    
  } else { 
    #Collect relevant information from the unfinished sample
    UnfinishedAlgoResults <- retrieve_UnfinishedAlgoResults_AMCMC(dir, oldtag, AlgoResults)
    AlgoResults <- UnfinishedAlgoResults$AlgoResults
    s_start <- UnfinishedAlgoResults$s_start
    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag)) 
  }
  
  for (s in s_start:AlgoParams$N_steps){
    if (s %% 20 == 0) { 
      cat('Saving results file\n')
      saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag))
    } else if(s %% 10 == 0){
      cat('Saving results file\n')
      saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag, "_backup"))
    }
    cat(paste0('Step: ', s,'\n'))
    
    #generate proposal
    Omega_prop <- multvarNormProp(xt=AlgoResults$Omega_sample[,s-1], propPars= c * AlgoResults$lambda_store[s-1] * AlgoResults$Sigma_store[,,s-1]) #perturb the proposal
    Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
    
    HP<- Model$HighLevelPriors(Omega_prop_phys %>% addTransfParams(), Model)
    if (HP> AlgoParams$ABC & Model$higherpriors){
      cat('   Rejected by higher level priors\n')
      #reject due to higher level priors
      AlgoResults$accprob_store[s] = 0
      AlgoResults$Omega_sample[,s] = AlgoResults$Omega_sample[,s-1]
      AlgoResults$Omega_sample_phys[,s] = AlgoResults$Omega_sample_phys[,s-1]
      AlgoResults$loss[s] = AlgoResults$loss[s-1]
      AlgoResults$u[,,,ifelse(s %% 3==0, 3, s%%3)] = AlgoResults$u[,,,ifelse((s-1) %% 3==0, 3, (s-1)%%3)]
      AlgoResults$u_selected[,s] = AlgoResults$u_selected[,s-1]
    } else {
      proposed = Omega_prop_phys %>% addTransfParams()
      
      #perturb local errors:
      perturbate_random_field <- function(event_i) {
        # Create raster grid
        r <- rast(crs = crs_all, extent = ext_all[[event_i]], resolution = res_all)
        grid <- as.data.frame(crds(r))
        names(grid) <- c("x", "y")
        
        # Define variogram model
        vgm_model <- vgm(psill = 1, model = "Mat", range = 0.5, kappa = 1)
        gstat_mod <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 3)
        
        # Generate random field perturbation
        RF_local <- suppressMessages(array(
          as.matrix(predict(gstat_mod, newdata = grid, nsim = AlgoParams$Np * AlgoParams$m_CRPS * 3)[, 3:(2 + AlgoParams$Np * AlgoParams$m_CRPS * 3)]),
          dim = c(ncell(r), AlgoParams$Np * AlgoParams$m_CRPS, 3)
        ))
        
        RF_perturbed <- AlgoParams$rho_local * AlgoResults$RF_current[[event_i]] + sqrt(1 - AlgoParams$rho_local^2) * RF_local
        return(RF_perturbed)
      }
      RF_prop <- mclapply(seq_along(ext_all), perturbate_random_field, mc.cores = AlgoParams$cores, mc.preschedule=F)
      proposed$u_local = RF_prop
      
      #perturb event-level error:
      cat('   Sampling local and event wide error terms\n')
      epsilon <- rnorm(length(AlgoResults$u[,,,ifelse((s-1) %% 3==0, 3, (s-1)%%3)]))
      u_prop <- AlgoParams$rho * c(AlgoResults$u[,,,ifelse((s-1) %% 3==0, 3, (s-1)%%3)]) + sqrt(1-AlgoParams$rho^2) * epsilon
      proposed$u = array(u_prop, dim=c(n_events, AlgoParams$m_CRPS*AlgoParams$Np, 3))
      
      #sample from model and calculate loss function
      cat('   Sampling from model given parameterisation and errors\n')
      impact_sample <- SampleImpact(dir = dir,Model = Model,
                                    proposed = proposed, 
                                    AlgoParams = AlgoParams,
                                    dat=dat)
      loss_all_prop <- CalcDist(impact_sample, AlgoParams) #, corr_noise = pnorm(u_rh_prop), corr_noise2 = pnorm(u_rh2_prop))
      loss_prop <- sum(loss_all_prop)
      cat(paste0('   Loss: ', loss_prop, '\n'))
      
      #calculate the acceptance probability:
      #print(modifyAcc(Omega_prop, Omega_sample_s[n,], Model))
      min_loss <- min(AlgoParams$learning_rate * loss_prop, AlgoParams$learning_rate * AlgoResults$loss[s-1])
      AlgoResults$accprob_store[s] <- min(1, exp(-AlgoParams$learning_rate * loss_prop + min_loss)/exp(-AlgoParams$learning_rate * AlgoResults$loss[s-1] + min_loss) * modifyAcc(Omega_prop, AlgoResults$Omega_sample[,s-1], Model, AlgoResults$lambda_store[s] * AlgoResults$Sigma_store[,,s]))
      
      u <- runif(1)
      if(u < AlgoResults$accprob_store[s]){
        #ACCEPT:
        cat('   Proposal accepted!\n')
        AlgoResults$Omega_sample[,s] = Omega_prop
        AlgoResults$Omega_sample_phys[,s] = unlist(Omega_prop_phys)
        AlgoResults$loss[s] = loss_prop
        AlgoResults$es_score[s] = loss_all_prop[1]
        AlgoResults$rh_score[s] = loss_all_prop[2]
        AlgoResults$u[,,,ifelse(s %% 3==0, 3, s%%3)] = proposed$u
        AlgoResults$u_selected[,s] =  c(proposed$u[1,1,1],proposed$u[2,2,2],proposed$u[3,3,3])
        AlgoResults$RF_current = RF_prop
        #AlgoResults$u_rh = u_rh_prop 
        #AlgoResults$u_rh2 = u_rh2_prop 
      }  else {
        cat('   Proposal rejected\n')
        AlgoResults$Omega_sample[,s] = AlgoResults$Omega_sample[,s-1]
        AlgoResults$Omega_sample_phys[,s] = AlgoResults$Omega_sample_phys[,s-1]
        AlgoResults$es_score[s] = AlgoResults$es_score[s-1]
        AlgoResults$rh_score[s] = AlgoResults$rh_score[s-1]
        AlgoResults$loss[s] = AlgoResults$loss[s-1]
        AlgoResults$u[,,,ifelse(s %% 3==0, 3, s%%3)] = AlgoResults$u[,,,ifelse((s-1) %% 3==0, 3, (s-1)%%3)]
        AlgoResults$u_selected[,s] = AlgoResults$u_selected[,s-1]
      }
    }
    
    #Update parameters of the adaptive MCMC algorithm:
    if (iter_func(s)==iter_func(s-1)){
      AlgoResults$mu_store[,s] = (s - iter_func(s))/(s-iter_func(s)+1) *AlgoResults$mu_store[,s-1] + 1/(s-iter_func(s)+1) * AlgoResults$Omega_sample[,s]
      AlgoResults$Sigma_store[,,s] = 1/(s-iter_func(s)+AlgoParams$v_0 + n_x + 2) * (
        (s - iter_func(s) + AlgoParams$v_0+n_x+1)*AlgoResults$Sigma_store[,,s-1] + 
          AlgoResults$Omega_sample[,s] %*% t(AlgoResults$Omega_sample[,s]) +
          (s-iter_func(s)) * AlgoResults$mu_store[,s-1] %*% t(AlgoResults$mu_store[,s-1]) -
          (s-iter_func(s)+1) * AlgoResults$mu_store[,s] %*% t(AlgoResults$mu_store[,s]))
    } else {
      AlgoResults$mu_store[,s] = AlgoResults$mu_store[,s-1] + 1/(s-iter_func(s)+1) * (AlgoResults$Omega_sample[,s]-AlgoResults$Omega_sample[,iter_func(s)-1])
      AlgoResults$Sigma_store[,,s] = AlgoResults$Sigma_store[,,s-1] + 1/(s-iter_func(s) + AlgoParams$v_0 + n_x+2) * (
        AlgoResults$Omega_sample[,s] %*% t(AlgoResults$Omega_sample[,s]) - 
          AlgoResults$Omega_sample[,iter_func(s)-1] %*% t(AlgoResults$Omega_sample[,iter_func(s)-1]) +
          (s-iter_func(s)+1) * (AlgoResults$mu_store[,s-1] %*% t(AlgoResults$mu_store[,s-1]) -  AlgoResults$mu_store[,s] %*% t(AlgoResults$mu_store[,s]))
      )
    }
    
    AlgoResults$lambda_store[s] = max(AlgoParams$lambda_min, AlgoResults$lambda_store[s-1] * exp(delta / (AlgoResults$s_restart + s) * (AlgoResults$accprob_store[s]-AlgoParams$a)))
    if (abs(log(AlgoResults$lambda_store[s]) - log(AlgoResults$lambda_start))>log(3)){
      AlgoResults$lambda_start = AlgoResults$lambda_store[s]
      AlgoResults$s_restart <- 5/((AlgoParams$a)*(1-AlgoParams$a)) - s
    }
    
  }
  saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag))
  cat('MCMC chain completed!\n')
}


#----------------------------------------------------------------------------------------------------------------------------
#--------------------------------Optimisation Approaches (not up to date)----------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

## Generate the Adaptive Metropolis Global Scaling Factor discount factor
# GenerateGamzy<-function(AlgoParams){
#   AlgoParams$gamzy0/(1:(AlgoParams$itermax+AlgoParams$cores))^
#     (seq(from=1/(1+AlgoParams$epsilon),to=1,length.out=(AlgoParams$itermax+AlgoParams$cores)))
# }
# 
# # Check that the proposed initial values, model and methodology parameters
# checkLTargs<-function(Model,iVals,AlgoParams){
#   if(length(unlist(Model$links))!=length(unlist(iVals$x0))) stop("Mismatching link functions and initial values of Omega space")
#   if(length(unlist(iVals$x0))!=nrow(iVals$COV) | length(unlist(iVals$x0))!=ncol(iVals$COV)) stop("Mismatching initial values and initial covariance matrix")
# }
# 
# NelderMeadOptim<-function(dir,Model,iVals,AlgoParams){
#   
#   # We don't need an initial guess of the covariance matrix for MLE optimisation
#   x0=Physical2Proposed(iVals$x0,Model)%>%unlist()
#   # Only optimise over the input iVals that are not NAs
#   if(is.null(AlgoParams$indices)) AlgoParams$indices<-length(x0)
#   # Cost function (note that this still includes the priors and ABC rejection, so isn't purely frequentist MLE)
#   Fopty<-function(vals){
#     #x0[!AlgoParams$indices]<-vals
#     x0 <- vals
#     # Convert proposal in to physical values ODDRIN understands
#     x0%<>%Proposed2Physical(Model)
#     # Trust me, it's nice to know the parameters tested, this optim algorithm can struggle with stochastic target distributions
#     print(unname(unlist(x0)))
#     # Posterior calculation
#     posterior<-logTarget(dir,Model,x0,AlgoParams,expLL = F)
#     print(posterior)
#     print("...")
#     return(posterior)
#   }
#   # Optimisation algorithm - Nelder & Mead, 1965, and includes outputting the Hessian
#   output<-optim(par=x0,
#                 fn = Fopty,control = list(maxit = AlgoParams$itermax,
#                                           fnscale=-1,
#                                           reltol=1.5e-3),
#                 hessian = T)
#   
#   
#   x0[!AlgoParams$indices]<-output$par
#   # Convert the optimised values into physical ones
#   x0%<>%Proposed2Physical(Model)
#   
#   return(list(PhysicalValues=x0,OptimisationOut=output,optimisedParams=(1:length(unlist(x0)))[!AlgoParams$indices]))
#   
# }
# 
# 
# # Using the bisection method to generate a new guess
# OptimMd<-function(FFF,x,LLs){
#   
#   Ldiff<-xdiff<-1
#   while(Ldiff>1e-4 & xdiff>1e-6){
#     xp<-c(mean(x[c(1,2)]),mean(x[c(2,3)]))
#     LLxp<-c(FFF(xp[1]),FFF(xp[2]))
#     LLt<-c(LLs[1],LLxp[1],LLs[2],LLxp[2],LLs[3]);  xt<-c(x[1],xp[1],x[2],xp[2],x[3])
#     Ldiff<-min(abs(diff(LLt))); xdiff<-min(abs(diff(xt)))
#     ix<-which.max(LLt)
#     if(any(is.na(x))) {print(c(x,LLs));print(c(xt,LLt))}
#     x<-xt[c(ix-1,ix,ix+1)];  LLs<-LLt[c(ix-1,ix,ix+1)]
#     
#   }
#   
#   ix<-which.max(LLs)
#   print(paste0("Modifier = ",x[ix],",    LLmd = ",LLs[ix]))
#   
#   return(x[ix])
#   
# }
# 
# GetModifier<-function(Fopty,init,lower,upper){
#   
#   skely<-rep(0,length(init))
#   
#   for(i in 1:length(init)){
#     
#     Ftmp<-function(x) {
#       y<-skely; y[i]<-x
#       as.numeric((Fopty(y))[i])
#     }
#     
#     if(init[i]>0){ # Calculate values between 0 and upper
#       xx<-(0:20)*upper/20
#       LLMd<-vapply(xx,Ftmp,numeric(1))
#       up<-upper
#       while(which.max(LLMd)==length(xx)){
#         up<-up*1.5
#         xx<-0.6*up+0.4*up*1:10/10
#         LLMd<-vapply(xx,Ftmp,numeric(1))
#       }
#       srtr<-sort(LLMd,decreasing = T,index.return = T); srtr$x<-srtr$x[1:3];  srtr$ix<-srtr$ix[1:3]
#       srtr2<-sort(xx[srtr$ix[1:3]],index.return = T);  xxx<-srtr2$x;  LLs<-srtr$x[srtr2$ix]
#       skely[i] <-OptimMd(Ftmp,xxx,LLs)
#     } else { # Calculate values between lower and 0
#       xx<-(0:20)*lower/20
#       LLMd<-vapply(xx,Ftmp,numeric(1))
#       lo<-lower
#       while(which.max(LLMd)==length(xx)){
#         lo<-lo*1.5
#         xx<-0.6*lo+0.4*lo*1:10/10
#         LLMd<-vapply(xx,Ftmp,numeric(1))
#       }
#       srtr<-sort(LLMd,decreasing = T,index.return = T); srtr$x<-srtr$x[1:3];  srtr$ix<-srtr$ix[1:3]
#       srtr2<-sort(xx[srtr$ix[1:3]],index.return = T);  xxx<-srtr2$x;  LLs<-srtr$x[srtr2$ix]
#       skely[i] <-OptimMd(Ftmp,xxx,LLs)
#     }
#     
#   }
#   
#   return(as.list(skely))
#   
# }
# 
# SingleEventsModifierCalc<-function(dir,Model,Omega,AlgoParams){
#   
#   warning("You need to set the subnational linear predictor variables to their average value before applying the modifiers")
#   
#   modifiers<-data.frame()
#   ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/ODDobjects"),pattern=Model$haz,recursive = T,ignore.case = T)
#   
#   if(Model$haz=="EQ"){
#     lower<--0.5
#     upper<-3
#   } else stop("Insert lower and upper bounds for the hazard studied, only have EQ so far")
#   
#   for(i in 1:length(ufiles)){
#     
#     print(ufiles[i])
#     
#     ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/",ufiles[i]))
#     ODDy@gmax%<>%arrange(desc(gmax))
#     
#     nisos<-length(unique(ODDy@gmax$iso3))
#     linp<-rep(list(0.),nisos)
#     names(linp)<-unique(ODDy@gmax$iso3)
#     ODDy@modifier<-linp
#       
#     tLL<-tryCatch(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams),
#                   error=function(e) NA)
#     
#     if(is.na(tLL)) stop(ufiles[i])
#     
#     init<-log(tLL@predictDisp$gmax/tLL@predictDisp$disp_predictor)*0.1
#     init[init<lower]<-0.5*lower; init[init>upper]<-0.5*upper; 
#     init%<>%as.list()
#     names(init)<-as.character(tLL@predictDisp$iso3)
#   
#     Fopty<-function(Md){
#       
#       Md%<>%as.list()
#       names(Md)<-as.character(ODDy@gmax$iso3)
#       ODDy@modifier<-Md
#       
#       lTarg<-logTargetSingle(ODDy,Model,Omega,AlgoParams)+
#         Model$HighLevelPriors(Omega,Model,modifier = list(Md))
#       
#       print(paste0(ufiles[i],"   :   lTarg={",strcat(as.character(lTarg),collapse = ","),"},    params= {",strcat(as.character(Md),collapse = ","),"}"))
#       return(lTarg)
#     }
#     
#     fin<-tryCatch(GetModifier(Fopty,init,lower,upper),error=function(e) NA)
#     if(is.na(fin)) {
#       print(paste0("Failed to calculate modifier for ",ufiles[i]))
#       print("...")
#       next
#     }
#     print(paste0("Final params= {",strcat(as.character(unlist(fin)),collapse = ","),"}"))
#     print("...")
#     names(fin)<-as.character(ODDy@gmax$iso3)
#     ODDy@modifier<-fin
#     
#     ODDy<-tryCatch(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams),
#                   error=function(e) NA)
#     if(is.na(ODDy)) stop(ufiles[i])
#     
#     modifiers%<>%rbind(data.frame(iso3=ODDy@predictDisp$iso3,modifier=as.numeric(unlist(fin)),
#                                   event=ufiles[i],gmax=ODDy@gmax$gmax,
#                                   eventid=as.numeric(strsplit(as.character(ufiles[i]),split = "_")[[1]][2]),
#                                   predictor=ODDy@predictDisp$disp_predictor))
#     
#     saveRDS(ODDy, paste0(dir,"IIDIPUS_Results/ODDobjects/",ufiles[i]))
#     saveRDS(modifiers,paste0(dir,"IIDIPUS_Results/ODDobjects/modifiers.Rdata"))
# 
#   }
# 
#   return(modifiers)
#   
# }

# CorrelateModifier<-function(modifiers,Model){
#   
#   # This function uses all the ODDobjects available to extract iso3, start date and event_id values for each hazard
#   DispData<-ExtractDispData_ODD(Model$haz)
#   # Now extract the vulnerability variables to be used in the correlation - these are defined in Model
#   val<-vulnerabilityVars(DispData,Model)
#   # Use only the variables we think are relevant
#   warning("Using only Hamish's pre-selected World Bank indicators for correlation")
#   indies<-read_csv("./IIDIPUS_Input/REDUCED_WB-WorldDevelopment_Indicators.csv")
#   val%<>%filter(indicator%in%indies$indicator_id)
#   # Merge the two data frames
#   modifiers$indicator<-"modifier"
#   modifiers<-val%>%group_by(eventid)%>%summarise(date=unique(date))%>%merge(modifiers,by="eventid")
#   # Housekeeping
#   modifiers%<>%transmute(eventid=eventid,nearval=modifier,indicator=indicator)%>%
#     dplyr::select(eventid,nearval); colnames(modifiers)[3]<-"modifier"
#   
#   # Get data frames in the correct structure
#   for(inds in unique(val$indicator[val$indicator!="modifier"])){
#     tmp<-filter(val,indicator==inds)
#     tn<-data.frame(value=tmp$nearval); colnames(tn)<-inds
#     modifiers%<>%cbind(tn)
#   }
#   rm(tn,ti,tmp)
#   
#   # Scale and center the vulnerability variables
#   modifiers[3:ncol(modifiers)]<-base::scale(modifiers[3:ncol(modifiers)])
#   # Remove all NA values
#   redmodies<-redmodies[!apply(redmodies,1,function(x) any(is.na(x))),]
#   # Remove all vulnerability variables that correlate strongly with one another - Variance Inflation Factor
#   ncor <- cor(redmodies[,4:ncol(redmodies)])
#   redmodies<-redmodies[,c(1,3,which(!1:ncol(redmodies)%in%caret::findCorrelation(ncor, cutoff=0.75)))]
#   redmodies<-redmodies[,c(1,3,2,4:ncol(redmodies))]
#   redmodies%<>%dplyr::select(-c(eventid))
#   
#   # Here we go!
#   vulncor<-LMFeatureSelection(final,nlim = 3,weights = weights)
#   # vulncor<-LMFeatureSelection(final,nlim = 3,weights = weights,fn = ":",intercept = T)
#   # vulncor<-LMFeatureSelection(final,nlim = 4,weights = weights,fn = ":",intercept = T)
#   
#   return(vulncor)
#   
# }


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


