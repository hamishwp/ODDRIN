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
                 minVar=1e-4 # Prevent certain parameters from being too sure of themselves
                 )
		 
if(is.null(AlgoParams$AllParallel)){
  if(AlgoParams$cores>4) { AlgoParams$AllParallel<-T
  } else AlgoParams$AllParallel<-F
}
# Choose the parameterisation algorithm - the string must match the function name exactly
Algorithm<-"AMCMC2" # "NelderMeadOptim", "AMCMC"

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
    proposed[i] <- match.fun(Model$links[[names(proposed)[i]]])(proposed[i])
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
    proposed[i] <- match.fun(Model$unlinks[[names(proposed)[i]]])(proposed[i])
  }
  # Reshape into desired structure
  proposed%>%relist(skeleton=Model$skeleton)
  
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

AMCMC <-function(dir,Model,iVals,AlgoParams){
  
  # Set Random Number Generator (RNG) initial seed
  set.seed(round(runif(1,0,100000))) 
  # Check no mistakes have been made in the model, methodology and initial values
  checkLTargs(Model,iVals,AlgoParams)
  
  t_0 <- 500
  
  xPrev<-propMu<-unlist(iVals$x0)
  n_x <- length(xPrev)
  xbar_tminus1 <- xPrev
  
  C_0 = diag(0.0001, nrow=n_x)
  
  s_d = (2.38)^2/n_x
  eps = diag(0.000000001, nrow=n_x)
  propCOV <- diag(n_x)
  epsilon=1
  
  output <- matrix(NA, nrow=AlgoParams$itermax, ncol=n_x+1)
  xNew <- rep(NA,n_x)
  lTargNew<-alpha<-c()
  # Create file names with unique names for storage reasons
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  
  # Find first log-target value using initial values
  output[1, ] <- c(logTarget(dir = dir,Model = Model,
                             proposed = xPrev %>%Proposed2Physical(Model),AlgoParams = AlgoParams, epsilon=epsilon),
                   xPrev)
  lTargOld<-output[1,1]
  
  # Start the iterations!
  it <- 2
  while (it <= AlgoParams$itermax){
    print(it)
    t <- it - 1
    epsilon = ifelse(t < t_0, 1 - t * (1-0.1)/t_0, 0.1)
    print(paste("epsilon", epsilon))
    # Parameter proposal
    xNew <- xPrev
    if (t > t_0){
      xNew <- multvarNormProp(xt=xPrev, propPars=propCOV)
    } else {
      xNew <- multvarNormProp(xt=xPrev, propPars=C_0)
    }
    
    # Check proposal is within the parameter space:
    if(any(xNew < Model$par_lb) | any(xNew > Model$par_ub) ){
      output[it,] <- c(lTargOld, xPrev)
      propCOV <- (t-1)/t * propCOV + s_d/(t+1) * (xPrev - xbar_tminus1) %*% t(xPrev - xbar_tminus1) + s_d /t * eps
      xbar_tminus1 <- (t * xbar_tminus1 + xPrev)/(t+1)
      it <- it + 1 
      next
    }
    
    # Convert parameters to physical/useable values
    if(!is.null(Model$links)) xProp<-xNew%>%Proposed2Physical(Model)
    
    # Calculate log-target value
    lTargNew <- tryCatch(logTarget(dir = dir,Model = Model,proposed = xProp,
                                   AlgoParams = AlgoParams, epsilon= epsilon), error=function(e) NA)
    
    # Check if we have a NaN
    if(is.na(lTargNew)|is.infinite(lTargNew)) {
      output[it,] <- c(lTargOld, xPrev)
      propCOV <- (t-1)/t * propCOV + s_d/(t+1) * (xPrev - xbar_tminus1) %*% t(xPrev - xbar_tminus1) + s_d /t * eps
      xbar_tminus1 <- (t * xbar_tminus1 + xPrev)/(t+1)
      it <- it + 1 
      next
    }
    
    #if proposed log likelihood is close to the old log likelihood, then resample the old log likelihood
    if (lTargNew > (lTargOld - 50)){
      lTargOld <- tryCatch(logTarget(dir = dir,Model = Model,proposed = xPrev %>%Proposed2Physical(Model),
                                     AlgoParams = AlgoParams, epsilon=epsilon), error=function(e) NA)
    }
    
    # Prepare for acceptance
    u <- runif(1)
    
    # Acceptance probability
    alpha <- min(c(exp(lTargNew - lTargOld),1))
    
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
    
    output[it,] <- c(lTargOld, xPrev)
    
    Intensity <- seq(0,10,0.1)
    Dfun<-function(I_ij, theta) h_0(I = I_ij,I0 = 4.5,theta = Omega$theta) 
    
    Omega_curr <- xPrev %>% 
      relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model)
    
    par(mfrow=c(4,4))
    
    # Plot S-curves for the actual and MAP parameterisation
    plot_S_curves(Omega, Omega_curr)
    
    #create trace plots
    for(i in 1:14){
      ylim=c(min(unlist(Omega)[i], output[1:it,i+1]), max(unlist(Omega)[i], output[1:it,i+1]))
      plot(output[1:it,i+1], type='l', ylab='', ylim=ylim)
      abline(h=unlist(Omega)[i], col='red')
    }
    
    
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

#Allows you to continue an unfinished MCMC chain that was interrupted
AMCMC_continueunfinished <-function(dir,Model,iVals,AlgoParams, filetag){
  
  # Set Random Number Generator (RNG) initial seed
  set.seed(round(runif(1,0,100000))) 
  # Check no mistakes have been made in the model, methodology and initial values
  checkLTargs(Model,iVals,AlgoParams)
  
  output <- readRDS(paste0(dir, '/IIDIPUS_Results/output_', filetag)) 
  xbar_tminus1 <- readRDS(paste0(dir, '/IIDIPUS_Results/xbar_tminus1_', filetag)) 
  propCOV <- readRDS(paste0(dir, '/IIDIPUS_Results/propCOV_', filetag)) 
  xPrev <- unlist(iVals$x0)
  interruption_point <- min(which(is.na(output[,1])))-1
  xPrev[1:14] <- output[interruption_point, 2:NCOL(output)]
  n_x <- length(xPrev)
  
  C_0 = diag(0.0001, nrow=n_x)
  
  s_d = (2.38)^2/n_x
  eps = diag(0.0000001, nrow=n_x)

  xNew <- rep(NA,n_x)
  lTargNew<-alpha<-c()
  # Create file names with unique names for storage reasons
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  
  lTargOld<-output[interruption_point, 1]
  
  # Start the iterations!
  it <- min(which(is.na(output[,1])))
  while (it <= AlgoParams$itermax){
    print(it)
    t <- it - 1
    
    # Parameter proposal
    xNew <- xPrev
    if (t > t_0){
      xNew <- multvarNormProp(xt=xPrev, propPars=propCOV)
    } else {
      xNew <- multvarNormProp(xt=xPrev, propPars=C_0)
    }
    #xNew[7:14] <- unlist(Omega)[7:14]
    
    # Check proposal is within the parameter space:
    if(any(xNew < Model$par_lb) | any(xNew > Model$par_ub) ){
      output[it,] <- c(lTargOld, xPrev)
      propCOV <- (t-1)/t * propCOV + s_d/(t+1) * (xPrev - xbar_tminus1) %*% t(xPrev - xbar_tminus1) + s_d /t * eps
      xbar_tminus1 <- (t * xbar_tminus1 + xPrev)/(t+1)
      it <- it + 1 
      next
    }
    
    # Convert parameters to physical/useable values
    if(!is.null(Model$links)) xProp<-xNew%>%Proposed2Physical(Model)
    
    # Calculate log-target value
    lTargNew <- tryCatch(logTarget(dir = dir,Model = Model,proposed = xProp,
                                   AlgoParams = AlgoParams), error=function(e) NA)
    
    # Check if we have a NaN
    if(is.na(lTargNew)|is.infinite(lTargNew)) {
      output[it,] <- c(lTargOld, xPrev)
      propCOV <- (t-1)/t * propCOV + s_d/(t+1) * (xPrev - xbar_tminus1) %*% t(xPrev - xbar_tminus1) + s_d /t * eps
      xbar_tminus1 <- (t * xbar_tminus1 + xPrev)/(t+1)
      it <- it + 1 
      next
    }
    
    if (lTargNew > (lTargOld - 50)){
      lTargOld <- tryCatch(logTarget(dir = dir,Model = Model,proposed = xPrev %>%Proposed2Physical(Model),
                                     AlgoParams = AlgoParams), error=function(e) NA)
    }
    
    # Prepare for acceptance
    u <- runif(1)
    
    # Acceptance probability
    alpha <- min(c(exp(lTargNew - lTargOld),1))
    
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
    
    output[it,] <- c(lTargOld, xPrev)
    
    Intensity <- seq(0,10,0.1)
    Dfun<-function(I_ij, theta) h_0(I = I_ij,I0 = 4.5,theta = Omega$theta) 
    
    Omega_curr <- xPrev %>% 
      relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model)
    
    par(mfrow=c(4,4))
    
    # Plot S-curves for the actual and MAP parameterisation
    plot_S_curves(Omega, Omega_curr)
    
    #create trace plots
    for(i in 1:14){
      ylim=c(min(unlist(Omega)[i], output[1:it,i+1]), max(unlist(Omega)[i], output[1:it,i+1]))
      plot(output[1:it,i+1], type='l', ylab='', ylim=ylim)
      abline(h=unlist(Omega)[i], col='red')
    }
    
    
    # Save log-target and parameters
    saveRDS(output,paste0(dir,"IIDIPUS_Results/output_",tag))
    # Save covariance matrix
    saveRDS(propCOV,paste0(dir,"IIDIPUS_Results/covariance_",tag))
    
    saveRDS(xbar_tminus1,paste0(dir,"IIDIPUS_Results/xbar_tminus1_",tag))
    
    
    it <- it + 1
  }
  
  return(list(PhysicalValues=output[which.max(output[,1]),2:ncol(output)] %>% 
                relist(skeleton=Model$skeleton) %>% unlist() %>% Proposed2Physical(Model), # MAP value 
              OptimisationOut=output))
  
}

Algorithm <- match.fun('AMCMC')


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






