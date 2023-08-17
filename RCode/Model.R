##########################################################################
##########################################################################
############################ Model formulation ###########################
### There are three sections to this file:                             ###
### 1) Define model variables, parameterisation and link functions     ###
###    Choose your important variables here, e.g. unemployment rate    ###
### 2) Linear predictor: damage*exp(linearpredictor)                   ###
###    This acts to modify the damage based on the country/area values ###
###    For example, local GDP or pop density decrease expected damage  ###
### 3) Log-likelihood, prior and posterior distribution definitions    ###
###    When optimising the parameters of the model based on historical ###
###    data, this section defines what function to minimise (cost fn)  ###
##########################################################################
##########################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Model variables and format, based on the specific hazard
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

Model<-list()

haz<-"EQ"
# haz<-"TC"
Model$haz<-haz
Model$I0 <- 4.3

WID_perc<- c("p0p10", # Bottom 10% share of Income Distribution
              "p10p20", # Income share held by 10th - 20th percentiles
              "p20p30", # Income share held by 20th - 30th percentiles
              "p30p40", # Income share held by 30th - 40th percentiles
              "p40p50", # Income share held by 40th - 50th percentiles
              "p50p60", # Income share held by 50th - 60th percentiles
              "p60p70", # Income share held by 60th - 70th percentiles
              "p70p80", # Income share held by 70th - 80th percentiles
              "p80p90", # Income share held by 80th - 90th percentiles
              "p90p100") # top 10% share of Income Distribution

Model%<>%c(list(WID_perc=WID_perc))

if(haz=='EQ'){
  Model$vuln_terms <- c('PDens', 'AveSchYrs','LifeExp', 'GNIc', 'Vs30', 'EQFreq')
}

ab_bounded <- function(x, a, b){
  return((b*exp(x)+a)/(exp(x)+1))
}

ab_bounded_inv <- function(x, a, b){
  return(log(x-a)-log(b-x))
}

ab_bounded_acc <- function(xnew, xold, a, b){
  return((b-xnew)/(b-xold)*(xnew-a)/(xold-a))
}

returnX <- function(x,a,b){
  return(x)
}

# Skeleton
Model$skeleton <- list(
  Lambda1=list(mu=NA,sigma=NA#,alpha=NA
               ), 
  Lambda2=list(mu=NA,sigma=NA),
  Lambda3=list(mu=NA,sigma=NA),
  Lambda4=list(mu=NA,sigma=NA),
  eps=list(local=NA,hazard=NA),
  vuln_coeff=list(PDens=NA, AveSchYrs=NA, LifeExp=NA, GNIc=NA, Vs30=NA, EQFreq=NA, Mag=NA),
  check=list(check=NA)
)

Model$Priors <- list( #All uniform so currently not included in the acceptance probability. 
  Lambda1=list(mu=list(dist='unif', min=6.5, max=10.5), 
               sigma=list(dist='unif', min=0.25, max=2) #, alpha=list(dist='unif', min=-0.1, max=0.5)
               ), 
  Lambda2=list(mu=list(dist='unif', min=9, max=12.5), 
               sigma=list(dist='unif', min=0.25, max=2)),
  Lambda3=list(mu=list(dist='unif', min=6.5, max=10), 
               sigma=list(dist='unif', min=0.25, max=2)),
  Lambda4=list(mu=list(dist='unif', min=8, max=12.5), 
               sigma=list(dist='unif', min=0.25, max=2.5)),
  eps=list(local=list(dist='unif', min=0, max=1),
           hazard=list(dist='unif', min=0, max=1)),
  vuln_coeff=list(PDens=list(dist='unif', min=-0.15, max=0.15),
                  EQFreq=list(dist='unif', min=-0.15, max=0.15),
                  AveSchYrs=list(dist='unif', min=-0.15, max=0.15),
                  LifeExp=list(dist='unif', min=-0.15, max=0.15),
                  GNIc=list(dist='unif', min=-0.15, max=0.15),
                  Vs30=list(dist='unif', min=-0.15, max=0.15), 
                  Mag=list(dist='unif', min=-0.15, max=0.15)),
  check=list(check=list(dist='unif', min=0, max=1))
)

#Set up the same structure to links, unlinks and acceptance transformations as Model$skeleton
# Links: transforms the parameters onto a more convenient domain (e.g. for MCMC proposals)
# Unlinks: untransforms the parameters
# AcceptTrans: function to apply to proposed and current parameters to modify acceptance probability to account for parameter transformations
#              (see here e.g. https://umbertopicchini.files.wordpress.com/2017/12/transformed-proposals2.pdf)

Model$links <-Model$unlinks <- Model$acceptTrans <- Model$skeleton

#Currently, all parameters use a lower and upper bound
for (i in 1:length(Model$links)){
  if (is.list(Model$links[[i]])){
    for (j in 1:length(Model$links[[i]])){
      Model$links[[i]][[j]] <- 'ab_bounded'
      Model$unlinks[[i]][[j]] <- 'ab_bounded_inv'
      Model$acceptTrans[[i]][[j]] <- 'ab_bounded_acc'
    }
  } else {
    Model$links[[i]] <- 'ab_bounded'
    Model$unlinks[[i]] <- 'ab_bounded_inv'
    Model$acceptTrans[[i]] <- 'ab_bounded_acc'
  }
}

#Set lower and upper bounds for the parameters, using the range of the uniform prior plaved over the parameters.
Model$par_lb <- c()
Model$par_ub <- c()
  
for (i in 1:length(Model$Priors)){
  if (is.list(Model$Priors[[i]])){
    for (j in 1:length(Model$Priors[[i]])){
      if(Model$Priors[[i]][[j]]$dist == 'unif'){
        Model$par_lb = c(Model$par_lb, Model$Priors[[i]][[j]]$min)
        Model$par_ub = c(Model$par_ub, Model$Priors[[i]][[j]]$max)
      } else {
        stop('Please update Method.R to adjust acceptance probability to account for non-uniform prior before continuing.')
      }
    }
  } else {
    if(Model$Priors[[i]]$dist == 'unif'){
      Model$par_lb = c(Model$par_lb, Model$Priors[[i]]$min)
      Model$par_ub = c(Model$par_ub, Model$Priors[[i]]$max)
    } else {
      stop('Please update Method.R to adjust acceptance probability to account for non-uniform prior before continuing.')
    }
  }
}

# Get the binary regression function (currently not implemented, only have option for normal cdf)
# Model$BinR<- "pnorm" #"weibull" # "gompertz" 

# Implement higher order Bayesian priors?
Model$higherpriors<-TRUE

Model$center<-ExtractCentering(dir,haz,T)

Model$impacts <- list(labels = c('mortality', 'displacement', 'buildDam', 'buildDest', 'buildDamDest'), 
                      qualifiers = c('qualifierMort', 'qualifierDisp', 'qualifierBuildDam', 'qualifierBuildDest', 'qualifierBuildDamDest'),
                      sampled = c('mort_sampled', 'disp_sampled', 'buildDam_sampled', 'buildDest_sampled'))

#Modifiers to capture change in probability of building damage from 1st to subsequent events
#e.g. We may expect P(Unaffected -> Damaged) is smaller in an aftershock as the building has been strong
# enough to survive the first shock. Alternatively, the first shock may weaken the building and make
# it more susceptible to damage. 
Model$DestDam_modifiers <- c(1,1,1) 
# Modifier 1 = change in P(Unaffected -> Damaged or Destroyed) from 1st to subsequent hzds
# Modifier 2 = change in P(Unaffected -> Destroyed) from 1st to subsequent hzds
# Modifier 3 = change from P(Unaffected -> Destroyed) in 1st hzd to P(Damaged -> Destroyed) in subsequent hzds
# Probability is placed to power of modifier. Therefore:
#     - Modifier > 1 means probability decreases, modifier < 1 means probability increases
#     - Modifier must be greater than 0 to ensure probabilities less than 1. 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Linear predictor calculations (act to modify the expected damage values)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

GetLP<-function(ODD,Omega,Params,Sinc,notnans, split_GNI=T){
  #if (split_GNI){return(array(1, dim=c(NROW(ODD@data), 8)))}
  #else {return(rep(1, NROW(ODD@data)))}
  
  LP_ij <- array(NA, dim=NROW(ODD@data))
  
  LP_ij[notnans] <- Omega$vuln_coeff$Mag * (max(ODD@hazinfo$magnitudes) - Params$center$Mag$mean) / Params$center$Mag$sd #0 # Omega$vuln_coeff$itc #intercept term
  
  #could perform all centering outside before model fitting? may allow a bit of speedup
  
  #Population density term:
  LP_ij[notnans] <- LP_ij[notnans] + Omega$vuln_coeff$PDens * ((log(ODD@data$PDens[notnans]+1) - Params$center$PDens$mean)/Params$center$PDens$sd)
  LP_ij[notnans] <- LP_ij[notnans] + Omega$vuln_coeff$EQFreq * ((log(ODD@data$EQFreq[notnans]+0.1) - Params$center$EQFreq$mean)/Params$center$EQFreq$sd)
  for (vuln_term in names(Omega$vuln_coeff)[!(names(Omega$vuln_coeff) %in%  c('itc', 'PDens', 'GNIc', 'EQFreq', 'Mag'))]){
    #All remaining terms except GNIc:
    LP_ij[notnans] <- LP_ij[notnans] + Omega$vuln_coeff[[vuln_term]] * ((ODD@data[notnans, vuln_term] - Params$center[[vuln_term]]$mean)/Params$center[[vuln_term]]$sd)
  }
  
  #GNIc:
  if (split_GNI==F){ #don't split into the eight GNIc deciles:  
    LP_ij[notnans] <-  LP_ij[notnans] + Omega$vuln_coeff$GNIc * (log(ODD@data$GNIc[notnans]) - Params$center$GNIc$mean)/Params$center$GNIc$sd
    return(LP_ij)
  }
  
  LP_ijs <- array(NA, dim=c(NROW(ODD@data),8))
  
  get_GNIc_vuln <- function(ij){
    vuln_GNIc_ij <- Omega$vuln_coeff$GNIc * (log(ODD@data$GNIc[ij] * Sinc[Sinc$iso3==ODD@data$ISO3C[ij],]$value * 12.5) - Params$center$GNIc$mean)/Params$center$GNIc$sd
    return(vuln_GNIc_ij)
  }

  LP_ijs[notnans,] <- sweep(t(vapply(t(notnans), get_GNIc_vuln, numeric(8))), 1, LP_ij[notnans], '+')
  
  return(LP_ijs)
}

# GetLP_single is the equivalent of GetLP for a single grid-cell ij 
# Used in higher-level prior to calculate linear predictor for a given set of vulnerability terms
GetLP_single <- function(Omega, center, vuln_terms){
  #return(1)
  LP_ij <- Omega$vuln_coeff$Mag * (vuln_terms[['Mag']] - center$Mag$mean) / center$Mag$sd #0 #Omega$vuln_coeff$itc 
  
  LP_ij <- LP_ij + Omega$vuln_coeff$PDens * ((log(vuln_terms[['PDens']]+1) - center$PDens$mean)/center$PDens$sd)
  LP_ij <- LP_ij + Omega$vuln_coeff$EQFreq * ((log(vuln_terms[['EQFreq']]+0.1) - center$EQFreq$mean)/center$EQFreq$sd)
  
  for (vuln_term in names(Omega$vuln_coeff)[!(names(Omega$vuln_coeff) %in%  c('itc', 'PDens', 'GNIc', 'EQFreq', 'Mag'))]){
    #All remaining terms except GNIc:
    LP_ij <- LP_ij + Omega$vuln_coeff[[vuln_term]] * ((vuln_terms[[vuln_term]] - center[[vuln_term]]$mean)/center[[vuln_term]]$sd)
  }
  
  LP_ij <- LP_ij + Omega$vuln_coeff$GNIc * (log(vuln_terms[['GNIc']]) - center$GNIc$mean)/center$GNIc$sd
  
  return(LP_ij)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Generalised Linear Models for building damage & displacement calcs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Stochastic damage function process
stochastic<-function(n,eps){
  return(rnorm(n=n, mean=0, sd=eps))
  #return(rgammaM(n = n,mu = 1, sig_percent = eps ))
}

# Baseline hazard function h_0
h_0<-function(I,I0, Omega){
  ind<-I>I0
  h<-rep(0,length(I))
  h[ind]<- I[ind]-I0
  h[ind] <- h[ind] #* exp(Omega$Lambda1$alpha*(I-6.5))
  return(h)
}

# Binomial displacement calculator function
rbiny<-function(size,p) rbinom(n = 1,size,p);
Fbdisp<-function(lPopS,Dprime) mapply(rbiny,lPopS,Dprime);

fBD<-function(nbuildings, D_BD) mapply(rbiny, nbuildings, D_BD)

# Calculate the unscaled damage function
fDamUnscaled<-function(I,Params,Omega){ 
  (h_0(I,Params$I0,Omega) +
     stochastic(Params$Np,Omega$eps$local)) %>%return()
  #(h_0(I,Params$I0,Omega$theta) + 
  #   stochastic(Params$Np,Omega$eps$local)) %>%return()
}

addTransfParams <- function(Omega, I0=4.5){
  Omega$Lambda1$loc <- Omega$Lambda1$mu - I0 #h_0(Omega$Lambda1$nu, I0, Omega$theta)
  Omega$Lambda2$loc <- Omega$Lambda2$mu - I0 #h_0(Omega$Lambda2$nu, I0, Omega$theta)
  Omega$Lambda3$loc <- Omega$Lambda3$mu - I0 #h_0(Omega$Lambda3$nu, I0, Omega$theta)
  Omega$Lambda4$loc <- Omega$Lambda4$mu - I0 #h_0(Omega$Lambda4$nu, I0, Omega$theta)
  return(Omega)
}

# Calculate Mortality and Displacement probabilities from the unscaled damage
D_MortDisp_calc <- function(Damage, Omega){
  D_Mort_and_Disp <- pnorm(Damage, mean = Omega$Lambda1$loc, sd = Omega$Lambda1$sigma)
  D_Mort <- pnorm(Damage, mean = Omega$Lambda2$loc, sd = Omega$Lambda2$sigma)
  D_Disp <- D_Mort_and_Disp - D_Mort
  D_Disp <- ifelse(D_Disp<0, 0, D_Disp)
  return(rbind(D_Mort, D_Disp))
}

# Calculate Mortality and Displacement probabilities from the unscaled damage
D_DestDam_calc <- function(Damage, Omega, first_haz=T, DestDam_modifiers = c(1,1,1), ind_dam=c()){
  if(first_haz == T){
    D_Dest_and_Dam <- pnorm(Damage, mean = Omega$Lambda3$loc, sd = Omega$Lambda3$sigma)
    D_Dest <- pnorm(Damage, mean = Omega$Lambda4$loc, sd = Omega$Lambda4$sigma)
    D_Dam <- pmax(D_Dest_and_Dam - D_Dest, 0)
  } else {
    if (length(ind_dam>0)){
      #D_Dam = 1 - D_Dest
      D_Dest_and_Dam <- D_Dam <- D_Dest <- rep(0, length(Damage))
      D_Dest_and_Dam[-ind_dam] <- pnorm(Damage[-ind_dam], mean = Omega$Lambda3$loc, sd = Omega$Lambda3$sigma) ^ DestDam_modifiers[1]
      D_Dest_and_Dam[ind_dam] <- 1
      D_Dest <- pnorm(Damage, mean = Omega$Lambda4$loc, sd = Omega$Lambda4$sigma)
      D_Dest[-ind_dam] <- D_Dest[-ind_dam] ^ DestDam_modifiers[2]
      D_Dest[ind_dam] <- D_Dest[ind_dam] ^ DestDam_modifiers[3]
      D_Dam <- pmax(D_Dest_and_Dam - D_Dest, 0)
    } else {
      D_Dest_and_Dam <- pnorm(Damage, mean = Omega$Lambda3$loc, sd = Omega$Lambda3$sigma) ^ DestDam_modifiers[1]
      D_Dest <- pnorm(Damage, mean = Omega$Lambda4$loc, sd = Omega$Lambda4$sigma) ^ DestDam_modifiers[2]
      D_Dam <- D_Dest_and_Dam - D_Dest
      D_Dam <- ifelse(D_Dam<0, 0, D_Dam)
    }
  }
  return(rbind(D_Dest, D_Dam))
}

#when working with buildings, D_Disp is equivalent to D_BuildDam and D_Mort is equivalent to D_BuildDest
rmultinomy<-function(size, D_Disp, D_Mort, D_Rem) rmultinom(n=1, size, c(D_Disp, D_Mort, D_Rem))
Fbdam<-function(PopRem, D_Disp, D_Mort, D_Rem) mapply(rmultinomy, PopRem, D_Disp, D_Mort, D_Rem)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Log likelihood, posterior and higher-level prior distribution calculations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# These high-level priors are to include expert opinion in the 
# model calculations, not in the individual parameters themselves (for this, see priors)
Model$HighLevelPriors <-function(Omega,Model,modifier=NULL){
  
  min_PDens <- 0; max_PDens <- 109043.5 #minimum and maximum population density from the training dataset
  min_AveSchYrs <- 0.342; max_AveSchYrs <- 18 #minimum and maximum from all regions in GDL dataset
  min_LifeExp <- 24.511; max_LifeExp <- 85.413 #minimum and maximum from all regions in GDL dataset
  min_GNIc <- exp(5.887395); max_GNIc <- exp(12.23771) #minimum and maximum from all regions in GDL dataset
  min_Vs30 <- 98; max_Vs30 <- 2197 #minimum and maximum from all regions in soil stiffness dataset
  min_EQFreq <- 0; max_EQFreq <- 949.4231 #minimum and maximum from all regions in PGA dataset
  min_Mag <- 4.5; max_Mag <- 8.2
  
  linp_min <- GetLP_single(Omega, Model$center, vuln_terms=list(PDens=ifelse(Omega$vuln_coeff$PDens>0, min_PDens, max_PDens), 
                                                                AveSchYrs=ifelse(Omega$vuln_coeff$AveSchYrs>0, min_AveSchYrs, max_AveSchYrs),
                                                                LifeExp=ifelse(Omega$vuln_coeff$LifeExp>0, min_LifeExp, max_LifeExp),
                                                                GNIc=ifelse(Omega$vuln_coeff$GNIc>0, min_GNIc, max_GNIc),
                                                                Vs30=ifelse(Omega$vuln_coeff$Vs30>0, min_Vs30, max_Vs30),
                                                                EQFreq=ifelse(Omega$vuln_coeff$EQFreq>0, min_EQFreq, max_EQFreq),
                                                                Mag=ifelse(Omega$vuln_coeff$Mag>0, min_Mag, max_Mag)))
  
  linp_max <- GetLP_single(Omega, Model$center, vuln_terms=list(PDens=ifelse(Omega$vuln_coeff$PDens<0, min_PDens, max_PDens), 
                                                                AveSchYrs=ifelse(Omega$vuln_coeff$AveSchYrs<0, min_AveSchYrs, max_AveSchYrs),
                                                                LifeExp=ifelse(Omega$vuln_coeff$LifeExp<0, min_LifeExp, max_LifeExp),
                                                                GNIc=ifelse(Omega$vuln_coeff$GNIc<0, min_GNIc, max_GNIc),
                                                                Vs30=ifelse(Omega$vuln_coeff$Vs30<0, min_Vs30, max_Vs30),
                                                                EQFreq=ifelse(Omega$vuln_coeff$EQFreq<0, min_EQFreq, max_EQFreq),
                                                                Mag=ifelse(Omega$vuln_coeff$Mag<0, min_Mag, max_Mag)))
  
  #if(!is.null(modifier)) lp<-exp(as.numeric(unlist(modifier))) else lp<-1.
  lp <- c(linp_min, linp_max) 
  if(Model$haz=="EQ"){
  
    # Lower and upper bounds on the impacts at I_ij = 4.6, 6, and 9
    # in the order (Mort, DispMort, Dest, DamDest),
    # where DispMort is the sum of the probabilities of displacement and mortality
    # and DamDest is the sum of the probabilities of building damage and destruction.
    
    Upp_bounds_4.6 <- c(0.03, 0.1, 0.03, 0.15)
    Low_bounds_6 <- c(0, 0.00001, 0, 0.00001)
    Upp_bounds_6 <- c(0.15, 0.6, 0.3, 0.75)
    Low_bounds_9 <- c(0.00005,0.01,0.001,0.15)
    Upp_bounds_9 <- c(0.7,0.995,0.99,0.995)
    
    Upp_bounds_4.6_zero_lp <- c(0.0005, 0.005, 0.001, 0.002)
    Low_bounds_6_zero_lp <- c(0, 0.001, 0, 0.001)
    Upp_bounds_6_zero_lp <- c(0.01, 0.1, 0.05, 0.1)
    Low_bounds_9_zero_lp <- c(0.0001,0.1,0.05,0.2)
    Upp_bounds_9_zero_lp <- c(0.2,0.9,0.5,0.8)
    
    HLP_impacts <- function(I_ij, lp, Omega){
      rbind(apply(D_MortDisp_calc(h_0(I_ij, I0=4.5, Omega=Omega) + lp, Omega),2,cumsum), 
            apply(D_DestDam_calc(h_0(I_ij, I0=4.5, Omega=Omega) + lp, Omega), 2,cumsum))
    }
    
    adder <- 0
    adder <- sum(apply(HLP_impacts(4.6, lp, Omega), 2, function (x) sum(x > Upp_bounds_4.6))) + 
      apply(HLP_impacts(4.6, 0, Omega), 2, function (x) sum(x > Upp_bounds_4.6_zero_lp))
    
    adder <- adder + sum(apply(HLP_impacts(6, lp, Omega), 2, 
                               function (x) sum(c(x > Upp_bounds_6, x<Low_bounds_6))))  + 
      apply(HLP_impacts(6, 0, Omega), 2, function (x) sum(x > Upp_bounds_6_zero_lp, x < Low_bounds_6_zero_lp))
    
    adder <- adder + sum(apply(HLP_impacts(9, lp, Omega), 2, 
                               function (x) sum(c(x > Upp_bounds_9, x<Low_bounds_9)))) + 
      apply(HLP_impacts(9, 0, Omega), 2, function (x) sum(x > Upp_bounds_9_zero_lp, x < Low_bounds_9_zero_lp))
    
    #check that at intensity 7, D_disp > D_mort and D_builddam > D_builddest
    impact_intens_7 <- HLP_impacts(7, lp, Omega)
    adder <- adder + sum(impact_intens_7[1,] > impact_intens_7[2,]) + sum(impact_intens_7[3,] > impact_intens_7[4,])

    return(adder) #looseend: need to address when including modifiers
    
  } else if(Model$haz=="TC"){
    # Would need to be udpated:
    Dfun<-function(I_ij) h_0(I = I_ij,I0 = 3,Omega = Omega) 
    Dispfun<-function(I_ij) c(BinR(Dfun(I_ij)*Dfun(I_ij)*Omega$Lambda1$kappa+Omega$Lambda1$nu*Dfun(I_ij) + Omega$Lambda1$omega,Omega$zeta)%o%lp)
    Damfun<-function(I_ij) c(BinR(Dfun(I_ij),Omega$zeta)%o%lp)
    
    # Add lower bound priors:
    adder<-sum(50*pweibull(c(Dispfun(3.05),Damfun(3.05)),3,0.001))
    # Middle range priors:
    # adder<-adder+sum(15*(1- pweibull(c(Dispfun(3.53),Damfun(3.53)),3,0.005)) +15*pweibull(c(Dispfun(3.53),Damfun(3.53)),15,0.8) )
    # Upper bound priors:
    adder<-adder+sum(50*(1-pweibull(c(Dispfun(45),Damfun(45)),30,0.85)))
    
    return(-adder)
    
  }
  # I<-seq(from=4.05,to=9.5,length.out = 200)
  # Intensity<-data.frame(I_ij=rep(I,2),value=c(vapply(I,Dispfun,numeric(1)),vapply(I,BDfun,numeric(1))),
  #                       term=c(rep("Displacement",200),rep("Building Damage",200)))
  
}

# Get the log-likelihood for the displacement data
CalcPolyDist <- function(Y,  kernel_sd, kernel, cap){
  
  if (any(c(is.nan(Y[,'observed']),is.nan(Y[,'sampled'])))) return(0)
  
  Dist <- 0
  k <- 10
  cap <- -100
  if (kernel=='log'){
    Dist = abs(log(Y[,'observed']+k) - log(Y[,'sampled']+k)) * unlist(kernel_sd)[Y[,'impact']]
  } else if (kernel == 'loglaplace'){ #use a laplace kernel 
    warning('Kernel_sd are currently set to act more as weights rather than standard deviations, so would need to be adjusted for this kernel.')
    Dist = log(dloglap(Y[,'observed']+k, location.ald = log(Y[,'sampled']+k), scale.ald = unlist(kernel_sd)[Y[,'impact']], tau = 0.5, log = FALSE)/
                        (1-ploglap(k, location.ald = log(Y[,'sampled']+k), scale.ald = unlist(kernel_sd)[Y[,'impact']], tau = 0.5, log = FALSE)))
  } else if (kernel == 'lognormal'){ #use a lognormal kernel 
    warning('Kernel_sd are currently set to act more as weights rather than standard deviations, so would need to be adjusted for this kernel.')
    Dist = log(dlnormTrunc(Y[,'observed']+k, log(Y[,'sampled']+k), sdlog=unlist(kernel_sd)[Y[,'impact']], min=k))
  } else {
    stop('Working with an unsupported distance kernel.')
  }
  
  Dist[which(is.na(Dist))] <- cap
  return(sum(Dist))
}

# Plot to compare normal and laplace kernels
# xrang <- seq(1000,2000,1)
# xobs <- 1500
# lapval <- dloglap(xobs+k, location.ald = log(xrang+k), scale.ald = AlgoParams$kernel_sd$displacement, tau = 0.5, log = FALSE)/
#   (1-ploglap(k, location.ald = log(xrang+k), scale.ald = AlgoParams$kernel_sd$displacement, tau = 0.5, log = FALSE))
# normval <- dlnormTrunc(xobs+k, log(xrang+k), sdlog=AlgoParams$kernel_sd$displacement, min=k)
# plot(xrang, lapval)
# points(xrang, normval, col='red')
# q <- 0.975
# xobs <- 1000
# qloglap(q, location.ald = log(xobs+k), scale.ald = 0.09, tau = 0.5, log = FALSE)/
#    (1-ploglap(k, location.ald = log(xobs+k), scale.ald = 0.1, tau = 0.5, log = FALSE))
# qlnormTrunc(q, log(xobs+k), sdlog=epsilon, min=k)

#Plot of bimodal kernel where we have conflicting data sources
# xrang <- seq(0,1000000,50)
# xobs1 <- 684800
# xobs2 <- 237655
# normval <- 0.5 * dlnormTrunc(xobs1+k, log(xrang+k), sdlog=AlgoParams$kernel_sd$displacement, min=k) + 0.5 * dlnormTrunc(xobs2+k, log(xrang+k), sdlog=AlgoParams$kernel_sd$displacement, min=k)
# plot(xrang, normval, col='red', type='l', xlab='Displacement', ylab="'Likelihood' assigned to simulated data")
# abline(v=xobs1)
# abline(v=xobs2)

# for (x in c(1,100,10000)){
#   dist_calc <- function(bound, x){
#     return(0.9*abs(log(bound+10)-log(x+10)))
#   }
#   print(uniroot(function(bound) dist_calc(bound,x=x)-1,c(-10,x))$root)
#   print(uniroot(function(bound) dist_calc(bound,x=x)-1,c(x,x*100))$root)
# }
# 
# plot(seq(0,1000,1), abs(log(300+k) - log(seq(0,1000,1)+k)))

# -------------------------------------------------------------------------------------------------------------------
# ------------------------------- Pulling out and breaking down the distances ---------------------------------------
# -------------------------------------------------------------------------------------------------------------------

SamplePolyImpact <-function(dir,Model,proposed,AlgoParams, dat='Train'){
  
  # Load ODD files
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  if(AlgoParams$kernel == 'crps' | AlgoParams$kernel == 'crps_with_mean'){
    #if using continuous ranked probability score as distance function, need multiple samples per particle
    AlgoParams$Np <- AlgoParams$Np * AlgoParams$m_CRPS
  }
  
  # Parallelise appropriately
  if(AlgoParams$AllParallel){
    # Task parallelism: this parallelisation calculates each event side-by-side, which is ideal if we have many CPU threads available and many ODD objects
    cores<-AlgoParams$cores
    AlgoParams$cores<-AlgoParams$NestedCores
    # When using task parallelisation, put the heaviest files first for optimisation reasons
    x <- file.info(paste0(folderin,ufiles))
    ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
    
    tmpFn<-function(filer){
      # Extract the ODD object
      ODDy<-readRDS(paste0(folderin,filer))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      #ODDy@fIndies<-Model$fIndies
      ODDy@impact%<>%as.data.frame.list()
      
      
      ODDy@impact$event_id = as.numeric(gsub(".*_(\\d+)$", "\\1", filer))
      
      
      # Apply DispX
      tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      #if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",filer));return(-Inf)}
      
      return(tLL) #SMC-CHANGE
      # Weight the likelihoods based on the number of events for that country
      cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==ODDy@gmax$iso3[1]]
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Return the average log-likelihood
      if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL))
      else return(cWeight*mean(tLL,na.rm=T))
    }
    return(pmap(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores), rbind)) # SMC-CHANGE
    #return(sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
  } 
  # } else {
  #   
  #   # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
  #   for(i in 1:length(ufiles)){
  #     # Extract the ODD object
  #     ODDy<-readRDS(paste0(folderin,ufiles[i]))
  #     # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
  #     ODDy@fIndies<-Model$fIndies
  #     ODDy@gmax%<>%as.data.frame.list()
  #     # Apply DispX
  #     tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams),
  #                   error=function(e) NA)
  #     # If all is good, add the LL to the total LL
  #     if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));return(-Inf)}
  #     # Weight the likelihoods based on the number of events for that country
  #     
  #     cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==ODDy@gmax$iso3[1]]
  #     # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
  #     maxLL<-max(tLL,na.rm = T)
  #     # Add the likelihood to the list of all events.
  #     if(expLL) {LL<-LL+cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)
  #     } else LL<-LL+cWeight*mean(tLL,na.rm=T)
  #   }
  #   return(LL)
  # }
}

SamplePointImpact <- function(dir,Model,proposed,AlgoParams,expLL=T, dat='Train'){
  # Load BD files
  folderin<-paste0(dir,"IIDIPUS_Input/BDobjects/")
  
  ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  # if(AlgoParams$kernel == 'crps'){
  #   #if using continuous ranked probability score as distance function, need multiple samples per particle
  #   AlgoParams$Np <- AlgoParams$Np * AlgoParams$m_CRPS
  # }
  
  # Parallelise appropriately
  if(AlgoParams$AllParallel){
    # Task Parallelisation: this parallelisation calculates each event side-by-side, which is ideal if we have many CPU threads available and many BD objects
    cores<-AlgoParams$cores
    AlgoParams$cores<-AlgoParams$NestedCores
    # When using task parallelisation, put the heaviest files first for optimisation reasons
    x <- file.info(paste0(folderin,ufiles))
    ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
    
    tmpFn<-function(filer){
      # Extract the BD object
      BDy<-readRDS(paste0(folderin,filer))
      
      if(nrow(BDy@data)==0){return()}
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      #BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=F),
                    error=function(e) NA)
      
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",filer));return(-Inf)}
      if (length(tLL > 1)){ return(tLL)}
      return(tLL) #SMC-CHANGE
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Return the average log-likelihood
      cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==BDy$ISO3C[which(!is.na(BDy$ISO3C))[1]]] #looseend
      if(is.na(cWeight)){return(-Inf)}                               
      if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)) 
      else return(cWeight*mean(tLL,na.rm=T))
      
    }
    return(do.call(abind, list(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores), along=1))) #SMC-CHANGE #d-change
    #return(LL + sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
  }  
  # } else {
  #   
  #   # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
  #   for(i in 1:length(ufiles)){
  #     # Extract the BD object
  #     BDy<-readRDS(paste0(folderin,ufiles[i]))
  #     # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
  #     BDy@fIndies<-Model$fIndies
  #     # Apply BDX
  #     tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=T),
  #                   error=function(e) NA)
  #     # If all is good, add the LL to the total LL
  #     if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));return(-Inf)}
  #     # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
  #     maxLL<-max(tLL,na.rm = T)
  #     # Add the likelihood to the list of all events.
  #     if(expLL) LL<-LL+log(mean(exp(tLL-maxLL),na.rm=T))+maxLL
  #     else LL<-LL+mean(tLL,na.rm=T)
  #   }
  #   
  #   return(LL)
  # }
}

SampleImpact <- function(dir,Model,proposed,AlgoParams,expLL=T, dat='Train'){

  impact_sample_poly<-SamplePolyImpact(dir,Model,proposed,AlgoParams, dat=dat)
  impact_sample_point<- SamplePointImpact(dir, Model, proposed, AlgoParams, expLL=T, dat=dat)

  return(list(poly=impact_sample_poly, point=impact_sample_point))
}

crps <- function(sample, obs){
  sample <- sort(sample)
  m <- length(sample)
  crps <- 0
  for (i in 1:m){
    crps <- crps + (sample[i]-obs)*(m*as.numeric(obs<sample[i]) -i + 0.5)
  }
  crps <- (crps*2)/(m^2)
  return(crps)
}

logTarget_CRPS <- function(impact_sample, AlgoParams, dist_poly_means=NULL){
  
  crps_eval <- function(n, obs, weights){
    samples_allocated <- 1:AlgoParams$m_CRPS * n
    samples_combined <- sapply(impact_sample$poly[samples_allocated], function(x){x$sampled})
    crps_vals <- sapply(1:NROW(samples_combined), function(i){crps(log(samples_combined[i,]+10), log(obs[i]+10))})*weights
    return(sum(crps_vals))
  }
  
  dist_poly <- unlist(mclapply(1:AlgoParams$Np, crps_eval, mc.cores=1, obs=impact_sample$poly[[1]]$observed, weights=unlist(AlgoParams$kernel_sd)[impact_sample$poly[[1]]$impact]))
  
  if(AlgoParams$kernel == 'crps_with_mean'){
    if (is.null(dist_poly_means)){
      k <- 10
      Y <- impact_sample$poly[[1]] %>% filter(impact %in% c('mortality', 'displacement'))
      dist_poly_means <- sum(abs(log(Y[,'observed']+k) - log(Y[,'mean']+k)) * unlist(AlgoParams$kernel_sd)[Y[,'impact']])
    }
  } else {
    dist_poly_means = 0
  }
  
  #is there a way to do this using a scoring rule as well? : 
  sumPointDat_dists <- function(PointDat_p){
    Dist_0.5 <- which(names(PointDat_p) %in% c('N12', 'N21', 'N23', 'N32'))
    Dist_1 <- which(names(PointDat_p) %in% c('N13', 'N31'))
    return(0.5*sum(PointDat_p[Dist_0.5])+sum(PointDat_p[Dist_1]))
  }
  

  if (length(impact_sample$point) > 0){
    dist_point <- apply(impact_sample$point, 2, sumPointDat_dists) * 0.1 #LOOSEEND: Currently using 0.1 as point data distance is around 10 times the poly data distance, but no real justification
  } else {
    dist_point <- 0
  }

  print(paste0('Dist_agg: ',dist_poly, ' Dist_agg_means: ', dist_poly_means,' Dist_sat: ', dist_point))
  dist_tot <- dist_poly + dist_poly_means + dist_point
  
  return(cbind(dist_poly, dist_poly_means, dist_point))
  #return(dist_tot)
}

CalcDist <- function(impact_sample, AlgoParams, dist_poly_means=NULL){
  
  if (AlgoParams$kernel == 'crps' | AlgoParams$kernel == 'crps_with_mean'){ return(logTarget_CRPS(impact_sample, AlgoParams, dist_poly_means))}
  
  sumDists <- function(poly_p){
    dist_p = 0 
    nrows <- NROW(poly_p)
    for (i in 1:nrows){
      dist_p = dist_p + CalcPolyDist(poly_p[i,], kernel_sd = AlgoParams$kernel_sd, kernel=AlgoParams$kernel, cap=-300)
    }
    return(dist_p)
  }
  dist_poly <- unlist(mclapply(impact_sample$poly,FUN = sumLLs,mc.cores = 1)) # sum log likelihoods
  
  sumPointDat_dists <- function(PointDat_p){
    Dist_0.5 <- which(names(PointDat_p) %in% c('N12', 'N21', 'N23', 'N32'))
    Dist_1 <- which(names(PointDat_p) %in% c('N13', 'N31'))
    return(0.5*sum(PointDat_p[Dist_0.5])+sum(PointDat_p[Dist_1]))
  }
  
  if (length(impact_sample$point) > 0){
    dist_point <- apply(impact_sample$point, 2, sumPointDat_dists) 
  } else {
    dist_point <- 0
  }

  print(paste0('Dist_agg: ',dist_poly, ' Dist_sat: ', dist_point))
  dist_tot <- dist_poly + dist_point 
  
  return(dist_tot)
}

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Currently not implemented:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Building damage baseline hazard function hBD_0
# hBD<-function(Ab,Population,rho,center){
#   exp(-rho$A*(log(Ab)-center$A) - rho$H*(log(Population)-center$H))
# }

# For building damage assessment data
# fDamageBuilding<-function(BD,I,Params,Omega,linp,Ik){
#   fDamUnscaled(I,Params,Omega)*linp #*hBD(BD$Ab,BD$Population,Omega$rho,Params$center[c("A","H")]),
# }

# qualifierDisp<-function(Disp,qualifier,mu) {
#   if(qualifier%in%c("total","approximately")) return(Disp)
#   else if(qualifier=="more than") {
#     return(vapply(Disp,function(Disp) rgammaM(n=1,mu = mu$muplus*Disp,
#                                               sig_percent = mu$sigplus),numeric(1)))
#   } else if(qualifier=="less than") {
#     return(vapply(Disp,function(Disp) rgammaM(n=1,mu = mu$muminus*Disp,
#                                               sig_percent = mu$sigminus),numeric(1)))
#   } else stop(paste0("qualifier is not recognised - ",qualifier))
# }

# lBD<-function(D_B, BD_params){
#   #input: probability of Building Damage D^B
#   #output: probability of Building Destruction D^BD
#   relative_probs = BDprob(D_B, BD_params) %>% mutate('damaged' = rowSums(.[1:4]))
#   D_BD = (relative_probs[['destroyed']] + relative_probs[['severe']])/rowSums(relative_probs)
#   return(D_BD)
# }

# Binary regression function (currently using normal distribution)
# if(Model$BinR=="weibull") {
#   BinR<-function(x,zeta) pweibull(x,shape=zeta$k,scale=zeta$lambda)
# } else if(Model$BinR=="gompertz") {
#   BinR<-function(x,zeta) flexsurv::pgompertz(x,shape=zeta$varrho,rate=zeta$eta)
# } else stop("Incorrect binary regression function name, try e.g. 'weibull'")

# GetIsoWeights<-function(dir){
#   Dispy<-readRDS(paste0(dir,"IIDIPUS_Input/DispData_EQ_V2.Rdata"))
#   WWW<-Dispy%>%group_by(iso3)%>%summarise(weights=1/length(gmax),.groups="drop_last")
# }
# 
# Model$IsoWeights<-GetIsoWeights(dir)
# Model$IsoWeights %<>% add_row(iso3='ABC', weights=1)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Old Code:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

############### Linear predictor function - ugly, but very fast!

# llinpred<-function(Params,beta,center,func) { # Note all of these are vectors, not single values
#   # if(is.null(func)) func<-replicate(length(beta),function(x) x)
#   linp<-list()
#   for(j in 1:length(Params)) {
#     vars<-Params[[j]]$var
#     # This part uses the name of the variable to search for the following:
#     # 1) link function - func, 2) centering value - center, 3) linear predictor parameterisation - beta
#     linp[[names(Params[j])]]<-prod(exp(vapply(1:nrow(vars),
#                                 function(i) beta[[vars$variable[i]]]*(do.call(func[[vars$variable[i]]],
#                                                                               list(vars$value[i]-center[[vars$variable[i]]] ))),
#                                 FUN.VALUE = numeric(1))))
#   }
#   return(linp)
# }
# 
# # Quicker version of llinpred for when only one vars$variable exists
# dlinpred<-function(vars,beta,center,func) { # Note all of these are vectors, not single values
#   # Quicker version of llinpred for when only one vars$variable exists
#   return(exp(beta[[vars$variable[1]]]*(do.call(func[[vars$variable[1]]],list(vars$value-center[[vars$variable[1]]] )))))
# }

# 
# WeibullScaleFromShape<-function(shape,center){
#   # Scale parameter
#   return(center/(log(2)^(1/shape)))
# }
# # Linear predictor for subnational variables
# locpred<-function(x,params){
#   exp(params$M*(pweibull(x,shape=params$shape,scale=params$scale)-0.5))
# }
# 
# GDPlinp<-function(ODD,Sinc,beta,center,notnans){
#   iGDP<-as.numeric(factor(interaction(ODD@data$GDP, ODD@data$ISO3C), levels=unique(interaction(ODD@data$GDP, ODD@data$ISO3C))))
#   dGDP<-data.frame(ind=iGDP[notnans],GDP=ODD@data$GDP[notnans],iso=ODD@data$ISO3C[notnans])
#   dGDP%<>%group_by(ind)%>%summarise(value=log(unique(GDP)*(Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]),"value"]/
#                                       Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]) & Sinc$variable=="p50p100","value"])),
#                               income=Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]),"variable"],
#                               variable="dollar",
#                               iso3c=unique(dGDP$iso[dGDP$ind==unique(ind)]),.groups = 'drop_last')
#   dGDP$linp<-rep(NA_real_,nrow(dGDP))
#   # Autocalculate the Weibull scale parameter using the shape and centering values
#   tp<-list(shape=beta$k,M=beta$M,
#         scale=WeibullScaleFromShape(shape = beta$k,center = center$dollar))
#   dGDP$linp[which(!is.na(dGDP$ind))]<-locpred(dGDP$value[!is.na(dGDP$ind)],tp)
#   # dGDP$linp[which(!is.na(dGDP$ind))]<-dlinpred(dGDP[!is.na(dGDP$ind),],beta,center,fIndies)
#   
#   return(list(dGDP=dplyr::select(dGDP,c(ind,income,linp)),iGDP=iGDP))
# }
# 
# Plinpred<-function(Pdens,beta,center,notnans){
#   
#   Plinp<-rep(NA_real_,length(Pdens))
#   # Plinp[notnans]<-dlinpred(data.frame(variable=rep("Pdens",length(Pdens[notnans])),value=log(Pdens[notnans]+1)),beta,center,fIndies)
#   # Autocalculate the Weibull scale parameter using the shape and centering values
#   tp<-list(shape=beta$k,M=beta$M,
#         scale=WeibullScaleFromShape(shape = beta$k,center = center$Pdens))
#   Plinp[notnans]<-locpred(log(Pdens[notnans]+1),tp)
#   return(Plinp)
#   
# }

#################### 'Likelihood' calculations:

# LL_beta_apply<-function(b,value,BD_params) do.call(BD_params$functions[[value]],as.list(c(x=b,unlist(BD_params$Params[[value]]))))
# 
# BDprob<-function(b,BD_params){
#   
#   lls<-array(dim=c(length(b),length(BD_params$Params)))
#   for(i in 1:length(BD_params$Params)){
#     value<-(names(BD_params$Params))[i]
#     lls[,i]<-vapply(b,FUN = LL_beta_apply,FUN.VALUE = numeric(1),
#                     value=value,BD_params=BD_params)
#   }
#   lls%<>%as.data.frame.array();colnames(lls)<-names(BD_params$Params)
#   return(lls)
# }
# 
# predBD<-function(b,BD_params){
#   lls<-BDprob(b,BD_params)
#   return(mean(apply(lls, 1, function(x) sample(1:5, 1, prob=x)))) #LOOSEEND: Mean?
# }
# 
# LL_BD<-function(b,classified,BD_params){
#   
#   lls<-BDprob(b,BD_params)
#   if(classified=="Damaged") {
#     tmp<-rowSums(lls)
#     # Sum of all rows for the classifications that predict at least some damage
#     return(log((tmp-lls[["notaffected"]])/tmp))
#   }
#   
#   out<-log(lls[[classified]]/rowSums(lls))
#   # machine level precision ~ -300
#   out[is.infinite(out)]<--300
#   return(out)
#   
# }
# 
# sampleBDdamage<-function(grading,n=10){
#   
#   vapply(1:length(grading),
#          FUN = function(i) {
#            
#            if(grading[i]=="Damaged") return(runif(n))
#            shapes<-Model$BD_params$Params[[grading[i]]]
#            return(rbeta(n,shape1 = shapes$shape1,shape2 = shapes$shape2))
#            
#          },
#          FUN.VALUE = numeric(n)
#   )
#   
# }

# Log-likelihood for displacement (ODD) objects
# LL_Displacement<-function(dir,Model,proposed,AlgoParams,expLL=T){
#   
#   # Load ODD files
#   folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
#   ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
#   
#   # Parallelise appropriately
#   if(AlgoParams$AllParallel){
#     # Task parallelism: this parallelisation calculates each event side-by-side, which is ideal if we have many CPU threads available and many ODD objects
#     cores<-AlgoParams$cores
#     AlgoParams$cores<-AlgoParams$NestedCores
#     # When using task parallelisation, put the heaviest files first for optimisation reasons
#     x <- file.info(paste0(folderin,ufiles))
#     ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))]) #looseend
#     
#     tmpFn<-function(filer){
#       # Extract the ODD object
#       ODDy<-readRDS(paste0(folderin,filer))
#       # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
#       ODDy@fIndies<-Model$fIndies
#       ODDy@impact%<>%as.data.frame.list()
#       # Apply DispX
#       tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams),
#                     error=function(e) NA)
#       # If all is good, add the LL to the total LL
#       if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",filer));return(-Inf)}
#       
#       return(tLL) #SMC-CHANGE
#       # Weight the likelihoods based on the number of events for that country
#       cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==ODDy@gmax$iso3[1]]
#       # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
#       maxLL<-max(tLL,na.rm = T)
#       # Return the average log-likelihood
#       if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL))
#       else return(cWeight*mean(tLL,na.rm=T))
#     }
#     return(do.call(rbind, mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))) # SMC-CHANGE
#     #return(sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
#     
#   } else {
#     # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
#     LL <- NULL
#     for(i in 1:length(ufiles)){
#       # Extract the ODD object
#       ODDy<-readRDS(paste0(folderin,ufiles[i]))
#       # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
#       ODDy@fIndies<-Model$fIndies
#       ODDy@impact%<>%as.data.frame.list()
#       # Apply DispX
#       tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams),
#                     error=function(e) NA)
#       # If all is good, add the LL to the total LL
#       if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));return(-Inf)}
#       
#       LL <- rbind(LL, tLL) #if want to return LL's for all events
#       next
#       # Weight the likelihoods based on the number of events for that country
#       cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==ODDy@gmax$iso3[1]]
#       # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
#       maxLL<-max(tLL,na.rm = T)
#       # Add the likelihood to the list of all events.
#       if(expLL) {LL<-LL+cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)
#       } else LL<-LL+cWeight*mean(tLL,na.rm=T)
#     }
#     return(LL)
#   }
# }
# 
# # Log-likelihood for building damage (BD) objects
# LL_Buildings<-function(dir,Model,proposed,AlgoParams,expLL=T){
#   # Load BD files
#   folderin<-paste0(dir,"IIDIPUS_Input/BDobjects/")
#   ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
#   
#   # Parallelise appropriately
#   if(AlgoParams$AllParallel){
#     # Task Parallelisation: this parallelisation calculates each event side-by-side, which is ideal if we have many CPU threads available and many BD objects
#     cores<-AlgoParams$cores
#     AlgoParams$cores<-AlgoParams$NestedCores
#     # When using task parallelisation, put the heaviest files first for optimisation reasons
#     x <- file.info(paste0(folderin,ufiles))
#     ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
#     
#     tmpFn<-function(filer){
#       # Extract the BD object
#       BDy<-readRDS(paste0(folderin,filer))
#       if(nrow(BDy@data)==0){return(0)}
#       # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
#       BDy@fIndies<-Model$fIndies
#       # Apply BDX
#       tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=T),
#                     error=function(e) NA)
#       
#       # If all is good, add the LL to the total LL
#       if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",filer));return(-Inf)}
#       return(tLL) #SMC-CHANGE
#       # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
#       maxLL<-max(tLL,na.rm = T)
#       # Return the average log-likelihood
#       cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==BDy$ISO3C[which(!is.na(BDy$ISO3C))[1]]] #looseend
#       if(is.na(cWeight)){return(-Inf)}                               
#       if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)) 
#       else return(cWeight*mean(tLL,na.rm=T))
#       
#     }
#     return(do.call(rbind, mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))) #SMC-CHANGE #d-change
#     #return(LL + sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
#     
#   } else {
#     LL <- NULL
#     # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
#     for(i in 1:length(ufiles)){
#       # Extract the BD object
#       BDy<-readRDS(paste0(folderin,ufiles[i]))
#       # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
#       BDy@fIndies<-Model$fIndies
#       # Apply BDX
#       tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=T),
#                     error=function(e) NA)
#       # If all is good, add the LL to the total LL
#       if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));return(-Inf)}
#       LL <- rbind(LL, tLL) #if want to return LL's for all events
#       next
#       # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
#       maxLL<-max(tLL,na.rm = T)
#       # Add the likelihood to the list of all events.
#       if(expLL) LL<-LL+log(mean(exp(tLL-maxLL),na.rm=T))+maxLL
#       else LL<-LL+mean(tLL,na.rm=T)
#     }
#     
#     return(LL)
#   }
# }
# 
# # Bayesian Posterior distribution: This is for a group of ODD objects with observed data
# logTarget<-function(dir,Model,proposed,AlgoParams,expLL=T){
#   
#   # Apply higher order priors
#   if(!is.null(Model$HighLevelPriors)){
#     HP<-Model$HighLevelPriors(proposed,Model)
#     print(paste0("Higher Level Priors = ",HP))
#   }
#   
#   
#   # Approximate Bayesian Computation rejection
#   #if(HP>AlgoParams$ABC) return(-5000000) # Cannot be infinite, but large (& negative) to not crash the Adaptive MCMC algorithm
#   if (HP> AlgoParams$ABC) return(-Inf)
#   
#   # Add the log-likelihood values from the ODD (displacement) objects
#   LL<-LL_Displacement(dir,Model,proposed,AlgoParams,expLL=T)
#   #print(paste0("LL Displacements = ",LL)) ; sLL<-LL
#   #if (any(LL == 0)){ #SMC-CHANGE
#   #  return(Inf)
#   #}
#   # Add the log-likelihood values from the BD (building damage) objects
#   #LL%<>%LL_Buildings(dir,Model,proposed,AlgoParams,expLL=T)
#   LL_Build<-LL_Buildings(dir, Model, proposed, AlgoParams, expLL=T)
#   #print(paste0("LL Building Damages = ",LL-sLL))
#   
#   posterior<-rbind(-LL, LL_Build)#LL #+HP
#   # Add Bayesian priors
#   if(!is.null(Model$priors)){
#     posterior<-posterior+sum(Priors(proposed,Model$priors),na.rm=T)
#   }
#   #print(paste0("Posterior = ",posterior))
#   
#   return(posterior)
#   
# }

# logTargetSingle<-function(ODDy,Model,Omega,AlgoParams){
#   
#   # HP<-0
#   # # Apply higher order priors
#   # if(!is.null(Model$HighLevelPriors)){
#   #   HP<-Model$HighLevelPriors(Omega,Model,modifier = tryCatch(ODDy@modifier,error=function(e) NULL))
#   #   # print(paste0("Higher Level Priors = ",HP))
#   # }
#   # 
#   # LL<-HP
#   
#   tLL<-tryCatch(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params,
#                       LL = T,Method = AlgoParams),
#                 error=function(e) NA)
#   
#   if(any(is.infinite(tLL)) | all(is.na(tLL))) stop("Failed to calculate LL")
#   
#   # LL<-LL+tLL
#   
#   return(tLL)
#   
# }

############## Parameter set-up (would be useful if links are not all the same)

# Link functions (MUST BE SAME LENGTH AS OMEGA)
# Model$links<-list(
#   Lambda1=list(nu='ab_bounded',omega='ab_bounded'),
#   Lambda2=list(nu='ab_bounded',omega='ab_bounded'),
#   Lambda3=list(nu='ab_bounded',omega='ab_bounded'),
#   Lambda4=list(nu='ab_bounded',omega='ab_bounded'), 
#   theta=list(e='ab_bounded'),
#   eps=list(eps='ab_bounded'),
#   vuln_coeff=list(itc='ab_bounded',PDens='ab_bounded', AveSchYrs='ab_bounded', 
#                   LifeExp='ab_bounded', GNIc='ab_bounded', Vs30='ab_bounded', EQFreq='ab_bounded') 
# )
# 
# Model$unlinks<-list(
#   Lambda1=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
#   Lambda2=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
#   Lambda3=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
#   Lambda4=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
#   theta=list(e='ab_bounded_inv'),
#   eps=list(eps='ab_bounded_inv'),
#   vuln_coeff=list(itc='ab_bounded_inv',PDens='ab_bounded_inv', AveSchYrs='ab_bounded_inv', 
#                   LifeExp='ab_bounded_inv', GNIc='ab_bounded_inv', Vs30='ab_bounded_inv', 
#                   EQFreq='ab_bounded_inv'))
# 
# Model$acceptTrans <- list(
#   Lambda1=list(nu='ab_bounded_acc', omega='ab_bounded_acc'),
#   Lambda2=list(nu='ab_bounded_acc', omega='ab_bounded_acc'),
#   Lambda3=list(nu='ab_bounded_acc', omega='ab_bounded_acc'),
#   Lambda4=list(nu='ab_bounded_acc', omega='ab_bounded_acc'),
#   theta=list(e='ab_bounded_acc'), 
#   eps=list(eps='ab_bounded_acc'),
#   vuln_coeff=list(itc='ab_bounded_acc',PDens='ab_bounded_acc', AveSchYrs='ab_bounded_acc', 
#                   LifeExp='ab_bounded_acc', GNIc='ab_bounded_acc', Vs30='ab_bounded_acc', 
#                   EQFreq='ab_bounded_acc') 
# )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Choose your vulnerability variables to be correlated with
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

vulnerabilityVars<-function(DispData,Model){
  
  vulny<-data.frame()
  if(!is.null(Model$modifiers$WB))      vulny%<>%rbind(GetWB_Vals(DispData))
  # These aren't ready yet, but are easy to do if you're willing!
  # if(!is.null(Model$modifiers$INFORM))  vulny%<>%rbind(GetINFORM_Vals(DispData))
  # if(!is.null(Model$modifiers$WID))  vulny%<>%rbind(GetWID_Vals(DispData))
  
  if(length(vulny)==0) stop("Please make sure modifier vulnerability databases are defined in Model")
  
  return(vulny)
  
}


# logTargetODDRIN<-function(dir,Model,proposed,AlgoParams,expLL=T){
#   
#   HP<-0
#   # Apply higher order priors
#   if(!is.null(Model$HighLevelPriors)){
#     HP<-Model$HighLevelPriors(proposed,Model)
#     print(paste0("Higher Level Priors = ",HP))
#   }
#   
#   if(HP<=-500) return(-5000)
#   
#   LL<-HP
#   # Load ODD files
#   ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/ODDobjects"),pattern=Model$haz,recursive = T,ignore.case = T)
#   for(i in 1:length(ufiles)){
#     ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/",ufiles[i]))
#     # print(paste0(dir,"IIDIPUS_Input/ODDobjects/",ufiles[i]))
#     ODDy@fIndies<-Model$fIndies
#     ODDy@gmax%<>%as.data.frame.list()
#     # if(max(ODDy@gmax$gmax)<3000) next
#     # Apply DispX
#     tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center,LL = T,Method = AlgoParams),
#                   error=function(e) NA)
#     # if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));next}
#     if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));return(-Inf)}
#     # Per country weighting of the LL, to prevent introducing a country bias
#     cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==ODDy@gmax$iso3[1]]
#     # print(paste0(ufiles[i]," Disp LL = ",mean(tLL,na.rm = T)))
#     # if(expLL) LL<-LL+cWeight*max(log(mean(exp(tLL),na.rm=T)),AlgoParams$cap,na.rm = T)
#     # else LL<-LL+cWeight*max(mean(tLL,na.rm=T),AlgoParams$cap,na.rm = T)
#     LL<-LL+cWeight*sum(tLL)
#   }
#   rm(ODDy)
#   
#   LL<-3*LL
#   print(paste0("LL Disp = ",LL-HP))
#   sLL<-LL
#   
#   # Load BD files
#   ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/BDobjects"),pattern=Model$haz,recursive = T,ignore.case = T)
#   for(i in 1:length(ufiles)){
#     BDy<-readRDS(paste0(dir,"IIDIPUS_Input/BDobjects/",ufiles[i]))
#     BDy@fIndies<-Model$fIndies
#     # Apply BDX
#     tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,center = Model$center,Method=AlgoParams),
#                   error=function(e) NA)
#     # tLL<-BDX(BD = BDy,Omega = proposed,center = Model$center,Method=AlgoParams)
#     # if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));next}
#     if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));return(-Inf)}
#     # print(paste0(ufiles[i]," BD LL = ",mean(tLL,na.rm = T)))
#     if(expLL) LL<-LL+max(log(mean(exp(tLL),na.rm=T)),AlgoParams$cap,na.rm = T)
#     else LL<-LL+max(mean(tLL,na.rm=T),AlgoParams$cap,na.rm = T)
#     
#     if(is.infinite(LL)) print(mean(tLL))
#   }
#   rm(BDy)
#   
#   print(paste0("LL BD = ",LL-sLL))
#   
#   # Add priors
#   if(!is.null(Model$priors)){
#     LL<-LL+sum(Priors(proposed,Model$priors),na.rm=T)
#   }
#   LL<-LL+HP
#   
#   return(LL)
#   
# }

# }
# vals<-seq(0.0001,0.99999,length.out = 50)
# plot(vals,dbeta(vals,Model$BD_params$Params$destroyed$shape1,
#                      Model$BD_params$Params$destroyed$shape2,log = T),
#      ylim=c(-50,5),col="black",ylab="Log Probability (beta)",xlab="Damage %",
#      type="l",lwd=3)
# lines(vals,dbeta(vals,Model$BD_params$Params$severe$shape1,
#                        Model$BD_params$Params$severe$shape2,log = T),col="red",
#       lwd=2,type="b",pch=0)
# lines(vals,dbeta(vals,Model$BD_params$Params$moderate$shape1,
#                   Model$BD_params$Params$moderate$shape2,log = T),col="orange",
#       lwd=2,type="b",pch=1,lty=2)
# lines(vals,dbeta(vals,Model$BD_params$Params$possible$shape1,
#                   Model$BD_params$Params$possible$shape2,log = T),col="green",
#       lwd=2,type="b",pch=2,lty=3)
# lines(vals,dbeta(vals,Model$BD_params$Params$notaffected$shape1,
#                   Model$BD_params$Params$notaffected$shape2,log = T),col="blue",
#       lwd=2,type="b",pch=3,lty=4)
# legend(x=0.45, y=-20, legend=c("Destroyed","Severe","Moderate","Possible","Unaffected"),
#        col=c("black","red","orange","green","blue"),lty=c(1,1,2,3,4),lwd = 2,pch = c(0,0:3))
