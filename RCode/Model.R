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

#library(DescTools)

Model<-list()

haz<-"EQ"
# haz<-"TC"
Model$haz<-haz
Model$I0 <- 4.3 #intensity threshold below which we assume no damage occurs

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

# parameter transformations: 

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

# Parameter set skeleton: 

Model$skeleton <- list(
  Lambda1=list(nu=NA,kappa=NA#,alpha=NA
               ), 
  Lambda2=list(nu=NA,kappa=NA),
  Lambda3=list(nu=NA,kappa=NA),
  eps=list(local=NA,hazard_mort=NA, hazard_disp=NA, hazard_bd=NA, hazard_cor=NA),
  vuln_coeff=list(PDens=NA, SHDI=NA, GNIc=NA, Vs30=NA, EQFreq=NA, #Mag=NA,
                  FirstHaz=NA, Night=NA, FirstHaz.Night=NA)
)

Model$Priors <- list( #All uniform so currently not included in the acceptance probability. 
  Lambda1=list(nu=list(dist='unif', min=6.5, max=10.5), 
               kappa=list(dist='unif', min=0.25, max=3) #, alpha=list(dist='unif', min=-0.1, max=0.5)
  ), 
  Lambda2=list(nu=list(dist='unif', min=9, max=13.5), 
               kappa=list(dist='unif', min=0.25, max=3)),
  Lambda3=list(nu=list(dist='unif', min=6.5, max=10), 
               kappa=list(dist='unif', min=0.25, max=3)),
  eps=list(local=list(dist='unif', min=0, max=2),
           hazard_mort=list(dist='unif', min=0.2, max=1.5),
           hazard_disp=list(dist='unif', min=0.2, max=1.5),
           hazard_bd=list(dist='unif', min=0.2, max=1.5),
           hazard_cor=list(dist='unif', min=0.2, max=1)),
  vuln_coeff=list(PDens=list(dist='laplace', location=0, scale=0.2), # 0.35
                  SHDI=list(dist='laplace', location=0, scale=0.2),
                  GNIc=list(dist='laplace', location=0, scale=0.2),
                  Vs30=list(dist='laplace', location=0, scale=0.2),
                  EQFreq=list(dist='laplace', location=0, scale=0.2),
                  FirstHaz=list(dist='laplace', location=0, scale=0.2),
                  Night=list(dist='laplace', location=0, scale=0.2),
                  FirstHaz.Night=list(dist='laplace', location=0, scale=0.2))
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

#Set lower and upper bounds for the parameters
Model$par_lb <- c()
Model$par_ub <- c()
  
for (i in 1:length(Model$Priors)){
  if (is.list(Model$Priors[[i]])){
    for (j in 1:length(Model$Priors[[i]])){
      if(Model$Priors[[i]][[j]]$dist == 'unif'){
        Model$par_lb = c(Model$par_lb, Model$Priors[[i]][[j]]$min)
        Model$par_ub = c(Model$par_ub, Model$Priors[[i]][[j]]$max)
      } else if (Model$Priors[[i]][[j]]$dist == 'norm'){
        Model$par_lb = c(Model$par_lb, Model$Priors[[i]][[j]]$mean - 6 * Model$Priors[[i]][[j]]$sd)
        Model$par_ub = c(Model$par_ub, Model$Priors[[i]][[j]]$mean + 6 * Model$Priors[[i]][[j]]$sd)
      } else if (Model$Priors[[i]][[j]]$dist == 'laplace'){
        Model$par_lb = c(Model$par_lb, Model$Priors[[i]][[j]]$location - 15 * Model$Priors[[i]][[j]]$scale)
        Model$par_ub = c(Model$par_ub, Model$Priors[[i]][[j]]$location + 15 * Model$Priors[[i]][[j]]$scale)
      } else {
        stop('Please update Method.R to adjust acceptance probability to account for other priors before continuing.')
      }
    }
  } else {
    if(Model$Priors[[i]]$dist == 'unif'){
      Model$par_lb = c(Model$par_lb, Model$Priors[[i]]$min)
      Model$par_ub = c(Model$par_ub, Model$Priors[[i]]$max)
    } else if (Model$Priors[[i]]$dist == 'norm'){
      Model$par_lb = c(Model$par_lb, Model$Priors[[i]]$mean - 6 * Model$Priors[[i]]$sd)
      Model$par_ub = c(Model$par_ub, Model$Priors[[i]]$mean + 6 * Model$Priors[[i]]$sd)
    } else if (Model$Priors[[i]]$dist == 'laplace'){
      Model$par_lb = c(Model$par_lb, Model$Priors[[i]]$location - 15 * Model$Priors[[i]]$scale)
      Model$par_ub = c(Model$par_ub, Model$Priors[[i]]$location + 15 * Model$Priors[[i]]$scale)
    } else {
      stop('Please update Method.R to adjust acceptance probability to account for other priors before continuing.')
    }
  }
}

# Get the binary regression function (currently not implemented, only have option for normal cdf)
# Model$BinR<- "pnorm" #"weibull" # "gompertz" 

# Implement higher order Bayesian priors?
Model$higherpriors<-TRUE

Model$center<-ExtractCentering(dir,haz,T)

Model$impacts <- list(labels = c('mortality', 'displacement', 'buildDam'), 
                      qualifiers = c('qualifierMort', 'qualifierDisp', 'qualifierBuildDam'),
                      sampled = c('mort_sampled', 'disp_sampled', 'buildDam_sampled'))

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

GetLP<-function(ODD_df,Omega,Params,Sinc,notnans, split_GNI=T){
  #if (split_GNI){return(array(1, dim=c(NROW(ODD@data), 8)))}
  #else {return(rep(1, NROW(ODD@data)))}
  
  LP_ij <- array(NA, dim=NROW(ODD_df))
  
  LP_ij[notnans] <- 0 #Omega$vuln_coeff_adj$Mag * (max(ODD@hazinfo$magnitudes) - Params$center$Mag$mean) / Params$center$Mag$sd # Omega$vuln_coeff_adj$itc #intercept term
  
  #Could perform all centering outside before model fitting? may allow a bit of speedup
  #   - All operations except the final sweep are actually very fast so little gains to be made compared with storing centred and uncentred forms
  #split(seq_along(ODD@data$SHDI), ODD@data$SHDI)
  
  #Population density term:
  LP_ij[notnans] <- LP_ij[notnans] + Omega$vuln_coeff_adj$PDens * ((log(ODD_df$PDens[notnans]+0.1) - Params$center$PDens$mean)/Params$center$PDens$sd)
  LP_ij[notnans] <- LP_ij[notnans] + Omega$vuln_coeff_adj$EQFreq * ((log(ODD_df$EQFreq[notnans]+0.001) - Params$center$EQFreq$mean)/Params$center$EQFreq$sd)
  LP_ij[notnans] <- LP_ij[notnans] + Omega$vuln_coeff_adj$Vs30 * ((log(ODD_df$Vs30[notnans]) - Params$center$Vs30$mean)/Params$center$Vs30$sd)

  for (vuln_term in names(Omega$vuln_coeff_adj)[!(names(Omega$vuln_coeff_adj) %in%  c('itc', 'PDens', 'GNIc', 'EQFreq', 'Mag', 'Vs30', 'FirstHaz', 'Night', 'FirstHaz.Night'))]){
    #All remaining terms except GNIc:
    LP_ij[notnans] <- LP_ij[notnans] + Omega$vuln_coeff_adj[[vuln_term]] * ((ODD_df[notnans, vuln_term] - Params$center[[vuln_term]]$mean)/Params$center[[vuln_term]]$sd)
  }
  
  #GNIc:
  if (split_GNI==F){ #don't split into the eight GNIc deciles:  
    LP_ij[notnans] <-  LP_ij[notnans] + Omega$vuln_coeff_adj$GNIc * log(ODD_df$GNIc[notnans])
    return(LP_ij)
  }
  
  LP_ijs <- array(NA, dim=c(NROW(ODD_df),8))
  
  # Slower:
  # get_GNIc_vuln <- function(ij){
  #   vuln_GNIc_ij <- Omega$vuln_coeff_adj$GNIc * (log(ODD@data$GNIc[ij] * Sinc[Sinc$iso3==ODD@data$ISO3C[ij],]$value * 12.5) - Params$center$GNIc$mean)/Params$center$GNIc$sd
  #   return(vuln_GNIc_ij)
  # }
  # LP_ijs[notnans,] <- sweep(t(vapply(t(notnans), get_GNIc_vuln, numeric(8))), 1, LP_ij[notnans], '+')

  GNIc_vuln <- Omega$vuln_coeff_adj$GNIc * (log(ODD_df$GNIc[notnans] * matrix(Sinc$value[unlist(split(seq_along(Sinc$iso3), Sinc$iso3)[ODD_df$ISO3C[notnans]])], ncol=8, byrow=T) * 12.5)- Params$center$GNIc$mean)/Params$center$GNIc$sd
  LP_ijs[notnans,] <- sweep(GNIc_vuln, 1, LP_ij[notnans], '+')
  
  return(LP_ijs)
}

getLP_event <- function(hazinfo, Omega, Params){
  coeffs <- rep(0, length(hazinfo$eventdates))
  coeffs <- Omega$vuln_coeff_adj$FirstHaz * (hazinfo$first_event - Params$center$FirstHaz$mean)/Params$center$FirstHaz$sd
  hour <- as.numeric(substr(hazinfo$eventtimes, 1, 2))
  night_flag <- ifelse(hour>=22 | hour < 6, 1, 0)
  coeffs <- coeffs + Omega$vuln_coeff_adj$Night * (night_flag - Params$center$Night$mean)/Params$center$Night$sd
  coeffs <- coeffs + Omega$vuln_coeff_adj$FirstHaz.Night * (night_flag*hazinfo$first_event - Params$center$FirstHaz.Night$mean)/Params$center$FirstHaz.Night$sd
  return(coeffs)
}

# GetLP_single is the equivalent of GetLP for a single grid-cell ij 
# Used in higher-level prior to calculate linear predictor for a given set of vulnerability terms
GetLP_single <- function(Omega, center, vuln_terms){
  #return(1)
  LP_ij <- 0 # Omega$vuln_coeff_adj$Mag * (vuln_terms[['Mag']] - center$Mag$mean) / center$Mag$sd #Omega$vuln_coeff_adj$itc 
  
  LP_ij <- LP_ij + Omega$vuln_coeff_adj$PDens * ((log(vuln_terms[['PDens']]+0.1) - center$PDens$mean)/center$PDens$sd)
  LP_ij <- LP_ij + Omega$vuln_coeff_adj$EQFreq * ((log(vuln_terms[['EQFreq']]+0.001) - center$EQFreq$mean)/center$EQFreq$sd)
  LP_ij <- LP_ij + Omega$vuln_coeff_adj$Vs30 * ((log(vuln_terms[['Vs30']]) - center$Vs30$mean)/center$Vs30$sd)
  
  for (vuln_term in names(Omega$vuln_coeff_adj)[!(names(Omega$vuln_coeff_adj) %in%  c('itc', 'PDens', 'GNIc', 'EQFreq', 'Mag', 'Vs30'))]){
    #All remaining terms except GNIc: (this includes the event vulnerabilities e.g. night and first event terms)
    LP_ij <- LP_ij + Omega$vuln_coeff_adj[[vuln_term]] * ((vuln_terms[[vuln_term]] - center[[vuln_term]]$mean)/center[[vuln_term]]$sd)
  }
  
  LP_ij <- LP_ij + Omega$vuln_coeff_adj$GNIc * (log(vuln_terms[['GNIc']]) - center$GNIc$mean)/center$GNIc$sd
  
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
  #h[ind]<- I[ind]-I0
  h[ind] <- I[ind] # exp(Omega$theta$theta1*I[ind]) #I[ind]#
  return(h)
}

# Binomial displacement calculator function
rbiny<-function(size,p) rbinom(n = 1,size,p);
Fbdisp<-function(lPopS,Dprime) mapply(rbiny,lPopS,Dprime);

fBD<-function(nbuildings, D_BD) mapply(rbiny, nbuildings, D_BD)

# Calculate the unscaled damage function
fDamUnscaled<-function(I,Params,Omega){ 
  return(h_0(I,Params$I0,Omega)) #+
     #stochastic(Params$Np,Omega$eps_adj$local)) %>%return()
  #(h_0(I,Params$I0,Omega$theta) + 
  #   stochastic(Params$Np,Omega$eps_adj$local)) %>%return()
}

fDamUnscaled_BD<-function(I,Params,Omega){ 
  return(h_0(I,Params$I0,Omega))
  #(h_0(I,Params$I0,Omega$theta) + 
  #   stochastic(Params$Np,Omega$eps_adj$local)) %>%return()
}

addTransfParams <- function(Omega, I0=Model$I0){
  #Omega$theta$theta1 <- 0.6
  Omega$Lambda1$loc <- Omega$Lambda1$nu #h_0(Omega$Lambda1$nu, I0, Omega)
  Omega$Lambda2$loc <- Omega$Lambda2$nu #h_0(Omega$Lambda2$nu, I0, Omega)
  Omega$Lambda3$loc <- Omega$Lambda3$nu #h_0(Omega$Lambda3$nu, I0, Omega)
  #Omega$Lambda4$loc <- Omega$Lambda4$nu #h_0(Omega$Lambda4$nu, I0, Omega)
  #h_10_minus_h_4.5 = h_0(10, I0, Omega) - h_0(4.5, I0, Omega)
  Omega$Lambda1$scale <- Omega$Lambda1$kappa #h_10_minus_h_4.5 / (6 * Omega$Lambda1$kappa)
  Omega$Lambda2$scale <- Omega$Lambda2$kappa #h_10_minus_h_4.5 / (6 * Omega$Lambda2$kappa)
  Omega$Lambda3$scale <- Omega$Lambda3$kappa #h_10_minus_h_4.5 / (6 * Omega$Lambda3$kappa)
  #Omega$Lambda4$scale <- Omega$Lambda4$kappa #h_10_minus_h_4.5 / (6 * Omega$Lambda4$kappa)
  if (Omega$eps$hazard_cor < 0){Omega$eps$hazard_cor = 0.001}
  Omega$vuln_coeff_adj <- Omega$vuln_coeff #lapply(Omega$vuln_coeff, function(x) x * Omega$Lambda2$scale)
  Omega$eps_adj <- Omega$eps #lapply(Omega$eps, function(x) x * Omega$Lambda2$scale)
  Omega$eps_adj$local <- Omega$eps$local / Omega$eps$hazard_mort # local mort / hazard wide mort = ratio
  return(Omega)
}

# Calculate Mortality and Displacement probabilities from the unscaled damage
D_MortDisp_calc <- function(Damage, Omega, stoch_event=rbind(0,0)){
  D_Mort_and_Disp <- pnorm(Damage + stoch_event[2,], mean = Omega$Lambda1$loc, sd = Omega$Lambda1$scale)
  D_Mort <- pnorm(Damage + stoch_event[1,], mean = Omega$Lambda2$loc, sd = Omega$Lambda2$scale)
  D_Disp <- D_Mort_and_Disp - D_Mort
  D_Disp <- ifelse(D_Disp<0, 10^(-60), D_Disp)
  return(rbind(D_Mort, D_Disp))
}

# Calculate Mortality and Displacement probabilities from the unscaled damage
D_Dam_calc <- function(Damage, Omega, stoch_bd=0){
  D_Dam <- pnorm(Damage + stoch_bd, mean = Omega$Lambda3$loc, sd = Omega$Lambda3$scale)
  return(D_Dam)
}

#when working with buildings, D_Disp is equivalent to D_BuildDam and D_Mort is equivalent to D_BuildDest
#rmultinomy<-function(size, D_Disp, D_Mort, D_Rem) return(round(size * c(D_Disp, D_Mort, D_Rem)))
# rmultinomy<-function(size, D_Disp, D_Mort, D_Rem) rmultinom(n=1, size, c(D_Disp, D_Mort, D_Rem))
# Fbdam<-function(PopRem, D_Disp, D_Mort, D_Rem) mapply(rmultinomy, PopRem, D_Disp, D_Mort, D_Rem)

Fbdam <- function(PopRem, D_Disp, D_Mort, D_Rem) {
   vapply(seq_along(PopRem), function(i) {
     rmultinom(1, PopRem[i], c(D_Disp[i], D_Mort[i], D_Rem[i]))
   }, numeric(3))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Log likelihood, posterior and higher-level prior distribution calculations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# These high-level priors are to include expert opinion in the 
# model calculations, not in the individual parameters themselves (for this, see priors)
Model$HighLevelPriors <-function(Omega,Model,modifier=NULL){
  # We begin by calculating the 1st and 99th percentiles of each vulnerability covariate
  # These values have been hard-coded in for speed, but the code left in comments
  
  ## Calculate the 1st and 99th quantiles of the population density from the training dataset:
  # path<-paste0("/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/")
  # ufiles<-list.files(path=path,pattern=haz,recursive = T,ignore.case = T)
  # ufiles<-ufiles[grepl(ufiles,pattern = haz)]
  # PDens<-c()
  # Vs30 <- c()
  # for(fff in ufiles){
  #   ODDy<-readRDS(paste0(path,fff))
  #   PDens<-append(PDens, ODDy@data$Population[!is.na(ODDy@data$Population)])
  #   #Vs30<-append(Vs30, ODDy@data$Vs30[!is.na(ODDy@data$Vs30)])
  # }
  ## min_PDens <- 0; max_PDens <- 102250.5 # range(PDens)
  min_PDens <- 0; max_PDens <- 10520.95 # quantile(PDens, c(0.01, 0.99))
  
  #GDLdata <- readGlobalDataLab()
  ##min_SHDI <- 0.172; max_SHDI <- 0.989 # range(GDLdata$SHDI, na.rm=T)
  min_SHDI <- 0.275; max_SHDI <- 0.941 #quantile(GDLdata$SHDI, c(0.01,0.99), na.rm=T)
  ##min_GNIc <- exp(5.873416); max_GNIc <- exp(12.253871) #log(range(GDLdata$GNIc, na.rm=T))
  min_GNIc <- exp(6.598482); max_GNIc <- exp(11.114198) #quantile(log(GDLdata$GNIc), c(0.01,0.99), na.rm=T)
  
  #stiff<-raster(paste0(dir,"Hazard_Data/global_vs30_tif/global_vs30.tif"))
  #stiffAgg <- aggregate(stiff, 5) #Need to aggregate a bit otherwise quantile crashes R
  #raster::quantile(stiffAgg, probs=c(0.01,0.99))
  min_Vs30 <- 180; max_Vs30 <- 880 #1st and 99th quantiles from all regions in soil stiffness dataset
  
  #pga <- raster(paste0(dir,"Hazard_Data/gdpga/pga_475y.tif"))
  #pga_vals <- values(pga)
  ##min_EQFreq <- 0; max_EQFreq <- 949.4231 #range(pga_vals)
  
  #pga <- rast(paste0(dir,"Hazard_Data/GEM-GSHM_PGA-475y-rock_v2023/v2023_1_pga_475_rock_3min.tif"))
  #pga_vals <- values(pga)
  min_EQFreq <- 0 ; max_EQFreq <- 0.5613953 #quantile(pga_vals, probs=c(0.01, 0.99), na.rm=T)
  #min_EQFreq <- Model$center$EQFreq$mean; max_EQFreq <- Model$center$EQFreq$mean
  
  
  min_FirstHaz <- 0; max_FirstHaz <- 1 #two possible values
  min_Night <- 0; max_Night <- 1 #two possible values
  min_FirstHaz.Night <- 0; max_FirstHaz.Night <- 1 #two possible values
  
  #should be around -3 and 3: 
  #print((c(min_SHDI,max_SHDI) - Model$center$SHDI$mean)/Model$center$SHDI$sd)
  #print((log(c(min_GNIc,max_GNIc)) - Model$center$GNIc$mean)/Model$center$GNIc$sd)
  #print((log(c(min_EQFreq,max_EQFreq)+1) - Model$center$EQFreq$mean)/Model$center$EQFreq$sd)
  #print((log(c(min_PDens,max_PDens)+0.1) - Model$center$PDens$mean)/Model$center$PDens$sd)
  #print((c(log(min_Vs30),log(max_Vs30)) - Model$center$Vs30$mean)/Model$center$Vs30$sd)
  
  linp_min <- GetLP_single(Omega, Model$center, vuln_terms=list(PDens=ifelse(Omega$vuln_coeff_adj$PDens>0, min_PDens, max_PDens), 
                                                                SHDI=ifelse(Omega$vuln_coeff_adj$SHDI>0, min_SHDI, max_SHDI),
                                                                GNIc=ifelse(Omega$vuln_coeff_adj$GNIc>0, min_GNIc, max_GNIc),
                                                                Vs30=ifelse(Omega$vuln_coeff_adj$Vs30>0, min_Vs30, max_Vs30),
                                                                EQFreq=ifelse(Omega$vuln_coeff_adj$EQFreq>0, min_EQFreq, max_EQFreq),
                                                                Mag=ifelse(Omega$vuln_coeff_adj$Mag>0, min_Mag, max_Mag),
                                                                FirstHaz=ifelse(Omega$vuln_coeff_adj$FirstHaz>0, min_FirstHaz, max_FirstHaz),
                                                                Night=ifelse(Omega$vuln_coeff_adj$Night>0, min_Night, max_Night),
                                                                FirstHaz.Night=ifelse(Omega$vuln_coeff_adj$FirstHaz.Night>0, min_FirstHaz.Night, max_FirstHaz.Night)))
  
  linp_max <- GetLP_single(Omega, Model$center, vuln_terms=list(PDens=ifelse(Omega$vuln_coeff_adj$PDens<0, min_PDens, max_PDens), 
                                                                SHDI=ifelse(Omega$vuln_coeff_adj$SHDI<0, min_SHDI, max_SHDI),
                                                                GNIc=ifelse(Omega$vuln_coeff_adj$GNIc<0, min_GNIc, max_GNIc),
                                                                Vs30=ifelse(Omega$vuln_coeff_adj$Vs30<0, min_Vs30, max_Vs30),
                                                                EQFreq=ifelse(Omega$vuln_coeff_adj$EQFreq<0, min_EQFreq, max_EQFreq),
                                                                Mag=ifelse(Omega$vuln_coeff_adj$Mag<0, min_Mag, max_Mag),
                                                                FirstHaz=ifelse(Omega$vuln_coeff_adj$FirstHaz<0, min_FirstHaz, max_FirstHaz),
                                                                Night=ifelse(Omega$vuln_coeff_adj$Night<0, min_Night, max_Night),
                                                                FirstHaz.Night=ifelse(Omega$vuln_coeff_adj$FirstHaz.Night<0, min_FirstHaz.Night, max_FirstHaz.Night)))
  
  #if(!is.null(modifier)) lp<-exp(as.numeric(unlist(modifier))) else lp<-1.
  lp <- c(0,0)#c(linp_min, linp_max) 
  if(Model$haz=="EQ"){
  
    # Lower and upper bounds on the impacts at I_ij = 4.6, 7, and 9.5
    # in the order (Mort, DispMort, BuildDam),
    # where DispMort is the sum of the probabilities of displacement and mortality
    # and BuildDam is the sum of the probabilities of building damage and destruction.
    
    Upp_bounds_4.6 <- c(10^(-6), 0.01, 0.05)
    Low_bounds_7 <- c(0, 0, 10^(-6))
    Upp_bounds_7 <- c(0.01, 0.2, 0.4)
    Low_bounds_9.5 <- c(10^(-6),0.2,0.3)
    
    HLP_impacts <- function(I_ij, lp, Omega){
      rbind(apply(D_MortDisp_calc(h_0(I_ij, Model$I0, Omega=Omega) + lp, Omega),2,cumsum), 
            D_Dam_calc(h_0(I_ij, Model$I0, Omega=Omega) + lp, Omega))
    }
    
    adder <- 0
    
    #Check upper bounds at I=4.6
    adder <- sum(apply(HLP_impacts(4.6, lp, Omega), 2, function (x) sum(x > Upp_bounds_4.6))) 
    
    #Check lower and upper bounds at I=7
    adder <- adder + sum(apply(HLP_impacts(7, lp, Omega), 2, 
                               function (x) sum(c(x > Upp_bounds_7, x<Low_bounds_7))))
    
    #Check lower bounds at I=9
    adder <- adder + sum(apply(HLP_impacts(9.5, lp, Omega), 2, 
                               function (x) sum(x<Low_bounds_9.5)))
    
    #check that at intensity 8, D_disp > D_mort
    impact_intens_8 <- HLP_impacts(8, lp, Omega)
    adder <- adder + sum(impact_intens_8[1,] > impact_intens_8[2,])
    
    #Print results (helpful for debugging / identifying issue with a proposed Omega)
    # print('Upper bounds at I=4.6:')
    # print(paste('Passed', c('Mort:', 'DispMort:', 'BuildDam:'), !HLP_impacts(4.6, lp[2], Omega)>Upp_bounds_4.6))
    # print('Lower bounds at I=7:')
    # print(paste('Passed', c('Mort:', 'DispMort:', 'BuildDam:'), !HLP_impacts(7, lp[1], Omega)<Low_bounds_7))
    # print('Upper bounds at I=7:')
    # print(paste('Passed', c('Mort:', 'DispMort:', 'BuildDam:'), !HLP_impacts(7, lp[2], Omega)>Upp_bounds_7))
    # print('Lower bounds at I=9.5:')
    # print(paste('Passed', c('Mort:', 'DispMort:', 'BuildDam:'), !HLP_impacts(9.5, lp[1], Omega)<Low_bounds_9.5))
    # print(paste('p(Disp)>p(Mort) at Intensity 8:', all(!(impact_intens_8[1,] > impact_intens_8[2,]))))
    

    return(adder) 
    
  } else if(Model$haz=="TC"){
    # Would need to be updated:
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
}


# -------------------------------------------------------------------------------------------------------------------
# ------------------------------- Pulling out and breaking down the distances ---------------------------------------
# -------------------------------------------------------------------------------------------------------------------

SamplePolyImpact <-function(dir,Model,proposed,AlgoParams, dat='Train', output='SampledAgg'){
  
  # Load ODD files
  folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects/")
  
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  x <- file.info(paste0(folderin,ufiles))
  ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
  
  #For Continuous Ranked Probability Score / Energy Score need multiple samples per particle
  AlgoParams$Np <- AlgoParams$Np * AlgoParams$m_CRPS
  
  tmpFn<-function(filer){
    # Extract the ODD object
    ODDy<-readODD(paste0(folderin,filer))
    
    # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
    #ODDy@fIndies<-Model$fIndies
    
    ODDy@impact%<>%as.data.frame.list()
    
    #remove inferred building damage and displacement observations:
    ODDy@impact <- ODDy@impact[!1:NROW(ODDy@impact) %in% which(ODDy@impact$impact %in% c('buildDam', 'displacement') & ODDy@impact$inferred == T),]
    
    ODDy@impact$event_id = as.numeric(gsub(".*_(-?\\d+)$", "\\1", filer))
   
    #DispX requires event_i if introducing correlation between error terms in subsequent model samples (i.e. for correlated MCMC)
    event_i = ifelse(is.null(proposed$u), NA, which(ufiles==filer))
    
    # Apply DispX
    impact_sample_event <- DispX(ODD = ODDy,Omega = proposed,center = Model$center, Method = AlgoParams, output=output, event_i=event_i)
    
    if (NROW(impact_sample_event[[1]]) > 0){
      train_flag = sub("/.*", "", filer)
      impact_sample_event <- lapply(impact_sample_event, function(x){x$train_flag = train_flag; return(x)})  
    }
    
    return(impact_sample_event) 
  }
  
  if (AlgoParams$AllParallel){
    # When using task parallelisation, put the heaviest files first for optimisation reasons
    
    return(pmap(mclapply(X = ufiles,FUN = tmpFn, mc.cores = AlgoParams$cores, mc.preschedule=F), rbind))
    
  } else { 
    #even when AlgoParams$cores = 1 the above will still work, but this can sometimes be useful for debugging:
    impact_sample_poly <- tmpFn(ufiles[1])
    for (file_i in 2:length(ufiles)){
      impact_sample_poly <- pmap(list(impact_sample_poly, tmpFn(ufiles[file_i])), rbind)
    }
  }
}


SamplePointImpact <- function(dir,Model,proposed,AlgoParams, dat='Train', output='LL'){
  #NOT CURRENTLY IMPLEMENTED:
  
  # Load BD files
  folderin<-paste0(dir,AlgoParams$input_folder,"BDobjects/")
  
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
      
      if(is.null(BDy)){return(rep(0, AlgoParams$Np))}
      if(nrow(BDy@data)==0){return(rep(0, AlgoParams$Np))}
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      #BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, output=output),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",filer));return(-Inf)}
      if (length(tLL > 1)){ return(tLL)}
      return(tLL) 
      # # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      # maxLL<-max(tLL,na.rm = T)
      # # Return the average log-likelihood
      # cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==BDy$ISO3C[which(!is.na(BDy$ISO3C))[1]]] #looseend
      # if(is.na(cWeight)){return(-Inf)}                               
      # if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)) 
      # else return(cWeight*mean(tLL,na.rm=T))
      
    }
    return(do.call(rbind, mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))) #SMC-CHANGE #d-change
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

SampleImpact <- function(dir,Model,proposed,AlgoParams, dat='Train', output='SampledAgg'){

  impact_sample_poly<-SamplePolyImpact(dir,Model,proposed,AlgoParams, dat=dat, output=output)
  if (AlgoParams$BD_weight > 0){
    impact_sample_point<- SamplePointImpact(dir, Model, proposed, AlgoParams, expLL=T, dat=dat)
  } else {
    impact_sample_point = NULL
  }
  return(list(poly=impact_sample_poly, point=impact_sample_point)) 
}

sample_quant <- function(x){
  indexes <- which(sort(as.numeric(x))==as.numeric(x[1]))
  return(ifelse(length(indexes)==1, indexes, sample(indexes,1)))
}

get_average_rank_single <- function(df, log=F){
 pre_ranks <- apply(apply(df, 1, rank), 1, mean)
 return(rank(pre_ranks,  ties.method ='random')[1])
}


#To use spantree function:
get_mst_rank_single <- function(df, log=F, noise=NA){
  # Calculates the minimum spanning tree pre-rank of the observation within the set of simulations.
  # Can correlate the noise between current step and proposal to reduce effect of randomness in how ties are broken.
  pre_ranks <- c()
  for (j in 1:NCOL(df)){
    if (NROW(df)>1){
      pre_ranks <- c(pre_ranks, sum(spantree(dist(t(df[,-j])))$dist))
    } else {
      pre_ranks <- c(pre_ranks, sum(spantree(dist(df[,-j]))$dist))
    }
  }
  if (all(!is.na(noise))){
    pre_ranks = pre_ranks + noise * min(abs(diff(unique(pre_ranks))))/2
  }
  return(rank(pre_ranks, ties.method ='random')[1])
}

get_banddepth_rank_single <- function(mat, noise = NA){
  # Calculates the band depth pre-rank of the observation within the set of simulations.
  # Can correlate the noise between current step and proposal to reduce effect of randomness in how ties are broken.
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
  return(rank(pre_ranks,  ties.method ='random')[1])
}

CalcDistPoly_EnergyScore <- function(impact_sample_poly, AlgoParams, corr_noise = NA, corr_noise2 = NA){
  #Calculate the mean energy score from impact_sample
  # Note that dist_poly has length 7, however, only the first element is currently used when calculating 
  # the distance. It has been kept with length greater than 1 in case we want to store up to 6 other values 
  # (e.g. the Anderson Darling Statistic of pre-ranks using the Minimum Spanning Tree), but these other
  # elements are not actually currently used when calculating the distance. 
  
  observed <- impact_sample_poly[[1]]$observed

  dist_poly <- array(NA, dim=c(AlgoParams$Np,7))
  impact_type <- impact_sample_poly[[1]]$impact
  impact_weightings <- unlist(AlgoParams$impact_weights[impact_type])
  event_id <- impact_sample_poly[[1]]$event_id
  grouped_events <- split(seq_along(event_id), event_id)
  
  for(n in 1:AlgoParams$Np){

    samples_allocated <- ((n-1)*AlgoParams$m_CRPS+1):(n*AlgoParams$m_CRPS)
    samples_combined <- sapply(impact_sample_poly[samples_allocated], function(x){x$sampled}) #doesn't work if samples_allocated is length 1
    
    dist_poly[n,1] <- 0
    
    es_store <- c()
    vs_store <- c()
    pre_ranks <- c()
    #pre_ranks_average <- c() #can also assess quantiles for uniformity based on the average pre-rank function
    #pre_ranks_mst <- c() #can also assess quantiles for uniformity based on the minimum spanning tree pre-rank function
    
    for (i in 1:length(grouped_events)){
      #For each event, compute the energy score of the observed data vs the 'prediction' (simulated data)
      #Each impact type is weighted differently, simply multiplying the observation and the simulations by this weight performs the weighting
      obs <- log(observed[grouped_events[[i]]]+AlgoParams$log_offset) *impact_weightings[grouped_events[[i]]]
      sims <- log(samples_combined[grouped_events[[i]],, drop=F]+AlgoParams$log_offset) * impact_weightings[grouped_events[[i]]]
      
      obs_rh <- log(observed[intersect(grouped_events[[i]], which(impact_type!='buildDam'))]+AlgoParams$log_offset) *impact_weightings[intersect(grouped_events[[i]], which(impact_type!='buildDam'))]
      sims_rh <- log(samples_combined[intersect(grouped_events[[i]], which(impact_type!='buildDam')),,drop=F]+AlgoParams$log_offset) * impact_weightings[intersect(grouped_events[[i]], which(impact_type!='buildDam'))]
      # get_mst_rank_single(cbind(obs_mst, sims_mst))
      
      #pre_ranks_mst <- c(pre_ranks_mst, get_mst_rank_single(cbind(obs_mst, sims_mst)))
      pre_ranks <- c(pre_ranks, get_banddepth_rank_single(cbind(obs_rh, sims_rh)))#, noise=corr_noise[i,]))
      #pre_ranks <- c(pre_ranks, get_mst_rank_single(cbind(obs_rh, sims_rh), noise=corr_noise[i,]))

      
      if (length(grouped_events[[i]])==1){
        #crps is the univariate case of the energy score, so use when dimension is 1. 
        es_store<- c(es_store, crps_sample(obs, sims))
        vs_store<- c(vs_store, 0)
        next
      } 
      es_store<- c(es_store, es_sample(obs, sims))
      w_vs <- sqrt(as.numeric(impact_weightings[grouped_events[[i]]]) %*% t(as.numeric(impact_weightings[grouped_events[[i]]]))) / (length(grouped_events[[i]])^2)
      if(vs_sample(obs, sims, w_vs = w_vs) > 5000){stop()}
      vs_store <- c(vs_store, vs_sample(obs, sims, w_vs = w_vs))
      #pre_ranks_average <- c(pre_ranks_average, get_average_rank_single(cbind(obs, sims)))
      #pre_ranks_mst <- c(pre_ranks_mst, get_mst_rank_single(cbind(obs,sims)))

    }
    
    ranks_std <- (pre_ranks-runif(length(pre_ranks),0,1))/(AlgoParams$m_CRPS + 1)
    #ranks_std <- (pre_ranks-corr_noise2)/(AlgoParams$m_CRPS + 1) #can also correlate the noise added here
  
    dist_poly[n,1] <- mean(es_store) #mean(crps_store[which(impact_type=='mortality')]) * unlist(AlgoParams$impact_weights['mortality'])
    dist_poly[n,2] <- AlgoParams$W_rankhist * (AndersonDarlingTest(ranks_std, null='punif')$statistic) #3*mean(vs_store) #mean(vs_store)/5
    dist_poly[n,3] <- 0
    
    ## Ways to assess uniformity of quantiles:
    #ranks_std_average <- (pre_ranks_average-runif(length(pre_ranks_average),0,1))/(AlgoParams$m_CRPS + 1)
    #ranks_std_mst <- (pre_ranks_mst-runif(length(pre_ranks_mst),0,1))/(AlgoParams$m_CRPS + 1)
    #bin_counts <- table(cut(ranks_std_mst, seq(0, 1, 0.1), include.lowest = TRUE, right = FALSE))
    #sum((bin_counts-length(ranks_std_mst)/10)^2)/(length(ranks_std_mst)/10)
    dist_poly[n,4] <- 0 #sum((bin_counts-length(ranks_std_mst)/10)^2)/(length(ranks_std_mst)/10) #0 #ks.test(ranks_std_average, y='punif')$p.value
    dist_poly[n,5] <- 0 
    dist_poly[n,6] <- 0 #0.2*(1 - AndersonDarlingTest(ranks_std_average, null='punif')$p.value) 
    dist_poly[n,7] <- 0 #0.2*(1 - AndersonDarlingTest(ranks_std_mst, null='punif')$p.value)
    
  }
  
  return(dist_poly)
  
}

CalcDistPoint <- function(impact_sample_point, AlgoParams){
  # NOT CURRENTLY IMPLEMENTED DUE TO ISSUES WITH SATELLITE DATA
  
  #is there a way to do this using a scoring rule as well? : 
  # sumPointDat_dists <- function(PointDat_p){
  #   Dist_0.5 <- which(names(PointDat_p) %in% c('N12', 'N21', 'N23', 'N32'))
  #   Dist_1 <- which(names(PointDat_p) %in% c('N13', 'N31'))
  #   return(0.5*sum(PointDat_p[Dist_0.5])+sum(PointDat_p[Dist_1]))
  # }
  
  #F1 score
  sumPointDat_dists <- function(PointDat_p){
    F1 = c()
    for (i in 1:(NROW(PointDat_p)/4)){
      event_dat <- PointDat_p[((i-1)*4+1):(i*4)]
      F1 = c(F1, (2 * event_dat[which(names(event_dat)=='N22')]) / (2 * event_dat[which(names(event_dat)=='N22')] + 2 * event_dat[which(names(event_dat)=='N12')] + 2 * event_dat[which(names(event_dat)=='N21')]))
    }
    return(sum(1-F1)*100) # LOOSEEND: don't have solid justification for this choice. 
  }
  
  if (length(impact_sample_point) > 0){
    dist_point <- apply(impact_sample_point, sumPointDat_dists)
  } else {
    dist_point <- 0
  }
  
  #print(paste0('Dist_agg: ',rowSums(dist_poly), ' Dist_agg_means: ', sum(dist_poly_means),' Dist_sat: ', dist_point))
  #dist_tot <- dist_poly + dist_poly_means + dist_point
  
  dist_tot = dist_poly[,which(names(dist_poly_raw[1:(length(dist_poly_raw)/AlgoParams$Np)]) %in% c('displacement', 'mortality'))] + rep(dist_poly_means, each=AlgoParams$Np)
  
  
  return(dist_tot)
}

CalcDist <- function(impact_sample, AlgoParams, corr_noise = NA, corr_noise2 = NA){
  if (AlgoParams$kernel == 'energy_score'){
    dist_poly <- CalcDistPoly_EnergyScore(impact_sample$poly, AlgoParams, corr_noise = corr_noise, corr_noise2 = corr_noise2)
    dist_point <- 0 #CalcDistPoint(impact_sample$point, AlgoParams)
    dist_tot <- cbind(dist_poly, dist_point)
  } else {
    stop('Currently no other distance functions implemented except the energy score')
  }
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
#   Lambda1=list(nu='ab_bounded_acc', omega='v'),
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

# obs <- c(100, 100, 100)
# sim1 <- rnorm(10000, 100, 0.00001)
# sim2 <- rnorm(10000, 100, 28)
# sim3 <- rnorm(10000,100, 0.000001)
# w_2 <- 0.125
# w_vs <- cbind(c(1,w_2,1),c(w_2,w_2,w_2), c(1,w_2,1))
# vs_sample(obs, rbind(sim1, sim2,sim3), w_vs=w_vs, p=0.5)
# 
# #Disp:
# obs <- c(100, 100, 100)
# sim1 <- rnorm(10000, 100, 7)
# sim2 <- rnorm(10000, 100, 7)
# sim3 <- rnorm(10000,100, 7)
# w_2 <- 1/7
# w_vs <- cbind(c(w_2,w_2,w_2),c(w_2,w_2,w_2), c(w_2,w_2,w_2))
# vs_sample(obs, rbind(sim1, sim2,sim3), w_vs=w_vs, p=0.5)
# 
# #Mort
# obs <- c(100, 100, 100)
# sim1 <- rnorm(10000, 100, 1)
# sim2 <- rnorm(10000, 100, 1)
# sim3 <- rnorm(10000,100, 1)
# w_2 <- 1
# w_vs <- cbind(c(w_2 ,w_2,w_2 ),c(w_2,w_2,w_2), c(w_2 ,w_2,w_2 ))
# vs_sample(obs, rbind(sim1, sim2,sim3), w_vs=w_vs, p=0.5)
# 
# #DispMort:
# obs <- c(100, 100, 100)
# sim1 <- rnorm(10000, 100, 1)
# sim2 <- rnorm(10000, 100, 3)
# sim3 <- rnorm(10000,100, 5)
# w_1 = 1; w_2 = 1/3; w_3 = 1/5
# w_12 <- sqrt(w_1*w_2); w_23=sqrt(w_2*w_3); w_13=sqrt(w_1*w_3)
# w_vs <- cbind(c(w_1 ,w_12, w_13),c(w_12,w_2,w_23), c(w_13,w_23,w_3))
# vs_sample(obs, rbind(sim1, sim2,sim3), w_vs=w_vs, p=0.5)
# 
# 
# xgr <- c(1,3,5,8,10)
# plot(xgr, c(1, 1/5, 1/13, 1/33, 1/50), type='l')
# lines(xgr,1/xgr, col='red')

