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

# Get linear predictor variables (currently the same for both building damage and human displacement)
if(haz%in%c("EQ","TC","FL")) {
  # Form a list of the different components used in the linear predictor
  INFORM_vars<-c("CC.INS.GOV.GE", # Government Effectiveness
                 "VU.SEV.AD", # Economic Dependency Vulnerability
                 "CC.INS.DRR", # Disaster Risk Reduction
                 "VU.SEV.PD", # Multi-dimensional Poverty
                 "CC.INF.PHY" # Physical Infrastructure
  )
  INFORM_vars%<>%c(grep(haz,c("HA.NAT.EQ","HA.NAT.TC","HA.NAT.FL"),value = T))
  # MUST BE SAME LENGTH AND CORRESPOND EXACTLY TO INFORM_vars + HAZ.NAT.EQ/TC/FL + Sinc/dollar,
  fIndies<-list(CC.INS.GOV.GE=returnX, # Government Effectiveness
                VU.SEV.AD=returnX, # Economic Dependency Vulnerability
                CC.INS.DRR=returnX, # Disaster Risk Reduction
                VU.SEV.PD=returnX, # Multi-dimensional Poverty
                CC.INF.PHY=returnX) # Physical Infrastructure
                # dollar=returnX, # IncomeDistribution*GDP
                # Pdens=returnX) # Population Density
  if(haz=="EQ") {fIndies%<>%c(HA.NAT.EQ=function(x) x*x) # Hazard Exposure)
  } else if (haz=="TC") {fIndies%<>%c(HA.NAT.TC=function(x) x*x) # Hazard Exposure)
  } else if (haz=="FL") {fIndies%<>%c(HA.NAT.EQ=function(x) x*x) }# Hazard Exposure)
  WID_perc<-   c("p10p100", # top 90% share of Income Distribution
                 "p20p100", # top 80% share of Income Distribution
                 "p30p100", # top 70% share of Income Distribution
                 "p40p100", # top 60% share of Income Distribution
                 "p50p100", # top 50% share of Income Distribution
                 "p60p100", # top 40% share of Income Distribution
                 "p70p100", # top 30% share of Income Distribution
                 "p80p100", # top 20% share of Income Distribution
                 "p90p100" # top 10% share of Income Distribution
  )
  
  Model%<>%c(list(INFORM_vars=INFORM_vars,WID_perc=WID_perc,fIndies=fIndies))
  
}

Model$modifiers$INFORM<- T
Model$modifiers$WB<-     T
Model$modifiers$WID<-    T

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

# Link functions (MUST BE SAME LENGTH AS OMEGA)
Model$links<-list(
  Lambda1=list(nu='ab_bounded',omega='ab_bounded'),
  Lambda2=list(nu='ab_bounded',omega='ab_bounded'),
  Lambda3=list(nu='ab_bounded',omega='ab_bounded'),
  Lambda4=list(nu='ab_bounded',omega='ab_bounded'), 
  # beta=list(xxx='exp',CC.INS.GOV.GE='exp',VU.SEV.AD='exp',CC.INS.DRR='exp',VU.SEV.PD='exp',CC.INF.PHY='exp'),
  Pdens=list(M='ab_bounded',k='ab_bounded'),
  dollar=list(M='ab_bounded',k='ab_bounded'),
  theta=list(e='ab_bounded'), #list(e=0.25),
  # rho=list(A='exp',H='exp'),
  eps=list(eps='ab_bounded')#,#,xi='exp')
  # mu=list(muplus='exp',muminus='exp',sigplus='exp',sigminus='exp')
  #lp=list(A='ab_bounded', B='ab_bounded', C='ab_bounded', D='ab_bounded', E='ab_bounded') # add linear predictor terms for testing
)

# Model$links<-list(
#   Lambda1=list(nu='ab_bounded',omega='ab_bounded'),
#   Lambda2=list(nu='ab_bounded',omega='ab_bounded'),
#   Lambda3=list(nu='ab_bounded',omega='ab_bounded'),
#   zeta=list(k='ab_bounded',lambda='ab_bounded'), # zeta=list(k=2.5,lambda=1.6),
#   #zeta1=list(k='exp',lambda='exp'),
#   #zeta2=list(k='exp',lambda='exp'),
#   #zeta3=list(k='exp',lambda='exp'),
#   # beta=list(xxx='exp',CC.INS.GOV.GE='exp',VU.SEV.AD='exp',CC.INS.DRR='exp',VU.SEV.PD='exp',CC.INF.PHY='exp'),
#   Pdens=list(M='ab_bounded',k='ab_bounded'),
#   dollar=list(M='ab_bounded',k='ab_bounded'),
#   theta=list(e='ab_bounded'), #list(e=0.25),
#   # rho=list(A='exp',H='exp'),
#   eps=list(eps='ab_bounded')#,xi='exp')
#   # mu=list(muplus='exp',muminus='exp',sigplus='exp',sigminus='exp')
# )

# And to go the other way....
# Model$unlinks<-list(
#   Lambda1=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
#   Lambda2=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
#   Lambda3=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
#   zeta=list(k='ab_bounded_inv',lambda='ab_bounded_inv'), # zeta=list(k=2.5,lambda=1.6),
#   #zeta1=list(k='log',lambda='log'),
#   #zeta2=list(k='log',lambda='log'),
#   #zeta3=list(k='log',lambda='log'),
#   # beta=list(xxx='log',CC.INS.GOV.GE='log',VU.SEV.AD='log',CC.INS.DRR='log',VU.SEV.PD='log',CC.INF.PHY='log'),
#   Pdens=list(M='ab_bounded_inv',k='ab_bounded_inv'),
#   dollar=list(M='ab_bounded_inv',k='ab_bounded_inv'),
#   theta=list(e='ab_bounded_inv'), #list(e=0.25),
#   # rho=list(A='log',H='log'),
#   eps=list(eps='ab_bounded_inv')#,xi='log')
#   # mu=list(muplus='exp',muminus='exp',sigplus='exp',sigminus='exp')
# )

Model$unlinks<-list(
  Lambda1=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
  Lambda2=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
  Lambda3=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
  Lambda4=list(nu='ab_bounded_inv',omega='ab_bounded_inv'),
  #zeta1=list(k='log',lambda='log'),
  #zeta2=list(k='log',lambda='log'),
  #zeta3=list(k='log',lambda='log'),
  # beta=list(xxx='log',CC.INS.GOV.GE='log',VU.SEV.AD='log',CC.INS.DRR='log',VU.SEV.PD='log',CC.INF.PHY='log'),
  Pdens=list(M='ab_bounded_inv',k='ab_bounded_inv'),
  dollar=list(M='ab_bounded_inv',k='ab_bounded_inv'),
  theta=list(e='ab_bounded_inv'), #list(e=0.25),
  # rho=list(A='log',H='log'),
  eps=list(eps='ab_bounded_inv')#,#,xi='log')
  # mu=list(muplus='exp',muminus='exp',sigplus='exp',sigminus='exp')
  #lp=list(A='ab_bounded_inv', B='ab_bounded_inv', C='ab_bounded_inv', D='ab_bounded_inv', E='ab_bounded_inv')
)

Model$acceptTrans <- list(
  Lambda1=list(nu='ab_bounded_acc', omega='ab_bounded_acc'),
  Lambda2=list(nu='ab_bounded_acc', omega='ab_bounded_acc'),
  Lambda3=list(nu='ab_bounded_acc', omega='ab_bounded_acc'),
  Lambda4=list(nu='ab_bounded_acc', omega='ab_bounded_acc'),
  #zeta1=list(k=NA,lambda=NA),
  #zeta2=list(k=NA,lambda=NA),
  #zeta3=list(k=NA,lambda=NA),
  # beta=list(xxx=NA,CC.INS.GOV.GE=NA,VU.SEV.AD=NA,CC.INS.DRR=NA,VU.SEV.PD=NA,CC.INF.PHY=NA),
  Pdens=list(M='ab_bounded_acc',k='ab_bounded_acc'),
  dollar=list(M='ab_bounded_acc',k='ab_bounded_acc'),
  theta=list(e='ab_bounded_acc'), #list(e=0.25),
  # rho=list(A=NA,H=NA),
  eps=list(eps='ab_bounded_acc')#,#,xi=NA)
  # mu=list(muplus=NA,muminus=NA,sigplus=NA,sigminus=NA)
  #lp=list(A='ab_bounded_acc', B='ab_bounded_acc', C='ab_bounded_acc', D='ab_bounded_acc', E='ab_bounded_acc')
)

# names(Model$unlinks$beta)[1]<-paste0("HA.NAT.",haz)
# Skeleton
Model$skeleton <- list(
  Lambda1=list(nu=NA,omega=NA),
  Lambda2=list(nu=NA,omega=NA),
  Lambda3=list(nu=NA,omega=NA),
  Lambda4=list(nu=NA,omega=NA),
  #zeta1=list(k=NA,lambda=NA),
  #zeta2=list(k=NA,lambda=NA),
  #zeta3=list(k=NA,lambda=NA),
  # beta=list(xxx=NA,CC.INS.GOV.GE=NA,VU.SEV.AD=NA,CC.INS.DRR=NA,VU.SEV.PD=NA,CC.INF.PHY=NA),
  Pdens=list(M=NA,k=NA),
  dollar=list(M=NA,k=NA),
  theta=list(e=NA), #list(e=0.25),
  # rho=list(A=NA,H=NA),
  eps=list(eps=NA)#,#,xi=NA)
  # mu=list(muplus=NA,muminus=NA,sigplus=NA,sigminus=NA)
  #lp=list(A=NA, B=NA, C=NA, D=NA, E=NA)
)
# names(Model$skeleton$beta)[1]<-paste0("HA.NAT.",haz)

#Set lower and upper bounds for the parameters
Model$par_lb <- c(6,  
                  0,
                  7.5, 
                  0,
                  6,  
                  0,  
                  0,  
                  0,  
                  0,  #PDens M
                  0,  #PDens k
                  -3,  #dollar M
                  0,   #dollar k
                  0.1, #theta_e
                  0)#, #epsilon
                  #-1, -1, -1, -1, -1) 

Model$par_ub <- c(9.5, 
                  6, 
                  10.5, 
                  6, 
                  9.5, 
                  6, 
                  10, 
                  6, 
                  3, #PDens M
                  10, #PDens k 
                  0,  #dollar M
                  10, #dollar k
                  1, #theta_e
                  1)#, #epsilon
                  #1, 1, 1, 1, 1) 


# Get the binary regression function
Model$BinR<-"weibull" # "gompertz"

# Implement higher order Bayesian priors?
# Model$higherpriors<-TRUE

# Get the building damage beta distribution parameters
Model$BD_params<-list(functions=list(
  destroyed="dbeta",
  severe="dbeta",
  moderate="dbeta",
  possible="dbeta",
  notaffected="dbeta"
  # damaged="alldamaged" # damaged is a combination of all categories except 'notaffected'
),
Params=list(
  destroyed=list(shape1=40,shape2=2),
  severe=list(shape1=40,shape2=25),
  moderate=list(shape1=15,shape2=30),
  possible=list(shape1=3,shape2=100),
  notaffected=list(shape1=0.05,shape2=100)
  # damaged=list() # damaged is a combination of all categories except 'notaffected'
)
)

Model$center<-ExtractCentering(dir,haz,T)

Model$impacts <- list(labels = c('mortality', 'displacement', 'buildDam', 'buildDest'), 
                      qualifiers = c('qualifierMort', 'qualifierDisp', 'qualifierBuildDam', 'qualifierBuildDest'),
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

# Linear predictor function - ugly, but very fast!
llinpred<-function(Params,beta,center,func) { # Note all of these are vectors, not single values
  # if(is.null(func)) func<-replicate(length(beta),function(x) x)
  linp<-list()
  for(j in 1:length(Params)) {
    vars<-Params[[j]]$var
    # This part uses the name of the variable to search for the following:
    # 1) link function - func, 2) centering value - center, 3) linear predictor parameterisation - beta
    linp[[names(Params[j])]]<-prod(exp(vapply(1:nrow(vars),
                                function(i) beta[[vars$variable[i]]]*(do.call(func[[vars$variable[i]]],
                                                                              list(vars$value[i]-center[[vars$variable[i]]] ))),
                                FUN.VALUE = numeric(1))))
  }
  return(linp)
}

# Quicker version of llinpred for when only one vars$variable exists
dlinpred<-function(vars,beta,center,func) { # Note all of these are vectors, not single values
  # Quicker version of llinpred for when only one vars$variable exists
  return(exp(beta[[vars$variable[1]]]*(do.call(func[[vars$variable[1]]],list(vars$value-center[[vars$variable[1]]] )))))
}

WeibullScaleFromShape<-function(shape,center){
  # Scale parameter
  return(center/(log(2)^(1/shape)))
}
# Linear predictor for subnational variables
locpred<-function(x,params){
  exp(params$M*(pweibull(x,shape=params$shape,scale=params$scale)-0.5))
}

GDPlinp<-function(ODD,Sinc,beta,center,notnans){
  iGDP<-as.numeric(factor(ODD@data$GDP,levels=unique(ODD@data$GDP)))
  dGDP<-data.frame(ind=iGDP[notnans],GDP=ODD@data$GDP[notnans],iso=ODD@data$ISO3C[notnans])
  dGDP%<>%group_by(ind)%>%summarise(value=log(unique(GDP)*(Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]),"value"]/
                                      Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]) & Sinc$variable=="p50p100","value"])),
                              income=Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]),"variable"],
                              variable="dollar",
                              iso3c=unique(dGDP$iso[dGDP$ind==unique(ind)]),.groups = 'drop_last')
  dGDP$linp<-rep(NA_real_,nrow(dGDP))
  # Autocalculate the Weibull scale parameter using the shape and centering values
  tp<-list(shape=beta$k,M=beta$M,
        scale=WeibullScaleFromShape(shape = beta$k,center = center$dollar))
  dGDP$linp[which(!is.na(dGDP$ind))]<-locpred(dGDP$value[!is.na(dGDP$ind)],tp)
  # dGDP$linp[which(!is.na(dGDP$ind))]<-dlinpred(dGDP[!is.na(dGDP$ind),],beta,center,fIndies)
  
  return(list(dGDP=dplyr::select(dGDP,c(ind,income,linp)),iGDP=iGDP))
}

Plinpred<-function(Pdens,beta,center,notnans){
  
  Plinp<-rep(NA_real_,length(Pdens))
  # Plinp[notnans]<-dlinpred(data.frame(variable=rep("Pdens",length(Pdens[notnans])),value=log(Pdens[notnans]+1)),beta,center,fIndies)
  # Autocalculate the Weibull scale parameter using the shape and centering values
  tp<-list(shape=beta$k,M=beta$M,
        scale=WeibullScaleFromShape(shape = beta$k,center = center$Pdens))
  Plinp[notnans]<-locpred(log(Pdens[notnans]+1),tp)
  return(Plinp)
  
}

GetLP<-function(ODD,Omega,Params,Sinc,notnans){
  
  # Calculate national vulnerability or apply modifier parameter
  if(!is.null(tryCatch(ODD@modifier,error=function(e) NULL))){ # Modifier parameter
    linp<-as.list(exp(as.numeric(unlist(ODD@modifier))))
    names(linp)<-names(ODD@modifier)
  } else { # National vulnerability
    linp<-rep(list(1.),length(unique(ODD@cIndies$iso3)))
    
    #testing convergence using dummy linear predictor terms:
    #for (iso3_filter in unique(ODD@cIndies$iso3)){
    #  linp_values <- (ODD@cIndies %>% dplyr::filter(iso3==iso3_filter))$value[1:5]
    #  linp[i] <- sum(exp((linp_values-5) * as.numeric(Omega$lp)))
    #}
    #names(linp)<-unique(ODD@cIndies$iso3)
    #llinpred(Params[c(unique(ODD@cIndies$iso3))],Omega$beta,Params$center,Params$fIndies)
  }
  # Calculate possible dollar (GDP*income_dist) linear predictor values
  GDP<-GDPlinp(ODD,Sinc,Omega$dollar,Params$center,notnans)
  # Calculate population density linear predictor values
  Plinp<-Plinpred(ODD@data$Population,
                  Omega$Pdens,
                  Params$center,
                  notnans)
  return(list(linp=linp,dGDP=GDP$dGDP,iGDP=GDP$iGDP,Plinp=Plinp))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Generalised Linear Models for building damage & displacement calcs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Binary regression function
if(Model$BinR=="weibull") {
  BinR<-function(x,zeta) pweibull(x,shape=zeta$k,scale=zeta$lambda)
} else if(Model$BinR=="gompertz") {
  BinR<-function(x,zeta) flexsurv::pgompertz(x,shape=zeta$varrho,rate=zeta$eta)
} else stop("Incorrect binary regression function name, try e.g. 'weibull'")
# Stochastic damage function process
stochastic<-function(n,eps){
  return(rgammaM(n = n,mu = 1, sig_percent = eps$eps ))
}

# Baseline hazard function h_0
h_0<-function(I,I0,theta){
  ind<-I>I0
  h<-rep(0,length(I))
  h[ind]<- exp( theta$e * (I[ind]-I0) ) -1
  return(h)
}

# Calculate the unscaled damage function
fDamUnscaled<-function(I,Params,Omega){ 
  (h_0(I,Params$I0,Omega$theta) * 
     stochastic(Params$Np,Omega$eps)) %>%return()
}

addTransfParams <- function(Omega, I0=4.5){
  Omega$Lambda1$loc <- h_0(Omega$Lambda1$nu, I0, Omega$theta)
  Omega$Lambda2$loc <- h_0(Omega$Lambda2$nu, I0, Omega$theta)
  Omega$Lambda3$loc <- h_0(Omega$Lambda3$nu, I0, Omega$theta)
  Omega$Lambda4$loc <- h_0(Omega$Lambda4$nu, I0, Omega$theta)
  Omega$Lambda1$sd <- exp(Omega$theta$e*Omega$Lambda1$omega)/6
  Omega$Lambda2$sd <- exp(Omega$theta$e*Omega$Lambda2$omega)/6
  Omega$Lambda3$sd <- exp(Omega$theta$e*Omega$Lambda3$omega)/6
  Omega$Lambda4$sd <- exp(Omega$theta$e*Omega$Lambda4$omega)/6
  return(Omega)
}

# Calculate Mortality and Displacement probabilities from the unscaled damage
D_MortDisp_calc <- function(Damage, Omega){
  D_Mort_and_Disp <- pnorm(Damage, mean = Omega$Lambda1$loc, sd = Omega$Lambda1$sd)
  D_Mort <- pnorm(Damage, mean = Omega$Lambda2$loc, sd = Omega$Lambda2$sd)
  D_Disp <- D_Mort_and_Disp - D_Mort
  D_Disp <- ifelse(D_Disp<0, 0, D_Disp)
  return(rbind(D_Mort, D_Disp))
}

# Calculate Mortality and Displacement probabilities from the unscaled damage
D_DestDam_calc <- function(Damage, Omega, first_haz=T, DestDam_modifiers = c(1,1,1), ind_dam=c()){
  if(first_haz == T){
    D_Dest_and_Dam <- pnorm(Damage, mean = Omega$Lambda3$loc, sd = Omega$Lambda3$sd)
    D_Dest <- pnorm(Damage, mean = Omega$Lambda4$loc, sd = Omega$Lambda4$sd)
    D_Dam <- D_Dest_and_Dam - D_Dest
    D_Dam <- ifelse(D_Dam<0, 0, D_Dam)
  } else {
    if (length(ind_dam>0)){
      D_Dest_and_Dam <- D_Dam <- D_Dest <- rep(0, length(Damage))
      D_Dest_and_Dam[-ind_dam] <- pnorm(Damage[-ind_dam], mean = Omega$Lambda3$loc, sd = Omega$Lambda3$sd) ^ DestDam_modifiers[1]
      D_Dest[-ind_dam] <- pnorm(Damage[-ind_dam], mean = Omega$Lambda4$loc, sd = Omega$Lambda4$sd) ^ DestDam_modifiers[2]
      D_Dest_and_Dam[ind_dam] <- 1
      D_Dest[ind_dam] <- pnorm(Damage[ind_dam], mean = Omega$Lambda4$loc, sd = Omega$Lambda4$sd) ^ DestDam_modifiers[3]
      D_Dam <- D_Dest_and_Dam - D_Dest
      D_Dam <- ifelse(D_Dam<0, 0, D_Dam)
    } else {
      D_Dest_and_Dam <- pnorm(Damage, mean = Omega$Lambda3$loc, sd = Omega$Lambda3$sd) ^ DestDam_modifiers[1]
      D_Dest <- pnorm(Damage, mean = Omega$Lambda4$loc, sd = Omega$Lambda4$sd) ^ DestDam_modifiers[2]
      D_Dam <- D_Dest_and_Dam - D_Dest
      D_Dam <- ifelse(D_Dam<0, 0, D_Dam)
    }
  }
  return(rbind(D_Dest, D_Dam))
}

# Building damage baseline hazard function hBD_0
hBD<-function(Ab,Population,rho,center){
  exp(-rho$A*(log(Ab)-center$A) - rho$H*(log(Population)-center$H))
}
# For building damage assessment data
fDamageBuilding<-function(BD,I,Params,Omega,linp,Ik){
  fDamUnscaled(I,Params,Omega)*linp#*hBD(BD$Ab,BD$Population,Omega$rho,Params$center[c("A","H")]),
}

qualifierDisp<-function(Disp,qualifier,mu) {
  if(qualifier%in%c("total","approximately")) return(Disp)
  else if(qualifier=="more than") {
    return(vapply(Disp,function(Disp) rgammaM(n=1,mu = mu$muplus*Disp,
                                              sig_percent = mu$sigplus),numeric(1)))
  } else if(qualifier=="less than") {
    return(vapply(Disp,function(Disp) rgammaM(n=1,mu = mu$muminus*Disp,
                                              sig_percent = mu$sigminus),numeric(1)))
  } else stop(paste0("qualifier is not recognised - ",qualifier))
}

# Binomial displacement calculator function
rbiny<-function(size,p) rbinom(n = 1,size,p); 
Fbdisp<-function(lPopS,Dprime) mapply(rbiny,lPopS,Dprime);

#when working with buildings, D_Disp is equivalent to D_BuildDam and D_Mort is equivalent to D_BuildDest
rmultinomy<-function(size, D_Disp, D_Mort, D_Rem) rmultinom(n=1, size, c(D_Disp, D_Mort, D_Rem))
Fbdam<-function(PopRem, D_Disp, D_Mort, D_Rem) mapply(rmultinomy, PopRem, D_Disp, D_Mort, D_Rem)

lBD<-function(D_B, BD_params){
  #input: probability of Building Damage D^B
  #output: probability of Building Destruction D^BD
  relative_probs = BDprob(D_B, BD_params) %>% mutate('damaged' = rowSums(.[1:4]))
  D_BD = (relative_probs[['destroyed']] + relative_probs[['severe']])/rowSums(relative_probs)
  return(D_BD)
}

fBD<-function(nbuildings, D_BD) mapply(rbiny, nbuildings, D_BD)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Log likelihood, posterior and prior distribution calculations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# These high-level priors are to include expert opinion in the 
# model calculations, not in the individual parameters themselves (for this, see priors)
Model$HighLevelPriors <-function(Omega,Model,modifier=NULL){
  
  min_gdp <- 2.200571; max_gdp <- 12.40876 #minimum and maximum log gdp from the dataset
  min_pdens <- 0; max_pdens <- 109043.5 #minimum and maximum population density from the dataset
  
  tp<-list(shape=Omega$dollar$k,M=Omega$dollar$M,
           scale=WeibullScaleFromShape(shape = Omega$dollar$k,center = Model$center$dollar))
  
  linp_min <- Plinpred(min_pdens, Omega$Pdens, Model$center, 1) * locpred(max_gdp,tp)
  linp_max <- Plinpred(max_pdens, Omega$Pdens, Model$center, 1) * locpred(min_gdp,tp)
  
  if(!is.null(modifier)) lp<-exp(as.numeric(unlist(modifier))) else lp<-1.
  lp <- c(linp_min, 1, linp_max) # lp_range - 0.361022 corresponds to lp for minimum GDP and PDens scaling, 2.861055 corresponds to maximum
  if(Model$haz=="EQ"){
  
    # Lower and upper bounds on the impacts at I_ij = 4.6, 6, and 9
    # in the order (DispMort, Mort, DamDest, Dest),
    # where DispMort is the sum of the probabilities of displacement and mortality
    # and DamDest is the sum of the probabilities of building damage and destruction.
                     #(Mort, DispMort, Dest, DamDest)
    Upp_bounds_4.6 <- c(0.01, 0.01, 0.01, 0.1)
    Low_bounds_6 <- c(0, 0.001, 0, 0.001)
    Upp_bounds_6 <- c(0.3, 0.5, 0.5, 0.8)
    Low_bounds_9 <- c(0.1,0.4,0.2,0.4)
    
    HLP_impacts <- function(I_ij, lp, Omega){
      rbind(apply(D_MortDisp_calc(h_0(I_ij, I0=4.5, theta=Omega$theta) * lp, Omega),2,cumsum), 
            apply(D_DestDam_calc(h_0(I_ij, I0=4.5, theta=Omega$theta) * lp, Omega), 2,cumsum))
    }
    adder <- 0
    adder <- sum(apply(HLP_impacts(4.6, lp, Omega), 2, function (x) sum(x > Upp_bounds_4.6)))
    
    adder <- adder + sum(apply(HLP_impacts(6, lp, Omega), 2, 
                               function (x) sum(c(x > Upp_bounds_6,x<Low_bounds_6))))
    
    adder <- adder + sum(apply(HLP_impacts(9, lp, Omega), 2, function (x) sum(x < Low_bounds_9)))

    return(adder) #looseend: need to address when including modifiers
    #return(adder)
    
  } else if(Model$haz=="TC"){
    
    Dfun<-function(I_ij) h_0(I = I_ij,I0 = 3,theta = Omega$theta) 
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
LL_IDP<-function(Y,  kernel_sd, kernel, cap){
  
  if (is.nan(Y[,'observed']) || is.nan(Y[,'sampled'])) return(0)
  
  LL <- 0
  k <- 10
  cap <- -100
  #impacts_observed <- intersect(which(Model$impacts$labels %in% colnames(Y)), which(predictions %in% colnames(Y)))

  if (kernel == 'loglaplace'){ #use a laplace kernel 
    LL_impact = log(dloglap(Y[,'observed']+k, location.ald = log(Y[,'sampled']+k), scale.ald = kernel_sd[[Y[1,'impact']]], tau = 0.5, log = FALSE)/
         (1-ploglap(k, location.ald = log(Y[,'sampled']+k), scale.ald = kernel_sd[[Y[1,'impact']]], tau = 0.5, log = FALSE)))
  } else if (kernel == 'lognormal'){ #use a lognormal kernel 
    LL_impact = log(dlnormTrunc(Y[,'observed']+k, log(Y[,'sampled']+k), sdlog=kernel_sd[[Y[1,'impact']]], min=k))
  } else {
    print(paste0("Failed to recognise kernel", AlgoParams$kernel))
    return(-Inf)
  }
  if(is.finite(LL_impact)){
    LL = LL + LL_impact
  } else {
    LL = LL + cap
  }
  return(LL)
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

LL_beta_apply<-function(b,value,BD_params) do.call(BD_params$functions[[value]],as.list(c(x=b,unlist(BD_params$Params[[value]]))))

BDprob<-function(b,BD_params){
  
  lls<-array(dim=c(length(b),length(BD_params$Params)))
  for(i in 1:length(BD_params$Params)){
    value<-(names(BD_params$Params))[i]
    lls[,i]<-vapply(b,FUN = LL_beta_apply,FUN.VALUE = numeric(1),
                    value=value,BD_params=BD_params)
  }
  lls%<>%as.data.frame.array();colnames(lls)<-names(BD_params$Params)
  return(lls)
}

predBD<-function(b,BD_params){
  lls<-BDprob(b,BD_params)
  return(mean(apply(lls, 1, function(x) sample(1:5, 1, prob=x)))) #LOOSEEND: Mean?
}

LL_BD<-function(b,classified,BD_params){
  
  lls<-BDprob(b,BD_params)
  if(classified=="Damaged") {
    tmp<-rowSums(lls)
    # Sum of all rows for the classifications that predict at least some damage
    return(log((tmp-lls[["notaffected"]])/tmp))
  }
  
  out<-log(lls[[classified]]/rowSums(lls))
  # machine level precision ~ -300
  out[is.infinite(out)]<--300
  return(out)
  
}

sampleBDdamage<-function(grading,n=10){
  
  vapply(1:length(grading),
         FUN = function(i) {
           
           if(grading[i]=="Damaged") return(runif(n))
           shapes<-Model$BD_params$Params[[grading[i]]]
           return(rbeta(n,shape1 = shapes$shape1,shape2 = shapes$shape2))
           
         },
         FUN.VALUE = numeric(n)
  )
  
}

GetIsoWeights<-function(dir){
  Dispy<-readRDS(paste0(dir,"IIDIPUS_Input/DispData_EQ_V2.Rdata"))
  WWW<-Dispy%>%group_by(iso3)%>%summarise(weights=1/length(gmax),.groups="drop_last")
}

Model$IsoWeights<-GetIsoWeights(dir)
Model$IsoWeights %<>% add_row(iso3='ABC', weights=1)

# Log-likelihood for displacement (ODD) objects
LL_Displacement<-function(dir,Model,proposed,AlgoParams,expLL=T){
  
  # Load ODD files
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  # Parallelise appropriately
  if(AlgoParams$AllParallel){
    # Task parallelism: this parallelisation calculates each event side-by-side, which is ideal if we have many CPU threads available and many ODD objects
    cores<-AlgoParams$cores
    AlgoParams$cores<-AlgoParams$NestedCores
    # When using task parallelisation, put the heaviest files first for optimisation reasons
    x <- file.info(paste0(folderin,ufiles))
    ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))]) #looseend
    
    tmpFn<-function(filer){
      # Extract the ODD object
      ODDy<-readRDS(paste0(folderin,filer))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      ODDy@fIndies<-Model$fIndies
      ODDy@impact%<>%as.data.frame.list()
      # Apply DispX
      tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",filer));return(-Inf)}
      
      return(tLL) #SMC-CHANGE
      # Weight the likelihoods based on the number of events for that country
      cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==ODDy@gmax$iso3[1]]
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Return the average log-likelihood
      if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL))
      else return(cWeight*mean(tLL,na.rm=T))
    }
    return(do.call(rbind, mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))) # SMC-CHANGE
    #return(sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
    
  } else {
    # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
    LL <- NULL
    for(i in 1:length(ufiles)){
      # Extract the ODD object
      ODDy<-readRDS(paste0(folderin,ufiles[i]))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      ODDy@fIndies<-Model$fIndies
      ODDy@impact%<>%as.data.frame.list()
      # Apply DispX
      tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));return(-Inf)}
      
      LL <- rbind(LL, tLL) #if want to return LL's for all events
      next
      # Weight the likelihoods based on the number of events for that country
      cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==ODDy@gmax$iso3[1]]
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Add the likelihood to the list of all events.
      if(expLL) {LL<-LL+cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)
      } else LL<-LL+cWeight*mean(tLL,na.rm=T)
    }
    
    return(LL)
    
  }
  
}

# Log-likelihood for building damage (BD) objects
LL_Buildings<-function(dir,Model,proposed,AlgoParams,expLL=T){
  # Load BD files
  folderin<-paste0(dir,"IIDIPUS_Input/BDobjects/")
  ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  
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
      if(nrow(BDy@data)==0){return(0)}
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=T),
                    error=function(e) NA)
      
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",filer));return(-Inf)}
      return(tLL) #SMC-CHANGE
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Return the average log-likelihood
      cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==BDy$ISO3C[which(!is.na(BDy$ISO3C))[1]]] #looseend
      if(is.na(cWeight)){return(-Inf)}                               
      if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)) 
      else return(cWeight*mean(tLL,na.rm=T))
      
    }
    return(do.call(rbind, mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))) #SMC-CHANGE #d-change
    #return(LL + sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
    
  } else {
    LL <- NULL
    # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
    for(i in 1:length(ufiles)){
      # Extract the BD object
      BDy<-readRDS(paste0(folderin,ufiles[i]))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=T),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));return(-Inf)}
      LL <- rbind(LL, tLL) #if want to return LL's for all events
      next
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Add the likelihood to the list of all events.
      if(expLL) LL<-LL+log(mean(exp(tLL-maxLL),na.rm=T))+maxLL
      else LL<-LL+mean(tLL,na.rm=T)
    }
    
    return(LL)
  }
}

# Bayesian Posterior distribution: This is for a group of ODD objects with observed data
logTarget<-function(dir,Model,proposed,AlgoParams,expLL=T){
  
  # Apply higher order priors
  if(!is.null(Model$HighLevelPriors)){
    HP<-Model$HighLevelPriors(proposed,Model)
    print(paste0("Higher Level Priors = ",HP))
  }

  
  # Approximate Bayesian Computation rejection
  #if(HP>AlgoParams$ABC) return(-5000000) # Cannot be infinite, but large (& negative) to not crash the Adaptive MCMC algorithm
  if (HP> AlgoParams$ABC) return(-Inf)
  
  # Add the log-likelihood values from the ODD (displacement) objects
  LL<-LL_Displacement(dir,Model,proposed,AlgoParams,expLL=T)
  #print(paste0("LL Displacements = ",LL)) ; sLL<-LL
  #if (any(LL == 0)){ #SMC-CHANGE
  #  return(Inf)
  #}
  # Add the log-likelihood values from the BD (building damage) objects
  #LL%<>%LL_Buildings(dir,Model,proposed,AlgoParams,expLL=T)
  LL_Build<-LL_Buildings(dir, Model, proposed, AlgoParams, expLL=T)
  #print(paste0("LL Building Damages = ",LL-sLL))
  
  posterior<-rbind(-LL, LL_Build)#LL #+HP
  # Add Bayesian priors
  if(!is.null(Model$priors)){
    posterior<-posterior+sum(Priors(proposed,Model$priors),na.rm=T)
  }
  #print(paste0("Posterior = ",posterior))
  
  return(posterior)
  
}

# -------------------------------------------------------------------------------------------------------------------
# ------------------------------- Pulling out and breaking down the distances ---------------------------------------
# -------------------------------------------------------------------------------------------------------------------

sampleDisps <-function(dir,Model,proposed,AlgoParams){
  
  # Load ODD files
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
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
      ODDy@fIndies<-Model$fIndies
      ODDy@impact%<>%as.data.frame.list()
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

sampleBDDist <- function(dir,Model,proposed,AlgoParams,expLL=T){
  # Load BD files
  folderin<-paste0(dir,"IIDIPUS_Input/BDobjects/")
  ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  
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
      BDy@fIndies<-Model$fIndies
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

sampleDist <- function(dir,Model,proposed,AlgoParams,expLL=T){

  Disps<-sampleDisps(dir,Model,proposed,AlgoParams)
  BDDists<- sampleBDDist(dir, Model, proposed, AlgoParams, expLL=T)

  return(list(Disps=Disps, BDDists=BDDists))
}

logTarget2 <- function(dist_sample, AlgoParams){
  
  sumLLs <- function(Disps_p){
    LL = 0 
    nrows <- NROW(Disps_p)
    for (i in 1:nrows){
      LL = LL + LL_IDP(Disps_p[i,], kernel_sd = AlgoParams$kernel_sd, kernel=AlgoParams$kernel, cap=-300)
    }
    return(LL)
  }
  LL_disps <- unlist(mclapply(dist_sample$Disps,FUN = sumLLs,mc.cores = 1)) # sum log likelihoods
  
  sumBD_dists <- function(BDDists_p){
    IC <- which(names(BDDists_p) == 'ICN' | names(BDDists_p) == 'ICD')
    sum(BDDists_p[IC])
  }
  
  LL_BD <- apply(dist_sample$BDDists, 2, sumBD_dists)
  
  dist_tot <- -LL_disps + LL_BD 
  
  return(dist_tot)
}

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

logTargetSingle<-function(ODDy,Model,Omega,AlgoParams){
  
  # HP<-0
  # # Apply higher order priors
  # if(!is.null(Model$HighLevelPriors)){
  #   HP<-Model$HighLevelPriors(Omega,Model,modifier = tryCatch(ODDy@modifier,error=function(e) NULL))
  #   # print(paste0("Higher Level Priors = ",HP))
  # }
  # 
  # LL<-HP
  
  tLL<-tryCatch(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params,
                      LL = T,Method = AlgoParams),
                error=function(e) NA)
  
  if(any(is.infinite(tLL)) | all(is.na(tLL))) stop("Failed to calculate LL")
  
  # LL<-LL+tLL
  
  return(tLL)
   
}


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
