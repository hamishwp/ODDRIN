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


# Link functions (MUST BE SAME LENGTH AS OMEGA)
Model$links<-list(
  Lambda1=list(nu='exp',omega='exp'),
  Lambda2=list(nu='exp',omega='exp'),
  Lambda3=list(nu='exp',omega='exp'),
  zeta=list(k='exp',lambda='exp'), # zeta=list(k=2.5,lambda=1.6),
  # beta=list(xxx='exp',CC.INS.GOV.GE='exp',VU.SEV.AD='exp',CC.INS.DRR='exp',VU.SEV.PD='exp',CC.INF.PHY='exp'),
  Pdens=list(M='exp',k='exp'),
  dollar=list(M='negexp',k='exp'),
  theta=list(e='exp'), #list(e=0.25),
  # rho=list(A='exp',H='exp'),
  eps=list(eps='exp')#,xi='exp')
  # mu=list(muplus='exp',muminus='exp',sigplus='exp',sigminus='exp')
)
# names(Model$links$beta)[1]<-paste0("HA.NAT.",haz)

# And to go the other way....
Model$unlinks<-list(
  Lambda1=list(nu='log',omega='log'),
  Lambda2=list(nu='log',omega='log'),
  Lambda3=list(nu='log',omega='log'),
  zeta=list(k='log',lambda='log'), # zeta=list(k=2.5,lambda=1.6),
  # beta=list(xxx='log',CC.INS.GOV.GE='log',VU.SEV.AD='log',CC.INS.DRR='log',VU.SEV.PD='log',CC.INF.PHY='log'),
  Pdens=list(M='log',k='log'),
  dollar=list(M='logneg',k='log'),
  theta=list(e='log'), #list(e=0.25),
  # rho=list(A='log',H='log'),
  eps=list(eps='log')#,xi='log')
  # mu=list(muplus='exp',muminus='exp',sigplus='exp',sigminus='exp')
)
# names(Model$unlinks$beta)[1]<-paste0("HA.NAT.",haz)
# Skeleton
Model$skeleton <- list(
  Lambda1=list(nu=NA,omega=NA),
  Lambda2=list(nu=NA,omega=NA),
  Lambda3=list(nu=NA,omega=NA),
  zeta=list(k=NA,lambda=NA), # zeta=list(k=2.5,lambda=1.6),
  # beta=list(xxx=NA,CC.INS.GOV.GE=NA,VU.SEV.AD=NA,CC.INS.DRR=NA,VU.SEV.PD=NA,CC.INF.PHY=NA),
  Pdens=list(M=NA,k=NA),
  dollar=list(M=NA,k=NA),
  theta=list(e=NA), #list(e=0.25),
  # rho=list(A=NA,H=NA),
  eps=list(eps=NA)#,xi=NA)
  # mu=list(muplus=NA,muminus=NA,sigplus=NA,sigminus=NA)
)
# names(Model$skeleton$beta)[1]<-paste0("HA.NAT.",haz)

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
    names(linp)<-unique(ODD@cIndies$iso3)
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
# Damage function
fDamUnscaled<-function(I,Params,Omega){
    (h_0(I,Params$I0,Omega$theta)*
       stochastic(Params$Np,Omega$eps))%>%return()
}
# Baseline hazard function h_0
h_0<-function(I,I0,theta){
  ind<-I>I0
  h<-rep(0,length(I))
  h[ind]<-exp( theta$e*(I[ind]-I0) ) -1
  return(h)
}
# Building damage baseline hazard function hBD_0
hBD<-function(Ab,Population,rho,center){
  exp(-rho$A*(log(Ab)-center$A) - rho$H*(log(Population)-center$H))
}
# For building damage assessment data
fDamageBuilding<-function(BD,I,Params,Omega,linp,Ik){
  BinR(fDamUnscaled(I,Params,Omega)*linp,#*hBD(BD$Ab,BD$Population,Omega$rho,Params$center[c("A","H")]),
       Omega$zeta)
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
Fbdisp<-function(lPopS,Dprime) mapply(rbiny,lPopS,Dprime)

rmultinomy<-function(size, D_Disp, D_Mort, D_Rem) rmultinom(n=1, size, c(D_Disp,D_Mort,D_Rem))
Fbdam<-function(PopRem, D_Disp, D_Mort, D_Rem) mapply(rmultinomy, PopRem, D_Disp, D_Mort, D_Rem)

lBD<-function(D_B, BD_params){
  #input: probability of Building Damage D^B
  #output: probability of Building Destruction D^BD
  relative_probs = BDprob(D_B, BD_params) %>% mutate('damaged' = rowSums(.[1:4]))
  D_BD = (relative_probs[['destroyed']] + relative_probs[['severe']])/rowSums(relative_probs)
  return(D_BD)
}

fBD<-function(nbuildings, D_BD) sapply(D_BD, function(i) rbiny(size=nbuildings, p=D_BD))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Log likelihood, posterior and prior distribution calculations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# These high-level priors are to include expert opinion in the 
# model calculations, not in the individual parameters themselves (for this, see priors)
Model$HighLevelPriors<-function(Omega,Model,modifier=NULL){
  
  if(!is.null(modifier)) lp<-exp(as.numeric(unlist(modifier))) else lp<-1.
  
  if(Model$haz=="EQ"){
    
    Dfun<-function(I_ij) h_0(I = I_ij,I0 = 4.5,theta = Omega$theta) 
    Damfun<-function(I_ij, type) {
      D <- Dfun(I_ij)
      if (type=='Damage'){
        return(BinR(D,Omega$zeta))
      }
      if (type=='Buildings Destroyed'){
        return(BinR(Omega$Lambda3$nu * D + Omega$Lambda3$omega, Omega$zeta))
      }
      DispMort <- BinR(Omega$Lambda1$nu * D + Omega$Lambda1$omega,Omega$zeta)
      Mort <- BinR(Omega$Lambda2$nu * D + Omega$Lambda2$omega,Omega$zeta) * DispMort
      if (type=='Displacement'){
        return(DispMort - Mort)
      }
      else if (type=='Mortality'){
        return(Mort)
      }
    }

    adder<-rep(0,length(lp))
    for(i in 1:length(adder)){
      # Add lower bound priors:
      adder[i]<-sum(pweibull(Damfun(4.6, type='Displacement')*lp[i],3,0.01), pweibull(Damfun(4.6, type='Mortality')*lp[i],3,0.01),
                    pweibull(Damfun(4.6, type='Buildings Destroyed')*lp[i],3,0.01), pweibull(Damfun(4.6, type='Damage')*lp[i],3,0.01))
                    pweibull(Damfun(4.6, type='Buildings Destroyed')*lp[i],3,0.01)
      # Middle range priors:
      adder[i]<-adder[i]+sum(pweibull(Damfun(6, type='Displacement')*lp[i],15,0.8), 1- pweibull(Damfun(6, type='Displacement')*lp[i],3,0.005),
                             pweibull(Damfun(6, type='Mortality')*lp[i],15,0.8), 1- pweibull(Damfun(6, type='Mortality')*lp[i],3,0.005),
                             pweibull(Damfun(6, type='Buildings Destroyed')*lp[i],15,0.8), 1- pweibull(Damfun(6, type='Buildings Destroyed')*lp[i],3,0.005),
                             pweibull(Damfun(6, type='Damage')*lp[i],15,0.8), 1- pweibull(Damfun(6, type='Damage')*lp[i],3,0.005)
                             )
      # Upper bound priors:
      adder[i]<-adder[i]+sum(1-pweibull(Damfun(9, type='Displacement')*lp[i],5,0.3), 1-pweibull(Damfun(9, type='Mortality')*lp[i],5,0.2),
                             1-pweibull(Damfun(9, type='Buildings Destroyed')*lp[i],5,0.5), 1-pweibull(Damfun(9, type='Damage')*lp[i],30,0.85)) 
                             1-pweibull(Damfun(9, type='Buildings Destroyed')*lp[i],5,0.5)
    }
    return(adder)
    
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
LL_IDP<-function(Y){ #LOOSEEND is it fair to marginalize over missing measurements? 
  LL = 0
  if(is.numeric(Y$gmax)){
    #LL = LL + dnorm(log(Y$gmax+1), log(Y$disp_predictor+1), 0.1, log=TRUE)
    #LL = LL + log(dLaplace((log(Y$gmax+1)-log(Y$disp_predictor+1))/log(Y$disp_predictor+1), 0, 0.0333)) #LOOSEEND +1 to avoid log issues
    LL = LL + log(dLaplace(log(Y$gmax+1), log(Y$disp_predictor+1), 0.1))
  }
  if(is.numeric(Y$mortality)){
    #LL = LL + dnorm(log(Y$mortality+1), log(Y$mort_predictor+1), 0.1, log=TRUE)
    LL = LL + log(dLaplace(log(Y$mortality+1), log(Y$mort_predictor+1), 0.1))
    #LL = LL + log(dLaplace((log(Y$mortality+1)-log(Y$mort_predictor+1))/log(Y$mort_predictor+1), 0, 0.00333))
  }
  if(is.numeric(Y$buildDestroyed)){
    #LL = LL + dnorm(Y$buildDestroyed, Y$nBD_predictor, 100, log=TRUE)
    #LL = LL + dnorm(log(Y$buildDestroyed+1), log(Y$nBD_predictor+1), 0.1, log=TRUE)
    LL = LL + log(dLaplace(log(Y$buildDestroyed+1), log(Y$nBD_predictor+1), 0.1))
    #LL = LL + log(dLaplace(log(Y$buildDestroyed+1)-log(Y$nBD_predictor+1), 0, 0.2)) 
  }
  return(LL)
  
  # -log(1+(Y$gmax-Y$predictor)^2/Y$gmax)
  # -((Y$gmax-Y$predictor)^2/Y$gmax)
  # dnorm((log(max(Ystar,1))-log(Y)), mean = 0, sd = log(Y), log = T)
  # dgammaM(Ystar,Y,0.1,log = T)
}

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
  return(mean(apply(lls,1,which.max)))
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
LL_Displacement<-function(LL,dir,Model,proposed,AlgoParams,expLL=T){
  
  # Load ODD files
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  
  # Parallelise appropriately
  if(AlgoParams$AllParallel){
    # Task parallelism: this parallelisation calculates each event side-by-side, which is ideal if we have many CPU threads available and many ODD objects
    cores<-AlgoParams$cores
    AlgoParams$cores<-AlgoParams$NestedCores
    # When using task parallelisation, put the heaviest files first for optimisation reasons
    x <- file.info(paste0(folderin,ufiles))
    ufiles<-ufiles[match(length(ufiles):1,rank(x$size))]
    
    tmpFn<-function(filer){
      # Extract the ODD object
      ODDy<-readRDS(paste0(folderin,filer))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      ODDy@fIndies<-Model$fIndies
      ODDy@gmax%<>%as.data.frame.list()
      # Apply DispX
      tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",filer));return(-Inf)}
      # Weight the likelihoods based on the number of events for that country
      cWeight<-Model$IsoWeights$weights[Model$IsoWeights$iso3==ODDy@gmax$iso3[1]]
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Return the average log-likelihood
      if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL))
      else return(cWeight*mean(tLL,na.rm=T))
    }
  
    return(sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
    
  } else {
  
    # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
    for(i in 1:length(ufiles)){
      # Extract the ODD object
      ODDy<-readRDS(paste0(folderin,ufiles[i]))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      ODDy@fIndies<-Model$fIndies
      ODDy@gmax%<>%as.data.frame.list()
      # Apply DispX
      tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));return(-Inf)}
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
LL_Buildings<-function(LL,dir,Model,proposed,AlgoParams,expLL=T){
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
    ufiles<-ufiles[match(length(ufiles):1,rank(x$size))]
    
    tmpFn<-function(filer){
      # Extract the BD object
      BDy<-readRDS(paste0(folderin,filer))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,center = Model$center,Method=AlgoParams),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",filer));return(-Inf)}
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Return the average log-likelihood
      if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL))
      else return(cWeight*mean(tLL,na.rm=T))
      
    }
    
    return(sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
    
  } else {
    
    # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
    for(i in 1:length(ufiles)){
      # Extract the BD object
      BDy<-readRDS(paste0(folderin,ufiles[i]))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,center = Model$center,Method=AlgoParams),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));return(-Inf)}
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
  if(HP>AlgoParams$ABC) return(-5000000) # Cannot be infinite, but large (& negative) to not crash the Adaptive MCMC algorithm
  
  # Add the log-likelihood values from the ODD (displacement) objects
  LL<-LL_Displacement(0,dir,Model,proposed,AlgoParams,expLL=T)
  print(paste0("LL Displacements = ",LL)) ; sLL<-LL
  
  # Add the log-likelihood values from the BD (building damage) objects
  #LL%<>%LL_Buildings(dir,Model,proposed,AlgoParams,expLL=T)
  #print(paste0("LL Building Damages = ",LL-sLL))
  
  posterior<-LL #+HP
  # Add Bayesian priors
  if(!is.null(Model$priors)){
    posterior<-posterior+sum(Priors(proposed,Model$priors),na.rm=T)
  }
  print(paste0("Posterior = ",posterior))
  
  return(posterior)
  
}

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