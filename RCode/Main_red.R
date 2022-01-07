#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%% IIDIPUS - Integrated Internal DIsplaced PopUlation Sampler %%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%% Coded by Hamish Patten %%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%% Collaboration between the University of Oxford %%%%%%%%%%#
#%%%%%%%%%% and the Internal Displacement Monitoring Centre %%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%% Started January 2020 %%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Extract Environment Variables
source('RCode/GetEnv.R'); packred<-T
# Set the working directory from your environment variables
setwd(directory)
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')

library(dplyr)
library(magrittr)

# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
# GetODDPackages()
GetODDPackages_red()
# Basic functions:
source('RCode/Functions_red.R')
source('RCode/GetInitialValues.R')
# S4 object classes required:
source('RCode/ODDobj.R')
source('RCode/BDobj.R')
source('RCode/HAZARDobj.R')
# Disaster related:
# source('RCode/GetGDACS.R')
# source('RCode/GetUSGS.R')
# source('RCode/GetDisaster.R')
# IDP estimate related:
# source('RCode/GetDisplacements.R')
# source('RCode/GetHelix.R')
# Demography & population related:
# source('RCode/GetPopDemo.R')
# source('RCode/GetSocioEconomic.R')
# source('RCode/GetINFORM.R')
# Damage estimate related:
# source('RCode/GetOSM.R')
# source('RCode/GetSatDamage.R')
# Sourcing the data:
# source('RCode/GetData.R')

# IIDIPUSModelValidation<-function(dir="/home/patten/Documents/Coding/Oxford/IIDIPUS/",
#                                  haz="Earthquake", extractedData=T){
#  
#   # Extract model functions and priors
#   source('RCode/Model.R')
#   # Extract Displacement & Hazard data and match & reduce
#   ODDpaths<-ExtractData(haz,dir,extractedData)
#   # Extract Monte Carlo algorithm
#   source('RCode/Method.R')
#   # Estimate initial values
#   iVals<-GetInitVals(dir,haz,Model,AlgoParams)
#   
#   # HERE WE GOOOOOOOOOOOOOOOO
#   output <- Algorithm(InputDataDir=ODDpaths,
#                              Model=Model,
#                              iVals=iVals,
#                              AlgoParams=AlgoParams)
#   
#   saveRDS(object = output,
#           file = paste0(dir,"IIDIPUS_Results/ForgeEngine/output_",
#                         DateTimeString(),".Rdata"))
#   
#   return(output)
#    
# }

source('RCode/Model.R')
# source('RCode/Method.R')

# if(Model$BinR=="weibull"){
  Omega<-list(
    Lambda=list(kappa=0.01786569,nu=0.03768294,omega=0.05954531),
    zeta=list(k=2.02265867,lambda=5.40390699), # zeta=list(k=2.5,lambda=1.6),
    beta=list(CC.INS.GOV.GE=0,VU.SEV.AD=0,CC.INS.DRR=0,VU.SEV.PD=0,CC.INF.PHY=0,HA.NAT.EQ=0,dollar=0,Pdens=0),
    theta=list(e=0.67431138), #list(e=0.25),
    rho=list(A=0,H=0),
    eps=list(eps=0.12709657,xi=3.52269924),
    mu=list(muplus=1,muminus=1,sigplus=0.001,sigminus=0.001)
  )
# } 
# if(Model$BinR=="gompertz"){
#   Omega<-list(
#     Lambda=list(kappa=0.01,nu=1.1,omega=0.2),
#     zeta=list(varrho=5,eta=0.3),
#     beta=list(CC.INS.GOV.GE=0,VU.SEV.AD=0,CC.INS.DRR=0,VU.SEV.PD=0,CC.INF.PHY=0,HA.NAT.EQ=0,dollar=0),
#     theta=list(e=0.12),
#     rho=list(A=0,H=0),
#     eps=list(eps=0.1,xi=3),
#     mu=list(muplus=1.1,muminus=0.8,sigplus=0.1,sigminus=0.1)
#   )
# }
center<-ExtractCentering(dir = ODDpath,haz = haz,Model = Model)
Model$center<-center

# LL<-logTarget(dir,Model,Omega,AlgoParams[c("Np","cores")],expLL = T)
# LL<-logTarget(dir,Model,Omega,AlgoParams[c("Np","cores")],expLL = F)

# First with Omega$beta <- Omega$eps$xi <- Omega$rho0 <-0
Fopty<-function(vals){
  
  # vals[1]%<>%match.fun(unname(unlist(Model$links[["Lambda"]]))[1])()
  # vals[2]%<>%match.fun(unname(unlist(Model$links[["Lambda"]]))[2])()
  # vals[3]%<>%match.fun(unname(unlist(Model$links[["zeta"]]))[1])()
  # vals[4]%<>%match.fun(unname(unlist(Model$links[["zeta"]]))[2])()
  # vals[5]%<>%match.fun(unname(unlist(Model$links[["theta"]]))[1])()
  # vals[6]%<>%match.fun(unname(unlist(Model$links[["eps"]]))[1])()
  # vals[7]%<>%match.fun(unname(unlist(Model$links[["eps"]]))[2])()
  
  vals%<>%exp()
  
  print(vals)
  
  Omega$Lambda[1:3]<-vals[1:3]
  Omega$zeta[1:2]<-vals[4:5]
  Omega$theta[1]<-vals[6]
  Omega$eps[1:2]<-vals[7:8]
  
  LL<-logTarget(dir,Model,Omega,AlgoParams[c("Np","cores")],expLL = F)
  print(LL)
  print("...")
  return(LL)
  
} 

ivals<-log(unname(unlist(Omega[c("Lambda","zeta","theta","eps")])))
# ivals<-c(ivals[1:2],log(ivals[3:5]),ivals[6],log(ivals[7]))

#####################################################################
output<-optim(par=ivals,
              fn = Fopty,control = list(maxit = 150,fnscale=-1,reltol=1.5e-4),hessian = T)
#####################################################################

vals<-output$par
vals[1]%<>%match.fun(unname(unlist(Model$links[["Lambda"]]))[1])()
vals[2]%<>%match.fun(unname(unlist(Model$links[["Lambda"]]))[2])()
vals[3]%<>%match.fun(unname(unlist(Model$links[["zeta"]]))[1])()
vals[4]%<>%match.fun(unname(unlist(Model$links[["zeta"]]))[2])()
vals[5]%<>%match.fun(unname(unlist(Model$links[["theta"]]))[1])()
vals[6]%<>%match.fun(unname(unlist(Model$links[["eps"]]))[1])()
vals[7]%<>%match.fun(unname(unlist(Model$links[["eps"]]))[2])()

Omega$Lambda[1:3]<-vals[1:3]
Omega$zeta[1:2]<-vals[4:5]
Omega$theta[1]<-vals[6]
Omega$eps[1:2]<-vals[7:8]

Fopty<-function(vals){
  
  for (i in 1:length(Omega$beta)) Omega$beta[i]<-match.fun(unname(unlist(Model$links$beta[names(Omega$beta[i])])))(vals[i])
  for (i in 1:2) Omega$eps[i]<-match.fun(unname(unlist(Model$links$eps[names(Omega$eps[i])])))(vals[i+length(Omega$beta)])
  
  print(unname(unlist(Omega[c("beta","eps")])))
  
  LL<-logTarget(dir,Model,Omega,AlgoParams[c("Np","cores")],expLL = F)
  print(LL)
  print("...")
  return(LL)
} 

# ivals<-unname(unlist(list(CC.INS.GOV.GE=5,VU.SEV.AD=5,CC.INS.DRR=5,VU.SEV.PD=5,CC.INF.PHY=5,HA.NAT.EQ=5,dollar=1e-1)))
ivals<-c(rep(log(2),6),log(0.1),0.1,log(unname(unlist(Omega$eps[1:2]))))

output<-optim(par=ivals,#method = "SANN",
              fn = Fopty,control = list(maxit = 200,fnscale=-1,reltol=1.5e-3),
              hessian = T)

vals<-output$par
for (i in 1:length(vals)) Omega$beta[i]<-match.fun(unname(unlist(Model$links$beta[names(Omega$beta[i])])))(vals[i])

Fopty<-function(vals){
  
  vals[1]%<>%match.fun(unname(unlist(Model$links[["Lambda"]]))[1])()
  vals[2]%<>%match.fun(unname(unlist(Model$links[["Lambda"]]))[2])()
  vals[3]%<>%match.fun(unname(unlist(Model$links[["zeta"]]))[1])()
  vals[4]%<>%match.fun(unname(unlist(Model$links[["zeta"]]))[2])()
  vals[5]%<>%match.fun(unname(unlist(Model$links[["theta"]]))[1])()
  vals[6]%<>%match.fun(unname(unlist(Model$links[["eps"]]))[1])()
  vals[7]%<>%match.fun(unname(unlist(Model$links[["eps"]]))[2])()
  
  for (i in 1:7) Omega$beta[i]<-match.fun(unname(unlist(Model$links$beta[names(Omega$beta[i])])))(vals[(i+7)])
  
  Omega$Lambda[1:2]<-vals[1:2]
  Omega$zeta[1:2]<-vals[3:4]
  Omega$theta[1]<-vals[5]
  Omega$eps[1:2]<-vals[6:7]
  
  print(unname(unlist(Omega[c("Lambda","zeta","theta","eps","beta")])))
  
  LL<-logTarget(dir,Model,Omega,AlgoParams[c("Np","cores")],expLL = F)
  print(LL)
  print("...")
  return(LL)
  
} 

ivals<-unname(unlist(Omega[c("Lambda","zeta","theta","eps","beta")]))
ivals<-c(ivals[1:2],log(ivals[3:5]),ivals[6],log(ivals[7]),log(ivals[8:13]),log(-ivals[14]))  

output<-optim(par=ivals,
              fn = Fopty,control = list(maxit = 200,fnscale=-1,reltol=1.5e-3),
              hessian = T)