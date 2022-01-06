#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@ FOR EARTHQUAKES @@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
Om_save<-Omega

vals<-c(0.03786569,0.19768294,0.09954531,2.02265867,5.40390699,0.67431138,0.12709657,3.52269924)
  # c(0.1601501,0.14,0.2053772,2.14934966,5.64216622,0.52240346,0.12552103,3.79169989)
  # c(0.02380763,0.03593319,0.05994856,1.61601944,4.71291825,0.59252792,0.11980378,3.79502670) # Not OLS
  # c(0.01601501,0.03270044,0.06053772,2.14934966,5.64216622,0.72240346,0.12552103,3.79169989) # LL=979
  # c(0.01786569,0.03768294,0.05954531,2.02265867,5.40390699,0.67431138,0.12709657,3.52269924) # LL_total=1159
  # c(0.01868896,0.03807269,0.05636423,2.06827208,5.46219578,0.68218565,0.12685656,3.35422034) # LL=990
  # c(0.01954243,0.03703514,0.05392705,2.08194530,5.32660921,0.67175072,0.12607799,3.74499350) #LL=989 2
  # c(0.01960936,0.03515436,0.05091546,2.17434362,5.12226427,0.68721752,0.13061263,4.25715482) # LL=983
  # c(0.01816326,0.03764612,0.05512244,2.09844350,5.12767547,0.66497198,0.12924055,3.93382976) # LL=989 1

Om_save$Lambda[1:3]<-vals[1:3]
Om_save$zeta[1:2]<-vals[4:5]
Om_save$theta[1]<-vals[6]
Om_save$eps[1:2]<-vals[7:8]
# Om_save$theta[1]<-vals[14]
# Om_save$eps[1:2]<-vals[15:16]

Dfun<-function(I_ij) mean(fDamUnscaled(I_ij,list(I0=4.5,Np=20),Om_save,0.032))
Dispfun<-function(I_ij) BinR(Dfun(I_ij)*Dfun(I_ij)*Om_save$Lambda$kappa+Om_save$Lambda$nu*Dfun(I_ij) + Om_save$Lambda$omega,Om_save$zeta)

Dispfun(4.6)
Dispfun(6)
Dispfun(9.5)

BDfun<-function(I_ij) BinR(Dfun(I_ij),Om_save$zeta)
I<-seq(from=4.05,to=9.5,length.out = 200)
Intensity<-data.frame(I_ij=rep(I,2),value=c(vapply(I,Dispfun,numeric(1)),vapply(I,BDfun,numeric(1))),
                      term=c(rep("Displacement",200),rep("Building Damage",200)))
Intensity$term%<>%as.factor()
ggplot(Intensity,aes(I_ij,value*100,group=term))+geom_point(aes(colour=term)) + 
  xlab("Earthquake Intensity [MMI]") + ylab("Value [%]")



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@ FOR CYCLONES @@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

Omega<-list(
  Lambda=list(kappa=0.1,nu=1.05,omega=0.004),
  zeta=list(k=9.5,lambda=0.23), # zeta=list(k=2.5,lambda=1.6),
  beta=list(CC.INS.GOV.GE=0,VU.SEV.AD=0,CC.INS.DRR=0,VU.SEV.PD=0,CC.INF.PHY=0,HA.NAT.TC=0,dollar=0,Pdens=0),
  theta=list(e=0.253965851), #list(e=0.25),
  # rho=list(A=0,H=0),
  eps=list(eps=0.060901098,xi=1.915351208)
  # mu=list(muplus=1,muminus=1,sigplus=0.001,sigminus=0.001)
)

Om_save<-Omega

vals<-c(0.17513286,0.66067661,0.00465764,5.58516690,0.21769465,0.43151790,0.08547787,2.03847584) # second round on cyclones, including extra data
  # c(0.03786569,0.19768294,0.09954531,2.02265867,5.40390699,0.67431138,0.12709657,3.52269924) # UNKNOWN, possibly first attempt

Om_save$Lambda[1:3]<-vals[1:3]
Om_save$zeta[1:2]<-vals[4:5]
# Om_save$theta[1]<-vals[6]
# Om_save$eps[1:2]<-vals[7:8]
Om_save$theta[1]<-vals[14]
Om_save$eps[1:2]<-vals[15:16]

Dfun<-function(I_ij) mean(fDamUnscaled(I_ij,list(I0=3,Np=20),Om_save,0.032))
Dispfun<-function(I_ij) BinR(Dfun(I_ij)*Dfun(I_ij)*Om_save$Lambda$kappa+Om_save$Lambda$nu*Dfun(I_ij) + Om_save$Lambda$omega,Om_save$zeta)
BDfun<-function(I_ij) BinR(Dfun(I_ij),Om_save$zeta)

BDfun(log(21))
BDfun(log(34))
BDfun(log(45))

I<-seq(from=3.04,to=4,length.out = 200)
Intensity<-data.frame(I_ij=rep(I,2),value=c(vapply(I,Dispfun,numeric(1)),vapply(I,BDfun,numeric(1))),
                      term=c(rep("Displacement",200),rep("Building Damage",200)))
Intensity$term%<>%as.factor()
ggplot(Intensity,aes(exp(I_ij)*3.6,value*100,group=term))+geom_point(aes(colour=term)) + 
  xlab("Cyclone Max. Wind Speed [km/h]") + ylab("Value [%]")

