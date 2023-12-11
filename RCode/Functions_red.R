# library(countrycode)
library(raster)
# library(geosphere)

returnX<-function(x,a=NULL,b=NULL) x
negexp <-function(x) -exp(x)
logneg <-function(x) log(-x)

# # 2D DENSITY PLOT WITH MODIFIED COLOUR AXIS
# ggplot(as.data.frame(buildings),aes(Longitude,Latitude))+
#   stat_density_2d(aes(fill = ..level..),breaks=c(0,0.5,1,3,5,10,30,50,100,300), 
#                   geom = "polygon")
# p1<-ggplot(as.data.frame(BDy),aes(Longitude,Latitude))+
#   stat_density_2d_filled(aes(fill=..level..),breaks=c(0,0.5,1,3,5,10,30,50,100,300)) +
#   ggtitle("Building Damage Data, EQ2015-04-25NPL")+theme(plot.title = element_text(hjust = 0.5))
# p2<-ggplot(as.data.frame(buildings),aes(Longitude,Latitude))+
#   stat_density_2d_filled(aes(fill=..level..),breaks=c(0,0.5,1,3,5,10,30,50,100,300)) +
#   ggtitle("OSM Building Data")+theme(plot.title = element_text(hjust = 0.5))
# gridExtra::grid.arrange(p1, p2, ncol=2)
# 
# tmp<-data.frame(PopDens=c(buildings$Pdensity,BDy$Population),ID=c(rep("OSM",nrow(buildings)),rep("Damaged",nrow(BDy))))
# ggplot(tmp,aes(x=PopDens,group=ID))+ stat_ecdf(aes(colour=ID),geom = "step",size=2) + scale_x_log10() + xlab("Population Density") +
#   ylab("Cumulative Distribution Function") + ggtitle("Nepal 2015 Earthquake") +theme(plot.title = element_text(hjust = 0.5))



# ggplot(tmp,aes(OSM,Damaged))+
#   stat_density_2d_filled(aes(fill=..level..),breaks=c(0,10,30,50,100,300,500,10000), geom="tile") +
#   ggtitle("Population Density of Damaged vs OSM Data")+theme(plot.title = element_text(hjust = 0.5))

AsYear<-function(date,red=F,limit=T){
  date%<>%as.Date
  if(!red) year<-as.numeric(format(date,"%Y"))
  else year<-as.numeric(format(date,"%y"))
  
  if(limit&any(year>as.numeric(format(Sys.Date(),"%Y")))) 
    year[year>as.numeric(format(Sys.Date(),"%Y"))]<-AsYear(Sys.Date())
  
  return(year)
}

AsMonth<-function(date){
  return(as.numeric(format(date,"%m")))
}
# Used to save output files by date and time for model validation comparison in time
DateTimeString<-function(){
  return(gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = ""))
}

# Assumes 8 income distribution percentiles
SplitSamplePop<-function(Pop,n=1){
  k<-length(Pop)
  return(array(vapply(Pop,function(tPop) rmultinom(n=n,
                                                   size=(tPop + rbinom(n=1,p=tPop%%1,size=1)), #LOOSEEND: same size for all Np
                                                   prob=rep(1/8,8)),FUN.VALUE = numeric(8L*n)),dim = c(8,k*n)))
}

rgammaM<-function(n,mu,sig_percent){
  # rgamma(n shape = alpha, scale = theta)
  # Note that the sig_percent is the constant coefficient of variation
  # Therefore, it is like a percentage of the mean
  ssq<-sig_percent*sig_percent
  rgamma(n,shape=1./ssq,scale=mu*ssq)
}

dgammaM<-function(x,mu,sig_percent,log=T){
  # rgamma(n shape = alpha, scale = theta)
  # Note that the sig_percent is the constant coefficient of variation
  # Therefore, it is like a percentage of the mean
  ssq<-sig_percent*sig_percent
  dgamma(x,shape=1./ssq,scale=mu*ssq,log=log)
}

extractnumbers<-function(str){
  return(as.numeric(unlist(regmatches(str,gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",str, perl=TRUE)))))
}