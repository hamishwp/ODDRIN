
sder<-sapply(1:6, function(i) sd(output[,names(Omega$beta)[i]]) )
p50<-sd(output[,"p50p100"])

