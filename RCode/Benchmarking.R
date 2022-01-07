
DispData<-as.data.frame(read_csv(paste0(directory,"IIDIPUS_Input/DispData_EQ_V2.csv")))
DispData$PAGER[DispData$PAGER=="green"]<-"Green";DispData$PAGER[DispData$PAGER=="yellow"]<-"Yellow";DispData$PAGER[DispData$PAGER=="red"]<-"Red";DispData$PAGER[DispData$PAGER=="orange"]<-"Orange"
DispData$PAGER<-factor(as.factor(DispData$PAGER), levels = c("Green","Yellow","Orange","Red",NA))
names(DispData)[1]<-"eventid"

tmp<-DispData%>%group_by(eventid)%>%summarise(IDP=sum(gmax),alert=unique(PAGER))
p<-ggplot(tmp)+geom_boxplot(aes(alert,IDP,fill=alert))+scale_y_log10() + geom_jitter(aes(alert,IDP),height = 0,width = 0.1) + 
  scale_fill_discrete(type = c("green","yellow2","orange","red","blue")) + ylab("Maximum Number of Displaced Persons") +
  xlab("USGS-PAGER Alertscore")
ggsave("PAGER_Benchmarking.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results'),width = 8,height = 5)

tmp<-DispData%>%group_by(eventid)%>%summarise(IDP=sum(gmax),alert=max(GDACS,na.rm=T),USGS=unique(PAGER))
model <- lm(IDP ~ exp(alert)+0, data = tmp)
tmp$pred<-model$fitted.values
tmp$gLL<-abs(tmp$pred-tmp$IDP)/tmp$IDP

p<-ggplot(tmp)+geom_point(aes(alert,IDP))+scale_y_log10() + 
  ylab("Maximum Number of Displaced Persons") + xlab("GDACS Alertscore")
ggsave("GDACS_Benchmarking.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results'),width = 8,height = 5)

p<-ggplot(tmp)+geom_point(aes(alert,IDP,colour=USGS),size=4)+geom_point(aes(alert,IDP),shape=4)+scale_y_log10() + 
  scale_color_discrete(type = c("green","yellow2","orange","red","blue")) +
  ylab("Maximum Number of Displaced Persons") + xlab("GDACS Alertscore") + ggtitle("Benchmark of GDACS vs USGS-PAGER")
ggsave("GDACS-PAGER_Benchmarking.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results'),width = 8,height = 5)

output<-ODDypreds(directory,"EQ")
output$iLL<-abs(output$gmax-output$predictor)/output$gmax

outy<-output%>%group_by(eventid)%>%summarise(gmax=sum(gmax,na.rm=T),predictor=sum(predictor,na.rm=T))
ttt<-merge(outy,tmp,by="eventid")
ttt%<>%filter(!is.infinite(alert))

USGS<-DispData%>%group_by(eventid)%>%summarise(IDP=sum(gmax),alert=unique(PAGER))
USGS%<>%filter(eventid%in%ttt$eventid & !is.na(alert))

USGS$factor<-as.numeric(USGS$alert)
umodel <- lm(IDP ~ factor, data = USGS)

tmp2<-ttt%>%filter(!is.na(alert))%>%group_by(eventid)%>%
  summarise(IDP=sum(gmax,na.rm = T),predictor=sum(predictor,na.rm = T),
            namer=unique(namer),alert=max(alert))
tmp2%<>%filter(eventid%in%USGS$eventid & !is.na(alert))

gmodel <- lm(IDP ~ alert+0, data = tmp2)
mmodel <- lm(IDP ~ predictor+0, data = tmp2)

sum(abs(gmodel$residuals))/sum(abs(mmodel$residuals))
sum(abs(umodel$residuals))/sum(abs(mmodel$residuals))

mean(100*abs(mmodel$residuals)/tmp2$IDP)
mean(100*abs(gmodel$residuals)/tmp2$IDP)
mean(100*abs(umodel$residuals)/USGS$IDP)

median(100*abs(mmodel$residuals)/tmp2$IDP)
median(100*abs(gmodel$residuals)/tmp2$IDP)
median(100*abs(umodel$residuals)/USGS$IDP)

benchmarker<-data.frame(Observed=rep(tmp2$IDP,3),
                        Prediction=c(mmodel$fitted.values,
                                     gmodel$fitted.values,
                                     umodel$fitted.values),
                        Model=c(rep("IIDIPUS",nrow(tmp2)),
                                rep("GDACS-Alertscore",nrow(tmp2)),
                                    rep("USGS-PAGER",nrow(tmp2))))

p<-ggplot(benchmarker,aes(Prediction,Observed,group=Model))+geom_point(aes(colour=Model)) + 
  scale_x_log10()+scale_y_log10()+geom_abline(slope=1);p
ggsave("ALLBenchmarking.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results'),width = 8,height = 5)

DispData$continent<-convIso3Continent(DispData$iso3)





Dfun<-function(I_ij) h_0(I = I_ij,I0 = 4.5,theta = Omega$theta) 
Dispfun<-function(I_ij) BinR(Dfun(I_ij)*Dfun(I_ij)*Omega$Lambda$kappa+Omega$Lambda$nu*Dfun(I_ij) + Omega$Lambda$omega,Omega$zeta)
Damfun<-function(I_ij) BinR(Dfun(I_ij),Omega$zeta)

noise<-as.numeric(Omega$eps$eps)

ODDy<-readRDS("IIDIPUS_Input/ODDobjects/EQ20150425NPL_-27")
# iso3<-"NPL"
notnans<-which(!(is.na(ODDy$Population) | is.na(ODDy$ISO3C) | is.na(ODDy$GDP)))
SincN<-seq(10,90,by = 10); Sinc<-ExtractCIndy(ODDy,var = paste0("p",SincN,"p100"))
GDP<-GDPlinp(ODDy,Sinc,Omega$dollar,Model$center,notnans)
NPLDollar<-max(GDP$dGDP$linp[GDP$dGDP$ind==GDP$iGDP[notnans]],na.rm=T)/min(GDP$dGDP$linp[GDP$dGDP$ind==GDP$iGDP[notnans]],na.rm=T)

max(1/(BinR(vapply(50:95/10,Dfun,numeric(1))*min(NPLDollar)*(1-noise),Omega$zeta)/
         BinR(vapply(50:95/10,Dfun,numeric(1))*max(NPLDollar)*(1+noise),Omega$zeta)))

NPLPlinp<-Plinpred(ODDy@data$Population,
                Omega$Pdens,
                Model$center,
                notnans)
# cIndies<-ODDy@cIndies[endsWith(ODDy@cIndies$variable,"p100") & ODDy@cIndies$iso3==iso3,]
# cIndies$value<-cIndies$value/cIndies$value[cIndies$variable=="p50p100"]
# GDP<-ODDy@data%>%group_by(ISO3C)%>%summarise(GDP=mean(GDP))
# 
# NPLdollarzzz<-exp((cIndies$value*log(GDP$GDP[GDP$ISO3C==iso3 & !is.na(GDP$ISO3C)])-Model$center$dollar)*Omega$dollar)

# noise<-0.08822774


# NPLGDP<-BinR(vapply(50:95/10,Dfun,numeric(1))*min(NPLdollarzzz),Omega$zeta)

ODDy<-readRDS("IIDIPUS_Input/ODDobjects/EQ20160824ITA_724")
notnans<-which(!(is.na(ODDy$Population) | is.na(ODDy$ISO3C) | is.na(ODDy$GDP)))
SincN<-seq(10,90,by = 10); Sinc<-ExtractCIndy(ODDy,var = paste0("p",SincN,"p100"))
GDP<-GDPlinp(ODDy,Sinc,Omega$dollar,Model$center,notnans)
ITADollar<-max(GDP$dGDP$linp[GDP$dGDP$ind==GDP$iGDP[notnans]],na.rm=T)/min(GDP$dGDP$linp[GDP$dGDP$ind==GDP$iGDP[notnans]],na.rm=T)

max(1/(BinR(vapply(50:95/10,Dfun,numeric(1))*min(ITADollar)*(1-noise),Omega$zeta)/
         BinR(vapply(50:95/10,Dfun,numeric(1))*max(ITADollar)*(1+noise),Omega$zeta)))


max(1/(BinR(vapply(50:95/10,Dfun,numeric(1))*min(NPLDollar)*(1-noise),Omega$zeta)/
         BinR(vapply(50:95/10,Dfun,numeric(1))*max(ITADollar)*(1+noise),Omega$zeta)))


ITAPlinp<-Plinpred(ODDy@data$Population,
                   Omega$Pdens,
                   Model$center,
                   notnans)

max(1/(BinR(vapply(50:95/10,Dfun,numeric(1))*(exp(Omega$Pdens$M))*(1-noise),Omega$zeta)/
         BinR(vapply(50:95/10,Dfun,numeric(1))*exp(-Omega$Pdens$M)*(1+noise),Omega$zeta)))

max(ITAPlinp,na.rm = T)/min(ITAPlinp,na.rm = T)
max(NPLPlinp,na.rm = T)/min(NPLPlinp,na.rm = T)






iso3<-"ITA"

cIndies<-ODDy@cIndies[endsWith(ODDy@cIndies$variable,"p100") & ODDy@cIndies$iso3==iso3,]
cIndies$value<-cIndies$value/cIndies$value[cIndies$variable=="p50p100"]
GDP<-ODDy@data%>%group_by(ISO3C)%>%summarise(GDP=mean(GDP))

ITAdollarzzz<-exp((cIndies$value*log(GDP$GDP[GDP$ISO3C==iso3 & !is.na(GDP$ISO3C)])-Model$center$dollar)*dollar)

max(1/(BinR(vapply(50:95/10,Dfun,numeric(1))*min(ITAdollarzzz)*(1-0.08822774),Omega$zeta)/
         BinR(vapply(50:95/10,Dfun,numeric(1))*max(ITAdollarzzz)*(1+0.08822774),Omega$zeta)))
# ITAGDP<-BinR(vapply(50:95/10,Dfun,numeric(1))*max(ITAdollarzzz)*(1+0.08822774),Omega$zeta)











tmp2<-outy%>%group_by(eventid)%>%
  summarise(IDP=sum(gmax,na.rm = T),predictor=sum(predmod,na.rm = T),
            namer=unique(namer))
tmp2$iLL<-abs(tmp2$IDP-tmp2$predictor)/tmp2$IDP

benchy<-merge(tmp2,tmp,by="eventid")
bdens<-rbind(cbind(unname(dplyr::select(benchy,IDP.x,iLL,predictor)),data.frame(Model=rep("IIDIPUS",nrow(benchy)))),
             cbind(unname(dplyr::select(benchy,IDP.x,gLL,pred)),data.frame(Model=rep("GDACS",nrow(benchy)))))
names(bdens)<-c("IDP","LL","Prediction","Model")

ggplot(bdens,aes((LL),group=Model))+geom_density(aes(fill=Model),alpha=0.5)+xlab("Log-Likelihood")+scale_x_log10()

p<-ggplot(tmp2)+geom_point(aes(predictor,IDP))+scale_y_log10() + scale_x_log10() +
  geom_abline(slope = 1,intercept = 0) +
  ylab("Maximum Number of Displaced Persons") + xlab("Predicted Max. Displaced")
ggsave("IIDIPUS_Benchmarking.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results'),width = 8,height = 5)

ttt<-cor.test(~IDP + exp(alert)+0,
         data=tmp,
         method = "pearson",
         conf.level = 0.95)
print(paste0("GDACS Correlation = ",unname(ttt$estimate)))
print(paste0("GDACS Correlation Confidence Interval = ",unname(ttt$conf.int)))

ttt<-cor.test(~ IDP + predictor+0,
         data=tmp2,
         method = "pearson",
         conf.level = 0.95)
print(paste0("IIDIPUS Correlation = ",unname(ttt$estimate)))
print(paste0("IIDIPUS Correlation Confidence Interval = ",unname(ttt$conf.int)))

tmp2$iALERT<-rep(NA,nrow(tmp2))
tmp2$iALERT[tmp2$predictor<2000]<-"Green"
tmp2$iALERT[tmp2$predictor>=2000 & tmp2$predictor<10000]<-"Yellow"
tmp2$iALERT[tmp2$predictor>=10000 & tmp2$predictor<100000]<-"Orange"
tmp2$iALERT[tmp2$predictor>=100000]<-"Red"
tmp2$iALERT<-factor(as.factor(tmp2$iALERT), c("Green","Yellow","Orange","Red",NA))

p1<-ggplot(tmp2)+geom_boxplot(aes(iALERT,IDP,fill=iALERT))+scale_y_log10() + geom_jitter(aes(iALERT,IDP),height = 0,width = 0.1) + 
  ylab("Maximum Number of Displaced Persons") +
  xlab("IIDIPUS Alertscore") + 
  scale_fill_manual(labels=c("Green (IDP<2,000)","Yellow (2,000<=IDP<10,000)",
                             "Orange (10,000<=IDP<100,000)","Red (IDP=>100,000)"),
                    values=c("green","yellow2","orange","red","blue"));p1

tPAGER<-DispData%>%group_by(eventid)%>%summarise(IDP=sum(gmax),alert=unique(PAGER))
p2<-ggplot(filter(tPAGER,!is.na(alert)))+geom_boxplot(aes(alert,IDP,fill=alert))+scale_y_log10() + geom_jitter(aes(alert,IDP),height = 0,width = 0.1) + 
  scale_fill_discrete(type = c("green","yellow2","orange","red","blue")) + ylab("Maximum Number of Displaced Persons") +
  xlab("USGS-PAGER Alertscore");p2

p<-cowplot::plot_grid(p1,p2,nrow=1,rel_widths = c(0.55,0.45))

ggsave("IIDIPUS-PAGER_Alertscore_Benchmarking.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results'),width = 12,height = 4)

summary(lm(IDP~exp(as.numeric(iALERT)),data = tmp2))
summary(lm(IDP~exp(as.numeric(alert)),data = tmp))

values<-1:1000/10
conf<-c()
for (i in values){
  
  ttt$low<-ttt$predmod/(2*i)
  ttt$up<-ttt$predmod*(i*2)
  tmp3<-ttt%>%group_by(eventid)%>%
    summarise(IDP=sum(gmax,na.rm = T),predictor=sum(predmod,na.rm = T),
              up=sum(up,na.rm = T),low=sum(low,na.rm = T),.groups = 'drop_last')
  conf%<>%c(sum(tmp3$up>tmp3$IDP & tmp3$low<tmp3$IDP)/length(tmp3))
  
}
View(sort(desc(conf)))
