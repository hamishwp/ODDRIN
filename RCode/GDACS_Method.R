folder<-"./IIDIPUS_Input/ODDobjects/"
filez<-list.files(folder)

out<-data.frame()
for(fff in filez){
  
  ODDy<-readRDS(paste0(folder,fff))
  for (iso in ODDy@gmax$iso3){
    inISO<-ODDy@data$ISO3C==iso
    
    INFORM<-InterpINFORMdata("INFORM",min(max(ODDy@hazdates[1],as.Date("2017-10-22")),as.Date("2020-12-30")),iso=iso)
    if(nrow(INFORM)==0) INFORM<-data.frame(value=NA)
    
    out%<>%rbind(
      data.frame(
        gmax=ODDy@gmax$gmax[ODDy@gmax$iso3==iso],
        iso3=iso,
        date=ODDy@hazdates[1],
        Exp5=sum(ODDy@data$Population[ODDy@data$hazMean1>=5 & inISO],na.rm = T),
        Exp5.5=sum(ODDy@data$Population[ODDy@data$hazMean1>=5.5 & inISO],na.rm = T),
        Exp6=sum(ODDy@data$Population[ODDy@data$hazMean1>=6 & inISO],na.rm = T),
        Exp6.5=sum(ODDy@data$Population[ODDy@data$hazMean1>=6.5 & inISO],na.rm = T),
        Exp7=sum(ODDy@data$Population[ODDy@data$hazMean1>=7 & inISO],na.rm = T),
        Exp7.5=sum(ODDy@data$Population[ODDy@data$hazMean1>=7.5 & inISO],na.rm = T),
        Exp8=sum(ODDy@data$Population[ODDy@data$hazMean1>=8 & inISO],na.rm = T),
        Exp8.5=sum(ODDy@data$Population[ODDy@data$hazMean1>=8.5 & inISO],na.rm = T),
        Exp9=sum(ODDy@data$Population[ODDy@data$hazMean1>=9 & inISO],na.rm = T),
        INFORM=INFORM$value
      )
    )
  }
}

tmp<-dplyr::select(out,-c("iso3","date"))
names(tmp)[1]<-"Y"
tmp[,c("Y",grep(names(tmp),pattern = "Exp",value = T))]<-log(tmp[,c("Y",grep(names(tmp),pattern = "Exp",value = T))]+1)
View(tmp)

# Function from the file CorrelateModifier.R:
predictions<-LMFeatureSelection(tmp,Nb=12,intercept=F,fn="+",nlim=10)
View(predictions)

predvals<-exp(predict(lm(formula = "Y ~ Exp5+Exp6+Exp8+Exp9+INFORM+0",data = tmp),tmp))
visu<-data.frame(predicted=predvals,observed=out$gmax)

ggplot(visu,aes(predicted,observed))+
  geom_point()+geom_abline(slope=1,intercept=0)+scale_y_log10()+scale_x_log10()+
  ggtitle("GDACS Method with Displacement Data")+xlab("Predicted Values") + ylab("Observed Values") +
  theme(plot.title = element_text(hjust = 0.5)) 

