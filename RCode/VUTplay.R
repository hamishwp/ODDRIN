country<-"Vanuatu"
namez<-c("Espiritu Santo","Aoba","Maewo","Pentecost Island","Malakula","Ambrym","Efate","Epi")

bbox<-c(166.29583,-16.99583,168.99583,-14.49583)
  
polyz<-list()
for(nn in namez){
  
  polyz%<>%c(getbb(paste0(nn,", ",country),
              viewbox=strcat(as.character(bbox),collapse = ","),
              format_out = 'sf_polygon'))
  
}

mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "toner-lite",zoom=9)

pp<-ggmap(mad_map) + xlim(c(166.3,169)) + ylim(c(-17,-14.5))
for(i in 1:length(namez)){
  pp<-pp+geom_sf(data=polyz[[i]],fill="red",alpha=0.2, inherit.aes = FALSE)
}
pp
ggsave("VUT_outline.png", plot=pp,path = paste0(directory,'Plots/IIDIPUS_Results/VUT'),width = 6,height = 5)
  
demonotone<-function(ni,p,q,nt){
  
  val<-rep(ni,nt)
  for(t in 1:nt){
    if(q>runif(1)) {val[t]<-ni; next}
    ni<-ni-rbinom(1,ni,p*(0.5+rbeta(1,20,20)))
    val[t]<-ni
  }
  val[val<0]<-0
  return(val)
  
}

holdback<-0.6
men<-0.1
women<-0.05
days<-40
total<-300
swomen<-rbinom(1,total,0.6)
smen<-total-swomen

plot(demonotone(swomen,women,holdback,days),col="blue",type="l",ylim=c(0,max(smen,swomen)))
lines(demonotone(smen,men,holdback,days),col="purple")

Insiders<-function(tLL,subpoly){
  
  if(length(class(subpoly[[1]][[1]]))==1) {
    subpoly<-st_coordinates(subpoly[[1]])
    if(length(unique(c(subpoly[,3:ncol(subpoly)])))!=1) print(unique(subpoly[,3:ncol(subpoly)]))
    return(sp::point.in.polygon(tLL@coords[,1],
                                tLL@coords[,2],
                                subpoly[,1],
                                subpoly[,2])>0)
  }
  
  subpoly<-subpoly[[1]][[1]][[1]]
  
  logicz<-rep(0,nrow(tLL@data))
  for(i in 1:length(subpoly)){
    
    subz<-subpoly[[i]]
    # subz<-st_coordinates(subpoly[[i]])
    # if(length(unique(c(subz[,3:ncol(subz)])))!=1) print(unique(subz[,3:ncol(subz)]))
    
    logicz<-logicz+as.numeric(sp::point.in.polygon(tLL@coords[,1],
                                tLL@coords[,2],
                                subz[,1],
                                subz[,2])>0)
  }
  
  return(logicz>0)

}

geoinsights<-data.frame()
indies<-data.frame()
for(i in 1:length(namez)){
  nn<-namez[i]
  
  inny<-Insiders(tLL,polyz[[i]])
  
  total<-as.integer(sum(tLL@data$Disp[inny],na.rm = T))
  swomen<-rbinom(1,total,0.55*(0.5+rbeta(1,20,20)))
  smen<-total-swomen
  
  indies%<>%rbind(data.frame(
    inny=which(inny),
    rep(nn,sum(inny))
  ))
  
  geoinsights%<>%rbind(
    data.frame(
      Displaced=c(demonotone(swomen,women,holdback,days),
                  demonotone(smen,men,holdback,days)),
      day=rep(1:days,2),
      Sex=c(rep("Female",days),rep("Male",days)),
      namez=rep(nn,2*days)
    )
  )
  
}

p<-ggplot(geoinsights,aes(day,Displaced,colour=Sex))+geom_line(aes(linetype=Sex),size=1.5)+
  scale_color_manual(values=c("Female"="darkolivegreen","Male"="darkorchid"))
p<-p+facet_wrap( ~ namez, scales = "free")+ theme(plot.title = element_text(hjust = 0.5));p
ggsave("FBgeoinsights_timeline_gender.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results/VUT'),width = 9,height = 5)

normy<-geoinsights%>%group_by(namez,Sex)%>%summarise(normy=Displaced/max(Displaced),day=day)

qq<-ggmap(mad_map) + xlim(c(166.3,169)) + ylim(c(-17,-14.5)) + xlab("Longitude") + ylab("Latitude")
maxy<-geoinsights%>%pull(Displaced)%>%max()
maxODD<-max(ODDy@data$Disp,na.rm = T)

for(j in 1:days){
  
  geosub<-filter(geoinsights,day==j)
  printers<-geosub%>%group_by(Sex)%>%summarise(sumy=sum(Displaced),.groups='drop_last')
  normsub<-filter(normy,day==j)
  
  pp1<-qq + ggtitle(paste0("Female Disp = ",printers$sumy[printers=="Female"])) + theme(plot.title = element_text(hjust = 0.5))
  pp2<-qq + ggtitle(paste0("Male Disp = ",printers$sumy[printers=="Male"])) + theme(plot.title = element_text(hjust = 0.5))
  rr1<-pp1+labs(fill = "Displaced")
  rr2<-pp2+labs(fill = "Displaced")
  
  for(i in 1:length(namez)){
    nn<-namez[i]
    colz<-geosub%>%filter(namez==nn)
    fdec<-pull(filter(colz,Sex=="Female"),Displaced)
    mdec<-pull(filter(colz,Sex=="Male"),Displaced)
    
    pp1<-pp1+geom_sf(data=polyz[[i]],fill="red",
                     alpha=fdec/maxy,
                                inherit.aes = FALSE)

    pp2<-pp2+geom_sf(data=polyz[[i]],fill="red",
                     alpha=mdec/maxy,
                                inherit.aes = FALSE)
    
    multy<-normsub%>%filter(namez==nn)
    
    tempy<-tLL
    inny<-Insiders(tLL,polyz[[i]])
    tempy@data$Disp[!inny]<-NA
    ftempy<-mtempy<-tempy
    ftempy@data$Disp[inny]<-tempy@data$Disp[inny]*as.numeric(multy[multy$Sex=="Female","normy"])
    mtempy@data$Disp[inny]<-tempy@data$Disp[inny]*as.numeric(multy[multy$Sex=="Male","normy"])
    
    rr1<-rr1 + geom_raster(data=as.data.frame(ftempy),aes(Longitude,Latitude,fill=Disp),
                           interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
      scale_fill_gradient2(low = "blue",mid="blue",high = "red",trans = "log",
                           breaks=c(0,1,10,100),
                           na.value = "transparent", limits=c(1,maxODD))

    rr2<-rr2 + geom_raster(data=as.data.frame(mtempy),aes(Longitude,Latitude,fill=Disp),
                           interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
      scale_fill_gradient2(low = "blue",mid="blue",high = "red",trans = "log",
                           breaks=c(0,1,10,100),
                           na.value = "transparent", limits=c(1,maxODD))
    
  }
  
  pp<-arrangeGrob(pp1,pp2,nrow=1)
  rr<-arrangeGrob(rr1,rr2,nrow=1)
  
  ggsave(paste0(formatC(j, width = 2, format = "d", flag = "0"),"VUT_FBoutline.png"), 
         plot=pp,path = paste0(directory,'Plots/IIDIPUS_Results/VUT/FB'),width = 10,height = 4)
  ggsave(paste0(formatC(j, width = 2, format = "d", flag = "0"),"VUT_IIDIPUS.png"), 
         plot=rr,path = paste0(directory,'Plots/IIDIPUS_Results/VUT/IIDIPUS'),width = 10,height = 4)
  
}




