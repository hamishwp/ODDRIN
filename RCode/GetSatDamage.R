library(sf)
library(ggplot2)
library("ggmap")
library(OpenStreetMap)
library(osmdata)
source("RCode/Functions.R")
library(gridExtra)
library(magrittr)

ExtractAllCOPERNICUS<-function(dir,haz="EQ"){
  
  print("WARNING: only point files and not polygons are currently read from COPERNICUS damage assessment")
  cfiles<-list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="crisis_information",recursive = T)
  cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="builtUpP",recursive = T))
  cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="built_up_p",recursive = T))
  cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="settlements_po",recursive = T))
  cfiles<-cfiles[grep(haz,cfiles)]; cfiles<-cfiles[grep(".shp",cfiles)] ; cfiles<-cfiles[grep(".shp.xml",cfiles,invert = T)]
  tmp<-strsplit(cfiles,"/",fixed = T) ; evname<-c()
  for (i in 1:length(tmp)) {
    evname%<>%c(tmp[[i]][1]) 
    cfiles[i]<-paste0(dir,"COPERNICUS_Damage/",cfiles[i])
  }
  
  print(paste0(length(unique(evname))," COPERNICUS building damage assessment files being loaded"))
  Damage<-data.frame()
  for(ev in unique(evname)){
    subfiles<-unique(grep(ev,cfiles,value = T,fixed = T))
    # I know that you can load all in at the same time, but let's save some RAM!!!
    for (f in unique(subfiles)){
      if(grepl(".shp.xml",f)) next
      tmp<-st_read(f,quiet = T)
      if((is.null(tmp$grading) & is.null(tmp$damage_gra))|is.null(tmp$geometry)) next
      if(!is.null(tmp$grading)) {tmp%<>%dplyr::select(grading,geometry)
      } else if (!is.null(tmp$damage_gra)) {
        tmp%<>%dplyr::select(damage_gra,geometry)
        colnames(tmp)<-c("grading","geometry")
      }
      if(!all(class(tmp$geometry)%in%c("sfc_POINT","sfc"))) tmp$geometry%<>%st_centroid
      tmp$geometry%<>%st_transform(crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
      # Remove unwanted grading terms
      tmp%<>%filter(!grading%in%
                         c("Null","Unknown","Not Applicable"))
      if(nrow(tmp)==0) next
      # convert point geometries into Longitude/Latitude
      ttt<-st_coordinates(tmp$geometry)
      tmp$Longitude<-ttt[,1]; tmp$Latitude<-ttt[,2]
      tmp$geometry<-NULL
      tmp$event<-rep(ev,nrow(tmp))
      Damage%<>%rbind(tmp)
      
    }
  }
  
  # Some more cleaning to do... ugly, right?
  Damage$grading[grepl("Completely Destroyed",Damage$grading,fixed = T) |
                 grepl("Destroyed",Damage$grading,fixed = T)]       <-"Completely Destroyed"
  Damage$grading[grepl("Highly Damaged",Damage$grading,fixed = T)]             <-"Highly Damaged"
  Damage$grading[grepl("Moderately Damaged",Damage$grading,fixed = T)]         <-"Moderately Damaged"
  Damage$grading[grepl("Negligible to slight damage",Damage$grading,fixed = T)|
                 grepl("Possibly damaged",Damage$grading,fixed = T)]<-"Negligible to slight damage"
  Damage$grading[grepl("Not Affected",Damage$grading,fixed = T) |
                 grepl("No visible damage",Damage$grading,fixed = T) ]               <-"Not Affected"
  
  Damage$grading%<>%as.factor()
  Damage%<>%droplevels()
  
  # levels(Damage$grading)<-unique(c(levels(Damage$grading),"Completely Destroyed","Highly Damaged",
  #                                  "Moderately Damaged","Negligible to slight damage",
  #                                  "Not Affected"))
  
  # Remove all other terms
  Damage%<>%filter(grading%in%
                  c("Completely Destroyed","Negligible to slight damage",
                    "Not Affected","Moderately Damaged","Highly Damaged","Damaged"))
  Damage$Confidence<-rep(NA,nrow(Damage))
  Damage<-Damage[!duplicated(Damage[,c("Longitude","Latitude","event")]),]
  
  print(table(Damage$grading))
  
  return(Damage)
  
}

ExtractAllUNOSAT<-function(dir,haz="EQ"){
  
  ufiles<-list.files(path=paste0(dir,"UNOSAT_Damage/"),pattern="damage_site",recursive = T,ignore.case = T)
  ufiles<-c(ufiles,list.files(path=paste0(dir,"UNOSAT_Damage/"),pattern="Damagepoint",recursive = T,ignore.case = T))
  ufiles<-c(ufiles,list.files(path=paste0(dir,"UNOSAT_Damage/"),pattern="damaged_structure",recursive = T,ignore.case = T))
  ufiles<-c(ufiles,list.files(path=paste0(dir,"UNOSAT_Damage/"),pattern="damage_points",recursive = T,ignore.case = T))
  ufiles<-c(ufiles,list.files(path=paste0(dir,"UNOSAT_Damage/"),pattern="damagesite",recursive = T,ignore.case = T))
  ufiles<-c(ufiles,list.files(path=paste0(dir,"UNOSAT_Damage/"),pattern="damagesettlement",recursive = T,ignore.case = T))
  ufiles<-c(ufiles,list.files(path=paste0(dir,"UNOSAT_Damage/"),pattern="damagedbuilding",recursive = T,ignore.case = T))
  ufiles<-ufiles[grep(haz,ufiles)]; ufiles<-ufiles[grep(".shp",ufiles)]; ufiles<-ufiles[grep(".shp.xml",ufiles,invert = T)]
  ufiles<-ufiles[endsWith(ufiles, ".shp")]
  
  tmp<-strsplit(ufiles,"/",fixed = T) ; evname<-c()
  for (i in 1:length(tmp)) {
    evname%<>%c(strsplit(tmp[[i]][1],"_",fixed = T)[[1]][1]) 
    ufiles[i]<-paste0(dir,"UNOSAT_Damage/",ufiles[i])
  }
  
  print(paste0(length(unique(evname))," UNOSAT building damage assessment files being loaded"))
  Damage<-data.frame()
  for(ev in unique(evname)){
    subfiles<-unique(grep(ev,ufiles,value = T,fixed = T))
    # I know that you can load all in at the same time, but let's save some RAM!!!
    for (f in unique(subfiles)){
      tmp<-st_read(f,quiet = T)
      if(is.null(tmp$geometry)|is.null(tmp$Grouped_Da)) next
      # Unify names with COPERNICUS dataframe:
      if(!is.null(tmp$DamageLeve)) {tmp$grading<-tmp$DamageLeve; tmp$DamageLeve<-NULL}
      else if(!is.null(tmp$Main_Damag)) {tmp$grading<-tmp$Main_Damag; tmp$Main_Damag<-NULL}
      else if(!is.null(tmp$Main_Dmg)) {tmp$grading<-tmp$Main_Dmg; tmp$Main_Dmg<-NULL}
      else {print(paste0("no damage data for: ",f));next}
      # Remove unwanted columns
      tmp%<>%dplyr::select(Confidence,grading,geometry)
      if(!all(class(tmp$geometry)%in%c("sfc_POINT","sfc"))) tmp$geometry%<>%st_centroid
      tmp$geometry%<>%st_transform(crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
      if(nrow(tmp)==0) next
      # convert point geometries into Longitude/Latitude
      ttt<-st_coordinates(tmp$geometry)
      tmp$Longitude<-ttt[,1]; tmp$Latitude<-ttt[,2]
      tmp$geometry<-NULL
      tmp$event<-rep(ev,nrow(tmp))
      Damage%<>%rbind(tmp)
      
    }
  }
  
  table(Damage$grading)
  return(Damage)
  
}

ExtractBDfiles<-function(dir="./",haz="EQ",org="all",rename=T){
  
  # COPERNICUS
  if(org%in%c("all","c"))       CopDamage<-ExtractAllCOPERNICUS(dir,haz)
  if(org%in%c("all","UNOSAT"))  UNODamage<-ExtractAllUNOSAT(dir,haz)
  # Merge dataframes together
  if(org=="all") {Damage<-rbind(CopDamage,UNODamage); rm(CopDamage,UNODamage)
  }else if(org=="c") {Damage<-CopDamage ; rm(CopDamage)
  }else {Damage<-UNODamage ; rm(UNODamage)}
  
  Damage$grading%<>%as.character()
  
  # Merge all naming conventions into one:
  Damage$grading[grep("Destroyed",Damage$grading,fixed = T)]                  <-"Completely Destroyed"
  Damage$grading[grep("Highly Damaged",Damage$grading,fixed = T)]             <-"Severe Damage"
  Damage$grading[grep("Moderately Damaged",Damage$grading,fixed = T)]         <-"Moderate Damage"
  Damage$grading[grep("Negligible to slight damage",Damage$grading,fixed = T)]<-"Possible Damage"
  Damage$grading[grep("No Visible Damage",Damage$grading,fixed = T)]          <-"Not Affected"
  
  # CHANGE NAMES FROM THE ORIGINAL EMS NAMING CONVENTION
  if(rename){
    Damage$grading[Damage$grading=="Completely Destroyed"]<-"destroyed"
    Damage$grading[Damage$grading=="Severe Damage"]<-"severe"
    Damage$grading[Damage$grading=="Moderate Damage"]<-"moderate"
    Damage$grading[Damage$grading=="Possible Damage"]<-"possible"
    Damage$grading[Damage$grading=="Not Affected"]<-"notaffected"
  }
  
  # Reduce name into haz,date,iso:
  Damage$hazard<-haz
  # Get dates first
  Damage$sdate<-as.Date(substring(Damage$event,3,10),tryFormats="%Y%m%d")
  Damage$iso3<-substring(Damage$event,11,13)
  # Damage$grading%<>%droplevels
  
  Damage%<>%distinct(.keep_all = T)
  
  return(Damage)
  
}

PlotBetaProbs<-function(){
  
  haz<-seq(0,1,length.out = 1000)
  # No damage
  plot(100*haz,dbeta(haz,0.065,100,log=T),ylim=c(-30,3),ylab="Log Probability (beta)",xlab="Damage %",type="l",lwd=4)  
  # Minor damage
  lines(100*haz,dbeta(haz,3,100,log=T),col="blue",lwd=4)
  # Major damage
  lines(100*haz,dbeta(haz,6,16,log = T),col="orange",lwd=4)  
  # Totally destroyed
  lines(100*haz,dbeta(haz,6,2,log = T),col="red",lwd=4)  
  # Add me a legend!
  legend(x = 30,y=-10,legend = c("No Damage", "Minor Damage","Major Damage", "Totally Destroyed"),
         col = c("black","blue","orange","red"),lty = 1,lwd = 4)
  
}

plotSAT<-function(Damage,country=NULL,zoomy=11){
  if(!is.null(country)) bbox<-c(getbb(country))
  else {
    var<-0.1*c(max(Damage$Longitude)-min(Damage$Longitude),max(Damage$Latitude)-min(Damage$Latitude))
    bbox<-c(min(Damage$Longitude-var[1]),min(Damage$Latitude-var[2]),max(Damage$Longitude+var[1]),max(Damage$Latitude+var[2]))
  }
  mad_map <- get_stamenmap(bbox,source = "stamen",zoom=zoomy)
  p<-ggmap(mad_map, maprange=FALSE) + coord_sf(crs = st_crs(4326))
  p<-p + geom_point(data=Damage,aes(Longitude,Latitude,colour=grading),size=2, inherit.aes = FALSE) +
    xlab("Longitude") + ylab("Latitude") +
    theme(plot.title=element_text(hjust=0.5));p
}

# 
# GetUNOSATDamage<-function(filer,plotty=T,titlz="UNOSAT Damage Assessment"){
#   # Damage<-sf::st_read(paste0(directory,"UNOSAT_Damage/all-damage-nepal-april-28th-2015/All_Damage_Nepal_April_28th_2015.shp"))
#   Damage<-st_read(filer)
#   Damage%<>%st_transform(crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
#   
#   if(plotty){
#     bbox<-unname(CheckBbox(st_bbox(Damage)))
#     
#     # bbox<-c(85,27.45,85.6,28)
#     mad_map <- get_stamenmap(bbox,source = "stamen",zoom=11)
#     
#     # names(bbox)<-c("xmin","ymin","xmax","ymax")
#     # Damage%<>%st_crop(bbox)
#     
#     p<-ggmap(mad_map, maprange=FALSE) + coord_sf(crs = st_crs(4326))
#     p<-p + geom_sf(data=Damage,aes(colour=Damage),size=2, inherit.aes = FALSE) +
#       xlab("Longitude") + ylab("Latitude") + ggtitle(titlz) +
#       theme(plot.title=element_text(hjust=0.5))
#     # if(!is.null(shakeid)){
#     #   # tp<-ShakeURL2Poly(1052901)
#     #   tp<-ShakeURL2Poly(shakeid)
#     #   tp%<>%filter(Longitude<bbox[3] & Longitude>bbox[1] & Latitude<bbox[4] & Latitude>bbox[2])
#     #   for (j in unique(tp$ncontour)){
#     #     p<-p+geom_polygon(data = filter(tp,ncontour==j),aes(x=Longitude,y=Latitude,group=Intensity,colour=Intensity),
#     #                       alpha=0,na.rm = T,size=2, inherit.aes = FALSE)
#     #   }
#     # }
#     p
#     ggsave(paste0(titlz,".png"), plot=p,path = paste0(directory,'/Plots/UNOSAT/'),width = 9,height = 7.)
#   }
#   
#   return(Damage)
# }
# 
# GetCOPERNICUSDamage<-function(directory,plotty=T){
#   
#   fff<-list.files(path=paste0(directory,"COPERNICUS_Damage/"),pattern="crisis_information_point_grading.shp",recursive = T)
#   
#   for (filer in fff){
#     Damage<-st_read(paste0(directory,"COPERNICUS_Damage/",filer))
#     # Damage<-st_read("./COPERNICUS_Damage/EMSR190_33SELLANO_GRADING_OVERVIEW_v1_vector/EMSR190_33SELLANO_02GRADING_v1_15000_settlements_poly_grading.shp")
#     Damage%<>%st_transform(crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
#     # st_combine(tmp,Damage)
#   }
#   
#   # bbox<-unname(CheckBbox(st_bbox(Damage)))
#   # ggplot() + geom_sf(data=Damage,aes(colour=grading),size=2, inherit.aes = FALSE)
#   
# }