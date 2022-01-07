library(ggplot2)
library(sf)
library("ggmap")
library(OpenStreetMap)
library(osmdata)
library(tidyverse)
library(sp)
library(gstat)
library(raster)
source("./RCode/Functions.R")
library(gridExtra)
library(dplyr)
library(magrittr)

dir<-directory<-"./"

# Saint Vincent boundary box
bbox<-c(-61.29443,13.10000,-61.09926,13.40346)

mad_map <- get_stamenmap(bbox,source = "stamen",zoom=12)
p<-ggmap(mad_map)

# Extract all damaged buildings from copernicus estimates
coperny<-ExtractAllCOPERNICUS(dir,haz="VC")

# Extract all OSM buildings on the island & sum per region, 
# including average RWI, average ground-floor area, etc
OSMbuilds<-ExtractOSMbuild(bbox)

# Extract FB population data
FBpop<-raster(x="/home/patten/Documents/Coding/Oxford/IIDIPUS/Demography_Data/Population/FB/population_vct_2018-10-01.tif")
FBpop%<>%convRaster2SPDF(name = "Population")
e <- as(raster::extent(c(bbox[c(1,3,2,4)])), 'SpatialPolygons')
crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
FBpop%<>%raster::crop(e)

# Extract SEDAC population data
SEDACpop<-GetPopulationBbox(directory = dir,bbox = bbox,density = FALSE,yr = "2020")

# Get FB Data for Good displacement estimates
files<-GetFBEstimates(paste0(directory,"FBdata/VC20210409VCT/"), name="2331841688749824", plotty=F)
files$name<-files$baseline_polygon_name
files$day<-files$cdates-min(files$cdates,na.rm = T)
files<-files[!files$day%in%73:76,]

# Get Network Coverage
network<-GetFBEstimates("./FBdata/VC20210409VCT/Network_Coverage/", 
                        name="2058533107746849", plotty=F)
names(network)<-c("date","Latitude","Longitude","quadkey","count_2G","count_3G",
                  "count_4G","country","zoom","time")

# Evacuation centres
evac<-read_csv("./VC20210409VCT/20210414-SVG-NEMO-Shelters-via DirectRelief.csv")
evac$event<-coperny$event[1]
evac$grading<-"Evacuation Shelter"
evac<-evac[!is.na(evac$Longitude)|!is.na(evac$Latitude),]

# IDMC estimates
helix<-readRDS("./Helix/FullDatabase-06-05-2021.Rdata")
helix%<>%filter(hazard_type=="Volcanic eruption")
helix%<>%filter(event_id=="10019")

household<-read_csv("./VC20210409VCT/Av-Household-Size.csv")
nbuilds<-read_csv("./VC20210409VCT/NumBuildings.csv")
ltable<-read_csv("./VC20210409VCT/LifeTable.csv")
ltable$Age<-as.numeric(ltable$Age)

household%>%xtable()
nbuilds%>%xtable()
ltable%>%xtable(digits = 4)

modhouse<-data.frame()
for (i in 2:ncol(household)){
  x<-household[,i]
  modhouse%<>%rbind(data.frame(Name=household$Name,
                               Year=as.numeric(rep(colnames(x),nrow(x))),
                               HSize=unname(x)))
}
mody<-lm(HSize~Year + Name, modhouse)
summary(mody)
newdata<-data.frame(Year=rep(2021,nrow(household)),Name=household$Name)
modhouse%<>%rbind(cbind(newdata,data.frame(HSize=predict(mody,newdata))))
modhouse%>%filter(Year==2021)%>%xtable()


q<-ggplot(ltable,aes(Age,Population))+geom_line()+ylim(c(0,11000))
ggsave("LtablePop.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 5.)  
q<-ggplot(ltable,aes(Age,Deaths))+geom_line()+ylim(c(0,150))
ggsave("LtableDeaths.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 5.)  
#ggplot(ltable,aes(Age,cumsum(Deaths)))+geom_line()+ylim(c(0,1000))

# boundaries<-lapply(unique(files$baseline_polygon_name),function(x) 
#   c(getbb(paste0(x),format_out = 'sf_polygon')))

inputs<-c("Charlotte",
          "Saint Andrew",
          "Saint David",
          "Saint George",
          "Saint Patrick")
boundaries<-lapply(inputs,function(x) 
  c(getbb(paste0(x,", Saint Vincent"),format_out = 'sf_polygon')))

# All of Saint Vincent
# bndVincent<-getbb("Saint Vincent",format_out = 'sf_polygon')
# coordz<-st_coordinates(boundaries[[i]]$polygon)
# if(ncol(coordz)>3) {polysum<-coordz[rowMeans(coordz[,3:ncol(coordz)])==1,1:2]
# } else polysum<-coordz[1:2]

network$region<-coperny$region<-evac$region<-OSMbuilds$region<-FBpop@data$region<-SEDACpop@data$region<-NA

for (i in 1:length(boundaries)){
  # for (i in 3:5){
  if(is.null((boundaries[[i]])$polygon)) {
    tmp<-st_cast((boundaries[[i]])$multipolygon,"POLYGON")
    tmp$geometry<-tmp$geometry[2]
    boundaries[[i]]$polygon<-tmp
  }
  
  coordz<-st_coordinates(boundaries[[i]]$polygon)
  if(ncol(coordz)>3) {polysum<-coordz[rowMeans(coordz[,3:ncol(coordz)])==1,1:2]
  } else polysum<-coordz[1:2]
  
  # Allocate damaged buildings to each region boundary
  insiders<-point.in.polygon(coperny$Longitude,
                                   coperny$Latitude,
                                   polysum[,1],
                                   polysum[,2])>0
  coperny$region[insiders]<-inputs[i]
  
  # Allocate evacuation shelters to each region boundary
  insiders<-point.in.polygon(evac$Longitude,
                             evac$Latitude,
                             polysum[,1],
                             polysum[,2])>0
  evac$region[insiders]<-inputs[i]
  
  # Allocate general OSM buildings to each region boundary
  insiders<-point.in.polygon(OSMbuilds$Longitude,
                             OSMbuilds$Latitude,
                             polysum[,1],
                             polysum[,2])>0
  OSMbuilds$region[insiders]<-inputs[i]
  
  # Allocate FB population data to each region boundary
  notnan<-!is.na(FBpop$Population)
  insiders<-point.in.polygon(FBpop@coords[notnan,1],
                             FBpop@coords[notnan,2],
                             polysum[,1],
                             polysum[,2])>0
  regy<-FBpop$region[notnan]
  regy[insiders]<-inputs[i]
  FBpop$region[notnan]<-regy
  
  # Allocate SEDACs population data to each region boundary
  notnan<-!is.na(SEDACpop$Population)
  insiders<-point.in.polygon(SEDACpop@coords[notnan,1],
                             SEDACpop@coords[notnan,2],
                             polysum[,1],
                             polysum[,2])>0
  regy<-SEDACpop$region[notnan]
  regy[insiders]<-inputs[i]
  SEDACpop$region[notnan]<-regy

  insiders<-point.in.polygon(network$Longitude,
                             network$Latitude,
                             polysum[,1],
                             polysum[,2])>0
  network$region[insiders]<-inputs[i]
    
}
rm(regy,insiders,notnan)

# Extract population data (compare SEDACS to FB) & sum people per region
stmp<-SEDACpop%>%as.data.frame()%>%filter(!is.na(region))%>%group_by(region)%>%summarise(sumy=sum(Population,na.rm = T), percentage=sum(Population,na.rm = T)/99147)
stmp%>%xtable()
FBpop%>%as.data.frame()%>%filter(!is.na(region))%>%group_by(region)%>%summarise(sumy=sum(Population,na.rm = T), percentage=sum(Population,na.rm = T)/99147)%>%xtable()



coperny%>%group_by(region,grading)%>%summarise(nbuilds=length(Longitude))

coperny%>%filter(!is.na(region))%>%group_by(region)%>%summarise(nbuilds=length(Longitude))%>%xtable()

evac%>%filter(!is.na(region))%>%group_by(region)%>%summarise(nbuilds=length(Longitude))%>%xtable()

# var<-0.1*c(max(coperny$Longitude)-min(coperny$Longitude),
#            max(coperny$Latitude)-min(coperny$Latitude))
# bbox<-c(min(coperny$Longitude-var[1]),
#         min(coperny$Latitude-var[2]),
#         max(coperny$Longitude+var[1]),
#         max(coperny$Latitude+var[2]))
# bbox[2]<-13.10

coperny%>%group_by(region,grading)%>%
  summarise(nbuilds=length(Longitude))%>%xtable()

files%>%filter((status=="Displaced" | status=="Displaced Abroad") &
                 cdates=="2021-04-25")%>%
  group_by(baseline_polygon_name,cdates,gender)%>%
  summarise(sumy=sum(n_people,na.rm = T))%>%xtable()

files%>%filter((status=="Displaced Abroad") &
                 cdates=="2021-04-25")%>%
  group_by(baseline_polygon_name,cdates,gender)%>%
  summarise(sumy=sum(n_people,na.rm = T))%>%xtable()

files%>%filter((status=="Displaced" | status=="Displaced Abroad") &
                 cdates=="2021-07-20")%>%
  group_by(baseline_polygon_name,cdates,gender)%>%
  summarise(sumy=sum(n_people,na.rm = T))%>%xtable()

files%>%filter((status=="Displaced Abroad") &
                 cdates=="2021-07-20")%>%
  group_by(baseline_polygon_name,cdates,gender)%>%
  summarise(sumy=sum(n_people,na.rm = T))%>%xtable()

tmp<-files%>%filter(cdates=="2021-06-25" & (status=="Displaced" | status=="Displaced Abroad"))%>%
  dplyr::select(baseline_polygon_name,gender,n_people,percent_of_total_population_normalized_to_overall)%>%
  group_by(baseline_polygon_name,gender)%>%summarise(total=sum(n_people/percent_of_total_population_normalized_to_overall,na.rm=TRUE),
                                                     totaldisp=sum(n_people,na.rm=TRUE),
                                                     avpc=mean(percent_of_total_population_normalized_to_overall,na.rm=TRUE))
colnames(tmp)[1]<-"region"
tmp%<>%merge(stmp,by="region")
tmp%<>%dplyr::select(region,gender,total,sumy,totaldisp,avpc)
tmp$scaleddisp<-tmp$totaldisp*tmp$sumy/tmp$total
tmp%>%xtable()

files%>%filter(day==5)%>%dplyr::select(baseline_polygon_name,status,n_people,gender,percent_of_total_population)

files%>%filter(day==5 & gender=="overall" & (status=="Displaced" | status=="Displaced Abroad"))%>%
  dplyr::select(baseline_polygon_name,n_people,percent_of_total_population_normalized_to_overall)%>%
  group_by(baseline_polygon_name)%>%summarise(totaldisp=sum(n_people,na.rm=TRUE),
                                               avpc=mean(percent_of_total_population_normalized_to_overall,na.rm=TRUE),
                                               total=sum(n_people/percent_of_total_population_normalized_to_overall,na.rm=TRUE))

files%>%filter(day==80 & gender=="overall" & (status=="Displaced" | status=="Displaced Abroad"))%>%
  dplyr::select(baseline_polygon_name,n_people,percent_of_total_population_normalized_to_overall)%>%
  group_by(baseline_polygon_name)%>%summarise(totaldisp=sum(n_people,na.rm=TRUE),
                                               avpc=mean(percent_of_total_population_normalized_to_overall,na.rm=TRUE),
                                               total=sum(n_people/percent_of_total_population_normalized_to_overall,na.rm=TRUE))

# RWI & krig onto buildings that were damaged
RWI<-read.csv("./Demography_Data/RWI/vct_relative_wealth_index.csv")
# coordinates(RWI) <- ~ Longitude + Latitude
names(RWI)<-c("Latitude","Longitude","Intensity","error")
# Add extra points to ensure a smooth interpolation at the centre of the island
extras<-data.frame(Longitude=c(-61.22,-61.175,-61.175,-61.175),
                   Latitude=c(13.24,13.25,13.29,13.32),
                   Intensity=rep(min(RWI$Intensity,na.rm = T),4),
                   error=RWI$error[1:4])
RWI%<>%rbind(extras)

# Krig the RWI onto the copernicus building locations
copRWI<-KrigMeUp(RWI,coperny[,c("Longitude","Latitude")])
names(copRWI)<-c("RWI","Error")
coperny$RWI<-copRWI$RWI
rm(copRWI)
# Krig the RWI onto the evacuation shelter locations
evaRWI<-KrigMeUp(RWI,evac[,c("Longitude","Latitude")])
names(evaRWI)<-c("RWI","Error")
evac$RWI<-evaRWI$RWI
rm(evaRWI)

osmRWI<-KrigMeUp(RWI,OSMbuilds[,c("Longitude","Latitude")])
names(osmRWI)<-c("RWI","Error")
OSMbuilds$RWI<-osmRWI$RWI
rm(osmRWI)

popRWI<-KrigMeUp(RWI,as.data.frame(SEDACpop)[,c("Longitude","Latitude")])
names(popRWI)<-c("RWI","Error")
SEDACpop$RWI<-popRWI$RWI
rm(popRWI)

netRWI<-KrigMeUp(RWI,network[,c("Longitude","Latitude")])
names(netRWI)<-c("RWI","Error")
network$RWI<-netRWI$RWI
rm(netRWI)

pop4G<-apply(SEDACpop@coords,1,function(x){
  disty<-geosphere::distHaversine(x,network[,c("Longitude","Latitude")])
  return(weighted.mean(network$count_4G,1./(disty^2),na.rm=TRUE))
})
SEDACpop$count_4G<-pop4G
rm(pop4G)

osm4G<-apply(OSMbuilds[,c("Longitude","Latitude")],1,function(x){
  disty<-geosphere::distHaversine(x,SEDACpop@coords)
  return(weighted.mean(SEDACpop$count_4G,1./(disty^2),na.rm=TRUE))
})
OSMbuilds$count_4G<-osm4G
rm(osm4G)

cop4G<-apply(coperny[,c("Longitude","Latitude")],1,function(x){
  disty<-geosphere::distHaversine(x,SEDACpop@coords)
  return(weighted.mean(SEDACpop$count_4G,1./(disty^2),na.rm=TRUE))
})
coperny$count_4G<-cop4G
rm(cop4G)

pop4G<-apply(evac[,c("Longitude","Latitude")],1,function(x){
  disty<-geosphere::distHaversine(x,network[,c("Longitude","Latitude")])
  return(weighted.mean(network$count_4G,1./(disty^2),na.rm=TRUE))
})
evac$count_4G<-pop4G
rm(pop4G)

SEDACpop$Pop4G<-log10(SEDACpop$Population*SEDACpop$count_4G-SEDACpop$Population)
SEDACpop$Pop4G<-SEDACpop$Pop4G/max(SEDACpop$Pop4G,na.rm = T)
SEDACpop$Pop4G[is.na(SEDACpop$Pop4G)]<-min(SEDACpop$Pop4G,na.rm = T)

coperny$OSMtot<-NA
tmp<-OSMbuilds%>%group_by(region)%>%
  summarise(totals=length(Longitude))
coperny%<>%merge(tmp,by="region")

coperny%>%group_by(region,grading)%>%summarise(total=mean(totals),damaged=length(Longitude),pctdam=100*length(Longitude)/mean(totals))%>%xtable()

# stmp<-as.data.frame(SEDACpop["Pop4G"]); names(stmp)[1]<-"Intensity"
# ObsProb<-KrigMeUp(stmp[!is.na(stmp$Intensity),],as.data.frame(SEDACpop@coords))
# SEDACpop$Pop4G<-ObsProb$var1.pred

# pop4G<-apply(FBpop@coords,1,function(x){
#   disty<-geosphere::distHaversine(x,network[,c("Longitude","Latitude")])
#   return(weighted.mean(network$count_4G,1./(disty^2),na.rm=TRUE))
# })
# FBpop$count_4G<-pop4G
# rm(pop4G)

# cop4G<-apply(coperny[,c("Longitude","Latitude")],1,function(x){
#   disty<-geosphere::distHaversine(x,network[,c("Longitude","Latitude")])
#   return(weighted.mean(network$count_4G,1./(disty^2),na.rm=TRUE))
# })
# coperny$count_4G<-cop4G
# rm(cop4G)

q<-ggplot()+geom_sf(data=bndVincent$multipolygon,
                    inherit.aes = FALSE,show.legend = TRUE,size=3)+
  geom_tile(data=as.data.frame(SEDACpop),mapping = aes(Longitude,Latitude,fill=Population),inherit.aes = FALSE,alpha=0.8,na.rm = F);q
ggsave("SEDACpop.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)  

SEDACpop$RWI[is.na(SEDACpop$Population)]<-NA
q<-ggplot()+geom_sf(data=bndVincent$multipolygon,
                    inherit.aes = FALSE,show.legend = TRUE,size=3)+
  geom_tile(data=as.data.frame(SEDACpop),mapping = aes(Longitude,Latitude,fill=RWI),inherit.aes = FALSE,alpha=0.8,na.rm = F);q
ggsave("RWIkriged.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)  

q<-ggplot()+geom_sf(data=bndVincent$multipolygon,
                    inherit.aes = FALSE,show.legend = TRUE,size=3)+
  geom_tile(data=as.data.frame(tSED),mapping = aes(Longitude,Latitude,fill=count_4G),inherit.aes = FALSE,alpha=0.8,na.rm = F);q
ggsave("4Gkriged.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)  

mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "terrain",zoom=12)
p<-ggmap(mad_map)+ xlab("Longitude") + ylab("Latitude");p

names(RWI)[3]<-"RWI"
q<-p + geom_point(data=RWI,mapping = aes(Longitude,Latitude,colour=RWI),inherit.aes = FALSE,na.rm = F)+
  scale_color_gradient2(low = "purple4",mid = "purple4",high = "red");q
ggsave("RWIkpoints.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)  

q<-p + geom_point(data=network,mapping = aes(Longitude,Latitude,colour=count_4G),inherit.aes = FALSE,na.rm = F)+
  scale_color_gradient2(low = "purple4",mid = "purple4",high = "red");q
ggsave("4Gkpoints.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)  


q<-p+geom_tile(data=as.data.frame(FBpop),mapping = aes(Longitude,Latitude,fill=Population),inherit.aes = FALSE)
ggsave("FBpop.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)  

q<-p+geom_point(data=coperny,
             mapping=aes(Longitude,Latitude,colour=RWI),alpha=0.5) + 
  xlab("Longitude") + ylab("Latitude") + 
  scale_color_gradient2(low = "purple4",mid = "purple4",high = "red")
ggsave("RWI-Krig_Copernicus.png", plot=q,path = "./VC20210409VCT/",width = 6,height = 8.)

coperny%>%group_by(region,grading)%>%
  summarise(meanRWI=mean(RWI,na.rm = T),
            sdRWI=sd(RWI,na.rm = T),
            lenny=length(RWI))%>%xtable()

OSMbuilds%>%group_by(region)%>%
  summarise(lenny=length(RWI),
            meanArea=mean(area,na.rm=T),
            sdArea=sd(area,na.rm=T),
            meanRWI=mean(RWI,na.rm = T),
            sdRWI=sd(RWI,na.rm = T))%>%
  xtable()

OSMbuilds%>%group_by(region)%>%summarise(numy=length(area))%>%xtable()

q<-ggplot(filter(OSMbuilds,!is.na(region)),aes(area,group=region))+
  geom_density(aes(fill=region),alpha=0.5)+scale_x_log10() +
  xlab("Building Surface Area") + ylab("Density")
ggsave("OSM_BuildArea-Region.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)

q<-ggplot(filter(OSMbuilds,!is.na(region)),aes(RWI,group=region))+
  geom_density(aes(fill=region),alpha=0.5)+
  xlab("Relative Wealth Index") + ylab("Density")
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "fixed") + theme(strip.text.x = element_text(size = 12));q
ggsave("OSM_RWI-Region.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)

q<-ggplot(filter(OSMbuilds,!is.na(region)),aes(area,RWI,group=region))+geom_point(aes(colour=region))+scale_x_log10();q
ggsave("OSM_RWI-Area-Region.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)

q<-p+geom_point(data = filter(OSMbuilds,!is.na(region)),aes(Longitude,Latitude),colour="red");q
ggsave("OSMpoints.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)

q<-ggplot(filter(OSMbuilds,!is.na(region)),aes(area,RWI,group=region)) +scale_x_log10(limits=c(10,1000))+
  geom_density_2d_filled(contour_var = "ndensity") + ylab("Relative Wealth Index") + xlab("Area")
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "fixed") + theme(strip.text.x = element_text(size = 12));q
ggsave("OSM_RWI-Area-Region_2D-density.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)  



q<-ggplot(filter(coperny,!is.na(region)),aes(RWI,group=grading))+
  geom_density(aes(fill=grading),alpha=0.5)+ xlim(c(0,1))+
  xlab("Relative Wealth Index") + ylab("Density");q
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "free_y") + theme(strip.text.x = element_text(size = 12));q
ggsave("CoperRWIRegion.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 8,height = 5.)  


# ggsave("Copernicus_4G-Region_density.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)  


# q<-ggplot(filter(as.data.frame(FBpop),!is.na(region)),aes(Population,group=region))+
#   geom_density(aes(fill=region),alpha=0.5)+ ylim(c(0,5)) +
#   xlab("Population Density") + ylab("Density");q
# q<-q+facet_wrap(. ~ region,nrow = 2,scales = "free_y") + theme(strip.text.x = element_text(size = 12));q

q<-ggplot(filter(as.data.frame(SEDACpop),!is.na(region)),aes(Population,group=region))+scale_x_log10(limits=c(50,1000))+
  geom_density(aes(fill=region),alpha=0.5)+ #geom_histogram(aes(fill=region),alpha=0.5)+
  xlab("Population Density") + ylab("Density")
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "free_y") + theme(strip.text.x = element_text(size = 12));q
ggsave("SEDACpop-Region_density_nohisty.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)  

q<-ggplot(filter(as.data.frame(SEDACpop),!is.na(region)),aes(Population,RWI,group=region)) +scale_x_log10(limits=c(50,1000))+
  geom_density_2d_filled(contour_var = "ndensity") + ylab("Relative Wealth Index") + xlab("Population")
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "fixed") + theme(strip.text.x = element_text(size = 12));q
ggsave("SEDACpop_RWI-Region_2Ddensity.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)  

q<-ggplot(filter(as.data.frame(SEDACpop),!is.na(region)),aes(Pop4G,RWI))+
  #scale_x_log10(limits=c(1,30))+
  geom_density_2d_filled(bins=30,show.legend = F,contour_var = "ndensity") + xlab("4G Coverage Per Capita") + ylab("Relative Wealth Index");q
ggsave("RWI4G2Ddens.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 5.)  

q<-ggplot(filter(network,!is.na(region)),aes(count_4G,RWI,group=region)) + scale_x_log10(limits=c(1,100))+
  geom_density_2d_filled(contour_var = "ndensity") + ylab("Relative Wealth Index") + xlab("4G Count")
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "fixed") + theme(strip.text.x = element_text(size = 12));q
ggsave("RWI-4Gcount-Region_2Ddensity.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)  

q<-ggplot(filter(as.data.frame(SEDACpop),!is.na(region)),aes(Population,count_4G,group=region)) + scale_x_log10(limits=c(50,1000))+
  geom_density_2d_filled(contour_var = "ndensity") + xlab("Population Count") + ylab("4G Count")
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "free") + theme(strip.text.x = element_text(size = 12));q
ggsave("Pop-4Gcount-Region_2Ddensity.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)  

q<-ggplot(filter(as.data.frame(SEDACpop),!is.na(region)),aes(count_4G/log10(Population),RWI,group=region))+
  scale_x_log10(limits=c(1,30))+
  geom_density_2d_filled(contour_var = "ndensity") + xlab("4G Coverage Per Capita") + ylab("Relative Wealth Index")
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "fixed") + theme(strip.text.x = element_text(size = 12));q
ggsave("RWI-PerCap-4Gcount-Region_2Ddensity.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)  

q<-ggplot(filter(evac,!is.na(region)),aes(count_4G,group=region)) + scale_x_log10(limits=c(1,100))+
  geom_density(aes(fill=region),alpha=0.5) + ylab("Density") + xlab("4G Count")
q<-q+facet_wrap(. ~ region,nrow = 2,scales = "fixed") + theme(strip.text.x = element_text(size = 12));q
ggsave("evac-4Gcount-Region_2Ddensity.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 5.)  


subnames<-c("Charlotte","Saint Andrew","Saint Patrick")

q1<-ggplot(filter(coperny,!is.na(region) & region%in%subnames),aes(count_4G,group=region)) + scale_x_log10(limits=c(3,50))+
  geom_density(aes(fill=region),alpha=0.5) + ylab("Density") + xlab("4G Count")
q1<-q1+facet_wrap(. ~ region,nrow = 1,scales = "fixed") + theme(strip.text.x = element_text(size = 12))+
  ggtitle("Copernicus Buildings")+theme(plot.title = element_text(hjust = 0.5)) 

q2<-ggplot(filter(evac,!is.na(region) & region%in%subnames),aes(count_4G,group=region)) + scale_x_log10(limits=c(3,50))+
  geom_density(aes(fill=region),alpha=0.5) + ylab("Density") + xlab("4G Count")
q2<-q2+facet_wrap(. ~ region,nrow = 1,scales = "fixed") + theme(strip.text.x = element_text(size = 12))+
  ggtitle("Evacuation Shelters")+theme(plot.title = element_text(hjust = 0.5)) 

q3<-ggplot(filter(OSMbuilds,region%in%subnames),aes(count_4G,group=region))+
  geom_density(aes(fill=region),alpha=0.5)+ scale_x_log10(limits=c(3,50))+
  xlab("4G Coverage") + ylab("Density");
q3<-q3+facet_wrap(. ~ region,nrow = 1,scales = "free_y") + theme(strip.text.x = element_text(size = 12)) +
  ggtitle("OSM Buildings")+theme(plot.title = element_text(hjust = 0.5)) ;q3

q<-ggpubr::as_ggplot(gridExtra::arrangeGrob(q1,q2,nrow=2))
ggsave("evac-coper_4Gcount-Region_density.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 6.)

q<-ggpubr::as_ggplot(gridExtra::arrangeGrob(q1,q3,nrow=2))
ggsave("OSM-coper_4Gcount-Region_density.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 8,height = 6.)  


q1<-ggplot(filter(coperny,!is.na(region) & region%in%subnames),aes(RWI,group=region)) + 
  geom_density(aes(fill=region),alpha=0.5) + ylab("Density") + xlab("Relative Wealth Index")
q1<-q1+facet_wrap(. ~ region,nrow = 1,scales = "fixed") + theme(strip.text.x = element_text(size = 12))+
  ggtitle("Copernicus Buildings")+theme(plot.title = element_text(hjust = 0.5)) 

q2<-ggplot(filter(evac,!is.na(region) & region%in%subnames),aes(RWI,group=region)) + 
  geom_density(aes(fill=region),alpha=0.5) + ylab("Density") + xlab("Relative Wealth Index")
q2<-q2+facet_wrap(. ~ region,nrow = 1,scales = "fixed") + theme(strip.text.x = element_text(size = 12))+
  ggtitle("Evacuation Shelters")+theme(plot.title = element_text(hjust = 0.5)) 

q<-ggpubr::as_ggplot(gridExtra::arrangeGrob(q1,q2,nrow=2))

ggsave("evac-coper_RWI-Region_density.png", plot=q,path = "./VC20210409VCT/",width = 8,height = 6.)  


mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "terrain",zoom=12)
p<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")

q1<-p+geom_point(evac,mapping = aes(Longitude,Latitude),colour="purple") + ggtitle("Evacuation Shelters")+theme(plot.title = element_text(hjust = 0.5)) 
q2<-p+geom_point(coperny,mapping = aes(Longitude,Latitude),colour="green4")+ ggtitle("Buildings Damaged")+theme(plot.title = element_text(hjust = 0.5)) 
q<-ggpubr::as_ggplot(gridExtra::arrangeGrob(q1,q2,nrow=1));q
ggsave("EvacVSCoperny.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 7,height = 5.)  

q1<-p+geom_point(evac,mapping = aes(Longitude,Latitude,colour=RWI))+scale_color_gradient2(low = "purple4",mid = "purple4",high = "red")
ggsave("evacRWIpoint.png", plot=q1,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)  
q2<-p+geom_point(evac,mapping = aes(Longitude,Latitude,colour=count_4G))+scale_color_gradient2(low = "purple4",mid = "purple4",high = "red")
ggsave("evac4Gpoint.png", plot=q2,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)  

q<-ggpubr::as_ggplot(gridExtra::arrangeGrob(q1,q2,nrow=1))
ggsave("evacRWI4Gpoint.png", plot=q,path = "./VC20210409VCT/TexTables/Figures/",width = 6,height = 8.)

q<-p+geom_point(rbind(dplyr::select(evac,Longitude,Latitude,grading),dplyr::select(coperny,Longitude,Latitude,grading)),
                mapping = aes(Longitude,Latitude,colour=grading));q
ggsave("coperny-evac_point-mapping.png", plot=q,path = "./VC20210409VCT/",width = 6,height = 8.)

ggsave("RoR.png", plot=p,path = "./VC20210409VCT/TexTables/Figures/",width = 11,height = 5.)

q<-ggplot()+geom_sf(data = bndVincent$multipolygon,size=2)
for(i in 1:length(boundaries)){
  q<-q+geom_sf(data = boundaries[[i]]$polygon,size=1)
}
q<-q+geom_contour_filled(data = as.data.frame(SEDACpop),
                     mapping = aes(Longitude,Latitude,z=Pop4G),alpha=0.7,bins = 30,show.legend = F,na.rm = T);q
ggsave("ObsProb_4GcountPop_2Ddensity.png", plot=q,path = "./VC20210409VCT/",width = 6,height = 8.)  

p+geom_contour_filled(data = as.data.frame(SEDACpop),
                      mapping = aes(Longitude,Latitude,z=1/(1+Pop4G)),alpha=0.7,bins=30)

mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "terrain",zoom=12)
qq<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")





m1<-wilcox.test( ~ | grading, data=coperny, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(m1)






maxy<-files%>%filter(gender!="overall" & (status=="Displaced"|status=="Displaced Abroad"))%>%
  group_by(name,cdates,gender)%>%
  summarise(sumy=sum(n_people,na.rm = T))%>%pull(sumy)%>%max()

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
sc<- scale_fill_gradientn(colours = myPalette(100), limits=c(0, maxy))

draw.curve<-function(j){
  
  if(j>86) {
    
    pp1<-qq+geom_sf(data=bndVincent$multipolygon,
                     inherit.aes = FALSE,show.legend = TRUE)+
      ggtitle("")
    pp<-ggpubr::as_ggplot(gridExtra::arrangeGrob(pp1,pp1,nrow=1))+
      ggtitle(paste0("Displaced Population on ",max(files$cdates)+(j-86)))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(pp)
    
  # } else if(j%in%73:76) {
  #   
  #   geosub<-files%>%filter(day==72)
  #   
  #   # pp1<-qq+geom_sf(data=bndVincent$multipolygon,
  #                   # inherit.aes = FALSE,show.legend = TRUE)+
  #     # ggtitle("")
  #   pp<-ggpubr::as_ggplot(gridExtra::arrangeGrob(pp1,pp2,nrow=1)) + 
  #     ggtitle(paste0("Displaced Population on ",geosub$cdates[1]))+theme(plot.title = element_text(hjust = 0.5)) 
  #   
  #   print(pp)
    
  } else {
    
    if(j%in%73:76) {
      geosub<-files%>%filter(day==72 & gender!="overall" &
                             (status=="Displaced"|status=="Displaced Abroad"))
      geosub$cdates<-geosub$cdates+(j-72)
    } else {
      geosub<-files%>%filter(day==j & gender!="overall" &
                               (status=="Displaced"|status=="Displaced Abroad"))
    }
      
    printers<-geosub%>%group_by(gender)%>%summarise(sumy=sum(n_people,na.rm = T),date=unique(cdates),.groups='drop_last')
    # normsub<-filter(normy,day==j)
    
    pp1<-qq + ggtitle(paste0("Displaced Women = ",printers$sumy[printers$gender=="female"])) +
      theme(plot.title = element_text(hjust = 0.5)) 
    pp2<-qq + ggtitle(paste0("Displaced Men = ",printers$sumy[printers$gender=="male"])) + 
      theme(plot.title = element_text(hjust = 0.5)) 
    
    for(i in 1:length(boundaries)){
      nn<-inputs[i]
      colz<-geosub%>%filter(name==nn)
      fdec<-pull(filter(colz,gender=="female"),n_people)%>%sum(na.rm = T)
      mdec<-pull(filter(colz,gender=="male"),n_people)%>%sum(na.rm = T)
      
      fval<-(fdec/maxy + 0.5)/1.5
      pp1<-pp1+geom_sf(data=boundaries[[i]]$polygon,fill=rgb(fval,0,(1-fval)),
                       # alpha=0.5,
                       inherit.aes = FALSE,show.legend = TRUE)
      fval<-(mdec/maxy + 0.5)/1.5
      pp2<-pp2+geom_sf(data=boundaries[[i]]$polygon,fill=rgb(fval,0,(1-fval)),
                       # alpha=0.5,
                       inherit.aes = FALSE,show.legend = TRUE)
      
      # print(max(fdec,mdec))
      
    }
    
    pp1<-pp1+labs(fill = paste0("Disp. Pop. ",geosub$cdates[1])) 
    pp2<-pp2+labs(fill = paste0("Disp. Pop. ",geosub$cdates[1])) 
    
    pp<-ggpubr::as_ggplot(gridExtra::arrangeGrob(pp1,pp2,nrow=1))+
      ggtitle(paste0("Displaced Population on ",geosub$cdates[1]))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(pp)
    # ggsave(paste0(formatC(j, width = 2, format = "d", flag = "0"),"VCT_FBoutline.png"), 
    #        plot=pp,path = './VC20210409VCT/FBplots/',width = 10,height = 4)
  }
}

draw.curve(1)
draw.curve(74)
draw.curve(99)


trace.animate <- function() {
  # lapply(unique(files$day), function(i) {
  lapply(1:100, function(i) {
    draw.curve(i)
  })
}
#save all iterations into one GIF
animation::saveGIF(trace.animate(), interval = .2, movie.name="trace_looped.gif", ani.height=380)

trace.animate <- function() {
  lapply(unique(files$day), function(i) {
    # lapply(1:100, function(i) {
    draw.curve(i)
  })
}
#save all iterations into one GIF
animation::saveGIF(trace.animate(), interval = .2, movie.name="trace.gif",loop=1, ani.height=380)







# colsy<-colorRampPalette(c("blue", "red"))( length(boundaries) )
colsy<-c("red","blue","purple","orange","green")
# colsy<-colorRampPalette(c("blue", "red"))( 3 )
for (i in 1:length(boundaries)){
  # for (i in 3:5){
  p<-p+geom_sf(data=(boundaries[[i]])$polygon,
               inherit.aes = FALSE,
               alpha=0.1,fill="grey",size=2)
}
p<-p+xlab("Longitude") + ylab("Latitude");p
ggsave("RegionBoundaries.png", plot=p,path = "./VC20210409VCT/",width = 6,height = 8.)  

# "Charlotte"     "Saint Andrew"  "Saint David"   "Saint George"  "Saint Patrick"







