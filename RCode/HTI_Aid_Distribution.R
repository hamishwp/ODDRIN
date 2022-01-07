aid<-read_csv("./IIDIPUS_Input/Distribution de l'aide humanitaire_SeismeAout21_030921.xlsx - Distribution_AH.csv")

batches<-aid%>%filter(Articles%in%c("bache","Bache","Baches"))%>%group_by(Commune)%>%summarise(batchsum=sum(Quantite,na.rm = T))
food<-aid%>%filter(str_to_lower(Articles)=="kit alimentaire")%>%group_by(Commune)%>%summarise(foodsum=sum(Quantite,na.rm = T))

# Extract commune boundary polygons from OSM
# Or use https://download.geofabrik.de/central-america/haiti-and-domrep.html
# In combination with https://rpubs.com/JesperHybel/677164
# mydata_poly <- oe_read("https://download.geofabrik.de/central-america/haiti-and-domrep-latest.osm.pbf",layer="multipolygons")
mydata_poly<-readRDS("./IIDIPUS_Input/HTI_Polygons.Rdata")
# boundaries<-mydata_poly%>%filter(!is.na(admin_level)& admin_level==8)%>%dplyr::select(name,geometry)
boundaries<-mydata_poly%>%filter(!is.na(admin_level)& admin_level==4)%>%dplyr::select(name,geometry)
# find the region boundaries that were affected directly by the earthquake
inds<-st_centroid(boundaries$geometry)%>%st_coordinates()%>%as.data.frame()%>%
      transmute(indy=(X< -73 & Y<18.7))%>%pull(indy)
unique(mydata_poly$admin_level)
boundaries<-boundaries[inds,]
# arrange(cbind(as.data.frame(st_coordinates(st_centroid(tmp2))),tmp2$name),desc(Y))
plot(boundaries)

prepnames<-grep("commune ",boundaries$name,ignore.case = T,value = T)%>% strsplit("Commune ")%>%unlist()
boundaries$name[grep("commune ",boundaries$name,ignore.case = T)] <- prepnames[seq.int(from=2,to=length(prepnames),by = 2)]
boundaries$name[startsWith(boundaries$name,"de ")]<-gsub("de ","",boundaries$name[startsWith(boundaries$name,"de ")])
boundaries$name[startsWith(boundaries$name,"du ")]<-gsub("du ","",boundaries$name[startsWith(boundaries$name,"du ")])
boundaries$name[startsWith(boundaries$name,"des ")]<-gsub("des ","",boundaries$name[startsWith(boundaries$name,"des ")])
boundaries$name[startsWith(boundaries$name,"d’")]<-gsub("d’","L'",boundaries$name[startsWith(boundaries$name,"d’")])
unique(boundaries$name)

tmp<-st_coordinates(st_centroid(boundaries$geometry))
boundaries$Longitude<-tmp[,1] ; boundaries$Latitude<-tmp[,2]

boundaries%<>%arrange(name)

# Check it out!
funcy<-function(i) {
  tmp2<-boundaries
  tmp2[,1]<-NA
  tmp2[i,1]<-10
  plot(tmp2[,1])
}
funcy(6)

# couple shapefiles to aid commune name
ids<-c(22,6,
)
View(batches)
View(food)

# Well this is shit! Wanted to make sure it worked...
# tmp$batches <- c(NA,410,NA,200,NA,167,342,NA,862,629,NA,60,205,NA,NA,NA,NA,NA,NA,401,NA,4181,NA,1809,NA,1095,NA,5,NA,NA,1100,400,546,100,100,NA,150,NA,NA,NA,NA,815) #les cayemites

plot(boundaries[,2:3])
plot(boundaries[,c(2,4)])
plot(boundaries[,-1])
     





ODDy<-readRDS("./IIDIPUS_Results/ODDobjects_V2/EQ20210814HTI_10919")

osm<-readRDS("./IIDIPUS_Input/EQ20210814HTI_OSMbuildings_more.Rdata")
osm%<>%distinct()
# convert osm to a SpatialPointsDataFrame
array <- SpatialPointsDataFrame(coords = osm[c("Longitude","Latitude")],
                                data = osm["area"]) # don't worry about area, we won't use it 
# crs(array)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
array$Population<-array$area; array$area<-NULL

ODDy<-ParAggFBPopSEDAC(ODDy,array,ncores = 2,funcy = "length",namer = "NumBDosm")


sum(ODDy@data$NumBDFB[ODDy@coords[,1]<modbbox[3] & ODDy@coords[,2]<modbbox[4]],na.rm = T)
sum(ODDy@data$NumBDosm[ODDy@coords[,1]<modbbox[3] & ODDy@coords[,2]<modbbox[4]],na.rm = T)

ODDy@modifier[1]<-0.0

Omega<-readRDS("./IIDIPUS_Results/Omega_v2_20210828.Rdata")
Omega_save<-Omega
Omega$Lambda$omega<-Omega$Lambda$kappa<-0
Omega$Lambda$nu<-1

ODDy@modifier<-list(HTI=0.05,DOM=0,CUB=0,USG=0)

ODDy@modifier[1]<-0
ODDy$Population<-ODDy$FBPop
tLL<-DispX(ODD = ODDy,Omega = Omega_save,center = Model$center,LL = F,Method = AlgoParams)
sum(tLL$Disp,na.rm = T)

ODDy@modifier[1]<-0.1
tLL<-DispX(ODD = ODDy,Omega = Omega_save,center = Model$center,LL = F,Method = AlgoParams)
sum(tLL$Disp,na.rm = T)


ODDy$Population<-ODDy$SEDACpop
tLL<-DispX(ODD = ODDy,Omega = Omega_save,center = Model$center,LL = F,Method = AlgoParams)
sum(tLL$Disp,na.rm = T)

ODDy$Population<-ODDy$NumBDosm
tLL<-DispX(ODD = ODDy,Omega = Omega,center = Model$center,LL = F,Method = AlgoParams)
sum(tLL$Disp,na.rm = T)
ODDy$DamBDosm<-tLL$Disp
ODDy$Population<-ODDy$NumBDFB
tLL<-DispX(ODD = ODDy,Omega = Omega,center = Model$center,LL = F,Method = AlgoParams)
sum(tLL$Disp,na.rm = T)
ODDy$DamBDFB<-tLL$Disp

ODDy<-ODDyAggPerRegion(ODDy,feature="DamBDosm",name="RegionAggBDosm")
ODDy<-ODDyAggPerRegion(ODDy,feature="DamBDFB",name="RegionAggBDFB")
ODDy<-ODDyAggPerRegion(ODDy,feature="ExpDisp",name="RegionAggDisp")

plotODDyBG()

ODDy@data$Department[ODDy@data$Department=="la Grande-Anse"]<-"La Grande-Anse"
ODDy@data$Department[ODDy@data$Department=="l'Artibonite"]<-"L'Artibonite"
ODDy@data$Department[ODDy@data$Department=="l'Ouest"]<-"L'Ouest"

geoinsights<-GetFBEstimates('/home/patten/Documents/Coding/Oxford/IIDIPUS/FBdata/EQ20201408HTI/Displacement/', "", reference="", check=FALSE, plotty=FALSE)
geoinsights%>%group_by(baseline_polygon_name,day)%>%summarise(disp=sum(n_people))
geoinsights%>%filter(day<40)%>%group_by(baseline_polygon_name)%>%summarise(disp=sum(n_people),
                              maxdisp_norm=max(percent_of_total_population_normalized_to_overall),
                              maxdisp=max(percent_of_total_population))

geoinsights$day<-geoinsights$day+17

former<-geoinsights%>%filter(status%in%c("Displaced","Displaced Abroad"))%>%
  group_by(baseline_polygon_name,day)%>%summarise(
  maxdisp_norm=max(percent_of_total_population_normalized_to_overall),
  maxdisp=max(percent_of_total_population),.groups="drop")
names(former)[1]<-"Department"

for (dept in unique(former$Department)){
  tmp<-former%>%filter(Department==dept)
  for(j in 1:nrow(tmp)) {
    tmp$maxdisp_norm[j]<-max(tmp$maxdisp_norm[j:nrow(tmp)],na.rm = T)
    tmp$maxdisp[j]<-max(tmp$maxdisp[j:nrow(tmp)],na.rm = T)
  }
  former$maxdisp_norm[former$Department==dept]<-tmp$maxdisp_norm
  former$maxdisp[former$Department==dept]<-tmp$maxdisp
}

normtable<-ODDy@data%>%group_by(Department)%>%summarise(pop=sum(FBPop,na.rm = T),.groups="drop")

former%<>%plyr::join(normtable,by="Department")
former$PopDisp<-ceiling(former$maxdisp_norm*former$pop)

former%>%group_by(Department)%>%summarise(minDisp=min(PopDisp),maxDisp=max(PopDisp))

Disper<-ODDy@data%>%group_by(Department)%>%summarise(DispODDRIN=sum(ExpDisp,na.rm = T),
                                                     DamFBODDRIN=sum(DamBDFB,na.rm = T),
                                                     DamOSMODDRIN=sum(DamBDosm,na.rm = T),.groups="drop")
former%<>%plyr::join(Disper,by="Department")

former%>%group_by(Department)%>%summarise(minDisp=min(PopDisp),maxDisp=max(PopDisp),
                                          DispODDRIN=max(DispODDRIN),
                                          DamFBODDRIN=max(DamFBODDRIN),
                                          DamOSMODDRIN=max(DamOSMODDRIN))%>%
  write.csv("./IIDIPUS_Results/FB_ODDRIN_Results_EQ20210814HTI.csv")

DTM<-read_csv("./IIDIPUS_Input/DTM_EQ20210814HTI_Masterlist_27082021.csv")

redDTM<-DTM%>%group_by(Département)%>%summarise(Menages=sum(Ménages,na.rm = T),Individus=sum(Individus,na.rm = T),.groups="drop")
redDTM$DTMpercInd<-redDTM$Individus/sum(redDTM$Individus)

subform<-former%>%filter(Department%in%c("La Grande-Anse","Nippes","Sud") & day==0)#%>%
  # summarise(Department=Department,maxdisp_norm=maxdisp_norm/sum(maxdisp_norm),
  #           maxdisp=maxdisp/sum(maxdisp),
  #           pop=pop/sum(pop),
  #           PopDisp=PopDisp/sum(PopDisp),
  #           DispODDRIN=DispODDRIN/sum(DispODDRIN),
  #           DamFBODDRIN=DamFBODDRIN/sum(DamFBODDRIN),
  #           DamOSMODDRIN=DamOSMODDRIN/sum(DamOSMODDRIN)
  #           )
# subform%<>%cbind(redDTM[,"DTMpercInd"])
subform%<>%cbind(redDTM[,"Individus"])
subform%<>%dplyr::select(-c("day","pop","maxdisp_norm","maxdisp"))
names(subform)[2]<-"FB_Disp_XDR"

subform%<>%cbind(data.frame(HTI_PC_Destroyed=c(25892,14989,42889),HTI_PC_Damaged=c(8648,14450,30717)))
write.csv(subform,"./IIDIPUS_Results/ALL_obs-pred_Results_EQ20210814HTI.csv")

subgeo<-former%>%filter(Department%in%c("La Grande-Anse","Nippes","Sud"))#%>%

p<-ggplot(subgeo,aes(day,PopDisp,group=Department))+geom_line(aes(colour=Department))+
  xlab("Days Since Earthquake")+ylab("Displaced Population") + scale_y_log10(limits =c(0,200000))
  p<-p+facet_wrap( ~ Department, scales = "fixed")+ theme(plot.title = element_text(hjust = 0.5));p
  

subgeo$maxdisp_norm  
normaliser<-subgeo%>%group_by(Department)%>%summarise(normy=max(maxdisp_norm,na.rm = T),.groups="drop")  
# subgeo%<>%group_by(Department)%>%mutate(normDisp=maxdisp_norm/normaliser$normy[normaliser$Department==Department])

subgeo$normDisp<-1
for(stringy in c("La Grande-Anse","Nippes","Sud")){
  subgeo$normDisp[subgeo$Department==stringy]<-subgeo$maxdisp_norm[subgeo$Department==stringy]/
    normaliser$normy[normaliser$Department==stringy]
}

mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "terrain-background",zoom=10)
q<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude") +theme(plot.title = element_text(hjust = 0.5))

# for (dd in unique(subgeo$day)){
for (dd in ddd){
  tmp<-ODDy
  for(stringy in c("La Grande-Anse","Nippes","Sud")){
    inds<-tmp@data$Department==stringy & !is.na(tmp@data$Department)
    
    normzy<-subgeo$normDisp[subgeo$day==dd & subgeo$Department==stringy]
    if(length(normzy)==0) normzy<-min(0.5,subgeo$normDisp[subgeo$day==dd & subgeo$Department=="Sud"])
    tmp@data$RegionAggDisp[inds]<-
      tmp@data$RegionAggDisp[inds]*normzy
  }
    
  charry<-as.character(as.Date("2021-08-14")+dd)
  
  p<-plotODDyBG(tmp,zoomy=10,var="RegionAggDisp",bbox = c(-74.5,18.,-72.8,18.7),
                breakings=c(0,2000,5000,17000,100000),alpha=0.7,
                map="terrain-background",hazardplot = F,p = q)+
    ggtitle(paste0("Haiti Earthquake Population Displaced, ",charry))+
    scale_fill_brewer(palette = "Oranges")+
    geom_text(data = boundaries, mapping = aes(Longitude, Latitude, label = name),size = 4,check_overlap = TRUE)
  # for(j in 1:nrow(boundaries)){
  #   
  #   nc <- as(boundaries$geometry[j], "Spatial")
  #   p<-p+geom_polygon(boun)
  #   
  # }
  
  ggsave(paste0("FEW_Disp_",charry,".png"), plot=p,path = paste0(directory,'Plots/IIDIPUS_Results/EQ20210814HTI/'),width = 13,height = 6.)
  
}
    
former2<-geoinsights%>%filter(status%in%c("Displaced","Displaced Abroad"))%>%
  group_by(baseline_polygon_name,day)%>%summarise(
    maxdisp_norm=max(percent_of_total_population_normalized_to_overall),
    maxdisp=max(percent_of_total_population),.groups="drop")
names(former2)[1]<-"Department"

max(former2$maxdisp_norm)
former2$displaced<-former2$maxdisp_norm*15000

xdr<-former2%>%group_by(day)%>%summarise(displaced=sum(displaced,na.rm = T))


oneDgeo<-subgeo%>%group_by(day)%>%summarise(displaced=sum(PopDisp,na.rm = T)/186719*sum(ODDy$ExpDisp),
                                            .groups="drop")  
oneDgeo$Source<-"True Displacement                      "

oneDgeo$displaced[oneDgeo$Source=="True Displacement                      "&oneDgeo$day>68]<-oneDgeo$displaced[oneDgeo$Source=="True Displacement                      "&oneDgeo$day>68]*0.7

cols <- c("True Displacement                      " = "red", "Mobile Phone - XDR" = "blue",
          "Emergency Shelter" = "green3", "Buildings Destroyed" = "orange",
          "Displaced Destroyed Homes" = "purple", "News/Gov Reports" = "grey")

p<-oneDgeo%>%
  ggplot(aes(day,displaced,group=Source))+geom_point(aes(colour=Source))+
  geom_line(aes(colour=Source))+scale_y_log10(limits = c(1e4,3.2e5)) +
  scale_colour_manual(values = cols[1]) + xlab("Days Since Event Onset") + ylab("Displaced Population")
ggsave("Timeline1.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results/EQ20210814HTI/'),width = 8,height = 4.)

oneDgeo<-rbind(oneDgeo,cbind(xdr,data.frame(Source=rep("Mobile Phone - XDR",nrow(xdr)))))
unique(oneDgeo$Source)  
oneDgeo$displaced[oneDgeo$Source=="Mobile Phone - XDR"&oneDgeo$day>70]<-oneDgeo$displaced[oneDgeo$Source=="Mobile Phone - XDR"&oneDgeo$day>70]*0.8

p<-oneDgeo%>%
  ggplot(aes(day,displaced,group=Source))+geom_point(aes(colour=Source))+
  geom_line(aes(colour=Source))+scale_y_log10(limits = c(1e4,3.2e5)) +
  scale_colour_manual(values = cols[1:2]) + xlab("Days Since Event Onset") + ylab("Displaced Population")
ggsave("Timeline2.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results/EQ20210814HTI/'),width = 8,height = 4.)

oneDgeo%<>%rbind(data.frame(day=c(28,54,68,75,85),displaced=c(1e5,5e4,4e4,2.5e4,2.5e4),
                            Source=rep("Emergency Shelter",5)))

p<-oneDgeo%>%
  ggplot(aes(day,displaced,group=Source))+geom_point(aes(colour=Source))+
  geom_line(aes(colour=Source))+scale_y_log10(limits = c(1e4,3.2e5)) +
  scale_colour_manual(values = cols[1:3]) + xlab("Days Since Event Onset") + ylab("Displaced Population")
ggsave("Timeline3.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results/EQ20210814HTI/'),width = 8,height = 4.)

oneDgeo%<>%rbind(data.frame(day=c(45),displaced=c(1.2e5/4.3),
                            Source="Buildings Destroyed"))

p<-oneDgeo%>%
  ggplot(aes(day,displaced,group=Source))+geom_point(aes(colour=Source))+
  geom_line(aes(colour=Source))+scale_y_log10(limits = c(1e4,3.2e5)) +
  scale_colour_manual(values = cols[1:4])  + xlab("Days Since Event Onset") + ylab("Displaced Population")
ggsave("Timeline4.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results/EQ20210814HTI/'),width = 8,height = 4.)

oneDgeo%<>%rbind(data.frame(day=c(45),displaced=c(1.2e5),
                            Source="Displaced Destroyed Homes"))

p<-oneDgeo%>%
  ggplot(aes(day,displaced,group=Source))+geom_point(aes(colour=Source))+
  geom_line(aes(colour=Source))+scale_y_log10(limits = c(1e4,3.2e5)) +
  scale_colour_manual(values = cols[1:5])  + xlab("Days Since Event Onset") + ylab("Displaced Population")
ggsave("Timeline5.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results/EQ20210814HTI/'),width = 8,height = 4.)

oneDgeo%<>%rbind(data.frame(day=c(23,50,77),displaced=c(2e5,3.1e5,4e4),
                            Source=rep("News/Gov Reports",3)))

p<-oneDgeo%>%
  ggplot(aes(day,displaced,group=Source))+geom_point(aes(colour=Source))+
  geom_line(aes(colour=Source))+scale_y_log10(limits = c(1e4,3.2e5)) +
  scale_colour_manual(values = cols) + xlab("Days Since Event Onset") + ylab("Displaced Population")
ggsave("Timeline6.png", plot=p,path = paste0(directory,'Plots/IIDIPUS_Results/EQ20210814HTI/'),width = 8,height = 4.)
# 
# p<-oneDgeo%>%
#   ggplot(aes(day,displaced,group=Source))+geom_point(aes(colour=Source))+
#   geom_line(aes(colour=Source))+scale_y_log10(limits = c(1e4,3.2e5)) + scale_colour_manual(values = cols)
# p
# geoinsights%>%filter(day==17 & status%in%c("Displaced","Displaced Abroad") & baseline_polygon_name%in%c("La Grande-Anse","Nippes","Sud"))%>%
  # dplyr::select(n_people,baseline_polygon_name)

# Find a few buildings from OSM that are in Copernicus/UNOSAT and plot the building shape and damage profile
# DO THIS REALLY AT THE LAST MINUTE - YOU MIGHT NOT HAVE TIME!!!!!!!!!!!



