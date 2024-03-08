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


getCopernicusSource <- function(event_i){
 
  event_i <- 39
  haz <- 'EQ'
  folderin_BD <- "/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/BDobjects/"
  ufiles_BD <- list.files(path=folderin_BD,pattern=Model$haz,recursive = T,ignore.case = T)
  file_match_BD <- ufiles_BD[which(as.numeric(sub(".+_(.+)$", "\\1", ufiles_BD))==event_i)]
  
  BDy <- readRDS(paste0(folderin_BD, file_match_BD))
  
  BDy_date <- BDy@hazdates[1]
  cfiles<-list.files(path=paste0(dir,"COPERNICUS_Damage/"))
  
  evname<-c()
  for (i in 1:length(list.files(path=paste0(dir,"COPERNICUS_Damage/")))) {
    evname%<>%c(list.files(path=paste0(dir,"COPERNICUS_Damage/"))[[i]][1])
  }
  shp_dates <- as.Date(gsub("\\D", "", evname), format = "%Y%m%d")
  ev <- evname[which.min(abs(shp_dates-BDy_date))]
  
  cfiles<-list.files(path=paste0(dir,"COPERNICUS_Damage/", ev, '/'),pattern="crisis_information",recursive = T)
  cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/", ev, '/'),pattern="builtUpP",recursive = T))
  cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/", ev, '/'),pattern="built_up_p",recursive = T))
  cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/", ev, '/'),pattern="settlements_po",recursive = T))
  cfiles<-cfiles[grep(".shp",cfiles)] ; cfiles<-cfiles[grep(".shp.xml",cfiles,invert = T)]
  
  source_file <- list.files(path=paste0(dir,"COPERNICUS_Damage/", ev, '/'),pattern="source", recursive=T)
  source_file <-source_file[grep(".dbf",source_file)] 
  source_file <-source_file[grep(".dbf.xml",source_file, invert=T)]
  source_file_AOI <- str_extract(source_file, "(?<=AOI)\\d+")
  
  for (f in cfiles){
    tmp <- st_read(paste0(dir,"COPERNICUS_Damage/",  ev, '/', f),quiet = T)
    AOI <- str_extract(f, "(?<=AOI)\\d+")
    or_tbl <- table(tmp$or_src_id)
    dmg_tbl <- table(tmp$dmg_src_id)
    grd_tbl <- table(tmp$damage_gra)
    print(paste0('Area of Interest: ', AOI,'.  ', paste(paste(names(grd_tbl), grd_tbl, sep = ": "), collapse=', ')))
    print('Origin Data:')
    source_file_f <- source_file[which(source_file_AOI==AOI)]
    source_dbf <- read.dbf(paste0(dir,"COPERNICUS_Damage/",  ev, '/', source_file_f))
    for (i in 1:length(dmg_tbl)){
      cat(paste0(source_dbf$source_nam[which(source_dbf$src_id==names(or_tbl)[i])], ' : ', or_tbl[i], '. '))
      print('')
    }
    print('Damage Data:')
    for (i in 1:length(dmg_tbl)){
      cat(paste0(source_dbf$source_nam[which(source_dbf$src_id==names(dmg_tbl)[i])], ' : ', dmg_tbl[i], '. '))
      print('')
    }
  }
}



library(data.table)
library(geojsonio)
library(jsonlite)
BD_increase_coverage_bing <- function(BDy, ODDy){
  
  #event_i <- 67
  
  #folderin_BD <- "/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/BDobjects/"
  #ufiles_BD <- list.files(path=folderin_BD,pattern=Model$haz,recursive = T,ignore.case = T)
  #file_match_BD <- ufiles_BD[which(as.numeric(sub(".+_(.+)$", "\\1", ufiles_BD))==event_i)]
  
  #BDy <- readRDS(paste0(folderin_BD, file_match_BD))
  
  BDy_date <- BDy@hazdates
  
  if (as.Date("2013-09-24") %in% BDy_date){ #UNOSAT
    shapefiles <- '/home/manderso/Documents/GitHub/ODDRIN/UNOSAT_Damage/EQ20130924PAK_UNOSAT/Damage_Sites.shp'
    ev <- 'EQ20130924PAK'
    buildings_in_MAR_polys <- which(BDy@coords[,1] > 65 & BDy@coords[,1] < 65.65 & BDy@coords[,2] < 27.2)
    BDy@data <- BDy@data[buildings_in_MAR_polys,]
    BDy@coords <- BDy@coords[buildings_in_MAR_polys,]
    return(BDy)
  } else {
    print("WARNING: only point files and not polygons are currently read from COPERNICUS damage assessment")
    cfiles<-list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="areaOfInterest",recursive = T)
    cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="area_of_interest",recursive = T))
    cfiles<-cfiles[grep(haz,cfiles)]; cfiles<-cfiles[grep(".shp",cfiles)] ; cfiles<-cfiles[grep(".shp.xml",cfiles,invert = T)]
    tmp<-strsplit(cfiles,"/",fixed = T) ;
    evname<-c()
    for (i in 1:length(tmp)) {
      evname%<>%c(tmp[[i]][1])
      cfiles[i]<-paste0(dir,"COPERNICUS_Damage/",cfiles[i])
    
      
    }
    shp_dates <- as.Date(gsub("\\D", "", evname), format = "%Y%m%d")
    ev <- evname[which.min(abs(shp_dates-BDy_date))]

    subfiles<-unique(grep(ev,cfiles,value = T,fixed = T))
    shapefiles <- subfiles<-unique(grep(ev,cfiles,value = T,fixed = T))
    
    contains_dam_data <- rep(T, length(shapefiles))
    for (j in 1:length(shapefiles)){
      subdir <- paste0(sub("/[^/]*$", "", shapefiles[j]), '/')
      files_subdir <- list.files(path=subdir,pattern="areaOfInterest",recursive = T)
      subdir_dam_files<-list.files(path=subdir,pattern="crisis_information",recursive = T)
      subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="builtUpP",recursive = T))
      subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="built_up_p",recursive = T))
      subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="settlements_po",recursive = T))
      subdir_dam_files<-subdir_dam_files[grep(".shp",subdir_dam_files)] ; 
      subdir_dam_files<-subdir_dam_files[grep(".shp.xml",subdir_dam_files,invert = T)]
      print(subdir_dam_files)
      if (length(subdir_dam_files) == 0){contains_dam_data[j] <- F}
    }
    
    #filter to only shapefiles that have relevant damage information:
    
    
  }
  
  shapefiles <- shapefiles[contains_dam_data]
  
  shape_data <- list()
  treat_as_MAR <- rep(NA, length(shape_data))
  for (i in 1:length(shapefiles)){
    shape_data[[i]] <- st_read(shapefiles[i])
    shape_data[[i]] <- st_transform(shape_data[[i]], crs = 4326)
    points_in_shape <- st_within(st_as_sf(as.data.frame(BDy@coords), coords = c("Longitude", "Latitude"), crs=st_crs(shape_data[[i]])), shape_data[[i]], sparse = FALSE)
    treat_as_MAR[i] <- ifelse(length(points_in_shape) < 10 | length(BDy$grading[points_in_shape]=='notaffected')/sum(points_in_shape) < 0.2, F, T)
  }
  
  shapefiles <- shapefiles[which(treat_as_MAR)]
  shape_data <- shape_data[which(treat_as_MAR)]
  
  combined_bbox <- st_bbox(shape_data[[1]])
  
  # Iterate through each element in shape_data starting from the second element
  for (i in seq_along(shape_data)[-1]) {
    
    # Update the combined bounding box to include the current object's bounding box
    combined_bbox["xmin"] <- min(combined_bbox["xmin"], st_bbox(shape_data[[i]])["xmin"])
    combined_bbox["ymin"] <- min(combined_bbox["ymin"], st_bbox(shape_data[[i]])["ymin"])
    combined_bbox["xmax"] <- max(combined_bbox["xmax"], st_bbox(shape_data[[i]])["xmax"])
    combined_bbox["ymax"] <- max(combined_bbox["ymax"], st_bbox(shape_data[[i]])["ymax"])
  }
  
  #BDy_coords_sf <- st_as_sf(data.table(BDy@coords), coords = c("Longitude", "Longitude"), crs = 4326) %>% st_transform("+proj=laea +lat_0=30 +lon_0=-95")
  
  #st_as_sf(data.table(BDy@coords), coords = c("Longitude", "Latitude"))
  
  col_values <- list('OpenBuildings'='white',
                     'notaffected'='blue',
                     'possible'='yellow', 
                     'moderate'='orange',
                     'Damaged'='red',
                     'severe' = 'hotpink',
                     'destroyed' = 'violetred')
  
  
  p <- ggplot() 
  for (i in 1:length(shape_data)){
    p <- p + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red")
  }
  p <- p + #geom_sf(data=st_as_sf(data.table(BDy@coords), coords = c("Longitude", "Latitude"), crs=4326), size=0.5) +#,  #color=BDy$grading) + 
    geom_point(data=data.table(BDy@coords), aes(x=Longitude, y=Latitude, color=BDy$grading),size=0.5) + theme_minimal() +
    scale_color_manual(values = col_values) + 
    xlim(combined_bbox[1], combined_bbox[3]) + ylim(combined_bbox[2], combined_bbox[4]) 
    #xlim(85.2, 85.4) + ylim(27.62, 27.75)
  p
  
  p <- list()
  for (i in 1:length(shape_data)){
    p[[i]] <- ggplot() + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red") + xlim(st_bbox(shape_data[[i]])[c(1,3)])+ylim(st_bbox(shape_data[[i]])[c(2,4)]) +
      geom_point(data=BDy@data[rev(order(BDy@data$grading)),], aes(x=BDy$Longitude, y=BDy$Latitude, color=BDy$grading), size=1) + 
      scale_color_manual(values = col_values) +  theme(legend.position = "none",  # Remove legend
                                                       axis.title.x = element_blank(),  # Remove x-axis title
                                                       axis.title.y = element_blank(),  # Remove y-axis title
                                                       axis.text.x = element_blank(),   # Remove x-axis labels
                                                       axis.text.y = element_blank()) 
  }
  do.call(grid.arrange,p)
  
  
  iso3_unique <- unique(BDy$ISO3C)[!is.na(unique(BDy$ISO3C))]
  lon_min <- combined_bbox[1]; lon_max <- combined_bbox[3]
  lat_min <- combined_bbox[2]; lat_max <- combined_bbox[4]
  
  if (all(iso3_unique %in% c('IDN', 'PHL'))){
    national_coverage = T
    open_buildings = T
  } else {
    open_buildings = F
  }
  
  if (open_buildings){ 
    
    #####################################################################################
    ################################# OPEN BUILDINGS ####################################
    #####################################################################################
    
    if (national_coverage){
      open_buildings_file <- paste0(dir, 'Demography_Data/Buildings/open_buildings_v2_points_your_own_wkt_polygon_',iso3_unique[1],'.csv')
      if(!file.exists(open_buildings_file)){
        stop('Download Open Buildings data from: https://colab.research.google.com/github/google-research/google-research/blob/master/building_detection/open_buildings_download_region_polygons.ipynb')
      }
      #When downloading, insert the polygon bounding the country of interest into your_own_wkt_polygon field
      # alternatively, can select the region_border_source and country, although this doesn't seem to work well when working with small islands in Philippines (and potentially other small regions)
      # select 'points' in the data_type field
    } else {
      open_buildings_files <- paste0(dir, 'Demography_Data/Buildings/open_buildings_v2_points_your_own_wkt_polygon_i',i,'.csv')
      if(!any(file.exists(open_buildings_files))){
        stop('Download Open Buildings data from: https://colab.research.google.com/github/google-research/google-research/blob/master/building_detection/open_buildings_download_region_polygons.ipynb')
      }
      open_buildings_file <- open_buildings_files[which(file.exists(open_buildings_files))]
    }
    
    # doesn't work once files exceed a certain size:
    # building_locs <- read.csv.sql(open_buildings_file,
    #                               paste0("select longitude, latitude from file where latitude > ", lat_min, ' AND longitude > ', lon_min,
    #                                      ' AND latitude < ',lat_max, ' AND longitude < ', lon_max))
    
    building_locs <- data.frame(latitude=double(), longitude=double())
    i <- 1
    nrows_file <- as.integer(strsplit(system(paste0('wc -l ', open_buildings_file), intern=T), ' ')[[1]][1])
    nrow_tmp <- 50000000
    nchunks <- ceiling(nrows_file/nrow_tmp)
    for (j in 1:nchunks){
      tmp <- fread(open_buildings_file,skip=i, nrows=nrow_tmp, select=c(1,2), col.names=c('latitude', 'longitude'))
      building_locs %<>% rbind(tmp %>% filter(latitude > lat_min, longitude > lon_min, latitude < lat_max, longitude < lon_max))
      i <- i + nrow_tmp
    }
  } else {
    
    #####################################################################################
    ################################# BING BUILDINGS ####################################
    #####################################################################################
    
    bbox <- combined_bbox
    zoom <- 9
    
    tiles <- getMercantileTiles(bbox[1], bbox[2], bbox[3], bbox[4], zoom)
    quad_keys <- list()
    for (i in 1:NROW(tiles)){
      quad_keys[[i]] <- getQuadKey(tiles[i,])
    }
    
    cat(sprintf("The input area spans %d tiles: %s\n", length(quad_keys), paste(quad_keys, collapse = ", ")))
    
    df <- read.csv("https://minedbuildings.blob.core.windows.net/global-buildings/dataset-links.csv")
    
    # Download the GeoJSON files for each tile that intersects the input geometry
    build_coords <- array(dim=c(0,2))
    
    missing_quadkeys <- unlist(quad_keys)[!unlist(quad_keys) %in% df$QuadKey]
    # if (length(missing_quadkeys)> (0.2 * length(quad_keys))){
    #   file_conn <- file(file_write, open = "a")
    #   writeLines(paste("Event:", event_id, ", Missing", length(missing_quadkeys)/length(quad_keys)*100, "percent of quad keys, not adding building data."), file_conn)
    #   close(file_conn) 
    #   return(ODD)
    # } else if (length(missing_quadkeys)> 0){
    #   file_conn <- file(file_write, open = "a")
    #   writeLines(paste("Event:", event_id, ", Missing", length(missing_quadkeys)/length(quad_keys)*100, "percent of quad keys, but still adding building data."), file_conn)
    #   close(file_conn) 
    # }
    if (length(missing_quadkeys)> 0){
      file_conn <- file(paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/BD_creation_notes'), open = "a")
      writeLines('Bing/Open Buildings building counts missing', file_conn)
      close(file_conn) 
      return(NULL)
      #stop(paste("Event:", event_i, ", Missing", length(missing_quadkeys)/length(quad_keys)*100, "percent of quad keys, not adding building data."))
    }
    
    missing_quadkeys_flag <- F
    
    for (quad_key in quad_keys) {
      rows <- df[df$QuadKey == quad_key, ]
      if (nrow(rows) >= 1) {
        for (url in rows$Url){
          tmp <- tempfile()
          download.file(url, destfile =tmp,quiet = FALSE, mode = "wb")
          out <- lapply(readLines(tmp), fromJSON)
          
          #convert building polygons to points by just taking the first in the coordinate 
          build_coords %<>% rbind(t(sapply(out, function(build) return(build$geometry$coordinates[1,1,]))))
        }
      } else {
        missing_quadkeys_flag <- T
        # file_conn <- file('IIDIPUS_Input_NonFinal/IIDIPUS_Input_June24_AggFactor5/BingGlobalNotes', open = "a")
        # writeLines(paste("               Event Date:", ODD@hazdates[1], ". No building data for quadkey", quad_key), file_conn)
        # close(file_conn) 
        # print(paste("QuadKey not found in dataset:", quad_key))
        #stop(paste("QuadKey not found in dataset:", quad_key))
      }
    }
    inside_bbox <- which(build_coords[,1] > bbox[1] & build_coords[,1] < bbox[3] & build_coords[,2] > bbox[2] & build_coords[,2] < bbox[4])
    building_locs <- build_coords[inside_bbox,]
    
    if (!missing_quadkeys_flag){
      print(paste("Complete Building Count from Global Bing Building Footprints"))
    }
    
    colnames(building_locs) <- c('longitude', 'latitude')
    building_locs <- data.table(building_locs)
  }
  
  building_locs_sf <- st_as_sf(building_locs, coords = c("longitude", "latitude"))
  
  # Initialize the poly_int column
  building_locs$poly_int <- FALSE
  
  # Loop through each shape in shape_data
  for (j in 1:length(shape_data)) {
    # Check if any points in building_locs_sf are within the current polygon
    shape_data[[j]] <- st_set_crs(shape_data[[j]], st_crs(building_locs_sf))
    within <- st_within(building_locs_sf, shape_data[[j]], sparse = FALSE)
    # Update poly_int for points within the polygon
    building_locs$poly_int[within] <- TRUE
  }
  building_locs <- building_locs %>% filter(poly_int==T)
  
  BDy_coords <- BDy@coords
  building_locs_miss <- building_locs[,1:2]
  building_locs_miss <- building_locs_miss[order(building_locs_miss$longitude),]
  building_locs_miss$rem <- T
  building_locs_miss$id <- 1:NROW(building_locs_miss)
  
  find_closest_loc <- function(lon1, lat1, building_locs_miss){
    filt <- building_locs_miss[latitude>(lat1-0.00025) & latitude < (lat1+0.00025) & longitude > (lon1-0.00025) & longitude < (lon1+0.00025) & rem==T,]
    #filt <- building_locs_miss %>% filter((latitude > (lat1-0.0005)) & (latitude < (lat1+0.0005)) &
    #                              (longitude > (lon1-0.0005)) & (longitude < (lon1+0.0005)))
    building_locs_miss$rem[filt$id[which.min(apply(filt[,1:2], 1, function(x) sum((x[1]-lat1)^2+(x[2]-lon1)^2)))]] = F
    return(building_locs_miss)
  }
  for (i in 1:NROW(BDy_coords)){
    building_locs_miss <- find_closest_loc(BDy_coords[i,1], BDy_coords[i,2], building_locs_miss)
    #if(length(i_match)==0) next
    #building_locs_miss <- building_locs_miss[-which(building_locs_miss$id==i_match), ]
  }
  building_locs_miss <- building_locs_miss[which(building_locs_miss$rem), ]
  
  
  # polygons_MAR <- 1:11 #c(1,2,4,5,7,8,9,12,13)
  # building_locs_miss_notMAR <- building_locs_miss
  # for (j in polygons_MAR){
  #   bbox_j <- st_bbox(shape_data[[j]])
  #   building_locs_miss_notMAR <- building_locs_miss_notMAR %>% filter(!(building_locs_miss_notMAR$longitude > bbox_j[1] & building_locs_miss_notMAR$longitude < bbox_j[3] &
  #                                                                  building_locs_miss_notMAR$latitude > bbox_j[2] & building_locs_miss_notMAR$latitude < bbox_j[4]))
  # }
  # building_locs_miss <- building_locs_miss_notMAR
  
  miniDam_old <- data.frame(grading=BDy$grading, Longitude=BDy$Longitude, Latitude=BDy$Latitude, event=ev, Confidence=BDy$Confidence, hazard='EQ', sdate=BDy@hazdates[1], iso3=BDy@hazdates[1])
  
  building_locs_old <- st_as_sf(miniDam_old, coords = c("Longitude", "Latitude"))
  building_locs_old$poly_int <- FALSE
  # Loop through each shape in shape_data
  for (j in 1:length(shape_data)) {
    # Check if any points in building_locs_sf are within the current polygon
    shape_data[[j]] <- st_set_crs(shape_data[[j]], st_crs(building_locs_old))
    within <- st_within(building_locs_old, shape_data[[j]], sparse = FALSE)
    # Update poly_int for points within the polygon
    building_locs_old$poly_int[within] <- TRUE
  }
  miniDam_old <- miniDam_old[which(building_locs_old$poly_int==T),]
  
  p <- ggplot()
  for (i in 1:length(shape_data)){
    p <- p + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red")
  }
  p <- p + geom_point(data=data.frame(building_locs_miss), aes(x=longitude, y=latitude, col='OpenBuildings'), size=0.1) + 
    geom_point(data=miniDam_old, aes(x=miniDam_old$Longitude, y=miniDam_old$Latitude, color=miniDam_old$grading), size=0.5) + 
    scale_color_manual(values = col_values)
  p
  
  i <- 1
  p <- ggplot() + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red") + xlim(st_bbox(shape_data[[i]])[c(1,3)])+ylim(st_bbox(shape_data[[i]])[c(2,4)]) +
    geom_point(data=miniDam_old, aes(x=miniDam_old$Longitude, y=miniDam_old$Latitude, color=miniDam_old$grading), size=1) +
    scale_color_manual(values = col_values)
  
  p_after <- p + geom_point(data=data.frame(building_locs_miss), aes(x=longitude, y=latitude, col='OpenBuildings'), size=0.1) + 
    geom_point(data=BDy@data[rev(order(BDy@data$grading)),], aes(x=BDy$Longitude, y=BDy$Latitude, color=BDy$grading), size=1)
  
  grid.arrange(p, p_after, nrow=2)
  if (NROW(building_locs_miss) > 0){
    miniDam <- rbind(miniDam_old, data.frame(grading='notaffected', Longitude=building_locs_miss$longitude, Latitude=building_locs_miss$latitude, 
                                             event=ev, Confidence=NA, hazard='EQ', sdate=BDy@hazdates[1], iso3=BDy@hazdates[1]))
  } else {
    miniDam <- miniDam_old
  }
  
  # folderin_ODD <- "/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/"
  # ufiles_ODD <- list.files(path=folderin_ODD,pattern=Model$haz,recursive = T,ignore.case = T)
  # file_match_ODD <- ufiles_ODD[which(as.numeric(sub(".+_(.+)$", "\\1", ufiles_ODD))==event_i)]
  # 
  # ODDy <- readRDS(paste0(folderin_ODD, file_match_ODD))
  
  BD_new <- new("BD",Damage=miniDam,ODD=ODDy)
  
  return(BD_new)
  
  # poly_3_i <- st_within(st_as_sf(as.data.frame(BD_new@coords), coords = c("Longitude", "Latitude")), shape_data[[3]], sparse = FALSE)
  # BD_new@data <- BD_new@data[poly_3_i,]
  # BD_new@coords <- BD_new@coords[poly_3_i,]
  
  plot_bd_df_int <- data.frame(intensity=numeric(), bd_new=logical(), obs_prop_bd=numeric(), tot_obs=integer())
  for(int in seq(4.25,9.25, 0.5)){
    plot_bd_df_int %<>% add_row(intensity=int, bd_new=T, 
                                obs_prop_bd=1-length(which(BD_new$grading[which(BD_new$hazMean1>int &BD_new$hazMean1<(int+0.2))]=='notaffected'))/length(which(BD_new$hazMean1>int & BD_new$hazMean1<(int+0.2))),
                                tot_obs=length(which(BD_new$hazMean1>int & BD_new$hazMean1<(int+0.2))))
    plot_bd_df_int %<>% add_row(intensity=int, bd_new=F, 
                                obs_prop_bd=1-length(which(BDy$grading[which(BDy$hazMean1>int &BDy$hazMean1<(int+0.2))]=='notaffected'))/length(which(BDy$hazMean1>int & BDy$hazMean1<(int+0.2))),
                                tot_obs=length(which(BDy$hazMean1>int & BDy$hazMean1<(int+0.2))))
  }
  
  plot_bd_df_int$std_dev <- sqrt(plot_bd_df_int$obs_prop_bd * (1 - plot_bd_df_int$obs_prop_bd) / plot_bd_df_int$tot_obs)
  # Create a ggplot with shaded region
  ggplot(plot_bd_df_int %>% filter(!is.na(plot_bd_df_int$obs_prop_bd)), aes(x = intensity, y = obs_prop_bd, group = bd_new, color = bd_new)) +
    geom_line() +
    geom_ribbon(aes(ymin = obs_prop_bd - 3 * std_dev, ymax = obs_prop_bd + 3 * std_dev, fill = bd_new), alpha=0.1,
                color = NA) +
    labs(title = "Line with Shaded Region for 3 Standard Deviations",
         x = "Intensity",
         y = "Observed Proportion")
  
  #saveRDS(plot_bd_df_int, paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/BD_propDam_vs_I/',ev))
}

# plot_bd_df_int_all <- rbind(plot_df_int_all, plot_bd_df_int )
# 
# folderin<-"/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/"
# ufiles_wFolder <- list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
# 
# 
# 
# for (file in ufiles_wFolder){
#   ODDy <- readRDS(paste0(folderin, file))
#   print(paste(file, 'contains nBuildings:', !is.null(ODDy$nBuildings)))
# }

plot_sat_damage_shps <- function(BDy, shape_data, ODDy=NULL, single_plot=T, plot_marks=rep(T, length(shape_data))){
  col_values <- list('OpenBuildings'='white',
                     'notaffected'='green',
                     'possible'='yellow', 
                     'moderate'='orange',
                     'Damaged'='red',
                     'severe' = 'purple',
                     'destroyed' = 'blue')
  
  if (single_plot){ # plot the whole region
    
    combined_bbox <- st_bbox(shape_data[[1]])
    
    # Iterate through each element in shape_data starting from the second element
    for (i in seq_along(shape_data)[-1]) {
      
      # Update the combined bounding box to include the current object's bounding box
      combined_bbox["xmin"] <- min(combined_bbox["xmin"], st_bbox(shape_data[[i]])["xmin"])
      combined_bbox["ymin"] <- min(combined_bbox["ymin"], st_bbox(shape_data[[i]])["ymin"])
      combined_bbox["xmax"] <- max(combined_bbox["xmax"], st_bbox(shape_data[[i]])["xmax"])
      combined_bbox["ymax"] <- max(combined_bbox["ymax"], st_bbox(shape_data[[i]])["ymax"])
    }
    
    p <- ggplot() 
    for (i in 1:length(shape_data)){
      p <- p + geom_sf(data = st_as_sf(shape_data[[i]]), color=ifelse(plot_marks[i],'black', "red") , fill = ifelse(plot_marks[i],'white', "red"), alpha=0.05, lwd=1)
    }
    p <- p + #geom_sf(data=st_as_sf(data.table(BDy@coords), coords = c("Longitude", "Latitude"), crs=4326), size=0.5) +#,  #color=BDy$grading) + 
      geom_point(data=data.table(BDy@coords), aes(x=Longitude, y=Latitude, color=BDy$grading),size=0.5) + theme_minimal() +
      scale_color_manual(values = col_values) + 
      xlim(combined_bbox[1], combined_bbox[3]) + ylim(combined_bbox[2], combined_bbox[4]) 
    #xlim(85.2, 85.4) + ylim(27.62, 27.75)
    if(!is.null(ODDy)){
      hazard<-rep(NA_real_,length(ODDy@data$hazMean1))
      for (variable in names(ODDy)[grepl("Mean",names(ODDy))]){
        tmp<-ODDy[variable]
        tmp$hazard<-hazard
        hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
      }
      ODDy@data$hazard<-hazard
      brks<-seq(9,ceiling(2*max(hazard,na.rm = T)),by=1)/2
      
      p<-p+geom_contour(data = as.data.frame(ODDy),
                        mapping = aes(Longitude,Latitude,z=ifelse(is.na(hazard), 0, hazard), alpha=after_stat(level)),
                        breaks = brks, lwd=0.7, col='black') +  labs(alpha = "Hazard Intensity")
    }
    return(p)
  } else {
    hazard<-rep(NA_real_,length(BDy@data$hazMean1))
    for (variable in names(BDy)[grepl("Mean",names(BDy))]){
      tmp<-BDy[variable]
      tmp$hazard<-hazard
      hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
    }
    
    p <- list()
    for (i in 1:length(shape_data)){
      BDy_in <- which(st_within(st_as_sf(as.data.frame(BDy@coords), coords = c("Longitude", "Latitude"), crs=st_crs(shape_data[[i]])), shape_data[[i]], sparse = FALSE))
      haz_range <- paste0('Haz Range: ', paste(round(range(hazard[BDy_in], na.rm=T), 1), collapse=' - '))
      p[[i]] <- ggplot() + geom_sf(data = st_as_sf(shape_data[[i]]), color=ifelse(plot_marks[i],'black', "red") , fill = ifelse(plot_marks[i],'white', "red"), alpha=0.05, lwd=1) + 
        xlim(st_bbox(shape_data[[i]])[c(1,3)])+ylim(st_bbox(shape_data[[i]])[c(2,4)]) +
        geom_point(data=BDy@data[rev(order(BDy@data$grading)),], aes(x=BDy$Longitude, y=BDy$Latitude, color=BDy$grading), size=1) + 
        scale_color_manual(values = col_values) +  theme(legend.position = "none",  # Remove legend
                                                         axis.title.y = element_blank(),
                                                         #axis.title.x = element_blank(),  # Remove x-axis title
                                                         axis.text.x = element_blank(),   # Remove x-axis labels
                                                         axis.text.y = element_blank()) + xlab(haz_range)
    }
    do.call("grid.arrange", c(p, ncol=4))
    return(do.call(grid.arrange,p))
  }
}
# 
# select_MAR_polys <- function(BDy, ODDy){
#   
#   #event_i <- 67
#   
#   #folderin_BD <- "/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/BDobjects/"
#   #ufiles_BD <- list.files(path=folderin_BD,pattern=Model$haz,recursive = T,ignore.case = T)
#   #file_match_BD <- ufiles_BD[which(as.numeric(sub(".+_(.+)$", "\\1", ufiles_BD))==event_i)]
#   
#   #BDy <- readRDS(paste0(folderin_BD, file_match_BD))
#   
#   BDy_date <- BDy@hazdates
#   
#   if (as.Date("2013-09-24") %in% BDy_date){ #UNOSAT
#     shapefiles <- '/home/manderso/Documents/GitHub/ODDRIN/UNOSAT_Damage/EQ20130924PAK_UNOSAT/Damage_Sites.shp'
#     ev <- 'EQ20130924PAK'
#     buildings_in_MAR_polys <- which(BDy@coords[,1] > 65 & BDy@coords[,1] < 65.65 & BDy@coords[,2] < 27.2)
#     BDy@data <- BDy@data[buildings_in_MAR_polys,]
#     BDy@coords <- BDy@coords[buildings_in_MAR_polys,]
#     return(BDy)
#   } else {
#     print("WARNING: only point files and not polygons are currently read from COPERNICUS damage assessment")
#     cfiles<-list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="areaOfInterest",recursive = T)
#     cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="area_of_interest",recursive = T))
#     cfiles<-cfiles[grep(haz,cfiles)]; cfiles<-cfiles[grep(".shp",cfiles)] ; cfiles<-cfiles[grep(".shp.xml",cfiles,invert = T)]
#     tmp<-strsplit(cfiles,"/",fixed = T) ;
#     evname<-c()
#     for (i in 1:length(tmp)) {
#       evname%<>%c(tmp[[i]][1])
#       cfiles[i]<-paste0(dir,"COPERNICUS_Damage/",cfiles[i])
#     }
#     shp_dates <- as.Date(gsub("\\D", "", evname), format = "%Y%m%d")
#     ev <- evname[which.min(abs(shp_dates-BDy_date))]
#     
#     subfiles<-unique(grep(ev,cfiles,value = T,fixed = T))
#     shapefiles <- subfiles<-unique(grep(ev,cfiles,value = T,fixed = T))
#     
#     contains_dam_data <- rep(T, length(shapefiles))
#     for (j in 1:length(shapefiles)){
#       subdir <- paste0(sub("/[^/]*$", "", shapefiles[j]), '/')
#       files_subdir <- list.files(path=subdir,pattern="areaOfInterest",recursive = T)
#       subdir_dam_files<-list.files(path=subdir,pattern="crisis_information",recursive = T)
#       subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="builtUpP",recursive = T))
#       subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="built_up_p",recursive = T))
#       subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="settlements_po",recursive = T))
#       subdir_dam_files<-subdir_dam_files[grep(".shp",subdir_dam_files)] ; 
#       subdir_dam_files<-subdir_dam_files[grep(".shp.xml",subdir_dam_files,invert = T)]
#       print(subdir_dam_files)
#       if (length(subdir_dam_files) == 0){contains_dam_data[j] <- F}
#     }
#     
#     #filter to only shapefiles that have relevant damage information:
#   }
#   
#   shapefiles <- shapefiles[contains_dam_data]
#   
#   shape_data <- list()
#   treat_as_MAR <- rep(NA, length(shape_data))
#   for (i in 1:length(shapefiles)){
#     shape_data[[i]] <- st_read(shapefiles[i])
#     shape_data[[i]] <- st_transform(shape_data[[i]], crs = 4326)
#     points_in_shape <- st_within(st_as_sf(as.data.frame(BDy@coords), coords = c("Longitude", "Latitude"), crs=st_crs(shape_data[[i]])), shape_data[[i]], sparse = FALSE)
#     treat_as_MAR[i] <- ifelse(sum(points_in_shape) < 10 | sum(BDy$grading[points_in_shape]=='notaffected')/sum(points_in_shape) < 0.2, F, T)
#   }
#   
#   file_conn <- file(paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/BD_creation_notes'), open = "a")
#   writeLines(paste0(c("Regions treated as MAR:", which(treat_as_MAR)), collapse=' '), file_conn)
#   writeLines(paste0(c("Regions removed:", which(!treat_as_MAR)), collapse=' '), file_conn)
#   close(file_conn) 
#   
#   if(sum(treat_as_MAR)==0){
#     return(NULL)
#   }
#   
#   print(plot_sat_damage_shps(BDy, shape_data, ODDy=ODDy, plot_marks=treat_as_MAR))
#   
#   shapefiles <- shapefiles[which(treat_as_MAR)]
#   shape_data <- shape_data[which(treat_as_MAR)]
#   
#   
#   miniDam_old <- data.frame(grading=BDy$grading, Longitude=BDy$Longitude, Latitude=BDy$Latitude, event=ev, Confidence=BDy$Confidence, hazard='EQ', sdate=BDy@hazdates[1], iso3=BDy@hazdates[1])
#   
#   building_locs_old <- st_as_sf(miniDam_old, coords = c("Longitude", "Latitude"))
#   building_locs_old$poly_int <- FALSE
#   # Loop through each shape in shape_data
#   for (j in 1:length(shape_data)){
#     # Check if any points in building_locs_sf are within the current polygon
#     shape_data[[j]] <- st_set_crs(shape_data[[j]], st_crs(building_locs_old))
#     within <- st_within(building_locs_old, shape_data[[j]], sparse = FALSE)
#     # Update poly_int for points within the polygon
#     building_locs_old$poly_int[within] <- TRUE
#     
#   }
#   miniDam_old <- miniDam_old[which(building_locs_old$poly_int==T),]
#   
#   # p <- ggplot()
#   # for (i in 1:length(shape_data)){
#   #   p <- p + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red")
#   # }
#   # p <- p + geom_point(data=data.frame(building_locs_miss), aes(x=longitude, y=latitude, col='OpenBuildings'), size=0.1) + 
#   #   geom_point(data=miniDam_old, aes(x=miniDam_old$Longitude, y=miniDam_old$Latitude, color=miniDam_old$grading), size=0.5) + 
#   #   scale_color_manual(values = col_values)
#   # p
#   # 
#   # i <- 1
#   # p <- ggplot() + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red") + xlim(st_bbox(shape_data[[i]])[c(1,3)])+ylim(st_bbox(shape_data[[i]])[c(2,4)]) +
#   #   geom_point(data=miniDam_old, aes(x=miniDam_old$Longitude, y=miniDam_old$Latitude, color=miniDam_old$grading), size=1) +
#   #   scale_color_manual(values = col_values)
#   # 
#   # p_after <- p + geom_point(data=data.frame(building_locs_miss), aes(x=longitude, y=latitude, col='OpenBuildings'), size=0.1) + 
#   #   geom_point(data=BDy@data[rev(order(BDy@data$grading)),], aes(x=BDy$Longitude, y=BDy$Latitude, color=BDy$grading), size=1)
#   # 
#   # grid.arrange(p, p_after, nrow=2)
#   
#   miniDam <- rbind(miniDam_old, data.frame(grading='notaffected', Longitude=building_locs_miss$longitude, Latitude=building_locs_miss$latitude, 
#                                            event=ev, Confidence=NA, hazard='EQ', sdate=BDy@hazdates[1], iso3=BDy@hazdates[1]))
#   
#   # folderin_ODD <- "/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/"
#   # ufiles_ODD <- list.files(path=folderin_ODD,pattern=Model$haz,recursive = T,ignore.case = T)
#   # file_match_ODD <- ufiles_ODD[which(as.numeric(sub(".+_(.+)$", "\\1", ufiles_ODD))==event_i)]
#   # 
#   # ODDy <- readRDS(paste0(folderin_ODD, file_match_ODD))
#   
#   BD_new <- new("BD",Damage=miniDam,ODD=ODDy)
# 
#   return(BD_new)
# }


#UNFINISHED:
select_MAR_polys <- function(BDy, ODDy, bySource=T){

  #event_i <- 67
  
  #folderin_BD <- "/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/BDobjects/"
  #ufiles_BD <- list.files(path=folderin_BD,pattern=Model$haz,recursive = T,ignore.case = T)
  #file_match_BD <- ufiles_BD[which(as.numeric(sub(".+_(.+)$", "\\1", ufiles_BD))==event_i)]
  
  #BDy <- readRDS(paste0(folderin_BD, file_match_BD))
  
  BDy_date <- BDy@hazdates[1]
  
  if (as.Date("2013-09-24") %in% BDy_date){ #UNOSAT
    shapefiles <- '/home/manderso/Documents/GitHub/ODDRIN/UNOSAT_Damage/EQ20130924PAK_UNOSAT/Damage_Sites.shp'
    ev <- 'EQ20130924PAK' 
    buildings_in_MAR_polys <- which(BDy@coords[,1] > 65 & BDy@coords[,1] < 65.65 & BDy@coords[,2] < 27.2)
    BDy@data <- BDy@data[buildings_in_MAR_polys,]
    BDy@coords <- BDy@coords[buildings_in_MAR_polys,]
    return(BDy)
  } else if (as.Date("2019-09-24") %in% BDy_date){
    #all buildings in unosat for event are possibly damaged
    return(NULL)
  } else {
    print("WARNING: only point files and not polygons are currently read from COPERNICUS damage assessment")
    cfiles<-list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="areaOfInterest",recursive = T)
    cfiles<-c(cfiles,list.files(path=paste0(dir,"COPERNICUS_Damage/"),pattern="area_of_interest",recursive = T))
    cfiles<-cfiles[grep(haz,cfiles)]; cfiles<-cfiles[grep(".shp",cfiles)] ; cfiles<-cfiles[grep(".shp.xml",cfiles,invert = T)]
    tmp<-strsplit(cfiles,"/",fixed = T) ;
    evname<-c()
    for (i in 1:length(tmp)) {
      evname%<>%c(tmp[[i]][1])
      cfiles[i]<-paste0(dir,"COPERNICUS_Damage/",cfiles[i])
    }
    shp_dates <- as.Date(gsub("\\D", "", evname), format = "%Y%m%d")
    ev <- evname[which.min(abs(shp_dates-BDy_date))]
    
    subfiles<-unique(grep(ev,cfiles,value = T,fixed = T))
    shapefiles <- subfiles<-unique(grep(ev,cfiles,value = T,fixed = T))
    
    contains_dam_data <- rep(T, length(shapefiles))
    for (j in 1:length(shapefiles)){
      subdir <- paste0(sub("/[^/]*$", "", shapefiles[j]), '/')
      files_subdir <- list.files(path=subdir,pattern="areaOfInterest",recursive = T)
      subdir_dam_files<-list.files(path=subdir,pattern="crisis_information",recursive = T)
      subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="builtUpP",recursive = T))
      subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="built_up_p",recursive = T))
      subdir_dam_files<-c(subdir_dam_files,list.files(path=subdir,pattern="settlements_po",recursive = T))
      subdir_dam_files<-subdir_dam_files[grep(".shp",subdir_dam_files)] ; 
      subdir_dam_files<-subdir_dam_files[grep(".shp.xml",subdir_dam_files,invert = T)]
      print(subdir_dam_files)
      if (length(subdir_dam_files) == 0){contains_dam_data[j] <- F}
    }
    if (bySource){
      Damage<-data.frame()
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
      
      subfiles<-unique(grep(ev,cfiles,value = T,fixed = T))
      # I know that you can load all in at the same time, but let's save some RAM!!!
      for (f in unique(subfiles)){
        if(grepl(".shp.xml",f)) next
        tmp<-st_read(f,quiet = T)
        if (is.null(tmp$or_src_id)){
          selected_cols <- c('geometry', 'source_nam', 'src_date', 'src_info')
        } else {
          selected_cols <- c('geometry', 'or_src_id', 'dmg_src_id', 'det_method', 'obj_type')
        }
        if((is.null(tmp$grading) & is.null(tmp$damage_gra))|is.null(tmp$geometry)) next
        if(!is.null(tmp$grading)) {
          tmp%<>%dplyr::select(grading,!!!syms(selected_cols))
        } else if (!is.null(tmp$damage_gra)) {
          tmp%<>%dplyr::select(damage_gra,!!!syms(selected_cols))
          colnames(tmp)[which(colnames(tmp)=='damage_gra')] <- 'grading'
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
        tmp$AOI <- as.numeric(str_extract(f, "(?<=AOI)\\d+"))
        
        #get dbs source file
        src_files <- list.files(path=paste0(sub("/[^/]*$", "", f), '/'),pattern="source",recursive = T)
        src_files <- src_files[grep(".dbf",src_files)]
        if (length(src_files) > 0){
          source_data <- read.dbf(paste0(paste0(sub("/[^/]*$", "", f), '/'), src_files))
          tmp_with_or_src <- merge(tmp, source_data, by.x='or_src_id', by.y='src_id', all.x=T, sort=F)
          names(tmp_with_or_src)[(NCOL(tmp)+1):NCOL(tmp_with_or_src)] <- paste0('or_', names(tmp_with_or_src)[(NCOL(tmp)+1):NCOL(tmp_with_or_src)] )
          tmp_with_src <- merge(tmp_with_or_src, source_data, by.x='dmg_src_id', by.y='src_id', all.x=T, sort=F)
          names(tmp_with_src)[(NCOL(tmp_with_or_src)+1):NCOL(tmp_with_src)] <- paste0('dmg_', names(tmp_with_src)[(NCOL(tmp_with_or_src)+1):NCOL(tmp_with_src)] )
        } else {
          tmp_with_src <- tmp
        }
  
        Damage%<>%rbind(tmp_with_src)
      }
    }
    #filter to only shapefiles that have relevant damage information:
    
  }
  
  Damage$grading[grepl("Completely Destroyed",Damage$grading,fixed = T) |
                  grepl("Destroyed",Damage$grading,fixed = T)]       <-"Completely Destroyed"
  Damage$grading[grepl("Highly Damaged",Damage$grading,fixed = T)]             <-"Highly Damaged"
  Damage$grading[grepl("Moderately Damaged",Damage$grading,fixed = T)]         <-"Moderately Damaged"
  Damage$grading[grepl("Negligible to slight damage",Damage$grading,fixed = T)] <-"Slight Damage"
  Damage$grading[grepl("Possibly damaged",Damage$grading,fixed = T) | grepl("Possible damage",Damage$grading,fixed = T)]<- "Possible damage"
  Damage$grading[grepl("Not Affected",Damage$grading,fixed = T) |
                    grepl("No visible damage",Damage$grading,fixed = T) ]               <-"Not Affected"
  
  Damage$grading[Damage$grading=="Completely Destroyed"]<-"destroyed"
  Damage$grading[Damage$grading=="Severe Damage"]<-"severe"
  Damage$grading[Damage$grading=="Moderate Damage"]<-"moderate"
  Damage$grading[Damage$grading=="Slight Damage"]<-"slight"
  Damage$grading[Damage$grading=="Possible damage" ]<-"possible"
  Damage$grading[Damage$grading=="Not Affected"]<-"notaffected"
  
  #group the data based on source information
  if (is.null(Damage$or_source_nam)){
    grouping_vars <- c('source_nam', 'src_date')
  } else {
    grouping_vars <- c('or_source_nam', 'or_src_date', 'dmg_source_nam', 'dmg_src_date')
  }
  groupedDamage <- Damage %>%
    group_by(!!!syms(grouping_vars)) %>%
    mutate(group_id = cur_group_id()) 
  
  #create BD object to interpolate hazard over the BD data
  miniDam <- data.frame(grading=groupedDamage$grading, Longitude=groupedDamage$Longitude, Latitude=groupedDamage$Latitude, 
                        event=ev, source=groupedDamage$group_id, hazard='EQ', sdate=min(BDy_date), iso3=NA)
  
  BD <- new("BD",Damage=miniDam,ODD=ODDy)
  BD_new <- BD
  
  # check which groups appear to be missing at random. Criteria:
  #    - More than 50% of buildings under intensity 7 are unaffected
  #    - More than 5% of buildings over intensity 7 are unaffected
  
  if (length(grep('hazMean', names(BD_new@data)))==1){
    haz_max = BD_new@data$hazMean1
  } else {
    haz_max <- apply(BD_new@data[,grep('hazMean', names(BD_new@data))], 1, max, na.rm=T)
  }
  remove_i <- c()
  for (source_unique in unique(BD_new@data$source)){
    source_i <- which(BD_new@data$source == source_unique)
    mean_notaffected <- mean(BD_new@data$grading[source_i] == 'notaffected')
    if (!is.na(mean_notaffected)){
      if (mean_notaffected < 0.1){
        remove_i <- c(remove_i, source_i)
      }
    }
    next
    mean_notaffected_lowI <- mean(BD_new@data$grading[intersect(source_i, which(haz_max < 7))] == 'notaffected')
    mean_notaffected_highI <- mean(BD_new@data$grading[intersect(source_i, which(haz_max > 7))] == 'notaffected')
    print(paste(mean_notaffected_lowI, mean_notaffected_highI))
    if(!is.na(mean_notaffected_lowI)){
      if (mean_notaffected_lowI < 0.3){
        remove_i <- c(remove_i, source_i)
        next
      }
    }
    if (!is.na(mean_notaffected_highI)){
      if (mean_notaffected_highI < 0.01){
        remove_i <- c(remove_i, source_i)
        next
      }
    }
  }
  remove_i <- union(remove_i, which(is.na(BD_new$grading)))
  if (length(remove_i) > 0){
    BD_new@data <- BD_new@data[-remove_i,] 
    BD_new@coords <- BD_new@coords[-remove_i,] 
  }
  
  # col_values <- list('OpenBuildings'='white',
  #                    'notaffected'='green',
  #                    'possible'='antiquewhite',
  #                    'slight'='yellow',
  #                    'moderate'='orange',
  #                    'Damaged'='red',
  #                    'severe' = 'purple',
  #                    'destroyed' = 'blue')
  # 
  # print(ggplot() + 
  #   geom_point(data=data.table(BD@coords), aes(x=Longitude, y=Latitude),color='black', size=0.5) + theme_minimal() +
  #   geom_point(data=data.table(BD_new@coords), aes(x=Longitude, y=Latitude, color=BD_new$grading),size=0.5) + theme_minimal() +
  #   scale_color_manual(values=col_values))
  #   
  
  return(BD_new)
  
  # BD_new@data %>% group_by(source) %>%  
  #   group_map(~summarize(.x, prop_notaffected = ifelse(mean(grading == 'notaffected') > 0.5, T, F))) %>% unlist()
  
  
  # plot the grouped data:
  
  # # Use confidence as placeholder for source information rather than adding a new slot to BDy:
  # groupedDamage$Confidence <- apply(groupedDamage[,grouping_vars], 1, paste, collapse=', ')
  # groupedDamage$iso3 <- NA
  # BD_new <- new("BD",Damage=groupedDamage,ODD=ODDy)
  # ggplot(groupedDamage, aes(x = Longitude, y = Latitude, shape=grading, col = as.factor(Confidence))) +
  #   geom_point() + labs(col = "Origin Source, Date")
  
  groupedDamage$source_full <- apply(groupedDamage[,grouping_vars], 1, paste, collapse=', ')
  unique(groupedDamage[,c('source_full', 'group_id')])
  BD_plot <- BD@data
  BD_plot %<>% merge(unique(groupedDamage[,c('source_full', 'group_id')]), by.x='source', by.y='group_id', sort=F, all.x=T)
  if (length(grep('hazMean', names(BD_new@data)))==1){
    BD_plot$haz_max = BD_plot$hazMean1
  } else {
    BD_plot$haz_max <- apply(BD_plot[,grep('hazMean', names(BD_plot))], 1, max, na.rm=T)
  }
  BD_plot$haz_max <- round(BD_plot$haz_max * 2) / 2
  BD_plot %<>% group_by(source_full, haz_max) %>% summarise(prop_unaff = mean(grading == 'notaffected'),
                                                            prop_possible = mean(grading == 'possible' | grading == 'notaffected'),
                                                            n_tot = n())
  ggplot(BD_plot, aes(x = haz_max, y = 1-prop_unaff, group=source_full, col = as.factor(source_full))) +
    geom_line() + geom_point(aes(size=n_tot)) + ylab('Proportion Damaged') +
    geom_ribbon(aes(x = haz_max, ymin = 1-prop_unaff, ymax = 1 - prop_possible,
                    group=source_full, fill = as.factor(source_full)), alpha=0.1, color = NA) +
    labs(col = "Origin Source, Date, Damage Source, Date") +
    guides(fill = FALSE)
  
  treat_as_MAR <- function(x, y){
    if(NROW(x) > 50 & (sum(x$grading==notaffected)/length(x$grading) < 0.8)){
      return(T)
    } else {
      return(F)
    }
  }
  
  # unlist(groupedDamage %>%
  #   group_map(~summarize(.x, prop_notaffected = mean(grading == 'notaffected') > 0.5)))
  
  plot(Damage$Longitude, Damage$Latitude)
  points(Damage$Longitude[which(Damage$or_src_id==1 & Damage$dmg_src_id==2 & Damage$det_method=='Photo-interpretation' & Damage$obj_type=='11-Residential Buildings')], Damage$Latitude[which(Damage$or_src_id==1 & Damage$dmg_src_id==2 & Damage$det_method=='Photo-interpretation' & Damage$obj_type=='11-Residential Buildings')], col='red')
  shapefiles <- shapefiles[contains_dam_data]
  
  shape_data <- list()
  treat_as_MAR <- rep(NA, length(shape_data))
  for (i in 1:length(shapefiles)){
    shape_data[[i]] <- st_read(shapefiles[i])
    shape_data[[i]] <- st_transform(shape_data[[i]], crs = 4326)
    points_in_shape <- st_within(st_as_sf(as.data.frame(BDy@coords), coords = c("Longitude", "Latitude"), crs=st_crs(shape_data[[i]])), shape_data[[i]], sparse = FALSE)
    treat_as_MAR[i] <- ifelse(sum(points_in_shape) < 10 | sum(BDy$grading[points_in_shape]=='notaffected')/sum(points_in_shape) < 0.2, F, T)
  }
  
  file_conn <- file(paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/BD_creation_notes'), open = "a")
  writeLines(paste0(c("Regions treated as MAR:", which(treat_as_MAR)), collapse=' '), file_conn)
  writeLines(paste0(c("Regions removed:", which(!treat_as_MAR)), collapse=' '), file_conn)
  close(file_conn) 
  
  if(sum(treat_as_MAR)==0){
    return(NULL)
  }
  
  print(plot_sat_damage_shps(BDy, shape_data, ODDy=ODDy, plot_marks=treat_as_MAR))
  
  shapefiles <- shapefiles[which(treat_as_MAR)]
  shape_data <- shape_data[which(treat_as_MAR)]
  
  
  miniDam_old <- data.frame(grading=BDy$grading, Longitude=BDy$Longitude, Latitude=BDy$Latitude, event=ev, Confidence=BDy$Confidence, hazard='EQ', sdate=BDy@hazdates[1], iso3=BDy@hazdates[1])
  
  building_locs_old <- st_as_sf(miniDam_old, coords = c("Longitude", "Latitude"))
  building_locs_old$poly_int <- FALSE
  # Loop through each shape in shape_data
  for (j in 1:length(shape_data)){
    # Check if any points in building_locs_sf are within the current polygon
    shape_data[[j]] <- st_set_crs(shape_data[[j]], st_crs(building_locs_old))
    within <- st_within(building_locs_old, shape_data[[j]], sparse = FALSE)
    # Update poly_int for points within the polygon
    building_locs_old$poly_int[within] <- TRUE
    
  }
  miniDam_old <- miniDam_old[which(building_locs_old$poly_int==T),]
  
  # p <- ggplot()
  # for (i in 1:length(shape_data)){
  #   p <- p + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red")
  # }
  # p <- p + geom_point(data=data.frame(building_locs_miss), aes(x=longitude, y=latitude, col='OpenBuildings'), size=0.1) + 
  #   geom_point(data=miniDam_old, aes(x=miniDam_old$Longitude, y=miniDam_old$Latitude, color=miniDam_old$grading), size=0.5) + 
  #   scale_color_manual(values = col_values)
  # p
  # 
  # i <- 1
  # p <- ggplot() + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red") + xlim(st_bbox(shape_data[[i]])[c(1,3)])+ylim(st_bbox(shape_data[[i]])[c(2,4)]) +
  #   geom_point(data=miniDam_old, aes(x=miniDam_old$Longitude, y=miniDam_old$Latitude, color=miniDam_old$grading), size=1) +
  #   scale_color_manual(values = col_values)
  # 
  # p_after <- p + geom_point(data=data.frame(building_locs_miss), aes(x=longitude, y=latitude, col='OpenBuildings'), size=0.1) + 
  #   geom_point(data=BDy@data[rev(order(BDy@data$grading)),], aes(x=BDy$Longitude, y=BDy$Latitude, color=BDy$grading), size=1)
  # 
  # grid.arrange(p, p_after, nrow=2)
  
  miniDam <- miniDam_old
  
  # folderin_ODD <- "/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/"
  # ufiles_ODD <- list.files(path=folderin_ODD,pattern=Model$haz,recursive = T,ignore.case = T)
  # file_match_ODD <- ufiles_ODD[which(as.numeric(sub(".+_(.+)$", "\\1", ufiles_ODD))==event_i)]
  # 
  # ODDy <- readRDS(paste0(folderin_ODD, file_match_ODD))
  
  BD_new <- new("BD",Damage=miniDam,ODD=ODDy)
  
  return(BD_new)
}  

