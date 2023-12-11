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
 
  event_i <- 67
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
  }
  
  shape_data <- list()
  for (i in 1:length(shapefiles)){
    shape_data[[i]] <- st_read(shapefiles[i])
    shape_data[[i]] <- st_transform(shape_data[[i]], crs = 4326)
  }
  
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
  
  p <- ggplot() 
  for (i in 1:length(shape_data)){
    p <- p + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red")
  }
  p <- p + #geom_sf(data=st_as_sf(data.table(BDy@coords), coords = c("Longitude", "Latitude"), crs=4326), size=0.5) +#,  #color=BDy$grading) + 
    geom_point(data=data.table(BDy@coords), aes(x=Longitude, y=Latitude, color=BDy$grading),size=0.5) + theme_minimal() +
    scale_color_manual(values = list('possible'='yellow', 'Damaged'='red', 'notaffected'='blue', 'destroyed'='purple')) + 
    xlim(combined_bbox[1], combined_bbox[3]) + ylim(combined_bbox[2], combined_bbox[4]) 
    #xlim(85.2, 85.4) + ylim(27.62, 27.75)
  p
  
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
      stop(paste("Event:", event_i, ", Missing", length(missing_quadkeys)/length(quad_keys)*100, "percent of quad keys, not adding building data."))
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
  
  p <- ggplot()
  for (i in 1:length(shape_data)){
    p <- p + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red")
  }
  p <- p + geom_point(data=data.frame(building_locs_miss), aes(x=longitude, y=latitude, col='OpenBuildings'), size=0.1) + 
    geom_point(data=BDy@data, aes(x=BDy$Longitude, y=BDy$Latitude, color=BDy$grading), size=0.5) + 
    scale_color_manual(values = list('possible'='yellow', 'Damaged'='red', 'notaffected'='blue', 'OpenBuildings'='green'))
  p
  
  miniDam <- data.frame(grading=BDy$grading, Longitude=BDy$Longitude, Latitude=BDy$Latitude, event=ev, Confidence=BDy$Confidence, hazard='EQ', sdate=BDy@hazdates[1], iso3=BDy@hazdates[1])
  miniDam <- rbind(miniDam, data.frame(grading='notaffected', Longitude=building_locs_miss$longitude, Latitude=building_locs_miss$latitude, 
                                       event=ev, Confidence=NA, hazard='EQ', sdate=BDy@hazdates[1], iso3=BDy@hazdates[1]))
  
  # folderin_ODD <- "/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/ODDobjects/"
  # ufiles_ODD <- list.files(path=folderin_ODD,pattern=Model$haz,recursive = T,ignore.case = T)
  # file_match_ODD <- ufiles_ODD[which(as.numeric(sub(".+_(.+)$", "\\1", ufiles_ODD))==event_i)]
  # 
  # ODDy <- readRDS(paste0(folderin_ODD, file_match_ODD))
  
  BD_new <- new("BD",Damage=miniDam,ODD=ODDy)
  
  return(BD_new)
  
  plot_bd_df_int <- data.frame(intensity=numeric(), event_i=integer(), bd_new=logical(), obs_prop_bd=numeric(), tot_obs=integer())
  for(int in seq(4.25,9.25, 0.5)){
    plot_bd_df_int %<>% add_row(intensity=int, event_i=event_i, bd_new=T, 
                                obs_prop_bd=1-length(which(BD_new$grading[which(BD_new$hazMean1>int &BD_new$hazMean1<(int+0.2))]=='notaffected'))/length(which(BD_new$hazMean1>int & BD_new$hazMean1<(int+0.2))),
                                tot_obs=length(which(BD_new$hazMean1>int & BD_new$hazMean1<(int+0.2))))
    plot_bd_df_int %<>% add_row(intensity=int, event_i=event_i, bd_new=F, 
                                obs_prop_bd=1-length(which(BDy$grading[which(BDy$hazMean1>int &BDy$hazMean1<(int+0.2))]=='notaffected'))/length(which(BDy$hazMean1>int & BDy$hazMean1<(int+0.2))),
                                tot_obs=length(which(BDy$hazMean1>int & BDy$hazMean1<(int+0.2))))
  }
  
  ggplot(plot_bd_df_int, aes(x=intensity, y=obs_prop_bd, group=bd_new, color=bd_new)) + 
    geom_line() + ylab('Proportion Buildings Damaged') + 
    geom_line(aes(x=intensity, y=sqrt(plot_bd_df_int$obs_prop_bd*(1-plot_bd_df_int$obs_prop_bd)/plot_bd_df_int$tot_obs)))
  
  plot_bd_df_int$std_dev <- sqrt(plot_bd_df_int$obs_prop_bd * (1 - plot_bd_df_int$obs_prop_bd) / plot_bd_df_int$tot_obs)
  
  # Create a ggplot with shaded region
  ggplot(plot_bd_df_int, aes(x = intensity, y = obs_prop_bd, group = bd_new, color = bd_new)) +
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




