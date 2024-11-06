
#------------------------------------------------------------------------------------------
#---------------------------Preliminary Functions------------------------------------------
#------------------------------------------------------------------------------------------

library(dplyr, include.only = c("revalue"))


extract_country_date <- function(filename) {
  # Extract date and country using regex
  date <- substr(filename, 3, 10)  # Extract date part (positions 3 to 10)
  country <- substr(filename, 11, 13)  # Extract country code (positions 11 to 13)
  
  # Convert date to desired format (YYYY-MM-DD)
  formatted_date <- format(as.Date(date, "%Y%m%d"), "%Y-%m-%d")
  
  # Combine formatted date and country
  result <- paste(formatted_date, country, sep = ", ")
  
  return(result)
}

input_folder <- 'IIDIPUS_Input_Alternatives/Aug24/BDobjects/'
BD_count_per_event <- function(input_folder){
  ufiles <- list.files(path=paste0(dir, input_folder),recursive = T,ignore.case = T)
  BD_dat_all <- data.frame(event_name=character(),
                           grading=character(),
                           pixel=integer())
  for (file in ufiles){
    BD <- readBD(paste0(dir, input_folder, file))
    BD_dat_all %<>% rbind(data.frame(event_name=extract_country_date(file), grading=BD$grading, pixel=BD$spatial_pixel))
  }
  
  BD_dat_all$grading <- factor(BD_dat_all$grading, levels = c(
    "destroyed", "severe", "Damaged","moderate","possible","notaffected"
  ))
  
  BD_dat_all$grading <- plyr::revalue(BD_dat_all$grading, c("destroyed" = "Destroyed",
                                                      'severe' = 'Severe damage',
                                                      'Damaged' = 'Damaged',
                                                      'moderate' = 'Moderate damage',
                                                      'possible' = 'Possibly damaged',
                                                      'notaffected' = 'Unaffected'))
  
  p0 <- ggplot(BD_dat_all, aes(x = as.factor(event_name), fill = grading)) +
    geom_bar(position = "stack") +
    labs(x = "Event", y = "Building Count", fill = "Grading") +
    theme_minimal() + 
    scale_fill_manual(
      values = c('Unaffected'='#99D097',
                 'Possibly damaged'='#ef9f00',
                 'Moderate damage'='#EF6200',
                 'Damaged'='#D50000',
                 'Severe damage' = '#b00041',
                 'Destroyed' = '#80002a',
                 'Missing building footprints'='#CECECE')
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Times New Roman", size=10),
          axis.title = element_text(family = "Times New Roman", size=12),  # Optional: Set specifically for axis titles
          plot.title = element_text(family = "Times New Roman"),   # Optional: Set specifically for plot titles
          legend.text = element_text(family = "Times New Roman", size=11),    # Legend text
          legend.title = element_text(family = "Times New Roman", size=12) )
  p0
  return(p0)
  #Proportion 'possibly damaged':
  # MAR_first_cull <- BD_dat_all %>% group_by(event_id) %>%
  #   filter(sum(grading == "notaffected") > 0)
  # NROW(MAR_first_cull) / MAR_first_cull %>% n_groups()
  # mean(MAR_first_cull$grading=='possible')
}

BD_count_per_event(input_folder)
#PointSummary.pdf, 10 x 6

collect_BD_shapefiles <- function(BDy){
  BDy_date <- BDy@hazdates
  
  if (as.Date("2013-09-24") %in% BDy_date){ #UNOSAT
    #cfiles<-list.files(path=paste0(dir,"UNOSAT_Damage/"),pattern="AnalysisExtent",recursive = T)
    
    shapefiles <- '/home/manderso/Documents/GitHub/ODDRIN/UNOSAT_Damage/EQ20130924PAK_UNOSAT/Damage_Sites.shp'
    ev <- 'EQ20130924PAK'
    
    stop('Handle UNOSAT EQ20130924PAK event.')
    #buildings_in_MAR_polys <- which(BDy@coords[,1] > 65 & BDy@coords[,1] < 65.65 & BDy@coords[,2] < 27.2)
    #BDy@data <- BDy@data[buildings_in_MAR_polys,]
    #BDy@coords <- BDy@coords[buildings_in_MAR_polys,]
    #return(BDy)
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
  
  shape_data <- list()
  #treat_as_MAR <- rep(NA, length(shape_data))
  for (i in 1:length(shapefiles)){
    shape_data[[i]] <- st_read(shapefiles[i])
    shape_data[[i]] <- st_transform(shape_data[[i]], crs = 4326)
  }
  return(shape_data)
}

plot_BD_with_shapefiles <- function(BDy){
  
  shape_data <- collect_BD_shapefiles(BDy)
  #shapefiles <- shapefiles[which(treat_as_MAR)]
  #shape_data <- shape_data[which(treat_as_MAR)]
  
  combined_bbox <- st_bbox(shape_data[[1]])
  
  # Iterate through each element in shape_data starting from the second element
  for (i in seq_along(shape_data)[-1]) {
    
    # Update the combined bounding box to include the current object's bounding box
    combined_bbox["xmin"] <- min(combined_bbox["xmin"], st_bbox(shape_data[[i]])["xmin"])
    combined_bbox["ymin"] <- min(combined_bbox["ymin"], st_bbox(shape_data[[i]])["ymin"])
    combined_bbox["xmax"] <- max(combined_bbox["xmax"], st_bbox(shape_data[[i]])["xmax"])
    combined_bbox["ymax"] <- max(combined_bbox["ymax"], st_bbox(shape_data[[i]])["ymax"])
  }
  
  col_values <- list('OpenBuildings'='white',
                     'notaffected'='blue',
                     'possible'='yellow', 
                     'moderate'='orange',
                     'Damaged'='red',
                     'severe' = 'hotpink',
                     'destroyed' = 'violetred')
  
  BDy_coords <- as.data.frame(crds(BDy))
  colnames(BDy_coords) <- c('Longitude', 'Latitude')
  
  p <- ggplot() 
  p <- p + #geom_sf(data=st_as_sf(data.table(BDy@coords), coords = c("Longitude", "Latitude"), crs=4326), size=0.5) +#,  #color=BDy$grading) + 
    geom_point(data=BDy_coords, aes(x=Longitude, y=Latitude, color=BDy$grading),size=0.5) + theme_minimal() +
    scale_color_manual(values = col_values) + 
    xlim(combined_bbox[1], combined_bbox[3]) + ylim(combined_bbox[2], combined_bbox[4]) 
  #xlim(85.2, 85.4) + ylim(27.62, 27.75)
  return(p)
  
  for (i in 1:length(shape_data)){
    p <- p + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red")
  }
  p <- p + #geom_sf(data=st_as_sf(data.table(BDy@coords), coords = c("Longitude", "Latitude"), crs=4326), size=0.5) +#,  #color=BDy$grading) + 
    geom_point(data=BDy_coords, aes(x=Longitude, y=Latitude, color=BDy$grading),size=0.5) + theme_minimal() +
    scale_color_manual(values = col_values) + 
    xlim(combined_bbox[1], combined_bbox[3]) + ylim(combined_bbox[2], combined_bbox[4]) 
  #xlim(85.2, 85.4) + ylim(27.62, 27.75)
  p
  
  # p <- list()
  # for (i in 1:length(shape_data)){
  #   p[[i]] <- ggplot() + geom_sf(data = st_as_sf(shape_data[[i]]), color = "red") + xlim(st_bbox(shape_data[[i]])[c(1,3)])+ylim(st_bbox(shape_data[[i]])[c(2,4)]) +
  #     geom_point(aes(x=BDy_coords$Longitude, y=BDy_coords$Latitude, color=BDy$grading), size=1) + 
  #     scale_color_manual(values = col_values) +  theme(legend.position = "none",  # Remove legend
  #                                                      axis.title.x = element_blank(),  # Remove x-axis title
  #                                                      axis.title.y = element_blank(),  # Remove y-axis title
  #                                                      axis.text.x = element_blank(),   # Remove x-axis labels
  #                                                      axis.text.y = element_blank()) 
  # }
  # do.call(grid.arrange,p)
}

library(foreign)

addBDsource <- function(BDy){
  BDy_date <- BDy@hazdates[1]
  
  
  if (as.Date("2013-09-24") %in% BDy_date){ #UNOSAT
    shapefiles <- '/home/manderso/Documents/GitHub/ODDRIN/UNOSAT_Damage/EQ20130924PAK_UNOSAT/Damage_Sites.shp'
    ev <- 'EQ20130924PAK'
    tmp<-st_read(shapefiles,quiet = T)
    
    tmp$geometry%<>%st_transform(crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
    # Remove unwanted grading terms
    tmp%<>%filter(!Main_Damag%in%
                    c("Null","Unknown","Not Applicable"))
    if(nrow(tmp)==0) next
    # convert point geometries into Longitude/Latitude
    ttt<-st_coordinates(tmp$geometry)
    tmp$Longitude<-ttt[,1]; tmp$Latitude<-ttt[,2]
    tmp$geometry<-NULL
    tmp <- tmp[,c('Longitude','Latitude', 'SensorID', 'Main_Damag')]
    names(tmp)[which(names(tmp)=='Main_Damag')] = 'grading'
    Damage <- tmp
    
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
      #print(subdir_dam_files)
      if (length(subdir_dam_files) == 0){contains_dam_data[j] <- F}
    }
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
    #filter to only shapefiles that have relevant damage information:
  }
  
  DamageNoDuplicates <- distinct(Damage[,-grep('grading', names(Damage))])
  plot_df <- as.data.frame(BDy, geom='XY',na.rm=F)
  plot_df$order <- 1:NROW(plot_df)
  names(plot_df)[which(names(plot_df)=='x')] = 'Longitude'
  names(plot_df)[which(names(plot_df)=='y')] = 'Latitude'
  plot_df_w_source <- merge(plot_df, DamageNoDuplicates, by=c('Longitude', 'Latitude'), all.x=T)
  plot_df_w_source <- plot_df_w_source[order(plot_df_w_source$order),]
  return(plot_df_w_source[, -grep('order', names(plot_df_w_source))])
}

  #------------------------------------------------------------------------------------------
  #----------------------------------Figures A and B-----------------------------------------
  #-------- Haiti 2018-10-06 earthquake without and with Bing Building Footprints------------
  #------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------

#Figure (a)
library(plyr)
BDy_A <- readBD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Aug24/BDobjects/EQ20181006HTI_97')
event_id_A=97
plot_df_A <- as.data.frame(BDy_A, geom='XY')
names(plot_df_A)[which(names(plot_df_A)=='x')] = 'Longitude'
names(plot_df_A)[which(names(plot_df_A)=='y')] = 'Latitude'

col_values <- list('Unaffected'='#99D097',
                   'Possibly damaged'='#EF8900',
                   'Moderate damage'='#EF6200',
                   'Damaged'='#D50000',
                   'Severe damage' = '#ce1256',
                   'Destroyed' = '#80002a',
                   'Missing building footprints'='#CECECE')

plot_df_A$grading <- factor(plot_df_A$grading, levels=c('notaffected', 'possible', 'moderate',
                                                        'Damaged', 'severe', 'destroyed', 'OpenBuildings'))
plot_df_A$grading <- revalue(plot_df_A$grading, c("notaffected" = "Unaffected",
                                                  'possible' = 'Possibly damaged',
                                                  'moderate' = 'Moderate damage',
                                                  'Damaged' = 'Damaged',
                                                  'severe' = 'Severe damage',
                                                  'destroyed' = 'Destroyed',
                                                  'OpenBuildings' = 'Missing building footprints'))

shape_data_A <- collect_BD_shapefiles(BDy_A)

points_sf_A <- st_as_sf(plot_df_A, coords = c("Longitude", "Latitude"), crs = 4326)

# Filter rows where points are contained in either of the two polygons
filtered_df_A <- plot_df_A[st_within(points_sf_A, st_as_sf(shape_data_A[[6]]), sparse = FALSE) |
                             st_within(points_sf_A, st_as_sf(shape_data_A[[7]]), sparse = FALSE),]

# If you want it as a dataframe again, drop the sf class
filtered_df_A <- st_drop_geometry(filtered_df_A)

p1 <- ggplot() +
  geom_sf(data = st_as_sf(shape_data_A[[6]]), color = "red", fill='grey', alpha=0.1) + 
  geom_sf(data = st_as_sf(shape_data_A[[7]]), color = "red", fill='grey', alpha=0.1) +
  geom_point(data=filtered_df_A, aes(x=Longitude, y=Latitude, color=grading),size=1) +
  theme_minimal() +
  scale_color_manual(values = col_values) + 
  scale_x_continuous(breaks = seq(-72.86,-72.76, by = 0.04), limits=c(-72.88, -72.75)) +
  ylim(19.9, 19.955) +
  labs(color = "Damage Classification") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3), 
        plot.title=element_text(hjust=-.2)) + 
  ggtitle('(a)') 
p1

# figure b:
event_id=97
buildingFootprints_B <- getBingBuildingsGlobal(BDy_A, event_id, aggregate=F)


plot_df_wBuildingFootprints_B <- as.data.frame(buildingFootprints_B)
colnames(plot_df_wBuildingFootprints_B) <- c('Longitude', 'Latitude')
plot_df_wBuildingFootprints_B$grading <- 'Missing building footprints'
plot_df_wBuildingFootprints_B %<>% rbind(plot_df_A[,c('Longitude','Latitude', 'grading')])

plot_df_wBuildingFootprints_B$grading <- factor(plot_df_wBuildingFootprints_B$grading, 
                                                levels=c('Unaffected', 'Possibly damaged', 'Moderate damage',
                                                         'Damaged', 'Severe damage', 'Destroyed', 'Missing building footprints'))

points_sf_B <- st_as_sf(plot_df_wBuildingFootprints_B, coords = c("Longitude", "Latitude"), crs = 4326)

# Filter rows where points are contained in either of the two polygons
filtered_df_B <- plot_df_wBuildingFootprints_B[st_within(points_sf_B, st_as_sf(shape_data_A[[6]]), sparse = FALSE) |
                                                 st_within(points_sf_B, st_as_sf(shape_data_A[[7]]), sparse = FALSE),]

# If you want it as a dataframe again, drop the sf class
filtered_df_B <- st_drop_geometry(filtered_df_B)

p2 <- ggplot() + 
  geom_sf(data = st_as_sf(shape_data_A[[6]]), color = "red", fill='grey', alpha=0.1) + 
  geom_sf(data = st_as_sf(shape_data_A[[7]]), color = "red", fill='grey', alpha=0.1) + 
  geom_point(aes(x=-72.81039, y=19.94458, color=as.factor('Unaffected')), size=0.1) +  # just so unaffected appears in legend. Should be hidden
  geom_point(data=filtered_df_B %>% filter(grading!='Missing building footprints'), aes(x=Longitude, y=Latitude, color=grading),size=1) + # show first in legend
  geom_point(data=filtered_df_B %>% filter(grading=='Missing building footprints'), aes(x=Longitude, y=Latitude, color=grading),size=0.5) + 
  geom_point(data=filtered_df_B %>% filter(grading!='Missing building footprints'), aes(x=Longitude, y=Latitude, color=grading),size=1) + #repeat so the points are on top
  theme_minimal() +
  scale_color_manual(values = unlist(col_values)) + 
  scale_x_continuous(breaks = seq(-72.86,-72.76, by = 0.04), limits=c(-72.88, -72.75)) +
  ylim(19.9, 19.955)  +
  labs(color = "Damage Classification") + 
  guides(color = guide_legend(override.aes = list(size = 1.5), title.position='top', title.hjust=0.5)) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3), 
        plot.title=element_text(hjust=-.2, family = "Liberation Serif"),
        #axis.text.x = element_text(angle = 45, hjust = 1, family = "Liberation Serif", size=10),
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=12) ) + 
  ggtitle('(a)      Haiti, 06-10-2018') + xlab('Longitude') + ylab('Latitude')

p2

#------------------------------------------------------------------------------------------
#----------------------------------Figures C and D-----------------------------------------
#-------- Nepal 2015-04-25 earthquake without and with Bing Building Footprints------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

BDy <- readBD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Aug24/BDobjects/EQ20150425NPL_31')

event_id=31
plot_df <- as.data.frame(BDy, geom='XY')
names(plot_df)[which(names(plot_df)=='x')] = 'Longitude'
names(plot_df)[which(names(plot_df)=='y')] = 'Latitude'

col_values <- list('Missing building footprints'='#CECECE',
                   'Unaffected'='#99D097',
                   'Possibly damaged'='#EF8900',
                   'Moderate damage'='#EF6200',
                   'Damaged'='#D50000',
                   'Severe damage' = '#ce1256',
                   'Destroyed' = '#80002a')

plot_df$grading <- factor(plot_df$grading, levels=c('OpenBuildings','notaffected', 'possible', 'moderate',
                                                    'Damaged', 'severe', 'destroyed'))
plot_df$grading <- revalue(plot_df$grading, c('OpenBuildings' = 'Missing building footprints',
                                              "notaffected" = "Unaffected",
                                              'possible' = 'Possibly damaged',
                                              'moderate' = 'Moderate damage',
                                              'Damaged' = 'Damaged',
                                              'severe' = 'Severe damage',
                                              'destroyed' = 'Destroyed'))

plot_df %<>% filter(Longitude > 83.95 & Longitude < 84.025 & Latitude > 28.2 & Latitude < 28.2425)

shape_data <- collect_BD_shapefiles(BDy)

p3 <- ggplot() +
  geom_sf(data = st_as_sf(shape_data[[3]]), color = "red", fill='grey', alpha=0.1) + 
  geom_point(data=plot_df %>% filter(grading=='Unaffected'), aes(x=Longitude, y=Latitude, color=grading),size=0.3) + theme_minimal() +
  geom_point(data=plot_df %>% filter(grading!='Unaffected'), aes(x=Longitude, y=Latitude, color=grading),size=1.2) + theme_minimal() +
  scale_color_manual(values = col_values) + 
  scale_x_continuous(breaks = seq(83.96,84.02, by = 0.02), limits=c(83.955, 84.025)) +
  ylim(28.196, 28.2425) +
  labs(color = "Damage Classification") + 
  guides(color = guide_legend(override.aes = list(size = 1.5))) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3), 
        plot.title=element_text(hjust=-.2))  + 
  ggtitle('(c)') 
p3

# figure (d):

buildingFootprints <- getBingBuildingsGlobal(BDy, event_id, aggregate=F)
buildingFootprints %<>% as.data.frame()
buildingFootprints %<>% filter(Longitude > 83.95 & Longitude < 84.025 & Latitude > 28.2 & Latitude < 28.2425)
#plot(buildingFootprints$Longitude, buildingFootprints$Latitude)

buildingFootprints_cleaned <- buildingFootprints
# Step 2: Loop over each point in plot_df to find the closest point in buildingFootprints
for (i in 1:nrow(plot_df)) {
  if (i %% 100 == 0) print(i)
  # Extract the current point in plot_df
  current_point <- plot_df[i, c("Longitude", "Latitude")]
  
  # Calculate distances between current point and all points in buildingFootprints
  distances <- distm(current_point, buildingFootprints_cleaned[, c("Longitude", "Latitude")])
  
  # Find the index of the closest point
  closest_index <- which.min(distances)
  
  buildingFootprints_cleaned[-closest_index,]
}

plot_df2 <- rbind(plot_df[,c('Longitude', 'Latitude', 'grading')], buildingFootprints_cleaned %>% add_column(grading='Missing building footprints'))

points_sf <- st_as_sf(plot_df2, coords = c("Longitude", "Latitude"), crs = 4326)

# Filter rows where points are contained in either of the two polygons
filtered_df <- plot_df2[st_within(points_sf, st_as_sf(shape_data[[3]]), sparse = FALSE),]

# If you want it as a dataframe again, drop the sf class
filtered_df <- st_drop_geometry(filtered_df)

p4 <- ggplot() +
  geom_sf(data = st_as_sf(shape_data[[3]]), color = "red", fill='grey', alpha=0.1) + 
  geom_point(data=filtered_df %>% filter(grading=='Missing building footprints'), aes(x=Longitude, y=Latitude, color=grading),size=0.3) + theme_minimal() +
  geom_point(data=filtered_df %>% filter(grading=='Unaffected'), aes(x=Longitude, y=Latitude, color=grading),size=0.3) + theme_minimal() +
  geom_point(data=filtered_df %>% filter(grading!='Unaffected' & grading != 'Missing building footprints'), aes(x=Longitude, y=Latitude, color=grading),size=1.2) + theme_minimal() +
  scale_color_manual(values = col_values) + 
  scale_x_continuous(breaks = seq(83.96,84.02, by = 0.02), limits=c(83.955, 84.025)) +
  ylim(28.196, 28.2425) +
  labs(color = "Damage Classification") + 
  guides(color = guide_legend(override.aes = list(size = 1.5))) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3), 
        plot.title=element_text(hjust=-.2, family = "Liberation Serif"),
        #axis.text.x = element_text(angle = 45, hjust = 1, family = "Liberation Serif", size=10),
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=12)) + 
  ggtitle('(b)   Nepal, 25-04-2015') 

p4

legend <- get_plot_component(p2 +
                               theme(legend.position="bottom"), 'guide-box', return_all=T)[[3]]


plot_grid( plot_grid( p1  + theme(legend.position="none") , p2  + theme(legend.position="none"), 
                      p3  + theme(legend.position="none"), p4 + theme(legend.position="none"),align = 'vh', nrow = 2),legend, ncol = 1, rel_heights=c(1,0.1))

plot_grid( plot_grid(p2  + theme(legend.position="none"), 
                     p4 + theme(legend.position="none"),
                     align = 'vh', nrow = 1),legend, ncol = 1, rel_heights=c(1,0.1))

#BuildPoint.pdf, 10 x 4.25

#------------------------------------------------------------------------------------------
#--------------------------------Figures E, F and G----------------------------------------
#---------Comparison of different satellite imagery sources, IDN event 2022-11-21----------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------


# figure (e):

BDy <- readBD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Aug24/BDobjects/EQ20221121IDN_168')
event_id = 168

shape_data <- collect_BD_shapefiles(BDy)

plot_df <- addBDsource(BDy)

col_values <- list('Missing building footprints'='#CECECE',
                   'Unaffected'='#99D097',
                   'Possibly damaged'='#EF8900',
                   'Moderate damage'='#EF6200',
                   'Damaged'='#D50000',
                   'Severe damage' = '#ce1256',
                   'Destroyed' = '#80002a')

plot_df$grading <- factor(plot_df$grading, levels=c('OpenBuildings','notaffected', 'possible', 'moderate',
                                                    'Damaged', 'severe', 'destroyed'))
plot_df$grading <- revalue(plot_df$grading, c('OpenBuildings' = 'Missing building footprints',
                                              "notaffected" = "Unaffected",
                                              'possible' = 'Possibly damaged',
                                              'moderate' = 'Moderate damage',
                                              'Damaged' = 'Damaged',
                                              'severe' = 'Severe damage',
                                              'destroyed' = 'Destroyed'))

#plot_df %<>% filter(Longitude > 83.95 & Longitude < 84.025 & Latitude > 28.2 & Latitude < 28.2425)

p5 <- ggplot() +
  geom_sf(data = st_as_sf(shape_data[[1]]), color = "red", fill='grey', alpha=0.1) + 
  geom_sf(data = st_as_sf(shape_data[[4]]), color = "red", fill='grey', alpha=0.1) + 
  geom_point(data=plot_df %>% filter(grading=='Unaffected'), aes(x=Longitude, y=Latitude, color=grading),size=1) + theme_minimal() +
  geom_point(data=plot_df %>% filter(grading!='Unaffected'), aes(x=Longitude, y=Latitude, color=grading),size=1) + theme(legend.position='bottom') +
  #geom_point(data=plot_df %>% filter(or_source_nam=='WorldView-2'), aes(x=Longitude, y=Latitude, color=grading),size=1) + theme_minimal() +
  scale_color_manual(values = col_values) + 
  labs(color = "Damage Classification") + 
  guides(color = guide_legend(override.aes = list(size = 1.5), title.position='top', title.hjust=0.5)) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        plot.title=element_text(hjust=-.2, family = "Liberation Serif"),
        #axis.text.x = element_text(angle = 45, hjust = 1, family = "Liberation Serif", size=10),
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=11)) + 
  geom_point(aes(x=107.0446, y=-6.82125), shape="\u2605", size=10, fill='blue', color='blue') +
  scale_x_continuous(breaks = seq(107.02,107.18, by = 0.04), limits=c(107.035, 107.17)) +
  scale_y_continuous(breaks = seq(-6.9,-6.74, by = 0.04), limits=c(-6.88, -6.75)) +
  #xlim(107.035, 107.17) + ylim(-6.88, -6.75) +
  ggtitle('(a) All Damage Data')


p6 <- ggplot() +
  geom_sf(data = st_as_sf(shape_data[[1]]), color = "red", fill='grey', alpha=0.1) + 
  geom_sf(data = st_as_sf(shape_data[[4]]), color = "red", fill='grey', alpha=0.1) + 
  geom_point(data=plot_df %>% filter(grading=='Unaffected' & or_source_nam=='Open Street Map'), aes(x=Longitude, y=Latitude, color=grading),size=1) + theme_minimal() +
  geom_point(data=plot_df %>% filter(grading!='Unaffected'& or_source_nam=='Open Street Map'), aes(x=Longitude, y=Latitude, color=grading),size=1) + theme_minimal() +
  scale_color_manual(values = col_values) + 
  scale_x_continuous(breaks = seq(107.02,107.18, by = 0.04), limits=c(107.035, 107.17)) +
  scale_y_continuous(breaks = seq(-6.9,-6.74, by = 0.04), limits=c(-6.88, -6.75)) +
  labs(color = "Damage Classification") + 
  guides(color = guide_legend(override.aes = list(size = 1.5))) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        plot.title=element_text(hjust=-1.8, family = "Liberation Serif"),
        #axis.text.x = element_text(angle = 45, hjust = 1, family = "Liberation Serif", size=10),
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=11)) + 
  geom_point(aes(x=107.0446, y=-6.82125), shape="\u2605", size=10, fill='blue', color='blue') +
  ggtitle('(b) Pre-event source: Open Street Map')


p7 <- ggplot() +
  geom_sf(data = st_as_sf(shape_data[[1]]), color = "red", fill='grey', alpha=0.1) + 
  geom_sf(data = st_as_sf(shape_data[[4]]), color = "red", fill='grey', alpha=0.1) + 
  geom_point(data=plot_df %>% filter(grading=='Unaffected' & or_source_nam=='WorldView-2'), aes(x=Longitude, y=Latitude, color=grading),size=1) + theme_minimal() +
  geom_point(data=plot_df %>% filter(grading!='Unaffected'& or_source_nam=='WorldView-2'), aes(x=Longitude, y=Latitude, color=grading),size=1) + 
  scale_color_manual(values = col_values) + 
  scale_x_continuous(breaks = seq(107.02,107.18, by = 0.04), limits=c(107.035, 107.17)) +
  scale_y_continuous(breaks = seq(-6.9,-6.74, by = 0.04), limits=c(-6.88, -6.75)) +
  labs(color = "Damage Classification") + 
  guides(color = guide_legend(override.aes = list(size = 1.5))) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        plot.title=element_text(hjust=-0.8, family = "Liberation Serif"),
        #axis.text.x = element_text(angle = 45, hjust = 1, family = "Liberation Serif", size=10),
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=11)) + 
  geom_point(aes(x=107.0446, y=-6.82125), shape="\u2605", size=10, fill='blue', color='blue')  +
  ggtitle('(c) Pre-event source: WorldView-2')



legend <- get_plot_component(p5 +
                               theme(legend.position="bottom"), 'guide-box', return_all=T)[[3]]

plot_grid( plot_grid( p5  + theme(legend.position="none"), p6  + theme(legend.position="none"), 
                      p7  + theme(legend.position="none"), align = 'vh', nrow = 1),legend, ncol = 1, rel_heights=c(1,0.1))


#12 x 5, IDNcopernicus.pdf


#------------------------------------------------------------------------------------------
#---------------------------------Figures H and I------------------------------------------
#--------Intensity vs Damage Proportion for original and filtered data---------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

input_folder <- 'IIDIPUS_Input_Alternatives/Aug24/BDobjects/'
#plot intensity vs observed proportion damaged for all BD objects in a file

plot_bd_df_int <- data.frame(event_name = character(), intensity=numeric(), obs_prop_bd=numeric(), tot_obs=integer())

ufiles_BD <- list.files(path=input_folder,pattern=Model$haz,recursive = T,ignore.case = T)
build_all <- data.frame(event_name=character(), grading=character())
for (i in 1:length(ufiles_BD)){
  file <- ufiles_BD[i]
  BDy <- readBD(paste0(input_folder, file))
  if (is.null(BDy) | NROW(BDy)==0){next}
  BDy_df <- as.data.frame(BDy)
  build_all %<>% rbind(data.frame(event_name=file, grading=BDy_df$grading))
  hazard<-rep(NA_real_,length(BDy_df$hazMean1))
  BDy_df$hazard <- apply(as.data.frame(BDy_df[,grepl("Mean",names(BDy_df))]), 1, max, na.rm=T)
  # for (variable in names(BDy)[grepl("Mean",names(BDy))]){
  #   tmp<-BDy[variable]
  #   tmp$hazard<-hazard
  #   hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
  # }
  # BDy$hazard <- hazard
  BDy_df$intensity <- round(BDy_df$hazard*2)/2

  plot_bd_df_int %<>% add_row((BDy_df %>% group_by(intensity) %>% dplyr::summarise(event_name=file,
                                                                                    obs_prop_bd = mean(grading !='notaffected'),
                                                                                    tot_obs=n())))
  
}

black_events <- c(
  "2018-02-26, PNG", "2018-07-29, IDN", "2018-10-06, HTI", 
  "2019-08-08, TUR", "2019-09-24, PAK", "2019-09-26, IDN", "2019-11-08, IRN", 
  "2019-11-26, ALB", "2020-01-24, TUR", "2020-10-30, TUR", "2020-12-29, HRV", 
  "2021-01-14, IDN", "2021-08-14, HTI"
)

# Create a color palette for the events
colour_events <- c("2013-09-24, PAK", "2015-04-25, NPL", "2016-04-16, ECU", 
                   "2016-08-24, ITA", "2016-10-26, ITA",  "2017-09-19, MEX", "2017-09-07, MEX", 
                   "2017-11-12, IRN", "2022-06-22, AFG", "2022-07-02, IRN", 
                   "2022-07-27, PHL", "2022-11-21, IDN", "2023-02-06, TUR"
)

soft_color_palette <- colorRampPalette(colors = c("#F0E442", "#E69F00", "#D55E00", "#9c1c13",
                                                  "#6a3499", "#CC79A7", "#56B4E9", "#0072B2", 
                                                  "#009E73", "#66CC99", "#B3B300", "#999999"),
                                       space = "Lab")(length(colour_events))

# Assign the generated soft color palette to the events, keeping black events black
color_palette <- setNames(soft_color_palette, colour_events )

for(event in black_events) {
  color_palette[event] <- "black"
}

plot_bd_df_int$event_name_fct <- extract_country_date(plot_bd_df_int$event_name)
plot_bd_df_int$event_name_fct <- factor(plot_bd_df_int$event_name_fct, 
                                        levels=c("2013-09-24, PAK", "2015-04-25, NPL", "2016-04-16, ECU", 
                                                 "2016-08-24, ITA", "2016-10-26, ITA", "2017-09-19, MEX",
                                                 "2017-09-07, MEX", "2017-11-12, IRN", "2022-06-22, AFG", 
                                                 "2022-07-02, IRN", "2022-07-27, PHL", "2022-11-21, IDN", 
                                                 "2023-02-06, TUR",
                                                 black_events))


p7 <- ggplot(plot_bd_df_int %>% filter(!is.na(plot_bd_df_int$obs_prop_bd)), 
             aes(x = intensity, y = obs_prop_bd, group = event_name_fct, color = event_name_fct)) +
  #geom_point() + 
  geom_line(lwd=1) +
  #geom_ribbon(aes(ymin = obs_prop_bd - 3 * std_dev, ymax = obs_prop_bd + 3 * std_dev, fill = event_name), alpha=0.1, color=NA) +
  labs(#title = "Line with Shaded Region for 3 Standard Deviations",
    x = "Exposed Hazard Intensity",
    y = "Observed Proportion Buildings Damaged") + ggtitle('') + ylim(0,1) +
  scale_color_manual(values=color_palette) + 
  labs(color = "Event") + 
  guides(color = guide_legend(override.aes = list(lwd = 1))) + 
  ggtitle('(a) All data') +
  theme_minimal() +
  theme(plot.title=element_text(hjust=-.07, family = "Liberation Serif"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        #axis.text.x = element_text(angle = 45, hjust = 1, family = "Liberation Serif", size=10),
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=12)) 


# Figure: Event vs Prop Damage MAR
input_folder <- 'IIDIPUS_Input_Alternatives/Aug24/BDobjects/'
#plot intensity vs observed proportion damaged for all BD objects in a file

plot_bd_df_int <- data.frame(event_name = character(), intensity=numeric(), obs_prop_bd=numeric(), tot_obs=integer())

ufiles_BD <- list.files(path=input_folder,pattern=Model$haz,recursive = T,ignore.case = T)
build_all <- data.frame(event_name=character(), grading=character())
save_BDMAR <- F
save_folder <- paste0('IIDIPUS_Input_Alternatives/Aug24/BDobjects_MAR/')
for (i in 1:length(ufiles_BD)){
  file <- ufiles_BD[i]
  BDy <- readBD(paste0(input_folder, file))
  if(all(BDy$grading=='possible')) next
  BDy_df <- addBDsource(BDy)
  BDy_df$index <- 1:NROW(BDy_df)
  
  BDy_df$hazard <- apply(as.data.frame(BDy_df[,grepl("hazMean",names(BDy_df))]), 1, max, na.rm=T)
  # for (variable in names(BDy)[grepl("Mean",names(BDy))]){
  #   tmp<-BDy[variable]
  #   tmp$hazard<-hazard
  #   hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
  # }
  # BDy$hazard <- hazard
  BDy_df$rounded_intensity <- round(BDy_df$hazard*2)/2
  
  BDy_df %<>% filter(grading !='possible')
  
  if (!is.null(BDy_df$SensorID)){
    grouping_vars <- c('SensorID')
  } else if (is.null(BDy_df$or_source_nam)){
    grouping_vars <- c('source_nam', 'src_date', 'src_info')
  } else {
    grouping_vars <- c('or_source_nam', 'or_src_date', 'dmg_source_nam', 'dmg_src_date')
  }
  groupedBDy_df <- BDy_df %>%
    group_by(!!!syms(grouping_vars)) %>%
    mutate(group_id = cur_group_id(), MAR=F) %>% group_split()
  
  for (i in 1:length(groupedBDy_df)){
    if (NROW(groupedBDy_df[[i]]) < 5) next
    prop_unaff_low_int <- mean(groupedBDy_df[[i]]$grading[groupedBDy_df[[i]]$hazard <= 7] == 'notaffected')
    prop_unaff_high_int <- mean(groupedBDy_df[[i]]$grading[groupedBDy_df[[i]]$hazard > 7] == 'notaffected')
    if (!is.na(prop_unaff_low_int)){
      if (prop_unaff_low_int<0.05) next
    }
    if (!is.na(prop_unaff_high_int)){
      if (prop_unaff_high_int<0.05) next
    }
    groupedBDy_df[[i]]$MAR=T
  }
  
  groupedBDy_df %<>% bind_rows()
  
  if(save_BDMAR){
    BDy <- BDy[groupedBDy_df$index[groupedBDy_df$MAR],]
    if(NROW(BDy)==0) next
    saveBD(BDy, paste0(dir, save_folder, file))
  }
  
  MAR_groupedBDy_df <- groupedBDy_df[groupedBDy_df$MAR,]
  
  # for (variable in names(BDy)[grepl("Mean",names(BDy))]){
  #   tmp<-BDy[variable]
  #   tmp$hazard<-hazard
  #   hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
  # }
  # BDy$hazard <- hazard
  MAR_groupedBDy_df$intensity <- MAR_groupedBDy_df$rounded_intensity
  plot_bd_df_int %<>% add_row((MAR_groupedBDy_df %>% group_by(intensity) %>% dplyr::summarise(event_name=file,
                                                                                       obs_prop_bd = mean(grading !='notaffected'),
                                                                                       tot_obs=n())))
  
}

BDy_HTI2021 <- readBD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Aug24/BDobjects_MAR/EQ20210814HTI_164')
BDy_df <- as.data.frame(BDy_HTI2021)
BDy_df %<>% filter(grading != 'possible')
BDy_df$intensity <- round(apply(BDy_df[,grep('hazMean', names(BDy_df))],1,max, na.rm=T)*2)/2
plot_bd_df_int %<>% add_row((BDy_df %>% group_by(intensity) %>% dplyr::summarise(event_name='EQ20210814HTI_164',
                                                                          obs_prop_bd = mean(grading !='notaffected'),
                                                                          tot_obs=n())))


color_palette['2021-08-14, HTI'] = '#f21818'

plot_bd_df_int$event_name_fct <- extract_country_date(plot_bd_df_int$event_name)
plot_bd_df_int$event_name_fct <- factor(plot_bd_df_int$event_name_fct, 
                                        levels=c("2013-09-24, PAK", "2015-04-25, NPL", "2016-04-16, ECU", 
                                                 "2016-08-24, ITA", "2016-10-26, ITA", "2017-09-19, MEX",
                                                 "2017-09-07, MEX", "2017-11-12, IRN", "2022-06-22, AFG", 
                                                 "2022-07-02, IRN", "2022-07-27, PHL", "2022-11-21, IDN", 
                                                 "2023-02-06, TUR",
                                                 black_events))

plot_bd_df_int_filtered <- plot_bd_df_int %>% filter(!is.na(plot_bd_df_int$obs_prop_bd))
linetype_palette <- rep("solid", length(unique(plot_bd_df_int_filtered$event_name_fct)))  # Start with all solid
linetype_palette[which(unique(plot_bd_df_int_filtered$event_name_fct) == "2021-08-14, HTI")] <- "dashed"  # Set the specific event to dashed

# Create a ggplot with shaded region
p8 <- ggplot(plot_bd_df_int_filtered, 
             aes(x = intensity, y = obs_prop_bd, group = event_name_fct, color = event_name_fct, 
                 linetype= event_name_fct)) +
  #geom_point() +
  geom_line(lwd=1) +
  #geom_ribbon(aes(ymin = obs_prop_bd - 3 * std_dev, ymax = obs_prop_bd + 3 * std_dev, fill = event_name), alpha=0.1, color=NA) +
  labs(#title = "Line with Shaded Region for 3 Standard Deviations",
    x = "Exposed Hazard Intensity",
    y = "Observed Proportion Buildings Damaged") + ggtitle('') + ylim(0,1) +
  scale_color_manual(values=color_palette) + 
  #scale_linetype_manual(values = linetype_palette) + 
  labs(color = "Event") + 
  scale_linetype_manual(values = linetype_palette) +
  guides(
    color = guide_legend("Event", override.aes = list(linetype=linetype_palette, lwd = 1)),  # Single legend with "Event" label
    linetype = "none"  # Remove the separate linetype legend
  ) +
  ggtitle('(b) Filtered Data') + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        plot.title=element_text(hjust=-.07, family = "Liberation Serif"),
        #axis.text.x = element_text(angle = 45, hjust = 1, family = "Liberation Serif", size=10),
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=12),
        legend.key.width=unit(2,"lines"))

p8
plot_grid( p7, p8, align = 'vh', nrow = 2)

#PointDatIntensity.pdf, 10 x 8

#------------------------------------------------------------------------------------------
#----------------------------------------Figure J------------------------------------------
#-------------------Precision-Recall curve for Simulated and Real data---------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

# Plot performance of fitted model
input_folder <- 'IIDIPUS_Input_Alternatives/Aug24/BDobjects_MAR/' #'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput4/BDobjects/'#
ufiles_BD <- list.files(path=input_folder,pattern=Model$haz,recursive = T,ignore.case = T)
#plot intensity vs observed proportion damaged for all BD objects in a file

nSamps <- 50
bd_df_int_real <- data.frame(event_name = character(), 
                             intensity=numeric(), 
                             obs_prop_bd=numeric(), 
                             tot_obs=integer())
for (i in 1:nSamps) {
  col_name <- paste0("pDamSamp.", i)
  bd_df_int_real[[col_name]] <- numeric()
}

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2024-10-31_164512_MCMC_RealAgg5_LR40_Rho0.9_15v0_adaptive_noHLP')
#AlgoResults$s_finish <- 160

set.seed(1)
for (i in 1:length(ufiles_BD)){
  file <- ufiles_BD[i]
  if (file == "EQ20210814HTI_164") next
  BDy <- readBD(paste0(input_folder, file))
  
  Omega <- relist(AlgoResults$Omega_sample_phys[,sample(950:1250,1)], skeleton=Model$skeleton)
  pDamSamp <- BDX(BDy,Omega %>% addTransfParams(),Model,Method = list(Np = 1, cores=2), output='results_analysis')
  pDamSamp$intensity <- apply(pDamSamp[,grep('hazMean', names(pDamSamp))], 1, max, na.rm=T)
  pDamSamp <- pDamSamp[, !names(pDamSamp) %in% c('LP', 'spatial_pixel', 'spatial_pixel_id', grep('hazMean', names(pDamSamp), value = TRUE))]
  pDamSamp$event_name = file
  
  for (j in 2:nSamps){
    Omega <- relist(AlgoResults$Omega_sample_phys[,sample(950:1250,1)], skeleton=Model$skeleton)
    pDamSamp_j <- BDX(BDy,Omega %>% addTransfParams(),Model,Method = list(Np = 1, cores=2), output='results_analysis')
    pDamSamp %<>% add_column(!!paste0('pDamSamp.', j) := pDamSamp_j$pDamSamp.1)
  }
  bd_df_int_real %<>% rbind(pDamSamp)
  
  # for (variable in names(BDy)[grepl("Mean",names(BDy))]){
  #   tmp<-BDy[variable]
  #   tmp$hazard<-hazard
  #   hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
  # }
  # BDy$hazard <- hazard
  
}

bd_df_int_real$rounded_intensity <- round(bd_df_int_real$intensity*2)/2
plot_bd_df_int_real = bd_df_int_real %>% group_by(event_name, rounded_intensity) %>% 
  dplyr::summarise(upd_mean=sum(tot_obs*obs_prop_bd)/sum(tot_obs),
            intensity=unique(rounded_intensity),
            across(starts_with("pDamSamp"), ~ sum(.x * tot_obs) / sum(tot_obs), .names = "{col}")
  )

#plot_bd_df_int_real$rounded_intensity <- round(  plot_bd_df_int_real$intensity)/2

black_events <- c(
  "2018-02-26, PNG", "2018-07-29, IDN", "2018-10-06, HTI", 
  "2019-08-08, TUR", "2019-09-24, PAK", "2019-09-26, IDN", "2019-11-08, IRN", 
  "2019-11-26, ALB", "2020-01-24, TUR", "2020-10-30, TUR", "2020-12-29, HRV", 
  "2021-01-14, IDN", "2021-08-14, HTI"
)

# Create a color palette for the events
colour_events <- c("2013-09-24, PAK", "2015-04-25, NPL", "2016-04-16, ECU", 
                   "2016-08-24, ITA", "2016-10-26, ITA",  "2017-09-19, MEX", "2017-09-07, MEX", 
                   "2017-11-12, IRN", "2022-06-22, AFG", "2022-07-02, IRN", 
                   "2022-07-27, PHL", "2022-11-21, IDN", "2023-02-06, TUR"
)

soft_color_palette <- color_palette <- colorRampPalette(colors = c("#F0E442", "#E69F00", "#D55E00", "#9c1c13",
                                                  "#6a3499", "#CC79A7", "#56B4E9", "#0072B2", 
                                                  "#009E73", "#66CC99", "#B3B300", "#999999"),
                                       space = "Lab")(length(colour_events))

plot_bd_df_int_real$event_name_fct <- extract_country_date(plot_bd_df_int_real$event_name)
plot_bd_df_int_real$event_name_fct <- factor(plot_bd_df_int_real$event_name_fct, 
                                             levels=c("2013-09-24, PAK", "2015-04-25, NPL", "2016-04-16, ECU", 
                                                      "2016-08-24, ITA", "2016-10-26, ITA", "2017-09-19, MEX",
                                                      "2017-09-07, MEX", "2017-11-12, IRN", "2022-06-22, AFG", 
                                                      "2022-07-02, IRN", "2022-07-27, PHL", "2022-11-21, IDN", 
                                                      "2023-02-06, TUR",
                                                      black_events))

event_names <- names(color_palette)
color_palette <- setNames(soft_color_palette, colour_events )
for(event in black_events) {
  color_palette[event] <- "black"
}

plot_bd_df_int_real$pDamSamp_median <- apply(plot_bd_df_int_real[,grep('pDamSamp.', names(plot_bd_df_int_real))],1 , median)
plot_bd_df_int_real$pDamSamp_q5 <- apply(plot_bd_df_int_real[,grep('pDamSamp.', names(plot_bd_df_int_real))],1 , quantile,0.05)
plot_bd_df_int_real$pDamSamp_q95 <- apply(plot_bd_df_int_real[,grep('pDamSamp.', names(plot_bd_df_int_real))],1 , quantile,0.95)

linetype_palette <- setNames(c('solid', 'dashed'), 
                             c('Observed Damage Proportion', 'Posterior Probability of Damage'))
# Create a ggplot with shaded region
ggplot(plot_bd_df_int_real[!is.na(plot_bd_df_int_real$upd_mean),], 
       aes(x = intensity, y = upd_mean, group = event_name_fct, color = event_name_fct, linetype='Observed Damage Proportion')) +
  geom_line(lwd=1) +
  geom_line(aes(x=intensity, y=pDamSamp_median, group=event_name_fct, color=event_name_fct, linetype='Posterior Probability of Damage')) +
  #geom_ribbon(aes(ymin = obs_prop_bd - 3 * std_dev, ymax = obs_prop_bd + 3 * std_dev, fill = event_name), alpha=0.1, color=NA) +
  labs(#title = "Line with Shaded Region for 3 Standard Deviations",
    x = "Exposed Hazard Intensity",
    y = "Observed Proportion Buildings Damaged") + ggtitle('') + ylim(0,1) +
  scale_color_manual(values=color_palette) + 
  scale_linetype_manual(values=linetype_palette) + 
  labs(color = "Event") + 
  guides(color = guide_legend(override.aes = list(lwd = 1)),
         linetype='none') +
  #ggtitle('(b) Filtered Data') + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        plot.title=element_text(hjust=-.07, family = "Liberation Serif"), 
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=12))

#PostMedianProbs.pdf, 6 x 10

plot_bd_df_int_real2 = plot_bd_df_int_real %>% filter(event_name != "EQ20210814HTI_164")
plot(plot_bd_df_int_real2$pDamSamp_median, plot_bd_df_int_real2$upd_mean)
abline(a=0, b=1)
#plot_grid( plot_grid( p7, p8, align = 'vh', nrow = 2),legend, ncol = 1, rel_heights=c(1,0.1))

bd_df_int_real$pDamMedian <- apply(bd_df_int_real[,grep('pDamSamp.', names(bd_df_int_real))], 1, median)

ROC_calc <- bd_df_int_real[, c('tot_obs', 'obs_prop_bd', 'pDamMedian', 'intensity')]
ROC_calc$obs_bd <- ROC_calc$tot_obs * ROC_calc$obs_prop_bd
#ROC_calc$samp_bd <- ROC_calc$tot_obs * ROC_calc$pDamMedian

# ROC_scores <- data.frame(threshold=numeric(),
#                          TP=numeric(), 
#                          FP=numeric())

# for(p_mult in seq(-1,1,0.01)){
#   ROC_calc$p_samp_bd <- pmax(0,pmin(1,ROC_calc$pDamMedian + p_mult))
#   TPR = sum(ROC_calc$p_samp_bd * ROC_calc$obs_bd) / sum(ROC_calc$obs_bd)
#   FPR = sum(ROC_calc$p_samp_bd * (ROC_calc$tot_obs - ROC_calc$obs_bd)) / sum(ROC_calc$tot_obs - ROC_calc$obs_bd)
#   ROC_scores %<>% add_row(
#     p_mult=p_mult,
#     TP = TPR,
#     FP = FPR
#   )
# }

# for(threshold in seq(0,1,0.0001)){
#   ROC_calc$p_samp_bd <- ifelse(ROC_calc$pDamMedian >threshold, 1, 0)
#   TPR = sum(ROC_calc$p_samp_bd * ROC_calc$obs_bd) / sum(ROC_calc$obs_bd)
#   FPR = sum(ROC_calc$p_samp_bd * (ROC_calc$tot_obs - ROC_calc$obs_bd)) / sum(ROC_calc$tot_obs - ROC_calc$obs_bd)
#   ROC_scores %<>% add_row(
#     threshold=threshold,
#     TP = TPR,
#     FP = FPR
#   )
# }

# Prec_Recall <- data.frame(threshold=numeric(),
#                           Precision=numeric(), 
#                           Recall=numeric())
# 
# exp(seq(-10, 0, 0.1))
# for(threshold in exp(seq(-20, 5, 0.1))){
#   ROC_calc$p_samp_bd <- pmax(0,pmin(ROC_calc$pDamMedian/threshold, 1))#ifelse(ROC_calc$pDamMedian >threshold, 1, 0)
#   TP = sum(ROC_calc$p_samp_bd * ROC_calc$obs_bd)
#   FP = sum(ROC_calc$p_samp_bd * (ROC_calc$tot_obs-ROC_calc$obs_bd))
#   FN = sum((1-ROC_calc$p_samp_bd) * ROC_calc$obs_bd)
#   Precision = TP / (TP + FP)
#   Recall = TP / (TP + FN)
#   Prec_Recall %<>% add_row(
#     threshold=threshold,
#     Precision = Precision,
#     Recall = Recall
#   )
# }
# lines(Prec_Recall$Recall, Prec_Recall$Precision, type='l')
# abline(h=sum(ROC_calc$obs_prop_bd * ROC_calc$tot_obs)/sum(ROC_calc$tot_obs),lty=2)


#ROC_calc$test_1 <- rbern(length(ROC_calc$obs_prop_bd),ROC_calc$obs_prop_bd)
#ROC_calc$pDamMedian <- runif(length(ROC_calc$pDamMedian))

Prec_Recall <- data.frame(threshold=numeric(),
                          Precision=numeric(), 
                          Recall=numeric(), 
                          FPR = numeric(),
                          TPR = numeric())

for(threshold in seq(0,1,0.0001)){
  # ROC_calc$p_samp_bd <- ifelse(ROC_calc$pDamMedian >threshold, 1, 0)
  # TP = sum(ROC_calc$p_samp_bd * ROC_calc$obs_bd)
  # FP = sum(ROC_calc$p_samp_bd * (ROC_calc$tot_obs-ROC_calc$obs_bd))
  # FN = sum((1-ROC_calc$p_samp_bd) * ROC_calc$obs_bd)
  ROC_calc$p_samp_bd <- ifelse(ROC_calc$pDamMedian >threshold, 1, 0)
  TP = sum(ROC_calc$p_samp_bd * ROC_calc$obs_bd)
  FP = sum(ROC_calc$p_samp_bd * (ROC_calc$tot_obs-ROC_calc$obs_bd))
  FN = sum((1-ROC_calc$p_samp_bd) * ROC_calc$obs_bd)
  TN = sum((1-ROC_calc$p_samp_bd) * (ROC_calc$tot_obs-ROC_calc$obs_bd))
  Precision = TP / (TP + FP)
  Recall = TP / (TP + FN)
  Prec_Recall %<>% add_row(
    threshold=threshold,
    Precision = Precision,
    Recall = Recall,
    TPR = TP/(TP + FN),
    FPR = FP/(FP+TN),
  )
}


# SIM DATA: 
input_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput4/BDobjects/'#
ufiles_BD <- list.files(path=input_folder,pattern=Model$haz,recursive = T,ignore.case = T)
#plot intensity vs observed proportion damaged for all BD objects in a file

nSamps <- 50
bd_df_int_sim <- data.frame(event_name = character(), 
                            intensity=numeric(), 
                            obs_prop_bd=numeric(), 
                            tot_obs=integer())
for (i in 1:nSamps) {
  col_name <- paste0("pDamSamp.", i)
  bd_df_int_sim[[col_name]] <- numeric()
}


Omega <- Omega_true <- list(Lambda1 = list(nu=8.75, kappa=0.6),
                            Lambda2 = list(nu=11.7, kappa=0.75), #list(nu=10.65, kappa=1.5), #
                            Lambda3 = list(nu=9.55, kappa=0.68),
                            Lambda4 = list(nu=9.9, kappa=1.6),
                            theta= list(theta1=0.6),
                            eps=list(local=0.8, hazard_mort=0.45, hazard_disp=0.6, hazard_bd=0.5, hazard_cor=0.55),
                            #eps = list(local=1.3, hazard_mort=0.8383464, hazard_disp=1, hazard_bd=0.9, hazard_cor=0.55),
                            vuln_coeff = list(PDens=0, SHDI=-0.08, GNIc=-0.02, Vs30=0.01, EQFreq=-0.02, FirstHaz=0.01, Night=0, FirstHaz.Night=0.05),
                            check = list(check=0.5))

for (i in 1:length(ufiles_BD)){
  file <- ufiles_BD[i]
  BDy <- readBD(paste0(input_folder, file))
  
  pDamSamp <- BDX(BDy,Omega %>% addTransfParams(),Model,Method = list(Np = nSamps, cores=2), output='results_analysis')
  pDamSamp$intensity <- apply(pDamSamp[,grep('hazMean', names(pDamSamp))], 1, max, na.rm=T)
  pDamSamp <- pDamSamp[, !names(pDamSamp) %in% c('LP', 'spatial_pixel', 'spatial_pixel_id', grep('hazMean', names(pDamSamp), value = TRUE))]
  pDamSamp$event_name = file
  
  bd_df_int_sim %<>% rbind(pDamSamp)
  
  
  # for (variable in names(BDy)[grepl("Mean",names(BDy))]){
  #   tmp<-BDy[variable]
  #   tmp$hazard<-hazard
  #   hazard<-apply(tmp@data,1,function(x) max(x,na.rm=T))
  # }
  # BDy$hazard <- hazard
  
}

bd_df_int_sim$rounded_intensity <- round(bd_df_int_sim$intensity*2)/2
plot_bd_df_int_sim = bd_df_int_sim %>% group_by(event_name, rounded_intensity) %>% 
  dplyr::summarise(upd_mean=sum(tot_obs*obs_prop_bd)/sum(tot_obs),
            intensity=unique(rounded_intensity),
            across(starts_with("pDamSamp"), ~ sum(.x * tot_obs) / sum(tot_obs), .names = "{col}")
  )

bd_df_int_sim$pDamMedian <- apply(bd_df_int_sim[,grep('pDamSamp.', names(bd_df_int_sim))], 1, median)


ROC_calc_sim <- bd_df_int_sim[, c('tot_obs', 'obs_prop_bd', 'pDamMedian', 'intensity')]
ROC_calc_sim$obs_bd <- ROC_calc_sim$tot_obs * ROC_calc_sim$obs_prop_bd


Prec_Recall2 <- data.frame(threshold=numeric(),
                           Precision=numeric(), 
                           Recall=numeric(), 
                           TPR=numeric(),
                           FPR=numeric())
for(threshold in exp(seq(-80, 0, 0.02))){
  # ROC_calc$p_samp_bd <- ifelse(ROC_calc$pDamMedian >threshold, 1, 0)
  # TP = sum(ROC_calc$p_samp_bd * ROC_calc$obs_bd)
  # FP = sum(ROC_calc$p_samp_bd * (ROC_calc$tot_obs-ROC_calc$obs_bd))
  # FN = sum((1-ROC_calc$p_samp_bd) * ROC_calc$obs_bd)
  ROC_calc_sim$p_samp_bd <- ifelse(ROC_calc_sim$pDamMedian >threshold, 1, 0)
  TP = sum(ROC_calc_sim$p_samp_bd * ROC_calc_sim$obs_bd)
  FP = sum(ROC_calc_sim$p_samp_bd * (ROC_calc_sim$tot_obs-ROC_calc_sim$obs_bd))
  FN = sum((1-ROC_calc_sim$p_samp_bd) * ROC_calc_sim$obs_bd)
  TN = sum((1-ROC_calc_sim$p_samp_bd) * (ROC_calc_sim$tot_obs-ROC_calc_sim$obs_bd))
  Precision = TP / (TP + FP)
  Recall = TP / (TP + FN)
  Prec_Recall2 %<>% add_row(
    threshold=threshold,
    Precision = Precision,
    Recall = Recall,
    TPR = TP/(TP + FN),
    FPR = FP/(FP+TN),
  )
}

plot(Prec_Recall$Recall, Prec_Recall$Precision, type='l', col='black', xlab='Recall', ylab='Precision')
lines(Prec_Recall2$Recall, Prec_Recall2$Precision, type='l', col='blue', xlab='Recall', ylab='Precision')
abline(h=sum(ROC_calc$obs_prop_bd * ROC_calc$tot_obs)/sum(ROC_calc$tot_obs),lty=2)
abline(h=sum(ROC_calc_sim$obs_prop_bd * ROC_calc_sim$tot_obs)/sum(ROC_calc_sim$tot_obs),col='blue',lty=2)


ggplot() + 
  geom_line(data=Prec_Recall, aes(x=FPR, y=TPR, col='Real Data')) +
  geom_line(data=Prec_Recall2, aes(x=FPR, y=TPR, col='Simulated Data')) + 
  theme_minimal() + xlab('False Positive Rate') + ylab('True Positive Rate') + 
  scale_color_manual(values=setNames(c('black', 'blue'), c('Real Data', 'Simulated Data')),
                     name = NULL  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        #axis.text.x = element_text(angle = 45, hjust = 1, family = "Liberation Serif", size=10),
        axis.title = element_text(family = "Liberation Serif", size=12),  
        legend.text = element_text(family = "Liberation Serif", size=11),    # Legend text
        legend.title = element_text(family = "Liberation Serif", size=12)) +
  geom_abline(slope=1,intercept=0, col='grey', linetype='dashed')

plot(Prec_Recall$FPR, Prec_Recall$TPR, type='l', col='black', xlab='False Positive Rate', ylab='True Positive Rate')
lines(Prec_Recall2$FPR, Prec_Recall2$TPR, type='l',col='blue', xlab='False Positive Rate', ylab='True Positive Rate')
abline(a=0, b=1, col='red')

#ROCcurves.pdf, 4 x 6 inches

# 
# plot(ROC_scores$FP, ROC_scores$TP, xlab='False Positive Rate', ylab='True Positive Rate')
# abline(0,1)
# points(ROC_scores$FP[which(ROC_scores$threshold==0.5)], ROC_scores$TP[which(ROC_scores$threshold==0.5)], col='red', pch=19)
# 
# plot(ROC_scores$FP, ROC_scores$TP, xlim=c(0,0.2), ylim=c(0,0.2),xlab='False Positive Rate', ylab='True Positive Rate')
# abline(0,1)


#------------------------------------------------------------------------------------------
#----------------------------------------Figure J------------------------------------------
#------------------------------Compare simulated and real data-----------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------


siminput_folder <- 'IIDIPUS_Input_Alternatives/IIDIPUS_SimInput4/BDobjects/'#
ufiles_BDsim <- list.files(path=siminput_folder,pattern=Model$haz,recursive = T,ignore.case = T)
#plot intensity vs observed proportion damaged for all BD objects in a file

input_folder <- 'IIDIPUS_Input_Alternatives/Aug24/BDobjects_MAR/'#
ufiles_BD <- list.files(path=input_folder,pattern=Model$haz,recursive = T,ignore.case = T)

bd_df_int_sim <- data.frame(event_name = character(), 
                             obs_prop_bd=numeric(), 
                             tot_obs=integer(), 
                             n_pixels=integer(), 
                             max_intensity=numeric(), 
                             n_buildings=integer())

bd_df_int_real <- bd_df_int_sim

for (i in 1:length(ufiles_BDsim)){
  file <- ufiles_BDsim[i]
  BDy <- readBD(paste0(siminput_folder, file))
  
  bd_df_int_sim %<>% rbind(data.frame(event_name=file,
                                       obs_prop_bd = mean(BDy$grading!='notaffected'),
                                       tot_obs = NROW(BDy),
                                       n_pixels = length(unique(BDy$spatial_pixel)),
                                       max_intensity=max(max(BDy[, grep('hazMean', names(BDy))], na.rm=T), na.rm=T),
                                       n_buildings=NROW(BDy)))

}

for (i in 1:length(ufiles_BD)){
  file <- ufiles_BD[i]
  BDy <- readBD(paste0(input_folder, file))
  
  bd_df_int_real%<>% rbind(data.frame(event_name=file,
                                      obs_prop_bd = mean(BDy$grading!='notaffected'),
                                      tot_obs = NROW(BDy),
                                      n_pixels = length(unique(BDy$spatial_pixel)),
                                      max_intensity=max(max(BDy[, grep('hazMean', names(BDy))], na.rm=T), na.rm=T),
                                      n_buildings=NROW(BDy)))
  
}

grid.arrange(ggplot(bd_df_int_real, aes(x=event_name, y=obs_prop_bd)) + geom_point() + ggtitle('Real'),
             ggplot(bd_df_int_sim, aes(x=event_name, y=obs_prop_bd)) + geom_point()+ ggtitle('Sim'),
             ggplot(bd_df_int_real, aes(x=event_name, y=n_buildings)) + geom_point() + ggtitle('Real'),
             ggplot(bd_df_int_sim, aes(x=event_name, y=n_buildings)) + geom_point()+ ggtitle('Sim'))
grid.arrange(ggplot(bd_df_int_real, aes(x=event_name, y=n_pixels)) + geom_point() + ggtitle('Real'),
       ggplot(bd_df_int_sim, aes(x=event_name, y=n_pixels)) + geom_point()+ ggtitle('Sim'))
grid.arrange(ggplot(bd_df_int_real, aes(x=event_name, y=max_intensity)) + geom_point() + ggtitle('Real'),
             ggplot(bd_df_int_sim, aes(x=event_name, y=max_intensity)) + geom_point()+ ggtitle('Sim'))
grid.arrange(ggplot(bd_df_int_real, aes(x=event_name, y=n_buildings)) + geom_point() + ggtitle('Real'),
             ggplot(bd_df_int_sim, aes(x=event_name, y=n_buildings)) + geom_point()+ ggtitle('Sim'))

mean(bd_df_int_real$obs_prop_bd*bd_df_int_real$n_buildings/sum(bd_df_int_real$n_buildings))
mean(bd_df_int_sim$obs_prop_bd*bd_df_int_sim$n_buildings/sum(bd_df_int_sim$n_buildings))
