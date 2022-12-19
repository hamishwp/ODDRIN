folder<-"./IIDIPUS_Input/HAZARDobjects/TC20200404VUT/"
bbox<-c(164.465332,-20.427013,173.474121,-12.640338)
# bbox<-c(153.65,-25.21,179.8,-2.69)

library(raster)
library(magrittr)
library(pracma)

GetEODISwind<-function(folder,bbox){
  
  e <- as(raster::extent(c(bbox[c(1,3,2,4)])), 'SpatialPolygons')
  crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  
  output<-data.frame()
  # for(i in c(strcat("0",as.character(1:9)),"10")){
  for(i in 3:8){
    
    SPEEDMAX<-brick(paste0(folder,"MERRA2_400.tavg1_2d_flx_Nx.2020040",i,".nc4"),varname="SPEEDMAX")
    SPEEDMAX%<>%raster::crop(e)
    
    tmp<-as.data.frame(getValues(SPEEDMAX))
    colnames(tmp)<-paste0("2020040",i,"_",1:24,"h")
    if(nrow(output)==0) {output<-tmp
    }else output%<>%cbind(tmp)
    
    print(paste0("2020-04-0",i," : ",max(SPEEDMAX[])))
    
  }
  rm(tmp)
  SPEEDMAX%<>%as("SpatialPixelsDataFrame")
  SPEEDMAX@data<-output
  colnames(SPEEDMAX@coords)<-rownames(SPEEDMAX@bbox)<-c("Longitude","Latitude")
  
  # library(ggmap)
  
  # mad_map <- get_stamenmap(SPEEDMAX@bbox,source = "stamen",maptype = "terrain",zoom=8)
  # q<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")
  
  # for(i in 1:ncol(SPEEDMAX@data)){
  #   
  #   tmp<-as.data.frame(SPEEDMAX[i])
  #   names(tmp)[1]<-"z"
  #   namer<-strsplit(names(SPEEDMAX[i]),split = "_")[[1]]
  #   
  #   p<-ggmap(mad_map,base_layer = ggplot(tmp,aes(x,y,z=z)))+
  #     geom_raster(aes(fill=z),alpha=0.7,interpolate = T) + coord_cartesian() +
  #   scale_fill_gradient2(low = "transparent",high = "red",
  #                        na.value = "transparent",limits=c(0,42)) +
  #     labs(fill = "Hourly Max. Wind [m/s]")+xlab("Longitude") + ylab("Latitude")+
  #     ggtitle(paste0(namer[2]," ",as.Date(namer[1],format = "%Y%m%d"))) +
  #     theme(plot.title = element_text(hjust = 0.5))
  #   
  #   ggsave(paste0(formatC(i, width = 4, format = "d", flag = "0"),"_TC_VUT.png"), plot=p,path = "./Plots/VUT/",width = 7,height = 5.)
  # }
  
  sdate<-"2020-04-03"
  fdate<-"2020-04-08"
  lenny<-ncol(SPEEDMAX)
  I0<-3
  lhazSDF<-c(list(bbox=bbox,sdate=sdate,fdate=fdate,NumEvents=lenny,hazard="TC",I0=I0))
  
  for(i in 1:ncol(SPEEDMAX)){
    # Create HAZARD object
    hazsdf<-SPEEDMAX[i]
    class(hazsdf)<-"HAZARD"
    hazsdf@hazard<-"TC"
    hazsdf@I0<-3
    hazsdf@eventdate<-as.Date(strsplit(names(hazsdf),split = "_")[[1]][1],format = "%Y%m%d")
    hazsdf@alertlevel<-NA_character_
    hazsdf@alertscore<-NA_real_
    names(hazsdf)<-"mean"
    ### REMEMBER TO TAKE THE LOG VALUE FOR TCs ###
    hazsdf@data[,1]%<>%log()
    # Add to the list of hazards
    lhazSDF%<>%c(hazsdf)
  }
  
  return(lhazSDF)
  
}

# bbox<-c(164.465332,-20.427013,173.474121,-12.640338)
# lhazSDF<-GetEODISwind(folder="./IIDIPUS_Input/HAZARDobjects/TC20200404VUT/",bbox=bbox)
# ODDy@data<-dplyr::select(ODDy@data,c("Population","GDP","ISO3C"))
# ODDy%<>%AddHazSDF(lhazSDF)


