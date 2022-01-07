library(ggplot2)
library(sf)
library("ggmap")
library(OpenStreetMap)
library(osmdata)
library(tidyverse)
library(geojsonR)
library(raster)
library(dplyr)
library(magrittr)

source('RCode/HAZARDobj.R')

ExtractUSGS<-function(url,namer,I0=NULL,plotty=F){
  
  # temp <- tempfile()
  temp<-paste0(namer,".zip")
  download.file(url,temp)
  unzip(paste0(temp),exdir = paste0(namer,"/"))
  # unzip(temp, exdir = namer)
  # a <- raster("~/Downloads/1241508_4_raster/mmi_mean.flt")
  meanhaz<-raster(file.path(namer,"mmi_mean.flt"))
  sdhaz<-raster(file.path(namer,"mmi_std.flt"))
  unlink(temp)
  
  sgdf <- as(meanhaz, 'SpatialPixelsDataFrame') ; rm(meanhaz)
  tmp<- as(sdhaz, 'SpatialPixelsDataFrame') ; rm(sdhaz)
  sgdf$mmi_std<-tmp$mmi_std ; rm(tmp)
  crs(sgdf)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
  
  colnames(sgdf@coords)<-rownames(sgdf@bbox)<-c("Longitude","Latitude")
  
  if(!is.null(I0)) sgdf<-sgdf[sgdf$mmi_mean>I0,]
  
  if(plotty){
    titlz<-strsplit(namer,"/")[[1]] ; titlz<-titlz[length(titlz)]
    
    mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "toner-lite",zoom=8)
    sdf<-as.data.frame(sgdf); colnames(sdf)<-c("Intensity","Longitude","Latitude")
    ggmap(mad_map) + geom_tile(data=sdf,mapping=aes(x=Longitude, y=Latitude,fill=Intensity),alpha=0.5) + coord_equal() +
      scale_fill_gradient(low = "yellow", high="red", limits=c(4.5,9), na.value = "transparent") +
      theme_bw() + xlab("Longitude") + ylab("Latitude") + ggtitle(titlz) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(sgdf)
  
}

SearchUSGSbbox<-function(bbox,sdate,fdate=NULL,minmag=5){
  
  if(is.null(fdate)) {
    fdate=min(Sys.Date(),(as.Date(sdate)+10))}
  else fdate=min(Sys.Date(),as.Date(fdate)+1)
  
  debut<-"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
  FROM_GeoJson(paste0(debut,"&starttime=",as.Date(sdate)-1,"&endtime=",as.Date(fdate),
                            "&minlongitude=",bbox[1],"&minlatitude=",bbox[2],
                            "&maxlongitude=",bbox[3],"&maxlatitude=",bbox[4],
                            "&minmagnitude=",minmag,"&orderby=magnitude",
                            "&producttype=shakemap"))%>%return
  
}

GetUSGS<-function(bbox,sdate,fdate=NULL,titlz="tmp",I0=4.5,minmag=5,dfGDACS=NULL){

  # Search through the USGS events
  USGS<-SearchUSGSbbox(bbox,sdate,fdate,minmag)
  lenny<-length(USGS$features)
  # if(lenny==0) USGS<-SearchUSGSbbox(bbox,sdate-2,fdate+2,minmag)
  # lenny<-length(USGS$features)
  if(lenny==0) return(NULL)
  
  # Find which countries the hazard occurred in
  pcoords<-array(NA,dim=c(lenny,2))
  for(l in 1:lenny) pcoords[l,]<-USGS$features[[l]]$geometry$coordinates
  iso3<-unique(coords2country(pcoords))
  
  # Get GDACS alertscore for benchmarking
  alertscores<-GetGDACSalertscore(haz = "EQ",bbox = bbox,sdater = sdate,
                                  fdater = fdate,isos = iso3,dfGDACS = dfGDACS)

  lhazdat<-list(hazard_info=list(bbox=bbox,sdate=sdate,fdate=fdate,NumEvents=lenny,hazard="EQ",I0=I0))
  tbbox<-rep(NA,4)
  for (i in 1:lenny){
    # Extract EQ raster image of hazard intensity
    tmp<-FROM_GeoJson(USGS$features[[i]]$properties$detail)
    hazsdf<-tryCatch(ExtractUSGS(url = tmp$properties$products$shakemap[[1]]$contents$`download/raster.zip`$url,
                        namer = paste0(directory,"Disaster_Data/USGS/",titlz,i)),
                     error=function(e) NULL)
    # Extra check (trust me, it's necessary... FML)
    if (is.null(hazsdf)){
      lhazdat$hazard_info$NumEvents<-lhazdat$hazard_info$NumEvents-1
      next
    }
    if (all(hazsdf@bbox[c(1,3)]<bbox[1]) | all(hazsdf@bbox[c(1,3)]>bbox[3]) |
        all(hazsdf@bbox[c(2,4)]<bbox[2]) | all(hazsdf@bbox[c(2,4)]>bbox[4]) |
        max(hazsdf@data$mmi_mean)<minmag) {
          lhazdat$hazard_info$NumEvents<-lhazdat$hazard_info$NumEvents-1
          next
    }
    # Create HAZARD object
    hazsdf<-new("HAZARD",
                obj=hazsdf,
                hazard="EQ",
                dater=sdate,
                I0=I0,
                alertlevel=ifelse(is.null(tmp$properties$alert),"green",tmp$properties$alert),
                alertscore=ifelse(i<=length(alertscores),alertscores[i],0))
    # Add to the list of hazards
    lhazdat[[length(lhazdat)+1]]<-hazsdf
    tbbox[c(1,2)]<-apply(cbind(tbbox[c(1,2)],hazsdf@bbox[c(1,2)]),1,min,na.rm=T)
    tbbox[c(3,4)]<-apply(cbind(tbbox[c(3,4)],hazsdf@bbox[c(3,4)]),1,max,na.rm=T)
  }
  if(any(is.na(tbbox))) return(NULL)
  
  lhazdat$hazard_info$bbox<-tbbox
  return(lhazdat)
  
}

ExtractShakeUSGS<-function(url,I0){
  shakemap<-FROM_GeoJson(url)
  for (i in 1:length(shakemap$features)){
    if(shakemap$features[[i]]$properties$value>=I0) break
  }
  
  if(i==1) stop("ERROR: Shakemap in USGS has all I>I0 hazard intensity")
  
  polylines<-shakemap$features[[(i-1)]]$geometry$coordinates[[1]]
  
  return(list(polylines=polylines,I0=shakemap$features[[(i-1)]]$properties$value))
}
