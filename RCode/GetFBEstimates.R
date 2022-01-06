library(reshape2)
library(tidyverse)
library(gridExtra)
library("ggmap")
library(OpenStreetMap)
library(osmdata)
library(ggplot2)
library(dplyr)
library(SparkR)
options(stringsAsFactors = FALSE)

### TO DO ###
# Output is a list of places, with raster polygon to map total area, 
# the number of people displaced per city in time
# and a list of places/city names to where people are displaced
# Write a function that takes the city names from-to and evaluates population density of both cities and distance & time between.
#############


# name<-'Flooding_India_10-19'
# reference<-"North-Eastern India"

# name<-'Earthquake'
# reference<-"Philippines"

name<-'Wildfires'
name<-'Typhoon'
# directory<-('/home/patten/Documents/IDMC/Facebook_Data/')
# FBdirectory<-('/home/patten/Documents/Coding/Oxford/IIDIPUS/FBdata/VC20210409VCT/')
# name<-c('Wildfires',"Adelaide")
# name<-c('Wildfires',"Eastern_New_South_Australia")
# name<-c('Wildfires',"East_Gippsland")
# name<-c('Wildfires',"South_Coast_NSW")
# reference<-"Australia"
reference<-"La_Soufriere_St_Vincent"
name<-"2331841688749824"

plotty<-TRUE
check=TRUE

GetFBEstimates(FBdirectory, name, check=TRUE, plotty=TRUE)

GetFBEstimates<-function(FBdirectory, name, reference="", check=FALSE, plotty=FALSE){
  
  tmp<-paste(name,collapse = '.*')
  namer=paste(name,collapse = '_')
  dirs <- list.files(path=FBdirectory,pattern=tmp)
  print(dirs)
  # Extract dates
  cdates<-t(as.data.frame(str_split(dirs,"_")))[,2]
  dates <- as.Date(cdates); ldays<-length(cdates)
  # Check that all dates are included in the decay rate
  if (check){
    dd<-as.numeric(format(dates,'%d'))
    mm<-as.numeric(format(dates,'%m'))
    yy<-as.numeric(format(dates,'%Y'))
    for (i in 2:ldays) {
      if (!(dd[i]==dd[i-1]+1)){
        if (!(mm[i]==mm[i-1]+1)){
          if (!(yy[i]==yy[i-1]+1)){stop(paste0("Not all days are included in the data, download the day ",toString(dates[i])," from GeoInsights"))}
        }
      }
    }
    print("All days are included! Data is good to go!")
  }
  # Extract files
  files <- lapply(paste0(FBdirectory,dirs),read.csv)
  names(files)<-cdates
  files<-bind_rows(files, .id = "cdates")
  files$cdates<-as.Date(files$ds); files$ds<-NULL
  
  # files<-files%>%mutate(time=as.integer(cdates-min(cdates)))
  files$day<-files$cdates-min(files$cdates,na.rm = T)
  
  return(files)
  
  #create dataframe with number of displaced, 
  disp<-filter(files,status=="Displaced"|status=="Displaced Abroad")
  dispA<-filter(files,status=="Displaced Abroad")
  notdisp<-filter(files,status=="Never Displaced")
  unknown<-filter(files,status=="Unknown")
  retn<-filter(files,status=="Returned")
  
  if(!plotty){return(list(disp,notdisp,retn,unknown,dispA))}
  
  sumvals<-retn%>%group_by(day)%>%summarise(max=max(n_people,na.rm = TRUE),min=min(n_people,na.rm = TRUE),sum=sum(n_people,na.rm = TRUE))%>%arrange(day)
  p<-ggplot(sumvals,aes(x=day,y=sum,group=gender))+geom_point(aes(colour=gender),size=2)+geom_line(aes(colour=gender))+#scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total Returned IDP Stock");p
  ggsave(paste0(namer,"_",reference,"_Total-Returned.png"), plot=p,path = paste0(directory,'Plots/'),width = 5,height = 4.)
  
  sumvals<-notdisp%>%group_by(day)%>%summarise(max=max(n_people,na.rm = TRUE),min=min(n_people,na.rm = TRUE),sum=sum(n_people,na.rm = TRUE))%>%arrange(day)
  p<-ggplot(sumvals,aes(x=day,y=sum))+geom_point(size=2)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total Never Displaced");p
  ggsave(paste0(namer,"_",reference,"_Total-NeverDisp.png"), plot=p,path = paste0(directory,'Plots/'),width = 5,height = 4.)
  
  sumvals<-disp%>%group_by(name,day,gender)%>%summarise(sum=sum(n_people,na.rm = TRUE))
  p<-ggplot(sumvals,aes(x=day,y=sum,group=gender))+geom_line(aes(colour=gender),size=1.5)+#scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total IDP Stock")+
    facet_wrap(. ~ name,nrow = 2,scales = "fixed") + theme(strip.text.x = element_text(size = 12))
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-IDPs.png"), plot=p,path = paste0(directory,'Plots/'),width = 5,height = 4.)
  
  mxd=sumvals$cdates[which.max(sumvals$sum)]
  mnd=sumvals$cdates[which.min(sumvals$sum[sumvals$cdates>mxd])]
  
  # mnlo<-min(c(min(files$baseline_polygon_centroid_longitude, na.rm = TRUE),min(files$today_polygon_centroid_longitude, na.rm = TRUE)))
  # mxlo<-max(c(max(files$baseline_polygon_centroid_longitude, na.rm = TRUE),max(files$today_polygon_centroid_longitude, na.rm = TRUE)))
  # mnla<-min(c(min(files$baseline_polygon_centroid_latitude,  na.rm = TRUE),min(files$today_polygon_centroid_latitude,  na.rm = TRUE)))
  # mxla<-max(c(max(files$baseline_polygon_centroid_latitude,  na.rm = TRUE),max(files$today_polygon_centroid_latitude,  na.rm = TRUE)))
  mxlo<-max(files$baseline_polygon_x, na.rm = TRUE)
  mnlo<-min(files$baseline_polygon_x, na.rm = TRUE)
  mxla<-max(files$baseline_polygon_y, na.rm = TRUE)
  mnla<-min(files$baseline_polygon_y, na.rm = TRUE)
  
  onelo<-onela<-FALSE
  if(mnlo==mxlo){mnlo<-mnlo-2.5;mxlo<-mxlo+2.5; print("Only one disaster location, resetting longitudinal map coordinate");onelo<-TRUE}
  if(mnla==mxla){mnla<-mnla-2.5;mxla<-mxla+2.5; print("Only one disaster location, resetting latitudinal map coordinate");onela<-TRUE}
  oneplace<-onelo&onela
  
  # bbox is bounding box in the form 'min lon, max lat, max lon, min lat'
  bbox<-c(mnlo,mxla,mxlo,mnla)
  
  elong<-2#0.05*(mxlo-mnla)
  elat<-2#0.05*(mxla-mnla)
  locGM<-c(left =as.integer(mnlo-elong), 
           bottom = as.integer(mnla-elat), 
           right = as.integer(mxlo+elong),
           top = as.integer(mxla+elat))
  
  annot<-TRUE
  cities<-maps::world.cities%>%filter(lat>bbox[2]&lat<bbox[4]&long>bbox[1]&long<bbox[3])%>%arrange(desc(pop))
  if(length(cities$name)==0){
    annot<-FALSE
    cities<-slice(cities,1:min(length(cities$name),5))
  }
  
  mad_map <- get_map(locGM,maptype = "toner-background")
  
  p<-ggmap(mad_map)+geom_point(filter(disp,cdates %in% c(mnd,mxd)),mapping=aes(x=baseline_polygon_x,y=baseline_polygon_y,size=n_people,colour=as.factor(cdates)), shape = 19,na.rm = TRUE)+
    xlab("Longitude") + ylab("Latitude")+
    scale_size_continuous(name="IDP Stock", breaks=c(400,800,1200,1600), range=c(1,15))+ scale_color_discrete(guide=FALSE)
  if(annot){p<-p+geom_label(data = cities, aes(long, lat, label = name), size = 3, fontface = "bold", nudge_x = 1.5)}
  p<-p+facet_grid(. ~ cdates) + scale_fill_viridis_c() + theme(strip.text.x = element_text(size = 16))
  p
  ggsave(paste0(namer,"_",reference,"_OSM_IDP-maxdayminday.png"), plot=p,path = paste0(directory,'Plots/'),width = 13,height = 5.)
  
  if(oneplace){return(list(disp,notdisp,retn,unknown,dispA))}
  
  mxdisp<-disp%>%filter(cdates==mxd:mxd+5)%>%group_by(baseline_polygon_name)%>%summarise(sum=sum(n_people))%>%arrange(desc(sum))
  # Choose top ten cities by number of people displaced
  topten<-as.character(mxdisp$baseline_polygon_name[1:10])
  
  sumvals<-filter(disp,baseline_polygon_name %in% topten)%>%
    group_by(cdates,baseline_polygon_name)%>%
    summarise(sum=sum(n_people,na.rm = TRUE))%>%arrange(cdates);
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total IDP Stock") 
  p<-p+facet_wrap(. ~ baseline_polygon_name,nrow = 2) + scale_fill_viridis_c() + theme(strip.text.x = element_text(size = 16))
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-IDPs_per-City.png"), plot=p,path = paste0(directory,'Plots/'),width = 13,height = 6.)
  
  sumvals<-filter(notdisp,baseline_polygon_name %in% topten)%>%
    group_by(cdates,baseline_polygon_name)%>%
    summarise(sum=sum(n_people,na.rm = TRUE))%>%arrange(cdates);
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total Not Displaced") 
  p<-p+facet_wrap(. ~ baseline_polygon_name,nrow = 2) + scale_fill_viridis_c() + theme(strip.text.x = element_text(size = 16))
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-Not-Displaced_per-City.png"), plot=p,path = paste0(directory,'Plots/'),width = 13,height = 6.)
  
  sumvals<-filter(retn,baseline_polygon_name %in% topten)%>%
    group_by(cdates,baseline_polygon_name)%>%
    summarise(sum=sum(n_people,na.rm = TRUE))%>%arrange(cdates);
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total Returned") 
  p<-p+facet_wrap(. ~ baseline_polygon_name,nrow = 2) + scale_fill_viridis_c() + theme(strip.text.x = element_text(size = 16))
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-Returned-IDPs_per-City.png"), plot=p,path = paste0(directory,'Plots/'),width = 13,height = 6.)
  
  boundaries<-lapply(unique(files$baseline_polygon_name),function(x) 
    c(getbb(paste0(x),format_out = 'sf_polygon')))
  
  inputs<-c("Charlotte, Saint Vincent",
            "St. Andrew, Saint Vincent",
            "Saint David, Saint Vincent",
            "Saint George, Saint Vincent",
            "Saint Patrick, Saint Vincent")
  boundaries<-lapply(inputs,function(x) 
    c(getbb(paste0(x),format_out = 'sf_polygon')))
  
  p<-ggmap(mad_map)
  # colsy<-colorRampPalette(c("blue", "red"))( length(boundaries) )
  colsy<-c("red","blue","purple","orange","green")
  # colsy<-colorRampPalette(c("blue", "red"))( 3 )
  for (i in 1:length(boundaries)){
  # for (i in 3:5){
    if(!is.null((boundaries[[i]])$polygon)) p<-p+geom_sf(data=(boundaries[[i]])$polygon,
                      inherit.aes = FALSE,
                      alpha=0.1,fill=colsy[i],colour=colsy[i],size=2)
    else  {
      tmp<-st_cast((boundaries[[i]])$multipolygon,"POLYGON")
      tmp$geometry<-tmp$geometry[2]
      p<-p+geom_sf(data=tmp,
                    inherit.aes = FALSE,
                 alpha=0.1,fill=colsy[i],colour=colsy[i],size=2)
    }

  }
  p<-p+xlab("Longitude") + ylab("Latitude");p
  
  return(list(disp,notdisp,retn,unknown,dispA))
  
  #General information on the area bounded by the data:
  
  # ggmap(mad_map)+geom_density_2d(data=disp,aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude),na.rm = TRUE,alpha=0.2)
  # 
  # xgrid <-  seq(mnlo, mxlo, length.out=200)
  # ygrid <-  seq(mnla, mxla, length.out=200)
  # 
  # # Plot three time values: beginning, middle and end values
  # p <- ggmap(mad_map)
  # for (fnn in c("min","median","max")){
  #   fn<-match.fun(fnn)
  #   tmp<-filter(disp,cdates==fn(cdates,na.rm = TRUE))
  #   data.loess <- loess(n_people ~ baseline_polygon_centroid_longitude * baseline_polygon_centroid_latitude, data = tmp)
  #   data.fit <-  expand.grid(baseline_polygon_centroid_longitude = xgrid, baseline_polygon_centroid_latitude = ygrid)
  #   dispy <-  predict(data.loess, newdata = data.fit)
  #   mtrx.melt<-melt(dispy,id.vars=c("baseline_polygon_centroid_longitude", "baseline_polygon_centroid_latitude"),measure.vars = "n_people")
  #   names(mtrx.melt)<-c("baseline_polygon_centroid_longitude", "baseline_polygon_centroid_latitude","n_people")
  #   mtrx.melt$baseline_polygon_centroid_longitude <- as.numeric(str_sub(mtrx.melt$baseline_polygon_centroid_longitude, str_locate(mtrx.melt$baseline_polygon_centroid_longitude, "=")[1,1] + 1))
  #   mtrx.melt$baseline_polygon_centroid_latitude <- as.numeric(str_sub(mtrx.melt$baseline_polygon_centroid_latitude, str_locate(mtrx.melt$baseline_polygon_centroid_latitude, "=")[1,1] + 1))
  #   p<-ggplot(mtrx.melt, aes(x = baseline_polygon_centroid_longitude, y = baseline_polygon_centroid_latitude, z = n_people)) +
  #     #geom_raster(aes(fill = n_people),alpha=0.3)+
  #     stat_contour(geom = "polygon", aes(fill = ..level..),bins=15,alpha=0.3,na.rm = TRUE) +
  #     geom_contour(colour = "white");p
  #   assign(x = paste0("p",fnn),value = 
  #            ggmap(mad_map,base_layer = p ))
  #            
  #            # p+stat_contour(mtrx.melt, mapping=aes(x = baseline_polygon_centroid_longitude, y = baseline_polygon_centroid_latitude, z = n_people,fill = n_people),geom = "polygon",bins=15,alpha=0.3,na.rm = TRUE) +
  #            # #geom_tile(aes(fill = n_people),alpha=0.3) +
  #            # xlab("Longitude")+ylab("Latitude")+ggtitle(name) +
  #            # guides(fill = guide_colorbar(title = "IDP Stock")))
  # }
  # grid.arrange(pmin,pmedian,pmax)
  # 
  # 
  # p<-p+facet_grid(. ~ cdates) + scale_fill_viridis_c();p
  #   
  #   p <- ggplot(mtrx.melt, aes(x = baseline_polygon_centroid_longitude, y = baseline_polygon_centroid_latitude, z = n_people)) +
  #         geom_raster(aes(fill = n_people),alpha=0.3)+
  #     geom_contour(colour = "white")
  #   #stat_contour(geom = "polygon", aes(fill = stat(level)),bins=15,alpha=0.3,na.rm = TRUE) +
  #     #geom_tile(aes(fill = n_people),alpha=0.3) +
  #     xlab("Longitude")+ylab("Latitude")+ggtitle(namer) +
  #     guides(fill = guide_colorbar(title = "IDP Stock"));p  
  
  # p<-ggmap(mad_map)+#geom_raster(data=tmp,aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude,z=n_people),na.rm = TRUE,bins=15,show.legend = FALSE)+
  #       #stat_density_2d(data=tmp,aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude,fill = stat(level)),alpha=0.2, geom = "polygon", na.rm = TRUE,bins=15) +
  #   ggplot(data=tmp,aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude,z=n_people))+stat_contour()
  #       xlab("Longitude")+ylab("Latitude")+ggtitle(name)
  # p<-p+facet_grid(. ~ cdates) + scale_fill_viridis_c();p
  
}




















GetFBEstimates_old<-function(FBdirectory, name, reference="", check=FALSE, plotty=FALSE){

  tmp<-paste(name,collapse = '.*')
  namer=paste(name,collapse = '_')
  dirs <- paste0(FBdirectory,list.files(path=FBdirectory,pattern=tmp))
  #dirs <- paste0(FBdirectory,list.files(path=FBdirectory,pattern=name))
  dfiles<-data.frame()
  for (dd in dirs){
    temp <- paste0(dd,'/',list.files(path=dd,pattern='.*csv'))
    files <- lapply(temp, read.csv)
    # The end of the file string is something like "...20191022.csv":
    cdates<-substr(temp,nchar(temp[1])-11,nchar(temp[1])-4)
    dates <- as.Date(cdates, "%Y%m%d")
    ldays<-length(cdates)
    
    # Check that all dates are included in the decay rate
    if (check){
      dd<-as.numeric(format(dates,'%d'))
      mm<-as.numeric(format(dates,'%m'))
      yy<-as.numeric(format(dates,'%Y'))
      for (i in 2:ldays) {
        if (!(dd[i]==dd[i-1]+1)){
          if (!(mm[i]==mm[i-1]+1)){
            if (!(yy[i]==yy[i-1]+1)){stop(paste0("Not all days are included in the data, download the day ",toString(dates[i])," from GeoInsights"))}
          }
        }
      }
      print("All days are included! Data is good to go!")
    }
    
    names(files)<-cdates
    files<-bind_rows(files, .id = "cdates")
    files$cdates<-as.Date(files$cdates,format = "%Y%m%d")
    dfiles<-rbind(dfiles,files)
  }
  rm(files)
  dfiles<-dfiles%>%mutate(time=as.integer(cdates-min(cdates)))

  #create dataframe with number of displaced, 
  disp<-filter(dfiles,status=="Displaced")
  dispA<-filter(dfiles,status=="Displaced Abroad")
  notdisp<-filter(dfiles,status=="Never Displaced")
  unknown<-filter(dfiles,status=="Unknown")
  retn<-filter(dfiles,status=="Returned")
  
  if(!plotty){return(list(disp,notdisp,retn,unknown,dispA))}
  
  sumvals<-retn%>%group_by(cdates)%>%summarise(max=max(n_people,na.rm = TRUE),min=min(n_people,na.rm = TRUE),sum=sum(n_people,na.rm = TRUE))%>%arrange(cdates)
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total Returned IDP Stock")
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-Returned.png"), plot=p,path = paste0(directory,'Plots/'),width = 5,height = 4.)
  sumvals<-notdisp%>%group_by(cdates)%>%summarise(max=max(n_people,na.rm = TRUE),min=min(n_people,na.rm = TRUE),sum=sum(n_people,na.rm = TRUE))%>%arrange(cdates)
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total Never Displaced")
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-NeverDisp.png"), plot=p,path = paste0(directory,'Plots/'),width = 5,height = 4.)
  sumvals<-disp%>%group_by(cdates)%>%summarise(sum=sum(n_people,na.rm = TRUE))
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total IDP Stock")
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-IDPs.png"), plot=p,path = paste0(directory,'Plots/'),width = 5,height = 4.)
  
  mxd=sumvals$cdates[which.max(sumvals$sum)]
  mnd=sumvals$cdates[which.min(sumvals$sum[sumvals$cdates>mxd])]
  
  # mnlo<-min(c(min(dfiles$baseline_polygon_centroid_longitude, na.rm = TRUE),min(dfiles$today_polygon_centroid_longitude, na.rm = TRUE)))
  # mxlo<-max(c(max(dfiles$baseline_polygon_centroid_longitude, na.rm = TRUE),max(dfiles$today_polygon_centroid_longitude, na.rm = TRUE)))
  # mnla<-min(c(min(dfiles$baseline_polygon_centroid_latitude,  na.rm = TRUE),min(dfiles$today_polygon_centroid_latitude,  na.rm = TRUE)))
  # mxla<-max(c(max(dfiles$baseline_polygon_centroid_latitude,  na.rm = TRUE),max(dfiles$today_polygon_centroid_latitude,  na.rm = TRUE)))
  mxlo<-max(dfiles$baseline_polygon_centroid_longitude, na.rm = TRUE)
  mnlo<-min(dfiles$baseline_polygon_centroid_longitude, na.rm = TRUE)
  mxla<-max(dfiles$baseline_polygon_centroid_latitude, na.rm = TRUE)
  mnla<-min(dfiles$baseline_polygon_centroid_latitude, na.rm = TRUE)
  
  onelo<-onela<-FALSE
  if(mnlo==mxlo){mnlo<-mnlo-2.5;mxlo<-mxlo+2.5; print("Only one disaster location, resetting longitudinal map coordinate");onelo<-TRUE}
  if(mnla==mxla){mnla<-mnla-2.5;mxla<-mxla+2.5; print("Only one disaster location, resetting latitudinal map coordinate");onela<-TRUE}
  oneplace<-onelo&onela
  
  # bbox is bounding box in the form 'min lon, max lat, max lon, min lat'
  bbox<-c(mnlo,mxla,mxlo,mnla)
  
  elong<-2#0.05*(mxlo-mnla)
  elat<-2#0.05*(mxla-mnla)
  locGM<-c(left =as.integer(mnlo-elong), 
           bottom = as.integer(mnla-elat), 
           right = as.integer(mxlo+elong),
           top = as.integer(mxla+elat))
  
  annot<-TRUE
  cities<-maps::world.cities%>%filter(lat>bbox[2]&lat<bbox[4]&long>bbox[1]&long<bbox[3])%>%arrange(desc(pop))
  if(length(cities$name)==0){
    annot<-FALSE
    cities<-slice(cities,1:min(length(cities$name),5))
  }
  
  mad_map <- get_map(locGM,maptype = "toner-background")
  
  p<-ggmap(mad_map)+geom_point(filter(disp,cdates %in% c(mnd,mxd)),mapping=aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude,size=n_people,colour=as.factor(cdates)), shape = 19,na.rm = TRUE)+
    xlab("Longitude") + ylab("Latitude")+
    scale_size_continuous(name="IDP Stock", breaks=c(400,800,1200,1600), range=c(1,15))+ scale_color_discrete(guide=FALSE)
  if(annot){p<-p+geom_label(data = cities, aes(long, lat, label = name), size = 3, fontface = "bold", nudge_x = 1.5)}
  p<-p+facet_grid(. ~ cdates) + scale_fill_viridis_c() + theme(strip.text.x = element_text(size = 16))
  print(p)
  ggsave(paste0(namer,"_",reference,"_OSM_IDP-maxdayminday.png"), plot=p,path = paste0(directory,'Plots/'),width = 13,height = 5.)
  
  if(oneplace){return(list(disp,notdisp,retn,unknown,dispA))}
  
  mxdisp<-disp%>%filter(cdates==mxd:mxd+5)%>%group_by(baseline_polygon_name)%>%summarise(sum=sum(n_people))%>%arrange(desc(sum))
  # Choose top ten cities by number of people displaced
  topten<-as.character(mxdisp$baseline_polygon_name[1:10])
  
  sumvals<-filter(disp,baseline_polygon_name %in% topten)%>%
    group_by(cdates,baseline_polygon_name)%>%
    summarise(sum=sum(n_people,na.rm = TRUE))%>%arrange(cdates);
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total IDP Stock") 
  p<-p+facet_wrap(. ~ baseline_polygon_name,nrow = 2) + scale_fill_viridis_c() + theme(strip.text.x = element_text(size = 16))
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-IDPs_per-City.png"), plot=p,path = paste0(directory,'Plots/'),width = 13,height = 6.)
  
  sumvals<-filter(notdisp,baseline_polygon_name %in% topten)%>%
    group_by(cdates,baseline_polygon_name)%>%
    summarise(sum=sum(n_people,na.rm = TRUE))%>%arrange(cdates);
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total Not Displaced") 
  p<-p+facet_wrap(. ~ baseline_polygon_name,nrow = 2) + scale_fill_viridis_c() + theme(strip.text.x = element_text(size = 16))
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-Not-Displaced_per-City.png"), plot=p,path = paste0(directory,'Plots/'),width = 13,height = 6.)
  
  sumvals<-filter(retn,baseline_polygon_name %in% topten)%>%
    group_by(cdates,baseline_polygon_name)%>%
    summarise(sum=sum(n_people,na.rm = TRUE))%>%arrange(cdates);
  p<-ggplot(sumvals,aes(x=cdates,y=sum))+geom_point(shape=4,size=3)+geom_line()+scale_x_date(labels = scales::date_format("%a %d-%m")) +
    xlab("Date")+ylab("Total Returned") 
  p<-p+facet_wrap(. ~ baseline_polygon_name,nrow = 2) + scale_fill_viridis_c() + theme(strip.text.x = element_text(size = 16))
  print(p)
  ggsave(paste0(namer,"_",reference,"_Total-Returned-IDPs_per-City.png"), plot=p,path = paste0(directory,'Plots/'),width = 13,height = 6.)

  return(list(disp,notdisp,retn,unknown,dispA))
  
  #General information on the area bounded by the data:
    
    # ggmap(mad_map)+geom_density_2d(data=disp,aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude),na.rm = TRUE,alpha=0.2)
    # 
    # xgrid <-  seq(mnlo, mxlo, length.out=200)
    # ygrid <-  seq(mnla, mxla, length.out=200)
    # 
    # # Plot three time values: beginning, middle and end values
    # p <- ggmap(mad_map)
    # for (fnn in c("min","median","max")){
    #   fn<-match.fun(fnn)
    #   tmp<-filter(disp,cdates==fn(cdates,na.rm = TRUE))
    #   data.loess <- loess(n_people ~ baseline_polygon_centroid_longitude * baseline_polygon_centroid_latitude, data = tmp)
    #   data.fit <-  expand.grid(baseline_polygon_centroid_longitude = xgrid, baseline_polygon_centroid_latitude = ygrid)
    #   dispy <-  predict(data.loess, newdata = data.fit)
    #   mtrx.melt<-melt(dispy,id.vars=c("baseline_polygon_centroid_longitude", "baseline_polygon_centroid_latitude"),measure.vars = "n_people")
    #   names(mtrx.melt)<-c("baseline_polygon_centroid_longitude", "baseline_polygon_centroid_latitude","n_people")
    #   mtrx.melt$baseline_polygon_centroid_longitude <- as.numeric(str_sub(mtrx.melt$baseline_polygon_centroid_longitude, str_locate(mtrx.melt$baseline_polygon_centroid_longitude, "=")[1,1] + 1))
    #   mtrx.melt$baseline_polygon_centroid_latitude <- as.numeric(str_sub(mtrx.melt$baseline_polygon_centroid_latitude, str_locate(mtrx.melt$baseline_polygon_centroid_latitude, "=")[1,1] + 1))
    #   p<-ggplot(mtrx.melt, aes(x = baseline_polygon_centroid_longitude, y = baseline_polygon_centroid_latitude, z = n_people)) +
    #     #geom_raster(aes(fill = n_people),alpha=0.3)+
    #     stat_contour(geom = "polygon", aes(fill = ..level..),bins=15,alpha=0.3,na.rm = TRUE) +
    #     geom_contour(colour = "white");p
    #   assign(x = paste0("p",fnn),value = 
    #            ggmap(mad_map,base_layer = p ))
    #            
    #            # p+stat_contour(mtrx.melt, mapping=aes(x = baseline_polygon_centroid_longitude, y = baseline_polygon_centroid_latitude, z = n_people,fill = n_people),geom = "polygon",bins=15,alpha=0.3,na.rm = TRUE) +
    #            # #geom_tile(aes(fill = n_people),alpha=0.3) +
    #            # xlab("Longitude")+ylab("Latitude")+ggtitle(name) +
    #            # guides(fill = guide_colorbar(title = "IDP Stock")))
    # }
    # grid.arrange(pmin,pmedian,pmax)
    # 
    # 
    # p<-p+facet_grid(. ~ cdates) + scale_fill_viridis_c();p
    #   
    #   p <- ggplot(mtrx.melt, aes(x = baseline_polygon_centroid_longitude, y = baseline_polygon_centroid_latitude, z = n_people)) +
    #         geom_raster(aes(fill = n_people),alpha=0.3)+
    #     geom_contour(colour = "white")
    #   #stat_contour(geom = "polygon", aes(fill = stat(level)),bins=15,alpha=0.3,na.rm = TRUE) +
    #     #geom_tile(aes(fill = n_people),alpha=0.3) +
    #     xlab("Longitude")+ylab("Latitude")+ggtitle(namer) +
    #     guides(fill = guide_colorbar(title = "IDP Stock"));p  
        
        # p<-ggmap(mad_map)+#geom_raster(data=tmp,aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude,z=n_people),na.rm = TRUE,bins=15,show.legend = FALSE)+
        #       #stat_density_2d(data=tmp,aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude,fill = stat(level)),alpha=0.2, geom = "polygon", na.rm = TRUE,bins=15) +
        #   ggplot(data=tmp,aes(x=baseline_polygon_centroid_longitude,y=baseline_polygon_centroid_latitude,z=n_people))+stat_contour()
        #       xlab("Longitude")+ylab("Latitude")+ggtitle(name)
        # p<-p+facet_grid(. ~ cdates) + scale_fill_viridis_c();p
  
}
  
  