library(ggplot2)
library(OpenStreetMap)
library(osmdata)
# Look at https://wiki.openstreetmap.org/wiki/Map_Features
library(sf)
library(ggmap)
# directory<-'/home/patten/Documents/IDMC/Facebook_Data/'
# folder<-c('Wildfire','East','New','South')

GetNetworkCoverage<-function(directory,folder,plotty=FALSE){
  
  name<-paste(folder,collapse = '.*')
  dirs <- paste0(directory,list.files(path=directory,pattern=name))
  temp <- paste0(dirs,'/Network_Coverage/',list.files(path=paste0(dirs,"/Network_Coverage"),pattern='.*csv'))
  NetCov<-read.csv(temp,header = TRUE,fill = TRUE,sep = ","); NetCov[is.na(NetCov)]<-0L
  
  if(plotty){
    
    bbox<-c(min(NetCov$lon),min(NetCov$lat),max(NetCov$lon),max(NetCov$lat))
    mad_map <- get_map(bbox,source = "stamen",maptype = "toner")
    
    p<-ggmap(mad_map)
    p+geom_contour(data=NetCov,mapping=aes(x = lon,y=lat,z=Count.4g,fill=Count.4g,color=Count.4g),bins=50,alpha=0.5)
    
    
    
    
    p<- p + stat_density2d(mapping = aes(x=lon, y=lat,fill=as.numeric(Count.4g)), 
                                       data=NetCov, bins=5, geom="polygon") + 
      #scale_fill_gradient2(low = "yellow", high = "red") +
      #scale_fill_continuous(name = "Number of Fires", limits = c(0.02,0.05), labels = seq(0.02,0.05,0.01), breaks = seq(0.02,0.05,0.01)) +
      #scale_alpha_continuous(range = c(0.02, 0.1)) + #, guide = FALSE) #+
      geom_density2d(colour="black", bins=10, mapping = aes(x=lon, y=lat,fill=Count.4g), data=NetCov) + 
      xlab("Longitude") + ylab("Latitude") + ggtitle("Network Coverage") + theme(plot.title = element_text(hjust = 0.5))
    print(p)
    
    m <- ggplot(NetCov, aes(x = lon, y = lat,z=Count.4g)) +
      xlab("Longitude") + ylab("Latitude") +
    #geom_raster(fill=NetCov$Count.4g) 
    geom_contour(colour = "white");m
    m <- m + stat_density_2d(aes(fill = Count.4g), geom = "polygon")
    print(m)
    ggsave(paste0(folder,"_NetworkCoverage4G.png"), plot=m,path = paste0(directory,'Plots/'),width = 5,height = 4.)
    
    m <- ggplot(NetCov, aes(x = lon, y = lat)) +
      xlim(0.5, 6) +  ylim(40, 110) + xlab("Longitude") + ylab("Latitude")
    # m + geom_density_2d()
    m <- m + stat_density_2d(aes(fill = stat(Count.3g)), geom = "polygon")
    print(m)
    ggsave(paste0(folder,"_NetworkCoverage3G.png"), plot=m,path = paste0(directory,'Plots/'),width = 5,height = 4.)
    
  }
  
  return(NetCov)
  
}