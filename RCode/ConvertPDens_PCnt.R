# Convert from pop density to pop count
Pdens<-c(Pdens=0)
nPdens<-0

ufiles<-list.files(path=paste0(directory,"IIDIPUS_Input/ODDobjects_full"),pattern=Model$haz,recursive = T,ignore.case = T)
for(i in 1:length(ufiles)){
  ODDy<-readRDS(paste0(directory,"IIDIPUS_Input/ODDobjects_full/",ufiles[i]))
  ODDy@fIndies<-Model$fIndies
  ODDy@gmax%<>%as.data.frame.list()
  
  bbox<-ODDy@bbox
  bbox[1]<-bbox[1]-0.5*ODDy@grid@cellsize[1]
  bbox[3]<-bbox[3]+0.5*ODDy@grid@cellsize[1]
  bbox[2]<-bbox[2]-0.5*ODDy@grid@cellsize[2]
  bbox[4]<-bbox[4]+0.5*ODDy@grid@cellsize[2]
  
  obj<-GetPopulationBbox(directory,bbox=bbox,density = F)
  
  inside<-obj@coords[,1]>=ODDy@bbox[1] & obj@coords[,1]<=ODDy@bbox[3] &
          obj@coords[,2]>=ODDy@bbox[2] & obj@coords[,2]<=ODDy@bbox[4]
  
  if(length(obj@data$Population[inside])!=nrow(ODDy)) stop(paste0("didn't work for - ",i))
  
  ODDy@data$Population<-obj@data$Population[inside]
  
  # FIND ALL NA HAZARD VALUES AND REPLACE 
  ind<-apply(ODDy@data[grep("hazMean",names(ODDy),value = T)],1,function(i) all(is.na(i)))
  ODDy@data$Population[ind]<-NA
  
  iso3c<-unique(ODDy@data$ISO3C[!ind]) ; iso3c<-iso3c[!is.na(iso3c)]
  Popfactors<-InterpPopWB(iso3c,min(ODDy@hazdates))
  
  for (iso in iso3c){
    indie<-ODDy@data$ISO3C==iso & !is.na(ODDy@data$ISO3C)
    ODDy@data$Population[indie]<-Popfactors$factor[Popfactors$iso3==iso]*ODDy@data$Population[indie]
  }
  
  Pdens<-Pdens+sum(log(ODDy@data$Population[ODDy@data$Population>0]),na.rm=T)
  nPdens<-nPdens+length(ODDy@data$Population[ODDy@data$Population>0 & !is.na(ODDy@data$Population)])
  
  saveRDS(ODDy,paste0("./","IIDIPUS_Input/ODDobjects_count/",ufiles[i]))
  
}

center<-readRDS(paste0(dir,"IIDIPUS_Input/centerings"))
center$Pdens<-Pdens/nPdens
saveRDS(center,paste0(dir,"IIDIPUS_Input/centerings"))

# BD DATASET POP DENS SWAPOUT


ufiles<-list.files(path=paste0(directory,"IIDIPUS_Input/BDobjects_save"),pattern=Model$haz,recursive = T,ignore.case = T)
for(i in 1:length(ufiles)){
  BDy<-readRDS(paste0(directory,"IIDIPUS_Input/BDobjects_save/",ufiles[i]))
  BDy@fIndies<-Model$fIndies
  
  ODDy<-readRDS(paste0(directory,"IIDIPUS_Input/ODDobjects_count/",ufiles[i]))
  
  tmp<-ODDy["Population"]%>%raster%>%raster::extract(BDy@coords)%>%data.frame()
  BDy@data$Population<-tmp[,1]
  tmp<-ODDy["GDP"]%>%raster%>%raster::extract(BDy@coords)%>%data.frame()
  BDy@data$GDP<-tmp[,1]
  
  BDy2<-readRDS(paste0(directory,"IIDIPUS_Input/BDobjects/",ufiles[i]))
  BDy@cIndies<-BDy2@cIndies
  BDy@coefs<-BDy2@coefs
  BDy@buildingsfile<-BDy2@buildingsfile
  
  saveRDS(BDy,paste0("./","IIDIPUS_Input/BDobjects_count/",ufiles[i]))
  
}

