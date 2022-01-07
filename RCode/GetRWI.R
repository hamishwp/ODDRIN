list.of.packages <- c("sf","tidyverse","gstat","dplyr","magrittr","sp")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gstat)
library(dplyr)
library(magrittr)
library(tidyverse)
library(sf)
library(sp)

#@@@@@@@@@@@@@@@@@@@ Extract and Interpolate the Relative Wealth Index @@@@@@@@@@@@@@@@@@@#
# Examples of the usage:
# rawRWI<-GetRWI(RWIfolder="./FB_RWI/",iso="ecu")
# interpolatedRWI<-GetRWI(RWIfolder="./FB_RWI/",iso="ecu",coords=SpatialPointsDataFrameObject,GPR=F)
# Where SpatialPointsObject is built in the following way:
#       SpatialPointsDataFrame(data.frame(Longitude=LongitudeValuesVector,Latitude=LatitudeValuesVector),
#                             data = data.frame(Intensity=rep(NA_real_,length(LongitudeValuesVector)),
#                                               Error=rep(NA_real_,length(LongitudeValuesVector))),
#                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

# Extract the RWI file from the folder and read in (assuming in .csv format)
ExtractRWI<-function(fRWI){
  
  if(length(fRWI)==0) {stop("no folder-file name provided to extract the RWI from")
  } else if (length(fRWI)>1){
    print("Multiple RWI files provided, concatenating into one data frame")
    RWI<-data.frame()
    for(i in 1:length(fRWI)) RWI%<>%rbind(as.data.frame(read.csv(fRWI[i])))
  } else{
    RWI<-read.csv(fRWI)
  }
  
  names(RWI)<-c("Latitude","Longitude","Intensity","Error")
  # array <- SpatialPointsDataFrame(coords = RWI[c("lon","lat")],data = RWI["Rwi"])
  RWI <- SpatialPointsDataFrame(coords = RWI[c("Longitude","Latitude")],data = RWI[c("Intensity","Error")])
  proj4string(RWI)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
  
  return(RWI)
  
}

# This function checks that the coordinates to interpolate on are in the correct format
ChangeCoords<-function(coords){
  
  # coords%<>%as.data.frame()
  if(any(dim(coords)==2)){
    if(which(dim(coords)==2)==1) rowz<-T else rowz<-F
    if(rowz) {namers<-rownames(coords)} else namers<-colnames(coords)
    
    if(is.null(namers)) {
      print("Warning: assuming first coordinate is Longitude and second is Latitude")
      namers<-c("Longitude","Latitude")
    }
    if(rowz) {
      lon<-coords[which(grepl("long",namers,ignore.case = T)),]
      lat<-coords[which(grepl("lat",namers,ignore.case = T)),]
    } else{
      lon<-coords[,which(grepl("long",namers,ignore.case = T))]
      lat<-coords[,which(grepl("lat",namers,ignore.case = T))]
    }
    # bbox<-c(min(lon,na.rm = T),min(lat,na.rm = T),max(lon,na.rm = T),max(lat,na.rm = T))  
    
    return(SpatialPointsDataFrame(data.frame(Longitude=lon,Latitude=lat),
                                  data = data.frame(Intensity=rep(NA_real_,length(lon)),Error=rep(NA_real_,length(lon))),
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")))
    
  } else if(class(coords)=="SpatialPoints") {
    
    return(SpatialPointsDataFrame(data.frame(Longitude=coords@coords[,1],Latitude=coords@coords[,2]),
                                  data = data.frame(Intensity=rep(NA_real_,nrow(coords@coords)),Error=rep(NA_real_,nrow(coords@coords))),
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")))
    
  } else stop("Please try coordinates that have 2 columns (or rows) as [Longitude,Latitude] instead of as a matrix")
  
}

# Instead of Gaussian Process Regression, using the weighted mean can cut down considerably on the compuation
WeightedMeanRWI<-function(interpRWI,RWI){
  
  # Note that the weighted mean could be based on a parameteric function of the distance, not just the squared distance
  # Alternatively, using X-fold cross validation, we could parameterise the function using a bespoke method
  # Also, we could use mclapply to speed things up, depending on the users computer infrastructure
  # but for simplicity, we use distance^2
  tmp<-apply(as.data.frame(interpRWI)[,c("Longitude","Latitude")],1,function(x){
    disty<-geosphere::distHaversine(x,coordinates(RWI))
    return(list(Exp=weighted.mean(RWI$Intensity,1./(disty^2),na.rm=TRUE),
                Error=weighted.mean(RWI$Error,1./(disty^2),na.rm=TRUE)))
  })
  for(i in 1:length(interpRWI)){
    interpRWI$Intensity[i]<-unlist(tmp[[i]][1])
    interpRWI$Error[i]<-unlist(tmp[[i]][2])
  }
  
  return(interpRWI)
}

# Using Gaussian Process Regression instead of the weighted mean
KrigMeUp<-function(interpRWI,RWI){

  if(class(RWI)!="SpatialPointsDataFrame") stop("wrong RWI points class, check the raw RWI data")
  # RWI%<>%as("SpatialPointsDataFrame")

  # Variogram for RWI Intensity
  print("Krigging the RWI Intensity")
  vario <- gstat::variogram(Intensity~1, RWI)
  fit <- gstat::fit.variogram(vario, model=gstat::vgm(model="Exp"))
  # Krig the values
  interpRWI$Intensity<-(gstat::krige(Intensity ~ 1, RWI, interpRWI, model=fit))@data$var1.pred
  # Variogram for RWI Error
  print("Krigging the RWI Error")
  vario <- gstat::variogram(Error~1, RWI)
  fit <- gstat::fit.variogram(vario, model=gstat::vgm(model="Exp"))
  # Krig the values
  interpRWI$Error<-(gstat::krige(Error ~ 1, RWI, interpRWI, model=fit))@data$var1.pred
  
  return(interpRWI)
  
}

intrpRWI<-function(RWI,coords,GPR=F){
  
  # What form are the coords in?
  interpRWI<-tryCatch({ChangeCoords(coords)}, error=function(e) NA)
  if(class(interpRWI)!="SpatialPointsDataFrame") stop("problem with the interpolation coordinates, try inputing as data.frame(Longitude=...,Latitude=...)")
  # What kind of interpolation technique?
  if(GPR){
    interpRWI%<>%KrigMeUp(RWI=RWI)
  } else{
    interpRWI%<>%WeightedMeanRWI(RWI=RWI)
  }
  
  return(interpRWI)
  
}

FormRWIfilename<-function(filer=NULL,RWIfolder=NULL,iso=NULL){
  
  # Find the folder location of RWI file
  if(is.null(RWIfolder)) RWIfolder<-getwd()
  # Housekeeping
  if(!endsWith(RWIfolder,"/")) RWIfolder<-paste0(RWIfolder,"/")
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Find the RWI file(s) @@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  # If country ISO3C is provided, use this to extract the appropriate RWI file
  if(!is.null(iso)) {
    if(nchar(iso)!=3) stop("wrong format of ISO3C code, needs to be three letters, e.g. 'ECU'")
    fRWI<-grep(list.files(RWIfolder),pattern = str_to_lower(iso),value = T)
    if(length(fRWI)==0) stop("ISO3C given instead of filename, but no file was found with name including the ISO3C provided")  
    fRWI<-paste0(RWIfolder,fRWI)
    # If filename is provided, use it to extract the RWI file
  } else if(!is.null(filer)) {
    fRWI<-paste0(RWIfolder,filer)
    # Otherwise, try searching for any .csv file in the folder name provided
  } else {
    fRWI<-grep(list.files(RWIfolder),pattern = ".csv",value = T)
    if(length(fRWI)==0) stop("no filename or ISO3C provided, and no other CSV files were found inside the folder provided")  
    fRWI<-paste0(RWIfolder,fRWI)
  }
  
  return(fRWI)
  
}

GetRWI<-function(filer=NULL,RWIfolder=NULL,iso=NULL,coords=NULL,GPR=F){
  
  # Get the filename and location to extract
  fRWI<-FormRWIfilename(filer=filer,RWIfolder=RWIfolder,iso=iso)
  # Extract RWI function - into spatialpointsdataframe
  RWI<-ExtractRWI(fRWI)
  # Check if we need to interpolate the RWI onto a provided grid
  if(is.null(coords)) return(RWI)
  # Interpolate function if coords are given
  return(intrpRWI(RWI,coords,GPR))
  
}
