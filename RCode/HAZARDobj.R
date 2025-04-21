# library(ggplot2)
# library(sf)
# library(tidyverse)
# library(sp)
# library(magrittr)
# library(dplyr)

ExtractParams<-function(haz="EQ"){
  if(haz=="EQ") return(list(I0=4.3,minmag=5))
  stop("Hazard type not recognised")
}

checkHAZARD<-function(object){

  if(!object@hazard%in%c("EQ","TC","FL")) stop("HAZARD object must be either TC, FL or EQ hazards")
  if(!object@alertlevel%in%c("red","orange","yellow","green", "null")) stop("HAZARD object must have either red, orange, yellow or green alertlevels")
    if(!(is.na(object@alertscore) | (object@alertscore>=0 & object@alertscore<10))) stop("HAZARD object must have 0<=alertscores<10")
  if(!(object@eventdate>=as.Date("2008-01-01") & object@eventdate<=Sys.Date())) stop("HAZARD object must have dates between now and 2008")
  
  TRUE
}

setClass("HAZARD", 
         slots = c(hazard="character",
                   alertscore="numeric",
                   alertlevel="character",
                   I0="numeric",
                   eventdate="Date", 
                   depth="numeric", 
                   magnitude="numeric", 
                   max_mmi="numeric",
                   USGS_id="character",
                   eventtime="character",
                   alertfull='list'),
         contains = "SpatRaster")

saveHAZ <- function(lhazSDF, path){
  for (i in 1:length(lhazSDF)){
    if(class(lhazSDF[[i]])=='HAZARD'){
      haz_list = list()
      slotnames = slotNames(lhazSDF[[i]])
      pointer_slot <- ifelse('ptr' %in% slotnames, 'ptr', ifelse('pnt' %in% slotnames, 'pnt', 'pntr'))
      for (slot in slotnames[slotnames!=pointer_slot]){
        haz_list[[slot]] = slot(lhazSDF[[i]], slot)
      }
      haz_list$spatrast <- wrap(lhazSDF[[i]])
      lhazSDF[[i]] = haz_list
    }
  }
  saveRDS(lhazSDF, path)
}

readHAZ <- function(path){
  lhazSDF <- readRDS(path)
  for (i in 1:length(lhazSDF)){
    if(!is.null(lhazSDF[[i]]$spatrast)){
      .Object <- new('HAZARD')
      slotnames <- slotNames(lhazSDF[[i]])
      pointer_slot <- ifelse('ptr' %in% slotnames, 'ptr', ifelse('pnt' %in% slotnames, 'pnt', 'pntr'))
      for (slot in slotnames[slotnames!=pointer_slot]){
        slot(.Object, slot) = lhazSDF[[i]][[slot]]
      }
      slot(.Object, pointer_slot) = slot(unwrap(lhazSDF[[i]]$spatrast), pointer_slot)
      lhazSDF[[i]] <- .Object
    }
  }
  return(lhazSDF)
}

setMethod(f="initialize", signature="HAZARD",
          # definition=function(.Object,bbox,hazSDF,dater=NULL,dir=directory,
          definition=function(.Object,obj=NULL,hazard=NULL,dater=NULL,I0=NULL,alertscore=NULL,alertlevel=NULL, 
                              depth=NULL, magnitude=NULL, max_mmi=NULL, eventtime=NULL, USGS_id=NULL, alertfull=NULL) {
            
            if(is.null(hazard)) {
              print("WARNING: no hazard type provided in HAZARD object initialisation, returning empty")
              return(.Object)
            }
            
            .Object@hazard<-hazard
            
            .Object@alertscore<-ifelse(is.null(alertscore),0,alertscore)
            .Object@alertlevel<-ifelse(is.null(alertlevel),"null",alertlevel)
            .Object@alertfull<- if(is.null(alertfull)){NULL} else {alertfull}
            .Object@I0<-        ifelse(is.null(I0),ExtractParams(hazard),I0)
            .Object@depth<-     ifelse(is.null(depth),NA, depth)
            .Object@magnitude<- ifelse(is.null(magnitude),NA, magnitude)
            .Object@max_mmi<-   ifelse(is.null(max_mmi),NA, max_mmi)
            .Object@eventtime<- ifelse(is.null(eventtime),NA, eventtime)
            .Object@USGS_id<-   ifelse(is.null(USGS_id),NA, USGS_id)
            if(!is.null(dater)) .Object@eventdate<-dater
            
            if(!is.null(obj)){

              inI0<-xyFromCell(obj, which(values(!is.na(obj[['mmi_mean']]) & obj[['mmi_mean']]>.Object@I0)))
              
              # Take a bounding box a little larger than the I>I0 object but inside original (the 5 is arbitrary)
              bbox<-c(max(ext(obj)[1],min(inI0[,1])-5*res(obj)[1]),
                      max(ext(obj)[3],min(inI0[,2])-5*res(obj)[2]),
                      min(ext(obj)[2],max(inI0[,1])+5*res(obj)[1]),
                      min(ext(obj)[4],max(inI0[,2])+5*res(obj)[2]))
              
              # Crop
              #e <- raster::extent(c(bbox[c(1,3,2,4)]))
              #proj4string(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"

              obj%<>%terra::crop(extent(c(bbox[c(1,3,2,4)])))
              .Object@pntr <- obj@pntr
              # .Object@file <- obj@file
              # .Object@data <- obj@data
              # .Object@legend <- obj@legend
              # .Object@title <- obj@title
              # .Object@extent <- obj@extent
              # .Object@rotated <- obj@rotated
              # .Object@rotation <- obj@rotation
              # .Object@ncols <- obj@ncols
              # .Object@nrows <- obj@nrows
              # .Object@crs <- obj@crs
              # .Object@srs <- obj@srs
              # .Object@history <- obj@history
              # .Object@z <- obj@z
            }
            #.Object@proj4string <-CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            names(.Object)<-c("mean","sd")
            
            print("Checking HAZARD values")
            checkHAZARD(.Object)
            
            return(.Object)
          }
)

# Convert from ISO3C to country ("GBR"="United Kingdom")
setGeneric("transIso", function(ODD) 
  standardGeneric("transIso") )
setMethod("transIso", "ODD", function(ODD)
  return(countrycode::countrycode(sourcevar = ODD@iso3,
                                  origin = "iso3c",
                                  destination = "country.name")))