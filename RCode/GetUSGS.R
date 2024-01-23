# library(ggplot2)
# library(sf)
# library("ggmap")
# library(OpenStreetMap)
# library(osmdata)
# library(tidyverse)
# library(geojsonR)
# library(raster)
# library(dplyr)
# library(magrittr)
library(xml2)
# source('RCode/HAZARDobj.R')

# Create an object of the required form from the USGS data
formUSGSobject<-function(meanhaz,sdhaz,I0=NULL){
  
  sgdf <- as(meanhaz, 'SpatialPixelsDataFrame') ; rm(meanhaz)
  tmp<- as(sdhaz, 'SpatialPixelsDataFrame') ; rm(sdhaz)
  sgdf$mmi_std<-tmp$mmi_std ; rm(tmp)
  proj4string(sgdf)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
  
  colnames(sgdf@coords)<-rownames(sgdf@bbox)<-c("Longitude","Latitude")
  
  if(!is.null(I0)) sgdf<-sgdf[sgdf$mmi_mean>I0,]
  
  return(sgdf)
  
}

# Extract the raster object from an API call to USGS and extract the zip folder
ExtractUSGS<-function(url,namer,I0=NULL,plotty=F){
  
  temp<-paste0(namer,".zip")
  # Download the raster file to the location 'temp'
  download.file(url,temp)
  # Unpack the files in the zip document
  unzip(paste0(temp),exdir = paste0(namer,"/"))
  # Extract the mean hazard intensity from raster
  meanhaz<-raster(file.path(namer,"mmi_mean.flt"))
  # Extract the variance of the hazard intensity from raster
  sdhaz<-raster(file.path(namer,"mmi_std.flt"))
  unlink(temp)
  
  # Form a standard USGS object
  sgdf<-formUSGSobject(meanhaz,sdhaz,I0) ; rm(meanhaz,sdhaz)
  
  return(sgdf)
  
}

# Using an API call, search through USGS database for a specific EQ & pre/aftershocks
SearchUSGSbbox<-function(bbox,sdate,fdate=NULL,minmag=5){
  
  debut<-"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
  geojsonR::FROM_GeoJson(paste0(debut,"&starttime=",as.Date(sdate)-1,"&endtime=",as.Date(fdate),
                                "&minlongitude=",bbox[1],"&minlatitude=",bbox[2],
                                "&maxlongitude=",bbox[3],"&maxlatitude=",bbox[4],
                                "&minmagnitude=",minmag,"&orderby=magnitude",
                                "&producttype=shakemap"))%>%return
  
}

check_preceding_hazards <- function(HAZobj){
  sdate <- HAZobj$hazard_info$sdate
  fdate <- HAZobj$hazard_info$fdate
  bbox <- HAZobj$hazard_info$bbox
  minmag=4
  debut<-"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
  haz_alls <- geojsonR::FROM_GeoJson(paste0(debut,"&starttime=",as.Date(sdate)-3,"&endtime=",as.Date(fdate)+2,
                                            "&minlongitude=",bbox[1],"&minlatitude=",bbox[2],
                                            "&maxlongitude=",bbox[3],"&maxlatitude=",bbox[4],
                                            "&minmagnitude=",minmag, "&producttype=shakemap"))
  
  getDateTimes <- function(x){
    coords <- x$geometry$coordinates
    timezone <- suppressWarnings(tz_lookup_coords(lat = coords[2], lon = coords[1], method = "fast"))
    eventtime <- as.POSIXct(x$properties$time / 1000, origin = "1970-01-01", tz = timezone)
    return( format(eventtime,  format = "%Y-%m-%d %H:%M:%S"))
  }
  
  eventtimes_all <- tryCatch(sort(unlist(lapply(haz_alls$features, function(x) GetUSGSdatetime(x$id)))),
                             error=function(e) NULL)
  
  if (is.null(eventtimes_all)){
    eventtimes_all <- sort(unlist(lapply(haz_alls$features, getDateTimes)))
  }
  
  eventtimes_modelled <- paste(HAZobj$hazard_info$eventdates, HAZobj$hazard_info$eventtimes)
  return(eventtimes_modelled==eventtimes_all[1])
}

# Make sure the USGS data wasn't empty or entirely outside of the specified boundary box
check_hazsdf<-function(hazsdf=NULL,minmag,bbox=NULL){
  if(is.null(hazsdf)) return(F)
  # Check if the bounding box of the hazard lies within the specified search area
  if(!is.null(bbox)){
    if (all(hazsdf@bbox[c(1,3)]<bbox[1]) | all(hazsdf@bbox[c(1,3)]>bbox[3]) |
        all(hazsdf@bbox[c(2,4)]<bbox[2]) | all(hazsdf@bbox[c(2,4)]>bbox[4])) {
      return(F)
    }
  }
  # Check that the maximum intensity of the earthquake is higher than our limit
  if(max(hazsdf@data$mmi_mean)<minmag) return(F)
  return(T)
}

  GetUSGS_id<-function(USGSid,titlz="tmp",I0=4.5,minmag=5){
  
  url<-paste0("https://earthquake.usgs.gov/fdsnws/event/1/query?eventid=",USGSid,"&format=geojson")
  tmp<-FROM_GeoJson(url)
  hazsdf<-tryCatch(ExtractUSGS(url = tmp$properties$products$shakemap[[1]]$contents$`download/raster.zip`$url,
                               namer = paste0(directory,"Disaster_Data/USGS/",titlz,"1"),
                               I0=I0),
                   error=function(e) NULL)
  if(is.null(hazsdf)) {
    print(paste0("Error extracting USGS id ",USGSid))
    return(NULL)
  } else if(!check_hazsdf(hazsdf,I0)) {
    print(paste0(max(hazsdf@data$mmi_mean,na.rm = T),
                 ": Either the hazard doesn't exist in USGS, or the magnitude is below 4.5 MMI"))
    return(NULL)
  }
  
  # Extract the date of the event
  sdate<-as.Date(tmp$properties$products$dyfi[[1]]$properties$eventtime)
  
  if(length(sdate)==0){
    sdate <- as.Date(tmp$properties$products$shakemap[[1]]$properties$eventtime)
  }
  
  print(sdate)
  
  return(new("HAZARD",
             obj=hazsdf,
             hazard="EQ",
             dater=as.Date(sdate),
             I0=I0,
             alertlevel=ifelse(is.null(tmp$properties$alert),"green",tmp$properties$alert),
             #alertscore=ifelse(i<=length(alertscores),alertscores[i],0))
             alertscore=NA_real_))
  
}

# WHY IS USGS SO DIFFICULT? GIVE ME A DATE!
GetUSGSdatetime<-function(USGSid){
  url<-paste0("https://earthquake.usgs.gov/fdsnws/event/1/query?eventid=",USGSid,"&format=geojson")
  tmp<-FROM_GeoJson(url)
  
  eventdatetime <- tmp$properties$products$dyfi[[1]]$properties$eventtime
  magnitude<- tmp$properties$products$dyfi[[1]]$properties$magnitude
  max_mmi <- tmp$properties$products$dyfi[[1]]$properties$maxmmi
  depth <- tmp$properties$products$dyfi[[1]]$properties$depth
  lat <- as.numeric(tmp$properties$products$dyfi[[1]]$properties$latitude)
  long <- as.numeric(tmp$properties$products$dyfi[[1]]$properties$longitude)
  
  if(is.null(eventdatetime)){
    eventdatetime <- tmp$properties$products$shakemap[[1]]$properties$eventtime
    magnitude<- tmp$properties$products$shakemap[[1]]$properties$magnitude
    max_mmi <- tmp$properties$products$shakemap[[1]]$properties$maxmmi
    depth <- tmp$properties$products$shakemap[[1]]$properties$depth
    lat <- as.numeric(tmp$properties$products$shakemap[[1]]$properties$latitude)
    long <- as.numeric(tmp$properties$products$shakemap[[1]]$properties$longitude)
  }
  
  
  # Latitude and longitude coordinates
  timezone <- suppressWarnings(tz_lookup_coords(lat = lat, lon = long, method = "fast"))
  
  # Get the timezone
  datetime_local <- as.POSIXct(with_tz(time = as_datetime(eventdatetime), tzone = timezone))
  
  datetime_reformatted <- format(datetime_local,  format = "%Y-%m-%d %H:%M:%S")
  return(datetime_reformatted)
}

GetUSGSfeatures<-function(USGSid){
  url<-paste0("https://earthquake.usgs.gov/fdsnws/event/1/query?eventid=",USGSid,"&format=geojson")
  tmp<-FROM_GeoJson(url)
  
  eventdatetime <- tmp$properties$products$dyfi[[1]]$properties$eventtime
  magnitude<- tmp$properties$products$dyfi[[1]]$properties$magnitude
  max_mmi <- tmp$properties$products$dyfi[[1]]$properties$maxmmi
  depth <- tmp$properties$products$dyfi[[1]]$properties$depth
  lat <- as.numeric(tmp$properties$products$dyfi[[1]]$properties$latitude)
  long <- as.numeric(tmp$properties$products$dyfi[[1]]$properties$longitude)
  alertlevel <- as.numeric(tmp$properties$alert)
  
  if(is.null(eventdatetime)){
    eventdatetime <- tmp$properties$products$shakemap[[1]]$properties$eventtime
    magnitude<- tmp$properties$products$shakemap[[1]]$properties$magnitude
    max_mmi <- tmp$properties$products$shakemap[[1]]$properties$maxmmi
    depth <- tmp$properties$products$shakemap[[1]]$properties$depth
    lat <- as.numeric(tmp$properties$products$shakemap[[1]]$properties$latitude)
    long <- as.numeric(tmp$properties$products$shakemap[[1]]$properties$longitude)
  }
  shakemap_max_mmis <- as.numeric(unlist(sapply(tmp$properties$products$shakemap, function(x) x$properties$maxmmi)))
  max_mmi_diff <- max(shakemap_max_mmis)-shakemap_max_mmis[1]
  
  # Latitude and longitude coordinates
  timezone <- suppressWarnings(tz_lookup_coords(lat = lat, lon = long, method = "fast"))
  
  # Get the timezone
  datetime_local <- as.POSIXct(with_tz(time = as_datetime(eventdatetime), tzone = timezone))
  
  eventtime <- format(datetime_local,  format = "%H:%M:%S")
  eventdate <- as.Date(format(datetime_local,  format = "%Y-%m-%d"))
  
  return(list(eventdate=eventdate, magnitude=as.numeric(magnitude), max_mmi=as.numeric(max_mmi_diff), depth=as.numeric(depth),
              eventtime=eventtime, alertlevel=alertlevel))
}

# Extract EQ data from USGS for a specified event
GetUSGS<-function(USGSid=NULL,bbox,sdate,fdate=NULL,titlz="tmp",I0=4.5,minmag=5){
  
  if(!is.null(USGSid)) {
    hazsdf<-GetUSGS_id(USGSid)
    if(is.null(hazsdf)) return(NULL)
    bbox<-hazsdf@bbox
    USGS<-SearchUSGSbbox(expandBbox(hazsdf@bbox,f = 200,scaling = F),
                         hazsdf@eventdate-7,hazsdf@eventdate+14,minmag)  
    sdate<-fdate<-hazsdf@eventdate
    lenny<-length(USGS$features)
    # Check that the original USGS id shakemap is contained inside the USGS extracted files
    ids<-unlist(sapply(USGS$features,function(x) x$id))
    if(lenny<=1 | !any(ids==USGSid)) {
      lhazdat<-list(hazard_info=list(bbox=bbox,sdate=sdate,fdate=fdate,NumEvents=1,
                                     hazard="EQ",I0=I0,eventdates=sdate),hazsdf)
      # lhazdat<-c(list(bbox=bbox,sdate=sdate,fdate=fdate,
      #               NumEvents=1,hazard="EQ",I0=I0,eventdates=sdate,hazsdf))
      return(lhazdat)
    } else{
      for (i in 1:lenny){
        USGSdate<-GetUSGSdatetime(USGS$features[[i]]$id)
        sdate<-min(sdate,USGSdate)
        fdate<-max(fdate,USGSdate)  
      }
    }
    
  } else {
    # Automatically assign end date if not specified (or badly specified)
    if(is.null(fdate)) {
      fdate=min(Sys.Date(),(as.Date(sdate)+10))
    } else fdate=min(Sys.Date(),as.Date(fdate)+1)
    # Search through the USGS events
    USGS<-SearchUSGSbbox(bbox,sdate,fdate,minmag)  
    lenny<-length(USGS$features)
  }
  
  # if no events are found:
  if(lenny==0) return(NULL)
  
  # lhazdat<-c(list(bbox=bbox,sdate=sdate,fdate=fdate,NumEvents=lenny,hazard="EQ",I0=I0,eventdates=NULL))
  lhazdat<-list(hazard_info=list(bbox=bbox,sdate=sdate,fdate=fdate,NumEvents=lenny,hazard="EQ",I0=I0,eventdates=c(), 
                                 depths=c(), magnitudes=c(), max_mmi=c(), eventtimes=c(), usgs_ids=c(), alertfull=list()))
  tbbox<-rep(NA,4)
  for (i in 1:lenny){
    #FOR ITALY i=43 only
    tmp<-FROM_GeoJson(USGS$features[[i]]$properties$detail)
    
    # Find the details of the raster file for the EQ in the USGS database
    #grid_mmi <- get_intensities_raw(tmp$properties$products$shakemap[[1]]$contents$`download/grid.xml`$url)
    #hazsdf <- ExtractUSGS_xml(grid_mmi, I0=I0)
    
    # Extract EQ raster of hazard intensity
    hazsdf<-tryCatch(ExtractUSGS(url = tmp$properties$products$shakemap[[1]]$contents$`download/raster.zip`$url,
                                 namer = paste0(directory,"Disaster_Data/USGS/",titlz,i),
                                 I0=I0),
                     error=function(e) NULL)
  
    # Check that this extracted event is in the correct form
    if(!check_hazsdf(hazsdf,minmag,bbox)){
      lhazdat$hazard_info$NumEvents<-lhazdat$hazard_info$NumEvents-1
      next
    }
    
    #check to see if there are any identical events stored in USGS and remove if so
    if (length(lhazdat) > 1){
      duplicate <- F
      for (j in 2:length(lhazdat)){
        if(length(hazsdf$mmi_mean)==length(lhazdat[[j]]$mean)){
          if(all(hazsdf$mmi_mean==lhazdat[[j]]$mean)) duplicate <- T
        }
      }
      if (duplicate) next
    }
    
    PAGER_alert <- GetPagerFatality(tmp$properties$products$losspager[[1]]$contents$pager.xml$url)
    
    # Create HAZARD object
    event_features <- GetUSGSfeatures(USGS$features[[i]]$id)
    hazsdf<-new("HAZARD",
                obj=hazsdf,
                hazard="EQ",
                dater=event_features$eventdate,
                I0=I0,
                alertlevel=tmp$properties$alert, # ifelse(is.null(tmp$properties$alert),"green",tmp$properties$alert),
                #alertscore=ifelse(i<=length(alertscores),alertscores[i],0))
                alertscore=0, 
                depth=event_features$depth, 
                magnitude=event_features$magnitude, 
                max_mmi=event_features$max_mmi,
                eventtime=event_features$eventtime,
                USGS_id=USGS$features[[i]]$id, 
                alertfull=PAGER_alert
    )
    # Add to the list of hazards
    lhazdat[[length(lhazdat)+1]]<-hazsdf
    # Extend the bounding box to account for this earthquake
    tbbox[c(1,2)]<-apply(cbind(tbbox[c(1,2)],hazsdf@bbox[c(1,2)]),1,min,na.rm=T)
    tbbox[c(3,4)]<-apply(cbind(tbbox[c(3,4)],hazsdf@bbox[c(3,4)]),1,max,na.rm=T)
    # Extract dates for each hazard event
    lhazdat$hazard_info$alertlevel%<>%c(hazsdf@alertlevel)
    lhazdat$hazard_info$eventdates%<>%c(as.character(hazsdf@eventdate))
    lhazdat$hazard_info$depths%<>%c(hazsdf@depth)
    lhazdat$hazard_info$magnitudes%<>%c(hazsdf@magnitude)
    lhazdat$hazard_info$max_mmi%<>%c(hazsdf@max_mmi)
    lhazdat$hazard_info$eventtimes%<>%c(hazsdf@eventtime)
    lhazdat$hazard_info$usgs_ids%<>%c(hazsdf@USGS_id)
    lhazdat$hazard_info$alertfull[[length(lhazdat$hazard_info$alertfull)+1]] <- hazsdf@alertfull
  }
  if(any(is.na(tbbox))) return(NULL)
  # Modify the bounding box to fit around all hazards extracted
  lhazdat$hazard_info$bbox<-tbbox
  lhazdat$hazard_info$eventdates%<>%as.Date()
  lhazdat$hazard_info$fdate<-max(lhazdat$hazard_info$eventdates)
  return(lhazdat)
  
}

GetPagerFatality <- function(url){
  if (is.null(url)) return(list())
  # Read the XML file from the URL
  xml_data <- read_xml(url)
  
  # Extract information from the XML
  alert_type <- xml_text(xml_find_first(xml_data, "//alert[@type='fatality']/@type"))
  alert_level <- xml_text(xml_find_first(xml_data, "//alert[@type='fatality']/@level"))
  summary <- xml_text(xml_find_first(xml_data, "//alert[@type='fatality']/@summary"))
  units <- xml_text(xml_find_first(xml_data, "//alert[@type='fatality']/@units"))
  
  # Extract bin information
  bins <- xml_find_all(xml_data, "//alert[@type='fatality']/bin")
  bin_list <- lapply(bins, function(bin) {
    min_val <- as.numeric(xml_attr(bin, "min"))
    max_val <- as.numeric(xml_attr(bin, "max"))
    probability <- as.numeric(xml_attr(bin, "probability"))
    color <- xml_attr(bin, "color")
    
    return(list(min = min_val, max = max_val, probability = probability, color = color))
  })
  
  # Create the list
  pager_alert_fatalities <- list(
    alert_type = alert_type,
    alert_level = alert_level,
    summary = summary,
    units = units,
    bins = bin_list
  )
  
  # Print the result list
  print(pager_alert_fatalities)
}

library(xml2)
get_intensities_raw <- function(url){
  # Read the XML data from the URL
  xml_file <- read_xml(url)
  grid <- xmlParse(xml_file)
  
  xml_data <- xmlToList(grid) 
  lines <- strsplit(xml_data[[20]], "\n")[[1]]
  
  intensities <- sapply(lines, function(x) as.numeric(strsplit(x, " ")[[1]][5]))
  longitude <- sapply(lines, function(x) as.numeric(strsplit(x, " ")[[1]][1]))
  latitude <- sapply(lines, function(x) as.numeric(strsplit(x, " ")[[1]][2]))
  plot_df <- data.frame(longitude=longitude, latitude=latitude, intensities=intensities)
  
  #ggplot(plot_df, aes(x=longitude, y=latitude, color=intensities)) + geom_point() +
  #  scale_color_gradient(low = "green", high = "red")
  
  return(plot_df[-1,])
}

ExtractUSGS_xml<-function(url,namer,I0=NULL,plotty=F){
  
  wide_df <- spread(grid_mmi, key = longitude, value = intensities)
  mean_haz <- raster(
    matrix(as.numeric(grid_mmi$intensities), nrow = length(unique(grid_mmi$latitude)), byrow=T),
    xmn = min(as.numeric(grid_mmi$longitude)),
    xmx = max(as.numeric(grid_mmi$longitude)),
    ymn = min(as.numeric(grid_mmi$latitude)),
    ymx = max(as.numeric(grid_mmi$latitude)),
    crs=crs(meanhaz)
  )
  names(mean_haz) <- 'mmi_mean'
  mean_haz <- resample(mean_haz, meanhaz, method='bilinear')
  sd_haz = mean_haz
  names(sd_haz) <- 'mmi_std'
  values(sd_haz) <- 0.1
  # Form a standard USGS object
  sgdf<-formUSGSobject(mean_haz,sd_haz,I0)
  
  return(sgdf)
  
}


