#--------------------------------------------------------------
#----------------------DEMOGRAPHIC DATA------------------------
#--------------------------------------------------------------

WorldPopURL<-function(iso3c,year,constrained=T,kmres=T,unadj=T,BGSM=T){
  
  #------------------CONSTRAINED DATASETS-------------------------------
  if(constrained & kmres & unadj & BGSM){
    base_url<-"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/"
    url<-paste0(base_url,year,"/BSGM/",str_to_upper(iso3c),
                "/",str_to_lower(iso3c),"_ppp_",year,"_UNadj_constrained.tif")  
  } else if(constrained & kmres & unadj & !BGSM){
    base_url<-"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/"
    url<-paste0(base_url,year,"/maxar_v1/",str_to_upper(iso3c),
                "/",str_to_lower(iso3c),"_ppp_",year,"_UNadj_constrained.tif")  
  } else if(constrained & kmres & !unadj){
    base_url<-"https://data.worldpop.org/GIS/Population/Individual_countries/"
    url<-paste0(base_url,year,"/",str_to_upper(iso3c),
                "/",str_to_lower(iso3c),"_ppp_",year,"_constrained.tif")  
  } else if(constrained & !kmres & !unadj){
    warning("100m WorldPop constrained dataset only available for 2010 and 2015")
    stop("please go and download by hand from https://data.worldpop.org/GIS/Population/Individual_countries/")
    year<-c(2010,2015)[which.min(abs(year-c(2010,2015)))]
    base_url<-"https://data.worldpop.org/GIS/Population/Individual_countries/"
    url<-paste0(base_url,str_to_upper(iso3c),"/",
                paste0(convIso3Country(iso3c), collapse = "_"),"_100m_Population/",
                paste0(convIso3Country(iso3c), collapse = "_"),"_100m_Population",".tif")  
  } else if(constrained & !kmres & !unadj){
    warning("100m WorldPop constrained dataset only available for 2010 and 2015")
    stop("please go and download by hand from https://data.worldpop.org/GIS/Population/Individual_countries/")
    year<-c(2010,2015)[which.min(abs(year-c(2010,2015)))]
    base_url<-"https://data.worldpop.org/GIS/Population/Individual_countries/"
    url<-paste0(base_url,str_to_upper(iso3c),
                "/",str_to_lower(iso3c),"_ppp_",year,"_constrained.tif")  
    
    #---------------------UNCONSTRAINED DATASETS-----------------------------
  } else if(!constrained & !kmres & unadj){    
    # UN-adjusted, 100m resolution and using the unconstrained method
    base_url<-"https://data.worldpop.org/GIS/Population/Global_2000_2020/"
    url<-paste0(base_url,year,"/",str_to_upper(iso3c),
                "/",str_to_lower(iso3c),"_ppp_",year,"_UNadj.tif")  
  } else if(!constrained & kmres & unadj){    
    # UN-adjusted, 1 km resolution and using the unconstrained method
    if (year <= 2020){
      base_url<-"https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/"
      url<-paste0(base_url,year,"/",str_to_upper(iso3c),
                  "/",str_to_lower(iso3c),"_ppp_",year,"_1km_Aggregated_UNadj.tif")
    } else {
      base_url<-"https://data.worldpop.org/GIS/Population/Global_2021_2022_1km_UNadj/unconstrained/"
      url<-paste0(base_url,year,"/",str_to_upper(iso3c),
                  "/",str_to_lower(iso3c),"_ppp_",year,"_1km_UNadj.tif")
    }

    #url<-paste0(base_url,year,"/",str_to_upper(iso3c),
    #            "/",str_to_lower(iso3c),"_ppp_",year,"_1km_Aggregated_UNadj.tif")
  } else if(!constrained & !kmres & !unadj){    
    # 100m resolution and using the unconstrained method
    base_url<-"https://data.worldpop.org/GIS/Population/Global_2000_2020/"
    url<-paste0(base_url,year,"/",str_to_upper(iso3c),
                "/",str_to_lower(iso3c),"_ppp_",year,".tif")  
  } else if(!constrained & kmres & !unadj){ 
    # 1 km resolution and using the unconstrained method
    base_url<-"https://data.worldpop.org/GIS/Population/Global_2000_2020_1km/"
    url<-paste0(base_url,year,"/",str_to_upper(iso3c),
                "/",str_to_lower(iso3c),"_ppp_",year,"_1km_Aggregated.tif")
  }
  return(url)
}

# Extract the required URL link to WorldPop and build the file names
CheckWPop_API<-function(iso3c,year,folder="./",constrained=T,kmres=T,unadj=T){
  # Main host URL
  url<-WorldPopURL(iso3c,year,constrained,kmres,unadj)
  # File name to be saved
  filer<-paste0("WorldPop_Population_UNadj_",ifelse(constrained,"constrained_",""),str_to_upper(iso3c),"_",year,".tif")
  locy<-paste0(folder,filer)
  # Attempt to extract file
  checker<-tryCatch(download.file(url,locy),error = function(e) NA)
  if(is.na(checker)) {
    print("No BGSM data for WorldPop, trying for Maxar datasets")
    # If error, try the maxar-based population data
    url<-WorldPopURL(iso3c,year,constrained,kmres,unadj,BGSM = F)
    checker<-tryCatch(download.file(url,locy, mode="wb"),error = function(e) NA)
    if(is.na(checker)) stop("Worldpop file not found, please check the WorldPop API address and your input iso3c & year")
    
  } 
  
  return(list(filer=filer,locy=locy))
}

#---------------------------------WORLDPOP POPULATION-----------------------------------------#
GetWorldPopISO3C<-function(iso3c,year=NULL,folder="./Data/Exposure/PopDemo/",constrained=T,kmres=T,unadj=T){
  # Try to download the most recent dataset
  if(is.null(year)) {year<-AsYear(Sys.Date()); mostrecent<-T} else mostrecent<-F
  # If we left the year blank, then let's search for the most recent CONSTRAINED dataset
  if(mostrecent){
    # Go through every year from now until 2020 until we find some data!
    yr<-year; extracter<-T
    # keep trying to extract most recent data, unless you go lower than 2020 which means WorldPop changed their API addresses
    while(extracter | yr>=2020){
      checker<-tryCatch(CheckWPop_API(iso3c,yr,folder,constrained,kmres,unadj),error = function(e) NA)
      # Check that something was returned!
      if(is.na(checker)) {yr<-yr-1; next} else extracter<-F
    }
  } else {
    checker<-tryCatch(CheckWPop_API(iso3c,year,folder,constrained,kmres,unadj),error = function(e) NA)
  }
  # If nothing happened... FML
  if(any(is.na(checker))) stop("something went wrong in extracting WorldPop data, check year, iso3c & output folder address")
  # Extract the mean hazard intensity from raster
  popy<-terra::rast(file.path(checker$locy))
  # Convert it to how we looove it!
  #popy%>%as("SpatialPixelsDataFrame")
}

#------------------------------------------------
#---------------EXTRACTION-----------------------
#------------------------------------------------

ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_All_2023May19/ODDobjects/EQ20170614GTM_64')

getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mode_non_na <- function(arr,...){
  non_zero <- which(arr!=0)
  if(length(non_zero)==0){ return(0)}
  else{return(getmode(arr[non_zero]))}
}

getWorldPop_ODD <- function(dir, year, bbox_vect, folder='Demography_Data/Population/WorldPop/'){
  
  bbox <- rbind(c(bbox_vect[1], bbox_vect[3]), c(bbox_vect[2],bbox_vect[4]))
  
  nations <- GetNationsBbox(dir, bbox)
  iso3c_all <- unique(nations$ISO3C)[which(!is.na(unique(nations$ISO3C)))]
  
  popy <- GetWorldPopISO3C(iso3c_all[1], year=year, constrained=F, folder=paste0(dir, folder)) %>% raster()
  popy_cropped <- crop(popy, extent(bbox), snap='out')
  if (length(iso3c_all)>1){
    for (iso3c in iso3c_all[2:length(iso3c_all)]){
      popy_add <- GetWorldPopISO3C(iso3c, year=year, constrained=F, folder=paste0(dir, folder)) %>% raster()
      popy_add_cropped <- crop(popy_add, extent(bbox), snap='out')
      popy_cropped %<>% merge(popy_add_cropped)
    }
  }
  spat_agg <- aggregate(popy_cropped, fact=2, fun=sum, expand=F) 
  
  names(spat_agg) <- 'Population'
  spat_agg$dummy <- 0 #need to add some variable so that pixels with 0 population are not lost when converting to SpatialPixelsDataFrame
  spdf <- as(spat_agg, 'SpatialPixelsDataFrame')
  spdf@data <- spdf@data[,-2, drop=F] #remove dummy
  
  iso3_lookup <- data.frame(ISO3C=c(NA, unique(iso3c_all)), id=0:length(iso3c_all))
  nations@data$order <- 1:NROW(nations@data)
  nations@data %<>% merge(iso3_lookup)
  nations@data <- nations@data[order(nations@data$order),]
    
  rastered_iso3 <- rasterize(nations, spat_agg, field='id', fun=mode_non_na)
  spdf_iso3 <- as(rastered_iso3, 'SpatialPixelsDataFrame')
  colnames(spdf_iso3@data) <- 'id'
  spdf_iso3@data$order <- 1:NROW(spdf_iso3@data)
  spdf_iso3@data %<>% merge(iso3_lookup, all.x=T)
  spdf_iso3@data <- spdf_iso3@data[order(spdf_iso3@data$order),'ISO3C', drop=F]
  
  spdf@data <- merge_rastered_spdf(spdf, spdf_iso3, 'ISO3C') 
 
  colnames(spdf@coords) <- c('Longitude', 'Latitude')
  rownames(spdf@bbox) <- c('Longitude', 'Latitude')
  return(spdf)
}

#----------------------GRIDDED WORLDPOP----------------------------
# Extract population data and put the conflict data onto the grid
GetPop<-function(ADM,ISO,constrained=T,ncores=2,outsiders=T,ext){
  # Get the WorldPop data
  if(length(list.files("./Data/Exposure/PopDemo/",paste0(ISO,"_")))==0){
    pop<-GetWorldPopISO3C(ISO,2020,folder="./Data/Exposure/PopDemo/",constrained=T,kmres=T,unadj=T)
    names(pop)[1]<-"POPULATION"
  } else {
    pop<-terra::rast(paste0("./Data/Exposure/PopDemo/",
                            list.files("./Data/Exposure/PopDemo/",paste0(ISO,"_"))))
    #%>%
    # as("SpatialPixelsDataFrame")
    names(pop)[1]<-"POPULATION"
  }
  # # Aggregate the population data to admin level 2
  # popvec<-Grid2ADM(pop,
  #                  ADM[ADM@data$ISO3CD==ISO,],
  #                  sumFn="sum",
  #                  index = which(names(pop)=="POPULATION"),
  #                  ncores = ncores,
  #                  outsiders=outsiders)
  # # Scale to make sure that the value is current
  # factor<-tryCatch(InterpPopWB(ISO,Sys.Date(),normdate=as.Date("2015-01-01"))$factor,error=function(e) NA)
  # popvec$all<-popvec$all*ifelse(is.na(factor),1,factor)
  # # Combine into one large data.frame
  # ADM@data%<>%cbind(data.frame(POPULATION=round(popvec$all)))
  
  #VErsion2: Aggregate population values per polygon:
  # ADM$Population <-Grid2ADM_v2(grid=pop,poly=ADM,Fn= "sum") %>%
  #   round(., 0) #make whole number
  pop%<>%terra::crop(ext)
  y<-pop%>%
    terra::extract(ADM,method='bilinear',fun=sum,na.rm=T,ID=FALSE)%>%
    unlist()%>%
    round(., 0)%>% #make whole number
    ifelse(is.nan(.) ==TRUE,NA,.)  #NaN to Na
  
  
  ADM<-cbind(ADM,y)
  names(ADM)[names(ADM) == "y"] <- "Population"
  
  return(ADM)
}



GetDemog<-function(Pop_totl,ADM,ISO){
  # Get World Bank data: Fem_prop, Und14,Ovr64
  indic<- c("SP.POP.TOTL.FE.ZS","SP.POP.0014.TO.ZS","SP.POP.65UP.TO.ZS")
  for(i in 1:4){ #becuase download always failing! repeat them!
    Pop_stats<-tryCatch(sapply(indic, function(x) WBcall(syear=AsYear(Sys.Date())-1 ,
                                                         indicator=x,ISO=ISO)[["Value"]]/100 * Pop_totl),
                        error=function(e) NA)
  }
  
  if(is.matrix(Pop_stats)==TRUE){
    Pop_stats %<>%
      as.data.frame()
  }else{
    Pop_stats %<>%
      dplyr::bind_rows() %>%
      as.data.frame()
  }
  
  colnames(Pop_stats) <-c("Fe_prop","Und14","Ovr64")
  
  ADM<-cbind(ADM,Pop_stats)
  return(ADM)
}



#Demographics from WorldPop
##----------------Get AgeSex_structure - Constrained UN adjusted data 1km.--------------- 
#highly simplified code: update later with more function parameters for the urls------
library(RCurl)
library(XML)

#Download data
GetWAgSx_L2<-function(iso,year,folder){
  base<-paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/constrained/",year)
  #From url to local dir
  iso_urls <-paste0(base,"/five_year_age_groups/",toupper(iso),"/")%>%
    getURL(.,verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>%
    getHTMLLinks(.,xpQuery = "//a/@href[contains(., '.tif')]")
  
  
  
  if(length(list.files(folder))!=42){
    print("Age_Sex structure .tif images not complete - redownloading from source")
    #clean folder then download:
    f <- list.files(folder, full.names = T)
    #remove the files
    file.remove(f)
    lapply(iso_urls, function(x) download.file(url=paste0(base,"/five_year_age_groups/",toupper(iso),"/",x),
                                               destfile = paste0(folder,"/",x), mode="wb"))
  }else{
    print("subdir has the images. no longer downloading from source")
  }
  
}  

#extract values to admin_polygons
GetAgeSexStruc<-function(ISO,ADM,year,path_to_files,create_subfolder=T){
  #check path dir
  if(create_subfolder == TRUE){
    subdir<-paste0(path_to_files,"/",ISO,"/")
    if(dir.exists(subdir)==TRUE){
      print("Subdirectoy exists, not creating a new one")
    }else{
      dir.create(subdir)
    }
  }else{
    if(dir.exists(subdir)==FALSE){
      print("Subdirectory does not exist, change create_subfolder to TRUE")
    }
  }
  #Download data is needed:
  GetWAgSx_L2(iso=ISO, year=year, folder=subdir)
  #Read
  fil<-list.files(path=subdir,
                  pattern = paste0(".*",tolower(ISO),".*\\.tif$"),full.names = TRUE)
  imgs<-terra::rast(fil)
  #extract values per ADMl2 unit
  # sx<-substr(names(imgs),5,5)
  # age<-sapply(names(imgs), function(x) {as.numeric(strsplit(x, "\\D+")[[1]][-1])%>%
  #     .[[1]]}) #first number is age
  # 
  # cat<-paste0(toupper(sx),age,"pop_UNAdj")
  #colnames:
  cat<-str_extract(names(imgs), paste0("(?<=",tolower(ISO),"_).*?(?=_",year,")"))
  
  imgs %<>%terra::crop(ADM)%>%
    terra::mask(ADM)
  
  AgSx<-imgs%>%terra::extract(ADM,method='bilinear',na.rm=T,fun=sum,ID=FALSE)%>%
    round(.,0)%>%
    as.data.frame()
  
  colnames(AgSx)<-cat
  
  # ADM<-cbind(ADM,AgSx)
  return(AgSx)
  
}