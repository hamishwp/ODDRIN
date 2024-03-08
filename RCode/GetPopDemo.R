# library(tiff)
# library(pracma)
# library(ggmap)
# library(OpenStreetMap)
# library(osmdata)
# library(ggplot2)
# library(geosphere)
# library(reshape2)
# library(tidyverse)
# library(maps)
# library(magrittr)

# bbox<-c(110.18,-44.93,115.01,-10.17)
# bbox<-opq("Australia")$bbox
# bbox<-as.double(unlist(strsplit(bbox,",")));bbox<-c(bbox[2],bbox[3],bbox[4],bbox[1])
# bbox<-c(110,-45,160,-10)

# directory<-"/home/patten/Documents/Coding/IIDIPUS/"
# bbox<-opq("Melbourne")$bbox
# bbox<-as.double(unlist(strsplit(bbox,",")));bbox<-c(bbox[2],bbox[3],bbox[4],bbox[1])
# namer<-"Melbourne_Australia_Population"
# country<-"Australia"
# 
# # London:
# # bbox - [mnlo,mnla,mxlo,mxla]
# bbox<-c(-1,51,1,52)
# country<-"UK"
# namer<-"London_UK"
# 
# bbox<-opq("East Gippsland")$bbox
# bbox<-as.double(unlist(strsplit(bbox,",")));bbox<-c(bbox[2],bbox[3],bbox[4],bbox[1])

# directory<-"/home/patten/Documents/Coding/Oxford/IIDIPUS/"
# bbox<-c(135,-45,155,-25)
# country<-"Australia"
# namer<-"NSW_Australia_Population"

###############################################################################################
# SortDemoData: filters out unwanted population data outside of bounding box provided by user
###############################################################################################
SortDemoData<-function(filer,bbox=NULL){
  
  info<-readLines(filer,n = 6); info<-strsplit(info, "\\s+")
  for (i in 1:6){ assign(info[[i]][1],as.numeric(info[[i]][2]))}; rm(info)
  
  #seems to be numerical discrepancies between nation and population data from SEDACs (on a very, very small order)
  #correct by rounding xllcorner, yllcorner, and cellsize to the nearest 1/10000th of an arcsecond
  one_tenthousandths_arcsecond <- 1/60/60/10000
  xllcorner <- round(xllcorner/one_tenthousandths_arcsecond)*one_tenthousandths_arcsecond
  yllcorner <- round(yllcorner/one_tenthousandths_arcsecond)*one_tenthousandths_arcsecond
  cellsize <- round(cellsize/one_tenthousandths_arcsecond)*one_tenthousandths_arcsecond
  
  popdemo<-read.csv(filer,header = FALSE,skip = 6,sep = " ",na.strings = NODATA_value,colClasses = "numeric")
  # dimensions of popdemo: [decreasing(latitude),longitude]
  sizer<-dim(popdemo)
  sizer2<-c(nrows,ncols+1)
  if(!all(length(sizer)==length(sizer2)) || !all(sizer==sizer2)){stop(paste0("ERROR! Incorrect dimensions: Check the population demography file for bounding box ",bbox," in file GetPopDemo.R"))}
  
  lat<-rev(seq(from=yllcorner+cellsize/2,by=cellsize, length.out=nrows))
  #lat<-seq(from=yllcorner+90-cellsize/2,by=-cellsize,length.out = nrows) #@@@ LATITUDE @@@# This doesn't work when we're not working with exact quadrants
  long<-seq(from=xllcorner+cellsize/2,by=cellsize,length.out = ncols)    #@@@ LONGITUDE @@@#
  colnames(popdemo)<-long
  row.names(popdemo)<-lat
  
  if(is.null(bbox)) return(popdemo %>% as.matrix() %>% pracma::rot90(-1))
  
  imnlo<-which.min(abs(bbox[1]-long))
  imxlo<-which.min(abs(bbox[3]-long))
  imnla<-which.min(abs(bbox[2]-lat))
  imxla<-which.min(abs(bbox[4]-lat))

  popdemo<-popdemo[imxla:imnla,imnlo:imxlo] %>% as.matrix() %>% pracma::rot90(-1)
  
  return(popdemo)
}

########################################################################################
# GetSEDACfnum: Get SEDAC file number that contains the demography data in the bounding box
########################################################################################
GetSEDACfnum<-function(long,lat){
  
  longCo<-c(-90L,0L,90L,180L)
  Llo<-longCo-long
  
  if(any(Llo< -270L)||any(Llo>360L)||(lat>90L)||(lat< -90L)){stop("Error in the longitude and latitude values for SEDACS data: see GetPopDemo.R\n")}
  
  SEDAC<-which.min(abs(Llo))
  llg<-longCo[SEDAC]
  # mLlo values have to be less than the box.
  if(Llo[SEDAC]<0L){SEDAC<-SEDAC+1L}
  if(lat<0L){SEDAC<-SEDAC+4L}
  
  return(as.integer(c(SEDAC,llg)))
}

######################################################################################
# Provided the bounding box, this function finds, filters and returns 
# the SEDAC population/demography data for that region.
######################################################################################
ExtractSEDACS<-function(strings, bbox){
  
  directory<-strings[1]
  loc<-strings[2]
  nom<-strings[3]
  
  LL<-GetSEDACfnum(bbox[1],bbox[2]) # Lower Left
  LR<-GetSEDACfnum(bbox[3],bbox[2]) # Lower right
  UL<-GetSEDACfnum(bbox[1],bbox[4]) # Upper Left
  UR<-GetSEDACfnum(bbox[3],bbox[4]) # Upper right
  
  if(abs(LL[1]-LR[1])>1L) {
    print("Bounding box of SEDAC population/demography data is too large, using lower resolution")
    stop("Haven't modified code to run lowres yet")
  }
  
  SEDAC  <-unique(c(UL[1],UR[1],LL[1],LR[1]))
  uniquer<- length(SEDAC)
  if(uniquer==1L){
    
    filer<-paste0(directory,loc,nom,SEDAC,".asc")
    population<-SortDemoData(filer,bbox)
    
  } else if (uniquer==2L){
    
    if(((SEDAC[1]<5L)&&(SEDAC[2]<5L))|((SEDAC[1]>4L)&&(SEDAC[2]>4L))){
      # LEFT AND RIGHT SEDACS POPULATION DATA PANELS
      # bbox - [mnlo,mnla,mxlo,mxla]
      bb1<-c(bbox[1],bbox[2],LL[2],bbox[4]) # LEFT PANEL
      bb2<-c(LL[2],bbox[2],bbox[3],bbox[4]) # RIGHT PANEL
      functionner<-"rbind"
    } else if (SEDAC[1]==(SEDAC[2]+4L)|SEDAC[1]==(SEDAC[2]-4L)) {
      # UPPER AND LOWER SEDACS POPULATION DATA PANELS
      # bbox - [mnlo,mnla,mxlo,mxla]
      bb1<-c(bbox[1],0,bbox[3],bbox[4]) # UPPER PANEL
      bb2<-c(bbox[1],bbox[2],bbox[3],0) # LOWER PANEL
      functionner<-"cbind"
    } else {stop("Something went wrong.... check SEDAC file ordering, see GetPopDemo.R\n")}
    
    filer<-paste0(directory,loc,nom,SEDAC[1],".asc")
    pop1<-SortDemoData(filer,bb1)
    filer<-paste0(directory,loc,nom,SEDAC[2],".asc")
    pop2<-SortDemoData(filer,bb2)    
    
    bindy<-match.fun(functionner)
    population<-bindy(pop1,pop2)
    
    rm(pop1,pop2)
    
  } else {
    
    print("Population/Demography bounding box is in 4 separate quadrants... please be patient while this loads. Bejinhos")
    # ALL BOUNDING BOX CORNERS ARE IN DIFFERENT QUADRANTS OF THE SEDAC FILE
    # bbox - [mnlo,mnla,mxlo,mxla]
    bb1 <- c(bbox[1],  0,        LL[2],    bbox[4])   # Upper Left
    bb2 <- c(LL[2],    0,        bbox[3],  bbox[4])   # Upper Right
    bb3 <- c(bbox[1],  bbox[2],  LL[2],    0      )   # Lower Left
    bb4 <- c(LL[2],    bbox[2],  bbox[3],  0      )   # Lower Right
    
    filer<-paste0(directory,loc,nom,UL[1],".asc")
    pop1<-SortDemoData(filer,bb1)
    filer<-paste0(directory,loc,nom,UR[1],".asc")
    pop2<-SortDemoData(filer,bb2)
    filer<-paste0(directory,loc,nom,LL[1],".asc")
    pop3<-SortDemoData(filer,bb3)
    filer<-paste0(directory,loc,nom,LR[1],".asc")
    pop4<-SortDemoData(filer,bb4)
    
    # LEFT AND RIGHT USE rbind, UPPER AND LOWER USE cbind
    population<-cbind(rbind(pop1,pop2),rbind(pop3,pop4))
    
    rm(pop1,pop2,pop3,pop4)
    
  }
  
  return(population)
  
}


# GetSEDACArea<-function(directory){
#   
#   cfiler<-paste0(directory,"Demography_Data/Population/gpw-v4-population-count-2015/",
#                  "gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_")
#   dfiler<-paste0(directory,"Demography_Data/Population/gpw-v4-population-density-2015/",
#                  "gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_")
#   saver<-paste0(directory,"Demography_Data/Population/gpw_v4_area_")
#   for (i in 1:8){
#     count<-SortDemoData(paste0(cfiler,i,".asc"))
#     dens<-SortDemoData(paste0(dfiler,i,".asc"))
#     area<-count/dens
#     save(area,file = paste0(saver,i,".Rdata"))
#     rm(count,dens,area)
#   }
#   
# }

######################################################################################
# Extracts the SEDAC population data for the bounding box and can also plot it
######################################################################################
GetPopulationBbox<-function(directory,bbox,density=F,lowres=FALSE,yr="2015",plotty=FALSE, namer="Population",ncity=1){
  # bbox is bounding box in the form 'min lon, max lat, max lon, min lat'
  # DATA: NASA - SEDAC
  
  if(density){ctds<-"density"} else {ctds<-"count"}
  
  if(abs(bbox[2])>90 | abs(bbox[4])>90 | abs(bbox[1])>180 | abs(bbox[3])>180) {stop("Error: non-physical bounding box values in GetPopulationBbox")}
  
  yr%<>%as.numeric()
  if(is.na(yr)) ("Please provide year as either numeric (2015) or character value ('2015')")
  
  # if(!is.null(date)) {
  #   date%<>%try(as.Date,silent=T)
  #   if(class(date)=="try-error") {
  #     print("WARNING: date badly specified in GetPopulationBbox (should be in format '2015-01-25'), no interpolation will be performed")
  #     date<-NULL
  #   }
  # }
  
  #   GetSEDACfnum(LONG,   LAT)
  LL<-GetSEDACfnum(bbox[1],bbox[2]) # Lower Left
  LR<-GetSEDACfnum(bbox[3],bbox[2]) # Lower right
  UL<-GetSEDACfnum(bbox[1],bbox[4]) # Upper Left
  UR<-GetSEDACfnum(bbox[3],bbox[4]) # Upper right
  
  # if(abs(LL[1]-LR[1])>1L) {
  #   print("Bounding box of SEDAC population data is too large, using lower resolution")
  #   lowres<-TRUE
  # }
  
  if (lowres){
    print("WARNING: using 2020 not 2015 lowres pop data ")
    filer<-paste0(directory,"Demography_Data/Population/",
                  "gpw_v4_population_",ctds,"_adjusted_to_2015_unwpp_country_totals_rev11_",yr,"_2pt5_min.asc")
    population<-SortDemoData(filer,bbox)
    
  } else {
    
    poploc<-paste0("Demography_Data/Population/gpw-v4-population-",ctds,"-",yr,"/")
    popnom<-paste0("gpw_v4_population_",ctds,"_adjusted_to_2015_unwpp_country_totals_rev11_",yr,"_30_sec_")
   
    strings<-c(directory,poploc,popnom)
    
    natloc<-paste0("Demography_Data/Population/gpw-v4-national-identifier-grid-rev11_30_sec_asc/")
    natnom<-paste0("gpw_v4_national_identifier_grid_rev11_30_sec_")
    strings_nat <- c(directory, natloc, natnom)
    
    population<-ExtractSEDACS(strings, bbox) 
    nations<-ExtractSEDACS(strings_nat,  bbox)
    nation_names_lookup <- read.csv(paste0(dir, natloc,'gpw_v4_national_identifier_grid_rev11_lookup.txt'), sep='\t')
    
  }
  
  if(plotty){
    
    longData<-melt(population)
    longData<-longData[longData$value!=0,]
    
    cities<-maps::world.cities%>%filter(lat>bbox[2]&lat<bbox[4]&long>bbox[1]&long<bbox[3])%>%arrange(desc(pop))
    if(ncity>1){wordcloud::wordcloud(words=cities$name,freq = cities$pop,max.words = 30,scale = c(2.5,0.2))}
    cities<-slice(cities,1:ncity)
    
    p<-ggplot(longData, aes(x = Var1, y = Var2)) + 
      geom_raster(aes(fill=log10(value))) + 
      scale_fill_gradient(low="yellow", high="red",name = "log_{10}(Population Count)") +
      labs(x="Longitude", y="Latitude" )+#, title="East Gippsland, Australia") +
      theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11)) +
      geom_label(data = cities, aes(long, lat, label = name), size = 4, fontface = "bold", nudge_x = 0.25,nudge_y = -0.25)
    print(p)
    ggsave(paste0(namer,"Population_",ctds,".eps"), plot=p,path = paste0(directory,'/'),width = 9,height = 7.)
    
    # mad_map <- get_map(bbox,maptype = "toner-background")
    # print(ggmap(mad_map))
    
  } 
  
  population %<>% convMat2SPDF(name="Population")
  nations %<>% convMat2SPDF(name='ISO_id')
  nations %<>% merge(nation_names_lookup %>% dplyr::select(Value, ISO3C = ISOCODE), by.x='ISO_id', by.y='Value', all.x=T, sort=F)
  nations$ISO_id <- NULL
  pop_with_iso3 <- cbind(population, nations)
  
  return(pop_with_iso3)
}

GetNationsBbox<-function(directory,bbox){
  # bbox is bounding box in the form 'min lon, max lat, max lon, min lat'
  # DATA: NASA - SEDAC
  
  if(abs(bbox[2])>90 | abs(bbox[4])>90 | abs(bbox[1])>180 | abs(bbox[3])>180) {stop("Error: non-physical bounding box values in GetPopulationBbox")}
  
  # if(!is.null(date)) {
  #   date%<>%try(as.Date,silent=T)
  #   if(class(date)=="try-error") {
  #     print("WARNING: date badly specified in GetPopulationBbox (should be in format '2015-01-25'), no interpolation will be performed")
  #     date<-NULL
  #   }
  # }
  
  #   GetSEDACfnum(LONG,   LAT)
  LL<-GetSEDACfnum(bbox[1],bbox[2]) # Lower Left
  LR<-GetSEDACfnum(bbox[3],bbox[2]) # Lower right
  UL<-GetSEDACfnum(bbox[1],bbox[4]) # Upper Left
  UR<-GetSEDACfnum(bbox[3],bbox[4]) # Upper right
  
  # if(abs(LL[1]-LR[1])>1L) {
  #   print("Bounding box of SEDAC population data is too large, using lower resolution")
  #   lowres<-TRUE
  # }
  
  natloc<-paste0("Demography_Data/Population/gpw-v4-national-identifier-grid-rev11_30_sec_asc/")
  natnom<-paste0("gpw_v4_national_identifier_grid_rev11_30_sec_")
  strings_nat <- c(directory, natloc, natnom)
  
  nations<-ExtractSEDACS(strings_nat,  bbox)
  nation_names_lookup <- read.csv(paste0(dir, natloc,'gpw_v4_national_identifier_grid_rev11_lookup.txt'), sep='\t')
  
  nations %<>% convMat2SPDF(name='ISO_id')
  nations %<>% merge(nation_names_lookup %>% dplyr::select(Value, ISO3C = ISOCODE), by.x='ISO_id', by.y='Value', all.x=T, sort=F)
  nations$ISO_id <- NULL
  
  return(nations)
}

# GetPopulationBboxSAFE<-function(bbox){
#   
#   population<-raster(x=paste0(directory,
#                               "Demography_Data/Population/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min.tif"))
#   
#   e <- as(raster::extent(c(bbox[c(1,3,2,4)])), 'SpatialPolygons')
#   crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
#   population%<>%raster::crop(e)
#   
#   rownames(population@bbox)<-c("Longitude","Latitude")
#   colnames(population@coords)<-c("Longitude","Latitude")
#   names(population)<-"Population"
#   
#   xo<-seq(from=bbox[1],to=bbox[3],by=population@grid@cellsize[1])
#   yo<-seq(from=bbox[2],to=bbox[4],by=population@grid@cellsize[2])
#   
#   tpop<-with(as.data.frame(population),akima::interp(x=Longitude,y=Latitude,z=Population,
#                                                      xo=xo,yo=yo,
#                                                      linear=F,extrap = F))$z
#   colnames(tpop)<-yo
#   rownames(tpop)<-xo
#   
#   tpop%>%convMat2SPDF(name="Population")%>%return()
#   
# }

######################################################################################
# Extracts the SEDAC aging data for the bounding box and can also plot it
######################################################################################
GetAgingBbox<-function(directory,bbox,lowres=FALSE,sex="b",plotty=FALSE, namer="Aging",ncity=1){
  # bbox is bounding box in the form 'min lon, max lat, max lon, min lat'
  # DATA: NASA - SEDAC
  
  # PERFORM BASIC CHECKS
  sex<-substr(sex,1,1) ;  sex<-tolower(sex)
  if(!(sex %in% c("b","m","f"))){stop("Error: choice of aging population demography 'sex' is either 'b'=both, 'f'=female or m'=male")}
  if(abs(bbox[2])>90 | abs(bbox[4])>90 | abs(bbox[1])>180 | abs(bbox[3])>180) {stop("Error: non-physical bounding box values in GetAgingBbox")}
  
  subfiler<-paste0("gpw_v4_basic_demographic_characteristics_rev11_a065plus",sex,"t_2010_dens_2pt5_min.asc")
  
  #   GetSEDACfnum(LONG,   LAT)
  LL<-GetSEDACfnum(bbox[1],bbox[2]) # Lower Left
  LR<-GetSEDACfnum(bbox[3],bbox[2]) # Lower right
  UL<-GetSEDACfnum(bbox[1],bbox[4]) # Upper Left
  UR<-GetSEDACfnum(bbox[3],bbox[4]) # Upper right
  
  if(abs(LL[1]-LR[1])>1L) {
    print("Bounding box of SEDAC aging population data is too large, using lower resolution")
    lowres<-TRUE
  }
  
  if (lowres){
    
    filer<-paste0(directory,"Demography_Data/Age/",
                  "gpw_v4_basic_demographic_characteristics_rev11_a065plus",sex,"t_2010_dens_2pt5_min.asc")
    aging<-SortDemoData(filer,bbox)
    
  } else {
    
    poploc<-paste0("Demography_Data/Age/gpw_v4_basic_demographic_characteristics_rev11_a065plus",sex,"t_2010_dens_30_sec")
    popnom<-paste0("gpw_v4_basic_demographic_characteristics_rev11_a065plus",sex,"t_2010_dens_30_sec_")
    strings<-c(directory,poploc,popnom)
    
    aging<-ExtractSEDACS(strings, bbox) 
    
  }
  
  if(plotty){
    
    longData<-melt(aging)
    longData<-longData[longData$value!=0,]
    
    cities<-maps::world.cities%>%filter(lat>bbox[2]&lat<bbox[4]&long>bbox[1]&long<bbox[3])%>%arrange(desc(pop))
    if(ncity>1){wordcloud::wordcloud(words=cities$name,freq = cities$pop,max.words = 30,scale = c(2.5,0.2))}
    cities<-slice(cities,1:ncity)
    
    p<-ggplot(longData, aes(x = Var1, y = Var2)) + 
      geom_raster(aes(fill=value)) + 
      scale_fill_gradient(low="grey90", high="red") +
      labs(x="Longitude", y="Latitude", title="Aging Demographic 65+") +
      theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11)) +
      geom_label(data = cities, aes(long, lat, label = name), size = 4, fontface = "bold", nudge_x = 1.5)#+
    print(p)
    ggsave(paste0(namer,"Aging60plusDemo.eps"), plot=p,path = paste0(directory,'/'),width = 9,height = 7.)
    
    # mad_map <- get_map(bbox,maptype = "toner-background")
    # print(ggmap(mad_map))
    
  }
  
  return(aging)
  
}

# place<-"East Gippsland"
# country<-"Australia"
# bbox<-opq("Australia")$bbox
# viewbox<-NULL#as.double(unlist(strsplit(bbox,",")))
# check<-FALSE

#############################################################################
# Check that place and country names can be used and don't generate rubbish
#############################################################################
PlaceCheck<-function(country,place,viewbox){
  # Place is village/town/city/region,etc but not country!
  if(strcmpi(place,country)){
    writeLines(paste0("Setting place=NULL for '",place,"' and using country '",country,"' in PlaceCheck function"))
    place=NULL
  }
  if(is.null(place)){bb <- getbb (country, viewbox = viewbox, silent = FALSE, format_out = 'polygon')
  }  else {bb <- getbb (paste0(place,", ",country), viewbox = viewbox, silent = FALSE, format_out = 'polygon')}
  if(length(bb)>1){
    for (i in 2:length(bb)){
      disty<-distm(c(mean(bb[[1]][,1]), mean(bb[[1]][,1])), c(mean(bb[[i]][,2]), mean(bb[[i]][,2])), fun = distHaversine) 
      if (disty>5000){    
        writeLines(paste0("
                       ######################################################################################### \n
                       Check the uniqueness of the place specified: '",place,", ",country,"' \n 
                       It appears that there are multiple places with that name, consider using 'viewbox' \n
                       Error location: function 'GetPopulationPlace' in file 'GetPopDemo.R' \n
                       ######################################################################################### \n"))
        stop()
      }
    }
  }
  
  return(bb[[1]])
  
}

######################################################################################
# Extracts the SEDAC population data for a given place and can also plots it
######################################################################################
GetPopulationPlace<-function(directory, country, place=NULL, lowres=FALSE,sex="b", viewbox=NULL,plotty=FALSE, namer=NULL,ncity=1){
  if(is.null(namer)){namer<-paste0(place,"_",country,"_population")}
  
  bb<-PlaceCheck(country,place,viewbox)
  
  # bbox - [mnlo,mnla,mxlo,mxla]
  bbox[1]<-min(bb[[1]][,1],na.rm = T)
  bbox[2]<-min(bb[[1]][,2],na.rm = T)
  bbox[3]<-max(bb[[1]][,1],na.rm = T)
  bbox[4]<-max(bb[[1]][,2],na.rm = T)
  # DATA: NASA - SEDAC
  population<-GetPopulationBbox(directory,bbox,lowres=FALSE,plotty=plotty,namer = namer,ncity=ncity)
  
  if(plotty){
    df<-data.frame(xx=bb[,1],yy=bb[,2])
    mad_map <- get_map(getbb(country),source="stamen",maptype = "toner-background")
    ggmap(mad_map) + coord_fixed() +
      geom_polygon(aes(x = xx, y = yy),data = df, colour = NA, fill = "red", alpha = .2)
  }
  
  stop("Error: code cannot yet integrate population in polygon")
  
  totalpop<-IntPoly(population,bb)
  
  return(totalpop)
  
}

######################################################################################
# Extracts the SEDAC population data for a given place and can also plot it
######################################################################################
GetAgingPlace<-function(directory, country, place=NULL, viewbox=NULL, plotty=FALSE, namer=NULL,ncity=1){
  if(is.null(namer)){namer<-paste0(place,"_",country,"_aging60plus")}
  
  bb<-PlaceCheck(country,place,viewbox)
  
  # bbox - [mnlo,mnla,mxlo,mxla]
  bbox[1]<-min(bb[[1]][,1],na.rm = T)
  bbox[2]<-min(bb[[1]][,2],na.rm = T)
  bbox[3]<-max(bb[[1]][,1],na.rm = T)
  bbox[4]<-max(bb[[1]][,2],na.rm = T)
  # DATA: NASA - SEDAC
  aging<-GetAgingBbox(directory,bbox,lowres=FALSE,plotty=plotty,namer = namer,ncity = ncity)
  
  if(plotty){
    df<-data.frame(xx=bb[[1]][,1],yy=bb[[1]][,2])
    mad_map <- get_map(getbb(country),source="stamen",maptype = "toner-background")
    ggmap(mad_map) + coord_fixed() +
      geom_polygon(aes(x = xx, y = yy),data = df, colour = NA, fill = "red", alpha = .2)
  }
  
  stop("Error: code cannot yet integrate aging population in polygon")
  
  avage<-Mode(aging)
  
  return(avage)
  
}

readFBpop<-function(bbox,saveloc="/tmp/tmp_hrsl.tif"){
  
  # Need bounding box in the form c(minlon,maxlat,maxlon,minlat)
  strbbox<-paste0(bbox[1]," ",bbox[4]," ",bbox[3]," ",bbox[2])
  if (bbox[3] < bbox[1]) strbbox<-paste0(bbox[1]," ",bbox[4]," 180.000 ",bbox[2])
  
  # Linux terminal call to extract the FB population data directly from AWS server
  liner<-paste0("gdal_translate  /vsicurl/https://dataforgood-fb-data.s3.amazonaws.com/hrsl-cogs/hrsl_general/hrsl_general-latest.vrt ", 
                saveloc," -projwin ",strbbox," -projwin_srs EPSG:4326")
  # Call it!
  system(liner)
  # Now read in the file
  FBpop<-raster(x=saveloc)
  # Convert it to how we need it
  FBpop%<>%convRaster2SP(name = "Population")
  # e <- as(raster::extent(c(bbox[c(1,3,2,4)])), 'SpatialPolygons')
  # crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # FBpop%<>%raster::crop(e)
  
  # If the bounding box doesn't cross over the longitude==+/-180 boundary line then return the value
  if (!(bbox[3] < bbox[1])) return(FBpop)
  
  print("Boundary box goes over the +/- 180 degree boundary line, please check values!")
  
  # Otherwise, calculate the rest of the data
  strbbox<-paste0("-180.000 ",bbox[4]," ",bbox[3]," ",bbox[2])
  liner<-paste0("gdal_translate  /vsicurl/https://dataforgood-fb-data.s3.amazonaws.com/hrsl-cogs/hrsl_general/hrsl_general-latest.vrt ", 
                saveloc," -projwin ",strbbox," -projwin_srs EPSG:4326")
  system(liner)
  FBpop2<-raster(x=saveloc)
  FBpop2%<>%convRaster2SP(name = "Population")
  
  fullFB<-FBpop
  fullFB@data$Population<-c(FBpop@data$Population,FBpop2@data$Population)
  fullFB@coords<-rbind(FBpop@coords,FBpop2@coords)
  fullFB@bbox<-c(min(fullFB@coords[,1]),min(fullFB@coords[,2]),max(fullFB@coords[,1]),max(fullFB@coords[,2]))
  
  return(fullFB)
  
}

getFBbuildings <- function(ODDy){
  
  #ODDy <- readRDS('/home/manderso/Documents/GitHub/IIDIPUS_InputRealwithMort/ODDobjects/EQ20110222NZL_-1')
  bbox <- ODDy@bbox
  FBpop <- readFBpop(bbox)
  rastered <- rasterize(FBpop, raster(ODDy), 'Population', fun='sum')
  rastered_spdf <- as(rastered, "SpatialPixelsDataFrame")
  
  data <- ODDy@data
  data$Longitude <-  round(ODDy@coords[,1], 8)
  data$Latitude <- round(ODDy@coords[,2], 8)
  data$id <- 1:NROW(data)
  data <- merge(data, data.frame(Longitude=round(rastered_spdf@coords[,1], 8), Latitude=round(rastered_spdf@coords[,2], 8), FBPop2 = rastered_spdf@data$layer), 
                by=c('Latitude', 'Longitude'), all.x = TRUE)
  data <- data[order(data$id),]
  data$FBPop2[which(is.na(data$ISO3C))] <- NA
  
  #need to retrieve aveHouseholdSize from Global Data Lab
  ODDy@data$FBPop2 <- data$FBPop2
  ODDy@data$nHouses <- data$FBPop2 / 3.17 #replace 3.17 with ODDy@data$aveHouseholdSize from global data lab
  
  return(ODDy)
  # Concerns:
  # - Granularity. Need number of houses to be an integer - rounding produces a fair bit of coarseness. Not so much an issue in dense areas, but e.g. in Canterbury 
  #   35% of cells have populations between 0 and 4 which is less than the average household size. 
  # - date of data vs date of event
  # - Meta Data for Good and SEDAC population data disagree a fair amount. 
  # - Assume Meta Data for Good is working with residential buildings
  # - Missing average household size data from Global Data Lab for some countries
  
}


######################################################################################
# Aggregate the FB population data onto SEDACs grid
######################################################################################
ParAggFBPopSEDAC<-function(ODDy,arrayz,ncores=8, funcy="sum",namer="FBPop", napop=T){
  
  funcy<-match.fun(funcy)
  
  grid<-as.data.frame(ODDy@grid)
  ODDy@data$array<-NA
  arrayz<-data.frame(Longitude=arrayz@coords[,1],
                    Latitude=arrayz@coords[,2],
                    Population=arrayz@data$Population)
  
  # Create a function that takes in input i in 1:ncores that reduces FBpop (rename as tmp) as it goes
  # and returns a list of N/ncores values
  parAGG<-function(kk){
    
    if(napop) ijs<-which(!is.na(ODDy@data$Population)) else ijs<-1:nrow(ODDy@data)
    iiis<-floor((kk-1L)*length(ijs)/ncores+1L):floor((kk)*length(ijs)/ncores)
    if(kk==ncores) iiis<-floor((kk-1L)*length(ijs)/ncores+1L):length(ijs)
    ijs<-ijs[iiis]
    
    inds<-arrayz$Longitude< (min(ODDy@coords[ijs,1]) - grid$cellsize[1]/2)&
      arrayz$Longitude>=(max(ODDy@coords[ijs,1]) + grid$cellsize[1]/2)&
      arrayz$Latitude< (min(ODDy@coords[ijs,2]) - grid$cellsize[2]/2)&
      arrayz$Latitude>=(max(ODDy@coords[ijs,2]) + grid$cellsize[2]/2)
    
    tmp<-arrayz[!inds,]
    
    output<-rep(NA,length(ijs))
    i<-1
    for (z in ijs){
      
      inds<-tmp$Longitude< (ODDy@coords[z,1] + grid$cellsize[1]/2)&
        tmp$Longitude>=(ODDy@coords[z,1] - grid$cellsize[1]/2)&
        tmp$Latitude< (ODDy@coords[z,2] + grid$cellsize[2]/2)&
        tmp$Latitude>=(ODDy@coords[z,2] - grid$cellsize[2]/2)
      output[i]<-funcy(tmp$Population[inds])
      
      tmp<-tmp[!inds,]
      i<-i+1
      
    }
    
    return(output)
    
  }
  
  if(napop) indies<-which(!is.na(ODDy@data$Population)) else indies<-1:nrow(ODDy@data)
  
  ODDy@data$array[indies]<-unlist(mclapply(1:ncores, FUN = parAGG, mc.cores = ncores))
  
  colnames(ODDy@data)[ncol(ODDy@data)]<-namer
  
  return(ODDy)
  
}

AggFBPopSEDAC<-function(ODDy,arrayz,iso3=NULL,funcy="sum",namer="FBPop", napop=T){
  
  funcy<-match.fun(funcy)
  
  grid<-as.data.frame(ODDy@grid)
  ODDy@data$array<-NA
  arrayz<-data.frame(Longitude=arrayz@coords[,1],
                     Latitude=arrayz@coords[,2],
                     Population=arrayz@data$Population)
  
  if(napop) {
    if(!is.null(iso3)) ijs<-which(!is.na(ODDy@data$Population) & ODDy@data$ISO3C==iso3) else ijs<-which(!is.na(ODDy@data$Population))
  } else {
    if(!is.null(iso3)) ijs<-which(ODDy@data$ISO3C==iso3) else ijs<-1:nrow(ODDy@data)
  }
  
  # Along the Longitude, parallelised over latitude arrays
  for (z in ijs){
    
    inds<-arrayz$Longitude< (ODDy@coords[z,1] + grid$cellsize[1]/2)&
      arrayz$Longitude>=(ODDy@coords[z,1] - grid$cellsize[1]/2)&
      arrayz$Latitude< (ODDy@coords[z,2] + grid$cellsize[2]/2)&
      arrayz$Latitude>=(ODDy@coords[z,2] - grid$cellsize[2]/2)
    ODDy@data$array[z]<-funcy(arrayz$Population[inds])
    
    arrayz<-arrayz[!inds,]
    
  }
  
  colnames(ODDy@data)[ncol(ODDy@data)]<-namer
  
  return(ODDy)
  
}

GridUpFBPop<-function(ODDy,ncores=2,funcy="sum",namer="FBPop"){
  
  indies<-!is.na(ODDy@data$Population)
  bbox<-c((min(ODDy@coords[indies,1]) - ODDy@grid@cellsize[1]/2),
          (min(ODDy@coords[indies,2]) - ODDy@grid@cellsize[2]/2),
          (max(ODDy@coords[indies,1]) + ODDy@grid@cellsize[1]/2),
          (max(ODDy@coords[indies,2]) + ODDy@grid@cellsize[2]/2))
  
  FBpop<-readFBpop(bbox)
  
  ODDy<-ParAggFBPopSEDAC(ODDy,FBpop,ncores = ncores,funcy = funcy,namer=namer)
  # ODDy<-AggFBPopSEDAC(ODDy,FBpop)
  rm(FBpop)
  
  if(funcy!="sum") return(ODDy)
  
  # Print the different values of population
  print(paste0("SEDACS: ",sum(ODDy@data$Population,na.rm = T),
               ", FB: ",sum(ODDy@data$FBPop,na.rm = T),
               ", DIFF(%) = ",100*(sum(ODDy@data$Population,na.rm = T)-sum(ODDy@data$FBPop,na.rm = T))/sum(ODDy@data$Population,na.rm = T)))
  print(".....")
  # Normalise values wrt World Bank population national estimate
  factor<-InterpPopWB(unique(na.omit(ODDy@data$ISO3C)),ODDy@hazdates[1],normdate=as.Date("2018-01-01"))  
  
  for(iso in unique(na.omit(ODDy@data$ISO3C))){
    ODDy@data$FBPop[!is.na(ODDy@data$ISO3C) & ODDy@data$ISO3C==iso]<-
      ODDy@data$FBPop[!is.na(ODDy@data$ISO3C) & ODDy@data$ISO3C==iso]*factor$factor[factor$iso3==iso]
  }
  
  return(ODDy)
}

ChangeAllSEDACS_FB<-function(folderin="./IIDIPUS_Results/ODDobjects/",folderout="./IIDIPUS_Results/ODDobjects_FB/"){
  
  filez<-list.files(folderin)
  x <- file.info(paste0(folderin,filez))
  filez<-filez[match(length(filez):1,rank(-x$size))]
  
  i<-136
  for (fff in filez[136:length(filez)]){
  # for (fff in filez){
    
    print(paste0(fff, ",   list number = ",i))
    ODDy<-readRDS(paste0(folderin,fff))
    ODDy<-GridUpFBPop(ODDy,ncores=2)
    saveRDS(ODDy,paste0(folderout,fff))
    i<-i+1
    
  }
  
  
}

# AggregateFBPopSEDAC<-function(ODDy,FBpop,cores=4){
#   
#   grid<-as.data.frame(ODDy@grid)
#   
#   # Along the Longitude, parallelised over latitude arrays
#   unlist(mclapply(X = 1:length(ODDy@data$Population),
#                   FUN = function(i) sum(FBpop@data$Population[FBpop@coords[,1]< (ODDy@coords[i,1] + grid$cellsize[1])&
#                                                                 FBpop@coords[,1]>=(ODDy@coords[i,1] - grid$cellsize[1])&
#                                                                 FBpop@coords[,2]< (ODDy@coords[i,2] + grid$cellsize[2])&
#                                                                 FBpop@coords[,2]>=(ODDy@coords[i,2] - grid$cellsize[2])]),
#                   mc.cores = cores))
#   
# }

# 
## compare OSM and Meta Data For Good - all very messy and will delete once I've sorted out the building count data

# ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/ODDy_NZL_OSMvsFB')
# which.max(ODDy@data$nHouses-ODDy@data$Buildings_OSM)
# ODDy@data$Buildings_OSM[8333]
# bbox_interest <- t(rbind(ODDy@coords[8333,]-ODDy@grid@cellsize/2, ODDy@coords[8333,]+ODDy@grid@cellsize/2))
# FBpop_interest<-readFBpop(bbox_interest)
# FBpop_interest@coords
# sum(FBpop_interest@data$Population)
# 
# 
# #Compare SEDAC and Meta Data for Good
# bbox_small <- ODDy@bbox
# bbox_small[1:2,1:2] <- rbind(c(124.9, 124.95), c(8.9,8.95))
# SEDAC_small <- GetPopulationBbox(dir, bbox=bbox_small, yr='2020')
# 
# buildings<-GetOSMbuildingsODD(ODDy, bbox = ODDy@bbox, timeout=60)
# rastered <- rasterize(cbind(buildings$Longitude, buildings$Latitude), raster(ODDy), fun='count')
# rastered_spdf <- as(rastered, "SpatialPixelsDataFrame")
# 
# data <- ODDy@data
# data$Longitude <-  round(ODDy@coords[,1], 8)
# data$Latitude <- round(ODDy@coords[,2], 8)
# data$id <- 1:NROW(data)
# data <- merge(data, data.frame(Longitude=round(rastered_spdf@coords[,1], 8), Latitude=round(rastered_spdf@coords[,2], 8), Buildings_OSM = rastered_spdf@data$layer),
#               by=c('Latitude', 'Longitude'), all.x = TRUE)
# data <- data[order(data$id),]
# data$Buildings_OSM[which(is.na(ODDy@data$Buildings_OSM))] <- 0
# data$Buildings_OSM[which(is.na(data$ISO3C))] <- NA
# 
# OSM_small <-
# FBpop_small<-readFBpop(SEDAC_small@bbox)
# plot(SEDAC_small)
# points(FBpop_small, col='white', pch=20, cex=FBpop_small$Population/5)
# 
# 
# #compare OSM with SEDAC
# 
# ODDy <- readRDS('/home/manderso/Documents/GitHub/IIDIPUS_InputRealUnedited/ODDobjects/EQ20170210PHL_816')
# bbox <- ODDy@bbox
# FBpop <- readFBpop(bbox)
# rastered <- rasterize(FBpop, raster(ODDy), 'Population', fun='sum')
# rastered_spdf <- as(rastered, "SpatialPixelsDataFrame")
# 
# data <- ODDy@data
# data$Longitude <-  round(ODDy@coords[,1], 8)
# data$Latitude <- round(ODDy@coords[,2], 8)
# data$id <- 1:NROW(data)
# data <- merge(data, data.frame(Longitude=round(rastered_spdf@coords[,1], 8), Latitude=round(rastered_spdf@coords[,2], 8), FBPop2 = rastered_spdf@data$layer), 
#               by=c('Latitude', 'Longitude'), all.x = TRUE)
# data <- data[order(data$id),]
# ODDy@data$FBPop2 <- data$FBPop2
# 
# bbox_small <- bbox
# bbox_small[1:2, 1:2] <- rbind(c(125.5, 125.6), c(9.6,9.71))
# 
# buildings<-GetOSMbuildingsODD(ODDy, bbox = bbox_small, timeout=60)
# rastered <- rasterize(cbind(buildings$Longitude, buildings$Latitude), raster(ODDy), fun='count')
# rastered_spdf <- as(rastered, "SpatialPixelsDataFrame")
# data <- ODDy@data
# data$Longitude <-  round(ODDy@coords[,1], 8)
# data$Latitude <- round(ODDy@coords[,2], 8)
# data$id <- 1:NROW(data)
# data <- merge(data, data.frame(Longitude=round(rastered_spdf@coords[,1], 8), Latitude=round(rastered_spdf@coords[,2], 8), Buildings_OSM = rastered_spdf@data$layer),
#               by=c('Latitude', 'Longitude'), all.x=T)
# data <- data[order(data$id),]
# data$Buildings_OSM[which(is.na(ODDy@data$Buildings_OSM))] <- 0
# data$Buildings_OSM[which(is.na(data$ISO3C))] <- NA
# 
# 
# 
# OSM_small <-
#   FBpop_small<-readFBpop(SEDAC_small@bbox)
# plot(SEDAC_small)
# points(FBpop_small, col='white', pch=20, cex=FBpop_small$Population/5)
# 
# 
# 
# 
# FBpop <- readFBpop(bbox)
# rastered <- rasterize(FBpop, raster(SEDAC_Pop), 'Population', fun='sum')
# rastered_spdf <- as(rastered, "SpatialPixelsDataFrame")
# 
# data <- SEDAC_Pop@data
# data$Longitude <-  round(SEDAC_Pop@coords[,1], 8)
# data$Latitude <- round(SEDAC_Pop@coords[,2], 8)
# data$id <- 1:NROW(data)
# data <- merge(data, data.frame(Longitude=round(rastered_spdf@coords[,1], 8), Latitude=round(rastered_spdf@coords[,2], 8), FBPop2 = rastered_spdf@data$layer),
#               by=c('Latitude', 'Longitude'), all.x = TRUE)
# data <- data[order(data$id),]
# data[which(is.na(ODDy@data$ISO3C)),'FBPop2'] <- NA
# par(mfrow=c(1,1))
# plot(data$Population, data$FBPop2, xlab='SEDAC Population', ylab='Meta Data for Good Population')
# abline(0,1)
# 
# plot(sort(data$Population/data$FBPop2), ylim=c(0,2))
# #need to retrieve aveHouseholdSize from Global Data Lab
# ODDy@data$FBPop2 <- data$FBPop2
# ODDy@data$nHouses <- data$FBPop2 / 3.17 #replace 3.17 with ODDy@data$aveHouseholdSize from global data lab
# 
# rbPal <- colorRampPalette(c('red','blue'))
# Col <- rbPal(10)[as.numeric(cut(FBpop@data$Population[141000:150000],breaks = 10))]
# plot(FBpop@coords[141000:150000,1], FBpop@coords[141000:150000,2], col=Col, pch=20, cex=0.1)
# plot_df <- data.frame(
#   x =  FBpop@coords[,1],
#   y = FBpop@coords[,2],
#   width = .08333333333333333333/30,
#   height = .08333333333333333333/30,
#   Pop = FBpop@data$Population
# )
# 
# cols = rainbow(26, s=.6, v=.9)[sample(1:26,26)]
# 
# ggplot(plot_df %>% filter(x>172.55 & x < 172.6 & y > -43.6 & y < -43.55), aes(x, y, width, height)) +
#   geom_tile(aes(fill=as.factor(Pop)))+scale_fill_manual(values=rainbow(length(Pop), s=.6, v=.9)[sample(1:length(Pop),length(Pop))])
# 
# SEDAC_Pop_aligned <- SEDAC_Pop@data$Population[1:23986]
# 
# plot(ODDy@data$FBPop2, SEDAC_Pop_aligned)
# SEDAC_Pop_aligned/ODDy@data$FBPop2
# plot(ODDy@data$Population, SEDAC_Pop_aligned)
# 
# plot(ODDy_OSM@data$FBPop2/3.14, ODDy_OSM@data$Buildings_OSM, ylab='Open Street Maps', xlab='Meta Data for good Pop / Ave Household Size')
# abline(0,1)
