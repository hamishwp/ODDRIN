
# # 2D DENSITY PLOT WITH MODIFIED COLOUR AXIS
# ggplot(as.data.frame(buildings),aes(Longitude,Latitude))+
#   stat_density_2d(aes(fill = ..level..),breaks=c(0,0.5,1,3,5,10,30,50,100,300), 
#                   geom = "polygon")
# p1<-ggplot(as.data.frame(BDy),aes(Longitude,Latitude))+
#   stat_density_2d_filled(aes(fill=..level..),breaks=c(0,0.5,1,3,5,10,30,50,100,300)) +
#   ggtitle("Building Damage Data, EQ2015-04-25NPL")+theme(plot.title = element_text(hjust = 0.5))
# p2<-ggplot(as.data.frame(buildings),aes(Longitude,Latitude))+
#   stat_density_2d_filled(aes(fill=..level..),breaks=c(0,0.5,1,3,5,10,30,50,100,300)) +
#   ggtitle("OSM Building Data")+theme(plot.title = element_text(hjust = 0.5))
# gridExtra::grid.arrange(p1, p2, ncol=2)
# 
# tmp<-data.frame(PopDens=c(buildings$Pdensity,BDy$Population),ID=c(rep("OSM",nrow(buildings)),rep("Damaged",nrow(BDy))))
# ggplot(tmp,aes(x=PopDens,group=ID))+ stat_ecdf(aes(colour=ID),geom = "step",size=2) + scale_x_log10() + xlab("Population Density") +
#   ylab("Cumulative Distribution Function") + ggtitle("Nepal 2015 Earthquake") +theme(plot.title = element_text(hjust = 0.5))



# ggplot(tmp,aes(OSM,Damaged))+
#   stat_density_2d_filled(aes(fill=..level..),breaks=c(0,10,30,50,100,300,500,10000), geom="tile") +
#   ggtitle("Population Density of Damaged vs OSM Data")+theme(plot.title = element_text(hjust = 0.5))


CheckBbox<-function(bbox,correct=TRUE){
  # bbox should be [min_longitude, min_latitude, max_longitude, max_latitude]
  if(abs(bbox[2])>90 | abs(bbox[4])>90 | abs(bbox[1])>180 | abs(bbox[3])>180) {
    stop("Error: non-physical bounding box values given to IIDIPUS")
  }
  
  if (bbox[1]>bbox[3]) {
    print("WARNING: bounding box range is long{-180,180}, lat{-90,90} with bbox[min_lon,min_lat,max_long,max_lat]")
    print("Found min_long>max_long")
    if(correct){
      tmp<-bbox[3]
      bbox[3]<-bbox[1]
      bbox[1]<-tmp
    }
  }
  if (bbox[2]>bbox[4]) {
    print("WARNING: bounding box range is long{-180,180}, lat{-90,90} with bbox[min_lon,min_lat,max_long,max_lat]")
    print("Found min_lat>max_lat")
    if(correct){
      tmp<-bbox[4]
      bbox[4]<-bbox[2]
      bbox[2]<-tmp
    }
  }
  
  return(bbox)
}

CheckArgs<-function(args){
  # INDEX - TYPE - NAME - DESCRIPTION
  # 1:4 - numeric   - bbox[min lon, max lat, max lon, min lat] - region bounding box
  # 5   - character - country - ISO format only
  # 6   - character - hazard  - using IDMC definition: {"Flood","Storm","Mass movement","Wildfire","Earthquake","Extreme temperature","Volcanic eruption","Drought"}
  # 7   - date (%Y%m%d) - sdate - start date of event
  # 8   - date (%Y%m%d) - fdate - start date of event
  bbox<-as.numeric(args[1:4])
  bbox<-CheckBbox(bbox)
  
  iso3<-as.character(args[5])
  if(nchar(args[5])>3) {
    print(paste0("Warning: detected country name instead of iso3: ",args[5]))
    args[5]%<>%countrycode(origin ='country.name', destination ='iso3c')
  } else if (nchar(args[5])==2){
    print(paste0("Warning: detected iso2 instead of iso3 for country: ",args[5]))
    args[5]%<>%countrycode(origin ='iso2c', destination ='iso3c')
  } else if (nchar(args[5])<2) stop("Error in country input, try using iso3 value")
  args[5]%<>%countrycode(origin ='iso3c', destination ='iso3c')
  
  # list of possible IDMC hazards
  hazard_type<-as.character(args[6])
  hazards<-c("Flood","Storm","Mass movement","Wildfire","Earthquake","Extreme temperature","Volcanic eruption","Drought")
  if(!(hazard_type %in% hazards)) stop("Error: hazard not found among possible IDMC hazard types")
  
  sdate<-args[7]%>%as.POSIXct()%>%as.Date(format = "%Y%m%d")
  fdate<-args[8]%>%as.POSIXct()%>%as.Date(format = "%Y%m%d")
  
  if(!any(grepl("20",c(sdate,fdate)))) stop("Error: date input requires full year e.g. 2019")
  if(fdate<sdate) stop("Error: hazard end date must preceed start date")
  if(any(c(sdate,fdate)>Sys.Date()+10)) stop("Error: hazard cannot be in the future")
  if(any(c(sdate,fdate)<"2017-01-01")) stop("Error: hazard cannot be before 2017")
  year<-format(sdate,"%Y")
  return(list(bbox=bbox,iso3=iso3,hazard_type=hazard_type,sdate=sdate,fdate=fdate,year=year))
  
}

areaBbox<-function(bbox){
  s1<-cbind(lon=c(bbox[1],bbox[3],bbox[1],bbox[3]),lat=c(bbox[2],bbox[2],bbox[4],bbox[4]))
  sp1 <- spPolygons(s1, crs="+proj=longlat +datum=WGS84")
  mp1 <- makePoly(sp1, interval=100000)
  return(areaPolygon(mp1)*1e-3)
  # R<-6378.137
  # return((pi/180)*R^2 *abs(sin(bbox[2])-sin(bbox[4]))*abs(bbox[1]-bbox[3]))
}

# convDF2SPDF<-function(DF,name=NULL,crs="WGS84"){
#   
#   
#   
#   if(is.null(name)) name<-"Value"
#   
#   0
#   
#   if(crs=="WGS84") {crs(DF)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
#   } else {stop("ERROR: Unknown coordinate system in convMat2SPDF, see Functions.R")}
#   
#   return(DF)
# }
checkMatlonglat<-function(array){
  
  long<-as.numeric(rownames(array))
  lat<-as.numeric(colnames(array))
  
  if(!(any(is.na(long)) | any(is.na(lat)))) {
    colnames(array)<-lat
    rownames(array)<-long
    
    array%<>%reshape2::melt()
    
    return(array)
  }
  
  ladiff<-median(diff(lat),na.rm = T)
  lodiff<-median(diff(long),na.rm = T)
  
  # Unfortunately sometimes the GetPopDemo function returns a character in col/rownames
  # Additionally, sometimes diff(long/lat) is not unique
  # Let's fix that!
  if(is.na(lat[1])) lat[1]<-lat[2]-ladiff
  if(is.na(lat[length(lat)])) lat[length(lat)]<-lat[length(lat)-1]+ladiff 
  if(is.na(long[1])) long[1]<-long[2]-lodiff 
  if(is.na(long[length(long)])) long[length(long)]<-long[length(long)-1]+lodiff 
  
  if(any(is.na(long)) | any(is.na(lat))) stop("nan values in longitude/latitude values of array col/row names")
  
  colnames(array)<-lat
  rownames(array)<-long
  
  array%<>%reshape2::melt()
  # if(array$Var1!=long | array$Var2!=lat) {
  #   
  # }
  
  return(array)
}

convSPDF2Array<-function(ODDy){
  if("ISO3C"%in%colnames(ODDy@data)) ODDy@data%<>%dplyr::select(-ISO3C)
  array(as.matrix(ODDy@data),
        dim = unname(c(ODDy@grid@cells.dim,
                       ncol(ODDy@data))))
}

convRaster2SPDF<-function(raster,name=NULL){
  
  raster%<>%as("SpatialPixelsDataFrame")
  colnames(raster@coords) <- c("Longitude","Latitude")
  rownames(raster@bbox) <- c("Longitude","Latitude")
  if(!is.null(name)) colnames(raster@data)<-name
  
  return(raster)
  
}

convRaster2SP<-function(raster,name=NULL){
  
  raster%<>%as("SpatialPointsDataFrame")
  colnames(raster@coords) <- c("Longitude","Latitude")
  rownames(raster@bbox) <- c("Longitude","Latitude")
  if(!is.null(name)) colnames(raster@data)<-name
  
  return(raster)
  
}

convMat2DF<-function(array,name=NULL){
  
  array%<>%checkMatlonglat()
  
  if(is.null(name)) name<-"Value"
  colnames(array)<-c("Longitude","Latitude",name)
  
  return(array)
  
}

convMat2raster<-function(array,name=NULL,crs="WGS84"){
  
  array%<>%checkMatlonglat()
  if(is.null(name)) name<-"Value"
  colnames(array)<-c("Longitude","Latitude",name)
  array %<>%raster
  
  if(crs=="WGS84") {crs(array)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
  } else {stop("ERROR: Unknown coordinate system in convMat2SPDF, see Functions.R")}
  
  return(array)
  
}

convMat2SPDF<-function(array,name=NULL,crs="WGS84"){
  
  array%<>%checkMatlonglat()
  
  if(is.null(name)) name<-"Value"
  colnames(array)<-c("Longitude","Latitude",name)
  array <- SpatialPixelsDataFrame(points = array[c("Longitude","Latitude")],
                                  data = array[name])
  
  if(crs=="WGS84") {crs(array)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
  } else {stop("ERROR: Unknown coordinate system in convMat2SPDF, see Functions.R")}
  
  ##### THIS SECTION IS TO ORDER THE VECTOR VALUES OF THE POPULATION MATRIX
  ##### OTHERWISE THE HAZARD INTERPOLATION IS SPLIT IN THE HORIZONTAL PLANE
  
  xo<-array@coords[1:array@grid@cells.dim[1],1]
  yo<-array@coords[1:array@grid@cells.dim[2]*array@grid@cells.dim[1]-array@grid@cells.dim[1]+1,2]
  
  if(!any(sort(yo)==yo)) {
    
    # find index to split data.frame
    ind<-which.min(array@coords[,2])-1L
    array@data[[name]]<-c(array@data[[name]][(ind+1):nrow(array)],array@data[[name]][1:ind])
    array@grid.index<-c(array@grid.index[(ind+1):nrow(array)],array@grid.index[1:ind])
    array@coords[,1]<-c(array@coords[(ind+1):nrow(array),1],array@coords[1:ind,1])
    array@coords[,2]<-c(array@coords[(ind+1):nrow(array),2],array@coords[1:ind,2])
    
  }
  # 
  # if(!any(sort(xo)==xo)) {
  # 
  #   # find index to split data.frame
  #   ind<-which.min(array@coords[1:array@grid@cells.dim[1],1])-1
  #   pop<-xco<-yco<-gind<-array(NA,nrow(array))
  #   for(i in 1:array@grid@cells.dim[2]){
  #     pop[(1:array@grid@cells.dim[1])*i]<-c(array@data[[name]][((ind+1):array@grid@cells.dim[1])*i],array@data[[name]][(1:ind)*i])
  #     gind[(1:array@grid@cells.dim[1])*i]<-c(array@grid.index[((ind+1):array@grid@cells.dim[1])*i],array@grid.index[(1:ind)*i])
  #     xco[(1:array@grid@cells.dim[1])*i]<-c(array@coords[((ind+1):array@grid@cells.dim[1])*i,1],array@coords[(1:ind)*i,1])
  #     yco[(1:array@grid@cells.dim[1])*i]<-c(array@coords[((ind+1):array@grid@cells.dim[1])*i,2],array@coords[(1:ind)*i,2])
  #   }
  #   array@data[[name]]<-pop
  #   array@grid.index<-gind
  #   array@coords[,1]<-xco
  #   array@coords[,2]<-yco
  # 
  # }
  
  return(array)
  
}

# Function to extract the coordinates of all the polygons
extractPolyCoords<-function(ADM){
  coords<-data.frame()
  for(i in 1:length(ADM@polygons)){
    for(j in 1:length(ADM@polygons[[i]]@Polygons)){
      coords%<>%rbind(data.frame(LONGITUDE=ADM@polygons[[1]]@Polygons[[1]]@coords[,1],
                                 LATITUDE=ADM@polygons[[1]]@Polygons[[1]]@coords[,2],
                                 i=i,j=j))
    }
  }
  return(coords)
}
# Find which coordinates of an S4 spatial object lie inside (or on the boundary of) a spatial polygon file
inPoly<-function(poly,pop,iii=1,sumFn="sum",reducer=NULL){
  
  Ifin<-rep(F,nrow(pop))
  if(is.null(reducer)) reducer<-!Ifin
  
  pop<-pop[reducer,]
  
  if(any(class(pop)=="SpatialPointsDataFrame") | any(class(pop)=="SpatialPixelsDataFrame")){
    coords<-pop@coords
    data<-pop@data
  } else {
    coords<-as.data.frame(pop[,c("Longitude","Latitude")])
    data<-as.data.frame(pop)
  }
  
  insidepoly<-rep(FALSE,nrow(pop))
  
  for (i in 1:length(poly@Polygons)){
    # Get rid of values outside the bounding box first
    minipoly<-rep(FALSE,length(insidepoly))
    indies<-coords[,1]>=min(poly@Polygons[[i]]@coords[,1]) &
      coords[,1]<=max(poly@Polygons[[i]]@coords[,1]) &
      coords[,2]>=min(poly@Polygons[[i]]@coords[,2]) &
      coords[,2]<=max(poly@Polygons[[i]]@coords[,2])
    # Now we only need to calculate a few points that lie inside the polygon!
    minipoly[indies]<-sp::point.in.polygon(coords[indies,1],
                                           coords[indies,2],
                                           poly@Polygons[[i]]@coords[,1],
                                           poly@Polygons[[i]]@coords[,2])>0
    # Add to the total
    insidepoly<- insidepoly | minipoly
  }
  
  outer<-match.fun(sumFn)(data[insidepoly,iii],na.rm=T)
  
  Ifin[reducer]<-insidepoly
  
  return(list(vals=outer,indies=Ifin))
}

expandBbox<-function(bbox,f,scaling=T){
  Dx<-bbox[3]-bbox[1]
  Dy<-bbox[4]-bbox[2]
  A<-Dx*Dy
  if(scaling) dx<-Dx*(-1+sqrt(f))
  else dx<-Dx*(-1+sqrt(1+(f-A)/A)) # use f as as the resulting area
  dy<-dx*Dy/Dx
  return(bbox+0.5*c(-dx,-dy,dx,dy))
}

convIso2Iso3<-function(iso2){
  countrycode::countrycode(sourcevar = iso2,
                           origin = "iso2c",
                           destination = "iso3c",warn = F)
}

convIso3Country<-function(iso3){
  countrycode::countrycode(sourcevar = iso3,
                           origin = "iso3c",
                           destination = "country.name",warn = F)
}

convIso3Continent<-function(iso3){
  # continents<-countrycode::countrycode(sourcevar = iso3,
  #                                      origin = "iso3c",
  #                                      destination = "continent",warn = F)
  left_join(data.frame(ISO3=iso3),raster::ccodes()[,c("ISO3","continent")],by="ISO3")$continent
}

InterpDay<-function(ndata,day){
  val<-data.frame()
  for (iso3c in unique(ndata$iso3)){
    nd<-filter(ndata,iso3==iso3c)
    if(all(is.na(nd$value))|length(nd$value)<=1) {
      print(paste0("Not enough data found for country ",iso3c," for normalisation spline for country indicators"))
      val%<>%rbind(data.frame(iso3=iso3c,value=NA))
      next
    }
    func = tryCatch(splinefun(x=nd$day,y=nd$value),error = function(e) NULL)
    if(is.null(func)) {
      print(paste0("No spline function possible for country ",iso3c," values: ",nd$value))
      value<-nd$value[which.min(abs(nd$day-day))]
      val%<>%rbind(data.frame(iso3=iso3c,value=value))
    } else val%<>%rbind(data.frame(iso3=iso3c,value=func(day)))
  }
  return(val)
}

coords2country = function(points,iso=T)
{  
  if(dim(points)[2]!=2) stop("coords2country Error: long/lat coords are invalid")
  if(dim(points)[1]==1) points=t(points[1,])
  
  countriesSP <- rworldmap::getMap(resolution='high')
  #setting CRS directly to that from rworldmap
  pointsSP = sp::SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = sp::over(pointsSP, countriesSP)
  
  # return the ISO3 of country
  if(iso) return(as.character(indices$ISO3))
  # return the ADMIN names of country
  return(as.character(indices$ADMIN))
}

countriesbbox<-function(iso3){
  
  countriesSP <- rworldmap::getMap(resolution='low')
  indies<-which(countriesSP$ISO3%in%iso3)
  
  mnlo<-mxlo<-mnla<-mxla<-c()
  for(c in 1:length(indies)){
    indy<-indies[c]
    for (i in 1:length(countriesSP@polygons[[indy]]@Polygons)){
      mnlo<-min(c(mnlo,countriesSP@polygons[[indy]]@Polygons[[i]]@coords[,1]))
      mxlo<-max(c(mxlo,countriesSP@polygons[[indy]]@Polygons[[i]]@coords[,1])) 
      mnla<-min(c(mnla,countriesSP@polygons[[indy]]@Polygons[[i]]@coords[,2]))
      mxla<-max(c(mxla,countriesSP@polygons[[indy]]@Polygons[[i]]@coords[,2]))
    }
  }
  
  return(c(mnlo,mnla,mxlo,mxla))
  
}

PlotDisaster<-function(pop,dfpoly,bbox=NULL,map=FALSE,ncity=1,namer="Disaster",filer="./"){
  
  if(is.null(bbox)) bbox<-as.numeric(c(min(rownames(pop)),min(colnames(pop)),max(rownames(pop)),max(colnames(pop))))
  longData<-reshape2::melt(pop)
  longData<-longData[longData$value!=0,]
  
  cities<-maps::world.cities%>%filter(lat>bbox[2]&lat<bbox[4]&long>bbox[1]&long<bbox[3])%>%arrange(desc(pop))
  if(ncity>1){wordcloud::wordcloud(words=cities$name,freq = cities$pop,max.words = 30,scale = c(2.5,0.2))}
  cities<-slice(cities,1:ncity)
  
  p<-ggplot(longData, aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low="gray80", high="black") +
    labs(x="Longitude", y="Latitude", title=namer) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  for (j in unique(dfpoly$ncontour)){
    tp<-filter(dfpoly,ncontour==j)
    p<-p+geom_polygon(data = tp,aes(x=Longitude,y=Latitude,group=Intensity,colour=Intensity),alpha=0,na.rm = T,size=2)+
      scale_color_gradient(low="mistyrose2", high="red")
  }
  p<-p+geom_label(data = cities, aes(long, lat, label = name), size = 4, fontface = "bold", nudge_x = 0.05*(bbox[3]-bbox[1]))
  
  print(p)
  if(!is.null(filer)) ggsave(paste0(namer,".eps"), plot=p,path = filer,width = 9,height = 7.)
  return(p)
}

# GetMapObj<-function(bbox,world=NULL){
#   
#   if(is.null(world)){
#     library("rnaturalearth")
#     library("rnaturalearthdata")
#     world <- ne_countries(scale = "medium", returnclass = "sf")
#   }
#   
#   aj<-c(abs(bbox[1]-bbox[3])*0.05,abs(bbox[4]-bbox[2])*0.05)
#   p<- ggplot(data = world) + geom_sf(fill= "antiquewhite") + 
#     coord_sf(xlim = c(bbox[1]-aj[1],bbox[3]+aj[1]), ylim = c(bbox[2]-aj[2],bbox[4]+aj[2]), expand = FALSE) + 
#     xlab("Longitude") + ylab("Latitude") +
#     theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))
#   
#   return(p)
# }

GetWeightedRoR<-function(x,y){
  
  len<-length(y)
  w<-rep(1,len)
  w[len]<-len
  
  fit<-lm(log(y) ~ x,weights = w)
  sumz<-summary(fit)
  # Extract gradient, intercept and p-value of gradient
  return(c(sumz$coefficients[2, 1],exp(sumz$coefficients[1, 1]),sumz$coefficients[2, 4]))
}

ggmap_bbox <- function(map,bbox) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  # map_bbox <- setNames(unlist(attr(map, "bb")), 
  # c("ymin", "xmin", "ymax", "xmax"))
  
  # Coonvert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  # bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox["ymin"]
  attr(map, "bb")$ll.lon <- bbox["xmin"]
  attr(map, "bb")$ur.lat <- bbox["ymax"]
  attr(map, "bb")$ur.lon <- bbox["xmax"]
  map
}

GetExtent<-function(ADM,expander=NULL){
  bbox<-ADM@bbox
  # Expand the bounding box, useful for the interpolation
  if(!is.null(expander)) bbox%<>%expandBbox(1.1)
  ext<-as(extent(bbox[c(1,3,2,4)]), 'SpatialPolygons')
  crs(ext) <- "+proj=longlat +datum=WGS84 +no_defs"  
  return(ext)
}

# Remove unnecessary columns that arise with the R st_multipolygon format
ExtractPolyCoords<-function(polys){
  coords<-st_coordinates(polys)
  ind<-c()
  for (j in 3:ncol(coords)){
    if(max(coords[,j])==1) ind<-c(ind,j)
  }
  if(!is.null(ind)) coords<-coords[,-ind]
  if(ncol(coords)==2) {
    coords<-cbind(coords,rep(1,dim(coords)[1]))
  } else if(ncol(coords)>3) {
    
    # filter per L3 column by L2 value. Add a value per unique value in L2
    coords%<>%as.data.frame()
    L1<-coords%>%filter(L1==unique(coords$L1)[1])
    L2<-max(L1$L2)
    for (j in unique(coords$L1)[-1]){
      L1<-coords%>%filter(L1==j)
      coords$L2[coords$L1==j]<-L1$L2+L2-min(L1$L2)+1
      L2<-max(coords$L2)
    }
    coords<-coords[,-3]
    coords%<>%arrange(L2)%>%as.matrix()
    
  } else if(ncol(coords)>5) {
    stop("Error in extracting polygon coordinates, see ExtractPolyCoords in GetGDACS.R")
  }
  
  return(coords)
  
}
