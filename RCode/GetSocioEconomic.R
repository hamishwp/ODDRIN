library(openxlsx)
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(dplyr)
# source('../../IIDIPUS/RCode/Functions.R')
library(magrittr)
library(wbstats)
library(wid)
library(reshape2)
library(raster)
library(sf)
library(sp)

FilterKummu<-function(GDP,bbox,melted=F){
  
  lat<-as.numeric(colnames(GDP))
  lon<-as.numeric(rownames(GDP))
  nlat<- length(lat)
  nlon<- length(lon)
  
  imnlon<-which.min(abs(lon-bbox[1]))
  if(imnlon>1) imnlon<-imnlon-1
  imxlon<-which.min(abs(lon-bbox[3]))
  if(imxlon<nlon) imxlon<-imxlon+1
  imnlat<-which.min(abs(lat-bbox[2]))
  if(imnlat>1) imnlat<-imnlat-1
  imxlat<-which.min(abs(lat-bbox[4]))
  if(imxlat<nlat) imxlat<-imxlat+1
  
  # GDP_PPP[longitude,latitude]
  if(!melted) return(GDP[imnlon:imxlon,imnlat:imxlat])
  GDP<-GDP[imnlon:imxlon,imnlat:imxlat]
  GDP<-melt(GDP);colnames(GDP)<-c("X","Y","data")
  return(GDP)
  
}

GetKummu<-function(dir,bbox=NULL,yr=2015L){
  
  iii<-yr-1989L
  
  # file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_PPP_30arcsec_v3.nc")
  file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc")
  GDP<-brick(file,varname="GDP_per_capita_PPP")
  GDP<-GDP[[iii]]
  
  if(!is.null(bbox)) {
    e <- as(raster::extent(c(bbox[c(1,3,2,4)])), 'SpatialPolygons')
    crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    GDP%<>%raster::crop(e)
  }
  
  GDP%<>%as('SpatialPixelsDataFrame')
  
  return(GDP)
  
}

plotGDP<-function(GDP,zoom=5){
  mad_map <- get_stamenmap(GDP@bbox,source = "stamen",maptype = "toner",zoom=zoom)
  p<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")
  p+geom_contour_filled(data = as.data.frame(GDP),
                        mapping = aes(x,y,z=X2015),
                        alpha=0.7)+ 
    labs(fill = "GDP-PPP [USD-2015]")
}

GetWID_natincome<-function(year=NULL,iso3c=NULL){
  
  if(is.null(year)) year<-"2015"
  
  WID<-download_wid(indicators = "aptinc",years=year, perc="p0p100", pop="j") %>%rename(value_lcu=value)
  ppp <-download_wid(indicators = "xlcusp",year = year) %>% # areas =unique(WID$country),
    rename(ppp=value) %>%select(-year, -percentile)
  WID %<>%merge(ppp, by="country") %>%mutate(value_ppp=value_lcu/ppp)
  
  
  WID%<>%mutate(iso3=convIso2Iso3(country))%>%
    dplyr::select(-c(variable.x,variable.y,country,year,percentile,value_lcu,ppp)) 
  if(!is.null(iso3c)) WID %<>% filter(iso3%in%iso3c)
  
  return(WID)
}

#https://data.worldbank.org/indicator/SI.DST.10TH.10?end=2012&locations=SB&start=2011 (data from 2012 but only source I can find)
SLB_WID <- data.frame(percentile = paste0('p',seq(0,90,10), 'p', seq(10,100,10)), 
                      value=c(0.028, 0.07-0.028,0.114/2,0.114/2,0.155/2, 0.155/2, 0.215/2, 0.215/2, .446-.292, 0.292), 
                      iso3='SLB')

#https://data.worldbank.org/indicator/SI.DST.10TH.10?end=2012&locations=SB&start=2011 (data from 2012 but only source I can find)
VUT_WID <- data.frame(percentile = paste0('p',seq(0,90,10), 'p', seq(10,100,10)), 
                      value=c(0.03, 0.075-0.03,0.124/2,0.124/2,0.172/2, 0.172/2, 0.23/2, 0.23/2, .399-.247, 0.247), 
                      iso3='VUT')

filter_WID_by_iso3c <- function(WID_all, iso3c){
  return(WID_all %>% filter(variable=='sptinc992j')%>%
                   mutate(iso3=ifelse(country=='KS', 'KOS', convIso2Iso3(country)))%>%
                   dplyr::select(-c(variable,country,year)) %>%
                   filter(iso3%in%iso3c))
}

getMissingWID <- function(missing_iso3c, WID, WID_all){
  if ('PRI' %in% missing_iso3c){ # Use WID data from Panama in place of Puerto Rico as it has an extremely similar Gini Coefficient 
              # according to: https://link.springer.com/article/10.1007/s11205-022-03010-8
    PRI_WID <- filter_WID_by_iso3c(WID_all, 'PAN')
    PRI_WID$iso3 = 'PRI'
    WID %<>% rbind(PRI_WID)
    missing_iso3c <- missing_iso3c[-which(missing_iso3c=='PRI')]
  }
  if ('CYM' %in% missing_iso3c){ #Population is small so just use the WID data from Cuba
    CYM_WID <- filter_WID_by_iso3c(WID_all, 'CUB')
    CYM_WID$iso3 = 'CYM'
    WID %<>% rbind(CYM_WID)
    missing_iso3c <- missing_iso3c[-which(missing_iso3c=='CYM')]
  }
  if ('GRD' %in% missing_iso3c){ #Population/exposed region is small so just use WID data from Venezuela 
    GRD_WID <- filter_WID_by_iso3c(WID_all, 'VEN')
    GRD_WID$iso3 = 'GRD'
    WID %<>% rbind(GRD_WID)
    missing_iso3c <- missing_iso3c[-which(missing_iso3c=='GRD')]
  }
  if ('VCT' %in% missing_iso3c){ #Population/exposed region is small so just use WID data from Venezuela 
    VCT_WID <- filter_WID_by_iso3c(WID_all, 'VEN')
    VCT_WID$iso3 = 'VCT'
    WID %<>% rbind(VCT_WID)
    missing_iso3c <- missing_iso3c[-which(missing_iso3c=='VCT')]
  }
  if ('SLB' %in%  missing_iso3c){
    WID %<>% rbind(SLB_WID)
    missing_iso3c <- missing_iso3c[-which(missing_iso3c=='SLB')]
  } 
  if ('VUT' %in%  missing_iso3c){
    WID %<>% rbind(VUT_WID)
    missing_iso3c <- missing_iso3c[-which(missing_iso3c=='VUT')]
  } 
  if (length(missing_iso3c) > 0){
    stop(paste('No WID data for', missing_iso3c, 'for', year))
  }
  return(WID)
}


GetWID_perc<-function(perc,iso3c,year){
  
  if (year > 2021) year <- "2021" #currently no data past 2021
  
  # Note that 'j' refers to the income divided equally between spouses 
  # (only chosen because it is the dataset with largest number of entries)
  perc <- paste0('p',seq(0,90,10), 'p', seq(10,100,10))
  WID_all <-download_wid(indicators = "sptinc",years=as.character(year),perc = as.character(perc), pop = "j")
  
  # Filter by most popular variable (usually 'sptinc992j')
  # WID%<>%filter(variable==names(which.max(table(WID$variable))))%>%
  WID <- WID_all %>% filter(variable=='sptinc992j')%>%
    mutate(iso3=ifelse(country=='KS', 'KOS', convIso2Iso3(country)))%>%
    dplyr::select(-c(variable,country,year)) %>%
    filter(iso3%in%iso3c)
  
  WID$value[which(WID$value==0)] <- 0.0001
  
  if (!all(iso3c %in% WID$iso3)){
    missing_iso3c <- iso3c[which(!iso3c %in% WID$iso3)]
    WID <- getMissingWID(missing_iso3c, WID, WID_all)
  }
  #WID$value<-1-WID$value
  #names(WID)[names(WID)=="percentile"]<-"variable"
  #mins<-WID%>%group_by(iso3)%>%summarise(mins=min(value),.groups = 'drop_last')
  #for(iso3c in mins$iso3) WID$value[WID$iso3==iso3c & WID$variable!="p10p100"]%<>%subtract(mins$mins[mins$iso3==iso3c])
  
  return(WID)
  
}

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
      
    } else stop("Please try coordinates that have 2 columns (or rows) as [Longitude,Latitude] instead of as a matrix")
  
}

WeightedMeanRWI<-function(interpRWI,RWI){
  
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

KrigMeUp<-function(interpRWI,RWI){

  if(class(RWI)!="SpatialPointsDataFrame") stop("wrong RWI points class, check the raw RWI data")
  # Create a default distance length scale for the variogram, based on the gridsize (bounding box area)
  bbox<-array(bbox(RWI),dim=4)
  dist<-0.5*sqrt((bbox[3]-bbox[1])*(bbox[4]-bbox[2]))
  
  # Variogram for RWI Intensity
  print("Krigging Intensity")
  vario <- gstat::variogram(Intensity~1, RWI)
  fit <- gstat::fit.variogram(vario, model=gstat::vgm(1, "Sph", dist, 1))
  interpRWI$Intensity<-gstat::krige(Intensity ~ 1, RWI, interpRWI@coords, model=fit)
  # Variogram for RWI Error
  print("Krigging Error")
  vario <- gstat::variogram(Error~1, RWI)
  fit <- gstat::fit.variogram(vario, model=gstat::vgm(1, "Sph", dist, 1))
  interpRWI$Error<-gstat::krige(Error ~ 1, RWI, interpRWI@coords, model=fit)
  
  return(interpRWI)

}

intrpRWI<-function(RWI,coords,GPR=F){
  
  # What form are the coords in?
  interpRWI<-tryCatch({ChangeCoords(coords)}, error=function(e) NA)
  if(is.na(interpRWI)) stop("problem with the interpolation coordinates, try inputing as data.frame(Longitude=...,Latitude=...)")
    
  if(GPR){
    interpRWI%<>%KrigMeUp(RWI=RWI)
  } else{
    interpRWI%<>%WeightedMeanRWI(RWI=RWI)
  }
  
  return(interpRWI)
  
}

GetRWI<-function(filer=NULL,RWIfolder=NULL,iso=NULL,intrp=F,coords=NULL,gridz=F){
  
  # Find the folder location of RWI file
  if(is.null(RWIfolder)) RWIfolder<-getwd()
  # Housekeeping
  if(!endsWith(RWIfolder,"/")) RWIfolder<-paste0(RWIfolder,"/")
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Find the RWI file(s) @@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  # If coutnry ISO3C is provided, use this to extract the appropriate RWI file
  if(!is.null(iso)) {
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
  
  # Extract RWI function - into spatialpointsdataframe
  RWI<-ExtractRWI(fRWI)
  
  # Check if we need to interpolate the RWI onto a provided grid
  if(is.null(coords)) return(RWI)
  # Interpolate function if coords are given
  return(intrpRWI(RWI,coords,gridz))
  
}

# 
# GetWB<-function(param,iso,Idate){
#   
#   param<-str_to_upper(gsub(" ", "", param, fixed = TRUE))
#   listy<-list("INFRASTRUCTURE"="GCI.2NDPILLAR.XQ","GDP"="NY.GDP.MKTP.PP.KD","POPULATION"="SP.POP.TOTL",
#               "POPDENSITY"="EN.POP.DNST")
#   indicator<-listy[[param]]
#   
#   if(is.null(indicator)) {
#     
#     ndata<-tryCatch(WBcall(as.character(2005),as.character(AsYear(Sys.Date())),param),error = function(e) NULL)
#     # ndata<-tryCatch(wb_data(indicator = param, start_date =as.character(2005), end_date = as.character(AsYear(Sys.Date())),date_as_class_date = T)%>%
#     # transmute(iso3=iso3c,date=date_ct,value=value),error = function(e) NULL)
#     if(is.null(ndata)) stop(paste0("ERROR: WB indicator not found '",param,"' see GetWB in GetSocioEconomic.R for examples"))
#     
#   } else {
#     
#     ndata<-tryCatch(WBcall(as.character(2005),as.character(AsYear(Sys.Date())),indicator),error = function(e) NULL)
#     # ndata<-wb_data(indicator = indicator, start_date = as.character(2005), end_date = as.character(AsYear(Sys.Date())),date_as_class_date = T)%>%
#     # transmute(iso3=iso3c,date=date_ct,value=value)
#   }
#   # ndata%<>%filter(iso3==iso&!is.na(value))%>%mutate(day=as.numeric(date-min(ndata$date)))
#   # tmp<-data.frame(iso3=iso,date=Idate)
#   if(length(iso))
#     
#     if(length(ndata$value)==0L) return(NA)
#   if(length(ndata$value)==1L) {print("WARNING: one value for WB ",indicator,", country ",iso) ;return(ndata$value)}
#   
#   func = splinefun(x=ndata$day,y=ndata$value, method="natural",  ties = mean)
#   return(func(as.numeric(Idate-min(ndata$date))))
#   
# }
# 
# GetINFORMinfo<-function(dir,year=2010){
#   
#   file<-paste0(dir,"Demography_Data/SocioEconomic/INFORM2020_TREND_2010_2020_v039_ALL.xls")
#   INFORMiso<-read.xlsx(file,sheetName = "Sheet1",header = TRUE)
#   
#   listy<-c("Population","Corruption Perception Index","Net ODA received (% of GNI)","Access to electricity",
#            "Government Effectiveness","Income Gini coefficient","Lack of Coping Capacity Index",
#            "Physical Infrastructure", "Disaster Risk Reduction", "INFORM Risk Index", "Vulnerability Index",
#            "Socio-Economic Vulnerability","Aid Dependency",
#            "FTS Current year","Estimated HDI from GDP per capita","Human Develpment Index",
#            "Physical exposure to earthquake MMI VI (relative) - raw","Physical exposure to flood (relative) - raw",
#            "Physical exposure to tropical cyclone of Saffir-Simpson category 1 (relative) - raw",
#            "Physical exposure to tsunami (relative) - raw","People affected by droughts (relative) - raw")
#   
#   INFORMiso%<>%filter(IndicatorName %in% listy & INFORMYear>year) %>% droplevels()
#   
#   drops<-c("IndicatorId","IndicatorType")
#   colnames(INFORMiso)[colnames(INFORMiso)=="Iso3"]<-"iso3"
#   INFORMiso$IndicatorName<-plyr::revalue(INFORMiso$IndicatorName,c("U5M"="Under 5 Mortality"))
#   INFORMiso<-INFORMiso[ , !(names(INFORMiso) %in% drops)]
#   
#   listy<-c("Population (total)"="POP","Corruption Perception Index"="CPI","Net ODA received (% of GNI)"="ODA",
#            "Estimated HDI from GDP per capita"="HDIGDP","Human Develpment Index"="HDI",
#            "Government Effectiveness"="GovEff","Lack of Coping Capacity Index"="CC",
#            "Physical Infrastructure"="PhysInf", "Disaster Risk Reduction"="DRR", "INFORM Risk Index"="INFORM", "Vulnerability Index"="Vuln",
#            "Socio-Economic Vulnerability"="SEVuln","Aid Dependency"="AidDep",
#            "Access to electricity"="ElecAcc","Income Gini coefficient"="GINI", "FTS Current year"="FTS",
#            "Physical exposure to earthquake MMI VI (relative) - raw"="EQexp","Physical exposure to flood (relative) - raw"="FLexp",
#            "Physical exposure to tropical cyclone of Saffir-Simpson category 1 (relative) - raw"="TCexp",
#            "Physical exposure to tsunami (relative) - raw"="TSexp","People affected by droughts (relative) - raw"="DRexp")
#   
#   INFORMiso$IndicatorName<-plyr::revalue(INFORMiso$IndicatorName, listy)
#   pc<-c("ElecAcc","CPI")
#   INFORMiso$IndicatorScore[INFORMiso$IndicatorName%in%pc]<-INFORMiso$IndicatorScore[INFORMiso$IndicatorName%in%pc]/100
#   
#   mINFORM<-data.frame()
#   for (ind in unique(INFORMiso$IndicatorName)){
#     t1<-INFORMiso%>%filter(IndicatorName==ind)
#     for (iso in unique(INFORMiso$iso3)){
#       tmp<-t1%>%filter(iso3==iso)%>%arrange(INFORMYear)
#       if(length(tmp$IndicatorScore)>0){
#         mINFORM<-rbind(mINFORM,data.frame(iso3=iso,IndicatorName=ind,dScore=max(tmp$IndicatorScore,na.rm = T)-min(tmp$IndicatorScore,na.rm = T),
#                                           wmScore=weighted.mean(x = tmp$IndicatorScore,w = 1:length(tmp$IndicatorScore))))
#       }
#     }
#     
#   }
#   for (ind in unique(mINFORM$IndicatorName)){
#     p<-ggplot(filter(mINFORM,IndicatorName==ind),aes(wmScore))+geom_density(fill="black",alpha=0.3)+
#       ylab("Density")+xlab(names(listy[listy==ind]))
#     ggsave(paste0("Density_",ind,'.png'), plot=p,path = paste0(directory,'Plots/WorldBank'),width = 7,height = 5)
#   }
#   
#   return(mINFORM)
#   
#   # corrupt<-xls %>% filter(IndicatorName=="Corruption Perception Index")
#   
# }
# 
# iso2WB<-function(iso,WBiso,ind="GINI",Score="wmScore"){
#   val<-WBiso%>%filter(iso3==iso & IndicatorName==ind)%>%pull(Score)
#   if(length(val)==0) return(NA)
#   return(val)
# }
# Viso2WB<-unname(Vectorize(iso2WB,vectorize.args = c("iso")))

# GetAllSocioEconomic<-function(dir){
#   
#   # group all data by country
#   WB<-GetWBinfo(dir)
#   WIB<-GetWIB(dir)
#   
#   filer<-paste0(dir,"Demography_Data/SocioEconomic/")
#   save(dfSE,filer)
#   
#   return(dfSE)
#   
# }

# GetWIB<-function(dir,iso3=NULL){
#   
#   file<-paste0(dir,"Demography_Data/SocioEconomic/")
#   
#   if(!is.null(iso3)) {
#     
#     return(WIB)
#   }
#   
#   # OUTPUT ALL DATA
#   
# }



# FillWBGaps<-function(listy,CM){
#   
#   nGDP<-listy[[1]]
#   nPop<-listy[[2]]
#   nPDens<-listy[[3]]
#   
#   interpy<-data.frame()
#   for (yr in unique(AsYear(CM$sdate))){
#     
#     isos<-CM%>%filter(AsYear(sdate)==yr)%>%pull(iso3)%>%unique()
#     
#     iGDP<-nGDP%>%filter(iso3%in%isos)
#     iGDP<-iGDP$iso3[is.null(iGDP[[as.character(yr)]])]
#     iPop<-nPop%>%filter(iso3%in%isos)
#     iPop<-iPop$iso3[is.null(iPop[[as.character(yr)]])]
#     iPDens<-nPDens%>%filter(iso3%in%isos)
#     iPDens<-iPDens$iso3[is.null(iPDens[[as.character(yr)]])]
#     
#     tmp<-data.frame()
#     if(length(iGDP)>0) tmp<-rbind(tmp,data.frame(iso3=iGDP,year=rep(yr,length(iGDP)),data=rep("GDP",length(iGDP))))
#     if(length(iPop)>0) tmp<-rbind(tmp,data.frame(iso3=iPop,year=rep(yr,length(iPop)),data=rep("Pop",length(iPop))))
#     if(length(iPDens)>0) tmp<-rbind(tmp,data.frame(iso3=iPDens,year=rep(yr,length(iPDens)),data=rep("PDens",length(iPDens))))
#     interpy%<>%rbind(tmp)
#     
#   }
#   rm(tmp)
#   
#   if(length(interpy)==0) return(listy)
#   
#   stop("Interpolation of WB data not yet done")
#   
# }




# GetKummu<-function(dir,bbox=NULL,yr=2015L){
#   
#   if(yr==2015L){iii<-3} else if (yr==2000L){iii<-2} else if(yr==1990L){iii<-1} else {stop("ERROR: incorrect year for GDP-PPP data")}
#   
#   # file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_PPP_30arcsec_v3.nc")
#   file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc")
#   if(!file.exists(file)) stop("no file found Kummu GDP")
#   ncin <- nc_open(file)
#   lat <- ncvar_get(ncin,"latitude")
#   nlat<- dim(lat)
#   lon <- ncvar_get(ncin,"longitude")
#   nlon<-dim(lon)
#   time <- ncvar_get(ncin,"time")
#   ntime<-dim(time)
#   
#   # GDP_PPP[longitude,latitude,time]
#   # GDP<-ncvar_get(ncin,"GDP_PPP",start = c(1,1,iii),count = c(nlon,nlat,1),collapse_degen = T)
#   GDP<-ncvar_get(ncin,"GDP_per_capita_PPP",start = c(1,1,iii),count = c(nlon,nlat,1),collapse_degen = T)
#   
#   nc_close(ncin)
#   
#   # GDP_PPP[longitude,latitude,time]
#   # GDP<-GDP[,,iii]
#   
#   colnames(GDP)<-lat
#   rownames(GDP)<-lon
#   
#   if(is.null(bbox)) {return(GDP)} else {GDP<-FilterKummu(GDP,bbox)}
#   
#   GDP%<>%convMat2SPDF(name="GDP")
#   # GDP_PPP[longitude,latitude]
#   return(GDP)
#   
# }