library(dplyr)
library(magrittr)

ExtractData<-function(haz="EQ",dir="./",extractedData=T){
  
  if(extractedData) return(paste0(dir,"IIDIPUS_Input/ODDobjects/"))
  
  # Get the human displacement data from IDMC Helix & GIDD databases and other sources filtered by hazard
  DamageData<-GetDisplacements(haz, saved=F, GIDD=F, EMDAT=T)
  # Extract GDACS database on the event for further validation & alertscore benchmarking
  dfGDACS<-FilterGDACS(haz=haz,syear=min(AsYear(DamageData$sdate)),fyear=max(AsYear(DamageData$sdate)),red=T)
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (ev in unique(DamageData$eventid)){
    # Subset displacement and disaster database objects
    miniDam<-DamageData%>%filter(eventid==ev)
    # Set some confining dates for the rest of the data to be assigned to this event
    maxdate<-miniDam$sdate-5
    if(is.na(miniDam$fdate)) mindate<-miniDam$sdate+3 else mindate<-miniDam$fdate+3
    # GDACS subset
    miniDACS<-dfGDACS%>%filter(iso3%in%unique(miniDam$iso3) & 
                                 sdate<mindate & sdate>maxdate)
    # Match displacement and hazard data and extract hazard maps
    # HazSDF includes SpatialPixelDataFrame object of hazmean & hazsd per date 
    # (list of each, bbox-cropped to remove M < minmag)
    #lhazSDF<-tryCatch(GetDisaster(miniDam,miniDACS),error=function(e) NULL)
    lhazSDF<-tryCatch(GetDisaster(miniDam),error=function(e) NULL)
    if(!is.null(lhazSDF)) {
      print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
                   " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
      next
    }
    
    # Create the ODD object:
    ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDam),error=function(e) NULL)
    if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDam$iso3)[1],
                  "_",ODDy@eventid)
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_Input/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    HAZARDpath<-paste0(dir,"IIDIPUS_Input/HAZARDobjects/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
    
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDam$iso3) & 
                               sdate<mindate & sdate>maxdate)
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input/BDobjects/",namer)
      # Save it out!
      saveRDS(BDy, BDpath)
    }
    
    # Path to file and ID
    path%<>%rbind(data.frame(ODDpath=ODDpath,
                             BDpath=BDpath,
                             eventid=ODDy@eventid))
    # Save some RAM
    rm(ODDy,BDy,miniDam)
  }
  
  return(path)
  
}

ExtractNewData<-function(isos,sdat,fdat=NULL,haz="EQ",dir="./"){
  # Make sure the event dates exist and are in the correct format
  sdat%<>%as.Date()
  if(is.null(fdat)) if(Sys.Date()-sdat<5) fdat<-Sys.Date() else fdat<-sdat+5
  fdat%<>%as.Date()
  # Form the right objects to then extract the hazard data
  DamageData<-data.frame(iso3=isos,
                         hazard=haz,
                         sdate=sdat,
                         fdate=fdat,
                         eventid=999)
  # Extract GDACS database on the event for further validation & alertscore benchmarking
  dfGDACS<-FilterGDACS(haz=haz,syear=AsYear(DamageData$sdate),fyear=AsYear(DamageData$sdate),red=T)%>%
    filter(sdate>=sdat-7 & fdate<=fdat & iso3%in%isos)
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (ev in unique(DamageData$eventid)){
    # Subset displacement and disaster database objects
    miniDam<-DamageData%>%filter(eventid==ev)
    # Set some confining dates for the rest of the data to be assigned to this event
    maxdate<-miniDam$sdate-5
    if(is.na(miniDam$fdate)) mindate<-miniDam$sdate+3 else mindate<-miniDam$fdate+3
    # GDACS subset
    miniDACS<-dfGDACS%>%filter(iso3%in%unique(miniDam$iso3) & 
                                 sdate<mindate & sdate>maxdate)
    # Match displacement and hazard data and extract hazard maps
    # HazSDF includes SpatialPixelDataFrame object of hazmean & hazsd per date 
    # (list of each, bbox-cropped to remove M < minmag)
    #lhazSDF<-tryCatch(GetDisaster(miniDam,miniDACS),error=function(e) NULL)
    lhazSDF<-tryCatch(GetDisaster(miniDam),error=function(e) NULL)
    if(is.null(lhazSDF)) {
      print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
                   " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
      next
    }
    
    # Create the ODD object:
    ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DamageData=miniDam),error=function(e) NULL)
    if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDam$iso3)[1],
                  "_",ODDy@eventid)
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_Input/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    HAZARDpath<-paste0(dir,"IIDIPUS_Input/HAZARDobjects/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
    
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDam$iso3) & 
                               sdate<mindate & sdate>maxdate)
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input/BDobjects/",namer)
      # Save it out!
      saveRDS(BDy, BDpath)
    }
    
    # Path to file and ID
    path%<>%rbind(data.frame(ODDpath=ODDpath,
                             BDpath=BDpath,
                             eventid=ODDy@eventid))
    # Save some RAM
    rm(ODDy,BDy,miniDam)
  }
  
  return(path)
  
}

# Calculate the predicted displaced populations per gridpoint and save out to a file
FormODDyOmega<-function(dir,Model,proposed,AlgoParams){
  
  output<-data.frame()
  # Load ODD files
  ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/ODDobjects"),pattern=Model$haz,recursive = T,ignore.case = T)
  for(i in 1:length(ufiles)){
    ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/",ufiles[i]))
    ODDy@fIndies<-Model$fIndies
    if(class(ODDy@gmax)=="list") ODDy@gmax%<>%as.data.frame.list()
    # Apply DispX
    ODDy<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams),
                   error=function(e) NA)
    # if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));next}
    if(is.na(ODDy)) stop(paste0("Failed to calculate Disp LL of ",ufiles[i]))
    
    saveRDS(ODDy,paste0(dir,"IIDIPUS_Results/ODDobjects/",ufiles[i]))
    
    print(ODDy@predictDisp)
    output%<>%rbind(ODDy@predictDisp)
  }
  return(output)
  
}

# # Total LL MAP:
# vals<-c(0.02337845,0.03644427,0.06072838,1.61965282,4.70430026,0.59366736,0.11988817,3.81639425)
# 
# # LL_Disp MAP:
# vals<-c(0.01080671,0.08229228,0.12590510,2.28739344,3.53698688,0.52015751,0.13003755,6.02350025)
# 
# Omega$Lambda[1:3]<-vals[1:3]
# Omega$zeta[1:2]<-vals[4:5]
# Omega$theta[1]<-vals[6]
# Omega$eps[1:2]<-vals[7:8]
# 
# FormODDyOmega(dir,Model,Omega,AlgoParams)
  

ODDypreds<-function(dir,haz,Model,Omega,AlgoParams){
  
  output<-data.frame()
  # Load ODD files
  ufiles<-list.files(path=paste0(dir,"IIDIPUS_Results/ODDobjects/"),pattern=haz,recursive = T,ignore.case = T)
  for(i in 1:length(ufiles)){
    # if(grepl(ufiles[i],pattern = "CHN")) next
    ODDy<-tryCatch(readRDS(paste0(dir,"IIDIPUS_Results/ODDobjects/",ufiles[i])),error=function(e) NA)
    ODDy@gmax%<>%arrange(desc(gmax))
    mody<-ODDy@modifier
    repmod<-rep(0,length(ODDy@modifier));names(repmod)<-ODDy@gmax$iso3; repmod%<>%as.list()
    ODDy@modifier<-repmod
    
    ttt<-data.frame(iso3=unique(ODDy@cIndies$iso3))
    # for (var in unique(ODDy@cIndies$variable)){
    var <-"p50p100"
    vt<-ODDy@cIndies[ODDy@cIndies$variable==var,c("iso3","value")]
    names(vt)[2]<-var
    ttt%<>%merge(vt,by="iso3")
    GDP<-ODDy@data%>%group_by(ISO3C)%>%summarise(meany=mean(GDP,na.rm=T),.groups="keep")%>%na.omit()%>%transmute(iso3=ISO3C,GDP=meany)
    ttt%<>%merge(GDP,by="iso3")
    
    ODDy<-tryCatch(DispX(ODD = ODDy,Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL = F,Method = AlgoParams),
                   error=function(e) NA)
    if(is.na(ODDy)|nrow(ODDy@predictDisp)<1) {print(paste0("no information found for ",ufiles[i])) ; next}
    
    tmp<-cbind(ODDy@predictDisp,data.frame(eventid=extractnumbers(ufiles[i])[2],namer=ufiles[i]))
    
    ttt$maxH<-max(ODDy$hazMean1,na.rm = T)
    ttt$aff6<-sum(ODDy$Population[ODDy$hazMean1>6],na.rm = T)
    
    ODDy@modifier<-mody
    saveRDS(ODDy,paste0(dir,"IIDIPUS_Results/ODDobjects/",ufiles[i]))
    
    output%<>%rbind(merge(tmp,ttt,by="iso3"))
    print(ODDy@predictDisp)
  }
  
  return(output)
  
}  
  
# ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/ODDobjects"),pattern=Model$haz,recursive = T,ignore.case = T)
GenerateBDs<-function(ufiles){
  for(f in ufiles){
    miniDam<-Damage%>%filter(grepl(strsplit(f,"_")[[1]][1],as.character(event)))
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      print(f)
      # Load the ODD object
      ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/",f))
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0(".........BD FAIL:     ",f)) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input/BDobjects/",f)
      BDy@cIndies<-WID%>%filter(year==AsYear(BDy@hazdates[1]) & 
                                  iso3%in%unique(BDy@data$ISO3C))%>%
        dplyr::select(-year)
      # Save it out!
      saveRDS(BDy, BDpath)
    }
    
  }
}
# 
# bfiles<-list.files(path=paste0(dir,"IIDIPUS_Input/BDobjects"),pattern=Model$haz,recursive = T,ignore.case = T)
# fff<-unlist(strsplit(bfiles,"_"))[seq.int(from=1,to=2L*length(bfiles)-1L,by = 2)]

# ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/BDobjects"),pattern=Model$haz,recursive = T,ignore.case = T)
SortBDout<-function(ufiles){
  for(f in ufiles){
    BDpath <-paste0(dir,"IIDIPUS_Input/BDobjects/",f)
    BDy<-readRDS(BDpath)
    BDy@cIndies<-WID%>%filter(year==AsYear(BDy@hazdates[1]) & 
                                iso3%in%unique(BDy@data$ISO3C))%>%
      dplyr::select(-year)
    print(paste0(unique(BDy@cIndies$iso3)," with ",unique(BDy@data$ISO3C)))
    # Save it out!
    saveRDS(BDy, BDpath)
    
  }
}

#update the ODD files to include mortality data from EMDAT
UpdateODD <- function(){
  
  startRow = 7 #ignore metadata stored on the first 6 rows
  EMDAT = openxlsx::read.xlsx(paste0(dir,"/Displacement_Data/emdat.xlsx"), startRow=startRow)
  EMDAT %<>% transmute(
    iso3 = ISO,
    mortality = ifelse(is.na(Total.Deaths), 0, Total.Deaths), #assume blank cells correspond to 0 deaths ??? 
    qualifierMort = 'total',
    hazard = ifelse(Disaster.Type == 'Earthquake', 'EQ', 'UK'), #Earthquake or Unknown
    sdate = as.Date(paste(Start.Day,Start.Month,Start.Year, sep='-'), format='%d-%m-%Y'),
    eventid = paste(Year, Seq, sep='-'), 
    fdate = as.Date(paste(Start.Day,Start.Month,Start.Year, sep='-'), format='%d-%m-%Y'), 
    inHelix = FALSE #tracks to avoid duplicates
  )
  
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
    
  for(i in 1:length(ufiles)){
    # Extract the ODD object
    ODDy<-readRDS(paste0(folderin,ufiles[i]))
    
    event_match <- which((EMDAT[,'iso3'] %in% ODDy@gmax$iso3 & 
            (EMDAT[,'sdate'] > (min(as.Date(ODDy@hazdates)) - 2)) & 
            (EMDAT[,'fdate'] < (max(as.Date(ODDy@hazdates)) + 2))))
    if (length(event_match) > 0){
      ODDy@gmax['mortality'] <- NA
      ODDy@gmax['qualifierMort'] <- 'Total'
      for (k in 1:length(ODDy@gmax$iso3)){
        country_match = which(EMDAT[event_match,'iso3'] == ODDy@gmax$iso3[k])
        if (length(country_match) > 1){
          print(paste('Matching EQ in', ODDy@gmax$iso3[k], 'from', min(as.Date(ODDy@hazdates)), 'to',  max(as.Date(ODDy@hazdates)), 
                      'with EQs in', EMDAT[event_match[country_match], 'iso3'], 'from', EMDAT[event_match[country_match], 'sdate'], 'to', EMDAT[event_match[country_match], 'fdate']))
          #if more than one match within the time period, take the sum of the mortality
          ODDy@gmax[k, 'mortality'] <- sum(EMDAT[event_match[country_match], 'mortality'])
        }
        else if (length(country_match) == 1){
          print(paste('Matching EQ in', ODDy@gmax$iso3[k], 'from', min(as.Date(ODDy@hazdates)), 'to',  max(as.Date(ODDy@hazdates)), 
                'with EQ in', EMDAT[event_match[country_match], 'iso3'], 'from', EMDAT[event_match[country_match], 'sdate'], 'to', EMDAT[event_match[country_match], 'fdate']))
          ODDy@gmax[k, 'mortality'] <- EMDAT[event_match[country_match], 'mortality']
        }
      }
      
      saveRDS(ODDy, paste0(folderin, ufiles[i]))
    }
  }
}
