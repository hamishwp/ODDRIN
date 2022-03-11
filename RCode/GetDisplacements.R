
ExtractCheckedDisp<-function(dir,haz="EQ"){
  
  DispData<-read.csv(paste0(dir,"IIDIPUS_Input/DispData_",haz,".csv"),na.strings = "-")
  DispData$sdate%<>%as.Date()
  DispData$fdate%<>%as.Date()
  return(DispData)
  
}

GetDisplacements<-function(haz, saved=T, reduce=T, GIDD=T, EMDAT=F, dir="./"){
  
  if(saved) return(ExtractCheckedDisp(dir))
  
  #CM<-GetHelix(haz=haz,reduce=reduce)
  
  Dispy <- readRDS(paste0(dir,"IIDIPUS_Input/DispData_EQ_V2.Rdata"))
  CM <- Dispy %>% transmute(
    iso3 = iso3, 
    gmax = gmax, 
    hazard = 'EQ',
    sdate = as.Date(sdate),
    fdate = as.Date(fdate),
    eventid = event_id,
    qualifierDisp = 'total')
  
  if(GIDD){
    # Ensure we don't access anything before 2018
    CM%<>%filter(AsYear(sdate)>=2017)
    # Extract older facts from GIDD (Global Internal Displacement Database)
    GIDD<-GetGIDD(dir,haz)
    # Combine Helix and GIDD data into one database
    CM<-MergeGIDD_Helix(GIDD,CM)
    # Ensure distinct events
    CM%<>%distinct
    # Remove what we have already extracted
    rm(GIDD)
  }
  if(EMDAT){
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
    
    CM = mergeEMDAT_Helix(EMDAT, CM)
    rm(EMDAT)
    
  }
  
  return(CM)
}

mergeEMDAT_Helix <- function(EMDAT, CM){
  #cannot figure out a more efficient way of doing this using dplyr/joins
  for (i in 1:NROW(CM)){
    #find a match between the CM and the EMDAT data, using a three day window?
    match_indexes = which((EMDAT[,'iso3'] == CM[i,'iso3']) & 
                            (EMDAT[,'sdate'] > (CM[i, 'sdate'] - 3)) & 
                            (EMDAT[,'fdate'] < (ifelse(is.na(CM[i, 'fdate']), CM[i, 'sdate'] + 3, CM[i, 'fdate']))))
    
    if (length(match_indexes) > 0){
      if (length(match_indexes) > 1){ #if more than one match, take the closest date
        #EMDAT rows 286 and 287 cannot be distinguished with the exception of eventid (China // 2019-09-08)
        match_index = match_indexes[which.min(abs(EMDAT[match_indexes,'sdate']-CM[i, 'sdate']))]
      } else match_index = match_indexes
      EMDAT[match_index,'inHelix'] = TRUE
      match_eventids = which(EMDAT[,'eventid'] == EMDAT[match_index,'eventid'])
      EMDAT[match_eventids, 'eventid'] = CM[i,'eventid'] #change event id for all events in EMDAT with the same id
      CM[i,'mortality'] = EMDAT[match_index, 'mortality']
      CM[i, 'qualifierMort'] = EMDAT[match_index, 'qualifierMort']
      #need to think about which start date we use if they clash.
      CM[i,'fdate'] = if_else(CM[i,'sdate'] == EMDAT[match_index, 'sdate'], EMDAT[match_index, 'fdate'], CM[i,'fdate'])
    }
  }
  
  CM %<>% rbind(EMDAT %>% filter(inHelix == FALSE) %>% transmute(
    iso3 = iso3, 
    gmax = NA,
    qualifierDisp = NULL,
    hazard = hazard, 
    sdate = sdate, 
    fdate = fdate, 
    mortality = mortality,
    eventid = gsub("-", "", eventid) %>% as.numeric(),
    qualifierMort = qualifierMort))
  
  return(CM)
}
