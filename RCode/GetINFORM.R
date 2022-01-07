library(httr)

GetINFORMid_year<-function(year){
  
  url<-paste0("https://drmkc.jrc.ec.europa.eu/inform-index/API/InformAPI/Workflows/GetByYear/",year)
  robj<-tryCatch(GET(url = url),error = function(e) NULL)
  if (robj$status_code!=200) stop(paste0("ERROR: INFORM id lookup table not found for ",year))
  
  return(content(robj,as = "parsed")[[1]][["WorkflowCompareId"]])
  
}

GetINFORMdata<-function(indicator,year,iso=NULL,normalise=F){
  # indicator= {INFORM, CC, VU, VU.SEV, CC.INF.PHY, } 
  # For more information, visit https://drmkc.jrc.ec.europa.eu/inform-index/In-depth/API-Demo
  
  id<-GetINFORMid_year(year)
  indicator<-str_to_upper(gsub(" ", "", indicator, fixed = TRUE))
  
  url<-paste0("https://drmkc.jrc.ec.europa.eu/inform-index/API/InformAPI/countries/Scores/?WorkflowId=",id,"&IndicatorId=",indicator)
  
  robj<-tryCatch(GET(url = url),error = function(e) NULL)
  if (robj$status_code!=200) stop(paste0("ERROR: INFORM id lookup table not found for ",year))
  
  tmp<-content(robj,as = "parsed")
  if(length(tmp)==0) stop(paste0("ERROR: WB indicator not found '",indicator,"' see GetINFORMdata in GetINFORM.R for examples"))
  
  data<-as.data.frame(t(matrix(unlist(tmp), nrow=length(unlist(tmp[[1]]))))); data<-data[,c(1,3)]; colnames(data)<-c("iso3",indicator)
  data[,2]<-as.numeric(data[,2])
  
  if(normalise) data[,2]<-data[,2]/10
  
  if(is.null(iso)) return(data)
  return(filter(data,iso3==iso))
  
}

INFORMmean<-function(indicator,iso=NULL){
  
  INFORM<-data.frame(iso3=character(),variable=character(),value=numeric())
  for (id in c(360,369,370,384,386)){
    for(indy in unique(indicator)){
      # Get 2020 value of the INFORM parameter
      url<-paste0("https://drmkc.jrc.ec.europa.eu/inform-index/API/InformAPI/countries/Scores/?WorkflowId=",id,"&IndicatorId=",indy)
      robj<-tryCatch(GET(url = url),error = function(e) NULL)
      if (robj$status_code!=200) stop(paste0("ERROR: INFORM id lookup table not found for 2019"))
      
      tmp<-content(robj,as = "parsed")
      if(length(tmp)==0) stop(paste0("ERROR: WB indicator not found '",indy,"' see GetINFORMdata in GetINFORM.R for examples"))
      
      # I know this is ugly... I can't be bothered to re-write it. Sorry!
      data<-as.data.frame(t(matrix(unlist(tmp), nrow=length(unlist(tmp[[1]]))))); data<-data[,c(1,3)]; colnames(data)<-c("iso3","value")
      data[,2]<-as.numeric(data[,2])
      
      if(!is.null(iso)) data%<>%filter(iso3%in%iso)
      INFORM%<>%rbind(cbind(data,data.frame(variable=rep(indy,nrow(data)))))
      
    }
    
    INFORM$value<-INFORM$value/10
  }
  
  INFORM%<>%group_by(iso3,variable)%>%summarise(value=mean(value,na.rm=T))
  INFORM<-INFORM[,c(1,3,2)]
  
  return(INFORM)
}

INFORM_LM<-function(iso,Indicators,DispData){
  # I ASSUME ONLY CC.INS.DRR is missing
  # I know... I'm sorry. So ugly
  ind<-Indicators$iso3%in%unique(DispData$iso3[DispData$iso3%in%unique(Indicators$iso3[Indicators$iso3%in%unique(Indicators$iso3[Indicators$variable=="CC.INS.DRR"])])])
  DF<-data.frame(DRR=Indicators$value[Indicators$variable=="CC.INS.DRR" & ind],
                 PHY=Indicators$value[Indicators$variable=="CC.INF.PHY" & ind],
                 GOV=Indicators$value[Indicators$variable=="CC.INS.GOV.GE" & ind],
                 EQ=Indicators$value[Indicators$variable=="HA.NAT.EQ" & ind],
                 VU=Indicators$value[Indicators$variable=="VU.SEV.AD" & ind],
                 PD=Indicators$value[Indicators$variable=="VU.SEV.PD" & ind])
  
  nind<-Indicators$iso3%in%unique(DispData$iso3[DispData$iso3%in%unique(Indicators$iso3[!Indicators$iso3%in%unique(Indicators$iso3[Indicators$variable=="CC.INS.DRR"])])])
  preddy<-data.frame(PHY=Indicators$value[Indicators$variable=="CC.INF.PHY" & nind],
                     GOV=Indicators$value[Indicators$variable=="CC.INS.GOV.GE" & nind],
                     EQ=Indicators$value[Indicators$variable=="HA.NAT.EQ" & nind],
                     VU=Indicators$value[Indicators$variable=="VU.SEV.AD" & nind],
                     PD=Indicators$value[Indicators$variable=="VU.SEV.PD" & nind])
  
  fit<-lm(data=DF,DRR~PHY+GOV+EQ+VU+PD)
  predictor<-predict.lm(fit,newdata = preddy)
  
  Indicators%<>%as.data.frame()
  Indicators%<>%rbind(data.frame(iso3=unique(Indicators$iso3[nind]),value=unname(predictor),variable=rep("CC.INS.DRR",length(predictor))))
  
  return(Indicators)
  
}

InterpINFORMdata<-function(indicator,dater,iso=NULL){
  # indicator= {INFORM, CC, VU, VU.SEV, CC.INF.PHY, ...} 
  # For more information, visit https://drmkc.jrc.ec.europa.eu/inform-index/In-depth/API-Demo
  
  # NOTE: 2016mid does not seem to exist...
  WorkflowId<-data.frame(year=c(2015,2015.5,2016,2017,2017.5,2018,2018.5,2019,2019.5,2020), # This is the year stated for the INFORM index
                         ID=c(193,211,258,261,354,360,369,370,384,386), # The Workflow ID used to access this particular year
                         Ddate=as.Date(c("2014-10-22","2015-03-05","2016-05-03","2016-08-25",
                                         "2017-03-13","2017-08-24","2018-03-28","2018-08-30",
                                         "2019-03-29","2019-06-21"))) # The actual dates the data was released
  mindate<-min(WorkflowId$Ddate)
  WorkflowId%<>%mutate(day=as.numeric((Ddate-mindate)))%>%dplyr::select(-c(Ddate))
  
  INFORM<-data.frame(iso3=character(),variable=character(),value=numeric())
  for(indy in unique(indicator)){
    IndData<-data.frame(iso3=character(),value=numeric(),day=numeric())
    for(id in unique(WorkflowId$ID)){
      
      url<-paste0("https://drmkc.jrc.ec.europa.eu/inform-index/API/InformAPI/countries/Scores/?WorkflowId=",id,"&IndicatorId=",indy)
      
      robj<-tryCatch(GET(url = url),error = function(e) NULL)
      if (robj$status_code!=200) stop(paste0("ERROR: INFORM id lookup table not found for ",year))
      
      tmp<-content(robj,as = "parsed")
      if(length(tmp)==0) stop(paste0("ERROR: WB indicator not found '",indy,"' see GetINFORMdata in GetINFORM.R for examples"))
      
      # I know this is ugly... I can't be bothered to re-write it. Sorry!
      data<-as.data.frame(t(matrix(unlist(tmp), nrow=length(unlist(tmp[[1]]))))); data<-data[,c(1,3)]; colnames(data)<-c("iso3","value")
      data[,2]<-as.numeric(data[,2])
      
      if(!is.null(iso)) data%<>%filter(iso3%in%iso)
      data$day=rep(WorkflowId$day[WorkflowId$ID==id],nrow(data))
      
      IndData%<>%rbind(data)
    }
    
    # Interpolate for actual day value
    intp<-InterpDay(ndata = IndData,day = as.numeric(dater-mindate))
    INFORM%<>%rbind(cbind(intp,data.frame(variable=rep(indy,nrow(intp)))))
    
  }
  
  return(INFORM)
  
}
