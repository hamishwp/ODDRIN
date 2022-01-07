#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%% IIDIPUS - Integrated Internal DIsplaced PopUlation Sampler %%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%% Coded by Hamish Patten %%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%% Collaboration between the University of Oxford %%%%%%%%%%#
#%%%%%%%%%% and the Internal Displacement Monitoring Centre %%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%% Started February 2020 %%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
source('RCode/Functions.R')
# Disaster related:
source('RCode/GetGDACS.R')
# source('RCode/GetDisaster.R')
# IDP estimate related:
source('RCode/GetHelix.R')
# Demography & population related:
source('RCode/GetPopDemo.R')
source('RCode/GetSocioEconomic.R')
source('RCode/GetINFORM.R')
# # Mobility related:
# source('RCode/GetOSM.R')
# source('RCode/GetFBconnectivity.R')
# # Damage estimate related:
# source('RCode/GetBuildingDamageEstimates.R')

library(tidyverse)
# detach(package:plyr)
library(dplyr)
library(magrittr)
library(reshape2)

IIDIPUS_Past<-function(dir,haz="Earthquake"){
  
  print("WARNING: Currently unable to handle event occuring in multiple countries")
  
  filer <- paste0(dir,"Statistical_Inference/past_data_calc",Sys.Date(),".Rdata")
  
  ############################ DISPLACEMENT DATA ############################
  # Helix<-GetHelix_F(dir) # Filtered Helix results for statistical inference
  CM<-GetHelix(haz=haz,reduce=T)
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
  # CM%<>% filter(!(eventid %in% c(-93,-7,,23,-43,-60,-74,-121,-126,-130,
  #                                    -131,6211,6588,1475,-59,-75,6754,-101,
  #                                    -108,-46,-32,-82,-103,-110,2195,2588,
  #                                    2587,3614,4506,6885,7263,7876,8870,-88,
  #                                    1479,1544,1549,1545,2216)))
  ################################################################
  
  ########################### DISASTER DATA ###########################
  # convert Helix hazard name into GDACS - e.g. "Earthquake" into "EQ"
  hazard<-GDACSconvDis(haz)
  # What years should we cover?
  syear<-min(AsYear(CM$sdate),na.rm = T)
  fyear<-max(AsYear(CM$sdate),na.rm = T)
  # Extract all relevant events from GDACS according to what is in IDMC database (GIDD & Helix)
  dfGDACS<-FilterGDACS(dir,haz=hazard,syear = syear,fyear = fyear)
  # Extract shakemaps/cyclone paths/flood zones, etc
  polyframe<-Match_HelixGDACS(dir,CM,hazard,dfGDACS)
  # Remove all Helix events that do not have shakemaps
  CM%<>%filter(eventid%in%unique(polyframe$helix_id))
  # Match the hazard severities in Helix DB with dfGDACS
  severity<-polyframe%>%filter(ncontour==0)%>%dplyr::select(!c(id,ncontour))%>%
    distinct()%>%distinct(across(helix_id),.keep_all = T)
  stmp<-CM%>%group_by(eventid)%>%transmute(hazard_severity=severity$Intensity[severity$helix_id%in%eventid],
                                           hazard_date=severity$date[severity$helix_id%in%eventid],
                                           hazard_alert=severity$alertscore[severity$helix_id%in%eventid])
  CM%<>%mutate(hazard_severity=stmp$hazard_severity,hazard_date=stmp$hazard_date,hazard_alert=stmp$hazard_alert)
  rm(severity,stmp)
  
  
  # Cheeky plot can't hurt!
  # p<-ggplot(CM,aes(x=(day+1),y=IDP/gmax,group=eventid))+geom_point()+scale_x_log10() +
  #   xlab("Hazard Occurance Date")+ylab("IDP Stock / max(count) ") + ylim(c(0,1));p
  # ggsave("CM_gmax_day.eps", plot=p,path = paste0(directory,'Plots/'),width = 6,height = 5)
  
  p<-ggplot(CM,aes(x=hazard_severity,y=gmax))+geom_point()+scale_y_log10()+
    xlab("Hazard Severity")+ylab("Maximum Population Displacement");p
  ggsave("CM_gmax_severity.eps", plot=p,path = paste0(directory,'Plots/'),width = 6,height = 5)
  
  p<-ggplot(CM,aes(x=hazard_alert,y=gmax))+geom_point()+scale_y_log10()+
    xlab("GDACS Alertscore")+ylab("Maximum Population Displacement");p         
  ggsave("CM_gmax_alertscore.eps", plot=p,path = paste0(directory,'Plots/'),width = 6,height = 5)
  
  
  # if(haz="Earthquake"){CM$hday<-(CM$sdate+CM$day-CM$hazard_date)}
  # CM$hday<-as.numeric(CM$sdate+CM$day-CM$hazard_date)
  # Choose intervals for polygon integration e.g. >4.5M, >5M, ..., Richter scale earthquakes
  IntMap<-GetIntMap(hazard)
  
  mnhaz<-4.5
  polyframe%<>%filter(Intensity>=mnhaz)
  # Reduce the IDMC and GDACS database into information that matches up between them!
  # Clean up / reduce CM
  # CM%<>%dplyr::select(-c(mxdate,Evac,Shelt,hazard_date))
  ################################################################
  
  ############## WORLD BANK DATA - NORMALISATION #################
  # nGDP<-GetWBGDP(dir,syear,fyear)
  nPop<-GetWBPop(dir,syear,fyear)
  nPDens<-GetWBPDens(dir,syear,fyear)
  ###############################################################
  
  ########################## GDP-PPP ############################
  GDP<-GetKummu(dir,yr=2015L)
  ###############################################################
  
  # SETUP TEMPORARY DATABASES
  lenIM<-length(IntMap)
  cIntMap<-as.character(IntMap)
  tmp <- data.frame(matrix(NA,ncol = lenIM, nrow = 1))
  
  PopCounts<-ZoneArea<-MaxGDP<-MaxPopDens<-data.frame(matrix(ncol = lenIM, nrow = 0))
  colnames(PopCounts)<-paste(cIntMap,"PopCounts",sep = "-")
  colnames(ZoneArea)<-paste(cIntMap,"ZoneArea",sep = "-")
  colnames(MaxGDP)<-paste(cIntMap,"MaxGDP",sep = "-")
  colnames(MaxPopDens)<-paste(cIntMap,"MaxPopDens",sep = "-")
  eventy<-c()
  # For machines with large amounts of RAM, we can load the population & pop-density data by country 
  # For GDP data, this extraction is through netCDF (long) and so we load it beforehand.
  # This will save on data extraction time (VERY SLOW FOR ~20GB DATASETS)
  for (iso in unique(CM$iso3)){
    
    subhel<-CM%>%filter(iso3==iso)
    # subdacs<-dfGDACS%>%filter(iso3==iso)
    polys<-polyframe%>%filter(helix_id%in%subhel$eventid)#|eventid%in%subdacs$eventid | )
    
    if(dim(polys)[1]==0L|max(polys$ncontour)==0) {
      
      event<-unique(subhel$eventid)
      cat(paste0("ncontour==0 only for event - ",event,"\n"))
      
      lenev<-length(event)
      
      ttt<-tmp %>% slice(rep(row_number(), lenev))
      
      PopCounts%<>%rbind(setNames(ttt,colnames(PopCounts)))
      ZoneArea%<>%rbind(setNames(ttt,colnames(ZoneArea)))
      MaxGDP%<>%rbind(setNames(ttt,colnames(MaxGDP)))
      MaxPopDens%<>%rbind(setNames(ttt,colnames(MaxPopDens)))
      eventy<-c(eventy,event)
      
      next
    }  
    
    # bbox - [mnlo,mnla,mxlo,mxla]
    bbox<-c(min(polys$Longitude,na.rm = T),min(polys$Latitude,na.rm = T),
            max(polys$Longitude,na.rm = T),max(polys$Latitude,na.rm = T))
    bbox<-CheckBbox(bbox)
    
    ############ POP-COUNT ############
    pCnt<-GetPopulationBbox(dir,bbox,density=FALSE,yr=2015L); pCnt<-melt(pCnt);colnames(pCnt)<-c("X","Y","data")
    ########### POP-DENSITY ###########
    pDens<-GetPopulationBbox(dir,bbox,density=TRUE,yr=2015L); pDens<-melt(pDens);colnames(pDens)<-c("X","Y","data")
    ############# GDP-PPP #############
    tGDP<-FilterKummu(GDP,bbox,melted = T)
    
    for (event in unique(subhel$eventid)){
      
      sshel<-subhel%>%filter(eventid==event)
      # Reduce shakemap polygons to overlap all activity per magnitude into one
      # Future work - separate id and ncontour values
      subpol<-tryCatch(polys%>%filter(helix_id==event)%>%ReducePolyEvent(),error = function(e) NULL)
      
      if(length(subpol)==0L|is.null(subpol)) {
        if(!is.null(subpol)) {print(paste0("No polygon info for event - ",event))
        } else print(paste0("Failed to gather polygon info for event - ",event))
        PopCounts%<>%rbind(rep(NA,lenIM))
        ZoneArea%<>%rbind(rep(NA,lenIM))
        MaxGDP%<>%rbind(rep(NA,lenIM))
        MaxPopDens%<>%rbind(rep(NA,lenIM))
        eventy<-c(eventy,event)
        
        next
      }
      
      ############ INTEGRATE DATA BY POLYGONS ############
      # Do mean and max of non-zero values only PER POLYGON (automatically done by PolyIntegrate)
      IntCounts<-PolyIntegrateData(pCnt,subpol,IntMap,sum)
      IntDens<-PolyIntegrateData(pDens,subpol,IntMap,max)
      IntGDP<-PolyIntegrateData(tGDP,subpol,IntMap,max)
      
      date<-min(sshel$sdate)
      
      ############ NORMALISE TO YEAR OF HAZARD ###########
      IntCounts$value<-ceiling(IntCounts$value*normaliseWB(nPop,iso,date))
      IntDens$value<-IntDens$value*normaliseWB(nPDens,iso,date)
      #IntGDP$value<-IntGDP$value*normaliseWB(nGDP,iso,date)
      
      ######### INTEGRATE INTO TEMPORARY DATABASE ########
      PopCounts<-rbind(PopCounts,IntCounts$value)
      ZoneArea<-rbind(ZoneArea,IntCounts$area)
      MaxGDP<-rbind(MaxGDP,IntGDP$value)
      MaxPopDens<-rbind(MaxPopDens,IntDens$value)
      eventy<-c(eventy,event)
      
    }
    
    saveRDS(object = list(PopCounts=PopCounts,ZoneArea=ZoneArea,MaxGDP=MaxGDP,MaxPopDens=MaxPopDens,eventy=eventy),file = paste0(dir,"IIDIPUS_Results/Historical_Analysis/OUTPUT_",hazard,"_",iso,".Rdata"))
    
  }
  rm(GDP,tGDP,pCnt,pDens)
  
  PopCounts<-listy$PopCounts
  ZoneArea<-listy$ZoneArea
  MaxGDP<-listy$MaxGDP
  MaxPopDens<-listy$MaxPopDens  
  # MaxGDP[is.infinite(MaxGDP)]<-0
  # MaxPopDens[is.infinite(MaxPopDens)]<-0
  eventy<-listy$eventy
  rm(listy)
  
  srt<-sort(eventy,decreasing = T,index.return=T)
  PopCounts<-PopCounts[srt$ix,];ZoneArea<-ZoneArea[srt$ix,];MaxGDP<-MaxGDP[srt$ix,];MaxPopDens<-MaxPopDens[srt$ix,]
  eventy<-srt$x; rm(srt)
  
  allnan<-apply(PopCounts, 1, function(x) all(is.na(x)))
  anynan<-apply(PopCounts, 1, function(x) any(is.na(x)))
  
  colnames(PopCounts)<-paste(cIntMap,"PopCounts",sep = "-")
  colnames(ZoneArea)<-paste(cIntMap,"ZoneArea",sep = "-")
  colnames(MaxGDP)<-paste(cIntMap,"MaxGDP",sep = "-")
  colnames(MaxPopDens)<-paste(cIntMap,"MaxPopDens",sep = "-")
  
  fCM<-CM%>%filter(eventid%in%eventy)
  subCM<-CM%>%filter(eventid%in%eventy[!anynan])
  
  ###################### WORLD BANK DATA #########################
  # Extract country specific 1-D information (e.g. physical exposure to hazard per country,
  # social capital, INFORM coping capacity index, etc)
  # WBiso<-GetWBinfo(dir)
  # PhysExp<-WBiso%>%filter(iso3==unique(subCM$iso3) & IndicatorName=="EQexp")%>%pull(wmScore)
  Ephys<-GetINFORMdata("HA.NAT.EQ",AsYear(Sys.Date()),normalise = T)
  print("WARNING: INFORM data based on current year, should be interpolated")
  Infra<-GetINFORMdata("CC.INF",AsYear(Sys.Date()),normalise = T)
  VU<-GetINFORMdata("VU",AsYear(Sys.Date()),normalise = T)
  INFORM<-GetINFORMdata("INFORM",AsYear(Sys.Date()),normalise = T)
  INCOME<-GetWB("NY.ADJ.NNTY.PC.CD",date=2015)
  
  ##################### GENERATE DATABASE ########################
  severity<-subCM%>%group_by(eventid)%>%summarise(haz=max(hazard_severity),gmax=max(gmax),IDP=max(IDP),iso3=unique(iso3)[1],alert=unique(hazard_alert)[1],hday=min(hday),day=min(day),sdate=min(sdate))
  severity<-Reduce(function(x,y) merge(x = x, y = y, by = "iso3"), list(severity,Ephys,VU,Infra,INFORM))
  
  severity%<>%arrange(desc(eventid))%>%filter(eventid%in%eventy)%>%cbind(PopCounts[!anynan,],MaxGDP[!anynan,4],MaxPopDens[!anynan,4])
  colnames(severity)[(dim(severity)[2]-1):dim(severity)[2]]<-c("MaxGDP_6M","MaxPopDens_6M")
  severity[["MaxPopDens_6M"]][is.infinite(severity[["MaxPopDens_6M"]])]<-0
  severity[["MaxGDP_6M"]][is.infinite(severity[["MaxGDP_6M"]])]<-0
  # severity%<>%filter(day<20)
  
  # severity%<>%filter(abs(hday)<15 & abs(day)<20)%>%select(-c(hday,day))
  if(hazard=="EQ") {
    severity%<>%select(-c("9-PopCounts"))
    cIntMap<-cIntMap[-lenIM]
    cIntMap_m<-as.character((IntMap-0.5)) ; cIntMap_m<-cIntMap_m[-lenIM]  
    IntMap<-IntMap[-lenIM]
  }
  # eventy<-unique(severity$eventid); iso<-unique(severity$iso3)
  # severity%<>%select(-c(eventid,iso3,IDP,haz,alert))
  
  library(GGally)
  
  subsev<-severity#%>%select(-c(IDP))
  
  # colnames(subsev)<-c(colnames(select(subsev,-c(paste(cIntMap,"PopCounts",sep = "-")))),
  # paste("Pop",cIntMap_m,sep="_"),paste0("ZoneArea_",cIntMap[1]))
  colnames(subsev)[14:21]<-paste("Pop",cIntMap_m,sep="_")
  subsev%<>%mutate(logmax=log(gmax+1))%>%
    select(-c(gmax))
  subsev%<>%filter(abs(hday)<30)%>%select(-c(hday,day))
  
  Imin<-4.5
  print("NOTE: we subtract 0.5 from EQ magnitude to account for the bizarre GDACS shakemap calculation")
  # tmp<-subsev%>%select(paste("Pop",cIntMap_m,sep="_"))%>%apply(1, function(x) log(1+x*exp((IntMap-0.5)/Imin)))%>%t()
  # tmp[tmp==0]<-NA
  # subsev[,paste("Pop",cIntMap_m,sep="_")]<-tmp
  subsev$haz<-exp(subsev$haz/Imin)
  subsev$MaxGDP_6M<-log(subsev$MaxGDP_6M+1)
  subsev$MaxPopDens_6M<-log(subsev$MaxPopDens_6M+1)
  subsev$MaxGDP_6M[subsev$MaxGDP_6M==0]<-NA
  subsev$MaxPopDens_6M[subsev$MaxPopDens_6M==0]<-NA
  subsev$IDP<-log(subsev$IDP)
  # colnames(subsev)<-c(colnames(select(subsev,-c(paste("Pop",cIntMap_m,sep="_")))),paste("Pop_exp",cIntMap_m,sep=""))
  # subsev$logPop[subsev$logPop==0]<-NA
  # subsev%<>%mutate(logPopDens5M=(logPop5M+logArea5M))%>%select(-logArea5M)
  subsev<-subsev[,c(dim(subsev)[2],1:(dim(subsev)[2]-1))]
  subsev%<>%mutate(Continent=countrycode::countrycode(iso3,"iso3c","continent"))#, INF_bool=as.factor(INFORM<0.4))
  
  
  
  # sevplot<-subsev%>%select(-c(paste("Pop",cIntMap_m,sep="_"),"iso3","sdate"))
  # sevplot<-subsev%>%select(-c("IDP","iso3","Continent","eventid","hday","day","sdate",paste("Pop",cIntMap_m[],sep="_")))
  
  # sevplot<-subsev%>%select(c(logmax,alert,INFORM,Pop_6))
  # sevplot$Pop_6[sevplot$Pop_6==0]<-NA
  # colnames(sevplot)<-c("log(Max. Disp.)","GDACS Alert","INFORM","log(Pop.>6M)")
  # 
  # # sevplot<-subsev%>%select(-c("IDP","iso3","Continent","eventid","hday","day","sdate","HA.NAT.EQ","VU","CC.INF","INFORM"))
  # # tmp<-subsev%>%select(paste("Pop",cIntMap_m,sep="_"))%>%apply(1, function(x) log(1+x*exp((IntMap-0.5)/Imin)))%>%t()
  # tmp<-subsev%>%select(paste("Pop",cIntMap_m,sep="_"))%>%apply(1, function(x) log(1+x))%>%t()
  # tmp[tmp==0]<-NA
  # sevplot[,paste("Pop",cIntMap_m,sep="_")]<-tmp
  # 
  # colnames(sevplot)<-c("log(Max Disp.)","exp(Mag.)","GDACS Alert","EQ Phys. Exp.","Vuln","Infra",
  #          "INFORM",paste("log(exp(Mag.')*Pop>",cIntMap_m,"M)",sep=""),"log(Pop.>4.5M)","log(Area>4.5M)","Bool(INFORM)")  
  # 
  # p<-ggpairs(sevplot, lower = list(continuous = "smooth"),mapping = ggplot2::aes(colour=subsev$INF_bool));p #+ theme(legend.position = "bottom");p
  # ggsave("EQ_CorrelationMatrix_mini.png", plot=p,path = paste0(directory,'IIDIPUS_Results'),width = 7,height = 5)
  # 
  # sevplot<-subsev%>%mutate(PopInfra_Exp=(Pop_5.5*CC.INF/HA.NAT.EQ),PopInfraExp=(Pop_5.5*CC.INF*HA.NAT.EQ),
  #                          PopInfraVU=(Pop_5.5*CC.INF*VU),PopVU_Exp=(Pop_5.5*VU/HA.NAT.EQ),PopVU_Exp=(Pop_5.5*VU*HA.NAT.EQ))
  # sevplot<-sevplot[,c(1,6:9,11:14,21:24)]
  # p<-ggpairs(sevplot, lower = list(continuous = "smooth"));p #+ theme(legend.position = "bottom");p
  
  
  
  tmp<-subsev%>%select(paste("Pop",cIntMap_m,sep="_"))%>%apply(1, function(x) log(1+x*exp((IntMap-0.5)/Imin)))%>%t()
  colnames(tmp)<-paste("expI_Pop",cIntMap_m,sep="_")
  # tmp[tmp==0]<-NA
  subsev%<>%cbind(tmp)
  tmp<-subsev%>%select(paste("Pop",cIntMap_m,sep="_"))%>%apply(1, function(x) log(1+x))%>%t()
  # tmp[tmp==0]<-NA
  subsev[,paste("Pop",cIntMap_m,sep="_")]<-tmp
  
  mExp<-median(subsev$HA.NAT.EQ)
  subsev%<>%mutate(modExp=(HA.NAT.EQ-mExp)*(HA.NAT.EQ-mExp)+1)
  # s_train<-filter(subsev,eventid!=6559)
  
  Ns = floor(0.85*nrow(subsev))
  set.seed(123)   # set seed to ensure you always have same random numbers generated
  Nb<-1000
  
  eqn<-c("logmax ~ Pop_5+Pop_5.5+Pop_6",
         "logmax ~ Pop_5.5+Pop_6",
         "logmax ~ Pop_5+Pop_5.5+Pop_6+MaxGDP_6M",
         "logmax ~ Pop_5.5+Pop_6+MaxGDP_6M",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6)*MaxGDP_6M",
         "logmax ~ (Pop_5.5+Pop_6)*MaxGDP_6M",
         "logmax ~ Pop_5+Pop_5.5+Pop_6 + 0",
         "logmax ~ Pop_5.5+Pop_6 + 0",
         "logmax ~ Pop_5+Pop_5.5+Pop_6+MaxGDP_6M + 0",
         "logmax ~ Pop_5.5+Pop_6+MaxGDP_6M + 0",
         "logmax ~ Pop_5+Pop_5.5+Pop_6+MaxGDP_6M + haz + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*VU + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*modExp + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)/modExp + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF/modExp + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*VU/modExp + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF*VU + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF*VU/modExp + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF*VU*modExp + 0",
         "logmax ~ (Pop_5+Pop_5.5+Pop_6+MaxGDP_6M)*INFORM + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*VU + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*modExp + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)/modExp + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF/modExp + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*VU/modExp + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF*VU + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF*VU/modExp + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*CC.INF*VU*modExp + 0",
         "logmax ~ (Pop_5.5+Pop_6+MaxGDP_6M)*INFORM + 0",
         "logmax ~ expI_Pop_5+expI_Pop_5.5+expI_Pop_6",
         "logmax ~ expI_Pop_5.5+expI_Pop_6",
         "logmax ~ expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M",
         "logmax ~ expI_Pop_5.5+expI_Pop_6+MaxGDP_6M",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6)*MaxGDP_6M",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6)*MaxGDP_6M",
         "logmax ~ expI_Pop_5+expI_Pop_5.5+expI_Pop_6 + 0",
         "logmax ~ expI_Pop_5.5+expI_Pop_6 + 0",
         "logmax ~ expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M + 0",
         "logmax ~ expI_Pop_5.5+expI_Pop_6+MaxGDP_6M + 0",
         "logmax ~ expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M + haz + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*VU + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*modExp + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)/modExp + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF/modExp + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*VU/modExp + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF*VU + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF*VU/modExp + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF*VU*modExp + 0",
         "logmax ~ (expI_Pop_5+expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*INFORM + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*VU + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*modExp + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)/modExp + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF/modExp + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*VU/modExp + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF*VU + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF*VU/modExp + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*CC.INF*VU*modExp + 0",
         "logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*INFORM + 0")
  
  prediction<-data.frame()         
  for(equation in eqn){
    predictors<-data.frame()
    for (i in 1:Nb){
      ind = sample(seq_len(nrow(subsev)),size = Ns)
      train<-subsev[ind,]
      test<-subsev[-ind,]
      # predy<-lm(formula = logmax ~ Pop_4.5+Pop_5+Pop_5.5+Pop_6+Pop_6.5, data = s_train)
      predy<-lm(formula = as.formula(equation),data = train)
      # summary(predy)
      # print(paste0(ceiling(AIC(predy)),", ",ceiling(BIC(predy))))
      predme<-tryCatch(predict(predy,test, interval = 'confidence'),error = function(e) rep(NA,15))
      # print(paste0("Prediction for PHL, 12-2019 : ",prediction[1],", with CI : [",prediction[2],"-",prediction[3],
      # "], Actual : ",filter(subsev,eventid==6559)%>%pull(logmax)%>%exp()))
      me<-tryCatch(colMeans(apply(predme,2,function(x) abs(exp(x)-exp(test$logmax))/exp(test$logmax)),na.rm = T),error = function(e) rep(NA,3))
      predictors%<>%rbind(data.frame(ME_E=me[1],ME_L=me[2],ME_U=me[3],AIC=AIC(predy),BIC=BIC(predy)))
    }
    prediction<-rbind(colMeans(predictors),prediction)
    # print(paste0("IIDIPUS Likelihood : ",round(prediction[1],digits = 2),", with CI : [",
    #              round(prediction[2],digits = 2),"-",round(prediction[3],digits = 2),"]"))
    # print(paste0("GDACS Likelihood : ",round(prediction[4],digits = 2),", with CI : [",
    #              round(prediction[5],digits = 2),"-",round(prediction[6],digits = 2),"]"))
  }
  colnames(prediction)<-colnames(predictors)
  prediction<-cbind(eqn,prediction)
  prediction%<>%mutate(CI=ME_L+ME_U)
  View(prediction)
  
  predG<-data.frame()
  for (i in 1:Nb){
    ind = sample(seq_len(nrow(subsev)),size = Ns)
    train<-subsev[ind,]
    test<-subsev[-ind,]
    
    p_GDACS<-lm(formula = logmax ~ alert+0,data = train)
    predGDACS<-predict(p_GDACS,test, interval = 'confidence')
    # print(paste0("GDACS Prediction for PHL, 12-2019 : ",predGDACS[1],", with CI : [",predGDACS[2],"-",predGDACS[3],
    # "], Actual : ",filter(subsev,eventid==6559)%>%pull(logmax)%>%exp()))
    gd<-colMeans(apply(predGDACS,2,function(x) abs(exp(x)-exp(test$logmax))/exp(test$logmax)),na.rm = T)
    predG%<>%rbind(data.frame(GD_E=gd[1],GD_L=gd[2],GD_U=gd[3]))
  }
  predG<-colMeans(predG)
  predG<-c(predG,predG[2]+predG[3])
  names(predG)<-c(names(predG)[1:3],"CI")
  # colnames(predG)<-c("GD_E","GD_L","GD_U")
  print(predG)
  
  
  # logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*modExp + 0
  
  fullmodel<-data.frame()
  for (i in 1:nrow(subsev)){
    train<-subsev[-i,]
    test<-subsev[i,]
    predy<-lm(formula = logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*modExp + 0 ,data = train)
    predme<-tryCatch(ceiling(exp(predict(predy,test, interval = 'confidence'))),error = function(e) rep(NA,15))
    # exp(x)-exp(test$logmax)
    p_GDACS<-lm(formula = logmax ~ alert+0,data = train)
    predGDACS<-ceiling(exp(predict(p_GDACS,test, interval = 'confidence')))
    
    fullmodel%<>%rbind(data.frame(ME_E=predme[1],ME_L=predme[2],ME_U=predme[3],pGDACS_E=predGDACS[1],pGDACS_L=predGDACS[2],pGDACS_U=predGDACS[3],ACTUAL=exp(test$logmax)))
  }
  fullmodel%<>%mutate(diff=abs(pGDACS_E-ACTUAL)/ACTUAL)
  
  finalsev<-cbind(subsev,select(fullmodel,ME_E,pGDACS_E,ACTUAL,diff,ME_L,ME_U,pGDACS_L,pGDACS_U))
  
  finalsev%<>%select(iso3,haz,sdate,ME_E,ACTUAL,diff,ME_L,ME_U,pGDACS_E,pGDACS_L,pGDACS_U)
  
  View(finalsev)
  finalsev%<>%arrange(diff)
  
  colnames(finalsev)<-c("Country (iso3C)","EQ Magnitude [R]","Date","LM Estimate","Actual Disp.","diff",
                        "LM Lower Bound","LM Higher Bound", "GDACS Estimate","GDACS Lower Bound", "GDACS Higher Bound")
  finalsev<-finalsev[,c(1:4,9,5:8,10,11)]
  finalsev[,2]<-4.5*log(finalsev[,2])
  write.csv(x = finalsev,file = paste0(directory,"IIDIPUS_Results/LM_Results.csv"))
  
  formula = logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*modExp + 0
  
  # 
  # predictors<-data.frame()
  # for (i in 1:Nb){
  #   ind = sample(seq_len(nrow(subsev)),size = Ns)
  #   train<-subsev[ind,]
  #   test<-subsev[-ind,]
  #   # predy<-lm(formula = logmax ~ Pop_4.5+Pop_5+Pop_5.5+Pop_6+Pop_6.5, data = s_train)
  #   predy<-lm(formula = logmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*modExp + 0 ,data = train)
  #   # summary(predy)
  #   # print(paste0(ceiling(AIC(predy)),", ",ceiling(BIC(predy))))
  #   predme<-tryCatch(predict(predy,test, interval = 'confidence'),error = function(e) rep(NA,15))
  #   # print(paste0("Prediction for PHL, 12-2019 : ",prediction[1],", with CI : [",prediction[2],"-",prediction[3],
  #   # "], Actual : ",filter(subsev,eventid==6559)%>%pull(logmax)%>%exp()))
  #   me<-colMeans(apply(predme,2,function(x) abs(exp(x)-exp(test$logmax))/exp(test$logmax)),na.rm = T)
  #   
  #   p_GDACS<-lm(formula = logmax ~ alert+0,data = train)
  #   predGDACS<-predict(p_GDACS,test, interval = 'confidence')
  #   gd<-colMeans(apply(predGDACS,2,function(x) abs(exp(x)-exp(test$logmax))/exp(test$logmax)),na.rm = T)
  #   
  #   predictors%<>%rbind(data.frame(ME_E=me[1],ME_L=me[2],ME_U=me[3],GD_E=gd[1],GD_L=gd[2],GD_U=gd[3]))
  # }
  
  tsev<-subsev
  tmp<-tsev%>%select(paste("Pop",cIntMap_m,sep="_"))%>%apply(1, function(x) x*exp((IntMap-0.5)/Imin))%>%t()
  colnames(tmp)<-paste("expI_Pop",cIntMap_m,sep="_")
  tsev[,paste("expI_Pop",cIntMap_m,sep="_")]<-tmp
  
  
  
  
  
  fullmodel<-data.frame()
  for (i in 1:nrow(tsev)){
    train<-tsev[-i,]
    test<-tsev[i,]
    predy<-lm(formula = gmax ~ (expI_Pop_5.5+expI_Pop_6+MaxGDP_6M)*modExp + 0 ,data = train)
    predme<-tryCatch(ceiling(predict(predy,test, interval = 'confidence')),error = function(e) rep(NA,15))
    # exp(x)-exp(test$logmax)
    p_GDACS<-lm(formula = gmax ~ alert,data = train)
    predGDACS<-ceiling(predict(p_GDACS,test, interval = 'confidence'))
    fullmodel%<>%rbind(data.frame(Expt=c(predme[1],predGDACS[1]),Low=c(predme[2],predGDACS[2]),Upp=c(predme[3],predGDACS[3]),Actual=c(test$gmax,test$gmax),Model=c("IIDIPUS","GDACS")))
    # fullmodel%<>%rbind(data.frame(ME_E=predme[1],ME_L=predme[2],ME_U=predme[3],pGDACS_E=predGDACS[1],pGDACS_L=predGDACS[2],pGDACS_U=predGDACS[3],ACTUAL=test$gmax))
  }
  
  fullmodel%<>%mutate(diff=abs(Expt-Actual)/Actual,SqE=(Expt-Actual)*(Expt-Actual)/Actual)
  ival<-ceiling(sqrt(mean(fullmodel$SqE[fullmodel$Model=="IIDIPUS"])))
  gval<-ceiling(sqrt(mean(fullmodel$SqE[fullmodel$Model=="GDACS"])))
  
  fullmodel%<>%mutate(CIdiff=abs(Low + Upp - 2*Actual)/Actual,CISqE=(Low + Upp - 2*Actual)*(Low + Upp - 2*Actual)/Actual)
  CIival<-ceiling(sqrt(mean(fullmodel$CISqE[fullmodel$Model=="IIDIPUS"])))
  CIgval<-ceiling(sqrt(mean(fullmodel$CISqE[fullmodel$Model=="GDACS"])))
  
  fullmodel[fullmodel<0]=0
  
  p<-ggplot(fullmodel,aes(diff,group=Model)) + geom_density(aes(fill=Model),alpha=0.4) + scale_x_log10(limits=c(1e-1,1e4)) + ylim(c(0,0.5)) +
    xlab("Relative Absolute Error") + ylab("Density") +
    annotate("text",label=paste0("Relative RMSE IIDIPUS: ",ival),x=5e2,y=0.35,size=6) +
    annotate("text",label=paste0("Relative RMSE GDACS: ",gval),x=5e2,y=0.32,size=6);p
  ggsave("LM_Model_ndiff_density.png", plot=p,path = paste0(directory,'IIDIPUS_Results/'),width = 8,height = 5)  
  
  p<-ggplot(fullmodel,aes(CIdiff,group=Model)) + geom_density(aes(fill=Model),alpha=0.4) + scale_x_log10(limits=c(1e-1,1e4)) +ylim(c(0,0.5)) +
    xlab("Relative Absolute Error of p.05 + p.95") + ylab("Density") + 
    annotate("text",label=paste0("CI Relative RMSE IIDIPUS: ",CIival),x=5e2,y=0.35,size=6) +
    annotate("text",label=paste0("CI Relative RMSE GDACS: ",CIgval),x=5e2,y=0.32,size=6);p
  ggsave("LM_Model_ndiff_density_CI.png", plot=p,path = paste0(directory,'IIDIPUS_Results/'),width = 8,height = 5)
  
  
  # 
  # p<-ggplot(data.frame(Residuals=c(abs(predy$fitted.values-tsev$gmax)/tsev$gmax,abs(p_GDACS$fitted.values-tsev$gmax)/tsev$gmax),
  #                      Model=c(rep("IIDIPUS",length(exp(predy$residuals))),rep("GDACS",length(exp(p_GDACS$residuals))))),
  #           aes(Residuals,group=Model)) + geom_density(aes(fill=Model),alpha=0.4) + scale_x_log10() +
  #   xlab("Predicted Value / Actual Value") + ylab("Density") + 
  #   annotate("text",label=paste0("IIDIPUS Residuals: ",ceiling(sum(predy$residuals))),x=5e2,y=0.35,size=6) +
  #   annotate("text",label=paste0("GDACS Residuals: ",ceiling(abs(sum(p_GDACS$residuals)))),x=7e2,y=0.32,size=6);p
  # ggsave("LM_Model_Residuals_density.png", plot=p,path = paste0(directory,'IIDIPUS_Results/'),width = 8,height = 5)  
  #  
  # p<-ggplot(data.frame(Residuals=c(abs(predy$fitted.values-tsev$gmax)/tsev$gmax,abs(p_GDACS$fitted.values-tsev$gmax)/tsev$gmax),
  #                   Model=c(rep("IIDIPUS",length(exp(predy$residuals))),rep("GDACS",length(exp(p_GDACS$residuals))))),
  #        aes(Residuals,group=Model)) + geom_density(aes(fill=Model),alpha=0.4) + scale_x_log10() +
  #   xlab("Predicted Value / Actual Value") + ylab("Density") + 
  #   annotate("text",label=paste0("IIDIPUS Residuals: ",ceiling(sum(predy$residuals))),x=5e2,y=0.35,size=6) +
  #   annotate("text",label=paste0("GDACS Residuals: ",ceiling(abs(sum(p_GDACS$residuals)))),x=7e2,y=0.32,size=6);p
  # ggsave("LM_Model_Residuals_density.png", plot=p,path = paste0(directory,'IIDIPUS_Results/'),width = 8,height = 5)  
  
  # finalsev%<>%mutate(diff=abs(`GDACS Estimate`-`Actual Disp.`)/`Actual Disp.`)
  # write.csv(x = finalsev,file = paste0(directory,"IIDIPUS_Results/GDACS_Results.csv"))
  
  # namer<-subsev%>%group_by(sdate)%>%mutate()
  # subsev$sdate<-paste0(AsMonth(subsev$sdate),AsYear(subsev$sdate,T)) #; subsev$sdate<-NULL
  # subsev%<>%group_by(iso3,sdate)%>%mutate(index=1:length(logmax))
  # 
  # library("factoextra")
  # library("FactoMineR")
  # library("corrplot")
  # 
  # subsev<-subsev[,c(2,1,3:19)]
  # subsev%<>%distinct()
  # subsev%<>%as.data.frame()
  # rownames(subsev)<-as.factor(paste0(subsev$iso3,subsev$sdate,"_",subsev$index))
  # subsev%<>%select(-c(index,sdate))
  # # subsev$iso3<-as.factor(subsev$iso3) ; subsev$Continent<-as.factor(subsev$Continent)
  # sev.pca <- PCA(subsev, scale.unit = TRUE, quali.sup = c(1,19), ncp = 6, graph = TRUE)
  # # sev.pca <- PCA(subsev, scale.unit = TRUE, ncp = 6, graph = TRUE)
  # var <- get_pca_var(sev.pca)
  # corrplot(var$cos2, is.corr=FALSE)
  # fviz_pca_var(sev.pca, col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
  # fviz_pca_ind(sev.pca, col.ind = "cos2", 
  #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  #              repel = TRUE)
  # 
  # fviz_pca_ind(sev.pca, col.ind = subsev$Continent,repel = TRUE)
  # 
  # fviz_pca_biplot(sev.pca)
  # 
  #gradient.cols = c("blue", "yellow", "red"),
  #legend.title = "Cont.Var")
  
  # ################### Population Count ###################
  # severity<-subCM%>%group_by(eventid)%>%summarise(haz=max(hazard_severity),gmax=max(gmax),IDP=max(IDP),iso=unique(iso3)[1],alert=unique(hazard_alert)[1],day=min(hday))
  # severity%<>%arrange(desc(eventid))%>%cbind(PopCounts)
  # png(filename = paste0(dir,"IIDIPUS_Results/",hazard,"_SumPopCount-NewDisp_6.0M.png"),width = 650,height = 500)
  # plot(severity[severity[,10]>0,10],severity[severity[,10]>0,3],xlab="No. People Exposed for each Intensity",
  #      ylab="Max No. Displaced People",log="xy",col=colsy[1],pch=19)
  # abline(a = 0,b=1)
  # for (i in 2:(lenIM-1)){
  #   par(new=TRUE)
  #   plot(severity[severity[,(8+i-1)]>0,(8+i-1)],severity[severity[,(8+i-1)]>0,3],log="xy",col=colsy[i],axes=F,pch=19,
  #        xlab="No. People Exposed for each Intensity",ylab="Max No. Displaced People")
  # }
  # legend("bottomright", legend=colnames(PopCounts)[-9],
  #        col=colsy, lty=1, cex=0.8)
  # dev.off()
  # 
  # ~################## Population Density ##################
  # severity<-fCM%>%group_by(eventid)%>%summarise(haz=max(hazard_severity),gmax=max(gmax),IDP=max(IDP),iso=unique(iso3)[1],alert=unique(hazard_alert)[1])
  # severity%<>%arrange(desc(eventid))%>%cbind(MaxPopDens)
  # 
  # png(filename = paste0(dir,"IIDIPUS_Results/",hazard,"_MaxPopDens-IDP.png"),width = 650,height = 500)
  # plot(severity[severity[,9]>0,8],severity[severity[,9]>0,4],xlab="Max Population Density Exposed for each Intensity",ylab="Max No. Displaced People",log="xy",col=colsy[1],pch=19)
  # abline(a = 0,b=1)
  # for (i in 2:(lenIM-1)){
  #   par(new=TRUE)
  #   plot(severity[severity[,(8+i-1)]>0,(8+i-1)],severity[severity[,(8+i-1)]>0,4],log="xy",col=colsy[i],axes=F,pch=19,xlab="Max Population Density Exposed for each Intensity",ylab="Max No. Displaced People")
  # }
  # legend("bottomright", legend=colnames(PopCounts)[-9],
  #        col=colsy, lty=1, cex=0.8)
  # dev.off()
  # 
  # ######################## GDP-PPP #######################
  # severity<-fCM%>%group_by(eventid)%>%summarise(haz=max(hazard_severity),gmax=max(gmax),IDP=max(IDP),iso=unique(iso3)[1],alert=unique(hazard_alert)[1])
  # severity%<>%arrange(desc(eventid))%>%cbind(MaxGDP)
  # 
  # png(filename = paste0(dir,"IIDIPUS_Results/",hazard,"_MaxGDP-IDP.png"),width = 650,height = 500)
  # plot(severity[severity[,9]>0,8],severity[severity[,9]>0,4],xlab="Max GDP-PPP (2011 USD) Exposed for each Intensity",ylab="Max No. Displaced People",log="xy",col=colsy[1],pch=19)
  # abline(a = 0,b=1)
  # for (i in 2:(lenIM-1)){
  #   par(new=TRUE)
  #   plot(severity[severity[,(8+i-1)]>0,(8+i-1)],severity[severity[,(8+i-1)]>0,4],log="xy",col=colsy[i],axes=F,pch=19,xlab="Max GDP-PPP (2011 USD) Exposed for each Intensity",ylab="Max No. Displaced People")
  # }
  # legend("bottomright", legend=colnames(PopCounts)[-9],
  #        col=colsy, lty=1, cex=0.8)
  # dev.off()
  
  ####################### ADD 1D COUNTRY-FACTORS #######################
  
  ##### Add columns for hazard 'shakemap' #####
  # POPULATION PER HAZARD ZONE:
  # tmp<-as.data.frame.array(array(NA,dim=c(length(CM$gmax),length(IntMap))));colnames(tmp)<-cIntMap
  # tmp2<-as.data.frame.array(array(NA,dim=c(length(CM$gmax),2)));colnames(tmp2)<-c("MaxPopDens","MaxGDP")
  # CM%<>%cbind(tmp,tmp2); rm(tmp,tmp2)
  # index<-1:dim(CM)[2]
  # minizone<-index[colnames(CM)==IntMap[1]]
  # maxizone<-index[colnames(CM)==IntMap[length(IntMap)]]
  # rm(index)
  # CM[CM$iso3==iso&CM$eventid==event,minizone:maxizone]<-FUCK YOU
  
}

