# library("ssh")
library(DBI)
library(tidyverse)
source(paste0(directory,"RCode/Functions.R"))
# detach(package:plyr)
library(dplyr)
library(latex2exp)
library(magrittr)

# GRID year for calculation:
#@@@ NOTE : Should be GRID year -1 !!!
myear<-2016L

# Helix Names : 
# unique(helix$hazard_type)
# [1] "Flood"               "Storm"               NA                   
# [4] "Mass Movement"       "Wildfire"            "Earthquake"         
# [7] "Extreme temperature" "Volcanic eruption"   "Drought"            
# [10] "Mass movement"  

# year number_events_recorded
# 2001      1
# 2004      1
# 2005      1
# 2008     12
# 2009     13
# 2010     11
# 2011     15
# 2012     14
# 2013     17
# 2014     19
# 2015     25
# 2016    612
# 2017    920
# 2018   1852
# 2019   2141
# 2020    326

### Frequency of sources used by IDMC FOR DISASTERS ###
# tmp<-helix%>%group_by(sources)%>%summarise(freq=length(figure))%>%arrange(desc(freq))
# sources                                                                                                    freq
# 1 "[\"Ministry of Home Affairs Disaster Management Division National Emergency Response Centre (NERC)\"]"  3102
# 2 "[\"Local Authorities\"]"                                                                                2287
# 3  NA                                                                                                      2228
# 4 "[\"National Disaster Relief Services Centre\"]"                                                         1945
# 5 "[\"International Organization for Migration (IOM)\"]"                                                   1390
# 6 "[\"National Agency for Disaster Management (BNPB)\"]"                                                   1110
# 7 "[\"Disaster Response Operations Monitoring and Information Center (DROMIC)\"]"                          1038
# 8 "[\"Federal Emergency Management Agency (FEMA)\"]"                                                        800
# 9 "[\"International Federation of the Red Cross (IFRC)\"]"                                                  528
# 10 "[\"Office for the Coordination of Humanitarian Affairs (OCHA)\"]"                                        517
# 11 "[\"IOM Displacement Tracking Matrix (IOM DTM)\"]"                                                        499
# 12 "[\"Nepal Disaster Risk Reduction Portal \"]"                                                             399
# 13 "[\"Red Cross/Red Crescent Movement\"]"                                                                   375
# 14 "[\"Government\"]"                                                                                        324
# 15 "[\"State Disaster Management Authority\"]"                                                               318
# 16 "[\"Cabinet Office\"]"                                                                                    280
# 17 "[\"National Disaster Management Authority (NDMA)\"]"                                                     264
# 18 "[\"Disaster Management Centre\"]"                                                                        251
# 19 "[\"Local Media\"]"                                                                                       221
# 20 "[\"Ministry of Civil Affairs\"]"                                                                         212


convHelixHaz<-function(d_type,helix=T){
  
  if(helix){
    d_choice<-c(
      "DR"="Drought",
      "TC"="Tropical Cyclones",
      "VW"="Violent Wind",
      "FL"="Flood",
      "VO"="Volcanic eruption",
      "EQ"="Earthquake"
  )} else {
  d_choice<-c(
    "Drought"="DR",
    "Tropical Cyclones"="TC",
    "Violent Wind"="VW",
    "Flood"="FL",
    "Volcanic eruption"="VO",
    "Earthquake"="EQ"
  )}
    
  namez<-d_choice[d_type]
  
  if (any(is.na(namez))){
    print(paste0("WARNING: Helix input disaster type ",d_type[is.na(namez)]," does not exist"))
  }  
  
  return(unname(namez))
  
}

GetGIDD<-function(dir,haz=NULL){
  library(openxlsx)
  file<-paste0(dir,"Helix/idmc_disaster_pre2017_dataset_hp.xlsx")
  # GIDD<-read.xlsx(file,sheetName = "Sheet1",header = TRUE,startRow = 2)
  GIDD<-read.xlsx(file,startRow = 2,sheet = 4)
  GIDD<-GIDD[,1:8]
  colnames(GIDD)<-c("iso3","country","year","start_date","event_name","hazard_category","hazard_type","figure")
  GIDD$hazard_type[GIDD$hazard_type %in% c("Dry mass movement","Wet mass movement","Wet Mass movement","Wet Mass Movement")]<-"Mass movement"
  GIDD$hazard_type[GIDD$hazard_type=="Volcanic activity"]<-"Volcanic eruption"
  GIDD%<>%dplyr::select(-hazard_category)
  GIDD%<>%filter(year<2017 & year>2012 & !is.na(figure))
    
  GIDD$hazard_type<-convHelixHaz(GIDD$hazard_type,helix = F)
  names(GIDD)[which(names(GIDD)=="hazard_type")]<-"hazard"
  if(!is.null(haz)) GIDD%<>%filter(hazard==haz)
  
  names(GIDD$hazard)<-NULL
  
  return(GIDD)
}

MergeGIDD_Helix<-function(GIDD,CM){
  len<-length(GIDD$iso3)
  evs<--1:-len
  qual<-rep("total",len)
  tCM<-GIDD%>%transmute(iso3=iso3,gmax=figure,hazard=hazard,sdate=start_date,eventid=evs,qualifier=qual)
  return(rbind(CM,tCM))
}

CleanHelix<-function(helix,reduce=F,haz="EQ"){
  
  # Remove conflict displacement
  helix<-helix[helix$displacement_type=="Disaster",]
  
  # Remove all Master facts: stock figures that represent an entire country or region for multiple disasters
  helix %<>% filter(!is.na(event_id) | !is.na(event_name) & !(class%in%c("Master","Sub")))
  # helix %<>% group_by(event_id)%>%filter(min(as.numeric(format(date,"%Y")),na.rm = TRUE)>=myear)
  
  #######################################################################
  # EXCEPTIONS : EUGHHHHHHHH
  helix<-helix %>% filter(!(id %in% c(34821,35943))) # Albania event_id=6472 remove housing damage so IDP stock is used.
  helix<-helix %>% filter(!(event_id %in% c(6611,6641,5331,5709))) 
  #######################################################################
  
  # CLEANING - TERM
  helix$hazard_type[helix$hazard_type %in% c("Dry mass movement","Wet mass movement")]<-"Mass movement"
  helix$term[helix$term=="relocated"]<-"Relocated"
  helix$term[helix$term=="displaced"]<-"Displaced"
  helix$term[helix$term=="evacuated"]<-"Evacuated"
  helix$term[helix$term=="evacuation"]<-"Evacuated"
  helix$term[helix$term=="evacuated, rescued"]<-"Evacuated"
  # helix$term[helix$term=="Failed Return / Returnee Displacement (Flow)"]<-
  helix$term[helix$term=="In Relief Camp"]<-"Sheltered"
  listy<-c("Affected","Forced to Flee","Failed Return / Returnee Displacement (Flow)")
  helix %<>% filter(!(term %in% listy & type!="IDPs (Stock)"))
  helix %<>% filter(term!="Homeless")
  
  # CLEANING - TYPE
  listy<-c("Partial (Flow)","Unverified (Flow)","Failed Local Integration (Flow)",
           "Unverified (Stock)","Cross-border Flight (Flow)",
           "Local Integration (Flow)","Cross-border Return (Flow)",
           "Relocation Elsewhere (Flow)","Provisional Solutions (Flow)",
           "Multiple Displacement (Flow)","Return (Flow)","Deaths (Flow)",
           "People Displaced Across Borders (Stock)","Partial (Stock)",
           "Failed Return / Returnee Displacement (Flow)","Returnees (Stock)")
  helix %<>% filter(!(type %in% listy))
  helix %<>% filter(!(term=="Relocated" & type %in% c("New Displacement (Flow)")))
  helix$type[helix$type=="Locally Integrated IDPs (Stock)"]<-"Returnees (Stock)"
  helix %<>% filter(event_role=="Recommended figure")
  
  names(helix)[which(names(helix)=="hazard_type")]<-"hazard"
  names(helix)[which(names(helix)=="event_id")]<-"eventid"
  
  if(reduce){
    
    print("WARNING: removing evacuation-related displacements from Helix (see CleanHelix)")
    helix%<>%filter(term!="Evacuated")%>%
      dplyr::select(c(event_name,iso3,type,figure,term,start_date,date,
                     eventid,hazard,event_start_date,qualifier))
    tmp<-helix%>%group_by(eventid)%>%
      summarise(gmax=max(figure,na.rm = T),
                sdate=min(c(event_start_date,start_date,date),na.rm = T),
                .groups = 'drop_last')
    
    helix%<>%merge(tmp,by="eventid")%>%
      dplyr::select(iso3,gmax,hazard,sdate,eventid,qualifier)
    
    if(haz=="EQ"){
      # ADD EXTRA DISPLACEMENT EVENTS FOR BUILDING DAMAGE
      isos<-c("PAK","NPL") #,"AFG","PAK","ECU")
      sdate<-c("2013-09-24","2015-04-25") #,"2015-10-26","2015-10-26","2016-04-17")
      lenny<-length(isos)
      eventid<-(max(helix$eventid)+1):(max(helix$eventid)+lenny)
      hazs<-rep("Earthquake",lenny)
      qual<-rep("approximately",lenny)
      
      gmax<-c(78000,2622733) #,54642,665812,259000) #,445343,54866,133631)
      # USGS<-c("b000jyiv","us20002926")
      source<-c("https://www.humanitarianresponse.info/sites/www.humanitarianresponse.info/files/documents/files/Balochistan%20Earthquake%202013%20Report%2023%20Oct_Final.pdf",
                "https://reliefweb.int/sites/reliefweb.int/files/resources/ETR_2020_web-1.pdf")
                # "GIDD","GIDD","GIDD")
      
      helix%<>%rbind(data.frame(iso3=isos,gmax=gmax,hazard=hazs,sdate=sdate,eventid=eventid,qualifier=qual))
    }
    
  }
  
  helix$hazard<-convHelixHaz(helix$hazard,helix = F)
  
  return(helix)
  
}

GetHelix<-function(clean=T,haz=NULL,reduce=F){
  # Connect to IDMC Helix database
  con <- dbConnect(RPostgres::Postgres(), dbname = "backend_production", host="localhost", port=5432, user="infographics", password="dummy")
  # dbListTables(con) 
  helix<-dbReadTable(conn = con,name = "facts_view")
  
  if(!is.null(haz)) helix%<>%filter(hazard_type==convHelixHaz(haz,helix = T))
  if(clean) helix%<>%CleanHelix(reduce=reduce,haz=haz)
  
  return(helix)
}





















GetCorrelationMatrix<-function(helix,gyear=2016,directory="/home/patten/Documents/Coding/IIDIPUS/"){
  
  helix%<>%CleanHelix()
  # What to do with "Partial (Stock)"?
  
  # Initialisations
  CM<-data.frame()
  
  ####################### LOOP OVER EVENTS #######################
  for (event in unique(helix$eventid)){
    # Extract event details from Helix
    Esub<-helix%>%filter(eventid==event & event_data_included)  
    hazard<-Esub$hazard_type
    for (iso in unique(Esub$iso3)){
      ECsub<-Esub%>%filter(iso3==iso)
      # Initialisations
      pHD<-fHD<-Settler<-Evac<-tEvac<-mxShelt<-mnShelt<-tShelt<-mx<-dmx<-RoR<-wRoR<-tweight<-nweight<-mpred<-NA    
      # If evacuations occured, this date could be different from IDP stock
      gmax<-max(ECsub$figure)
      # Sometimes max values are repeated, making it seem like there are more data points
      sdate<-max(ECsub$date[ECsub$figure==gmax])
      
      # Partial Housing Damage
      tmp<-ECsub%>%filter(type=="New Displacement (Flow)" & term=="Partially Destroyed Housing")%>%pull(figure)
      if(length(tmp)>0) pHD<-max(tmp,na.rm = T)
      # Fully Destroyed Housing Damage
      tmp<-ECsub%>%filter(type=="New Displacement (Flow)" & term=="Destroyed Housing")%>%pull(figure)
      if(length(tmp)>0) fHD<-max(tmp,na.rm = T)
      # IDPs settled elsewhere
      tmp<-ECsub%>%filter(type=="IDPs Settled Elsewhere (Stock)")%>%pull(figure)
      if(length(tmp)>0) Settler<-max(tmp,na.rm = T)
      # Evacuated people (short term)
      EECsub<-ECsub%>%filter(term=="Evacuated")
      tmp<-EECsub%>%pull(figure)
      if(length(tmp)==1) {
        Evac<-max(tmp,na.rm = T)
      } else if (length(tmp)>1){
        Evac<-max(tmp,na.rm = T)
        tEvac<-as.numeric(max(EECsub$date,na.rm = T)-min(EECsub$date,na.rm = T))
        if(is.na(tEvac)|tEvac==0) tEvac<-NA
      }
      # People in shelters (long & short term)
      SECsub<-ECsub%>%filter(term=="Sheltered") #%>%group_by(date)%>%summarise(stock=max(figure,na.rm = T))
      tmp<-SECsub%>%pull(figure)
      if(length(tmp)==1) {
        mxShelt<-max(tmp,na.rm = T)
      } else if (length(tmp)>1){
        tmp<-SECsub%>%group_by(date)%>%summarise(figure=max(figure,na.rm = T))
        mxShelt<-max(tmp$figure)
        mnShelt<-min(tmp$figure)
        tShelt<-as.numeric(max(tmp$date)-min(tmp$date))
        if(tShelt==0) tShelt<-mnShelt<-NA
      }
      
      ##### STOCK FIGURE EVALUATION #####
      Isum<-ECsub%>%filter(type=="IDPs (Stock)")
      if (length(Isum$figure)<=1) {
        if (length(Isum$figure)==1) {
          # Find peak value of event
          mx<-max(Isum$figure)
          # What is the start date of the stock?
          dmx<-max(Isum$date[Isum$figure==mx])
        }
        CM<-rbind(CM,data.frame(pHD=pHD,fHD=fHD,Settler=Settler,Evac=Evac,tEvac=tEvac,
                                mxShelt=mxShelt,mnShelt=mnShelt,tShelt=tShelt,mx=mx,
                                dmx=dmx,RoR=RoR,wRoR=wRoR,tweight=tweight,nweight=nweight,
                                mpred=mpred,iso3=iso,sdate=sdate,eventid=event,gmax=gmax,hazard=hazard))
        next
      }
      
      # Find peak value of event
      mx<-max(Isum$figure)
      dmx<-min(Isum$date[Isum$figure==mx])
      Isum$day<-as.numeric(Isum$date-dmx)
      Isum%<>%group_by(day)%>%summarise(figure=max(figure,na.rm = T))
      # Isum%<>%group_by(day,qualifier)%>%summarise(figure=max(figure))
      # IDP stock start day
      Isum%<>%filter(day>=14)
      
      if(length(unique(Isum$day))<=1) {
        CM<-rbind(CM,data.frame(pHD=pHD,fHD=fHD,Settler=Settler,Evac=Evac,tEvac=tEvac,
                                mxShelt=mxShelt,mnShelt=mnShelt,tShelt=tShelt,mx=mx,
                                dmx=dmx,RoR=RoR,wRoR=wRoR,tweight=tweight,nweight=nweight,
                                mpred=mpred,iso3=iso,sdate=sdate,eventid=event,gmax=gmax,hazard=hazard))
        next
      }
      
      # Take all values>14days after dmx
      # check how many unique stock figures exist
      # log-linear fit, extract p-value (no normalisation used)
      # For log, can't have log(0) so set to 1
      Isum$figure[Isum$figure<1]<-1
      tlist<-GetWeightedRoR(Isum$day,Isum$figure)
      # tlist<-GetWeightedRoR(Isum$day,Isum$figure,Isum$qualifier)
      
      # temporary allocation of decay rate for testing
      tRR<-tlist[1]
      # Using decay rate, predict how many people were displaced at start date, if more than mx set all to NA
      if (tRR<=0 & tlist[2]<mx) {
        RoR<-tlist[1]
        # 2 weightings: fit significance and days between first and last fitted points
        wRoR<-tlist[3]
        tweight<-as.numeric(max(Isum$day))
        nweight<-length(unique(Isum$figure))
        mpred<-tlist[2]
      }
      
      CM<-rbind(CM,data.frame(pHD=pHD,fHD=fHD,Settler=Settler,Evac=Evac,tEvac=tEvac,
                              mxShelt=mxShelt,mnShelt=mnShelt,tShelt=tShelt,mx=mx,
                              dmx=dmx,RoR=RoR,wRoR=wRoR,tweight=tweight,nweight=nweight,
                              mpred=mpred,iso3=iso,sdate=sdate,eventid=event,gmax=gmax,hazard=hazard))
    }
    
  }
  
  # Normalise the weighting for tweight and nweight?
  # CM$tweight=CM$tweight/max(CM$tweight,na.rm = T)
  # CM$nweight=CM$nweight/max(CM$nweight,na.rm = T)
  
  save(CM,paste0(directory,"Helix/CorrelationMatrix.Rdata"))
  
  return(CM)
  
}

###### Plots of correlation matrix #####
PlotCorrelationMatrix<-function(CM,WBiso=NULL,directory=directory<-"/home/patten/Documents/Coding/IIDIPUS/",haz=NULL){
  
  for (i in 1:length(CM)){
    print(paste0("Number of non-NA's in ",colnames(CM)[i]," = ",length(CM[!is.na(CM[,i]),i])))
  }
  
  length(unique(CM[,11]))
  
  # if(!is.null(hazard)) CM%<>%filter(hazard==haz)
  
  p<-ggplot(CM,aes(x=RoR,y=wRoR))+
    geom_point()+ggtitle("Significance of Fit") +ylab("P-value") +xlab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$"));p
  ggsave('RoR_pvalue.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 5,height = 5)
  
  p<-ggplot(CM,aes(x=RoR,y=nweight))+
    geom_point()+ggtitle("Significance of Fit") +ylab("Number of points for fit") +xlab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$"));p
  ggsave('RoR_nweight.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 5,height = 5)
  
  p<-ggplot(CM,aes(x=RoR,y=tweight))+geom_point()+ggtitle("Significance of Fit") +ylab("Time span for fit") +
    xlab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$")) + scale_y_log10();p
  ggsave('RoR_tweight.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 5,height = 5)
  
  p<-ggplot(data=CM,aes(x=as.factor(hazard),y=RoR))+
    geom_violin(draw_quantiles = ,na.rm = T,scale="width",fill="red")+xlab("Hazard") + ggtitle("IDP Rate of Return") +ylab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$"));p
  ggsave('RoR_violin.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 9,height = 5)
  
  rCM<-CM%>%filter(!(hazard%in%c("Wildfire","Mass movement","Drought","Volcanic eruption","Extreme temperature",NA)))
  
  p<-ggplot(rCM,aes(x=mpred,y=RoR,colour=tweight))+
    geom_point()+ggtitle("IDP Rate of Return") +xlab(TeX("IDP Stock Intercept $N_0$")) +ylab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$")) +
    scale_x_log10() + scale_color_gradient(name = "Max(Days after Event)", trans = "log10")
  p<-p+facet_wrap( ~ hazard, scales = "free")+ theme(plot.title = element_text(hjust = 0.5));p
  ggsave('MpredvsRoR_tweight.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 8,height = 5)
  
  p<-ggplot(rCM,aes(x=mx,y=RoR,colour=tweight))+
    geom_point()+ggtitle("IDP Rate of Return") +xlab("Maximum IDP stock") +ylab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$")) +
    scale_x_log10() + scale_color_gradient(name = "Max(Days after Event)", trans = "log10")
  p<-p+facet_wrap( ~ hazard, scales = "free")+ theme(plot.title = element_text(hjust = 0.5));p
  ggsave('IDPvsRoR.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 8,height = 5)
  
  p<-ggplot(rCM,aes(x=gmax,y=RoR,colour=tweight))+
    geom_point()+ggtitle("IDP Rate of Return") +xlab("Maximum Displaced") +ylab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$")) +
    scale_x_log10() + scale_color_gradient(name = "Max(Days after Event)", trans = "log10")
  p<-p+facet_wrap( ~ hazard, scales = "free")+ theme(plot.title = element_text(hjust = 0.5));p
  ggsave('gMAXvsRoR.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 8,height = 5)
  
  p<-ggplot(rCM,aes(x=fHD,y=RoR,colour=tweight))+
    geom_point()+ggtitle("IDP Rate of Return") +xlab("Number of people displaced due to house destruction") +ylab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$")) +
    scale_x_log10() + scale_color_gradient(name = "Max(Days after Event)", trans = "log10")
  p<-p+facet_wrap( ~ hazard, scales = "free")+ theme(plot.title = element_text(hjust = 0.5));p
  ggsave('fHDvsRoR.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 8,height = 5)
  
  p<-ggplot(rCM,aes(x=pHD,y=RoR,colour=tweight))+
    geom_point()+ggtitle("IDP Rate of Return") +xlab("Number of people displaced due to partial house destruction") +ylab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$")) +
    scale_x_log10() + scale_color_gradient(name = "Max(Days after Event)", trans = "log10")
  p<-p+facet_wrap( ~ hazard, scales = "free")+ theme(plot.title = element_text(hjust = 0.5));p
  ggsave('pHDvsRoR.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 8,height = 5)
  
  if(is.null(WBiso)) tWBiso<-GetWBinfo(directory)
  
  FTS<-Viso2WB(as.character(CM$iso3),WBiso,ind="FTS",Score="wmScore")
  dFTS<-Viso2WB(as.character(CM$iso3),WBiso,ind="FTS",Score="dScore")
  GINI<-Viso2WB(as.character(CM$iso3),WBiso,ind="GINI",Score="wmScore")
  U5M<-Viso2WB(as.character(CM$iso3),WBiso,ind="U5M",Score="wmScore")
  ELEC<-Viso2WB(as.character(CM$iso3),WBiso,ind="ElecAcc",Score="wmScore")
  ODA<-Viso2WB(as.character(CM$iso3),WBiso,ind="ODA",Score="wmScore")
  dODA<-Viso2WB(as.character(CM$iso3),WBiso,ind="ODA",Score="dScore")
  HDI<-Viso2WB(as.character(CM$iso3),WBiso,ind="HDI",Score="wmScore")
  dHDI<-Viso2WB(as.character(CM$iso3),WBiso,ind="HDI",Score="dScore")
  HDIGDP<-Viso2WB(as.character(CM$iso3),WBiso,ind="HDIGDP",Score="wmScore")
  dHDIGDP<-Viso2WB(as.character(CM$iso3),WBiso,ind="HDIGDP",Score="dScore")
  
  # Evacuated/(fHD+pHD) vs socioeconomic
  p<-ggplot(CM,aes(x=Evac/mx,y=HDI))+
    geom_point(aes(colour=dHDI))+
    xlab("Number of people Evacuated / IDP stock after") + ylab("HDI - Human Development Index")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Variation in HDI since 2010")
  p
  ggsave('EvacIDP_HDI.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=Evac/mx,y=HDIGDP))+
    geom_point(aes(colour=dHDIGDP))+
    xlab("Number of people Evacuated / IDP stock after") + ylab("HDI-GDP - Human Development Index (GDP based)")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Variation in HDIGDP since 2010")
  p
  ggsave('EvacIDP_HDIGDP.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=Evac/mx,y=FTS))+
    geom_point(aes(colour=dFTS))+
    xlab("Number of people Evacuated / IDP stock after") + ylab("Humanitarian Aid Contributions [dollars]")+
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Variation in HAC since 2010 [dollars]")
  p
  ggsave('EvacIDP_FTS.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=Evac/mx,y=GINI))+
    geom_point(aes(colour=U5M))+
    xlab("Number of people Evacuated / IDP stock after") + ylab("GINI Index")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Under 5 Mortality")
  p
  ggsave('EvacIDP_GINI_U5M.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=Evac/mx,y=U5M))+
    geom_point(aes(colour=GINI))+
    xlab("Number of people Evacuated / IDP stock after") + ylab("Under 5 Mortality")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "GINI Index")
  p
  ggsave('EvacIDP_U5M_GINI.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=Evac/mx,y=ODA))+
    geom_point(aes(colour=dODA))+
    xlab("Number of people Evacuated / IDP stock after") + ylab("ODA - Official Development Assistance")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Variation in ODA since 2010")
  p
  ggsave('EvacIDP_ODA.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=Evac/mx,y=U5M))+
    geom_point(aes(colour=GINI))+
    xlab("Number of people Evacuated / IDP stock after") + ylab("Under 5 Mortality")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "GINI Index")
  p
  ggsave('EvacIDP_U5M_GINI.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  # max(IDP)/Sheltered vs socioeconomic
  p<-ggplot(CM,aes(x=mx/mxShelt,y=HDI))+
    geom_point(aes(colour=dHDI))+
    xlab("max (IDP stock) / max(Sheltered)") + ylab("HDI - Human Development Index")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Variation in HDI since 2010")
  p
  ggsave('IDPShelt_HDI.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=mx/mxShelt,y=HDIGDP))+
    geom_point(aes(colour=dHDIGDP))+
    xlab("max (IDP stock) / max(Sheltered)") + ylab("HDI-GDP - Human Development Index (GDP based)")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Variation in HDIGDP since 2010")
  p
  ggsave('IDPShelt_HDIGDP.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=mx/mxShelt,y=FTS))+
    geom_point(aes(colour=tShelt))+
    xlab("max (IDP stock) / max(Sheltered)") + ylab("Humanitarian Aid Contributions [dollars]")+
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Shelter Data Time Coverage",trans="log")
  p
  ggsave('IDPShelt_FTS.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=mx/mxShelt,y=GINI))+
    geom_point(aes(colour=tShelt))+
    xlab("max (IDP stock) / max(Sheltered)") + ylab("GINI Index")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Shelter Data Time Coverage",trans="log")
  p
  ggsave('IDPShelt_GINI.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=mx/mxShelt,y=ELEC))+
    geom_point(aes(colour=tShelt))+
    xlab("max (IDP stock) / max(Sheltered)") + ylab("Access to Electricity")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Shelter Data Time Coverage",trans="log")
  p
  ggsave('IDPShelt_ELEC.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=mx/mxShelt,y=ODA))+
    geom_point(aes(colour=tShelt))+
    xlab("max (IDP stock) / max(Sheltered)") + ylab("ODA - Official Development Assistance")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Shelter Data Time Coverage",trans="log")
  p
  ggsave('IDPShelt_ODA.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  p<-ggplot(CM,aes(x=mx/mxShelt,y=U5M))+
    geom_point(aes(colour=tShelt))+
    xlab("max (IDP stock) / max(Sheltered)") + ylab("Under 5 Mortality")+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_gradient(name = "Shelter Data Time Coverage",trans="log")
  p
  ggsave('IDPShelt_U5M.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 7,height = 5)
  
  # Do same for Shelter and tShelter?
  
  # +geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2)
  # +ggtitle("IDP Rate of Return") +xlab("Number of people displaced due to partial house destruction") +ylab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$")) +
  #   scale_x_log10() + scale_color_gradient(name = "Max(Days after Event)", trans = "log10")
  #
  
  
  
  p<-ggplot(rCM,aes(x=pHD,y=RoR,colour=tweight))+
    geom_point()+ggtitle("IDP Rate of Return") +xlab("Number of people displaced due to partial house destruction") +ylab(TeX("$\\alpha$ from $N=N_0e^{\\alpha t}$")) +
    scale_x_log10() + scale_color_gradient(name = "Max(Days after Event)", trans = "log10")
  
  # p<-ggplot(,aes())+geom_point(aes())+ggtitle()+ylab()+xlab();
  # p<-p+facet_wrap( ~ DS, scales = "free_y",labeller = as_labeller(TeX,default = label_parsed), nrow = 2)+ theme(plot.title = element_text(hjust = 0.5)) + ggtitle();p
  # ggsave('.png', plot=p,path = paste0(directory,'Plots/Helix'),width = 6,height = 5)
  
  
}

NormaliseHelix<-function(helix){
  
  normy<-helix%>%group_by(eventid)%>%mutate(maxy=max(figure)) %>% mutate(figure=figure/maxy)
  return(normy)
  
}

GetIDETECT<-function(){
  # Connect to IDMC Helix database
  con <- dbConnect(RPostgres::Postgres(), dbname = "idetect", 
                   host="localhost", port=5433, user="postgres")
  dbListTables(con)
  #   PHL<-dbGetQuery(con = con,
  #                   "SELECT F.*, A.* FROM idetect_facts AS F LEFT JOIN 
  # idetect_analysis_facts AS
  #                    AF ON AF.fact = F.id LEFT JOIN idetect_analyses AS A 
  # ON A.gkg_id = AF.analysis
  #                    WHERE UPPER(iso3) = 'PHL' AND EXTRACT(YEAR FROM 
  # A.publication_date) IN (2019, 2020)")
  ALL<-dbGetQuery(con = con,
                  "SELECT F.*, A.* FROM idetect_facts AS F LEFT JOIN 
idetect_analysis_facts AS
                   AF ON AF.fact = F.id LEFT JOIN idetect_analyses AS A 
ON A.gkg_id = AF.analysis
                   WHERE EXTRACT(YEAR FROM A.publication_date) IN (2019, 
2020)")
  gkg<-dbReadTable(conn = con,name = "gkg")
  
  return(helix)
}


# CM%>%group_by(hazard)%>%summarise(freq=length(RoR))%>%arrange(desc(freq))
# # A tibble: 9 x 2
# hazard                 freq
# 1 Flood               13008
# 2 Storm               11605
# 3 Mass movement        1692
# 4 Wildfire             1592
# 5 Earthquake           1077
# 6 Drought               401
# 7 Volcanic eruption     315
# 8 Extreme temperature   232
# 9 NA                     25

# "Number of non-NA's in pHD = 8166"
# "Number of non-NA's in fHD = 15128"
# "Number of non-NA's in Settler = 48"
# "Number of non-NA's in Evac = 19381"
# "Number of non-NA's in tEvac = 12857"
# "Number of non-NA's in mxShelt = 11426"
# "Number of non-NA's in mnShelt = 8262"
# "Number of non-NA's in tShelt = 8262"
# "Number of non-NA's in mx = 20398"
# "Number of non-NA's in dmx = 20391"
# "Number of non-NA's in RoR = 3448"
# "Number of non-NA's in wRoR = 3174"
# "Number of non-NA's in tweight = 3448"
# "Number of non-NA's in nweight = 3448"
# "Number of non-NA's in mpred = 3448"
# "Number of non-NA's in iso3 = 29947"
# "Number of non-NA's in sdate = 29821"
# "Number of non-NA's in eventid = 29947"
# "Number of non-NA's in gmax = 29838"
# "Number of non-NA's in hazard = 29922"
