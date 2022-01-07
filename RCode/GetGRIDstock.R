library(ssh)
library(DBI)
library(tidyverse)
source("RCode/Functions.R")
library(GGally)
library(GauPro)
detach(package:plyr)
library(dplyr)

# GRID year for calculation:
#@@@ NOTE : Should be GRID year -1 !!!
gyear<-2019L
# Do we calculate using events for year > gyear or year = gyear?
moregyear<-TRUE

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

GetGRIDyear<-function(helix,gyear){
  
  # Remove conflict displacement
  helix<-helix[helix$displacement_type=="Disaster",]
  
  #######################################################################
  # EXCEPTIONS : EUGHHHHHHHH
  helix<-helix %>% filter(!(id %in% c(34821,35943))) # Albania event_id=6472 remove housing damage so IDP stock is used.
  helix<-helix %>% filter(!(event_id %in% c(6611,6641,5331,5709))) 
  
  #######################################################################
  
  # Group mass movement
  ind<-helix$hazard_type %in% c("Dry mass movement","Wet mass movement")
  helix$hazard_type[ind]<-"Mass movement"
  
  # Initialisations
  hazard<-hazy<-data.frame()
  
  for (cntry in unique(helix$country)) {
    # Initialisations
    masters<-0
    hazard<-data.frame()
    
    # Split helix database by country
    Csubs<-filter(helix,country==cntry )
    
    print(paste0("Number of events in ",cntry," = ",length(unique(Csubs$event_id))))
    
    #@@@@@ SCENARIO 1 : OCCURING BEFORE GYEAR YET ONGOING @@@@@#
    # Find events that started before gyear and are ongoing through and after gyear
    scen1<-Csubs%>%group_by(event_id)%>%filter(min(as.numeric(format(date,"%Y")),na.rm = TRUE)<gyear & 
                                                 max(as.numeric(format(date,"%Y")),na.rm = TRUE)>=gyear)
    # Extract stock estimate for last value produced in gyear 
    if(length(scen1$event_id)>0){
      if(!moregyear){scen1<-scen1%>%filter(date<paste0(gyear,"-12-31"))}
      tmp<-scen1 %>% group_by(event_id) %>% 
        filter(as.POSIXct(date)==as.POSIXct(max(date)) & grepl("IDPs (",type,fixed = TRUE) & !str_detect(class,"Master"))%>%
        arrange(desc(figure))%>%slice(1L)
      hazard<-rbind(hazard,data.frame(stock=tmp$figure, scenario=rep("sc1",length(tmp$figure)), hazard=tmp$hazard_type, event=tmp$event_id))
      Csubs<-Csubs%>%filter(!(event_id %in% tmp$event_id & country==cntry))
      rm(tmp,scen1)
    }
    # Remove all events created before gyear 
    if(moregyear){Csubs<-Csubs%>%filter(year>=gyear)} else {Csubs<- Csubs%>%filter(year==gyear)}
    # Remove all events that started after gyear
    Csubs<-Csubs%>%group_by(event_id)%>%filter(!min(as.numeric(format(date,"%Y")),na.rm = TRUE)>gyear)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
    
    #@@@@@ SCENARIO 6 : MASTER FACTS @@@@@#
    # Find number of NA event_id's and event_names
    #mastsub<-Csubs%>%filter(is.na(event_id) & is.na(event_name))
    mastsub<-Csubs%>%filter(str_detect(class,"Master") & grepl("IDPs (",type,fixed = TRUE) )# & term %in% c("Displaced","displaced"))
    # Do any of these have master facts?
    if(length(mastsub$figure)>0){
      masters<-mastsub%>%filter(as.POSIXct(date)==as.POSIXct(max(date)))%>%summarise(IDP=sum(figure))%>%pull(IDP)
      if(sum(masters,na.rm = TRUE)>0){print(cntry)}
      if(length(mastsub$figure)>0){Csubs<-Csubs%>%filter(!(event_id %in% mastsub$event_id & country==cntry))}
    }
    rm(mastsub)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
    
    Csubs<-Csubs%>%filter(!is.na(event_id) | !is.na(event_name) & !str_detect(class,"Master"))
    
    
    
    ####################### LOOP OVER EVENTS #######################
    for (event in unique(Csubs$event_id)){
      ECsub<-Csubs%>%filter(event_id==event & event_data_included)
      #if(length(unique(ECsub$event_id))==0){next}
      # for (loc in unique(ECsub$centroid)){
      lECsub<-ECsub#%>%filter(centroid==loc)
      IlECsub<-lECsub%>%filter(grepl("IDPs (",type,fixed = TRUE) )#& term %in% c("Displaced","displaced"))
      HlECsub<-lECsub%>%filter(term=="Destroyed Housing")
      
      #@@@@@ SCENARIO 7 : IGNORE FOR NO OR ONE IDP STOCK VALUE & NO HOUSING DAMAGE INFO @@@@@#
      if(length(IlECsub$figure)<=1 & length(HlECsub$figure)==0){
        
        if(length(IlECsub$groups)>0 & length(IlECsub$event_role)>0){
          if(length(IlECsub$figure)==1 & grepl(IlECsub$groups,pattern = paste0("GRID ",gyear+1)) & str_detect(IlECsub$event_role,"Recommended figure")){
            hazard<-rbind(hazard,data.frame(stock=IlECsub$figure, scenario="sc7*", hazard=IlECsub$hazard_type[1],event=IlECsub$event_id[1]))
              next
            }
          }
        
          hazard<-rbind(hazard,data.frame(stock=NA, scenario="sc7", hazard=IlECsub$hazard_type[1],event=NA))
          next
          
        }
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        
        #@@@@@ SCENARIO 2 & 4 : HAS HOUSING DAMAGE INFO BUT NO IDP STOCK @@@@@#
        if(length(IlECsub$figure)<=1){
          sc2_4<-HlECsub%>%filter(grepl("New Displacement",type,fixed = TRUE))%>%
                                    arrange(desc(figure))%>%slice(1L)%>%pull(figure)
                                    # figure==max(figure,na.rm = TRUE))%>%pull(figure)
          if(sum(sc2_4,na.rm = TRUE)>0){hazard<-rbind(hazard,data.frame(stock=sc2_4, scenario="sc2_4", hazard=HlECsub$hazard_type[1],event=HlECsub$event_id[1]))}
          next
        }
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        
        #@@@@@ SCENARIO 5 : MULTIPLE IDP STOCK VALUES BUT AT START OF YEAR - USE HOUSING DAMAGE INFO @@@@@#
        if(length(HlECsub$figure)>0 & (max(IlECsub$date,na.rm = T) - max(HlECsub$date,na.rm = T) <60) & !(cntry=="Philippines")){ 
          sc5<-HlECsub%>%filter(grepl("New Displacement",type,fixed = TRUE))%>%
            arrange(desc(figure))%>%slice(1L)%>%pull(figure)
          # sc5<-HlECsub%>%filter(grepl("New Displacement",type,fixed = TRUE) & 
          #         as.POSIXct(date)==as.POSIXct(max(date)))%>%arrange(desc(figure))%>%slice(1L)%>% pull(figure)
          if(sum(sc5,na.rm = TRUE)>0){hazard<-rbind(hazard,data.frame(stock=sc5, scenario="sc5", hazard=HlECsub$hazard_type[1],event=HlECsub$event_id[1]))}
          next
        }
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        
        #@@@@@ SCENARIO 3 : MULTIPLE IDP STOCK VALUES @@@@@#
        sc3<-IlECsub%>%filter(as.POSIXct(date)==as.POSIXct(max(date)) & !term=="Homeless")%>%
          summarise(stock=sum(figure,na.rm = TRUE))%>%arrange(desc(stock))%>%pull(stock)%>%sum(na.rm = TRUE)
        # sc3<-IlECsub%>%filter(as.POSIXct(date)==as.POSIXct(max(date)) & grepl(IlECsub$groups,pattern = paste0("GRID ",gyear+1)))%>%
        #   summarise(stock=sum(figure,na.rm = TRUE))%>%arrange(desc(stock))%>%pull(stock)%>%sum(na.rm = TRUE)
        # sc3<-IlECsub%>%filter(as.POSIXct(date)==as.POSIXct(max(date)) & !str_detect(IlECsub$event_role,"Triangulation"))%>%
        #   summarise(stock=sum(figure,na.rm = TRUE))%>%arrange(desc(stock))%>%pull(stock)%>%sum(na.rm = TRUE)
        if(sum(sc3,na.rm = TRUE)>0){hazard<-rbind(hazard,data.frame(stock=sc3, scenario="sc3", hazard=IlECsub$hazard_type[1],event=IlECsub$event_id[1]))}
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
      # }
    }
    ############################################################
    
    
    
    stock<-sum(hazard$stock,na.rm = TRUE)
    
    # Finish scenario 6!
    if(which.max(c(stock,masters))==1){
      hazy<-rbind(hazy,data.frame(country=rep(paste0(Csubs$iso3[!is.na(Csubs$iso3)][1],", ",cntry),length(hazard$stock)),IDP=hazard$stock,hazard=hazard$hazard,scenario=hazard$scenario,
                                  event=hazard$event)) # paste0(Csubs$iso3,", ",cntry)
      print(paste0("Total IDP stock in ",cntry," = ",sum(hazard$stock,na.rm = TRUE)))
    } else {
      hazy<-rbind(hazy,data.frame(country=paste0(Csubs$iso3[!is.na(Csubs$iso3)][1],", ",cntry),IDP=masters,hazard=NA,scenario="sc6",event=NA))
      print(" ")
      print(paste0("SCENARIO 6: Total IDP stock in ",cntry," sc6 = ",masters))
      print(" ")
    }
    
  }
  
  scen<-hazy%>%group_by(country)%>%summarise(IDP=sum(IDP,na.rm = T))
  write_csv(arrange(scen,desc(IDP)),path = paste0(directory,"GRID/IDPstock_GRID",gyear+1,"_country_scenarios.csv"))
  
  scen6<-hazy%>%filter(str_detect(scenario,"sc6"))%>%group_by(country)%>%summarise(IDP=sum(IDP,na.rm = T))
  # scen6$ISO<-sapply(strsplit(as.character(hazy$country), split=',', fixed=TRUE), `[`, 1)
  write_csv(arrange(scen6,desc(IDP)),path = paste0(directory,"GRID/IDPstock_GRID",gyear+1,"_country_scenario-6.csv"))
  
  hazy<-hazy%>%group_by(country,hazard,event)%>%summarise(IDP=sum(IDP,na.rm = TRUE))
  cIDP<-hazy%>%group_by(country,event)%>%summarise(IDP=sum(IDP,na.rm = TRUE))
  write_csv(arrange(cIDP,desc(IDP)),path = paste0(directory,"GRID/IDPstock_GRID",gyear+1,"_country_event.csv"))
  
  hazy<-hazy%>%group_by(country,hazard)%>%summarise(IDP=sum(IDP,na.rm = TRUE))
  
  write_csv(arrange(hazy,desc(IDP)),path = paste0(directory,"GRID/IDPstock_GRID",gyear+1,"_country_hazard.csv"))
  hIDP<-hazy%>%group_by(hazard)%>%summarise(IDP=sum(IDP,na.rm = TRUE))
  write_csv(arrange(hIDP,desc(IDP)),path = paste0(directory,"GRID/IDPstock_GRID",gyear+1,"_hazard.csv"))
  cIDP<-hazy%>%group_by(country)%>%summarise(IDP=sum(IDP,na.rm = TRUE))
  write_csv(arrange(cIDP,desc(IDP)),path = paste0(directory,"GRID/IDPstock_GRID",gyear+1,"_country.csv"))
  
  IDP<-sum(hazy$IDP,na.rm = TRUE)
  
  print(paste0("Total IDP stock for ",gyear+1," = ",IDP))
  
  return(hazy)
  
}

NormaliseHelix<-function(helix){
  
  normy<-helix%>%group_by(event_id)%>%mutate(maxy=max(figure)) %>% mutate(figure=figure/maxy)
  return(normy)
  
}

GetHelix<-function(){
  # Connect to IDMC Helix database
  con <- dbConnect(RPostgres::Postgres(), dbname = "backend_production", host="localhost", port=5432, user="infographics", password="dummy")
  dbListTables(con) 
  helix<-dbReadTable(conn = con,name = "facts_view")
  
  return(helix)
}

# Get data from Helix
helix<-GetHelix()
# Clean that database!
helix<-GetGRIDyear(helix,gyear)  

