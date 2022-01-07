# Extract Environment Variables
source('RCode/GetEnv.R')
# Set the working directory from your environment variables
setwd(directory)

Dispy<-read.csv(paste0(directory,"/Helix/EQs20210808_DispData.csv"))
# Convert the retrieve column into a boolean
Dispy$Retrieve<-as.logical(as.integer(Dispy$Retrieve)-1L)
Dispy$iso3<-as.character(Dispy$iso3)
fails<-c()
for (ev in unique(Dispy$event_id)){
# for (ev in unique(Dispy$event_id)[54:length(unique(Dispy$event_id))]){
  sDisp<-Dispy%>%filter(event_id==ev)%>%arrange(desc(gmax))
  # Skip if the ODD object already exists
  if(!sDisp$Retrieve[1]) next
  # Extract the name of the event
  namer<-paste0("EQ",str_remove_all(sDisp$sdate[1],pattern="-"),sDisp$iso3[1],"_",ev)
  print(namer)
  # Extract Earthquake
  EQ<-tryCatch(GetUSGS(USGSid=as.character(sDisp$USGSid[1])),error=function(e) NULL)
  if(is.null(EQ)) {
    print(paste0("HAZARD fail ",namer))
    fails%<>%c(namer)
    next
  }
  saveRDS(EQ,paste0("./IIDIPUS_Input/HAZARDobjects/",namer))
  # Form an ODD object from the hazard(s)
  ODDy<-tryCatch(new("ODD",lhazSDF=EQ, dir=dir, Model=Model),error=function(e) NULL)
  if(is.null(ODDy)) {
    print(paste0("ODD fail ",namer))
    fails%<>%c(namer)
    next
  }
  # Add observed displacement estimate (provide a default qualifier for now)
  sDisp$qualifier<-"total"
  ODDy@gmax<-sDisp%>%dplyr::select(gmax,qualifier,iso3)
  ODDy@eventid<-ev
  ODDy@cIndies<-WID%>%filter(year==AsYear(ODDy@hazdates[1]) & 
                               iso3%in%unique(ODDy@data$ISO3C))%>%
                dplyr::select(-year)
  
  saveRDS(ODDy,paste0("./IIDIPUS_Input/ODDobjects/",namer))
}

fff<-as.integer(unlist(strsplit(fails,"_"))[seq.int(from=2,to=2L*length(fails),by = 2)])
Dispy%<>%filter(!(event_id%in%fff))
saveRDS(Dispy,"./IIDIPUS_Input/DispData_EQ_V2.Rdata")

WID_perc=   c("p10p100", # top 90% share of Income Distribution
              "p20p100", # top 80% share of Income Distribution
              "p30p100", # top 70% share of Income Distribution
              "p40p100", # top 60% share of Income Distribution
              "p50p100", # top 50% share of Income Distribution
              "p60p100", # top 40% share of Income Distribution
              "p70p100", # top 30% share of Income Distribution
              "p80p100", # top 20% share of Income Distribution
              "p90p100" # top 10% share of Income Distribution
)

WID<-GetWID_perc(WID_perc,Dispy$iso3,AsYear(Dispy$sdate))

# Add on the stragglers - Puerto Rico and Soloman Islands
# Distributions from:
# PRI - https://datausa.io/profile/geo/puerto-rico/
# SLB - https://data.worldbank.org/country/solomon-islands
PRI<-read.csv("./Demography_Data/SocioEconomic/PRI_IncomeDistribution.csv")
out<-interp1(PRI$percentile,PRI$income,xi=seq.int(from=0.1,to=0.9,by=0.1),method="spline")
PRI<-data.frame(variable=unique(WID$variable),
                value=out)
PRI$iso3<-"PRI"
PRI$year<-2020
WID%<>%rbind(PRI)

SLB<-data.frame(quintile=c(0,0.20,0.40,0.60,0.80,1.),income=c(0,0.07,0.114,0.155,0.215,0.446))
SLB$cum<-cumsum(SLB$income)
out<-interp1(SLB$quintile,SLB$cum,xi=seq.int(from=0.10,to=0.90,by=0.10),method="spline")
SLB<-data.frame(variable=unique(WID$variable),
                value=out)
SLB$iso3<-"SLB"
SLB$year<-2016
WID%<>%rbind(SLB)

# Go through each country and add the WID data
path<-"./IIDIPUS_Input/ODDobjects/"
tpath<-"./IIDIPUS_Input/tODDobjects/"
files <- list.files(path=path)
for (fev in files){
  sDisp<-Dispy%>%filter(event_id==as.integer(strsplit(fev,"_")[[1]][2]))%>%
    arrange(desc(gmax))
  ODDy<-readRDS(paste0(path,fev))
  ODDy@eventid<-sDisp$event_id[1]
  namer<-paste0("EQ",str_remove_all(sDisp$sdate[1],pattern="-"),sDisp$iso3[1],"_",sDisp$event_id[1])
  saveRDS(ODDy,paste0(tpath,namer))
}

files <- list.files(path=path)
tfiles <- list.files(path=tpath)

ids<-extractnumbers(files)[seq.int(from=2,to=2L*length(files),by = 2)]
# Sort both by the order of event ids
ix<-sort(ids,index.return = T)$ix
files<-files[ix]
files<-as.Date(as.character((extractnumbers(files)[seq.int(from=1,to=2L*length(files)-1,by = 2)])),format="%Y%m%d")

ids<-extractnumbers(tfiles)[seq.int(from=2,to=2L*length(tfiles),by = 2)]
# Sort both by the order of event ids
ix<-sort(ids,index.return = T)$ix
tfiles<-tfiles[ix]
tfiles<-as.Date(as.character((extractnumbers(tfiles)[seq.int(from=1,to=2L*length(tfiles)-1,by = 2)])),format="%Y%m%d")
(files-tfiles)/12

for(namer in list.files(path=path)){
  ODDy<-readRDS(paste0(path,namer))
  print(extractnumbers(namer)[2])
  ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
}

for(namer in list.files(path=path)){
  ODDy<-readRDS(paste0(path,namer))
  isos<-Dispy%>%filter(event_id==extractnumbers(namer)[2])%>%arrange(desc(gmax))%>%
    pull(iso3)
  if(!any(isos%in%unique(ODDy$ISO3C))) print(namer)
}

for(namer in list.files(path=path)){
  ODDy<-readRDS(paste0(path,namer))
  if(nrow(ODDy@cIndies)<1) {
    print(namer)
    dater<-as.Date(as.character(extractnumbers(namer)[1]),format="%Y%m%d")
    ODDy@cIndies<-WID%>%filter(year==AsYear(dater) & 
                               iso3%in%unique(ODDy@data$ISO3C))%>%
    dplyr::select(-year)
    saveRDS(ODDy,paste0(path,namer))
  }
}
