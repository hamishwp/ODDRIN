
ExtractCheckedDisp<-function(dir,haz="EQ"){
  
  DispData<-read.csv(paste0(dir,"IIDIPUS_Input/DispData_",haz,".csv"),na.strings = "-")
  DispData$sdate%<>%as.Date()
  DispData$fdate%<>%as.Date()
  return(DispData)
  
}

GetDisplacements<-function(haz, saved=T, reduce=T, GIDD=T, EMDAT=F, dir="./"){
  
  if(saved) return(ExtractCheckedDisp(dir))
  
  CM<-GetHelix(haz=haz,reduce=reduce)
  
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
    openxlsx::read.xlsx(paste0(dir,"/Displacement_Data/emdat.xls"))
    stop("Not ready to filter EMDAT data yet")
  }
  
  return(CM)
  
}
