
folder<-c('Wildfire','East','Gippsland')

GetNetworkCoverage<-function(FBdirectory,folder){
  
  name<-paste(folder,collapse = '.*')
  dirs <- paste0(FBdirectory,list.files(path=FBdirectory,pattern=name))
  temp <- paste0(dirs,'/Network_Coverage/',list.files(path=paste0(dirs,"/Network_Coverage"),pattern='.*csv'))
  print(temp)
  NetCov<-read.csv(temp,header = TRUE,fill = TRUE,sep = ","); NetCov[is.na(NetCov)]<-0L

  return(NetCov)
  
}
