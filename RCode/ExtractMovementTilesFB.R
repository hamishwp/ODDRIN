filez<-list.files("./IIDIPUS_Input/Haiti_MovementTiles_2021-08-15_2021-08-28_csv/")

out<-data.frame()
for(fff in filez){
  
  movement<-read_csv(paste0("./IIDIPUS_Input/Haiti_MovementTiles_2021-08-15_2021-08-28_csv/",fff))
  out<-movement%>%group_by(start_polygon_name)%>%summarise(displaced=sum(n_difference,na.rm = T))%>%rbind(out)
  
}

sumzy<-out%>%group_by(start_polygon_name)%>%summarise(sder=sd(displaced),displaced=sum(displaced))

sumzy$lowbound<-sumzy$displaced-3*sumzy$sder
sumzy$highbound<-sumzy$displaced+3*sumzy$sder
View(sumzy)
