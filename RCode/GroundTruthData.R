subnat_file <- '~/Downloads/EQ_SubNational.xlsx'
SubNatData <- openxlsx::read.xlsx(subnat_file, colNames = TRUE , na.strings = c("","NA"))
cleanSubNatData <- function(SubNatData){
  # Clean data from xlsx file by converting data to the appropriate formats/data types
  SubNatData$source_date <- openxlsx::convertToDate(SubNatData$source_date)
  SubNatData$sdate <- openxlsx::convertToDate(SubNatData$sdate)
  SubNatData$fdate <- openxlsx::convertToDate(SubNatData$fdate)
  SubNatData$mortality <- as.integer(SubNatData$mortality)
  SubNatData$displacement <- as.integer(SubNatData$displacement)
  SubNatData$buildDam <- as.integer(SubNatData$buildDam)
  SubNatData$buildDest <- as.integer(SubNatData$buildDest)
  return(SubNatData)
}
SubNatData %<>% cleanSubNatData()
# Identify events by name and sdate (make sure all rows corresponding to the same event have the same name and sdate!)
SubNatDataByEvent <- SubNatData %>% group_by(event_name, sdate) %>% group_split()
mort_list <- list()
disp_list <- list()
buildDam_list <- list()
buildDest_list <- list()
j_mort <- 1
j_disp <- 1
j_buildDam <- 1
j_buildDest <- 1
for (i in 1:length(SubNatDataByEvent)){
  # NatDataByEvent <- SubNatDataByEvent[[i]] %>% filter(is.na(Region) & is.na(Subregion))
  # morts_recorded <- NatDataByEvent %>% group_by(iso3) %>% dplyr::select('mortality') %>% na.omit() %>% group_split()
  # disps_recorded <- NatDataByEvent %>% group_by(iso3) %>% dplyr::select('displacement') %>% na.omit() %>% group_split()
  # buildDam_recorded <- NatDataByEvent %>% group_by(iso3) %>% dplyr::select('buildDam') %>% na.omit() %>% group_split()
  # buildDest_recorded <- NatDataByEvent %>% group_by(iso3) %>% dplyr::select('buildDest') %>% na.omit() %>% group_split()
  PolyDataByEvent <- SubNatDataByEvent[[i]] %>% group_by(iso3, Region, Subregion)
  morts_recorded <- PolyDataByEvent %>% dplyr::select('mortality') %>% dplyr::filter(!is.na(mortality)) %>% group_split()
  disps_recorded <- PolyDataByEvent %>% dplyr::select('displacement') %>% dplyr::filter(!is.na(displacement)) %>% group_split()
  buildDam_recorded <- PolyDataByEvent %>% dplyr::select('buildDam') %>% dplyr::filter(!is.na(buildDam)) %>% group_split()
  buildDest_recorded <- PolyDataByEvent %>% dplyr::select('buildDest') %>% dplyr::filter(!is.na(buildDest)) %>% group_split()
  if (length(morts_recorded) > 0 ){
    for (k in 1:length(morts_recorded)){
      mort_list[[j_mort]] <- morts_recorded[[k]]$mortality
      j_mort <- j_mort + 1
    }
  }
  if (length(disps_recorded) > 0 ){
    for (k in 1:length(disps_recorded)){
      disp_list[[j_disp]] <- disps_recorded[[k]]$displacement
      j_disp <- j_disp + 1
    }
  }
  if (length(buildDam_recorded) > 0 ){
    for (k in 1:length(buildDam_recorded)){
      buildDam_list[[j_buildDam]] <- buildDam_recorded[[k]]$buildDam
      j_buildDam <- j_buildDam + 1
    }
  }
  if (length(buildDest_recorded) > 0 ){
    for (k in 1:length(buildDest_recorded)){
      buildDest_list[[j_buildDest]] <- buildDest_recorded[[k]]$buildDest
      j_buildDest <- j_buildDest + 1
    }
  }
}
k <- 10
fun_log_dif <- function(array,mn=F){
  if(length(array)==1 | sd(array)==0){
    return(NA)
  }
  if(mn) return(mean(array))
  log_difs <- sd(abs(array+1))/mean(array+1)
}

impies<-data.frame(diffSD=sapply(mort_list, fun_log_dif),
                   meddie=sapply(mort_list, fun_log_dif,mn=T),
                   impact="mortality")
impies%<>%rbind(data.frame(diffSD=sapply(disp_list, fun_log_dif),
                   meddie=sapply(disp_list, fun_log_dif,mn=T),
                   impact="displacement"))
impies%<>%rbind(data.frame(diffSD=sapply(buildDam_list, fun_log_dif),
                           meddie=sapply(buildDam_list, fun_log_dif,mn=T),
                           impact="buildDam"))
impies%<>%rbind(data.frame(diffSD=sapply(buildDest_list, fun_log_dif),
                           meddie=sapply(buildDest_list, fun_log_dif,mn=T),
                           impact="buildDest"))

impies%<>%filter(!is.na(diffSD) & !is.na(meddie))

p<-impies%>%ggplot(aes(meddie,diffSD,group=impact))+geom_point(aes(colour=impact))+
  scale_x_log10()+scale_y_log10() + xlab("Mean Impact [log-scale]") + ylab("Standard Deviation / Mean Impact [log-scale]")+
  labs(colour="Impact Type")
ggsave("GroundTruth-Impacts.eps",p,path="./Plots/IIDIPUS_Results/",width=9,height=4,device = grDevices::cairo_ps)  
impies%>%ggplot(aes(diffSD,group=impact))+geom_density(aes(fill=impact),alpha=0.3)+
  scale_x_log10()
#





mort_sd <- sapply(mort_list, fun_log_dif)
mort_sd <- mort_sd[which(!is.na(mort_sd))]
disp_sd <- sapply(disp_list, fun_log_dif)
disp_sd <- disp_sd[which(!is.na(disp_sd))]
buildDam_sd <- sapply(buildDam_list, fun_log_dif)
buildDam_sd <- buildDam_sd[which(!is.na(buildDam_sd))]
buildDest_sd <- sapply(buildDest_list, fun_log_dif)
buildDest_sd <- buildDest_sd[which(!is.na(buildDest_sd))]
mean(mort_sd)/mean(mort_sd)
mean(disp_sd)/mean(mort_sd)
mean(buildDam_sd)/mean(mort_sd)
mean(buildDest_sd)/mean(mort_sd)
(1/mean(mort_sd))/(1/mean(disp_sd))
(1/mean(disp_sd))/(1/mean(disp_sd))
(1/mean(buildDam_sd))/(1/mean(disp_sd))
(1/mean(buildDest_sd))/(1/mean(disp_sd))
plot(x=1,y=1, xlim=c(1,6), ylim=c(0,2))
for (i in 1:length(disp_list)){
  if(length(disp_list[[i]])>1){
    obs <- disp_list[[i]]
    for (j in 1:length(obs)){
      points(log10(median(obs)), obs[j]/median(obs))
    }
  }
}
abline(0,1)















