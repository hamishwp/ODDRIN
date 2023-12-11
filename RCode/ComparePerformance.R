
# GDACS:
for (i in )

# PAGER:
Hazy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/HAZARDobjects_additionalInfo3/EQ20230206TUR_169')
Hazy



ufiles <- ODD

addPAGER <- function(folder_in='IIDIPUS_Input'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  haz_files_dir <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/HAZARDobjects_PAGERfull/'
  ufiles_haz <- list.files(path=haz_files_dir,pattern=Model$haz,recursive = T,ignore.case = T)
  haz_ids <- as.numeric(sapply(ufiles_haz, function(x) as.numeric(strsplit(x, "_")[[1]][2])))
  for (file in ufiles){
    event_id <- as.numeric(strsplit(file, "_")[[1]][2])
    ODDy <- readRDS(paste0(ODD_folderin, file))
    HAZy <- readRDS(paste0(haz_files_dir, ufiles_haz[which(haz_ids==event_id)]))
    ODDy_withPAGER <- ODDy
    ODDy_withPAGER@hazinfo$alertfull <- HAZy$hazard_info$alertfull
    if(is.null(ODDy_withPAGER)){
      print(event_id)
      next
    }
    saveRDS(ODDy_withPAGER, paste0(ODD_folderout, file))
  }
}

makeNational <- function(folder_in='IIDIPUS_Input'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  for (file in ufiles){
    event_id <- as.numeric(strsplit(file, "_")[[1]][2])
    ODDy <- readRDS(paste0(ODD_folderin, file))
    impact_summ <- data.frame(impact=character(), observed=integer(), iso3=character(), polygon=integer())
    for (impact_type in c('mortality', 'displacement')){
      impact_filt <- ODDy@impact[which(ODDy@impact$impact==impact_type),]
      if (NROW(impact_filt)==0) next
      if (length(Reduce(intersect,sapply(ODDy@polygons[impact_filt$polygon], function(x) x$indexes)))> min(0.5 * sapply(ODDy@polygons[impact_filt$polygon], function(x) length(x$indexes)))){
        stop()
      } else {
        impact_summ %<>% add_row(impact=impact_type, observed=sum(impact_filt$observed), iso3='TOT', polygon=1)
      }
    }
    ODDy@impact <- impact_summ
    ODDy@hazinfo$alertlevel_fatalities <- sapply(ODDy@hazinfo$alertfull, function(x) x$alert_level)
    alertmax_fatalities <- ifelse('red' %in% ODDy@hazinfo$alertlevel_fatalities, 'red', ifelse('orange' %in% ODDy@hazinfo$alertlevel_fatalities, 'orange', ifelse('yellow' %in% ODDy@hazinfo$alertlevel_fatalities, 'yellow', ifelse('green' %in% ODDy@hazinfo$alertlevel_fatalities, 'green', 'null'))))
    alertnull_fatalities <- 'null' %in% ODDy@hazinfo$alertlevel_fatalities
    ODDy@impact$alertlevel_fatalities = alertmax_fatalities
    ODDy@impact$alertnull_fatalities = alertnull_fatalities
    ODDy@polygons <- list(list(name='Total', indexes=1:NROW(ODDy@data)))
    if(is.null(ODDy@impact)){
      print(event_id)
      next
    }
    saveRDS(ODDy, paste0(ODD_folderout, file))
  }
}
