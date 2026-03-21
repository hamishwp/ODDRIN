library(feather)

# --------------------------------------------------------------------------------
#----------- Prepare data frame of neural network inputs/outputs -----------------
# --------------------------------------------------------------------------------

#Get total impact for each event (e.g. total mortality):
retrieveTotalImpact <- function(ODD){
  if (ODD@hazdates[1] == "2011-02-22"){
    ODD@impact$polygon = 1 #assume reflects total impact
  }
  ODD_df <- as.data.frame(ODD, na.rm=F)
  df_SampledTot <- list()
  impact_types <- unique(ODD@impact$impact)
  
  for (impact_type in impact_types){
    polygon_names <- unlist(lapply(ODD@polygons[ODD@impact$polygon], function(x) x$name))
    if (any(tolower(polygon_names[which(ODD@impact$impact==impact_type)]) %in% c('tot', 'total'))){
      nonmatch <- which(!tolower(polygon_names[which(ODD@impact$impact==impact_type)]) %in% c('tot', 'total'))
      if (length(nonmatch)>0){
        ODD@impact <- ODD@impact[-which(ODD@impact$impact==impact_type)[nonmatch],] # in the case of total and subnational data, remove the subnational
      }
    }
  }
  
  #Many ODD objects don't contain 'total' impact values (to avoid double counting), so we need to obtain these. 
  #We combine the polygons with observations so long as the proportion of overlapping pixels is less than 10%
  #and the total coverage of the exposed area is greater than 90%. Then sum the observations across these polygons.
  observed_total=rep(NA, length(impact_types))
  exposed_haz <- which(apply(ODD_df[,grep('hazMean', names(ODD_df)), drop=F], 1, function(row) any(!is.na(row))) & !is.na(ODD_df$ISO3C) & (ODD_df$Population>0))
  get_overlap_coverage <- function(impact_type){
    indexes_list <- lapply(ODD@polygons[ODD@impact$polygon[which(ODD@impact$impact==impact_type)]], function(x) x$indexes)
    universal_set <- Reduce(union, indexes_list)
    overlap <- Reduce(intersect, indexes_list)
    overlap <- intersect(overlap, exposed_haz)
    prop_overlap <- ifelse(length(indexes_list)>1,length(overlap)/length(universal_set),0)
    prop_coverage <- length(intersect(unique(universal_set), exposed_haz))/length(exposed_haz)#sum(!is.na(ODD_df$ISO3C))
    return(c(prop_overlap, prop_coverage))
  }
  
  overlap_coverage <- sapply( impact_types,get_overlap_coverage)
  
  for (i in 1:length(impact_types)){
    observed_total[i] = sum(ODD@impact$observed[which(ODD@impact$impact==impact_types[i])])
  }
  
  
  df_SampledTot <- data.frame(polygon=0,
                                   iso3 = paste0(unique(ODD@impact$iso3), collapse=' '),
                                   impact=impact_types, 
                                   observed=observed_total,
                                   qualifier='total',
                                   overlap = overlap_coverage[1,],
                                   coverage = overlap_coverage[2,],
                                   inferred=F)
  
  df_SampledTot %<>% filter(overlap < 0.1 & coverage > 0.9)
  return(df_SampledTot)
}

#input_folder <- "IIDIPUS_Input_Alternatives/Nov24Agg/"
input_folder <- "IIDIPUS_Input_Alternatives/Sept25/"
collect_data <- function(input_folder, dat='all'){
  folderin<-paste0(dir,input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  ufiles = rev(ufiles)
  
  impact_all = data.frame(event_id = integer(), iso3=character(), sdate=as.Date(character()), impact=character(), observed=numeric(),
                          Pop4=numeric(), Pop5=numeric(), Pop6=numeric(), Pop7=numeric(), Pop8=numeric(), Pop9=numeric(),
                          Build4=numeric(), Build5=numeric(), Build6=numeric(), Build7=numeric(), Build8=numeric(), Build9=numeric())
  for (filer in ufiles){
    ODD <- readODD(paste0(folderin, filer))
    
    impact_total <- retrieveTotalImpact(ODD)
    if (NROW(impact_total %>% filter(impact=='mortality'))==0){ #incomplete coverage, assess manually whether full impact is included in subnational:
      if( as.numeric(sub(".*_([0-9]+).*", "\\1", filer))==89){
        impact_total = ODD@impact %>% group_by(impact) %>% summarise( polygon=0, observed = sum(observed), iso3='IDN', qualifier='total', overlap=0, coverage=1, inferred=F)
      } else {
        stop()
      }
    }  else {
      print(filer)
    } 
    impact_total$event_id = as.numeric(sub(".*_([0-9]+).*", "\\1", filer))
    #store the number of people exposed to intensity 7 or larger in each pixel
    
    impact_total[paste0('Pop', 4:9)] <- NA
    impact_total[paste0('Build', 4:9)] <- NA
    if (length(grep('hazMean', names(ODD))) > 1){
      hazMax <- apply(values(ODD[[grep('hazMean', names(ODD))]]),1 ,max, na.rm=T)
    } else {
      hazMax <- values(ODD[[grep('hazMean', names(ODD))]])
    }

    #if (ODD@impact$observed[i]==2152) stop()
    impact_total$Pop4 <- sum(ODD$Population[which(hazMax>4 & hazMax < 5)], na.rm=T)
    impact_total$Pop5 <- sum(ODD$Population[which(hazMax>5 & hazMax < 6)], na.rm=T)
    impact_total$Pop6 <- sum(ODD$Population[which(hazMax>6 & hazMax < 7)], na.rm=T)
    impact_total$Pop7 <- sum(ODD$Population[which(hazMax>7 & hazMax < 8)], na.rm=T)
    impact_total$Pop8 <- sum(ODD$Population[which(hazMax>8 & hazMax < 9)], na.rm=T)
    impact_total$Pop9 <- sum(ODD$Population[which(hazMax>9)], na.rm=T)
    
    if('nBuildings' %in% names(ODD)){
      impact_total$Build4 <- sum(ODD$nBuildings[which(hazMax>4 & hazMax < 5)], na.rm=T)
      impact_total$Build5 <- sum(ODD$nBuildings[which(hazMax>5 & hazMax < 6)], na.rm=T)
      impact_total$Build6 <- sum(ODD$nBuildings[which(hazMax>6 & hazMax < 7)], na.rm=T)
      impact_total$Build7 <- sum(ODD$nBuildings[which(hazMax>7 & hazMax < 8)], na.rm=T)
      impact_total$Build8 <- sum(ODD$nBuildings[which(hazMax>8 & hazMax < 9)], na.rm=T)
      impact_total$Build9 <- sum(ODD$nBuildings[which(hazMax>9)], na.rm=T)
    }

    impact_total$event_id = as.numeric(sub(".*_([0-9]+).*", "\\1", filer))
    impact_total$polygon = NULL
    impact_total[,c('polygon', 'qualifier', 'overlap', 'coverage', 'inferred')] = NULL
    impact_total$sdate = min(ODD@hazdates)
    impact_all %<>% rbind(impact_total)
  }
  
  write_feather(impact_all, 'IIDIPUS_Input_Alternatives/Sept25/event_summaries')
  
}

input_folder <- "IIDIPUS_Input_Alternatives/Sept25/"
collect_data_covar <- function(input_folder, dat='all'){
  folderin<-paste0(dir,input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  ufiles = rev(ufiles)
  
  impact_all = data.frame(event_id = integer(), iso3=character(), sdate=as.Date(character()), impact=character(), observed=numeric(),
                          Pop4=numeric(), Pop5=numeric(), Pop6=numeric(), Pop7=numeric(), Pop8=numeric(), Pop9=numeric(),
                          Build4=numeric(), Build5=numeric(), Build6=numeric(), Build7=numeric(), Build8=numeric(), Build9=numeric())
  for (filer in ufiles){
    ODD <- readODD(paste0(folderin, filer))
    
    impact_total <- retrieveTotalImpact(ODD)
    if (NROW(impact_total %>% filter(impact=='mortality'))==0){ #incomplete coverage, assess manually whether full impact is included in subnational:
      if( as.numeric(sub(".*_([0-9]+).*", "\\1", filer))==89){
        impact_total = ODD@impact %>% group_by(impact) %>% summarise( polygon=0, observed = sum(observed), iso3='IDN', qualifier='total', overlap=0, coverage=1, inferred=F)
      } else {
        stop()
      }
    }  else {
      print(filer)
    } 
    impact_total$event_id = as.numeric(sub(".*_([0-9]+).*", "\\1", filer))
    #store the number of people exposed to intensity 7 or larger in each pixel
    
    impact_total[paste0('Pop', 4:9)] <- NA
    impact_total[paste0('Build', 4:9)] <- NA
    if (length(grep('hazMean', names(ODD))) > 1){
      hazMax <- apply(values(ODD[[grep('hazMean', names(ODD))]]),1 ,max, na.rm=T)
    } else {
      hazMax <- values(ODD[[grep('hazMean', names(ODD))]])
    }
    
    #if (ODD@impact$observed[i]==2152) stop()
    impact_total$Pop4 <- sum(ODD$Population[which(hazMax>4 & hazMax < 5)], na.rm=T)
    impact_total$Pop5 <- sum(ODD$Population[which(hazMax>5 & hazMax < 6)], na.rm=T)
    impact_total$Pop6 <- sum(ODD$Population[which(hazMax>6 & hazMax < 7)], na.rm=T)
    impact_total$Pop7 <- sum(ODD$Population[which(hazMax>7 & hazMax < 8)], na.rm=T)
    impact_total$Pop8 <- sum(ODD$Population[which(hazMax>8 & hazMax < 9)], na.rm=T)
    impact_total$Pop9 <- sum(ODD$Population[which(hazMax>9)], na.rm=T)
    
    if('nBuildings' %in% names(ODD)){
      impact_total$Build4 <- sum(ODD$nBuildings[which(hazMax>4 & hazMax < 5)], na.rm=T)
      impact_total$Build5 <- sum(ODD$nBuildings[which(hazMax>5 & hazMax < 6)], na.rm=T)
      impact_total$Build6 <- sum(ODD$nBuildings[which(hazMax>6 & hazMax < 7)], na.rm=T)
      impact_total$Build7 <- sum(ODD$nBuildings[which(hazMax>7 & hazMax < 8)], na.rm=T)
      impact_total$Build8 <- sum(ODD$nBuildings[which(hazMax>8 & hazMax < 9)], na.rm=T)
      impact_total$Build9 <- sum(ODD$nBuildings[which(hazMax>9)], na.rm=T)
    }
    
    impact_total$SHDI = (mean(pull(ODD$SHDI[which(hazMax>4)]), na.rm=T) - Model$center$SHDI$mean)/Model$center$SHDI$sd
    impact_total$EQFreq = (log(mean(pull(ODD$EQFreq[which(hazMax>4)]), na.rm=T)+0.001) - Model$center$EQFreq$mean)/Model$center$EQFreq$sd
    impact_total$Vs30 = (log(mean(pull(ODD$Vs30[which(hazMax>4)]), na.rm=T))  - Model$center$Vs30$mean)/Model$center$Vs30$sd
      
    impact_total$event_id = as.numeric(sub(".*_([0-9]+).*", "\\1", filer))
    impact_total$polygon = NULL
    impact_total[,c('polygon', 'qualifier', 'overlap', 'coverage', 'inferred')] = NULL
    impact_total$sdate = min(ODD@hazdates)
    impact_all %<>% rbind(impact_total)
  }
  
  write_feather(impact_all, 'IIDIPUS_Input_Alternatives/Sept25/event_summaries_covar')
  
}

input_folder <- "IIDIPUS_Input_Alternatives/Sept25/"
collect_data_0.5Intervals <- function(input_folder, dat='all'){
  folderin <- paste0(dir, input_folder, "ODDobjects/")
  ufiles <- na.omit(list.files(path = folderin, pattern = Model$haz, recursive = TRUE, ignore.case = TRUE))
  if (tolower(dat) == 'train') {
    ufiles <- grep('^Train/', ufiles, value = TRUE)
  } else if (tolower(dat) == 'test') {
    ufiles <- grep('^Test/', ufiles, value = TRUE)
  }
  ufiles <- rev(ufiles)
  
  # helper: nice lower-bound labels (4, 4.5, ..., 8.5) + final 9plus
  fmt <- function(x) sub("\\.0$", "", format(round(x, 1), nsmall = 1))
  edges  <- c(seq(4, 9, by = 0.5), Inf)
  lowers <- edges[-length(edges)]
  uppers <- edges[-1]
  pop_cols   <- c(paste0("Pop", fmt(lowers[is.finite(uppers)])), "Pop9plus")
  build_cols <- c(paste0("Build", fmt(lowers[is.finite(uppers)])), "Build9plus")
  
  impact_all <- data.frame(
    event_id = integer(), iso3 = character(), sdate = as.Date(character()),
    impact = character(), observed = numeric(),
    # columns created later
    stringsAsFactors = FALSE
  )
  
  for (filer in ufiles){
    ODD <- readODD(paste0(folderin, filer))
    
    impact_total <- retrieveTotalImpact(ODD)
    if (NROW(dplyr::filter(impact_total, impact == 'mortality')) == 0) {  # incomplete coverage
      if (as.numeric(sub(".*_([0-9]+).*", "\\1", filer)) == 89){
        impact_total <- ODD@impact %>%
          dplyr::group_by(impact) %>%
          dplyr::summarise(polygon = 0, observed = sum(observed), iso3 = 'IDN',
                           qualifier = 'total', overlap = 0, coverage = 1, inferred = FALSE,
                           .groups = "drop")
      } else {
        stop("Missing total mortality impact; encountered non-special case.")
      }
    } else {
      message(filer)
    }
    
    impact_total$event_id <- as.numeric(sub(".*_([0-9]+).*", "\\1", filer))
    
    # Hazard maxima
    haz_idx <- grep('hazMean', names(ODD))
    if (length(haz_idx) > 1){
      hazMax <- apply(values(ODD[[haz_idx]]), 1, max, na.rm = TRUE)
    } else {
      hazMax <- values(ODD[[haz_idx]])
    }
    
    impact_total[paste0('Pop', c(seq(4, 8.5, 0.5), '9plus'))] <- NA
    impact_total[paste0('Build', c(seq(4, 8.5, 0.5), '9plus'))] <- NA
    
    # Compute population in 0.5 bins
    pop_vals <- vapply(seq_along(lowers), function(i){
      if (is.finite(uppers[i])) {
        sum(ODD$Population[which(hazMax > lowers[i] & hazMax <= uppers[i])], na.rm = TRUE)
      } else {
        sum(ODD$Population[which(hazMax > lowers[i])], na.rm = TRUE)  # >9
      }
    }, numeric(1))
    
    # Assign population columns
    impact_total[, pop_cols] <- as.list(pop_vals)
    
    # Buildings (if available)
    if ('nBuildings' %in% names(ODD)){
      build_vals <- vapply(seq_along(lowers), function(i){
        if (is.finite(uppers[i])) {
          sum(ODD$nBuildings[which(hazMax > lowers[i] & hazMax <= uppers[i])], na.rm = TRUE)
        } else {
          sum(ODD$nBuildings[which(hazMax > lowers[i])], na.rm = TRUE)
        }
      }, numeric(1))
      impact_total[, build_cols] <- as.list(build_vals)
    }
    
    # Clean & keep metadata
    impact_total$polygon <- NULL
    impact_total[, c('polygon', 'qualifier', 'overlap', 'coverage', 'inferred')] <- NULL
    impact_total$sdate <- min(ODD@hazdates)
    
    impact_all %<>% rbind(impact_total)
  }
  
  write_feather(impact_all, 'IIDIPUS_Input_Alternatives/Sept25/event_summaries0.5Intervals')
}


input_folder <- "IIDIPUS_Input_Alternatives/Sept25/"
collect_data_0.5Intervals <- function(input_folder, dat='all'){
  folderin <- paste0(dir, input_folder, "ODDobjects/")
  ufiles <- na.omit(list.files(path = folderin, pattern = Model$haz, recursive = TRUE, ignore.case = TRUE))
  if (tolower(dat) == 'train') {
    ufiles <- grep('^Train/', ufiles, value = TRUE)
  } else if (tolower(dat) == 'test') {
    ufiles <- grep('^Test/', ufiles, value = TRUE)
  }
  ufiles <- rev(ufiles)
  
  # helper: nice lower-bound labels (4, 4.5, ..., 8.5) + final 9plus
  fmt <- function(x) sub("\\.0$", "", format(round(x, 1), nsmall = 1))
  edges  <- c(seq(4, 9, by = 0.5), Inf)
  lowers <- edges[-length(edges)]
  uppers <- edges[-1]
  pop_cols   <- c(paste0("Pop", fmt(lowers[is.finite(uppers)])), "Pop9plus")
  build_cols <- c(paste0("Build", fmt(lowers[is.finite(uppers)])), "Build9plus")
  
  impact_all <- data.frame(
    event_id = integer(), iso3 = character(), sdate = as.Date(character()),
    impact = character(), observed = numeric(),
    # columns created later
    stringsAsFactors = FALSE
  )
  
  for (filer in ufiles){
    ODD <- readODD(paste0(folderin, filer))
    
    impact_total <- retrieveTotalImpact(ODD)
    if (NROW(dplyr::filter(impact_total, impact == 'mortality')) == 0) {  # incomplete coverage
      if (as.numeric(sub(".*_([0-9]+).*", "\\1", filer)) == 89){
        impact_total <- ODD@impact %>%
          dplyr::group_by(impact) %>%
          dplyr::summarise(polygon = 0, observed = sum(observed), iso3 = 'IDN',
                           qualifier = 'total', overlap = 0, coverage = 1, inferred = FALSE,
                           .groups = "drop")
      } else {
        stop("Missing total mortality impact; encountered non-special case.")
      }
    } else {
      message(filer)
    }
    
    impact_total$event_id <- as.numeric(sub(".*_([0-9]+).*", "\\1", filer))
    
    # Hazard maxima
    haz_idx <- grep('hazMean', names(ODD))
    if (length(haz_idx) > 1){
      hazMax <- apply(values(ODD[[haz_idx]]), 1, max, na.rm = TRUE)
    } else {
      hazMax <- values(ODD[[haz_idx]])
    }
    
    impact_total[paste0('Pop', c(seq(4, 8.5, 0.5), '9plus'))] <- NA
    impact_total[paste0('Build', c(seq(4, 8.5, 0.5), '9plus'))] <- NA
    
    # Compute population in 0.5 bins
    pop_vals <- vapply(seq_along(lowers), function(i){
      if (is.finite(uppers[i])) {
        sum(ODD$Population[which(hazMax > lowers[i] & hazMax <= uppers[i])], na.rm = TRUE)
      } else {
        sum(ODD$Population[which(hazMax > lowers[i])], na.rm = TRUE)  # >9
      }
    }, numeric(1))
    
    # Assign population columns
    impact_total[, pop_cols] <- as.list(pop_vals)
    
    # Buildings (if available)
    if ('nBuildings' %in% names(ODD)){
      build_vals <- vapply(seq_along(lowers), function(i){
        if (is.finite(uppers[i])) {
          sum(ODD$nBuildings[which(hazMax > lowers[i] & hazMax <= uppers[i])], na.rm = TRUE)
        } else {
          sum(ODD$nBuildings[which(hazMax > lowers[i])], na.rm = TRUE)
        }
      }, numeric(1))
      impact_total[, build_cols] <- as.list(build_vals)
    }
    
    # Clean & keep metadata
    impact_total$polygon <- NULL
    impact_total[, c('polygon', 'qualifier', 'overlap', 'coverage', 'inferred')] <- NULL
    impact_total$sdate <- min(ODD@hazdates)
    
    impact_all %<>% rbind(impact_total)
  }
  
  write_feather(impact_all, 'IIDIPUS_Input_Alternatives/Sept25/event_summaries0.5Intervals')
}

input_folder <- "IIDIPUS_Input_Alternatives/Sept25/"
collect_data_disagg <- function(input_folder, dat='all'){
  folderin<-paste0(dir,input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  ufiles = rev(ufiles)
  
  impact_all = data.frame(event_id = integer(), iso3=character(), sdate=as.Date(character()), impact=character(), observed=numeric(),
                          Pop4=numeric(), Pop5=numeric(), Pop6=numeric(), Pop7=numeric(), Pop8=numeric(), Pop9=numeric(),
                          Build4=numeric(), Build5=numeric(), Build6=numeric(), Build7=numeric(), Build8=numeric(), Build9=numeric())
  for (filer in ufiles){
    ODD <- readODD(paste0(folderin, filer))
    
    impact_mort = ODD@impact %>% filter(impact == 'mortality')
    # impact_total <- retrieveTotalImpact(ODD)
    # if (NROW(impact_total %>% filter(impact=='mortality'))==0){ #incomplete coverage, assess manually whether full impact is included in subnational:
    #   if( as.numeric(sub(".*_([0-9]+).*", "\\1", filer))==89){
    #     impact_total = ODD@impact %>% group_by(impact) %>% summarise( polygon=0, observed = sum(observed), iso3='IDN', qualifier='total', overlap=0, coverage=1, inferred=F)
    #   } else {
    #     stop()
    #   }
    # }  else {
    #   print(filer)
    # } 
    impact_mort$event_id = as.numeric(sub(".*_([0-9]+).*", "\\1", filer))
    #store the number of people exposed to intensity 7 or larger in each pixel
    
    impact_mort[paste0('Pop', 4:9)] <- NA
    impact_mort[paste0('Build', 4:9)] <- NA
    if (length(grep('hazMean', names(ODD))) > 1){
      hazMax <- apply(values(ODD[[grep('hazMean', names(ODD))]]),1 ,max, na.rm=T)
    } else {
      hazMax <- values(ODD[[grep('hazMean', names(ODD))]])
    }
    
    for (i in 1:NROW(impact_mort)){
      poly_indexes = ODD@polygons[[impact_mort$polygon[i]]]$indexes
      impact_mort$Pop4[i] <- sum(ODD$Population[intersect(poly_indexes, which(hazMax>4 & hazMax < 5))], na.rm=T)
      impact_mort$Pop5[i] <- sum(ODD$Population[intersect(poly_indexes,which(hazMax>5 & hazMax < 6))], na.rm=T)
      impact_mort$Pop6[i] <- sum(ODD$Population[intersect(poly_indexes,which(hazMax>6 & hazMax < 7))], na.rm=T)
      impact_mort$Pop7[i] <- sum(ODD$Population[intersect(poly_indexes,which(hazMax>7 & hazMax < 8))], na.rm=T)
      impact_mort$Pop8[i] <- sum(ODD$Population[intersect(poly_indexes,which(hazMax>8 & hazMax < 9))], na.rm=T)
      impact_mort$Pop9[i] <- sum(ODD$Population[intersect(poly_indexes,which(hazMax>9))], na.rm=T)
    }
    
    #if (ODD@impact$observed[i]==2152) stop()
    
    # if('nBuildings' %in% names(ODD)){
    #   impact_total$Build4 <- sum(ODD$nBuildings[which(hazMax>4 & hazMax < 5)], na.rm=T)
    #   impact_total$Build5 <- sum(ODD$nBuildings[which(hazMax>5 & hazMax < 6)], na.rm=T)
    #   impact_total$Build6 <- sum(ODD$nBuildings[which(hazMax>6 & hazMax < 7)], na.rm=T)
    #   impact_total$Build7 <- sum(ODD$nBuildings[which(hazMax>7 & hazMax < 8)], na.rm=T)
    #   impact_total$Build8 <- sum(ODD$nBuildings[which(hazMax>8 & hazMax < 9)], na.rm=T)
    #   impact_total$Build9 <- sum(ODD$nBuildings[which(hazMax>9)], na.rm=T)
    # }
    
    impact_mort$event_id = as.numeric(sub(".*_([0-9]+).*", "\\1", filer))
    #impact_total$polygon = NULL
    impact_mort[,c('qualifier', 'overlap', 'coverage', 'inferred')] = NULL
    impact_mort$sdate = min(ODD@hazdates)
    impact_all %<>% rbind(impact_mort)
  }
  
  write_feather(impact_all, 'IIDIPUS_Input_Alternatives/Sept25/event_summaries_disagg')
  
}

# --------------------------------------------------------------------------------
#-------- For saving ODD objects as data frames using feather: -------------------
# --------------------------------------------------------------------------------
saveFeather <- function(){
  input_folder <- "IIDIPUS_Input_Alternatives/Sept25/"
  folderin<-paste0(dir,input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
  
  event_ids_all <- as.numeric(sub(".*_(\\d+)$", "\\1", ufiles))
  
  impact_observed = data.frame(event_i = numeric(), mort=integer())
  polygon_impacts = list()
  
  n <- 171
  
  data_grids <- data.frame(
    xmin = numeric(n),
    xmax = numeric(n),
    ymin = numeric(n),
    ymax = numeric(n),
    ncol = integer(n),
    nrow = integer(n),
    res_x = numeric(n),
    res_y = numeric(n),
    crs = rep(NA_character_, n),  # character column with NAs
    stringsAsFactors = FALSE
  )
  
  exposed_cells_list = list()
  
  for (i in event_ids_all){
    filer <- ufiles[which(event_ids_all==i)]
    ODD <- readODD(paste0(folderin, filer))
    
    if( as.numeric(sub(".*_([0-9]+).*", "\\1", filer))==32){
      ODD@impact = ODD@impact[1:3,]
    }
    
    impact_total = retrieveTotalImpact(ODD)
      
    if (NROW(impact_total %>% filter(impact=='mortality'))==0){ #incomplete coverage, assess manually whether full impact is included in subnational:
      if( as.numeric(sub(".*_([0-9]+).*", "\\1", filer))==89){
        impact_total = ODD@impact %>% group_by(impact) %>% summarise( polygon=0, observed = sum(observed), iso3='IDN', qualifier='total', overlap=0, coverage=1, inferred=F)
      }  else {
        stop()
      }#else if( as.numeric(sub(".*_([0-9]+).*", "\\1", filer))==32){
        #impact_total = ODD@impact %>% group_by(impact) %>% summarise( polygon=0, observed = 0, iso3='CHN', qualifier='total', overlap=0, coverage=1, inferred=F)
        #impact_total$observed = c(76615, 9000,3)
      #}
    } 
    
    impact_observed %<>% add_row(event_i = i,
                                 mort = (impact_total %>% filter(impact=='mortality'))$observed)
    
    impact_mort = ODD@impact %>% filter(impact=='mortality')    
    polygon_impacts[[i]] = list(event_id = i, 
                                polygons = list())
    
    
    ODD_df <- as.data.frame(ODD, na.rm=F)
    ODD_df$PDens = (log(ODD_df$PDens) - Model$center$PDens$mean)/Model$center$PDens$sd
    ODD_df$GNIc = (log(ODD_df$GNIc) - Model$center$GNIc$mean)/Model$center$GNIc$sd
    ODD_df$SHDI = (ODD_df$SHDI - Model$center$SHDI$mean)/Model$center$SHDI$sd
    ODD_df$Vs30 = (log(ODD_df$Vs30) - Model$center$Vs30$mean)/Model$center$Vs30$sd
    ODD_df$EQFreq = (log(ODD_df$EQFreq) - Model$center$EQFreq$mean)/Model$center$EQFreq$sd
    exposed_flag = !apply(ODD_df[, grep('hazMean', names(ODD_df)), drop=F],1, function(row) all(is.na(row)))
    
    #ODD_df$hazMax <- apply(ODD_df[,grep('hazMean', names(ODD_df)), drop=F], 1, max, na.rm=T)
    rows_keep = which(!is.na(ODD_df$Population) & !is.na(ODD_df$ISO3C) & !is.na(ODD_df$SHDI) & !is.na(ODD_df$EQFreq) & exposed_flag)
    
    exposed_cells_list[[i]] = rows_keep
    
    ODD_df$updated_id=NA
    ODD_df$updated_id[rows_keep] = 1:length(rows_keep)
    
    #get updated indexes corresponding to old polygons
    for (j in 1:NROW(impact_mort)){
      if (is.null(ODD@polygons[[impact_mort$polygon[j]]]$indexes)){
        indexes = 1:length(rows_keep)
      } else {
        indexes = ODD_df$updated_id[ODD@polygons[[impact_mort$polygon[j]]]$indexes]
        indexes = indexes[!is.na(indexes)]
        if (length(indexes)==0 & impact_mort$observed[j] > 0) print(event_id)
      }
      polygon_impacts[[i]]$polygons[[j]] = list(
        polygon_id = impact_mort$polygon[j],
        observed = impact_mort$observed[j],
        polygon_indexes = indexes
      )
    }
    
    
    ODD_df %<>% filter(!is.na(Population) & !is.na(ISO3C) & !is.na(SHDI) & !is.na(EQFreq) & exposed_flag)
    ODD_df = ODD_df[, - grep('hazSD', names(ODD_df))]
    
    write_feather(ODD_df, paste0(dir, input_folder, 'ODDobjects_feather/', filer))
    
    hour <- as.numeric(substr(ODD@hazinfo$eventtimes, 1, 2))
    night_flag <- ifelse(hour>=22 | hour < 6, 1, 0)
    haz_info_df = data.frame(haz_number=1:length(ODD@hazdates),
                             night=(night_flag - Model$center$Night$mean)/Model$center$Night$sd,
                             first_haz=(ifelse(ODD@hazinfo$first_event,1,0) - Model$center$FirstHaz$mean)/Model$center$FirstHaz$sd,
                             first_haz_night=(ODD@hazinfo$first_event*night_flag - Model$center$FirstHaz.Night$mean)/Model$center$FirstHaz.Night$sd)
    
    write_feather(haz_info_df, paste0(dir, input_folder, 'ODDobjects_feather/hazinfo/', filer))
    
    data_grids[i, ] = data.frame(
      xmin = ext(ODD)[1],
      xmax = ext(ODD)[2],
      ymin = ext(ODD)[3],
      ymax = ext(ODD)[4],
      ncol = ncol(ODD),
      nrow = nrow(ODD),
      res_x = res(ODD)[1],
      res_y = res(ODD)[2],
      crs = crs(ODD, proj=TRUE)
    )
    
  }
  write_feather(data_grids, paste0(dir, input_folder, "ODDobjects_feather/grid_metadata.feather"))
  write_json(exposed_cells_list, paste0(dir, input_folder, "ODDobjects_feather/exposed_cells.json"))
  
  write_feather(impact_observed, paste0(dir, input_folder, 'ODDobjects_feather/impact_all'))
  write_json(polygon_impacts, paste0(dir, input_folder, 'ODDobjects_feather/subnatImpacts.json'), pretty = TRUE, auto_unbox = TRUE)
  
}

ODD_all <- read_feather('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Nov24Agg/event_summaries')

ODD_all %<>% filter(impact=='mortality')
theta0 = -35
theta1 = 3.8
ev_err = 0
ODD_all$impact_pred = exp(theta0 + theta1 * 5 + ev_err)* ODD_all$Pop5 +
                          exp(theta0 + theta1 * 6 + ev_err)* ODD_all$Pop6 +
                          exp(theta0 + theta1 * 7 + ev_err)* ODD_all$Pop7 +
                          exp(theta0 + theta1 * 8 + ev_err)* ODD_all$Pop8 +
                          exp(theta0 + theta1 * 9 + ev_err)* ODD_all$Pop9

plot(log(ODD_all$observed+10), log(ODD_all$impact_pred+10))
abline(0,1)


#-------- Posterior Predictive Samples:

library(arrow)

df <- read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/posterior_samples_2025-09-04_05-32-49.pt')

posterior_samples = matrix(unlist(df), nrow = nrow(df), byrow=T)

AlgoResults = list(
  Omega_sample_phys = array(0, dim = c(nrow(posterior_samples), length(unlist(Model$skeleton)),2)),
  W = array(1, dim=c(nrow(posterior_samples), 2))
)

#add prior samples
n_x = NCOL(AlgoResults$Omega_sample_phys)
s <- 1
for (i in 1:length(Model$Priors)){
  for (j in 1:length(Model$Priors[[i]])){
    prior_dist <- Model$Priors[[i]][[j]]
    AlgoResults$Omega_sample_phys[,s,1] <- do.call(match.fun(paste0('r', prior_dist$dist)), c(list(n=NROW(AlgoResults$Omega_sample_phys)), prior_dist[2:length(prior_dist)]))
    s <- s + 1
  }
}

AlgoResults$Omega_sample_phys[, c(8, 7, 3, 4, 12, 13, 14, 15, 16, 17, 18, 19),2] = posterior_samples

source('RCode/ResultsAnalysisFunctions2.R')
p1 = plot_correlated_posteriors_ggplot(AlgoResults, pairings = rbind(c(3,4), c(8,7), c(15,12), c(13, 14), c(16, 17), c(18, 19)), s_finish=2, plots_ncol = 2)
p1

df_pred <- read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/predictive_samples_2025-09-04_05-32-49.pt')

predictive_samples = exp(matrix(unlist(df_pred), nrow = nrow(df_pred), byrow=T)) - 10

impact_data = read_feather('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/impact_all')

df_postpredictive = data.frame(event_id=NA,
                               iso3 = 'AAA',
                               polygon=0, 
                               impact='mortality',
                               observed = impact_data$mort)

y_all_ordered = unlist(read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/y_all_ordered'))
df_postpredictive$observed = exp(y_all_ordered)-10

sample_column_names <- paste0('sampled.', 1:nrow(predictive_samples))
df_postpredictive[sample_column_names] <- t(predictive_samples)
df_postpredictive$train_flag = c(rep('TRAIN', 113), rep('TEST', 57))

input_folder <- "IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/"
train_dir <- file.path(dir, input_folder, "Train")
test_dir <- file.path(dir, input_folder, "Test")

# List all files in Train and Test directories
train_files <- list.files(train_dir, full.names = TRUE)
test_files <- list.files(test_dir, full.names = TRUE)

# Combine Train and Test files
files <- c(train_files, test_files)

# Extract event IDs from filenames
event_ids <- as.integer(str_extract(basename(files), "\\d+$"))


p2 = plot_df_postpredictive(df_postpredictive %>% filter(train_flag=='TEST'),'mortality')  + guides(color="none") 

legend <- get_plot_component(p1[[1]] +
                               theme(legend.position="bottom",
                                     plot.margin = margin(5, 5, 5, 25)), 'guide-box', return_all=T)[[3]]

plot_grid(plot_grid(plot_grid(p1[[1]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"), 
                    p1[[2]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"),
                    p1[[3]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"), 
                    p1[[4]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"),
                    p1[[5]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"), 
                    p1[[6]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"), 
                    ncol=2,
                    labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'),
                    label_fontfamily = "Times New Roman",
                    label_fontface = "plain"),plot_grid(plot_spacer(),legend, nrow=1, rel_widths=c(0.1, 0.8)), nrow=2, ncol=1,rel_heights=c(0.95, 0.05)),
          p2+ theme(plot.margin = margin(10, 20, 10, 10)), labels=c('', '(g)'), ncol=2,
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

#SBI_National.pdf, 10 x 5


#----------------------------------------------------------------------------------------------------
#----------------------------- SUBNATIONAL POST PREDICTIVE PLOTS ------------------------------------
#----------------------------------------------------------------------------------------------------

library(arrow)

df <- read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/posterior_samples_subnat_2025-08-19_06-20-31.pt')

posterior_samples = matrix(unlist(df), nrow = nrow(df), byrow=T)

AlgoResults = list(
  Omega_sample_phys = array(0, dim = c(nrow(posterior_samples), length(unlist(Model$skeleton)),2)),
  W = array(1, dim=c(nrow(posterior_samples), 2))
)

#add prior samples
n_x = NCOL(AlgoResults$Omega_sample_phys)
s <- 1
for (i in 1:length(Model$Priors)){
  for (j in 1:length(Model$Priors[[i]])){
    prior_dist <- Model$Priors[[i]][[j]]
    AlgoResults$Omega_sample_phys[,s,1] <- do.call(match.fun(paste0('r', prior_dist$dist)), c(list(n=NROW(AlgoResults$Omega_sample_phys)), prior_dist[2:length(prior_dist)]))
    s <- s + 1
  }
}

AlgoResults$Omega_sample_phys[, c(8, 7, 3, 4, 12, 13, 14, 15, 16, 17, 18, 19),2] = posterior_samples

library(patchwork)
p1 = plot_correlated_posteriors_ggplot(AlgoResults, pairings = rbind(c(3,4), c(8,7), c(15,12), c(13, 14), c(16, 17), c(18, 19)), s_finish=2, plots_ncol = 2)
p1

df_pred <- read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/predictive_samples_subnat_2025-08-19_06-20-31.pt')

predictive_samples = exp(matrix(unlist(df_pred), nrow = nrow(df_pred), byrow=T)) - 10

y_all_ordered = unlist(read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/y_all_ordered_subnat'))
df_postpredictive = data.frame(event_id=NA,
                               iso3 = 'AAA',
                               polygon=0, 
                               impact='mortality',
                               observed = exp(y_all_ordered)-10)

#df_postpredictive$observed = exp(y_all_ordered)-10

sample_column_names <- paste0('sampled.', 1:nrow(predictive_samples))
df_postpredictive[sample_column_names] <- t(predictive_samples)
df_postpredictive$train_flag = c(rep('TRAIN', 576), rep('TEST', 649))

input_folder <- "IIDIPUS_Input_Alternatives/Sept25/ODDobjects_feather/"
train_dir <- file.path(dir, input_folder, "Train")
test_dir <- file.path(dir, input_folder, "Test")

# List all files in Train and Test directories
train_files <- list.files(train_dir, full.names = TRUE)
test_files <- list.files(test_dir, full.names = TRUE)

# Combine Train and Test files
files <- c(train_files, test_files)

# Extract event IDs from filenames
event_ids <- as.integer(str_extract(basename(files), "\\d+$"))


p2 = plot_df_postpredictive(df_postpredictive %>% filter(train_flag=='TEST'),'mortality')  + guides(color="none") 

p2

legend <- get_plot_component(p1[[1]] +
                               theme(legend.position="bottom",
                                     plot.margin = margin(5, 5, 5, 25)), 'guide-box', return_all=T)[[3]]

plot_grid(plot_grid(plot_grid(p1[[1]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"), 
                              p1[[2]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"),
                              p1[[3]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"), 
                              p1[[4]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"),
                              p1[[5]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"), 
                              p1[[6]]+ theme(plot.margin = margin(5, 5, 5, 25), legend.position = "none"), 
                              ncol=2,
                              labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'),
                              label_fontfamily = "Times New Roman",
                              label_fontface = "plain"),plot_grid(plot_spacer(),legend, nrow=1, rel_widths=c(0.1, 0.8)), nrow=2, ncol=1,rel_heights=c(0.95, 0.05)),
          p2+ theme(plot.margin = margin(10, 20, 10, 10)), labels=c('', '(g)'), ncol=2,
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

# plot_grid(plot_grid(p1[[1]]+ theme(plot.margin = margin(5, 5, 5, 25)), 
#                     p1[[2]]+ theme(plot.margin = margin(5, 5, 5, 25)),
#                     p1[[3]]+ theme(plot.margin = margin(5, 5, 5, 25)), 
#                     p1[[4]]+ theme(plot.margin = margin(5, 5, 5, 25)),
#                     p1[[5]]+ theme(plot.margin = margin(5, 5, 5, 25)), 
#                     p1[[6]]+ theme(plot.margin = margin(5, 5, 5, 25)), ncol=2,
#                     labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'),
#                     label_fontfamily = "Times New Roman",
#                     label_fontface = "plain"),
#           p2+ theme(plot.margin = margin(10, 20, 10, 10)), labels=c('', '(g)'), ncol=2,
#           label_fontfamily = "Times New Roman",
#           label_fontface = "plain")

#SBI_SubNational.pdf, 10 x 5
  
  