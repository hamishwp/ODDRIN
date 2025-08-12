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

#input_folder <- "IIDIPUS_Input_Alternatives/July25Agg/"
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
    if (NROW(impact_total)==0) stop()
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
  
  write_feather(impact_all, 'IIDIPUS_Input_Alternatives/July25Agg/event_summaries')
  
}

# --------------------------------------------------------------------------------
#-------- For saving ODD objects as data frames using feather: -------------------
# --------------------------------------------------------------------------------
saveFeather <- function(){
  input_folder <- "IIDIPUS_Input_Alternatives/July25Agg/"
  folderin<-paste0(dir,input_folder, "ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
  
  event_ids_all <- as.numeric(sub(".*_(\\d+)$", "\\1", ufiles))
  
  impact_observed = data.frame(event_i = numeric(), mort=integer())
  SHDI_df = data.frame(event_i = numeric(), SHDI=numeric())
  
  data_grids = data.frame(xmin = numeric(),
                          xmax = numeric(),
                          ymin = numeric(),
                          ymax = numeric(),
                          ncol = integer(),
                          nrow = integer(),
                          res_x = numeric(),
                          res_y = numeric(),
                          crs = character())
  
  exposed_cells_list = list()
  
  for (i in event_ids_all){
    filer <- ufiles[which(event_ids_all==i)]
    ODD <- readODD(paste0(folderin, filer))
    
    impact_total = retrieveTotalImpact(ODD)
    
    impact_observed %<>% add_row(event_i = i,
                                 mort = (impact_total %>% filter(impact=='mortality'))$observed)
    
    ODD_df <- as.data.frame(ODD, na.rm=F)
    ODD_df$PDens = (log(ODD_df$PDens) - Model$center$PDens$mean)/Model$center$PDens$sd
    ODD_df$GNIc = (log(ODD_df$GNIc) - Model$center$GNIc$mean)/Model$center$GNIc$sd
    ODD_df$SHDI = (ODD_df$SHDI - Model$center$SHDI$mean)/Model$center$SHDI$sd
    ODD_df$Vs30 = (log(ODD_df$Vs30) - Model$center$Vs30$mean)/Model$center$Vs30$sd
    ODD_df$EQFreq = (log(ODD_df$EQFreq) - Model$center$EQFreq$mean)/Model$center$EQFreq$sd
    exposed_flag = !apply(ODD_df[, grep('hazMean', names(ODD_df)), drop=F],1, function(row) all(is.na(row)))
    
    #ODD_df$hazMax <- apply(ODD_df[,grep('hazMean', names(ODD_df)), drop=F], 1, max, na.rm=T)
    
    exposed_cells = which(!is.na(ODD_df$Population) & !is.na(ODD_df$ISO3C) & !is.na(ODD_df$SHDI) & !is.na(ODD_df$EQFreq) & exposed_flag)
    exposed_cells_list[[i]] = exposed_cells
    ODD_df = ODD_df[exposed_cells, ]
    ODD_df = ODD_df[, - grep('hazSD', names(ODD_df))]
    
    write_feather(ODD_df, paste0(dir, input_folder, 'ODDobjects_feather/', filer))
    
    hour <- as.numeric(substr(ODD@hazinfo$eventtimes, 1, 2))
    night_flag <- ifelse(hour>=22 | hour < 6, 1, 0)
    haz_info_df = data.frame(haz_number=1:length(ODD@hazdates),
                             night=(night_flag - Model$center$Night$mean)/Model$center$Night$sd,
                             first_haz=(ifelse(ODD@hazinfo$first_event,1,0) - Model$center$FirstHaz$mean)/Model$center$FirstHaz$sd,
                             first_haz_night=(ODD@hazinfo$first_event*night_flag - Model$center$FirstHaz.Night$mean)/Model$center$FirstHaz.Night$sd)
    
    write_feather(haz_info_df, paste0(dir, input_folder, 'ODDobjects_feather/hazinfo/', filer))
    
    data_grids %<>% rbind(data.frame(
      xmin = ext(ODD)[1],
      xmax = ext(ODD)[2],
      ymin = ext(ODD)[3],
      ymax = ext(ODD)[4],
      ncol = ncol(ODD),
      nrow = nrow(ODD),
      res_x = res(ODD)[1],
      res_y = res(ODD)[2],
      crs = crs(ODD, proj=TRUE)
    ))
    
    
  }
  # Save as Feather
  write_feather(data_grids, paste0(dir, input_folder, "grid_metadata.feather"))
  
  write_json(exposed_cells_list, paste0(dir, input_folder, "exposed_cells.json"))
  
  write_feather(impact_observed, paste0(dir, input_folder, 'ODDobjects_feather/impact_all'))
  
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



#-----------------------------------------------------------------------------------------------------------------
#-------------------------------------------- PLOT SBI RESULTS ---------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

library(arrow)

df <- read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/July25Agg/ODDobjects_feather/posterior_samples_2025-08-07_11-59-08.pt')

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

AlgoResults$Omega_sample_phys[, c(8, 7, 3, 4, 12, 13, 15, 16, 17, 18),2] = posterior_samples

p1 = plot_correlated_posteriors_ggplot(AlgoResults, pairings = rbind(c(3,4), c(8,7), c(15,12), c(13, 14), c(16, 17), c(18, 19)), s_finish=2, plots_ncol = 2)
p1

#-------- Posterior Predictive Samples:

df_pred <- read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/July25Agg/ODDobjects_feather/predictive_samples_2025-08-07_14-11-21.pt')

predictive_samples = exp(matrix(unlist(df_pred), nrow = nrow(df_pred), byrow=T)) - 10

impact_data = read_feather('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/July25Agg/ODDobjects_feather/impact_all')

df_postpredictive = data.frame(event_id=NA,
                               iso3 = 'AAA',
                               polygon=0, 
                               impact='mortality',
                               observed = impact_data$mort)

y_all_ordered = unlist(read_parquet('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/July25Agg/ODDobjects_feather/y_all_ordered'))
df_postpredictive$observed = exp(y_all_ordered)-10

sample_column_names <- paste0('sampled.', 1:nrow(predictive_samples))
df_postpredictive[sample_column_names] <- t(predictive_samples)
df_postpredictive$train_flag = c(rep('TRAIN', 112), rep('TEST', 57))

input_folder <- "IIDIPUS_Input_Alternatives/July25Agg/ODDobjects_feather/"
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

plot_grid(plot_grid(p1[[1]]+ theme(plot.margin = margin(5, 5, 5, 25)), 
                    p1[[2]]+ theme(plot.margin = margin(5, 5, 5, 25)),
                     p1[[3]]+ theme(plot.margin = margin(5, 5, 5, 25)), 
                     p1[[4]]+ theme(plot.margin = margin(5, 5, 5, 25)),
                     p1[[5]]+ theme(plot.margin = margin(5, 5, 5, 25)), 
                     p1[[6]]+ theme(plot.margin = margin(5, 5, 5, 25)), ncol=2,
             labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'),
             label_fontfamily = "Times New Roman",
            label_fontface = "plain"),
             p2+ theme(plot.margin = margin(10, 20, 10, 10)), labels=c('', '(g)'), ncol=2,
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

#SBI_National.pdf, 10 x 5