
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-07-30_185759'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-07-22_193115'
#results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-08-04_134932'
results_file <- '/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/abcsmc_2023-08-11_164714'

AlgoResults <- readRDS(results_file)

addAlgoParams <- function(AlgoResults){
  AlgoResults$s_finish <- which(is.na(AlgoResults$Omega_sample_phys[1,1,]))[1]-1
  AlgoResults$n_x <- length(unlist(Model$skeleton))
  AlgoResults$Npart <- NROW(AlgoResults$W)
  return(AlgoResults)
}

plot_density_vs_step = function(AlgoResults, Omega=NULL){
  AlgoResults %<>% addAlgoParams()
  par(mfrow=c(5,4))
  plot_titles <- names(unlist(Model$skeleton))
  plot_titles[1:8] <- paste0(ifelse(1:8 %% 2 == 0, 'sigma_', 'mu_'), rep(1:4,each=2))
  for (i in 1:length(unlist(Model$skeleton))){
    ymin= min(AlgoResults$Omega_sample_phys[,i,1:AlgoResults$s_finish])
    ymax= max(AlgoResults$Omega_sample_phys[,i,1:AlgoResults$s_finish])
    if(!is.null(Omega)){
      ymin <- min(ymin, unlist(Omega)[i])
      ymax <- max(ymax, unlist(Omega)[i])
    }
    plot(rep(1:AlgoResults$s_finish,each=AlgoResults$Npart), AlgoResults$Omega_sample_phys[,i,1:AlgoResults$s_finish], ylim=c(ymin, ymax), main=plot_titles[i], xlab='Step', ylab='')
    if (!is.null(Omega)){
      abline(h=unlist(Omega)[i], col='red')
    }
  }
  par(mfrow=c(1,1)) # 800 x 800
}

plot_d_vs_step = function(AlgoResults){
  AlgoResults %<>% addAlgoParams()
  
  par(mfrow=c(1,1))
  ymin=min(AlgoResults$d, na.rm=T)
  ymax=max(AlgoResults$d[which(is.finite(AlgoResults$d))], na.rm=T)
  plot(rep(1, AlgoResults$Npart), apply(adrop(AlgoResults$d[,,1, drop=F], drop=3), 1, median), xlim=c(1, AlgoResults$s_finish), ylim=c(ymin, ymax), xlab='Step', ylab='Median distance for each particle')
  for (s in 2:AlgoResults$s_finish){
    nonzero_weights <- which(AlgoResults$W[,s] != 0)
    points(rep(s, length(nonzero_weights)), apply(adrop(AlgoResults$d[nonzero_weights,,s, drop=F], drop=3), 1, median))
  }
}

plot_correlated_posteriors = function(AlgoResults, include_priors=T, Omega=NULL,
                                      pairings=rbind(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12))){
  AlgoResults %<>% addAlgoParams()
  post_samples <- AlgoResults$Omega_sample_phys[,,AlgoResults$s_finish]
  if (include_priors) prior_samples <- AlgoResults$Omega_sample_phys[,,82]
  
  par(mfrow=c(2,3))
  for (p in 1:NROW(pairings)){
    xmin= min(post_samples[,pairings[p,1]]); xmax= max(post_samples[,pairings[p,1]])
    ymin= min(post_samples[,pairings[p,2]]); ymax= max(post_samples[,pairings[p,2]])
    if(include_priors){
      xmin <- min(xmin, prior_samples[,pairings[p,1]]); xmax <- max(xmax, prior_samples[,pairings[p,1]])
      ymin <- min(ymin, prior_samples[,pairings[p,2]]); ymax <- max(ymax, prior_samples[,pairings[p,2]])
    }
    plot(post_samples[,pairings[p,1]], post_samples[,pairings[p,2]], 
         xlab=names(unlist(Model$skeleton))[pairings[p,1]], xlim=c(xmin, xmax),
         ylab=names(unlist(Model$skeleton))[pairings[p,2]], ylim=c(ymin, ymax))
    if (include_priors) points(prior_samples[,pairings[p,1]], prior_samples[,pairings[p,2]], col='blue')
    if (!is.null(Omega)){
      points(unlist(Omega)[pairings[p,1]], unlist(Omega)[pairings[p,2]], col='red', pch=4, cex=2, lwd=4)
    }
  }
  
}

sample_post_predictive <- function(AlgoResults, M, s, dat='Train', single_particle=F, particle_i = NULL, Omega=NULL, 
                                   return_type='Poly'){
  if (!single_particle){
    sampled_part <- sample(1:AlgoResults$Npart, M, prob=AlgoResults$W[, s], replace=T)
    impact_sample <- SampleImpact(dir = dir,Model = Model,
                                  proposed = AlgoResults$Omega_sample_phys[sampled_part[1],,s] %>% relist(Model$skeleton) %>% addTransfParams(), 
                                  AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                  dat=dat)
    poly_sampled <- impact_sample$poly[[1]][,c('event_id', 'iso3', 'sdate', 'polygon', 'impact', 'observed', 'sampled')]
    point_sampled <- impact_sample$point
    
    if (M>1){
      for (m in 2:M){
        impact_sample <- SampleImpact(dir = dir,Model = Model,
                                      proposed = AlgoResults$Omega_sample_phys[sampled_part[m],,s] %>% relist(Model$skeleton) %>% addTransfParams(), 
                                      AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 1) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                      dat=dat)
        poly_sampled <- cbind(poly_sampled, impact_sample$poly[[1]]$sampled)
        point_sampled <- cbind(point_sampled, impact_sample$point)
      }
    }
  } else {
    if(is.null(Omega)){
      sampled_part <- ifelse(is.null(particle_i), sample(1:AlgoResults$Npart, 1, prob=AlgoResults$W[, s],replace=T), particle_i)
      Omega <- AlgoResults$Omega_sample_phys[sampled_part[1],,s] %>% relist(Model$skeleton) 
    } 
    impact_sample <- SampleImpact(dir = dir,Model = Model,
                                  proposed = Omega %>% addTransfParams(), 
                                  AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), M) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                  dat=dat)
    poly_sampled <- impact_sample$poly[[1]][,c('event_id', 'iso3', 'sdate', 'polygon', 'impact', 'inferred', 'observed', 'sampled')]
    point_sampled <- impact_sample$point
    
    if (M>1){
      for (m in 2:M){
        poly_sampled <- cbind(poly_sampled, impact_sample$poly[[m]]$sampled)
      }
    }
  }
  
  df_poly <- poly_sampled
  names(df_poly)[grep("sampled", names(df_poly))] = paste0('sampled.', 1:M)
  df_poly$obs_id[order(df_poly$observed)] <- 1:NROW(df_poly)
  df_poly$crps <- NA
  
  for (i in 1:NROW(df_poly)){
    df_poly$crps[i] <- crps(log(as.numeric(df_poly[i, grep("sampled", names(df_poly))])+10), log(df_poly$observed[i]+10)) * 
      unlist(AlgoParams$kernel_sd)[df_poly$impact[i]]
  }
  
  if (tolower(dat)=='train'){
    df_poly$train_flag = 'TRAIN'
  } else if (tolower(dat)=='test'){
    df_poly$train_flag = 'TEST'
  } else {
    folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
    ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
    ufiles_train <- grep('^Train/' , ufiles, value = TRUE)
    i_train <- as.numeric(gsub(".*_(\\d+)", "\\1", ufiles_train))
    df_poly$train_flag <- ifelse(df_poly$event_id %in% i_train, 'TRAIN', 'TEST')
  }
  
  if(tolower(return_type=='poly')){
    return(df_poly)
  } else {
    df_point = point_sampled
    colnames(df_point) = 'count'
    df_point <- cbind(df_point, 0)
    
    df_point[,2] = ifelse(rownames(df_point) %in% c('N12', 'N21', 'N23', 'N32'), 0.5 * df_point[,1], ifelse(rownames(df_point) %in% c('N13', 'N31'), df_point[,1], 0)) * 0.1
    return(list(poly=df_poly, point=df_point))
  }
}

compare_CRPS_breakdown_at_different_s <- function(AlgoResults, M){
  M <- 1
  AlgoResults %<>% addAlgoParams()
  df_poly1 <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, single_particle=T, particle_i = 1)
  df_poly2 <- sample_post_predictive(AlgoResults, M, 20, single_particle=T, particle_i = 1)
  
  
  ix <- which(df_poly1$impact=='mortality')
  ymin=min(df_poly1$crps[ix], df_poly2$crps[ix])
  ymax=max(df_poly1$crps[ix], df_poly2$crps[ix])
  plot(df_poly1$crps[ix], xlim=c(310, 350))
  points(df_poly2$crps[ix], col='red')
  
}

plot_predictive <- function(AlgoResults, dat='Train'){
  
  M <- 1
  AlgoResults %<>% addAlgoParams()
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat=dat, single_particle=T, particle=particle_min.d)
  
  impact_type <- 'displacement'
  k <- 10
  ix <- which(df_poly$impact==impact_type)
  inferred_cols <- ifelse(df_poly$inferred[ix], 'red', 'black')
  plot(log(df_poly$observed[ix]+k), log(df_poly$sampled.1[ix]+k), 
       xlab=paste0('log(observed+',k,')'), ylab=paste0('log(sampled+',k,')'), col=inferred_cols)
  abline(0,1)
  
  #large_discrepancies
  df_poly[which(abs(log(df_poly$sampled.1+1)-log(df_poly$observed+1)) >7),] %>% filter(impact==impact_type)
  
  #plot covariates against discrepancy between sampled and observed
  df_poly$sampled_median <- apply(df_poly[grep("sampled", names(df_poly))], 1, median)
  df_poly <- add_covar(df_poly, covar='EQFreq', dat=dat)
  plot(df_poly$EQFreq, df_poly$sampled_median)
  
  abline(0,1)
}

plot_predictive_train_vs_test <- function(AlgoResults){
  
  M <- 1
  AlgoResults %<>% addAlgoParams()
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  df_poly_train <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat='Train', single_particle=T, particle=particle_min.d)
  df_poly_test <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat='Test', single_particle=T, particle=particle_min.d)
  
  impact_type <- 'buildDest'
  k <- 10
  ix_train <- which(df_poly_train$impact==impact_type)
  df_poly <- rbind(df_poly_train, df_poly_test)
  ix <- which(df_poly$impact==impact_type)
  cols = c(rep('black', length(ix_train)), rep('red', length(ix)-length(ix_train)))
  plot(log(df_poly$observed[ix]+k), log(df_poly$sampled.1[ix]+k), xlab=paste0('log(observed+',k,')'), ylab=paste0('log(sampled+',k,')'), col=cols)
  abline(0,1)
  
  #large_discrepancies
  df_poly[which(abs(log(df_poly$sampled.1+1)-log(df_poly$observed+1)) >7),] %>% filter(impact==impact_type)
}

add_covar <- function(df_poly, covars='EQFreq', dat='all'){

  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  for (covar in covars){
    df_poly[[covar]] <- NA
  }
  
  for(i in 1:NROW(df_poly)){
    file_match <- grep(paste0("_", df_poly$event_id[i], "\\b"),  ufiles, value = TRUE)
    if (length(file_match) > 1) stop('Multiple ODD files for event with the same ID.')
    ODDy <- readRDS(paste0(folderin,file_match))
    poly_indexes <- ODDy@polygons[[df_poly$polygon[i]]]$indexes
    if ('hazMean' %in% covars){
      poly_pop_restricted <- poly_indexes[which(ODDy@data$Population[poly_indexes] > 1000)]
      df_poly[['hazMean']][i] <- max(ODDy@data[poly_pop_restricted,grep("hazMean",names(ODDy),value = T)], na.rm=T)
    } 
    #df_poly[[covar]][i] <- median(ODDy@data[[covar]][poly_indexes])
    #find mode:
    for (covar in covars[which(covars != 'hazMean')]){
      tbl <- table(ODDy@data[[covar]][poly_indexes])
      if (length(tbl)==0){df_poly[[covar]][i] <- NA; next }
      df_poly[[covar]][i] <- as.numeric(names(tbl)[which.max(tbl)])
    }
  }
  return(df_poly)
}


# ------------- Regionalise countries according to PAGER proposed regionalisation --------------

region_list <- list(
  list(name = "Australia, USA, and Canada", countries = c("Australia", "Canada", "United States (without California)", "Mexico", "Saint Pierre and Miquelon", "United States Minor Outlying Islands")),
  list(name = "New Zealand and California", countries = c("New Zealand", "California (USA)")),
  list(name = "Central America", countries = c("Costa Rica", "Panama")),
  list(name = "South Central America", countries = c("Dominican Republic", "Jamaica", "Guadeloupe", "El Salvador")),
  list(name = "Caribbean and Central America", countries = c("Guatemala", "Belize", "Honduras", "Nicaragua", "Haiti", "Puerto Rico", "Cayman Islands", "Turks and Caicos Islands", "Anguilla", "Montserrat", "Cuba", "Bahamas", "Saint Kitts and Nevis", "Saint Lucia", "Antigua and Barbuda", "Trinidad and Tobago", "Aruba", "Netherlands Antilles", "Dominica", "Grenada", "Saint Vincent and the Grenadines", "Martinique", "British Virgin Islands", "U.S. Virgin Islands", "Barbados", "Saint Barthelemy", "Saint Martin (France)")),
  list(name = "Western South America", countries = c("Colombia", "Ecuador", "Peru", "Chile", "Argentina", "Uruguay", "Brazil", "Paraguay")),
  list(name = "Eastern South America", countries = c("Venezuela", "Bolivia", "Brazil", "Uruguay", "Guyana", "Suriname", "Paraguay", "French Guiana")),
  list(name = "North Africa", countries = c("Algeria", "Egypt", "Tunisia", "Western Sahara")),
  list(name = "South-central Africa", countries = c("Botswana", "Namibia", "South Africa", "Swaziland", "Zimbabwe", "Morocco", "Sudan", "Chad", "Central African Republic", "Cameroon", "Congo", "DRP Congo", "Gabon", "Equatorial Guinea", "Sao Tome and Principe", "Angola", "Mauritania", "Senegal", "Gambia", "Guinea-Bissau", "Sierra Leone", "Liberia", "Cote d'Ivoire", "Ghana", "Togo", "Benin", "Niger", "Nigeria", "Mali", "Burkina Faso", "Guinea", "Yemen", "Eritrea", "Djibouti", "Ethiopia", "Somalia", "Kenya", "Uganda", "Rwanda", "Burundi", "United Republic of Tanzania", "Malawi", "Madagascar", "Mozambique", "Zambia", "Lesotho")),
  list(name = "Italy", countries = c("Italy", "Holy See", "Malta", "San Marino")),
  list(name = "Northern Europe", countries = c("Norway", "Sweden", "Finland", "Denmark", "Germany", "Belgium", "France", "Austria", "Switzerland", "Aland Islands", "Monaco", "Poland", "Bouvet Island", "United Kingdom", "Ireland", "Guernsey", "Isle of Man", "Jersey", "Falkland Islands (Malvinas)", "Saint Helena", "South Georgia and the South Sandwich Islands", "Iceland", "Faroe Islands", "Greenland", "Svalbard and Jan Mayen", "Liechtenstein", "Luxembourg", "Netherlands", "Greece", "Spain", "Portugal", "Gibraltar", "Cape Verde", "Andorra")),
  list(name = "Eastern Europe", countries = c("Czech Republic", "Slovenia", "Slovakia", "Hungary", "Bosnia and Herzegovina", "Croatia", "Serbia", "Montenegro", "Romania", "Albania", "Former Yugoslav Republic of Macedonia", "Bulgaria", "Republic of Moldova")),
  list(name = "Baltic States and Russia", countries = c("Estonia", "Latvia", "Lithuania", "Belarus", "Ukraine", "Russian Federation", "Georgia", "Armenia", "Azerbaijan")),
  list(name = "Central Asia", countries = c("Kazakhstan", "Uzbekistan", "Turkmenistan", "Kyrgyzstan", "Tajikistan")),
  list(name = "Arabian Peninsula", countries = c("Turkey", "Oman", "United Arab Emirates", "Qatar", "Saudi Arabia", "Bahrain", "Kuwait", "Lebanon", "Jordan", "Palestinian Territory", "Syrian Arab Republic", "Israel", "Cyprus", "Libyan Arab Jamahiriya")),
  list(name = "Iran & Iraq", countries = c("Iran", "Iraq", "Afghanistan", "Pakistan")),
  list(name = "Chinese Peninsula", countries = c("Brunei Darussalam", "China", "North Korea", "South Korea", "Macao", "Mongolia")),
  list(name = "Philippines and Malaysian Peninsula", countries = c("Singapore", "Thailand", "Hong Kong", "Malaysia", "Philippines")),
  list(name = "Indian Peninsula", countries = c("India", "Sri Lanka", "Bangladesh", "Nepal", "Bhutan", "Myanmar")),
  list(name = "Indonesian Peninsula", countries = c("Cambodia", "Laos", "Timor-Leste", "Vietnam", "Papua New Guinea", "American Samoa", "Samoa", "Tokelau", "Tuvalu", "Fiji", "Tonga", "Vanuatu", "Wallis and Futuna", "Niue", "Nauru", "New Caledonia", "Solomon Islands", "Palau", "Guam", "Northern Mariana Islands", "Marshall Islands", "Federated States of Micronesia", "Kiribati", "Cook Islands", "French Polynesia", "Norfolk Island", "Pitcairn", "British Indian Ocean Territory", "Christmas Island", "Cocos (Keeling) Islands", "French Southern Territories", "Heard Island and McDonald Islands", "Maldives", "Comoros", "Mauritius", "Mayotte", "Reunion", "Seychelles", "Bermuda", "Antarctica", "Indonesia")),
  list(name = "Japan & South Korea", countries = c("Japan", "Taiwan"))
)

#simplify region list given that we just know country, and so that it matches output of countrycode()
region_list[[1]]$countries[3] <- 'United States'
region_list[[2]]$countries <- region_list[[2]]$countries[1]
region_list[[9]]$countries[which(region_list[[9]]$countries=='United Republic of Tanzania')] = 'Tanzania'
region_list[[9]]$countries[which(region_list[[9]]$countries=='DRP Congo')] = 'Congo - Kinshasa'
region_list[[12]]$countries[which(region_list[[12]]$countries=='Former Yugoslav Republic of Macedonia')] = 'North Macedonia'
region_list[[15]]$countries[11] <- 'Syria'
region_list[[19]]$countries[which(region_list[[19]]$countries=='Myanmar')] = 'Myanmar (Burma)'


iso3_to_region <- function(iso3, region_list) {
  #having some issues with getting some country name matches with countrycode() but this seems to work ok for now:
  if (iso3 %in% c('TTO', 'BIH')){
    country <- countrycode(iso3, origin='iso3c', destination='iso.name.en')
  } else {
    country <- ifelse(iso3=='KOS', 'Serbia', countrycode(iso3, origin='iso3c', destination='country.name'))
  }
  for (i in seq_along(region_list)) {
    sublist <- region_list[[i]]$countries
    if (country %in% sublist) {
      return(region_list[[i]]$name)
    }
  }
  return(NA)  # Country not found in the list
}

add_regions <- function(df_poly){
  df_poly$region <- NA
  for (i in 1:NROW(df_poly)){
    iso3=df_poly$iso3[i]
    if (df_poly$iso3[i] == 'TOT'){
      df_poly_event <- filter(df_poly, event_id==df_poly$event_id[i] & iso3 != 'TOT')
      if (NROW(df_poly_event)==0){
        if (df_poly$event_id[i]==81){ iso3='PNG'}
        if (df_poly$event_id[i]==90){ iso3='CRI'}
        if (iso3=='TOT') {stop(df_poly$event_id[i])}
      } else {
        iso3=df_poly_event$iso3[which.max(df_poly_event$observed)]
      }
    }
    df_poly$region[i] <- iso3_to_region(iso3, region_list)
  }
  return(df_poly)
}

add_landslide_flag <- function(poly_df){
  EQIL <- read.csv('/home/manderso/Downloads/EQIL_Database_2022.csv', stringsAsFactors=FALSE, fileEncoding="latin1")[, c('Event', 'Country', 'Year', 'Month', 'Day', 'Mw', 'Fault.Type', 'Depth..km.', 'Ms', 'Landslide.Fatalities')]
  df_poly$landslide_flag <- F
  df_poly$Landslide.Fatalities <- NA
  df_poly$Mw <- NA
  df_poly$Fault.Type <- NA
  df_poly$Depth <- NA
  df_poly$Ms <- NA
  
  EQIL$date <- as.Date(paste(EQIL$Year, EQIL$Month, EQIL$Day, sep = "-"))
  EQIL %<>% filter(date > '2007-01-01')
  for (i in unique(df_poly$event_id)){
    EQIL_match <- c()
    corr_rows <- which(df_poly$event_id==i)
    date_match <- which(EQIL$date > df_poly$sdate[corr_rows[1]] - 1 & EQIL$date < df_poly$sdate[corr_rows[1]] + 1)
    for(j in date_match){
      if(any(countrycode(df_poly$iso3[corr_rows], origin='iso3c', destination='country.name') == EQIL$Country[j])){
        EQIL_match <- c(EQIL_match, j)
      }
    }
    if (length(EQIL_match)>1){
      stop()
    }
    if (length(EQIL_match)==1){
      df_poly$landslide_flag[corr_rows] <- T
      df_poly$Landslide.Fatalities[corr_rows] <- EQIL$Landslide.Fatalities[EQIL_match]
      df_poly$Mw[corr_rows] <- EQIL$Mw[EQIL_match]
      df_poly$Fault.Type[corr_rows] <- EQIL$Fault.Type[EQIL_match]
      df_poly$Depth[corr_rows] <- EQIL$Depth..km.[EQIL_match]
      df_poly$Ms[corr_rows] <- EQIL$Ms[EQIL_match]
    }
  }
  return(df_poly)
}

add_hazard_info <- function(df_poly){
  folderin <- paste0(dir,"IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12/HAZARDobjects_additionalInfo/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
  
  df_poly$depth <- NA
  df_poly$max_mmi <- NA
  df_poly$magnitude <- NA
  
  for(i in 1:NROW(df_poly)){
    file_match <- grep(paste0("_", df_poly$event_id[i], "\\b"),  ufiles, value = TRUE)
    HAZARDobj <- readRDS(paste0(folderin, file_match))
    max_mmi_i <- which.max(HAZARDobj$hazard_info$max_mmi)
    df_poly$max_mmi[i] <- HAZARDobj$hazard_info$max_mmi[max_mmi_i]
    df_poly$depth[i] <- HAZARDobj$hazard_info$depth[max_mmi_i]
    df_poly$magnitude[i] <- max(HAZARDobj$hazard_info$magnitude)
  }
  
  return(df_poly)
}

addSHDI <- function(df_poly){ 
  GDLdata <- readGlobalDataLab()
  df_poly$shdi <- NA
  df_poly$healthindex <- NA
  df_poly$ExpSchYrs <- NA
  df_poly$sgdi <- NA
  for (i in 1:NROW(df_poly)){
    GDLdata_i <- which(GDLdata$AveSchYrs==df_poly$AveSchYrs[i] & GDLdata$LifeExp==df_poly$LifeExp[i] & GDLdata$GNIc==df_poly$GNIc[i])
    if (length(GDLdata_i)> 1){
      stop(paste('Duplicate Matches:', i))
    } 
    if (length(GDLdata_i)==0){
      print(paste('No Match', df_poly$iso3[i]))
      next
    }
    df_poly$shdi[i] <- GDLdata[GDLdata_i, 'shdi']
    df_poly$healthindex[i] <- GDLdata[GDLdata_i, 'healthindex']
    df_poly$ExpSchYrs[i] <- GDLdata[GDLdata_i, 'ExpSchYrs']
    df_poly$sgdi[i] <- GDLdata[GDLdata_i, 'sgdi']
  }
  
}



# -----------------------------------------------------------------------------------------------

manual_modify_params = function(AlgoResults){
  M=3
  AlgoResults %<>% addAlgoParams()
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  particle_min.d <- 40
  omega1 <- AlgoResults$Omega_sample_phys[particle_min.d, , AlgoResults$s_finish] %>% relist(Model$skeleton)
  omega1$eps$local = 0
  omega1$eps$hazard = 0
  omega2 <- omega1
  omega2$vuln_coeff$AveSchYrs = -0.1
  omega2$vuln_coeff$GNIc <- 0
  omega2$Lambda2$mu <- 11.68
  #omega2$vuln_coeff$GNIc = -0.05
  Model$HighLevelPriors(omega1 %>% addTransfParams(), Model)
  Model$HighLevelPriors(omega2 %>% addTransfParams(), Model)
  
  df_all1 <- sample_post_predictive(AlgoResults, 3, AlgoResults$s_finish, dat='train', single_particle=T, Omega = omega1)
  df_all2 <- sample_post_predictive(AlgoResults, 3, AlgoResults$s_finish, dat='train', single_particle=T, Omega = omega2)
  
  df_poly1 <- df_all1$poly
  df_poly2 <- df_all2$poly
  sum(df_poly1$crps); sum(df_poly2$crps)
  
  df_point1 <- df_all1$point
  df_point2 <- df_all2$point
  sum(df_point1[,2]); sum(df_point2[,2])
  
  impact_sample1 <- SampleImpact(dir = dir,Model = Model,
                                proposed = omega1 %>% addTransfParams(), 
                                AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 3) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                dat='Train')
  
  impact_sample2 <- SampleImpact(dir = dir,Model = Model,
                                proposed = omega2 %>% addTransfParams(), 
                                AlgoParams = AlgoParams %>% replace(which(names(AlgoParams)==c('m_CRPS')), 3) %>% replace(which(names(AlgoParams)==c('Np')), 1),
                                dat='Train')
  
  sum(df_poly1$crps)
  sum(df_poly2$crps)
  
  dist1 <- CalcDist(impact_sample1, AlgoParams)
  dist2 <- CalcDist(impact_sample2, AlgoParams)
  hist(AlgoResults$d[,1,AlgoResults$s_finish])
  Model$HighLevelPriors(omega2 %>% addTransfParams(), Model, AlgoParams)
  
  ix <- which(impact_sample2$poly[[1]]$impact=='mortality')
  plot(log(impact_sample1$poly[[1]]$observed[ix] + k), log(impact_sample1$poly[[1]]$sampled[ix]+k) - (log(impact_sample1$poly[[1]]$observed[ix] + k)))
  points(log(impact_sample1$poly[[1]]$observed[ix] + k), log(impact_sample2$poly[[1]]$sampled[ix]+k) - (log(impact_sample2$poly[[1]]$observed[ix] + k)), col='red')
  
  
}



plot_covar_vs_error = function(AlgoResults, covar='EQFreq', dat='all'){
  #plot covariates against discrepancy between sampled and observed
  
  M <- 5
  AlgoResults %<>% addAlgoParams()
  particle_min.d <- which(AlgoResults$d[,,AlgoResults$s_finish] == min(AlgoResults$d[which(AlgoResults$W[, AlgoResults$s_finish] > 0),,AlgoResults$s_finish]), arr.ind=T)
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish, dat=dat, single_particle=T, particle_i = particle_min.d)
  df_poly$sampled_median <- apply(df_poly[grep("sampled", names(df_poly))], 1, median)
  df_poly$sampled_mean <- apply(df_poly[grep("sampled", names(df_poly))], 1, mean)
  
  impact_type <- 'mortality'
  k <- 10
  
  ggplot(df_poly %>% filter(impact==impact_type), aes(x= log(observed+k), y=log(sampled_mean+k))) +
    geom_point(aes(color=train_flag)) +theme_minimal()

  # df_poly %>% filter(magnitude> 7.5 & ((log(sampled_median+k)-log(observed+k))>2))
  df_poly %>% filter(impact==impact_type)
  
  df_poly %<>% add_hazard_info()
  df_poly %<>% add_covar(covars=c('hazMean', 'EQFreq', 'GNIc', 'Vs30', 'AveSchYrs', 'LifeExp'), dat=dat)
  
  ggplot(df_poly %>% filter(impact==impact_type), aes(x= log(observed+k), y=log(sampled_median+k))) + 
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color=AveSchYrs))  +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  #df_poly %>% filter(magnitude < 6.5 & (log(sampled_median+k) - log(observed+k) < -2.5) & impact==impact_type)
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = log(EQFreq+0.1), y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = log(EQFreq))) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type & train_flag=='TRAIN'), aes(x = log(GNIc), y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = log(GNIc))) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type & train_flag=='TRAIN'), aes(x = AveSchYrs, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = AveSchYrs)) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type & train_flag=='TRAIN'), aes(x = LifeExp, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = LifeExp)) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = Vs30, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = Vs30)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn"))) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = LifeExp, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point(aes(color = train_flag, shape=train_flag))
  
  
  
  ggplot(df_poly %>% filter(impact == impact_type & train_flag=='TRAIN'), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(color = "black", size = 1, stroke = 1) +
    geom_point(aes(color = AveSchYrs)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn"))) + theme_minimal()
  
  
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = log(depth), y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(aes(color = magnitude)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn"))) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(aes(color = hazMean)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlGn"))) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) + geom_point() +
    geom_point(aes(color = log(depth))) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = magnitude, y = log(sampled_median + k) - log(observed + k))) + geom_point()
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = max_mmi, y = log(sampled_median + k) - log(observed + k))) + geom_point()
  
  df_poly %<>% add_covar(covar='hazMean', dat=dat)
  df_poly %<>% add_covar(covar='AveSchYrs', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k))) +
    geom_point(aes(color = AveSchYrs)) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdYlGn")) + theme_minimal()
  
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed + k)))
  
  df_poly %<>% add_covar(covar='GNIc', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = log(GNIc), y = log(sampled_median + k) - log(observed + k))) + 
    geom_point()
  
  df_poly %<>% add_covar(covar='LifeExp', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = LifeExp, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point()
  
  df_poly %<>% add_covar(covar='AveSchYrs', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = AveSchYrs, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point()
  
  df_poly %<>% add_covar(covar='Vs30', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = Vs30, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point(aes(color=as.factor(event_id))) + theme(legend.position='none')
  
  df_poly %<>% add_covar(covar='EQFreq', dat=dat)
  ggplot(df_poly %>% filter(impact == impact_type), aes(x = EQFreq, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point(aes(color=as.factor(event_id))) + theme(legend.position='none')
  
  #df_poly[ix[which(7.7 < df_poly$hazMean[ix] & df_poly$hazMean[ix] < 8.3)],]
  #which(log(df_poly$sampled_median[ix]+k) - log(df_poly$observed[ix]+k) > 2.5)
  
  #how reliable is landslide data?
  df_poly_landslide = add_landslide_flag(df_poly)
  ggplot(df_poly_landslide %>% filter(impact == impact_type), aes(x = Depth, y = log(sampled_median + k) - log(observed + k))) + 
    geom_point(aes(color=landslide_flag))
  
  # ggplot(df_poly_landslide %>% filter(impact == impact_type), aes(x = hazMean, y = log(sampled_median + k) - log(observed -ifelse(is.na(Landslide.Fatalities), 0,Landslide.Fatalities) + k))) + 
  #   geom_point(aes(color=landslide_flag))
  
  
  
  df_poly %<>% add_regions()
  ggplot(df_poly %>% filter(impact==impact_type), aes(x=region, y=log(observed+k)-log(sampled_median+k))) +
    geom_point(aes(color=as.factor(event_id))) + ggtitle(impact_type) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='none')
  
  + 
    geom_point(aes(shape=as.factor(event_id), color=as.factor(event_id)))
  
  plot(df_poly$EQFreq[ix], log(df_poly$sampled_median[ix]+k)-log(df_poly$observed[ix]+k), xlab=covar, ylab='log(Sampled+10) - log(Observed+10)')
  
  ggplot(df_poly, aes(x=GNIc, y=log(observed+k)-log(sampled_median+k))) + geom_point() 
  #ggplot(df_poly, aes(x=countrycode(iso3, origin='iso3c', destination='continent'), y=log(observed+k)-log(sampled_median+k))) + geom_point()
  
  
  
  
  
  abline(0,0)
}

plot_impact_curves = function(AlgoResults){
  #plot covariates against discrepancy between sampled and observed
  
  M <- 5
  AlgoResults %<>% addAlgoParams()
  
  n_post_samples <- 500
  particle_sample <- sample(1:AlgoResults$Npart, n_post_samples, prob=AlgoResults$W[,AlgoResults$s_finish], replace=T)
  I <- seq(4.5, 5.5, 0.01)
  Omega <- relist(AlgoResults$Omega_sample_phys[particle_sample[1],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()
  #fDamUnscaled(I,list(I0=Model$I0, Np=length(I)),relist(AlgoResults$Omega_sample_phys[particle_sample[1],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()) 
  Damage <- h_0(I,Model$I0, Omega)
  plot(I, D_MortDisp_calc(Damage, Omega)[2,], type='l', xlim=c(4.5, 5.5), ylim=c(0,0.3))
  lines(I, D_MortDisp_calc(Damage, Omega)[1,], col='red')
  for (i in 2:n_post_samples){
    Omega <- relist(AlgoResults$Omega_sample_phys[particle_sample[i],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()
    #fDamUnscaled(I,list(I0=Model$I0, Np=length(I)),relist(AlgoResults$Omega_sample_phys[particle_sample[1],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()) 
    Damage <- h_0(I,Model$I0, Omega)
    lines(I, D_MortDisp_calc(Damage, Omega)[2,])
    lines(I, D_MortDisp_calc(Damage, Omega)[1,], col='red')
    lines(I, D_DestDam_calc(Damage, Omega)[2,], col='blue')
    lines(I, D_DestDam_calc(Damage, Omega)[1,], col='green')
  }
}

plot_vuln = function(AlgoResults, dat='all'){
  
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  AlgoResults %<>% addAlgoParams()
  n_post_samples <- 5
  df_vuln <- data.frame(
    post_sample = integer(),
    iso3 = character(),
    vuln_overall = numeric(),
    date = as.Date(character()),  # Convert the date column to Date type
    PDens = numeric(),
    AveSchYrs = numeric(),
    LifeExp = numeric(),
    GNIc = numeric(),
    Vs30 = numeric(),
    EQFreq = numeric(),
    mag = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion of character to factors
  )
  for (i in 1:n_post_samples){
    particle_sample <- sample(1:AlgoResults$Npart, 1, prob=AlgoResults$W[,AlgoResults$s_finish], replace=T)
    Omega <- relist(AlgoResults$Omega_sample_phys[particle_sample[1],,AlgoResults$s_finish], Model$skeleton) %>% addTransfParams()
    EQfreqs <- c()
    for (file in ufiles){
      ODDy <- readRDS(paste0(folderin, file))
      pop_restricted = which(ODDy@data$Population > 1000)
      max_haz <- max(ODDy@data[pop_restricted, grep("hazMean",names(ODDy),value = T)], na.rm=T)
      max_haz_i <- which(ODDy@data[pop_restricted, grep("hazMean",names(ODDy),value = T), drop=F] == max_haz, arr.ind=T)[1]
      EQfreqs <- c(EQfreqs, log(ODDy@data[pop_restricted[max_haz_i], 'EQFreq']+0.1))
      vuln_overall <- GetLP_single(Omega, Model$center, vuln_terms=list(PDens=ODDy@data[pop_restricted[max_haz_i], 'PDens'], 
                                                        AveSchYrs=ODDy@data[pop_restricted[max_haz_i], 'AveSchYrs'],
                                                        LifeExp=ODDy@data[pop_restricted[max_haz_i], 'LifeExp'],
                                                        GNIc=ODDy@data[pop_restricted[max_haz_i], 'GNIc'],
                                                        Vs30=ODDy@data[pop_restricted[max_haz_i], 'Vs30'],
                                                        EQFreq=ODDy@data[pop_restricted[max_haz_i], 'EQFreq'],
                                                        Mag=6.187425))
      df_vuln %<>% add_row(post_sample=i, 
                           iso3=ODDy@data[pop_restricted[max_haz_i], 'ISO3C'], 
                           vuln_overall=vuln_overall, 
                           date=ODDy@hazdates[1],
                           PDens=Omega$vuln_coeff$PDens * ((log(ODDy@data[pop_restricted[max_haz_i], 'PDens']+1) - Model$center$PDens$mean)/Model$center$PDens$sd),
                           AveSchYrs=Omega$vuln_coeff$AveSchYrs * ((ODDy@data[pop_restricted[max_haz_i], 'AveSchYrs'] - Model$center$AveSchYrs$mean)/Model$center$AveSchYrs$sd),
                           LifeExp=Omega$vuln_coeff$LifeExp * ((ODDy@data[pop_restricted[max_haz_i], 'LifeExp'] - Model$center$LifeExp$mean)/Model$center$LifeExp$sd),
                           GNIc=Omega$vuln_coeff$GNIc * ((log(ODDy@data[pop_restricted[max_haz_i], 'GNIc']) - Model$center$GNIc$mean)/Model$center$GNIc$sd),
                           Vs30=Omega$vuln_coeff$Vs30 * ((ODDy@data[pop_restricted[max_haz_i], 'Vs30'] - Model$center$Vs30$mean)/Model$center$Vs30$sd),
                           EQFreq=Omega$vuln_coeff$EQFreq * ((log(ODDy@data[pop_restricted[max_haz_i], 'EQFreq']+0.1) - Model$center$EQFreq$mean)/Model$center$EQFreq$sd), 
                           mag=Omega$vuln_coeff$Mag * ((6.187425 - Model$center$Mag$mean)/Model$center$Mag$sd))
    }
  }
  plot(df_vuln$date, df_vuln$EQFreq, col='blue',  ylim=c(-0.3,0.3))
  points(df_vuln$date, df_vuln$GNIc, col='green')
  points(df_vuln$date, df_vuln$AveSchYrs)
  points(df_vuln$date, df_vuln$PDens, col='red')
  points(df_vuln$date, df_vuln$Vs30, col='yellow')
  points(df_vuln$date, df_vuln$LifeExp, col='orange')
  
}

plot_post_predictive = function(AlgoResults, M){
  #mapply(crps, df_plot[grep("^sampled\\.", names(df_sampled))], df_plot$observed)
  AlgoResults %<>% addAlgoParams()
  
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish)
  
  df_long_sampled <-   df_poly %>% pivot_longer(paste0('sampled.', 1:M))
  df_long_sampled$log_observed <- log(df_long_sampled$observed+10)
  df_long_sampled$log_sampled <- log(df_long_sampled$value+10)
  
  # Create the ggplot
  ggplot(df_long_sampled %>% filter(impact=='mortality', obs_id %in% 900:1000), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
  
  ggplot(df_long_sampled %>% filter(impact=='mortality', observed != 0, obs_id >1000), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
  
  geom_density(alpha = 0.5) + geom_vline(aes(xintercept = observed), col='red') + theme_minimal()
  
  ggplot(df_sampled, aes(x=)) + + geom_density()
  
}

model_deepdive = function(AlgoResults, M){
  #mapply(crps, df_plot[grep("^sampled\\.", names(df_sampled))], df_plot$observed)
  AlgoResults %<>% addAlgoParams()
  
  df_poly <- sample_post_predictive(AlgoResults, M, AlgoResults$s_finish)
  
  df_long_sampled <-   df_poly %>% pivot_longer(paste0('sampled.', 1:M))
  df_long_sampled$log_observed <- log(df_long_sampled$observed+10)
  df_long_sampled$log_sampled <- log(df_long_sampled$value+10)
  
  # Create the ggplot
  ggplot(df_long_sampled %>% filter(impact=='mortality', obs_id %in% 900:1000), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
  
  ggplot(df_long_sampled %>% filter(impact=='mortality', observed != 0, obs_id >1000), aes(x = as.factor(obs_id), y=log_sampled)) + 
    geom_violin(color = "darkgray", trim = FALSE, scale='width') + geom_point(aes(x=as.factor(obs_id), y=log_observed))
  
  geom_density(alpha = 0.5) + geom_vline(aes(xintercept = observed), col='red') + theme_minimal()
  
  ggplot(df_sampled, aes(x=)) + + geom_density()
  
}


plot_correlated_posteriors(AlgoResults, include_priors=T)
plot_density_vs_step(AlgoResults, Omega)



