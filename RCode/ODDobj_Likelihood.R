
setGeneric("DispX_Likelihood", function(ODD,Omega,center, Method, latent_var_event, event_i=NA)
  standardGeneric("DispX_Likelihood") )
# Code that calculates/predicts the total human displacement 
setMethod("DispX_Likelihood", "ODD", function(ODD,Omega,center,
                                   Method=list(Np=20,cores=8,cap=-300, 
                                               impact_weights=list(displacement=1,mortality=7,buildDam=0.6), 
                                               kernel='energy_score'), latent_var_event,
                                   event_i = NA
                                
){
  # Samples from the model to simulate the impact of an event.
  # Output type options:
  #     - SampledFull: Simulates and returns the impact in every pixel
  #     - SampledAgg: Simulates the impact in every pixel, but returns a data frame only containing
  #                   simulations aggregated at the same level as the true observations
  #     - SampledTotal: Simulates the impact in every pixel, but returns a data frame containing
  #                   only the total for each impact type (aggregated over all pixels)
  #     - ODDwithSampled: Returns the ODD object with columns Mort, Disp, and BuildDam for the simulated impact
  
  #elapsed_time <- c()
  #start_time <- Sys.time()
  
  # Extract 0D parameters & speed up loop
  Params<-FormParams(ODD,list(Np=Method$Np,center=center))
  Params$I0 <- Model$I0 #some objects have different I0 but don't want this to affect the model
  
  # Income distribution percentiles & extract income percentile  
  SincN<-paste0('p',seq(10,80,10), 'p', seq(20,90,10))
  Sinc<-ExtractCIndy(ODD,var = SincN)
  
  ODD_df <- as.data.frame(ODD, na.rm=F)
  
  
  #ODD_df$ISO3C <- levels(ODD[['ISO3C']])[[1]]$VALUE[ODD_df$ISO3C]
  
  # Speed-up calculation (through accurate cpu-work distribution) to only values that are not NA
  notnans<-which(!(is.na(ODD_df$Population) | is.na(ODD_df$ISO3C) | is.na(ODD_df$SHDI)))
  
  # Calculate non-local linear predictor values
  LP<-GetLP(ODD_df,Omega,Params,Sinc,notnans, split_GNI=T)
  LP_buildings <- GetLP(ODD_df,Omega,Params,Sinc,notnans, split_GNI=F)
  
  #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, GetLP = finish_time-start_time); start_time <- Sys.time()
  
  BD_data_present <- !is.null(ODD_df$nBuildings)
  hrange<-grep("hazMean",names(ODD_df),value = T)
  hrange_order <- order(paste(ODD@hazinfo$eventdates, ODD@hazinfo$eventtimes))
  
  event_lp <- getLP_event(ODD@hazinfo, Omega, Params)
  
  cov_mort_disp = Omega$eps$hazard_cor * Omega$eps_adj$hazard_mort * Omega$eps_adj$hazard_disp
  cov_mort_bd = Omega$eps$hazard_cor * Omega$eps_adj$hazard_mort * Omega$eps_adj$hazard_bd
  cov_disp_bd = Omega$eps$hazard_cor * Omega$eps_adj$hazard_disp * Omega$eps_adj$hazard_bd
  covar_matrix = cbind(c(Omega$eps_adj$hazard_mort^2, cov_mort_disp, cov_mort_bd), c(0, Omega$eps_adj$hazard_disp^2, cov_disp_bd), c(0, 0, Omega$eps_adj$hazard_bd^2))
  covar_matrix[upper.tri(covar_matrix)] = covar_matrix[lower.tri(covar_matrix)]
  
  
  covar_matrix_local = covar_matrix * Omega$eps_adj$local
  # covar_matrix_local[1,2] = covar_matrix_local[2,1]  = 1 * sqrt(covar_matrix_local[1,1] * covar_matrix_local[2,2])
  # covar_matrix_local[1,3] = covar_matrix_local[3,1]  =  1 * sqrt(covar_matrix_local[1,1] * covar_matrix_local[3,3])
  # covar_matrix_local[2,3] = covar_matrix_local[3,2]  =  1 * sqrt(covar_matrix_local[2,2] * covar_matrix_local[3,3])
  # 
  # eps_local_ij <- array(0, dim=c(length(hrange), 3, Method$Np))
  # for (i in 1:Method$Np){
  #   eps_local_ij[,,i] <- rmvnorm(length(hrange), rep(0,3), sigma=covar_matrix_local)
  # }
  
  #eps_event <- array(0, dim=c(3, Method$Np))
  # if (is.na(event_i)){
  #   eps_event <- t(rmvnorm(Method$Np, rep(0, 3), sigma=covar_matrix))
  # } else {
  #   eps_event <-  chol(covar_matrix) %*% t(Omega$u[event_i,,])
  # }
  
  eps_event =  t(chol(covar_matrix)) %*% latent_var_event$eps_event

  eps_local = aperm(apply(latent_var_event$eps_local, c(1, 2), function(x) t(chol(covar_matrix_local)) %*% x), c(2,3,1))
  
  
  #REMOVE ONCE THIS IS SPECIFIED ELSEWHERE:
  #latent_var_event$obs <- array(0, dim=c(nrow(ODD_df), length(hrange), 3)) # (grid-cell, hazard, impact type)
  #latent_var_event$eps_local_ij <- rmvnorm(length(hrange)*nrow(ODD_df), rep(0,3), sigma=covar_matrix_local)
  
  #eps_local_ij <- aperm(array(eps_local_long, dim=c(length(hrange), Method$Np, 3)), c(1,3,2))
  ODD_df$Population <- round(ODD_df$Population)
  LL_Mort_sum <- 0
  LL_Disp_sum <- 0
  LL_BuildDam_sum <- 0
  for (h_i in 1:length(hrange)){
    h <- hrange[h_i]
    I_ij<-ODD_df[,h]
    ind <- intersect(notnans, which(!is.na(I_ij)))
    D <- fDamUnscaled(I_ij[ind],list(I0=Params$I0, Np=1),Omega) + LP_buildings[ind] + event_lp[h_i]
    p_MortDisp <- pmin(D_MortDisp_calc(D, Omega, t(sweep(eps_local[ind,h_i,1:2], 2, eps_event[1:2,], '+'))), 0.999) #First row of D_MortDisp is D_Mort, second row is D_Disp
    if (h_i == 1){
      LL_Mort = dbinom(latent_var_event$obs[ind,h_i, 1], ODD_df$Population[ind], p_MortDisp[1,], log=T)
      LL_Disp = dbinom(latent_var_event$obs[ind,h_i, 2], ODD_df$Population[ind] - latent_var_event$obs[ind,h_i, 1], p_MortDisp[2,], log=T)
    } else {
      LL_Mort = dbinom(latent_var_event$obs[ind,h_i, 1], ODD_df$Population[ind] - rowSums(latent_var_event$obs[ind, 1:(h_i-1), 1:2], dim=1), p_MortDisp[1,], log=T)
      LL_Disp = dbinom(latent_var_event$obs[ind,h_i, 2], ODD_df$Population[ind] - rowSums(latent_var_event$obs[ind, 1:(h_i-1), 1:2], dim=1) - latent_var_event$obs[ind,h_i, 1], p_MortDisp[2,], log=T)
    }
    if(BD_data_present){
      
      p_BuildDam <- pmin(D_Dam_calc(D, Omega, eps_local[ind,h_i,3] + eps_event[3,]), 0.999)
      
      if (h_i==1){
        LL_BuildDam <- dbinom(latent_var_event$obs[ind,h_i, 3], ODD_df$nBuildings[ind], p_BuildDam, log=T)
      } else {
        LL_BuildDam <- dbinom(latent_var_event$obs[ind,h_i, 3], ODD_df$nBuildings[ind] - rowSums(latent_var_event$obs[ind,1:(h_i-1), 3, drop=F]), p_BuildDam, log=T)
      }
    } else {
      LL_BuildDam=0
    }
    LL_Mort_sum <- LL_Mort_sum + sum(LL_Mort, na.rm=T)
    LL_Disp_sum <- LL_Disp_sum + sum(LL_Disp, na.rm=T)
    LL_BuildDam_sum <- LL_BuildDam_sum + sum(LL_BuildDam, na.rm=T)
  }
  return(sum(LL_Mort_sum, LL_Disp_sum, LL_BuildDam_sum, na.rm=T))
})
  
CalcLL <-function(dir,Model,proposed,AlgoParams, latent_var, dat='Train', LL_old = NULL, event_updated = NA){
  
  # Load ODD files
  folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects_Unweighted/")
  
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  # When using task parallelisation, put the heaviest files first for optimisation reasons
  x <- file.info(paste0(folderin,ufiles))
  ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
  
  ufiles <- ufiles[3]
  
  #For Continuous Ranked Probability Score / Energy Score need multiple samples per particle
  AlgoParams$Np <- 1
  
  tmpFn<-function(filer){
    # Extract the ODD object
    ODDy<-readODD(paste0(folderin,filer))
    
    # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
    #ODDy@fIndies<-Model$fIndies
    
    #ODDy@impact%<>%as.data.frame.list()
    
    #remove inferred building damage and displacement observations:
    #ODDy@impact <- ODDy@impact[!1:NROW(ODDy@impact) %in% which(ODDy@impact$impact %in% c('buildDam', 'displacement') & ODDy@impact$inferred == T),]
    
    ODDy@impact$event_id = as.numeric(gsub(".*_(-?\\d+)$", "\\1", filer))
    
    #DispX requires event_i if introducing correlation between error terms in subsequent model samples (i.e. for correlated MCMC)
    event_i = which(ufiles==filer) #ifelse(is.null(proposed$u), NA, which(ufiles==filer))
    
    latent_var_event = list(
      eps_event = latent_var$eps_event[event_i,],
      eps_local = latent_var$eps_local[[event_i]],
      obs = latent_var$obs[[event_i]]
    )
    
    # Apply DispX
    LL <- DispX_Likelihood(ODD = ODDy,Omega = proposed,center = Model$center, Method = AlgoParams, latent_var_event=latent_var_event, event_i=event_i)
    
    return(LL) 
  }
  
  if (AlgoParams$AllParallel){
    
    return(pmap(mclapply(X = ufiles,FUN = tmpFn,mc.cores = AlgoParams$cores), c)[[1]])
    
  } else { 
    #even when AlgoParams$cores = 1 the above will still work, but this can sometimes be useful for debugging:
    LL <- tmpFn(ufiles[1])
    for (file_i in 2:length(ufiles)){
      LL <- c(LL, tmpFn(ufiles[file_i]))
    }
  }
}

library(mc2d)

init_latent_var <- function(AlgoParams, dat='Train'){
  folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects_Unweighted/")
  
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  # When using task parallelisation, put the heaviest files first for optimisation reasons
  x <- file.info(paste0(folderin,ufiles))
  ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  ufiles <- ufiles[3]
  
  latent_var <- list()
  latent_var$eps_event <- matrix(rnorm(length(ufiles)*3), ncol=3)
  latent_var$eps_local <- vector('list', length(ufiles))
  n_pixels <- vector('list', length(ufiles))
  n_haz <- vector('list', length(ufiles))
  
  for (file_i in 1:length(ufiles)){
    ODD <- readODD(paste0(folderin, ufiles[file_i]))
    ODD_df <- as.data.frame(ODD, na.rm=F)
    n_pixels[[file_i]] = nrow(ODD_df)
    n_haz[[file_i]] = length(grep("hazMean",names(ODD_df),value = T))
    latent_var$eps_local[[file_i]] = array(rnorm(n_pixels[[file_i]] * n_haz[[file_i]] * 3), dim=c(n_pixels[[file_i]],  n_haz[[file_i]], 3))
  }
  
  for (file_i in 1:length(ufiles)){
    ODD <- readODD(paste0(folderin, ufiles[file_i]))
    #ODD <- removeWeights(ODD)
    ODD_df <- as.data.frame(ODD, na.rm=F)
    latent_var$obs[[file_i]] <- array(NA, dim=c(nrow(ODD_df), length(grep("hazMean",names(ODD_df),value = T)), 3))
    
    for (impact_type_i in seq_along(c('mortality', 'displacement', 'buildDam'))){
      if (impact_type_i==3 & is.null(ODD_df$nBuildings)) next
      impact_type = c('mortality', 'displacement', 'buildDam')[impact_type_i]
      impact_filt <- ODD@impact %>% filter(impact==impact_type)
      
      pixel_fixed <- rep(F, NROW(ODD_df))
      
      if(NROW(impact_filt)==0) next
      impact_filt = impact_filt[order(impact_filt$observed),]
      if (impact_type_i==3){ n_unaff = ODD_df$nBuildings
      } else if (impact_type_i == 2){
        if ( n_haz[[file_i]]==1){
          n_unaff= round(ODD_df$Population) - latent_var$obs[[file_i]][,,1]
          n_unaff[is.na(n_unaff)] = 0
        } else {
          n_unaff= round(ODD_df$Population) - rowSums(latent_var$obs[[file_i]][,,1], na.rm=T)
          n_unaff[is.na(n_unaff)] = 0
        }
      } else {
        n_unaff = round(ODD_df$Population)
        n_unaff[is.na(n_unaff)] = 0
      } 
      for (j in 1:NROW(impact_filt)){
        
        #if (impact_filt$observed[j]==560) stop()
        region_indexes_all = ODD@polygons[[impact_filt$polygon[j]]]$indexes
        region_indexes <- intersect(region_indexes_all, which(!pixel_fixed)) # unallocated
        if (length(region_indexes)==0) next
        if (impact_filt$observed[j] == 0){
          latent_var$obs[[file_i]][region_indexes,,impact_type_i] = 0
          pixel_fixed[region_indexes] = T
        } else {
          #print(sum(latent_var$obs[[file_i]][region_indexes,,impact_type_i], na.rm=T))
          obs_rem = impact_filt$observed[j] - sum(latent_var$obs[[file_i]][region_indexes_all,,impact_type_i], na.rm=T) # some observations already allocated to pixels via subnational
          people = unlist(lapply(1:length(region_indexes), function(x) rep(x, n_unaff[region_indexes][x])))
          if (length(people) < obs_rem){
            obs_rem = length(people)
            print(paste('Observed value exceeds possible pixel pops. Event = ',ufiles[file_i], '. Impact type =', impact_type))
            #warning(paste('Observed value exceeds possible pixel pops. Event = ',ufiles[file_i], '. Impact type =', impact_type))
          } 
          
          people_aff = sample(people, size=obs_rem, replace=F)
          tot_people_aff <- table(people_aff)
          latent_var$obs[[file_i]][region_indexes,1,impact_type_i] <- 0
          latent_var$obs[[file_i]][region_indexes[as.numeric(names(tot_people_aff))],1,impact_type_i] = as.numeric(tot_people_aff)
          n_unaff[region_indexes[as.numeric(names(tot_people_aff))]] = n_unaff[region_indexes[as.numeric(names(tot_people_aff))]] - as.numeric(tot_people_aff)
          pixel_fixed[region_indexes] = T
          
          if (sum(latent_var$obs[[file_i]][region_indexes_all,,impact_type_i],na.rm=T) > impact_filt$observed[j]){
            print(paste('More people/buildings affected than observed. Event = ',ufiles[file_i], '. Impact type =', impact_type))
            #warning(paste('More people/buildings affected than observed. Event = ',ufiles[file_i], '. Impact type =', impact_type))
          } else if (sum(latent_var$obs[[file_i]][region_indexes_all,,impact_type_i], na.rm=T) < impact_filt$observed[j]){
            print(paste('Less people/buildings affected than observed. Event = ',ufiles[file_i], '. Impact type =', impact_type))
            #warning(paste('Less people/buildings affected than observed. Event = ',ufiles[file_i], '. Impact type =', impact_type))
          }
        }
      }
      if (n_haz[[file_i]] > 1){
        impact_across_haz <- latent_var$obs[[file_i]][,1,impact_type_i]
        inds_nonzero <- which(!is.na(impact_across_haz) & impact_across_haz > 0)
        alloc_probs <- ODD_df[,grep('hazMean', names(ODD_df))]
        alloc_probs[is.na(alloc_probs)] <- 0.001
        latent_var$obs[[file_i]][inds_nonzero,,impact_type_i] = rmultinomial(length(inds_nonzero), impact_across_haz[inds_nonzero], exp(alloc_probs[inds_nonzero,]))
        latent_var$obs[[file_i]][which(impact_across_haz==0),,impact_type_i]=0
      }
    }
  }
  
  return(latent_var)
  #latent_var$eps_local <- rnorm(sum(unlist(n_pixels) * unlist(n_haz) * 3))
  
}

#remove weights
removeWeightsAll <- function(){
  folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects/")

  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  # When using task parallelisation, put the heaviest files first for optimisation reasons
  x <- file.info(paste0(folderin,ufiles))
  ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
  
  for (file_i in 1:length(ufiles)){
    filer <- ufiles[file_i]
    ODD <- readODD(paste0(folderin, filer))
    ODDunweighted <- tryCatch(removeWeights(ODD),error=function(e) filer)
    saveODD(ODDunweighted, paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Nov24Agg/ODDobjects_Unweighted2/',filer ))
  }
}

move_latent_obs <- function(obs, AlgoParams, dat='Train'){
  folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects_Unweighted/")
  
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  # When using task parallelisation, put the heaviest files first for optimisation reasons
  x <- file.info(paste0(folderin,ufiles))
  ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
  
  if (tolower(dat)=='train'){
    ufiles <- grep('^Train/' , ufiles, value = TRUE)
  } else if (tolower(dat)=='test'){
    ufiles <- grep('^Test/' , ufiles, value = TRUE)
  }
  
  ufiles <- ufiles[3]
  
  seed <- round(runif(1, 0, 100000))
  print(seed)
  set.seed(seed)
  
  obs_new = obs
  updated_file_i <- 1#sample(1:length(ufiles),20)
  for (file_i in updated_file_i){
    ODD <- readODD(paste0(folderin, ufiles[file_i]))
    #ODD <- removeWeights(ODD)
    ODD_df <- as.data.frame(ODD, na.rm=F)
    #obs_new[[file_i]] <- array(NA, dim=c(nrow(ODD_df), length(grep("hazMean",names(ODD_df),value = T)), 3))
    
    for (impact_type_i in seq_along(c('mortality', 'displacement', 'buildDam'))){ # does this order need to be random for it to be reversible?
      if (impact_type_i==3 & is.null(ODD_df$nBuildings)) next
      impact_type = c('mortality', 'displacement', 'buildDam')[impact_type_i]
      impact_filt <- ODD@impact %>% filter(impact==impact_type)
      if(NROW(impact_filt)==0) next
      impact_filt = impact_filt[order(impact_filt$observed),]
      
      pixel_fixed = rep(F, NROW(ODD_df))
      
      if (impact_type_i==3){ 
        n_unaff = ODD_df$nBuildings - rowSums(obs_new[[file_i]][,,3, drop=F])
      } else { #else if (impact_type_i == 2){
        # if ( length(grep("hazMean",names(ODD_df),value = T))==1){
        #   n_unaff = round(ODD_df$Population) - obs_new[[file_i]][,,1]
        #   n_unaff[is.na(n_unaff)] = 0
        # } else {
        #   n_unaff = round(ODD_df$Population) - rowSums(obs_new[[file_i]][,,1], na.rm=T)
        #   n_unaff[is.na(n_unaff)] = 0
        # }
      #} 
        n_unaff = round(ODD_df$Population) - rowSums(obs_new[[file_i]][,,1:2, drop=F], na.rm=T)
        n_unaff[is.na(n_unaff)] = 0
      } 
      for (j in 1:NROW(impact_filt)){
        region_indexes_all = ODD@polygons[[impact_filt$polygon[j]]]$indexes
        #region_indexes <- region_indexes_all
        #region_indexes <- intersect(region_indexes_all, which(is.na(obs_new[[file_i]][,1,impact_type_i]))) # unallocated
        region_indexes <- intersect(region_indexes_all, which(!pixel_fixed))
        if (length(region_indexes)==0) next
        if (impact_filt$observed[j] == 0){ #| impact_filt$observed[j] == sum(obs[[file_i]][region_indexes_all,,impact_type_i], na.rm=T)){
          obs_new[[file_i]][region_indexes,,impact_type_i] = 0
          pixel_fixed[region_indexes] = T
        } else {
          #if (impact_type_i==2) stop()
          #print(sum(latent_var$obs[[file_i]][region_indexes,,impact_type_i], na.rm=T))
          #indices <- which(is.numeric(obs[[file_i]][region_indexes,,impact_type_i, drop=F]), arr.ind=T)
          pixels_haz_long <- c(obs_new[[file_i]][region_indexes,,impact_type_i, drop=F])
          if (length(grep("hazMean",names(ODD_df),value = T)) > 1){
            pixels_long <- region_indexes[c(row(adrop(obs_new[[file_i]][region_indexes,,impact_type_i, drop=F], drop=3)))]
            haz_long <-c(col(adrop(obs_new[[file_i]][region_indexes,,impact_type_i, drop=F], drop=3)))
          } else {
            pixels_long <- region_indexes[1:length(obs_new[[file_i]][region_indexes,,impact_type_i])]
            haz_long <- rep(1, length(region_indexes))
          }
          
          #current_allocs = which(obs[[file_i]][region_indexes,,impact_type_i, drop=F] != 0, arr.ind=T)
          current_allocs = which(pixels_haz_long != 0)
          
          if (length(current_allocs)==0){
            pixel_fixed[region_indexes] = T  
            print(paste0('Event: ', ufiles[file_i], '. Requesting to move observations but no observations found'))
            next
          }
          #current_allocs_regions = pixels_long[current_allocs]
          #current_allocs_haz = haz_long[current_allocs]
          current_obs = pixels_haz_long[current_allocs]
          
          #sample 1/10th of the current obs:
          current_people = current_allocs[unlist(lapply(1:length(current_allocs), function(x) rep(x, current_obs[x])))]
          if (length(current_people)==ceil(length(current_people)/100)){
            moved_people = current_people
          } else {
            moved_people = sample(current_people, ceil(length(current_people)/100), replace=F)
          }
          
          # select where the moved people go:
          pop_range <- range(ODD_df[region_indexes, 'Population'], na.rm=T)
          haz_range <- range(ODD_df[region_indexes,grep("hazMean",names(ODD_df),value = T)], na.rm=T)
          pop_haz_mix <- (ODD_df[region_indexes, 'Population'] - pop_range[1])/diff(pop_range) + (ODD_df[region_indexes,grep("hazMean",names(ODD_df),value = T)] - haz_range[1])/diff(haz_range)
          
          #sorted_hazint <- order(unlist(ODD_df[region_indexes, grep("hazMean",names(ODD_df),value = T)])) # sort pixels by intensity # check this
          sorted_hazint <- order(pop_haz_mix) # sort pixels by intensity # check this
          
          #sorted_pixels <- intersect(pixels_long[sorted_pixels], region_indexes)
          
          if (impact_type_i == 3){
            sorted_hazint <- sorted_hazint[which(round(ODD_df$nBuildings)[pixels_long[sorted_hazint]]>0)]
          } else {
            sorted_hazint <- sorted_hazint[which(round(ODD_df$Population)[pixels_long[sorted_hazint]]>0)]
          }
          
          
          for (moved_person in moved_people){
            index_match = which(sorted_hazint == moved_person)
            if (length(sorted_hazint) < 41){
              possible_moves <- sorted_hazint
            } else if ((index_match-20)<1){
              possible_moves <- sorted_hazint[1:(index_match+20)]
              possible_moves <- c(possible_moves, rep(sorted_hazint[1],41-length(possible_moves)))
            } else if ((index_match+20)>length(sorted_hazint)){
              possible_moves <- sorted_hazint[(index_match-20):length(sorted_hazint)]
              possible_moves <- c(possible_moves, rep(sorted_hazint[length(sorted_hazint)],41-length(possible_moves)))
            } else {
              possible_moves <- sorted_hazint[(index_match-20):(index_match+20)]
            }
            #remove from possible moves all those that are at capacity:
            possible_moves <- possible_moves[which(n_unaff[pixels_long[possible_moves]] != 0)]
            if (length(possible_moves)==0){
              next
            }
            
            #print(range(ODD_df$hazMean1[pixels_long[possible_moves]]))
            
            if(impact_type_i == 3){
              reallocated_pixel <- sample(possible_moves, 1)
              #reallocated_pixel <- sample(possible_moves, 1, prob = ODD_df$nBuilding[pixels_long[possible_moves]])
            } else {
              reallocated_pixel <- sample(possible_moves, 1)
              #reallocated_pixel <- sample(possible_moves, 1, prob = ODD_df$Population[pixels_long[possible_moves]])
            }
            
            obs_new[[file_i]][cbind(pixels_long[reallocated_pixel],haz_long[reallocated_pixel],impact_type_i)] = obs_new[[file_i]][cbind(pixels_long[reallocated_pixel],haz_long[reallocated_pixel],impact_type_i)] + 1
            obs_new[[file_i]][cbind(pixels_long[moved_person],haz_long[moved_person],impact_type_i)] = obs_new[[file_i]][cbind(pixels_long[moved_person],haz_long[moved_person],impact_type_i)] - 1
            n_unaff[pixels_long[moved_person]] = n_unaff[pixels_long[moved_person]] + 1
            n_unaff[pixels_long[reallocated_pixel]] = n_unaff[pixels_long[reallocated_pixel]] - 1
            #if(any(round(ODD_df$Population[ind]) - rowSums(obs_new[[file_i]][ind, 1:(h_i-1), 1:2], dim=1) < 0, na.rm=T)) stop('Hi')
          }
          if (any(obs_new[[file_i]]<0, na.rm=T)) stop()
          #now just adjust for if there is too many / too few left. This is due to overlap between regions
          if (sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) > impact_filt$observed[j]){
            #stop()
            n_excess = sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) - impact_filt$observed[j]
            #obs_rem = impact_filt$observed[j] - sum(latent_var$obs[[file_i]][region_indexes_all,,impact_type_i], na.rm=T) # some observations already allocated to pixels via subnational
            
            pixels_haz_long <- c(obs_new[[file_i]][region_indexes,,impact_type_i, drop=F])
            if (length(grep("hazMean",names(ODD_df),value = T)) > 1){
              pixels_long <- region_indexes[c(row(adrop(obs_new[[file_i]][region_indexes,,impact_type_i, drop=F],drop=3)))]
              haz_long <-c(col(adrop(obs_new[[file_i]][region_indexes,,impact_type_i, drop=F], drop=3)))
            } else {
              pixels_long <- region_indexes[1:length(obs_new[[file_i]][region_indexes,,impact_type_i])]
              haz_long <- rep(1, length(region_indexes))
            }
            
            current_allocs = which(pixels_haz_long != 0)
            #current_allocs_regions = pixels_long[current_allocs]
            #current_allocs_haz = haz_long[current_allocs]
            current_obs = pixels_haz_long[current_allocs]
            
            people = current_allocs[unlist(lapply(1:length(current_allocs), function(x) rep(x, current_obs[x])))]
            
            if (length(people)<=n_excess){
              people_removed = people
            } else {
              people_removed = sample(people, n_excess, replace=F)
            }
            tot_people_removed <- table(people_removed)
            haz_pixel_removed <- as.numeric(names(tot_people_removed))
            obs_new[[file_i]][cbind(pixels_long[haz_pixel_removed],haz_long[haz_pixel_removed],impact_type_i)] = obs_new[[file_i]][cbind(pixels_long[haz_pixel_removed],haz_long[haz_pixel_removed],impact_type_i)] - as.numeric(tot_people_removed)
            n_unaff[pixels_long[haz_pixel_removed]] = n_unaff[pixels_long[haz_pixel_removed]] +  as.numeric(tot_people_removed)
            if (sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) != impact_filt$observed[j]){
              print(paste0('Event:', ufiles[file_i], '. Sum of new observed (', sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]),') does not equal true observed (', impact_filt$observed[j], ')'))
              #warning(paste0('Event:', ufiles, '. Sum of new observed does not equal true observed'))
            }
            if (any(obs_new[[file_i]]<0, na.rm=T)) stop('HI')
          } else if (sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) < impact_filt$observed[j]){
            #stop()
            n_short = impact_filt$observed[j] - sum(obs_new[[file_i]][region_indexes_all,,impact_type_i])
            #obs_rem = impact_filt$observed[j] - sum(latent_var$obs[[file_i]][region_indexes_all,,impact_type_i], na.rm=T) # some observations already allocated to pixels via subnational
            people = region_indexes[unlist(lapply(1:length(region_indexes), function(x) rep(x, n_unaff[region_indexes][x])))]
            if (length(people)<=n_short){
              people_added = people
            } else {
              people_added = sample(people, n_short, replace=F)
            }
            tot_people_added <- table(people_added)
            if (length(grep("hazMean",names(ODD_df),value = T)) > 1){
              haz_allocs = sample(length(grep("hazMean",names(ODD_df),value = T)), length(as.numeric(names(tot_people_added))), replace=T)
            } else {
              haz_allocs = 1
            }
            obs_new[[file_i]][cbind(as.numeric(names(tot_people_added)),haz_allocs,rep(impact_type_i, length(haz_allocs)))] = obs_new[[file_i]][cbind(as.numeric(names(tot_people_added)),haz_allocs,rep(impact_type_i, length(haz_allocs)))] + as.numeric(tot_people_added)
            n_unaff[as.numeric(names(tot_people_added))] = n_unaff[as.numeric(names(tot_people_added))] -  as.numeric(tot_people_added)
            
            if (sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) != impact_filt$observed[j]){
              print(paste0('Event:', ufiles[file_i], '. Sum of new observed (', sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]),') does not equal true observed (', impact_filt$observed[j], ')'))
              #warning(paste0('Event:', ufiles, '. New observed does not equal true observation'))
            }
            #stop('Have removed people')
          }
          
          pixel_fixed[region_indexes] = T  
          #if(any(ODD_df$Population[ind] - rowSums(obs_new[[file_i]][ind, 1:(h_i-1), 1:2], dim=1) < 0, na.rm=T)) stop()
          # obs_rem = impact_filt$observed[j] - sum(latent_var$obs[[file_i]][region_indexes_all,,impact_type_i], na.rm=T) # some observations already allocated to pixels via subnational
          # people = unlist(lapply(1:length(region_indexes), function(x) rep(x, max_capacity[region_indexes][x])))
          # if (length(people) < obs_rem){
          #   obs_rem = length(people)
          #   warning(paste('Observed value exceeds possible pixel pops. Event = ',ufiles[file_i], '. Impact type =', impact_type))
          # } 
          # 
          # people_aff = sample(people, size=obs_rem, replace=F)
          # tot_people_aff <- table(people_aff)
          # latent_var$obs[[file_i]][region_indexes,1,impact_type_i][is.na(latent_var$obs[[file_i]][region_indexes,1,impact_type_i])] <- 0
          # latent_var$obs[[file_i]][region_indexes[as.numeric(names(tot_people_aff))],1,impact_type_i] = latent_var$obs[[file_i]][region_indexes[as.numeric(names(tot_people_aff))],1,impact_type_i] + as.numeric(tot_people_aff)
          # max_capacity[region_indexes[as.numeric(names(tot_people_aff))]] = max_capacity[region_indexes[as.numeric(names(tot_people_aff))]] - as.numeric(tot_people_aff)
        }
      }
      
      # if (length(grep("hazMean",names(ODD_df),value = T)) > 1){
      #   impact_across_haz <- latent_var$obs[[file_i]][,1,impact_type_i]
      #   inds_nonzero <- which(!is.na(impact_across_haz) & impact_across_haz > 0)
      #   alloc_probs <- ODD_df[,grep('hazMean', names(ODD_df))]
      #   alloc_probs[is.na(alloc_probs)] <- 0.001
      #   latent_var$obs[[file_i]][inds_nonzero,,impact_type_i] = rmultinomial(length(inds_nonzero), impact_across_haz[inds_nonzero], exp(alloc_probs[inds_nonzero,]))
      #   latent_var$obs[[file_i]][which(impact_across_haz==0),,impact_type_i]=0
      # }
    }
  }
  
  return(obs_new)
}
  
#==========================================================================
# 
# if h length is 1:
#   region_indexes_all = ODD@polygons[[impact_filt$polygon[j]]]$indexes
# #region_indexes <- region_indexes_all
# #region_indexes <- intersect(region_indexes_all, which(is.na(obs_new[[file_i]][,1,impact_type_i]))) # unallocated
# region_indexes <- intersect(region_indexes_all, which(!pixel_fixed))
# if (length(region_indexes)==0) next
# if (impact_filt$observed[j] == 0){
#   latent_var$obs[[file_i]][region_indexes,,impact_type_i] = 0
#   pixel_fixed[region_indexes] = T
# } else {
#   #if (impact_type_i==2) stop()
#   #print(sum(latent_var$obs[[file_i]][region_indexes,,impact_type_i], na.rm=T))
#   current_allocs = region_indexes[which(obs[[file_i]][region_indexes,,impact_type_i, drop=F] != 0)]
#   current_obs = obs[[file_i]][current_allocs,, impact_type_i]
#   
#   #sample 1/10th of the current obs:
#   current_people = current_allocs[unlist(lapply(1:length(current_allocs), function(x) rep(x, current_obs[x])))]
#   if (length(current_people)==ceil(length(current_people)/10)){
#     moved_people = current_people
#   } else {
#     moved_people = sample(current_people, ceil(length(current_people)/10), replace=F)
#   }
#   
#   # select where the moved people go:
#   sorted_pixels <- order(ODD_df[, grep("hazMean",names(ODD_df),value = T)]) # sort pixels by intensity
#   sorted_pixels <- intersect(sorted_pixels, which(round(ODD_df$Population) > 0))
#   sorted_pixels <- intersect(sorted_pixels, region_indexes)
#   for (moved_person in moved_people){
#     index_match = which(sorted_pixels == moved_person)
#     if (length(sorted_pixels) < 21){
#       possible_moves <- sorted_pixels
#     } else if ((index_match-10)<1){
#       possible_moves <- sorted_pixels[1:(index_match+10)]
#       possible_moves <- c(possible_moves, rep(sorted_pixels[1],21-length(possible_moves)))
#     } else if ((index_match+10)>length(sorted_pixels)){
#       possible_moves <- sorted_pixels[(index_match-10):length(sorted_pixels)]
#       possible_moves <- c(possible_moves, rep(sorted_pixels[length(sorted_pixels)],21-length(possible_moves)))
#     } else {
#       possible_moves <- sorted_pixels[(index_match-10):(index_match+10)]
#     }
#     #remove from possible moves all those that are at capacity:
#     possible_moves <- possible_moves[which(n_unaff[possible_moves] != 0)]
#     reallocated_pixel <- sample(possible_moves, 1)
#     obs_new[[file_i]][reallocated_pixel,,impact_type_i] = obs_new[[file_i]][reallocated_pixel,,impact_type_i] + 1
#     obs_new[[file_i]][moved_person,,impact_type_i] = obs_new[[file_i]][moved_person,,impact_type_i] - 1
#     n_unaff[moved_person] = n_unaff[moved_person] + 1
#     n_unaff[reallocated_pixel] = n_unaff[reallocated_pixel] - 1
#   }
#   
#   #now just adjust for if there is too many / too few left. This is due to overlap between regions
#   if (sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) > impact_filt$observed[j]){
#     #stop()
#     n_excess = sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) - impact_filt$observed[j]
#     #obs_rem = impact_filt$observed[j] - sum(latent_var$obs[[file_i]][region_indexes_all,,impact_type_i], na.rm=T) # some observations already allocated to pixels via subnational
#     current_allocs = region_indexes[which(obs[[file_i]][region_indexes,,impact_type_i, drop=F] != 0)]
#     current_obs = obs[[file_i]][current_allocs]
#     people = current_allocs[unlist(lapply(1:length(current_allocs), function(x) rep(x, current_obs[x])))]
#     
#     if (length(people)==n_excess){
#       people_removed = people
#     } else {
#       people_removed = sample(people, n_excess, replace=F)
#     }
#     tot_people_removed <- table(people_removed)
#     obs_new[[file_i]][as.numeric(names(tot_people_removed)),,impact_type_i] = obs_new[[file_i]][as.numeric(names(tot_people_removed)),,impact_type_i] - as.numeric(tot_people_removed)
#     n_unaff[as.numeric(names(tot_people_removed))] = n_unaff[as.numeric(names(tot_people_removed))] +  as.numeric(tot_people_removed)
#     if (sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) != impact_filt$observed[j]){
#       warning('still an issue')
#     }
#   } else if (sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) < impact_filt$observed[j]){
#     #stop()
#     n_short = impact_filt$observed[j] - sum(obs_new[[file_i]][region_indexes_all,,impact_type_i])
#     #obs_rem = impact_filt$observed[j] - sum(latent_var$obs[[file_i]][region_indexes_all,,impact_type_i], na.rm=T) # some observations already allocated to pixels via subnational
#     people = region_indexes[unlist(lapply(1:length(region_indexes), function(x) rep(x, n_unaff[region_indexes][x])))]
#     if (length(people)==n_short){
#       people_added = people
#     } else {
#       people_added = sample(people, n_short, replace=F)
#     }
#     tot_people_added <- table(people_added)
#     obs_new[[file_i]][as.numeric(names(tot_people_added)),,impact_type_i] = obs_new[[file_i]][as.numeric(names(tot_people_added)),,impact_type_i] - as.numeric(tot_people_added)
#     n_unaff[as.numeric(names(tot_people_added))] = n_unaff[as.numeric(names(tot_people_added))] +  as.numeric(tot_people_added)
#     
#     if (sum(obs_new[[file_i]][region_indexes_all,,impact_type_i]) != impact_filt$observed[j]){
#       warning('still an issue')
#     }
#     stop('Have removed people')
#   }
#   
#   pixel_fixed[region_indexes] = T  
  
  
  #==========================================================================

HLPrior_samples <- readRDS(paste0(dir, 'IIDIPUS_Input/HLPriorSamples_SMCOut'))
propCOV <- cov(HLPrior_samples)/5
init_val_phys <- Proposed2Physical(HLPrior_samples[1,] %>% relist(skeleton=Model$skeleton) %>% unlist(), Model)

#init_val_phys <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/abcsmc_2024-12-10_060355_alphaAdaptive_M100_Npart1000NovAgg5_propCOVmult0.2_further')$Omega_sample_phys[1,,170]



#init_val_phys = unlist()
AlgoParams$input_folder <- 'IIDIPUS_Input_Alternatives/Nov24Agg/'
AlgoParams$rho <- 0.9
AlgoParams$lambda_rate <- 0.5
AlgoParams$alpha_star <- 0.1
AlgoParams$a <- 0.234
AlgoParams$v_0 <- 15
AlgoParams$lambda_min <- 0.05
AlgoParams$N_steps <- 10000

Model$HighLevelPriors <- function(Omega, Model){
  return(0)
}

options(warn=2)

#Method:
correlated_AMCMC_LL <- function(AlgoParams, Model, propCOV = NULL, init_val_phys = NULL, unfinished=F, oldtag=NULL, tag_notes=NULL){
  #Input: 
  # - AlgoParams: Parameters describing the ABC-MCMC Algorithm (e.g. the ABC rejection threshold for higher level priors)
  # - Model: Describes the data simulation and calculation of the distance measure 
  # - Unfinished: If TRUE, then include oldtag - the tag (end of the filename) of an unfinished ABC-SMC run to be completed.
  # Output:
  # - 
  # Details:
  # - 
  
  n_x <- length(Model$par_lb) #n_x = number of parameters
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  tag<-ifelse(is.null(tag_notes), tag, paste0(tag, '_', tag_notes))
  folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects_Unweighted/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
  n_events <- length(ufiles)
  
  c = 2.38^2 / n_x
  A = pnorm(AlgoParams$a/2)
  delta = (1-1/n_x)*(sqrt(2*pi)*exp((A^2)/2)/(2*A)) + 1/(n_x*AlgoParams$a*(1-AlgoParams$a))
  
  iter_func = function(s){
    return(floor(s/2))
  }
  
  AlgoResults <- list(
    #Omega_sample = array(NA, dim=c(AlgoParams$smc_Npart, n_x, AlgoParams$smc_steps)), #store sampled parameters on the transformed space
    input_folder = AlgoParams$input_folder,
    Omega_sample = array(NA, dim=c(n_x, AlgoParams$N_steps)),
    Omega_sample_phys = array(NA, dim=c(n_x, AlgoParams$N_steps)), #store sampled parameters on the untransformed space
    loss = array(Inf, dim=c( AlgoParams$N_steps)), #Distances
    #u_event = array(NA, dim=c(n_events, AlgoParams$m_CRPS*AlgoParams$Np, 3, 3)),
    #u_event_selected = array(NA, dim=c(3, AlgoParams$N_steps)),
    lambda_store=array(NA, AlgoParams$N_steps), 
    mu_store=array(NA, dim=c(n_x, AlgoParams$N_steps)), 
    Sigma_store=array(NA, dim=c(n_x, n_x, AlgoParams$N_steps)),
    accprob_store=array(NA, AlgoParams$N_steps),
    sampled_full = NULL,
    tolerancestore=array(NA, AlgoParams$smc_steps),
    essstore=array(NA, AlgoParams$smc_steps),
    accrate_store=array(NA, AlgoParams$smc_steps),
    propCOV=array(NA, dim=c(n_x, n_x, AlgoParams$smc_steps)),
    propCOV_multiplier=array(NA, AlgoParams$smc_steps),
    s_restart = round(5/((AlgoParams$a)*(1-AlgoParams$a))),
    lambda_start = 1
  )
  
  if (is.null(propCOV)){
    stop('Please provide initial covariance of perturbation kernel')
  }
  
  if(unfinished==F){ 
    #Initialize and perform sampling for s=1
    if (is.null(init_val_phys)){
      stop('Please provide initial value until we provide functionality to sample from prior')
    }
    # STEP 1:
    AlgoResults$Omega_sample_phys[,1] = unlist(init_val_phys)
    AlgoResults$Omega_sample[,1] = init_val_phys %>% Physical2Proposed(Model) %>% unlist()
    # AlgoResults$u[,,,1] = rnorm(length(AlgoResults$u[,,,1]))
    proposed =  AlgoResults$Omega_sample_phys[,1] %>% relist(skeleton=Model$skeleton)
    # proposed$u = AlgoResults$u[,,,1]
    
    #initialise latent var:
    latent_var = init_latent_var(AlgoParams)
    latent_var_local_dim = lapply(latent_var$eps_local, dim)
    latent_var_obs_store <- latent_var$obs
    
    AlgoResults$loss[1] = CalcLL(dir, Model, proposed %>% addTransfParams(), AlgoParams, latent_var) #SampleImpact(dir, Model, proposed %>% addTransfParams(), AlgoParams)
    #print(impact_sample)
    #impact_sample = CalcLL(dir, Model, proposed %>% addTransfParams, AlgoParams, latent_var) #SampleImpact(dir, Model, proposed %>% addTransfParams(), AlgoParams)
    #print(impact_sample)
    #AlgoResults$loss[1] = CalcDist(impact_sample, AlgoParams)[1]
    AlgoResults$lambda_store[1] = AlgoResults$lambda_start
    AlgoResults$mu_store[,1] = AlgoResults$Omega_sample[,1]
    AlgoResults$Sigma_store[,,1] = propCOV
    AlgoResults$accprob_store[1] = AlgoParams$alpha_star
    #AlgoResults$u_selected[,1] = c(AlgoResults$u[1,1,1,1],AlgoResults$u[2,2,2,1],AlgoResults$u[3,3,3,1])
    
    # STEP 2:
    
    upd_params_vs_latentvar <- runif(1)
    if (upd_params_vs_latentvar < 0.1){
      #update the parameters and but keep latent variables the same
      Omega_prop <- AlgoResults$Omega_sample[,1]#multvarNormProp(xt=AlgoResults$Omega_sample[,1], propPars= c * AlgoResults$lambda_store[1] * AlgoResults$Sigma_store[,,1]) #perturb the proposal
      Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
      latent_var_prop <- latent_var
    } else {
      #update the latent variables but keep the parameters the same 
      Omega_prop <- AlgoResults$Omega_sample[,1]
      Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
      latent_var_prop <- list()
      latent_var_prop$eps_event <- AlgoParams$rho * latent_var$eps_event + sqrt(1-AlgoParams$rho^2) * rnorm(length(latent_var$eps_event))
      latent_var_loc_prop <- AlgoParams$rho * unlist(latent_var$eps_local) + sqrt(1-AlgoParams$rho^2) * rnorm(length(unlist(latent_var$eps_local)))
      latent_var_prop$eps_local = mapply(function(vec, d) array(vec, dim = c(d[1], d[2], d[3])), 
                                         split(latent_var_loc_prop, rep(1:length(latent_var_local_dim), sapply(latent_var_local_dim, prod))), 
                                         latent_var_local_dim, 
                                         SIMPLIFY = FALSE)
      latent_var_prop$obs <- move_latent_obs(latent_var$obs, AlgoParams)
      loss_prop = CalcLL(dir, Model, proposed %>% addTransfParams(), AlgoParams, latent_var_prop)
      print(loss_prop>AlgoResults$loss[1])
      print(paste(2, loss_prop))
    }
    #epsilon <- rnorm(length(unlist(latent_var)))
    #u_prop <- AlgoParams$rho * c(AlgoResults$u[,,,1]) + sqrt(1-AlgoParams$rho^2) * epsilon
    
    
    HP<- Model$HighLevelPriors(Omega_prop_phys %>% addTransfParams(), Model)
    if (HP> AlgoParams$ABC & Model$higherpriors){
      #REJECT DUE TO HIGHER LEVEL PRIOR
      AlgoResults$accprob_store[2] = 0
      AlgoResults$Omega_sample[,2] = AlgoResults$Omega_sample[,1]
      AlgoResults$Omega_sample_phys[,2] = AlgoResults$Omega_sample_phys[,1]
      AlgoResults$loss[2] = AlgoResults$loss[1]
      latent_var = latent_var
      #AlgoResults$u[,,,2] = AlgoResults$u[,,,1]
      #AlgoResults$u_selected[,2] = AlgoResults$u_selected[,1]
    }  else {
      proposed = Omega_prop_phys
      #proposed$u = array(u_prop, dim=c(n_events, AlgoParams$m_CRPS*AlgoParams$Np, 3))
      # impact_sample <- SampleImpact(dir = dir,Model = Model,
      #                               proposed = proposed, 
      #                               AlgoParams = AlgoParams)
      # loss_prop <- CalcDist(impact_sample, AlgoParams)[1]
      loss_prop = CalcLL(dir, Model, proposed %>% addTransfParams(), AlgoParams, latent_var_prop)
      print(paste(2, loss_prop))
      
      #calculate the acceptance probability:
      #print(modifyAcc(Omega_prop, Omega_sample_s[n,], Model))
      min_loss <- min(loss_prop, AlgoResults$loss[1])
      AlgoResults$accprob_store[2] <- min(1, exp(loss_prop - min_loss)/exp(AlgoResults$loss[1] - min_loss) * modifyAcc(Omega_prop, AlgoResults$Omega_sample[,1], Model, AlgoResults$lambda_store[1] * AlgoResults$Sigma_store[,,1]))
      
      u <- runif(1)
      if(u < AlgoResults$accprob_store[2]){
        #ACCEPT
        print('Accepted')
        AlgoResults$Omega_sample[,2] = Omega_prop
        AlgoResults$Omega_sample_phys[,2] = unlist(Omega_prop_phys)
        AlgoResults$loss[2] = loss_prop
        latent_var = latent_var_prop
        #AlgoResults$u[,,,2] = proposed$u
        #AlgoResults$u_selected[,2] = c(proposed$u[1,1,1],proposed$u[2,2,2],proposed$u[3,3,3])
      }  else {
        #REJECT
        AlgoResults$Omega_sample[,2] = AlgoResults$Omega_sample[,1]
        AlgoResults$Omega_sample_phys[,2] = AlgoResults$Omega_sample_phys[,1]
        AlgoResults$loss[2] = AlgoResults$loss[1]
        #AlgoResults$u[,,,2] = AlgoResults$u[,,,1]
        #AlgoResults$u_selected[,2] = AlgoResults$u_selected[,1]
      }
    }
    
    AlgoResults$mu_store[,2] = (AlgoResults$Omega_sample[,2] +  AlgoResults$Omega_sample[,1])/2
    AlgoResults$Sigma_store[,,2] = 1 / (AlgoParams$v_0 + n_x + 3) * 
      ((AlgoResults$mu_store[,2] %*% t(AlgoResults$mu_store[,2]) + AlgoResults$mu_store[,1] %*% t(AlgoResults$mu_store[,1])) - 
         2*AlgoResults$mu_store[,2] %*% t(AlgoResults$mu_store[,2]) +
         (AlgoParams$v_0 + n_x + 1)* AlgoResults$Sigma_store[,,1] )
    AlgoResults$lambda_store[2] = max(AlgoParams$lambda_min, AlgoResults$lambda_store[1] * exp(delta / (2) * (AlgoResults$accprob_store[2]-AlgoParams$a)))
    
    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag))  
    
    s_start = 3
    
  } else { 
    #Collect relevant information from the unfinished sample
    UnfinishedAlgoResults <- retrieve_UnfinishedAlgoResults_AMCMC3(dir, oldtag, AlgoResults)
    AlgoResults <- UnfinishedAlgoResults$AlgoResults
    s_start <- UnfinishedAlgoResults$s_start
    saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag)) 
  }
  
  for (s in s_start:AlgoParams$N_steps){
    if (s %% 20 == 0) { 
      for (i in 1:length(latent_var_obs_store)){
        latent_var_obs_store[[i]] <- abind(latent_var_obs_store[[i]], latent_var$obs[[i]], along=4)
      }
      saveRDS(latent_var_obs_store, paste0(dir,"IIDIPUS_Results/latent_obs_",tag))
      saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag))
    } else if(s %% 10 == 0){
      saveRDS(AlgoResults, paste0(dir,"IIDIPUS_Results/mcmc_",tag, "_backup"))
    }
    
    if (runif(1) < -0.1){
      Omega_prop <- AlgoResults$Omega_sample[,s-1] #multvarNormProp(xt=AlgoResults$Omega_sample[,s-1], propPars= c * AlgoResults$lambda_store[s-1] * AlgoResults$Sigma_store[,,s-1]) #perturb the proposal
      Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
      latent_var_prop = latent_var
    } else {
      #Omega_prop <-  AlgoResults$Omega_sample[,s-1] #
      Omega_prop <- multvarNormProp(xt=AlgoResults$Omega_sample[,s-1], propPars= c * AlgoResults$lambda_store[s-1] * AlgoResults$Sigma_store[,,s-1]) #perturb the proposal
#AlgoResults$Omega_sample[,s-1]
      Omega_prop_phys <- Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model)
      #epsilon <- rnorm(length(AlgoResults$u[,,,ifelse((s-1) %% 3==0, 3, (s-1)%%3)]))
      #u_prop <- AlgoParams$rho * c(AlgoResults$u[,,,ifelse((s-1) %% 3==0, 3, (s-1)%%3)]) + sqrt(1-AlgoParams$rho^2) * epsilon
      
      #epsilon <- rnorm(length(unlist(latent_var)))
      #latent_var_prop <- relist(AlgoParams$rho * unlist(latent_var) + sqrt(1-AlgoParams$rho^2) * epsilon)
      latent_var_prop <- list()
      latent_var_prop$eps_event <- AlgoParams$rho * latent_var$eps_event + sqrt(1-AlgoParams$rho^2) * rnorm(length(latent_var$eps_event))
      latent_var_loc_prop <- AlgoParams$rho * unlist(latent_var$eps_local) + sqrt(1-AlgoParams$rho^2) * rnorm(length(unlist(latent_var$eps_local)))
      latent_var_prop$eps_local = mapply(function(vec, d) array(vec, dim = c(d[1], d[2], d[3])), 
                                         split(latent_var_loc_prop, rep(1:length(latent_var_local_dim), sapply(latent_var_local_dim, prod))), 
                                         latent_var_local_dim, 
                                         SIMPLIFY = FALSE)
      latent_var_prop$eps_local <- latent_var$eps_local #REMOVE
      latent_var_prop$eps_event <- latent_var$eps_event #REMOVE
      latent_var_prop$obs <- move_latent_obs(latent_var$obs, AlgoParams)
    }
    
    
    HP<- Model$HighLevelPriors(Omega_prop_phys %>% addTransfParams(), Model)
    if (HP> AlgoParams$ABC & Model$higherpriors){
      AlgoResults$accprob_store[s] = 0
      AlgoResults$Omega_sample[,s] = AlgoResults$Omega_sample[,s-1]
      AlgoResults$Omega_sample_phys[,s] = AlgoResults$Omega_sample_phys[,s-1]
      AlgoResults$loss[s] = AlgoResults$loss[s-1]
      #AlgoResults$u[,,,ifelse(s %% 3==0, 3, s%%3)] = AlgoResults$u[,,,ifelse((s-1) %% 3==0, 3, (s-1)%%3)]
      #AlgoResults$u_selected[,s] = AlgoResults$u_selected[,s-1]
      latent_var=latent_var
    } else {
      proposed = Omega_prop_phys %>% addTransfParams()
      #proposed$u = array(u_prop, dim=c(n_events, AlgoParams$m_CRPS*AlgoParams$Np, 3))
      # impact_sample <- SampleImpact(dir = dir,Model = Model,
      #                               proposed = proposed, 
      #                               AlgoParams = AlgoParams)
      # loss_prop <- CalcDist(impact_sample, AlgoParams)[1]
      loss_prop = CalcLL(dir, Model, proposed %>% addTransfParams(), AlgoParams, latent_var_prop)
      print(paste(s, loss_prop))
      
      #calculate the acceptance probability:
      #print(modifyAcc(Omega_prop, Omega_sample_s[n,], Model))
      min_loss <- min(loss_prop, AlgoResults$loss[s-1])
      #AlgoResults$accprob_store[s] <- min(1, exp(-AlgoParams$learning_rate * loss_prop + min_loss)/exp(-AlgoParams$learning_rate * AlgoResults$loss[s-1] + min_loss) * modifyAcc(Omega_prop, AlgoResults$Omega_sample[,s-1], Model, AlgoResults$lambda_store[s] * AlgoResults$Sigma_store[,,s]))
      AlgoResults$accprob_store[s] <- min(1, exp(loss_prop - min_loss)/exp(AlgoResults$loss[s-1] - min_loss) * modifyAcc(Omega_prop, AlgoResults$Omega_sample[,s-1], Model, AlgoResults$lambda_store[s] * AlgoResults$Sigma_store[,,s]))
      
      u <- runif(1)
      if(u < AlgoResults$accprob_store[s]){
        print('Accepted')
        AlgoResults$Omega_sample[,s] = Omega_prop
        AlgoResults$Omega_sample_phys[,s] = unlist(Omega_prop_phys)
        AlgoResults$loss[s] = loss_prop
        #AlgoResults$u[,,,ifelse(s %% 3==0, 3, s%%3)] = proposed$u
        #AlgoResults$u_selected[,s] =  c(proposed$u[1,1,1],proposed$u[2,2,2],proposed$u[3,3,3])
        latent_var= latent_var_prop
      }  else {
        AlgoResults$Omega_sample[,s] = AlgoResults$Omega_sample[,s-1]
        AlgoResults$Omega_sample_phys[,s] = AlgoResults$Omega_sample_phys[,s-1]
        AlgoResults$loss[s] = AlgoResults$loss[s-1]
        # AlgoResults$u[,,,ifelse(s %% 3==0, 3, s%%3)] = AlgoResults$u[,,,ifelse((s-1) %% 3==0, 3, (s-1)%%3)]
        # AlgoResults$u_selected[,s] = AlgoResults$u_selected[,s-1]
        latent_var=latent_var
      }
    }
    
    plot(AlgoResults$Omega_sample_phys[3,])
    
    if (iter_func(s)==iter_func(s-1)){
      AlgoResults$mu_store[,s] = (s - iter_func(s))/(s-iter_func(s)+1) *AlgoResults$mu_store[,s-1] + 1/(s-iter_func(s)+1) * AlgoResults$Omega_sample[,s]
      AlgoResults$Sigma_store[,,s] = 1/(s-iter_func(s)+AlgoParams$v_0 + n_x + 2) * (
        (s - iter_func(s) + AlgoParams$v_0+n_x+1)*AlgoResults$Sigma_store[,,s-1] + 
          AlgoResults$Omega_sample[,s] %*% t(AlgoResults$Omega_sample[,s]) +
          (s-iter_func(s)) * AlgoResults$mu_store[,s-1] %*% t(AlgoResults$mu_store[,s-1]) -
          (s-iter_func(s)+1) * AlgoResults$mu_store[,s] %*% t(AlgoResults$mu_store[,s]))
    } else {
      AlgoResults$mu_store[,s] = AlgoResults$mu_store[,s-1] + 1/(s-iter_func(s)+1) * (AlgoResults$Omega_sample[,s]-AlgoResults$Omega_sample[,iter_func(s)-1])
      AlgoResults$Sigma_store[,,s] = AlgoResults$Sigma_store[,,s-1] + 1/(s-iter_func(s) + AlgoParams$v_0 + n_x+2) * (
        AlgoResults$Omega_sample[,s] %*% t(AlgoResults$Omega_sample[,s]) - 
          AlgoResults$Omega_sample[,iter_func(s)-1] %*% t(AlgoResults$Omega_sample[,iter_func(s)-1]) +
          (s-iter_func(s)+1) * (AlgoResults$mu_store[,s-1] %*% t(AlgoResults$mu_store[,s-1]) -  AlgoResults$mu_store[,s] %*% t(AlgoResults$mu_store[,s]))
      )
    }
    
    AlgoResults$lambda_store[s] = max(AlgoParams$lambda_min, AlgoResults$lambda_store[s-1] * exp(delta / (AlgoResults$s_restart + s) * (AlgoResults$accprob_store[s]-AlgoParams$a)))
    if (abs(log(AlgoResults$lambda_store[s]) - log(AlgoResults$lambda_start))>log(3)){
      AlgoResults$lambda_start = AlgoResults$lambda_store[s]
      AlgoResults$s_restart <- 5/((AlgoParams$a)*(1-AlgoParams$a)) - s
    }
    
  }
}

#plot results:
latent_obs_store <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/latent_obs_2025-01-16_162110.15832')
latent_obs_final <- latent_obs_store[[1]][,,,dim(latent_obs_store[[1]])[4]]

folderin<-paste0(dir,AlgoParams$input_folder, "ODDobjects_Unweighted/")

ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) 
ufiles <- grep('^Train/' , ufiles, value = TRUE)
x <- file.info(paste0(folderin,ufiles))
ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
ODDy<-readODD(paste0(folderin,ufiles[3]))
plot(ODDy)
ODDy$mortality <- latent_obs_final[,1]
ODDy$displacement <- latent_obs_final[,2]  
ODDy$buildDam <- latent_obs_final[,3]
plot(ODDy)
ODD_df <- as.data.frame(ODDy, na.rm=F)

plot(ODD_df$hazMean1[ODDy@polygons[[1]]$indexes], ODD_df$mortality[ODDy@polygons[[1]]$indexes])
plot(ODD_df$hazMean1[ODDy@polygons[[1]]$indexes], ODD_df$Population[ODDy@polygons[[1]]$indexes])
points(ODD_df$hazMean1[ODDy@polygons[[1]]$indexes], ODD_df$Population[ODDy@polygons[[1]]$indexes], 
       col=ODD_df$mortality[ODDy@polygons[[1]]$indexes])

df_subset <- ODD_df[ODDy@polygons[[1]]$indexes, ]
df_subset <- df_subset[!is.na(df_subset$hazMean1) & !is.na(df_subset$Population) & !is.infinite(df_subset$hazMean1) & !is.infinite(df_subset$Population), ]
ggplot(df_subset, aes(x=hazMean1, y=Population, col=as.factor(mortality))) + geom_point()

# df_subset <- ODD_df[ODDy@polygons[[158]]$indexes, ]
# df_subset <- df_subset[!is.na(df_subset$hazMean1) & !is.na(df_subset$Population) & !is.infinite(df_subset$hazMean1) & !is.infinite(df_subset$Population), ]
# ggplot(as.data.frame(df_subset), aes(x=hazMean1, y=Population, col=displacement)) + geom_point()
       
indexes <- ODDy@polygons[[100]]$indexes
#indexes <- indexes[!(indexes %in% ODDy@polygons[[100]]$indexes)]
plot(ODD_df[indexes, 'hazMean1'], ODD_df[indexes, 'displacement'])
#points(ODD_df[ODDy@polygons[[100]]$indexes, 'hazMean1'], ODD_df[ODDy@polygons[[100]]$indexes, 'displacement'], col='red')

impact_with_poly_names <- merge(ODDy@impact, 
      data.frame(id=1:length(ODDy@polygons),
           polygon_name=unlist(lapply(ODDy@polygons, function(x) x$name))),
      by.x='polygon', by.y='id')
impact_with_poly_names[112:190,]


# vs_sample(log(c(1000, 1000)+10), log(matrix(rnorm(200, 1000, 90)+10, nrow=2)))
# vs_sample(log(c(100, 100)+10), log(matrix(rlnorm(200, 100, 10), nrow=2)))
# vs_sample(log(c(10, 10)+10), log(matrix(rlnorm(200, 10, 10)+10, nrow=2)))
# 
# vs_sample(log(c(1000,1000)+10), log(matrix(rnorm(1000, 1000, 100)+10, nrow=2)))/5
# vs_sample(log(c(100,100)+10), log(matrix(rnorm(1000, 100, 10)+10, nrow=2)))/5
# vs_sample(log(c(10,10)+10), log(matrix(rnorm(1000, 10, 2)+10, nrow=2)))/5
# vs_sample(log(c(1,1)+10), log(matrix(rnorm(1000, 1, 1)+10, nrow=2)))/5
# 
# es_sample(log(c(1000,1000)+10), log(matrix(rnorm(1000, 1000, 100)+10, nrow=2)))
# es_sample(log(c(100,100)+10), log(matrix(rnorm(1000, 100, 10)+10, nrow=2)))
# es_sample(log(c(10,10)+10), log(matrix(rnorm(1000, 10, 2)+10, nrow=2)))
# es_sample(log(c(1,1)+10), log(matrix(rnorm(1000, 1, 1)+10, nrow=2)))
# 
# plot(t(matrix(rnorm(1000, 100, 10), nrow=2)), xlim=c(60, 240), ylim=c(60, 240))
# points(t(matrix(rnorm(1000, 220, 5), nrow=2)), col='cyan')
# points(t(matrix(rnorm(1000, 200, 10), nrow=2)), col='blue')
# points(cbind(rnorm(500, 200, 10), rnorm(500, 100, 10)), col='green')
# points(rmvnorm(1000, c(150, 150), cbind(c(500, 450), c(450,500))), col='magenta')
# points(t(c(100,100)), col='red', pch=19)
# vs_sample(c(100,100), matrix(rnorm(1000, 100, 10), nrow=2))
# vs_sample(c(100,100), matrix(rnorm(1000, 200, 10), nrow=2))
# vs_sample(c(100,100),t(cbind(rnorm(500, 200, 10), rnorm(500, 100, 10))))
# vs_sample(c(100,100), t(rmvnorm(1000, c(150, 150), cbind(c(500, 450), c(450,500)))))
# vs_sample(c(100,100), matrix(rnorm(1000, 220, 5), nrow=2))
# 
# vs_sample(log(c(100,100)+10), log(matrix(rnorm(1000, 100, 10)+10, nrow=2)))/5 + es_sample(log(c(100,100)+10), log(matrix(rnorm(1000, 100, 10)+10, nrow=2)))
# vs_sample(log(c(100,100)+10), log(matrix(rnorm(1000, 200, 10)+10, nrow=2)))/5 + es_sample(log(c(100,100)+10), log(matrix(rnorm(1000, 200, 10)+10, nrow=2)))
# 
# vs_sample(log(c(exp(7.5), exp(7.5))+10), log(matrix(rlnorm(200, 7.5, 0.2)+10, nrow=2)))
# vs_sample(log(c(exp(5.5), exp(5.5))+10), log(matrix(rlnorm(200, 5.5, 0.2)+10, nrow=2)))
# vs_sample(log(c(exp(3.5), exp(3.5))+10), log(matrix(rlnorm(200, 3.5, 0.2)+10, nrow=2)))
# vs_sample(log(c(exp(0), exp(0))+10), log(matrix(rlnorm(200, 0, 0.5)+10, nrow=2)))
# 
# vs_sample(c(exp(7.5), exp(7.5)), matrix(rlnorm(200, 7.5, 0.2), nrow=2))
# vs_sample(c(exp(3.5), exp(3.5)), matrix(rlnorm(200, 3.5, 0.2), nrow=2))
# vs_sample(c(exp(0), exp(0)), matrix(rlnorm(200, 0, 0.5), nrow=2))
# 
# vs_sample(c(100, 100), matrix(rnorm(200, 100, 10)+10, nrow=2))
# vs_sample(c(10, 10),matrix(rnorm(200, 10, 10)+10, nrow=2))
# 
# obs <- c(100, 1000)
# weights <- log(obs) + 10
# weights <- weights %*% t(weights)
# ?vs_sample

#finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); start_time <- Sys.time()
  
#   
#   #Function to predict damage per gridpoint
#   CalcDam<-function(ij){
#     
#     # elapsed_time <- c()
#     # start_time <- Sys.time()
#     
#     #finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochasticregen = finish_time-start_time); start_time <- Sys.time()
#     
#     eps_local_long <- rmvnorm(length(hrange)*Method$Np, rep(0,3), sigma=covar_matrix_local)
#     eps_local_ij <- aperm(array(eps_local_long, dim=c(length(hrange), Method$Np, 3)), c(1,3,2))
#     
#     ##SLOWER:
#     # for (i in 1:Method$Np){
#     #   eps_local_ij[,,i] <- rmvnorm(length(hrange), rep(0,3), sigma=covar_matrix_local)
#     # }
#     # notnans_ij <- which(notnans==ij)
#     #eps_local_ij <- adrop(eps_local_transf[,,,notnans_ij, drop=F], drop=4)
#     
#     locallinp<- LP[ij,] # LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]]*LP$Plinp[ij]*LP$linp[[iso3c]] 
#     
#     # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochasticreshape = finish_time-start_time); start_time <- Sys.time()
#     
#     # Sample population per income distribution (Assumes 8 percentiles):
#     # Population is split evenly between the income quantiles, with remainders randomly allocated between
#     lPopS <- matrix(ODD_df$Population[ij] %/% 8, nrow=8, ncol=Method$Np) + rmultinom(Method$Np,ODD_df$Population[ij] %% 8,rep(1/8,8)) 
#     #lPopS <- SplitSamplePop(Pop=ODD@data$Population[ij],Method$Np) #matrix(round(ODD@data$Population[ij]/length(locallinp)), nrow=length(locallinp), ncol = Method$Np)
#     #lPopS <- matrix(round(ODD@data$Population[ij]/ 8), nrow=8, ncol=Method$Np)
#     
#     # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); start_time <- Sys.time()
#     
#     
#     lPopDisp <- array(0, dim=c(length(locallinp), Method$Np))
#     lPopMort <- array(0, dim=c(length(locallinp), Method$Np))
#     tPop <-array(0,c(3, Method$Np)) #row 1 = tDisp, #row 2 = tMort, #row 3 = tRem
#     tPop[3,]=colSums(lPopS)
#     
#     # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); 
#     
#     for(h_i in hrange_order){
#       start_time <- Sys.time()
#       h <- hrange[h_i]
#       
#       if(is.na(ODD_df[ij,h])) next
#       
#       nonzero_pop <- which(lPopS != 0, arr.ind=T)
#       if (length(nonzero_pop)==0) next
#       
#       
#       # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); start_time <- Sys.time()
#       
#       ## Sample hazard Intensity 
#       ## Doesn't work very well as don't have the covariance of the errors, so must instead assume they are independent
#       ## So instead just use the mean value
#       # I_ij<-rnorm(n = Method$Np,
#       #             mean = ODD@data[ij,paste0("hazMean",h)],
#       #             sd = ODD@data[ij,paste0("hazSD",h)]/10)
#       I_ij<-ODD_df[ij,h]
#       Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=NROW(nonzero_pop)),Omega) + locallinp[nonzero_pop[,1]] + event_lp[h_i], error=function(e) NA) #+ rep(eps_local[h_i,], each=8), error=function(e) NA)
#       
#       # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, damcalc = finish_time-start_time); start_time <- Sys.time()
#       
#       D_MortDisp <- D_MortDisp_calc(Damage, Omega, eps_event[1:2, nonzero_pop[,2], drop=F] + adrop(eps_local_ij[h_i,1:2,nonzero_pop[,2], drop=F], drop = 1)) #First row of D_MortDisp is D_Mort, second row is D_Disp
#       
#       # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, mortdisp = finish_time-start_time); start_time <- Sys.time()
#       
#       D_Rem <- pmax(0, 1 - D_MortDisp[2,] - D_MortDisp[1,]) #probability of neither death nor displacement. Use pmax to avoid errors caused by numerical accuracy.
#       
#       # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, d_rem = finish_time-start_time); start_time <- Sys.time()
#       
#       Dam <- Fbdam(lPopS[nonzero_pop], D_MortDisp[2,], D_MortDisp[1,], D_Rem)
#       
#       # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, fbdam = finish_time-start_time); start_time <- Sys.time()
#       
#       lPopS[nonzero_pop] <- Dam[3,]
#       lPopDisp[nonzero_pop] <- lPopDisp[nonzero_pop] + Dam[1,]
#       lPopMort[nonzero_pop] <- lPopMort[nonzero_pop] + Dam[2,]
#       tPop[3,] <- colSums(lPopS)
#       
#       # finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, sum_pops = finish_time-start_time); 
#       
#       # # This is a bit clearer than the above but slower:
#       # Accumulate the number of people displaced/deceased, but don't accumulate the remaining population
#       # tPop[3,ind]<-0
#       # for (s in 1:length(SincN)){ #Separate into income distributions (as each have 10% of population, order doesn't matter)
#       #   if(all(lPopS[s,]==0)) next
#       #   # Predict damage at coordinate {i,j} (vector with MC particles)
#       #   Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=sum(ind)),Omega) + locallinp[s] + eps_event[h_i,ind], error=function(e) NA)
#       #   if(any(is.na(Damage))) print(ij)
#       # 
#       #   #LOOSEEND: Include [ind] here
#       #   D_MortDisp <- D_MortDisp_calc(Damage, Omega) #First row of D_MortDisp is D_Mort, second row is D_Disp
#       #   D_Rem <- pmax(0, 1 - D_MortDisp[2,] - D_MortDisp[1,]) #probability of neither death nor displacement. Use pmax to avoid errors caused by numerical accuracy.
#       # 
#       #   tPop[,ind]<-tPop[,ind] + Fbdam(lPopS[s,ind],D_MortDisp[2,], D_MortDisp[1,], D_Rem)
#       # }
#     }
#     
#     # start_time <- Sys.time()
#     
#     tPop[1,] <- colSums(lPopDisp)
#     tPop[2,] <- colSums(lPopMort)
#     
#     #ensure the total displaced, deceased or remaining does not exceed total population
#     tPop[tPop>ODD_df$Population[ij]] <- floor(ODD_df$Population[ij])
#     
#     # if (length(elapsed_time)==9) {return(elapsed_time)
#     # } else {(return(rep(0,9)))}
#     
#     #if no building destruction data:
#     if(!BD_data_present) return(list(samples=rbind(tPop[1:2,, drop=FALSE], rep(NA, Method$Np))))#, 
#     
#     locallinp_buildings <- LP_buildings[ij]
#     
#     nUnaff = rep(ODD_df$nBuildings[ij], Method$Np)
#     nDam = rep(0, Method$Np)
#     
#     #finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time); start_time <- Sys.time()
#     
#     for (h_i in hrange_order){
#       h <- hrange[h_i]
#       if(is.na(ODD_df[ij,h])) next
#       if(all(nUnaff==0)) break #if no remaining buildings, skip modelling
#       
#       I_ij<-ODD_df[ij,h]
#       Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Method$Np),Omega) + locallinp_buildings + event_lp[h_i], error=function(e) NA) #+ eps_local[h_i,], error=function(e) NA) #calculate unscaled damage (excluding GDP)
#       
#       D_Dam <- D_Dam_calc(Damage, Omega, eps_event[3,] + eps_local_ij[h_i,3,]) 
#       
#       # Accumulate the number of buildings damaged/destroyed, but not the number of buildings remaining
#       nDam_new <- rbinom(Method$Np, nUnaff, D_Dam)
#       nUnaff <- nUnaff - nDam_new
#       nDam <- nDam + nDam_new
#       
#     }
#     
#     #finish_time <- Sys.time(); elapsed_time <- c(elapsed_time, stochastic_sample = finish_time-start_time);
#     
#     #return(elapsed_time)
#     
#     
#     return(list(samples = rbind(tPop[1:2,,drop=FALSE], nDam[1:Method$Np])))
#   }
#   
#   Dam<-array(0,c(nrow(ODD_df),Method$Np,3)) # Dam[,,1] = Displacement, Dam[,,2] = Mortality, Dam[,,3] = Buildings Damaged
#   Dam_means<-array(0,c(nrow(ODD_df),3))
#   
#   if(Method$NestedCores>1) { 
#     CalcDam_out <- mclapply(X = notnans,FUN = CalcDam,mc.cores = Method$NestedCores)
#   } else  {CalcDam_out <- lapply(X = notnans,FUN = CalcDam)}
#   
#   Dam[notnans,,]<-aperm(simplify2array(lapply(CalcDam_out, function(x) x$samples)), perm=c(3,2,1))
#   
#   if (output=='SampledFull'){
#     return(Dam)
#   } else if (output == 'SampledTotal'){
#     SampledTot <- colSums(Dam)
#     df_SampledTot <- list()
#     impact_types <- unique(ODD@impact$impact)
#     
#     for (impact_type in impact_types){
#       polygon_names <- unlist(lapply(ODD@polygons[ODD@impact$polygon], function(x) x$name))
#       if (any(tolower(polygon_names[which(ODD@impact$impact==impact_type)]) %in% c('tot', 'total'))){
#         nonmatch <- which(!tolower(polygon_names[which(ODD@impact$impact==impact_type)]) %in% c('tot', 'total'))
#         if (length(nonmatch)>0){
#           ODD@impact <- ODD@impact[-which(ODD@impact$impact==impact_type)[nonmatch],] # in the case of total and subnational data, remove the subnational
#         }
#       }
#     }
#     
#     #Many ODD objects don't contain 'total' impact values (to avoid double counting), so we need to obtain these. 
#     #We combine the polygons with observations so long as the proportion of overlapping pixels is less than 10%
#     #and the total coverage of the exposed area is greater than 90%. Then sum the observations across these polygons.
#     observed_total=rep(NA, length(impact_types))
#     exposed_haz <- which(apply(ODD_df[,grep('hazMean', names(ODD_df)), drop=F], 1, function(row) any(!is.na(row))) & !is.na(ODD_df$ISO3C))
#     get_overlap_coverage <- function(impact_type){
#       indexes_list <- lapply(ODD@polygons[ODD@impact$polygon[which(ODD@impact$impact==impact_type)]], function(x) x$indexes)
#       universal_set <- Reduce(union, indexes_list)
#       overlap <- Reduce(intersect, indexes_list)
#       overlap <- intersect(overlap, exposed_haz)
#       prop_overlap <- ifelse(length(indexes_list)>1,length(overlap)/length(universal_set),0)
#       prop_coverage <- length(intersect(unique(universal_set), exposed_haz))/length(exposed_haz)#sum(!is.na(ODD_df$ISO3C))
#       return(c(prop_overlap, prop_coverage))
#     }
#     
#     overlap_coverage <- sapply( impact_types,get_overlap_coverage)
#     
#     for (i in 1:length(impact_types)){
#       observed_total[i] = sum(ODD@impact$observed[which(ODD@impact$impact==impact_types[i])])
#     }
#     
#     
#     for (i in 1:NROW(SampledTot)){
#       df_SampledTot[[i]] <- data.frame(event_id = ODD@impact$event_id[1],
#                                        polygon=0,
#                                        iso3 = paste0(unique(ODD@impact$iso3), collapse=' '),
#                                        impact=impact_types, 
#                                        sampled=SampledTot[i,match(impact_types, c('displacement', 'mortality','buildDam'))], 
#                                        observed=observed_total,
#                                        qualifier='total',
#                                        overlap = overlap_coverage[1,],
#                                        coverage = overlap_coverage[2,],
#                                        inferred=F)
#       df_SampledTot[[i]] %<>% filter(overlap < 0.1 & coverage > 0.9)
#     }
#     
#     
#     return(df_SampledTot)
#   }
#   
#   
#   #finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, dam_sample = finish_time-start_time); start_time <- Sys.time()
#   
#   funcy<-function(i, impact_weights=Method$impact_weights, kernel=Method$kernel, cap=Method$cap, polygons_indexes=ODD@polygons) {
#     
#     # funcy() aggregates the simulate impact (for each impact type) across each polygon, when
#     # there is a matching observation for that impact type and polygon
#     
#     tmp<-data.frame(iso3=ODD_df$ISO3C, displacement=Dam[,i,1], mortality=Dam[,i,2], buildDam=Dam[,i,3])#, 
#     #mort_mean=Dam_means[,1], disp_mean=Dam_means[,2], buildDam_mean=Dam_means[,3])
#     impact_sampled<-data.frame(polygon = numeric(), impact = character(), sampled = numeric(), mean=numeric())
#     
#     for (polygon_id in unique(ODD@impact$polygon)){
#       polygon_impacts <- ODD@impact$impact[which(ODD@impact$polygon==polygon_id)]
#       for (impact in polygon_impacts){
#         if (impact == 'mortality'){
#           impact_sampled <- rbind(impact_sampled, data.frame(polygon=polygon_id,
#                                                              impact=impact,
#                                                              sampled=floor(sum(tmp[polygons_indexes[[polygon_id]]$indexes,impact]  * polygons_indexes[[polygon_id]]$weights, na.rm=T))))
#         } else if (impact=='displacement'){
#           impact_sampled <- rbind(impact_sampled, data.frame(polygon=polygon_id,
#                                                              impact=impact,
#                                                              sampled=floor(sum(tmp[polygons_indexes[[polygon_id]]$indexes,impact]  * polygons_indexes[[polygon_id]]$weights, na.rm=T))))
#         } else {
#           impact_sampled <- rbind(impact_sampled, data.frame(polygon=polygon_id,
#                                                              impact=impact,
#                                                              sampled=floor(sum(tmp[polygons_indexes[[polygon_id]]$indexes,impact]  * polygons_indexes[[polygon_id]]$weights, na.rm=T))))
#         }
#       }
#     }
#     impact_obs_sampled <- arrange(merge(impact_sampled, ODD@impact, by=c("polygon", "impact")),desc(observed)) 
#     return(impact_obs_sampled)
#   }
#   
#   if (output == 'ODDwithSampled'){ #usually used for generating simulated data
#     ODD[['Disp']]<- Dam[,1,1]  
#     ODD[['Mort']]<- Dam[,1,2]
#     ODD[['BuildDam']] <- Dam[,1,3]  
#     ODD@predictDisp<-funcy(1) 
#     return(ODD)
#   }
#   
#   for (i in 1:length(ODD@polygons)){ #weights can be used to downweight pixels that are only partly in a polygon
#     if(is.null(ODD@polygons[[i]]$weights)){
#       ODD@polygons[[i]]$weights <- rep(1, length=length(ODD@polygons[[i]]$indexes))
#     }
#   }
#   
#   ## IF TIMING: --------------------------------
#   # dummy <- lapply(1:Method$Np, funcy, LLout=F)
#   # finish_time <-  Sys.time(); elapsed_time <- c(elapsed_time, agg_dam = finish_time-start_time); start_time <- Sys.time()
#   # return(elapsed_time)
#   # --------------------------------------------
#   
#   if(is.null(ODD_df$nBuildings)){
#     ODD@impact <- ODD@impact[!ODD@impact$impact %in% c('buildDam', 'buildDest', 'buildDamDest'),]
#   }
#   
#   if(output == 'SampledAgg'){
#     return(lapply(1:Method$Np, funcy))
#   }
#   
#   stop('Output type selected as "output" not recognised, please select either SampledAgg, SampledTotal,
#        SampledFull, or ODDwithSampled')
#   
#   return(ODD)
#   
# })