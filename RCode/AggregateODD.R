

#----------------------------------------------------------------------------------------------------------------------------------
# ---------------- AGGREGATE BY REGIONS WITH SIMILAR VULNERABILITY / HAZARD INTENSITY / POPULATION DENSITY ------------------------
#----------------------------------------------------------------------------------------------------------------------------------

# Out of date since move from spdf to raster:

# increaseAggregation_all <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12'){
#   ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
#   ODD_folderout<-paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_Aug31_Agg', '/ODDobjects/')
#   ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
#   for (file in ufiles[1:length(ufiles)]){
#     event_id <- as.numeric(strsplit(file, "_")[[1]][2])
#     ODDy <- readRDS(paste0(ODD_folderin, file))
#     saveRDS(ODDy, paste0(dir, folder_in, '/ODDobjects/', file))
#     ODDyAgg <- tryCatch(increaseAggregation(removeWeights(ODDy)),error=function(e) NULL)
#     if(is.null(ODDyAgg)){
#       print(event_id)
#       next
#     }
#     saveRDS(ODDyAgg, paste0(ODD_folderout, file))
#   }
# }
# 
# 
# removeWeights <- function(ODD){
#   pixels_of_interest <- which(!is.na(ODD$ISO3C))
#   for (i in pixels_of_interest){
#     weight_sum = c(0, 0, 0)
#     polygon_matches = list(c(), c(), c())
#     for (p in 1:length(ODD@polygons)){
#       poly <- ODD@polygons[[p]]
#       match_i <- which(poly$indexes == i)
#       if (length(match_i) > 0){
#         gadm_level <- str_count(poly$name, ',')
#         polygon_matches[[gadm_level+1]] %<>% rbind(c(p, match_i))
#         weight_sum[gadm_level+1] <- weight_sum[gadm_level+1] + poly$weights[match_i]
#       }
#     }
#     for (g in 2:3){
#       if (is.null(polygon_matches[[g]])) next
#       if (NROW(polygon_matches[[g]])==1) next
#       if (round(weight_sum[g],3)>1) stop()
#       weights <- rep(0, NROW(polygon_matches[[g]]))
#       for (j in 1:NROW(polygon_matches[[g]])){
#         weights[j] <- ODD@polygons[[polygon_matches[[g]][j,1]]]$weights[polygon_matches[[g]][j,2]]
#       }
#       j_max_weight <- which.max(weights)
#       for (j in 1:NROW(polygon_matches[[g]])){
#         if (j == j_max_weight){
#           ODD@polygons[[polygon_matches[[g]][j,1]]]$weights[polygon_matches[[g]][j,2]] = 1
#         } else {
#           ODD@polygons[[polygon_matches[[g]][j,1]]]$indexes <- ODD@polygons[[polygon_matches[[g]][j,1]]]$indexes[-polygon_matches[[g]][j,2]]
#           ODD@polygons[[polygon_matches[[g]][j,1]]]$weights <- ODD@polygons[[polygon_matches[[g]][j,1]]]$weights[-polygon_matches[[g]][j,2]]
#         }
#       }
#     }
#   }
#   return(ODD)
# }
# 
# increaseAggregation <- function(ODD){
#   ODD$PDens <- exp(round(log(ODD$PDens)))
#   ODD$Vs30 <- round(ODD$Vs30, -2)
#   ODD$EQFreq <- exp(round(log(ODD$EQFreq+0.1)/0.2, 0)*0.2)-0.1
#   hrange<-grep("hazMean",names(ODD),value = T)
#   ODD@data[,hrange] <- log(round(1.3^(1.5*ODD@data[,hrange])),base=1.3)/1.5  #round(ODD@data[,hrange]/0.2)*0.2 #
#   
#   polyMatch <- list()
#   polyMatch[1:NROW(ODD@data)] <- ''
#   for (i in 1:length(ODD@polygons)){
#     polyMatch[ODD@polygons[[i]]$indexes] <- paste0(polyMatch[ODD@polygons[[i]]$indexes], i,',')
#   }
#   ODD@data$polyMatch <- unlist(polyMatch)
#   ODD@data$index_original <- 1:NROW(ODD@data)
#   grouped_by_covar <- ODD@data %>% group_by(across(all_of(c(hrange, 'ISO3C',  'PDens', 'Vs30', 'EQFreq', 'AveSchYrs', 'LifeExp', 'GNIc', 'SHDI', 'polyMatch'))))
#   
#   if (!is.null(ODD@data$nBuildings)){
#     summarised <- grouped_by_covar %>% summarize(Population=sum(Population, na.rm=T), nBuildings=sum(nBuildings, na.rm=T), disAgg_indexes=list(index_original))
#   } else {
#     summarised <- grouped_by_covar %>% summarize(Population=sum(Population, na.rm=T), disAgg_indexes=list(index_original))
#   }
#   
#   ODDagg <- ODD
#   ODDagg@data <- as.data.frame(summarised)
#   for (i in 1:length(ODDagg@polygons)){
#     ODDagg@polygons[[i]]$indexes <- c()
#   }
#   
#   for (i in 1:NROW(ODDagg@data)){
#     match_polys <- as.numeric(unlist(str_extract_all(ODDagg@data$polyMatch[i], "\\d+")))
#     for (j in match_polys){
#       ODDagg@polygons[[j]]$indexes <- c(ODDagg@polygons[[j]]$indexes, i)
#     }
#   }
#   for (i in 1:length(ODDagg@polygons)){
#     ODDagg@polygons[[i]]$weights <- rep(1, length(ODDagg@polygons[[i]]$indexes))
#   }
#   return(ODDagg)
# }
# 
# increaseAggregation_all <- function(folder_in='IIDIPUS_Input_NonFinal/IIDIPUS_Input_July12'){
#   ODD_folderin<-paste0(dir, folder_in, '/ODDobjects/')
#   ODD_folderout<-paste0(dir, 'IIDIPUS_Input_NonFinal/IIDIPUS_Input_Aug31_Agg', '/ODDobjects/')
#   ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
#   for (file in ufiles[1:length(ufiles)]){
#     event_id <- as.numeric(strsplit(file, "_")[[1]][2])
#     ODDy <- readRDS(paste0(ODD_folderin, file))
#     saveRDS(ODDy, paste0(dir, folder_in, '/ODDobjects/', file))
#     ODDyAgg <- tryCatch(increaseAggregation(removeWeights(ODDy)),error=function(e) NULL)
#     if(is.null(ODDyAgg)){
#       print(event_id)
#       next
#     }
#     saveRDS(ODDyAgg, paste0(ODD_folderout, file))
#   }
# }
# 
# plot_aggregated_regions <- function(ODDyAgg, ODD){
#   plot(ODD@coords)
#   available_colors <- colors()
#   for (i in 1:NROW(ODDyAgg)){
#     points(ODD@coords[ODDyAgg$disAgg_indexes[[i]],], col=sample(available_colors, 1), pch=19)
#   }
# }

#----------------------------------------------------------------------------------------------------------------------------------
# -----------------------------Aggregate by combining cells with small population -------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

# Out of date since move from spdf to raster:

# aggregateODDbyPop <- function(ODD){
#   coords_df <- data.frame(ODD@coords)
#   coords_df$id <- 1:NROW(ODD@coords)
#   reshape(coords_df, timevar = "Latitude", idvar = "Longitude", direction = "wide")
#   
#   coords_wide <- as.matrix(coords_df %>% pivot_wider(names_from = Longitude, values_from = id, values_fill = NA))
#   rownames(coords_wide) <- coords_wide[,1]
#   coords_wide <- coords_wide[,-1]
#   1:NCOL(coords_wide)
#   col_chunks <- split(1:NCOL(coords_wide), ceiling(seq_along(1:NCOL(coords_wide))/4))
#   row_chunks <- split(1:NROW(coords_wide), ceiling(seq_along(1:NROW(coords_wide))/4))
#   agged <- c()
#   pop <- c()
#   ODDyAgg_groupings <- list()
#   for (col_chunk in col_chunks){
#     for (row_chunk in row_chunks){
#       indexes <- c(coords_wide[row_chunk, col_chunk])
#       pop_sum <- sum(ODD$Population[indexes], na.rm=T)
#       if (pop_sum == 0) next
#       if (pop_sum < 500){
#         ODDyAgg_groupings[[length(ODDyAgg_groupings)+1]] <- indexes
#         agged <- c(agged,2)
#         pop <- c(pop, pop_sum)
#       }
#       else {
#         for (row_chunk_half in split(row_chunk,ceiling(seq_along(row_chunk)/2))){
#           for (col_chunk_half in split(col_chunk,ceiling(seq_along(col_chunk)/2))){
#             indexes <- c(coords_wide[row_chunk_half, col_chunk_half])
#             pop_sum <- sum(ODD$Population[indexes], na.rm=T)
#             if (pop_sum == 0) next
#             if (pop_sum < 500){
#               ODDyAgg_groupings[[length(ODDyAgg_groupings)+1]] <- indexes
#               agged <- c(agged,1)
#               pop <- c(pop, pop_sum)
#             } else {
#               ODDyAgg_groupings <- append(ODDyAgg_groupings, indexes) #append each individually
#               agged <- c(agged, 0)
#               pop <- c(pop, ODD$Population[indexes])
#             }
#           }
#         }
#       }
#     }
#   }
#   length(ODDyAgg_groupings)/sum(!is.na(ODD$Population) &  ODD$Population>0)
# }
# Doesn't actually seem to reduce computation as much as would be expected
# e.g. Nepal 31 only reduces number of pixels by about half

#----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------Aggregate by doubling/tripling/... size of every pixel -----------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

# aggregateISO3C <- function(cell_values, na.rm=F){
#   if (all(is.na(cell_values))) return(NA)
#   cell_values <- cell_values[!is.na(cell_values)] 
#   return(as.numeric(names(sort(-table(cell_values)))[1]))
# }

# aggregateIndex <- function(cell_values, na.rm=T){
#   return(as.factor(paste(cell_values, collapse=' '))) #factor is annoying but aggregate function won't allow characters or arrays
# }

aggregateODDbyX <- function(ODD, aggFactor){
  
  #ODDAgg <- new('ODD')
  ODDAgg <- aggregate(ODD[['Population']], aggFactor, fun='sum')
  # ODDAggPop <- #brick(stack(ODDAggPop, aggregate(ODD[['ISO3C']], aggFactor, fun='aggregateISO3C')))
  # ODDAgg@file <- ODDAggPop@file
  # ODDAgg@data <- ODDAggPop@data
  # ODDAgg@legend <- ODDAggPop@legend
  # ODDAgg@title <- ODDAggPop@title
  # ODDAgg@extent <-ODDAggPop@extent
  # ODDAgg@rotated <- ODDAggPop@rotated
  # ODDAgg@rotation <- ODDAggPop@rotation
  # ODDAgg@ncols <- ODDAggPop@ncols
  # ODDAgg@nrows <- ODDAggPop@nrows
  # ODDAgg@crs <- ODDAggPop@crs
  # ODDAgg@srs <- ODDAggPop@srs
  # ODDAgg@history <- ODDAggPop@history
  # ODDAgg@z <- ODDAggPop@z
  
  ODDAgg@impact = ODD@impact
  ODDAgg@polygons = ODD@polygons
  for (i in 1:length(ODDAgg@polygons)){
    ODDAgg@polygons[[i]]$indexes <- c()
    ODDAgg@polygons[[i]]$weights <- c()
  }
  
  sum_vars <- c('PDens')
  if ('nBuildings' %in% names(ODD)) sum_vars <- c(sum_vars, 'nBuildings')
  
  for (var_name in sum_vars){
    ODDAgg[[var_name]] <- aggregate(ODD[[var_name]], aggFactor, fun='sum')
  }
  for (var_name in names(ODD)){
    print(var_name)
    if (var_name %in% sum_vars | var_name == 'Population' | var_name=='ISO3C') next
    ODDAgg[[var_name]] <- aggregate(ODD[[var_name]], aggFactor, fun='mean')
  }
  
  ODDAgg[['ISO3C']] <- aggregate(ODD[['ISO3C']], aggFactor, fun='modal', na.rm=T)
  
  #list with length = number of new cells, which each list element containing a vector of the corresponding old cells
  pixels_as_poly <- as.polygons(ODDAgg, aggregate=F, na.rm=F)
  coords_ODD <- vect(xyFromCell(ODD, 1:ncell(ODD)))
  intersect_oldnew <- relate(coords_ODD, pixels_as_poly, 'intersects')
  new_old_cell_ref <- apply(intersect_oldnew, 2, function(x) which(x)) 
  
  # index_agg <- aggregate(ODD[['index']],aggFactor, fun='aggregateIndex')
  # new_old_cell_ref <- list()
  # for (new_cell in levels(index_agg)[[1]]$ID){
  #   old_cells <- strsplit(levels(index_agg)[[1]]$VALUE[new_cell], ' ')[[1]]
  #   new_old_cell_ref[[new_cell]] = as.numeric(old_cells[old_cells != 'NaN'])
  # }
  
  pop_values <- values(ODD[['Population']])
  ncells = length(new_old_cell_ref)
  for (new_cell in 1:ncells){
    pop_weights <- pop_values[new_old_cell_ref[[new_cell]]]/sum(pop_values[new_old_cell_ref[[new_cell]]], na.rm=T)
    pop_weights[is.na(pop_weights)] = 0
    
    match_mat_weights <- do.call(rbind,lapply(ODD@polygons, function(x) x$weights[match(new_old_cell_ref[[new_cell]], x$indexes)]))
    match_mat_weights[which(is.na(match_mat_weights))] = 0
    
    if (NROW(match_mat_weights)>1){
      #weight of new pixel for polygon = sum over (proportion of new pixel's population that lies in that old pixel times
      #                                              old pixel's weight for that polygon)
      matched_weights <- rowSums(sweep(match_mat_weights, 2, pop_weights,'*'))
    } else {
      matched_weights <- sum(match_mat_weights * pop_weights)
    }
    
    for (p in which(matched_weights>0)){
      ODDAgg@polygons[[p]]$indexes = c(ODDAgg@polygons[[p]]$indexes, new_cell)
      ODDAgg@polygons[[p]]$weights = c(ODDAgg@polygons[[p]]$weights,matched_weights[p])
    }
  }

  return(ODDAgg)
}


increaseAggregation_all <- function(folder_in='IIDIPUS_Input'){
  ODD_folderin<-paste0(dir, folder_in, '/ODDobjects_RealFull/')
  ODD_folderout<-paste0(dir, folder_in, '/ODDobjects_RealFullAgg5/')
  ufiles<-list.files(path=ODD_folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  for (file in ufiles[1:length(ufiles)]){
    event_id <- as.numeric(strsplit(file, "_")[[1]][2])
    print(event_id)
    ODDy <- readRDS(paste0(ODD_folderin, file))
    ODDyAgg <- tryCatch(aggregateODDbyX(ODDy, 5),error=function(e) NULL)
    if(is.null(ODDyAgg)){
      print(paste('FAIL', event_id))
      next
    }
    saveRDS(ODDyAgg, paste0(ODD_folderout, file))
  }
}


#----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------Compare aggregated vs disaggregated impact ----------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
# AlgoParams$Np <- 500
# obs_index <- 1
# 
# #ODDyAgg <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects/Train/EQ20191215PHL_135')
# event <- 'EQ20180907ECU_95' #'EQ20191215PHL_135'#
# ODDyAgg <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects_RealAgg5/Train/', event))
# AggImpact <- DispX(ODDyAgg, Omega %>% addTransfParams(), Model$center, AlgoParams, output='SampledAgg')
# impAgg_sampled <- matrix(unlist(lapply(AggImpact, function(x) x$sampled[c(1,2)])), ncol=2, byrow=T)
# 
# #ODDy <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects_RealFull/Train/EQ20180216MEX_80')
# ODDy <- readRDS(paste0('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input/ODDobjects_RealFull/Train/', event))
# DisaggImpact <- DispX(ODDy, Omega %>% addTransfParams(), Model$center, AlgoParams, output='SampledAgg')
# impDisagg_sampled <- matrix(unlist(lapply(DisaggImpact, function(x) x$sampled[c(1,2)])), ncol=2, byrow=T)
# 
# qqplot(impAgg_sampled[,2], impDisagg_sampled[,2])
# abline(0,1)
# plot(log(impAgg_sampled+10))
# points(impDisagg_sampled, col='blue')
# 
# ggplot()  + 
#   #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000), labels = label_comma())  + 
#   #scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000), labels = label_comma()) +
#   geom_histogram(data=data.frame(sampled=log(impAgg_sampled[,1]+10)), aes(x=sampled,y=after_stat(count)), alpha=0.4, col="blue", lwd=0.2, fill='blue') +
#   geom_histogram(data=data.frame(sampled=log(impDisagg_sampled[,1]+10)), aes(x=sampled,y=after_stat(count)), alpha=0.4, col='red', lwd=0.2, fill='red') +
#   theme_bw() + ylab('Count') + xlab('log(Sampled Impact+10)')
