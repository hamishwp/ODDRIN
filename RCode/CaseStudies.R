
source('RCode/AutoQuake.R')
source('RCode/AutoQuake_functions.R')
folder_write='IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/'
library(ggridges)
library(GGally, include.only = c("ggpairs"))
#plotting functions

plot_predictive_dist_subnat <- function(impact_mort){
  #plots predictive distributions over subnational regions
  
  impact_mort$polygon_name = sub("(^([^,]*),[^,]*),.*", "\\1", impact_mort$polygon_name)
  long_df <- impact_mort %>%
    pivot_longer(cols = starts_with("sampled."),
                 names_to = "sample",
                 values_to = "value") %>%
    mutate(sample = as.integer(str_remove(sample, "sampled\\.")))  # Optional: just cleans sample id
  
  long_df$value %<>% as.numeric()
  
  long_df$polygon_name <- factor(long_df$polygon_name, levels = rev(impact_mort$polygon_name))
  impact_mort$polygon_name <- factor(impact_mort$polygon_name, levels = rev(impact_mort$polygon_name))
  
  # Plot using ggplot
  #more options: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
  ggplot(long_df, aes(x = value, y = polygon_name)) +
    geom_density_ridges(
      scale = 0.9,
      alpha=0.6,
      fill = '#440154',
      color = 'black',
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 2, point_alpha = 0.8, alpha = 0.7,
    ) +
    geom_point(
      data = impact_mort,
      aes(x = observed, y = polygon_name),
      color = "red",
      size = 3,
      shape = 19,  # Solid diamond, or use shape = 16 for filled circle
      inherit.aes = FALSE
    ) +
    labs(
      x = "Mortality",
      y = "Municipality"#,
      #title = "Posterior predictive distribution for mortality by municipality"
    ) +
    theme_minimal() + scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.5, base = 100),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      labels = scales::comma_format(),
      minor_breaks = NULL
    ) + theme(
      axis.title = element_text(family = "Liberation Serif", size=12),
      axis.text = element_text(family = "Liberation Serif", size=12)
    )
}

plot_predictive_dist_subnat <- function(impact_mort, n_cols = 1) {
  # plots predictive distributions over subnational regions
  # n_cols = 1 (default) for a single column; set to 2 for two columns, etc.
  
  # Keep only the first two comma-separated parts of the polygon name
  impact_mort$polygon_name <- sub("(^([^,]*),[^,]*),.*", "\\1", impact_mort$polygon_name)
  
  # Long format for the posterior samples
  long_df <- impact_mort |>
    tidyr::pivot_longer(
      cols = starts_with("sampled."),
      names_to = "sample",
      values_to = "value"
    ) |>
    dplyr::mutate(sample = as.integer(stringr::str_remove(sample, "sampled\\.")))
  
  long_df$value <- as.numeric(long_df$value)
  
  # Consistent ordering (reverse so the first appears at top)
  impact_mort$polygon_name <- factor(impact_mort$polygon_name,
                                     levels = rev(impact_mort$polygon_name))
  long_df$polygon_name     <- factor(long_df$polygon_name,
                                     levels = levels(impact_mort$polygon_name))
  
  # If requested, split the y categories across n_cols facets
  if (n_cols > 1) {
    levs <- levels(impact_mort$polygon_name)
    n <- length(levs)
    # Split by index into roughly equal groups
    idx <- seq_len(n)
    # create breaks for n_cols chunks
    brks <- round(seq(0, n, length.out = n_cols + 1))
    brks[1] <- 0
    brks[length(brks)] <- n
    grp <- cut(idx, breaks = brks, include.lowest = TRUE,
               labels = paste0("Panel ", seq_len(n_cols)))
    map_df <- data.frame(polygon_name = levs, facet_col = grp)
    
    impact_mort <- merge(impact_mort, map_df, by = "polygon_name", all.x = TRUE)
    long_df     <- merge(long_df,     map_df, by = "polygon_name", all.x = TRUE)
  }
  
  
  p <- ggplot(long_df, aes(x = value, y = polygon_name)) +
    ggridges::geom_density_ridges(
      scale = 1,
      fill = "#440154",
      color = "black",
      jittered_points = TRUE,
      position = ggridges::position_points_jitter(width = 0.05, height = 0),
      point_shape = "|", point_size = 2, point_alpha = 0.8,
      alpha = 0.6
    ) +
    # Observed point
    geom_point(
      data = impact_mort,
      aes(x = observed, y = polygon_name),
      color = "red",
      size = 3,
      shape = 19,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Mortality",
      y = "Municipality"
    ) +
    theme_minimal() +
    scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.5, base = 100),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      labels =  expression(phantom(1) * 0^phantom(0),phantom(1) * 1^phantom(0), phantom(1) * 10^phantom(1),10^2,10^3, 10^4),
      minor_breaks = NULL
    ) +
    theme(
      axis.title = element_text(family = "Liberation Serif", size = 12),
      axis.text  = element_text(family = "Liberation Serif", size = 12),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.text = element_blank(), strip.background = element_blank()
    )
  
  if (n_cols > 1) {
    p <- p + facet_wrap(~ facet_col, ncol = n_cols, scales = "free_y")
  }
  
  p
}

plot_GADM_impact_polygons <- function(ODDy, polygons_list, gadm_level = 2, impact_type = 'mortality', summary = 'mean', plot_bbox = NULL){
  #GADM_level 0 = national
  #GADM_level 1 = admin level 1
  #GADM_level 2 = admin level 2
  #summary = mean, median, q5, q10, q20, q80, q90, q95
  if (is.null(plot_bbox)){plot_bbox = ext(ODDy)}
  
  ODDy$polygon_id = 0
  for (i in 1:length(ODDy@polygons)){
    if (str_count(ODDy@polygons[[i]]$name, ',') != gadm_level | ODDy@polygons[[i]]$name == 'TOTAL') next
    ODDy$polygon_id[ODDy@polygons[[i]]$indexes] = i
  }
  
  impact_filt = ODDy@impact[ODDy@impact$impact == impact_type,]
  impact_filt$median = apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) median(as.numeric(x)))
  impact_filt$mean = apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) mean(as.numeric(x)))
  impact_filt[paste0('q', c(05, 10, 20, 80, 90, 95))] = t(apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) quantile(as.numeric(x), c(.05, 0.10, 0.20, 0.80, 0.90, 0.95))))
  
  # poly_ids = as.data.frame(ODDy$polygon_id)
  # poly_ids$order = 1:nrow(poly_ids)
  # poly_ids_impact = merge(poly_ids,impact_filt[,c('polygon', 'median')], by.x='polygon_id', by.y='polygon', all.x=T)
  # poly_ids_impact = poly_ids_impact[order(poly_ids_impact$order),]
  # ODDy[[paste0(impact_type, '_median')]] = poly_ids_impact['median']
  # plot(ODDy[[paste0(impact_type, '_median')]])
  
  # Convert SpatialPolygonsDataFrame to a tidy format for ggplot2
  polygon_list_plot <- list()
  for (i in 1:length(polygons_list)){
    if (is.null(polygons_list[[i]]$sf_polygon$NAME_2)) next #polygons_list$polygons_list[[i]]$sf_polygon$NAME_2 = ''
    #if(grep("^[^,]+,[^,]+$", polygons_list[[i]]$polygon_name)) next
    polygon_list_plot[[i]] = st_as_sf(polygons_list[[i]]$sf_polygon)
    poly_match = which(tolower(gsub("\\s+", "", impact_filt$polygon_name)) == tolower(gsub("\\s+", "", polygons_list[[i]]$polygon_name)))
    if (length(poly_match)>0){
      polygon_list_plot[[i]]$impact = pull(impact_filt[summary])[poly_match]
    } else {
      polygon_list_plot[[i]]$impact = 0
    }
  }
  combined_polygons <- do.call(rbind, polygon_list_plot)
  
  p_impact <- ggplot(data = combined_polygons) +
    geom_sf(aes(fill = impact), color = "black", size = 0.2) +
    scale_fill_gradientn(colors = c("#FFF5E1", "#FDBB84", "#D73027"), trans = "sqrt")+  # Color scale
    theme_minimal() +
    geom_sf_text(aes(label = NAME_2), size = 3, color = "black") + 
    labs(title = "Impact by Region", fill = paste(impact_type, summary)) + 
    xlim(plot_bbox[c(1,2)]) + ylim(plot_bbox[c(3,4)])
  
  hazard = apply(values(ODDy[grep("hazMean",names(ODDy),value = T), drop=F]),1,max, na.rm=T)
  ODDy$hazard<-hazard
  brks<-seq(9,ceiling(2*max(hazard,na.rm = T)),by=1)/2
  
  ODDy_df <- as.data.frame(ODDy, na.rm=F, xy=T)
  names(ODDy_df)[which(names(ODDy_df)=='x')] = 'Longitude'
  names(ODDy_df)[which(names(ODDy_df)=='y')] = 'Latitude'
  
  p<-p_impact+geom_contour(data = ODDy_df,
                           mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
                           alpha=0.5,breaks = brks, size=1) +
    scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
    labs(colour = "Hazard Intensity")
  
  return(p)
}

plot_GADM_impact_polygons_pixellated <- function(ODDy, polygons_list, gadm_level = 2, impact_type = 'mortality', summary = 'mean', plot_bbox = NULL, max_break = 500){
  #GADM_level 0 = national
  #GADM_level 1 = admin level 1
  #GADM_level 2 = admin level 2
  #summary = mean, median, q5, q10, q20, q80, q90, q95
  if (is.null(plot_bbox)){plot_bbox = ext(ODDy)}
  
  
  
  
  plot_df <- as.data.frame(ODDy, xy=T, na.rm=F)
  #plot_df[,grep(impact_type_samples, names(plot_df))] = t(apply(as.matrix(plot_df[, grep(impact_type_samples, names(plot_df)), drop = FALSE]), MARGIN=1, sort))
  
  var = 'mortality.mean'
  #plot_df[which(is.na(plot_df$ISO3C)), var] <- NA
  names(plot_df)[which(names(plot_df)=='x')] = 'Longitude'
  names(plot_df)[which(names(plot_df)=='y')] = 'Latitude'
  
  bbox <- matrix(ODDy@hazinfo$bbox, nrow=2)
  #gadm_iso <- getData("GADM", country="NZL", level=2)
  iso3_unique <- unique(ODDy$ISO3C)
  iso3_unique <- iso3_unique[!is.na(iso3_unique)]
  gadm_iso <- as(geodata::gadm(country=iso3_unique[1], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial')
  if (length(iso3_unique) > 1){
    for (i in 2:length(iso3_unique)){
      gadm_iso %<>% rbind(as(geodata::gadm(country=iso3_unique[i], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial'))
    }
  }
  
  
  #gadm_iso <- gSimplify(gadm_iso, 0.01)
  gadm_iso <- intersect(gadm_iso, bbox)
  gadm_map <- fortify(gadm_iso)
  
  gadm_map_cropped <- gadm_map %>%
    filter(long >= plot_bbox[1], long <= plot_bbox[2],
           lat  >= plot_bbox[3], lat  <= plot_bbox[4])
  
  p <- ggplot() + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
  
  p <- p + geom_map(map=gadm_map_cropped, data=gadm_map_cropped, aes(map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
  p <- p + coord_map() + geom_polygon(data=gadm_map_cropped, aes(x=long, y=lat, group=group), fill='white', color='black')
  
  if (impact_type == 'mortality'){
    breaks = c(-Inf, 1, 3, 5, 10, 20, 50, 100, 500, Inf)
    if (max_break < 500 & max_break > 100){
      breaks = c(-Inf, 1, 3, 5, 10, 20, 50, 100, Inf)
    }
    if (max_break < 100){
      breaks = c(-Inf, 1, 3, 5, 10, 20, 50, Inf)
    }
    if (max_break < 50){
      breaks = c(-Inf, 1, 5, 10, 20, 30, Inf)
    }
    if (max_break < 25){
      breaks = c(-Inf, 1, 5, 10, 15, 20, Inf)
    }
    if (max_break < 17){
      breaks = c(-Inf, 1, 4, 7, 10, 13, 16, Inf)
    }
    
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
  } else {
    breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, 10000, 50000, Inf)
    if (max_break < 50000 & max_break > 10000){
      breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, 10000, 50000, Inf)
    }
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
  }
  
  bin_colors <- c(viridis::viridis(length(levels(plot_df$bin_category))-2), 'orange', 'red')
  
  
  # Modify the plot
  map_plot <- p + 
    geom_raster(data = plot_df, aes(x = Longitude, y = Latitude, fill = bin_category), alpha = 0.75) + 
    coord_equal() +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal() + 
    theme(
      axis.title = element_text(family = "Liberation Serif", size = 12),  
      legend.text = element_text(family = "Liberation Serif", size = 11),
      legend.title = element_text(family = "Liberation Serif", size = 12)
    ) + 
    scale_fill_manual(name = paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type)), ' mean'),#paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type)),' ', quantile,'% Quantile'), 
                      values = bin_colors, na.value = "transparent") + # Discrete color scale 
    xlim(plot_bbox[c(1,2)]) + ylim(plot_bbox[c(3,4)])
  
  
  # Define the bounding box
  full_extent <- st_as_sfc(st_bbox(c(
    xmin = plot_bbox[1],
    xmax = plot_bbox[2],
    ymin = plot_bbox[3],
    ymax = plot_bbox[4]
  ), crs = st_crs(gadm_iso)))
  
  mask_poly <- st_difference(full_extent, st_union( st_as_sf(gadm_iso)))
  
  polygon_list_plot <- list()
  for (i in 1:length(polygons_list)){
    if (is.null(polygons_list[[i]]$sf_polygon$NAME_2)) next #polygons_list$polygons_list[[i]]$sf_polygon$NAME_2 = ''
    polygon_list_plot[[i]] = st_as_sf(polygons_list[[i]]$sf_polygon)
  }
  combined_polygons <- do.call(rbind, polygon_list_plot)
  
  p = map_plot + geom_sf(data=mask_poly,fill='lightgrey', color=NA) + 
    geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='white', alpha=0,color='black') +
    geom_contour(data = plot_df, aes(Longitude, Latitude, z = hazMean1, colour = ..level..),
                 alpha = 0.7, lwd = 0.8) +
    scale_color_gradientn(colors = c("transparent", "#fc9272", "#ef3b2c")) +
    labs(colour = "Hazard Intensity") +
    geom_sf_text(data=combined_polygons, aes(label = NAME_2), size = 3, color = "white", fontface = "bold") + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + coord_sf(
      xlim = plot_bbox[c(1,2)],
      ylim = plot_bbox[c(3,4)],
      expand = FALSE
    )
  
  return(p)
}

#Pairplot 
pairplot_regions <- function(impact_mort, cross_size = 3){
  
  sampled_df <- (impact_mort %>%
                   dplyr::select(polygon_name, starts_with("sampled.")) %>%
                   pivot_longer(cols = -polygon_name, names_to = "sample", values_to = "value") %>%
                   pivot_wider(names_from = polygon_name, values_from = value) %>%
                   mutate(type = "Sampled"))[,-1]
  
  sampled_df[1:(ncol(sampled_df)-1)] <- lapply(sampled_df[1:(ncol(sampled_df)-1)], as.numeric)
  
  # Step 2: Get observed values, same format
  observed_row <- impact_mort %>%
    dplyr::select(polygon_name, observed) %>%
    pivot_wider(names_from = polygon_name, values_from = observed) %>%
    mutate(type = "Observed")
  
  # Step 3: Ensure column types match
  # Make sure all except `type` are numeric
  observed_row[1:(ncol(observed_row)-1)] <- lapply(observed_row[1:(ncol(observed_row)-1)], as.numeric)
  
  # Step 4: Combine
  plot_df <- bind_rows(sampled_df, observed_row)
  
  obs_vals <- filter(plot_df, type == "Observed")
  
  # Function to overlay observed points
  overlay_observed <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(alpha = 0.7, size = 0.5, color = "#440154", ...) +
      geom_point(data = obs_vals, mapping = mapping, color = "red", size = cross_size, shape=4, stroke = 1.5)
  }
  # overlay_observed <- function(data, mapping, ...) {
  #   ggplot(data = data, mapping = mapping) +
  #     geom_density_2d_filled(alpha = 0.8, contour_var = "density") +
  #     geom_point(data = obs_vals, mapping = mapping, color = "red", size = 3, shape = 4, stroke = 1.5)
  # }
  
  gpairs_lower <- function(g){
    g$plots <- g$plots[-(1:g$nrow)]
    g$yAxisLabels <- g$yAxisLabels[-1]
    g$nrow <- g$nrow -1
    
    g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
    g$xAxisLabels <- g$xAxisLabels[-g$ncol]
    g$ncol <- g$ncol - 1
    
    g
  }
  
  # Pair plot with custom overlay
  g <- ggpairs(plot_df %>% dplyr::select(-type),
               upper = list(continuous = "blank"),
               lower = list(continuous = overlay_observed),
               diag = list(continuous = "blankDiag"), 
               showStrips=T
  ) +
    theme_minimal(base_family = "Times") +
    theme(strip.text = element_text(size = 10.5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) + 
    scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 2, base = 2),
      breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
      labels = expression(phantom(1) * 0^phantom(0),phantom(1) * 10^phantom(1),10^2,10^3, 10^4, 10^5, 10^6, 10^7),
      minor_breaks = NULL
    ) + 
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 2, base = 2),
      breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
      labels = expression(phantom(1) * 0^phantom(0),phantom(1) * 10^phantom(1),10^2,10^3, 10^4, 10^5, 10^6, 10^7),
      minor_breaks = NULL
    )
  
  gpairs_lower(g)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------------------------------------------------------------------
#==========     Case Study 0: Philippines Earthquake, October 29th 2019  ============
#------------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

input<-list(
  sdate=as.Date("2019-10-28"),
  fdate=as.Date("2019-10-30"),
  iso3="PHL",
  datadir=dir,
  plotdir="Plots/"
)

input%<>%append(list(Model=Model,
                     PosteriorFileLoc='IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
                     Method=AlgoParams))

options(timeout = 500)
ODDy_with_namer <- prepareODD(dir, input, getGADMregions=T, folder_write=folder_write)
ODDy = ODDy_with_namer$ODDy
namer = ODDy_with_namer$namer

#ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20191029PHL_-1')

values(ODDy$nBuildings)[which(is.na(values(ODDy$Population)))] = NA
values(ODDy$Vs30)[which(is.na(values(ODDy$Population)))] = NA

#Paper:
par(mfrow = c(2, 4))
par(family = "Times")
# Plot first figure
plot(ODDy$Population, main = "Population", font.main = 1)
mtext("(a)", side = 3, line = 2.2, adj = -0.05, cex = 1)

plot(ODDy$nBuildings, main = "Building Count", font.main = 1)
mtext("(b)", side = 3, line = 2.2, adj = -0.05, cex = 1)

plot(ODDy$hazMean1, main = "Hazard 1 Mean", font.main = 1)
mtext("(c)", side = 3, line = 2.2, adj = -0.05, cex = 1)

plot(ODDy$hazMean2, main = "Hazard 2 Mean", font.main = 1)
mtext("(d)", side = 3, line = 2.2, adj = -0.05, cex = 1)

plot(ODDy$SHDI, main = "SHDI", type = 'continuous', font.main = 1)
mtext("(e)", side = 3, line = 2.2, adj = -0.05, cex = 1)

plot(ODDy$GNIc, main = "GNIc", type = 'continuous', font.main = 1)
mtext("(f)", side = 3, line =  2.2, adj = -0.05, cex = 1)

plot(ODDy$EQFreq, main = "EQFreq", font.main = 1)
mtext("(h)", side = 3, line = 2.2, adj = -0.05, cex = 1)

plot(ODDy$Vs30, main = "Vs30", font.main = 1)
mtext("(i)", side = 3, line = 2.2, adj = -0.05, cex = 1)

# plot(ODDy$PDens, main = "PDens")
# mtext("(j)", side = 3, line = 1, adj = 0, cex = 1)

par(mfrow=c(1,1))
#PHL_Data.pdf, 11 x 5 (or 1300 x 680)


#Thesis:

dev.off()
# windows(width=8, height=8); quartz(width=8, height=8); x11(width=8, height=8)


par(family = "Times")
par(
  mfrow = c(3, 3),
  family = "Times",
  cex = 1.05,
  # much tighter inner margins: (bottom, left, top, right)
  mar = c(2, 2, 1.2, 0.6),
  # no outer margins
  oma = c(0, 0, 0, 0),
  # pull axes in a bit and tighten label spacing
  mgp = c(1.4, 0.4, 0),
  xaxs = "i", yaxs = "i"
)

plot(ODDy$Population, font.main = 1, cex.main = 1)
mtext("(a) Population", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$nBuildings, font.main = 1, cex.main = 1)
mtext("(b) Building Count", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$hazMean1, font.main = 1, cex.main = 1)
mtext("(c) Hazard 1 MMI", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$hazMean2, font.main = 1, cex.main = 1)
mtext("(d) Hazard 2 MMI", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$SHDI, type = 'continuous', font.main = 1, cex.main = 1)
mtext("(e) SHDI", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$GNIc, type = 'continuous', font.main = 1, cex.main = 1)
mtext("(f) GNIc", side = 3, line =  2.45, adj = -0.05, cex = 1)

plot(ODDy$EQFreq, font.main = 1, cex.main = 1)
mtext("(h) EQFreq", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$Vs30, font.main = 1, cex.main = 1)
mtext("(i) Vs30", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$PDens, font.main = 1, cex.main = 1)
mtext("(j) PopDens", side = 3, line = 2.45, adj = -0.05, cex = 1)

#PHL_Data_all.pdf, 8.7 x 8 (or 1300 x 680)

# Plot first figure
ODDy_df <- as.data.frame(ODDy)
coords  <- as.data.frame(crds(ODDy, na.rm = FALSE))
names(coords) <- c("x", "y")
ODDy_df$x <- coords$x
ODDy_df$y <- coords$y

# Helper to make a single raster plot (legend off to keep spacing tight)
library(scales)   # for label_scientific
library(grid)     # for unit()

make_panel <- function(var, title) {
  ggplot(ODDy_df, aes(x = x, y = y, fill = .data[[var]])) +
    geom_raster() +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    scale_fill_viridis_c(na.value = "white",
                         labels = scales::label_comma()) +
    coord_fixed(expand = FALSE) +
    theme(
      text = element_text(family = "Times New Roman"),
      axis.text = element_text(size = 11, color="black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      panel.background = element_blank(),
      plot.title = element_text(size = 12.5, hjust = 0),
      plot.margin = margin(1, 50, 20, 1),
      
      legend.position = c(1.07, 0.5),              # (x, y) in NPC
      legend.justification = c(0, 0.5), 
      #legend.direction = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 5, 0, 0),
      legend.spacing = unit(0, "pt"),
      legend.key.height = unit(10, "pt"),
      legend.key.width  = unit(6, "pt"),
      legend.text = element_text(size = 11),
      legend.title = element_blank()
    ) +
    guides(fill = guide_colorbar(
      title = title,
      #direction = "horizontal",
      barheight = unit(130, "pt"),
      barwidth  = unit(6,  "pt"),
      ticks = TRUE
    ))
}

# Build the nine panels
p1 <- make_panel("Population",  "Population")
p2 <- make_panel("nBuildings",  "Building Count")
p3 <- make_panel("hazMean1",    "Hazard 1 MMI")
p4 <- make_panel("hazMean2",    "Hazard 2 MMI")
p5 <- make_panel("SHDI",        "SHDI")
p6 <- make_panel("GNIc",        "GNIc")
p7 <- make_panel("EQFreq",      "EQFreq")
p8 <- make_panel("Vs30",        "Vs30")
p9 <- make_panel("PDens",       "PopDens")

# Arrange in a tight 3x3 grid
dev.off()
# windows(8, 8); quartz(8, 8); x11(8, 8)  # optionally open a fresh square device
plot_grid(p1, p2, p3,
          p4, p5, p6,
          p7, p8, p9,
          ncol = 3, align = "hv", 
          labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain", 
          label_size = 12.5, label_x = 0.01, label_y = 1.015)

#8.5 x 8.25, PHL_Data_all.pdf

label_descaled <- function(scale_factor = 1, accuracy = NULL) {
  function(x) {
    vals <- x / scale_factor
    label_number(accuracy = accuracy, big.mark = ",", trim = TRUE)(vals)
  }
}

make_panel <- function(var, title, scale_factor = 1, accuracy = NULL) {
  # Build the legend title, e.g. "Population × 10^5" if scale_factor == 1e5
  title_expr <- if (scale_factor == 1) {
    title
  } else {
    # use plotmath so it renders as "title × 10^k"
    # k = log10(scale_factor)
    k <- round(log10(scale_factor))
    bquote(.(title) ~ "\u00D7" ~ 10^.(k))
  }
  
  ggplot(ODDy_df, aes(x = x, y = y, fill = .data[[var]])) +
    geom_raster() +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(
      na.value = "white",
      labels = label_descaled(scale_factor, accuracy)
    ) +
    coord_fixed(expand = FALSE) +
    #ggtitle(title_expr) +
    theme(
      text = element_text(family = "Liberation Serif"),
      axis.text = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      panel.background = element_blank(),
      plot.margin = margin(1, 1, 1, 1),
      
      # legend at bottom; title above the bar
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = title_expr,
      legend.text  = element_text(size = 8),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(10, 0, 0, 0),
      legend.spacing = unit(0, "pt"),
      legend.key.height = unit(10, "pt"),
      legend.key.width  = unit(6, "pt")
    ) +
    guides(fill = guide_colorbar(
      #title = element_blank(),              # title above
      title.position = "top",
      title.hjust = 0.5,                 # left-align title
      #direction = "horizontal",
      barheight = unit(6, "pt"),
      barwidth  = unit(100, "pt"),
      ticks = TRUE
    ))
}

p1 <- make_panel("Population", "Population", scale_factor = 1e5, accuracy = 0.1)
p2 <- make_panel("nBuildings", "Building Count", scale_factor = 1e4, accuracy = 1)
p3 <- make_panel("hazMean1",   "Hazard 1 MMI",  scale_factor = 1)
p4 <- make_panel("hazMean2",   "Hazard 2 MMI",  scale_factor = 1)
p5 <- make_panel("SHDI",       "SHDI",          scale_factor = 1,  accuracy = 0.01)
p6 <- make_panel("GNIc",       "GNIc",          scale_factor = 1e3, accuracy = 1)
p7 <- make_panel("EQFreq",     "EQFreq",        scale_factor = 1,  accuracy = 0.01)
p8 <- make_panel("Vs30",       "Vs30",          scale_factor = 1)
p9 <- make_panel("PDens",      "PopDens",       scale_factor = 1e3, accuracy = 1)

#dev.off()
plot_grid(p1, p2, p3,
          p4, p5, p6,
          p7, p8, p9,
          ncol = 3, align = "hv")


ODDyWithImpact <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Apr25Agg/ODDobjects/Train/EQ20191029PHL_125')
plot_GADM_impact_observed <- function(ODDy, polygons_list, gadm_level = 2, impact_type = 'mortality', plot_bbox = NULL, plot_zeros=T, haz_legend=F){
  #GADM_level 0 = national
  #GADM_level 1 = admin level 1
  #GADM_level 2 = admin level 2
  if (is.null(plot_bbox)){plot_bbox = ext(ODDy)}
  summary='observed'
  
  ODDy$polygon_id = 0
  for (i in 1:length(ODDy@polygons)){
    if (str_count(ODDy@polygons[[i]]$name, ',') != gadm_level | ODDy@polygons[[i]]$name == 'TOTAL') next
    ODDy$polygon_id[ODDy@polygons[[i]]$indexes] = i
  }
  
  impact_filt = ODDy@impact[ODDy@impact$impact == impact_type,]
  impact_filt$median = apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) median(as.numeric(x)))
  impact_filt$mean = apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) mean(as.numeric(x)))
  impact_filt[paste0('q', c(05, 10, 20, 80, 90, 95))] = t(apply(impact_filt[,grep('sampled', names(impact_filt))],1 , function(x) quantile(as.numeric(x), c(.05, 0.10, 0.20, 0.80, 0.90, 0.95))))
  
  # poly_ids = as.data.frame(ODDy$polygon_id)
  # poly_ids$order = 1:nrow(poly_ids)
  # poly_ids_impact = merge(poly_ids,impact_filt[,c('polygon', 'median')], by.x='polygon_id', by.y='polygon', all.x=T)
  # poly_ids_impact = poly_ids_impact[order(poly_ids_impact$order),]
  # ODDy[[paste0(impact_type, '_median')]] = poly_ids_impact['median']
  # plot(ODDy[[paste0(impact_type, '_median')]])
  
  # Convert SpatialPolygonsDataFrame to a tidy format for ggplot2
  polygon_list_plot <- list()
  for (i in 1:length(polygons_list)){
    if (gadm_level == 2){
      if (is.null(polygons_list[[i]]$sf_polygon$NAME_2)) next 
    } else if (gadm_level == 1){
      if (!is.null(polygons_list[[i]]$sf_polygon$NAME_2)) next 
    }
    #if (is.null(polygons_list[[i]]$sf_polygon$NAME_2) & gadm_level == 2) next #polygons_list$polygons_list[[i]]$sf_polygon$NAME_2 = ''
    #if(grep("^[^,]+,[^,]+$", polygons_list[[i]]$polygon_name)) next
    polygon_list_plot[[i]] = st_as_sf(polygons_list[[i]]$sf_polygon)
    poly_match = which(tolower(gsub("\\s+", "", impact_filt$polygon_name)) == tolower(gsub("\\s+", "", polygons_list[[i]]$polygon_name)))
    if (length(poly_match)>0){
      polygon_list_plot[[i]]$impact = pull(impact_filt[summary])[poly_match]
    } else {
      polygon_list_plot[[i]]$impact = NA
    }
  }
  combined_polygons <- do.call(rbind, polygon_list_plot)
  if (!plot_zeros){
    combined_polygons$impact[which(combined_polygons$impact==0)] = NA
  }
  if (impact_type=='displacement'){
    #combined_polygons = st_union(combined_polygons)
    combined_polygons$impact = impact_filt$observed
  }
  
  if (impact_type != 'displacement'){
    p_impact <- ggplot(data = combined_polygons) +
      geom_sf(aes(fill = impact), color = "white", size = 0.2) +
      scale_fill_viridis(limits = c(0, NA)) +  # Color scale
      theme_minimal() +
      #geom_sf_text(aes(label = NAME_2), size = 3, color = "black") + 
      labs(fill = ifelse(impact_type=='buildDam', 'Building Damage', ifelse(impact_type=='mortality', 'Mortality', 'Displacement'))) + 
      coord_sf(xlim = plot_bbox[c(1,2)],
               ylim = plot_bbox[c(3,4)],
               expand = FALSE) +
      theme(
        plot.margin = unit(c(5,5,5,5), "pt"),  # same margins across panels
        axis.title.y = element_text(family = "Times New Roman", size = 12),
        axis.text.x = element_text(family = "Times New Roman", size = 12),
        axis.text.y = element_text(family = "Times New Roman", size = 12),
        axis.title.x = element_text(family = "Times New Roman", size = 12),
        plot.title = element_text(family = "Times New Roman", size = 14),
        legend.text = element_text(family = "Times New Roman", size = 12),
        legend.title = element_text(family = "Times New Roman", size = 12)
      )
    
    p_impact <- p_impact +
      guides(fill = guide_colorbar(barheight = unit(4, "mm"),
                                   barwidth  = unit(35, "mm"),
                                   title.position = "top")) +
      theme(legend.spacing.y = unit(0, "mm"),
            legend.box.margin = margin(0,0,0,0))
    
    if (impact_type=='buildDam'){
      p_impact = p_impact + scale_fill_viridis(breaks=c(1000,10000,20000))
    }
  } else {
    p_impact <- ggplot(data = combined_polygons) +
      geom_sf(aes(fill = as.factor(impact)), color = '#440154', size = 0) +#color = '#440154', size = 0) +
      scale_fill_viridis_d() +  # Color scale
      theme_minimal() +
      #geom_sf_text(aes(label = NAME_2), size = 3, color = "black") + 
      labs(fill = ifelse(impact_type=='buildDam', 'Building Damage', ifelse(impact_type=='mortality', 'Mortality', 'Displacement'))) + 
      coord_sf(xlim = plot_bbox[c(1,2)],
               ylim = plot_bbox[c(3,4)],
               expand = FALSE) +
      theme(
        plot.margin = unit(c(5,5,5,5), "pt"),  # same margins across panels
        element_text(family = "Times New Roman"),
        axis.title.y = element_text(family = "Times New Roman", size = 12),
        axis.text.x = element_text(family = "Times New Roman", size = 12),
        axis.text.y = element_text(family = "Times New Roman", size = 12),
        axis.title.x = element_text(family = "Times New Roman", size = 12),
        plot.title = element_text(family = "Times New Roman", size = 14),
        legend.text = element_text(family = "Times New Roman", size = 12),
        legend.title = element_text(family = "Times New Roman", size = 12))
    
    p_impact <- p_impact +
      guides(fill = guide_legend(nrow = 1,
                                 keyheight = unit(4, "mm"),
                                 keywidth  = unit(6, "mm"),
                                 title.position = "top")) +
      theme(legend.spacing.y = unit(0, "mm"),
            legend.box.margin = margin(0,0,0,0))
  }
  
  #hazard = apply(values(ODDy[grep("hazMean",names(ODDy),value = T), drop=F]),1,max, na.rm=T)
  hazard = values(ODDy$hazMean1)
  hazard[hazard < 4.98] = NA
  ODDy$hazard<-hazard
  brks<-seq(9,ceiling(2*max(hazard,na.rm = T)),by=1)/2
  
  ODDy_df <- as.data.frame(ODDy, na.rm=F, xy=T)
  names(ODDy_df)[which(names(ODDy_df)=='x')] = 'Longitude'
  names(ODDy_df)[which(names(ODDy_df)=='y')] = 'Latitude'
  
  x_breaks <- seq( floor(plot_bbox[1] * 2) / 2,
                   ceiling(plot_bbox[2] * 2) / 2, by = 0.5 )
  y_breaks <- seq( floor(plot_bbox[3] * 2) / 2,
                   ceiling(plot_bbox[4] * 2) / 2, by = 0.5 )
  
  p<-p_impact+geom_contour(data = ODDy_df,
                           mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
                           alpha=0.8,breaks = brks, size=1, 
                           show.legend = ifelse(haz_legend, T, F)) +
    scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
    labs(colour = "MMI") + 
    scale_x_continuous(breaks = x_breaks, labels = scales::number_format(accuracy = 0.5), limits = plot_bbox[c(1,2)], expand = c(0,0)) +
    scale_y_continuous(breaks = y_breaks, labels = scales::number_format(accuracy = 0.5), limits = plot_bbox[c(3,4)], expand = c(0,0)) +
    theme(
      panel.grid = element_blank(),  # Remove gridlines
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = element_line(),   # Ensure ticks are drawn
      axis.text = element_text() ,    # Ensure axis labels are visible
      element_text(family = "Times New Roman"),
    )
  
  return(p)
}



polygons_list = readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp_polygons/EQ20191029PHL_-1polygons')$polygons_list

ODDyWithImpact@impact$polygon_name = unlist(lapply(ODDyWithImpact@polygons[ODDyWithImpact@impact$polygon], function(x) x$name))

ODDyWithImpact@impact$observed[ODDyWithImpact@impact$observed==0 & ODDyWithImpact@impact$impact=='mortality'] = NA
p1 <- plot_GADM_impact_observed(ODDyWithImpact, polygons_list, impact_type = 'mortality', gadm_level = 1, plot_zeros=T)
p2 <- plot_GADM_impact_observed(ODDyWithImpact, polygons_list, impact_type = 'buildDam', gadm_level = 2, plot_zeros=F)
p3 <- plot_GADM_impact_observed(ODDyWithImpact, polygons_list, impact_type = 'displacement', gadm_level = 1, plot_zeros=T)
p3_with_hazlegend =  plot_GADM_impact_observed(ODDyWithImpact, polygons_list, impact_type = 'displacement', gadm_level = 1, plot_zeros=T, haz_legend=T)

legend_contour <- get_legend(
  p3_with_hazlegend +
    guides(fill = "none") +  # remove fill legend
    theme(legend.position = "right")
)

plot_grid(p1 + theme(legend.position="bottom",
                     legend.spacing = unit(0.5, "cm")), 
          p2 + theme(legend.position="bottom",
                     legend.spacing = unit(0.5, "cm")), 
          p3+ theme(legend.position="bottom",
                     legend.spacing = unit(0.5, "cm")), plot_grid(legend_contour, NULL, nrow=2, rel_heights=c(0.99, 0.01)), ncol=4, rel_widths=c(1,1,1,0.35), labels=c('(a)', '(b)', '(c)', ''),
          label_fontfamily = "Times New Roman",
          align = "hv", axis = "tblr",
          label_fontface = "plain")
#PHL_impact.pdf, 9 x 4

# Omega = list(Lambda1 = list(nu=9.6, kappa=1.9049508),
#                       Lambda2 = list(nu=9.8626452, kappa=1.9049508),
#                       Lambda3 = list(nu=9.3, kappa=0.9),
#                       Lambda4 = list(nu=9.9, kappa=1.9),
#                       theta= list(theta1=0.6),
#                       eps = list(local=1.2, hazard_mort=0.8383464, hazard_disp=0.9, hazard_bd=0.9, hazard_cor=0.55),
#                       vuln_coeff = list(PDens=0, SHDI=-0.5, GNIc=-0.1, Vs30=0.1, EQFreq=-0.1, FirstHaz=0.05, Night=0.05, FirstHaz.Night=0.1),
#                       check = list(check=0.5))

AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-04-13_120740.076549_MCMC_BDScore.05nocorr_M100_Npart1000NovAgg5_RandomFieldThree_rfNoTot')
Omega = AlgoResults$Omega_sample_phys[,1000] %>% relist(skeleton=Model$skeleton)

source('RCode/Model Variations/ODDobj_NoTot.R')

sampled_full = DispX(ODDy, Omega %>% addTransfParams(), Model$center, Method = list(Np=1, m_CRPS=1, cores=1, NestedCores=1), output='SampledFull')

ODDy$Displacement = sampled_full[,1,1]
values(ODDy$Displacement)[which(is.na(values(ODDy$Population)))] = NA
ODDy$Mortality = sampled_full[,1,2]
values(ODDy$Mortality)[which(is.na(values(ODDy$Population)))] = NA
ODDy$buildDam = sampled_full[,1,3]
values(ODDy$buildDam)[which(is.na(values(ODDy$Population)))] = NA

par(mfrow=c(1,3))
plot(ODDy$Displacement, main = "Displacement", type = 'continuous', font.main = 1)
mtext("(a)", side = 3, line =  2.2, adj = -0.05, cex = 1)

plot(ODDy$Mortality, main = "Mortality", type = 'continuous', font.main = 1)
mtext("(b)", side = 3, line =  2.2, adj = -0.05, cex = 1)

plot(ODDy$buildDam, main = "Building Damage", type = 'continuous', font.main = 1)
mtext("(c)", side = 3, line =  2.2, adj = -0.05, cex = 1)
par(mfrow=c(1,1))

displacement_df <- as.data.frame(ODDy$Displacement, xy = TRUE, na.rm=F)
mortality_df    <- as.data.frame(ODDy$Mortality, xy = TRUE, na.rm=F)
buildDam_df     <- as.data.frame(ODDy$buildDam, xy = TRUE, na.rm=F)
names(displacement_df)[which(names(displacement_df)=='x')] = 'Longitude'
names(displacement_df)[which(names(displacement_df)=='y')] = 'Latitude'
names(mortality_df)[which(names(mortality_df)=='x')] = 'Longitude'
names(mortality_df)[which(names(mortality_df)=='y')] = 'Latitude'
names(buildDam_df)[which(names(buildDam_df)=='x')] = 'Longitude'
names(buildDam_df)[which(names(buildDam_df)=='y')] = 'Latitude'

# Rename the value column for consistency
names(displacement_df)[3] <- "value"
names(mortality_df)[3] <- "value"
names(buildDam_df)[3] <- "value"


hazard = values(ODDy$hazMean1)
hazard[hazard < 4.98] = NA
ODDy$hazard<-hazard
brks<-seq(9,ceiling(2*max(hazard,na.rm = T)),by=1)/2

displacement_df$hazard = ODDy_df$hazard
mortality_df$hazard = ODDy_df$hazard
buildDam_df$hazard = buildDam_df$hazard

# Create each plot with legend at bottom and labeled
p4 <-  ggplot(mortality_df, aes(x = Longitude, y = Latitude, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Mortality", na.value = "transparent", breaks=c(0,1,2)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_fixed() 
  # guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5)) +
  # geom_contour(aes(Longitude,Latitude,z=hazard,colour=..level..),
  #              alpha=0.8,breaks = brks, size=1, show.legend = F) +
  # scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
  # labs(colour = "MMI")

p5 <- ggplot(displacement_df, aes(x = Longitude, y = Latitude, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Displacement", na.value = "transparent", breaks=c(0,1000,2000)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_fixed() 
  # guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5)) +
  # geom_contour(aes(Longitude,Latitude,z=hazard,colour=..level..),
  #              alpha=0.8,breaks = brks, size=1, show.legend = F) +
  # scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
  # labs(colour = "MMI")


p6 <- ggplot(buildDam_df, aes(x = Longitude, y = Latitude, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Building Damage", na.value = "transparent") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_fixed() 
  # guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5)) +
  # geom_contour(aes(Longitude,Latitude,z=hazard,colour=..level..),
  #              alpha=0.8,breaks = brks, size=1, show.legend = F) +
  # scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
  # labs(colour = "MMI")


# plot_grid(p1 + theme(legend.position="bottom",
#                      legend.spacing = unit(0.5, "cm")), 
#           p2 + theme(legend.position="bottom",
#                      legend.spacing = unit(0.5, "cm")), 
#           p3 + theme(legend.position="bottom",
#                      legend.spacing = unit(0.5, "cm")), 
#           legend_contour,
#           p4, 
#           p5, 
#           p6,  ncol=4, rel_widths=c(1,1,1,0.3), labels=c('(a)', '(b)', '(c)', '', '(d)', '(e)', '(f)'),
#           label_fontfamily = "Times New Roman",
#           label_fontface = "plain")
plot_grid(p4, 
          p5, 
          p6,  ncol=3, labels=c('(a)', '(b)', '(c)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

#PHL_SampledImpact.pdf, 9 x 3


#Prior Predicitive Checks

# AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-04-13_120740.076549_MCMC_BDScore.05nocorr_M100_Npart1000NovAgg5_RandomFieldThree_rfNoTot')
# s_finish = which(is.na(AlgoResults$Omega_sample_phys[1,]))[1]-1
# HLP = rep(0, s_finish)
# Model$higherpriors = T
# for (j in 1:s_finish){
#   HLP[j] = Model$HighLevelPriors(AlgoResults$Omega_sample_phys[,j] %>% relist(skeleton=Model$skeleton) %>% addTransfParams() , Model)
# }
# plot(HLP, type='l')


# Prior Predictive Checks:
plot_priors = function(n_samples){
  I_ij <- seq(4.5, 10, 0.01)
  n_samples = 50
  # Collect results into a tidy data frame
  plot_data <- map_dfr(1:n_samples, function(i) {
    params <- HLPrior_sample(Model, AlgoParams)
    Omega <- params %>%
      Array2Physical(Model) %>%
      addTransfParams()
    
    MortDisp <- D_MortDisp_calc(h_0(I_ij, Model$I0, Omega = Omega), Omega)
    BuildDam <- D_Dam_calc(h_0(I_ij, Model$I0, Omega = Omega), Omega)
    
    tibble(
      I = I_ij,
      Mort = MortDisp[1,],
      Disp = MortDisp[2,],
      BuildDam = BuildDam,
      sample = i
    )
  })
  
  ggplot() +
    geom_line(data = plot_data, aes(x = I, y = Mort, group = sample, color = "Mortality"), alpha = 0.5) +
    geom_line(data = plot_data, aes(x = I, y = Disp, group = sample, color = "Displacement"), alpha = 0.5) +
    geom_line(data = plot_data, aes(x = I, y = BuildDam, group = sample, color = "Building Damage"), alpha = 0.5) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Mortality" = "#440154",
        "Displacement" = "#1FA187",
        "Building Damage" = "#FDE725"
      ),
      guide = guide_legend(override.aes = list(size = 4, alpha=1, lwd=1))  # Thicker lines in legend
    ) +
    labs(x = "MMI", y = "Impact Probability") +
    coord_cartesian(xlim = c(4.5, 10), ylim = c(0, 1)) +
    theme_minimal(base_family = "Times New Roman") +
    scale_x_continuous(limits=c(4.5, 10), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    theme(
      legend.position = "bottom",
      text = element_text(family = "Times New Roman"),
      axis.title = element_text(family = "Liberation Serif", size=12),
      axis.text = element_text(family = "Liberation Serif", size=12),
      legend.text = element_text(family = "Liberation Serif", size=12),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
}
p1 = plot_priors(50)
p1

#create an AlgoResults object with prior samples:
AlgoResults <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Results/HPC/mcmc_2025-07-31_101619.999807__July25Agg_Normal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_Chain3_backup')
AlgoResultsPrior = AlgoResults
for (i in 1:(which(is.na(AlgoResultsPrior$Omega_sample_phys[1,]))[1]-1)){
  AlgoResultsPrior$Omega_sample[,i] = as.numeric(HLPrior_sample(Model, AlgoParams))
  AlgoResultsPrior$Omega_sample_phys[,i] = unlist(Array2Physical(AlgoResultsPrior$Omega_sample[,i], Model))
}
AlgoResultsPrior$RF_current = NULL
saveRDS(AlgoResultsPrior, paste0(dir, 'IIDIPUS_Results/AlgoResultsPriors'))

ODDy = readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20191029PHL_-1')
namer = 'EQ20191029PHL_-1'

ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename='IIDIPUS_Results/AlgoResultsPriors',
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200)

ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

mortality_true = c('Davao del Sur' = 2, 'South Cotabato' = 2, 'North Cotabato' = 8, 'Bukidnon' = 0)

displacement_true = c('Philippines' = 223105)

buildDam_true = c(
  #'Bansalan, Davao del Sur' = 3675,
  #'Hagonoy, Davao del Sur' = 842,
  #'Kiblawan, Davao del Sur' = 52,
  'Magsaysay, Davao del Sur' = 7402,
  #'Matanao, Davao del Sur' = 725,
  #'Padada, Davao del Sur' = 375,
  #'Santa Cruz, Davao del Sur' = 11,
  #'Sulop, Davao del Sur' = 150,
  'Malita, Davao del Sur' = 1,
  'Santa Maria, Davao del Sur' = 200,
  #'Antipas, North Cotabato' = 312,
  #'Arakan, North Cotabato' = 192,
  #'Kabacan, North Cotabato' = 261,
  'Kidapawan City, North Cotabato' = 5697,
  #'Magpet, North Cotabato' = 974,
  'Makilala, North Cotabato' = 20704,
  #'Matalam, North Cotabato' = 534,
  #'MLang, North Cotabato' = 1991,
  #'President Roxas, North Cotabato' = 259,
  #'Tulunan, North Cotabato' = 3063,
  #'Tupi, South Cotabato' = 18,
  'Columbio, Sultan Kudarat' = 38
)



#names(mortality_true) = paste0(names(mortality_true), ', Philippines')

observed_df <- data.frame(
  polygon_name2 = names(mortality_true),
  observed = as.vector(mortality_true),
  stringsAsFactors = FALSE,
  polygon_name = paste0(gsub("\\s+", "", names(mortality_true)), ", Philippines"),
  impact = 'mortality'
)

observed_df %<>% rbind(data.frame(
  polygon_name2 = paste0(names(buildDam_true), ', Philippines'),
  observed = as.vector(buildDam_true),
  stringsAsFactors = FALSE,
  polygon_name = paste0(
    sub(",", ", ", gsub("\\s+", "", names(buildDam_true)), fixed = TRUE),
    ", Philippines"
  ),
  impact='buildDam'
))

observed_df %<>% rbind(data.frame(
  polygon_name2 = names(displacement_true),
  observed = as.vector(displacement_true),
  stringsAsFactors = FALSE,
  polygon_name = names(displacement_true),
  impact='displacement'
))

impact_mort = ODDy@impact #%>% filter(impact=='mortality')
impact_mort$observed = NULL


impact_mort %<>% merge(observed_df, by=c('polygon_name', 'impact'))
impact_mort$observed[which(is.na(impact_mort$observed))] = 0
impact_mort <- impact_mort[
  order(factor(impact_mort$impact, levels = c("displacement", "buildDam", "mortality")),
        -impact_mort$observed),
]

impact_mort$polygon_name = impact_mort$polygon_name2
impact_mort$polygon_name2 = NULL

#impact_mort %<>% filter(polygon_name %in% paste0(names(mortality_true), ', Philippines'))

impact_mort_shortened_names = impact_mort
impact_mort_shortened_names$polygon_name <- sub(",.*", "", impact_mort_shortened_names$polygon_name)
#pairplot_regions(impact_mort_shortened_names)


plot_predictive_dist_subnat_PHL <- function(impact_mort){
  #plots predictive distributions over subnational regions
  
  impact_mort$polygon_name = sub("(^([^,]*),[^,]*),.*", "\\1", impact_mort$polygon_name)
  long_df <- impact_mort %>%
    pivot_longer(cols = starts_with("sampled."),
                 names_to = "sample",
                 values_to = "value") %>%
    mutate(sample = as.integer(str_remove(sample, "sampled\\.")))  # Optional: just cleans sample id
  
  long_df$value %<>% as.numeric()
  
  impact_mort$polygon_name <- sub(",\\s*", ",\n", impact_mort$polygon_name)
  long_df$polygon_name <- sub(",\\s*", ",\n", long_df$polygon_name)
  
  long_df$polygon_name <- factor(long_df$polygon_name, levels = rev(impact_mort$polygon_name))
  impact_mort$polygon_name <- factor(impact_mort$polygon_name, levels = rev(impact_mort$polygon_name))
  
  # Plot using ggplot
  #more options: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
  ggplot(long_df, aes(x = value, y = polygon_name, fill=impact)) +
    geom_density_ridges(
      scale = 0.9,
      alpha=0.6,
      #fill = impact,
      color = 'black',
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 2, point_alpha = 0.8, alpha = 0.7,
    ) +
    scale_fill_manual(
      name = 'Impact Type',
      labels = c(
        "mortality" = 'Mortality',
        "displacement" = 'Displacement',
        "buildDam" = 'Building Damage'
      ),
      values = c(
        "mortality" = "#440154",
        "displacement" = "#1FA187",
        "buildDam" = "#FDE725"
      )) +
    geom_point(
      data = impact_mort,
      aes(x = observed, y = polygon_name),
      color = "red",
      size = 3,
      shape = 19,  # Solid diamond, or use shape = 16 for filled circle
      inherit.aes = FALSE
    ) +
    labs(
      x = "Impact Size",
      y = "Region"#,
      #title = "Posterior predictive distribution for mortality by municipality"
    ) +
    theme_minimal() + scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 10, base = 10),
      breaks = c(0, 100, 1000, 10000, 100000, 1000000, 10000000),
      labels = function(x) {
        ifelse(x == 0, "0", parse(text = paste0("10^", log10(x))))
      },
      minor_breaks = NULL
    ) + theme(
      axis.title = element_text(family = "Liberation Serif", size=12),
      axis.text = element_text(family = "Liberation Serif", size=12),
      legend.position = "bottom",
      text = element_text(family = "Times New Roman"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) 
}

subnat_predDist = plot_predictive_dist_subnat_PHL(impact_mort)
p2 <- subnat_predDist

plot_grid(p1, p2, ncol=2, nrow=1,rel_widths=c(0.4, 0.6), labels=c('(a)', '(b)'), label_x = c(-0.02, -0.02), label_y = c(1.01,1.01),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")
#PriorPredCheck.pdf, 9.2 x 5.5


#-------------------------- RSS Prior Plot -----------------------

plot_predictive_dist_subnat_PHL <- function(impact_mort){
  #plots predictive distributions over subnational regions
  
  impact_mort$polygon_name = sub("(^([^,]*),[^,]*),.*", "\\1", impact_mort$polygon_name)
  long_df <- impact_mort %>%
    pivot_longer(cols = starts_with("sampled."),
                 names_to = "sample",
                 values_to = "value") %>%
    mutate(sample = as.integer(str_remove(sample, "sampled\\.")))  # Optional: just cleans sample id
  
  long_df$value %<>% as.numeric()
  
  #impact_mort$polygon_name <- sub(",\\s*", ",\n", impact_mort$polygon_name)
  #long_df$polygon_name <- sub(",\\s*", ",\n", long_df$polygon_name)
  
  long_df$polygon_name <- factor(long_df$polygon_name, levels = rev(impact_mort$polygon_name))
  impact_mort$polygon_name <- factor(impact_mort$polygon_name, levels = rev(impact_mort$polygon_name))
  
  # Plot using ggplot
  #more options: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
  ggplot(long_df, aes(x = value, y = polygon_name, fill=impact)) +
    geom_density_ridges(
      scale = 0.9,
      alpha=0.6,
      #fill = impact,
      color = 'black',
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 2, point_alpha = 0.8, alpha = 0.7,
      from = 0
    ) +
    scale_fill_manual(
      name = 'Impact Type',
      labels = c(
        "mortality" = 'Mortality',
        "displacement" = 'Displacement',
        "buildDam" = 'Building Damage'
      ),
      values = c(
        "mortality" = "#440154",
        "displacement" = "#1FA187",
        "buildDam" = "#FDE725"
      )) +
    geom_point(
      data = impact_mort,
      aes(x = observed, y = polygon_name),
      color = "red",
      size = 3,
      shape = 19,  # Solid diamond, or use shape = 16 for filled circle
      inherit.aes = FALSE
    ) +
    labs(
      x = "Impact Size",
      y = "Region"#,
      #title = "Posterior predictive distribution for mortality by municipality"
    ) +
    theme_minimal() + scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 10, base = 10),
      breaks = c(0, 100, 1000, 10000, 100000, 1000000, 10000000),
      labels = function(x) {
        ifelse(x == 0, "0", parse(text = paste0("10^", log10(x))))
      },
      minor_breaks = NULL
    ) + theme(
      axis.title = element_text(family = "Liberation Serif", size=12),
      axis.text = element_text(family = "Liberation Serif", size=12),
      legend.position = "bottom",
      text = element_text(family = "Times New Roman"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) 
}

subnat_predDist = plot_predictive_dist_subnat_PHL(impact_mort)
p2 <- subnat_predDist

plot_grid(p1 + theme(panel.border = element_blank()), 
          p2 + theme(panel.border = element_blank()), ncol=2, nrow=1,rel_widths=c(0.4, 0.6), labels=c('(a)', '(b)'), label_x = c(-0.02, -0.02), label_y = c(1.01,1.01),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")
#PriorPredCheck_zeroTrunc.pdf, 12 x 4.5

#------------------------------------------------------

# folderin <- paste0(dir, 'IIDIPUS_Input_Alternatives/Sept25/ODDobjects/')
# ufiles <- na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T))
# for (file in ufiles) {
#   print(file)
#   ODDy <- readODD(file.path(folderin, file))
#   plot(ODDy)
#   ans <- readline(prompt = "Press [Enter] for next, or type q to quit: ")
#   if (tolower(ans) == "q") break
# }

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects/Train/EQ20191029PHL_125')
plot(ODDy$hazMean1) # 350 x 350
plot(ODDy$SHDI + 0.00005*ODDy$Vs30) # 350 x 350
plot(ODDy$EQFreq)
ODDy$Fixed = 0
plot(ODDy$Fixed) # yellow

ODDy$Fixed = 0.5
values(ODDy$Fixed)[1] = 1
values(ODDy$Fixed)[2] = 0
plot(ODDy$Fixed) # green

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects/Test/EQ20171112IRN_70')
plot(ODDy$hazMean1) # 350 x 350
plot(ODDy$GNIc)
plot(ODDy$SHDI + 0.00005*ODDy$GNIc + 0.001*ODDy$Vs30 + 0.3 * ODDy$EQFreq)
# Error
grid <- as.data.frame(xyFromCell(ODDy, 1:ncell(ODDy)))  # Extract grid coordinates
names(grid) <- c("x", "y")  
vgm_model <- vgm(psill = 1, model = "Exp", range = 0.5)
set.seed(1)
gstat_mod = gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 3)
ODDy$Error = predict(gstat_mod, newdata = grid, nsim = 1)[,3]
plot(ODDy$Error)
ODDy$xx = ODDy$hazMean1 + 0.2*(ODDy$SHDI + 0.00005*ODDy$GNIc + 0.001*ODDy$Vs30 + 0.3 * ODDy$EQFreq) + ODDy$Error
plot(sqrt(exp(ODDy$xx)))
values(ODDy$xx)[which(is.na(values(ODDy$xx)))] = -3
plot(ODDy$xx)

sampledImpact = DispX(ODDy, Omega %>% addTransfParams(), center=Model$center, Method=AlgoParams, output='SampledFull')

ODDy$Mortality = sampledImpact[,1,1]
plot((ODDy$Mortality/100)^(1/3))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------------------------------------------------------------------
#==========     Case Study 1: Japan Earthquake, January 1st 2024  ===================
#------------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Prepare the ODD object, as in AutoQuake.R
input<-list(
  sdate=as.Date("2024-01-01"), 
  fdate=as.Date("2024-01-01"), 
  iso3="JPN", 
  datadir=dir,
  plotdir="Plots/"
)
# source('RCode/Model Variations/ODDobj_NoTot.R')
# #source('RCode/Model.R')
# source('RCode/Model Variations/Model_MVRrank.R')

#source('RCode/Model Variations/ODDobj_NoTot.R')

# input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
#                      PosteriorFileLoc='IIDIPUS_Results/mcmc_2025-05-05_225121.802448__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
#                      Method=AlgoParams))
# 
# input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
#                      PosteriorFileLoc='IIDIPUS_Results/HPC/mcmc_2025-05-15_211354.026585__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5_withHLPriors_flexHLP_chain2',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
#                      Method=AlgoParams))

# input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
#                      PosteriorFileLoc='IIDIPUS_Results/HPC/mcmc_2025-06-17_151849.678089__Apr25Agg_LogNormal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_withHLPriors_flexHLP',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_014135.935892__Apr25Agg_Lognormal_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
#                      Method=AlgoParams))

# #RSS:
# input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
#                      PosteriorFileLoc='IIDIPUS_Results/HPC/abcsmc_2025-08-03_093559.602748__July25Agg_NormalCDF_ESplus0.05MVR_RF_kappa0.5_range0.5',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_014135.935892__Apr25Agg_Lognormal_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
#                      Method=AlgoParams))

#Thesis:
input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
                     PosteriorFileLoc='IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_014135.935892__Apr25Agg_Lognormal_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
                     Method=AlgoParams))




#source('RCode/Model Variations/Model_Logistic.R')
#source('RCode/ODDobj.R')

options(timeout = 500)
ODDy_with_namer <- prepareODD(dir, input, getGADMregions=T, folder_write=folder_write)
ODDy = ODDy_with_namer$ODDy
namer = ODDy_with_namer$namer
#ODDy <- readODD(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240101JPN_-1'))


# Sample impact:
namer = 'EQ20240101JPN_-1'
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 500)

ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

# poly_names = unlist(lapply(ODDy@polygons, function(x) x$name))
# poly_names_df = data.frame(polygon_id = 1:length(poly_names), polygon_name = poly_names)

#ODDy <- readODD(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/EQ20240101JPN_-1_WithImpact200'))

# Manually load true data:
# mortality_true = c('Jōetsu, Niigata' = 12, 'Wajima, Ishikawa' = 180, 'Suzu, Ishikawa' = 104, 'Shika, Ishikawa' = 75,
#                    'Noto, Ishikawa' = 110, 'Nanao, Ishikawa' =  94, 'Nakanoto, Ishikawa' = 17,  'Kanazawa, Ishikawa' = 23, 
#                    'Hakui, Ishikawa' = 23, 'Anamizu, Ishikawa' = 87)
#Direct casualties: https://en.wikipedia.org/wiki/2024_Noto_earthquake
# mortality_true = c('Wajima, Ishikawa' = 101, 'Suzu, Ishikawa' = 97, 'Shika, Ishikawa' = 2,
#                     'Noto, Ishikawa' = 2, 'Nanao, Ishikawa' =  5,
#                     'Hakui, Ishikawa' = 1, 'Anamizu, Ishikawa' = 20)
# https://www.adrc.asia/view_disaster_en.php?Lang=en&Key=2653
mortality_true = c('Wajima, Ishikawa' = 80, 'Suzu, Ishikawa' = 54, 'Shika, Ishikawa' = 17,
                   'Noto, Ishikawa' = 49, 'Nanao, Ishikawa' =  37,
                   'Hakui, Ishikawa' = 3, 'Anamizu, Ishikawa' = 22,
                   'Komatsu, Ishikawa'=1, 'Hakusan, Ishikawa' = 1,
                   'Kahoku, Ishikawa'=5,'Nakanoto, Ishikawa'=1)

names(mortality_true) = paste0(names(mortality_true), ', Japan')

observed_df <- data.frame(
  polygon_name = names(mortality_true),
  observed = as.vector(mortality_true),
  stringsAsFactors = FALSE
)

impact_mort = ODDy@impact %>% filter(impact=='mortality')
impact_mort$observed = NULL


impact_mort %<>% merge(observed_df, by='polygon_name', all.x=T)
impact_mort$observed[which(is.na(impact_mort$observed))] = 0
impact_mort = impact_mort[order(impact_mort$observed, decreasing=T),]

impact_mort %<>% filter(polygon_name %in% c("Nakanoto, Ishikawa, Japan","Nanao, Ishikawa, Japan", "Noto, Ishikawa, Japan", "Shika, Ishikawa, Japan", "Suzu, Ishikawa, Japan",
                                            "Anamizu, Ishikawa, Japan", "Hakui, Ishikawa, Japan" , "Wajima, Ishikawa, Japan", "Hakui, Ishikawa, Japan", "Hodatsushimizu, Ishikawa, Japan",
                                            "Takaoka, Toyama, Japan", "Imizu, Toyama, Japan", "Komatsu, Ishikawa, Japan", "Hakusan, Ishikawa, Japan", "Kahoku, Ishikawa, Japan", "Himi, Toyama, Japan"))

impact_mort_shortened_names = impact_mort
impact_mort_shortened_names$polygon_name <- sub(",.*", "", impact_mort_shortened_names$polygon_name)
impact_mort_shortened_names %<>% filter(polygon_name != 'Hodatsushimizu')
pairplot_regions(impact_mort_shortened_names, cross_size=1.5)
#JPN_Pairplot.png, 1300 x 880
#JPN_Pairplot.pdf, 8.6 x 8.6, pdf

subnat_predDist = plot_predictive_dist_subnat(impact_mort)
#JPN_subnat_predDist.pdf, 8 x 4

# Load GADM regions:
polygons_list <- readRDS(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp_polygons/', namer,'polygons'))$polygons_list
plot_bbox = c(136, 138, 36.5, 38)
#plot_GADM_impact_polygons(ODDy, polygons_list, plot_bbox=plot_bbox)

p_pixellated <- plot_GADM_impact_polygons_pixellated(ODDy, polygons_list, plot_bbox=plot_bbox, max_break = 15)
#JPN_pixellated_mortPred.pdf, 11 x 9


grid.arrange(p_pixellated + theme(legend.position="bottom",
                     legend.spacing = unit(2, "cm")),
             subnat_predDist, ncol=2, widths=c(0.65, 0.35)
)



plot_grid(p_pixellated + theme(legend.position="bottom",
                               legend.spacing = unit(2, "cm")),
          subnat_predDist, ncol=2, nrow=1,rel_widths=c(0.5, 0.5), labels=c('(a)', '(b)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

#JPN_subnatPred, 14 x 8


subnat_predDist = plot_predictive_dist_subnat(impact_mort[1:12,])

plot_grid(p_pixellated +
            theme(
              legend.position      = "bottom",
              legend.box           = "vertical",   # stack multiple guides vertically
              legend.box.just      = "center",
              legend.spacing.y     = unit(6, "pt"),
              legend.spacing.x     = unit(0, "pt"),
              panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              legend.title.align = 0.5 
            ) +
            guides(
              # adjust to your aesthetics:
              fill   = guide_legend(direction = "horizontal", title.position = "top"),
              colour = guide_colorbar(direction = "horizontal", title.position = "top")
              # for discrete legends, use guide_legend() instead
            ), rel_widths =c(1,0.9),
          subnat_predDist + ylab(''), ncol=2, nrow=1, labels=c('(a)', '(b)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

#JPN_subnatPred, 9.2 x 7



obs1_i = 212
obs2_i = 220

joint_pred = data.frame(
  samples_imp1 = as.numeric(ODDy@impact[obs1_i, grep('sampled.', names(ODDy@impact))]),
  samples_imp2 = as.numeric(ODDy@impact[obs2_i, grep('sampled.', names(ODDy@impact))]),
  class = 'Initial Prediction'
)
joint_pred_updated = joint_pred %>% filter(samples_imp1 > 40 & samples_imp1 < 60)
joint_pred_updated$class = 'Updated Prediction'
joint_pred %<>% rbind(joint_pred_updated)

p1 <- ggplot(joint_pred %>% filter(class=='Initial Prediction'), aes(x = samples_imp1, y = samples_imp2), alpha=0.3) +
  geom_point(alpha = 0.8, color='#440154') +
  theme_minimal() +
  geom_point(aes(x = 49, y = 80, color = "Observed"), size = 4, shape = 4, stroke=1.5) +  # Crimson red point
  xlab("Mortality in Noto") +
  ylab("Mortality in Wajima") +
  theme(
    axis.title.y = element_text(family = "Times New Roman", size = 12),
    axis.text.x = element_text(family = "Times New Roman", size = 12),
    axis.text.y = element_text(family = "Times New Roman", size = 12),
    axis.title.x = element_text(family = "Times New Roman", size = 12),
    plot.title = element_text(family = "Times New Roman", size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0, 20, 0, 15), "pt"),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank()  # Remove legend title
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  scale_color_manual(
    values = c(# Viridis purple
      "Observed" = "#D62728"                      # Crimson red last
    ))

p2 <- ggplot(joint_pred, aes(x = samples_imp2, color = class, fill = class)) +
  geom_density(alpha = 0.2, adjust = 1) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  xlab("Predicted Mortality in Wajima") +
  ylab("Density") +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "Initial Prediction" = "#1FA187",
      "Updated Prediction" = "#440154"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Initial Prediction" = "#1FA187",
      "Updated Prediction" = "#440154"
    )
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank(),
    legend.text = element_text(family = "Times New Roman", size = 12)
  ) + 
  geom_vline(xintercept = 80, color = "#D62728", size = 1)


plot_grid(p1 + theme(legend.position = 'none', plot.margin = margin(t = 3, b = 10, l = 2, r = 10)), 
                 p2 + theme(legend.position='bottom', panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                                                                                  plot.margin = margin(t = 3, b = 10, l = 2, r = 10)), 
                 ncol=2, nrow=1,rel_widths=c(0.5, 0.5), labels=c('(a)', '(b)'),
                 label_fontfamily = "Times New Roman",
                 label_fontface = "plain",
                 label_x = c(-0.01, -0.02), align = "v",   # align by height
       axis  = "t" )

#JPN_subnatUpdate.png, 8 x 3.8

#-------------------------------------------------
#----------------- RSS Version -------------------
#-------------------------------------------------

plot_GADM_impact_polygons_pixellated <- function(ODDy, polygons_list, gadm_level = 2, impact_type = 'mortality', summary = 'mean', plot_bbox = NULL, max_break = 500){
  #GADM_level 0 = national
  #GADM_level 1 = admin level 1
  #GADM_level 2 = admin level 2
  #summary = mean, median, q5, q10, q20, q80, q90, q95
  if (is.null(plot_bbox)){plot_bbox = ext(ODDy)}
  
  
  
  
  plot_df <- as.data.frame(ODDy, xy=T, na.rm=F)
  #plot_df[,grep(impact_type_samples, names(plot_df))] = t(apply(as.matrix(plot_df[, grep(impact_type_samples, names(plot_df)), drop = FALSE]), MARGIN=1, sort))
  
  var = 'mortality.mean'
  #plot_df[which(is.na(plot_df$ISO3C)), var] <- NA
  names(plot_df)[which(names(plot_df)=='x')] = 'Longitude'
  names(plot_df)[which(names(plot_df)=='y')] = 'Latitude'
  
  bbox <- matrix(ODDy@hazinfo$bbox, nrow=2)
  #gadm_iso <- getData("GADM", country="NZL", level=2)
  iso3_unique <- unique(ODDy$ISO3C)
  iso3_unique <- iso3_unique[!is.na(iso3_unique)]
  gadm_iso <- as(geodata::gadm(country=iso3_unique[1], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial')
  if (length(iso3_unique) > 1){
    for (i in 2:length(iso3_unique)){
      gadm_iso %<>% rbind(as(geodata::gadm(country=iso3_unique[i], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial'))
    }
  }
  
  
  #gadm_iso <- gSimplify(gadm_iso, 0.01)
  gadm_iso <- intersect(gadm_iso, bbox)
  gadm_map <- fortify(gadm_iso)
  
  gadm_map_cropped <- gadm_map %>%
    filter(long >= plot_bbox[1], long <= plot_bbox[2],
           lat  >= plot_bbox[3], lat  <= plot_bbox[4])
  
  p <- ggplot() + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
  
  p <- p + geom_map(map=gadm_map_cropped, data=gadm_map_cropped, aes(map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
  p <- p + coord_map() + geom_polygon(data=gadm_map_cropped, aes(x=long, y=lat, group=group), fill='white', color='black')
  
  if (impact_type == 'mortality'){
    breaks = c(-Inf, 1, 3, 5, 10, 20, 50, 100, 500, Inf)
    if (max_break < 500 & max_break > 100){
      breaks = c(-Inf, 1, 3, 5, 10, 20, 50, 100, Inf)
    }
    if (max_break < 100){
      breaks = c(-Inf, 1, 3, 5, 10, 20, 50, Inf)
    }
    if (max_break < 50){
      breaks = c(-Inf, 1, 5, 10, 20, 30, Inf)
    }
    if (max_break < 25){
      breaks = c(-Inf, 1, 5, 10, 15, 20, Inf)
    }
    if (max_break < 17){
      breaks = c(-Inf, 1, 4, 7, 10, 13, 16, Inf)
    }
    
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
  } else {
    breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, 10000, 50000, Inf)
    if (max_break < 50000 & max_break > 10000){
      breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, 10000, 50000, Inf)
    }
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
  }
  
  bin_colors <- c(viridis::viridis(length(levels(plot_df$bin_category))))
  
  
  # Modify the plot
  map_plot <- p + 
    geom_raster(data = plot_df, aes(x = Longitude, y = Latitude, fill = bin_category), alpha = 0.75) + 
    coord_equal() +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal() + 
    theme(
      axis.title = element_text(family = "Liberation Serif", size = 12),  
      legend.text = element_text(family = "Liberation Serif", size = 11),
      legend.title = element_text(family = "Liberation Serif", size = 12)
    ) + 
    scale_fill_manual(name = paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type)), ' mean'),#paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type)),' ', quantile,'% Quantile'), 
                      values = bin_colors, na.value = "transparent") + # Discrete color scale 
    xlim(plot_bbox[c(1,2)]) + ylim(plot_bbox[c(3,4)])
  
  
  # Define the bounding box
  full_extent <- st_as_sfc(st_bbox(c(
    xmin = plot_bbox[1],
    xmax = plot_bbox[2],
    ymin = plot_bbox[3],
    ymax = plot_bbox[4]
  ), crs = st_crs(gadm_iso)))
  
  mask_poly <- st_difference(full_extent, st_union( st_as_sf(gadm_iso)))
  
  polygon_list_plot <- list()
  for (i in 1:length(polygons_list)){
    if (is.null(polygons_list[[i]]$sf_polygon$NAME_2)) next #polygons_list$polygons_list[[i]]$sf_polygon$NAME_2 = ''
    polygon_list_plot[[i]] = st_as_sf(polygons_list[[i]]$sf_polygon)
  }
  combined_polygons <- do.call(rbind, polygon_list_plot)
  
  p = map_plot + geom_sf(data=mask_poly,fill='lightgrey', color=NA) + 
    geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='white', alpha=0,color='black') +
    geom_contour(data = plot_df, aes(Longitude, Latitude, z = hazMean1, colour = ..level..),
                 alpha = 0.7, lwd = 0.8) +
    scale_color_gradientn(colors = c("transparent", "#fc9272", "#ef3b2c")) +
    labs(colour = "Hazard Intensity") +
    geom_sf_text(data=combined_polygons, aes(label = NAME_2), size = 3, color = "white", fontface = "bold") + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + coord_sf(
      xlim = plot_bbox[c(1,2)],
      ylim = plot_bbox[c(3,4)],
      expand = FALSE
    )
  
  return(p)
}


plot_predictive_dist_subnat <- function(impact_mort){
  #plots predictive distributions over subnational regions
  
  impact_mort$polygon_name = sub("(^([^,]*),[^,]*),.*", "\\1", impact_mort$polygon_name)
  long_df <- impact_mort %>%
    pivot_longer(cols = starts_with("sampled."),
                 names_to = "sample",
                 values_to = "value") %>%
    mutate(sample = as.integer(str_remove(sample, "sampled\\.")))  # Optional: just cleans sample id
  
  long_df$value %<>% as.numeric()
  
  long_df$polygon_name <- factor(long_df$polygon_name, levels = rev(impact_mort$polygon_name))
  impact_mort$polygon_name <- factor(impact_mort$polygon_name, levels = rev(impact_mort$polygon_name))
  
  long_df <- long_df %>%
    mutate(height = ifelse(polygon_name %in% levels(polygon_name)[11:15], 1.2,0.9))
  
  # Plot using ggplot
  #more options: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
  
  ggplot(long_df, aes(x = value, y = polygon_name, scale=height)) +
    geom_density_ridges(
      alpha=0.55,
      fill = '#440154',
      color = 'black',
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.04, height = 0),
      point_shape = '|', point_size = 2, point_alpha = 0.8, alpha = 0.7,
      from = 0
    ) +
    geom_point(
      data = impact_mort,
      aes(x = observed, y = polygon_name),
      color = "red",
      size = 3,
      shape = 19,  # Solid diamond, or use shape = 16 for filled circle
      inherit.aes = FALSE
    ) +
    labs(
      x = "Mortality",
      y = "Municipality"#,
      #title = "Posterior predictive distribution for mortality by municipality"
    ) +
    theme_minimal() + scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.5, base = 100),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      labels = scales::comma_format(),
      minor_breaks = NULL
    ) + theme(
      axis.title = element_text(family = "Liberation Serif", size=12),
      axis.text = element_text(family = "Liberation Serif", size=12)
    )
}

subnat_predDist = plot_predictive_dist_subnat(impact_mort)

p_pixellated <- plot_GADM_impact_polygons_pixellated(ODDy, polygons_list, plot_bbox=plot_bbox, max_break = 15)

plot_grid(p_pixellated + theme(legend.position="bottom",
                               legend.spacing = unit(1, "cm")) +
                         guides(fill = guide_legend(nrow = 2, byrow = TRUE)),
          subnat_predDist, ncol=2, nrow=1,rel_widths=c(0.65, 0.35), labels=c('(a)', '(b)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

#JPN_subnatPred_zeroTrunc, 14 x 8

#-------------------------------------------------

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------------------------------------------------------------------
#==========     Case Study 2: Morocco Earthquake, September 2023  ===================
#------------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#source('RCode/Model Variations/ODDobj_NoTot.R')

input<-list(
  sdate=as.Date("2023-09-07"), 
  fdate=as.Date("2023-09-09"), 
  iso3="MAR", 
  datadir=dir,
  plotdir="Plots/"
)

input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
                     PosteriorFileLoc='IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_014135.935892__Apr25Agg_Lognormal_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
                     Method=AlgoParams))

# input%<>%append(list(Model=Model, 
#                      PosteriorFileLoc='IIDIPUS_Results/mcmc_2025-05-05_225121.802448__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot',  # 'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2', #'IIDIPUS_Results/mcmc_2025-04-21_203713.178972_backup_NoTot',
#                      #PosteriorFileLoc='IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
#                      Method=AlgoParams))

options(timeout = 200)
ODDy_with_namer <- prepareODD(dir, input, getGADMregions=T, folder_write=folder_write)
ODDy = ODDy_with_namer$ODDy
namer = ODDy_with_namer$namer

#source('RCode/Model Variations/ODDobj_NoTot.R')
#source('RCode/Model Variations/Model_Logistic.R')
ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20230908MAR_-1')

# Sample impact:
namer = 'EQ20230908MAR_-1'
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

# ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/EQ20230908MAR_-1_WithImpact200')
# sampled_full <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/sampled_full/EQ20230908MAR_-1_fullSampledImpact')

mortality_true = c('Marrakech, Marrakech-Tensift-AlHaouz, Morocco' = 15,
                   'AlHaouz, Marrakech-Tensift-AlHaouz, Morocco' = 1684,
                   'Taroudannt, Souss-Massa-Draâ, Morocco' = 980, 
                   'Chichaoua, Marrakech-Tensift-AlHaouz, Morocco' = 202,
                   'Azilal, Tadla-Azilal, Morocco' = 11,
                   'Ouarzazate, Souss-Massa-Draâ, Morocco' = 41,
                   'Morocco' = 2946)

observed_df <- data.frame(
  polygon_name = names(mortality_true),
  observed = as.vector(mortality_true),
  stringsAsFactors = FALSE
)

impact_mort = ODDy@impact %>% filter(impact=='mortality')
impact_mort$observed = NULL

impact_mort %<>% merge(observed_df, by='polygon_name', all.x=T)
impact_mort$observed[which(is.na(impact_mort$observed))] = 0
impact_mort = impact_mort[order(impact_mort$observed, decreasing=T),]

impact_mort %<>% filter(polygon_name %in% c("Marrakech, Marrakech-Tensift-AlHaouz, Morocco","AlHaouz, Marrakech-Tensift-AlHaouz, Morocco", "Taroudannt, Souss-Massa-Draâ, Morocco", "Chichaoua, Marrakech-Tensift-AlHaouz, Morocco", 
                                            "Azilal, Tadla-Azilal, Morocco", "Ouarzazate, Souss-Massa-Draâ, Morocco"))

subnat_predDist = plot_predictive_dist_subnat(impact_mort)

impact_mort$polygon_name = c('Al Haouz', 'Taroudant', 'Chichaoua', 'Ouarzazate', 'Marrakech', 'Azilal')

pairplot_regions(impact_mort)
#MAR_pairplot.pdf, 12 x 12

plot(log(as.numeric(impact_mort[2,9:108])+10), log(as.numeric(impact_mort[3,9:108])+10))
points(log(impact_mort[2,109]+10), log(impact_mort[3,109]+10), col='red')
dat <- matrix(log(as.numeric(as.matrix(impact_mort[1:5, 9:109])) + 10), nrow=5,byrow=F)

mmdsss <- c()
esss <- c()
for (i in 1:101){
  mmdsss <- c(mmdsss, mmds_sample(dat[,i], dat[,-i]))
  esss <- c(esss, es_sample(dat[,i], dat[,-i]))
}
plot(mmdsss, esss)
points(mmdsss[101], esss[101], col='blue', pch=19)

polygons_list = readRDS(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp_polygons/EQ20230908MAR_-1polygons'))$polygons_list
plot_GADM_impact_polygons(ODDy, polygons_list, gadm_level = 2, impact_type = 'mortality', summary = 'mean')
  
plotODDy_GADM(ODDy, var='displacement', gadm_level=2)
# r <- rast(ODDy)
# grid <- as.data.frame(xyFromCell(ODDy, 1:ncell(ODDy)))  # Extract grid coordinates
# names(grid) <- c("x", "y")  # Name the columns
# vgm_model <- vgm(psill = 1, model = "Mat", range = 0.2, kappa = 0.2)
# gstat_mod = gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 3)
# 
# # set.seed(1)
# #RF_local <- as.matrix(predict(gstat_mod, newdata = grid, nsim = Method$Np)[, 3:(2+Method$Np)])
# RF_mort <- Omega$eps$local * as.matrix(predict(gstat_mod, newdata = grid, nsim = 1)[, 3:3])
# r$rf_mort = RF_mort
# plot(r)

  
displacement_true = c('TOTAL' = 500000)

obs1_i = 74
obs2_i = 73

joint_pred = data.frame(
  samples_imp1 = as.numeric(ODDy@impact[obs1_i, grep('sampled.', names(ODDy@impact))]),
  samples_imp2 = as.numeric(ODDy@impact[obs2_i, grep('sampled.', names(ODDy@impact))]),
  class = 'Initial Prediction'
)
joint_pred_updated = joint_pred %>% filter(samples_imp1 > 2000 & samples_imp1 < 4000)
joint_pred_updated$class = 'Updated Prediction'
joint_pred %<>% rbind(joint_pred_updated)

p1 <- ggplot(joint_pred %>% filter(class=='Initial Prediction'), aes(x = samples_imp1, y = samples_imp2), alpha=0.3) +
  geom_point(alpha = 0.8, color='#440154') +
  theme_minimal() +
  geom_point(aes(x = 2946, y = 500000, color = "Observed"), size = 4, shape = 4, stroke=1.5) +  # Crimson red point
  xlab("Total Mortality") +
  ylab("Total Population Displacement") +
  theme(
    axis.title.y = element_text(family = "Times New Roman", size = 12),
    axis.text.x = element_text(family = "Times New Roman", size = 12),
    axis.text.y = element_text(family = "Times New Roman", size = 12),
    axis.title.x = element_text(family = "Times New Roman", size = 12),
    plot.title = element_text(family = "Times New Roman", size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0, 20, 0, 15), "pt"),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank()  # Remove legend title
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  scale_color_manual(
    values = c(# Viridis purple
      "Observed" = "#D62728"                      # Crimson red last
    ))

p2 <- ggplot(joint_pred, aes(x = samples_imp2, color = class, fill = class)) +
  geom_density(alpha = 0.2, adjust = 1) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 100, base = 10),
    breaks = c(0, 1000, 10000, 100000, 1000000, 10000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  xlab("Predicted Population Displacement") +
  ylab("Density") +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "Initial Prediction" = "#1FA187",
      "Updated Prediction" = "#440154"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Initial Prediction" = "#1FA187",
      "Updated Prediction" = "#440154"
    )
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank(),
    legend.text = element_text(family = "Times New Roman", size = 12)
  ) + 
  geom_vline(xintercept = 500000, color = "#D62728", size = 1)


ggdraw(plot_grid(p1, p2, ncol=2, nrow=1,rel_widths=c(0.5, 0.5), labels=c('(a)', '(b)'),
                 label_fontfamily = "Times New Roman",
                 label_fontface = "plain")) + 
  theme(
    plot.margin = margin(t = 10, b = 10, l = 0, r = 0) # 20pt top and bottom margin
  )


obs1_i = 26
obs2_i = 74

joint_pred = data.frame(
  samples_imp1 = as.numeric(ODDy@impact[obs1_i, grep('sampled.', names(ODDy@impact))]),
  samples_imp2 = as.numeric(ODDy@impact[obs2_i, grep('sampled.', names(ODDy@impact))]),
  class = 'Initial Prediction'
)
joint_pred_updated = joint_pred %>% filter(samples_imp1 > 1000 & samples_imp1 < 2000)
joint_pred_updated$class = 'Updated Prediction'
joint_pred %<>% rbind(joint_pred_updated)

p1 <- ggplot(joint_pred %>% filter(class=='Initial Prediction'), aes(x = samples_imp1, y = samples_imp2), alpha=0.3) +
  geom_point(alpha = 0.8, color='#440154') +
  theme_minimal() +
  geom_point(aes(x = 1684, y = 2946, color = "Observed"), size = 4, shape = 4, stroke=1.5) +  # Crimson red point
  xlab("Mortality in Al Haouz") +
  ylab("Total Mortality") +
  theme(
    axis.title.y = element_text(family = "Times New Roman", size = 12),
    axis.text.x = element_text(family = "Times New Roman", size = 12),
    axis.text.y = element_text(family = "Times New Roman", size = 12),
    axis.title.x = element_text(family = "Times New Roman", size = 12),
    plot.title = element_text(family = "Times New Roman", size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0, 20, 0, 15), "pt"),
    legend.text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank()  # Remove legend title
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  scale_color_manual(
    values = c(# Viridis purple
      "Observed" = "#D62728"                      # Crimson red last
    ))

p2 <- ggplot(joint_pred, aes(x = samples_imp2, color = class, fill = class)) +
  geom_density(alpha = 0.2, adjust = 1) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 100, base = 10),
    breaks = c(0, 1000, 10000, 100000, 1000000, 10000000),
    labels = scales::comma_format(),
    minor_breaks = NULL
  ) +
  xlab("Predicted Total Mortality") +
  ylab("Density") +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "Initial Prediction" = "#1FA187",
      "Updated Prediction" = "#440154"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Initial Prediction" = "#1FA187",
      "Updated Prediction" = "#440154"
    )
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    legend.title = element_blank(),
    legend.text = element_text(family = "Times New Roman", size = 12)
  ) + 
  geom_vline(xintercept = 2946, color = "#D62728", size = 1)


ggdraw(plot_grid(p1, p2, ncol=2, nrow=1,rel_widths=c(0.5, 0.5), labels=c('(a)', '(b)'),
                 label_fontfamily = "Times New Roman",
                 label_fontface = "plain")) + 
  theme(
    plot.margin = margin(t = 10, b = 10, l = 0, r = 0) # 20pt top and bottom margin
  )

#MAR_MortDisp.png, 1350 x 450







#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------------------------------------------------------------------
#==========     Case Study 3: Turkey Earthquake, February 2023  =====================
#------------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Prepare the ODD object, as in AutoQuake.R
input<-list(
  sdate=as.Date("2023-02-05"), 
  fdate=as.Date("2023-02-07"), 
  iso3="TUR", 
  datadir=dir,
  plotdir="Plots/"
)

input%<>%append(list(Model=Model, 
                     PosteriorFileLoc='IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
                     Method=AlgoParams))

options(timeout = 500)
ODDy_with_namer <- prepareODD(dir, input, getGADMregions=T, folder_write=folder_write)
ODDy = ODDy_with_namer$ODDy
namer = ODDy_with_namer$namer

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200)

ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------------------------------------------------------------------
#==========     Case Study 3: Myanmar Earthquake, March 2025  =====================
#------------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Prepare the ODD object, as in AutoQuake.R
input<-list(
  sdate=as.Date("2025-03-28"), 
  fdate=as.Date("2025-03-28"), 
  iso3="MMR", 
  datadir=dir,
  plotdir="Plots/"
)

input%<>%append(list(Model=Model, 
                     PosteriorFileLoc='IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1',
                     Method=AlgoParams))

options(timeout = 500)
ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
ODDy = ODDy_with_namer$ODDy
namer = ODDy_with_namer$namer

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/EQ20250328MMR_-1_WithBuildings_WithImpactSummaries')
ODDy@impact = ODDy@impact[,1:8]
ODDy@impact$sdate = NULL
namer = 'EQ20250328MMR_-1'
ODDy@impact %<>% add_row(data.frame(impact='buildDam', iso3='MMR', observed=NA, qualifier=NA, build_type=NA, inferred=NA, polygon=476))


# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200)
sampled_full = ODDy_with_sampled$sampled_full
ODDy = ODDy_with_sampled$ODDyWithImpact

impact_plot = ODDy@impact %>% filter(polygon_name == 'Myanmar (Burma)')
impact_plot$observed = c(157000, 200000, 10000)


impact_plot$polygon_name = paste0(paste0(
  toupper(substring(impact_plot$impact, 1, 1)),
  substring(impact_plot$impact, 2)
), ':\n', impact_plot$polygon_name)

impact_plot = impact_plot[order(impact_plot$observed, decreasing=T),]

pairplot_regions(impact_plot)




mandalay_dam = st_read("/home/manderso/Downloads/skysat_20250329_080228_predictions_merged.gpkg")
bbox_obsdam = project(ext(mandalay_dam), from=crs(mandalay_dam), to="EPSG:4326")

overlapping_cells = cells(ODDy, bbox_obsdam)
plot(t(sampled_full[overlapping_cells[c(1,2)],3,]))

#mandalay_dam_to_ODDy <- st_transform(mandalay_dam, st_crs(crs(ODDy)))
#mandalay_dam_to_ODDy_v <- project(vect(mandalay_dam), crs(ODDy))

ODD_cropped <- crop(ODDy, bbox_obsdam)

ODDy_df = as.data.frame(ODD_cropped, na.rm=F)
var = 'nBuildings'
ODDy_df$x = crds(ODD_cropped, na.rm=F)[,1]
ODDy_df$y = crds(ODD_cropped, na.rm=F)[,2]

dam_pts_v <- centroids(mandalay_dam_v, inside = TRUE)

dam_pts <- cbind(
  as.data.frame(dam_pts_v, na.rm=F),           # attributes, incl. `damaged`
  crds(dam_pts_v)[, c("x","y")]  # coordinates
)

ggplot(ODDy_df) +
  geom_raster(aes(x = x, y = y, fill = .data[[var]])) +
  labs(title = '', x = NULL, y = NULL, fill = NULL) +
  scale_fill_viridis_c(na.value = "white",
                       labels = scales::label_comma()) +
  coord_fixed(expand = FALSE) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 11, color="black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    panel.background = element_blank(),
    plot.title = element_text(size = 12.5, hjust = 0),
    plot.margin = margin(1, 50, 20, 1),
    
    legend.position = c(1.07, 0.5),              # (x, y) in NPC
    legend.justification = c(0, 0.5), 
    #legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 5, 0, 0),
    legend.spacing = unit(0, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.key.width  = unit(6, "pt"),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white", colour = NA)  # keeps light grey visible
  ) +
  guides(fill = guide_colorbar(
    title = title,
    #direction = "horizontal",
    barheight = unit(130, "pt"),
    barwidth  = unit(6,  "pt"),
    ticks = TRUE
  ),
  colour = guide_legend(
    override.aes = list(size = 3.5),  # bigger points in legend only
    title.position = "top"
  )) + 
  geom_point(
    data = subset(dam_pts, damaged == 0),
    aes(x = x, y = y, colour = "Undamaged"),
    shape = 16, size = 0.5
  ) +
  geom_point(
    data = subset(dam_pts, damaged == 1),
    aes(x = x, y = y, colour = "Damaged"),
    shape = 16, size = 0.5
  ) +
  scale_colour_manual(name = "Damage status", values = c("Undamaged" = "lightgrey", "Damaged" = "red"))

#project(mandalay_dam, from=crs(mandalay_dam), to="EPSG:4326")
mandalay_dam_v <- project(vect(mandalay_dam),crs(ODD_cropped))
mandalay_dam_v$cond_1 <- as.numeric(mandalay_dam_v$damage_pct_0m > 0 & mandalay_dam_v$unknown_pct < .50)
mandalay_dam_v$cond_2 <- as.numeric(mandalay_dam_v$unknown_pct > .50)

mean(mandalay_dam_v$damage_pct_0m > 0)
mean(mandalay_dam_v$damaged > 0)


plot(ODD_cropped$nBuildings)
points(mandalay_dam_v)

dam_agg = data.frame(
  cell = numeric(),
  n_buildings = numeric(),
  damage_count = numeric(),
  unknown_count = numeric()
)

for (cell_num in overlapping_cells){
  
  xy <- xyFromCell(ODDy, cell=cell_num)
  
  # Cell resolution
  res_x <- res(ODDy)[1]
  res_y <- res(ODDy)[2]
  
  # Calculate extent: xmin, xmax, ymin, ymax
  cell_ext <- ext(
    xy[1] - res_x / 2,
    xy[1] + res_x / 2,
    xy[2] - res_y / 2,
    xy[2] + res_y / 2
  )
  
  cropped_builds = crop(mandalay_dam_v, cell_ext)
  dam_agg %<>% add_row(
    cell = cell_num,
    n_buildings = length(cropped_builds$cond_1),
    damage_count = sum(cropped_builds$cond_1, na.rm=T),
    unknown_count = sum(cropped_builds$cond_2, na.rm=T))
  
}

plot(colSums(sampled_full[overlapping_cells, 3, ])/sum(values(ODDy$nBuildings)[overlapping_cells], na.rm=T), 
     log10(colSums(sampled_full[,3,])),
     xlab='Prop Dam Mandalay', ylab='Total Buildings Damaged')
abline(v=sum(dam_agg$damage_count)/sum(dam_agg$n_buildings), col='red')
abline(h=log10(150000), col='red')

df_plot <- tibble::tibble(
  prop_mandalay = colSums(sampled_full[overlapping_cells, 3, ], na.rm = TRUE) /
    sum(values(ODDy$nBuildings)[overlapping_cells], na.rm = TRUE),
  total_log10   = colSums(sampled_full[, 3, ], na.rm = TRUE)
)

obs_prop   <- sum(dam_agg$damage_count, na.rm = TRUE) / sum(dam_agg$n_buildings, na.rm = TRUE)
true_total <- 150000

p <- ggplot(df_plot, aes(x = prop_mandalay, y = total_log10)) +
  geom_point(alpha = 0.7, colour = "#440154") +
  geom_point(aes(x=obs_prop, y=true_total), col='red', size = 3, shape = 4, stroke = 1.5) +
  labs(
    x = "Damage Proportion Mandalay",
    y = "Total Buildings Damaged"
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.margin  = margin(t = 10, r = 10, b = 20, l = 10)
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 2, base = 2),
    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
    labels = expression(phantom(1) * 0^phantom(0),phantom(1) * 10^phantom(1),10^2,10^3, 10^4, 10^5, 10^6, 10^7),
    minor_breaks = NULL
  ) + scale_x_continuous(expand=c(0,0))

p



cells_plot = c(15,6)
plot(t(sampled_full[overlapping_cells[cells_plot],3,]/values(ODDy$nBuildings)[overlapping_cells[cells_plot]]),
     xlab='Prop Dam Cell 1',
     ylab='Prop Dam Cell 2')
true_prop = dam_agg[cells_plot,'damage_count']/dam_agg[cells_plot,'n_buildings']
points(ifelse(is.na(true_prop[1]), 0, true_prop[1]), ifelse(is.na(true_prop[2]), 0, true_prop[2]), col='red', pch=19)

plot(colSums(sampled_full[overlapping_cells, 3, ])/sum(values(ODDy$nBuildings)[overlapping_cells], na.rm=T), 
     log10(colSums(sampled_full[,3,])),
     xlab='Prop Dam Mandalay', ylab='Total Buildings Damaged')
abline(v=sum(dam_agg$damage_count)/sum(dam_agg$n_buildings), col='red')
abline(h=log10(150000), col='red')

plot(colSums(sampled_full[overlapping_cells, 3, ])/sum(values(ODDy$nBuildings)[overlapping_cells], na.rm=T), 
     colSums(sampled_full[,1,]),
     xlab='Prop Dam Mandalay', ylab='Total Population Displacement')
abline(v=sum(dam_agg$damage_count)/sum(dam_agg$n_buildings), col='red')

plot(colSums(sampled_full[overlapping_cells, 3, ])/sum(values(ODDy$nBuildings)[overlapping_cells], na.rm=T), 
     log10(colSums(sampled_full[,2,])),
     xlab='Prop Dam Mandalay', ylab='Total Population Displacement')
abline(v=sum(dam_agg$damage_count)/sum(dam_agg$n_buildings), col='red')

ODDyWithImpact <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/EQ20250328MMR_-1_WithBuildings_WithImpactSummaries')
sampled_full <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp/EQ20250328MMR_-1_WithBuildings_fullSampledImpact')

#-----------------------------------

naypyidaw_dam = st_read("/home/manderso/Downloads/naypyidaw_myanmar_results_03313025_skysat.gpkg")
bbox_loc2 = project(ext(naypyidaw_dam), from=crs(naypyidaw_dam), to="EPSG:4326")

overlapping_cells2 = cells(ODDy, bbox_loc2)
plot(t(sampled_full[overlapping_cells2[c(1,2)],3,]))

ODD_cropped2 <- crop(ODDy, bbox_loc2)

#project(naypyidaw_dam, from=crs(naypyidaw_dam), to="EPSG:4326")
naypyidaw_dam_v <- project(vect(naypyidaw_dam),crs(ODD_cropped))
naypyidaw_dam_v$cond_1 <- as.numeric(naypyidaw_dam_v$damage_pct_0m > 0 & naypyidaw_dam_v$unknown_pct < .50)
naypyidaw_dam_v$cond_2 <- as.numeric(naypyidaw_dam_v$unknown_pct > .50)

dam_agg2 = data.frame(
  cell = numeric(),
  n_buildings = numeric(),
  damage_count = numeric(),
  unknown_count = numeric()
)

for (cell_num in overlapping_cells2){
  
  xy <- xyFromCell(ODDy, cell=cell_num)
  
  # Cell resolution
  res_x <- res(ODDy)[1]
  res_y <- res(ODDy)[2]
  
  # Calculate extent: xmin, xmax, ymin, ymax
  cell_ext <- ext(
    xy[1] - res_x / 2,
    xy[1] + res_x / 2,
    xy[2] - res_y / 2,
    xy[2] + res_y / 2
  )
  
  cropped_builds = crop(naypyidaw_dam_v, cell_ext)
  dam_agg2 %<>% add_row(
    cell = cell_num,
    n_buildings = length(cropped_builds$cond_1),
    damage_count = sum(cropped_builds$cond_1, na.rm=T),
    unknown_count = sum(cropped_builds$cond_2, na.rm=T))
  
}

plot(colSums(sampled_full[overlapping_cells, 3,])/sum(values(ODDy$nBuildings)[overlapping_cells]), 
     colSums(sampled_full[overlapping_cells2, 3,])/sum(values(ODDy$nBuildings)[overlapping_cells2]), 
     xlab='Mandalay Dam Prop', ylab='Naypyidaw Dam Prop')
points(sum(dam_agg$damage_count)/sum(dam_agg$n_buildings), sum(dam_agg2$damage_count)/sum(dam_agg2$n_buildings),
       col='red', pch=19)

plot(colSums(sampled_full[overlapping_cells2, 3, ])/sum(values(ODDy$nBuildings)[overlapping_cells2], na.rm=T), 
     log10(colSums(sampled_full[,3,])),
     xlab='Prop Dam Mandalay', ylab='Total Buildings Damaged')
abline(v=sum(dam_agg2$damage_count)/sum(dam_agg2$n_buildings), col='red')
abline(h=log10(150000), col='red')


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------------------------------------------------------------------
#=======================     ODDRIN vs PAGER Comparison  ============================
#------------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ODDRIN_vs_PAGER_hists <- function(ODDy, sampled_full, namer, impact_type = 'mortality'){
  #impact_type_samples = paste0(impact_type, '.s')
  
  plot_df <- as.data.frame(ODDy, xy=T, na.rm=F)
  #tot_impact = colSums(as.matrix(plot_df[,grep(impact_type_samples, names(plot_df))]))
  impact_i = which(c('displacement', 'mortality', 'buildDam') == impact_type)
  tot_impact = colSums(sampled_full[,impact_i,])
  folderin_haz <- paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/HAZARDobjects/')
  ufiles_haz <- na.omit(list.files(path=folderin_haz,pattern=Model$haz,recursive = T,ignore.case = T))
  
  
  file_match <- grep(namer,  ufiles_haz, value = TRUE)[1]
  if (length(file_match)==0){ stop()}
  HAZy <- readHAZ(paste0(folderin_haz, file_match ))
  which.max.mmi <- which.max(sapply(HAZy[2:length(HAZy)], function(x) max(values(x$mean), na.rm=T)))
  if (!all(HAZy$hazard_info$bbox == ODDy@hazinfo$bbox)){stop('Mismatched HAZARD and ODD objects')}
  
  binned_preds = data.frame(bin_lower=numeric(),
                            bin_upper=numeric(), 
                            ODDRIN_prob = numeric(),
                            PAGER_prob = numeric())
  
  bins = HAZy$hazard_info$alertfull[[which.max.mmi]]$bins
  
  if (is.null(bins) | length(bins)==0){
    bins = list(list(min=0, max=1, probability=0),
                list(min=1, max=10, probability=0),
                list(min=10, max=100, probability=0),
                list(min=100, max=1000, probability=0),
                list(min=1000, max=10000, probability=0),
                list(min=10000, max=100000, probability=0),
                list(min=100000, max=10000000, probability=0))
    prob_mult = 1
  } else {
    prob_mult <- 1/sum(unlist(lapply(bins, function(x) return(x$probability)))) # some sum to 1 and some to 100
  }
  
  for (j in 1:7){
    oddrin_preds <- tot_impact
    oddrin_prob = mean((oddrin_preds >= bins[[j]]$min) & (oddrin_preds <= bins[[j]]$max))
    
    binned_preds %<>% add_row(
      bin_lower = bins[[j]]$min,
      bin_upper = bins[[j]]$max,
      ODDRIN_prob =  mean((oddrin_preds >= bins[[j]]$min) & (oddrin_preds <= bins[[j]]$max)),
      PAGER_prob = bins[[j]]$probability * prob_mult
    )
  }
  binned_preds$bin_upper[1] = 0
  
  binned_preds$color <- cut(binned_preds$bin_lower, 
                            breaks = c(-Inf, 1, 100, 1000, Inf), 
                            labels = c("green", "yellow", "orange", "red"), 
                            right = FALSE)
  
  binned_preds$bin_label <- paste(binned_preds$bin_lower, "-", binned_preds$bin_upper)
  binned_preds$bin_expr = binned_preds$bin_label 

  binned_preds$bin_expr[1] <- parse(text = 'phantom("1"^1) * "0" * phantom("1"^1)')
  binned_preds$bin_expr[2] <- parse(text = 'phantom("1"^1) * "1 - 10" * phantom("1"^1)')
  binned_preds$bin_expr[3] = parse(text = "10 - 10^2")
  binned_preds$bin_expr[4] = parse(text = "10^2 - 10^3")
  binned_preds$bin_expr[5] = parse(text = "10^3 - 10^4")
  binned_preds$bin_expr[6] = parse(text = "10^4 - 10^5")
  binned_preds$bin_expr[7] = expression("10"^5 * "+")
  
  
  # Create the histogram
  # hist_ODDRIN <- ggplot(binned_preds, aes(x = factor(bin_label, levels = unique(bin_label)), y = ODDRIN_prob, fill = color)) +
  #   geom_bar(width = 1, stat = "identity", color = "grey") +
  #   scale_fill_identity() +
  #   scale_x_discrete(labels = binned_preds$bin_expr) +  # Use math expressions as labels
  #   labs(title = "PAGER Predicted Fatalities", y = "Probability", x="Fatalities") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
  #         plot.title = element_text(hjust = 0.5),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())
  # 
  # hist_PAGER <- ggplot(binned_preds, aes(x = factor(paste(bin_lower, ' - ', bin_upper)), y = PAGER_prob, fill = color)) +
  #   geom_bar(width=1,stat = "identity", color = "grey") +  # Black border for clarity
  #   scale_fill_identity() +  # Use pre-defined colors
  #   labs(title = "PAGER Predicted Fatalities", y = "Probability", x="Fatalities") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5),
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank()) 
  
  binned_long <- binned_preds %>%
    pivot_longer(cols = c(ODDRIN_prob, PAGER_prob),
                 names_to = "model",
                 values_to = "probability") %>%
    mutate(bin_label = paste(bin_lower, "-", bin_upper),
           model = ifelse(model == "ODDRIN_prob", "ODDRIN", "PAGER"))

  plot_compare = ggplot(binned_long, aes(x = factor(bin_label), y = probability, fill = model)) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.5, color = "black") +
    scale_fill_manual(values = c("ODDRIN" = "#1FA187", "PAGER" = "#440154")) +
    labs(x = "Predicted Fatalities", y = "Probability", fill = "Model") +
    theme_minimal(base_family = "serif") +
    scale_x_discrete(labels = binned_preds$bin_expr) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=2),
          plot.title = element_text(hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) #+ geom_vline(xintercept=5.3, color='#D73027', size=1)
  
  return(plot_compare)
  
  #plot_df[,grep(impact_type_samples, names(plot_df))] = t(apply(as.matrix(plot_df[, grep(impact_type_samples, names(plot_df)), drop = FALSE]), MARGIN=1, sort))
  
}

HAZy <- readHAZ('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/HAZARDobjects/EQ20250328MMR_-1_Sunday30th')
ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20250328MMR_-1')
sampled_full <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/sampled_full/EQ20250328MMR_-1_fullSampledImpact')

namer = 'EQ20250328MMR_-1'
hist1 = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist1 +  geom_vline(xintercept=5.3, color='#D73027', size=1) + xlab("Predicted Fatalities (MMR 2025-03-28)")




# ----------- SETUP ------------
source('RCode/ODDobj.R')
input$PosteriorFileLoc = 'IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1'
posterior_file_loc = 'IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1'

folder_write = 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/'

#---------- JPN 2024-01-01 ---------
#sampled_full <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/sampled_full/EQ20240101JPN_-1_fullSampledImpact')
ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240101JPN_-1')
namer = 'EQ20240101JPN_-1'

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200, 
                                        s_finish = 250)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

hist_JPN = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist_JPN +  geom_vline(xintercept=4, color='#D73027', size=1) + xlab("Predicted Fatalities (JPN 2024-01-01)")



# ------ IDN 2024-01-01 --------
# 
# input<-list(
#   sdate=as.Date("2024-01-08"), 
#   fdate=as.Date("2024-01-08"), 
#   iso3="PHL", 
#   datadir=dir,
#   plotdir="Plots/"
# )

input<-list(USGSid="us6000m2jp",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

input%<>%append(list(Model=Model, 
                     PosteriorFileLoc=posterior_file_loc,
                     Method=AlgoParams))

# options(timeout = 200)
# ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# ODDy = ODDy_with_namer$ODDy
# namer = ODDy_with_namer$namer

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240109PHL_-1')
namer <- 'EQ20240109PHL_-1'

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200,
                                        s_finish = 250)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

hist_PHL = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist_PHL+  geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (PHL 2024-01-08)")

# #---------------------------------------------------------------------
# # USA 2024-01-13
# #---------------------------------------------------------------------
# 
# input<-list(USGSid="ok2024awfn",
#             datadir=dir, # Location of the main folder to access the data 
#             plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
# )
# 
# input%<>%append(list(Model=Model, 
#                      PosteriorFileLoc=posterior_file_loc,
#                      Method=AlgoParams))
# 
# # options(timeout = 500)
# # ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# # ODDy = ODDy_with_namer$ODDy
# # namer = ODDy_with_namer$namer
# 
# ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240113USA_-1')
# namer = 'EQ20240113USA_-1'
# 
# # Sample impact:
# ODDy_with_sampled = PosteriorImpactPred(ODDy, 
#                                         AlgoResultsFilename=input$PosteriorFileLoc,
#                                         folder_write=folder_write,
#                                         namer=namer,
#                                         N_samples = 200)
# ODDy <- ODDy_with_sampled$ODDyWithImpact
# sampled_full <- ODDy_with_sampled$sampled_full
# 
# hist_USA = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
# hist_USA +  geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (USA 2024-01-13)")

#---------------------------------------------------------------------
# COL 2024-01-19
#---------------------------------------------------------------------

input<-list(USGSid="us6000m4n2",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

input%<>%append(list(Model=Model, 
                     PosteriorFileLoc=posterior_file_loc,
                     Method=AlgoParams))

# options(timeout = 200)
# ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# ODDy = ODDy_with_namer$ODDy
# namer = ODDy_with_namer$namer

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240119COL_-1')
namer <- 'EQ20240119COL_-1'

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200,
                                        s_finish = 250)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

hist_COL = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist_COL +  geom_vline(xintercept=1.7, color='#D73027', size=1) + xlab("Predicted Fatalities (COL 2024-01-19)")

# #---------------------------------------------------------------------
# # USA 2024-01-19
# #---------------------------------------------------------------------
# 
# input<-list(USGSid="ak024vrx24o",
#             datadir=dir, # Location of the main folder to access the data 
#             plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
# )
# 
# input%<>%append(list(Model=Model, 
#                      PosteriorFileLoc=posterior_file_loc,
#                      Method=AlgoParams))
# 
# # options(timeout = 200)
# # ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# # ODDy = ODDy_with_namer$ODDy
# # namer = ODDy_with_namer$namer
# 
# ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240119USA_-1')
# namer = 'EQ20240119USA_-1'
# 
# # Sample impact:
# ODDy_with_sampled = PosteriorImpactPred(ODDy, 
#                                         AlgoResultsFilename=input$PosteriorFileLoc,
#                                         folder_write=folder_write,
#                                         namer=namer,
#                                         N_samples = 200)
# ODDy <- ODDy_with_sampled$ODDyWithImpact
# sampled_full <- ODDy_with_sampled$sampled_full
# 
# hist_ALA = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
# hist_ALA +  geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (USA 2024-01-19)")

# #---------------------------------------------------------------------
# # CHN 2024-01-22
# #---------------------------------------------------------------------
# 
input<-list(USGSid="us7000lsze",
            datadir=dir, # Location of the main folder to access the data
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

input%<>%append(list(Model=Model,
                     PosteriorFileLoc=posterior_file_loc,
                     Method=AlgoParams))

# options(timeout = 200)
# ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# ODDy = ODDy_with_namer$ODDy
# namer = ODDy_with_namer$namer

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240123CHN_-1')
namer = 'EQ20240123CHN_-1'

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy,
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200,
                                        s_finish = 250)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

hist_CHN = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist_CHN +  geom_vline(xintercept=1.9, color='#D73027', size=1) + xlab("Predicted Fatalities (CHN 2024-01-22)")

#---------------------------------------------------------------------
# RUS 2024-01-15
#---------------------------------------------------------------------

input<-list(USGSid="us6000m3sa",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

input%<>%append(list(Model=Model, 
                     PosteriorFileLoc=posterior_file_loc,
                     Method=AlgoParams))

# options(timeout = 500)
# ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# ODDy = ODDy_with_namer$ODDy
# namer = ODDy_with_namer$namer

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240123CHN_-1')
namer = 'EQ20240123CHN_-1'

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200,
                                        s_finish = 250)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

hist_RUS = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist_RUS +  geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (RUS 2024-01-15)")

# #---------------------------------------------------------------------
# #AFG 2023-10-07
# #---------------------------------------------------------------------
# 
# input<-list(USGSid="us6000ldpm",
#             datadir=dir, # Location of the main folder to access the data 
#             plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
# )
# 
# input%<>%append(list(Model=Model, 
#                      PosteriorFileLoc=posterior_file_loc,
#                      Method=AlgoParams))
# 
# # options(timeout = 500)
# # ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# # ODDy = ODDy_with_namer$ODDy
# # namer = ODDy_with_namer$namer
# 
# ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20231007AFG_-1')
# namer = 'EQ20231007AFG_-1'
# 
# # Sample impact:
# ODDy_with_sampled = PosteriorImpactPred(ODDy, 
#                                         AlgoResultsFilename=input$PosteriorFileLoc,
#                                         folder_write=folder_write,
#                                         namer=namer,
#                                         N_samples = 200)
# ODDy <- ODDy_with_sampled$ODDyWithImpact
# sampled_full <- ODDy_with_sampled$sampled_full
# 
# hist_AFG = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
# hist_AFG +  geom_vline(xintercept=4.7, color='#D73027', size=1) + xlab("Predicted Fatalities (AFG 2023-10-07)")


#---------------------------------------------------------------------
# TWN 2024-04-02
#---------------------------------------------------------------------

input<-list(USGSid="us7000m9g4",
            datadir=dir, # Location of the main folder to access the data
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

input%<>%append(list(Model=Model,
                     PosteriorFileLoc=posterior_file_loc,
                     Method=AlgoParams))
# 
# options(timeout = 500)
# ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# ODDy = ODDy_with_namer$ODDy
# namer = ODDy_with_namer$namer

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20240403CHN_-1')
namer = 'EQ20240403CHN_-1'

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200,
                                        s_finish = 250)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

hist_TWN = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist_TWN +  geom_vline(xintercept=2.5+9/90+1/20, color='#D73027', size=1) + xlab("Predicted Fatalities (TWN 2024-04-03)")


# #---------------------------------------------------------------------
# #PNG 2018-02-26
# #---------------------------------------------------------------------
# 
# input<-list(USGSid="us2000d7q6",
#             datadir=dir, # Location of the main folder to access the data 
#             plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
# )
# 
# input%<>%append(list(Model=Model, 
#                      PosteriorFileLoc=posterior_file_loc,
#                      Method=AlgoParams))
# 
# options(timeout = 500)
# ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# ODDy = ODDy_with_namer$ODDy
# namer = ODDy_with_namer$namer
# 
# ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20180226PNG_-1')
# namer = 'EQ20180226PNG_-1'
# 
# # Sample impact:
# ODDy_with_sampled = PosteriorImpactPred(ODDy, 
#                                         AlgoResultsFilename=input$PosteriorFileLoc,
#                                         folder_write=folder_write,
#                                         namer=namer,
#                                         N_samples = 200)
# ODDy <- ODDy_with_sampled$ODDyWithImpact
# sampled_full <- ODDy_with_sampled$sampled_full
# 
# hist_PNG = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
# hist_PNG +  geom_vline(xintercept=4.7, color='#D73027', size=1) + xlab("Predicted Fatalities (PNG 2018-02-26)")

#---------------------------------------------------------------------
#MMR 2025-03-28
#---------------------------------------------------------------------

input<-list(USGSid="us7000pn9s",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

input%<>%append(list(Model=Model, 
                     PosteriorFileLoc=posterior_file_loc,
                     Method=AlgoParams))

# options(timeout = 500)
# ODDy_with_namer <- prepareODD(dir, input, getGADMregions=F, folder_write=folder_write)
# ODDy = ODDy_with_namer$ODDy
# namer = ODDy_with_namer$namer

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20250328MMR_-1')
namer = 'EQ20250328MMR_-1'

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200,
                                        s_finish = 250)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

hist_MMR = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist_MMR +  geom_vline(xintercept=4.5, color='#D73027', size=1) + xlab("Predicted Fatalities (MMR 2025-03-28)")


# -------------------------- PLOT: ---------------------------

legend <- get_legend(hist_MMR + theme(legend.position = "right"))

# Combine plots without legends
plots <- plot_grid(
  hist_JPN + theme(legend.position = "none", 
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                   plot.margin = margin(t=2, l=8, b=8, r=4)) + 
                   scale_y_continuous(expand=expansion(mult = c(0, 0.05))) +
                   geom_vline(xintercept=3.5 + 270/1000, color='#D73027', size=1) + xlab("Predicted Fatalities (JPN 2024-01-01)"),
  hist_PHL + theme(legend.position = "none", 
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                   plot.margin = margin(t=2, l=8, b=8, r=4)) + 
    scale_y_continuous(expand=expansion(mult = c(0, 0.05))) +
    geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (PHL 2024-01-08)"),
  hist_RUS + theme(legend.position = "none", 
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                   plot.margin = margin(t=2, l=8, b=8, r=4)) + 
    scale_y_continuous(expand=expansion(mult = c(0, 0.05))) +
    geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (RUS 2024-01-15)"),
  hist_CHN + theme(legend.position = "none", 
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                   plot.margin = margin(t=2, l=8, b=8, r=4)) + 
    scale_y_continuous(expand=expansion(mult = c(0, 0.05))) +
    geom_vline(xintercept=1.9, color='#D73027', size=1) + xlab("Predicted Fatalities (CHN 2024-01-22)"),
  #hist_COL + theme(legend.position = "none") +  geom_vline(xintercept=1.5 + 1.5/10, color='#D73027', size=1) + xlab("Predicted Fatalities (COL 2024-01-19)"),
  #hist_ALA + theme(legend.position = "none") +  geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (USA 2024-01-19)"),
  #hist_CHN + theme(legend.position = "none") +  geom_vline(xintercept=1.5 + 3/10, color='#D73027', size=1) + xlab("Predicted Fatalities (CHN 2024-01-22)"),
  hist_TWN + theme(legend.position = "none", 
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                   plot.margin = margin(t=2, l=8, b=8, r=4)) + 
    scale_y_continuous(expand=expansion(mult = c(0, 0.05))) +
    geom_vline(xintercept=2.5 + 9/90+1/10, color='#D73027', size=1) + xlab("Predicted Fatalities (TWN 2024-04-03)"),
  hist_MMR + theme(legend.position = "none", 
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                   plot.margin = margin(t=2, l=8, b=8, r=4)) + 
    scale_y_continuous(expand=expansion(mult = c(0, 0.05))) +
    geom_vline(xintercept=4.5 + 2804/9000+1/10, color='#D73027', size=1) + xlab("Predicted Fatalities (MMR 2025-03-28)"),
  ncol = 2,
  labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'),
  label_fontfamily = "Times New Roman",
  label_fontface = "plain",
  label_x = -0.02,  # Move label to the left
  label_y = 1.01
)

# Combine with legend on the right
final_plot <- plot_grid(
  plots,
  legend,
  rel_widths = c(1, 0.15),
  nrow = 1
)

final_plot
#ODDRINvsPAGER.pdf, 16 x 5
#ODDRINvsPAGER.pdf, 9.5 x 7




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------------------------------------------------------------------------
#========================     Misc Performance Evaluations  =========================
#------------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# ------- Mexico -----------------

source('RCode/Model Variations/Model_MVRrank.R')
source('RCode/ODDobj.R')

input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
                     PosteriorFileLoc='IIDIPUS_Results/HPC/mcmc_2025-06-17_151849.678089__Apr25Agg_LogNormal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_withHLPriors_flexHLP',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_014135.935892__Apr25Agg_Lognormal_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
                     Method=AlgoParams))

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Apr25Agg/ODDobjects/Train/EQ20170919MEX_68')
namer = 'EQ20170919MEX_68'
  
# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 100)



ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

names(mortality_true) = paste0(names(mortality_true), ', Mexico')

observed_df <- data.frame(
  polygon_name = names(mortality_true),
  observed = as.vector(mortality_true),
  stringsAsFactors = FALSE
)

impact_mort = ODDy@impact %>% filter(impact=='mortality')

impact_mort = impact_mort[order(impact_mort$observed, decreasing=T),]
impact_mort_shortened_names = impact_mort
impact_mort_shortened_names$polygon_name <- sub(",.*", "", impact_mort_shortened_names$polygon_name)
pairplot_regions(impact_mort_shortened_names)


subnat_predDist = plot_predictive_dist_subnat(impact_mort)


# ------- Japan -----------------

source('RCode/Model Variations/Model_MVRrank.R')
source('RCode/ODDobj.R')

input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
                     PosteriorFileLoc='IIDIPUS_Results/HPC/mcmc_2025-06-17_151849.678089__Apr25Agg_LogNormal_ESplus.05MVRrank_RFwithTot_range.5_kappa.5_withHLPriors_flexHLP',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_014135.935892__Apr25Agg_Lognormal_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
                     Method=AlgoParams))

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/EQ20240101JPN_-1_WithImpact200')
namer = 'EQ20240101JPN_-1'

ODDy@impact = ODDy@impact[, -grep('sampled', names(ODDy@impact))]
ODDy@impact$polygon_name = NULL

# Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy, 
                                        AlgoResultsFilename=input$PosteriorFileLoc,
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 100)



ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

impact_mort = ODDy@impact %>% filter(impact=='mortality')

impact_mort = impact_mort[order(impact_mort$observed, decreasing=T),]
impact_mort_shortened_names = impact_mort
impact_mort_shortened_names$polygon_name <- sub(",.*", "", impact_mort_shortened_names$polygon_name)
pairplot_regions(impact_mort_shortened_names)


subnat_predDist = plot_predictive_dist_subnat(impact_mort)






#------------------------------------------------------------------------------------
#---------------------------------Haiti 2021-08-14 ----------------------------------
#------------------------------------------------------------------------------------

#
input<-list(USGSid="us6000f65h",
            datadir=dir, # Location of the main folder to access the data
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)
# 
input%<>%append(list(Model=Model,
                     PosteriorFileLoc='IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1',
                     Method=AlgoParams))
# 
# # options(timeout = 500)
ODDy_with_namer <- prepareODD(dir, input, getGADMregions=T, folder_write=folder_write)
ODDy = ODDy_with_namer$ODDy
namer = ODDy_with_namer$namer
# 
#ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20210814HTI_-1')
#namer = 'EQ20210814HTI_164'

folder_write = 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/'
# 
# # Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy,
                                        AlgoResultsFilename='IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1',
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 500)

sampled_impact = ODDy_with_sampled$sampled_full

mean(colSums(sampled_impact[,1,])>100000)

impact_full = ODDy_with_sampled$ODDyWithImpact@impact %>% as.data.frame()
impact_full$sdate = NULL
impact_mort = impact_full %>% dplyr::filter(impact=='mortality')
#impact_mort$observed = NULL

ODDy2 <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects/Test/EQ20210814HTI_164')

ODDy2@impact$polygon_name = ''
for (i in 1:NROW(ODDy2@impact)){
  ODDy2@impact$polygon_name[i] = ODDy2@polygons[[ODDy2@impact$polygon[i]]]$name
}

mortality_true = c("Grand'Anse, Haiti" = 256, 
                   #'Haiti' = 219745, 
                   'Nippes, Haiti' = 167,
                   'Nord-Ouest, Haiti' = 2, 
                   'Sud, Haiti' =  2152,
                   'Centre, Haiti' = 0, 
                   'Sud-Est, Haiti' = 0,
                   'Ouest, Haiti'=0, 
                   'Nord-Est, Haiti' = 0,
                   'Nord, Haiti'=0,
                   "L'Artibonite, Haiti"=1)

#names(mortality_true) = paste0(names(mortality_true), ', Japan')

observed_df <- data.frame(
  polygon_name = names(mortality_true),
  observed = as.vector(mortality_true),
  stringsAsFactors = FALSE
)

#impact_mort = ODDy@impact %>% filter(impact=='mortality')
#impact_mort$observed = NULL


impact_mort_merged = impact_mort %>% merge(observed_df, by='polygon_name')
#impact_mort_merged$observed[which(is.na(impact_mort_merged$observed))] = 0
impact_mort_merged$observed = impact_mort_merged$observed.y
impact_mort_merged$observed.y = NULL; impact_mort_merged$observed.x = NULL

impact_mort_merged %<>% rbind(impact_full %>% filter(polygon_name == 'Haiti' & impact=='displacement'))
impact_mort_merged$polygon_name[nrow(impact_mort_merged)] = 'Haiti Total'
impact_mort_merged$observed[nrow(impact_mort_merged)] = 219745
impact_mort_merged %<>% rbind(impact_full %>% filter(polygon_name == 'Haiti' & impact=='mortality'))
impact_mort_merged$polygon_name[nrow(impact_mort_merged)] = 'Haiti Total'
impact_mort_merged$observed[nrow(impact_mort_merged)] = 2577



impact_mort_merged = impact_mort_merged[order(impact_mort_merged$observed, decreasing=T),]
impact_mort_shortened_names = impact_mort_merged
impact_mort_shortened_names$polygon_name <- paste0(paste0(
  toupper(substring(impact_mort_shortened_names$impact, 1, 1)),
  substring(impact_mort_shortened_names$impact, 2)
), ':\n',sub(",.*", "", impact_mort_shortened_names$polygon_name))
# impact_mort_shortened_names %<>% filter(polygon_name != 'Hodatsushimizu')
impact_mort_shortened_names = impact_mort_shortened_names[1:6,]
pairplot_regions(impact_mort_shortened_names, cross_size = 1.5)
# HTI_pairplot.pdf, 7 x 7



plot(log10(colSums(sampled_impact[,1,])), log10(colSums(sampled_impact[,2,])))
points(log10(219745), log10(2248), col='red', pch=19)

ODDy[['displacement.median']] = apply(sampled_impact[,1,],1, quantile, 0.5)
plot(ODDy[['displacement.median']])



plot_GADM_impact_polygons_pixellated_HTI <- function(ODDy, polygons_list, var, gadm_level = 2, impact_type = 'displacement', summary = 'median', plot_bbox = NULL, max_break = 500){
  #GADM_level 0 = national
  #GADM_level 1 = admin level 1
  #GADM_level 2 = admin level 2
  #summary = mean, median, q5, q10, q20, q80, q90, q95
  
  if (is.null(plot_bbox)){plot_bbox = c(-75, -72, 17.95, 20)}
  
  
  plot_df <- as.data.frame(ODDy, xy=T, na.rm=F)
  #plot_df[,grep(impact_type_samples, names(plot_df))] = t(apply(as.matrix(plot_df[, grep(impact_type_samples, names(plot_df)), drop = FALSE]), MARGIN=1, sort))

  #plot_df[which(is.na(plot_df$ISO3C)), var] <- NA
  names(plot_df)[which(names(plot_df)=='x')] = 'Longitude'
  names(plot_df)[which(names(plot_df)=='y')] = 'Latitude'
  
  bbox <- matrix(ODDy@hazinfo$bbox, nrow=2)
  #gadm_iso <- getData("GADM", country="NZL", level=2)
  iso3_unique <- 'HTI'
  gadm_iso <- as(geodata::gadm(country=iso3_unique[1], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial')
  if (length(iso3_unique) > 1){
    for (i in 2:length(iso3_unique)){
      gadm_iso %<>% rbind(as(geodata::gadm(country=iso3_unique[i], level=gadm_level, path=paste0(dir, 'Demography_Data/GADM/')), 'Spatial'))
    }
  }
  
  
  #gadm_iso <- gSimplify(gadm_iso, 0.01)
  gadm_iso <- intersect(gadm_iso, bbox)
  gadm_map <- fortify(gadm_iso)
  
  gadm_map_cropped <- gadm_map %>%
    filter(long >= plot_bbox[1], long <= plot_bbox[2],
           lat  >= plot_bbox[3], lat  <= plot_bbox[4])
  
  p <- ggplot() + xlab("Longitude") + ylab("Latitude")#+ theme(legend.position = "none")
  
  p <- p + geom_map(map=gadm_map_cropped, data=gadm_map_cropped, aes(map_id=id, group=id)) + xlim(bbox[1,1],bbox[1,2]) + ylim(bbox[2,1], bbox[2,2])
  p <- p + coord_map() + geom_polygon(data=gadm_map_cropped, aes(x=long, y=lat, group=group), fill='white', color='black')
  
  if (impact_type == 'mortality'){
    breaks = c(-Inf, 1, 3, 5, 10, 20, 50, 100, 500, Inf)
    if (max_break < 500 & max_break > 100){
      breaks = c(-Inf, 1, 3, 5, 10, 20, 50, 100, Inf)
    }
    if (max_break < 100){
      breaks = c(-Inf, 1, 3, 5, 10, 20, 50, Inf)
    }
    if (max_break < 50){
      breaks = c(-Inf, 1, 5, 10, 20, 30, Inf)
    }
    if (max_break < 25){
      breaks = c(-Inf, 1, 5, 10, 15, 20, Inf)
    }
    if (max_break < 17){
      breaks = c(-Inf, 1, 4, 7, 10, 13, 16, Inf)
    }
    
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
    title = ifelse(var=='mortality.q10', 'Mortality 10% Quantile', ifelse(var=='mortality.q90', 'Mortality 90% Quantile',
                                                                                'Median Mortality'))
    legend_title = 'Mortality'
  } else if (impact_type == 'Population'){
    breaks = c(-Inf, 1, 100, 1000, 5000, 10000, 50000, 100000, Inf)
    # if (max_break < 50000 & max_break > 10000){
    #   breaks = c(-Inf, 1, 50, 100, 500, 1000, 10000, 50000, Inf)
    # }
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
    
    title = 'Population'
    legend_title = 'Population'
  } else {
    breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, 5000, 10000, Inf)
    if (max_break < 50000 & max_break > 10000){
      breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, 10000, 50000, Inf)
    }
    plot_df$bin_category <- cut(
      plot_df[[var]], 
      breaks = breaks, 
      labels = c(0, paste0(breaks[2:(length(breaks)-2)], ' - ', breaks[3:(length(breaks)-1)]-1), paste0(breaks[length(breaks)-1], '+')), #c("0", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-500","500+"),
      right = FALSE
    )
    
    title = ifelse(var=='displacement.q10', 'Displacement 10% Quantile', ifelse(var=='displacement.q90', 'Displacement 90% Quantile',
                                                                                'Median Displacement'))
    legend_title = 'Displacement'
  }
  
  bin_colors <- viridis::viridis(length(levels(plot_df$bin_category))) # c(viridis::viridis(length(levels(plot_df$bin_category))-2), 'orange', 'red')

  # Modify the plot
  map_plot <- p + 
    geom_raster(data = plot_df, aes(x = Longitude, y = Latitude, fill = bin_category), alpha = 0.75) + 
    coord_equal() +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal() + 
    labs(title=title) +
    theme(
      axis.title = element_text(family = "Liberation Serif", size = 12),  
      legend.text = element_text(family = "Liberation Serif", size = 11),
      legend.title = element_text(family = "Liberation Serif", size = 12),
      plot.title = element_text(
        family = "Liberation Serif",
        hjust = 0,           
        size = 13           
      ),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.margin = margin(t=2, l=0, b=8, r=0)
    ) + 
    scale_fill_manual(name = legend_title,#paste0(toupper(substr(impact_type, 1, 1)), substr(impact_type, 2, nchar(impact_type)),' ', quantile,'% Quantile'), 
                      values = bin_colors, na.value = "transparent") + # Discrete color scale 
    xlim(plot_bbox[c(1,2)]) + ylim(plot_bbox[c(3,4)])
  
  
  # Define the bounding box
  full_extent <- st_as_sfc(st_bbox(c( xmin = plot_bbox[1], xmax = plot_bbox[2], ymin = plot_bbox[3], ymax = plot_bbox[4] ), crs = st_crs(gadm_iso)))
  
  mask_poly <- st_difference(full_extent, st_union( st_as_sf(gadm_iso)))
  
  polygon_list_plot <- list()
  for (i in 1:length(polygons_list[[2]])){
    if (is.null(polygons_list[[2]][[i]]$sf_polygon$NAME_2)) next #polygons_list$polygons_list[[i]]$sf_polygon$NAME_2 = ''
    polygon_list_plot[[i]] = st_as_sf(polygons_list[[2]][[i]]$sf_polygon)
  }
  combined_polygons <- do.call(rbind, polygon_list_plot)
  
  p = map_plot + geom_sf(data=mask_poly,fill='lightgrey', color=NA) + 
    geom_polygon(data=gadm_map, aes(x=long, y=lat, group=group), fill='white', alpha=0,color='black') +
    geom_contour(data = plot_df, aes(Longitude, Latitude, z = hazMean1, colour = ..level..),
                 alpha = 0.7, lwd = 0.8) +
    scale_color_gradientn(colors = c("transparent", "#fc9272", "#ef3b2c")) +
    labs(colour = "Hazard Intensity") +
    #geom_sf_text(data=combined_polygons, aes(label = NAME_2), size = 3, color = "white", fontface = "bold") + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + coord_sf(
      xlim = plot_bbox[c(1,2)],
      ylim = plot_bbox[c(3,4)],
      expand = FALSE
    )
  
  p
  
  return(p)
}

ODDy[['displacement.q10']] = apply(sampled_impact[,1,],1, quantile, 0.1)
ODDy[['displacement.median']] = apply(sampled_impact[,1,],1, quantile, 0.5)
ODDy[['displacement.q90']] = apply(sampled_impact[,1,],1, quantile, 0.9)

polygons_list <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp_polygons/EQ20210814HTI_-1polygons')


p10 <- plot_GADM_impact_polygons_pixellated_HTI(ODDy, polygons_list, var = 'displacement.q10')
p50 <- plot_GADM_impact_polygons_pixellated_HTI(ODDy, polygons_list, var = 'displacement.median')
p90 <- plot_GADM_impact_polygons_pixellated_HTI(ODDy, polygons_list, var = 'displacement.q90')

legend_contour <- get_legend(
  p90 +
    #guides(fill = "none") +  # remove fill legend
    #guides(fill = guide_colorbar(direction = "horizontal", title.position = "top"),
    #       colour = guide_colorbar(direction = "horizontal", title.position = "top")) +
    theme(legend.position = "right",
          legend.box = "horizontal")
)



plot_grid(p10 + theme(legend.position="none"), 
          p50  + theme(legend.position="none"), 
          p90  + theme(legend.position="none"), plot_grid( NULL, legend_contour,rel_widths=c(0.2,1)), nrow=2,
          labels=c('(a)', '(b)', '(c)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain", 
          label_x = 0,
          label_y = 1.02
          )
#HTI_Disp_Quantiles, 10 x 7


#--------- Dashboard example for slides ----------------

plot_GADM_impact_polygons_pixellated_HTI(ODDy, polygons_list, var = 'Population', impact_type = 'Population', plot_bbox = c(-74.6, -72, 17.95, 20))
plot_GADM_impact_polygons_pixellated_HTI(ODDy, polygons_list, var= 'displacement.median', plot_bbox = c(-74.6, -72, 17.95, 20))


ODDy[['mortality.q10']] = apply(sampled_impact[,2,],1, quantile, 0.1)
ODDy[['mortality.median']] = apply(sampled_impact[,2,],1, quantile, 0.5)
ODDy[['mortality.q90']] = apply(sampled_impact[,2,],1, quantile, 0.9)

polygons_list <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp_polygons/EQ20210814HTI_-1polygons')


p_mort10 <- plot_GADM_impact_polygons_pixellated_HTI(ODDy, polygons_list, impact_type = 'mortality', var = 'mortality.q10')
p_mort50 <- plot_GADM_impact_polygons_pixellated_HTI(ODDy, polygons_list, impact_type = 'mortality', var = 'mortality.median')
p_mort90 <- plot_GADM_impact_polygons_pixellated_HTI(ODDy, polygons_list, impact_type = 'mortality', var = 'mortality.q90')

legend_contour_mort <- get_legend(
  p_mort90 +
    #guides(fill = "none") +  # remove fill legend
    #guides(fill = guide_colorbar(direction = "horizontal", title.position = "top"),
    #       colour = guide_colorbar(direction = "horizontal", title.position = "top")) +
    theme(legend.position = "right",
          legend.box = "horizontal")
)



plot_grid(plot_grid(p_mort10 + theme(legend.position="none"), 
          p_mort50  + theme(legend.position="none"), 
          p_mort90  + theme(legend.position="none"), plot_grid( NULL, legend_contour_mort,rel_widths=c(0.2,1)),nrow=1, rel_widths=c(1,1,1,0.5)),
          plot_grid(p10 + theme(legend.position="none"), 
          p50  + theme(legend.position="none"), 
          p90  + theme(legend.position="none"), plot_grid( NULL, legend_contour,rel_widths=c(0.3,2)),nrow=1, rel_widths=c(1,1,1,0.5)),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain", 
          label_x = 0,
          label_y = 1.02, nrow=2
)

p50

ODDRIN_hist2(sampled_full, impact_type='mortality') + ggtitle('')
ODDRIN_hist2(sampled_full, impact_type='displacement') + ggtitle('')

#-------------------------



polygons_list <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp_polygons/EQ20210814HTI_-1polygons')

#plot_GADM_impact_polygons_pixellated(ODDy, polygons_list, summary='median')



impact_agg = ODDy_with_sampled$ODDyWithImpact@impact
impact_agg = impact_agg[order(impact_agg$observed, decreasing=T),]

i_1 = 5; i_2 = 3;
plot(log10(as.numeric(impact_agg[i_1, 11:110])), log10(as.numeric(impact_agg[i_2, 11:110])))
points(log10(as.numeric(impact_agg[i_1, 5])),log10(as.numeric(impact_agg[i_2, 5])), col='red', pch=19)


admin3_shp = read_sf('/home/manderso/Downloads/hti_adm_cnigs_20181129/hti_admbndl_admALL_cnigs_itos_20181129.shp')
plot(ODDy$Population, xlim=c(-74, -73.5), ylim=c(18, 18.4))
lines(admin3_shp, col='red')

coords = crds(ODDy, na.rm=F)
coords %<>% as.data.frame()

#---------------------------- Streamlined Haiti Displacement -------------------------------

library(ggrepel)    # for nicer labels
# library(patchwork)  # only if you want to assemble separate plots

# --- 1) Define the 9 sections with their bounding boxes -----------------------
# NOTE: Replace any bbox if you want to tweak selection windows.
sections <- tibble::tribble(
  ~name,                  ~xmin,    ~xmax,    ~ymin,   ~ymax,
  "Bourdet",              -73.82,   -73.70,   18.15,   18.27,
  "Fond Rouge Torbeck",   -74.15,   -74.10,   18.60,   18.65,
  "Fond Rouge Dahere",    -74.15,   -74.10,   18.60,   18.65,  # adjust if needed
  "Bariadelle",           -74.45,   -74.40,   18.50,   18.55,
  "Frangipane",           -73.10,   -73.00,   18.20,   18.30,
  "Champlois",            -73.90,   -73.85,   18.33,   18.35,
  "Levy Marsan",          -73.90,   -73.85,   18.26,   18.31,
  "Brodequin",            -73.43,   -73.40,   18.26,   18.31,
  "Barbois",              -73.90,   -73.85,   18.12,   18.15
)

# --- 2) Find cell indices per section -----------------------------------------
# Assumes you already have:
#   coords <- as.data.frame(crds(ODDy, na.rm = FALSE)) with columns x, y
#   ODDy_df <- as.data.frame(ODDy, na.rm = FALSE); ODDy_df$x <- coords$x; ODDy_df$y <- coords$y
#   admin3_shp is sf
#   sampled_impact is [cells x ? x samples], where [,1,] selects event 1 across samples

cell_index_for_bbox <- function(xmin, xmax, ymin, ymax, coords_df) {
  which(coords_df$x > xmin & coords_df$x < xmax &
          coords_df$y > ymin & coords_df$y < ymax)
}

sections <- sections %>%
  mutate(
    cell_idx = pmap(
      list(xmin, xmax, ymin, ymax),
      ~ cell_index_for_bbox(..1, ..2, ..3, ..4, coords)
    )
  )

# Label positions: use mean of the selected cells (fallback to bbox center if empty)
sections <- sections %>%
  dplyr::mutate(
    lab_x = map_dbl(cell_idx, ~ if(length(.x)) mean(coords$x[.x], na.rm = TRUE) else (xmin + xmax)/2),
    lab_y = map_dbl(cell_idx, ~ if(length(.x)) mean(coords$y[.x], na.rm = TRUE) else (ymin + ymax)/2)
  )

sections[which(sections$name=='Fond Rouge Torbeck'), c('lab_x', 'lab_y')] = list(-74.14, 18.63)
sections[which(sections$name=='Fond Rouge Dahere'), c('lab_x', 'lab_y')] = list(-74.14, 18.61)
sections[which(sections$name=='Bariadelle'), c('lab_x', 'lab_y')] = list(-74.43, 18.52)
sections[which(sections$name=='Levy Marsan'), c('lab_x', 'lab_y')] = list(-73.85, 18.3)
sections[which(sections$name=='Champlois'), c('lab_x', 'lab_y')] = list(-73.85, 18.36)
sections[which(sections$name=='Barbois'), c('lab_x', 'lab_y')] = list(-73.9, 18.10)
sections[which(sections$name=='Bourdet'), c('lab_x', 'lab_y')] = list(-73.78, 18.2)

ODDy_df = as.data.frame(ODDy, na.rm=F)
ODDy_df$x = crds(ODDy, na.rm=F)[,1]
ODDy_df$y = crds(ODDy, na.rm=F)[,2]
#ODDy_df$Displacement = apply(sampled_impact[,1,],1, quantile, 0.5)

# --- 3) Full map with points + labels -----------------------------------------
map_full <- ggplot() +
  geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
  geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.5) +
  coord_sf(xlim = c(-74.55, -73.00), ylim = c(17.95, 18.72), expand = FALSE) +
  scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor, na.value = "darkgrey") +
  geom_point(data = sections, aes(x = lab_x, y = lab_y), colour = "white", size = 2) +
  ggrepel::geom_text_repel(
    data = sections,
    aes(x = lab_x, y = lab_y, label = name),
    colour = "white", fontface = "bold", box.padding = 0.25, min.segment.length = 1
  ) +
  labs(x = NULL, y = NULL, fill = "Population") +
  theme_minimal() +
  theme(
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    legend.title   = element_text(family = "Times New Roman"),
    legend.text    = element_text(family = "Times New Roman"),
    legend.position      = c(0.02, 0.03),     # bottom-left inside panel (x,y in [0,1])
    legend.justification = c(0, 0),           # anchor the legend’s bottom-left corner
    legend.background    = element_rect(fill = alpha("white", 0.6), colour = "black"),
    legend.key           = element_rect(fill = "white", colour = NA),
    legend.margin        = margin(6, 6, 6, 6)
  )

map_full
# print(map_full)

# --- 4) Build data for 9 scatter plots (prop vs total) -------------------------
# For event 1 across samples:
imp_evt1 <- sampled_impact[, 1, , drop = FALSE] # dims: [cells x 1 x samples]
# Collapse the middle dim:
imp_mat <- drop(imp_evt1)                         # dims: [cells x samples]

# Total displacement per sample (across all cells)
total_disp_per_sample <- colSums(imp_mat, na.rm = TRUE)

# For each section, compute proportion displaced (per sample) within that section
pop_vec <- values(ODDy$Population)

section_summaries <- sections %>%
  mutate(
    df = map(cell_idx, function(idx) {
      if (length(idx) == 0) return(tibble(prop = NA_real_, total = total_disp_per_sample))
      # numerator: total displaced in section per sample
      num <- colSums(imp_mat[idx, , drop = FALSE], na.rm = TRUE)
      # denominator: population of section
      den <- sum(pop_vec[idx], na.rm = TRUE)
      tibble(
        prop  = if (is.finite(den) && den > 0) num / den else NA_real_,
        total = total_disp_per_sample
      )
    })
  ) %>%
  dplyr::select(name, df) %>%
  unnest(df) %>%
  mutate(total_log10 = log10(total))

# --- 4a) Faceted grid of the 9 scatter plots (clean + compact) ----------------
scatters_faceted <- ggplot(section_summaries,
                           aes(x = prop, y = total_log10)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ name, ncol = 3, scales = "free_x") +
  labs(
    x = "Proportion displaced (section)",
    y = "Total displacement (log10)"
  ) +
  theme_minimal()

print(scatters_faceted)

# --- 4b) (Optional) Produce 9 separate plots as a named list -------------------
true_total <- 219745

observed_prop <- c(
  "Bourdet"              = 0.32,
  "Fond Rouge Torbeck"   = 0.25,
  "Fond Rouge Dahere"    = 0.28,
  "Bariadelle"           = 0.96,
  "Frangipane"           = 0.49,
  "Champlois"            = 0.80,
  "Levy Marsan"          = 0.37,
  "Brodequin"            = 0.15,
  "Barbois"              = 0.32
)

# 3) Build the per-section plots with transformed y and the red point
scatter_list <- section_summaries %>%
  group_split(name) %>%
  setNames(unique(section_summaries$name)) %>%
  lapply(function(df_sec) {
    sec <- unique(df_sec$name)
    obs_x <- unname(observed_prop[sec])
    
    ggplot(df_sec, aes(x = prop, y = total)) +
      geom_point(alpha = 0.7, col = '#440154') +
      geom_point(aes(x = obs_x, y = true_total),
                 inherit.aes = FALSE, colour = "red", size = 3, shape = 4, stroke = 1.5) +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 10, base = 100),
        breaks = c(0, 100, 1000, 10000, 100000, 1000000),
        labels = expression(phantom(1) * 0^phantom(0), 10^2, 10^3, 10^4, 10^5, 10^6),
        minor_breaks = NULL,
        expand=c(0,0)
      ) +
      scale_x_continuous(expand=c(0,0)) +
      labs(
        x = paste0("Proportion displaced: ", sec),
        y = "Total displacement"
      ) +
      theme_minimal(base_family = "Times New Roman") +
      theme(
        # border around the **panel**
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        plot.margin = margin(t=10, l=10, b=10, r=10)
      ) 
  })

# Preview
scatter_list[[1]]

plot_grid(map_full, 
          plot_grid(scatter_list[[1]], scatter_list[[2]], scatter_list[[3]],
          scatter_list[[4]], scatter_list[[5]], scatter_list[[6]],
          scatter_list[[7]], scatter_list[[8]], scatter_list[[9]],
          ncol = 3, align = "hv", 
          labels=c('(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain", 
          label_size = 12.5, label_x = 0.005, label_y = 0.99), nrow=2, rel_heights=c(0.8, 1.2),
          labels=c('(a)', ''),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain", 
          label_x = -0.007)
# HTI_Disp_Joint, 9.5 x 10.5
# HTI_Disp_Joint_Margs, 9.5 x 10.5 (portrait)

#----------------------------------------------------
#------------- Marginal Densities -------------------

library(ggExtra)

scatter_list <- section_summaries %>%
  group_split(name) %>%
  setNames(unique(section_summaries$name)) %>%
  lapply(function(df_sec) {
    
    sec <- unique(df_sec$name)
    obs_x <- unname(observed_prop[sec])
    
    base_plot <- ggplot(df_sec, aes(x = prop, y = total)) +
      geom_point(alpha = 0.6, colour = '#440154') +
      geom_point(aes(x = obs_x, y = true_total),
                 inherit.aes = FALSE, colour = "red",
                 size = 3, shape = 4, stroke = 1.3) +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 10, base = 100),
        breaks = c(0, 100, 1000, 10000, 100000, 1000000),
        labels = expression(0, 10^2, 10^3, 10^4, 10^5, 10^6),
        minor_breaks = NULL
      ) +
      labs(
        x = paste0("Proportion displaced: ", sec),
        y = "Total displacement"
      ) +
      theme_minimal(base_family = "Times New Roman") +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        plot.margin = margin(t=10, l=20, b=10, r=10)
      )
    
    # ---- NEW: add marginal density plots ----
    ggMarginal(
      base_plot,
      type = "histogram",      # or "histogram"
      fill = "#440154",      # optional, can remove for auto colour
      alpha = 0.5,
      size = 4
    ) 
    
  })

#----------------------------------------------------


pairplot_communes <- function(impact_prop){
  
  # 1) Samples to wide: one column per commune, rows = samples
  sampled_df <- (impact_prop %>%
                   dplyr::select(polygon_name, dplyr::starts_with("sampled.")) %>%
                   tidyr::pivot_longer(cols = -polygon_name, names_to = "sample", values_to = "value") %>%
                   tidyr::pivot_wider(names_from = polygon_name, values_from = value) %>%
                   dplyr::mutate(type = "Sampled"))[, -1]
  
  # ensure numeric
  sampled_df[1:(ncol(sampled_df)-1)] <- lapply(sampled_df[1:(ncol(sampled_df)-1)], as.numeric)
  
  # 2) Observed row in the same wide format
  observed_row <- impact_prop %>%
    dplyr::select(polygon_name, observed) %>%
    tidyr::pivot_wider(names_from = polygon_name, values_from = observed) %>%
    dplyr::mutate(type = "Observed")
  
  observed_row[1:(ncol(observed_row)-1)] <- lapply(observed_row[1:(ncol(observed_row)-1)], as.numeric)
  
  # 3) Combine
  plot_df <- dplyr::bind_rows(sampled_df, observed_row)
  obs_vals <- dplyr::filter(plot_df, type == "Observed")
  
  # Lower-panel painter: sampled points + observed red cross
  overlay_observed <- function(data, mapping, ...) {
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::geom_point(alpha = 0.7, size = 0.5, color = "#440154", ...) +
      ggplot2::geom_point(data = obs_vals, mapping = mapping,
                          color = "red", size = 3, shape = 4, stroke = 1.5)
  }
  
  # Utility to drop the upper triangle & first row/col (like your gpairs_lower)
  gpairs_lower <- function(g){
    g$plots <- g$plots[-(1:g$nrow)]
    g$yAxisLabels <- g$yAxisLabels[-1]
    g$nrow <- g$nrow - 1
    
    g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
    g$xAxisLabels <- g$xAxisLabels[-g$ncol]
    g$ncol <- g$ncol - 1
    g
  }
  
  # 4) Pairplot (proportions, so keep 0–1 scales)
  g <- GGally::ggpairs(
    plot_df %>% dplyr::select(-type),
    upper = list(continuous = "blank"),
    lower = list(continuous = overlay_observed),
    diag  = list(continuous = "blankDiag"),
    showStrips = TRUE
  ) +
    ggplot2::theme_minimal(base_family = "Times") +
    ggplot2::theme(
      strip.text   = ggplot2::element_text(size = 10.5),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2), minor_breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), minor_breaks = NULL)
  
  gpairs_lower(g)
}

pop_vec <- values(ODDy$Population)

# 2) Helper: return per-sample *proportion displaced* for a given section index set
section_props <- function(idx) {
  if (length(idx) == 0) return(rep(NA_real_, ncol(imp_mat)))
  den <- sum(pop_vec[idx], na.rm = TRUE)
  if (!is.finite(den) || den <= 0) return(rep(NA_real_, ncol(imp_mat)))
  num <- colSums(imp_mat[idx, , drop = FALSE], na.rm = TRUE)
  num / den
}

# 3) Build wide data.frame: one row per commune/section
n_samp <- ncol(imp_mat)
sample_cols <- paste0("sampled.", seq_len(n_samp))

impact_prop <- purrr::map2_dfr(
  sections$cell_idx,
  sections$name,
  ~ {
    props <- section_props(.x)
    tibble::tibble(
      polygon_name = .y,
      observed     = unname(observed_prop[.y]),
      !!!rlang::set_names(as.list(as.numeric(props)), sample_cols)
    )
  }
)

# (Optional) sanity checks
# head(impact_prop)
# colSums(!is.na(impact_prop[sample_cols]))  # how many non-NA per sample
# summary(impact_prop$observed)

# 4) Make the pairplot (using the function you have)
g <- pairplot_communes(impact_prop)
print(g)











# bourdet_prop_disp = as.numeric(colSums(sampled_impact[bourdet_cells,1, , drop=F])/sum(values(ODDy$Population)[bourdet_cells], na.rm=T))
# 
# plot(density(log10(colSums(sampled_impact[,1,]))))
# lines(density(log10(colSums(sampled_impact[,1,])[bourdet_prop_disp > 0.1 & bourdet_prop_disp < 0.3])), col='red')
# abline(v=log10(219745))



# ODDy <- ODDy_with_sampled$ODDyWithImpact
# sampled_full <- ODDy_with_sampled$sampled_full
# 
# hist_AFG = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
# hist_AFG +  geom_vline(xintercept=4.7, color='#D73027', size=1) + xlab("Predicted Fatalities (AFG 2023-10-07)")

#------------------------------------------------------------------------------------
#--------------------------------- Turkey 2023-02-06 --------------------------------
#------------------------------------------------------------------------------------


ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjects/EQ20230206TUR_169_withBuildings')
namer = 'EQ20230206TUR_169'

folder_write = 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/'
# 
# # Sample impact:
ODDy_with_sampled = PosteriorImpactPred(ODDy,
                                        AlgoResultsFilename='IIDIPUS_Results/HPC/abcsmc_start_step_2025-10-11_174442.924329_Sept25_NormalCDF_ESplus0.1AVE_RFwithTot_range1',
                                        folder_write=folder_write,
                                        namer=namer,
                                        N_samples = 200)

sampled_impact = ODDy_with_sampled$sampled_full
ODDy = ODDy_with_sampled$ODDyWithImpact

par(mfrow=c(2,2))
plot(log10(colSums(sampled_impact[,1,])), log10(colSums(sampled_impact[,2,])))
points(log10(4500000), log10(51882), col='red', pch=19)

plot(log10(colSums(sampled_impact[,2,])), log10(colSums(sampled_impact[,3,])))
points(log10(51882),log10(736760), col='red', pch=19)

plot(log10(colSums(sampled_impact[,1,])), log10(colSums(sampled_impact[,3,])))
points(log10(4500000),log10(736760), col='red', pch=19)

impact_df = ODDy_with_sampled$ODDyWithImpact@impact
impact_df$polygon_name = paste0(paste0(
  toupper(substring(impact_df$impact, 1, 1)),
  substring(impact_df$impact, 2)
), ':\n', impact_df$polygon_name)

impact_df = impact_df[order(impact_df$observed, decreasing=T),]

pairplot_regions(impact_df, cross_size = 1.5)
#TUR_Pairplot, 8.5 x 8.5

# 1) Align CRS (work in the raster CRS)
kah_shp <- read_sf('/home/manderso/Documents/USGS/For Max/1New_Turkey/1New_Turkey/1Kahramanmaras/Kahramanmaras_GT_BDarea_adminID.shp')
kah_sf <- st_transform(kah_shp, crs(ODDy))
kah_v  <- vect(kah_sf)                # sf -> SpatVector

ant_shp <- read_sf('/home/manderso/Documents/USGS/For Max/1New_Turkey/1New_Turkey/1Antakya/Antakya_GT_BDarea_adminID.shp')
ant_sf <- st_transform(ant_shp, crs(ODDy))
ant_v  <- vect(ant_sf)                # sf -> SpatVector

osm_shp <- read_sf('/home/manderso/Documents/USGS/For Max/1New_Turkey/1New_Turkey/1Osmaniye/Osmaniye_GT_BDarea_adminID.shp')
osm_sf <- st_transform(osm_shp, crs(ODDy))
osm_v  <- vect(osm_sf)                # sf -> SpatVector
osm_v <- crop(osm_v, ext(36.2, 36.3, 37.04, 37.12))
mean(!is.na(osm_v$GT))

gaz_shp <- read_sf('/home/manderso/Documents/USGS/For Max/1New_Turkey/1New_Turkey/1Gaziantep/Gaziantep_GT_BDarea_adminID.shp')
gaz_sf <- st_transform(gaz_shp, crs(ODDy))
gaz_v  <- vect(gaz_sf)                # sf -> SpatVector
mean(!is.na(gaz_v$GT))
gaz_v <- crop(gaz_v, ext(37.3,37.47, 37, 37.12 ))
mean(!is.na(gaz_v$GT))

make_city_plot_nomarg <- function(city_v, city_name, obs_total = NA_real_, show_legend = FALSE) {
  
  # cells whose centers fall inside the city's extent
  ex <- terra::extract(ODDy, ext(city_v), cells = TRUE)
  stopifnot(nrow(ex) > 0)
  
  # x: per-sample proportion in this city's extent
  x_prop <- as.numeric(colSums(sampled_impact[ex$cell, 3, , drop = FALSE])) /
    sum(values(ODDy[["nBuildings"]])[ex$cell])
  
  # y: total across all cells (same for every city) -- replace with your set if needed
  y_tot <- colSums(sampled_impact[, 3, ], na.rm = TRUE)
  
  df <- tibble(prop = x_prop, total = y_tot)
  
  # observed x from GT flag (your earlier definition)
  obs_x <- mean(!is.na(city_v$GT))
  
  # provider points for legend + overlay
  provider_pts <- tibble(x = numeric(0), y = numeric(0), source = character(0))
  
  # red cross = MEUCC (keep it at your existing location)
  if (is.finite(obs_x) && !is.na(obs_total)) {
    provider_pts <- bind_rows(provider_pts,
                              tibble(x = obs_x, y = obs_total, source = "MEUCC Field Survey"))
  }
  
  # yellow cross = Microsoft AI for Good (Kahramanmaraş only)
  if (identical(city_name, "Kahramanmaraş")) {
    provider_pts <- bind_rows(provider_pts,
                              tibble(x = 0.0744, y = obs_total, source = "Microsoft AI for Good"))
  }
  
  # green crosses = Copernicus (Kahramanmaraş and Gaziantep)
  if (identical(city_name, "Kahramanmaraş")) {
    provider_pts <- bind_rows(provider_pts,
                              tibble(x = 0.1593238, y = obs_total, source = "Copernicus"))
  }
  if (identical(city_name, "Gaziantep")) {
    provider_pts <- bind_rows(provider_pts,
                              tibble(x = 0.004234298, y = obs_total, source = "Copernicus"))
  }
  
  # base scatter in your house style
  p <- ggplot(df, aes(x = prop, y = total)) +
    geom_point(alpha = 0.7, colour = "#440154") +
    # overlay provider crosses (all as shape 4; colors via scale below)
    { if (nrow(provider_pts) > 0)
      geom_point(
        data = provider_pts,
        aes(x = x, y = y, colour = source),
        inherit.aes = FALSE,
        shape = 4, size = 2, stroke = 2.2
      )
      else NULL } +
    scale_colour_manual(
      name   = "Source of City-Level Damage Data",
      values = c("MEUCC Field Survey"     = "red",
                 "Microsoft AI for Good"  = "#FDE725",
                 "Copernicus"             = "#1FA187"),
      breaks = c("MEUCC Field Survey", "Microsoft AI for Good", "Copernicus")
    ) +
    guides(
      colour = guide_legend(override.aes = list(shape = 4, size = 2, stroke = 2.2))
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 10, base = 100),
      breaks = c(0, 100, 1000, 10000, 100000, 1e6),
      labels = expression(phantom(1) * 0^phantom(0), 10^2, 10^3, 10^4, 10^5, 10^6),
      minor_breaks = NULL, expand = c(0, 0)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(
      x = paste0("Proportion damaged: ", city_name),
      y = "Total buildings damaged"
    ) +
    # scale_colour_manual(
    #   name   = 'Source of City-Level Damage Data',
    #   values = c("MEUCC Field Survey" = "red",
    #              "Microsoft AI for Good" = "gold",
    #              "Copernicus" = "forestgreen"),
    #   breaks = c("MEUCC Field Survey", "Microsoft AI for Good", "Copernicus")
    # ) +
    # scale_shape_manual(
    #   name   = NULL,
    #   values = c("MEUCC Field Survey" = 4,
    #              "Microsoft AI for Good" = 4,
    #              "Copernicus" = 4),
    #   breaks = c("MEUCC Field Survey", "Microsoft AI for Good", "Copernicus")
    # ) +
    # guides(
    #   colour = guide_legend(override.aes = list(size = 4, stroke = 1.5)),
    #   shape  = guide_legend(override.aes = list(size = 4, stroke = 1.5))
    # ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.margin    = margin(t = 10, l = 10, b = 20, r = 10),
      legend.position = if (show_legend) "right" else "none",
      legend.text    = element_text(size = 11, family = "Times New Roman"),
      legend.title   = element_text(size = 11, family = "Times New Roman"),
      axis.text = element_text(size = 11, family = "Times New Roman")
    )
  
  p
}

make_city_plot <- function(city_v, city_name, obs_total = NA_real_, show_legend = FALSE) {
  
  # cells whose centers fall inside the city's extent
  ex <- terra::extract(ODDy, ext(city_v), cells = TRUE)
  stopifnot(nrow(ex) > 0)
  
  x_prop <- as.numeric(colSums(sampled_impact[ex$cell, 3, , drop = FALSE])) /
    sum(values(ODDy[["nBuildings"]])[ex$cell])
  
  y_tot <- colSums(sampled_impact[, 3, ], na.rm = TRUE)
  
  df <- tibble(prop = x_prop, total = y_tot)
  
  obs_x <- mean(!is.na(city_v$GT))
  
  provider_pts <- tibble(x = numeric(0), y = numeric(0), source = character(0))
  
  if (is.finite(obs_x) && !is.na(obs_total)) {
    provider_pts <- bind_rows(provider_pts,
                              tibble(x = obs_x, y = obs_total, source = "MEUCC Field Survey"))
  }
  if (identical(city_name, "Kahramanmaraş")) {
    provider_pts <- bind_rows(provider_pts,
                              tibble(x = 0.0744, y = obs_total, source = "Microsoft AI for Good"),
                              tibble(x = 0.1593238, y = obs_total, source = "Copernicus"))
  }
  if (identical(city_name, "Gaziantep")) {
    provider_pts <- bind_rows(provider_pts,
                              tibble(x = 0.004234298, y = obs_total, source = "Copernicus"))
  }
  
  # ------------------ BUFFER TO PREVENT CLIPPING -----------------------
  xr <- range(df$prop, provider_pts$x, na.rm = TRUE)
  yr <- range(df$total, provider_pts$y, na.rm = TRUE)
  xb <- diff(xr) * 0.03   # 3% buffer
  yb <- diff(yr) * 0.03
  
  # ------------------ MAIN SCATTER PLOT ------------------------------
  p <- ggplot(df, aes(x = prop, y = total)) +
    geom_point(alpha = 0.7, colour = "#440154") +
    { if (nrow(provider_pts) > 0)
      geom_point(
        data = provider_pts,
        aes(x = x, y = y, colour = source),
        inherit.aes = FALSE,
        shape = 4, size = 2, stroke = 2.2
      )
      else NULL } +
    scale_colour_manual(
      name   = "Source of City-Level Damage Data",
      values = c(
        "MEUCC Field Survey"    = "red",
        "Microsoft AI for Good" = "#FDE725",
        "Copernicus"            = "#1FA187"
      )
    ) +
    guides(
      colour = guide_legend(override.aes = list(shape = 4, size = 2, stroke = 2.2))
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 10, base = 100),
      breaks = c(0, 100, 1000, 10000, 100000, 1e6),
      labels = expression(phantom(1) * 0^phantom(0), 10^2, 10^3, 10^4, 10^5, 10^6),
      minor_breaks = NULL
    ) +
    # coord_cartesian(
    #   xlim = c(xr[1] - xb, xr[2] + xb),
    #   ylim = c(yr[1] - yb, yr[2] + yb)
    # ) +
    labs(
      x = paste0("Proportion damaged: ", city_name),
      y = "Total buildings damaged"
    ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.margin    = margin(t = 20, r = 20, b = 25, l = 25),
      legend.position = if (show_legend) "right" else "none",
      legend.text    = element_text(size = 11, family = "Times New Roman"),
      legend.title   = element_text(size = 11, family = "Times New Roman"),
      axis.text      = element_text(size = 11, family = "Times New Roman")
    )
  
  # ------------------ ADD MARGINAL DENSITIES ---------------------------
  ggMarginal(
    p,
    type = "histogram",
    margins = "both",
    fill = "#440154",
    alpha = 0.5,
    size = 4
  )
}

# --- observed total used for all crosses (per your snippet) --------------------
obs_total <- 736760

# --- build panels (no legends inside) -----------------------------------------
p_kah <- make_city_plot(kah_v, "Kahramanmaraş", obs_total, show_legend = FALSE)
p_ant <- make_city_plot(ant_v, "Antakya",       obs_total, show_legend = FALSE)
p_osm <- make_city_plot(osm_v, "Osmaniye",      obs_total, show_legend = FALSE)
p_gaz <- make_city_plot(gaz_v, "Gaziantep",     obs_total, show_legend = FALSE)

# --- pull a clean legend from a copy of one plot -------------------------------
legend_plot <- make_city_plot_nomarg(kah_v, "Kahramanmaraş", obs_total, show_legend = TRUE) + theme(legend.position = 'bottom')
legend_grob <- get_plot_component(legend_plot, 'guide-box', return_all=T)[[3]]

# --- arrange 2x2 grid + legend to the right -----------------------------------
grid_2x2 <- plot_grid(
  p_kah, p_ant, p_osm, p_gaz,
  ncol = 2, align = "hv",
  labels = c("(a)", "(b)", "(c)", "(d)"),
  label_fontfamily = "Times New Roman",
  label_fontface   = "plain",
  label_size       = 12.5,
  label_x          = 0, label_y = 0.98
)

plot_grid(grid_2x2, legend_grob, ncol = 1, rel_heights = c(1, 0.1))
# TUR_BuildingUpdates.pdf, 8 x 6





BDy_TUR <- readBD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Aug24/BDobjects_MAR/EQ20230206TUR_169')
plot(ODDy$Population, xlim=c(36, 38), ylim=c(36.5, 38))
points(BDy_TUR, col='red', lwd=0.01)

kah_BDy = crop(BDy_TUR, ext(kah_sf))
sum(kah_BDy$grading != 'notaffected')/sum(kah_BDy$grading != 'possible')

gaz_BDy = crop(BDy_TUR, ext(gaz_sf))
sum(gaz_BDy$grading != 'notaffected')/sum(gaz_BDy$grading != 'possible')

# 
# rastered_BDy = rasterize(BDy_TUR, ODDy, fun='sum')
# 
# BDy_TUR_nonna = BDy_TUR[BDy_TUR$grading != 'possible',]
# rastered_BDy_count = rasterize(BDy_TUR_nonna, ODDy, fun='sum')
# 
# BDy_TUR$dam_col = BDy_TUR$grading
# BDy_TUR$dam_col = ifelse(BDy_TUR$dam_col == 'notaffected', 0, 1)
# 
# rastered_BDy <- terra::rasterize(
#   BDy_TUR,
#   ODDy,
#   field = "dam_col",
#   fun   = mean,
#   na.rm = TRUE
# )
# 
# ODDy$nSatBuildings = rastered_BDy_count$sum
# ODDy$propDam = rastered_BDy$mean
# 
# ODDy_df = as.data.frame(ODDy, na.rm=F)
# ODDy_df %>% filter(!is.na(nSatBuildings))
# ggplot(ODDy_df, aes(x=hazMean1, y=propDam, size=nSatBuildings)) + 
#   geom_point()
# 
# non_na_ids = which(!is.na(ODDy_df$nSatBuildings))
# plot(colSums(sampled_impact[non_na_ids,3,])/sum(ODDy$nBuildings[non_na_ids]), log10(colSums(sampled_impact[,3,])))
# abline(h=log10(736760))
# abline(v=0.07)
# abline(v=sum(ODDy_df$propDam * ODDy_df$nSatBuildings, na.rm=T))
# 
# coords = crds(ODDy, na.rm=F)
# ids_kahramanmaras <- which(coords[,1]> 36.78372 & coords[,1] < 37.03483 & coords[,2] > 37.50237 & coords[,2] < 37.64031)
# plot(colSums(sampled_impact[ids_kahramanmaras,3,])/sum(ODDy$nBuildings[ids_kahramanmaras]), log10(colSums(sampled_impact[,3,])))
# abline(h=log10(736760)); abline(v=0.4)
# 
# BDy <- readBD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/BDobjects/EQ20230206TUR_169')
# BDy_df = as.data.frame(BDy, na.rm=F)
# BDy_df$longitude = crds(BDy)[,1]
# BDy_df$latitude = crds(BDy)[,2]
# BDy_df %<>% filter(grading != 'possible')
# ggplot(BDy_df, aes(x=longitude, y=latitude, col=grading)) + 
#   geom_point() + xlim(36.75, 37) + ylim(37.5, 37.7)
#   + xlim(37.3, 37.5) + ylim(37, 37.2) #

# IGNORE: 
# 
# bourdet_cells = which(coords$x > -73.82 & coords$x < -73.7 & coords$y > 18.15 & coords$y < 18.27)
# bourdet_cells
# points(coords[bourdet_cells,])
# 
# ODDy_df = as.data.frame(ODDy, na.rm=F)
# ODDy_df$x = coords$x
# ODDy_df$y = coords$y
# 
# ODDy_df$hazMean1[which(is.na(ODDy_df$Population))] = NA
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = hazMean1)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.5) +
#   geom_sf(data = admin3_shp[1,], fill = 'green', colour = "red", linewidth = 0.5) +
#   coord_sf(xlim = c(-74.55, -73), ylim = c(17.95, 18.72), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(5, 8.5), oob = scales::censor) +
#   geom_text(aes(x=-73.76, y=18.21,label='Bourdet'),col='white', fontface = "bold") +
#   geom_text(aes(x=-74.18, y=18.65,label='Fond Rouge Torbeck'),col='white', fontface = "bold") + geom_point(aes(x=-74.18, y=18.663), col='white')
# 
# 
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(xlim = c(-74, -73.5), ylim = c(18, 18.4), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[bourdet_cells,], aes(x=x, y=y))
# 
# #bourdet_cells = bourdet_cells[5]
# propDisplaced = sweep(sampled_impact[,1,], 1, values(ODDy$Population), "/")
# plot(colSums(sampled_impact[bourdet_cells,1, , drop=F])/sum(values(ODDy$Population)[bourdet_cells], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Bourdet Proportion displaced', 
#      ylab= 'Total Displacement')
# points(0.32, log10(219745), col='red', pch=19)
# 
# 
# rouge_torbeck_cells = which(coords$x > -74.15 & coords$x < -74.1 & coords$y > 18.6 & coords$y < 18.65)
# 
# # Fond Rouge torbeck
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(ylim = c(18.3, 18.75), xlim = c( -74.4, -73.7), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[rouge_torbeck_cells,], aes(x=x, y=y))
# 
# plot(colSums(sampled_impact[rouge_torbeck_cells ,1, , drop=F])/sum(values(ODDy$Population)[rouge_torbeck_cells ], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Bourdet Proportion displaced', 
#      ylab= 'Total Displacement')
# points(0.25, log10(219745) , col='red', pch=19)
# 
# 
# rouge_dahere_cells = which(coords$x > -74.15 & coords$x < -74.1 & coords$y > 18.6 & coords$y < 18.65)
# 
# # Fond Rouge Dahere
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(ylim = c(18.3, 18.75), xlim = c( -74.4, -73.7), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[rouge_dahere_cells,], aes(x=x, y=y))
# 
# plot(colSums(sampled_impact[rouge_dahere_cells,1, , drop=F])/sum(values(ODDy$Population)[rouge_dahere_cells], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Bourdet Proportion displaced', 
#      ylab= 'Total Displacement')
# points(0.28, log10(219745) , col='red', pch=19)
# 
# 
# plot(colSums(sampled_impact[rouge_dahere_cells,1, , drop=F])/sum(values(ODDy$Population)[rouge_dahere_cells], na.rm=T),
#      colSums(sampled_impact[bourdet_cells,1, , drop=F])/sum(values(ODDy$Population)[bourdet_cells], na.rm=T))
# points(0.25, 0.32, col='red', pch=19)
# 
# 
# bariadelle_cells = which(coords$x > -74.45 & coords$x < -74.4 & coords$y > 18.5 & coords$y < 18.55)
# 
# # Bariadelle
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(xlim = c( -74.7, -74), ylim = c(18.4, 18.8), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[bariadelle_cells,], aes(x=x, y=y))
# 
# plot(colSums(sampled_impact[bariadelle_cells,1, , drop=F])/sum(values(ODDy$Population)[bariadelle_cells], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Bariadelle Proportion displaced', 
#      ylab= 'Total Displacement')
# points(c(0.96), log10(219745), col='red', pch=19)
# 
# 
# frangipane_cells = which(coords$x > -73.1 & coords$x < -73 & coords$y > 18.2 & coords$y < 18.3)
# 
# # Frangipane
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(xlim = c(-74.55, -73), ylim = c(17.95, 18.72), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[frangipane_cells,], aes(x=x, y=y))
# 
# plot(colSums(sampled_impact[frangipane_cells,1, , drop=F])/sum(values(ODDy$Population)[frangipane_cells], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Bariadelle Proportion displaced', 
#      ylab= 'Total Displacement')
# points(c(0.49), log10(219745), col='red', pch=19)
# 
# 
# # Champlois
# 
# champlois_cells = which(coords$x > -73.9 & coords$x < -73.85 & coords$y > 18.33 & coords$y < 18.35)
# 
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(xlim = c(-74.55, -73), ylim = c(17.95, 18.72), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[champlois_cells,], aes(x=x, y=y))
# 
# plot(colSums(sampled_impact[champlois_cells,1, , drop=F])/sum(values(ODDy$Population)[champlois_cells], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Bariadelle Proportion displaced', 
#      ylab= 'Total Displacement')
# points(c(0.8), log10(219745), col='red', pch=19)
# 
# # Levy Marsan
# levymarsan_cells = which(coords$x > -73.9 & coords$x < -73.85 & coords$y > 18.26 & coords$y < 18.31)
# 
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(xlim = c(-74.55, -73), ylim = c(17.95, 18.72), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[levymarsan_cells,], aes(x=x, y=y))
# 
# plot(colSums(sampled_impact[levymarsan_cells,1, , drop=F])/sum(values(ODDy$Population)[levymarsan_cells], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Bariadelle Proportion displaced', 
#      ylab= 'Total Displacement')
# points(c(0.37), log10(219745), col='red', pch=19)
# 
# 
# brodequin_cells = which(coords$x > -73.43 & coords$x < -73.4 & coords$y > 18.26 & coords$y < 18.31)
# 
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(xlim = c(-74.55, -73), ylim = c(17.95, 18.72), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[brodequin_cells,], aes(x=x, y=y))
# 
# plot(colSums(sampled_impact[brodequin_cells,1, , drop=F])/sum(values(ODDy$Population)[brodequin_cells], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Bariadelle Proportion displaced', 
#      ylab= 'Total Displacement')
# points(c(0.15), log10(219745), col='red', pch=19)
# 
# 
# barbois_cells = which(coords$x > -73.9 & coords$x < -73.85 & coords$y > 18.12 & coords$y < 18.15)
# 
# ggplot() +
#   geom_tile(data = ODDy_df, aes(x = x, y = y, fill = Population)) +
#   geom_sf(data = admin3_shp, fill = NA, colour = "red", linewidth = 0.3) +
#   coord_sf(xlim = c(-74.55, -73), ylim = c(17.95, 18.72), expand = FALSE) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   scale_fill_viridis_c(limits = c(0, 50000), oob = scales::censor) + 
#   geom_point(data=coords[barbois_cells,], aes(x=x, y=y))
# 
# plot(colSums(sampled_impact[barbois_cells,1, , drop=F])/sum(values(ODDy$Population)[barbois_cells], na.rm=T), log10(colSums(sampled_impact[,1,])), 
#      xlab= 'Barbois Proportion displaced', 
#      ylab= 'Total Displacement')
# points(c(0.32), log10(219745), col='red', pch=19)

