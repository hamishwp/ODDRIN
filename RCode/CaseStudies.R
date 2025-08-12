
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
  
  bin_colors <- viridis::viridis(length(levels(plot_df$bin_category)))
  
  
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
pairplot_regions <- function(impact_mort){
  
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
      geom_point(data = obs_vals, mapping = mapping, color = "red", size = 3, shape=4, stroke = 1.5)
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
    theme(strip.text = element_text(size = 10),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) + 
    scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 1, base = 10),
      breaks = c(0, 10, 100, 1000),
      labels = expression(phantom(1) * 0^phantom(0),phantom(1) * 10^phantom(1),10^2,10^3),
      minor_breaks = NULL
    ) + 
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 1, base = 10),
      breaks = c(0, 10, 100, 1000),
      labels = expression(phantom(1) * 0^phantom(0),phantom(1) * 10^phantom(1),10^2,10^3),
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
par(mfrow = c(3, 3), cex = 1.05,   # reduce inner margins
    oma = c(0, 0, 0, 0) )
par(family = "Times")
# Plot first figure
plot(ODDy$Population, main = "Population", font.main = 1, cex.main = 1)
mtext("(a)", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$nBuildings, main = "Building Count", font.main = 1, cex.main = 1)
mtext("(b)", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$hazMean1, main = "Hazard 1 MMI", font.main = 1, cex.main = 1)
mtext("(c)", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$hazMean2, main = "Hazard 2 MMI", font.main = 1, cex.main = 1)
mtext("(d)", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$SHDI, main = "SHDI", type = 'continuous', font.main = 1, cex.main = 1)
mtext("(e)", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$GNIc, main = "GNIc", type = 'continuous', font.main = 1, cex.main = 1)
mtext("(f)", side = 3, line =  2.45, adj = -0.05, cex = 1)

plot(ODDy$EQFreq, main = "EQFreq", font.main = 1, cex.main = 1)
mtext("(h)", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$Vs30, main = "Vs30", font.main = 1, cex.main = 1)
mtext("(i)", side = 3, line = 2.45, adj = -0.05, cex = 1)

plot(ODDy$PDens, main = "PopDens", font.main = 1, cex.main = 1)
mtext("(j)", side = 3, line = 2.45, adj = -0.05, cex = 1)

#PHL_Data_all.pdf, 8.7 x 8 (or 1300 x 680)

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
      polygon_list_plot[[i]]$impact = 0
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
      scale_fill_viridis() +  # Color scale
      theme_minimal() +
      #geom_sf_text(aes(label = NAME_2), size = 3, color = "black") + 
      labs(fill = ifelse(impact_type=='buildDam', 'Building Damage', ifelse(impact_type=='mortality', 'Mortality', 'Displacement'))) + 
      xlim(plot_bbox[c(1,2)]) + ylim(plot_bbox[c(3,4)]) #+
      # theme(
      #   axis.title.y = element_text(family = "Times New Roman", size = 12),
      #   axis.text.x = element_text(family = "Times New Roman", size = 12),
      #   axis.text.y = element_text(family = "Times New Roman", size = 12),
      #   axis.title.x = element_text(family = "Times New Roman", size = 12),
      #   plot.title = element_text(family = "Times New Roman", size = 14),
      #   legend.text = element_text(family = "Times New Roman", size = 12),
      #   legend.title = element_text(family = "Times New Roman", size = 12)
      # ) 
    if (impact_type=='buildDam'){
      p_impact = p_impact + scale_fill_viridis(breaks=c(1000,10000,20000))
    }
  } else {
    p_impact <- ggplot(data = combined_polygons) +
      geom_sf(aes(fill = as.factor(impact)), color = '#440154', size = 0) +
      scale_fill_viridis_d() +  # Color scale
      theme_minimal() +
      #geom_sf_text(aes(label = NAME_2), size = 3, color = "black") + 
      labs(fill = ifelse(impact_type=='buildDam', 'Building Damage', ifelse(impact_type=='mortality', 'Mortality', 'Displacement'))) + 
      xlim(plot_bbox[c(1,2)]) + ylim(plot_bbox[c(3,4)]) #+ 
      # theme(
      #   axis.title.y = element_text(family = "serif", size = 12),
      #   axis.text.x = element_text(family = "serif", size = 12),
      #   axis.text.y = element_text(family = "serif", size = 12),
      #   axis.title.x = element_text(family = "serif", size = 12),
      #   plot.title = element_text(family = "serif", size = 14),
      #   legend.text = element_text(family = "serif", size = 12),
      #   legend.title = element_text(family = "serif", size = 12)
      # ) 
  }
  
  #hazard = apply(values(ODDy[grep("hazMean",names(ODDy),value = T), drop=F]),1,max, na.rm=T)
  hazard = values(ODDy$hazMean1)
  hazard[hazard < 4.98] = NA
  ODDy$hazard<-hazard
  brks<-seq(9,ceiling(2*max(hazard,na.rm = T)),by=1)/2
  
  ODDy_df <- as.data.frame(ODDy, na.rm=F, xy=T)
  names(ODDy_df)[which(names(ODDy_df)=='x')] = 'Longitude'
  names(ODDy_df)[which(names(ODDy_df)=='y')] = 'Latitude'
  
  p<-p_impact+geom_contour(data = ODDy_df,
                           mapping = aes(Longitude,Latitude,z=hazard,colour=..level..),
                           alpha=0.8,breaks = brks, size=1, 
                           show.legend = ifelse(haz_legend, T, F)) +
    scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
    labs(colour = "MMI") + 
    scale_x_continuous(labels = scales::number_format(accuracy = 1), limits=c(plot_bbox[c(1,2)]), expand=c(0,0)) +  # Remove degrees, keep plain numbers
    scale_y_continuous(labels = scales::number_format(accuracy = 1), limits=c(plot_bbox[c(3,4)]), expand=c(0,0)) +
    theme(
      panel.grid = element_blank(),  # Remove gridlines
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = element_line(),   # Ensure ticks are drawn
      axis.text = element_text()     # Ensure axis labels are visible
    )
  
  return(p)
}



polygons_list = readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp_polygons/EQ20191029PHL_-1polygons')$polygons_list

ODDyWithImpact@impact$polygon_name = unlist(lapply(ODDyWithImpact@polygons[ODDyWithImpact@impact$polygon], function(x) x$name))

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
          p3 + theme(legend.position="bottom",
                     legend.spacing = unit(0.5, "cm")), legend_contour, ncol=4, rel_widths=c(1,1,1,0.3), labels=c('(a)', '(b)', '(c)', ''),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")
#PHL_impact.pdf, 15 x 6

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
      guide = guide_legend(override.aes = list(size = 4, alpha=1))  # Thicker lines in legend
    ) +
    labs(x = "MMI", y = "Impact Probability") +
    coord_cartesian(xlim = c(4.5, 10), ylim = c(0, 1)) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      legend.position = "bottom",
      text = element_text(family = "Times New Roman")
    )
}
p1 = plot_priors(50)

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
      text = element_text(family = "Times New Roman")
    ) 
}

subnat_predDist = plot_predictive_dist_subnat_PHL(impact_mort)
p2 <- subnat_predDist

plot_grid(p1, p2, ncol=2, nrow=1,rel_widths=c(0.4, 0.6), labels=c('(a)', '(b)'), label_x = c(-0.01, -0.01),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

#PriorPredCheck.pdf, 15 x 6

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

input%<>%append(list(Model=Model, #'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup', #
                     PosteriorFileLoc='IIDIPUS_Results/HPC/abcsmc_2025-08-03_093559.602748__July25Agg_NormalCDF_ESplus0.05MVR_RF_kappa0.5_range0.5',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_014135.935892__Apr25Agg_Lognormal_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot', #IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2',#IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
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
pairplot_regions(impact_mort_shortened_names)
#JPN_Pairplot.png, 1300 x 880

subnat_predDist = plot_predictive_dist_subnat(impact_mort)
#JPN_subnat_predDist.pdf, 8 x 4

# Load GADM regions:
polygons_list <- readRDS(paste0(dir, 'IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/tmp_polygons/', namer,'polygons'))$polygons_list
plot_bbox = c(136, 138, 36.5, 38)
plot_GADM_impact_polygons(ODDy, polygons_list, plot_bbox=plot_bbox)

p_pixellated <- plot_GADM_impact_polygons_pixellated(ODDy, polygons_list, plot_bbox=plot_bbox, max_break = 24)
#JPN_pixellated_mortPred.pdf, 11 x 9


grid.arrange(p_pixellated + theme(legend.position="bottom",
                     legend.spacing = unit(2, "cm")),
             subnat_predDist, ncol=2, widths=c(0.65, 0.35)
)

plot_grid(p_pixellated + theme(legend.position="bottom",
                               legend.spacing = unit(2, "cm")),
          subnat_predDist, ncol=2, nrow=1,rel_widths=c(0.65, 0.4), labels=c('(a)', '(b)'),
          label_fontfamily = "Times New Roman",
          label_fontface = "plain")

#JPN_subnatPred, 14 x 8



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


ggdraw(plot_grid(p1, p2, ncol=2, nrow=1,rel_widths=c(0.5, 0.5), labels=c('(a)', '(b)'),
                 label_fontfamily = "Times New Roman",
                 label_fontface = "plain",
                 label_x = c(-0.01, -0.02))) + 
  theme(
    plot.margin = margin(t = 10, b = 10, l = 0, r = 0) # 20pt top and bottom margin
  )

#JPN_subnatUpdate.png, 1300 x 365



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

input%<>%append(list(Model=Model, 
                     PosteriorFileLoc='IIDIPUS_Results/mcmc_2025-05-05_225121.802448__Apr25Agg_NormalCDF_ESplus0.05BD_RFNoTot_range.5_kappa.5',#'IIDIPUS_Results/HPC/mcmc_2025-04-27_015007.755409__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot',  # 'IIDIPUS_Results/HPC/mcmc_2025-04-27_014331.085601__Apr25Agg_Logistic_ESplus0.05BD_RFwithTot_backup',#'IIDIPUS_Results/HPC/mcmc_2025-04-28_203305.949572__Apr25Agg_NormalCDF_ESplus0.05BD_RFwithTot_RFsmooth0.2', #'IIDIPUS_Results/mcmc_2025-04-21_203713.178972_backup_NoTot',
                     #PosteriorFileLoc='IIDIPUS_Results/mcmc_2025-04-21_202717.886201_backup',
                     Method=AlgoParams))

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

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/EQ20230908MAR_-1_WithImpact200')
sampled_full <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/sampled_full/EQ20230908MAR_-1_fullSampledImpact')

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
plot_GADM_impact_polygons(ODDy, polygons_list, gadm_level = 2, impact_type = 'mortality', summary = 'mean', plot_bbox = NULL)
  
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


#ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/EQ20250328MMR_-1_WithImpact200')
#^ BAD. Each sampled impact is duplicated
sampled_full <- readRDS('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/IIDIPUS_Input_NewEvents/ODDobjectsWithImpact/sampled_full/EQ20250328MMR_-1_fullSampledImpact')



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
posterior_file_loc = 'IIDIPUS_Results/HPC/abcsmc_2025-08-03_093559.602748__July25Agg_NormalCDF_ESplus0.05MVR_RF_kappa0.5_range0.5'
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
                                        N_samples = 200)
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
                                        N_samples = 200)
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
                                        N_samples = 200)
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
                                        N_samples = 200)
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
                                        N_samples = 200)
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
                                        N_samples = 100)
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
                                        N_samples = 200)
ODDy <- ODDy_with_sampled$ODDyWithImpact
sampled_full <- ODDy_with_sampled$sampled_full

hist_MMR = ODDRIN_vs_PAGER_hists(ODDy, sampled_full, namer)
hist_MMR +  geom_vline(xintercept=4.5, color='#D73027', size=1) + xlab("Predicted Fatalities (MMR 2025-03-28)")


# -------------------------- PLOT: ---------------------------

legend <- get_legend(hist_MMR + theme(legend.position = "right"))

# Combine plots without legends
plots <- plot_grid(
  hist_JPN + theme(legend.position = "none") +  geom_vline(xintercept=3.5 + 270/1000, color='#D73027', size=1) + xlab("Predicted Fatalities (JPN 2024-01-01)"),
  hist_PHL + theme(legend.position = "none") +  geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (PHL 2024-01-08)"),
  hist_RUS + theme(legend.position = "none") +  geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (RUS 2024-01-15)"),
  hist_CHN +  geom_vline(xintercept=1.9, color='#D73027', size=1) + xlab("Predicted Fatalities (CHN 2024-01-22)"),
  #hist_COL + theme(legend.position = "none") +  geom_vline(xintercept=1.5 + 1.5/10, color='#D73027', size=1) + xlab("Predicted Fatalities (COL 2024-01-19)"),
  #hist_ALA + theme(legend.position = "none") +  geom_vline(xintercept=1, color='#D73027', size=1) + xlab("Predicted Fatalities (USA 2024-01-19)"),
  #hist_CHN + theme(legend.position = "none") +  geom_vline(xintercept=1.5 + 3/10, color='#D73027', size=1) + xlab("Predicted Fatalities (CHN 2024-01-22)"),
  hist_TWN + theme(legend.position = "none") +  geom_vline(xintercept=2.5 + 9/90+1/10, color='#D73027', size=1) + xlab("Predicted Fatalities (TWN 2024-04-03)"),
  hist_MMR + theme(legend.position = "none") +  geom_vline(xintercept=4.5 + 2804/9000+1/10, color='#D73027', size=1) + xlab("Predicted Fatalities (MMR 2025-03-28)"),
  ncol = 3,
  labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'),
  label_fontfamily = "Times New Roman",
  label_fontface = "plain",
  label_x = -0.02  # Move label to the left
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
