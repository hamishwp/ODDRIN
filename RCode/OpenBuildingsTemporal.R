

ODDy <- readODD('/home/manderso/Documents/GitHub/ODDRIN/IIDIPUS_Input_Alternatives/Sept25/ODDobjects/Train/EQ20170529IDN_62')
ext_ODD = ext(ODDy)
wkd_ODD = paste0('POLYGON((',ext_ODD[1], ' ', ext_ODD[3],', ', ext_ODD[2],' ',ext_ODD[3],', ', ext_ODD[2],' ',ext_ODD[4],', ', ext_ODD[1],' ',ext_ODD[4],', ', ext_ODD[1],' ',ext_ODD[3],'))')

#Run colab, 
#then run:
#/snap/bin/gsutil -m cp -I PolyBuildings/ < downloadable_urls.txt

library(terra)

# 1) list your tiles
files <- list.files("/home/manderso/Downloads/PolyBuildings/", pattern = "\\.tif$", full.names = TRUE)

# 2) raster -> polygons where building_presence >= 0.7
tile_to_polys <- function(f, threshold = 0.7) {
  r <- rast(f)                         # SpatRaster
  crs(r) <- "EPSG:32723" 
  r_ll <- project(r, "EPSG:4326", method = "bilinear")
  projectRaster(r$building_fractional_count,
                crs = crs(pop))
  
  rb <- rast(f)  # your .tif
  b  <- if ("building_fractional_count" %in% names(rb)) rb[["building_fractional_count"]] else rb[[1]]
  
  # 2) Make sure target is a SpatRaster too
  pop_tgt <- rast(pop)
  
  xx = project(b, pop_tgt, method = "bilinear")
  
  # pick the presence band (fallback to first band if names absent)
  bname <- "building_fractional_count"
  pres  <- r[[bname]]
  r$building_fractional_count = 1
  r_on_ODDy <- project(r$building_fractional_count, pop, method = "sum")  
  
  pop = rast(raster(ODDy[['Population']]))
  resample(r,pop, method='sum')
  
  # TRUE/FALSE raster of “likely building”
  mask <- pres >= threshold
  
  # polygonize TRUE pixels; drop FALSE/NA
  v <- as.polygons(mask, dissolve = TRUE, values = TRUE, na.rm = TRUE)
  v <- subset(v, v$lyr.1 == 1)
  v$tile <- basename(f)
  v
}

# 3) apply to all tiles and merge
vec_list <- lapply(files, tile_to_polys)
buildings_vect <- do.call(rbind, vec_list)   # SpatVector
buildings_vect



