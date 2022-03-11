
GetSourceFiles<-function(packred){
  
  #@@@@@ SOURCE FILES @@@@@#
  # Basic functions:
  source('RCode/Functions_red.R')
  source('RCode/GetInitialValues.R')
  # S4 object classes required:
  source('RCode/ODDobj.R')
  source('RCode/BDobj.R')
  source('RCode/HAZARDobj.R')
  
  if(!packred){
  
    # Basic functions:
    source('RCode/Functions.R')
    # Disaster related:
    source('RCode/GetGDACS.R')
    source('RCode/GetUSGS.R')
    source('RCode/GetDisaster.R')
    # IDP estimate related:
    source('RCode/GetDisplacements.R')
    source('RCode/GetHelix.R')
    # Demography & population related:
    source('RCode/GetPopDemo.R')
    source('RCode/GetSocioEconomic.R')
    source('RCode/GetINFORM.R')
    source('RCode/GetWorldBank.R')
    # Damage estimate related:
    source('RCode/GetOSM.R')
    source('RCode/GetSatDamage.R')  
    
  }
  
}

LoadLibraries<-function(packred){
  
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(ggplot2)
  library(sf)
  library(ggmap)
  library(geojsonR)
  library(countrycode)
  library(stringr)
  library(pracma)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(ExtDist)
  
  if(!packred) {
    library(codetools)
    library(osmdata)
    library(OpenStreetMap)
    library(osmdata)
  }
  
}

GetODDPackages<-function(packred){

  list.of.packages <- c("ggplot2","sf","tidyverse","openxlsx","pracma",
                        "geojsonR", "tiff", "gstat", "mvtnorm",
                        "RColorBrewer", "geosphere","GGally", "wbstats",
                        "countrycode","rworldmap","rworldxtra","chron","ncdf4",
                        "GADMTools","akima","adehabitatMA","flexsurv", "ExtDist")
  
  if(!packred) list.of.packages<-c(list.of.packages,
                                   "codetools","latex2exp",
                                   "rJava","devtools","OpenStreetMap","osmdata",
                                   "tidyRSS","geojsonR", "tiff", "gstat",
                                   "FactoMineR","factoextra","xtable",
                                   "gsubfn","mapsapi","leaflet", "ssh","RPostgres",
                                   "GADMTools")
  
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  # devtools::install_github('daroczig/fbRads')
  # if(length(list.of.packages[!("openrouteservice" %in% installed.packages()[,"Package"])])){devtools::install_github("rCarto/osrm")}
  if(length(list.of.packages[!("ggmap" %in% installed.packages()[,"Package"])])){devtools::install_github("dkahle/ggmap")}
  # if(length(list.of.packages[!("countrycodes" %in% installed.packages()[,"Package"])])){devtools::install_github("vincentarelbundock/countrycode")}
  if(length(list.of.packages[!("wbstats" %in% installed.packages()[,"Package"])])){devtools::install_github('nset-ornl/wbstats')}
  if(length(list.of.packages[!("wid" %in% installed.packages()[,"Package"])])){devtools::install_github("WIDworld/wid-r-tool")}
  
  LoadLibraries(packred)
  GetSourceFiles(packred)
  
}

GetODDPackages(packred)

filers<-c(paste0(dir,"Demography_Data"),
          paste0(dir,"Demography_Data/SocioEconomic"),
          paste0(dir,"Demography_Data/SocioEconomic/KUMMU"),
          paste0(dir,"Demography_Data/Population"),
          paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2015"),
          paste0(dir,"Plots"),
          paste0(dir,"Plots/IIDIPUS_Results"),
          paste0(dir,"Plots/IIDIPUS_BG"),
          paste0(dir,"IIDIPUS_Input"),
          paste0(dir,"IIDIPUS_Input/ODDobjects"),
          paste0(dir,"IIDIPUS_Input/BDobjects"),
          paste0(dir,"IIDIPUS_Input/HAZARDobjects"),
          paste0(dir,"IIDIPUS_Results"),
          paste0(dir,"IIDIPUS_Results/ODDobjects"),
          paste0(dir,"COPERNICUS_Damage"),
          paste0(dir,"UNOSAT_Damage"),
          paste0(dir,"tmp"))

tmp<-vapply(filers, function(fff) dir.create(fff, showWarnings = FALSE),numeric(1)) ; rm(tmp)

filers<-c(paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc"),
          paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2015/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc"))

# Make some checks
if(!file.exists(filers[1])){
  
  print("Please download manually the Kummu GDP dataset from")
  print("https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0")
  print("then extract the file and make sure it can be found at the location")
  print(filers[1])
  print("For more information, please look at the `Installation` section of the README.md file on the Github landing page ")
  stop("GetEnv.R problem finding GDP data")
  
}

if(!file.exists(filers[1])){
  
  print("Please download manually the CIESIN population count dataset from")
  print("https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11/data-download")
  print("then extract the file and make sure it can be found at the location")
  print(filers[2])
  print("For more information, please look at the `Installation` section of the README.md file on the Github landing page ")
  stop("GetEnv.R problem finding CIESIN population data")
  
}