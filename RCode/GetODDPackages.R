
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
                        "GADMTools","akima","adehabitatMA","flexsurv")
  
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