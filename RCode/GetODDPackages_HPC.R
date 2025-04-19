

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
    source('RCode/GetSubNationalData.R')
    source('RCode/AggregateODD.R')
    
    # IDP estimate related:
    # source('RCode/GetDisplacements.R')
    # source('RCode/GetHelix.R')
    
    # Demography & population related:
    source('RCode/AddVulnerability.R')
    source('RCode/GetSocioEconomic.R')
    source('RCode/GetWorldPop.R')
    source('RCode/GetPopDemo.R')
    #source('RCode/GetINFORM.R')
    #source('RCode/GetWorldBank.R')
    
    # Damage estimate related:
    source('RCode/GetBuildingCounts.R')
    source('RCode/GetSatDamage.R') 
    #source('RCode/GetOSM.R')
    
    #Simulated data Related:
    source('RCode/Simulate.R')
    
    #Results Analysis:
    source('RCode/ResultsAnalysisFunctions.R')
    
  }
  
}

LoadLibraries<-function(packred, loadRmpi=F){
  
  library(tidyverse)
  library(magrittr)
  #library(ggmap)
  library(pracma)
  library(parallel)
  library(mvtnorm)
  library(abind)
  library(countrycode)
  library(raster)
  library(scoringRules)
  library(viridis)
  library(gridExtra)
  library(devtools)
  library(shiny)
  library(usethis)
  library(LaplacesDemon)
  library(cowplot)
  library(terra) 
  library(DescTools)
  library(gstat)
  if (loadRmpi){library(Rmpi)}
  
  if(!packred) {
    library(sf)  
    library(geojsonR)
    library(geosphere)
    library(rworldmap)
    library(lutz)
    library(xml2)
    library(gstat)
    library(wid)
    library(geojsonio)
    library(jsonlite)
    library(reshape2)
    library(ecochange)
    library(openxlsx)
    library(osmdata)
    library(tidyxl)
    
    # library(doParallel)
    # library(foreach)
    # library(EnvStats)
    # library(posterior)
    # library(VGAM)
    # library(DescTools)
    # library(codetools)
    # library(OpenStreetMap)
    # library(ecochange)
    # library(geodata)
    # library(spatialEco)
    # library(ExtDist)
  }
}

GetODDPackages<-function(packred, loadRmpi=F){
  
  # list.of.packages <- c("dplyr", "ggplot2","sf","tidyverse","openxlsx","pracma",
  #                        "tiff", "gstat", "mvtnorm", 'DescTools', #'rgeos', 
  #                       "RColorBrewer", "geosphere","GGally", "wbstats",
  #                       "countrycode","rworldmap","rworldxtra","chron","ncdf4",
  #                       "GADMTools","akima","adehabitatMA","flexsurv", "ExtDist", 
  #                       'EnvStats', 'posterior', 'doParallel', 'VGAM', 'abind',
  #                       'Rmpi', 'openxlsx', 'ecochange', 'lutz', 'scoringRules',
  #                       'vegan')
  
  list.of.packages <- c('tidyverse', 'magrittr', 'pracma', 'parallel',
                        'mvtnorm', 'abind', 'countrycode', 'raster', 'scoringRules',
                        'viridis', 'gridExtra', 'devtools', 'shiny', 'usethis',
                        'LaplacesDemon', 'cowplot', 'DescTools', 'gstat')
  if (loadRmpi){ list.of.packages <- c(list.of.packages,'Rmpi')}
  
  if(!packred) list.of.packages<-c(list.of.packages,
                                   'sf', 'geojsonR', 'geosphere', 'rworldmap', 'lutz',
                                   'xml2', 'gstat', 'terra', 'geojsonio', 'jsonlite', 'reshape2',
                                   'ecochange', 'openxlsx', 'osmdata', 'tidyxl')
                                   # "codetools","latex2exp", "geojsonR", 'geojsonio',
                                   # "rJava","devtools","OpenStreetMap","osmdata",
                                   # "tidyRSS","geojsonR", "tiff", "gstat",
                                   # "FactoMineR","factoextra","xtable",
                                   # "gsubfn","mapsapi","leaflet", "ssh","RPostgres",
                                   # "GADMTools", "pscl","multiColl", 'lutz', 'ecochange',
                                   # 'geodata', 'spatialEco')
  
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) install.packages(new.packages, repos='http://cran.us.r-project.org')
  
  # This package makes sure you're uptodate on all packages, including things like tidyverse 
  #if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos='http://cran.us.r-project.org')
  #BiocManager::install("Biobase", version = "3.19")
  
  # devtools::install_github('daroczig/fbRads')
  # if(length(list.of.packages[!("openrouteservice" %in% installed.packages()[,"Package"])])){devtools::install_github("rCarto/osrm")}
  #if(length(list.of.packages[!("ggmap" %in% installed.packages()[,"Package"])])){devtools::install_github("dkahle/ggmap")}
  #if(length(list.of.packages[!("ggmap" %in% installed.packages()[,"Package"])])){devtools::install_github("stadiamaps/ggmap")}
  # if(length(list.of.packages[!("countrycodes" %in% installed.packages()[,"Package"])])){devtools::install_github("vincentarelbundock/countrycode")}
  #if(length(list.of.packages[!("wbstats" %in% installed.packages()[,"Package"])])){devtools::install_github('nset-ornl/wbstats')}
  #if(length(list.of.packages[!("parallelsugar" %in% installed.packages()[,"Package"])])){devtools::install_github('nathanvan/parallelsugar')}

  if (!packred){
    if(length(list.of.packages[!("wid" %in% installed.packages()[,"Package"])])){devtools::install_github("WIDworld/wid-r-tool")}
  }
  
  LoadLibraries(packred, loadRmpi)
  GetSourceFiles(packred)
  
  #see what packages different files are using
  #NCmisc::list.functions.in.file(filename = paste0(dir,"RCode/GetSatDamage.R"), alphabetic = FALSE) |> print()
}

GetODDPackages(packred, loadRmpi)

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

#Check that required data is downloaded:

if (!packred){
  #Check that required data files have been downloaded:
  
  # #KUMMU GDP Data (currently not being used):
  # if(!file.exists(paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc"))){
  #   print("Please download manually the Kummu GDP dataset from")
  #   print("https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0")
  #   print("then extract the file and make sure it can be found at the location")
  #   print(filers[1])
  #   print("For more information, please look at the `Installation` section of the README.md file on the Github landing page ")
  #   stop("GetEnv.R problem finding GDP data")
  # }
  
  #CIESIN Population Data:
  #even though we are now using WorldPop rather than CIESIN, we are still using the nation identifiers from the CIESIN data (national_identifier_grid files)
  if(!file.exists(paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2015/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc"))){
    print("Please download manually the CIESIN population count dataset from")
    print("https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11/data-download")
    print("then extract the file and make sure it can be found at the location")
    print(filers[2])
    print("For more information, please look at the `Installation` section of the README.md file on the Github landing page ")
    stop("GetEnv.R problem finding CIESIN population data")
  }
  
  #GDL data:
  if(!file.exists(paste0(dir,"Demography_Data/SocioEconomic/GlobalDataLab/GDL Shapefiles V6/shdi2022_World_large.shp"))){
    print("Please download manually the vulnerability data shapefiles here (you'll have to create a free account with GDL first):")
    print("https://globaldatalab.org/asset/378/GDL%20Shapefiles%20V7.zip")
    print("then extract the zipped files and place them in: ./Demography_Data/SocioEconomic/GlobalDataLab/GDL Shapefiles V6/")
    print("For more information, please look at the `Installation` section of the README.md file on the Github landing page ")
    stop("GetODDPackages.R problem finding Global Data Lab Vulnerability data")
  }
  if(!file.exists(paste0(dir,"Demography_Data/SocioEconomic/GlobalDataLab/SHDI-SGDI-Total 7.0.csv"))){
    print("Please download manually the vulnerability data here (you'll have to create an account with GDL first):")
    print("https://globaldatalab.org/asset/348/SHDI-SGDI-Total%207.0.csv")
    print("then name the file 'SHDI-SGDI-Total 7.0.csv' and place here: ./Demography_Data/SocioEconomic/GlobalDataLab/")
    print("For more information, please look at the `Installation` section of the README.md file on the Github landing page ")
    stop("GetODDPackages.R problem finding Global Data Lab Vulnerability data")
  }
  
  #Vs30 data
  if(!file.exists(paste0(dir,"Hazard_Data/global_vs30_tif/global_vs30.tif"))) {
    print("Please download the VS30 dataset (geotiff and auxiliary files) from here https://www.usgs.gov/programs/earthquake-hazards/science/vs30-models-and-data")
    print('Then place in ./Hazard_Data/global_vs30_tif/')
    stop("GetODDPackages.R problem finding Vs30 data")
  }
  
  #PGA data:
  # if(!file.exists(paste0(dir,"Hazard_Data/gdpga/pga_475y.tif"))){
  #   print("Please download the hazard frequency data from here: https://www.geonode-gfdrrlab.org/layers/hazard:pga_475y")
  #   print("Select Download Layer -> Data -> Original Dataset")
  #   print("Then extract the file pga_475y.tif and move to ./Hazard_Data/gdpga/")
  #   stop("GetODDPackages.R problem finding Earthquake Frequency data")
  # } 
  if (!file.exists(paste0(dir, 'Hazard_Data/GEM-GSHM_PGA-475y-rock_v2023/v2023_1_pga_475_rock_3min.tif'))){
      print("Please download the hazard frequency data from here: https://www.globalquakemodel.org/product/global-seismic-hazard-map") #
      print('Extract the zip file and move to ./HazardData/')
      stop('GetODDPackages.R problem finding Earthquake Frequency data')
  }
}
