# Where is the main folder with all the code and data
dir<-directory<-"/home/patten/Documents/Coding/Oxford/IIDIPUS/"
# Set the working directory from your environment variables
setwd(directory)
# Directory of the Data for Good data, e.g. Disaster Mapping, 4G connectivity, etc
FBdirectory<-'/home/patten/Documents/IDMC/Facebook_Data/'
# Do you want only the reduced packages or all? Choose via packred
packred<-F

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

