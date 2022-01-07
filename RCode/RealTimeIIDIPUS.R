#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%% IIDIPUS - Integrated Internal DIsplaced PopUlation Sampler %%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%% Coded by Hamish Patten %%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%% Collaboration between the University of Oxford %%%%%%%%%%#
#%%%%%%%%%% and the Internal Displacement Monitoring Centre %%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%% Started January 2020 %%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Extract Environment Variables
source('RCode/GetEnv.R')
# Set the working directory from your environment variables
setwd(directory)
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')

#@@@@@ INPUT @@@@@#
# INDEX - TYPE - NAME - DESCRIPTION
# 1:4 - numeric   - bbox[min lon, max lat, max lon, min lat] - region bounding box
# 5   - character - country - ISO format only
# 6   - character - hazard  - using IDMC definition: {"Flood","Storm","Mass movement","Wildfire","Earthquake","Extreme temperature","Volcanic eruption","Drought"}
# 7   - date (%Y%m%d) - sdate - start date of event
# 8   - date (%Y%m%d) - fdate - start date of event
args <- commandArgs(trailingOnly = TRUE)
args<-CheckArgs(args)
bbox<-args$bbox
iso3<-args$iso3
hazard_type<-args$hazard_type
sdate<-args$sdate
fdate<-args$fdate
year<-args$year

###### PAST DATA FOR STATISTICAL INFERENCE #####
if(length(list.files(paste0(directory,"Statistical_Inference/*")))==0) {
  print("WARNING: IIDIPUS is going to recalculate all past data for statistical inference... this will take a while")
  IIDIPUS_Past(directory)
}
filer <- file.info(list.files(paste0(directory,"Statistical_Inference"), full.names = T))
# Get most recent file
filer <- rownames(filer)[which.max(filer$mtime)]

############ DISASTER ###########
# OUTPUT: List of Functions - 2D Cubic Spline (long,lat) FOR EACH DAY
disaster<-GetDisaster()

############ IDP COUNT ############
# OUTPUT: List - Place Name, polygon or point, type (polygon or point), estimated IDP per day, displaced-to/from location
IDP_FB<-GetFBEstimates()
IDP_tot<-GetTotalEstimates()
IDP_shelter<-GetShelterEstimates()
IDP_CDR<-GetCDREstimates()

############ MOBILITY ############
# OUTPUT: Function - 2D Cubic Spline (long,lat)
SmobF<-GetOSMMobility()
BmobF<-GetOSMMobility()

############ CONNECTIVITY ############
# OUTPUT:
NetCov<-GetNetworkCoverage(directory,fnamer)
FBcon<-GetFBconnectivity()

############ POPULATION ############
# OUTPUT: Monte Carlo sample of population, including demography: elderly & income
population<-GetPopulation() %>% GetDemography()
# WBiso<-GetWBinfo(directory)

############ MODEL ############
### PARAMETERS: ###
# Mobility: Force effect like gravity F = M(1) * M(2) * G / D  +  connectivity  +  
# Where D is some function of distance via road mapping and time to target via public transport and car
# M is the population density+road density and socio-economic indicator of each of two places
# connectivity is FB connectivity between two points
# Then convert from 


#grid mapping of likelihood of displacement to (long,lat) coordinate. 
# The probability is calculated based on population+road environement
# Then the probability of being displaced to area is exponential decay


# Pmob_g, Pmob_i - physical mobililty as a continuous function of income per age grouping (old, not old)
# Smob_g, Smob_i - length 2 : physical mobililty as a continuous function of income, per age grouping (old, not old)


########################## Code Outline ##############################
# Likelihood of being displaced to shelter
# Likelihood of returning based on a) displaced location (shelter?) b) risk of further/continuing hazard and c) magnitude of damage

# How to integrate shelter, global, FB and Flowkit stock estimates together?
# Flowkit and FB should be easy: complementary for given region. FB has surveys too.
# Shelter links into Flowkit & FB via social connectivity levels, and evacuation period.
# 

# Calculate evacuation period: disaster info, twitter & insta posts.
# Find population demographic
# Physical mobility and length-scales: how far to 'safety' (safe town or shelter)
# 

# Correlate FB RWI and GDP to social network connectivity to provide an estimate for rural areas.

# Population demographic
#        - Ageing populations (SEDACS)
#        - Income (FB RWI)
#        - GDP PPP (Matti KUMMU)
# Create mobility mapping using 
#        - FB social network index 
#        - Distance & time to shelter (mapsapi)
#        - Distance & time to local city (second derivatives of SEDAC population density)
# Disaster risk zone
#        - Damage & magnitude 2D mapping (GDACS)
#        - Building damage (IOM???)
#        - Social media presence (instaR & Rtweet, searching for translated event e.g. 'earthquake' in Philipinnes language)
#          (THIS COULD BE GLOBAL SEARCH AND NOT LOCATED TO THE REGION) - 
#          using Hawkes processes for return, depending on social mobility value of individual.
#          If significantly more insta posts & tweets (& retweets) referring to disaster after event occurance are produced after the evacuation period or more than previous days, search for words like 'safe' which will imply a time-to-return.



# PEOPLE ALSO DISPLACED TO HOTELS AND NOT JUST TO SOCIAL CONNECTIONS



# Main variables in model:
#        - evacuation probability
#        - likelihood of evacuation by car or public transport
#        - likelihood to have a car? (Search OSM for nearest garages?)
#        - probability to have home destroyed or damaged (from IOM, etc, data, or train networks...)
#        - end-of-disaster date
#        - 
#        - probability of long-term displacement
#        - probability of evacuating to shelter
#        - probability of being undetected (displacement to family/friend w/o FB)
#        - 
#        - 


#@@@@@@@@@@@@@@@@@@@@@@@ IDP Estimate Data @@@@@@@@@@@@@@@@@@@@@@@@@@#
############### Government/News/NGO global estimates #################
# - Bayes sample values from Gamma dist with SD
######################### Shelter estimates ##########################

########################### FB estimates #############################

###################### FlowKit/CDR Estimates #########################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


#@@@@@@@@@@@@@@@@@@@@@@@@@ Disaster Data @@@@@@@@@@@@@@@@@@@@@@@@@@@@#
############################## GDACS #################################
#
# 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


#@@@@@@@@@@@@@@@@@@@@@@@ Demographic Data @@@@@@@@@@@@@@@@@@@@@@@@@@@#
#################### Ageing populations - SEDAC ######################
#
######################### Income - FB RWI ############################
#
##################### Income - KUMMU GDP PPP #########################
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


#@@@@@@@@@@@@@@@@@@@@@@@ Social Media Data @@@@@@@@@@@@@@@@@@@@@@@@@@#
# Increase in social media use inside the disaster bounding box:
# Twitter, FB, etc as a function of days after event.
############################ Twitter #################################

########################### Instagram ################################
# instaR pacakge: 'searchInstagram' function, including mindate & maxdate.
# First, search hashtags e.g. 'earthquake'
# Second, search all posts close to certain location
#
# Output parameters: 
#     - posts per population density
#     - Number of posts out of the norm
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


