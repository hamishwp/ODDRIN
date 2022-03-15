
# ODDRIN 

### Oxford Disaster Displacement Real-time Information Network

![*Figure 1: ODDRIN interactive data visualisation platform - ODD-Mapping*](https://i.ibb.co/kx1mqtk/ODDRIN.png)

Welcome to the ODDRIN code, comprised of both a front-end visualisation component, known as ODD-Mapping, and the back-end statistical engine and real-time updating software, IIDIPUS - Integrated Internal DIsplaced PopUlation Sampler.

The aim of this software is to predict the number of people displaced and the number of buildings damaged in the early phases of rapid-onset natural (and human-generated climate-change related) hazards. Predictions are made as accurate as possible by training the model on hundreds of events and hundreds of thousands of damaged buildings, across a broad demographic of countries and hazard severities.

### Authors

ODDRIN was designed, developed, made operational and is managed by Dr. Hamish Patten at the Department of Statistics, University of Oxford [@HamishPatten](https://www.linkedin.com/in/hamish-patten), with support and formal supervision from Professor David Steinsaltz [@DavidSteinsaltz](https://www.steinsaltz.me.uk/).

## Documentation

[IDMC GRID 2021 Background Paper](https://www.internal-displacement.org/global-report/grid2021/downloads/background_papers/background_paper-analytics.pdf) - This was a non-peer reviewed article written for the Internal Displacement Monitoring Centre (IDMC) early 2021 to be included with the Global Report on Internal Displacement (GRID).

[arXiv](https://arxiv.org) - Unfortunately, the article intended for peer-review is not yet ready to be published, nor placed on arXiv. Watch this space!

___
## Code Layout

**The code can be decomposed into several sections:**

  1. **Main**
  2. **Model**
  3. **Method**
  4. **Data**
  5. **Object Orientated Programming (OOP) class formation**
  6. **Additional functions**

The majority of the programming was spent on number 4: extracting all the necessary data. Here we try to explain without too much detail the most important files from each section.

![*Figure 2: an example of ODDRIN output - a surface plot of the predicted displaced population from the Haiti earthquake on the 14/08/2021, including contour lines of the earthquake shakemap intensity (MMI)*](https://i.ibb.co/cvm25ZF/Exp-Disp-EQ20210814-HTI-10919.png)

### Main

##### Key Files & Roles:

1. `Main.R` - This is where the ODDRIN model parameterisation occurs, for each hazard. Here we extract the pre-formatted data, model formulas and structuring, and then run the model-training algorithm to parameterise the model.
2. `AutoQuake.R` - This file allows, requiring minimal input, an automated extraction of everything necessary to predict the spatial distribution and magnitude of the displaced population in the immediate aftermath of earthquakes on a population, including fore-shocks and after-shocks.
3. `RealTimeIIDIPUS.R` (not yet included) - Real-time tracking of the occurrence of rapid-onset hazards, including predicting the magnitude and spatial distribution of the displaced population, then broadcasting this to partners, such as the IFRC GO Platform.

##### Key Functions:

  - `IIDIPUSModelTraining` - Extract data, model and methodology, then train the model using a specific Adaptive Markov Chain Monte Carlo (MCMC). Note that initial values are given by the Nelder and Mead optimisation algorithm.
  - `AutoQuake` - Extracts the earthquake intensity data, creates an object that automatically extracts the relevant exposure and vulnerability information, then makes a prediction on the number of people estimated to be displaced, per gridpoint.

### Model

##### Key Files & Roles:

  1. `Model.R` - Here we can find everything model-related. This includes damage function equation definitions, declaring the chosen imported vulnerability indicators required, and, finally, the pseudo-marginal log-likelihood, prior and posterior distribution equations for population displacement and also satellite building damage estimations. Also includes the linear predictor terms that parameterise the systemic vulnerability.
  
##### Key Functions:

  - `logTarget` - The posterior distribution (which for likelihood-free pseudo-marginal approaches is known as the target distribution), this function tells you; given the observed number of people displaced and manually-classified satellite building damage assessment classifications, the Bayesian priors and the Approximate Bayesian Computing (ABC) rejection classification; how accurately the specified model parameterisation fits the data.
  - `LL_BD` - Log-likelihood calculation for the satellite image-based building damage assessment data
  - `LL_IDP` - Log-likelihood calculation for the observed population displacement data
  - `HighLevelPriors` - Approximate Bayesian Computing (ABC) method of rejection
  - `BinR` - Binary regression function used for the damage function of the model, which depends on the hazard type.
  - `llinpred` - Nationally-aggregated variable *linear predictor* (this refers to calculating the amount by which a variables contributes to systemic-vulnerability)
  - `GDPlinp` - Barely-sub-nationally-aggregated variable *linear predictor*, such as the Kummu, et al, GDP dataset. Speeds up the process.
  - `Plinpred` - Fully sub-nationally aggregated variable *linear predictor*, such as population density.

### Method

##### Key Files & Roles:

  1. `Method.R` - Define the Adaptive MCMC algorithm to be used to parameterise the model via Bayesian statistics. The specific algorithm is described under [Andrieu and Thoms, Stats and Computing, 2008](https://doi.org/10.1007/s11222-008-9110-y) under `Algorithm 4` in section 5.1.2. Furthermore, this file contains the routines used to parameterise the nationally aggregated systemic-vulnerability variables against the difference between the observed and predicted displacement values.
  2. `GetInitialValues.R` - This file both allows the user to input initial guess values for the model parameterisation and to use a standard optimisation technique to improve on this guess based on optimising the Maximum Likelihood Estimate (MLE). The optimisation algorithm is based on [Nelder and Mead, 1965](https://doi.org/10.1093/comjnl/7.4.308).
  3. `GDACS_Method.R` - Following the GDACS method of a basic Generalised Linear Model (GLM) which predicts the displaced population based on the number of people exposed to a hazard of at least a given intensity, for increasing hazard intensities, and includes a systemic vulnerability variable. The method is outlined on the [GDACS website](https://www.gdacs.org/Knowledge/models_eq.aspx).

##### Key Functions:

  - `Algorithm` - Adaptive MCMC algorithm
  - `SingleEventsModifierCalc` - Algorithm that calculates the `modifier` parameter that is used in this research to correlate the difference between the predicted displacement and the observed displacement with national indicators that infer systemic disaster vulnerability, e.g. population expansion or government corruption index. The value of this modifier parameter comes from a convex function, and therefore, we chose to use the bisection method to parameterise it based on MLE. This bisection method is described in the `OptimMd` function.
  - `Proposed2Physical` and `Physical2Proposed` - link functions for the model parameters, to ensure that the model parameters sampled by the MCMC proposal distribution are on the real line with infinite support.
  - `multvarNormProp` - proposal parameter set generation function (multivariate normal distribution)

### Data

##### Key Files & Roles:

  - `GetData.R` - In order to train the ODDRIN model parameters, the data must be extracted and unified. `GetData` does the following:
      - Takes displacement data from IDMC
      - Matches it to a specific hazard from the GDACS database
      - Finds if any satellite image-based building data is also on your computer
      - Extracts the hazard intensity data (data source depends on hazard)
      - Builds `ODD` and/or `BD` objects from all the available data
      - This therefore includes automated extraction of socio-economic and systemic vulnerability data
  - `GetPopDemo.R` - Population and demography data extraction, mostly built around the CIESIN data, but now includes the Facebook Data for Good Population Mapping data.
  - `GetSatDamage.R` - Given the location of UNOSAT-UNITAR or COPERNICUS building damage assessment data files, this extracts the buildings and harmonises the format of the data to be used by ODDRIN later on, when provided to initialise a `BD` object
  - `GetDisaster.R` - Extract the hazard intensity data, the source of which depends on the hazard type. For example, earthquakes rely on USGS
  - `GetGDACS.R` - This file is really hideous, I apologise, I was learning to code in R at the time. This file is to access the Global Disaster Alert and Coordination System (GDACS) database, this is key to the real-time component of ODDRIN
  - `GetUSGS.R` - Access earthquake shakemaps and other information automatically from the United States Geological Survey (USGS).
  - `GetSocioEconomic.R` - Extract socio-economic data, with the use of APIs and extracting from files that the user must download manually. This includes the aforementioned Kummu GDP dataset, and also the World Inequality Database set of income distribution indicators. Note that this is non-exhaustive in terms of all the socioeconomic data included in this research.
  - `GetINFORM.R` - Access the national indicator dataset produced by the Joint Research Centre (European Commission), including access to past values.
  - `GetWorldBank.R` - All things World Bank, as national aggregated values but with temporal trends. For example, it is easy to access temporal trends of population count for most countries around the world, and even have access to when the last time the data was updated via national surveys.
  - `GetOSM.R` - The file that accesses OpenStreetMaps, including downloading buildings and roads located within a certain bounding box, country or region polygon. Be careful what you wish for here... if your search is too broad, you'll never be able to access anything! Go in small chunks and slowly cover the area you want.

##### Key Functions:

  - `ExtractData` - This is the function that extracts all that we need for ODDRIN for a given hazard occurence. However, currently, this is only automated for earthquakes. Provided a collection of estimates of the **maximum** observed displaced population of an event (e.g. from IDMC), including the date of the event and the country, this function will find the matching value in GDACS, then create ODD classes, then do the same for the satellite image-based building assessment data.
  - `GetPopulationBbox` - Provided only a bounding box and the folder name for the population data, this function extracts the population data from CIESIN in a memory efficient way, also ensuring things like continuity across the longitude=0 plane.
  - `InterpPopWB` and `InterpGDPWB` - These functions extract the population and GDP nationally aggregated values, respecitively, from the World Bank exactly on the date provided (through interpolation/extrapolation techniques), which are used to make sure that the CIESIN and Kummu pop & GDP values are updated to reflect the value on the day of the hazard.
  - `ExtractBDfiles` - Provided the location of the folder where the satellite image-based building damage assessment data is kept, this function will extract all UNOSAT and all Copernicus data ready to create an instance of the `BD` class.
  - `GetDisaster` - This is the function that, provided with only minimal input (bounding box, start and end date, hazard type), can extract hazard intensity raster data, and output HAZARD objects made from the data.
  - `GetEarthquake` - Automated extraction for earthquakes from USGS, forming a list of `HAZARD` objects.  
  - `FilterGDACS` - Extract data from GDACS through their API, then filter it to get only what we need.
  - `GetUSGS` - For earthquakes, `GetDisaster` depends entirely on this function, which accesses the USGS database, given minimal input, and extracts the important data and forms a HAZARD class instance from it.
  
### Object Orientated Programming (OOP) class formation

##### Key Files & Roles:

  - `ODDobj.R` - The principle ODDRIN class, whereby hazard intensities (from all hazards included), exposed population and exposed buildings, as well as vulnerability information are all included as fields/attributes. The methods of the class greatly facilitate automating the initialisation of objects with only minimal data provided as input, whereby, for example, the interpolation of hazard intensities onto the population grid is automated.
  - `BDobj.R` - This class is for the satellite image-based building damage assessment data. The main difference from the `ODD` class is that this is that the data is not on a grid but can be considered as a list of points in space whereby the hazard intensities and other information is interpolated.
  - `HAZARDobj.r` - Hazard intensity data is read in and then structured into the correct form to be provided to the `ODD` or `BD` classes. This class is used heavily by `GetUSGS.R`, for example.

##### Key Functions:

  - `initialize` - Initialises the objects, with a unique initialisation function per object mentioned above
  - `DispX` - Predict the number of people displaced, based on the hazard, the model and the parameterisation
  - `BDX` - Predict the building damage classification level, based on the hazard, the model and the parameterisation
  - `BDinterpODD` - Creating an instance of the `BD` class is facilitated by providing the instance of the `ODD` class that corresponds with the same hazard(s).

### Additional Functions

##### Key Files & Roles:
  
  - `Functions.R` - Includes all the miscellaneous functions that are required by the ODDRIN code.
  
##### Key Functions:

  - `convRaster2SPDF` - Converts rasters that are imported from a `.tif` file into the `SpatialPixelsDataFrame` format which ODDRIN relies on.
  - `convMat2SPDF` - Same but from matrix format to `SpatialPixelsDataFrame` format
  - `coords2country` Given a longitude and latitude coordinate, which country does it belong to?
  - `countriesbbox` Given the ISO3C code for a given country, what is the countries bounding box?

___________

## Installation

Before doing anything, please change the directory location environment variable `dir` and `directory` (these two must be equal) in the `GetEnv.R` file. The simplest installation of ODDRIN is to download and load only the most fundamentally important packages and to source only the files that you will need. To do this, change `packred<-T` in the file `GetEnv.R`. Installation in RStudio is simple, just run the following:

```R
  
  # Extract Environment Variables
  source('RCode/GetEnv.R')
  # Download and install the necessary packages:
  source('RCode/GetODDPackages.R')  


```

**This will error if you have not already installed the necessary data** (see the end of installation instructions - *'Data to be downloaded manually'*). In order to install the full ODDRIN package, I can only guide you through Linux and possibly Mac distributions (sorry!). For Macs, this should, in theory, be the same, but I can't say for sure. The problem is getting the R package `rJava` to work. For Linux and Mac distributions, open a terminal and run the following:

```bash

    sudo apt-get install libcurl4-openssl-dev libxml2-dev libjq-dev libprotobuf-dev
                         libv8-dev protobuf-compiler openjdk-8-jdk libssh-dev libssl-dev
                         libgdal-dev libudunits2-dev
                         
```

This installs all sorts of important software, not just `rJava`, but this one is what causes the problems - note the `openjdk-8-jdk` package is the difficult one. Next part is to make sure that you have an enviroment variable for the location of your java libraries. In your `/etc/environment` file add the java libraries environment variable (using `sudo nano /etc/environment`):
    
```bash

    LD_LIBRARY_PATH=/usr/lib/jvm/java-8-openjdk-amd64/lib/amd64/:/usr/lib/jvm/java-8-openjdk-amd64/lib/amd64/
    
```

Please check that the folder `/usr/lib/jvm/java-8-openjdk-amd64/lib/amd64/` actually exists! Otherwise, insert the folder location you find (another example could be `LD_LIBRARY_PATH=/usr/lib/jvm/jre/lib/amd64:/usr/lib/jvm/jre/lib/amd64/default`). *NOW RESTART YOUR COMPUTER!* Follow this up with:

```bash

    source /etc/environment
    sudo R CMD javareconf JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/bin/jar
    
```

Finally, with `packred<-F` in the file `GetEnv.R`, run the following:

```R
  
  # Extract Environment Variables
  source('RCode/GetEnv.R')
  # Download and install the necessary packages:
  source('RCode/GetODDPackages.R')  
  

```

##### Data to be Downloaded Manually

Furthermore, you will also be required to download the following datasets:

  - Download the [sub-national GDP dataset](https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0) - Taken from the article by [Kummu, et al, 2020](https://www.nature.com/articles/sdata20184). Once downloaded, place the file called `GDP_per_capita_PPP_1990_2015_v2.nc` in the folder `Demography_Data/SocioEconomic/KUMMU/` which you must create inside the ODDRIN folder such that, from the console, `file.exists(paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc"))==TRUE`.
  - Then download the [high resolution population count dataset](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11/data-download) - Taken from the article produced by [Center for International Earth Science Information Network - CIESIN](https://doi.org/10.7927/H4F47M65). Please download **all years via the *Single Year* option, in ASCII format, at a resolution of 30 arc-seconds** as shown in the screenshot below. From what I remember, you may need to download each year at a time, otherwise the CIESIN servers hang. In all cases, you need the 2015 file as this is used to normalise the ODDRIN population values. The ASCII format is to ensure quick computation, and also to ensure that no problems are found when crossing the longitude=0 plane. Once downloaded, place the files in folders separately split by year, such that folders such as `paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2015/")` or `paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2005/")` exist. The files should then exist inside the correct folder under, for example, the name `gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc` (for the 2015 year). Please check that there should be 8 files in each folder finishing with the same extension `.asc` but as `...30_sec_2.asc`, `...30_sec_3.asc`, `...30_sec_4.asc` etc.
  - Download the [EM-DAT database](https://public.emdat.be/) using the query tool and selecting `Natural -> Geophysical -> Earthquake` as the Disaster Classification and filtering from 2008 - 2021. Name the file `emdat.xlsx` and save in the folder `Displacement_Data` such that the folder `paste0(dir,"Displacement Data/emdat.xlsx")` exists. 
  
![*Figure 3: screenshot of the CIESIN website download options for the population count data used in the ODDRIN software*](https://i.ibb.co/723h9BX/Screenshot-from-2022-01-05-16-20-30.png)

## Environment Variables

To run this software, you will need to add the following environment variables to the `GetEnv.R` file:

`directory` = `dir` (for the lazy writers), `FBdirectory`, `packred`

Note that `FBdirectory` is optional, but is important when using data extracted from the Facebook Data for Good platform.

## Installation Checks

Please run the `InstallationChecks.R` file:

```R

  source('RCode/InstallationChecks.R')
  

```


## Usage/Examples

In the folder `IIDIPUS_Input`, there are a few example files, which contain the three different objects used by ODDRIN: `HAZARD`, `ODD` and `BD`, corresponding to the object containing hazard information and raster data, then the principal ODDRIN object that is used to predict population displacement, and the building damage object, respectively. You can have a basic explore using the following:

```R

# Extract Environment Variables
source('RCode/GetEnv.R')
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
# Extract model functions and priors
source('RCode/Model.R')
# Extract the model parameterisation algorithm, default = Adaptive MCMC
source('RCode/Method.R')

# This file holds all the required data - hazard, exposure, vulnerability, as well as event information and observed displacement estimates in the 'ODD' class format
ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/EQ20210814HTI_10919_example"))
# This is the model parameterisation, currently trained only on earthquakes on a global level
Omega<-readRDS(paste0(dir,"IIDIPUS_Results/Omega_v2_20210828.Rdata"))
# Test to see if the displacement prediction calculations are working
ODDy%<>%DispX(Omega = Omega,center = Model$center, BD_params = Model$BD_params, LL=F,Method = AlgoParams)
# Plot the ODD object using base R functions:
plot(ODDy) # default plots the CIESIN population data
plot(ODDy["hazMean1"]) # plots the principle hazard intensity
# Or you can use some of the ODD class methods
p<-MakeODDPlots(ODDy) # plots hazard and population side-by-side
p<-plotODDy(ODDy) # plots hazard contour lines ontop of displaced population surface plot
# Then we can also tune the plot to our greatest desires
p<-plotODDy(ODDy,breakings = c(0,10,50,100,500,1000,5000,10000),bbox=c(-74.5,17.9,-72.5,19),zoomy = 9)


```

Additionally, you can try extracting earthquake information directly from USGS through their API via the in-built ODDRIN functions written up in `AutoQuake.R`, where an example is given. Note that the USGS ID of an earthquake can be seen in the website name, e.g. [https://earthquake.usgs.gov/earthquakes/eventpage/at00qxtxcn/executive](https://earthquake.usgs.gov/earthquakes/eventpage/at00qxtxcn/executive) is for the HTI 2021 earthquake, and the USGS ID is `at00qxtxcn`.


```R

# Download and install the necessary packages, and load source files & data:
source('RCode/GetODDPackages.R')
# Search for an earthquake by providing a date range and longitude-latitude bounding box
input<-list(
  sdate=as.Date("2019-12-13"), # "YYYY-MM-DD"
  fdate=as.Date("2019-12-17"), # "YYYY-MM-DD"
  iso3="PHL", # Country code in ISO-3C form
  datadir=dir, # Location of the main folder to access the data 
  plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)
# Or extract the data purely based on the USGS id number
input<-list(USGSid="at00qxtxcn",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)
#%%%%%%%%%%%%% Variables and functions from ODDRIN files %%%%%%%%%%%%%%%#
input%<>%append(list(Model=Model, # Model parameters - equations, parameters, ... (Model.R)
                     Omega=Omega, # Parameterisation of the model (GetInitialValues.R)
                     Method=AlgoParams)) # Number of CPUs, particles in SMC, etc. (Method.R)

ODDy<-AutoQuake(input)


```

To run the model parameter training algorithm, use the `Main.R` file:


```R

source('RCode/Main.R')
output<-IIDIPUSModelTraining()


```

