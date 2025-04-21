
# ODDRIN

### Oxford Disaster Displacement Real-time Information Network

![*Figure 1: ODDRIN interactive data visualisation platform - ODD-Mapping*](https://i.ibb.co/kx1mqtk/ODDRIN.png)

Welcome to the ODDRIN code, comprised of both a front-end visualisation component, known as ODD-Mapping, and the back-end statistical engine and real-time updating software, IIDIPUS - Integrated Internal DIsplaced PopUlation Sampler.

The aim of this software is to predict the number of people displaced, the number of fatalities, and the number of buildings damaged in the early phases of rapid-onset natural (and human-generated climate-change related) hazards. Predictions are made as accurate as possible by training the model on hundreds of events and hundreds of thousands of damaged buildings, across a broad demographic of countries and hazard severities.

### Authors

ODDRIN was designed, developed, made operational and is managed by Dr. Hamish Patten [@HamishPatten](https://www.linkedin.com/in/hamishpatten/),
Max Anderson Loake [@MaxLoake](www.linkedin.com/in/max-anderson-loake-943045195)
and Professor David Steinsaltz [@DavidSteinsaltz](https://www.stats.ox.ac.uk/people/david-steinsaltz), as part of a project developed at the Department of Statistics, University of Oxford.

## Documentation
[A Bayesian Approach to Disaster Impact Modelling](https://arxiv.org/abs/2412.15791) - This pre-print, submitted to the Royal Statistical Society (RSS), details the model, method, and results.

[Data-Driven Earthquake Multi-impact Modeling: A Comparison of Models](https://link.springer.com/article/10.1007/s13753-024-00567-5) - Published in the International Journal of Disaster Risk Science, this paper compares a range of machine-learning approaches to earthquake impact modelling.

[IDMC GRID 2021 Background Paper](https://www.internal-displacement.org/global-report/grid2021/downloads/background_papers/background_paper-analytics.pdf) - This was a non-peer reviewed article written for the Internal Displacement Monitoring Centre (IDMC) early 2021 to be included with the Global Report on Internal Displacement (GRID).
___
## Code Layout

**The code can be decomposed into several sections:**

  1. **Main**
  2. **Model**
  3. **Method**
  4. **Data**
  5. **Object Orientated Programming (OOP) class formation**
  6. **Additional functions**

Here we try to explain without too much detail the most important files from each section.

![*Figure 2: an example of ODDRIN output - a surface plot of the predicted displaced population from the Haiti earthquake on the 14/08/2021, including contour lines of the earthquake shakemap intensity (MMI)*](https://i.ibb.co/cvm25ZF/Exp-Disp-EQ20210814-HTI-10919.png)

### Main

##### Key Files & Roles:

1. `Main.R` - This is where the ODDRIN model parameterisation occurs. Here we extract the pre-formatted data, model formulas and structuring, and then run the model-training algorithm to parameterise the model.
2. `AutoQuake.R` - This file allows, requiring minimal input, an automated extraction of everything necessary to predict the spatial distribution and magnitude of the mortality, displaced population, and building damage in the immediate aftermath of earthquakes, including fore-shocks and after-shocks.
3. `RealTimeIIDIPUS.R` (not yet included) - Real-time tracking of the occurrence of rapid-onset hazards, including predicting the magnitude and spatial distribution of the displaced population, then broadcasting this to partners, such as the IFRC GO Platform.

##### Key Functions:

  - `IIDIPUSModelTraining` - Extract data, model and methodology, then train the model using an Adaptive Markov Chain Monte Carlo (AMCMC) algorithm or Sequential Monte Carlo (SMC) algorithm.
  - `AutoQuake` - Extracts the earthquake intensity data, creates an object that automatically extracts the relevant exposure and vulnerability information, then makes a prediction on the fatalities, population displacement and building damage, per gridpoint.

### Model

##### Key Files & Roles:

  1. `Model.R` - Here we can find everything model-related. This includes damage function equation definitions, declaring the chosen imported vulnerability indicators required, and, finally, the pseudo-marginal log-likelihood, prior and posterior distribution equations for population displacement and also satellite building damage estimations. Also includes the linear predictor terms that parameterise the systemic vulnerability.

##### Key Functions:
- `HighLevelPriors` - Approximate Bayesian Computing (ABC) method of rejection
- `GetLP` - Calculate the exposure-related component of the vulnerability over all grid-cells (e.g. using the SHDI data)
- `GetLP_single` - Calculate the exposure-related component of the vulnerability for a single grid cell
- `getLP_event` - Calculate the hazard-related component of the vulnerability (e.g. using the night time indicator)
- `addTransfParams` - Transforms parameters to reduce correlation between parameters
- `SamplePolyImpact` - Sample the impact for each event in the provided event set
- `SamplePointImpact` - Sample the impact for each building in the point building dataset
- `CalcDist` - Calculates the loss function comparing the sampled and observed data

<!--  - `logTarget` - The posterior distribution (which for likelihood-free pseudo-marginal approaches is known as the target distribution), this function tells you; given the observed number of people displaced and manually-classified satellite building damage assessment classifications, the Bayesian priors and the Approximate Bayesian Computing (ABC) rejection classification; how accurately the specified model parameterisation fits the data.
  - `LL_BD` - Log-likelihood calculation for the satellite image-based building damage assessment data
  - `LL_IDP` - Log-likelihood calculation for the observed population displacement data
  - `HighLevelPriors` - Approximate Bayesian Computing (ABC) method of rejection
  - `BinR` - Binary regression function used for the damage function of the model, which depends on the hazard type.
  - `llinpred` - Nationally-aggregated variable *linear predictor* (this refers to calculating the amount by which a variables contributes to systemic-vulnerability)
  - `GDPlinp` - Barely-sub-nationally-aggregated variable *linear predictor*, such as the Kummu, et al, GDP dataset. Speeds up the process.
  - `Plinpred` - Fully sub-nationally aggregated variable *linear predictor*, such as population density. -->

### Method

##### Key Files & Roles:

  1. `Method.R` - Define the two algorithms that are used to parameterise the model via likelihood-free Bayesian statistics. The options are the Adaptive MCMC algorithm described in [Del Moral, Doucet and Jasra, 2012](https://link.springer.com/article/10.1007/s11222-011-9271-y), and the ABC-SMC algorithm described in [Spencer, 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12344).
  2. `GetInitialValues.R` - This file allows the initialisation of the AMCMC algorithm, either by using samples from past model runs or samples from the prior to estimate an appropriate proposal covariance.
<!--This file both allows the user to input initial guess values for the model parameterisation and to use a standard optimisation technique to improve on this guess based on optimising the Maximum Likelihood Estimate (MLE). The optimisation algorithm is based on [Nelder and Mead, 1965](https://doi.org/10.1093/comjnl/7.4.308).
  3. `GDACS_Method.R` - Following the GDACS method of a basic Generalised Linear Model (GLM) which predicts the displaced population based on the number of people exposed to a hazard of at least a given intensity, for increasing hazard intensities, and includes a systemic vulnerability variable. The method is outlined on the [GDACS website](https://www.gdacs.org/Knowledge/models_eq.aspx). -->

##### Key Functions:

<!--- - `Algorithm` - Adaptive MCMC algorithm
  - `SingleEventsModifierCalc` - Algorithm that calculates the `modifier` parameter that is used in this research to correlate the difference between the predicted displacement and the observed displacement with national indicators that infer systemic disaster vulnerability, e.g. population expansion or government corruption index. The value of this modifier parameter comes from a convex function, and therefore, we chose to use the bisection method to parameterise it based on MLE. This bisection method is described in the `OptimMd` function. -->
  - `Proposed2Physical` and `Physical2Proposed` - link functions for the model parameters, to ensure that the model parameters sampled by the MCMC proposal distribution are on the real line with infinite support.
  - `multvarNormProp` - proposal parameter set generation function (multivariate normal distribution)
  - `AMCMC` - Runs the adaptive MCMC algorithm
  - `ABCSMC` - Runs the ABC-SMC algorithm

### Data

##### Key Files & Roles:

  <!-- - `GetData.R` - In order to train the ODDRIN model parameters, the data must be extracted and unified. `GetData` does the following:
      - Takes displacement data from IDMC
      - Matches it to a specific hazard from the GDACS database
      - Finds if any satellite image-based building data is also on your computer
      - Extracts the hazard intensity data (data source depends on hazard)
      - Builds `ODD` and/or `BD` objects from all the available data
      - This therefore includes automated extraction of socio-economic and systemic vulnerability data -->

  - `GetPopDemo.R` - Population and demography data extraction, mostly built around the CIESIN data, but now includes the Facebook Data for Good Population Mapping data.
  - `GetSatDamage.R` - Given the location of UNOSAT-UNITAR or COPERNICUS building damage assessment data files, this extracts the buildings and harmonises the format of the data to be used by ODDRIN later on, when provided to initialise a `BD` object
  - `GetDisaster.R` - Extract the hazard intensity data, the source of which depends on the hazard type. For example, earthquakes rely on USGS
  - `GetGDACS.R` - This file is really hideous, I apologise, I was learning to code in R at the time. This file is to access the Global Disaster Alert and Coordination System (GDACS) database, this is key to the real-time component of ODDRIN
  - `GetUSGS.R` - Access earthquake shakemaps and other information automatically from the United States Geological Survey (USGS).
  - `AddVulnerability.R` - Extract and add the vulnerability data related to the exposure, such as the Vs30 data (from USGS) and the SHDI data (from Global Data Lab)
  - `GetWorldBank.R` - All things World Bank, as national aggregated values but with temporal trends. For example, it is easy to access temporal trends of population count for most countries around the world, and even have access to when the last time the data was updated via national surveys.
  - `GetOSM.R` - The file that accesses OpenStreetMaps, including downloading buildings and roads located within a certain bounding box, country or region polygon. Be careful what you wish for here... if your search is too broad, you'll never be able to access anything! Go in small chunks and slowly cover the area you want.
  - `GetBuildingCounts.R` - The file that accesses building footprint data from Microsoft/Bing Building Footprint datasets.

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
  - `HAZARDobj.R` - Hazard intensity data is read in and then structured into the correct form to be provided to the `ODD` or `BD` classes. This class is used heavily by `GetUSGS.R`, for example.

##### Key Functions:

  - `initialize` - Initialises the objects, with a unique initialisation function per object mentioned above
  - `DispX` - Predict the number of people displaced, based on the hazard, the model and the parameterisation
  - `readODD` - read in a saved ODD file (stored as RDS). For BD objects and HAZARD objects, the relevant functions are readBD() and readHAZ() respectively.
  - `saveODD` - save an ODD file (stored as RDS). For BD objects and HAZARD objects, the relevant functions are saveBD() and saveHAZ() respectively.
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

**This will error if you have not already installed the necessary data** (see the end of installation instructions - *'Data to be downloaded manually'*).

For a full installation of ODDRIN, the problem is getting the R package `rJava` to work. For Linux and Mac distributions, follow the installation instructions in *'Linux and Mac Installation'* below. For Windows, follow the instructions under *'Windows Installation'*.

##### Linux and Mac Installation

In order to install the full ODDRIN package for Linux and Mac distributions, open a terminal and run the following:

```bash

    sudo apt-get install libcurl4-openssl-dev libxml2-dev libjq-dev libprotobuf-dev
                         libv8-dev protobuf-compiler openjdk-8-jdk libssh-dev libssl-dev
                         libgdal-dev libudunits2-dev libopenmpi-dev

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

##### Windows Installation

To install the full ODDRIN package on Windows, first install Ubuntu from the Microsoft Store. This allows Linux command syntax to be run on a Windows machine. Follow all instructions as per the *'Linux and Mac Installation'* section above *EXCEPT* the instruction to set the LD_LIBRARY_PATH. Ubuntu does not permit the user to set LD_LIBRARY_PATH in the '/etc/environment' file, so run the following in Ubuntu instead:

```bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-8-openjdk-amd64/lib/amd64/:/usr/lib/jvm/java-8-openjdk-amd64/lib/amd64/
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-8-openjdk-amd64/lib/amd64/:/usr/lib/jvm/java-8-openjdk-amd64/lib/amd64/' >> ~/.bashrc

```

Follow all the remaining steps after the setting of the LD_LIBRARY_PATH instruction as detailed in the *'Linux and Mac Installation'* section above, i.e. check that the folder actually exists, restart your computer, etc. For the sudo R CMD, please ensure that R is installed for use on Ubuntu, otherwise you will get a command not found error.

##### Data to be Downloaded Manually
<!--
Furthermore, you will also be required to download the following datasets:

  - Download the [sub-national GDP dataset](https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0) - Taken from the article by [Kummu, et al, 2020](https://www.nature.com/articles/sdata20184). Once downloaded, place the file called `GDP_per_capita_PPP_1990_2015_v2.nc` in the folder `Demography_Data/SocioEconomic/KUMMU/` which you must create inside the ODDRIN folder such that, from the console, `file.exists(paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc"))==TRUE`.
  - Then download the [high resolution population count dataset](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11/data-download) - Taken from the article produced by [Center for International Earth Science Information Network - CIESIN](https://doi.org/10.7927/H4F47M65). Please download **all years via the *Single Year* option, in ASCII format, at a resolution of 30 arc-seconds** as shown in the screenshot below. From what I remember, you may need to download each year at a time, otherwise the CIESIN servers hang. In all cases, you need the 2015 file as this is used to normalise the ODDRIN population values. The ASCII format is to ensure quick computation, and also to ensure that no problems are found when crossing the longitude=0 plane. Once downloaded, place the files in folders separately split by year, such that folders such as `paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2015/")` or `paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2005/")` exist. The files should then exist inside the correct folder under, for example, the name `gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc` (for the 2015 year). Please check that there should be 8 files in each folder finishing with the same extension `.asc` but as `...30_sec_2.asc`, `...30_sec_3.asc`, `...30_sec_4.asc` etc.
  - Download the [EM-DAT database](https://public.emdat.be/) using the query tool and selecting `Natural -> Geophysical -> Earthquake` as the Disaster Classification and filtering from 2008 - 2021. Name the file `emdat.xlsx` and save in the folder `Displacement_Data` such that the folder `paste0(dir,"Displacement_Data/emdat.xlsx")` exists.

  ### Required Datasets -->

In addition to installing the necessary packages, you are also required to **manually download** several datasets due to licensing and access restrictions. Please follow the instructions below carefully:
<!--
- **Sub-national GDP dataset**  
  Download the dataset from [Dryad Repository](https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0), as described in [Kummu et al., 2020](https://www.nature.com/articles/sdata20184). After downloading, place the file named `GDP_per_capita_PPP_1990_2015_v2.nc` inside the folder `Demography_Data/SocioEconomic/KUMMU/` (you will need to create this folder). You should be able to confirm the file has been correctly placed using the R check:  
  `file.exists(paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc")) == TRUE`.

- **High-resolution population count dataset**  
  Download the [GPWv4 population count dataset](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11/data-download) from [CIESIN](https://doi.org/10.7927/H4F47M65). Use the *Single Year* option and choose the **ASCII format at 30 arc-second resolution**. You will need to download files **year by year** (e.g., 2000, 2005, 2010, 2015), ensuring that each year has its own folder (e.g., `Demography_Data/Population/gpw-v4-population-count-2015/`). For the model, only the 2015 dataset is required, particularly the file named:  
  `gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc`  
  Make sure this file exists at:  
  `file.exists(paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2015/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc")) == TRUE`.  
  Each folder should contain 8 `.asc` files (`...30_sec_1.asc` to `...30_sec_8.asc`).
-->
- **Global Data Lab Vulnerability Data (SHDI/SGDI)**  
  You will need to download two sets of files:
  1. The **GDL shapefiles** from [this link](https://globaldatalab.org/asset/378/GDL%20Shapefiles%20V7.zip) (requires free account). Extract and place all files in:  
     `Demography_Data/SocioEconomic/GlobalDataLab/GDL Shapefiles V6/`. Ensure that the file `shdi2022_World_large.shp` is in that folder.
  2. The **CSV data** containing SHDI/SGDI values from [this link](https://globaldatalab.org/asset/348/SHDI-SGDI-Total%207.0.csv) (also requires free account). Name the file exactly:  
     `SHDI-SGDI-Total 7.0.csv` and place it in:  
     `Demography_Data/SocioEconomic/GlobalDataLab/`.

- **VS30 dataset (soil shear wave velocity)**  
  Download the dataset from the [USGS VS30 page](https://www.usgs.gov/programs/earthquake-hazards/science/vs30-models-and-data). Place the extracted `global_vs30.tif` and any auxiliary files in:  
  `Hazard_Data/global_vs30_tif/`.

- **Global Earthquake Hazard Frequency Data (PGA)**  
  Download the most recent version of the PGA hazard data from the [Global Earthquake Model's GSHM page](https://www.globalquakemodel.org/product/global-seismic-hazard-map). Extract the file named:  
  `v2023_1_pga_475_rock_3min.tif`  
  and place it in:  
  `Hazard_Data/GEM-GSHM_PGA-475y-rock_v2023/`.

- **High-resolution population count dataset**  
  Note that the CIESIN data does not currently seem to be available online. This dataset is not necessary as, if not downloaded, coords2country() is instead used to label the country of each grid cell. If the data does become available again, the download instructions are as follows:

  Download the [GPWv4 population count dataset](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11/data-download) from [CIESIN](https://doi.org/10.7927/H4F47M65). Use the *Single Year* option and choose the **ASCII format at 30 arc-second resolution**. You will need to download files **year by year** (e.g., 2000, 2005, 2010, 2015), ensuring that each year has its own folder (e.g., `Demography_Data/Population/gpw-v4-population-count-2015/`). For the model, only the 2015 dataset is required, particularly the file named:  
  `gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc`  
  Make sure this file exists at:  
  `file.exists(paste0(dir,"Demography_Data/Population/gpw-v4-population-count-2015/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc")) == TRUE`.  
  Each folder should contain 8 `.asc` files (`...30_sec_1.asc` to `...30_sec_8.asc`).

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
The files InstallationChecks.R, Main.R, and Autoquake.R provide usage examples.

<!--
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
-->
