#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%% AutoQuake - automatically generate object from earthquake data %%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%% Developed by Hamish Patten %%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%% Collaboration between the University of Oxford %%%%%%%%%%%%#
#%%%%%%%%%%%% and the Internal Displacement Monitoring Centre %%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%% Funded by the EPSRC - Impact Acceleration Account %%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%% Started March 2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Set the working directory
dir<-directory<-"/home/patten/Documents/Coding/Oxford/AutoQuake/" ; setwd(dir)
# Download and install the necessary packages, and load source files & data:
source('RCode/GetODDPackages.R')
# Check that background data exists - Population & GDP (GetData.R)
CheckBGData(input$datadir)

#%%%%%%%%%%%%% User defined input - bare minimum required %%%%%%%%%%%%%%%#
input<-list(
  sdate=as.Date("2019-12-13"), # "YYYY-MM-DD"
  fdate=as.Date("2019-12-17"), # "YYYY-MM-DD"
  iso3="PHL", # Country code in ISO-3C form
  datadir=dir, # Location of the main folder to access the data 
  plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)
# Or extract the data purely based on the USGS id number
input<-list(USGSid="usp000huvq",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

input<-list(USGSid="at00qxtxcn",
            datadir=dir, # Location of the main folder to access the data 
            plotdir="Plots/" # Location for plots as paste0(datadir,plotdir)
)

#%%%%%%%%%%%%% Variables and functions from IIDIPUS files %%%%%%%%%%%%%%%#
input%<>%append(list(Model=Model, # Model parameters - equations, parameters, ... (Model.R)
                     Omega=Omega, # Parameterisation of the model (GetInitialValues.R)
                     Method=AlgoParams)) # Number of CPUs, particles in SMC, etc. (Method.R)

#%%%%%%% Function to extract data, predict displacement & plot %%%%%%%%%%#
AutoQuake<-function(input,extras=T){
  # Find the earthquake data via APIs (GetDisaster.R & GetUSGS.R)
  EQ<-GetEarthquake(input)
  # Create ODD object from EQ, Population, GDP, ... datasets (ODDobj.R)
  ODDy<-new("ODD",lhazSDF=EQ, dir=input$datadir, Model=input$Model) ; rm(EQ) # Model defined in Model.R
  
  if(!extras) return(ODDy)
  
  # Predict displacement (ODDobj.R)
  ODDy%<>%DispX(Omega = input$Omega, 
                center = input$Model$center, 
                LL = F, # Do not return the log-likelihood, but the displacement prediction
                Method = input$AlgoParams) 
  # Make plots and store them in a specific folder (ODDobj.R)
  pout<-MakeODDPlots(ODDy, input)
  
  return(ODDy)
}

# RUN IT! (Cross your fingers that the EQ exist in USGS)
ODDy<-AutoQuake(input)


