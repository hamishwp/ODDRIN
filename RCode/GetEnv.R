# Where is the main folder with all the code and data
dir<-directory<- "/home/manderso/Documents/GitHub/ODDRIN/" #"./" #"/home/hamishwp/Documents/BEAST/Coding/Oxford/ODDRIN/"
# Set the working directory from your environment variables
setwd(directory)
# Directory of the Data for Good data, e.g. Disaster Mapping, 4G connectivity, etc
FBdirectory<-'/home/patten/Documents/IDMC/Facebook_Data/'
# Do you want only the reduced packages or all? Choose via packred
packred<-F
# Do you want to load the Rmpi package? This is used for MPI (parallelising across nodes) - only required for the parallel implementation of the ABC-SMC algorithm
loadRmpi <- T

