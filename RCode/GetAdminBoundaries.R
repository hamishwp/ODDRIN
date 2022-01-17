# Make sure the necessary packages are installed
list.of.packages <- c("jsonlite","geojsonsf")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Function that extracts the boundary data of a country, on many different admin levels
GetNatBnd<-function(iso="ALL",admLev=1,sf_out=T){
  
  # Find the metadata
  url<-paste0("https://www.geoboundaries.org/api/current/gbOpen/",iso,"/ADM",admLev,"/")
  # Extract metadata (note that this includes the boundary data source location)
  infos<-jsonlite::fromJSON(url)
  # Extract only the shapefile data, preferring sf format as output
  if(sf_out) return(geojsonsf::geojson_sf(infos$gjDownloadURL))
  return(jsonlite::fromJSON(infos$gjDownloadURL))
  
}

