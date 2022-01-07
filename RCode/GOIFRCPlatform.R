library(httr)

ODDy<-readRDS("./IIDIPUS_Results/EQ20210814HTI_10919_5")
ODDy@data%<>%dplyr::select(Population,hazMean1,ISO3C,NumBDFB,FBPop,NumBDosm,ExpDisp,DamBDosm,DamBDFB,RWI)
colnames(ODDy@data)[1]<-"SEDACpop"

# Unique GLIDE ID
glide<- "EQ-2021-000018-IRN"
imainhaz<-grep(names(ODDy@data),pattern = "hazMean")[1]
hazard<-"EQ";haz<-"EQ"
# Define the ISO3C as the country which has the largest number of datapoints exposed to at least I0 intensity hazard
exposure<-ODDy@data[,imainhaz]>ODDy@I0 & !is.na(ODDy@data[,imainhaz])
iso3<-names(which.max(table(ODDy$ISO3C[exposure])))
# Centroid of the hazard intensity
loc<-ODDy@coords[which.max(ODDy@data[,imainhaz]),]
Longitude<-unname(loc[1])
Latitude<-unname(loc[2])
# Exposed population
PopExp<-sum(ODDy@data$FBPop,na.rm=T)
# Displaced population
Disp<-c(sum(ODDy@data$LowCIDisp,na.rm = T),
        sum(ODDy@data$ExpDisp,na.rm = T),
        sum(ODDy@data$HiCIDisp,na.rm = T))
Disp[1]<-269821
Disp[3]<-533393
modbbox<-c(min(ODDy@coords[exposure,1]),
           min(ODDy@coords[exposure,2]),
           max(ODDy@coords[exposure,1]),
           max(ODDy@coords[exposure,2]))

hazard_title<-paste0(haz,gsub("-","",ODDy@hazdates[1]),iso3,"_1")

colnames(ODDy@data)[2]<-"hazMean20210814_1"

# modbbox<-c(-74.7,-17.25,-73.2,18.7)
# tbuildings<-tryCatch(ExtractOSMbuild(modbbox),error=function(e) NULL)
BdExp<-sum(ODDy@data$NumBDFB,na.rm = T)

filepather<-"./IFRC_tmpODDRIN.tiff"
filepather<-paste0(getwd(),"/IFRC_tmpODDRIN.tiff")
# ODDy$ISO3C<-NULL
# class(ODDy)<-"SpatialPixelsDataFrame"
# ODDy%<>%brick()
# writeRaster(ODDy, filename=filepather, 
#             options="INTERLEAVE=BAND", overwrite=TRUE)

# tRast<-convODDy2raster(ODDy,filepather)

url<-"https://risk-module-api.togglecorp.com/api/v1/global-exposure-data/"
body<-list(hazard_title=hazard_title,
           hazard_type=hazard,
           glide_number=glide,
           source_type="oddrin",
           latitude=Latitude,
           longitude=Longitude,
           iso3=iso3,
           people_exposed=as.integer(PopExp),
           buildings_exposed=as.integer(BdExp),
           people_displaced=as.integer(Disp[2]),
           lower_displacement=Disp[1],
           higher_displacement=Disp[3],
           file=httr::upload_file(filepather),
           file_type="raster")

headers = c('Content-Type' = 'multipart/form-data; boundary=----WebKitFormBoundary7MA4YWxkTrZu0gW')
r<-httr::POST(url = url, httr::add_headers(.headers=headers), body = body,verbose())






# rebody <- jsonlite::toJSON(body, auto_unbox = TRUE, force = TRUE)
# headers = c('Content-Type' = 'application/json; charset=UTF-8')
# 
# r<-httr::POST(url = url, httr::add_headers(.headers=headers), body = rebody,encode="json", verbose())
# r<-httr::POST(url = url, httr::add_headers(headers), body = rebody,encode="json", verbose())
# r<-httr::POST(url = url, body = rebody,encode="json",multipart=F, verbose())
# 
# r<-httr::POST(url = url, httr::add_headers(.headers=headers), body = body,encode="json", verbose())
# httr::content(r)
# content(r, as="text") %>%
#   jsonlite::fromJSON(flatten=TRUE) %>%
#   glimpse()
# 
# rebofile<-jsonlite::toJSON(list(file=httr::upload_file(filepather)),
#                            auto_unbox = TRUE, force = TRUE)
# # file=httr::upload_file(filepather, type = "raster/tiff")
# 
# # fullbody<-list(data=body,files=list(files=httr::upload_file(filepather)))
# # rebo <- jsonlite::toJSON(fullbody, auto_unbox = TRUE, force = TRUE)
# 
# r<-httr::POST(url = url, httr::add_headers(.headers=headers), body = rebody,verbose())
# httr::content(r)
# 
# # rebody<-'{"GLIDE":"EQ-2021-000018-IRN"}'
# 
# # headers = c('Content-Type' = 'application/json; charset=UTF-8')
# headers = c('Content-Type' = 'multipart/form-data')
# r<-httr::POST(url = url, httr::add_headers(.headers=headers), body = rebody,verbose())
# r<-httr::POST(url = url, body = rebody,verbose())
# httr::content(r)
# 
# 
# 
# py_install("requests")
# library(reticulate)
# 
# reqy<-import("requests")
# reqy$Request(url = url,files = rebofile, data = rebody,
#              headers = httr::add_headers(.headers=headers))
# 
# respy<-reqy$(url = url,data = rebody,
#              headers = httr::add_headers(.headers=headers))
# 
# 




