bbox=c(166.32,-18.5,169.2,-14)

bbox<-c(ODDy@bbox)

mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "terrain",zoom=8)

p<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")

# p+geom_contour_filled(data = as.data.frame(ODDy),
                         # mapping = aes(Longitude,Latitude,z=log10(ODDy@data[["Population"]])))

q<- p+ geom_raster(data=as.data.frame(ODDy),aes(Longitude,Latitude,fill=Population),
                   alpha=0.5,interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
  scale_fill_gradient2(low = "blue",mid="blue",high = "red",trans = "log",
                       # breaks=c(0,1,10,100),
                       na.value = "transparent")+#,limits=c(0,100)) +
  # labs(fill = "Displaced Pop. / Total Pop.")+xlab("Longitude") + ylab("Latitude");p
  labs(fill = "Population Density")+xlab("Longitude") + ylab("Latitude")

tmpy<-as(polys, "SpatialPolygonsDataFrame")
tmpy@data$id<-1:nrow(tmpy@data)
tmpy%<>%fortify(region='id')

outer<-q+geom_sf(data=st_as_sf(as(polys, "SpatialPolygonsDataFrame")),alpha=0,size=1.,colour="orange",inherit.aes = F)



ggsave(outer,filename = "~/Documents/VUT.png")
