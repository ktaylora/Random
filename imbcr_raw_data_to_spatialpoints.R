#
# build a consistent sample grid from raw IMBCR point count data and validate in a GIS
#

require(raster)
require(rgdal)

setwd("~/Downloads")

zones <- list()
    t <- read.csv("BCR17_NE-BCR18_CO-BCR18_raw_bird_counts.csv")

for(zone in unique(na.omit(t$PtVisitZone))){
  zones[[length(zones)+1]] <- na.omit(t[t$PtVisitZone == zone,])
  zones[[length(zones)]] <- SpatialPointsDataFrame(coords=data.frame(x=zones[[length(zones)]]$PtVisitEasting,
                                                   y=zones[[length(zones)]]$PtVisitNorthing),
                                                   data=zones[[length(zones)]],
                                                   proj4string=CRS(projection(paste("+init=epsg:269",zone,sep=""))))
}

zones <- lapply(zones,FUN=spTransform,CRS(projection(zones[[1]])))
  zones <- do.call(raster::merge,zones)
    zones$FID <- 1:nrow(zones)

 grid <- raster(resolution=1000,crs=CRS(projection(zones)),ext=extent(zones))
   grid <- rasterize(zones,grid,field="FID")
    #  grid[is.na(grid)] <- 0
    #    grid <- as(grid,'SpatialPolygonsDataFrame')

writeOGR(grid,".","grid",driver="ESRI Shapefile")
writeOGR(zones,".","imbcr_sample_pts",driver="ESRI Shapefile")
