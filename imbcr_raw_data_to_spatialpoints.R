#
# build a consistent sample grid from raw IMBCR point count data and validate in a GIS
#

require(raster)
require(rgdal)

recursiveFindFile <- function(name=NULL,root=Sys.getenv("HOME")){
  return(list.files(root,pattern=name,recursive=T,full.names=T))
}

#
# Parse ALL BCR raw data tables into a single table
#
t <- recursiveFindFile(name="BCR1.*._raw_data.csv")
  t <- lapply(t,read.csv)

for(i in 1:length(t)){
  names(t[[i]]) <- tolower(names(t[[i]]))
}

t <- do.call(what=rbind,t)

#
# iterate over each UTM zone in the table, creating SpatialPoints
# projected to a focal UTM.  Then merge all of the zones together into
# a single shapfile with an arbitrary CRS.
#

zones <- list()
for (zone in unique(na.omit(t$ptvisitzone))){
  zones[[length(zones) + 1]] <- na.omit(t[t$ptvisitzone == zone,])
  zones[[length(zones)]] <- SpatialPointsDataFrame(
                              coords = data.frame(x = zones[[length(zones)]]$ptvisiteasting,
                              y = zones[[length(zones)]]$ptvisitnorthing),
                              data = zones[[length(zones)]],
                              proj4string = CRS(projection(paste("+init=epsg:269", zone, sep = "")))
                            )
}

zones <- lapply(zones,FUN=spTransform,CRS(projection(zones[[1]])))
  zones <- do.call(raster::merge,zones)
    zones$FID <- 1:nrow(zones)

#
# make a 1 kilometer grid consistent with the extent of our SpatialPoints
# data -- this will be used later for re-gridding LANDFIRE variables to
# project occupancy / abundance estimates.
#
grid <- raster(resolution=1000,crs=CRS(projection(zones)),ext=extent(zones))
  grid <- rasterize(zones,grid,field="FID")
  #  grid[is.na(grid)] <- 0
  #    grid <- as(grid,'SpatialPolygonsDataFrame')

writeOGR(grid,".","grid",driver="ESRI Shapefile")
writeOGR(zones,".","imbcr_sample_pts",driver="ESRI Shapefile",overwrite=T)
