require(raster)
require(rgdal)
require(utils)

# Set our Project Directory
setwd("Downloads")
# Decompress our source Lek point data and extent buffer zipfiles
unlink("/tmp/lek_pts",recursive=T,force=T); utils::unzip("All_LekPoints_NotTemporal.zip",exdir="/tmp/lek_pts")
unlink("/tmp/lek_project_boundaries",recursive=T,force=T); utils::unzip("data.zip",exdir="/tmp/lek_project_boundaries")
# read and process the source shapefiles, dropping lek points that are outside of the extent of the study region
lek_pts <- readOGR("/tmp/lek_pts/","All_LekPoints_NotTemporal",verbose=F)
  lek_pts <- lek_pts[lek_pts$Year > 2002,]
b <- spTransform(readOGR("/tmp/lek_project_boundaries","LEPC_EOR10_11082013_AlbersEAC",verbose=F),CRS(projection(lek_pts)))
  lek_pts <- lek_pts[!is.na(sp::over(lek_pts,b)[,1]),]

# Calculate a standard UNIX time variable to inform our hierarchical clustering
lek_pts$uTime <- as.numeric(as.POSIXct(lek_pts$observatio))
# Calculate a distance metric that we can believe in using X/Y and our time metric (in seconds)
lek_dist_metric <- dist(lek_pts@data[,c('POINT_X','POINT_Y','uTime')])
# Build our hierarchical cluster
out<-hclust(lek_dist_metric,method='complete')
#
lek_pts$cMemb <- cutree(out, k=500)
