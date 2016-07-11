require(raster)
require(rgdal)
require(rgeos)

calcLekCentroids <- function(x=NULL){
  leks <- as.vector(unique(x$LekGridID))

  f <- function(y) gCentroid(x[x$LekGridID == y,])
    centroids <- lapply(leks,FUN=f)
      centroids <- do.call(rbind,centroids)

  return(SpatialPointsDataFrame(centroids,
    data=data.frame(LekGridID=leks))
}

grid <- readOGR("Vector","WAFWA_LekGrid_NM_TX")

# attribute lek points with grid ID's
lek_pts <- spTransform(readOGR("Vector", "All_LekPoints_NotTemporal", verbose=F), CRS(projection(grid)))
  lek_pts$LekGridID <- sp::over(lek_pts, grid[,"LekGridID"])$LekGridID
    lek_pts <- lek_pts[!is.na(lek_pts$LekGridID),]

# calculate centroids
centroids <- calcLekCentroids(lek_pts)
