require(raster)
require(rgdal)
require(rgeos)

calcLekCentroids <- function(x=NULL){
  leks <- as.vector(unique(x$LekGridID))

  f <- function(y) gCentroid(x[x$LekGridID == y,])
    centroids <- lapply(leks,FUN=f)
      centroids <- do.call(rbind,centroids)

  return(SpatialPointsDataFrame(centroids,
    data=data.frame(LekGridID=leks)))
}


grid <- readOGR("Vector","WAFWA_LekGrid_NM_TX")

# attribute lek points with grid ID's
lek_pts <- spTransform(readOGR("Vector", "All_LekPoints_NotTemporal", verbose=F), CRS(projection(grid)))
  lek_pts$LekGridID <- sp::over(lek_pts, grid[,"LekGridID"])$LekGridID
    lek_pts <- lek_pts[!is.na(lek_pts$LekGridID),]
       lek_pts$active <- as.numeric(lek_pts$Count_ > 0 & lek_pts$Year >= 2005) # Active within the last 7 years (last obs year = 2012)

# Calculate the sum of our lek counts for each grid unit, aggregating by ID
t <- lek_pts@data[,c("LekGridID","active","Year")]
  t <- aggregate(active~LekGridID, FUN=sum, data=t, na.rm=T)
    t$active <- as.numeric(t$active>0)
      t <- merge(grid,t, by="LekGridID", all=F)
          t <- t[,c("LekGridID","active.x")]
            names(t) <- c("LekGridID","active")
              grid <- t

writeOGR(grid,".","wafwa_grid_attributed",driver="ESRI Shapefile", overwrite=T)
