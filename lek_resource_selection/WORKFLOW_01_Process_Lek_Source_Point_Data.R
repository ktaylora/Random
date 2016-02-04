require(raster)
require(rgdal)
require(utils)
require(FNN)
require(dbscan)

# Define a mode function, just in case we want to use it
Mode <- function(x,na.rm=T){
  if(na.rm) x <- na.omit(x)
  ux <- unique(x)
  ux[which.max(tabulate(match(x,ux)))]
}
# Solve for k using Moors' Kurtosis
mKurtosis <- function(x,na.rm=T){
  if(na.rm) x <- na.omit(x);
    z <- (x-mean(x))/sd(x)
      var(z^2)+1
}

#
# MAIN
#

# Set our Project Directory
setwd("~/Downloads")
cat(" -- decompressing our source lek point data and extent buffer zipfiles\n ");
unlink("/tmp/lek_pts",recursive=T,force=T); utils::unzip("All_LekPoints_NotTemporal.zip",exdir="/tmp/lek_pts")
unlink("/tmp/lek_project_boundaries",recursive=T,force=T); utils::unzip("data.zip",exdir="/tmp/lek_project_boundaries")
cat(" -- reading and processing the source shapefiles, dropping lek points that are outside of the extent of the study region\n")
lek_pts <- readOGR("/tmp/lek_pts/","All_LekPoints_NotTemporal",verbose=F)
  lek_pts <- lek_pts[lek_pts$Year > 2002,]
# christian's active vs not-active rules
keep <- unique(lek_pts[lek_pts$Year == 2011 & lek_pts$Count_ > 0,]$SFs_sfid) # we had at least one chicken observed last year?

# overlapping points?
lek_pts$overlaps <- sqrt(lek_pts$SFs_sfid*lek_pts$POINT_X*lek_pts$POINT_Y)
  lek_pts$overlaps <- duplicated(lek_pts$overlaps) | duplicated(lek_pts$overlaps, fromLast=TRUE)
b <- spTransform(readOGR("/tmp/lek_project_boundaries","LEPC_EOR10_11082013_AlbersEAC",verbose=F),CRS(projection(lek_pts)))
  lek_pts <- lek_pts[!is.na(sp::over(lek_pts,b)[,1]),]
cat(" -- calculating a standard UNIX time variable to inform our hierarchical clustering\n")
lek_pts$uTime <- as.numeric(as.POSIXct(lek_pts$observatio))
cat(" -- rasterizing to merge overlapping points at 30m resolution, then converting back to a points shapefile\n")
lek_pts_r <- raster(lek_pts,res=c(30,30),crs=CRS(projection(lek_pts)))
  lek_pts_r <- rasterize(lek_pts,field=c('SFs_sfid','POINT_X','POINT_Y','Year','uTime','Count_'),fun=median,y=lek_pts_r)
    lek_pts <- rasterToPoints(lek_pts_r,spatial=T)
      names(lek_pts) <- c('SFs_sfid','POINT_X','POINT_Y','year','uTime','count')

cat(" -- calculating cluster densities from spatial data (knn=1,5) and Moors' Kurtosis\n")
lek_pts$knn1 <- rowMeans(get.knn(lek_pts@data[,c('POINT_X','POINT_Y')],k=1)$nn.dist)
lek_pts$knn5 <- rowMeans(get.knn(lek_pts@data[,c('POINT_X','POINT_Y')],k=5)$nn.dist)

# Moors' Kurtosis
f <- get.knn(lek_pts@data[,c('POINT_X','POINT_Y')],k=30)$nn.dist
  lek_pts$mKurt30 <- apply(f,1,mKurtosis); rm(f)
cat(" -- calculating a distance metric that we can believe in using X/Y and our time metric (in seconds)\n")
lek_dist_metric <- dist(lek_pts@data[,c('POINT_X','POINT_Y','knn1','knn5','mKurt30')],method='manhattan')
cat(" -- building our hierarchical cluster and cut our cluster trees (using rule of thumb for k)\n")
lek_pts$cMemb_hk <- cutree(hclust(lek_dist_metric,method='centroid'), k=round(nrow(lek_pts)/2))
cat(" -- building our k-means cluster (using rule of thumb for k)\n")
lek_pts$cMemb_km <- kmeans(lek_pts@data[,c('POINT_X','POINT_Y','knn1','knn5','mKurt30')],centers=round(nrow(lek_pts@data)/2),nstart=50,iter=500,trace=F)$cluster

# evaluation
# for(i in 1:500){ plot(lek_pts[lek_pts$cMemb_hk == i,],pch=15, main=paste("cluster",i)); scalebar(type="bar",lwd=4); Sys.sleep(1.5); }
# hist(lek_pts$cMemb_hk,plot=F,breaks=250)$counts
# hist(lek_pts$cMemb_km,plot=F,breaks=250)$counts

cat(" -- writing our output to disk\n")
writeOGR(lek_pts,"/tmp/lek_pts/","lek_pts_processed",driver="ESRI Shapefile",overwrite=T)
