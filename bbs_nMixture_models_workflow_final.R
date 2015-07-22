#
# Workflow for Implementing N-Mixture Model with Breeding Bird Survey Data
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

#
# calculate site-level landscape metrics for the focal route [parallelized]
#
processFocalRoute <- function(route=NULL){
  require(landscapeAnalysis)
  require(habitatWorkbench)
  require(rgdal)
  require(raster)
  landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")
  t_routes_bcr1819 <- read.csv("/home/ktaylora/PLJV/species_data/bbs_data/bbs_routes_bcr1819.csv")
  points  <- habitatWorkbench::sampleAroundVertices(s=route,maxDist=2200,n=250) # generate a number of sampling points around each route to derive our site-level landscape metrics
  buffers <- landscapeAnalysis::subsampleSurface(x=landcover,pts=points, width=750) # sample buffers from our source landcover dataset with the points derived from the current route
  # parse our buffers out into focal landcover types
  buffers_grassland     <- landscapeAnalysis::lReclass(buffers,inValues=c(31,37,39,71,75))
  buffers_agriculture   <- landscapeAnalysis::lReclass(buffers,inValues=c(38,201,202,203,205,206,207,208,209,210,211,212))
  buffers_shrubland     <- landscapeAnalysis::lReclass(buffers,inValues=c(83,85,87,81,82))
  # calculate the requisite landscape metrics for each cover type
  metrics_grassland <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_grassland)[[1]] # need: total area and mean patch area
  grassland_total_area <- try(metricsListToVector(metrics_grassland,'total.area'))
    if(class(grassland_total_area)=="try-error"){
      print(str(metrics_grassland))
      grassland_total_area <- 0
    } else {
      grassland_total_area <- ifelse(is.null(grassland_total_area),0,grassland_total_area)
        grassland_total_area[is.na(grassland_total_area)] <- 0
          grassland_total_area <- mean(grassland_total_area)
    }
  grassland_mean_patch_area <- try(metricsListToVector(metrics_grassland,'mean.patch.area'))
    if(class(grassland_mean_patch_area)=="try-error"){
      print(str(metrics_grassland))
      grassland_mean_patch_area <- 0
    } else {
      grassland_mean_patch_area <- ifelse(is.null(grassland_mean_patch_area),0,grassland_mean_patch_area)
        grassland_mean_patch_area[is.na(grassland_mean_patch_area)] <- 0
          grassland_mean_patch_area <- mean(grassland_mean_patch_area)
    }
  metrics_agriculture   <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_agriculture) # need: total area
  agriculture_total_area <- try(metricsListToVector(metrics_agriculture,'total.area'))
    if(class(agriculture_total_area)=="try-error"){
      print(str(metrics_agriculture))
      agriculture_total_area <- 0
    } else {
      agriculture_total_area <- ifelse(is.null(agriculture_total_area),0,agriculture_total_area)
        agriculture_total_area[is.na(agriculture_total_area)] <- 0
          agriculture_total_area <- mean(agriculture_total_area)
    }
  metrics_shrubland     <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_shrubland) # need: total area
  shrubland_total_area <- try(metricsListToVector(metrics_shrubland,'total.area'))
    if(class(shrubland_total_area)=="try-error"){
      print(str(metrics_shrubland))
      shrubland_total_area <- 0
    } else {
      shrubland_total_area <- ifelse(is.null(shrubland_total_area),0,shrubland_total_area)
        shrubland_total_area[is.na(shrubland_total_area)] <- 0
          shrubland_total_area <- mean(shrubland_total_area)
    }
  # write to output table
  return(data.frame(route=route$RTENO,grass_total_area=grassland_total_area,grass_mean_patch_area=grassland_mean_patch_area,
    ag_total_area=agriculture_total_area,shrub_total_area=shrubland_total_area))
}

#
# MAIN
#

require(rgdal)
require(raster)
require(habitatWorkbench)
require(parallel); cl <- makeCluster(getOption("cl.cores", 7),outfile='outfile.log')

# read-in our local (study area) routes
t_routes_bcr1819 <<- read.csv("/home/ktaylora/PLJV/species_data/bbs_data/bbs_routes_bcr1819.csv")

# read-in the national bbs routes data and parse accordingly
data(bbsRoutes); s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_routes_bcr1819$RTENO,]

# ensure we have adequate land cover data for our bbs routes
landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")
  landcover <- extract(landcover,as(spTransform(s_bbsRoutes,CRS(projection(landcover))),'SpatialPointsDataFrame'),sp=T)
s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% unique(landcover[!is.na(landcover$Final_LC_8bit),]$RTENO),] # make sure that our route data has landcover data available
  s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% unique(landcover[landcover$Final_LC_8bit != 0,]$RTENO),]

# split our routes so that they can be processed in parallel on a multi-core machine and then calculate some site level statistics for landcover at each route
s_bbsRoutes <- split(s_bbsRoutes,f=1:nrow(s_bbsRoutes))
if(!file.exists("site_level_parameters.csv")){
  cat(" -- sampling and processing BBS routes\n");
    out  <- parLapply(cl=cl,fun=processFocalRoute,X=s_bbsRoutes)
      out <- do.call(rbind,out) # bind our list into a single data.frame
        write.csv(out,"site_level_parameters.csv",row.names=F)
} else {
  out <- read.csv("site_level_parameters.csv")
}

## calculate our isolation metrics (at 3.3 and 30 km2 scales)
if(!grepl(names(out),pattern="isolation")){
    landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")
  s_bbsRoutes <-lapply(s_bbsRoutes,spTransform,CRSobj=CRS(projection(landcover)))
    centroids <- parLapply(cl=cl,s_bbsRoutes,fun=getBbsRouteLocations,centroid=T) # calculate the centroid of each route
buffers_3.3km <- parLapply(cl=cl,X=centroids,fun=rgeos::gBuffer,width=3300/2)      # calculate our 3.3 km regions
  m <- mcmapply(buffers_3.3km, FUN=raster::crop, MoreArgs=list(x=landcover))
    buffers_3.3km <- mcmapply(FUN=raster::mask, x=m, mask=buffers_3.3km); rm(m);
      buffers_3.3km <- lReclass(buffers_3.3km,inValues=c(31,37,39,71,75))
        buffers_3.3km <- parLapply(cl=cl,buffers_3.3km,fun=landscapeAnalysis::rasterToPolygons)
          # check for NA values
          na_values<-as.vector(unlist(lapply(X=buffers_3.3km,FUN=is.na)))
          if(sum(na_values)>0){
            buffers_3.3km[na_values] <- buffers_3.3km[which(!na_values)[1]] # overwrite our NA values with something valid
            buffers_3.3km<-lapply(X=buffers_3.3km,FUN=getSpPPolygonsLabptSlots)

          } else {
            buffers_3.3km<-lapply(X=buffers_3.3km,FUN=getSpPPolygonsLabptSlots)
          }
          # do our NN assessment
          d <- function(x,na.rm=F){ o<-try(FNN::knn.dist(x,k=1)); if(class(o) != "try-error") { x <- o; } else { x[[i]] <- NA }; return(x)}
          buffers_3.3km <- lapply(buffers_3.3km,FUN=d);rm(d);
            buffers_3.3km <- lapply(buffers_3.3km, FUN=mean)
              buffers_3.3km[na_values] <- NA # restore our NA values
                buffers_3.3km
                  out$isolation_3.3km<-as.vector(unlist(buffers_3.3km))
}
# buffers_30km  <- parLapply(cl=cl,X=centroids,fun=rgeos::buffer,width=30000/2)     # calculate our 30 km regions
#   m <- parLapply(cl=cl, buffers_30km, fun=raster::crop, x=raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif"))
#     buffers_30km <- mcmapply(FUN=raster::mask, x=m, mask=buffers_30km); rm(m);
#       buffers_30km <- lReclass(buffers_30km,inValues=c(31,37,39,71,75))
