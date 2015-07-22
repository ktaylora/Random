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
  t_routes_bcr1819 <<- read.csv("/home/ktaylora/PLJV/species_data/bbs_data/bbs_routes_bcr1819.csv")
  points  <- habitatWorkbench::sampleAroundVertices(s=route,maxDist=2200,n=250) # generate a number of sampling points around each route to derive our site-level landscape metrics
  buffers <- landscapeAnalysis::subsampleSurface(x=landcover,pts=points, width=750) # sample buffers from our source landcover dataset with the points derived from the current route
  # parse our buffers out into focal landcover types
  buffers_grassland     <- landscapeAnalysis::lReclass(buffers,inValues=c(31,37,39,71,75))
  buffers_agriculture   <- landscapeAnalysis::lReclass(buffers,inValues=c(38,201,202,203,205,206,207,208,209,210,211,212))
  buffers_shrubland     <- landscapeAnalysis::lReclass(buffers,inValues=c(83,85,87,81,82))
  # calculate the requisite landscape metrics for each cover type
  metrics_grassland <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_grassland)[[1]] # need: total area and mean patch area
  grassland_total_area <- unlist(lapply(metrics_grassland, function(x){ x$total.area }))
    grassland_total_area <- ifelse(is.null(grassland_total_area),0,grassland_total_area)
      grassland_total_area[is.na(grassland_total_area)] <- 0
        grassland_total_area <- mean(grassland_total_area)
  grassland_mean_patch_area <- unlist(lapply(metrics_grassland, function(x){ x$mean.patch.area }))
    grassland_mean_patch_area[is.na(grassland_mean_patch_area)] <- 0
      grassland_mean_patch_area <- ifelse(is.null(grassland_mean_patch_area),0,grassland_mean_patch_area)
        grassland_mean_patch_area[is.na(grassland_mean_patch_area)] <- 0
          grassland_mean_patch_area <- mean(grassland_mean_patch_area)
  metrics_agriculture   <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_agriculture) # need: total area
    agriculture_total_area <- unlist(lapply(metrics_agriculture, function(x){ x$total.area }))
      agriculture_total_area <- ifelse(is.null(agriculture_total_area),0,agriculture_total_area)
        agriculture_total_area[is.na(agriculture_total_area)] <- 0
          agriculture_total_area <- mean(agriculture_total_area)
  metrics_shrubland     <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_shrubland) # need: total area
    shrubland_total_area <- unlist(lapply(metrics_shrubland, function(x){ x$total.area }))
      shrubland_total_area <- ifelse(is.null(shrubland_total_area),0,shrubland_total_area)
        shrubland_total_area[is.na(shrubland_total_area)] <- 0
          shrubland_total_area <- mean(shrubland_total_area)
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
require(parallel); cl <- makeCluster(getOption("cl.cores", 7))

# read-in our local (study area) routes
t_routes_bcr1819 <<- read.csv("/home/ktaylora/PLJV/species_data/bbs_data/bbs_routes_bcr1819.csv")

# read-in the national bbs routes data and parse accordingly
data(bbsRoutes); s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_routes_bcr1819$RTENO,]

# ensure we have adequate land cover data for our bbs routes
landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")
  landcover <- extract(landcover,as(spTransform(s_bbsRoutes,CRS(projection(landcover))),'SpatialPointsDataFrame'),sp=T)
s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% unique(landcover[!is.na(landcover$Final_LC_8bit),]$RTENO),] # make sure that our route data has landcover data available

# split our routes so that they can be processed in parallel on a multi-core machine
s_bbsRoutes <- split(s_bbsRoutes,f=1:nrow(s_bbsRoutes))
cat(" -- sampling and processing BBS routes\n");
out  <- parLapply(cl=cl,fun=processFocalRoute,X=s_bbsRoutes)
  out <- do.call(rbind,out) # bind our list into a single data.frame
