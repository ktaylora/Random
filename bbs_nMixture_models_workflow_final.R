#
# Workflow for Implementing N-Mixture Model with Breeding Bird Survey Data
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

require(landscapeAnalysis)
require(habitatWorkbench)
require(rgdal)
require(raster)

landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")

t_routes_bcr1819 <- read.csv("/home/ktaylora/PLJV/species_data/bbs_data/bbs_routes_bcr1819.csv")  # read-in our local (study area) routes
  data(bbsRoutes)                                                                                 # read-in the national bbs routes data and parse accordingly
    s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_routes_bcr1819$RTENO,]

# split and process all routes for the study area one-at-a-time
s_bbsRoutes <- split(s_bbsRoutes,f=1:nrow(s_bbsRoutes))
masterTable <- data.frame(route=rep(NA,length(s_bbsRoutes)),grass_total_area=NA,grass_mean_patch_area=NA,ag_total_area=NA,shrub_total_area=NA)

require(parallel)
cl <- makeCluster(getOption("cl.cores", 6))

processFocal <- function(route){
  line <- which(t_routes_bcr1819$RTENO == route)
  points  <- habitatWorkbench::sampleAroundVertices(s=route,maxDist=2200,n=250) # generate a number of sampling points around each route to derive our site-level landscape metrics
  buffers <- landscapeAnalysis::subsampleSurface(x=landcover,pts=points, width=750) # sample buffers from our source landcover dataset with the points derived from the current route
  # parse our buffers out into focal landcover types
  buffers_grassland     <- lReclass(buffers,inValues=c(31,39,71,75))
  buffers_agriculture   <- lReclass(buffers,inValues=c(38,201,202,203,205,206,207,208,209,210,211,212))
  buffers_shrubland     <- lReclass(buffers,inValues=c(83,85,87,81,82))
  # calculate the requisite landscape metrics for each cover type
  cat(".")
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
  cat(".")
  metrics_agriculture   <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_agriculture) # need: total area
    agriculture_total_area <- unlist(lapply(metrics_agriculture, function(x){ x$total.area }))
      agriculture_total_area <- ifelse(is.null(agriculture_total_area),0,agriculture_total_area)
        agriculture_total_area[is.na(agriculture_total_area)] <- 0
          agriculture_total_area <- mean(agriculture_total_area)
  cat(".")
  metrics_shrubland     <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_shrubland) # need: total area
    shrubland_total_area <- unlist(lapply(metrics_shrubland, function(x){ x$total.area }))
      shrubland_total_area <- ifelse(is.null(shrubland_total_area),0,shrubland_total_area)
        shrubland_total_area[is.na(shrubland_total_area)] <- 0
          shrubland_total_area <- mean(shrubland_total_area)
  # write to output table
  masterTable[line,] <- c(route$RTENO,grassland_total_area,grassland_mean_patch_area,agriculture_total_area,shrubland_total_area)
  cat("+")
}

cat(" -- sampling and processing BBS routes: ");
parLapply(cl,fun=processFocal,X=s_bbsRoutes)
