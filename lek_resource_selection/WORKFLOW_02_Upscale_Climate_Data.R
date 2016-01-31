#
# WORKFLOW 2 : "Upscale" climate data using a high-resolution DEM
#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2016]
#

# load our default libraries
require(raster)
require(rgdal)
require(parallel)

# by default, let's assume the user passed us a directory path at runtime
argv <- commandArgs(trailingOnly=T)
path <- if(file.exists(argv),argv,".")

# read-in our source raster data from a directory path provided by the user
elevation  <- raster(paste(path,"elevation.tif",sep="//"))
cliRasters <- ifelse(file.exists(argv[1]), argv[1], file.choose())
  cliRasters <- list.files(cliRasters, pattern=paste(c("temp.img$","precip.img$","yada.img$"),collapse="|")) # filter files by expected names
    cliRasters <- lapply(as.list(cliRasters),fun=raster)

# crop our rasters to the extent of the study area
studyAreaExtent <- readOGR(path,"studyAreaExtent",verbose=F)
  studyAreaExtent <- spTransform(studyAreaExtent,CRS(projection(elevation))) # re-project our shapefile to the CRS of our elevation DEM

cliRasters <- lapply(cliRasters,FUN=crop,studyAreaExtent)

# perform a bi-linear interpolation of our climate rasters so they are a consistent resolution with our elevation DEM

# generate a large sample of random points across the extent of our study area

# extract point values of elevation and interpolated climate data from our climate rasters

# overfit GLMs of climate ~ f(elevation) for the current landscape (one for each climate variable)

# project our overfit models across the extent of the current landscape

# save our output rasters
