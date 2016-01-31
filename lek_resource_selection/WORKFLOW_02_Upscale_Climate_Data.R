# load our default libraries

require(raster)
require(rgdal)

argv <- commandArgs(trailingOnly=T)

# read-in our source climate raster data from the path provided by the user
rasters <- ifelse(file.exists(argv[1]),argv[1],file.choose())
  rasters <- list.files(rasters,pattern=paste(c("temp,precip,")))
