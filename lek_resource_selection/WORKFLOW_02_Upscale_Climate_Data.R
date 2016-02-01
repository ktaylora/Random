#
# WORKFLOW 2 : "Upscale" climate data using a high-resolution DEM
#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2016]
#

#
# LOAD OUR DEFAULT LIBRARIES
#

require(raster)   # for manipulating georasters
require(rgdal)    # for working with ESRI Shapefiles
require(parallel) # for parallelizing our workflow across multiple CPU cores

#
# DEFINE LOCAL FUNCTIONS FOR THE WORKFLOW
#

#
# fetchClimateData()
# Local function that downloads a zipfile containing the Reihfeldt climate variable of
# interest from the USFS FTP server.  Does some redundant checking to make sure we aren't
# duplicating downloads if the file already exists in the CWD.  Make this parallelized by
# including the necessary 'R' libraries inside the function call.
#

fetchClimateData <- function(var=NULL,dest="."){
  require(raster)
  # is the focal variable fetchable from the Reihfeldt dataset?
  fetchable <- c("map","mat_tenths","ffp") # hardcoded values from USFS FTP server
  if(!sum(grepl(var,fetchable))){
    warning(paste(dest,"/",fetchable[var],".zip isn't a known Reihfeldt variable. Skipping.",sep=""))
    return(NULL)
  } else {
    var <- which(grepl(var,fetchable))
    if(file.exists(paste(dest,"/",fetchable[var],".zip",sep=""))){
      warning(paste(dest,"/",fetchable[var],".zip already exists. Skipping.",sep=""))
    } else {
      cat(" -- downloading.\n")
      download.file(url=paste("http://forest.moscowfsl.wsu.edu/climate/current/allNA/derivedGrids/",fetchable[var],".zip",sep=""),
                    destfile=paste(dest,"/",fetchable[var],".zip",sep=""))
    }
  }
  # Decompress and process the raster file, returning to user as a raster object
  if(!file.exists(paste(fetchable[var],"asc",sep="."))){
    cat(" -- decompressing.\n")
    utils::unzip(paste(dest,"/",fetchable[var],".zip",sep=""),exdir=dest)
    cat(" -- processing.\n")
    file.rename(paste(dest,"/",fetchable[var],".txt",sep=""),paste(dest,"/",fetchable[var],".asc",sep=""))
  } else {
    warning(paste(dest,"/",fetchable[var],".asc already exists. Did nothing",sep=""))
  }
  cat(" -- done.\n")
  var <- raster(paste(dest,"/",fetchable[var],".asc",sep=""),crs=CRS(projection("+init=epsg:4326"))) # read raster from file with a forced EPSG resolution
  return(var)
}

#
# MAIN WORKFLOW
#

# prepare our local 'cluster'
cl <- parallel::makeCluster(4) # use 4 CPU cores

# by default, let's assume the user passed us a directory path at runtime
argv <- commandArgs(trailingOnly=T)
path <- ifelse(file.exists(argv),argv,".")

# read-in our source raster data from a directory path provided by the user
elevation  <- raster(paste(path,"elevation.img",sep="/"))
cliRasters <- ifelse(file.exists(path), path, file.choose())
  cliRasters <- parLapply(cl,as.list(c("mat_tenths","map","ffp")),fun=fetchClimateData, dest=cliRasters)
if(projection(elevation) != projection(cliRasters[[1]])){
  cat(" -- reprojecting elevation DEM to the CRS of our climate data\n")
  elevation <- projectRaster(elevation,cliRasters[[1]]) # note: projectRaster will use our cl object by default -- no need to code a special wrapper function for it.
    writeRaster(elevation,paste(path,"elevation.img",sep="/"),overwrite=T)
}

# crop our rasters to the extent of the study area
studyAreaExtent <- readOGR(path,"Study_area",verbose=F)
  studyAreaExtent <- spTransform(studyAreaExtent,CRS(projection(cliRasters[[1]]))) # re-project our shapefile to the CRS of our elevation DEM

cat(" -- cropping our climate rasters to the extent of our study region\n")
cliRasters <- parLapply(cl,cliRasters,fun=crop,studyAreaExtent)

# perform a bilinear interpolation of our climate rasters so they are a consistent resolution with our elevation DEM
cliRasters <- parLapply(cl,cliRasters,fun=raster::resample,y=elevation,method='bilinear')

# generate a large sample of random points across the extent of our study area

# extract point values of elevation and interpolated climate data from our climate rasters

# overfit GLMs of climate ~ f(elevation) for the current landscape (one for each climate variable)

# project our overfit models across the extent of the current landscape

# save our output rasters
