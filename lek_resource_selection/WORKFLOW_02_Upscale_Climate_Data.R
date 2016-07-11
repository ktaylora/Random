#
# WORKFLOW 2 : Upsample climate data using a high-resolution DEM
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

if(file.exists(argv)){
  path <- argv
} else {
  path <- "."
}

# read-in our source raster data from a directory path provided by the user
elevation  <- raster(paste(path,"elevation.tif",sep="/"))
cliRasters <- ifelse(file.exists(path), path, file.choose())
  cliRasters <- parLapply(cl,as.list(c("mat_tenths","map","ffp")),fun=fetchClimateData, dest=cliRasters)

# crop our rasters to the extent of the study area
studyAreaExtent <- readOGR(path,"Study_area",verbose=F)
  studyAreaExtent <- spTransform(studyAreaExtent,CRS(projection(cliRasters[[1]]))) # re-project our shapefile to the CRS of our elevation DEM

cat(" -- cropping and trimming our climate rasters to the extent of our study region\n")
cliRasters <- parLapply(cl,cliRasters,fun=crop,studyAreaExtent)

if(projection(elevation) != projection(cliRasters[[1]])){
  cat(" -- reprojecting climate rasters to the CRS of our elevation DEM\n")
    cliRasters <- parLapply(cl,cliRasters,fun=projectRaster,crs=CRS(projection(elevation)))
}

cat(" -- performing a bilinear interpolation on source climate data so they are spatially consistent with our DEM\n")
cliRasters <- parLapply(cl,cliRasters,fun=raster::resample,y=elevation,method='bilinear')
cat(" -- generating random point samples for GLM training data across elevation and climate variable surfaces\n")
cliRasterSamplePts <- parLapply(cl,cliRasters,fun=raster::sampleRandom,size=ncell(cliRasters[[1]])*0.00010,sp=TRUE,ext=extent(spTransform(studyAreaExtent,CRS(projection(cliRasters[[1]]))))*0.95) # grab a %0.01.5 random sample of points from our grid, returning as SpatialPoints*

# extract point values of elevation and interpolated climate data from our climate rasters and
# overfit a GLM of climate ~ f(elevation) for the current landscape (one for each climate variable)
# note: this operation doesn't fit well into the built-in 'R' apply methods... doing it the ugly way with for() loops

elevSamplePts <- list();
 trainingData <- list();
       models <- list();

cat(" -- fitting regression models:\n")

for(i in 1:length(cliRasterSamplePts)){
  elevSamplePts[[length(elevSamplePts)+1]] <- raster::extract(y=cliRasterSamplePts[[i]],x=elevation,df=T)
  trainingData[[length(trainingData)+1]] <- cbind(cliRasterSamplePts[[i]]@data,elevSamplePts[[i]])
    trainingData[[i]] <- trainingData[[i]][,names(trainingData[[i]])!="ID"] # strip out a looming "ID" field from our explanatory data
  # fit a simple linear regression and report the fit to our user
  models[[length(models)+1]] <- suppressWarnings(rlm(formula(paste(names(trainingData[[i]]),collapse="~")),
                                    scale.est="Huber", psi=psi.hampel, init="lts", na.rm=T,
                                    data=trainingData[[i]])) # model formula : climate variable ~ f(elevation)
  cat(paste("      -- formula: ",paste(names(trainingData[[i]]),collapse="~"),"\n",sep=""))
}

# project our overfit models across the extent of the current landscape
cat(" -- projecting models across current landscape DEM:\n")
fitSurfaces <- parLapply(cl,fun=raster::predict,X=models,object=elevation)
# save our output rasters
cat(" -- writing fit raster surfaces to disk:\n")
mapply(fitSurfaces,FUN=writeRaster,filename=as.list(paste(names(stack(cliRasters)),".upsampled.tif",sep="")),overwrite=T)
