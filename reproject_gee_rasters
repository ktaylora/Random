#!/usr/bin/Rscript

argv <- commandArgs(trailingOnly=T)

if(length(argv)<1) stop("requires one argument: ",
  "[full path to gee zipfile with raster to reproject] or a directory of ",
  "tifs specified with the -d argument")

require(snow)
require(raster)

beginCluster(n=2)
rasterOptions(tmpdir="/home/ktaylora/r_raster_tmp")

#
# local function declarations
#

reproject_raster <- function(r=NULL,template=NULL){
  if(!inherits(r,"Raster")){
    r <- raster::raster(r)
  }
  # fixed raster template object we are projecting "to"
  if(is.null(template)){
    template <- raster("/global_workspace/ring_necked_pheasant_imbcr_models/
                        raster/2016_crp_107x107.tif")
  }
  cat(" -- reprojecting:",names(r),"\n")
  reprojected <- raster::projectRaster(r, crs=CRS(projection(template)),progress='text')
    final <- raster::projectRaster(reprojected,to=template, progress='text')
  return(final)
}


#
# MAIN
#

using_directory <- which(grepl(tolower(argv),pattern="-d"))

if(sum(using_directory)>0){
  rasters <- list.files(argv[using_directory+1], pattern="tif$", fullnames=T)
  if(length(rasters)<1){
    stop("-d directory argument didn't contain any .tif images")
  }
# if there wasn't a directory specification we will treat this as one image
} else {
  if(!file.exists(argv)){
    stop("file",argv,"doesn't exist")
  } else {
    if(grepl(argv,pattern="zip$")){
      if(dir.exists("unpack")){
        unlink("unpack")
      } 
      dir.create("unpack")
      utils::unzip(argv, exdir="unpack", overwrite=T)
      rasters <- list.files("unpack", pattern="tif$", full.names=T)
    } else { # assume it's a raster : AFNP 
        rasters <- argv
    }
  }
}

for(r in rasters){
  target <- <- gsub(r, pattern="tif$", replacement="reprojected.tif")
  if(file.exists(target)){
    cat("target file",target,"already exists... skipping...\n")
  } else {
    r <- reproject_raster(r)  
    writeRaster(r, target, overwrite=T, progress='text')
    rm(r);
    raster::removeTmpFiles();
    gc()
  } 
}

endCluster()
