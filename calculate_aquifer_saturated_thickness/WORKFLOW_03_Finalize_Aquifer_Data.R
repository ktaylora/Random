#
# Clean up the output from the ArcGIS Workflow and write averages for training/evaluation data to disk
#

require(raster)
require(rgdal)
require(parallel)

cl <- makeCluster(4)
argv <- commandArgs(trailingOnly=T)
  if(length(argv)<1) argv <- "." # by default, let's assume that the user is working in the CD. 

if(!file.exists(argv[1])) stop("usage: R < script [directory with rasters and maskfile]")

# mask our rasters against the boundaries of the HP aquifer
s <- readOGR(argv[1],"hpbedrock",verbose=F)
r <- list.files(argv[1],pattern="tif$")
  r <- lapply(as.list(r),FUN=raster)

s <- spTransform(s,CRS(projection(r[[1]])))
  r <- parLapply(cl,r,mask,mask=s)

satThick_10_14 <- stackApply(stack(r[grepl(unlist(lapply(r,names)),pattern="2010|2011|2012|2013|2014")]),fun=mean,indices=1)

# write our output to disk
writeRaster(satThick_10_14,"satThick_10_14.tif",overwrite=T)
