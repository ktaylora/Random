#
# Clean up the output from the ArcGIS Workflow and write averages for training/evaluation data to disk
#

require(raster)
require(rgdal)
require(MASS)
require(parallel)

cl <- makeCluster(4)
argv <- commandArgs(trailingOnly=T)
  if(length(argv)<1) argv <- "." # by default, let's assume that the user is working in the CD.

if(!file.exists(argv[1])) stop("usage: R --no-save --vanilla --slave [directory with aquifer sat. thickness rasters] < script.R")

# mask our rasters against the boundaries of the HP aquifer
s <- readOGR(argv[1],"hpbedrock",verbose=F)
r <- list.files(argv[1],pattern="tif$")
  r <- r[!grepl(r,pattern="satThick")]
    r <- lapply(as.list(r),FUN=raster)

s <- spTransform(s,CRS(projection(r[[1]])))
  r <- parLapply(cl,r,mask,mask=s)

# calculate mean saturated thickness for a time-series consistent with NASS crop data
satThick_10_14 <- stackApply(stack(r[grepl(unlist(lapply(r,names)),pattern="2010|2011|2012|2013|2014")]),fun=mean,indices=1)

# build a table of trainable values and calculate a robust time-series regression for our aquifer depletion time-series
pts <- raster::sampleRandom(r[[1]],size=ncell(r[[1]])*0.35,sp=T,na.rm=T)

trainingData <- data.frame(sat.thickness=extract(r[[1]],pts),
                           yr=as.numeric(substr(names(r[[1]]),2,5)),
                           yr.sqrd=as.numeric(substr(names(r[[1]]),2,5))^2,
                           lat=pts@coords[,2],
                           lon=pts@coords[,1])

for(i in 2:length(r)){
  trainingData <- rbind(trainingData,
                        data.frame(sat.thickness=extract(r[[i]],pts),
                                   yr=as.numeric(substr(names(r[[i]]),2,5)),
                                   yr.sqrd=as.numeric(substr(names(r[[i]]),2,5))^2,
                                   lat=pts@coords[,2],
                                   lon=pts@coords[,1])

                       )
}

cat(" -- fitting robust regression model\n")
rl.model <- rlm(as.formula(paste(names(trainingData[1]), ".", sep=" ~ ")),scale.est="Huber", psi=psi.hampel, init="lts",data=trainingData)



# write our output to disk
writeRaster(satThick_10_14,"satThick_10_14.tif",overwrite=T)
