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

# clean-up previous workspace
if(sum(file.remove(list.files(pattern="resid|coef|pred")))>0) cat(" -- flushing previous raster products from working directory\n")

# mask our rasters against the boundaries of the HP aquifer
s <- readOGR(argv[1],"hpbedrock",verbose=F)
r <- list.files(argv[1],pattern="tif$")
  r <- r[!grepl(r,pattern="satThick|pred|coef")]
    r <- lapply(as.list(r),FUN=raster)

s <- spTransform(s,CRS(projection(r[[1]])))
  r <- parLapply(cl,r,mask,mask=s)

# calculate mean saturated thickness for a time-series consistent with NASS crop data
satThick_10_14 <- stackApply(stack(r[grepl(unlist(lapply(r,names)),pattern="2010|2011|2012|2013|2014")]),fun=mean,indices=1)

# build a table of trainable values and calculate a robust time-series regression for our aquifer depletion time-series
pts <- raster::sampleRandom(r[[1]],size=ncell(r[[1]])*0.35,sp=T,na.rm=T)

trainingData <- data.frame(sat.thickness=extract(r[[1]],pts),
                           yr=as.numeric(substr(names(r[[1]]),2,5)),
                           lat=pts@coords[,2],
                           lon=pts@coords[,1])

for(i in 2:length(r)){
  trainingData <- rbind(trainingData,
                        data.frame(sat.thickness=extract(r[[i]],pts),
                                   yr=as.numeric(substr(names(r[[i]]),2,5)),
                                   lat=pts@coords[,2],
                                   lon=pts@coords[,1])

                       )
}

cat(" -- fitting robust regression model to randomized samples (for evaluating variable importance and influence of space)\n")
n <- names(trainingData)
  rl.model <- rlm(as.formula(paste(names(trainingData[1]), "~I(", n[2], ")^2+I(",n[2],")+",n[3],"+",n[4],sep="")),scale.est="Huber", psi=psi.hampel, init="lts",data=trainingData)
    print(summary(rl.model)); cat("\n");

cat(" -- fitting implicit, non-spatial regression to raster time-series\n")
time <- 1:length(r);

calcCoefficients <- function(y) { # stolen and modified from R. Hijmans --
    if(all(is.na(y))) {
      c(NA, NA)
    } else {
      lm(y ~ time)$coefficients
    }
}

# Fit a model to our raster time-series, keeping the coefficient values so we can project a raster
r_coef <- calc(stack(r), calcCoefficients) # will return a stack as intercept[1],slope[2](time)

# Residuals = observed - expected (masked against NAs)
r_resid <- unlist(lapply(r,FUN=names))
  r_resid <- round(median(as.numeric(substr(r_resid,2,5))))
    r_resid <- (stackApply(stack(r),indices=1,fun=median)-(r_coef[[2]]*setValues(r[[1]],r_resid)+r_coef[[1]]))*!is.na(r[[1]])

# Extrapolate out over a short time-period consistent with the range of our training data
cat(" -- extrapolating to a novel time-period\n")
target_year <- diff(range(as.numeric(substr(unlist(lapply(r,FUN=names)),2,5))))+max(as.numeric(substr(unlist(lapply(r,FUN=names)),2,5)))

r_predicted <- (r_coef[[2]]*setValues(r[[1]],target_year)+r_coef[[1]])*!is.na(r[[1]])
  #r_predicted[r_predicted<0] <- 0 # let's assume projected negative thickness is ~0

# Write our output to disk
cat(" -- writing rasters to disk\n")
writeRaster(satThick_10_14,"satThick_10_14.tif",overwrite=T)
  writeRaster(r_coef,paste(paste(range(as.numeric(substr(unlist(lapply(r,FUN=names)),2,5))),collapse="_"),"coef_time.tif",sep="_"),overwrite=T)
    writeRaster(r_resid,paste(paste(range(as.numeric(substr(unlist(lapply(r,FUN=names)),2,5))),collapse="_"),"resid_error.tif",sep="_"),overwrite=T)
      writeRaster(r_predicted,paste(target_year,"predicted.tif",sep="_"),overwrite=T)
