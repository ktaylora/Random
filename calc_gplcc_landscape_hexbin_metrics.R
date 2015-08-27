#
# calculate landscape metrics per hexagonal unit across the GPLCC pilot region
#
# Author : Kyle Taylor (kyle.taylor@pljv.org)
# currently takes around 20 seconds per unit, per core
#

require(rgdal)
require(raster)

require(landscapeAnalysis)

#
# LOCAL FUNCTIONS 
#

calculateLMetricsForFocalCover <- function(sampleUnits=NULL,habitat=NULL,metrics=NULL,bufferDistance=750){
    t <- data.frame(matrix(nrow=length(sampleUnits), ncol=length(metrics)+1))
    # iterate over units, sampling as we go
    cat("[sampling",length(sampleUnits),"units]\n");
    for(i in 1:length(sampleUnits)){
      out <- subsampleSurface(x=sampleUnits[[i]], n=100, width=bufferDistance) # keep the same buffer width that was used for our route-level analysis
        out <- lReclass(out, inValues=habitat)  
      # iterate over our samples, calculating statistics as we go
      row       <- rep(NA,length(metrics));
      lm_output <- lCalculateLandscapeMetrics(x=out, metric=metrics, DEBUG=F)
        names(lm_output) <- metrics
      
      for(j in 1:length(metrics)){ row[j] <- median(metricsListToVector(lm_output, metrics[j]),na.rm=T) }
      row <- c(i,row);
      
      # contribute to table
      t[i,] <- row;
    }; cat("\n");
    names(t) <- c("id",metrics)
    return(t);
}

parseStoredMetricsForRanges <- function(path=".") {                                                               
  f        <- list.files(path, pattern="rdata$",full.names=T)
  ranges   <- unlist(strsplit(f,split="-"))
    ranges <- as.numeric(gsub("[^0-9]", "", ranges))
    j      <- seq(1,length(ranges),2) # jumps between range numbers
    
  getEnv   <- function(x,var=NULL,na.rm=F) { 
    stored <- new.env(); 
    load(x,envir=stored); 
    seq<-as.numeric(gsub(x=unlist(strsplit(x,split="-")),pattern="\\D", ""));
    t <- get(var,envir=stored);
    t$id <- seq(seq[1],seq[2],1);
    return(t) 
  }
  
  grass    <- lapply(as.list(f), FUN=getEnv, var="t_landcover_grass")
  ag       <- lapply(as.list(f), FUN=getEnv, var="t_landcover_ag")
  shrub    <- lapply(as.list(f), FUN=getEnv, var="t_landcover_shrub")
  
  # re-assign id's based on the names associated with run data (.rdata) in the CWD.  They were not recorded properly internally
  grass <- do.call(rbind,grass)
    names(grass) <- paste("grass_",names(grass),sep="");
  
  ag <- do.call(rbind,ag)
   names(ag) <- paste("ag_",names(ag),sep="");
  
  shrub <- do.call(rbind,shrub)
   names(shrub) <- paste("shrub_",names(shrub),sep="");
  
  t<-cbind(ag[,'ag_id'],ag[,names(ag)!='ag_id'],grass[,names(grass)!="grass_id"],shrub[,names(shrub) != "shrub_id"])
    names(t) <- c("id",names(ag)[names(ag)!='ag_id'],names(grass)[names(grass)!="grass_id"],names(shrub)[names(shrub)!="shrub_id"])
      return(t)
}

#
# MAIN
#

s    <- readOGR(paste(Sys.getenv("HOME"),"/PLJV/boundaries/GPLCC Pilot/GPLCC_Pilot_Region/",sep=""),"pilot_region_hexagonal_units",verbose=F)
r    <- raster(paste(Sys.getenv("HOME"),"/PLJV/landcover/Final_LC_8bit",".tif",sep=""))

argv <- commandArgs(trailingOnly=T)

stored <- new.env();
if(file.exists(paste(argv[1],"-",argv[2],"_so.far.rdata",sep=""))) { 
  cat(" -- loading previous session data\n");
  load(paste(argv[1],"-",argv[2],"_so.far.rdata",sep=""), envir=stored)
}

if(as.numeric(argv[2]) > nrow(s)){ argv[2] <- nrow(s) }

unitRange           <- seq(as.numeric(argv[1]),as.numeric(argv[2]),1)
hexSampleUnits      <- cropRasterByPolygons(r=r, s=s[unitRange,])
sampleUnitCentroids <- NULL
ids                 <- as.vector(s$id)

# define our cover types 
grass_habitat <- c(39, # CRP
                   37, # Pasture
                   71, # Mixedgrass
                   75) # Shortgrass

agricultural_habitat <- c(38,202,201,203,204,205,206,207,208,209,210) # Cropland

shrubland_habitat <- c(82, # Shrubland
                       83, # Mesquite
                       81  # Savannah
                       85, # Shinnery
                       87) # Sand sage

if(!exists("t_landcover_grass",envir=stored)){
  cat(" -- calculating grassland fragstats\n");
  # calculate landscape metrics for hexagonal units 
  t_landcover_grass <- calculateLMetricsForFocalCover(hexSampleUnits,habitat=grass_habitat, metrics=c("total.area","mean.patch.area","shape.index"), bufferDistance=750);

  assign(x="t_landcover_grass",value=t_landcover_grass, envir=stored);  
  save(list=ls(stored),file=paste(argv[1],"-",argv[2],"_so.far.rdata",sep=""),envir=stored,compress=T);
}

#if(!exists("t_connectivity_grass",envir=stored)){
#cat(" -- calculating grassland connectivity\n");
#  sampleUnitCentroids <- SpatialPoints(coords=getSpPPolygonsLabptSlots(s[unitRange,]),proj4string=CRS(projection(s)));
#  t_connectivity_grass <- subsampleSurface(x=r, pts=sampleUnitCentroids,width=20000);
#    t_connectivity_grass <- lReclass(t_connectivity_grass, inValues=grass_habitat)
#      t_connectivity_grass <- lCalculateLandscapeMetrics(x=t_connectivity_grass, metric="patch.issolation",DEBUG=F);
#
#  assign(x="t_connectivity_grass",value=t_connectivity_grass, envir=stored);
#  save(list=ls(stored),file=paste(argv[1],"-",argv[2],"_so.far.rdata",sep=""),envir=stored,compress=T);
#}

if(!exists("t_landcover_ag",envir=stored)){
  cat(" -- calculating ag fragstats\n")
  # calculate landscape metrics for hexagonal units 
  t_landcover_ag <- calculateLMetricsForFocalCover(hexSampleUnits,habitat=agricultural_habitat, metrics=c("total.area","mean.patch.area","shape.index"), bufferDistance=750);

  assign(x="t_landcover_ag",value=t_landcover_ag, envir=stored);
  save(list=ls(stored),file=paste(argv[1],"-",argv[2],"_so.far.rdata",sep=""),envir=stored, compress=T);
}

#if(!exists("t_connectivity_ag",envir=stored)){
#cat(" -- calculating agriculture connectivity\n");
#  sampleUnitCentroids <- SpatialPoints(coords=getSpPPolygonsLabptSlots(s[unitRange,]),proj4string=CRS(projection(s)));
#  t_connectivity_ag <- subsampleSurface(x=r, pts=sampleUnitCentroids,width=20000);
#    t_connectivity_ag <- lReclass(t_connectivity_ag, inValues=agricultural_habitat)
#      t_connectivity_ag <- lCalculateLandscapeMetrics(x=t_connectivity_ag, metric="patch.issolation",DEBUG=F);
#
#  assign(x="t_connectivity_ag",value=t_connectivity_ag, envir=stored);
#  save(list=ls(stored),file=paste(argv[1],"-",argv[2],"_so.far.rdata",sep=""),envir=stored,compress=T);
#}

if(!exists("t_landcover_shrub",envir=stored)){
  cat(" -- calculating shrub fragstats\n");
  # calculate landscape metrics for hexagonal units 
  t_landcover_shrub <- calculateLMetricsForFocalCover(hexSampleUnits,habitat=shrubland_habitat, metrics=c("total.area","mean.patch.area","shape.index"), bufferDistance=750);

  assign(x="t_landcover_shrub",value=t_landcover_shrub, envir=stored);
  save(list=ls(stored),file=paste(argv[1],"-",argv[2],"_so.far.rdata",sep=""),envir=stored,compress=T);
}

#if(!exists("t_connectivity_shrub",envir=stored)){
#cat(" -- calculating shrub connectivity\n");
#  sampleUnitCentroids <- SpatialPoints(coords=getSpPPolygonsLabptSlots(s[unitRange,]),proj4string=CRS(projection(s)));
#  t_connectivity_shrub <- subsampleSurface(x=r, pts=sampleUnitCentroids,width=20000);
#    t_connectivity_shrub <- lReclass(t_connectivity_shrub, inValues=shrubland_habitat)
#      t_connectivity_shrub <- lCalculateLandscapeMetrics(x=t_connectivity_shrub, metric="patch.issolation",DEBUG=F);
#
#  assign(x="t_connectivity_shrub",value=t_connectivity_shrub, envir=stored);
#  save(list=ls(stored),file=paste(argv[1],"-",argv[2],"_so.far.rdata",sep=""),envir=stored,compress=T);
#}
