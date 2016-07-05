require(raster)
require(landscapeAnalysis)
require(rgdal)
require(rgeos)

argv <- commandArgs(trailingOnly=T)

s    <- readOGR(paste(Sys.getenv("HOME"),"/Incoming/",sep=""),"bcr_18_units",verbose=F)
r    <- raster(paste(Sys.getenv("HOME"),"/Incoming/Final_LC_8bit",".tif",sep=""))
#r <- raster("/home/ktaylora/Projects/large_block_sliding_window_analysis/plv_landcover_gplcc_ag_forecast_2025.tif")

s <- spTransform(s[seq(as.numeric(argv[1]),as.numeric(argv[2]),1),],CRS(projection(r)))

pts <- isolation_1650 <- isolation_16500 <- vector('list',nrow(s))

for(i in 1:length(pts)){
  pts[[i]] <- rgeos::gBuffer(spsample(s[i,],type="random",n=100),width=750,byid=T)
    pts[[i]] <- as(pts[[i]],'SpatialPolygonsDataFrame')
  isolation_1650[[i]] <- rgeos::gBuffer(gCentroid(s[i,]),width=1650,byid=T)
    isolation_1650[[i]] <- as(isolation_1650[[i]],'SpatialPolygonsDataFrame')
  isolation_16500[[i]] <- rgeos::gBuffer(gCentroid(s[i,]),width=16500,byid=T)
    isolation_16500[[i]] <- as(isolation_16500[[i]],'SpatialPolygonsDataFrame')
}

e <- as(landscapeAnalysis::multiplyExtent(extent(s),1.05),'SpatialPolygons')
  projection(e) <- projection(s)
    e <- spTransform(e,CRS(projection(r)))

agricultural_habitat <- c(38,202,201,203,204,205,206,207,208,209,210) # Cropland

shrubland_habitat <- c(82, # Shrubland
                       83, # Mesquite
                       81, # Savannah
                       85, # Shinnery
                       87) # Sand sage

grass_habitat <- c(39, # CRP
                   31, # CRP - GRASS
                   37, # Pasture
                   71, # Mixedgrass
                   75) # Shortgrass

t <- data.frame()

for(i in 1:length(pts)){
  cat(" -- cropping unit samples:")
  o <- landscapeAnalysis::cropRasterByPolygons(r,pts[[i]])
  cat(" -- cropping isolation sample (1650m):")
  o_1650 <- landscapeAnalysis::cropRasterByPolygons(r,isolation_1650[[i]])
    o_1650 <- lReclass(o_1650,inValues=grass_habitat)
  cat(" -- cropping isolation sample (16500m):")
  o_16500 <- landscapeAnalysis::cropRasterByPolygons(r,isolation_16500[[i]])
    o_16500 <- lReclass(o_16500,inValues=grass_habitat)

  grass_units <- lReclass(o,inValues=grass_habitat)
    grass_stats <- lCalculateLandscapeMetrics(grass_units,metric=c("total.area","mean.patch.area"))
      grass_stats[,2:3] <- grass_stats[,2:3] * (30*30) # num. cells -> meters-squared
        means <- apply(grass_stats[,2:3],2,FUN=mean)
      grass_stats <- data.frame(cbind(id=s$id[i],total.area=means[1],mean.patch.area=means[2]))

  ag_units <- lReclass(o,inValues=agricultural_habitat)
    ag_stats <- lCalculateLandscapeMetrics(ag_units,metric=c("total.area","mean.patch.area"))
      ag_stats[,2:3] <- ag_stats[,2:3] * (30*30)
        means <- apply(grass_stats[,2:3],2,FUN=mean)
      ag_stats <- data.frame(cbind(id=s$id[i],total.area=means[1],mean.patch.area=means[2]))

  sh_units <- lReclass(o,inValues=shrubland_habitat)
    sh_stats <- lCalculateLandscapeMetrics(sh_units,metric=c("total.area","mean.patch.area"))
      sh_stats[,2:3] <- sh_stats[,2:3] * (30*30)
        means <- apply(grass_stats[,2:3],2,FUN=mean)
      sh_stats <- data.frame(cbind(id=s$id[i],total.area=means[1],mean.patch.area=means[2]))

  o_1650 <- try(landscapeAnalysis::calcPatchIsolation(o_1650))
    if(class(o_1650) == "try-error"){
      o_1650 <- NA;
    }
  o_16500 <- try(landscapeAnalysis::calcPatchIsolation(o_16500))
    if(class(o_16500) == "try-error"){
      o_16500 <- NA;
    }

  is <- data.frame(cbind(is_1650=o_1650, is_16500=o_16500))

  t <- rbind(t,
             cbind(grass_stats,
             ag_stats[,!grepl(names(ag_stats),pattern="id")],
             sh_stats[,!grepl(names(sh_stats),pattern="id")],
             is[,!grepl(names(is),pattern="id")])
            );
}


rownames(t) <- c()
  names(t) <- <- c("id","tAr_gr","pAr_gr","tAr_ag","pAr_ag","tAr_sh","pAr_sh","is_1650","is_16500")

t$id <- s$id

write.csv(t,paste("results_",argv[1],"-",argv[2],".csv",sep=""),row.names=F)
