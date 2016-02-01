#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2015]
# Find source data at: http://ne.water.usgs.gov/ogw/hpwlms/data.html
#
# Assuming +init=epsg:4269 for CRS.  No meta data available and "NAD83" with unit degrees is
# vague.  In the first build, I designed a really hideous function-of-functions (observed below) to ease parallelization.
# Future Kyle, if you happen across this and need to edit it -- I am sorry.

#
# WORKFLOW : PARSE RAW AQUIFER WELL DATA AND INTERPOLATE / PROCESS A FIXED RASTER FOR THE HIGH PLAINS AQUIFER
#

processFocalYear <- function(x,write=F){
  require(spatstat)
  require(raster)
  require(rgdal)
  require(utils)
  # LOCAL FUNCTIONS
  unpackZip <- function(x){
    f <- utils::unzip(x,list=T)[,1];
      utils::unzip(x);
    f <-read.csv(f,sep="\t",comment.char="#",stringsAsFactors=FALSE)
      f[f<=-999] <- NA
         f[f>=9999] <- NA
    return(f)
  }

  multiplyExtent <- function(x,extentMultiplier=1.1){
    e <- extent(x)

    if(!is.null(extentMultiplier)) {
      e@xmin <- e@xmin*extentMultiplier
      e@xmax <- e@xmax+abs(e@xmax*(extentMultiplier-1))
      e@ymin <- e@ymin-abs(e@ymin*(extentMultiplier-1))
      e@ymax <- e@ymax*extentMultiplier
    }

    return(e)
  }

  as.owin <- function(x,extentMultiplier=1.1){
    e <- multiplyExtent(x,extentMultiplier=extentMultiplier)
    return(owin(xrange=c(e@xmin,e@xmax), yrange=c(e@ymin,e@ymax)))
  }

  spatialPointsToPPP <- function(x,extentMultiplier=1.1,field=NULL){
    # attribute 'data' to 'marks' for our PPP
    if(grepl(class(x),pattern="SpatialPoints")){
      d <- x@data
      c <- x@coords
      if(grepl(class(x),pattern="SpatialPointsDataFrame")){
        if(!is.null(field)){
          x <- ppp(x=c[,1], y=c[,2], window=as.owin(x,extentMultiplier=extentMultiplier), marks=d[,field])
        } else {
          x <- ppp(x=c[,1], y=c[,2], window=as.owin(x,extentMultiplier=extentMultiplier), marks=d)
        }
      } else {
        x <- ppp(x=c[,1], y=c[,2], window=as.owin(x,extentMultiplier=extentMultiplier))
      }
    }
    return(x)
  }

  calcSaturatedThickness <- function(x){
    # censor non-sense values

    if(sum(x$lev_va_ft<0,na.rm=T)>1){ x$lev_va_ft[x$lev_va_ft<0] <- 0; warning("censored non-sense depth to water observations.") }
    if(sum(x$bedrock<0,na.rm=T)>1){ x$bedrock[x$bedrock<0] <- 0; warning("censored non-sense bedrock observations.") }
    x$saturated_thickness <- x$bedrock-x$lev_va_ft
    # censor non-sense values
    if(sum(x$saturated_thickness<0,na.rm=T)>1){ x$saturated_thickness[x$saturated_thickness<0] <- 0; warning("censored non-sense saturated-thickness observations.") }

    return(x[!is.na(x$saturated_thickness),])
  }

  spatialPts <- function(x,usingPPP=F){
    x <- x[!is.na(x$long_dd_NAD83),]
      x <- x[!is.na(x$lat_dd_NAD83),]
    pts <- SpatialPointsDataFrame(coords=data.frame(x=x$long_dd_NAD83,y=x$lat_dd_NAD83),data=x[,grepl(names(x),pattern="lev|depth|saturated|bed")])
      projection(pts) <- projection("+init=epsg:4269")

    if(!usingPPP) {
      return(pts)
    } else {
      return(spatialPointsToPPP(pts))
    }
  }

  interpolateBedrockDepth <- function(){
    require(spatstat)
    require(raster)
    require(rgdal)
    cat(" -- extracting bedrock depths at well point locations\n")
      s <- spatialPts(unpackZip("WL_ALL_2009.zip"))
    sat <- raster("hp_satthk09")

    o <- extract(sat,s,df=T)[,2]
    o <- (o-(s$well_depth_ft-s$lev_va_ft))+s$well_depth_ft

    t <- unpackZip("WL_ALL_2009.zip")
      t <- t[!is.na(t$long_dd_NAD83),]
        t <- t[!is.na(t$lat_dd_NAD83),]
          t$bedrock <- o

    t <- spatialPts(t)
      t <- t[!is.na(t$bedrock),]

    p <- spatialPointsToPPP(t,extentMultiplier=1.05,field="bedrock")
    cat(" -- building a raster template based on parsed point data\n")
    rasterTemplate <- raster(ext=multiplyExtent(t,extentMultiplier=1.05),crs=CRS(projection(t)),resolution=(500/111319.9))
            coords <- xyFromCell(rasterTemplate,cell=1:ncell(rasterTemplate))

    cat(" -- calculating inverse-distance surface\n")
    i <- raster(idw(p,at="pixels", eps=(500/111319.9),dimyx=c(rasterTemplate@nrows,rasterTemplate@ncols),xy=list(x=coords[,1],y=coords[,2])))
     projection(i) <- projection(rasterTemplate)
       i <- focal(i, w=matrix(1, 51, 51), mean)
         writeRaster(i,"bedrock.tif")
  }

  focalYear <- gsub(unlist(strsplit(x,split="_"))[3],pattern="[.]zip",replacement="")
  if(file.exists(paste("saturated_thickness.",focalYear,".tif",sep=""))){ cat(" -- found existing raster for year in CWD.  Quitting.\n"); return(FALSE) }

  cat(" -- parsing zip data for:",focalYear,"\n")
  t <- unpackZip(x)
    t <- t[!is.na(t$long_dd_NAD83),]
      t <- t[!is.na(t$lat_dd_NAD83),]
  cat(" -- extracting bedrock depth\n")
  if(!file.exists("bedrock.tif")){ interpolateBedrockDepth() }
  t$bedrock <- extract(raster("bedrock.tif"),spatialPts(t),df=T)[,2]
  t <- calcSaturatedThickness(t)
    t <- spatialPts(t)
      writeOGR(t,".",focalYear,driver="ESRI Shapefile",overwrite=T)

  # EXPORT SHAPEFILES TO ARC
  #p <- spatialPointsToPPP(t,extentMultiplier=1.05,field="saturated_thickness")

  #cat(" -- building a raster template based on parsed point data\n")
  #rasterTemplate <- raster(ext=multiplyExtent(t,extentMultiplier=1.05),crs=CRS(projection(t)),resolution=(500/111319.9))
  #        coords <- xyFromCell(rasterTemplate,cell=1:ncell(rasterTemplate))

  #cat(" -- calculating inverse-distance surface\n")
  #i <- raster(idw(p,at="pixels", eps=(500/111319.9),dimyx=c(rasterTemplate@nrows,rasterTemplate@ncols),xy=list(x=coords[,1],y=coords[,2])))
  # projection(i) <- projection(rasterTemplate)

  # Area-weighting step, similar to V.L. McGuire (2013).
  # Should validate this approach again Thiessen polygonization using cross-validation to see which approach works better.
  #cat(" -- smothing with a 50x50 moving window\n")
  #i <- focal(i, w=matrix(1, 51, 51), mean)

  #cat(" -- cropping/masking\n")
  #b <- spTransform(readOGR("ds543/","hp_bound2010"),CRS(projection("+init=epsg:4269")))
  #  i <- mask(i,b)
  #if(write){
  #  writeRaster(i,paste("saturated_thickness.",focalYear,".tif",sep=""),overwrite=T)
  #} else {
  #  return(i)
  #}
}

require(parallel)
cl <- makeCluster(8)
dataZips <- list.files(pattern="zip$")
  dataZips <- dataZips[grepl(dataZips,pattern="WL_ALL")]
parLapply(cl,as.list(dataZips),fun=processFocalYear,write=T)
