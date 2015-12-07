#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
# Find source data at: http://ne.water.usgs.gov/ogw/hpwlms/data.html
#
# Assuming +init=epsg:4269 for CRS.  No meta data available and "NAD83" with unit degrees is
# vague.
#



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
    if(sum(x$well_depth_ft<0,na.rm=T)>1){ x$well_depth_ft[x$well_depth_ft<0] <- 0; warning("censored non-sense well depth observations.") }
    x$saturated_thickness <- x$well_depth_ft-x$lev_va_ft
    # censor non-sense values
    if(sum(x$saturated_thickness<0,na.rm=T)>1){ x$saturated_thickness[x$saturated_thickness<0] <- 0; warning("censored non-sense saturated-thickness observations.") }

    return(x[!is.na(x$saturated_thickness),])
  }

  spatialPts <- function(x,usingPPP=F){
    pts <- SpatialPointsDataFrame(coords=data.frame(x=x$long_dd_NAD83,y=x$lat_dd_NAD83),data=x[,grepl(names(x),pattern="lev|depth|saturated")])
      projection(pts) <- projection("+init=epsg:4269")

    if(!usingPPP) {
      return(pts)
    } else {
      return(spatialPointsToPPP(pts))
    }
  }

  focalYear <- gsub(unlist(strsplit(x,split="_"))[3],pattern="[.]zip",replacement="")
  if(file.exists(paste("saturated_thickness.",focalYear,".tif",sep=""))){ cat(" -- found existing raster for year in CWD.  Quitting.\n"); return(FALSE) }

  cat(" -- parsing zip data for:",focalYear,"\n")
  t <- unpackZip(x)
    t <- calcSaturatedThickness(t)
      t <- spatialPts(t)
  p <- spatialPointsToPPP(t,extentMultiplier=1.05,field="saturated_thickness")

  cat(" -- building a raster template based on parsed point data\n")
  rasterTemplate <- raster(ext=multiplyExtent(t,extentMultiplier=1.05),crs=CRS(projection(t)),resolution=(500/111319.9))
          coords <- xyFromCell(rasterTemplate,cell=1:ncell(rasterTemplate))

  cat(" -- calculating inverse-distance surface\n")
  i <- raster(idw(p,at="pixels", eps=(500/111319.9),dimyx=c(rasterTemplate@nrows,rasterTemplate@ncols),xy=list(x=coords[,1],y=coords[,2])))
   projection(i) <- projection(rasterTemplate)
   
  # Area-weighting step, similar to V.L. McGuire (2013). 
  # Should validate this approach again Thiessen polygonization using cross-validation to see which approach works better.
  cat(" -- smothing with a 50x50 moving window\n")
  i <- focal(i, w=matrix(1, 51, 51), mean) 

  cat(" -- cropping/masking\n")
  b <- spTransform(readOGR("ds543/","hp_bound2010"),CRS(projection("+init=epsg:4269")))
    i <- mask(i,b)
  if(write){
    writeRaster(i,paste("saturated_thickness.",focalYear,".tif",sep=""),overwrite=T)
  } else {
    return(i)
  }
}

require(parallel)
cl <- makeCluster(8)
dataZips <- list.files(pattern="zip$")
parLapply(cl,as.list(dataZips),fun=processFocalYear,write=T)

