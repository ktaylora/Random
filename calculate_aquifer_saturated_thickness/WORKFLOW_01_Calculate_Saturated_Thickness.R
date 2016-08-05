#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2015]
# Find source data at: http://ne.water.usgs.gov/ogw/hpwlms/data.html
#
# Assuming +init=epsg:4269 for CRS.  No meta data available and "NAD83" with unit degrees is
# vague.  In the first build, I designed a really hideous function-of-functions (observed below) to ease parallelization.
# Future Kyle, if you happen across this and need to edit it -- I am sorry.

#
# LOCAL INCLUDES
#
require(parallel)           # For parallization
require(landscapeAnalysis)  # Kyle's handy package for various spatial things

#
# LOCAL FUNCTIONS
#
aquiferAnalysis <- new.env()
#
# fetchWellPointData()
#
fetchWellPointData <- function(){
  baseUrl <- "http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_"
    years <- 1995:2013
    toGet <- zips  <- paste("WL_ALL_",years,".zip",sep="")

  nExistingZips <- length(list.files(pattern="^WL_ALL_.*.zip$"))
  if(nExistingZips<length(years)){
    if(nExistingZips > 0){
      toGet <- list.files(pattern="^WL_ALL_.*.zip$")
        years <- years[which(!grepl(zips,pattern=paste(toGet,collapse="|")))]
          toGet <- zips[!grepl(zips, pattern=paste(toGet,collapse="|"))]
    }
    fetchURLs <- paste(baseUrl,years,"/",toGet,sep="")
      for(i in 1:length(fetchURLs)){ utils::download.file(fetchURLs[i],destfile=toGet[i]) }
  }
  return(list.files(pattern="^WL_ALL_.*.zip$"))
}
#
# unpackZip()
#
aquiferAnalysis$unpackZip <- function(x){
  f <- utils::unzip(x,list=T)[,1];
    utils::unzip(x);
  f <-read.csv(f,sep="\t",comment.char="#",stringsAsFactors=FALSE)
    f[f<=-999] <- NA # looming NA values in the source CSV have inconsistently large values. Treat as NA's
       f[f>=9999] <- NA
  return(f)
}
#
# calcSaturatedThickness()
#
# this is outdated and could use some refactorin' so that it is actually useful.
#
aquiferAnalysis$calcSaturatedThickness <- function(x){
  # censor non-sense values

  if(sum(x$lev_va_ft<0,na.rm=T)>1){ x$lev_va_ft[x$lev_va_ft<0] <- 0; warning("censored non-sense depth to water observations.") }
  if(sum(x$bedrock<0,na.rm=T)>1){ x$bedrock[x$bedrock<0] <- 0; warning("censored non-sense bedrock observations.") }
  x$saturated_thickness <- x$bedrock-x$lev_va_ft
  # censor non-sense values
  if(sum(x$saturated_thickness<0,na.rm=T)>1){ x$saturated_thickness[x$saturated_thickness<0] <- 0; warning("censored non-sense saturated-thickness observations.") }

  return(x[!is.na(x$saturated_thickness),])
}
#
# spatialPts()
#
aquiferAnalysis$spatialPts <- function(x,usingPPP=F){
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
#
# interpolateBedrockDepth()
#
aquiferAnalysis$interpolateBedrockDepth <- function(){
  require(spatstat)
  require(raster)
  require(rgdal)
  cat(" -- extracting bedrock depths at well point locations\n")
    s <- spatialPts(unpackZip("WL_ALL_2009.zip")) # use well data for 2009 for our calculation
  sat <- raster("hp_satthk09") # Use the published saturated thickness (2009) product from V. McGuire to back-calculate bedrock depth for the region

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
#' Use kriging and back-filling to interpolate a saturated thickness surface from well point data
#' and a historical (reference) surface representing saturated thickness in areas currently unsampled
aquiferAnalysis$interpolateSaturatedThickness <- function(s){
  # generate a large number of random points within the spatial extent of our s= well points
  background_pts <- raster::sampleRandom(raster(extent(s),vals=1,nrow=20000,ncol=20000),size=90000,sp=T)
    projection(background_pts) <- projection(s)
  # select only those random points that are within the boundary of the HP Aquifer
  background_pts <- background_pts[!is.na(sp::over(background_pts,
    spTransform(rgdal::readOGR("hp_bound2010","hp_bound2010",verbose=F),
      CRS(projection(background_pts)))))[,1],]
  # back-fill missing areas with values from a historical depletion surface : generate random background points with historical sat_thickness
  background <- raster::raster("/home/ktaylora/Incoming/aquifer_saturated_thickness_products/Raster/satThick_10_14.tif")
  background_pts <- sp::spTransform(background_pts,CRS(projection(background)))
    background_pts$strtd_t <- raster::extract(background,background_pts)
  # buffer the initial well points and ensure that our background points
  # are outside of that buffer distance
  s_buffered <- suppressWarnings(rgeos::gBuffer(s,width=0.035))
  doesnt_overlap <- is.na(sp::over(background_pts,
                                   spTransform(s_buffered,CRS(projection(background_pts))))
                       )
  background_pts <- background_pts[doesnt_overlap,]
  # merge the two geometries together, keeping only the saturated thickness field
  s@data <- data.frame(strtd_t=s@data[,'strtd_t'])
  background_pts@data <- data.frame(strtd_t=background_pts@data[,'strtd_t'])
  pts <- raster::bind(s,sp::spTransform(background_pts,CRS(projection(s))))
    pts$x <- pts@coords[,1]
    pts$y <- pts@coords[,2]
  pts <- pts[which(!is.na(rowSums(pts@data))),] # drop NA values, or gstat::variogram will complain
  cat(" -- computing the (quadratic) variogram of our points")
  f <- as.formula(strtd_t ~ x + y + I(x*x)+I(y*y) + I(x*y))
  variogram  <- gstat::variogram(f, data=pts[sample(1:nrow(pts),size=round(0.25*nrow(pts))),], cloud = TRUE)
}
#
# calcSaturatedThickness_byYear()
# Parse raw USGS aquifer well data into a shapefile that we can hand-off to ArcGIS. This will generate spatial points
# with saturated thickness indicated for each year.
#
calcSaturatedThickness_byYear <- function(x,write=F,env=aquiferAnalysis){

  require(spatstat)
  require(raster)
  require(rgdal)
  require(utils)
  require(landscapeAnalysis)

  attach(env)

  # unpack the zipfile passed by the user for the focal year
  focalYear <- gsub(unlist(strsplit(x,split="_"))[3],pattern="[.]zip",replacement="")
    if(focalYear == "pd") focalYear <- "historical" # Catch 'PD' string from historical well point data
  if(file.exists(paste("saturated_thickness.",focalYear,".tif",sep=""))){ cat(" -- found existing raster for year in CWD.  Quitting.\n"); return(FALSE) }
  cat(" -- parsing zip data for:",focalYear,"\n")
  t <- unpackZip(x)
    t <- t[!is.na(t$long_dd_NAD83),]
      t <- t[!is.na(t$lat_dd_NAD83),]
  # do a back-of-the-envelope calculation of bedrock depth using the published saturated thickness product from 2008 from V. McGuire
  cat(" -- extracting bedrock depth\n")
  if(!file.exists("bedrock.tif")){ interpolateBedrockDepth() }
  t$bedrock <- extract(raster("bedrock.tif"),spatialPts(t),df=T)[,2]
  # calculate saturated thickness for our spatial points of well depth and write to disk
  t <- calcSaturatedThickness(t)
    t <- spatialPts(t)
      writeOGR(t,".",focalYear,driver="ESRI Shapefile",overwrite=T)
}

#
# MAIN WORKFLOW
#

cl <- makeCluster(parallel::detectCores()-1)
  dataZips <- fetchWellPointData()
    s <- parLapply(cl,as.list(dataZips),fun=calcSaturatedThickness_byYear,write=T,env=aquiferAnalysis)
