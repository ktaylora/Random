require(raster)
require(rgdal)

normalize <- function(t,fun=median) {
  mean <- apply(t,MARGIN=2, FUN=fun, na.rm=T)
    sd <- apply(t,MARGIN=2, FUN=sd, na.rm=T)

  n <- function(x) (x-mean)/sd
    return(t(apply(t,1,n)))
}

applyWeights <- function(t,weights=c(1,1,1,1)){
  t <- t(apply(t,1,function(x) x*weights))
    return(t)
}

calculateScores <- function(s){
  t <- s[,!grepl(names(s),pattern="score|id|sc")]@data
    t$dstToRm <- (1/sqrt(t$dstToRm+1)) # inverse-distance to Rio Mora
      t <- normalize(t)

  # scores, weighted by 1/distance from Rio Mora
  scores <- applyWeights(t,c(1,0.5,0.5,0.5))
    s$scDstRm <- as.vector(rowSums(scores))
  # scores, weighted by area
  scores <- applyWeights(t,c(0.5,1,0.5,0.5))
    s$scArea <- as.vector(rowSums(scores))
  # scores, weighted by percent grass cover
  scores <- applyWeights(t,c(0.5,0.5,1,0.5))
    s$scPrcGr <- as.vector(rowSums(scores))
  # scores, weighted by stream fetch length
  scores <- applyWeights(t,c(0.5,0.5,0.5,1))
    s$scStFch <- as.vector(rowSums(scores))

  return(s)
}

vectorToRaster <- function(s,county_prefix=NULL){
  # convert each score to a raster object
  template <- raster(extent(s),resolution=30,crs=CRS(projection(s)))

  scDstRm <- rasterize(s,template,field='scDstRm',progress='text')
    writeRaster(scDstRm, paste("score_distance_to_rio_mora_-",county_prefix,".tif",sep=""), overwrite=T)
  scArea <- rasterize(s,template,field='scArea',progress='text')
    writeRaster(scArea, paste("score_parcel_area_-",county_prefix,".tif",sep=""), overwrite=T)
  scPrcGr <- rasterize(s,template,field='scPrcGr',progress='text')
    writeRaster(scPrcGr, paste("score_percent_grass_-",county_prefix,".tif",sep=""), overwrite=T)
  scStFch <- rasterize(s,template,field='scStFch',progress='text')
    writeRaster(scStFch, paste("score_stream_fetch_-",county_prefix,".tif",sep=""), overwrite=T)
}

s <- readOGR(".","mora_county_parcels")
  s <- calculateScores(s)
    vectorToRaster(s, county_prefix="mora_county")

s <- readOGR(".","san_miguel_county_parcels")
  s <- calculateScores(s)
    vectorToRaster(scounty_prefix="san_miguel_county")

s <- readOGR(".","colfax_county_parcels")
  s <- calculateScores(s)
    vectorToRaster(scounty_prefix="colfax_county")

# writeOGR(s,layer="tri_county_area_landowner_boundaries_v2",
#   "tri_county_area_landowner_boundaries_v2.geoJson", driver="GeoJSON",
#      overwrite=T)
