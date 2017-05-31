#
# CP1/2 Grassland Productivity Summary Analysis
# This workflow details a pre-cursor analysis to looking at the feasibility of
# using CP1/2 treatments in a region-wide IMBCR analysis
#
# Author : Kyle Taylor (kyle.taylor@pljv.org) [2017]
#

require(raster)
require(rgdal)
require(rgdal)
require(ggplot2)

RASTER_DIR <- "/home/ktaylora/Workspace/tpw_imbcr_grassland_birds_workflow/Raster"

calc_seasonal_ndvi_by <- function(s=NULL,r=NULL){
  s = sp::spTransform(s,sp::CRS(raster::projection(r)))
  return(raster::extract(r,s,na.rm=T,fun=mean,progress='text',df=T))
}

sort_ndvi_rasters_by_season <- function(path=NULL, pattern=NULL){
  f <- list.files(path, pattern=pattern,full.names=T)
  spring <- f[grep(x=f,pattern="spring")[1]]
  summer <- f[grep(x=f,pattern="summer")[1]]
  fall   <- f[grep(x=f,pattern="fall")[1]]
  winter <- f[grep(x=f,pattern="winter")[1]]
  f <- raster::stack(c(spring,summer,fall,winter))
    names(f) <- c("spring","summer","fall","winter")
  return(f)
}

calc_transect_centroids <- function(s=NULL){
  centroids <- s[!duplicated(s$transectnum),]
  centroids <- rgeos::gCentroid(centroids,byid=T,id=centroids$transectnum)
  centroids <- sp::SpatialPointsDataFrame(centroids,data=data.frame(transectnum=s[!duplicated(s$transectnum),'transectnum']))
  centroids@data <- data.frame(transectnum=centroids@data[,1])
  return(centroids)
}

calc_distance_to_features <- function(s=NULL, to=NULL, converter=6.21371e-4){
  transect_dist_to <- rgeos::gDistance(calc_transect_centroids(s),
    to, byid=T)
  transect_dist_to <- data.frame(distance=apply(transect_dist_to, MARGIN=2, FUN=function(x){
      round(min(x*converter),2) # meters -> miles
    }))
  return(transect_dist_to)
}
#' we have to talk about your function names...
calc_ndvi_at_treatment_sample_pts <- function(pattern="cp1.*.ndvi.*.tif$", write=NULL){
  cp_seasons <- sort_ndvi_rasters_by_season(RASTER_DIR, pattern=pattern)
    cp_seasonal_points <- sampleRandom(cp_seasons,size=99999,sp=T,na.rm=T)
  if(!is.null(write)){
    writeOGR(cp_seasonal_points, ".", layer=write, driver="ESRI Shapefile")
  }
  return(cp_seasonal_points)
}
#
# Grab IMBCR data
#

s <- spTransform(OpenIMBCR:::imbcrTableToShapefile(OpenIMBCR:::recursiveFindFile(
  "master_imbcr_table.csv")), CRS(projection("+init=epsg:2163")))

#
# Read-in treatment vector datasets
#

s_cp_1 <- sp::spTransform(readOGR("Vector","cp_1", verbose=F),
  CRS(projection("+init=epsg:2163")))
s_cp_2 <- sp::spTransform(readOGR("Vector","cp_2", verbose=F),
  CRS(projection("+init=epsg:2163")))

#
# Calculate distance to treatment vectors
#

transect_dist_to_cp_1 <- calc_distance_to_features(s=s, to=s_cp_1)
transect_dist_to_cp_2 <- calc_distance_to_features(s=s, to=s_cp_2)

dev.new()
ggplot(transect_dist_to_cp_1, aes(distance)) +
  geom_histogram(aes(y=..density..), color="grey",fill="#879EAD", bins=50) +
  geom_density(color="#1F78B4") +
  xlab("Distance from Transect to Nearest CP2 Field (miles)") +
  theme(plot.title=element_text(size=16),
        axis.title.x=element_text(size = 16),
        axis.title.y=element_text(size = 16),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size = 14))

dev.new()
ggplot(transect_dist_to_cp_2, aes(distance)) +
  geom_histogram(aes(y=..density..), color="grey",fill="#879EAD", bins=50) +
  geom_density(color="#1F78B4") +
  ylab("Sample Density") +
  xlab("Distance from Transect to Nearest CP2 Field (miles)") +
  theme(plot.title=element_text(size=16),
        axis.title.x=element_text(size = 16),
        axis.title.y=element_text(size = 16),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size = 14))

#
# Calculate NDVI summary statistics
#

# cp1_seasons <- sort_ndvi_rasters_by_season(RASTER_DIR, pattern="cp1.*.ndvi.*.tif$")
#   cp1_seasonal_points <- sampleRandom(cp1_seasons,size=99999,sp=T,na.rm=T)
#     writeOGR(cp1_seasonal_points, ".", "cp1_seasonal_points", driver="ESRI Shapefile")

calc_ndvi_at_treatment_sample_pts(pattern="cp1.*.ndvi.*.tif$", write="cp1_seasonal_points")
  #cp1_seasons <- calc_seasonal_ndvi_by(s_cp_1,r=cp1_seasons)
calc_ndvi_at_treatment_sample_pts(pattern="cp2.*.ndvi.*.tif$", write="cp2_seasonal_points")
  #cp2_seasons <- calc_seasonal_ndvi_by(s_cp_2,r=cp2_seasons)
