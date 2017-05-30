#
# CP1/2 Grassland Productivity Summary Analysis
# This workflow details a pre-cursor analysis to looking at the feasibility of
# using CP1/2 treatments in a region-wide IMBCR analysis
#
# Author : Kyle Taylor (kyle.taylor@pljv.org) [2017]
#

require(raster)
require(rgdal)

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

#
# Grab IMBCR data
#

#
# Read-in treatment vector datasets
#

s_cp_1 <- readOGR("Vector","cp_1", verbose=F)
s_cp_2 <- readOGR("Vector","cp_2", verbose=F)

#
# Calculate NDVI summary statistics
#

cp1_seasons <- sort_ndvi_rasters_by_season(RASTER_DIR, pattern="cp1.*.ndvi.*.tif$")
  cp1_seasons <- calc_seasonal_ndvi_by(s_cp_1,r=cp1_seasons)

cp2_seasons <- sort_ndvi_rasters_by_season(path=RASTER_DIR, pattern="cp2.*.ndvi.*.tif$")
  cp2_seasons <- calc_seasonal_ndvi_by(s_cp_2,r=cp2_seasons)
