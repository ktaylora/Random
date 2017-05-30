require(raster)
require(rgdal)

calc_seasonal_ndvi_by <- function(s=NULL,r=NULL){
  s = sp::spTransform(s,sp::CRS(raster::projection(r)))
  return(raster::extract(r,s,na.rm=T,fun=mean,progress='text',df=T))
}

cp1_seasons <- raster::stack(list.files(pattern="cp1.*.ndvi.*.tif$"))
  names(cp1_seasons) <- c("fall","spring","summer","winter")
    cp1_seasons <- cp1_seasons[[c(2,3,1,4)]]

cp1_seasons <- calc_seasonal_ndvi_by(s_cp_1,r=cp1_seasons)

cp2_seasons <- raster::stack(list.files(pattern="cp2.*.ndvi.*.tif$"))
  names(cp2_seasons) <- c("fall","spring","summer","winter")
    cp2_seasons <- cp2_seasons[[c(2,3,1,4)]]

cp2_seasons <- calc_seasonal_ndvi_by(s_cp_1,r=cp1_seasons)
