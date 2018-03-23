reprojectToTarget <- function(x=NULL){
  require(raster)
  return(projectRaster(
      from=raster::raster(x),
      to=raster::raster(paste(
        "/gis_data/Elevation/30m_D3P",
        "_elevation_pljv_region.tif",
        sep=""))
    ))
}

#
# MAIN
#

require(parallel)
require(raster)
require(rgdal)

files <- list.files(pattern="250m.tif$")

# don't try and do files that are already finished
files <- files[
    !( files %in% gsub(files, pattern="250m.tif", replacement="30m.tif") )
  ]

if(length(files)<12){
  splits <- 0:length(files)
} else {
  splits <- seq(0, length(files), by=12)
}


for(i in 1:(length(splits)-1)){
  focal <- files[(splits[i]+1):(splits[i+1])]
  focal <- focal[!
        (gsub(focal, pattern="_250m.tif", replacement="_30m.tif") %in%
          list.files(
            ".",
            recursive=T,
            pattern="_30m.tif$",
            full.names=F
          ))
      ]
  if(length(focal)>0){
    pCl <- parallel::makeCluster(12)
    # reproject
    rasters <- parallel::parLapply(
        cl=pCl,
        X=focal,
        fun=reprojectToTarget
      )
     # save rasters
     mapply(rasters,
         FUN=writeRaster,
         gsub(paste(unlist(lapply(
             rasters,
             names
           )),"_30m.tif",sep=""),
           pattern = "_250m",
           replacement = ""
         ),
         progress='text',
         overwrite=T
       )
     # clean-up
     unlink(list.files(dirname(
            raster::rasterTmpFile()),
            recursive = T,
            full.names = T,
            include.dirs = T
          ),
          recursive = T,
          force = T
       )
     # drop our cluster and free-up swap
     parallel::stopCluster(pCl)
     file.remove(
         dirname(rasterTmpFile()),
         recursive=T
       )
  }
}
