#
# Workflow for generating sampling grids for ground-based, mid-winter
# waterfowl surveys
#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

require(rgdal)
require(raster)
require(rgeos)

merge_playas_by_id <- function(s=NULL){
  dissolved_playas <- rgeos::gUnaryUnion(
      s,
      id=s$PPv4_ID
    )
  join_attributes <-
      sp::over(
        rgeos::gBuffer(dissolved_playas, width=2, byid=T),
        s,
        returnList=T
      )
  join_attributes <- do.call(
          rbind,
          # arbitrarily take last playa id's values
          lapply(join_attributes, FUN=function(x) { x[1,] })
        )
  rownames(join_attributes) <- names(dissolved_playas)
  return(
    SpatialPolygonsDataFrame(
      dissolved_playas,
      data=join_attributes
    )
  )
}

calc_basin_station_angles <- function(playa_basins=NULL){
  # calculate the centroid from our SpatialLines features to accomodate a
  # Northing angle for our camera positioning
  # centroids <- rgeos::gCentroid(
  #     playa_basin_samples[(splits[i]+1):splits[i + 1], ],
  #     byid=T
  #   )

  # playa_basins[[length(playa_basins)]]$angle <- apply(
  #   X=playa_basins[[length(playa_basins)]]@coords,
  #   MARGIN=1,
  #   FUN=function(x){
  #       x2 <- centroid@coords[,1]
  #       y2 <- centroid@coords[,2]
  #       x1 <- x[1]
  #       y1 <- x[2]
  #       a <- atan2(y2-y1,x2-x1) * 180 / pi
  #       if(a<0) a <- a+360 # add 2 pi radians if we have a negative angle
  #     }
  #   )

}

generate_playa_upland_pts <- function(playas = NULL,
                                      chunk_size = 1000,
                                      n = 20,
                                      buffer_size=200){
  # buffer an arbitrary dist away from the known
  # playa boundary to keep your feet dry
  playa_upland_samples <- rgeos::gBuffer(
      playas,
      width=buffer_size,
      byid=T
    )
  splits       <- seq(0, nrow(playa_upland_samples), by=chunk_size)
  playa_uplands <- list()
  cat(" -- chunking and processing (this could take awhile): ")
  for(i in 1:(length(splits)-1)){
    playa_uplands[[length(playa_uplands)+1]] <-
      rgeos::gDifference(
          playa_upland_samples[(splits[i]+1):splits[i + 1], ],
          playas[(splits[i]+1):splits[i + 1], ],
          byid = F,
        )
    # generate n random points around each playa
    sample_pts <- raster::rasterToPoints(raster::rasterize(
        playa_uplands[[length(playa_uplands)]],
        raster::raster(extent(playa_uplands[[length(playa_uplands)]]),
        resolution=50
      )),
      sp=T)
    # build a data.frame containing attributes for our generated sample points
    join_attributes <-
      sp::over(
        rgeos::gBuffer(sample_pts, width=buffer_size, byid=T),
        playas[(splits[i]+1):splits[i + 1], ],
        returnList=T
      )
    # merge SpatialPoints against our data.frame
    playa_uplands[[length(playa_uplands)]] <- SpatialPointsDataFrame(
        sample_pts,
        data=do.call(
            rbind,
            # arbitrarily take last playa id's values
            lapply(join_attributes, FUN=function(x) { x[1,] })
          )
      )
    cat("[", i, "/", length(splits)-1, "]", sep = "")
  }; cat("\n")
# bind our and return to the user
return(
    do.call(
        rbind,
        playa_uplands
      )
  )
}

generate_playa_basin_pts <- function(playas=NULL,
                                     chunk_size = 1000,
                                     n = 20,
                                     buffer_size=5){

  # buffer an arbitrary dist away from the known
  # playa boundary to keep your feet dry
  playa_basin_samples <- rgeos::gBuffer(
     focal_playas,
     width=buffer_size,
     byid=T
   )

  splits       <- seq(0, nrow(playa_basin_samples), by=chunk_size)
  playa_basins <- list()

  cat(" -- chunking and processing (this could take awhile): ")
  for(i in 1:(length(splits)-1)){
    playa_basins[[length(playa_basins)+1]] <-
      as(
          playa_basin_samples[(splits[i]+1):splits[i + 1], ],
          'SpatialLinesDataFrame'
        )
    # generate n regular points around each playa
    sample_pts <- sp::spsample(
        playa_basins[[length(playa_basins)]],
        n=chunk_size*n,
        type="regular"
      )
    # build a data.frame containing attributes for our generated sample points
    join_attributes <-
      sp::over(
        rgeos::gBuffer(sample_pts, width=6, byid=T),
        playa_basins[[length(playa_basins)]],
        returnList=T
      )
    # merge SpatialPoints against our data.frame
    playa_basins[[length(playa_basins)]] <- SpatialPointsDataFrame(
        sample_pts,
        data=do.call(
            rbind,
            # arbitrarily take first playa id's values
            lapply(join_attributes, FUN=function(x) { x[1,] })
          )
      )

    cat("[", i, "/", length(splits)-1, "]", sep = "")
  }; cat("\n")
  # bind our and return to the user
  return(
      do.call(
          rbind,
          playa_basins
        )
    )
}

#
# MAIN
#

cat(" -- applying a quantile filter to hydroperiod for project area playas\n")

fall_hydroperiod   <- raster(
    "~/Incoming/playas/fall_years_wet_normalized.tif"
  )

winter_hydroperiod <- raster(
    "~/Incoming/playas/winter_years_wet_normalized.tif"
  )

cat(" -- merging regional playa segments (by id)\n")

focal_playas <- merge_playas_by_id(readOGR(
    dsn="/home/ktaylora/Incoming/playas/",
    layer="playa_cluster_nm_tx_ks_tri_cluster_subset",
    verbose=F
  ))

focal_playas$fall_hydroperiod <- as.vector(raster::extract(
    fall_hydroperiod,
    spTransform(focal_playas,CRS(projection(fall_hydroperiod))),
    fun=mean,
    na.rm=T
  ))

focal_playas$winter_hydroperiod <- as.vector(raster::extract(
    winter_hydroperiod,
    spTransform(focal_playas,CRS(projection(winter_hydroperiod))),
    fun=mean,
    na.rm=T
  ))

keep <-
  focal_playas$fall_hydroperiod >
    quantile(focal_playas$fall_hydroperiod, p=0.6, na.rm=T) &
  focal_playas$winter_hydroperiod >
    quantile(focal_playas$winter_hydroperiod, p=0.6, na.rm=T)

keep[is.na(keep)] <- FALSE # our bone dry playas will have NA values

focal_playas <- focal_playas[keep,]; rm(keep)

writeOGR(
    focal_playas,
    ".",
    paste("playa_cluster_nm_tx_ks_tri_cluster_subset_",
          "joined_by_common_playa_ids", sep = ""),
    driver = "ESRI Shapefile",
    overwrite = T
  )

cat(" -- look for a chunksize that makes sense:\n")

print(data.frame(chunksize=nrow(focal_playas)/1:20, chunks=1:20))

basin_sampling_pts <- generate_playa_basin_pts(
    focal_playas,
    chunk_size=581
  )

writeOGR(
    basin_sampling_pts,
    ".",
    "playa_cluster_nm_tx_ks_tri_cluster_subset_basin_sampling_pts",
    driver="ESRI Shapefile",
    overwrite=T
  )

upland_sampling_pts <- generate_playa_upland_pts(
      focal_playas,
      chunk_size=581
    )
    
writeOGR(
    upland_sampling_pts,
    ".",
    "playa_cluster_nm_tx_ks_tri_cluster_subset_upland_sampling_pts",
    driver="ESRI Shapefile",
    overwrite=T
  )
