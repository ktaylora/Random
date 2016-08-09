require(raster)
require(rgdal)
require(landscapeAnalysis)

processLandfireZip <- function(x){
  zip_name <- unlist(strsplit(x,"/"))
    zip_name <- zip_name[length(zip_name)]
  zip_folders <- unlist(strsplit(zip_name,split="_"))
    zip_folders <- gsub(paste(zip_folders[2:3],collapse="_"),pattern=".zip",replacement="")

  system(paste("unzip -o ",x," -d /tmp",sep=""));

  if(grepl(x,pattern="EVH")){ # is this VEG HEIGHT data?
    cat(" -- calculating vegetation height for associations\n")

    veg_height <- raster(paste("/tmp",toupper(zip_folders),tolower(zip_folders),sep="/"));

    grass_height <- veg_height
      grass_height[grass_height<101] <- 0
      grass_height[grass_height>103] <- 0
      grass_height <- ((2*(grass_height-100)*0.25))-0.25 # units are meters

    shrub_height <- veg_height
      shrub_height[shrub_height<104] <- 0
      shrub_height[shrub_height>107] <- 0
      shrub_height[shrub_height==104] <- 0.25
      shrub_height[shrub_height==105] <- 0.75
      shrub_height[shrub_height==106] <- 2
      shrub_height[shrub_height==107] <- 3.25

    tree_height <- veg_height
      tree_height[tree_height<108] <- 0
      tree_height[tree_height>111] <- 0
      tree_height[tree_height==108] <- 2.5
      tree_height[tree_height==109] <- 7.5
      tree_height[tree_height==110] <- 17.5
      tree_height[tree_height==111] <- 37.5

    return(list(grass_height,shrub_height,tree_height))

  } else if(grepl(x,pattern="EVC")) { # is this veg cover data?
    cat(" -- calculating vegetation % cover for associations\n")

    veg_cover <- raster(paste("/tmp",toupper(zip_folders),tolower(zip_folders),sep="/"));

    grass_perc_cover <- veg_cover
      grass_perc_cover[grass_perc_cover < 121] <- 0
      grass_perc_cover[grass_perc_cover > 129] <- 0
      grass_perc_cover <- ((grass_perc_cover-120)*10)+5

    shrub_perc_cover <- veg_cover
      shrub_perc_cover[shrub_perc_cover < 111] <- 0
      shrub_perc_cover[shrub_perc_cover > 119] <- 0
      shrub_perc_cover <- ((shrub_perc_cover-110)*10)+5


    tree_perc_cover <- veg_cover
      tree_perc_cover[tree_perc_cover < 101] <- 0
      tree_perc_cover[tree_perc_cover > 109] <- 0
      tree_perc_cover <- ((tree_perc_cover-100)*10)+5

    return(list(grass_perc_cover,shrub_perc_cover,tree_perc_cover))
  }
  return(NULL)
}
# merge an Nx3 dimensional list of Landfire raster pieces
merge_by <- function(x, column=NULL, filename=NULL){
  # fetch the i'th item from each list in our list-of-lists
  # e.g., x[1:4,i]
  focal <- lapply(x, FUN=function(x,i=NULL){ return(x[[i]]) }, i=column)
  writeRaster(do.call(raster::merge,focal),filename=filename,overwrite=T)
}

# process explanatory data for our grid points
if(length(list.files(pattern="grass_|shrub_|tree_"))!=6){
  cover_zips <- list.files("Raster/LANDFIRE/broader_regional_landscape/",pattern="EVC",full.names=T)
  height_zips <- list.files("Raster/LANDFIRE/broader_regional_landscape/",pattern="EVH",full.names=T)

  cover <- vector('list',length(cover_zips))
  for(i in 1:length(cover_zips)){
    cover[[i]] <- processLandfireZip(x=cover_zips[i])
  }

  # i.e., processLandfireZip() returns : grass_perc_cover, shrub_perc_cover, tree_perc_cover
  merge_by(cover,column=1,filename="grass_perc_cover.tif")
    merge_by(cover,column=2,filename="shrub_perc_cover.tif")
      merge_by(cover,column=3,filename="tree_perc_cover.tif")

  # clean-up
  rm(cover); gc()

  height <- vector('list',length(height_zips))
  for(i in 1:length(height_zips)){
    height[[i]] <- processLandfireZip(x=height_zips[i])
  }

  merge_by(height,column=1,filename="grass_height.tif")
    merge_by(height,column=2,filename="shrub_height.tif")
      merge_by(height,column=3,filename="tree_height.tif")

  rm(height); gc()

}

# read our explanatory habitat variables and re-grid
# so they are spatially consistent with gridded eBird response data

grid_pts <- list.files(Sys.getenv("HOME"),recursive=T,full.names=T,pattern="spp_grid_points_processed.shp$")[1]

if(file.exists(grid_pts)){
  grid_pts <- landscapeAnalysis:::.readOGRfromPath(grid_pts)
  if(file.exists("grid.tif")){
    grid <- raster("grid.tif")
  } else {
    s <- list.files(Sys.getenv("HOME"),recursive=T,full.names=T,pattern="pts_2002_2012.shp$")[1]
      s <- landscapeAnalysis:::.readOGRfromPath(s)
    grid <- raster(ext=extent(s),nrows=200,ncols=200,crs=CRS(projection("+init=epsg:4326")))
    rm(s);
  }
} else {
  stop("recursive search didn't find a grid_pts_processed.shp to work with")
}

height <- list()
cover  <- list()

height[[1]] <- raster("grass_height.tif")
height[[2]] <- raster("shrub_height.tif")
height[[3]] <- raster("tree_height.tif")

 cover[[1]] <- raster("grass_perc_cover.tif")
 cover[[2]] <- raster("shrub_perc_cover.tif")
 cover[[3]] <- raster("tree_perc_cover.tif")

height <- landscapeAnalysis::snapTo(height,to=grid)
 cover <-  landscapeAnalysis::snapTo(cover,to=grid)
