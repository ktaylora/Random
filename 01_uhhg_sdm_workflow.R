#
# WORKFLOW : Big Sagebrush Subspecies and Invasive Brome Climate SDMs
#
# This is a full-run generating content outlined in chapters 2/3 of my thesis. I wrote it
# on a unix machine. You may need to switch "/" with "\" in paths where I slip-up if
# you are running this on windows. Your workstation should have a good amount of ram (30 gigs)
# and about 500 gigs of temp file space (on unix, check with `df -h | grep /tmp`). If
# you don't have that much, place raster::rasterOptions(tmpdir="/home/my_home_dir") that
# points to a partition with enought space somewhere after the includes, below. Email/github
# me if you have any questions.
#
# Author : Kyle Taylor (kyle.taylor@uwyo.edu) [2014] (https://github.com/ktaylora)
#

include <- function(x,from="cran",repo=NULL){
  if(from == "cran"){
    if(!do.call(require,as.list(x))) install.packages(x, repos=c("http://cran.revolutionanalytics.com","http://cran.us.r-project.org"));
    if(!do.call(require,as.list(x))) stop("auto installation of package ",x," failed.\n")
  } else if(from == "github"){
    if(!do.call(require,as.list(x))){
      if(!do.call(require,as.list('devtools'))) install.packages('devtools', repos=c("http://cran.revolutionanalytics.com","http://cran.us.r-project.org"));
      require('devtools');
      install_github(paste(repo,x,sep="/"));
    }
  } else{
    stop(paste("couldn't find package:", x))
  }
}

include("foreign")
include("raster")
include("rgdal")
include("rgeos")
include("randomForest")
include("rfUtilities")
include("rvest")

#' attempts to download and unpack a gbif herbarium dataset, by fixed session ID
gbifFetchAndUnpack <- function(gbifID=NULL){
  if(!file.exists(file.path(paste("herbarium_records/",gbifID,".csv",sep="")))){
    download.file(paste("http://api.gbif.org/v1/occurrence/download/request/",gbifID,".zip",sep=""), destfile=paste(gbifID,".zip",sep="")) # this URL should persist indefinitely
    unzip(paste(gbifID,".zip",sep=""), exdir="herbarium_records")
  }
}
#' converts post-processed herbarium data to spatial points dataframe
gbifDf_toSpatialPoints <- function(t=NULL) {
  sp::SpatialPointsDataFrame(coords = data.frame(x=t$decimallongitude,y=t$decimallatitude),
                             proj4string=sp::CRS("+init=epsg:4326"), # assume WGS84 for all GBIF records
                             data=data.frame(response=rep(1,length(t$decimallatitude)))
                            )
}
#' gbifParseDf splits raw gbif herbarium data into US/Non-US training/hold-out datasets
#' that are later merged into the boundary regions of our GAP dataset (non-us) and used to evaluate
#' later model predictions (US)
gbifParseDf <- function(gbifID=NULL,t=NULL){
  # simultaneous update -- sample boundaries outside of the US, keep the US records for model validation
  t <- read.csv(file.path(paste("herbarium_records/",gbifID,".csv",sep="")),sep="\t")
  t <- t[!is.na(t$decimallatitude),] # kludge : if we can't get a reliable 'decimallatitude', throw it out

  holdout  <- t[t$countrycode == "US",]
  training <- t[t$countrycode != "US",]

  ret <- list(training,holdout)
    names(ret) <- c("training","holdout")
      return(ret)
}
#' wrapper function that processes (gbif) herbarium records for a given scientific name (matched by pattern=)
#' @export
fetchAndProcessGbifRecords <- function(gbifID=NULL, pattern=NULL){
 gbifFetchAndUnpack(gbifID=gbifID)
   a <- gbifParseDf(gbifID=gbifID) # returns list of non-us (training) and us (holdout) records

 a_training <- try(gbifDf_toSpatialPoints(a$training[grepl(tolower(a$training$scientificname),pattern=pattern),]))
 a_holdout  <- try(gbifDf_toSpatialPoints(a$holdout[grepl(tolower(a$holdout$scientificname),pattern=pattern),]))

 ret <- list(a_training,a_holdout)
   names(ret) <- names(a)
     return(ret)
}
#' an iterative and faster implementation of a stratified random sampling
#' algorithm for large rasters
rasterSampleStratified <- function(r=NULL, strata=c(1,0), blockSize=20000, nCell=round(raster::ncell(r)*0.000002),e=1.1){
  require(sp)
  cat(" -- blockwise processing input raster data:")
  blocks <- raster::sampleRandom(r,size=blockSize,sp=T)
  cat(paste("[pos:",sum(blocks$layer==1),";neg:",sum(blocks$layer==0),"]",sep=""))
  while(sum(blocks$layer == 0) < nCell*e | sum(blocks$layer == 1) < nCell*e ){
    block <- raster::sampleRandom(r,size=blockSize,sp=T)
    if(sum(blocks$layer==1)< nCell*e){
      blocks <- rbind(blocks, block[block$layer==1,])
    }
    if(sum(blocks$layer==0)< nCell*e){
      blocks <- rbind(blocks, block[block$layer==0,])
    }
    cat(paste("[pos:",sum(blocks$layer==1),";neg:",sum(blocks$layer==0),"]",sep=""))
  }
  # downsample to our target density
  blocks <- rbind(blocks[sample(which(blocks$layer==1),size=nCell),],
                  blocks[sample(which(blocks$layer==0),size=nCell),])

  names(blocks) <- "response"
  return(blocks)
}
#' wrapper for raster::sampleStratified that estimates a sample-size to draw
#' random samples from a gap raster before converting to stratified random
#' points
gapRaster_toPts <- function(r, size=NULL) {
  if(is.null(size)){
    size = round(raster::ncell(r)*0.000002)
  }
  return(rasterSampleStratified(r, blockSize=20000, nCell=size))
  # sloooow :
  # r <- raster::sampleStratified(r, size=size,
  #                               strata=c(1,0),
  #                               sp=T,
  #                               exp=5)
  # r@data <- data.frame(response=r@data$layer)
}

#' grep-by pattern
grep_by <- function(x, pattern=NULL){
  x[grepl(as.character(x), pattern = pattern)]
}
#' grab all zips by default and parse through what the user requested
parse_hrefs <- function(url=NULL){
  a <- rvest::html(url)
    return(as.character(rvest::html_nodes(a,"a")))
}
#' grab the full string within a set of parantheses
parse_url <- function(x=NULL){
  unlist(lapply(strsplit(x, split = "\""),
                FUN = function(x) x[which(grepl(x, pattern = "zip$"))]
                ))
}
#' scrape EDDMaps counties for latest bromus tectorum data
scrapeEddmapsData <- function(base_url="https://goo.gl/375Nb7"){
  # hack ...
  return("http://bigbeatbox.duckdns.org/b_tectorum_eddmaps.zip")
}
fetchAndProcessEddmap <- function(url){
  download.file(url, destfil=file.path("boundaries/b_tectorum_edds_counties.zip"))
  unzip("boundaries/b_tectorum_edds_counties.zip",exdir="boundaries/eddmaps")
}
#'
bTectorumEddmap <- function(){
  if(!dir.exists("boundaries/eddmaps")){
    fetchAndProcessEddmap(scrapeEddmapsData())
  }
  if(!dir.exists("boundaries/counties")){
    dir.create("boundaries/counties")
    counties <- scrapeUsBoundaries(pattern="us_county_500k")
      download.file(counties,destfile=file.path("boundaries/counties/counties.zip"))
        unzip(file.path("boundaries/counties/counties.zip"),exdir="boundaries/counties")
  }
  # point in polygon
  counties  <- readOGR("boundaries/counties/","cb_2015_us_county_500k",verbose=F)
    cg_points <- sp::spTransform(rgdal::readOGR("boundaries/eddmaps", "points",verbose=F), CRS(projection(counties)))
      cg_points <- rgeos::gBuffer(cg_points,width=0.6,byid=T) #we need to buffer to re-create the map
  over <- sp::over(counties,cg_points)[,1]
  # punt
  counties@data <- data.frame(response=!is.na(over))
  return(counties)
}
#' scrape TIGER for us boundary
scrapeUsBoundaries <- function(base_url="http://www2.census.gov/geo/tiger/GENZ2015/shp/", pattern="us_state_500k"){

  hrefs <- parse_hrefs(url = base_url)
  hrefs <- grep_by(hrefs,pattern=pattern)
  hrefs <- unlist(lapply(hrefs, FUN = parse_url))[1]

  return(paste(base_url,hrefs,sep=""))
}
#' scrape Robert Heijman's worldclim data
scrapeWorldclimCmip5 <- function(base_url="http://worldclim.org/cmip5_30s",
                                 models=NULL,            # vector of two-letter codes for corresponding model (see: website)
                                 climate=FALSE,
                                 bioclim=FALSE,
                                 start_of_century=FALSE,
                                 mid_century=FALSE,
                                 end_of_century=FALSE,
                                 scen_45=FALSE,
                                 scen_85=FALSE
                                 ){
  # fetch and parse index.html
  hrefs <- parse_hrefs(url = base_url)
  # grep hack through user parameters
  time_periods <- vector()
  climate_sets <- vector()
  scenarios    <- vector()
  # grep across time aggregations (handle "current" last, it is it's own sep query)
  if (mid_century){
    time_periods <- append(time_periods,"50[.]zip")
  }
  if (end_of_century){
    time_periods <- append(time_periods,"70[.]zip")
  }
  hrefs <- grep_by(hrefs, pattern = paste(time_periods, collapse = "|") )
  # grep across raw climate / bioclim
  if (bioclim){
    climate_sets <- append(climate_sets, "bi.[0-9]")
  }
  if (climate){
    # not implemented -- do some digging
    #hrefs <- grep_by(hrefs, pattern = ".bi..[.]zip")
  }
  hrefs <- grep_by(hrefs, pattern = paste(climate_sets, collapse = "|") )
  # grep across emissions scenarios
  if (scen_45){
    scenarios <- append(scenarios,"/..45")
  }
  if (scen_85){
    scenarios <- append(scenarios,"/..85")
  }
  hrefs <- grep_by(hrefs, pattern = paste(scenarios, collapse = "|") )
  # grep across GCMs
  if(!is.null(models[1])){
    hrefs <- grep_by(hrefs, pattern = paste(paste("/",models,sep=""), collapse = "[0-9]|") )
  }
  # append "current" climate if it was requested
  if (start_of_century){
    # current climate data are stashed on another webpage...
    hrefs_current <- parse_hrefs("http://worldclim.org/current")
    hrefs_current <- grep_by(hrefs_current, pattern = "cur" )
    hrefs_current <- grep_by(hrefs_current, pattern = "30s" )
    hrefs_current <- grep_by(hrefs_current, pattern = "bil.zip" )
    # bioclim or original climate variables?
    if(bioclim){
      hrefs_current <- grep_by(hrefs_current, pattern = "/bi.[0-9]" )
    }
    time_periods <- append(time_periods,"cur")
    hrefs <- append(hrefs, hrefs_current)
  }
  # return URLs to user
  return(unlist(lapply(hrefs, FUN = parse_url)))
}
#' unpack climate data
worldclimFetchAndUnpack <- function(urls=NULL){
  climate_zips <- unlist(lapply(strsplit(urls,split="/"),FUN=function(x){return(x[length(x)])}))
  if( sum( climate_zips %in% list.files("climate_data") ) != length(climate_zips) ){
    cat(" -- fetching missing climate data\n")
    missing <- which(!(climate_zips %in% list.files("climate_data")))
    for(i in missing){
      download.file(urls[i],destfile=file.path(paste("climate_data",climate_zips[i],sep="/")))
    }
  }
  if(sum( list.files( file.path(paste("climate_data/",gsub(climate_zips, pattern=".zip",replacement="*.tif"))) ))){
    cat(" -- unpacking climate data\n")
    for(i in climate_zips){
      unzip(file.path(paste("climate_data",i,sep="/")),exdir="climate_data")
    }
  }
}
#' return a RasterStack of worldclim data, parsed by flags and/or pattern
worldclimStack <- function(vars=NULL, time=NULL, bioclim=T, climate=F, pattern=NULL){
  all_vars <- list.files("climate_data",pattern="bil$|tif$", full.names=T)
  climate_flags <- vector()
  if(grepl(tolower(time), pattern="cur")){
    all_vars <- grep_by(all_vars, pattern="bil$") # current climate data are always .bil files
  } else if(grepl(as.character(time),pattern="50")){
    all_vars <- grep_by(all_vars, pattern=".50.")
  } else if(grepl(as.character(time),pattern="70")){
    all_vars <- grep_by(all_vars, pattern=".70.")
  }
  if(bioclim){
    climate_flags <- append(climate_flags,"bi")
  }
  if(climate){
    climate_flags <- append(climate_flags,"temp|prec")
  }

  all_vars <- grep_by(all_vars, pattern=paste(climate_flags,collapse="|"))

  if(!is.null(vars)){
    # parse out weird worldclim var formating (differs by time scenario)
    if(grepl(tolower(time), pattern="cur")){
      v <- unlist(lapply(strsplit(all_vars,split="/"), FUN=function(x){x[length(x)]}))
      v <- lapply(strsplit(v,split="_"),FUN=function(x){x[length(x)]})
      v <- as.numeric(gsub(v, pattern=".bil", replacement=""))
        all_vars <- all_vars[ v %in% vars ]
    } else if(grepl(tolower(time), pattern="50|70")){
      split <- gsub(as.character(time), pattern="20", replacement="")
      v <- unlist(lapply(strsplit(all_vars,split="/"), FUN=function(x){x[length(x)]}))
      v <- lapply(strsplit(v,split=split),FUN=function(x){x[length(x)]})
      v <- as.numeric(gsub(v, pattern=".tif", replacement=""))
        all_vars <- all_vars[ v %in% vars ]
    }
  }
  # any final custom grep strings passed?
  if(!is.null(pattern)){
    all_vars <- grep_by(all_vars, pattern=pattern)
  }

  raster::stack(all_vars)
}
#' will reclassify and downsample a gap dataset to a reasonable density (arbitrary total landscape * 0.002%)
#' and make a randomized stratification vector of community-level presence/absence for a focal
#' species.
parse_gap_training_data <- function(r=NULL, gap=NULL, gap_rat=NULL){
  r <- gap_rat[which(as.character(gap_rat[,"ECOLSYS_LU"]) == r), ]$Value
  return(gapRaster_toPts(gap == r))
}
#' mask out records that don't occur within a given boundary (e.g., US + Omernik)
mask_by_boundary <- function(s=NULL, region_boundary=NULL){
  over <- sp::over(s, sp::spTransform(region_boundary,CRS(projection(s))))
  if(is.null(ncol(over))){
    return(s[!is.na(over),])
  } else {
    return(s[!is.na(over[,1]),])
  }
}
#' downsample an over-abundant class to the density of the smaller class. This
#' only supports presence-absence (1/0).
balanceClasses <- function(s){
  if ( sum(s$response == 1) > sum(s$response == 0) ){
    big_class <- 1
    small_class <- 0
  } else if ( sum(s$response == 0) > sum(s$response == 1) ){
    big_class <- 0
    small_class <- 1
  } else {
    return(s)
  }
  # downsample to the density of our "small class"
  pos <- sample(which(s$response == big_class), size=sum(s$response == small_class))
    pos <- s[pos,]
  return(rbind(pos,s[s$response == small_class,]))
}
#' generate random pseudo-absences (i.e, in Canada and Mexico) at a density consistent
#' with the area of the western US (US + Omernik) boundary. Default metohds fail
#' for large ,irregular polygons like ours.
#' @export
generate_pseudoabsences <- function(s=NULL, src_region_boundary=NULL, target_region_boundary=NULL){
  require(sp)
  s <- sp::spTransform(s,sp::CRS(raster::projection(target_region_boundary)))
  # point in polygon extent operation
  randPts <- function(s=NULL,n=NULL){
    extent <- raster::extent(s)
    proj   <- sp::CRS(raster::projection(s))
    x_pts  <- runif(min = extent@xmin, max=extent@xmax, n = n)
    y_pts  <- runif(min = extent@ymin, max=extent@ymax, n = n)
    SpatialPointsDataFrame(coords = data.frame(x=x_pts,y=y_pts),
                           data = data.frame(response=rep(0,n)),
                           proj4string = proj
                          )
  }
  cat(" -- generating random pseudo-absences for Canada/Mexico:")
  pts_to_generate <- sum(s$response == 1) # we want our # absences == # presences to satisfy the balanced-class requirements of RF
    pts_to_generate <- pts_to_generate - sum(s$response == 0) # account for the absences we already have from US GAP

  if(pts_to_generate==0) stop("it looks like you already have balanced 1/0's in your training set. quitting.")

  # total area (km2) of our source region (e.g., US)
  src_area_km <- rgeos::gArea(sp::spTransform(src_region_boundary,
                              sp::CRS(raster::projection("+init=epsg:2163")))) * 10^-6
  # total area (km2) of our target region (e.g., non-US)
  target_area_km <- rgeos::gArea(sp::spTransform(target_region_boundary,
                                 sp::CRS(raster::projection("+init=epsg:2163")))) * 10^-6
  # re-scale our pts_to_generate number so that we sample a density consistent
  # with our GAP sampling effort
  target_pts_upsampled <- round((sum(s$response == 0)/src_area_km)*(target_area_km))
  pts_generated <- randPts(target_region_boundary,n=target_pts_upsampled)
  pts_generated <- pts_generated[is.na(sp::over(pts_generated, target_region_boundary)),]
  cat(".");
  while(nrow(pts_generated) < target_pts_upsampled){
    focal <- randPts(target_region_boundary, n=target_pts_upsampled)
    pts_generated <- rbind(pts_generated, focal)
    pts_generated <- pts_generated[!is.na(sp::over(pts_generated,target_region_boundary)),]
    cat(".")
  };
  cat("\n");
  # downsample back to our target number, merge, and return to user
  if(target_pts_upsampled < pts_to_generate){
    warning("To keep consistent sampling density we had to under-sample in our target region. Check for class imbalances.")
  } else {
    pts_generated <- pts_generated[sample(1:nrow(pts_generated),size=pts_to_generate),]
  }

  return(rbind(s,pts_generated))
}
#' generate in/outsample (rows) from a data.frame, t=
genInOutSample <- function(t=NULL,size=0.3){
  ret <- vector("list",2)
  ret[[2]] <- sample(1:nrow(t),size=nrow(t)*size)
  ret[[1]] <- which(!(1:nrow(t) %in% ret[[2]]))
  names(ret) <- c("out_sample","in_sample")
  return(ret)
}
#' build a training dataset for our random forest model using a combination
#' of GAP and herbarium (gbif) records. This function is a kludge, lo siento.
buildTrainingDataset <- function(spp="spp_tridentata", gap=NULL,
                                 gap_rat=NULL, ecolsys_lu=NULL,
                                 gbif_df=NULL,
                                 region_boundary_us=NULL,
                                 region_boundary_non_us=NULL){

  cat(paste(" -- ",spp,"...\n",sep=""))
  s <- parse_gap_training_data(ecolsys_lu, gap=gap, gap_rat=gap_rat)
    writeOGR(s,"training_data",paste(spp,"_gap_source_data",sep=""),driver="ESRI Shapefile", overwrite=T)

  s <- mask_by_boundary(s, region_boundary=region_boundary_us) # (e.g., US + Omernik)

  # merge-in any GBIF training records that are available
  if(class(gbif_df$training)!="try-error"){
    ssp_tridentata <- rbind(sp::spTransform(gbif_df$training,
                            CRS(projection(s))), s)
  }

  # generate pseudo-absences in areas unsampled by GAP (Canada + Mexico + Omernik)
  s <- generate_pseudoabsences(s,
                               src_region_boundary=region_boundary_us,
                               target_region_boundary=region_boundary_non_us)
  # dupe check and balance classes
  s <- balanceClasses(s[!duplicated(s@coords),])

  # write to disk for validation
  writeOGR(s,"training_data",paste(spp,"_gap_training_data",sep=""),driver="ESRI Shapefile",overwrite=T)

  return(s)
}
#' accepts spatial point data and re-grids observations to the scale of raster r=
regridTo <- function(s=NULL, r=NULL){
  cat(" -- re-gridding:\n")
  s <- sp::spTransform(s,sp::CRS(raster::projection(r)))
  r <- raster::rasterize(x=s,y=r,field='response',progress='text')
  s <- raster::rasterToPoints(r,spatial=T,progress='text')
  names(s) <- "response"
  return(s)
}
#' accepts a Random Forest model object, reports it's internal (OOB) validation statistics
#' and then evaluates using an external validation dataset of your choosing
reportValidation <- function(m,t){
  cat(" -- confusion matrix (OOB):\n")
  print(m$confusion)
  if(length(t$response)/length(m$y) < 0.2 ){
    cat(" -- warning holdout dataset is less-than 20% of the size of the training dataset\n")
  }
  cat(" -- accuracy vs holdout : ",sum(predict(m,t) == t$response,na.rm=T)/length(t$response),"\n")
}
#' fit a Random Forest and validate using herbarium presences and GAP holdout
fitAndValidate <- function(main=NULL, s=NULL, s_holdout=NULL){
  if(!is.null(main)) cat(paste(" -- model fitting for ",main,":\n",sep=""))
  sample <- genInOutSample(t=s@data)
  m <- randomForest(as.factor(response)~.,data=na.omit(s@data[sample$in_sample,]),do.trace=F)
    reportValidation(m,s_holdout@data)
    reportValidation(m,s@data[sample$out_sample,])
  return(m)
}
#' contour feature extraction for raster predictions that gives a broad (~100 km) probability density
#' for exploratory analysis
featureExtractionContour <- function(r=NULL){
  suitability_contour <- focal(r,fun=mean,w=matrix(1,nrow=3,ncol=3))
  suitability_contour <- focal(suitability_contour,fun=mean,w=matrix(1,nrow=11,ncol=11))
  suitability_contour <- focal(suitability_contour,fun=mean,w=matrix(1,nrow=33,ncol=33))
  suitability_contour <- focal(suitability_contour,fun=mean,w=matrix(1,nrow=99,ncol=99))
    raster::rasterToContour(suitability_contour)
}
#' re-scale (normalize) proportional votes to "probability" based on predictions made across the entire region. Usually
#' not needed (i.e., when our proportional votes have a min/max of 0/1), but sometimes it is (see also: 'rfUtilities'). We
#' optionally take the inverse of the current prediction, because raster::predict doesn't support selecting focal class
#' when working with raster objects, which means that all of our model predictions will default to probability of absence, rather than probability
#' of occurrence.
scale.dist <- function(d, invert=T) {
  if(invert){
    d <- cellStats(d,"max",na.rm=T) - d
  }
  return( ( d - cellStats(d,"min",na.rm=T) )/( cellStats(d,"max",na.rm=T) - cellStats(d,"min",na.rm=T) ) )
}
#' rename raster stack vars from e.g. "ipbio504" to the "bio_4" strings models "see" under current climate
#' conditions worldclim rasters
namesToBioclim <- function(r, prefix="bio_", split_pattern="bio50"){
  n <- paste(prefix,unlist(lapply(strsplit(names(r),split=split_pattern),FUN=function(x){ return(x[2]) })),sep="")
  names(r) <- n
    return(r)
}
#' use grep splitting to slice and project across a worldclim raster stack
projectModelAcrossGCMs <- function(m=NULL,r=NULL, split_pattern="bi50",rcp="45"){
  gcms <- vector("list",length(gcm_abbrevs))
    names(gcms) <- gcm_abbrevs
  cat(paste(" -- projecting across GCM's [",paste(gcm_abbrevs,collapse=","),"]\n",sep=""))
  for(i in 1:length(gcms)){
    focal <- subset(r,names(r)[grepl(names(r),pattern=gcm_abbrevs[i])])
    focal <- subset(focal,names(focal)[grepl(names(focal),pattern=rcp)])
    focal <- namesToBioclim(focal, split_pattern=split_pattern)
    gcms[[i]] <- scale.dist(predict(focal, m, type="prob", progress='text'))
  }
  return(gcms)
}
#' this is getting a little hackish. returns a list-of-lists
#' with aggregate statistics calculated across our predictive surfaces
projectModelByYearAndScenario <- function(spp=NULL,m=NULL,year=c("50","70"),rcps=c("45","85")){
  ret <- vector("list", (2*length(year)*length(rcps)))
  j <- 1
  for(y in year){
    for(pathway in rcps){
      mean_sd <- vector("list",2);
        names(mean_sd) <- c("mean","sd")
      i  <- projectModelAcrossGCMs(m=m, r=get(paste("climate_conditions_20",y,sep="")), split_pattern=y, rcp=pathway)
      ret[[j]]   <- stackApply(stack(i),fun=mean,na.rm=T,indices=1)
        names(ret[[j]]) <- paste(y,pathway,"mean",sep="_")
      ret[[j+1]] <- stackApply(stack(i),fun=sd,na.rm=T,indices=1)
        names(ret[[j+1]]) <- paste(y,pathway,"sd",sep="_")
      j <- j+2
    }
  }
  return(ret)
}
#' wrapper for writeRaster that can comprehend our raster hodge-podge
saveRasters <- function(prefix=NULL, x=NULL){
  cat(" -- saving prediction rasters:\n")
  dir.create("rasters")
  names <- gsub(unlist(lapply(x,FUN=names)),pattern="X",replacement="20")
    names <- gsub(names,pattern="0_",replacement="0_rcp")
      names <- file.path(paste("rasters/",prefix,"_",names,".tif",sep=""))
  mapply(x,FUN=raster::writeRaster,filename=names,progress='text',overwrite=T)
}

################################################################################
##
#   MAIN
##
################################################################################

#
# fetch and read-in omernik ecoregions, national boundaries, and county-level bromus tectorum
# data to 1.) mask relavent GAP observations and 2.) define our absense space (minding the naughty-naughts)
# [https://www.jstor.org/stable/3683828]
#

# pre-define climate and topographic conditions thought to be limiting to big sagebrush,
# but still inclusive of our occurrence data (SEMIARID HIGHLANDS)

raster::rasterOptions(tmpdir="/home/ktaylora/raster_tmp")

study_region_levels <- c("NORTH AMERICAN DESERTS",
                         "MEDITERRANEAN CALIFORNIA",
                         "SOUTHERN SEMIARID HIGHLANDS",
                         "TEMPERATE SIERRAS",
                         "NORTHWESTERN FORESTED MOUNTAINS",
                         "MARINE WEST COAST FOREST",
                         "GREAT PLAINS")

cat(" -- fetching / reading omernik ecoregions data for absence record generation and study region delineation\n")
if(length(list.files("boundaries", pattern="NA_CEC_Eco_Level1"))==0){
  download.file("ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/na_cec_eco_l1.zip",destfile="na_cec_eco_l1.zip")
    unzip("na_cec_eco_l1.zip",exdir="boundaries")
}

region_boundary <- rgdal::readOGR("boundaries","NA_CEC_Eco_Level1",verbose=F)
  region_boundary <- region_boundary[region_boundary@data[,"NA_L1NAME"] %in% study_region_levels,]
    region_boundary <- rgeos::gUnionCascaded(region_boundary)

if(length(list.files("boundaries", pattern="us_state_500k"))==0){
  url <- scrapeUsBoundaries()
    download.file(url,destfile="cb_2015_us_state_500k.zip")
      unzip("cb_2015_us_state_500k.zip",exdir="boundaries")
}
# for ensuring GAP raster observations are within the US
us_boundary <- rgdal::readOGR("boundaries","cb_2015_us_state_500k",verbose=F)
  us_boundary <- us_boundary[!tolower(us_boundary$NAME) %in% c("alaska","puerto rico","hawaii"),]
    us_boundary <- rgeos::gUnionCascaded(us_boundary)

# for masking GAP records to extent of ~ western US
region_boundary_us <- rgeos::gIntersection(us_boundary,
                                           sp::spTransform(region_boundary,
                                                           sp::CRS(raster::projection(us_boundary))
                                                          ))
# for generating non-us pseudo-absences at ecological boundaries (~mexico/canada)
region_boundary_non_us <- rgeos::gDifference(region_boundary, sp::spTransform(us_boundary, CRS(projection(region_boundary)))) # mask out our gap sample space (us boundary)
full_study_region      <- rgeos::gUnion(region_boundary_us, sp::spTransform(region_boundary, CRS(projection(region_boundary_us))))
  writeOGR(SpatialPolygonsDataFrame(full_study_region, data=data.frame(id=1)),
                                    "boundaries", "full_study_region_extent",
                                    driver="ESRI Shapefile", overwrite=T)

#
# fetch our explanatory climate data
#

cat(" -- fetching / cropping / reading worldclim data\n")

# GCM's selected to maximise the range of predictions for temperature/precipitation
# made across all available scenarios/models for the inter-mountain region.
# See VisTrails-SAHM documentation for an overview.

gcm_abbrevs  <- c("ac","gf","in","ip")
bioclim_vars <- c(3,4,11,15,18) # variables selected from the optimal variable set outlined in Schlaepfer et al., 2012. Not interested in variable comparison here
time_periods <- c("current","2050","2070")

dir.create("climate_data")

urls <- scrapeWorldclimCmip5(models=gcm_abbrevs,
                             bioclim=T,
                             start_of_century=sum(grepl(time_periods,pattern="current")),
                             mid_century=sum(grepl(time_periods,pattern="50")),
                             end_of_century=sum(grepl(time_periods,pattern="70")),
                             scen_45=TRUE,
                             scen_85=TRUE
                            )

worldclimFetchAndUnpack(urls)

climate_conditions_current <- worldclimStack(vars=bioclim_vars, bioclim=T, time="current")
  climate_conditions_current <- raster::crop(climate_conditions_current,sp::spTransform(full_study_region,CRS(projection(climate_conditions_current))))
climate_conditions_2050    <- worldclimStack(vars=bioclim_vars, bioclim=T, time=2050, pattern=paste(gcm_abbrevs,collapse="|"))
  climate_conditions_2050 <- raster::crop(climate_conditions_2050,sp::spTransform(full_study_region,CRS(projection(climate_conditions_2050))))
climate_conditions_2070    <- worldclimStack(vars=bioclim_vars, bioclim=T, time=2070, pattern=paste(gcm_abbrevs,collapse="|"))
  climate_conditions_2070 <- raster::crop(climate_conditions_2070,sp::spTransform(full_study_region,CRS(projection(climate_conditions_2070))))

#
# Read-in gbif herbarium records
# GAP lacks coverage for Canada / Mexico, but there are important boundary
# conditions worth sampling in both geographies. We are going to fill-in with
# presence records with georeferenced herbarium records publised via GBIF
# and keep a larger presence holdout dataset for evaluation
#
# See query info here : https://goo.gl/4nF37C
#

ssp_tridentata_gbif   <- fetchAndProcessGbifRecords(gbifID="0075776-160910150852091", pattern=".*nutt*.|.*subsp.*.tridentata.*")
ssp_wyomingensis_gbif <- fetchAndProcessGbifRecords(gbifID="0075776-160910150852091", pattern=".*subsp.*.wyomingensis.*")
ssp_vaseyana_gbif     <- fetchAndProcessGbifRecords(gbifID="0075776-160910150852091", pattern=".*subsp.*.vaseyana.*")

#
# Let's do the same for our b. tectorum records
#
# See query info here : https://goo.gl/HX0169
#

b_tectorum_gbif   <- fetchAndProcessGbifRecords(gbifID="0075763-160910150852091", pattern=".") # this is a complete dataset -- no pattern needed

#
# fetch and read-in GAP raster training data
#
cat(" -- fetching / reading national GAP data\n")

if(length(list.files("national_gap", pattern="natgaplandcov_v2_2_1"))==0){
  dir.create("national_gap")
    download.file("https://usgs-gap-data.s3.amazonaws.com/NAT_LC/Nat_GAP_LandCover.zip",destfile="Nat_GAP_LandCover.zip")
      unzip("Nat_GAP_LandCover.zip",exdir="national_gap")
} else if(length(list.files("national_gap",pattern="natgaplandcov_v2_2_1"))==1){
  unzip("Nat_GAP_LandCover.zip",exdir="national_gap")
}

gap <- raster::raster(file.path("national_gap/natgaplandcov_v2_2_1.img"))
gap_rat <- foreign::read.dbf(file.path("national_gap/natgaplandcov_v2_2_1.img.vat.dbf"))

cat(" -- ensuring GAP is cropped to the extent of the united states -- this may take a bit\n")
gap <- raster::crop(gap, sp::spTransform(us_boundary,raster::projection(gap)))


#
# The 'ECOLSYS_LU' field allows us to crosswalk with big sagebrush subspecies.
# The folks at GAP note that it is difficult distinguising wyomingensis/tridentata
# subspecies. But there are important community-level differences noted in GAP
# See LANDFIRE: https://www.landfire.gov/documents/LF-GAPMapUnitDescriptions.pdf
# search by:
#
# 'Dominants : Artemisia tridentata'
# 'dominated by Artemisia tridentata ssp. wyomingensis'
#
# On Data : there are positives and negatives associated with using GAP as species occurrence data vs.
# herbarium records. On the positive side, GAP gives us a systematic representation
# of species presence-absence across the US. Herbarium records only offer presence information, forcing
# us to fake absences. The GAP dataset is large (hundreds of thousands of records),
# while herbarium records are limited (300-800 records) and often of questionable spatial,
# temporal, and thematic accuracy. On the negative side, GAP represents vegetation as a
# community; not as individual species. Artemesia subspeices presence-absence in GAP is confounded by
# biotic interactions not represented in model space and we have to ignore these effects in our work.
# But that being said, herbarium records collected in the field are confounded by unmeasured biotic interactions
# as well. GAP at-least represents what community a subspecies is in, while herbarium
# records often do not record the plant associations that a species was observed in,
# leaving us guessing as to what interactions in a broader plant community look like. It is for this
# reason that I have chosen to primarily fit GAP data for this modeling effort, while only using
# herbarium records to augment data at regional boundaries (Canada/Mexico) and as a validation
# dataset to test predictive accuracy. I assume that if our GAP trained model is a poor predictor
# of subspecies occurrence/absence, that it will evaluate poorly when validated against our herbarium
# record hold-out data. In the future, we might consider record weighting (e.g., Random Forest), or
# model averaging / data augmentation (such as with HB) across both datasets. The current approach
# is intuitive and will work well for now.
#

cat(" -- performing a binary reclassification of GAP explanatory data for our (sub)species of interest\n")
if(sum(grepl(names(gap_rat),pattern="ECOLSYS_LU"))==0) stop("couldn't find ECOLSYS_LU in GAP attribute table. This shouldn't happen.")

cat(" -- building training datasets for our species of interest:\n")
ssp_tridentata <- buildTrainingDataset(spp="spp_tridentata", gap=gap,
                                       gap_rat=gap_rat, ecolsys_lu="Inter-Mountain Basins Big Sagebrush Shrubland",
                                       gbif_df=ssp_tridentata_gbif,
                                       region_boundary_us=region_boundary_us,
                                       region_boundary_non_us=region_boundary_non_us)

ssp_tridentata <- regridTo(ssp_tridentata, climate_conditions_current)

ssp_wyomingensis <- buildTrainingDataset(spp="ssp_wyomingensis", gap=gap,
                                         gap_rat=gap_rat, ecolsys_lu="Inter-Mountain Basins Big Sagebrush Steppe",
                                         gbif_df=ssp_wyomingensis_gbif,
                                         region_boundary_us=region_boundary_us,
                                         region_boundary_non_us=region_boundary_non_us)

ssp_wyomingensis <- regridTo(ssp_wyomingensis, climate_conditions_current)

ssp_vaseyana <- buildTrainingDataset(spp="ssp_vaseyana", gap=gap,
                                     gap_rat=gap_rat, ecolsys_lu="Inter-Mountain Basins Montane Sagebrush Steppe",
                                     gbif_df=ssp_vaseyana_gbif,
                                     region_boundary_us=region_boundary_us,
                                     region_boundary_non_us=region_boundary_non_us)

ssp_vaseyana <- regridTo(ssp_vaseyana, climate_conditions_current)

# bromus tectorum requires some censoring to mask out observations of annual bromes that are not cheatgrass (e.g., red brome)
# to do this, we will use the most recent EDDMaps county dataset that is available. We will roll our own, here. This will look messy.

cat(" -- b. tectorum...\n")
invaded_counties <- bTectorumEddmap()
  invaded_counties <- invaded_counties[invaded_counties$response ==1,]

bromus_tectorum <- parse_gap_training_data("Introduced Upland Vegetation - Annual Grassland", gap=gap, gap_rat=gap_rat)
  writeOGR(bromus_tectorum,"training_data","bromus_tectorum_gap_source_data",driver="ESRI Shapefile", overwrite=T)

presences <- bromus_tectorum[bromus_tectorum$response==1,]
  presences <- mask_by_boundary(presences, region_boundary=invaded_counties) # mask GAP presence locations containing non-cg invasive annuals
bromus_tectorum <- rbind(presences,bromus_tectorum[bromus_tectorum$response==0,])
  rm(presences)

# merge-in any GBIF training records that are available
if(class(b_tectorum_gbif$training)!="try-error"){
  bromus_tectorum <- rbind(sp::spTransform(b_tectorum_gbif$training,
                          CRS(projection(bromus_tectorum))), bromus_tectorum)
}

# mask out observations made outside of western north america -- e.g., south america (from gbif).
# we are interested in capturing it's climate niche in the invaded inter-mountain region

over <- !is.na(sp::over(bromus_tectorum,full_study_region))
  bromus_tectorum <- bromus_tectorum[over,]

# generate pseudo-absences in areas unsampled by GAP (Canada + Mexico + Omernik)
bromus_tectorum <- generate_pseudoabsences(bromus_tectorum,
                                            src_region_boundary=region_boundary_us,
                                            target_region_boundary=region_boundary_non_us)

# re-grid to scale of our explanatory climate data
bromus_tectorum <- regridTo(bromus_tectorum, climate_conditions_current)

# balance our classes
bromus_tectorum <- balanceClasses(bromus_tectorum)

# write to disk
writeOGR(bromus_tectorum,"training_data","bromus_tectorum_gap_training_data",driver="ESRI Shapefile", overwrite=T)

# extract our climate data
        ssp_tridentata <- raster::extract(climate_conditions_current,ssp_tridentata,sp=T)
ssp_tridentata_holdout <- raster::extract(climate_conditions_current,ssp_tridentata_gbif$holdout,sp=T)

        ssp_wyomingensis <- raster::extract(climate_conditions_current,ssp_wyomingensis,sp=T)
ssp_wyomingensis_holdout <- raster::extract(climate_conditions_current,ssp_wyomingensis_gbif$holdout,sp=T)

        ssp_vaseyana <- raster::extract(climate_conditions_current,ssp_vaseyana,sp=T)
ssp_vaseyana_holdout <- raster::extract(climate_conditions_current,ssp_vaseyana_gbif$holdout,sp=T)

        bromus_tectorum <- raster::extract(climate_conditions_current,bromus_tectorum,sp=T)
bromus_tectorum_holdout <- raster::extract(climate_conditions_current,b_tectorum_gbif$holdout,sp=T)

# fit our Random Forest models
m_tridentata      <- fitAndValidate(main="ssp. tridentata", s=ssp_tridentata, s_holdout=ssp_tridentata_holdout)
m_wyomingensis    <- fitAndValidate(main="ssp. wyomingensis", s=ssp_wyomingensis, s_holdout=ssp_wyomingensis_holdout)
m_vaseyana        <- fitAndValidate(main="ssp. vaseyana", s=ssp_vaseyana, s_holdout=ssp_vaseyana_holdout)
m_bromus_tectorum <- fitAndValidate(main="bromus tectorum")

# project out as rasters
dir.create("rasters")
cat(" -- projecting models across current climate conditions:\n")
ssp_tridentata_suitability_current   <- scale.dist(predict(climate_conditions_current, m_tridentata, type="prob", progress='text'))
  raster::writeRaster(ssp_tridentata_suitability_current,filename=file.path("rasters/ssp_tridentata_current.tif"),progress='text', overwrite=T)
ssp_wyomingensis_suitability_current <- scale.dist(predict(climate_conditions_current, m_wyomingensis, type="prob", progress='text'))
  raster::writeRaster(ssp_wyomingensis_suitability_current,filename=file.path("rasters/ssp_wyomingensis_current.tif"),progress='text', overwrite=T)
ssp_vaseyana_suitability_current     <- scale.dist(predict(climate_conditions_current, m_vaseyana, type="prob", progress='text'))
  raster::writeRaster(ssp_vaseyana_suitability_current,filename=file.path("rasters/ssp_vaseyana_current.tif"),progress='text', overwrite=T)
bromus_tectorum_suitability_current  <- scale.dist(predict(climate_conditions_current, m_bromus_tectorum, type="prob", progress='text'))
  raster::writeRaster(bromus_tectorum_suitability_current,filename=file.path("rasters/bromus_tectorum_current.tif"),progress='text', overwrite=T)

# project out under future climate conditions
cat(" -- projecting models across future climate conditions:\n")
ssp_tridentata_suitability_future    <- projectModelByYearAndScenario(spp="ssp_tridentata", m=m_tridentata)
ssp_wyomingensis_suitability_future  <- projectModelByYearAndScenario(spp="ssp_wyomingensis", m=m_wyomingensis)
ssp_vaseyana_suitability_future      <- projectModelByYearAndScenario(spp="ssp_vaseyana", m=m_vaseyana)
bromus_tectorum_suitability_future   <- projectModelByYearAndScenario(spp="bromus_tectorum", m=m_bromus_tectorum)

# write our rasters to disk
saveRasters(prefix="ssp_tridentata",   ssp_tridentata_suitability_future)
saveRasters(prefix="ssp_wyomingensis", ssp_wyomingensis_suitability_future)
saveRasters(prefix="ssp_vaseyana",     ssp_vaseyana_suitability_future)
saveRasters(prefix="bromus_tectorum",  bromus_tectorum_suitability_future)

# save and exit
dir.create("models")
save(list=ls(pattern="^m_|^ssp_|^bromus_"),file="models/ktaylor_modeling_workspace.rdata",compress=T) # check out response curves for exp. variables, etc...
