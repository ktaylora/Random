#
# WORKFLOW : Big Sagebrush Subspecies and Invasive Brome Climate SDMs
#
# This is a full-run generating content outlined in chapters 2/3 of my thesis. I wrote
# on a unix machine. You may need to switch "/" with "\" where I slip-up if
# you are running this on windows.
#
# Author : Kyle Taylor (kyle.taylor@uwyo.edu) [2014] (https://github.com/ktaylora)
#

library("foreign")

if (!require("raster")){
  install.packages("raster", repos = "https://cran.revolutionanalytics.com")
  library("raster")
}
if (!require("rgdal")){
  install.packages("rgdal", repos = "https://cran.revolutionanalytics.com")
  library("rgdal")
}
if (!require("randomForest")){
  install.packages("randomForest", repos = "https://cran.revolutionanalytics.com")
  library("randomForest")
}
if (!require("rfUtilities")){
  install.packages("rfUtilities", repos = "https://cran.revolutionanalytics.com")
  library("rfUtilities")
}
if (!require("rvest")){
  install.packages("rvest", repos = "https://cran.revolutionanalytics.com")
  library("rvest")
}
if (!require("rgeos")){
  install.packages("rvest", repos = "https://cran.revolutionanalytics.com")
  library("rgeos")
}
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
scrapeEddmapsData <- function(x){
  return(NA)
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
    target_pts_upsampled <- pts_generated[sample(1:nrow(pts_generated),size=pts_to_generate),]
  }

  return(rbind(s,target_pts_upsampled))
}

################################################################################
##
#   MAIN
##
################################################################################

#
# fetch our explanatory climate data
#

# GCM's selected to maximise the range of predictions for temperature/precipitation
# made across all available scenarios/models for the inter-mountain region.
# See VisTrails-SAHM documentation for an overview.
gcm_abbrevs  <- c("ac","gf","in","ip")
time_periods <- c("current","50","70")

dir.create("climate_data")

urls <- scrapeWorldclimCmip5(models=gcm_abbrevs,
                             bioclim=T,
                             start_of_century=sum(grepl(time_periods,pattern="current")),
                             mid_century=sum(grepl(time_periods,pattern="50")),
                             end_of_century=sum(grepl(time_periods,pattern="70")),
                             scen_45=TRUE,
                             scen_85=TRUE
                            )


#
# fetch and read-in omernik ecoregions, national boundaries, and county-level bromus tectorum
# data to 1.) mask relavent GAP observations and 2.) define our absense space (minding the naughty-naughts)
# [https://www.jstor.org/stable/3683828]
#

# pre-define climate and topographic conditions thought to be limiting to big sagebrush,
# but still inclusive of our occurrence data (SEMIARID HIGHLANDS)
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
full_study_region      <- rgeos::gUnion(region_boundary_us, sp::spTransform(region_boundary, CRS(projection(region_boundary))))
  writeOGR(SpatialPolygonsDataFrame(full_study_region, data=data.frame(id=1)),
                                    "boundaries", "full_study_region_extent",
                                    driver="ESRI Shapefile", overwrite=T)

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
# record hold-out data.
#

cat(" -- performing a binary reclassification of GAP explanatory data for our (sub)species of interest\n")
if(sum(grepl(names(gap_rat),pattern="ECOLSYS_LU"))==0) stop("couldn't find ECOLSYS_LU in GAP attribute table")

cat(" -- a. tridentata ssp tridentata...\n")
ssp_tridentata <- parse_gap_training_data("Inter-Mountain Basins Big Sagebrush Shrubland", gap=gap, gap_rat=gap_rat)
  ssp_tridentata <- mask_by_boundary(ssp_tridentata, region_boundary=region_boundary_us) # (US + Omernik)
# merge-in any GBIF training records that are available
if(class(ssp_tridentata_gbif$training)!="try-error"){
  ssp_tridentata <- rbind(sp::spTransform(ssp_tridentata_gbif$training,
                          CRS(projection(ssp_tridentata))), ssp_tridentata)
}
# generate pseudo-absences in areas unsampled by GAP (Canada + Mexico + Omernik)
ssp_tridentata <- generate_pseudoabsences(ssp_tridentata,
                                          src_region_boundary=region_boundary_us,
                                          target_region_boundary=region_boundary_non_us)
# dupe check
ssp_tridentata <- ssp_tridentata[!duplicated(ssp_tridentata@coords),]
# write to disk for validation
writeOGR(ssp_tridentata,"training_data","ssp_tridentata_gap_training_data",driver="ESRI Shapefile",overwrite=T)

cat(" -- a. tridentata ssp wyomingensis...\n")
ssp_wyomingensis <- parse_gap_training_data("Inter-Mountain Basins Big Sagebrush Steppe", gap=gap, gap_rat=gap_rat)
  ssp_wyomingensis <- mask_by_boundary(ssp_wyomingensis, region_boundary=region_boundary_us)
# merge-in any GBIF training records that are available
if(class(ssp_wyomingensis_gbif$training)!="try-error"){
  ssp_wyomingensis <- rbind(sp::spTransform(ssp_wyomingensis_gbif$training,
                          CRS(projection(ssp_wyomingensis))), ssp_wyomingensis)
}
# generate pseudo-absences in areas unsampled by GAP (Canada + Mexico + Omernik)
ssp_wyomingensis <- generate_pseudoabsences(ssp_wyomingensis,
                                            src_region_boundary=region_boundary_us,
                                            target_region_boundary=region_boundary_non_us)
# dupe check
ssp_wyomingensis <- ssp_wyomingensis[!duplicated(ssp_wyomingensis@coords),]
# write to disk for validation
writeOGR(ssp_wyomingensis,"training_data","ssp_wyomingensis_gap_training_data",driver="ESRI Shapefile", overwrite=T)

cat(" -- a. tridentata ssp vaseyana...\n")
ssp_vaseyana <- parse_gap_training_data("Inter-Mountain Basins Montane Sagebrush Steppe", gap=gap, gap_rat=gap_rat)
  ssp_vaseyana <- mask_by_boundary(ssp_vaseyana, region_boundary=region_boundary_us)
# merge-in any GBIF training records that are available
if(class(ssp_vaseyana_gbif$training)!="try-error"){
  ssp_vaseyana <- rbind(sp::spTransform(ssp_vaseyana_gbif$training,
                          CRS(projection(ssp_vaseyana))), ssp_vaseyana)
}
# generate pseudo-absences in areas unsampled by GAP (Canada + Mexico + Omernik)
ssp_vaseyana <- generate_pseudoabsences(ssp_vaseyana,
                                        src_region_boundary=region_boundary_us,
                                        target_region_boundary=region_boundary_non_us)
# dupe check
ssp_vaseyana <- ssp_vaseyana[!duplicated(ssp_vaseyana@coords),]
# write to disk for validation
writeOGR(ssp_vaseyana,"training_data","ssp_vaseyana_gap_training_data",driver="ESRI Shapefile", overwrite=T)

cat(" -- b. tectorum...\n")
bromus_tectorum <- parse_gap_training_data("Introduced Upland Vegetation - Annual Grassland", gap=gap, gap_rat=gap_rat)
  bromus_tectorum <- mask_by_boundary(bromus_tectorum, region_boundary=region_boundary_us)
  # mask our GAP data against county-level EDDSmap data
  #writeOGR(bromus_tectorum,"training_data","bromus_tectorum_gap_training_data",driver="ESRI Shapefile")




# process our herbarium records (for augmentation at us boundary and validation)
