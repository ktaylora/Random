#
# Workflow for Implementing N-Mixture Model with Breeding Bird Survey Data
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

HOME <- Sys.getenv("HOME")

#
# calculate site-level landscape metrics for the focal route [parallelized]
#

processFocalRoute <- function(route=NULL){
  require(landscapeAnalysis)
  require(habitatWorkbench)
  require(rgdal)
  require(raster)
  landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")
  t_routes_bcr1819 <- read.csv("/home/ktaylora/PLJV/species_data/bbs_data/bbs_routes_bcr1819.csv")
  points  <- habitatWorkbench::sampleAroundVertices(s=route,maxDist=2200,n=250) # generate a number of sampling points around each route to derive our site-level landscape metrics
  buffers <- landscapeAnalysis::subsampleSurface(x=landcover,pts=points, width=750) # sample buffers from our source landcover dataset with the points derived from the current route
  # parse our buffers out into focal landcover types
  buffers_grassland     <- landscapeAnalysis::lReclass(buffers,inValues=c(31,37,39,71,75))
  buffers_agriculture   <- landscapeAnalysis::lReclass(buffers,inValues=c(38,201,202,203,205,206,207,208,209,210,211,212))
  buffers_shrubland     <- landscapeAnalysis::lReclass(buffers,inValues=c(83,85,87,81,82))
  # calculate the requisite landscape metrics for each cover type
  metrics_grassland <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_grassland)[[1]] # need: total area and mean patch area
  grassland_total_area <- try(metricsListToVector(metrics_grassland,'total.area'))
    if(class(grassland_total_area)=="try-error"){
      print(str(metrics_grassland))
      grassland_total_area <- 0
    } else {
      grassland_total_area <- ifelse(is.null(grassland_total_area),0,grassland_total_area)
        grassland_total_area[is.na(grassland_total_area)] <- 0
          grassland_total_area <- mean(grassland_total_area)
    }
  grassland_mean_patch_area <- try(metricsListToVector(metrics_grassland,'mean.patch.area'))
    if(class(grassland_mean_patch_area)=="try-error"){
      print(str(metrics_grassland))
      grassland_mean_patch_area <- 0
    } else {
      grassland_mean_patch_area <- ifelse(is.null(grassland_mean_patch_area),0,grassland_mean_patch_area)
        grassland_mean_patch_area[is.na(grassland_mean_patch_area)] <- 0
          grassland_mean_patch_area <- mean(grassland_mean_patch_area)
    }
  metrics_agriculture   <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_agriculture) # need: total area
  agriculture_total_area <- try(metricsListToVector(metrics_agriculture,'total.area'))
    if(class(agriculture_total_area)=="try-error"){
      print(str(metrics_agriculture))
      agriculture_total_area <- 0
    } else {
      agriculture_total_area <- ifelse(is.null(agriculture_total_area),0,agriculture_total_area)
        agriculture_total_area[is.na(agriculture_total_area)] <- 0
          agriculture_total_area <- mean(agriculture_total_area)
    }
  metrics_shrubland     <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_shrubland) # need: total area
  shrubland_total_area <- try(metricsListToVector(metrics_shrubland,'total.area'))
    if(class(shrubland_total_area)=="try-error"){
      print(str(metrics_shrubland))
      shrubland_total_area <- 0
    } else {
      shrubland_total_area <- ifelse(is.null(shrubland_total_area),0,shrubland_total_area)
        shrubland_total_area[is.na(shrubland_total_area)] <- 0
          shrubland_total_area <- mean(shrubland_total_area)
    }
  # write to output table
  return(data.frame(route=route$RTENO,grass_total_area=grassland_total_area,grass_mean_patch_area=grassland_mean_patch_area,
    ag_total_area=agriculture_total_area,shrub_total_area=shrubland_total_area))
}

#
# MAIN
#

require(rgdal)
require(raster)
require(habitatWorkbench)
require(parallel); cl <- makeCluster(getOption("cl.cores", 7),outfile='outfile.log')

# read-in our local (study area) routes
t_routes_bcr1819 <<- read.csv("/home/ktaylora/PLJV/species_data/bbs_data/bbs_routes_bcr1819.csv") # pre-define the route numbers needed in an external GIS

# read-in the national bbs routes data and parse accordingly
data(bbsRoutes); s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_routes_bcr1819$RTENO,] # RTENO : [STATE(2-digit)][ROUTE(3-digit)]

# ensure we have adequate land cover data for our bbs routes
landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")
  landcover <- extract(landcover,as(spTransform(s_bbsRoutes,CRS(projection(landcover))),'SpatialPointsDataFrame'),sp=T)
s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% unique(landcover[!is.na(landcover$Final_LC_8bit),]$RTENO),] # make sure that our route data has landcover data available
  s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% unique(landcover[landcover$Final_LC_8bit != 0,]$RTENO),]

# split our routes so that they can be processed in parallel on a multi-core machine and then calculate some site level statistics for landcover at each route
s_bbsRoutes <- split(s_bbsRoutes,f=1:nrow(s_bbsRoutes))
if(!file.exists("site_level_parameters.csv")){
  cat(" -- sampling and processing BBS routes\n");
    out  <- parLapply(cl=cl,fun=processFocalRoute,X=s_bbsRoutes)
      out <- do.call(rbind,out) # bind our list into a single data.frame
        write.csv(out,"site_level_parameters.csv",row.names=F)
} else {
  out <- read.csv("site_level_parameters.csv")
}

## CALCULATE OUR ISOLATION METRICS AT 1650 AND 30000 METER BUFFER DISTANCES

if(!grepl(names(out),pattern="isolation_3.3")){
      landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")
    s_bbsRoutes <-lapply(s_bbsRoutes,spTransform,CRSobj=CRS(projection(landcover)))
      centroids <- parLapply(cl=cl,s_bbsRoutes,fun=getBbsRouteLocations,centroid=T) # calculate the centroid of each route
  buffers_3.3km <- parLapply(cl=cl,X=centroids,fun=rgeos::gBuffer,width=3300/2)      # calculate our 3.3 km regions
    m <- mcmapply(buffers_3.3km, FUN=raster::crop, MoreArgs=list(x=landcover))
      buffers_3.3km <- mcmapply(FUN=raster::mask, x=m, mask=buffers_3.3km); rm(m);
        buffers_3.3km <- lReclass(buffers_3.3km,inValues=c(31,37,39,71,75))
          buffers_3.3km <- parLapply(cl=cl,buffers_3.3km,fun=landscapeAnalysis::rasterToPolygons)
            # check for NA values
            na_values<-as.vector(unlist(lapply(X=buffers_3.3km,FUN=is.na)))
            if(sum(na_values)>0){
              buffers_3.3km[na_values] <- buffers_3.3km[which(!na_values)[1]] # overwrite our NA values with something valid
              buffers_3.3km<-lapply(X=buffers_3.3km,FUN=getSpPPolygonsLabptSlots)

            } else {
              buffers_3.3km<-lapply(X=buffers_3.3km,FUN=getSpPPolygonsLabptSlots)
            }
            # do our NN assessment
            d <- function(x,na.rm=F){ o<-try(FNN::knn.dist(x,k=1)); if(class(o) != "try-error") { x <- o; } else { x[[i]] <- NA }; return(x)}
            buffers_3.3km <- lapply(buffers_3.3km,FUN=d);rm(d);
              buffers_3.3km <- lapply(buffers_3.3km, FUN=mean)
                buffers_3.3km[na_values] <- NA # restore our NA values
                  buffers_3.3km
                    out$isolation_3.3km<-as.vector(unlist(buffers_3.3km))
}

if(!grepl(names(out),pattern="isolation_30km")){
      landcover <- raster("/home/ktaylora/PLJV/landcover/orig/Final_LC_8bit.tif")
    s_bbsRoutes <-lapply(s_bbsRoutes,spTransform,CRSobj=CRS(projection(landcover)))
      centroids <- parLapply(cl=cl,s_bbsRoutes,fun=getBbsRouteLocations,centroid=T) # calculate the centroid of each route
  buffers_30km <- parLapply(cl=cl,X=centroids,fun=rgeos::gBuffer,width=30000)      # calculate our 3.3 km regions
    m <- mcmapply(buffers_30km, FUN=raster::crop, MoreArgs=list(x=landcover))
      buffers_30km <- mcmapply(FUN=raster::mask, x=m, mask=buffers_30km); rm(m);
        buffers_30km <- lReclass(buffers_30km,inValues=c(31,37,39,71,75))
          buffers_30km <- parLapply(cl=cl,buffers_30km,fun=landscapeAnalysis::rasterToPolygons)
            # check for NA values
            na_values<-as.vector(unlist(lapply(X=buffers_30km,FUN=is.na)))
            if(sum(na_values)>0){
              buffers_30km[na_values] <- buffers_30km[which(!na_values)[1]] # overwrite our NA values with something valid
              buffers_30km<-lapply(X=buffers_30km,FUN=getSpPPolygonsLabptSlots)
            } else {
              buffers_30km<-lapply(X=buffers_30km,FUN=getSpPPolygonsLabptSlots)
            }
            # do our NN assessment
            d <- function(x,na.rm=F){ o<-try(FNN::knn.dist(x,k=1)); if(class(o) != "try-error") { x <- o; } else { x[[i]] <- NA }; return(x)}
            buffers_30km <- lapply(buffers_30km,FUN=d);rm(d);
              buffers_30km <- lapply(buffers_30km, FUN=mean)
                buffers_30km[na_values] <- NA # restore our NA values
                    out$isolation_30km<-as.vector(unlist(buffers_30km))
}

write.csv(out,"site_level_parameters.csv",row.names=F)

# FETCH BBS COUNT DATA AND ASSOCIATE IT WITH THE FOCAL ROUTES

habitatWorkbench:::.fetchBbsStopData();
t_counts <- habitatWorkbench:::.unpackToCSV(list.files(pattern="F.*.zip")); # ROUTES here are designated by their 3-digit code only. See the accompanying statenum column
  t_counts <- t_counts[t_counts$year %in% 1999:2014,]

route_code <- as.numeric(substring(s_bbsRoutes$RTENO, nchar(s_bbsRoutes$RTENO)-2,nchar(s_bbsRoutes$RTENO))) # 3-digit route code
state_code <- as.numeric(substring(s_bbsRoutes$RTENO, 1,nchar(s_bbsRoutes$RTENO)-3)) # 2-digit state code

t_counts <- t_counts[t_counts$statenum %in% state_code & t_counts$Route %in% route_code,] # parse our route counts to routes within the focal region
  t_counts$RTENO <- paste(sprintf("%02d",t_counts$statenum),sprintf("%03d",t_counts$Route),sep="")

s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_counts$RTENO,] # parse out the 'in region' routes we have count data for
   t_counts <- t_counts[t_counts$RTENO %in% s_bbsRoutes$RTENO,]    # reciprocate

# find the routes that have non-zero counts for our focal species, keeping the zero routes.
n <- names(t_counts)
nonZero <- t_counts[which(t_counts$AOU == as.numeric(habitatWorkbench:::.sppToAou("long-billed curlew"))),]
   zero <- t_counts[!(t_counts$RTENO %in% nonZero$RTENO),]
     zero[,n[grepl(n,pattern="Stop")]] <- 0
       zero <- zero[!duplicated(zero$RTENO),] # remove duplicate zero values for routes
         zero$AOU <- as.numeric(habitatWorkbench:::.sppToAou("long-billed curlew"))

   t_counts <- rbind(nonZero,zero);rm(list=c('nonZero','zero'))
     t_counts <- t_counts[!duplicated(t_counts),] # let's really make sure we don't have duplicate entries in the table, eh?
s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_counts$RTENO,]

# transpose counts by year
output <- data.frame()
routes <- unique(t_counts$RTENO)
for(y in sort(unique(t_counts$year))){
  # parse stops for the focal year
  cnts_focal <- data.frame(RTENO=routes,NA)
      counts <- as.numeric(rowSums(t_counts[t_counts$year == y,n[grepl(n,pattern="Stop")]],na.rm=F))
        counts <- counts[!duplicated(t_counts[t_counts$year == y,"RTENO"])] # make sure we don't have a looming duplicate route
  stop_routes <- t_counts[t_counts$year == y,"RTENO"]
      cnts_focal[match(cnts_focal$RTENO,stop_routes,nomatch=F),2] <- counts
        names(cnts_focal) <- c("RTENO",paste("cnt",y,sep="."))
          cnts_focal <- cnts_focal[!duplicated(cnts_focal$RTENO),]
  # merge columns
  if(!nrow(output)){
    output <- cnts_focal
  } else {
    output <- cbind(output,unique(cnts_focal))
  }
}; names <- names(output); output <- cbind(RTENO=routes,output[,names[grepl(names,pattern="cnt[.]")]])
# BIND TO OUR LANDSCAPE METRICS TO MAKE A FULL TABLE OF SITE-LEVEL COVARIATES
out <- read.csv("site_level_parameters.csv")
  out <- out[!duplicated(out$route),]
    out <- out[out$route %in% output$RTENO,]
      out <- out[order(out$route),]
      output <- output[output$RTENO %in% out$route,]
        output <- output[order(output$RTENO),]
          out <- cbind(out[match(out$route,output$RTENO),],output)
            out<-out[,!duplicated(names(out))]
              if(sum(out$route != out$RTENO) > 0){ cat("error -- failed to match routes between tables.")}
                out <- out[,names(out) != "RTENO"]

# for show, we can merge our tabular data into our spatial data
s<-s_bbsRoutes; s@data <- out[which(out$route %in% s_bbsRoutes$RTENO),]
# GRAB ACCOMPANYING TABULAR DATA FOR FITTING COEFFICIENTS FOR DETECTION
detection_covariates <- list()
# calculate noise and cars present
habitatWorkbench:::.fetchVehicleData()
unzip("VehicleSummary.zip")
t <- read.csv("VehicleSummary.csv")
  t$RTENO <- paste(t$state,sprintf("%03d", t$Route),sep="")
    t <- t[t$RTENO %in% out$route,]
      t <- t[,c("RTENO","Year","StopTotal","NoiseTotal")]
# transpose by year
t_data_full <- matrix()
for(y in 1999:2014){
  # create a matrix with NAs for all routes
  t_data <- matrix(rep(NA,2*length(out$route)),ncol=2)
    rownames(t_data) <- out$route
  # parse the noise and car presences for THIS year
  t_focal <- as.matrix(t[t$Year == y,c(3,4)])
    colnames(t_focal) <- paste(names(t)[3:4],y,sep=".")
    rownames(t_focal) <- t[t$Year == y,]$RTENO
    t_focal <- t_focal[!duplicated(rownames(t_focal)),] # ensure we don't have duplicate routes
  # assign data to those routes that actually have observations this year
  t_data[which(rownames(t_data) %in% rownames(t_focal)),] <- t_focal[which(rownames(t_focal) %in% rownames(t_data)),]
    colnames(t_data) <- paste(names(t)[3:4],y,sep=".")
  if(ncol(t_data_full)>=2){
    t_data_full <- cbind(t_data_full,t_data)
  } else {
    t_data_full <- t_data
  }
}
# assign to a named list that 'unmarked' will understand
detection_covariates[[1]] <- as.matrix(data.frame(t_data_full)[,grepl(colnames(t_data_full),pattern="Noise")]) # Noise Observed at Routes
detection_covariates[[2]] <- as.matrix(data.frame(t_data_full)[,grepl(colnames(t_data_full),pattern="Stop")])  # Number of Cars Observed at Routes
# calculate distance to transmission lines along each route
# -- removed this because it doesn't vary across time.  Consider using it as a site-level covariate, perhaps
# require(raster)
# distanceToTrans <- raster(paste(HOME,"/PLJV/infrastructure/products/bcr_18_19_distanceToTransmissionLines.tif",sep=""))
#               t <- raster::extract(distanceToTrans,spTransform(s_bbsRoutes,CRS(projection(distanceToTrans)),df=T,progress='text')
#                 t <- unlist(lapply(t,FUN=mean))
# NOTE: An 'observer' variable is conspicuously missing here on probability of detection.  I don't want to make a flat assumption in extrapolating observer bias
# in geographic space.  Let's see how well the model does without it.

# calculate a dummy variable representing years in the time-series (one for each route)
detection_covariates[[3]] <- matrix(rep(letters[seq(1,length(years))],nrow(t_data_full)),nrow=nrow(t_data_full))
names(detection_covariates) <- c("noise","cars","timeTrend") # name our list for compatibility with 'unmarked'

# build a model in 'unmarked'
require(unmarked)
t_pco <- unmarked::unmarkedFramePCO(y=out[,8:ncol(out)],siteCovs=out[,2:7], obsCovs=detection_covariates, numPrimary=length(1999:2014))

m_lbcu <- pcountOpen(lambdaformula=~grass_total_area+grass_mean_patch_area+ag_total_area+shrub_total_area+isolation_3.3km+isolation_30km,
                     gammaformula=~1,
                     omegaformula=~1,
                     pformula=~noise+cars+timeTrend,
                     data=t_pco,mixture='P',se=T, K=floor(max(max(out[,grepl(names(out),pattern="cnt")],na.rm=T))*1.3))
