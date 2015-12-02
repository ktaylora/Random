#
# Workflow for Implementing N-Mixture Model with Breeding Bird Survey Data
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

HOME <<- Sys.getenv("HOME")

#
# calculate site-level landscape metrics for the focal route [parallelized -- verified working KT]
#

processFocalRoute <- function(route=NULL){
  require(landscapeAnalysis)
  require(habitatWorkbench)
  require(rgdal)
  require(raster)

  HOME <- Sys.getenv("HOME")

  names(route@data) <- tolower(names(route@data)) #
  cat(" -- reading landcover raster\n")
  landcover <- raster(paste(sep="",HOME,"/PLJV/landcover/orig/Final_LC_8bit.tif"))
  cat(" -- sampling around vertices for route:",unique(route$rteno),"\n")
  points  <- habitatWorkbench::sampleAroundVertices(s=route,maxDist=2200,n=350) # generate a number of sampling points around each route to derive our site-level landscape metrics
  cat(" -- n=",length(points),"points generated\n")
  cat(" -- generating buffered raster surface around sample points for route:",unique(route$rteno),"\n")
  buffers <- landscapeAnalysis::subsampleSurface(x=landcover,pts=points, width=750) # sample buffers from our source landcover dataset with the points derived from the current route
  # parse our buffers out into focal landcover types
  cat(" -- reclassifying landcover buffers to binary for route:",unique(route$rteno),"\n")
  buffers_grassland     <- landscapeAnalysis::lReclass(buffers,inValues=c(31,37,39,71,75))
  buffers_agriculture   <- landscapeAnalysis::lReclass(buffers,inValues=c(38,201,202,203,205,206,207,208,209,210,211,212))
  buffers_shrubland     <- landscapeAnalysis::lReclass(buffers,inValues=c(83,85,87,81,82))
  # calculate the requisite landscape metrics for each cover type
  cat(" -- calculating grassland total area and mean patch area metrics for route:",unique(route$rteno),"\n")
  metrics_grassland <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_grassland)[[1]] # need: total area and mean patch area
  grassland_total_area <- try(metricsListToVector(metrics_grassland,'total.area'))
    if(class(grassland_total_area)=="try-error"){
      print(str(metrics_grassland))
      grassland_total_area <- 0
    } else {
      grassland_total_area <- ifelse(is.null(grassland_total_area),0,mean(grassland_total_area,na.rm=T))
        if(is.na(grassland_total_area)) grassland_total_area <- 0
    }
  grassland_mean_patch_area <- try(metricsListToVector(metrics_grassland,'mean.patch.area'))
    if(class(grassland_mean_patch_area)=="try-error"){
      print(str(metrics_grassland))
      grassland_mean_patch_area <- 0
    } else {
      grassland_mean_patch_area <- ifelse(is.null(grassland_mean_patch_area),0,mean(grassland_mean_patch_area,na.rm=T))
        if(is.na(grassland_mean_patch_area)) grassland_mean_patch_area <- 0
    }
  cat(" -- calculating agriculture total area for route:",unique(route$rteno),"\n")
  metrics_agriculture   <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_agriculture) # need: total area
  agriculture_total_area <- try(metricsListToVector(metrics_agriculture,'total.area'))
    if(class(agriculture_total_area)=="try-error"){
      print(str(metrics_agriculture))
      agriculture_total_area <- 0
    } else {
      agriculture_total_area <- ifelse(is.null(agriculture_total_area),0,mean(agriculture_total_area,na.rm=T))
        if(is.na(agriculture_total_area)) agriculture_total_area <- 0
    }
  cat(" -- calculating shrubland total area for route:",unique(route$rteno),"\n")
  metrics_shrubland     <- landscapeAnalysis::lCalculateLandscapeMetrics(buffers_shrubland) # need: total area
  shrubland_total_area <- try(metricsListToVector(metrics_shrubland,'total.area'))
    if(class(shrubland_total_area)=="try-error"){
      print(str(metrics_shrubland))
      shrubland_total_area <- 0
    } else {
      shrubland_total_area <- ifelse(is.null(shrubland_total_area),0,mean(shrubland_total_area,na.rm=T))
        if(is.na(shrubland_total_area)) shrubland_total_area <- 0
    }
  # write to output table
  cat(" -- returning data.frame from thread for route:",unique(route$rteno),"\n")
  return(data.frame(route=unique(route$rteno),grass_total_area=grassland_total_area,grass_mean_patch_area=grassland_mean_patch_area,
    ag_total_area=agriculture_total_area,shrub_total_area=shrubland_total_area))
}

#
# MAIN
#

require(rgdal)
require(raster)
require(habitatWorkbench)
require(parallel); cl <- makeCluster(getOption("cl.cores", parallel::detectCores()-1),outfile='outfile.log')

# read-in our local (study area) routes
# t_routes_bcr1819 <<- read.csv(paste(sep="",HOME,"PLJV/species_data/bbs_data/bbs_routes_bcr1819.csv")) # pre-define the route numbers needed in an external GIS

# THE PUBLISHED ROUTE NUMBERS WITH THE 1999 BBS SHAPEFILE IS OUT OF DATE -- DO NOT USE IT
# read-in the national bbs routes data and parse accordingly
# data(bbsRoutes); s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_routes_bcr1819$RTENO,] # RTENO : [STATE(2-digit)][ROUTE(3-digit)]

# ensure we have adequate land cover data for our bbs routes
lc <- raster(paste(sep="",HOME,"/PLJV/landcover/orig/Final_LC_8bit.tif"))
 s_bbsRoutes <- spTransform(readOGR(paste(sep="",HOME,"/PLJV/species_data/bbs_data/"),"Breeding_Bird_Surveys_1998to2013",verbose=F),CRS(projection(lc)))
    lc <- extract(lc,as(spTransform(s_bbsRoutes,CRS(projection(lc))),'SpatialPointsDataFrame'),sp=T)

 s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$rteno %in% unique(lc[!is.na(lc$Final_LC_8bit),]$rteno),] # make sure that our route data has landcover data available
   s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$rteno %in% unique(lc[lc$Final_LC_8bit != 0,]$rteno),]

 s_bbsRoutes_sampling <- s_bbsRoutes[!duplicated(s_bbsRoutes$rteno),] # there are duplicate points for each year.  For sampling LMs, let's just trim this to a single year.
# split our routes so that they can be processed in parallel on a multi-core machine and then calculate some site level statistics for landcover at each route
s_bbsRoutes_sampling <- split(s_bbsRoutes_sampling,f=as.factor(s_bbsRoutes_sampling$rteno))
if(!file.exists("site_level_parameters.csv")){
  cat(" -- sampling and processing BBS routes\n");
    out  <- parLapply(cl=cl,fun=processFocalRoute,X=s_bbsRoutes_sampling)
      out <- do.call(rbind,out) # bind our list into a single data.frame
        write.csv(out,"site_level_parameters.csv",row.names=F)
} else {
  out <- read.csv("site_level_parameters.csv")
}

## CALCULATE OUR ISOLATION METRICS AT 1650 AND 30000 METER BUFFER DISTANCES

lc <- raster(paste(sep="",HOME,"/PLJV/landcover/orig/Final_LC_8bit.tif"))
centroids <- parLapply(cl=cl,s_bbsRoutes_sampling,fun=getBbsRouteCoords,centroid=T)
  buffers <- subsampleSurface(x=lc,
                              pts=SpatialPoints(do.call(rbind,centroids),proj4string=CRS("+proj=utm +zone=14 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")),
                              width=3300/2)
    buffers <- lReclass(buffers,inValues=c(31,37,39,71,75))
t_isolation <- calcPatchIsolation(buffers,parallel=T)
out$isolation_1650m <- as.vector(unlist(t_isolation))

  buffers <- subsampleSurface(x=lc,
                              pts=SpatialPoints(do.call(rbind,centroids),proj4string=CRS("+proj=utm +zone=14 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")),
                              width=30000)
    buffers <- lReclass(buffers,inValues=c(31,37,39,71,75))
t_isolation <- calcPatchIsolation(buffers,parallel=T)
out$isolation_30000m <- as.vector(unlist(t_isolation))

write.csv(out,"site_level_parameters.csv",row.names=F)

# FETCH BBS COUNT DATA AND ASSOCIATE IT WITH THE FOCAL ROUTES
# NOTE: The default shapefile published by USGS is horribly out-of-date.  Don't use it.
#
# habitatWorkbench:::.fetchBbsStopData();
# t_counts <- habitatWorkbench:::.unpackToCSV(list.files(pattern="F.*.zip")); # ROUTES here are designated by their 3-digit code only. See the accompanying statenum column
#   t_counts <- t_counts[t_counts$year %in% 1999:2014,]
#
# route_code <- as.numeric(substring(s_bbsRoutes$RTENO, nchar(s_bbsRoutes$RTENO)-2,nchar(s_bbsRoutes$RTENO))) # 3-digit route code
# state_code <- as.numeric(substring(s_bbsRoutes$RTENO, 1,nchar(s_bbsRoutes$RTENO)-3)) # 2-digit state code
#
# t_counts <- t_counts[t_counts$statenum %in% state_code & t_counts$Route %in% route_code,] # parse our route counts to routes within the focal region
#   t_counts$RTENO <- paste(sprintf("%02d",t_counts$statenum),sprintf("%03d",t_counts$Route),sep="")
#
# s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_counts$RTENO,] # parse out the 'in region' routes we have count data for
#    t_counts <- t_counts[t_counts$RTENO %in% s_bbsRoutes$RTENO,]    # reciprocate
#
# # find the routes that have non-zero counts for our focal species, keeping the zero routes.
# n <- names(t_counts)
# nonZero <- t_counts[which(t_counts$AOU == as.numeric(habitatWorkbench:::.sppToAou("long-billed curlew"))),]
#    zero <- t_counts[!(t_counts$RTENO %in% nonZero$RTENO),]
#      zero[,n[grepl(n,pattern="Stop")]] <- 0
#        zero <- zero[!duplicated(zero$RTENO),] # remove duplicate zero values for routes
#          zero$AOU <- as.numeric(habitatWorkbench:::.sppToAou("long-billed curlew"))
#
#    t_counts <- rbind(nonZero,zero);rm(list=c('nonZero','zero'))
#      t_counts <- t_counts[!duplicated(t_counts),] # let's really make sure we don't have duplicate entries in the table, eh?
# s_bbsRoutes <- s_bbsRoutes[s_bbsRoutes$RTENO %in% t_counts$RTENO,]
#
# # transpose counts by year
# output <- data.frame()
# routes <- unique(t_counts$RTENO)
# for(y in sort(unique(t_counts$year))){
#   # parse stops for the focal year
#   cnts_focal <- data.frame(RTENO=routes,NA)
#       counts <- as.numeric(rowSums(t_counts[t_counts$year == y,n[grepl(n,pattern="Stop")]],na.rm=F))
#         counts <- counts[!duplicated(t_counts[t_counts$year == y,"RTENO"])] # make sure we don't have a looming duplicate route
#   stop_routes <- t_counts[t_counts$year == y,"RTENO"]
#       cnts_focal[match(cnts_focal$RTENO,stop_routes,nomatch=F),2] <- counts
#         names(cnts_focal) <- c("RTENO",paste("cnt",y,sep="."))
#           cnts_focal <- cnts_focal[!duplicated(cnts_focal$RTENO),]
#   # merge columns
#   if(!nrow(output)){
#     output <- cnts_focal
#   } else {
#     output <- cbind(output,unique(cnts_focal))
#   }
# }; names <- names(output); output <- cbind(RTENO=routes,output[,names[grepl(names,pattern="cnt[.]")]])

# BIND TO OUR LANDSCAPE METRICS TO MAKE A FULL TABLE OF SITE-LEVEL COVARIATES
# out <- read.csv("site_level_parameters.csv")
#   out <- out[!duplicated(out$route),]
#     out <- out[out$route %in% output$RTENO,]
#       out <- out[order(out$route),]
#       output <- output[output$RTENO %in% out$route,]
#         output <- output[order(output$RTENO),]
#           out <- cbind(out[match(out$route,output$RTENO),],output)
#             out<-out[,!duplicated(names(out))]
#               if(sum(out$route != out$RTENO) > 0){ cat("error -- failed to match routes between tables.")}
#                 out <- out[,names(out) != "RTENO"]
#
# # for show, we can merge our tabular data into our spatial data
# s<-s_bbsRoutes[!duplicated(s_bbsRoutes@data$RTENO),]; s@data <- out[which(out$route %in% s_bbsRoutes$RTENO),]

out <- out[order(out$route,decreasing=F),]
t_counts <- s_bbsRoutes[,c('LBCU','Year','rteno')]@data
  t_counts <- t_counts[order(t_counts$rteno,decreasing=F),]
    t_counts <- split(t_counts,f=as.factor(t_counts$rteno))
counts_table <- matrix(rep(NA,47*length(1998:2013)),ncol=length(1998:2013))
  colnames(counts_table) <- paste("cnts",1998:2013,sep=".")
  rownames(counts_table) <- out$route
for(i in 1:47){
  focal_count <- do.call(rbind,lapply(split(t_counts[[i]],f=as.factor(t_counts[[i]]$Year)),FUN=colSums))
  counts_table[i,grepl(x=colnames(counts_table),pattern=paste(rownames(focal_count),collapse="|"))] <- focal_count[,1]
}; rm(t_counts);
counts_table <- data.frame(counts_table)
  counts_table$route <- rownames(counts_table)


out <- cbind(out,counts_table[,1:(ncol(counts_table)-1)])


# GRAB ACCOMPANYING TABULAR DATA FOR FITTING COEFFICIENTS FOR DETECTION
detection_covariates <- list()
# calculate noise and cars present
habitatWorkbench:::.fetchVehicleData()
unzip("VehicleSummary.zip")
t <- read.csv("VehicleSummary.csv")
  t$RTENO <- paste(t$state,sprintf("%03d", t$Route),sep="")
    t <- t[t$RTENO %in% out$route,]
      t <- t[,c("RTENO","Year","StopTotal","NoiseTotal")]
# transpose noise and cars by year
t_data_full <- matrix()
for(y in 1998:2013){
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
detection_covariates[[3]] <- matrix(rep(letters[seq(1,length(1998:2013))],nrow(t_data_full)),nrow=nrow(t_data_full))
names(detection_covariates) <- c("noise","cars","timeTrend_") # name our list for compatibility with 'unmarked'

# build a model in 'unmarked'
require(unmarked)
t_pco <- unmarked::unmarkedFramePCO(y=out[,8:ncol(out)],siteCovs=out[,2:7], obsCovs=detection_covariates, numPrimary=length(1998:2013))

# m_lbcu_zip <- pcountOpen(lambdaformula=~grass_total_area+grass_mean_patch_area+ag_total_area+shrub_total_area+isolation_1650m+isolation_30000m,
#                      gammaformula=~1,
#                      omegaformula=~1,
#                      pformula=~noise+cars+timeTrend_,
#                      data=t_pco,mixture='ZIP',se=T, K=floor(max(max(out[,grepl(names(out),pattern="cnt")],na.rm=T))*2.1))
#
# m_lbcu_nb <- pcountOpen(lambdaformula=~grass_total_area+grass_mean_patch_area+ag_total_area+shrub_total_area+isolation_1650m+isolation_30000m,
#                     gammaformula=~1,
#                     omegaformula=~1,
#                     pformula=~noise+cars+timeTrend_,
#                     data=t_pco,mixture='NB',se=T, K=floor(max(max(out[,grepl(names(out),pattern="cnt")],na.rm=T))*2.1))

# pois model demonstrates ~7% improvement over null model in explaining residual error.  ZIP demonstrates similar performance, and negative binomial
# doesn't work for this question.  Perhaps if we took the ratio of counts to number of BBS stops, NB would be useful. -Kyle
m_lbcu_pois <- pcountOpen(lambdaformula=~grass_total_area+grass_mean_patch_area+ag_total_area+shrub_total_area+isolation_1650m+isolation_30000m,
                    gammaformula=~1,
                    omegaformula=~1,
                    pformula=~noise+cars+timeTrend_,
                    data=t_pco,mixture='P',se=T, K=floor(max(max(out[,grepl(names(out),pattern="cnt")],na.rm=T))*2.1))

# worse than null model
# m_lbcu <- pcountOpen(lambdaformula=~grass_total_area+grass_mean_patch_area+ag_total_area+shrub_total_area+isolation_3.3km+isolation_30km,
#                     gammaformula=~1,
#                     omegaformula=~1,
#                     pformula=~grass_total_area+grass_mean_patch_area+ag_total_area+shrub_total_area+isolation_3.3km+isolation_30km,
#                     data=t_pco,mixture='P',se=T, K=floor(max(max(out[,grepl(names(out),pattern="cnt")],na.rm=T))*1.3))

# % error variance explained by model vs null
model_sse <- sum(as.vector(m_lbcu@data@y-fitted(m_lbcu))^2, na.rm=T)
  null=ceiling(median(as.numeric(unlist(na.omit(out[,grepl(names(out),pattern="cnt")]))))) # just take the median as a null count model
    null_sse <- sum((out[,grepl(names(out),pattern="cnt")]-null)^2,na.rm=T)
cat(" -- pseudo-rsquared:",1-(model_sse/null_sse),"\n")

## CALCULATE LANDSCAPE METRICS WITHIN GPLCC CD UNITS FOR BBS MODELS
## CALCULATE LANDSCAPE METRICS FOR NON-PARAMETRIC MODELS
## EXTRAPOLATE BBS MODELS

ratioCU_to_BBS_Route <- (4.253/((3.1415926*(39.4289/50/2)^2)*50))

m_buow <- function(x,type="normalized"){
  x <- rowSums(sweep(x[,c('gr_tAr','gr_patAr','ag_tAr','sh_tAr','is_1650','is_30km')],MARGIN=2,as.numeric(c(0.000717,0.001151,0.001088,0.001379,-0.000491,-0.004840)),`*`))
    x <- as.vector(x+1.997769) # intercept
      x[x<0] <- 0 # don't return irrational counts
    # normalize to regional counts?
    if(grepl(tolower(type),pattern="normal")){
      x <- x/mean(x) # normalize to the mean count
        x <- x/max(x) # scale to values consistent with RF
    # scale our counts from the area of a BBS route to the area of our conservation delivery units
    } else if (grepl(tolower(type),pattern="counts")){
      x <- x*ratioCU_to_BBS_Route
    }
  return(x)
}

m_feha<- function(x,type="normalized"){
  x <- rowSums(sweep(x[,c('gr_tAr','gr_patAr','ag_tAr','sh_tAr','is_1650','is_30km')], MARGIN=2,as.numeric(c(-0.001144,0.000341,-0.002519,-0.001254,-0.000543,0.005354)),`*`))
    x <- as.vector(x+1.164134) # intercept
      x[x<0] <- 0 # don't return irrational counts
    # normalize to regional counts?
    if(grepl(tolower(type),pattern="normal")){
      x <- x/mean(x) # normalize to the mean count
        x <- x/max(x) # scale to values consistent with RF
    # scale our counts from the area of a BBS route to the area of our conservation delivery units
    } else if (grepl(tolower(type),pattern="counts")){
      x <- x*ratioCU_to_BBS_Route
    }
  return(x)
}

m_lbcu<- function(x,type="normalized"){
  x <- rowSums(sweep(x[,c('gr_tAr','gr_patAr','ag_tAr','sh_tAr','is_1650','is_30km')], MARGIN=2,as.numeric(c(2.81e-03,-4.58e-03,3.21e-04,6.46e-05,3.31e-03,-4.07e-04)),`*`))
    x <- as.vector(x+7.36e-02) # intercept
      x[x<0] <- 0 # don't return irrational counts
    # normalize to regional counts?
    if(grepl(tolower(type),pattern="normal")){
      x <- x/mean(x) # normalize to the mean count
        x <- x/max(x) # scale to values consistent with RF
    # scale our counts from the area of a BBS route to the area of our conservation delivery units
    } else if (grepl(tolower(type),pattern="counts")){
      x <- x*ratioCU_to_BBS_Route
    }
  return(x)
}

## EXTRAPOLATE NON-PARAMETRIC MODELS
