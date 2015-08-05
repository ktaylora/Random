#
# This is a kitchen sink. I am sorry. Kyle
#

require(randomForest)
require(raster)
require(rgdal)

#
# spatialSmoothing()
# Implements a gaussian smoothing window as implemented by Jeff Evans 2014 (see: http://evansmurphy.wix.com/evansspatial#!spatial-smoothing/ch1)
#

spatialSmoothing <- function(x, s=1, d=5, filename=FALSE, ...) {
    if (!require(sp)) stop("sp PACKAGE MISSING")
    if (!require(raster)) stop("raster PACKAGE MISSING")
    if (!require(rgdal)) stop("rgdal PACKAGE MISSING")
    if (!inherits(x, "RasterLayer")) stop("MUST BE RasterLayer OBJECT")
       GaussianKernel <- function(sigma=s, n=d) {
          m <- matrix(nc=n, nr=n)
            col <- rep(1:n, n)
            row <- rep(1:n, each=n)
          x <- col - ceiling(n/2)
          y <- row - ceiling(n/2)
         m[cbind(row, col)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
        m / sum(m)
       }
   if (filename != FALSE) {
      focal(x, w=GaussianKernel(sigma=s, n=d), filename=filename, ...)
      print(paste("RASTER WRITTEN TO", filename, sep=": "))
        } else {
      return(focal(x, w=GaussianKernel(sigma=s, n=d), ...))
  }
}

#
# Stole this fire from GDAL
#

rasterToPolygons <- function(r=NULL, method='gdal'){
  r_name=deparse(substitute(r))
  writeRaster(r,paste(r_name,"tif",sep="."),overwrite=T);
  unlink(paste(r_name,c("shp","xml","shx","prj","dbf"),sep="."))
  if(grepl(Sys.getenv("HOME"), pattern="home")){
    system(paste("gdal_polygonize.py -8",paste(r_name,"tif",sep="."),"-f \"ESRI Shapefile\"",paste(r_name,"shp",sep="."),sep=" "))
  } else {
    system(paste("/opt/local/bin/python2.7 /opt/local/share/doc/py27-gdal/examples/scripts/gdal_polygonize.py -8",paste(r_name,"tif",sep="."),"-f \"ESRI Shapefile\"",paste(r_name,"shp",sep="."),sep=" "))
  }
  return(rgdal::readOGR(".",r_name,verbose=F));
}

#
# Extract polygon contours for high and low suitability predictions from a raster surface
#

extractHighLowDensities <- function(x,s=5,d=15, p=c("0.50","0.90")){
  # extract  range contours for raster surface x
  q <- as.character(seq(0.05,0.95,0.05))
    if(sum(as.character(p) %in% as.character(q)) != length(p)) stop("quantiles are typically extracted in 0.05 interval steps")
  # smooth
  smoothed <- SpatialSmoothing(x,s=s)
    h <- hist(smoothed, plot=F)
      smoothed[smoothed<=h$mids[2]] <- NA
        quantiles <- as.vector(quantile(smoothed,probs=as.numeric(q)))

  out <- list()
  for(focal in p){
    smoothed_focal <- smoothed>=quantiles[which(as.character(q) == as.character(focal))]
      smoothed_focal <- match(smoothed_focal,1,nomatch=NA)
        out[[length(out)+1]] <- rasterToPolygons(smoothed_focal);
  }

  return(out)
}

extractPointsFromContour <- function(contour,r){
  rast <- raster(res=0.005)
    rast <- crop(rast,extent(contour))
  pts  <- rasterize(contour,rast)
    pts <- rasterToPoints(pts,spatial=T)

    t <- data.frame(val=extract(r,pts))
    pts@data <- t
    return(pts)
}

#
# returns area in kilometers squared
#
envelopeToSurfaceArea <- function(x) sum(unlist(lapply(spTransform(x, CRS(projection("+init=epsg:2163")))@polygons, FUN=slot, "area")))/1000000

HOME      <- Sys.getenv("HOME");

f<-list.files(paste(HOME,"/Products/weather/worldclim/30_sec/bioclim/",sep=""),pattern="bil$",full.names=T)
  f<-f[grepl(f, pattern="_11|18|12|_1[.]|15")]
    climate_data_current <- raster::stack(f)

f<-list.files(paste(HOME,"/Products/weather/worldclim/future/2050_RCP45",sep=""),pattern="tif$", full.names=T)
  f<-f[grepl(f, pattern="[.]11[.]|[.]18[.]|[.]12[.]|[.]01[.]|[.]15[.]")]
    climate_data_2050_RCP45 <- raster::stack(f)

f<-list.files(paste(HOME,"/Products/weather/worldclim/future/2050_RCP85",sep=""),pattern="tif$", full.names=T)
  f<-f[grepl(f, pattern="[.]11[.]|[.]18[.]|[.]12[.]|[.]01[.]|[.]15[.]")]
    climate_data_2050_RCP85 <- raster::stack(f)

f<-list.files(paste(HOME,"/Products/weather/worldclim/future/2070_RCP45",sep=""),pattern="tif$", full.names=T)
  f<-f[grepl(f, pattern="[.]11[.]|[.]18[.]|[.]12[.]|[.]01[.]|[.]15[.]")]
    climate_data_2070_RCP45 <- raster::stack(f)

f<-list.files(paste(HOME,"/Products/weather/worldclim/future/2070_RCP85",sep=""),pattern="tif$", full.names=T)
  f<-f[grepl(f, pattern="[.]11[.]|[.]18[.]|[.]12[.]|[.]01[.]|[.]15[.]")]
    climate_data_2070_RCP85 <- raster::stack(f)

pts       <- readOGR(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/vector/gap_abs_pts/",sep=""), "pres_abs_records_gap", verbose=F)
presences <- pts[pts$resp==1,]
absences  <- pts[pts$resp==0,]

b <- spTransform(readOGR(paste(HOME,"/Products/boundaries/",sep=""),"western_north_american_boundaries"), CRS(projection(climate_data_current)))

# build a suitability model under current conditions

t_pres <- data.frame(resp=1, raster::extract(climate_data_current,presences))
t_abs  <- data.frame(resp=0, raster::extract(climate_data_current,absences))

t     <- na.omit(rbind(t_pres,t_abs))
m_rf  <- randomForest(as.factor(resp) ~., do.trace=T, ntree=600, importance=T, data=t)
m_glm <- glm(formula=resp~I(bio_1^4)+I(bio_11^4)+I(bio_12^2)+I(bio_15^4)+I(bio_18^2), family=binomial(link="logit"),data=t) # I fit the orders to this polynomial using stepwise AIC selection

# project across the entire range

climate_data_current <- raster::crop(climate_data_current,b)

if(file.exists("/Volumes/NO NAME/a.tridentata.GLM.regap.random_AllData_Full_GLM.tif")){
  current_glm <- raster("/Volumes/NO NAME/a.tridentata.GLM.regap.random_AllData_Full_GLM.tif")
} else {
  current_glm <- predict(climate_data_current,m_glm,type="response",progress='text')
}

if(file.exists("/Volumes/NO NAME/a.tridentata.RF.regap.random_AllData_Full_RF.tif")){
  current_rf <- raster("/Volumes/NO NAME/a.tridentata.RF.regap.random_AllData_Full_RF.tif")
} else {
  current_rf <- predict(climate_data_current,m_rf,type="prob",progress='text')
    current_rf <- abs(current_rf-1)
}

# ggplot2 plots of current conditions models

current_glm_matrix<-rasterToPoints(current_glm)
  current_glm_matrix <- data.frame(current_glm_matrix)
    colnames(current_glm_matrix) <- c("Longitude","Latitude","Value")
current_rf_matrix <- rasterToPoints(current_rf)
  current_rf_matrix <- data.frame(current_rf_matrix)
    current_rf_matrix[,3] <- abs(current_rf_matrix[,3])
    colnames(current_rf_matrix) <- c("Longitude","Latitude","Value")

p1<-ggplot(data=current_glm_matrix, aes(y=Latitude, x=Longitude)) +
    geom_raster(aes(fill=Value)) +
    geom_polygon(data=fortify(b), aes(x = long, y = lat, group = group), fill=NA, color = "grey", size = 0.25) +
    scale_colour_brewer(name = "Probability of Occurrence") +
    coord_map(ylim=c(32.5,55), xlim=c(-128,-100)) +
    annotate("text", x = -100, y = 54, label = "A", color="white") +
    theme_bw() +
    coord_equal()

p2<-ggplot(data=current_rf_matrix, aes(y=Latitude, x=Longitude)) +
    geom_raster(aes(fill=Value)) +
    geom_polygon(data=fortify(b), aes(x = long, y = lat, group = group), fill=NA, color = "grey", size = 0.25) +
    scale_colour_brewer(name = "Probability of Occurrence") +
    coord_map(ylim=c(32.5,55), xlim=c(-128,-100)) +
    annotate("text", x = -100, y = 54, label = "B", color="white") +
    theme_bw() +
    coord_equal()

g<-arrangeGrob(p1, p2, ncol=2)

# "ensemble"
# current_ensemble <- (current_rf+current_glm)/2

# project into future CMIP5 conditions

# names(climate_data_2050_RCP45) <-c("bio_1","bio_11","bio_12","bio_15","bio_18")
# names(climate_data_2050_RCP85) <-c("bio_1","bio_11","bio_12","bio_15","bio_18")
# names(climate_data_2070_RCP45) <-c("bio_1","bio_11","bio_12","bio_15","bio_18")
# names(climate_data_2070_RCP85) <-c("bio_1","bio_11","bio_12","bio_15","bio_18")

# climate_data_2050_RCP45 <- raster::crop(climate_data_2050_RCP45, b)
# climate_data_2050_RCP85 <- raster::crop(climate_data_2050_RCP85, b)
# climate_data_2070_RCP45 <- raster::crop(climate_data_2070_RCP45, b)
# climate_data_2070_RCP85 <- raster::crop(climate_data_2070_RCP85, b)

# RCP45_rf_2050  <- predict(climate_data_2050_RCP45,m_rf,type="prob")
#   RCP45_rf_2050 <- abs(RCP45_rf_2050-1)
# RCP85_rf_2050  <- predict(climate_data_2050_RCP85,m_rf,type="prob")
#   RCP85_rf_2050 <- abs(RCP85_rf_2050-1)
# RCP45_rf_2070  <- predict(climate_data_2070_RCP45,m_rf,type="prob")
#   RCP45_rf_2070 <- abs(RCP45_rf_2070-1)
# RCP85_rf_2070  <- predict(climate_data_2070_RCP85,m_rf,type="prob")
#   RCP85_rf_2070 <- abs(RCP85_rf_2070-1)

# RCP45_glm_2050  <- predict(climate_data_2050_RCP45,m_glm,type="response")
# RCP85_glm_2050  <- predict(climate_data_2050_RCP85,m_glm,type="response")
# RCP45_glm_2070  <- predict(climate_data_2070_RCP45,m_glm,type="response")
# RCP85_glm_2070  <- predict(climate_data_2070_RCP85,m_glm,type="response")

rasterOption(timer=T,maxmemory=10000)

h <- extractHighLowDensities(RCP45_rf_2050)
  RCP45_rf_2050_90 <- h[[1]]
    RCP45_rf_2050_50 <- h[[2]]

h <- extractHighLowDensities(RCP45_rf_2070)
  RCP45_rf_2070_90 <- h[[1]]
    RCP45_rf_2070_50 <- h[[2]]

h <- extractHighLowDensities(RCP45_glm_2070)
  RCP45_glm_2070_90 <- h[[1]]
    RCP45_glm_2070_50 <- h[[2]]

h <- extractHighLowDensities(RCP45_glm_2050)
  RCP45_glm_2050_90 <- h[[1]]
    RCP45_glm_2050_50 <- h[[2]]

h <- extractHighLowDensities(RCP85_glm_2050)
  RCP85_glm_2050_90 <- h[[1]]
    RCP85_glm_2050_50 <- h[[2]]

h <- extractHighLowDensities(RCP85_glm_2070)
  RCP85_glm_2070_90 <- h[[1]]
    RCP85_glm_2070_50 <- h[[2]]

h <- extractHighLowDensities(RCP85_rf_2050)
  RCP85_rf_2050_90 <- h[[1]]
    RCP85_rf_2050_50 <- h[[2]]

h <- extractHighLowDensities(RCP85_rf_2070)
  RCP85_rf_2070_90 <- h[[1]]
    RCP85_rf_2070_50 <- h[[2]]

save.image(paste(HOME,"/a_tridentata_current_vs_future_ensemble.rdata",sep=""),compress=T)

##
## sample points within each p50 and p90 contour
##

elev <- raster("/Volumes/NO NAME/lf_elev1km/w001001.adf")

lat_elev_current_glm_50    <- extractPointsFromContour(contour=p_current_glm_60, r=elev)
lat_elev_current_glm_90    <- extractPointsFromContour(contour=p_current_glm_90, r=elev)
lat_elev_RCP85_glm_2050_50 <- extractPointsFromContour(contour=RCP85_glm_2050_50, r=elev)
lat_elev_RCP85_glm_2050_90 <- extractPointsFromContour(contour=RCP85_glm_2050_90, r=elev)

lat_elev_presence_bias_50_glm  <- extractPointsFromContour(contour=p_glm_presence_bias_50,r=elev)
  lat_elev_presence_bias_50_glm <- lat_elev_presence_bias_50_glm[!is.na(lat_elev_presence_bias_50_glm@data[,1]),]
 lat_elev_absence_bias_50_glm  <- extractPointsFromContour(contour=p_glm_absence_bias_50,r=elev)
   lat_elev_absence_bias_50_glm <- lat_elev_absence_bias_50_glm[!is.na(lat_elev_absence_bias_50_glm@data[,1]),]
 lat_elev_climate_bias_50_glm  <- extractPointsFromContour(contour=p_glm_absence_bias_50,r=elev)
   lat_elev_climate_bias_50_glm <- lat_elev_climate_bias_50_glm[!is.na(lat_elev_climate_bias_50_glm@data[,1]),]

lat_elev_presence_bias_90_glm  <- extractPointsFromContour(contour=p_glm_presence_bias_90_glm,r=elev)
  lat_elev_presence_bias_90_glm <- lat_elev_presence_bias_90_glm[!is.na(lat_elev_presence_bias_90_glm@data[,1]),]
 lat_elev_absence_bias_90_glm  <- extractPointsFromContour(contour=p_glm_absence_bias_90,r=elev)
   lat_elev_absence_bias_90_glm <- lat_elev_absence_bias_90_glm[!is.na(lat_elev_absence_bias_90_glm@data[,1]),]
 lat_elev_climate_bias_90_glm  <- extractPointsFromContour(contour=p_glm_climate_bias_90,r=elev)
   lat_elev_climate_bias_90_glm <- lat_elev_climate_bias_90_glm[!is.na(lat_elev_climate_bias_90_glm@data[,1]),]

lat_elev_current_rf_50    <- extractPointsFromContour(contour=p_current_rf_60, r=elev)
lat_elev_current_rf_90    <- extractPointsFromContour(contour=p_current_rf_90, r=elev)
lat_elev_presence_bias_50_rf  <- extractPointsFromContour(contour=p_rf_presence_bias_50,r=elev)
  lat_elev_presence_bias_50_rf <- lat_elev_presence_bias_50_rf[!is.na(lat_elev_presence_bias_50_rf@data[,1]),]
 lat_elev_absence_bias_50_rf  <- extractPointsFromContour(contour=p_rf_absence_bias_50,r=elev)
   lat_elev_absence_bias_50_rf <- lat_elev_absence_bias_50_rf[!is.na(lat_elev_absence_bias_50_rf@data[,1]),]
 lat_elev_climate_bias_50_rf  <- extractPointsFromContour(contour=p_rf_absence_bias_50,r=elev)
   lat_elev_climate_bias_50_rf <- lat_elev_climate_bias_50_rf[!is.na(lat_elev_climate_bias_50_rf@data[,1]),]

lat_elev_presence_bias_90_rf  <- extractPointsFromContour(contour=p_rf_presence_bias_90,r=elev)
  lat_elev_presence_bias_90_rf <- lat_elev_presence_bias_90_rf[!is.na(lat_elev_presence_bias_90_rf@data[,1]),]
 lat_elev_absence_bias_90_rf  <- extractPointsFromContour(contour=p_rf_absence_bias_90,r=elev)
   lat_elev_absence_bias_90_rf <- lat_elev_absence_bias_90_rf[!is.na(lat_elev_absence_bias_90_rf@data[,1]),]
 lat_elev_climate_bias_90_rf  <- extractPointsFromContour(contour=p_rf_climate_bias_90,r=elev)
   lat_elev_climate_bias_90_rf <- lat_elev_climate_bias_90_rf[!is.na(lat_elev_climate_bias_90_rf@data[,1]),]

require(gridExtra)
require(reshape2)
require(ggplot2)
require(grid)

grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}

make_subplot <- function(df=NULL,label=NULL,sub="elevation",xlab="Source of Uncertainty",ylab="Elevation (m)"){
  if(sub=="elevation"){
    ymax <- 3950;
  } else if(sub=="latitude"){
    ymax <- 50
  }
  p <- ggplot(df, aes(factor(uncertainty), get(sub))) + geom_boxplot(outlier.shape=NA) +
               xlab(xlab) + ylab(ylab) +
               annotate("text", x = 1.25, y = ymax, label = label, size=3.85) +
               theme_bw() + theme(legend.title=element_blank());
  if(sub=="elevation") { p <- p + ylim(0,4000) }
  if(sub=="latitude") { p <- p + ylim(35,50)}
  return(p)
}

sub <- "elevation"
# GLM P50 Interval (Elevation)
df <- data.frame(presence=sample(lat_elev_presence_bias_50_glm@data[,1],500),
                 absence=sample(lat_elev_absence_bias_50_glm@data[,1],500),
                 climate=sample(lat_elev_climate_bias_50_glm@data[,1],500),
                 control=sample(lat_elev_current_glm_50@data[,1],500))

df <- melt(df,variable.name="uncertainty",value.name="elevation")

p1 <- make_subplot(df,y="elevation",label="p=0.50 (elevation) [GLM]")

# GLM P90 Interval (Elevation)
df <- data.frame(presence=sample(lat_elev_presence_bias_90_glm@data[,1],500),
                 absence=sample(lat_elev_absence_bias_90_glm@data[,1],500),
                 climate=sample(lat_elev_climate_bias_90_glm@data[,1],500),
                 control=sample(lat_elev_current_glm_90@data[,1],500))

df <- melt(df,variable.name="uncertainty",value.name="elevation")

p2 <- make_subplot(df,y="elevation",label="p=0.90 (elevation) [GLM]")

# RF P50 Interval (Elevation)
df <- data.frame(presence=sample(lat_elev_presence_bias_50_rf@data[,1],500),
                 absence=sample(lat_elev_absence_bias_50_rf@data[,1],500),
                 climate=sample(lat_elev_climate_bias_50_rf@data[,1],500),
                 control=sample(lat_elev_current_rf_50@data[,1],500),na.rm=T)

df <- melt(df,variable.name="uncertainty",value.name="elevation")

p3 <- make_subplot(df,label="p=0.50 (elevation) [RF]")

# RF P90 Interval (Elevation)
df <- data.frame(presence=sample(lat_elev_presence_bias_90_rf@data[,1],500),
                 absence=sample(lat_elev_absence_bias_90_rf@data[,1],500),
                 climate=sample(lat_elev_climate_bias_90_rf@data[,1],500),
                 control=sample(lat_elev_current_rf_90@data[,1],500),na.rm=T)

df <- melt(df,variable.name="uncertainty",value.name="elevation")

p4 <- make_subplot(df,label="p=0.90 (elevation) [RF]")

grid.arrange(p1,p2,p3,p4)

require(raster)
require(rgdal)

# GLM P50 Interval (Latitude)

df <- data.frame(presence=sample(spTransform(lat_elev_presence_bias_50_glm,CRS(projection("+init=epsg:4326")))@coords[,2],500),
                 absence=sample(spTransform(lat_elev_absence_bias_50_glm,CRS(projection("+init=epsg:4326")))@coords[,2],500),
                 climate=sample(spTransform(lat_elev_climate_bias_50_glm,CRS(projection("+init=epsg:4326")))@coords[,2],500),
                 control=sample(spTransform(lat_elev_current_glm_50,CRS(projection("+init=epsg:4326")))@coords[,2],500))

df <- melt(df,variable.name="uncertainty",value.name="latitude")

sub <- "latitude" # why is this happening?
p1 <- make_subplot(df,sub='latitude',label="p=0.50 (latitude) [GLM]",ylab="Latitude (degrees)")

# GLM P90 Interval (Latitude)
df <- data.frame(presence=sample(spTransform(lat_elev_presence_bias_90_glm,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 absence=sample(spTransform(lat_elev_absence_bias_90_glm,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 climate=sample(spTransform(lat_elev_climate_bias_90_glm,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 control=sample(spTransform(lat_elev_current_glm_90,CRS(projection("+init=epsg:4326")))@coords[,2],500))

df <- melt(df,variable.name="uncertainty",value.name="latitude")
p2 <- make_subplot(df,sub='latitude',label="p=0.90 (latitude) [GLM]",ylab="Latitude (degrees)")

# RF P50 Interval (Latitude)

df <- data.frame(presence=sample(spTransform(lat_elev_presence_bias_50_rf,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 absence=sample(spTransform(lat_elev_absence_bias_50_rf,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 climate=sample(spTransform(lat_elev_climate_bias_50_rf,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 control=sample(spTransform(lat_elev_current_rf_50,CRS(projection("+init=epsg:4326")))@coords[,2],500))

df <- melt(df,variable.name="uncertainty",value.name="latitude")
p3 <- make_subplot(df,sub='latitude',label="p=0.50 (latitude) [RF]",ylab="Latitude (degrees)")

# RF P90 Interval (Latitude)

df <- data.frame(presence=sample(spTransform(lat_elev_presence_bias_90_rf,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 absence=sample(spTransform(lat_elev_absence_bias_90_rf,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 climate=sample(spTransform(lat_elev_climate_bias_90_rf,CRS(projection("+init=epsg:4326")))@coords[,2],500,replace=T),
                 control=sample(spTransform(lat_elev_current_rf_90,CRS(projection("+init=epsg:4326")))@coords[,2],500))

df <- melt(df,variable.name="uncertainty",value.name="latitude")
p4 <- make_subplot(df,sub='latitude',label="p=0.90 (latitude) [RF]",ylab="Latitude (degrees)")

grid.arrange(p1, p2, p3, p4)


# "ensemble"
# RCP45_ens_2050 <-(RCP45_rf_2050+RCP45_glm_2050)/2
# RCP85_ens_2050 <-(RCP85_rf_2050+RCP85_glm_2050)/2
# RCP45_ens_2070 <-(RCP45_rf_2070+RCP45_glm_2070)/2
# RCP85_ens_2070 <-(RCP85_rf_2070+RCP85_glm_2070)/2

# difference_RCP85_glm_2070 <- -1*((current_glm*(current_glm>0.2))-(RCP85_glm_2070*(RCP85_glm_2070>0.2)))
# difference_RCP45_glm_2050 <- -1*((current_glm*(current_glm>0.2))-(RCP45_glm_2050*(RCP45_glm_2050>0.2)))
# difference_RCP85_glm_2050 <- -1*((current_glm*(current_glm>0.2))-(RCP85_glm_2050*(RCP85_glm_2050>0.2)))
# difference_RCP45_glm_2070 <- -1*((current_glm*(current_glm>0.2))-(RCP45_glm_2070*(RCP45_glm_2070>0.2)))

# difference_RCP85_rf_2070 <- -1*((current_rf*(current_rf>0.2))-(RCP85_rf_2070*(RCP85_rf_2070>0.2)))
# difference_RCP45_rf_2050 <- -1*((current_rf*(current_rf>0.2))-(RCP45_rf_2050*(RCP45_rf_2050>0.2)))
# difference_RCP85_rf_2050 <- -1*((current_rf*(current_rf>0.2))-(RCP85_rf_2050*(RCP85_rf_2050>0.2)))
# difference_RCP45_rf_2070 <- -1*((current_rf*(current_rf>0.2))-(RCP45_rf_2070*(RCP45_rf_2070>0.2)))

# for(i in ls(pattern="difference_RCP")){ f<-get(i); f[f==0]<-NA; assign(i,f); }

save.image(paste(HOME,"/a_tridentata_current_vs_future_ensemble.rdata",sep=""),compress=T)

#
# do plots
#

## Plot Current VS Future Conditions

# COLORS <- colorRampPalette(colors=c("red","yellow","green"))

# current_glm_contours <- rasterToContour(current_glm)
# current_rf_contours  <- rasterToContour(current_rf)

# png(filename="climate_current_vs_2050_CMIP5.png", units="px", height=800, width=1200);
# par(mfrow=c(2,2))
# plot(difference_RCP45_glm_2050, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(GLM) RCP 4.5 (2050)",cex=0.8);
# plot(b, ylim=c(33,53), xlim=c(-130,-90),add=T);
# plot(current_glm_contours[seq(1,nrow(current_glm_contours),2),], col="#34343499",add=T, ylim=c(33,53), xlim=c(-130,-90))
# plot(difference_RCP85_glm_2050, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(GLM) RCP 8.5 (2050)",cex=0.8);
# plot(b, ylim=c(33,53), xlim=c(-130,-90),add=T);
# plot(current_glm_contours[seq(1,nrow(current_glm_contours),2),], col="#34343499",add=T, ylim=c(33,53), xlim=c(-130,-90))

# plot(difference_RCP45_rf_2050, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(RF) RCP 4.5 (2050)",cex=0.8);
# plot(b, ylim=c(33,53), xlim=c(-130,-90),add=T);
# plot(current_rf_contours[seq(1,nrow(current_rf_contours),3),], col="#34343499",add=T, ylim=c(33,53), xlim=c(-130,-90))
# plot(difference_RCP85_rf_2050, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(RF) RCP 8.5 (2050)",cex=0.8);
# plot(b, ylim=c(33,53), xlim=c(-130,-90),add=T);
# plot(current_rf_contours[seq(1,nrow(current_rf_contours),3),], col="#34343499",add=T, ylim=c(33,53), xlim=c(-130,-90))

# graphics.off();

# png(filename="climate_current_vs_2070_CMIP5.png", height=800, width=1200);
# par(mfrow=c(2,2))
# plot(difference_RCP45_glm_2070, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(GLM) RCP 4.5 (2070)",cex=0.8);
# plot(b, ylim=c(33,53), xlim=c(-130,-90),add=T);
# plot(current_glm_contours[seq(1,nrow(current_glm_contours),2),], col="#34343499",add=T, ylim=c(33,53), xlim=c(-130,-90))
# plot(difference_RCP85_glm_2070, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(GLM) RCP 8.5 (2070)",cex=0.8);
# plot(b, ylim=c(33,53), xlim=c(-130,-90),add=T);#34343499
# plot(current_glm_contours[seq(1,nrow(current_glm_contours),2),], col="#34343499",add=T, ylim=c(33,53), xlim=c(-130,-90))

# plot(difference_RCP45_rf_2070, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(RF) RCP 4.5 (2070)",cex=0.8);
# plot(b, ylim=c(33,53), xlim=c(-130,-90),add=T);
# plot(current_rf_contours[seq(1,nrow(current_rf_contours),3),], col="#34343499",add=T, ylim=c(33,53), xlim=c(-130,-90))
# plot(difference_RCP85_rf_2070, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(RF) RCP 8.5 (2070)",cex=0.8);
# plot(b, ylim=c(33,53), xlim=c(-130,-90),add=T);
# plot(current_rf_contours[seq(1,nrow(current_rf_contours),3),], col="#34343499",add=T, ylim=c(33,53), xlim=c(-130,-90))

# graphics.off()

## Plot Current and Future Conditions w/ Uncertainty Surfaces

herbarium_pts <- readOGR(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/vector/herbarium_presence_pts",sep=""), "presences")

# absence bias surfaces
glm_abs_bias         <- raster(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/raster/products/a.tridentata.GLM.absence.bias.ensemble.tif",sep=""))
rf_abs_bias          <- raster(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/raster/products/a.tridentata.RF.absence.bias.ensemble.tif",sep=""))

# presence bias surfaces
glm_pres_bias         <- raster(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/raster/products/a.tridentata.GLM.presence.bias.exp.ensemble.tif",sep=""))
rf_pres_bias          <- raster(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/raster/products/a.tridentata.RF.presence.bias.exp.ensemble.tif",sep=""))

# climate bias surfaces
glm_climate_bias         <- raster(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/raster/products/a.tridentata.glm.experimental.climate.bias.ensemble.tif",sep=""))
rf_climate_bias          <- raster(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/raster/products/a.tridentata.rf.experimental.climate.bias.ensemble.tif",sep=""))

# presence bias surface calculation
# glm_pres_bias_control <- current_glm # use the unbiased current model as a control
#   o <- -1*((glm_pres_bias_control*max(glm_pres_bias))-glm_pres_bias) # rescale the binomial response to the log-difference units of our biased run series. Take the difference of the control and the biased ensemble
#     samplePts<-sampleRandom(current_glm >= as.numeric(quantile(current_glm,probs=0.90)),sp=T,size=500)
#     cutoff <- mean(extract(o,samplePts[samplePts$layer==1,]))
#       o[o<=cutoff] <- NA
#   glm_pres_bias_final <- o; rm(o,samplePts,cutoff);
# # absence bias surface calculation
# glm_abs_bias_control <- current_glm
#   o <- -1*((glm_abs_bias_control*max(glm_abs_bias))-glm_abs_bias)
#     samplePts<-sampleRandom(current_glm >= as.numeric(quantile(current_glm,probs=0.90)),sp=T,size=500)
#     cutoff <- mean(extract(o,samplePts[samplePts$layer==1,]))
#       o[o<=cutoff] <- NA
#   glm_abs_bias_final <- o; rm(o,samplePts,cutoff);
# # climate bias surface calculation
# glm_climate_bias_control <- current_glm
#   o <- -1*((glm_climate_bias_control*max(glm_climate_bias))-glm_climate_bias)
#     samplePts<-sampleRandom(current_glm >= as.numeric(quantile(current_glm,probs=0.90)),sp=T,size=500)
#     cutoff <- mean(extract(o,samplePts[samplePts$layer==1,]))
#       o[o<=cutoff] <- NA
#   glm_climate_bias_final <- o; rm(o,samplePts,cutoff);


rasterOption(timer=T,maxmemory=10000)

# extract hith and low range contours for GLM
# current_glm_smoothed<-SpatialSmoothing(current_glm,s=5)
# h <- hist(current_glm_smoothed, plot=F)
#   current_glm_smoothed[current_glm_smoothed<=h$mids[2]] <- NA
#     quantiles <- as.vector(quantile(current_glm_smoothed,probs=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)))

# r_current_glm_90 <- current_glm_smoothed>=quantiles[10]
#   r_current_glm_90 <- match(r_current_glm_90,1,nomatch=NA)
#     p_current_glm_90 <- rasterToPolygons(r_current_glm_90);

# r_current_glm_50 <- current_glm_smoothed>=quantiles[6]
#   r_current_glm_50 <- match(r_current_glm_50,1,nomatch=NA)
#     p_current_glm_50 <- rasterToPolygons(r_current_glm_50);

# # extract high and low range contours for RF
# current_rf_smoothed<-SpatialSmoothing(current_rf,d=25,s=5)

# h <- hist(current_rf_smoothed,plot=F)
#   current_rf_smoothed[current_rf_smoothed<=h$mids[1]] <- NA
#     quantiles <- as.vector(quantile(current_rf_smoothed,probs=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)))

# r_current_rf_90 <- current_rf_smoothed>=quantiles[10]
#   r_current_rf_90 <- match(r_current_rf_90,1,nomatch=NA)
#     #p_current_rf_90 <- rasterToPolygons(r_current_rf_90,na.rm=T,dissolve=T)
#     p_current_rf_90 <- rasterToPolygons(r_current_rf_90);

# r_current_rf_50 <- current_rf_smoothed>=quantiles[6]
#   r_current_rf_50 <- match(r_current_rf_50,1,nomatch=NA)
#     #p_current_rf_50 <- rasterToPolygons(r_current_rf_60,na.rm=T,dissolve=T)
#     p_current_rf_50 <- rasterToPolygons(r_current_rf_50)

# extract contours for presence bias runs

glm_pres_bias          <- mask(crop(glm_pres_bias,b),b)
glm_pres_bias_smoothed <- SpatialSmoothing(glm_pres_bias,s=5)

o <- extractHighLowDensities(glm_pres_bias_smoothed,p=c(0.90,0.10))
  p_glm_presence_bias_90 <- o[[1]]
  p_glm_presence_bias_10 <- o[[2]]

rf_pres_bias          <- mask(crop(rf_pres_bias,b),b)
rf_pres_bias_smoothed <- SpatialSmoothing(rf_pres_bias,s=5)

o <- extractHighLowDensities(rf_pres_bias_smoothed,p=c(0.90,0.10))
  p_rf_presence_bias_90 <- o[[1]]
  p_rf_presence_bias_10 <- o[[2]]

# absences
glm_abs_bias <- mask(crop(glm_abs_bias,b),b)
glm_abs_bias_smoothed <- SpatialSmoothing(glm_abs_bias,s=5)

o <- extractHighLowDensities(glm_abs_bias_smoothed,p=c(0.90,0.10))
  p_glm_absence_bias_90 <- o[[1]]
  p_glm_absence_bias_10 <- o[[2]]

rf_abs_bias <- mask(crop(rf_abs_bias,b),b)
rf_abs_bias_smoothed <- SpatialSmoothing(rf_abs_bias,s=5)

o <- extractHighLowDensities(rf_abs_bias_smoothed,p=c(0.90,0.10))
  p_rf_absence_bias_90 <- o[[1]]
  p_rf_absence_bias_10 <- o[[2]]

# climate
glm_climate_bias <- mask(crop(glm_climate_bias,b),b)
glm_climate_bias_smoothed <- SpatialSmoothing(glm_climate_bias,s=5)

o <- extractHighLowDensities(glm_climate_bias_smoothed,p=c(0.90,0.10))
  p_glm_climate_bias_90 <- o[[1]]
  p_glm_climate_bias_10 <- o[[2]]

rf_climate_bias_smoothed <- mask(crop(rf_climate_bias,b),b)
rf_climate_bias_smoothed <- SpatialSmoothing(rf_climate_bias_smoothed,d=25, s=5)

o <- extractHighLowDensities(rf_climate_bias_smoothed,p=c(0.90,0.10))
  p_rf_climate_bias_90 <- o[[1]]
  p_rf_climate_bias_10 <- o[[2]]

## plot them

require(gridExtra)

#RF, current conditions, with envelopes

layers <- list();
layers[[1]] <- fortify(b)
layers[[2]] <- fortify(rgeos::gSymdifference(p_rf_climate_bias_90,p_rf_climate_bias_10))   # climate bias envelope
layers[[3]] <- fortify(rgeos::gSymdifference(p_rf_presence_bias_90,p_rf_presence_bias_10)) # presence bias envelope
layers[[4]] <- fortify(rgeos::gSymdifference(p_rf_absence_bias_90,p_rf_absence_bias_10))   # absence bias envelope
layers[[5]] <- fortify(p_rf_current_tss_threshold)
layers[[6]] <- cbind(data.frame(herbarium_pts@coords),id=herbarium_pts$id)                 # points

# climate, presence, and absence bias (GLM)

dev.new(height=4.4, width=11)
par(mfrow=c(1,4))
par(mar=par()$mar/2)

plot(p_glm_climate_bias_10,border=NA,add=F,col="white",axes=T, xlim=c(-134,-100), ylim=c(27,60))
  grid(lwd=0.8,lty=1)
    plot(p_glm_current_tss_threshold, border=NA, col="#009900", add=T)
      plot(p_glm_climate_bias_90, border=NA, col="#FFFF00", add=T)
        plot(rgeos::gSymdifference(p_glm_climate_bias_90,p_glm_climate_bias_10), col="#C0C0C099", border=NA, add=T)
          plot(b, add=T)
            text(x=-106.5,y=62, "A) Climate", cex=1.5,col="#000000")

plot(p_glm_presence_bias_10,border=NA,add=F,col="white",axes=T, xlim=c(-134,-100), ylim=c(27,60))
  grid(lwd=0.8,lty=1)
    plot(p_glm_current_tss_threshold, border=NA, col="#009900", add=T)
      plot(p_glm_presence_bias_90, border=NA, col="#FFFF00", add=T)
        plot(rgeos::gSymdifference(p_glm_presence_bias_90,p_glm_presence_bias_10), col="#C0C0C099", border=NA, add=T)
          plot(b, add=T)
            text(x=-109,y=62, "B) Presence", cex=1.5,col="#000000")

plot(p_glm_absence_bias_10,border=NA,add=F,col="white",axes=T, xlim=c(-134,-100), ylim=c(27,60))
  grid(lwd=0.8,lty=1)
    plot(p_glm_current_tss_threshold, border=NA, col="#009900", add=T)
      plot(p_glm_absence_bias_90, border=NA, col="#FFFF00", add=T)
        plot(rgeos::gSymdifference(p_glm_absence_bias_90,p_glm_absence_bias_10), col="#C0C0C099", border=NA, add=T)
          plot(b, add=T)
            text(x=-109,y=62, "C) Absence", cex=1.5,col="#000000")

o <- rgeos::gIntersection(p_glm_climate_bias_90,p_glm_presence_bias_90)
  o <- rgeos::gIntersection(o@polyobj,p_glm_absence_bias_90)@polyobj

plot(o, add=F,col="white",axes=T, border=NA, xlim=c(-134,-100), ylim=c(27,60))
  grid(lwd=0.8,lty=1)
    plot(o, col="#0000FF",border=NA,add=T)
      plot(b, add=T)
        text(x=-109,y=62, "D) Consensus", cex=1.5,col="#000000",border=NA)

# climate, presence, and absence bias (RF)

dev.new(height=4.4, width=11)
par(mfrow=c(1,4))
par(mar=par()$mar/2)

plot(p_rf_climate_bias_10,border=NA,add=F,col="white",axes=T, xlim=c(-134,-100), ylim=c(27,60))
  grid(lwd=0.8,lty=1)
    plot(p_rf_current_tss_threshold, border=NA, col="#009900", add=T)
      plot(p_rf_climate_bias_90, border=NA, col="#FFFF00", add=T)
        plot(rgeos::gSymdifference(p_rf_climate_bias_90,p_rf_climate_bias_10), col="#C0C0C099", border=NA, add=T)
          plot(b, add=T)
            text(x=-106.5,y=62, "E) Climate", cex=1.5,col="#000000")

plot(p_rf_presence_bias_10,border=NA,add=F,col="white",axes=T, xlim=c(-134,-100), ylim=c(27,60))
  grid(lwd=0.8,lty=1)
    plot(p_rf_current_tss_threshold, border=NA, col="#009900", add=T)
      plot(p_rf_presence_bias_90, border=NA, col="#FFFF00", add=T)
        plot(rgeos::gSymdifference(p_rf_presence_bias_90,p_rf_presence_bias_10), col="#C0C0C099", border=NA, add=T)
          plot(b, add=T)
            text(x=-109,y=62, "F) Presence", cex=1.5,col="#000000")

plot(p_rf_absence_bias_10,border=NA,add=F,col="white",axes=T, xlim=c(-134,-100), ylim=c(27,60))
  grid(lwd=0.8,lty=1)
    plot(p_rf_current_tss_threshold, border=NA, col="#009900", add=T)
      plot(p_rf_absence_bias_90, border=NA, col="#FFFF00", add=T)
        plot(rgeos::gSymdifference(p_rf_absence_bias_90,p_rf_absence_bias_10), col="#C0C0C099", border=NA, add=T)
          plot(b, add=T)
            text(x=-109,y=62, "G) Absence", cex=1.5,col="#000000")

o <- rgeos::gIntersection(p_rf_climate_bias_90,p_rf_presence_bias_90)
  o <- rgeos::gIntersection(o@polyobj,p_rf_absence_bias_90)@polyobj

plot(o, add=F,col="white",axes=T, border=NA, xlim=c(-134,-100), ylim=c(27,60))
  grid(lwd=0.8,lty=1)
    plot(o, col="#0000FF", border=NA, add=T)
      plot(b, add=T)
        text(x=-109,y=62, "H) Consensus", cex=1.5,col="#000000")

# try to get this to work with ggplot

p2<-ggplot()+ geom_polygon(data=layers[[5]], aes(x = long, y = lat, group = group), fill="#009900", alpha=1, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[3]], aes(x = long, y = lat, group = group), fill="blue", alpha=0.5, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[1]], aes(x = long, y = lat, group = group), fill=NA, color = "black", size = 0.25) +
              geom_point(data=layers[[6]], aes(x=coords.x1,y=coords.x2), size=0.3, shape=5, alpha=0.35) +
              coord_map(ylim=c(30,58), xlim=c(-130,-97)) +
              theme(plot.margin=unit(c(0.5,0.5,-0.5,0.5), "cm")) +
              ylab("") + xlab("") +
              annotate("text", x = -102.5, y = 54, label = "B") +
              theme_bw()

p3<-ggplot()+ geom_polygon(data=layers[[5]], aes(x = long, y = lat, group = group), fill="#009900", alpha=1, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[4]], aes(x = long, y = lat, group = group), fill="yellow", alpha=0.5, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[1]], aes(x = long, y = lat, group = group), fill=NA, color = "black", size = 0.25) +
              geom_point(data=layers[[6]], aes(x=coords.x1,y=coords.x2), size=0.3, shape=5, alpha=0.35) +
              coord_map(ylim=c(30,58), xlim=c(-130,-97)) +
              theme(plot.margin=unit(c(0.5,0.5,-0.5,0.5), "cm")) +
              ylab("") + xlab("") +
              annotate("text", x = -102.5, y = 54, label = "A") +
              theme_bw()

layers <- list();
layers[[1]] <- fortify(b)
layers[[2]] <- fortify(p_rf_climate_bias_90)
layers[[3]] <- fortify(p_rf_presence_bias_90)
layers[[4]] <- fortify(p_rf_absence_bias_90)
layers[[5]] <- cbind(data.frame(herbarium_pts@coords),id=herbarium_pts$id) # points

p2<-ggplot()+ geom_polygon(data=layers[[2]], aes(x = long, y = lat, group = group), fill="#003380", alpha=1, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[3]], aes(x = long, y = lat, group = group), fill="#0066FF", alpha=0.75, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[4]], aes(x = long, y = lat, group = group), fill="#6685E0", alpha=0.65, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[1]], aes(x = long, y = lat, group = group), fill=NA, color = "black", size = 0.25) +
              geom_point(data=layers[[5]], aes(x=coords.x1,y=coords.x2), size=0.2, shape=8, alpha=0.35) +
              coord_map(ylim=c(32.5,55), xlim=c(-128,-100)) +
              theme(plot.margin=unit(c(0.5,0.5,-0.5,0.5), "cm")) +
              ylab("") + xlab("") +
              annotate("text", x = -102.5, y = 54, label = "B") +
              theme_bw()


#GLM, current conditions, "highly suitable climate" vs "moderately suitable climate"

layers <- list();
layers[[1]] <- fortify(b)
layers[[2]] <- fortify(p_glm_climate_bias_50)
layers[[3]] <- fortify(p_glm_presence_bias_50)
layers[[4]] <- fortify(p_glm_absence_bias_50)
layers[[5]] <- cbind(data.frame(herbarium_pts@coords),id=herbarium_pts$id) # points

p3<-ggplot()+ geom_polygon(data=layers[[2]], aes(x = long, y = lat, group = group), fill="#003380", alpha=1, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[3]], aes(x = long, y = lat, group = group), fill="#0066FF", alpha=0.75, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[4]], aes(x = long, y = lat, group = group), fill="#6685E0", alpha=0.65, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[1]], aes(x = long, y = lat, group = group), fill=NA, color = "black", size = 0.25) +
              geom_point(data=layers[[5]], aes(x=coords.x1,y=coords.x2), size=0.2, shape=8, alpha=0.35) +
              coord_map(ylim=c(32.5,55), xlim=c(-128,-100)) +
              theme(plot.margin=unit(c(0.5,0.5,-0.5,0.5), "cm")) +
              ylab("") + xlab("") +
              annotate("text", x = -102.5, y = 54, label = "C") +
              theme_bw()

layers[[2]] <- fortify(p_glm_climate_bias_90)
layers[[3]] <- fortify(p_glm_presence_bias_90)
layers[[4]] <- fortify(p_glm_absence_bias_90)

p4<-ggplot()+ geom_polygon(data=layers[[2]], aes(x = long, y = lat, group = group), fill="#003380", alpha=1, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[3]], aes(x = long, y = lat, group = group), fill="#0066FF", alpha=0.75, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[4]], aes(x = long, y = lat, group = group), fill="#6685E0", alpha=0.65, color = "grey", size = 0.05) +
              geom_polygon(data=layers[[1]], aes(x = long, y = lat, group = group), fill=NA, color = "black", size = 0.25) +
              geom_point(data=layers[[5]], aes(x=coords.x1,y=coords.x2), size=0.2, shape=8, alpha=0.35) +
              coord_map(ylim=c(32.5,55), xlim=c(-128,-100)) +
              theme(plot.margin=unit(c(0.5,0.5,-0.5,0.5), "cm")) +
              ylab("") + xlab("") +
              annotate("text", x = -102.5, y = 54, label = "D") +
              theme_bw()

g<-arrangeGrob(p1, p2, p3, p4, ncol=2)
ggsave(file="/Users/ktaylora/Desktop/rf_glm_high_vs_moderate_suit_w_bias_map.png",height=222, width=180, unit="mm", dpi=800,g)

s <- function(x) spTransform(x,CRS(projection("+init=epsg:2163")))

# calculate surface areas of envelopes

surface_areas <- data.frame(
glm=
c((gArea(rgeos::gSymdifference(s(p_glm_absence_bias_90),s(p_glm_absence_bias_10)))/gArea(s(p_glm_current_tss_threshold))),
  (gArea(rgeos::gSymdifference(s(p_glm_presence_bias_90),s(p_glm_presence_bias_10)))/gArea(s(p_glm_current_tss_threshold))),
  (gArea(rgeos::gSymdifference(s(p_glm_climate_bias_90),s(p_glm_climate_bias_10)))/gArea(s(p_glm_current_tss_threshold)))),
rf=
c((gArea(rgeos::gSymdifference(s(p_rf_absence_bias_90),s(p_rf_absence_bias_10)))/gArea(s(p_rf_current_tss_threshold))),
  (gArea(rgeos::gSymdifference(s(p_rf_presence_bias_90),s(p_rf_presence_bias_10)))/gArea(s(p_rf_current_tss_threshold))),
  (gArea(rgeos::gSymdifference(s(p_rf_climate_bias_90),s(p_rf_climate_bias_10)))/gArea(s(p_rf_current_tss_threshold))))
); row.names(surface_areas) <- c("absence","presence","climate");

surface_areas<-round(surface_areas,2)

o_glm <- rgeos::gIntersection(p_glm_climate_bias_90,p_glm_presence_bias_90)
  o_glm <- rgeos::gIntersection(o_glm@polyobj,p_glm_absence_bias_90)@polyobj
o_rf <- rgeos::gIntersection(p_rf_climate_bias_90,p_rf_presence_bias_90)
  o_rf <- rgeos::gIntersection(o_rf@polyobj,p_rf_absence_bias_90)@polyobj

consensus_areas <- data.frame(
  glm=gArea(s(o_glm)),
  rf=gArea(s(o_rf))
)
rownames(consensus_areas) <- "consensus area"
rm(o_glm,o_rf)

perimeters <- data.frame(
glm_p50=
c((perimeter(p_glm_absence_bias_50)/perimeter(p_current_glm_60))-1,
  (perimeter(p_glm_presence_bias_50)/perimeter(p_current_glm_60))-1,
  (perimeter(p_glm_climate_bias_50)/perimeter(p_current_glm_60))-1),
rf_p50=
c((perimeter(p_rf_absence_bias_50)/perimeter(p_current_rf_60))-1,
  (perimeter(p_rf_presence_bias_50)/perimeter(p_current_rf_60))-1,
  (perimeter(p_rf_climate_bias_50)/perimeter(p_current_rf_60))-1),
glm_p90=
c((perimeter(p_glm_absence_bias_90)/perimeter(p_current_glm_90))-1,
  (perimeter(p_glm_presence_bias_90)/perimeter(p_current_glm_90))-1,
  (perimeter(p_glm_climate_bias_90)/perimeter(p_current_glm_90))-1),
rf_p90=
c((perimeter(p_rf_absence_bias_90)/perimeter(p_current_rf_90))-1,
  (perimeter(p_rf_presence_bias_90)/perimeter(p_current_rf_90))-1,
  (perimeter(p_rf_climate_bias_90)/perimeter(p_current_rf_90))-1)
); row.names(perimeters) <- c("absence","presence","climate");

# make violin plots

par(mfrow=c(2,1))
par(mar=par()$mar/2)

t<-data.frame(glm_abs_exp=sample(glm_abs_record_experimental, size=500, replace=T),
              glm_abs_con=sample(glm_abs_record_control, size=500, replace=T),
              glm_pres_exp=sample(glm_pres_record_experimental, size=500, replace=T),
              glm_pres_cont=sample(glm_pres_record_control, size=500, replace=T),
              glm_cli_exp=sample(glm_climate_record_experimental, size=500, replace=T),
              glm_cli_con=sample(glm_climate_record_control, size=500, replace=T),
              rf_abs_exp=sample(rf_abs_record_experimental, size=500, replace=T),
              rf_abs_con=sample(rf_abs_record_control, size=500, replace=T),
              rf_pres_exp=sample(rf_pres_record_experimental, size=500, replace=T),
              rf_pres_cont=sample(rf_pres_record_control, size=500, replace=T),
              rf_cli_exp=sample(rf_climate_record_experimental, size=500, replace=T),
              rf_cli_con=sample(rf_climate_record_control, size=500, replace=T))

df.m <- reshape2::melt(t, id.vars = NULL)
require(gridExtra)
p1<-ggplot(df.m[grep(df.m$variable,pattern="glm_abs_exp|glm_abs_con"),], aes(x = variable, y = value)) + ylim(0.86,1) +
    scale_x_discrete(breaks=c("glm_abs_exp", "glm_abs_con"),labels=c("GLM abs", "GLM abs (control)")) + geom_violin() +
    xlab("") + ylab("Mean AUC/TSS") + theme_bw() + theme(axis.title.y = element_text(size = rel(0.75)))
p2<-ggplot(df.m[grep(df.m$variable,pattern="rf_abs_exp|rf_abs_con"),], aes(x = variable, y = value)) + ylim(0.86,1) +
    scale_x_discrete(breaks=c("rf_abs_exp", "rf_abs_con"),labels=c("RF abs", "RF abs (control)")) + geom_violin() +
    xlab("") + ylab("") + theme_bw()+ theme(axis.title.y = element_text(size = rel(0.75)))
p3<-ggplot(df.m[grep(df.m$variable,pattern="glm_pres_exp|glm_pres_cont"),], aes(x = variable, y = value)) + geom_violin() +
    scale_x_discrete(breaks=c("glm_pres_exp", "glm_pres_cont"),labels=c("GLM pres", "GLM pres (control)")) + geom_violin() +
    xlab("") + ylab("Mean AUC/TSS") + theme_bw()+ theme(axis.title.y = element_text(size = rel(0.75)))
p4<-ggplot(df.m[grep(df.m$variable,pattern="rf_pres_exp|rf_pres_con"),], aes(x = variable, y = value)) + geom_violin() +
    scale_x_discrete(breaks=c("rf_pres_exp", "rf_pres_cont"),labels=c("RF pres", "RF pres (control)")) + geom_violin() +
    xlab("") + ylab("") + theme_bw()+ theme(axis.title.y = element_text(size = rel(0.75)))
p5<-ggplot(df.m[grep(df.m$variable,pattern="glm_cli_exp|glm_cli_con"),], aes(x = variable, y = value)) + ylim(0.85,1) + geom_violin() +
    scale_x_discrete(breaks=c("glm_cli_exp", "glm_cli_con"),labels=c("GLM climate", "GLM climate (control)")) + geom_violin() +
    xlab("") + ylab("Mean AUC/TSS") + theme_bw()+ theme(axis.title.y = element_text(size = rel(0.75)))
p6<-ggplot(df.m[grep(df.m$variable,pattern="rf_cli_exp|rf_cli_con"),], aes(x = variable, y = value)) + ylim(0.85,1) + geom_violin() +
    scale_x_discrete(breaks=c("rf_cli_exp", "rf_cli_con"),labels=c("RF climate", "RF climate (control)")) + geom_violin() +
    xlab("") + ylab("") + theme_bw()+ theme(axis.title.y = element_text(size = rel(0.75)))
grid.arrange(p1, p2, p3,p4,p5,p6,ncol=2,nrow=3)


#plot(difference_RCP85_glm_2070, ylim=c(30,55), xlim=c(-130,-90),col=COLORS(100), main="(GLM) RCP 8.5 (2070)",cex=0.8);
#plot(b, add=T, border="#34343460")

#
# Boxplots
#

require(ggplot2)
require(gridExtra)
setwd(paste(HOME,"/Products/uw/big_sagebrush_sensitivity_analysis/sdm-evaluation/sdm-evaluation-scores/sdm-evaluation-scores",sep=""))
source("process.R")

names <- ls(pattern="glm.")
  names <- names[grepl(names,pattern="Runs")]
    d_min <- min(unlist(lapply(lapply(as.list(names),FUN=get),FUN=length)))
    d_out <- lapply(lapply(as.list(names),FUN=get),FUN=sample, size=d_min)
t <- data.frame(cont=d_out[which(grepl(names,pattern="AbsBiasControlRuns"))],
                pres=d_out[which(grepl(names,pattern="ccurrenceBiasRuns"))],
                abs=d_out[which(grepl(names,pattern="bsBiasRuns"))],
                climate=d_out[which(grepl(names,pattern="imateBiasRuns"))])
names(t) <- c("control","pres","abs","climate")

p1 <- ggplot(melt(t)) + geom_boxplot(aes(y=value,x=variable),outlier.shape = NA,notch=F) +
                        xlab("") + ylab("TSS") + ylim(c(0.8,1)) + theme_bw() +
                        annotate("text", y = 0.995, x=0.75, label = "A")

names <- ls(pattern="rf.")
  names <- names[grepl(names,pattern="Runs")]
    d_min <- min(unlist(lapply(lapply(as.list(names),FUN=get),FUN=length)))
    d_out <- lapply(lapply(as.list(names),FUN=get),FUN=sample, size=d_min)
t <- data.frame(cont=d_out[which(grepl(names,pattern="AbsBiasControlRuns"))],
                pres=d_out[which(grepl(names,pattern="ccurrenceBiasRuns"))],
                abs=d_out[which(grepl(names,pattern="bsBiasRuns"))],
                climate=d_out[which(grepl(names,pattern="imateBiasRuns"))])

names(t) <- c("control","pres","abs","climate")

p2 <- ggplot(melt(t)) + geom_boxplot(aes(y=value,x=variable),outlier.shape = NA, notch=F) +
                        xlab("") + ylab("TSS") + ylim(c(0.8,1)) + theme_bw() +
                        annotate("text", y = 0.995, x=0.75, label = "B")

grid.arrange(p1, p2, ncol=1)
