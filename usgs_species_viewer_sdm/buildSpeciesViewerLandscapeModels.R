#
# Quick and Dirty non-parametric models of species occurrence/non-occurrence from USGS range data and grassland habitat data
#

require(rgdal)
require(raster,quietly=TRUE)
require(randomForest,quietly=TRUE)
require(landscapeAnalysis)

HOME <- Sys.getenv("HOME")

argv <- commandArgs(trailingOnly=T)

# local functions
calculateLMetricsForFocalCover <- function(sampleUnits=NULL,metrics=NULL){
  t <- data.frame(matrix(nrow=length(sampleUnits), ncol=length(metrics)+1))
  # iterate over units, sampling as we go # iterate over our samples, calculating statistics as we go
  for(i in 1:length(sampleUnits)){
    row       <- rep(NA,length(metrics));
    lm_output <- lCalculateLandscapeMetrics(x=list(sampleUnits[[i]]), metric=metrics, DEBUG=F)
    for(j in 1:length(metrics)){ row[j] <- median(metricsListToVector(lm_output, metrics[j])) }
  	  row <- c(i,row);
    # contribute to table
    t[i,] <- row;
  }; cat("\n");
  names(t) <- c("id",metrics)
  return(t);
}

# check arguments

if(!is.na(argv[1])) {
	sessions <- paste(HOME,"/",argv[1],".presence_absence_surfaces.rdata",sep="")
} else {
  sessions <- list.files(path=HOME,pattern="pres_abs_pts_usgs_species_viewer",full.names=T)
}

# define our cover types
grass_habitat <- c(39, # CRP
                   37, # Pasture
                   71, # Mixedgrass
                   75, # Shortgrass
                   85, # Shinnery
                   87) # Sand sage

agricultural_habitat <- c(38,202,201,203,204,205,206,207,208,209,210) # Cropland

shrubland_habitat <- c(82, # Shrubland
                       83, # Mesquite
                       81) # Savannah

# read-in our source explanatory data
r                         <- raster(paste(HOME,"/PLJV/landcover/orig/Final_LC_8bit.tif",sep=""))
buffered_gplcc_boundary   <- readOGR(paste(HOME,"/PLJV/boundaries/GPLCC Pilot/GPLCC_Pilot_Region/",sep=""),"GPLCC_pilot_region_boundary_buffered",verbose=F)
topo_roughness_3x3        <- raster(paste(HOME,"/PLJV/DEM/rough3.tif",sep=""))
topo_roughness_27x27      <- raster(paste(HOME,"/PLJV/DEM/rough27.tif",sep=""))
mean_t_coldest_quarter    <- raster(paste(HOME,"/Products/weather/worldclim/30_sec/bioclim/bio_11.bil",sep=""))
prec_driest_quarter       <- raster(paste(HOME,"/Products/weather/worldclim/30_sec/bioclim/bio_17.bil",sep=""))

# crop 700m buffer regions around our presence/absence sample points
for(s in sessions){
  if(!file.exists(s)){
    save.image(file=s,compress=T)
  }
  # don't let a workspace overwrite the sessions or s variables
  SESSIONS <- sessions; S<-s; load(s);
  sessions <- SESSIONS; s<-S; rm(S,SESSIONS);
  if(length(ls(pattern="absences")) == 0){
     cat(" -- sampling landcover for focal species\n");
     #parse those points that are within some buffered distance of the GPLCC pilot region
     if(is.na(argv[1])){
       n <- unlist(strsplit(unlist(strsplit(s,split="_"))[1],split="/"));
         n <- n[length(n)]
           argv[1] <- paste(n,"_pres_abs_pts_usgs_species_viewer",sep="")
     }
     spp <- readOGR(paste(HOME,"/PLJV/species_data/usgs_gap_species_ranges/vector/",sep=""),argv[1], verbose=F)
	     spp <- spTransform(spp,CRS(projection(buffered_gplcc_boundary)))
	       spp <- spp[as.vector(!is.na(sp::over(spp,buffered_gplcc_boundary))),]

      presences <- SpatialPoints(spp[spp$resp==1,])
        projection(presences) <- projection(spp)

      t_3   <- extract(topo_roughness_3x3,presences)
      t_27  <- extract(topo_roughness_27x27,presences)
      bio11 <- extract(mean_t_coldest_quarter,presences)
      bio17 <- extract(prec_driest_quarter,presences)

      presences_climate_topography <- data.frame(rough3=t_3,rough27=t_27,bio11=bio11,bio17=bio17)

      absences <- SpatialPoints(spp[spp$resp==0,])
	      projection(absences) <- projection(spp)

      t_3   <- extract(topo_roughness_3x3,absences)
      t_27  <- extract(topo_roughness_27x27,absences)
      bio11 <- extract(mean_t_coldest_quarter,absences)
      bio17 <- extract(prec_driest_quarter,absences)

      absences_climate_topography <- data.frame(rough3=t_3,rough27=t_27,bio11=bio11,bio17=bio17)

		presences <- subsampleSurface(r,pts=presences, width=750)
		absences <- subsampleSurface(r,pts=absences, width=750)

		rm(t_3,t_27,bio11,bio17);
		save.image(file=s,compress=T)
	};

# classify our landcover into various cover types
  session_data <- new.env(); load(s,envir=session_data);
  if(length(ls(session_data,pattern="sh_absences"))==0){
	  cat(" -- reclassing presence/absence landcover for:",s,"\n");
	  if(length(ls(session_data,pattern="grass_presences"))==0){ # assumes if presences aren't in session_data, that absences aren't either.
	    assign(x="grass_presences", value=lReclass(get("presences",session_data), inValues=grass_habitat), envir=session_data)
	    assign(x="grass_absences", value=lReclass(get("absences",session_data), inValues=grass_habitat), envir=session_data)
	  }
	  if(length(ls(session_data,pattern="ag_presences"))==0){
	    assign(x="ag_presences", value=lReclass(get("presences",session_data), inValues=agricultural_habitat), envir=session_data)
	    assign(x="ag_absences", value=lReclass(get("absences",session_data), inValues=agricultural_habitat), envir=session_data)
	  }
	  if(length(ls(session_data,pattern="sh_presences"))==0){
	    assign(x="sh_presences", value=lReclass(get("presences",session_data), inValues=shrubland_habitat), envir=session_data)
	    assign(x="sh_absences", value=lReclass(get("absences",session_data), inValues=shrubland_habitat), envir=session_data)
	  }
	  save(list=ls(session_data), envir=session_data, file=s, compress=T);
	}

# Calculate and aggregate our landscape metrics and other explanatory variables into a data.frame
	session_data <- new.env(); load(s,envir=session_data);
	if(sum(ls(session_data)=="t")==0){
		cat(" -- calculating landscape metrics for:",s,"\n");
	  t_landcover_grass_presences <- cbind(resp=1,calculateLMetricsForFocalCover(get("grass_presences",session_data), metrics=c("total.area","mean.patch.area","shape.index")));
	  t_landcover_grass_absences  <- cbind(resp=0,calculateLMetricsForFocalCover(get("grass_absences",session_data), metrics=c("total.area","mean.patch.area","shape.index")));
	  t_landcover_ag_presences    <- cbind(resp=1,calculateLMetricsForFocalCover(get("ag_presences",session_data), metrics=c("total.area","mean.patch.area","shape.index")));
	  t_landcover_ag_absences     <- cbind(resp=0,calculateLMetricsForFocalCover(get("ag_absences",session_data), metrics=c("total.area","mean.patch.area","shape.index")));
	  t_landcover_sh_presences    <- cbind(resp=1,calculateLMetricsForFocalCover(get("sh_presences",session_data), metrics=c("total.area","mean.patch.area","shape.index")));
	  t_landcover_sh_absences     <- cbind(resp=0,calculateLMetricsForFocalCover(get("sh_absences",session_data), metrics=c("total.area","mean.patch.area","shape.index")));

	  t_1 <- cbind(t_landcover_grass_presences[,1:5], t_landcover_ag_presences[,3:5], t_landcover_sh_presences[,3:5])
	  t_2 <- cbind(t_landcover_grass_absences[,1:5], t_landcover_ag_absences[,3:5], t_landcover_sh_absences[,3:5])

	  t<-rbind(t_1,t_2);

	  names(t) <-
	  c(names(t_landcover_grass_presences)[1:2],paste("grass_",names(t_landcover_grass_presences)[3:5],sep=""),
	  	                                        paste("ag_",names(t_landcover_ag_presences)[3:5],sep=""),
	  	                                        paste("sh_",names(t_landcover_sh_presences)[3:5],sep=""));

	  t<-cbind(t,rbind(presences_climate_topography, absences_climate_topography))
	  assign(x="t", value=t, envir=session_data);
	  save(list=ls(session_data), envir=session_data, file=s, compress=T);
	}

# Build our initial RF models
  session_data <- new.env(); load(s,envir=session_data);
  if(length(ls(session_data,pattern="m_rf"))==0){
		cat(" -- building RF model for:",s,"\n")
		t<-t[,names(t)!="id"]
		t[is.na(t)] <- 0
	  assign(x="m_rf", value=randomForest(as.factor(resp)~., data=t, ntree=1000, do.trace=T,importance=T), envir=session_data);
	  save(list=ls(session_data), envir=session_data, file=s, compress=T);
  }

# Apply to GPLCC landscape
  load(s)
  units <- readOGR(paste(HOME,"/PLJV/landscape_metrics/products/",sep=""),"gplcc_landscape_current_4326",verbose=F)
  t_units <- units@data
    t_units[is.na(t_units)] <- 0
      t_units<-t_units[,names(t_units)!='id']
  # re-name our columns so they are consistent with what the RF models were trained with
  names(t_units) <- c("ag_total.area","ag_mean.patch.area","ag_shape.index", "grass_total.area","grass_mean.patch.area","grass_shape.index","sh_total.area","sh_mean.patch.area","sh_shape.index","bio11","bio17","rough27","rough3")
  # predict
  units@data <- data.frame(predicted=as.numeric(as.character(predict(m_rf, newdata=t_units, type='prob')[,2])))
  # figure out name and write to CWD
  n <- unlist(strsplit(unlist(strsplit(s,split="_"))[1],split="/"));
    n <- n[length(n)]
  writeOGR(units, ".",paste(n,"_predicted_gplcc_units",sep=""), driver="ESRI Shapefile", overwrite=T)
}
