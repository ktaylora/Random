#
# Build non-parametric models of species occurrence/non-occurrence from USGS range data, landcover data, and more
#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2015]
#
#

require(rgdal)
require(raster,quietly=TRUE)
require(randomForest,quietly=TRUE)

HOME <- Sys.getenv("HOME")

require(landscapeAnalysis)

argv <- commandArgs(trailingOnly=T)

# check runtime arguments

if(!is.na(argv[1])) { 
	sessions <- paste(HOME,"/",argv[1],".presence_absence_surfaces.rdata",sep="") 
} else {
  sessions <- list.files(path=HOME,pattern="pres_abs_pts_usgs_species_viewer",full.names=T)
}

# read in our regional boundaries, species data, and land cover data

r                         <- raster(paste(HOME,"/PLJV/landcover/orig/Final_LC_8bit.tif",sep=""))                                                             # PLJV Composite land cover
buffered_gplcc_boundary   <- readOGR(paste(HOME,"/PLJV/boundaries/GPLCC Pilot/GPLCC_Pilot_Region/",sep=""),"GPLCC_pilot_region_boundary_buffered",verbose=F) # GPLCC Pilot Region Boundary

spp <- readOGR(paste(HOME,"/PLJV/species_data/usgs_gap_species_ranges/vector/",sep=""),argv[1], verbose=F)
  spp <- spTransform(spp,CRS(projection(buffered_gplcc_boundary)))
    spp <- spp[as.vector(!is.na(sp::over(spp,buffered_gplcc_boundary))),] 

presences <- SpatialPoints(spp[spp$resp==1,])
  projection(presences) <- projection(spp)
absences <- SpatialPoints(spp[spp$resp==0,])
  projection(absences) <- projection(spp)

# define local functions

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

# define our 30m land cover types at the association level

grass_habitat <- c(39, # CRP
                   37, # Pasture
                   71, # Mixedgrass
                   75) # Shortgrass

agricultural_habitat <- c(38,202,201,203,204,205,206,207,208,209,210) # Cropland (Wheat, Corn, Soybeans, et al.)

shrubland_habitat <- c(82, # Shrubland
                       83, # Mesquite
                       81, # Savannah
                       85, # Shinnery
                       87) # Sand sage
                       
# read-in our presence-absence and explanatory data
# there are a lot of variables here -- for our first round of EDA, we are going to use them all across a very large forest and
# let the presence/absence data define which are the most informative.  To limit complexity, we will select the most important / biologically
# meaningful variables from each group and use them in the final model.
#

# SPATIAL / TOPOGRAPHIC VARIABLES
topo_roughness_3x3        <- raster(paste(HOME,"/PLJV/DEM/rough3.tif",sep=""))                                                                               # Topographic roughness (3x3)
topo_roughness_27x27      <- raster(paste(HOME,"/PLJV/DEM/rough27.tif",sep=""))                                                                              # Topographic roughness (27x27)
# CLIMATE VARIABLES
bioclim                   <- raster::stack(list.files(paste(HOME,"/Products/weather/worldclim/30_sec/bioclim/",sep=""),full.names=T,pattern="bil$"))                # 19 BIOCLIM variables
# HUMAN DIMENSIONS VARIABLES
population_density        <- raster(paste(HOME,"/PLJV/human_dimensions/usa_gpwv3_pdens_bil_25/usads00ag.bil",sep=""),crs=CRS(projection("+init=epsg:4326")))
distance_to_transmission  <- raster(paste(HOME,"/PLJV/infrastructure/products/bcr_18_19_distanceToTransmissionLines.tif",sep=""))
distance_to_roads         <- raster()
# SOIL TEXTURE VARIABLES 
sand_fraction_0_5         <- raster()
sand_fraction_5_10        <- raster()
sand_fraction_15_20       <- raster()
silt_fraction_0_5         <- raster()
silt_fraction_5_10        <- raster()
silt_fraction_15_20       <- raster()
clay_fraction_0_5         <- raster()
clay_fraction_5_10        <- raster()
clay_fraction_15_20       <- raster()

#soil_texture              <- raster::stack(list.files(paste(HOME,"/Products/weather/worldclim/30_sec/bioclim/",full.names=T,pattern="bil$")))

## PROCESS INPUT DATA REDUNDANTLY, CHECKING RDATA FILES FOR OUR PREVIOUS CONTENT FOR FOCAL SPECIES BEFORE REPEATING EXTRACTIONS/ANALYSES 

# crop 700m buffer regions around our presence/absence sample points and sample land cover and topography
for(s in sessions){
  if(!file.exists(s)){
    save.image(file=s,compress=T)
  }
  # don't let a workspace overwrite the sessions or s variables 
  SESSIONS <- sessions; S<-s; load(s); 
  sessions <- SESSIONS; s<-S; rm(S,SESSIONS);
  if(length(ls(pattern="absences_topography")) == 0){
     cat(" -- sampling landcover for focal species\n");
     #parse those points that are within some buffered distance of the GPLCC pilot region
     if(is.na(argv[1])){ 
       n <- unlist(strsplit(unlist(strsplit(s,split="_"))[1],split="/"));
         n <- n[length(n)]
           argv[1] <- paste(n,"_pres_abs_pts_usgs_species_viewer",sep="")
     }
                
    t_3   <- extract(topo_roughness_3x3,presences)
    t_27  <- extract(topo_roughness_27x27,presences)
                
    presences_topography <- data.frame(rough3=t_3,rough27=t_27)
		  
    t_3   <- extract(topo_roughness_3x3,absences)
    t_27  <- extract(topo_roughness_27x27,absences)
     
    absences_topography <- data.frame(rough3=t_3,rough27=t_27)
                
	  rm(t_3,t_27);
	  save.image(file=s,compress=T)
  }; 
}
# extract climate data
for(s in sessions){
  session_data <- new.env(); load(s,envir=session_data);
  if(length(ls(pattern="bioclim_absences")) == 0){  
    assign(x="bioclim_presences", value=data.frame(extract(bioclim,presences)), envir=session_data)
    assign(x="bioclim_absences", value=data.frame(extract(bioclim,absences)), envir=session_data)
    save(list=ls(session_data), envir=session_data, file=s, compress=T);
  }
# classify our landcover into various cover types
  session_data        <- new.env(); load(s,envir=session_data);
  presence_landcover  <- subsampleSurface(r,pts=presences, width=750)
  absence_landcover   <- subsampleSurface(r,pts=absences, width=750)

  if(length(ls(session_data,pattern="sh_absences"))==0){
	  cat(" -- reclassing presence/absence landcover for:",s,"\n");
	  if(length(ls(session_data,pattern="grass_presences"))==0){ # assumes if presences aren't in session_data, that absences aren't either.
	    assign(x="grass_presences", value=lReclass(presence_landcover, inValues=grass_habitat), envir=session_data)
	    assign(x="grass_absences", value=lReclass(absence_landcover, inValues=grass_habitat), envir=session_data)
	  }
	  if(length(ls(session_data,pattern="ag_presences"))==0){
	    assign(x="ag_presences", value=lReclass(presence_landcover, inValues=agricultural_habitat), envir=session_data)
	    assign(x="ag_absences", value=lReclass(absence_landcover, inValues=agricultural_habitat), envir=session_data)
	  }
	  if(length(ls(session_data,pattern="sh_presences"))==0){
	    assign(x="sh_presences", value=lReclass(presence_landcover, inValues=shrubland_habitat), envir=session_data)
	    assign(x="sh_absences", value=lReclass(absence_landcover, inValues=shrubland_habitat), envir=session_data)
	  }
	  save(list=ls(session_data), envir=session_data, file=s, compress=T);
  }    
# extract human dimensions variables
for(s in sessions){
  session_data <- new.env(); load(s,envir=session_data);
  if(length(ls(pattern="population_density_absences")) == 0){  
    assign(x="popDensity_presences", value=data.frame(popDensity=extract(population_density,presences)), envir=session_data)
    assign(x="popDensity_absences", value=data.frame(popDensity=extract(population_density,absences)), envir=session_data)
    assign(x="distTransmission_presences", value=data.frame(distTrans=extract(distance_to_transmission,presences)), envir=session_data)
    assign(x="distTransmission_absences", value=data.frame(distTrans=extract(distance_to_transmission,absences)), envir=session_data)
    assign(x="distRoads_presences", value=data.frame(distRoads=extract(distance_to_roads,presences)), envir=session_data)
    assign(x="distRoads_absences", value=data.frame(distRoads=extract(distance_to_roads,absences)), envir=session_data)
    save(list=ls(session_data), envir=session_data, file=s, compress=T);
  }
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
	  
	t <-cbind(t,rbind(get('presences_topography', envir=session_data), get('absences_topography', envir=session_data)))
    t <- cbind(t,rbind(get('bioclim_presences',envir=session_data), get('bioclim_absences', envir=session_data)))
	assign(x="t", value=t, envir=session_data);
	save(list=ls(session_data), envir=session_data, file=s, compress=T);
  }

# Build our initial RF models
  session_data <- new.env(); load(s,envir=session_data);
  if(length(ls(session_data,pattern="m_rf"))==0){
		cat(" -- building RF model for:",s,"\n")
		t<-t[,names(t)!="id"]
		t[is.na(t)] <- 0
	  assign(x="m_rf", value=randomForest(as.factor(resp)~., data=t, ntree=1000, do.trace=T,importance=T), envir=session_data); # save the full model with diagnostics so we can evaluate it / predict with it
	  save(list=ls(session_data), envir=session_data, file=s, compress=T);
  }

# Apply extracted to GPLCC pilot region conservation delivery units 
for(s in sessions){
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
