require(rgdal)
require(raster,quietly=TRUE)
require(randomForest,quietly=TRUE)

HOME <- Sys.getenv("HOME")

require(landscapeAnalysis)

argv <- "mass_pres_abs_pts_usgs_species_viewer"

# check runtime arguments

if(!is.na(argv[1])) {
	sessions <- paste(HOME,"/",argv[1],".presence_absence_surfaces.rdata",sep="")
} else {
  sessions <- list.files(path=HOME,pattern="pres_abs_pts_usgs_species_viewer",full.names=T)
}

# read in our regional boundaries, species data, and land cover data

r                         <- raster(paste(HOME,"/PLJV/landcover/orig/Final_LC_8bit.tif",sep=""))                                                             # PLJV Composite land cover
buffered_gplcc_boundary   <- readOGR(paste(HOME,"/PLJV/boundaries/GPLCC Pilot/GPLCC_Pilot_Region/",sep=""),"GPLCC_pilot_region_boundary_buffered",verbose=F) # GPLCC Pilot Region Boundary

# CLIMATE VARIABLES
bioclim                   <- raster::stack(list.files(paste(HOME,"/Products/weather/worldclim/30_sec/bioclim/",sep=""),full.names=T,pattern="bio_18[.]bil$"))

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

## PROCESS INPUT DATA REDUNDANTLY, CHECKING RDATA FILES FOR OUR PREVIOUS CONTENT FOR FOCAL SPECIES BEFORE REPEATING EXTRACTIONS/ANALYSES

# crop 700m buffer regions around our presence/absence sample points and sample land cover and topography
s=sessions;
if(!file.exists(s)){
 save.image(file=s,compress=T)
}
# don't let a workspace overwrite the sessions or s variables
SESSIONS <- sessions; S<-s; load(s);
sessions <- SESSIONS; s<-S; rm(S,SESSIONS);
# extract climate data
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

if(length(ls(session_data,pattern="ag_absences"))==0){
  cat(" -- reclassing presence/absence landcover for:",s,"\n");
  if(length(ls(session_data,pattern="grass_presences"))==0){ # assumes if presences aren't in session_data, that absences aren't either.
    assign(x="shrub_presences", value=lReclass(presence_landcover, inValues=shrubland_habitat), envir=session_data)
    assign(x="shrub_absences", value=lReclass(absence_landcover, inValues=shrubland_habitat), envir=session_data)
  }
  if(length(ls(session_data,pattern="ag_presences"))==0){
    assign(x="ag_presences", value=lReclass(presence_landcover, inValues=agricultural_habitat), envir=session_data)
    assign(x="ag_absences", value=lReclass(absence_landcover, inValues=agricultural_habitat), envir=session_data)
  }
  save(list=ls(session_data), envir=session_data, file=s, compress=T);
}
# Extract topographic variables
if(length(ls(pattern="absences_topography")) == 0){
     cat(" -- sampling landcover for focal species\n");
     #parse those points that are within some buffered distance of the GPLCC pilot region
     if(is.na(argv[1])){ 
       n <- unlist(strsplit(unlist(strsplit(s,split="_"))[1],split="/"));
         n <- n[length(n)]
           argv[1] <- paste(n,"_pres_abs_pts_usgs_species_viewer",sep="")
     }
                
    t_27  <- extract(topo_roughness_27x27,presences)
                
    presences_topography <- data.frame(rough27=t_27)
		  
    t_27  <- extract(topo_roughness_27x27,absences)
     
    absences_topography <- data.frame(rough27=t_27)
                
    rm(t_27);
    save.image(file=s,compress=T)
}; 
# Calculate human dimensions variables
assign(x="distRoads_presences", value=data.frame(distRoads=extract(distance_to_roads,presences)), envir=session_data)
assign(x="distRoads_absences", value=data.frame(distRoads=extract(distance_to_roads,absences)), envir=session_data)
save(list=ls(session_data), envir=session_data, file=s, compress=T);

# Calculate and aggregate our landscape metrics and other explanatory variables into a data.frame
 session_data <- new.env(); load(s,envir=session_data);
 if(sum(ls(session_data)=="t")==0){
	cat(" -- calculating landscape metrics for:",s,"\n");
        t_landcover_shrub_presences <- cbind(resp=1,calculateLMetricsForFocalCover(get("shrub_presences",session_data), metrics=c("total.area")));
	t_landcover_shrub_absences  <- cbind(resp=0,calculateLMetricsForFocalCover(get("shrubs_absences",session_data), metrics=c("total.area")));
	t_landcover_ag_presences    <- cbind(resp=1,calculateLMetricsForFocalCover(get("ag_presences",session_data), metrics=c("total.area")));
	t_landcover_ag_absences     <- cbind(resp=0,calculateLMetricsForFocalCover(get("ag_absences",session_data), metrics=c("total.area")));

	t_1 <- cbind(t_landcover_shrub_presences[,1:2], t_landcover_ag_presences[,3])
        t_2 <- cbind(t_landcover_shrub_absences[,1:2], t_landcover_ag_absences[,3])

	 t <- rbind(t_1,t_2);

	names(t) <-
	c(names(t_landcover_grass_presences)[1:2],paste("grass_",names(t_landcover_grass_presences)[3],sep=""),
		                                        paste("ag_",names(t_landcover_ag_presences)[3],sep=""));

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
