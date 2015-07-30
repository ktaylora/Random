#
# Big sagebrush subspecies SDMs
# A custom implementation of random forest and GLMs, along with spatial aggregation of records
#
# Author: Kyle Taylor (kyle.taylor@uwyo.edu)
#

require(raster)
require(spatstat)
require(rgdal)
require(rgeos)

require(landscapeAnalysis)

HOME <- Sys.getenv("HOME")


#
#  responsePlot()
#

responsePlot <- function(x,var=NULL,plot=T){
  m <- x[[1]][[1]]
  if(is.null(var)) stop("var= argument undefined");
  if(sum(names(m$data) %in% var)==0) stop("var= not found in model object")
  # define names
  names <- names(m$data);
    names <- names[names != "resp"]

  d        <- x[[1]][[1]]$data
  focal    <- d[,var]
  notFocal <- d[,names[names!=var]]

  x <- apply(notFocal,2,FUN=median, na.rm=T)
    x <- data.frame(matrix(x,nrow=1))
      names(x) <- names[names!=var]

  spread <- sort(runif(min=min(focal,na.rm=T),max=max(focal,na.rm=T),n=nrow(notFocal)))
    x <- cbind(spread,x)
      names(x) <- c(var,names[names!=var])

  prob <- as.vector(predict(m,x,type="resp"))
    prob[prob < 0] <- 0; prob[prob>1] <- 1
      prob <- data.frame(cbind(prob,x[,var]))
        names(prob) <- c("prob",var)
  if(plot){
    plot(prob~get(var),type="l",data=prob, ylim=c(0,1),xlab=var,ylab="p(occ)", col="white",cex=1.8)
      grid(lwd=1.2); lines(prob~get(var),lwd=1.2,col="red",data=prob)
  }
}

#
# spatialPointsToPPP()
# Convert a standard SpatialPoints* object to a ppp for use in spatstat
#

spatialPointsToPPP <- function(x,extentMultiplier=1.1){
  # default includes
  require(rgdal)
  require(raster)
  require(spatstat)

  e <- extent(x)

  if(!is.null(extentMultiplier)) {
    e@xmin <- e@xmin*extentMultiplier
    e@xmax <- e@xmax+abs(e@xmax*(extentMultiplier-1))
    e@ymin <- e@ymin-abs(e@ymin*(extentMultiplier-1))
    e@ymax <- e@ymax*extentMultiplier
  }

  if(class(x) == "SpatialPointsDataFrame"){
    x <- x@coords
    x <- ppp(x=x[,1], y=x[,2], window=owin(xrange=c(e@xmin,e@xmax), yrange=c(e@ymin,e@ymax)))
  }

  return(x)
}

#
# buildTrainingEvaluationSets()
# split a spatial points data frame into training or evaluation datasets based on spatial segregation of data
#

buildTrainingEvaluationSets <- function(p_focal, type="spatially-uniform", debug=T){
 # generate pseudo-absences using an 8-degrees from presences method
  r_template <- raster(res=0.008333333) # consistent with the CRS and resolution of our climate data
  p_focal    <- spTransform(p_focal, CRS(projection(r_template)))
  r_focal    <- crop(r_template, extent(p_focal)*2)
    cat(" -- rasterizing species presence points to a grid that is consistent with our climate data\n")
      r_focal <- rasterize(p_focal[p_focal$resp ==1,],field='resp', fun=min, r_focal,background=0);
  # generate a sample of potential absence points (i.e., anything that isn't in a presence cell)
  abs_pts <- sampleRandom(r_focal,size=ceiling(20*nrow(p_focal[p_focal$resp==1,])),sp=T) # let's make a large pool of absence records to start with
    abs_pts <- abs_pts[abs_pts$layer==0,]
  # throw out absences whose minimum distance is > 8 degrees from a presence record
  cat(" -- calculating distances between absence points and presence points\n")
  boundaries <- spTransform(boundaries,CRS(projection(abs_pts)))
  t_l <- rgeos::gDistance(abs_pts,p_focal[p_focal$resp==1,],byid=T)
    t_l <- as.vector(apply(t_l,2,FUN=min) < 8) # are we 8 degrees or less from a presence record?
      abs_pts <- abs_pts[t_l,]
        abs_pts <- rasterToPoints(rasterize(abs_pts,field='layer', fun=min, r_focal),sp=T) # re-grid our absence points to the resolution of our climate data so there are no duplicates
          abs_pts <- abs_pts[as.vector(!is.na(sp::over(abs_pts,boundaries))),]
            abs_pts <- abs_pts[sample(1:nrow(abs_pts),size=nrow(p_focal[p_focal$resp==1,])),] # downsample to a count consistent with our presence records
  # generate a sampling grid for k-fold crossvalidation from our presence data

  ## implement a spatially uniform split?
  if(grepl(type,pattern="uniform")){
    cat(" -- performing spatial uniform sampling\n")
    x <- spatialPointsToPPP(p_focal[p_focal$resp==1,])
      x <- envelope(x, r=seq(0,0.8,0.001), fun=Jest, 1000)
    cat(" -- values for r (degrees) that intersect with the same level of clustering observed in presence records:\n");
    print(seq(0,0.8,0.001)[which(x$lo < x$obs)])
    intersection <- seq(0,0.8,0.001)[which(x$lo < x$obs)][5] # treat the first 5 values as burn-in and ignore them
    if(debug){
      dev.new(height=6, width=8)
        plot(x,col="white", main=paste("J-Function For ",deparse(substitute(p_focal)),sep=""))
          grid();
            plot(x, add=T)
              abline(v=intersection)
    }
    # build a rasterized checkerboard sampling grid
    samplingGrid <- crop(raster(res=intersection),boundaries)
      samplingMatrix <- matrix(c(0,1), nrow=nrow(samplingGrid), ncol=ifelse(!(ncol(samplingGrid)%%2),ncol(samplingGrid)+1,ncol(samplingGrid)), byrow=T)
        samplingGrid <- setValues(samplingGrid,samplingMatrix[1:nrow(samplingGrid),1:ncol(samplingGrid)])
    # plot our records space
    if(debug){
      dev.new(height=6,width=8)
      plot(samplingGrid, add=F, legend=F)
      plot(abs_pts,pch=15, cex=0.5, col="red", add=T)
      plot(p_focal[p_focal$resp==1,], col="DarkBlue", cex=0.5, pch=15,add=T)
      plot(boundaries, add=T)
    }
    # extract training and evaluation data
    climate_variables <- crop(climate_variables,boundaries,progress='text')
    cat(" -- extracting climate data for model training\n")
    training <- raster::extract(samplingGrid,p_focal[p_focal$resp==1,],sp=T)
      training_presence     <- training[training$layer == 1,]
        training_presence@data <- data.frame(resp=rep(1,nrow(training_presence@data)))
          training_presence@data <- cbind(training_presence@data,extract(climate_variables,training_presence))
      evaluation_presence  <- training[training$layer == 0,]
        evaluation_presence@data <- data.frame(resp=rep(1,nrow(evaluation_presence@data)))
          evaluation_presence@data <- cbind(evaluation_presence@data,extract(climate_variables,evaluation_presence))

    training <- raster::extract(samplingGrid,abs_pts,sp=T)
      training_abs     <- training[as.vector(training[,2]@data == 1),]
        training_abs@data <- data.frame(resp=rep(0,nrow(training_abs@data)))
          training_abs@data <- cbind(training_abs@data,extract(climate_variables,training_abs))
      evaluation_abs   <- training[as.vector(training[,2]@data == 0),]
        evaluation_abs@data <- data.frame(resp=rep(0,nrow(evaluation_abs@data)))
          evaluation_abs@data <- cbind(evaluation_abs@data,extract(climate_variables,evaluation_abs))

    training   <- rbind(training_presence,training_abs)
    evaluation <- rbind(evaluation_presence,evaluation_abs)

  ## implement the longitudinal strips of Bahn, V., and B. J. McGill. 2013?
  } else if(grepl(type,pattern="longitudinal")){
    cat(" -- performing longitudinal strip sampling\n")
    # identify the longitudinal quantiles or our presence records
    lon_quantiles <- as.vector(quantile(p_focal[p_focal$resp==1,]@coords[,1],p=c(0.25,0.5,0.75)))
    presence <- p_focal[p_focal$resp == 1,];
       training_presence <- presence[presence@coords[,1] < lon_quantiles[1],] # 1st quantile
      training_presence <- rbind(training_presence,
                        presence[presence@coords[,1] > lon_quantiles[2] & # 3rd quantile
                                 presence@coords[,1] < lon_quantiles[3],])
    training_abs <- abs_pts[abs_pts@coords[,1] < lon_quantiles[1],] # 1st quantile
      training_abs <- rbind(training_abs,
                            abs_pts[abs_pts@coords[,1] > lon_quantiles[2] & # 3rd quantile
                                    abs_pts@coords[,1] < lon_quantiles[3],])

    evaluation_presence <- presence[presence@coords[,1] > lon_quantiles[1] & # 2nd quantile
                                   presence@coords[,1] < lon_quantiles[2],]
      evaluation_presence <- rbind(evaluation_presence,
                                   presence[presence@coords[,1] > lon_quantiles[3],]) # 4th quantile
    evaluation_abs <- abs_pts[abs_pts@coords[,1] > lon_quantiles[1] & # 2nd quantile
                              abs_pts@coords[,1] < lon_quantiles[2],]
      evaluation_abs <- rbind(evaluation_abs,
                              abs_pts[abs_pts@coords[,1] > lon_quantiles[3],]) # 4th quantile
    # extract training and evaluation data
    climate_variables <- crop(climate_variables,boundaries,progress='text')
    cat(" -- extracting climate data for model training")
    training_presence@data <- data.frame(resp=rep(1,nrow(training_presence@data)))
      training_presence@data <- cbind(training_presence@data,extract(climate_variables,training_presence))
    evaluation_presence@data <- data.frame(resp=rep(1,nrow(evaluation_presence@data)))
      evaluation_presence@data <- cbind(evaluation_presence@data,extract(climate_variables,evaluation_presence))
    training_abs@data <- data.frame(resp=rep(0,nrow(training_abs@data)))
      training_abs@data <- cbind(training_abs@data,extract(climate_variables,training_abs))
    evaluation_abs@data <- data.frame(resp=rep(0,nrow(evaluation_abs@data)))
      evaluation_abs@data <- cbind(evaluation_abs@data,extract(climate_variables,evaluation_abs))

    training   <- rbind(training_presence,training_abs)
    evaluation <- rbind(evaluation_presence,evaluation_abs)
  }

  return(list(training,evaluation))
}


#
# buildSDMS_ssp()
# build random forest and glms for focal sagebrush ssp
#

build_GLM <- function(training,evaluation,formula=NULL,debug=F){
  # build some GLMs and evaluate the best model using hold-out data
  orders   <- expand.grid(rep(list(1:3), 5)) # all possible combinations of GLM polynomial orders (from 1->3) for each of our explanatory variables
  formulas <- list();
  models   <- list();
  cat(" -- building GLM series:")
  for(i in 1:nrow(orders)){
    # build a formula object for this step of the 'orders' matrix
    j <- paste(names(training)[2:6],orders[i,],sep="^")
      j<-paste("I(",j,sep="")
        j<-paste(j,")",sep="")
          j<-paste(j,collapse="+")
            # tack-on a string of all potential variable interactions
            j<-paste(j,paste(apply(apply(combn(names(training)[2:6],m=2),MARGIN=1,FUN=c),1,FUN=paste,collapse=":"),collapse="+"),sep="+")
              formulas[[length(formulas)+1]] <- formula(paste(names(training)[1],j,sep="~"))
    # fit a GLM to our formula object and store in a massively inefficient list
    models[[length(models)+1]] <- glm(formula=formulas[[length(formulas)]],family=binomial,data=training); cat(".");
  };cat("\n");

  # Make Some Q/D GLM AIC plots
  if(debug){
    dev.new(height=6,width=8)
    hist(unlist(lapply(models, FUN=AIC)),breaks=70)
  }

  # identify which model was the best out of the series by minimizing AIC
  best  <- which(unlist(lapply(models, FUN=AIC)) == min(unlist(lapply(models, FUN=AIC))))
    m_glm <- models[best]

  cutoffs <- seq(0.01,0.99,0.01)
  performance <- vector()
  # find a threshold that optimizes overall model performance against the evaluation dataset and use it to find a cut-off for binary image
  for(c in cutoffs){
    test <- as.numeric(round(as.numeric(predict(m_glm[[1]], newdata=evaluation[,2:ncol(evaluation)],type="resp")),3) > c) == evaluation[,1]
    glm_eval <-
      data.frame(
                  overall=sum(test,na.rm=T)/length(evaluation[,1]),
                  sensitivity=sum(test[1:max(which(evaluation$resp == 1))],na.rm=T)/length(test[1:max(which(evaluation$resp == 1))]),
                  specificity=sum(test[max(which(evaluation$resp == 1)):max(which(evaluation$resp == 0))],na.rm=T)/length(test[max(which(evaluation$resp == 1)):max(which(evaluation$resp == 0))])
                )
    performance <- append(performance,glm_eval[,1])
  }
  # run against our optimal cut-off value and report the value
  c=cutoffs[which(performance == max(performance))];
  test <- as.numeric(round(as.numeric(predict(m_glm[[1]], newdata=evaluation[,2:ncol(evaluation)],type="resp")),3) > c) == evaluation[,1]
  glm_eval <-
      data.frame(
                  overall=sum(test,na.rm=T)/length(evaluation[,1]),
                  sensitivity=sum(test[1:max(which(evaluation$resp == 1))],na.rm=T)/length(test[1:max(which(evaluation$resp == 1))]),
                  specificity=sum(test[max(which(evaluation$resp == 1)):max(which(evaluation$resp == 0))],na.rm=T)/length(test[max(which(evaluation$resp == 1)):max(which(evaluation$resp == 0))]),
                  cutoff=c
                )
  return(list(m_glm,glm_eval))
}

#
# buildRF()
#

build_RF <- function(training,evaluation,debug=F){
# build and evaluate a random forest using hold-out data
  require(randomForest)
  m_rf <- randomForest(as.factor(resp)~., data=na.omit(training),importance=T,ntree=2500,do.trace=T)
    test <- as.numeric(as.vector(predict(m_rf,newdata=evaluation[,2:ncol(evaluation)]))) == evaluation[,1]
  rf_eval <-
    data.frame(
                overall=sum(test,na.rm=T)/length(evaluation[,1]),
                sensitivity=sum(test[1:max(which(evaluation$resp == 1))],na.rm=T)/length(test[1:max(which(evaluation$resp == 1))]),
                specificity=sum(test[max(which(evaluation$resp == 1)):max(which(evaluation$resp == 0))],na.rm=T)/length(test[max(which(evaluation$resp == 1)):max(which(evaluation$resp == 0))])
              )
  return(list(m_rf,rf_eval))
}

#
# MAIN
#

boundaries <- readOGR(paste(HOME,"/Products/boundaries/",sep=""), "western_north_american_boundaries",verbose=F)
# read-in and parse our climate variables
climate_variables <- list.files(paste(HOME,"/Products/weather/worldclim/30_sec/bioclim",sep=""),pattern="[.]bil$", full.names=T)
  climate_variables <- raster::stack(climate_variables[grepl(climate_variables, pattern="bio_3|bio_4|bio_11|bio_15|bio_18")])
    climate_variables <- crop(climate_variables,spTransform(boundaries,CRS(projection(climate_variables))),progress='text')
# read-in our previously extracted SpatialPoints shapefiles indicating subspecies presence-absence from GAP
p_vaseyana     <- readOGR(paste(HOME,"/Products/uw/big_sagebrush_subspp_analysis/vectors",sep=""), "vaseyana_gap_records.1",verbose=F)
p_wyomingensis <- readOGR(paste(HOME,"/Products/uw/big_sagebrush_subspp_analysis/vectors",sep=""), "wyomingensis_gap_records",verbose=F)
p_tridentata   <- readOGR(paste(HOME,"/Products/uw/big_sagebrush_subspp_analysis/vectors",sep=""), "tridentata_gap_records",verbose=F)

# build our models for evaluation via longitudinal strips.
o <- buildTrainingEvaluationSets(p_vaseyana,type="longitudinal",debug=T)
   training <- o[[1]]@data; evaluation <- o[[2]]@data; rm(o)

vaseyana_glm_lon <- build_GLM(training,evaluation)
vaseyana_rf_lon  <- build_RF(training,evaluation)

# now build a real model using a not-crazy sampling approach that will actually produce
# meaningful predictions
o <- buildTrainingEvaluationSets(p_vaseyana,type="uniform",debug=F)
  training <- o[[1]]; evaluation <- o[[2]]; rm(o)

vaseyana_glm_unif <- build_GLM(training,evaluation)
vaseyana_rf_unif  <- build_RF(training,evaluation)

cat(" -- projecting model rasters\n")
r_glm_vaseyana_current <- predict(climate_variables, vaseyana_glm_unif[[1]][[1]], type='resp', progress='text')
r_rf_vaseyana_current  <- 1-predict(climate_variables, vaseyana_rf_unif[[1]], type='prob', progress='text')
ens_vaseyana_current   <- stackApply(stack(r_glm_vaseyana_current,r_rf_vaseyana_current),fun=mean,indices=1,progress='text')

o <- buildTrainingEvaluationSets(p_tridentata,type="longitudinal",debug=T)
  training <- o[[1]]@data; evaluation <- o[[2]]@data; rm(o)

tridentata_glm_lon <- build_GLM(training,evaluation)
tridentata_rf_lon  <- build_RF(training,evaluation)

o <- buildTrainingEvaluationSets(p_tridentata,type="uniform",debug=F)
  training <- o[[1]]@data; evaluation <- o[[2]]@data; rm(o)

tridentata_glm_unif <- build_GLM(training,evaluation)
tridentata_rf_unif  <- build_RF(training,evaluation)

cat(" -- projecting model rasters\n")
r_glm_tridentata_current <- predict(climate_variables, tridentata_glm_unif[[1]][[1]], type='resp', progress='text')
r_rf_tridentata_current  <- 1-predict(climate_variables, tridentata_rf_unif[[1]], type='prob', progress='text')
ens_tridentata_current   <- stackApply(stack(r_glm_tridentata_current,r_rf_tridentata_current),fun=mean,indices=1,progress='text')

o <- buildTrainingEvaluationSets(p_wyomingensis,type="longitudinal",debug=T)
  training <- o[[1]]@data; evaluation <- o[[2]]@data; rm(o)

wyomingensis_glm_lon <- build_GLM(training,evaluation)
wyomingensis_rf_lon  <- build_RF(training,evaluation)

o <- buildTrainingEvaluationSets(p_wyomingensis,type="uniform",debug=F)
  training <- o[[1]]@data; evaluation <- o[[2]]@data; rm(o)

wyomingensis_glm_unif <- build_GLM(training,evaluation)
wyomingensis_rf_unif  <- build_RF(training,evaluation)

cat(" -- projecting model rasters\n")
r_glm_wyomingensis_current <- predict(climate_variables, wyomingensis_glm_unif[[1]][[1]], type='resp', progress='text')
r_rf_wyomingensis_current  <- 1-predict(climate_variables, wyomingensis_rf_unif[[1]], type='prob', progress='text')
ens_wyomingensis_current   <- stackApply(stack(r_glm_wyomingensis_current,r_rf_wyomingensis_current),fun=mean,indices=1,progress='text')

dev.new(height=5,width=10)
par(mfrow=c(1,3))

plot(ens_tridentata_current,main="ssp. tridentata");
plot(boundaries, border=rgb(0, 0, 0, 0.5),add=T);
plot(ens_wyomingensis_current,main="ssp. wyomingensis");
plot(boundaries, border=rgb(0, 0, 0, 0.5),add=T);
plot(ens_vaseyana_current,main="ssp. vaseyana");
plot(boundaries, border=rgb(0, 0, 0, 0.5),add=T);

# magnitude of difference surfaces
# r_mag_dif_tridentata<-round(r_rf_tridentata_current/r_glm_tridentata_current,3)
#   r_mag_dif_tridentata[r_mag_dif_tridentata>3] <- 3
#     plot(r_mag_dif_tridentata, main="Tridentata. (Magnitude of difference [RF/GLM])"); plot(boundaries,add=T)
# r_mag_dif_wyomingensis<-round(r_rf_wyomingensis_current/r_glm_wyomingensis_current,3)
#   r_mag_dif_wyomingensis[r_mag_dif_wyomingensis>3] <- 3
#     plot(r_mag_dif_wyomingensis, main="Wyomingensis. (Magnitude of difference [RF/GLM])"); plot(boundaries,add=T)
# r_mag_dif_vaseyana<-round(r_rf_vaseyana_current/r_glm_vaseyana_current,3)
#   r_mag_dif_vaseyana[r_mag_dif_vaseyana>3] <- 3
#     plot(r_mag_dif_vaseyana, main="Vaseyana. (Magnitude of difference [RF/GLM])"); plot(boundaries,add=T)

# Calculate Agreement / Disagreement Vector Surfaces
tridentata_current_quantiles   <- extractDensities(ens_tridentata_current,p=c(0.5,0.75,0.95))
wyomingensis_current_quantiles <- extractDensities(ens_wyomingensis_current,p=c(0.5,0.75,0.95))
vaseyana_current_quantiles     <- extractDensities(ens_vaseyana_current,p=c(0.5,0.75,0.95))

intersect_p50_current <- rgeos::gIntersection(rgeos::gIntersection(tridentata_current_quantiles[[1]],wyomingensis_current_quantiles[[1]])@polyobj,vaseyana_current_quantiles[[1]])@polyobj
intersect_p50_wyo_tri_current <- rgeos::gIntersection(tridentata_current_quantiles[[1]],wyomingensis_current_quantiles[[1]])@polyobj
#intersect_p75_current <- rgeos::gIntersection(rgeos::gIntersection(tridentata_current_quantiles[[2]],wyomingensis_current_quantiles[[2]])@polyobj,vaseyana_current_quantiles[[2]])@polyobj

unique_p50_current_wyomingensis <- rgeos::gDifference(wyomingensis_current_quantiles[[1]],tridentata_current_quantiles[[1]]) # wyomingensis has greater range
  unique_p50_current_wyomingensis <- rgeos::gDifference(unique_p50_current_wyomingensis, vaseyana_current_quantiles[[1]])

unique_p50_current_tridentata <- rgeos::gDifference(tridentata_current_quantiles[[1]],wyomingensis_current_quantiles[[1]])
    unique_p50_current_tridentata <- rgeos::gDifference(unique_p50_current_tridentata, vaseyana_current_quantiles[[1]])

unique_p50_current_vaseyana <- rgeos::gDifference(vaseyana_current_quantiles[[1]],wyomingensis_current_quantiles[[1]])
    unique_p50_current_vaseyana <- rgeos::gDifference(unique_p50_current_vaseyana, tridentata_current_quantiles[[1]])

# plot unique regions and consensus amoung ssp
dev.new(width=9.807017,height=5.385965)
par(mfrow=c(1,2))
par(mar=par()$mar/2)
plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.8)
  plot(spTransform(intersect_p50_current,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
  plot(spTransform(unique_p50_current_tridentata,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
  plot(spTransform(unique_p50_current_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#4D94DB",border=NA,add=T)
  plot(spTransform(unique_p50_current_vaseyana,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
  plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","tridentata","wyomingensis","vaseyana"), cex=0.8, fill=c("#003D7A","#005CB8","#4D94DB","#B2D1F0"),bg = "white");

plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.8)
  plot(spTransform(intersect_p50_wyo_tri_current,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
  plot(spTransform(unique_p50_current_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
  plot(spTransform(unique_p50_current_tridentata,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
  plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), border="#000000B3",add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","wyomingensis","tridentata"), cex=0.8, fill=c("#003D7A","#005CB8","#B2D1F0"),bg = "white");

# plot overlapping envelopes
dev.new(height=5.359649,width=6.219298)

plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.8)
  plot(spTransform(wyomingensis_current_quantiles[[1]],CRS(projection("+init=epsg:2163"))), col="#00330099",border=NA,add=T)
  plot(spTransform(tridentata_current_quantiles[[1]],CRS(projection("+init=epsg:2163"))), col="#19751980",border=NA,add=T)
  plot(spTransform(vaseyana_current_quantiles[[1]],CRS(projection("+init=epsg:2163"))), col="#80B280B3",border=NA,add=T)
  plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), border="#000000B3",add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("wyomingensis","tridentata","vaseyana"), cex=0.8, fill=c("#003300","#197519","#80B280"),bg = "white");

# response plots for GLMs
png(paste(sep="/",HOME,"/Desktop/wyo_response_plots_glm.png"),height=1200,width=800)
  par(mfrow=c(5,1),cex.lab=2.8,cex.axis=2.8)
    for(v in vars){
      responsePlot(wyomingensis_glm_unif,var=v)
      out <- partialPlot(wyomingensis_rf_unif[[1]],x.var=as.character(v),pred.data=na.omit(wyomingensis_glm_unif[[1]][[1]]$data),which.class=1,plot=F)
      out$y <- exp(out$y);
      out$y <- (out$y/max(out$y))
      lines(y=out$y, x=out$x,lwd=1.8,col="blue",main="",xlab=as.character(v))
    }
      graphics.off()
png(paste(sep="/",HOME,"/Desktop/tri_response_plots_glm.png"),height=1200,width=800)
    par(mfrow=c(5,1),cex.lab=2.8,cex.axis=2.8)
      for(v in vars){ responsePlot(tridentata_glm_unif,var=v) }
        graphics.off()
png(paste(sep="/",HOME,"/Desktop/vas_response_plots_glm.png"),height=1200,width=800)
  par(mfrow=c(5,1),cex.lab=2.8,cex.axis=2.8)
    for(v in vars){ responsePlot(vaseyana_glm_unif,var=v) }
      graphics.off()

# response plots for RFs
png(paste(sep="/",HOME,"/Desktop/wyomingnesis_response_plots_rf.png"),height=1200,width=800)
  par(mfrow=c(5,1),cex.lab=2.8,cex.axis=2.8)
    for(v in vars){

    }
graphics.off()
png(paste(sep="/",HOME,"/Desktop/tridentata_response_plots_rf.png"),height=1200,width=800)
  par(mfrow=c(5,1),cex.lab=2.8,cex.axis=2.8)
    for(v in vars){
      partialPlot(tridentata_rf_unif[[1]],x.var=as.character(v),pred.data=na.omit(tridentata_glm_unif[[1]][[1]]$data),which.class=1,lwd=2,col="red",main="",xlab=as.character(v))
    }
graphics.off()
png(paste(sep="/",HOME,"/Desktop/vaseyana_response_plots_rf.png"),height=1200,width=800)
  par(mfrow=c(5,1),cex.lab=2.8,cex.axis=2.8)
    for(v in vars){
      partialPlot(vaseyana_rf_unif[[1]],x.var=as.character(v),pred.data=na.omit(vaseyana_glm_unif[[1]][[1]]$data),which.class=1,lwd=2,col="red",main="",xlab=as.character(v))
    }
graphics.off()

#source("do_2050_projections.R")
#source("do_2070_projections.R")

# cat(" -- projecting model rasters\n")
#   r_glm_wyomingensis_current <- predict(climate_variables, wyomingensis_glm[[1]][[1]], type='resp', progress='text')
#   r_rf_wyomingensis_current  <- predict(climate_variables, wyomingensis_rf[[1]][[1]], type='resp', progress='text')

# make some raster projections
# dev.new(height=5,width=8.5)
# par(mfrow=c(1,3))
# par(mar=par()$mar/1.5)
# mask <- r_glm_tridentata_current >= tridentata_glm[[2]][,4] # mask with binary image
#   plot(r_glm_tridentata_current*mask,main="ssp tridentata (current climate -- glm)"); plot(boundaries, add=T, border="DarkGrey")
# mask <- r_glm_wyomingensis_current >= wyomingensis_glm[[2]][,4]
#   plot(r_glm_wyomingensis_current*mask,main="ssp wyomingensis (current climate -- glm)"); plot(boundaries, add=T, border="DarkGrey")
# mask <- r_glm_vaseyana_current >= vaseyana_glm[[2]][,4]
#   plot(r_glm_vaseyana_current*mask,main="ssp vaseyana (current climate -- glm)"); plot(boundaries, add=T, border="DarkGrey")

# dev.new(height=5,width=8.5)
# par(mfrow=c(1,3))
# par(mar=par()$mar/1.5)
#   plot(r_rf_tridentata_current,main="ssp tridentata (current climate -- rf)"); plot(boundaries, add=T, border="DarkGrey")
#   plot(r_rf_wyomingensis_current,main="ssp wyomingensis (current climate -- rf)"); plot(boundaries, add=T, border="DarkGrey")
#   plot(r_rf_vaseyana_current,main="ssp vaseyana (current climate -- rf)"); plot(boundaries, add=T, border="DarkGrey")

cat(" -- done\n")
