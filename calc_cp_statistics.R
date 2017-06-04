#
# CP1/2 Grassland Productivity Summary Analysis
# This workflow details a pre-cursor analysis to looking at the feasibility of
# using CP1/2 treatments in a region-wide IMBCR analysis
#
# Author : Kyle Taylor (kyle.taylor@pljv.org) [2017]
#

require(raster)
require(rgdal)
require(rgdal)
require(ggplot2)

RASTER_DIR <- "/home/ktaylora/Workspace/tpw_imbcr_grassland_birds_workflow/Raster"

#' code for re-scaling proportional votes to the range of predicted values observed
#' in vector d
#' @export
scale.dist <- function(d) (exp(d) - min(exp(d)))/(max(exp(d)) - min(exp(d)))
logRescalePorpVotes <- function(pr=NULL, class=1, w=rep(1,nrow(pr)), focus=NULL) {
  focus <- which(grepl(colnames(pr), pattern=class))
    pr[pr == 0] <- .Machine$double.eps
      return( stats::weighted.mean(log(pr[,focus]) - rowMeans(log(pr), w, na.rm = TRUE)) )
}
#' generate a partial plot for a Random Forests model object, given some input table, a focal predictor variable, and a target class
#' @export
partialPredict <- function(m=NULL,t=NULL,var=NULL,class=1){
  # re-train a new model for each unique value of x, averaging the predictions across all non-focal variables as we go
  x <- seq((min(t[,var])-sd(t[,var])),(max(t[,var])+sd(t[,var])),length.out=100)
  y <- rep(0,length(x))
  for(i in 1:length(x)){
    run_table <- t
      run_table[,var] <- x[i]
    p <- predict(m,run_table,type="prob")
      y[i] <- logRescalePorpVotes(p)
  }
  y <- scale.dist(y)
  return(data.frame(x=x,y=y))
}
#' given a smoother data.frame (returned by rfLoessSmoother), make a
#' pretty plot with ggplot2
#' @export
loessPartialPlot <- function(smoother_df=NULL,
                             median_pred=NULL,
                             focal_var=NULL, sourceTable=NULL, plot=T,
                             xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL){
  require(ggplot2)

  xlim <- if(is.null(xlim)) range(smoother_df$x, na.rm=T) else xlim
  ylim <- if(is.null(ylim)) range(smoother_df$y, na.rm=T) else ylim
  xlab <- if(is.null(xlab)) focal_var else xlab
  ylab <- if(is.null(ylab)) "Probability Lek is Active" else ylab

  # build a rug compatible with ggplot by randomly sampling
  # n=number of rows in a typical source table from the full range
  # of x-values considered across the span of random forest models
  rug <- as.vector(
           unlist(do.call(cbind,
                          lapply(sourceTable,FUN=function(x) x[,focal_var])
                )))
  rug <- data.frame(fit=0,x=sample(rug,size=nrow(sourceTable[[1]])))
  # build a data.frame for our plot data out of out_l
  t <- as.data.frame(smoother_df)
    t$lower <- t$fit-(t$fit*qnorm(0.8))
    t$upper <- t$fit+(t$fit*qnorm(0.8))
  # plot specification
  ggplot(data=t, aes(x,fit)) +
    xlab(xlab) +
    ylab(ylab) +
    geom_line(size=0.5, alpha=0.7) +
    geom_line(aes(x, y), data=median_pred, linetype=2, alpha=0.3) +
    geom_rug(aes(x), sides="b", alpha=0.1, size=0.5, data=rug) +
    geom_ribbon(data=t, aes(ymin=lower,ymax=upper),alpha=0.15) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    theme(axis.text = element_text(face="bold",size=15),
          text = element_text(face="bold",size=15)) +
    xlim(xlim)
}
#' hack that adds a loess smoother to a random forest partial plot.  Need to add in a step that iteritively rebuilds
#' the input dataset and trains a new forest for each iteration, using model predictions to construct a predictive interval,
#' ala Baruch-Mordoche et al.  Consider this an MVP.
#' @export
rfLoessSmoother <- function(m=NULL, sourceTable=NULL, focal_var=NULL, class=1,
                             main=NULL, xlab=NULL, xlim=NULL, ylim=NULL, plot=T,
                             degree=2, nCores=NULL){
  require(parallel)
  require(randomForest)
  nCores <- if (is.null(nCores)) parallel::detectCores()-1 else nCores
    cl <- parallel::makeCluster(nCores)
      clusterExport(cl=cl,varlist=c("randomForest","logRescalePorpVotes","scale.dist"))
  #' sum of slopes calculates the sum of slopes for each predicted value of pred
  sum_of_slopes <- function(pred) sum(lm(as.formula("y~x"),data=pred)$coefficients)
  # have we built a series of RF models that we are leveraging for our loess parial plot?
  if(class(m) == "list" & class(sourceTable) == "list"){
    cat(" -- predicting across variable series for all models:")
    predictions <- clusterMap(cl, m, sourceTable, fun=partialPredict, var=focal_var)
    parallel::stopCluster(cl)
    # find the median prediction across our model series using a simple coefficient of 'x' determination
    median <- unlist(lapply(predictions, FUN=sum_of_slopes))
    median <- which(abs(median-median(median,na.rm=T))==min(abs(median-(median(median,na.rm=T)))))
      median <- median[length(median)] # in case there are two-or-more identical preds
        median <- predictions[[median]]
    # bind all of our predictions together to estimate a predictive interval
    predictions <- do.call(rbind,predictions)
    # bug fix : don't try to use loess on an extremely large table
    if(nrow(predictions)>30000){
      size <- 30000
    } else {
      size <- nrow(predictions)
    }
    # run a loess smoother across our predictions to make some envelopes
    x <- seq(min(predictions$x),max(predictions$x),length.out=500)
    m_loess <- smoother_df <- loess(y~x,data=predictions[sample(1:nrow(predictions),
                                                         size=size),],degree=degree)
      smoother_df <- predict(smoother_df,newdata=x,se=T)
        smoother_df$x <- x
    smoother_df$sd.fit <- 2*(sqrt(length(m))*smoother_df$se.fit) # 2(SD) = 2(se*sqrt(n))
      smoother_df <- as.data.frame(smoother_df)

  # are we fitting a loess curve to a single model/table/variable combination?
  } else {
    if(sum(names(sourceTable) %in% focal_var)==0){
      stop(paste("focal_var=",focal_var," not found in sourceTable\n",sep=""))
    }
    predictions <- sourceTable
    median <- sourceTable
    out <- partialPredict(m,t=predictions,var=focal_var[1],class=class,main=main,plot=F)
    smoother_df <- loess("y~x",data=out,degree=degree)
      smoother_df <- predict(smoother_df,se=T)
        smoother_df$x <- x
    smoother_df$sd.fit <- 2*(sqrt(length(m))*se.fit)
      smoother_df <- as.data.frame(smoother_df)
  }
  # if asked, build our plots with ggplot2
  if(plot){
    loessPartialPlot(smoother_df=smoother_df, median_pred=median,
                     focal_var=focal_var,
                     sourceTable=sourceTable, xlim=xlim, ylim=ylim,
                     xlab=xlab, ylab=ylab)
  }
  # clean-up
  if(!plot){
    ret <- vector('list', 3)
    ret[[1]] <- m_loess
    ret[[2]] <- median
    ret[[3]] <- smoother_df
    names(ret) <- c("loess_model","median_prediction","smoother_df")
    return(ret)
  }
}
#' rebag our input data table -- and balance the presence/absences by randomly-downsampling the over-abundant class
#' @export
rebag <- function(t, y.var="response"){
  final <- t[sample(which(t[, y.var] == 1),
    size=sum(t[, y.var] == 0)), ]
  final <- rbind(final, t[t[, y.var] == 0,])
  # drop our latent variables used for testing
  final <- final[, ! grepl(names(final), pattern="ID$|X$")]

  final <- final[sample(1:nrow(final),size=nrow(final),replace=F),] # randomise the order of our input data
}
#' bootstrap the rf model selection process by re-sampling (with replacement) the input training
#' dataset. Returns a table with summary importance metrics averaged across n boot-strapped model runs.
#'
bs_model_select <- function(training_data=NULL,
                          n=100, nCores=NULL){
  require(parallel)
  nCores <- if (is.null(nCores)) parallel::detectCores()-1 else nCores
  cl <- parallel::makeCluster(nCores)
    clusterExport(cl=cl,varlist=c("rebag"))
  # pre-allocate what will become an aggregate importance table
  importance_table <- vector('list', n)
  # define local function for model_selection with apply comprehension
  model_select <- function(t=NULL, x=NULL){
    require(rfUtilities)
    t <- rebag(t)
    m <- rfUtilities::rf.modelSel( xdata=t[,!grepl(names(t),
                                                   pattern="response")],
                                   ydata=as.factor(t$response),
                                   imp.scale="mir",
                                   r=seq(0,1,0.015),
                                   ntree=800,
                                   seed=25
                                  )
    i <- m$importance
      i <- data.frame(var=rownames(i),imp=i[, 1])
        return(i)
  }
  cat(" -- parallelizing 'rfUtilities' model selection across our cluster\n")
  importance_table <- parallel::parLapply(cl,
                                          importance_table,
                                          fun=model_select,
                                          t=training_data)
  # pretty-up
  var_names <- as.character(importance_table[[1]]$var)
  importance_table <- do.call(cbind,importance_table)
    importance_table <- as.data.frame(
                          cbind(var_names,
                            round(rowMeans(importance_table[,!grepl(colnames(importance_table),pattern="var")]),3)
                        ))
      rownames(importance_table) <- NULL
        colnames(importance_table) <- c("vars","imp")
  # re-order our table and return to user
  return(importance_table[order(importance_table[,2],decreasing=T),])
}
#' bootstrap random forest
#'
bs_random_forest <- function(training_data=NULL,
                             vars=NULL,
                             n=500, nCores=NULL){
  require(parallel)
  nCores <- if (is.null(nCores)) parallel::detectCores()-1 else nCores
      cl <- parallel::makeCluster(nCores)
  clusterExport(cl=cl,varlist=c("rebag"))
  # train our final model(s) for fancy-dancey partial plots with confidence intervals, similar to Baruch-Mordo et al. (2013)
  random_forest <- function(x=NULL,t=NULL,vars=vars){
    require(randomForest)
    t <- rebag(t)
    m <- randomForest(y=as.factor(t$response),
                      x=t[,vars],
                      do.trace=T,
                      ntree=800,
                      importance=T)
    ret <- list(m,t)
      names(ret) <- c("models","tables")
        return(ret)
  }
  #' get importance from a random forest object with apply comprehension
  get_importance <- function(x){
    i <- x$models$importance
      return(data.frame(var=rownames(i),imp=i[,'1'])) # we want importance for our "active" class
  }
  cat(" -- parallelizing 'randomForest' across cluster\n")
  runs <- parallel::parLapply(cl,vector('list',n), fun=random_forest, t=training_data, vars=vars)
    names(runs) <- "modeling"
  parallel::stopCluster(cl)
  cat(" -- building rf importance tables\n")
  importance_table <- do.call(cbind,lapply(runs, FUN=get_importance))
    importance_table <- rowMeans(importance_table[,!grepl(colnames(importance_table), pattern="var")])
  # pretty-up the output and normalize importance (rel to most important variable)
  importance_table <- data.frame(var=names(importance_table),imp=round(importance_table,3))
    importance_table[,2] <- round(importance_table[,2]/max(importance_table[,2]),3)
      rownames(importance_table) <- NULL

  runs$overall_importance <- importance_table[order(importance_table[,2],
                                                    decreasing=T),]
  return(runs)
}
#' uses a cluster to run random predict across a series variables
#' using values returned by our loess smoother (above) to estimate
#' the range of y values taken across all models
#' @export
calc_effect_size <- function(models=NULL,
                             vars=NULL, nCores=NULL){
  require(parallel)
  # pre-allocated result frame containing variable name and effect size
  res <- data.frame(matrix(0,nrow=length(vars),ncol=2))
    res[,1] <- vars
      colnames(res) <- c("variable","effect")

  nCores <- if (is.null(nCores)) parallel::detectCores()-1 else nCores

  cl <- parallel::makeCluster(nCores)
  clusterExport(cl=cl,
         varlist=c("scale.dist","logRescalePorpVotes","partialPredict","rebag"))
  # predict (without plotting) across our variable list using
  # our cluster
  effect_size <-
  parallel::parLapply(cl, X=as.list(vars), fun=loessPartialPlot,
            m=l_unpack(models,'models'),
            sourceTable=l_unpack(models,'tables'),
            ylim=c(0,1),plot=F)
  # clean-up
  parallel::stopCluster(cl)
  res[,2] <- unlist(lapply(effect_size,
                           FUN=function(e) {
                             e$fit[e$fit<0] <- 0
                             e$fit[e$fit>1] <- 1
                             diff(range(e$fit,na.rm=T))
                           }
                    ))
  # sort and return to user
  return(res[order(res[,2],decreasing=T),])
}
#' shorthand fetch list from a 'list-of-lists' by list name
l_unpack <- function(x=NULL,var=NULL) {
  x <- lapply(x, FUN=function(x) x[[var]])
  valid <- as.vector(!unlist(lapply(x,is.null)))
  if(sum(!valid)>0){
    warning("null values detected in l_unpack. dropping from output.")
  }
  return(x[valid])
}
#' downsample our CP1/2 sample-space across seasons
post_stratify_across_seasons <- function(s=NULL, write=NULL,
                                         method="quantile", ...){
  seasons <- c('spring','summer','fall','winter')
  for (season in seasons){
    cat(" -- season:",season,"\n")
    # testing quantiles vs. normal dist for downsampling
    if(grepl(tolower(method),pattern="norm")){
      s = s[OpenIMBCR::downsample_by_normal_dist(s@data[,season],
        use_mean=T, byid=T, ...), ]
    } else if(grepl(tolower(method),pattern="quant")){
      s = s[OpenIMBCR::downsample_by_quantile(s@data[,season], byid=T,
            ...), ]
    }
  }
  if (!is.null(write)){
    rgdal::writeOGR(s,
        ".", layer=write,
        driver="ESRI Shapefile",
        overwrite=T
      )
  }
  return(s)
}

simple_hist_plot <- function(x,xlab=NULL,main=NULL, xlim=NULL){
  hist(x, main=main, xlab=xlab, breaks=50,
       col="#4A8EBC", border="#D6EEFF",
       ylab=NULL, xlim=xlim
    )
  grid(lty=2,lwd=0.5,col="LightGrey")
  abline(lwd=1.5, v=mean(x)+(sd(x)/sqrt(length(x))),col="DarkGrey")
  abline(lwd=1.5, v=mean(x)-(sd(x)/sqrt(length(x))),col="DarkGrey")
  abline(lwd=1.5, v=mean(x),col="red")
  abline(lwd=1.5, v=mean(x)+sd(x),col="Grey")
  abline(lwd=1.5, v=mean(x)-sd(x),col="Grey")
  abline(lwd=1.5, v=mean(x)+(2*sd(x)),col="LightGrey")
  abline(lwd=1.5, v=mean(x)-(2*sd(x)),col="LightGrey")
}

seasonal_cp1_cp2_ndvi_histograms <- function(s_cp1=NULL,
                                             s_cp2=NULL,
                                             xlim=c(-0.15,0.30)){
  dev.new()
  par(mfrow=c(2,4))
  simple_hist_plot(s_cp1$spring,
      xlab="CP1 (Spring)", xlim=xlim
    )
  simple_hist_plot(s_cp1$summer,
      xlab="CP1 (Summer)", xlim=xlim
    )
  simple_hist_plot(s_cp1$fall,
      xlab="CP1 (Fall)", xlim=xlim
    )
  simple_hist_plot(s_cp1$winter,
      xlab="CP1 (Winter)", xlim=xlim
    )
  simple_hist_plot(s_cp2$spring,
      xlab="CP2 (Spring)", xlim=xlim
    )
  simple_hist_plot(s_cp2$summer,
      xlab="CP2 (Summer)", xlim=xlim
    )
  simple_hist_plot(s_cp2$fall,
      xlab="CP2 (Fall)", xlim=xlim
    )
  simple_hist_plot(s_cp2$winter,
      xlab="CP2 (Winter)", xlim=xlim
    )
}

calc_seasonal_ndvi_by <- function(s=NULL,r=NULL){
  s = sp::spTransform(s,sp::CRS(raster::projection(r)))
  return(raster::extract(r,s,na.rm=T,fun=mean,progress='text',df=T))
}

sort_ndvi_rasters_by_season <- function(path=NULL, pattern=NULL){
  f <- list.files(path, pattern=pattern,full.names=T)
  spring <- f[grep(x=f,pattern="spring")[1]]
  summer <- f[grep(x=f,pattern="summer")[1]]
  fall   <- f[grep(x=f,pattern="fall")[1]]
  winter <- f[grep(x=f,pattern="winter")[1]]
  f <- raster::stack(c(spring,summer,fall,winter))
    names(f) <- c("spring","summer","fall","winter")
  return(f)
}

calc_transect_centroids <- function(s=NULL){
  centroids <- s[!duplicated(s$transectnum),]
  centroids <- rgeos::gCentroid(centroids,byid=T,id=centroids$transectnum)
  centroids <- sp::SpatialPointsDataFrame(centroids,
      data=data.frame(transectnum=s[!duplicated(s$transectnum),'transectnum'])
    )
  centroids@data <- data.frame(transectnum=centroids@data[,1])
  return(centroids)
}

calc_distance_to_features <- function(s=NULL, to=NULL, converter=6.21371e-4){
  transect_dist_to <- rgeos::gDistance(calc_transect_centroids(s),
    to, byid=T)
  transect_dist_to <- data.frame(distance=apply(transect_dist_to, MARGIN=2, FUN=function(x){
      round(min(x*converter),2) # meters -> miles
    }))
  return(transect_dist_to)
}
#' we have to talk about your function names...
calc_ndvi_at_treatment_sample_pts <- function(pattern="cp1.*.ndvi.*.tif$", write=NULL){
  cp_seasons <- sort_ndvi_rasters_by_season(RASTER_DIR, pattern=pattern)
    cp_seasonal_points <- raster::rasterToPoints(cp_seasons,spatial=T,na.rm=T)
  if(!is.null(write)){
    writeOGR(cp_seasonal_points, ".",
        layer=write, driver="ESRI Shapefile", overwrite=T
      )
  }
  return(cp_seasonal_points)
}

#
# Grab IMBCR data
#

s <- spTransform(OpenIMBCR:::imbcrTableToShapefile(OpenIMBCR:::
    recursiveFindFile("master_imbcr_table.csv")),
    CRS(projection("+init=epsg:2163"))
  )

#
# Read-in treatment vector datasets
#

s_cp_1 <- sp::spTransform(readOGR("Vector","cp_1", verbose=F),
    CRS(projection("+init=epsg:2163"))
  )
s_cp_2 <- sp::spTransform(readOGR("Vector","cp_2", verbose=F),
    CRS(projection("+init=epsg:2163"))
  )

#
# Calculate distance to treatment vectors
#

transect_dist_to_cp_1 <- calc_distance_to_features(s=s, to=s_cp_1)
transect_dist_to_cp_2 <- calc_distance_to_features(s=s, to=s_cp_2)

dev.new()
ggplot(transect_dist_to_cp_1, aes(distance)) +
  geom_histogram(aes(y=..density..), color="grey",fill="#879EAD", bins=50) +
  geom_density(color="#1F78B4") +
  xlab("Distance from Transect to Nearest CP2 Field (miles)") +
  theme(plot.title=element_text(size=16),
        axis.title.x=element_text(size = 16),
        axis.title.y=element_text(size = 16),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size = 14))

dev.new()
ggplot(transect_dist_to_cp_2, aes(distance)) +
  geom_histogram(aes(y=..density..), color="grey",fill="#879EAD", bins=50) +
  geom_density(color="#1F78B4") +
  ylab("Sample Density") +
  xlab("Distance from Transect to Nearest CP2 Field (miles)") +
  theme(plot.title=element_text(size=16),
        axis.title.x=element_text(size = 16),
        axis.title.y=element_text(size = 16),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size = 14))

#
# Calculate NDVI summary statistics
#

# cp1_seasons <- sort_ndvi_rasters_by_season(RASTER_DIR, pattern="cp1.*.ndvi.*.tif$")
#   cp1_seasonal_points <- sampleRandom(cp1_seasons,size=99999,sp=T,na.rm=T)
#     writeOGR(cp1_seasonal_points, ".", "cp1_seasonal_points", driver="ESRI Shapefile")

system.time(
  cp1_ndvi_samples <- calc_ndvi_at_treatment_sample_pts(
      pattern="cp1.*.ndvi.*.tif$",
      write="cp1_seasonal_points"
    )
)

# system.time(
#   cp1_seasonal_points_post_stratified <- post_stratify_across_seasons(
#       cp1_ndvi_samples,
#       write="cp1_seasonal_points_post_stratified",
#       method="normal",
#       bins=11,
#       shape=1.2
#     )
#   )

system.time(
    cp1_seasonal_points_post_stratified <- post_stratify_across_seasons(
        cp1_ndvi_samples,
        write="cp1_seasonal_points_post_stratified",
        gte=T,
        p=0.07
      )
  )

system.time(
    cp2_ndvi_samples <- calc_ndvi_at_treatment_sample_pts(
        pattern="cp2.*.ndvi.*.tif$",
        write="cp2_seasonal_points",
        gte=T,
        p=0.07
      )
  )

system.time(
    cp2_seasonal_points_post_stratified <- post_stratify_across_seasons(
        cp2_ndvi_samples,
        write="cp2_seasonal_points_post_stratified",
        gte=T,
        p=0.07
      )
  )

#
# Take a peek at our sample means
#

# Means comparison
seasonal_cp1_cp2_ndvi_histograms(
    cp1_ndvi_samples,
    cp2_ndvi_samples
  )

seasonal_cp1_cp2_ndvi_histograms(
    cp1_seasonal_points_post_stratified,
    cp2_seasonal_points_post_stratified,
    xlim=c(0.15,0.30)
  )

# ML classifier comparison

rf_full_dataset <- rbind(cp1@data,cp2@data)
rf_full_dataset <- cbind(treatment=c(rep("cp1",nrow(cp1@data)),
    rep("cp2",nrow(cp2@data))),rf_full_dataset
  )

m <- randomForest(as.factor(treatment)~.,
    data=rf_testframe[sample(1:nrow(rf_testframe),size=20000), ],
    ntree=1000, do.trace=T
  )
