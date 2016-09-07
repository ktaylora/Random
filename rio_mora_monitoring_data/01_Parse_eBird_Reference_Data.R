#
# Parse eBird Reference Data
#
# Attempt to parse an eBird ERD file, accepting all records found within a
# spatial extent specified in-file.  Will also perform sensoring of records
# within urban areas, if a population density raster is found after a recursive
# search of the user's home directory.
#
# Assumes:
# (1) ERD File: /Incoming/ebird_number_crunching/count_data/erd_us48_data_grouped_by_year_v5.0.tar.gz exists.
# (2) That somewhere in the user's home directory, a population density raster named gpw-v4-population-density_2015.tif exists
require(tcltk)
#

require(sqldf)
require(rgdal)
require(raster)
require(parallel)
require(landscapeAnalysis)

# priority species for refuge planning, as mined from Rio Mora's LPP
priority_spp_four_letter_codes <- c(
  # "LOSH",  # loggerhead shrike
  # "LASP",  # lark sparrow
  "PRFA",  # prairie falcon
  # "PIJA",  # pinyon jay
  # "LBCU",  # long-billed curlew
  # "BUOW",  # burrowing owl
  # "CASP",  # cassin's sparrow
  # "CCLO",  # chestnut-collared longspur
  # "BEVI",  # bell's vireo
  "GOEA",  # golden eagle
  # "FEHA",  # ferruginous hawk
  # "GRSP",  # grasshopper sparrow
  # "LEWO",  # lewis's woodpecker
  # "MOPL",  # mountain plover
  # "NOHA",  # northern harrier
  # "SPOW",  # spotted owl (Mexican)
  "PEFA"  # peregrin falcon
  # "SWFL",  # southwestern willow flycatcher
  # "SWHA",  # swainson’s hawk
  # "YEWA"   # yellow warbler
)

priority_spp_binomials <- c(
  # "laniusludovicianus",       # loggerhead shrike
  # "chondestesgrammacus",      # lark sparrow
  "falcomexicanus",           # prairie falcon
  # "gymnorhinuscyanocephalus", # pinyon jay
  # "numeniusamericanus",       # long-billed curlew
  # "athenecunicularia",        # burrowing owl
  # "peucaeacassinii",          # cassin's sparrow
  # "calcariusornatus",         # chestnut-collared longspur
  # "vireobellii",              # bell's vireo
  "aquilachrysaetos",         # golden eagle
  # "buteoregalis",             # ferruginous hawk
  # "ammodramussavannarum",     # grasshopper sparrow
  # "melanerpeslewis",          # lewis's woodpecker
  # "charadriusmontanus",       # mountain plover
  # "circuscyaneus",            # northern harrier
  # "strixoccidentalislucida",  # spotted owl (Mexican)
  "falcoperegrinus"          # peregrin falcon
  # "empidonaxtrailliiextimus", # southwestern willow flycatcher
  # "buteoswainsoni",           # swainson’s hawk
  # "setophagapetechia"         # yellow warbler
)

#
# MAIN
#

# if we have an existing dataset, let's use it
if(file.exists(paste(Sys.getenv("HOME"),"/Incoming/ebird_number_crunching/pts_2002_2012.csv",sep=""))){
  s <- rgdal::readOGR(paste(Sys.getenv("HOME"),"/Incoming/ebird_number_crunching/",sep=""),"pts_2002_2012",verbose=F)
    s@data <- read.csv(paste(Sys.getenv("HOME"),"/Incoming/ebird_number_crunching/pts_2002_2012.csv",sep=""))
  raster::projection(s) <- CRS(raster::projection("+init=epsg:4326"))
} else {
  # unpack our eBird Reference Data
  if(!dir.exists("~/erd")){
    dir.create("~/erd")
    system("tar xzf ~/Incoming/ebird_number_crunching/count_data/erd_us48_data_grouped_by_year_v5.0.tar.gz -C ~/erd");
  }
  # Treat each year as an independent observation within the time-series
  merged <- list()
  cat(" -- processing:");
  for(i in 2002:2012){
    cat(paste(i,",",sep=""))
    f <- file(paste("~/erd",i,"checklists.csv",sep="/"))
    # Bounding box for the Rio Mora (tri-county area): -105.735,-103.36; 35.04,37.011
    t <- sqldf("SELECT * from f WHERE LONGITUDE <= -100.6 and LONGITUDE >= -108.735
                and LATITUDE <= 40.011 and LATITUDE >= 32.04",
                # and COUNT_TYPE in ('P22','P34')", # P22 and P34 are eBird's "Traveling" designators
               dbname = tempfile(),
               file.format = list(header = T, row.names = F)
               )
    # parse species names into something greppable
    counts <- t[,18:ncol(t)] # first 17 columns are junk
    names <- names(counts)
      names <- tolower(gsub(names(counts),pattern="[.]|_",replacement=""))
        names(counts) <- names
          counts <- counts[,nchar(names)<200] # bug-fix: throw-out ridiculously long species names as junk
    # get spatial and temporal information
    t <- t[,grepl(names(t),pattern="LATITUDE|LONGITUDE|YEAR|DAY")]
      t <- t[,c(1,2,3,4)]
    merged[[length(merged)+1]] <- cbind(t,counts)
  }; cat("\n");
  # merge and create spatial points with a year field to check for clustering in space and time in outside analysis
  merged <- do.call(rbind,merged)
  s <- SpatialPointsDataFrame(data.frame(x=merged$LONGITUDE, y=merged$LATITUDE),
         data=data.frame(year=merged$YEAR,day=merged$DAY))
    writeOGR(s,paste(Sys.getenv("HOME"),"/Incoming/ebird_number_crunching",sep=""),
      "pts_2002_2012",driver="ESRI Shapefile", overwrite=T)
    write.csv(merged,paste(Sys.getenv("HOME"),
      "/Incoming/ebird_number_crunching/pts_2002_2012.csv",sep=""))
  # parse presence absence for priority species
  n <- names(merged)
  for(i in 1:length(priority_spp_binomials)){
    if(is.null(ncol(merged[,grepl(n,pattern=priority_spp_binomials[i])]))){
      v <- merged[,grepl(n,pattern=priority_spp_binomials[i])] >= 1
    } else {
      v <- as.numeric(as.numeric(apply(merged[,grepl(n,pattern=priority_spp_binomials[i])],MARGIN=1,FUN=max)) >= 1)
    }
    v[is.na(v)] <- 1 # "X" is noted on some submissions -- counted as a presence
      cat(paste(priority_spp_four_letter_codes[i],": ",sum(v),"/",nrow(merged),"\n",sep=""))
  }
}
# ERD marks some presences with "checked" or other non-numeric values -- swap them out with 1's
for(i in 1:ncol(s@data)){
  x <- as.numeric(s@data[,i]);
    x[is.na(x)] <- 1;
  s@data[,i] <- x
}
# re-grid our data onto a 100 x 100 unit surface (or use a grid specified by the user in the CWD)
if(file.exists("grid.tif")){
  grid <- raster("grid.tif")
} else {
  grid <- raster(ext=extent(s),nrows=200,ncols=200,crs=CRS(projection("+init=epsg:4326"))) # arbitrary grid
  # grid <- raster::raster(spTransform(s,CRS(projection("+init=epsg:2163"))),resolution=632,crs="+init=epsg:2163") # grid-cell resolution = mean parcel-size in Mora County
  #   grid <- raster::projectRaster(grid,CRS(raster::projection(s)))
  #     grid <- suppressWarnings(raster::projectRaster(grid,crs=CRS(projection(s))))
}
# create a stack for each of our priority species that we have data for
if(sum(!priority_spp_binomials %in% names(s@data))>0){
  warning("not all species binomial names were found in data frame --
    only re-gridding those species we have data for")
  have <- priority_spp_binomials %in% names(s@data)
  priority_spp_binomials <- priority_spp_binomials[have]
  priority_spp_four_letter_codes <- priority_spp_four_letter_codes[have]
}

# sanity-check : make sure there are absences in the ERD to work with -- warn if not
for(i in 1:length(priority_spp_binomials)){
  if(min(s[,priority_spp_binomials[i]]@data) != 0){
    cat(" -- warning: no zero's found for:",priority_spp_four_letter_codes[i], "; adjusting values accordingly. Please check the source ERD for inconsistencies.\n")
    s@data[,priority_spp_binomials[i]] <- s@data[,priority_spp_binomials[i]]-min(s@data[,priority_spp_binomials[i]])
  }
}

# re-grid our original ERD points into a raster stack, then parse the stack back to points,
# treating the data as presence(1)/absence(0)

grid_pts <- rasterize(s,grid,field=priority_spp_binomials, fun=function(x,na.rm=T){ as.numeric(sum(x,na.rm=T)>0) }, na.rm=T)
  grid_pts <- rasterToPoints(grid_pts[[1:nlayers(grid_pts)]],sp=T)
    names(grid_pts) <- priority_spp_four_letter_codes

# calculate species richness and append to our table
grid_pts$RICHNESS <- rowSums(grid_pts@data)

# sensor eBird observations that occur in urban areas (i.e., population densities > 150 sq-km)
cat(" -- searching for population density raster: gpw-v4-population-density_2015.tif\n")
population <- raster(list.files(Sys.getenv("HOME"),pattern="gpw-v4-population-density_2015.tif$",recursive=T,full.names=T)[1])
  population <- landscapeAnalysis::snapTo(population,to=grid)

keep <- extract(population,grid_pts)
  keep[is.na(keep)] <- 151

grid_pts <- grid_pts[keep <= 150,]

# Downsample over-abundant richness "bins" to a consistent density to deal with value inflation
# In effect, we are sensoring average values of richness and allowing more extreme values to have a greater
# influence on model fit.
#
# sanity-check : make sure we have some heterogeneity in species richness in the first place
#

if(length(unique(abs(grid_pts$RICHNESS-mean(grid_pts$RICHNESS)))) > 1){
  # convert RICHNESS to units of standard deviation
  unit_sd <- abs(grid_pts$RICHNESS-mean(grid_pts$RICHNESS))/sd(grid_pts$RICHNESS)
     keep <- vector() # source data values we will keep
  # build a histogram that we can use to bin our standardized data
  h <- hist(abs(unit_sd),plot=F,breaks=30)
  # downsample over-abundant bins to a density consistent with the mean count across bins
  target_count <- round(mean(h$counts[h$counts>0]))
          bins <- h$mids[h$counts > target_count]
     bin_width <- h$breaks[2]-h$mids[1]
  for(i in 1:length(bins)){
    keep <- append(keep,sample(which(unit_sd < bins[i]+bin_width & unit_sd > bins[i]-bin_width),size=target_count))
  }
  # now restore all values from the tail of the distribution and parse our grid_pts
  keep <- append(keep,which(unit_sd > bins[length(bins)]+bin_width))
  grid_pts <- grid_pts[keep,]
} else {
  cat(" -- warning: all locations had the same species richness value.  This doesn't bode well for modeling.\n")
}

# write to disk
cat(" -- writing grid_pts to disk: spp_grid_points_processed.shp\n")
writeOGR(grid_pts,".","spp_grid_points_processed", driver="ESRI Shapefile",overwrite=T)
