require(sqldf)
require(tcltk)
require(rgdal)
require(raster)

# priority species for refuge planning, as mined from Rio Mora's LPP
priority_spp_four_letter_codes <- c(
  "LOSH",  # loggerhead shrike
  "LASP",  # lark sparrow
  "PRFA",  # prairie falcon
  "PIJA",  # pinyon jay
  "LBCU",  # long-billed curlew
  "BUOW",  # burrowing owl
  "CASP",  # cassin's sparrow
  "CCLO",  # chestnut-collared longspur
  "BEVI",  # bell's vireo
  "GOEA",  # golden eagle
  "FEHA",  # ferruginous hawk
  "GRSP",  # grasshopper sparrow
  "LEWO",  # lewis's woodpecker
  "MOPL",  # mountain plover
  "NOHA",  # northern harrier
  "SPOW",  # spotted owl (Mexican)
  "PEFA",  # peregrin falcon
  "SWFL",  # southwestern willow flycatcher
  "SWHA",  # swainson’s hawk
  "YEWA"   # yellow warbler
)

priority_spp_binomials <- c(
  "laniusludovicianus",       # loggerhead shrike
  "chondestesgrammacus",      # lark sparrow
  "falcomexicanus",           # prairie falcon
  "gymnorhinuscyanocephalus", # pinyon jay
  "numeniusamericanus",       # long-billed curlew
  "athenecunicularia",        # burrowing owl
  "peucaeacassinii",          # cassin's sparrow
  "calcariusornatus",         # chestnut-collared longspur
  "vireobellii",              # bell's vireo
  "aquilachrysaetos",         # golden eagle
  "buteoregalis",             # ferruginous hawk
  "ammodramussavannarum",     # grasshopper sparrow
  "melanerpeslewis",          # lewis's woodpecker
  "charadriusmontanus",       # mountain plover
  "circuscyaneus",            # northern harrier
  "strixoccidentalislucida",  # spotted owl (Mexican)
  "falcoperegrinus",          # peregrin falcon
  "empidonaxtrailliiextimus", # southwestern willow flycatcher
  "buteoswainsoni",           # swainson’s hawk
  "setophagapetechia"         # yellow warbler
)
# unpack our eBird Reference Data
if(!dir.exists("~/erd")){
  dir.create("~/erd")
  system("tar xzf ~/Incoming/ebird_number_crunching/count_data/erd_us48_data_grouped_by_year_v5.0.tar.gz -C ~/erd");
}

# if we have an existing dataset, let's use it
if(file.exists(paste(Sys.getenv("HOME"),"/Incoming/ebird_number_crunching/pts_2002_2012.csv"))){
  s <- rgdal::readOGR(paste(Sys.getenv("HOME"),"/Incoming/ebird_number_crunching/",sep=""),"pts_2002_2012")
    s@data <- read.csv(paste(Sys.getenv("HOME"),"/Incoming/ebird_number_crunching/pts_2002_2012.csv",sep=""))
} else {
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
    counts <- t[,18:ncol(t)]
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
  # parse presence absence for species
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
  grid <- raster(ext=extent(s),nrows=100,ncols=100,crs=CRS(projection("+init=epsg:4326")))
}
# create a stack for each of our priority species that we have data for
if(sum(!priority_spp_binomials %in% names(s@data))>0){
  warning("not all species binomial names were found in data frame --
    only re-gridding those species we have data for")
  in <- priority_spp_binomials %in% names(s@data)
  priority_spp_binomials <- priority_spp_binomials[in]
  priority_spp_four_letter_codes <- priority_spp_four_letter_codes[in]
}
# re-grid our original ERD points into a raster stack, then parse the stack back to points,
# treating the data as presence(1)/absence(0)
grid <- rasterize(s, field=priority_spp_binomials, fun=function(x,na.rm=T){ as.numeric(sum(x,na.rm=T)>0) }, na.rm=T)
grid_pts <- rasterToPoints(grid_attributed[[1:nlayers(grid_attributed)]],sp=T)
names(grid_pts) <- priority_spp_four_letter_codes
writeOGR(grid_pts,".","spp_grid_points_processed", driver="ESRI Shapefile",overwrite=T)
