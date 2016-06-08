#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2016]
#

require(raster)
require(rgdal)
require(rgeos)
require(landscapeAnalysis)
require(randomForest)

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

## CLEAN-UP MONITORING DATA AND REPORT OBSERVATIONS TO USER

t_monitoring_data <- read.csv(
                           "~/Incoming/rio_mora_ld/rio_mora_transect_data.csv")
      t_locations <- read.csv(
                      "~/Incoming/rio_mora_ld/rio_mora_transect_locations.csv")

      t_locations$Area <- toupper(gsub(as.vector(t_locations$Area),pattern=" ",replacement=""))

# strip out observations not taken during a minute interval
suppressWarnings( t_monitoring_data <- t_monitoring_data[!
  is.na(as.numeric(as.vector(t_monitoring_data$Minute.Interval))), ] )
# don't trust observations that haven't been coded with a 4-letter species code
t_monitoring_data <- t_monitoring_data[nchar(as.vector(t_monitoring_data$X4.letter.code)) == 4, ]
# build a consistent string for our site names
t_monitoring_data$Array.Site <- gsub(toupper(
                                t_monitoring_data$Array.Site),
                                pattern = " ",
                                replacement = ""
                                )
# determine season and year from date string and add to attribute table
date <- strsplit(as.vector(t_monitoring_data$Date), split = "/")
for (i in 1:length(date)){
  season <- vector()
    year <- as.numeric(date[[i]][3]) + 2000
  if (as.numeric(date[[i]][1]) < 4 || as.numeric(date[[1]][1]) == 12){
    season <- "winter"
  } else if (as.numeric(date[[i]][1]) < 7){
    season <- "spring"
  } else if (as.numeric(date[[i]][1]) < 9){
    season <- "summer"
  } else {
    season <- "fall"
  }
  t_monitoring_data[i, "Year"] <- year
  t_monitoring_data[i, "Season"] <- season
}
# report seasonality of detections to user
dev.new()
barplot(table(t_monitoring_data$Season),
        ylab = "detections",
        main = "seasonality of detections for all species")
# report summary statistics for number of priority species detections
cat(" -- number of priority species detected in monitoring data: ",
    length(
      unique(
        as.vector(
          t_monitoring_data$X4.letter.code)[
            (toupper(as.vector(t_monitoring_data$X4.letter.code))
            %in%
            priority_spp_four_letter_codes)]
      )
    ),
    "/",
    length(priority_spp_four_letter_codes),
    "\n", sep = ""
)

# determine how many of our priority species were observed in the raw
# monitoring data and report to the user
t <- sort(
      table(
        as.vector(t_monitoring_data$X4.letter.code)
        [
          (toupper(as.vector(t_monitoring_data$X4.letter.code))
            %in%
          priority_spp_four_letter_codes)
        ]
      )
     )

t_missing <- priority_spp_four_letter_codes[
             !(
               priority_spp_four_letter_codes
               %in%
               (toupper(as.vector(t_monitoring_data$X4.letter.code)))
              )
             ]

a <- as.table(rep(0, length(t_missing)))
  names(a) <- t_missing
    t_missing <- a
      t <- sort(c(t, t_missing), decreasing = T)
        rm(t_missing, a);

dev.new()
p <- barplot(t, xaxt = "n", ylab = "detections", xlab = "priority species", main = "Rio Mora Bird Monitoring Data : 2012-14")
  text(cex = 0.75, x = p + 0.65,  y = -0.45,  names(t),  xpd = TRUE,  srt = 45,  pos = 2)
    grid(lwd = 1.2)

## PARSE OUR MONITORING DATA SO WE ARE ONLY OBSERVING DATA FOR PRIORITY SPECIES
## AND TABULATE DETECTION / ABUNDANCE
# drop species that are not on our refuge monitoring list
t_monitoring_data <- t_monitoring_data[
  ((toupper(as.vector(t_monitoring_data$X4.letter.code))))
  %in%
  priority_spp_four_letter_codes, ]
# determine detections and abundance for each relevant species at each site/year/season/day
# please forgive this monstrosity of nested for-loops.  refactor me later.
detections <- list();
for (spp in unique(as.vector(t_monitoring_data$X4.letter.code))){
  focal <- t_monitoring_data[t_monitoring_data$X4.letter.code == spp,]
  for (y in unique(focal$Year)){
    year <- focal[focal$Year == y, ]
    if (nrow(year) > 0){
      for (S in unique(year$Season)){
        season <- year[year$Season == S, ]
        if (nrow(season) > 0){
          for (l in unique(season$Array.Site)){
            site <- season[season$Array.Site == l, ]
            if (nrow(site) > 0){
              for(m in unique(as.vector(site$Date))){
                day <- site[site$Date == m, ]
                if(nrow(day)>0){
                  detection <- rep(0,10) # 10-minute intervals
                  detection[as.numeric(as.vector(day$Minute.Interval))] <- 1;
                  # convert distance string to something tractable
                  distances <- as.vector(day$Distance..m.)
                  distances <- gsub(distances, pattern = "-", replacement = "+")
                  distances <- gsub(distances, pattern = ">", replacement = "")
                  for (i in 1:length(distances)){
                    if (grepl(distances[i], pattern = "[+]")){
                      distances[i] <- (eval(parse(text = distances[i])) / 2)
                    }
                  }
                  # build a table for focal site
                  detections[[length(detections)+1]] <-
                  data.frame(spp = as.vector(day$X4.letter.code),
                             site = day$Array.Site[1],
                             date = day$Date[1],
                             year = day$Year[1],
                             season = day$Season[1],
                             wind_speed = day$Wind.Speed[1],
                             temp = mean(day$Temp),
                             time = round(median(day$Time)),
                             abundance = sum(as.numeric(day$Number), na.rm = T),
                             distance = mean(as.numeric(distances)),
                             det_hist = paste(as.character(detection), collapse = "")
                             )
                }
              }
            }
          }
        }
      }
    }
  }
}

# merge our tables with point coordinates for analysis with LANDFIRE
detections <- do.call(rbind, detections) # make this into a table
  detections$Area <- detections$site
    detections <- merge(t_locations[, c("Area", "x", "y")], detections, by = "Area")
      detections <- detections[,names(detections) != "Area"]

s <- SpatialPointsDataFrame(coords=data.frame(x=detections$x,y=detections$y),data=detections[,3:ncol(detections)])
  projection(s) <- projection("+init=epsg:4326")
    s <- s[!duplicated(s@data), ] # make sure there aren't any superflous entries

# decompress our LANDFIRE data

system("unzip -o ~/Incoming/rio_mora_ld/Raster/LANDFIRE/*US_130EVH.zip -d /tmp");
system("unzip -o ~/Incoming/rio_mora_ld/Raster/LANDFIRE/*US_130EVT.zip -d /tmp");
system("unzip -o ~/Incoming/rio_mora_ld/Raster/LANDFIRE/*US_130EVC.zip -d /tmp");

cat(" -- calculating vegetation height for associations\n")

veg_height <- raster("/tmp/US_130EVH/us_130evh");
  if(!file.exists("grass_height.tif")){
    grass_height <- veg_height
      grass_height[grass_height<101] <- NA
      grass_height[grass_height>103] <- NA
        grass_height <- ((2*(grass_height-100)*0.25))-0.25 # units are meters
    writeRaster(grass_height,"grass_height.tif",overwrite=T)
  } else {
    grass_height <- raster("grass_height.tif")
  }
  if(!file.exists("shrub_height.tif")){
    shrub_height <- veg_height
      shrub_height[shrub_height<104] <- NA
      shrub_height[shrub_height>107] <- NA
        shrub_height[shrub_height==104] <- 0.25
        shrub_height[shrub_height==105] <- 0.75
        shrub_height[shrub_height==106] <- 2
        shrub_height[shrub_height==107] <- 3.25
    writeRaster(shrub_height,"shrub_height.tif",overwrite=T)
  } else {
    shrub_height <- raster("shrub_height.tif")
  }
  if(!file.exists("tree_height.tif")){
    tree_height <- veg_height
      tree_height[tree_height<108] <- NA
      tree_height[tree_height>111] <- NA
        tree_height[tree_height==108] <- 2.5
        tree_height[tree_height==109] <- 7.5
        tree_height[tree_height==110] <- 17.5
        tree_height[tree_height==111] <- 37.5
    writeRaster(tree_height,"tree_height.tif",overwrite=T)
  } else {
    tree_height <- raster("tree_height.tif")
  }

cat(" -- calculating vegetation % cover for associations\n")

veg_cover <- raster("/tmp/US_130EVC/us_130evc");
  if(!file.exists("grass_perc_cover.tif")){
    grass_perc_cover <- veg_cover
     grass_perc_cover[grass_perc_cover < 121] <- NA
     grass_perc_cover[grass_perc_cover > 129] <- NA
       grass_perc_cover <- ((grass_perc_cover-120)*10)+5
    writeRaster(grass_perc_cover,"grass_perc_cover.tif",overwrite=T)
  } else {
    grass_perc_cover <- raster("grass_perc_cover.tif")
  }
  if(!file.exists("shrub_perc_cover.tif")){
    shrub_perc_cover <- veg_cover
     shrub_perc_cover[shrub_perc_cover < 111] <- NA
     shrub_perc_cover[shrub_perc_cover > 119] <- NA
       shrub_perc_cover <- ((shrub_perc_cover-110)*10)+5
     writeRaster(shrub_perc_cover,"shrub_perc_cover.tif",overwrite=T)
  } else {
    shrub_perc_cover <- raster("shrub_perc_cover.tif")
  }
  if(!file.exists("tree_perc_cover.tif")){
    tree_perc_cover <- veg_cover
     tree_perc_cover[tree_perc_cover < 101] <- NA
     tree_perc_cover[tree_perc_cover > 109] <- NA
       tree_perc_cover <- ((tree_perc_cover-100)*10)+5
    writeRaster(tree_perc_cover,"tree_perc_cover.tif",overwrite=T)
  } else {
    tree_perc_cover <- raster("tree_perc_cover.tif")
  }

# DON'T BOTHER WITH PATCH METRICS FOR grasslands.
# Most of the region is "grass" and the vertical and horizontal structral characteristics are better captured
# in the veg height and % cover calculated above for associations.
#
# cat(" -- calculating vegetation patch metrics for associations\n")
# veg_type <- raster("/tmp/US_130EVT/us_130evt");
#   if(!file.exists("grass.tif")){
#     t <- veg_type@data@attributes
#       grass <- veg_type %in% t[t$EVT_PHYS == "Grassland",]$ID
#     writeRaster(grass,"grass.tif")
#   } else {
#     grass <- raster("grass.tif")
#   }

if(!file.exists("rio_mora_monitoring_data_processed.csv")){
  # prepare buffers around sites
  buffers <- vector('list', 30) # create enough space for a range of buffers from 100m -> 3000m
        s <- spTransform(s, CRS(projection(veg_height))); # use a consistent CRS in meters

  for(i in 1:length(buffers)){
    buffers[[i]] <- gBuffer(s, byid = T, width = i*100, capStyle = "square");
  }
  # extract
  for(i in 1:length(buffers)){
    buffers[[i]]$tree_perc_cover  <- suppressWarnings(raster::extract(tree_perc_cover,buffers[[i]],fun=mean,na.rm=T)); cat(".")
    buffers[[i]]$shrub_perc_cover <- suppressWarnings(raster::extract(shrub_perc_cover,buffers[[i]],fun=mean,na.rm=T)); cat(".")
    buffers[[i]]$grass_perc_cover <- suppressWarnings(raster::extract(grass_perc_cover,buffers[[i]],fun=mean,na.rm=T)); cat(".")

    buffers[[i]]$tree_height  <- suppressWarnings(raster::extract(tree_height,buffers[[i]],fun=mean,na.rm=T)); cat(".")
    buffers[[i]]$shrub_height <- suppressWarnings(raster::extract(shrub_height,buffers[[i]],fun=mean,na.rm=T)); cat(".")
    buffers[[i]]$grass_height <- suppressWarnings(raster::extract(grass_height,buffers[[i]],fun=mean,na.rm=T)); cat(".")

    names(buffers[[i]]@data) <- c(names(buffers[[i]]@data[,1:11]),paste(names(buffers[[i]]@data[,12:17]),i,sep="_"))
    cat(paste("(",i,"/",length(buffers),")",sep=""))
  }; cat("\n");
  # merge into a single table
  merged <- buffers[[1]]@data
    for(i in 2:length(buffers)){
      merged <- cbind(merged,buffers[[i]]@data[,12:17])
    }
      s@data <- merged
  # write to disk
  write.csv(s@data,"rio_mora_monitoring_data_processed.csv",row.names=F) # too many variables to store as a spatial object
} else {
  t <- read.csv("rio_mora_monitoring_data_processed.csv")
    s@data <- t
}

# train some preliminary random forests to each species
t <- t[,!grepl(names(t),pattern="site|date|det_hist|distance")]

models <- list()
for(spp in unique(t$spp)){
  training <- t[t$spp == spp,] # figure out our non-zero observations
  absences <- t[t$spp != spp,]
      absences$abundance <- 0
  training <- rbind(training,absences)
    training <- training[,!grepl(names(training),pattern="spp")]
      training[is.na(training)] <- 0

  models[[length(models)+1]] <- randomForest(abundance~.,data=training,ntree=8000)
}

names(models) <- as.character(unique(t$spp))

# Beyond summary statistics, these data don't look that promising for fitting any models
# sort(importance(models$LASP)[,1]/max(importance(models$LASP)[,1]),decreasing=T)
# sort(importance(models$PIJA)[,1]/max(importance(models$PIJA)[,1]),decreasing=T)
