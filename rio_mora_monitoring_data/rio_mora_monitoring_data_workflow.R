#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2016]
#

require(raster)
require(rgdal)
require(rgeos)
require(landscapeAnalysis)

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
  "YEWA"  # yellow warbler
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

system("unzip -o ~/Incoming/rio_mora_ld/LANDFIRE/*US_130EVH.zip -d /tmp");
system("unzip -o ~/Incoming/rio_mora_ld/LANDFIRE/*US_130EVT.zip -d /tmp");
system("unzip -o ~/Incoming/rio_mora_ld/LANDFIRE/*US_130EVC.zip -d /tmp");

veg_height <- raster("/tmp/US_130EVH/us_130evh");
  veg_type <- raster("/tmp/US_130EVT/us_130evt");
 veg_cover <- raster("/tmp/US_130EVC/us_130evc");

# prepare buffers around sites
      s <- spTransform(s, CRS(projection(veg_height))); # use a consistent CRS in meters

  s_3x3 <- gBuffer(s, byid = T, width = 45, capStyle = "square");
  s_6x6 <- gBuffer(s, byid = T, width = 90, capStyle = "square");
s_12x12 <- gBuffer(s, byid = T, width = 180, capStyle = "square");
s_24x24 <- gBuffer(s, byid = T, width = 360, capStyle = "square");
