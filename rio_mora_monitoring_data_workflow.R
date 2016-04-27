require(raster)
require(rgdal)
require(rgeos)
require(landscapeAnalysis)

priority_spp_four_letter_codes <- c(
  "LOSH", # loggerhead shrike
  "LASP", # lark sparrow
  "PRFA", # prairie falcon
  "PIJA", # pinyon jay
  "LBCU", # long-billed curlew
  "BUOW", # burrowing owl
  "CASP", # cassin's sparrow
  "CCLO", # chestnut-collared longspur
  "BEVI", # bell's vireo
  "GOEA", # golden eagle
  "FEHA", # ferruginous hawk
  "GRSP", # grasshopper sparrow
  "LEWO", # lewis's woodpecker
  "MOPL", # mountain plover
  "NOHA", # northern harrier
  "SPOW", # spotted owl (Mexican)
  "PEFA", # peregrin falcon
  "SWFL", # southwestern willow flycatcher
  "SWHA", # swainson’s hawk
  "YEWA"  # yellow warbler
  )

t_monitoring_data <- read.csv("~/Incoming/rio_mora_ld/rio_mora_transect_data.csv")
      t_locations <- read.csv("~/Incoming/rio_mora_ld/rio_mora_transect_locations.csv")

s <- SpatialPointsDataFrame(coords=data.frame(x=t_locations$x,y=t_locations$y), data=data.frame(site=as.vector(t_locations$New.Name)))
  projection(s) <- projection("+init=epsg:4326")

t_monitoring_data <- t_monitoring_data[!is.na(as.numeric(as.vector(t_monitoring_data$Minute.Interval))),] # strip out observations not taken during a minute interval

veg_height <- raster("/tmp/US_130EVH/us_130evh")
  veg_type <- raster("/tmp/US_130EVT/us_130evt")
 veg_cover <- raster("/tmp/US_130EVC/us_130evc")

  s_3x3 <- gBuffer(s,byid=T,width=45, capStyle="square")
  s_6x6 <- gBuffer(s,byid=T,width=90, capStyle="square")
s_12x12 <- gBuffer(s,byid=T,width=180, capStyle="square")
s_24x24 <- gBuffer(s,byid=T,width=360, capStyle="square")

# report summary statistics for number of priority species detections
cat(" -- number of priority species detected in monitoring data: ",
    length(unique(as.vector(t_monitoring_data$X4.letter.code)[(toupper(as.vector(t_monitoring_data$X4.letter.code)) %in% priority_spp_four_letter_codes)])),
    "/",
    length(priority_spp_four_letter_codes),
    "\n",sep="")

t <- sort(table(as.vector(t_monitoring_data$X4.letter.code)[(toupper(as.vector(t_monitoring_data$X4.letter.code)) %in% priority_spp_four_letter_codes)]))
t_missing <- priority_spp_four_letter_codes[!(priority_spp_four_letter_codes %in% (toupper(as.vector(t_monitoring_data$X4.letter.code))))]
  a <- as.table(rep(0,length(t_missing)))
    names(a) <- t_missing
      t_missing <- a
        t <- c(t,t_missing);
        
rm(t_missing,a);

p <- barplot(t,xaxt="n",ylab="detections",xlab="priority species",main="Rio Mora Bird Monitoring Data : 2012-14")
  text(cex=0.75, x=p+0.65, y=-0.45, names(t), xpd=TRUE, srt=45, pos=2)
    grid(lwd=1.2)

# parse our monitoring data so we are only observing data for priority species
