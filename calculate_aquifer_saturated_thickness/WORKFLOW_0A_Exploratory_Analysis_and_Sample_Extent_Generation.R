require(raster)
require(rgdal)

plot.saturatedThickness <- function(s, m, lat_quantile=NULL, lon_quantile=NULL, main=""){
  plot(s$strtd_t ~ s$year, pch="-", cex=0.85, col="white", xlab="Year", ylab="Saturated Thickness (ft)", main = main)
    grid(); grid();
      points(s$lv_v_ft ~ s$year, pch="-", cex=0.85, col="DarkGrey")
      years <- seq(min(s$year), max(s$year), 1)
      demo <- data.frame(year = years,
                          lon = rep(quantile(s$lon, ifelse(is.null(lon_quantile), 0.5, lon_quantile)), length(years)),
                          lat = rep(quantile(s$lat, ifelse(is.null(lat_quantile), 0.5, lat_quantile)), length(years))
              )
      p <- predict(m, newdata = demo)
        lines(y = p, x = years, col = "DarkRed")
}

results <- new.env()

file.copy("processed_well_pts_all_years.zip", to = "/tmp", overwrite = T)
setwd("/tmp")
s <- list.files(pattern = "shp$")
  s <- s[grepl(s, pattern = "19|20")]
    s <- gsub(s, pattern="[.]shp", replacement = "")

years <- as.numeric(s)
    s <- lapply(s,readOGR,dsn=".",verbose=F)

# assign a year to points in well dataset, then merge into a single shapefile
for (i in 1:length(years)){
  s[[i]]$year <- years[i]
}
s <- do.call(rbind, s)
  s$lat <- s@coords[, 2]
    s$lon <- s@coords[, 1]

grid <- raster(extent(spTransform(s,CRS(projection("+init=epsg:2163")))), crs=CRS(projection("+init=epsg:2163")),res=90)

training <- data.frame(thickness = s$strtd_t, year = s$year, lat = s$lat, lon = s$lon)
       m <- loess(thickness ~ year + lat + lon, data=training)
# plot some trends for northern and southern extremes of the highplains region
par(mfrow=c(2,2))

plot.saturatedThickness(s=s, m=m, lat_quantile=0.9, main="Northern Region")
plot.saturatedThickness(s=s, m=m, lat_quantile=0.1, main="Southern Region")
plot.saturatedThickness(s=s, m=m, lon_quantile=0.1, main="Western Region")
plot.saturatedThickness(s=s, m=m, lon_quantile=0.9, main="Eastern Region")

# dump a copy of the local regression model -- because it takes forever to fit
assign(x = "m_loess", value = m, envir = results)
save(file="session.Rdata", list=ls(results), envir=results)
# grid overlapping points to define a base layer for deriving the extent of our study region
raster()
