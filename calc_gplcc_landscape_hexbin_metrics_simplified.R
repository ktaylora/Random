require(raster)
require(landscapeAnalysis)
require(rgdal)

HOME <- "/home/ktaylora"
argv <- commandArgs(trailingOnly=T)

s    <- readOGR(paste(Sys.getenv("HOME"),"/PLJV/boundaries/GPLCC Pilot/GPLCC_Pilot_Region/",sep=""),"pilot_region_hexagonal_units",verbose=F)
r    <- raster(paste(Sys.getenv("HOME"),"/PLJV/landcover/orig/Final_LC_8bit",".tif",sep=""))
#r <- raster("/home/ktaylora/Projects/large_block_sliding_window_analysis/plv_landcover_gplcc_ag_forecast_2025.tif")

s <- s[seq(as.numeric(argv[1]),as.numeric(argv[2]),1),]
o<-landscapeAnalysis::cropRasterByPolygons(r,s,parallel=F)

agricultural_habitat <- c(38,202,201,203,204,205,206,207,208,209,210) # Cropland

shrubland_habitat <- c(82, # Shrubland
                       83, # Mesquite
                       81, # Savannah
                       85, # Shinnery
                       87) # Sand sage

grass_habitat <- c(39, # CRP
                   31, # CRP - GRASS
                   37, # Pasture
                   71, # Mixedgrass
                   75) # Shortgrass

grass_units <- lReclass(o,inValues=grass_habitat)
  grass_stats <- lCalculateLandscapeMetrics(grass_units,metric=c("total.area","mean.patch.area"))

ag_units <- lReclass(o,inValues=agricultural_habitat)
  ag_stats <- lCalculateLandscapeMetrics(ag_units,metric=c("total.area","mean.patch.area"))

sh_units <- lReclass(o,inValues=shrubland_habitat)
  sh_stats <- lCalculateLandscapeMetrics(sh_units,metric=c("total.area","mean.patch.area"))

t <- cbind(grass_stats,
           ag_stats[,!grepl(names(ag_stats),pattern="id")],
           sh_stats[,!grepl(names(sh_stats),pattern="id")])

names(t) <- c("id","total.area_grass","mean.patch.area_grass","total.area_ag","mean.patch.area_ag","total.area_shrub","mean.patch.area_shrub")

t$id <- s$id

write.csv(t,paste("results_",argv[1],"-",argv[2],".csv",sep=""),row.names=F)
