require(rgdal)
require(rgeos)
require(ggplot2)
require(reshape2)
require(gridExtra)
require(grid)
require(raster)

HOME <- Sys.getenv("HOME")

suitability_projections <- list.files('/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/raster/suitability_projections',
                                      pattern="tif$", full.names=T)
             boundaries <- readOGR(paste(HOME,"/Products/boundaries/",sep=""), "western_north_american_boundaries",verbose=F)

# read-in our contours for current conditions
ens_tridentata_current   <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "tridentata_current_climate_ens",verbose=F)
ens_wyomingensis_current <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "wyomingensis_current_climate_ens",verbose=F)
ens_vaseyana_current     <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "vaseyana_current_climate_ens",verbose=F)


# Figure 1 - Climate Suitabilty for subspecies under current conditions
# dev.new(height=5,width=10)
# par(mfrow=c(1,3))
#
# plot(ens_tridentata_current,col="#00990050", border=NA, main="ssp. tridentata");
#   plot(boundaries, border=rgb(0, 0, 0, 0.5),add=T);
# box(); grid(lty=1,col="#00000030")
# plot(ens_wyomingensis_current,col="#00990050", border=NA,main="ssp. wyomingensis");
#   plot(boundaries, border=rgb(0, 0, 0, 0.5),add=T);
# box(); grid(lty=1,col="#00000030")
# plot(ens_vaseyana_current,col="#00990050", border=NA,main="ssp. vaseyana");
#   plot(boundaries, border=rgb(0, 0, 0, 0.5),add=T);
# box(); grid(lty=1,col="#00000030")

# Figure 2 - Plot unique regions and consensus amoung ssp

# Calculate Agreement / Disagreement Vector Surfaces

intersect_p50_current <- rgeos::gIntersection(rgeos::gIntersection(ens_tridentata_current,ens_wyomingensis_current)@polyobj,ens_vaseyana_current)@polyobj
intersect_p50_wyo_tri_current <- rgeos::gIntersection(ens_tridentata_current,ens_wyomingensis_current)@polyobj

unique_p50_current_wyomingensis <- rgeos::gDifference(ens_wyomingensis_current,ens_tridentata_current) # wyomingensis has greater range
  unique_p50_current_wyomingensis <- rgeos::gDifference(unique_p50_current_wyomingensis, ens_vaseyana_current)

unique_p50_current_tridentata <- rgeos::gDifference(ens_tridentata_current,ens_wyomingensis_current)
    unique_p50_current_tridentata <- rgeos::gDifference(unique_p50_current_tridentata, ens_vaseyana_current)

unique_p50_current_vaseyana <- rgeos::gDifference(ens_vaseyana_current,ens_wyomingensis_current)
    unique_p50_current_vaseyana <- rgeos::gDifference(unique_p50_current_vaseyana, ens_tridentata_current)


png(filename="Fig_2.png",width=941.473632,height=517.05264)
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
text("topleft", "A")

plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.8)
  plot(spTransform(intersect_p50_wyo_tri_current,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
  plot(spTransform(unique_p50_current_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
  plot(spTransform(unique_p50_current_tridentata,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
  plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), border="#000000B3",add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","wyomingensis","tridentata"), cex=0.8, fill=c("#003D7A","#005CB8","#B2D1F0"),bg = "white");
text("topleft", "B")
graphics.off()

# Figure 3 - Climate Suitabilty for subspecies under future conditions

# read-in our contours for future conditions

# 2050 RCP 4.5
ens_tridentata_2050_rcp_45   <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "2050_45_tri_ens",verbose=F)
ens_wyomingensis_2050_rcp_45 <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "2050_45_wyo_ens",verbose=F)
ens_vaseyana_2050_rcp_45     <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "2050_45_vas_ens",verbose=F)

intersect_p50_2050_rcp_45 <- rgeos::gIntersection(rgeos::gIntersection(ens_tridentata_2050_rcp_45,ens_wyomingensis_2050_rcp_45)@polyobj,ens_vaseyana_2050_rcp_45)@polyobj
intersect_p50_wyo_tri_2050_rcp_45 <- rgeos::gIntersection(ens_tridentata_2050_rcp_45,ens_wyomingensis_2050_rcp_45)@polyobj

unique_p50_2050_rcp_45_wyomingensis <- rgeos::gDifference(ens_wyomingensis_2050_rcp_45,ens_tridentata_2050_rcp_45) # wyomingensis has greater range
 unique_p50_2050_rcp_45_wyomingensis <- rgeos::gDifference(unique_p50_2050_rcp_45_wyomingensis, ens_vaseyana_2050_rcp_45)

unique_p50_2050_rcp_45_tridentata <- rgeos::gDifference(ens_tridentata_2050_rcp_45,ens_wyomingensis_2050_rcp_45)
   unique_p50_2050_rcp_45_tridentata <- rgeos::gDifference(unique_p50_2050_rcp_45_tridentata, ens_vaseyana_2050_rcp_45)

unique_p50_2050_rcp_45_vaseyana <- rgeos::gDifference(ens_vaseyana_2050_rcp_45,ens_wyomingensis_2050_rcp_45)
   unique_p50_2050_rcp_45_vaseyana <- rgeos::gDifference(unique_p50_2050_rcp_45_vaseyana, ens_tridentata_2050_rcp_45)

png(filename="Fig_3.png",width=941.473632,height=517.05264)
par(mfrow=c(1,2))
par(mar=par()$mar/2)
plot(main=NA,spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
 plot(spTransform(intersect_p50_2050_rcp_45,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
 plot(spTransform(unique_p50_2050_rcp_45_tridentata,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
 plot(spTransform(unique_p50_2050_rcp_45_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#4D94DB",border=NA,add=T)
 plot(spTransform(unique_p50_2050_rcp_45_vaseyana,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
 plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","tridentata","wyomingensis","vaseyana"), cex=0.8, fill=c("#003D7A","#005CB8","#4D94DB","#B2D1F0"),bg = "white");
text(x=-1500000,y=1300000,"2050 (RCP 4.5)",cex=0.85)
plot(main=NA,spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
  plot(spTransform(intersect_p50_wyo_tri_2050_rcp_45,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
  plot(spTransform(unique_p50_2050_rcp_45_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
  plot(spTransform(unique_p50_2050_rcp_45_tridentata,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
  plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), border="#000000B3",add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","wyomingensis","tridentata"), cex=0.8, fill=c("#003D7A","#005CB8","#B2D1F0"),bg = "white");
text(x=-1500000,y=1300000,"2050 (RCP 4.5)",cex=0.85)
graphics.off()

# 2070 RCP 4.5
ens_tridentata_2070_rcp_45   <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                "2070_45_tri_ens",verbose=F)
ens_wyomingensis_2070_rcp_45 <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                "2070_45_wyo_ens",verbose=F)
ens_vaseyana_2070_rcp_45     <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                "2070_45_vas_ens",verbose=F)

intersect_p50_2070_rcp_45 <- rgeos::gIntersection(rgeos::gIntersection(ens_tridentata_2070_rcp_45,ens_wyomingensis_2070_rcp_45)@polyobj,ens_vaseyana_2070_rcp_45)@polyobj
intersect_p50_wyo_tri_2070_rcp_45 <- rgeos::gIntersection(ens_tridentata_2070_rcp_45,ens_wyomingensis_2070_rcp_45)@polyobj

unique_p50_2070_rcp_45_wyomingensis <- rgeos::gDifference(ens_wyomingensis_2070_rcp_45,ens_tridentata_2070_rcp_45) # wyomingensis has greater range
 unique_p50_2070_rcp_45_wyomingensis <- rgeos::gDifference(unique_p50_2070_rcp_45_wyomingensis, ens_vaseyana_2070_rcp_45)

unique_p50_2070_rcp_45_tridentata <- rgeos::gDifference(ens_tridentata_2070_rcp_45,ens_wyomingensis_2070_rcp_45)
   unique_p50_2070_rcp_45_tridentata <- rgeos::gDifference(unique_p50_2070_rcp_45_tridentata, ens_vaseyana_2070_rcp_45)

unique_p50_2070_rcp_45_vaseyana <- rgeos::gDifference(ens_vaseyana_2070_rcp_45,ens_wyomingensis_2070_rcp_45)
   unique_p50_2070_rcp_45_vaseyana <- rgeos::gDifference(unique_p50_2070_rcp_45_vaseyana, ens_tridentata_2070_rcp_45)

# Figure 4
png(filename="Fig_4.png",width=941.473632,height=517.05264)
par(mfrow=c(1,2))
par(mar=par()$mar/2)
plot(main=NA,spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
 plot(spTransform(intersect_p50_2070_rcp_45,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
 plot(spTransform(unique_p50_2070_rcp_45_tridentata,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
 plot(spTransform(unique_p50_2070_rcp_45_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#4D94DB",border=NA,add=T)
 plot(spTransform(unique_p50_2070_rcp_45_vaseyana,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
 plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","tridentata","wyomingensis","vaseyana"), cex=0.8, fill=c("#003D7A","#005CB8","#4D94DB","#B2D1F0"),bg = "white");
text(x=-1500000,y=1300000,"2070 (RCP 4.5)",cex=0.85)
plot(main=NA,spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
  plot(spTransform(intersect_p50_wyo_tri_2070_rcp_45,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
  plot(spTransform(unique_p50_2070_rcp_45_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
  plot(spTransform(unique_p50_2070_rcp_45_tridentata,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
  plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), border="#000000B3",add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","wyomingensis","tridentata"), cex=0.8, fill=c("#003D7A","#005CB8","#B2D1F0"),bg = "white");
text(x=-1500000,y=1300000,"2070 (RCP 4.5)",cex=0.85)
graphics.off()

# 2050 RCP 8.5
ens_tridentata_2050_rcp_85   <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "2050_85_tri_ens",verbose=F)
ens_wyomingensis_2050_rcp_85 <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "2050_85_wyo_ens",verbose=F)
ens_vaseyana_2050_rcp_85     <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                 "2050_85_vas_ens",verbose=F)

intersect_p50_2050_rcp_85 <- rgeos::gIntersection(rgeos::gIntersection(ens_tridentata_2050_rcp_85,ens_wyomingensis_2050_rcp_85)@polyobj,ens_vaseyana_2050_rcp_85)@polyobj
intersect_p50_wyo_tri_2050_rcp_85 <- rgeos::gIntersection(ens_tridentata_2050_rcp_85,ens_wyomingensis_2050_rcp_85)@polyobj

unique_p50_2050_rcp_85_wyomingensis <- rgeos::gDifference(ens_wyomingensis_2050_rcp_85,ens_tridentata_2050_rcp_85) # wyomingensis has greater range
 unique_p50_2050_rcp_85_wyomingensis <- rgeos::gDifference(unique_p50_2050_rcp_85_wyomingensis, ens_vaseyana_2050_rcp_85)

unique_p50_2050_rcp_85_tridentata <- rgeos::gDifference(ens_tridentata_2050_rcp_85,ens_wyomingensis_2050_rcp_85)
   unique_p50_2050_rcp_85_tridentata <- rgeos::gDifference(unique_p50_2050_rcp_85_tridentata, ens_vaseyana_2050_rcp_85)

unique_p50_2050_rcp_85_vaseyana <- rgeos::gDifference(ens_vaseyana_2050_rcp_85,ens_wyomingensis_2050_rcp_85)
   unique_p50_2050_rcp_85_vaseyana <- rgeos::gDifference(unique_p50_2050_rcp_85_vaseyana, ens_tridentata_2050_rcp_85)

# Figure 5
png(filename="Fig_5.png",width=941.473632,height=517.05264)
par(mfrow=c(1,2))
par(mar=par()$mar/2)
plot(main=NA,spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
 plot(spTransform(intersect_p50_2050_rcp_85,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
 plot(spTransform(unique_p50_2050_rcp_85_tridentata,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
 plot(spTransform(unique_p50_2050_rcp_85_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#4D94DB",border=NA,add=T)
 plot(spTransform(unique_p50_2050_rcp_85_vaseyana,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
 plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","tridentata","wyomingensis","vaseyana"), cex=0.8, fill=c("#003D7A","#005CB8","#4D94DB","#B2D1F0"),bg = "white");
text(x=-1500000,y=1300000,"2050 (RCP 8.5)",cex=0.85)
plot(main=NA,spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
  plot(spTransform(intersect_p50_wyo_tri_2050_rcp_85,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
  plot(spTransform(unique_p50_2050_rcp_85_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
  plot(spTransform(unique_p50_2050_rcp_85_tridentata,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
  plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), border="#000000B3",add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","wyomingensis","tridentata"), cex=0.8, fill=c("#003D7A","#005CB8","#B2D1F0"),bg = "white");
text(x=-1500000,y=1300000,"2050 (RCP 8.5)",cex=0.85)
graphics.off()

# 2070 RCP 8.5
ens_tridentata_2070_rcp_85   <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                "2070_85_tri_ens",verbose=F)
ens_wyomingensis_2070_rcp_85 <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                "2070_85_wyo_ens",verbose=F)
ens_vaseyana_2070_rcp_85     <- readOGR("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/vector/suitability_contours",
                                "2070_85_vas_ens",verbose=F)

intersect_p50_2070_rcp_85 <- rgeos::gIntersection(rgeos::gIntersection(ens_tridentata_2070_rcp_85,ens_wyomingensis_2070_rcp_85)@polyobj,ens_vaseyana_2070_rcp_85)@polyobj
intersect_p50_wyo_tri_2070_rcp_85 <- rgeos::gIntersection(ens_tridentata_2070_rcp_85,ens_wyomingensis_2070_rcp_85)@polyobj

unique_p50_2070_rcp_85_wyomingensis <- rgeos::gDifference(ens_wyomingensis_2070_rcp_85,ens_tridentata_2070_rcp_85) # wyomingensis has greater range
 unique_p50_2070_rcp_85_wyomingensis <- rgeos::gDifference(unique_p50_2070_rcp_85_wyomingensis, ens_vaseyana_2070_rcp_85)

unique_p50_2070_rcp_85_tridentata <- rgeos::gDifference(ens_tridentata_2070_rcp_85,ens_wyomingensis_2070_rcp_85)
   unique_p50_2070_rcp_85_tridentata <- rgeos::gDifference(unique_p50_2070_rcp_85_tridentata, ens_vaseyana_2070_rcp_85)

unique_p50_2070_rcp_85_vaseyana <- rgeos::gDifference(ens_vaseyana_2070_rcp_85,ens_wyomingensis_2070_rcp_85)
   unique_p50_2070_rcp_85_vaseyana <- rgeos::gDifference(unique_p50_2070_rcp_85_vaseyana, ens_tridentata_2070_rcp_85)

# Figure 6
png(filename="Fig_6.png",width=941.473632,height=517.05264)
par(mfrow=c(1,2))
par(mar=par()$mar/2)
plot(main=NA,spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
 plot(spTransform(intersect_p50_2070_rcp_85,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
 plot(spTransform(unique_p50_2070_rcp_85_tridentata,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
 plot(spTransform(unique_p50_2070_rcp_85_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#4D94DB",border=NA,add=T)
 plot(spTransform(unique_p50_2070_rcp_85_vaseyana,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
 plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","tridentata","wyomingensis","vaseyana"), cex=0.8, fill=c("#003D7A","#005CB8","#4D94DB","#B2D1F0"),bg = "white");
text(x=-1500000,y=1300000,"2070 (RCP 8.5)",cex=0.85)
plot(main=NA,spTransform(boundaries,CRS(projection("+init=epsg:2163"))),col="white",border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
  plot(spTransform(intersect_p50_wyo_tri_2070_rcp_85,CRS(projection("+init=epsg:2163"))), col="#003D7A",border=NA,add=T)
  plot(spTransform(unique_p50_2070_rcp_85_wyomingensis,CRS(projection("+init=epsg:2163"))), col="#005CB8",border=NA,add=T)
  plot(spTransform(unique_p50_2070_rcp_85_tridentata,CRS(projection("+init=epsg:2163"))), col="#B2D1F0",border=NA,add=T)
  plot(spTransform(boundaries,CRS(projection("+init=epsg:2163"))), border="#000000B3",add=T)
box(); grid(lty=1,col="#00000030")
legend("topright", c("consensus","wyomingensis","tridentata"), cex=0.8, fill=c("#003D7A","#005CB8","#B2D1F0"),bg = "white");
text(x=-1500000,y=1300000,"2070 (RCP 8.5)",cex=0.85)
graphics.off()

# Figure 7 -- step plots for latitude / elevation under current climate conditions

png(filename="Fig_7.png",width=941.473632,height=517.05264)

grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}

elevation <- raster("/media/ktaylora/80C0-4F29/lf_elev1km/w001001x.adf")

p_current_elev_agreement <- rasterToPoints(rasterize(spTransform(spsample(intersect_p50_current,n=15000,type="random"),CRS(projection(elevation))),elevation),spatial=T)
  p_current_elev_agreement <- extract(elevation,p_current_elev_agreement,df=T)[,2]
p_current_elev_tri <- rasterToPoints(rasterize(spTransform(spsample(ens_tridentata_current,n=15000,type="random"),CRS(projection(elevation))),elevation),spatial=T)
  p_current_elev_tri <- extract(elevation,p_current_elev_tri,df=T)[,2]
p_current_elev_wyo <- rasterToPoints(rasterize(spTransform(spsample(ens_wyomingensis_current,n=15000,type="random"),CRS(projection(elevation))),elevation),spatial=T)
  p_current_elev_wyo <- extract(elevation,p_current_elev_wyo,df=T)[,2]
p_current_elev_vas <- rasterToPoints(rasterize(spTransform(spsample(ens_vaseyana_current,n=15000,type="random"),CRS(projection(elevation))),elevation),spatial=T)
  p_current_elev_vas <- extract(elevation,p_current_elev_vas,df=T)[,2]

df <- data.frame(agreement=sample(p_current_elev_agreement,5000,replace=T),
                 tridentata=sample(p_current_elev_tri,5000,replace=T),
                 vaseyana=sample(p_current_elev_vas,5000,replace=T),
                 wyomingensis=sample(p_current_elev_wyo,5000))
df <- melt(df,variable.name="subspecies",value.name="elevation")
p1 <- ggplot(df, aes(factor(subspecies), elevation)) +
             geom_boxplot(aes(x='agreement', y=elevation, color="Agreement")) +
             geom_boxplot(aes(x='tridentata', y=elevation, color="Tridentata")) +
             geom_boxplot(aes(x='vaseyana', y=elevation,color="Vaseyana")) +
             geom_boxplot(aes(x='wyomingensis', y=elevation,color="Wyomingensis")) +
             #geom_vline(xintercept=mean(sample(lat_elev_current_glm_90@coords[,2],500)), colour="grey", linetype = "longdash") +
             xlab("Subspecies") + ylab("Elevation") +
             theme_bw() + theme(legend.title=element_blank())

p_current_lat_agreement <- spsample(intersect_p50_current,n=15000,type="random")@coords[,2]
p_current_lat_tri <- spsample(ens_tridentata_current,n=15000,type="random")@coords[,2]
p_current_lat_wyo <- spsample(ens_wyomingensis_current,n=15000,type="random")@coords[,2]
p_current_lat_vas <- spsample(ens_vaseyana_current,n=15000,type="random")@coords[,2]

df <- data.frame(agreement=sample(p_current_lat_agreement,5000,replace=T),
                 tridentata=sample(p_current_lat_tri,5000,replace=T),
                 vaseyana=sample(p_current_lat_wyo,5000,replace=T),
                 wyomingensis=sample(p_current_lat_wyo,5000))

df <- melt(df,variable.name="subspecies",value.name="latitude")

p2 <- ggplot(df, aes(factor(subspecies), latitude)) +
            geom_boxplot(aes(x='agreement', y=latitude, color="Agreement")) +
            geom_boxplot(aes(x='tridentata', y=latitude, color="Tridentata")) +
            geom_boxplot(aes(x='vaseyana', y=latitude,color="Vaseyana")) +
            geom_boxplot(aes(x='wyomingensis', y=latitude,color="Wyomingensis")) +
            #geom_vline(xintercept=mean(sample(lat_elev_current_glm_90@coords[,2],500)), colour="grey", linetype = "longdash") +
            xlab("Subspecies") + ylab("Latitude") +
            theme_bw() + theme(legend.title=element_blank())

grid_arrange_shared_legend(p1, p2)
graphics.off()

## Prepare current climate conditions for boxplots

current_climate <- list.files("/media/ktaylora/big_black/intermediates/weather/worldclim/current",pattern="bil$",full.names=T)
  current_climate <- raster::stack(current_climate[grepl(current_climate,pattern="bio_1[.]|bio_11[.]|bio_12[.]|bio_15[.]|bio_18[.]")])

# ex: scenario="RCP_45_2050", subset="gain", value="25"
data.frame(scenario=NA,subset=NA,value=NA)

rcp_45_2050_climate_gf
rcp_45_2050_climate_in
rcp_45_2050_climate_ip
rcp_85_2050_climate <- list.files("/media/ktaylora/big_black/intermediates/weather/sagebrush_subspp_future_conditions/focal_for_ssp_manuscript",pattern="bil$",full.names=T)
rcp_45_2070_climate <- list.files("/media/ktaylora/big_black/intermediates/weather/sagebrush_subspp_future_conditions/focal_for_ssp_manuscript",pattern="bil$",full.names=T)
rcp_85_2070_climate <- list.files("/media/ktaylora/big_black/intermediates/weather/sagebrush_subspp_future_conditions/focal_for_ssp_manuscript",pattern="bil$",full.names=T)

## Figure 8 -- difference plots for each subspecies
t_crs <- CRS(projection("+init=epsg:2163"))

# ssp. tridentata 2050 (4.5)

gain <- rgeos::gDifference(ens_tridentata_2050_rcp_45,ens_tridentata_current)
loss <- rgeos::gDifference(ens_tridentata_current,ens_tridentata_2050_rcp_45)

png(file="Desktop/2050_45_difference_plots.png",width=1273.488192,height=469.516512)
par(mfrow=c(1,3))

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_tridentata_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030")
text("A",cex=1.4,y=1300000,x=-1700000)

# extract for boxplots
#
# focalVars <- c("01[.]tif|11[.]tif|12[.]tif|15[.]tif|18[.]tif")
# zips_path <- list.files("/media/ktaylora/big_black/intermediates/weather/sagebrush_subspp_future_conditions/focal_for_ssp_manuscript",pattern="zip$",full.names=T)
# rcp_45_2050_climate_ac <- zips_path[grepl(zips_path,pattern="/ac*.*45bi50*")]
#   names <- unzip(rcp_45_2050_climate_ac,list=T)$Name;
#     names <- names[grepl(names,pattern=focalVars)];
#   unlink("/tmp/focal_zip",force=T,recursive=T);
#   unzip(rcp_45_2050_climate_ac,files=names,overwrite=T,exdir="/tmp/focal_zip")
#     rcp_45_2050_climate_ac <- raster::stack(list.files("/tmp/focal_zip",pattern="tif$",full.names=T));


# species GAP records
wyo_gap <- spTransform(readOGR("/home/ktaylora/Products/uw/big_sagebrush_subspp_analysis/vectors","wyomingensis_gap_records",verbose=F),t_crs)
  wyo_gap <- wyo_gap[wyo_gap$resp==1,]
tri_gap <- spTransform(readOGR("/home/ktaylora/Products/uw/big_sagebrush_subspp_analysis/vectors","tridentata_gap_records",verbose=F),t_crs)
  tri_gap <- tri_gap[tri_gap$resp==1,]
vas_gap <- spTransform(readOGR("/home/ktaylora/Products/uw/big_sagebrush_subspp_analysis/vectors","vaseyana_gap_records",verbose=F),t_crs)
  vas_gap <- vas_gap[vas_gap$resp==1,]

# ssp. wyomingensis 2050 (4.5)

gain <- rgeos::gDifference(ens_wyomingensis_2050_rcp_45,ens_wyomingensis_current)
loss <- rgeos::gDifference(ens_wyomingensis_current,ens_wyomingensis_2050_rcp_45)

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_wyomingensis_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(wyo_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("B",cex=1.4,y=1300000,x=-1700000)

# ssp. vaseyana 2050 (4.5)

gain <- rgeos::gDifference(ens_vaseyana_2050_rcp_45,ens_vaseyana_current)
loss <- rgeos::gDifference(ens_vaseyana_current,ens_vaseyana_2050_rcp_45)

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_vaseyana_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(vas_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("C",cex=1.4,y=1300000,x=-1700000)

graphics.off()

# ssp. tridentata 2050 (8.5)

gain <- rgeos::gDifference(ens_tridentata_2050_rcp_85,ens_tridentata_current)
loss <- rgeos::gDifference(ens_tridentata_current,ens_tridentata_2050_rcp_85)

png(file="Desktop/2050_85_difference_plots.png",width=1273.488192,height=469.516512)
par(mfrow=c(1,3))

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_tridentata_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(tri_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("D",cex=1.4,y=1300000,x=-1700000)

# ssp. wyomingensis 2050 (8.5)

gain <- rgeos::gDifference(ens_wyomingensis_2050_rcp_85,ens_wyomingensis_current)
loss <- rgeos::gDifference(ens_wyomingensis_current,ens_wyomingensis_2050_rcp_85)

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_wyomingensis_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(wyo_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("E",cex=1.4,y=1300000,x=-1700000)

# ssp. vaseyana 2050 (8.5)

gain <- rgeos::gDifference(ens_vaseyana_2050_rcp_85,ens_vaseyana_current)
loss <- rgeos::gDifference(ens_vaseyana_current,ens_vaseyana_2050_rcp_85)

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_vaseyana_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(vas_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("F",cex=1.4,y=1300000,x=-1700000)

graphics.off()

# ssp. tridentata 2070 (4.5)

gain <- rgeos::gDifference(ens_tridentata_2070_rcp_45,ens_tridentata_current)
loss <- rgeos::gDifference(ens_tridentata_current,ens_tridentata_2070_rcp_45)

png(file="Desktop/2070_45_difference_plots.png",width=1273.488192,height=469.516512)
par(mfrow=c(1,3))

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_tridentata_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(tri_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("G",cex=1.4,y=1300000,x=-1700000)

# ssp. wyomingensis 2070 (4.5)

gain <- rgeos::gDifference(ens_wyomingensis_2070_rcp_45,ens_wyomingensis_current)
loss <- rgeos::gDifference(ens_wyomingensis_current,ens_wyomingensis_2070_rcp_45)

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_wyomingensis_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(wyo_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("H",cex=1.4,y=1300000,x=-1700000)

# ssp. vaseyana 2070 (4.5)

gain <- rgeos::gDifference(ens_vaseyana_2070_rcp_45,ens_vaseyana_current)
loss <- rgeos::gDifference(ens_vaseyana_current,ens_vaseyana_2070_rcp_45)

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_vaseyana_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(vas_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("I",cex=1.4,y=1300000,x=-1700000)

graphics.off()

# ssp. tridentata 2070 (8.5)

gain <- rgeos::gDifference(ens_tridentata_2070_rcp_85,ens_tridentata_current)
loss <- rgeos::gDifference(ens_tridentata_current,ens_tridentata_2070_rcp_85)

png(file="Desktop/2070_85_difference_plots.png",width=1273.488192,height=469.516512)
par(mfrow=c(1,3))

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_tridentata_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(tri_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("J",cex=1.4,y=1300000,x=-1700000)

# ssp. wyomingensis 2070 (8.5)

gain <- rgeos::gDifference(ens_wyomingensis_2070_rcp_85,ens_wyomingensis_current)
loss <- rgeos::gDifference(ens_wyomingensis_current,ens_wyomingensis_2070_rcp_85)

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_wyomingensis_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(wyo_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("K",cex=1.4,y=1300000,x=-1700000)

# ssp. vaseyana 2070 (8.5)

gain <- rgeos::gDifference(ens_vaseyana_2070_rcp_85,ens_vaseyana_current)
loss <- rgeos::gDifference(ens_vaseyana_current,ens_vaseyana_2070_rcp_85)

plot(main=NA,spTransform(boundaries,t_crs),col="white",
     border=NA,xlim=c(-1449624,103145),ylim=c(-1258277,1363081), axes=T,cex=0.65)
plot(spTransform(ens_vaseyana_current,t_crs), border=NA, col="#006600",add=T)
plot(spTransform(gain,t_crs),border=NA,col="#0000FF",add=T)
plot(spTransform(loss,t_crs),border=NA,col="#CC0000",add=T)
plot(vas_gap,pch=15,cex=0.65,add=T)
plot(spTransform(boundaries,t_crs), border=rgb(0, 0, 0, 0.5),add=T);
box(); grid(lty=1,col="#00000030");legend("topright", c("maintain","gain","loss"), cex=0.8, fill=c("#006600","#0000FF","#CC0000"),bg = "white");
text("L",cex=1.4,y=1300000,x=-1700000)

graphics.off()

# Figure 8 -- Response plots with histogram overlay

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
      grid(lwd=1.4); lines(prob~get(var),lwd=2.5,col="red",data=prob)
  }
}

respPlotDensityOverlay <- function(m_glm=NULL,m_rf=NULL,name="wyo"){
  vars <<- sort(names(m_glm[[1]][[1]]$data)[grepl(names(m_glm[[1]][[1]]$data),pattern="bio")]) # this has to be global because of some weird bug in random forest
  # response plots for wyomingensis
  png(paste(sep="",HOME,"/Desktop/",name,"_response_plots_glm.png"),height=1250,width=850)
    par(mfrow=c(5,1),cex.lab=1.8,cex.axis=1.8)
      for(i in 1:length(vars)){
        responsePlot(m_glm,var=vars[i])
        out <- partialPlot(m_rf[[1]],x.var=vars[i],pred.data=na.omit(m_glm[[1]][[1]]$data),which.class=1,plot=F)
        out$y <- exp(out$y);
        out$y <- (out$y/max(out$y))
        lines(y=out$y, x=out$x,lwd=2.5,col="blue",main="",xlab=as.character(vars[i]))
        h1 <- density(na.omit(m_glm[[1]][[1]]$data[m_glm[[1]][[1]]$data$resp==1,vars[i]]))
          h1$y <- h1$y/max(h1$y)
        lines(h1,col="#99996699",lwd=2.5,lty=15)
        h2 <- density(na.omit(m_glm[[1]][[1]]$data[m_glm[[1]][[1]]$data$resp==0,vars[i]]))
          h2$y <- h2$y/max(h2$y)
        lines(h2,col="#CCCC0099",lwd=2.5,lty=15)
      };
  graphics.off();
}

load("/media/ktaylora/big_black/products/uw/big_sagebrush_subspp_analysis/models.Rdata")

# response plots for wyomingensis
respPlotDensityOverlay(m_glm=wyomingensis_glm_unif,m_rf=wyomingensis_rf_unif,name="wyo")
respPlotDensityOverlay(m_glm=tridentata_glm_unif,m_rf=tridentata_rf_unif,name="tri")
respPlotDensityOverlay(m_glm=vaseyana_glm_unif,m_rf=vaseyana_rf_unif,name="vas")

# png(paste(sep="/",HOME,"/Desktop/wyo_response_plots_glm.png"),height=1250,width=850)
#   par(mfrow=c(5,1),cex.lab=1.8,cex.axis=1.8)
#     for(v in vars){
#       responsePlot(wyomingensis_glm_unif,var=v)
#       out <- partialPlot(wyomingensis_rf_unif[[1]],x.var=as.character(v),pred.data=na.omit(wyomingensis_glm_unif[[1]][[1]]$data),which.class=1,plot=F)
#       out$y <- exp(out$y);
#       out$y <- (out$y/max(out$y))
#       lines(y=out$y, x=out$x,lwd=2.5,col="blue",main="",xlab=as.character(v))
#       h1 <- density(na.omit(wyomingensis_glm_unif[[1]][[1]]$data[wyomingensis_glm_unif[[1]][[1]]$data$resp==1,v]))
#         h1$y <- h1$y/max(h1$y)
#       lines(h1,col="#99996699",lwd=2.5,lty=15)
#       h2 <- density(na.omit(wyomingensis_glm_unif[[1]][[1]]$data[wyomingensis_glm_unif[[1]][[1]]$data$resp==0,v]))
#         h2$y <- h2$y/max(h2$y)
#       lines(h2,col="#CCCC0099",lwd=2.5,lty=15)
#     };
# graphics.off();
