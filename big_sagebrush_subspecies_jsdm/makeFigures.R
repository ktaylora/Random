require(rgdal)
require(rgeos)
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
elevation <- raster("/media/ktaylora/80C0-4F29/lf_elev1km/w001001x.adf")

p_current_elev_agreement <- rasterToPoints(rasterize(spTransform(spsample(intersect_p50_current,n=15000,type="random"),CRS(projection(elevation))),elevation),spatial=T)
  p_current_elev_agreement <- extract(p_current_elev_agreement,elevation,df=T)[,2]
p_current_elev_tri <- rasterToPoints(rasterize(spTransform(spsample(ens_tridentata_current,n=15000,type="random"),CRS(projection(elevation))),elevation),spatial=T)
  p_current_elev_tri <- extract(elevation,p_current_elev_tri,df=T)[,2]
p_current_elev_wyo <- rasterToPoints(rasterize(spTransform(spsample(ens_wyomingensis_current,n=15000,type="random"),CRS(projection(elevation))),elevation),spatial=T)
  p_current_elev_wyo <- extract(elevation,p_current_elev_wyo,df=T)[,2]
p_current_elev_vas <- rasterToPoints(rasterize(spTransform(spsample(ens_vaseyana_current,n=15000,type="random"),CRS(projection(elevation))),elevation),spatial=T)
  p_current_elev_vas <- extract(elevation,p_current_elev_vas,df=T)[,2]
