require(raster)
require(rgdal)
require(randomForest)

# project into future climate conditions
setwd("/Volumes/big_black/intermediates/weather/sagebrush_subspp_future_conditions/focal_for_ssp_manuscript")
focal_vars <- paste(paste("*.",c("03.tif$","04.tif$","11.tif$","15.tif$","18.tif$"),sep=""),collapse="|")

## 2050 (Wyomingensis)
# gf45
utils::unzip("gf45bi70.zip", exdir="/tmp/gf45bi70/")
gfdl_cm3_2070_45_wyo <- list.files("/tmp/gf45bi70",pattern="tif$",full.names=T)
  gfdl_cm3_2070_45_wyo <- raster::stack(gfdl_cm3_2070_45_wyo[grepl(gfdl_cm3_2070_45_wyo,pattern=focal_vars)])
    gfdl_cm3_2070_45_wyo <- raster::crop(gfdl_cm3_2070_45_wyo,climate_variables,progress='text')
      names(gfdl_cm3_2070_45_wyo) <- names(climate_variables) # the order of the files are in agreement

gfdl_cm3_2070_45_wyo_glm <- predict(gfdl_cm3_2070_45_wyo,wyomingensis_glm_unif[[1]][[1]],type='resp',progress='text')
gfdl_cm3_2070_45_wyo_rf  <- 1-predict(gfdl_cm3_2070_45_wyo,wyomingensis_rf_unif[[1]],type='prob',progress='text')
gfdl_cm3_2070_45_wyo_ens <- stackApply(stack(gfdl_cm3_2070_45_wyo_glm,gfdl_cm3_2070_45_wyo_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/gf45bi70",recursive=T,force=T)
# gf85
utils::unzip("gf85bi70.zip", exdir="/tmp/gf85bi70/")
gfdl_cm3_2070_85_wyo <- list.files("/tmp/gf85bi70",pattern="tif$",full.names=T)
  gfdl_cm3_2070_85_wyo <- raster::stack(gfdl_cm3_2070_85_wyo[grepl(gfdl_cm3_2070_85_wyo,pattern=focal_vars)])
    gfdl_cm3_2070_85_wyo <- raster::crop(gfdl_cm3_2070_85_wyo,climate_variables,progress='text')
      names(gfdl_cm3_2070_85_wyo) <- names(climate_variables)

gfdl_cm3_2070_85_wyo_glm <- predict(gfdl_cm3_2070_85_wyo,wyomingensis_glm_unif[[1]][[1]],type='resp',progress='text')
gfdl_cm3_2070_85_wyo_rf  <- 1-predict(gfdl_cm3_2070_85_wyo,wyomingensis_rf_unif[[1]],type='prob',progress='text')
gfdl_cm3_2070_85_wyo_ens <- stackApply(stack(gfdl_cm3_2070_85_wyo_glm,gfdl_cm3_2070_85_wyo_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/gf85bi70",recursive=T,force=T)
# ip45
utils::unzip("ip45bi70.zip", exdir="/tmp/ip45bi70/")
ipsl_cm5b_lr_2070_45_wyo <- list.files("/tmp/ip45bi70",pattern="tif$",full.names=T)
  ipsl_cm5b_lr_2070_45_wyo <- raster::stack(ipsl_cm5b_lr_2070_45_wyo[grepl(ipsl_cm5b_lr_2070_45_wyo,pattern=focal_vars)])
    ipsl_cm5b_lr_2070_45_wyo <- raster::crop(ipsl_cm5b_lr_2070_45_wyo,climate_variables,progress='text')
      names(ipsl_cm5b_lr_2070_45_wyo) <- names(climate_variables) # the order of the files are in agreement

ipsl_cm5b_lr_2070_45_wyo_glm <- predict(ipsl_cm5b_lr_2070_45_wyo,wyomingensis_glm_unif[[1]][[1]],type='resp',progress='text')
ipsl_cm5b_lr_2070_45_wyo_rf  <- 1-predict(ipsl_cm5b_lr_2070_45_wyo,wyomingensis_rf_unif[[1]],type='prob',progress='text')
ipsl_cm5b_lr_2070_45_wyo_ens <- stackApply(stack(ipsl_cm5b_lr_2070_45_wyo_glm,ipsl_cm5b_lr_2070_45_wyo_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ip45bi70",recursive=T,force=T)
# ip85
utils::unzip("ip85bi70.zip", exdir="/tmp/ip85bi70/")
ipsl_cm5b_lr_2070_85_wyo <- list.files("/tmp/ip85bi70",pattern="tif$",full.names=T)
  ipsl_cm5b_lr_2070_85_wyo <- raster::stack(ipsl_cm5b_lr_2070_85_wyo[grepl(ipsl_cm5b_lr_2070_85_wyo,pattern=focal_vars)])
    ipsl_cm5b_lr_2070_85_wyo <- raster::crop(ipsl_cm5b_lr_2070_85_wyo,climate_variables,progress='text')
      names(ipsl_cm5b_lr_2070_85_wyo) <- names(climate_variables) # the order of the files are in agreement

ipsl_cm5b_lr_2070_85_wyo_glm <- predict(ipsl_cm5b_lr_2070_85_wyo,wyomingensis_glm_unif[[1]][[1]],type='resp',progress='text')
ipsl_cm5b_lr_2070_85_wyo_rf  <- 1-predict(ipsl_cm5b_lr_2070_85_wyo,wyomingensis_rf_unif[[1]],type='prob',progress='text')
ipsl_cm5b_lr_2070_85_wyo_ens <- stackApply(stack(ipsl_cm5b_lr_2070_85_wyo_glm,ipsl_cm5b_lr_2070_85_wyo_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ip85bi70",recursive=T,force=T)
# ac45
utils::unzip("ac45bi70.zip", exdir="/tmp/ac45bi70/")
access_1_0_2070_45_wyo <- list.files("/tmp/ac45bi70",pattern="tif$",full.names=T)
  access_1_0_2070_45_wyo <- raster::stack(access_1_0_2070_45_wyo[grepl(access_1_0_2070_45_wyo,pattern=focal_vars)])
    access_1_0_2070_45_wyo <- raster::crop(access_1_0_2070_45_wyo,climate_variables,progress='text')
      names(access_1_0_2070_45_wyo) <- names(climate_variables) # the order of the files are in agreement

access_1_0_2070_45_wyo_glm <- predict(access_1_0_2070_45_wyo,wyomingensis_glm_unif[[1]][[1]],type='resp',progress='text')
access_1_0_2070_45_wyo_rf  <- 1-predict(access_1_0_2070_45_wyo,wyomingensis_rf_unif[[1]],type='prob',progress='text')
access_1_0_2070_45_wyo_ens <- stackApply(stack(access_1_0_2070_45_wyo_glm,access_1_0_2070_45_wyo_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ac45bi70",recursive=T,force=T)
# ac85
utils::unzip("ac85bi70.zip", exdir="/tmp/ac85bi70/")
access_1_0_2070_85_wyo <- list.files("/tmp/ac85bi70",pattern="tif$",full.names=T)
  access_1_0_2070_85_wyo <- raster::stack(access_1_0_2070_85_wyo[grepl(access_1_0_2070_85_wyo,pattern=focal_vars)])
    access_1_0_2070_85_wyo <- raster::crop(access_1_0_2070_85_wyo,climate_variables,progress='text')
      names(access_1_0_2070_85_wyo) <- names(climate_variables) # the order of the files are in agreement

access_1_0_2070_85_wyo_glm <- predict(access_1_0_2070_85_wyo,wyomingensis_glm_unif[[1]][[1]],type='resp',progress='text')
access_1_0_2070_85_wyo_rf  <- 1-predict(access_1_0_2070_85_wyo,wyomingensis_rf_unif[[1]],type='prob',progress='text')
access_1_0_2070_85_wyo_ens <- stackApply(stack(access_1_0_2070_85_wyo_glm,access_1_0_2070_85_wyo_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ac85bi70",recursive=T,force=T)
# in45
utils::unzip("in45bi70.zip", exdir="/tmp/in45bi70/")
inmcm4_2070_45 <- list.files("/tmp/in45bi70",pattern="tif$",full.names=T)
  inmcm4_2070_45 <- raster::stack(inmcm4_2070_45[grepl(inmcm4_2070_45,pattern=focal_vars)])
    inmcm4_2070_45 <- raster::crop(inmcm4_2070_45,climate_variables,progress='text')
      names(inmcm4_2070_45) <- names(climate_variables) # the order of the files are in agreement

inmcm4_2070_45_glm <- predict(inmcm4_2070_45,wyomingensis_glm_unif[[1]][[1]],type='resp',progress='text')
inmcm4_2070_45_rf  <- 1-predict(inmcm4_2070_45,wyomingensis_rf_unif[[1]],type='prob',progress='text')
inmcm4_2070_45_ens <- stackApply(stack(inmcm4_2070_45_glm,inmcm4_2070_45_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/in45bi70",recursive=T,force=T)
# in85
utils::unzip("in85bi70.zip", exdir="/tmp/in85bi70/")
inmcm4_2070_85_wyo <- list.files("/tmp/in85bi70",pattern="tif$",full.names=T)
  inmcm4_2070_85_wyo <- raster::stack(inmcm4_2070_85_wyo[grepl(inmcm4_2070_85_wyo,pattern=focal_vars)])
    inmcm4_2070_85_wyo <- raster::crop(inmcm4_2070_85_wyo,climate_variables,progress='text')
      names(inmcm4_2070_85_wyo) <- names(climate_variables) # the order of the files are in agreement

inmcm4_2070_85_wyo_glm <- predict(inmcm4_2070_85_wyo,wyomingensis_glm_unif[[1]][[1]],type='resp',progress='text')
inmcm4_2070_85_wyo_rf  <- 1-predict(inmcm4_2070_85_wyo,wyomingensis_rf_unif[[1]],type='prob',progress='text')
inmcm4_2070_85_wyo_ens <- stackApply(stack(inmcm4_2070_85_wyo_glm,inmcm4_2070_85_wyo_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/in85bi70",recursive=T,force=T)

## 2050 (Tridentata)
# gf45
utils::unzip("gf45bi70.zip", exdir="/tmp/gf45bi70/")
gfdl_cm3_2070_45_tri <- list.files("/tmp/gf45bi70",pattern="tif$",full.names=T)
  gfdl_cm3_2070_45_tri <- raster::stack(gfdl_cm3_2070_45_tri[grepl(gfdl_cm3_2070_45_tri,pattern=focal_vars)])
    gfdl_cm3_2070_45_tri <- raster::crop(gfdl_cm3_2070_45_tri,climate_variables,progress='text')
      names(gfdl_cm3_2070_45_tri) <- names(climate_variables) # the order of the files are in agreement

gfdl_cm3_2070_45_tri_glm <- predict(gfdl_cm3_2070_45_tri,tridentata_glm_unif[[1]][[1]],type='resp',progress='text')
gfdl_cm3_2070_45_tri_rf  <- 1-predict(gfdl_cm3_2070_45_tri,tridentata_rf_unif[[1]],type='prob',progress='text')
gfdl_cm3_2070_45_tri_ens <- stackApply(stack(gfdl_cm3_2070_45_tri_glm,gfdl_cm3_2070_45_tri_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/gf45bi70",recursive=T,force=T)
# gf85
utils::unzip("gf85bi70.zip", exdir="/tmp/gf85bi70/")
gfdl_cm3_2070_85_tri <- list.files("/tmp/gf85bi70",pattern="tif$",full.names=T)
  gfdl_cm3_2070_85_tri <- raster::stack(gfdl_cm3_2070_85_tri[grepl(gfdl_cm3_2070_85_tri,pattern=focal_vars)])
    gfdl_cm3_2070_85_tri <- raster::crop(gfdl_cm3_2070_85_tri,climate_variables,progress='text')
      names(gfdl_cm3_2070_85_tri) <- names(climate_variables)

gfdl_cm3_2070_85_tri_glm <- predict(gfdl_cm3_2070_85_tri,tridentata_glm_unif[[1]][[1]],type='resp',progress='text')
gfdl_cm3_2070_85_tri_rf  <- 1-predict(gfdl_cm3_2070_85_tri,tridentata_rf_unif[[1]],type='prob',progress='text')
gfdl_cm3_2070_85_tri_ens <- stackApply(stack(gfdl_cm3_2070_85_tri_glm,gfdl_cm3_2070_85_tri_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/gf85bi70",recursive=T,force=T)
# ip45
utils::unzip("ip45bi70.zip", exdir="/tmp/ip45bi70/")
ipsl_cm5b_lr_2070_45_tri <- list.files("/tmp/ip45bi70",pattern="tif$",full.names=T)
  ipsl_cm5b_lr_2070_45_tri <- raster::stack(ipsl_cm5b_lr_2070_45_tri[grepl(ipsl_cm5b_lr_2070_45_tri,pattern=focal_vars)])
    ipsl_cm5b_lr_2070_45_tri <- raster::crop(ipsl_cm5b_lr_2070_45_tri,climate_variables,progress='text')
      names(ipsl_cm5b_lr_2070_45_tri) <- names(climate_variables) # the order of the files are in agreement

ipsl_cm5b_lr_2070_45_tri_glm <- predict(ipsl_cm5b_lr_2070_45_tri,tridentata_glm_unif[[1]][[1]],type='resp',progress='text')
ipsl_cm5b_lr_2070_45_tri_rf  <- 1-predict(ipsl_cm5b_lr_2070_45_tri,tridentata_rf_unif[[1]],type='prob',progress='text')
ipsl_cm5b_lr_2070_45_tri_ens <- stackApply(stack(ipsl_cm5b_lr_2070_45_tri_glm,ipsl_cm5b_lr_2070_45_tri_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ip45bi70",recursive=T,force=T)
# ip85
utils::unzip("ip85bi70.zip", exdir="/tmp/ip85bi70/")
ipsl_cm5b_lr_2070_85_tri <- list.files("/tmp/ip85bi70",pattern="tif$",full.names=T)
  ipsl_cm5b_lr_2070_85_tri <- raster::stack(ipsl_cm5b_lr_2070_85_tri[grepl(ipsl_cm5b_lr_2070_85_tri,pattern=focal_vars)])
    ipsl_cm5b_lr_2070_85_tri <- raster::crop(ipsl_cm5b_lr_2070_85_tri,climate_variables,progress='text')
      names(ipsl_cm5b_lr_2070_85_tri) <- names(climate_variables) # the order of the files are in agreement

ipsl_cm5b_lr_2070_85_tri_glm <- predict(ipsl_cm5b_lr_2070_85_tri,tridentata_glm_unif[[1]][[1]],type='resp',progress='text')
ipsl_cm5b_lr_2070_85_tri_rf  <- 1-predict(ipsl_cm5b_lr_2070_85_tri,tridentata_rf_unif[[1]],type='prob',progress='text')
ipsl_cm5b_lr_2070_85_tri_ens <- stackApply(stack(ipsl_cm5b_lr_2070_85_tri_glm,ipsl_cm5b_lr_2070_85_tri_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/gf85bi70",recursive=T,force=T)
# ac45
utils::unzip("ac45bi70.zip", exdir="/tmp/ac45bi70/")
access_1_0_2070_45_tri <- list.files("/tmp/ac45bi70",pattern="tif$",full.names=T)
  access_1_0_2070_45_tri <- raster::stack(access_1_0_2070_45_tri[grepl(access_1_0_2070_45_tri,pattern=focal_vars)])
    access_1_0_2070_45_tri <- raster::crop(access_1_0_2070_45_tri,climate_variables,progress='text')
      names(access_1_0_2070_45_tri) <- names(climate_variables) # the order of the files are in agreement

access_1_0_2070_45_tri_glm <- predict(access_1_0_2070_45_tri,tridentata_glm_unif[[1]][[1]],type='resp',progress='text')
access_1_0_2070_45_tri_rf  <- 1-predict(access_1_0_2070_45_tri,tridentata_rf_unif[[1]],type='prob',progress='text')
access_1_0_2070_45_tri_ens <- stackApply(stack(access_1_0_2070_45_tri_glm,access_1_0_2070_45_tri_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ac45bi70",recursive=T,force=T)
# ac85
utils::unzip("ac85bi70.zip", exdir="/tmp/ac85bi70/")
access_1_0_2070_85_tri <- list.files("/tmp/ac85bi70",pattern="tif$",full.names=T)
  access_1_0_2070_85_tri <- raster::stack(access_1_0_2070_85_tri[grepl(access_1_0_2070_85_tri,pattern=focal_vars)])
    access_1_0_2070_85_tri <- raster::crop(access_1_0_2070_85_tri,climate_variables,progress='text')
      names(access_1_0_2070_85_tri) <- names(climate_variables) # the order of the files are in agreement

access_1_0_2070_85_tri_glm <- predict(access_1_0_2070_85_tri,tridentata_glm_unif[[1]][[1]],type='resp',progress='text')
access_1_0_2070_85_tri_rf  <- 1-predict(access_1_0_2070_85_tri,tridentata_rf_unif[[1]],type='prob',progress='text')
access_1_0_2070_85_tri_ens <- stackApply(stack(access_1_0_2070_85_tri_glm,access_1_0_2070_85_tri_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ac85bi70",recursive=T,force=T)
# in45
utils::unzip("in45bi70.zip", exdir="/tmp/in45bi70/")
inmcm4_2070_45 <- list.files("/tmp/in45bi70",pattern="tif$",full.names=T)
  inmcm4_2070_45 <- raster::stack(inmcm4_2070_45[grepl(inmcm4_2070_45,pattern=focal_vars)])
    inmcm4_2070_45 <- raster::crop(inmcm4_2070_45,climate_variables,progress='text')
      names(inmcm4_2070_45) <- names(climate_variables) # the order of the files are in agreement

inmcm4_2070_45_glm <- predict(inmcm4_2070_45,tridentata_glm_unif[[1]][[1]],type='resp',progress='text')
inmcm4_2070_45_rf  <- 1-predict(inmcm4_2070_45,tridentata_rf_unif[[1]],type='prob',progress='text')
inmcm4_2070_45_ens <- stackApply(stack(inmcm4_2070_45_glm,inmcm4_2070_45_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ac45bi70",recursive=T,force=T)
# in85
utils::unzip("in85bi70.zip", exdir="/tmp/in85bi70/")
inmcm4_2070_85_tri <- list.files("/tmp/in85bi70",pattern="tif$",full.names=T)
  inmcm4_2070_85_tri <- raster::stack(inmcm4_2070_85_tri[grepl(inmcm4_2070_85_tri,pattern=focal_vars)])
    inmcm4_2070_85_tri <- raster::crop(inmcm4_2070_85_tri,climate_variables,progress='text')
      names(inmcm4_2070_85_tri) <- names(climate_variables) # the order of the files are in agreement

inmcm4_2070_85_tri_glm <- predict(inmcm4_2070_85_tri,tridentata_glm_unif[[1]][[1]],type='resp',progress='text')
inmcm4_2070_85_tri_rf  <- 1-predict(inmcm4_2070_85_tri,tridentata_rf_unif[[1]],type='prob',progress='text')
inmcm4_2070_85_tri_ens <- stackApply(stack(inmcm4_2070_85_tri_glm,inmcm4_2070_85_tri_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/in85bi70",recursive=T,force=T)

## 2050 (Vaseyana)
# gf45
utils::unzip("gf45bi70.zip", exdir="/tmp/gf45bi70/")
gfdl_cm3_2070_45_vas <- list.files("/tmp/gf45bi70",pattern="tif$",full.names=T)
  gfdl_cm3_2070_45_vas <- raster::stack(gfdl_cm3_2070_45_vas[grepl(gfdl_cm3_2070_45_vas,pattern=focal_vars)])
    gfdl_cm3_2070_45_vas <- raster::crop(gfdl_cm3_2070_45_vas,climate_variables,progress='text')
      names(gfdl_cm3_2070_45_vas) <- names(climate_variables) # the order of the files are in agreement

gfdl_cm3_2070_45_vas_glm <- predict(gfdl_cm3_2070_45_vas,vaseyana_glm_unif[[1]][[1]],type='resp',progress='text')
gfdl_cm3_2070_45_vas_rf  <- 1-predict(gfdl_cm3_2070_45_vas,vaseyana_rf_unif[[1]],type='prob',progress='text')
gfdl_cm3_2070_45_vas_ens <- stackApply(stack(gfdl_cm3_2070_45_vas_glm,gfdl_cm3_2070_45_vas_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/gf45bi70",recursive=T,force=T)
# gf85
utils::unzip("gf85bi70.zip", exdir="/tmp/gf85bi70/")
gfdl_cm3_2070_85_vas <- list.files("/tmp/gf85bi70",pattern="tif$",full.names=T)
  gfdl_cm3_2070_85_vas <- raster::stack(gfdl_cm3_2070_85_vas[grepl(gfdl_cm3_2070_85_vas,pattern=focal_vars)])
    gfdl_cm3_2070_85_vas <- raster::crop(gfdl_cm3_2070_85_vas,climate_variables,progress='text')
      names(gfdl_cm3_2070_85_vas) <- names(climate_variables)

gfdl_cm3_2070_85_vas_glm <- predict(gfdl_cm3_2070_85_vas,vaseyana_glm_unif[[1]][[1]],type='resp',progress='text')
gfdl_cm3_2070_85_vas_rf  <- 1-predict(gfdl_cm3_2070_85_vas,vaseyana_rf_unif[[1]],type='prob',progress='text')
gfdl_cm3_2070_85_vas_ens <- stackApply(stack(gfdl_cm3_2070_85_vas_glm,gfdl_cm3_2070_85_vas_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/gf85bi70",recursive=T,force=T)
# ip45
utils::unzip("ip45bi70.zip", exdir="/tmp/ip45bi70/")
ipsl_cm5b_lr_2070_45_vas <- list.files("/tmp/gf45bi70",pattern="tif$",full.names=T)
  ipsl_cm5b_lr_2070_45_vas <- raster::stack(ipsl_cm5b_lr_2070_45_vas[grepl(ipsl_cm5b_lr_2070_45_vas,pattern=focal_vars)])
    ipsl_cm5b_lr_2070_45_vas <- raster::crop(ipsl_cm5b_lr_2070_45_vas,climate_variables,progress='text')
      names(ipsl_cm5b_lr_2070_45_vas) <- names(climate_variables) # the order of the files are in agreement

ipsl_cm5b_lr_2070_45_vas_glm <- predict(ipsl_cm5b_lr_2070_45_vas,vaseyana_glm_unif[[1]][[1]],type='resp',progress='text')
ipsl_cm5b_lr_2070_45_vas_rf  <- 1-predict(ipsl_cm5b_lr_2070_45_vas,vaseyana_rf_unif[[1]],type='prob',progress='text')
ipsl_cm5b_lr_2070_45_vas_ens <- stackApply(stack(ipsl_cm5b_lr_2070_45_vas_glm,ipsl_cm5b_lr_2070_45_vas_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ip45bi70",recursive=T,force=T)
# ip85
utils::unzip("ip85bi70.zip", exdir="/tmp/ip85bi70/")
ipsl_cm5b_lr_2070_85_vas <- list.files("/tmp/ip85bi70",pattern="tif$",full.names=T)
  ipsl_cm5b_lr_2070_85_vas <- raster::stack(ipsl_cm5b_lr_2070_85_vas[grepl(ipsl_cm5b_lr_2070_85_vas,pattern=focal_vars)])
    ipsl_cm5b_lr_2070_85_vas <- raster::crop(ipsl_cm5b_lr_2070_85_vas,climate_variables,progress='text')
      names(ipsl_cm5b_lr_2070_85_vas) <- names(climate_variables) # the order of the files are in agreement

ipsl_cm5b_lr_2070_85_vas_glm <- predict(ipsl_cm5b_lr_2070_85_vas,vaseyana_glm_unif[[1]][[1]],type='resp',progress='text')
ipsl_cm5b_lr_2070_85_vas_rf  <- 1-predict(ipsl_cm5b_lr_2070_85_vas,vaseyana_rf_unif[[1]],type='prob',progress='text')
ipsl_cm5b_lr_2070_85_vas_ens <- stackApply(stack(ipsl_cm5b_lr_2070_85_vas_glm,ipsl_cm5b_lr_2070_85_vas_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/gf85bi70",recursive=T,force=T)
# ac45
utils::unzip("ac45bi70.zip", exdir="/tmp/ac45bi70/")
access_1_0_2070_45_vas <- list.files("/tmp/ac45bi70",pattern="tif$",full.names=T)
  access_1_0_2070_45_vas <- raster::stack(access_1_0_2070_45_vas[grepl(access_1_0_2070_45_vas,pattern=focal_vars)])
    access_1_0_2070_45_vas <- raster::crop(access_1_0_2070_45_vas,climate_variables,progress='text')
      names(access_1_0_2070_45_vas) <- names(climate_variables) # the order of the files are in agreement

access_1_0_2070_45_vas_glm <- predict(access_1_0_2070_45_vas,vaseyana_glm_unif[[1]][[1]],type='resp',progress='text')
access_1_0_2070_45_vas_rf  <- 1-predict(access_1_0_2070_45_vas,vaseyana_rf_unif[[1]],type='prob',progress='text')
access_1_0_2070_45_vas_ens <- stackApply(stack(access_1_0_2070_45_vas_glm,access_1_0_2070_45_vas_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ac45bi70",recursive=T,force=T)
# ac85
utils::unzip("ac85bi70.zip", exdir="/tmp/ac85bi70/")
access_1_0_2070_85_vas <- list.files("/tmp/ac85bi70",pattern="tif$",full.names=T)
  access_1_0_2070_85_vas <- raster::stack(access_1_0_2070_85_vas[grepl(access_1_0_2070_85_vas,pattern=focal_vars)])
    access_1_0_2070_85_vas <- raster::crop(access_1_0_2070_85_vas,climate_variables,progress='text')
      names(access_1_0_2070_85_vas) <- names(climate_variables) # the order of the files are in agreement

access_1_0_2070_85_vas_glm <- predict(access_1_0_2070_85_vas,vaseyana_glm_unif[[1]][[1]],type='resp',progress='text')
access_1_0_2070_85_vas_rf  <- 1-predict(access_1_0_2070_85_vas,vaseyana_rf_unif[[1]],type='prob',progress='text')
access_1_0_2070_85_vas_ens <- stackApply(stack(access_1_0_2070_85_vas_glm,access_1_0_2070_85_vas_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ac85bi70",recursive=T,force=T)
# in45
utils::unzip("in45bi70.zip", exdir="/tmp/in45bi70/")
inmcm4_2070_45 <- list.files("/tmp/in45bi70",pattern="tif$",full.names=T)
  inmcm4_2070_45 <- raster::stack(inmcm4_2070_45[grepl(inmcm4_2070_45,pattern=focal_vars)])
    inmcm4_2070_45 <- raster::crop(inmcm4_2070_45,climate_variables,progress='text')
      names(inmcm4_2070_45) <- names(climate_variables) # the order of the files are in agreement

inmcm4_2070_45_glm <- predict(inmcm4_2070_45,vaseyana_glm_unif[[1]][[1]],type='resp',progress='text')
inmcm4_2070_45_rf  <- 1-predict(inmcm4_2070_45,vaseyana_rf_unif[[1]],type='prob',progress='text')
inmcm4_2070_45_ens <- stackApply(stack(inmcm4_2070_45_glm,inmcm4_2070_45_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/ac45bi70",recursive=T,force=T)
# in85
utils::unzip("in85bi70.zip", exdir="/tmp/in85bi70/")
inmcm4_2070_85_vas <- list.files("/tmp/in85bi70",pattern="tif$",full.names=T)
  inmcm4_2070_85_vas <- raster::stack(inmcm4_2070_85_vas[grepl(inmcm4_2070_85_vas,pattern=focal_vars)])
    inmcm4_2070_85_vas <- raster::crop(inmcm4_2070_85_vas,climate_variables,progress='text')
      names(inmcm4_2070_85_vas) <- names(climate_variables) # the order of the files are in agreement

inmcm4_2070_85_vas_glm <- predict(inmcm4_2070_85_vas,vaseyana_glm_unif[[1]][[1]],type='resp',progress='text')
inmcm4_2070_85_vas_rf  <- 1-predict(inmcm4_2070_85_vas,vaseyana_rf_unif[[1]],type='prob',progress='text')
inmcm4_2070_85_vas_ens <- stackApply(stack(inmcm4_2070_85_vas_glm,inmcm4_2070_85_vas_rf),fun=mean,indices=1,progress='text')
unlink("/tmp/in85bi70",recursive=T,force=T)
