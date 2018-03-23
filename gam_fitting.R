var_response_plot <- function(m=NULL, var=NULL, rug=T){
  dev.new();
  plot(
    m, 
    which=var, 
    type="l", 
    col="red", 
    lwd=1.5,
    rug=T
  );
  grid();
  grid();
}

m_gam <- gamboost(
  as.factor(response)~
    bmono(raster_115)+
    bbs(windspeed_100m_snapTo)+
    btree(raster_230)+
    btree(raster_345)+
    bmono(raster_138)+
    btree(slope_wgs)+
    btree(rgh_11_wgs)+
    btree(raster_69), 
  control = boost_control(center=T), 
  family=Binomial(),
  data=training_data
)

test <- predict(m_gam, newdata=training_data, type="response")
sum(as.numeric( test > 0.5 ) == training_data$response) / nrow(training_data)
