# author: Kyle Taylor (kyle.taylor@pljv.org)
# see inspiration: http://www.petrkeil.com/?p=1050

require(ncf)
require(quantreg)

boot_correlogram <- function(s=NULL,nboot=50){
  x <- vector()
  y <- vector()
  for (i in 1:nboot){
    cat(paste("(", i, "/", nboot, ")", sep = ""));
    s_down <- s[sample(1:nrow(s), 0.65*nrow(s), replace=T), ]
    o <- correlog(x = s_down$x, y = s_down$y, z = s_down$near_dist, resamp = 0, increment = 9000, quiet = T)
    x <- append(x, o$mean.of.class)
    y <- append(y, o$correlation)
  }; cat("\n");
  return(list(x,y))
}

s_control <- readOGR(".","randombcr18pts")
  s_control <- data.frame(x=s_control@coords[,1],y=s_control@coords[,2],near_dist=s_control$NEAR_DIST)

o <- boot_correlogram(s_control)
x_control <- o[[1]]
y_control <- o[[2]]

plot(y_control~x_control, cex = 0.25, pch = 21, xlab = "Lag Distance (m)", ylab = "Correlation", col="DarkRed",ylim=c(-6,6))

s_playas <- readOGR(".","randomplayaclusterpts")
  s_playas <- data.frame(x=s_playas@coords[,1],y=s_playas@coords[,2],near_dist=s_playas$NEAR_DIST)

o <- boot_correlogram(s_playas)
x_experimental <- o[[1]]
y_experimental <- o[[2]]

points(y_experimental~x_experimental, cex = 0.25, pch = 21, col="DarkBlue")
grid()
abline(h=0)
abline(v=0)
