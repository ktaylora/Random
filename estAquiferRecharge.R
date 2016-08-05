require(raster)

r <- stack(paste(c(2004:2011,2013:2014),".tif",sep=""))
out <- stackApply(r,fun=sum,indices=1)

year_sums <- vector('list',nlayers(out))

for(i in 1:nlayers(r)){ 
  year_sums[[i]] <- sum(values(subset(r,i)),na.rm=T) 
}

names(year_sums) <- names(r)

cat(" -- mean:",mean(as.numeric(year_sums)),"\n")
cat(" -- sd:",sd(as.numeric(year_sums)),"\n")

hist(as.numeric(year_sums),plot=T,prob=T,xlim=c(-200,800),col="grey",main="",xlab="Saturated Ft of Recharge (per year)")
lines(density(as.numeric(year_sums)), col="blue", lwd=2)
lines(density(as.numeric(year_sums), adjust=2), lty="dotted", col="darkgreen", lwd=2)
abline(v=mean(as.numeric(year_sums),lwd=14,col="red"))
