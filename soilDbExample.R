#
# Do something with this later
#

require(rgdal)
require(soilDB)

# get map unit keys and polygons for this bbox
b <- c(-120.54,38.61,-120.41,38.70)
x <- mapunit_geom_by_ll_bbox(b)
m <- as.character(x@data$mukey)

# format and submit SQL "in" statement
in.statement <- format_SQL_in_statement(m)

q <- paste("SELECT component.cokey, mukey, compname, comppct_r, hzdept_r, hzdepb_r,
                   hzname, sandtotal_r, silttotal_r, claytotal_r
                   FROM component JOIN chorizon ON component.cokey = chorizon.cokey
                   WHERE majcompflag = 'Yes' AND mukey IN ", in.statement,
                   "ORDER BY mukey, comppct_r DESC, hzdept_r ASC", sep="")

# now get component and horizon-level data for these map unit keys
res <- SDA_query(q)
# merge to our shapefile
x@data <- merge(x@data,res[,c(2,4:10)],by="mukey")
# plot with some pretty colors
plot(x,col=x$claytotal_r)
