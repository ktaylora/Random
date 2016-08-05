#!/bin/bash

#
# GDAL-based GIS Workflow for Estimating Volume of Water Contributed by Playas to the Aquifer
#
# Assumes that the CWD contains (1) a time-series of LANDSAT "wetness" surfaces as .img files, 
# (2) focal counties as a shapefile, and (3) the PLJV probable playas (v4) shapefile.
#
# Will produce a series of .tif files (one for each year) that indicate volume of water (in acre-feet)
# contributed to ground-water (i.e., the aquifer) using expert opinion collected from the PLJV playa symposium.
#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2016]
#

# search the CWD for all filenames matching wet*.img -- establish the min and max numbers of our 'years' for wetness
years=(`ls -1 wet*.img | awk '{ sub("wet","",$1); print $0 }' | awk '{ sub(".img",""); print $0 }' | sort -n | sed -n '1p;$p'`)

# pre-run clean-up
rm -rf `ls -1 *.tif | grep -v "playa"`
# crop wetness to our focal counties to speed up resampling
seq ${years[0]} ${years[1]} | awk '{ print "gdalwarp -cutline counties.shp -crop_to_cutline wet"$1".img wet"$1".tif" }' | /bin/bash
# resample our wetness rasters to a target CRS consistent with our playas raster
for i in $(seq 2004 2014); do
  gdalwarp -to playas_curry_quay_counties.tif -tr 30.00149 30.00019 -t_srs \
    "+proj=utm +zone=14 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0" \
    wet$i.tif wet$i.resamp.tif
done;

# crop resampled wetness surfaces to extent of our focal counties for a final time
seq ${years[0]} ${years[1]} | awk '{ print "rm -rf wet" $1".tif" }' | /bin/bash
seq ${years[0]} ${years[1]}| awk '{ print "gdalwarp -cutline counties.shp -crop_to_cutline wet"$1".resamp.tif wet"$1".tif" }' | /bin/bash
# calculate saturated feet contributed by wet playas
# Per Cell: Wet[1,0] * Playa [1,0] * { (1181.1610236 * 1181.1098425) in^2 * (3 inches) = 4,185,242.73165 in^3 = 0.05560190720488647 acre-feet }
seq ${years[0]} ${years[1]} | awk '{ print "gdal_calc.py --type=Float64 -A wet"$1".tif -B playas_curry_quay_counties.tif --calc=A*B*0.05560190720488647 --outfile="$1"_1.tif" }' | /bin/bash
# calculate saturated feet contributed by interspace (CURRENLTY ZERO'd)
seq ${years[0]} ${years[1]} | awk '{ print "gdal_calc.py --type=Float64 -A wet"$1".tif -B not_playas_curry_quay_counties.tif --calc=A*B*0.000231660367*0 --outfile="$1"_2.tif" }' | /bin/bash
# combine playas and interspace into a joint raster surface
seq ${years[0]} ${years[1]} | awk '{ print "gdal_calc.py --type=Float64 -A "$1"_1.tif -B "$1"_2.tif --calc=A+B --outfile="$1"_3.tif" }' | /bin/bash
# clean-up
rm -rf *_1*.tif *_2*.tif
seq ${years[0]} ${years[1]}| awk '{ print "mv "$1"_3.tif "$1".tif; rm -rf "$1"_3.tif" }' | /bin/bash
