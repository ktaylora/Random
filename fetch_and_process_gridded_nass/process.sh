#!/bin/bash

#bash get.sh

for z in `ls *.zip`; do
  7za x $z
done

#rm -rf *.zip

for z in `ls *.img`; do
  echo " -- converting" $z "to GeoTIFF"
  f=`echo $z | awk '{ print substr($1,1,index($1,".")) "tif" }'`
  gdalwarp -of GTIFF $z $f
  if [ -r $1 ]; then
    gdalwarp -cutline $1 -crop_to_cutline $f `echo $f | awk '{ print substr($1,1,index($1,".")) "cut.tif" }'` # append with .cut.tif
    mv `echo $f | awk '{ print substr($1,1,index($1,".")) "cut.tif" }'` $f # over-write original
  fi
done

rm -rf *.img *.ige *.rrd *.rde
