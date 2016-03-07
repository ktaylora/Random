bash get.sh

for z in `ls *.zip`; do
  7za x $z
done

rm -rf *.zip

for z in `ls *.img`; do
  echo " -- converting" $z "to GeoTIFF"
  f=`echo $z | awk '{ print substr($1,1,index($1,".")) "tif" }'`
  gdalwarp -of GTIFF $z $f
done

rm -rf *.img *.ige *.rrd *.rde
