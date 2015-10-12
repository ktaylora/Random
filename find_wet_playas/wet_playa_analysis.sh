#!/bin/bash

#
# Quick and efficient Unix BASH interface to do a simple landscape wetness assessment from LANDSAT imagery using GDAL
#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

# Hack to get this working with the Mac OS X ports collection
if [[ `uname -a | awk '{ print $1 }'` == "Darwin" ]]; then
  PATH=$PATH:/opt/local/share/doc/py27-gdal/examples/scripts/
  PYTHON="python2.7"
fi
# SANITY CHECK FOR INPUT ARGUMENTS
if [ $# -eq 0 ]; then
    echo "usage: -p(rocess), -c(onvert to polygon shapefile) -(m)erge slices in CWD"
    exit;
# PROCESS SOURCE IMAGERY
elif [[ $1 == "-p" ]]; then 
    echo " -- processing source landsat gzip files"
    for zip in `ls -1 *.gz`; do
        gunzip -q -c $zip | tar xvf -;
        if [[ `ls -1 *.TIF | grep -v "_OUTPUT" | head -n1 | awk '{ print substr($1,1,3) }'` == "LT8" ]]; then # for Landsat 8, the wetness calculation is 6<4
          echo " -- Landsat 8 imagery detected.";
          tifs=`ls -1 *.TIF | grep -E "B4|B6"`;
          echo -n " -- processing: ";
          echo $tifs | awk '{ print "gdal_calc.py -A "$1" -B " $2 " --outfile=" substr($1,1,length($1)-6)"OUTPUT.TIF --calc=\"(B<A)\"" }' | /bin/bash;
        elif [[ `ls -1 *.TIF | grep -v "_OUTPUT" | head -n1 | awk '{ print substr($1,1,3) }'` == "LE7" ]]; then # for Landsat 7, the wetness calculation is 5<3
          echo " -- Landsat 7 imagery detected.";
          tifs=`ls -1 *.TIF | grep -E "B3|B5"`;
          echo -n " -- processing: ";
          echo $tifs | awk '{ print "gdal_calc.py -A "$1" -B " $2 " --outfile=" substr($1,1,length($1)-6)"OUTPUT.TIF --calc=\"(B<A)\"" }' | /bin/bash;
        elif [[ `ls -1 *.TIF | grep -v "_OUTPUT" | head -n1 | awk '{ print substr($1,1,3) }'` == "LT5" ]]; then # for Landsat 5, the wetness calculation is 5<3
          echo " -- Landsat 5 imagery detected.";
          tifs=`ls -1 *.TIF | grep -E "B3|B5"`;
          echo -n " -- processing: ";
          echo $tifs | awk '{ print "gdal_calc.py -A "$1" -B " $2 " --outfile=" substr($1,1,length($1)-6)"OUTPUT.TIF --calc=\"(B<A)\"" }' | /bin/bash
        else 
          echo " -- Unknown landsat imagery prefix. Quitting";
          exit;
        fi
        # remove intermediate raster imagery
        ls *.TIF | grep -v "OUTPUT" | awk '{ print "rm -rf " $1 }' | /bin/bash;
        rm -rf *.txt *.jpg *.GTF;
    done;
    # mask landsat scene periphery to minimize noise
    echo -n " -- masking landsat scene periphery:"
    rasters=`ls -1 *_OUTPUT.TIF`
      paths=`ls -1 *_OUTPUT.TIF | awk '{ print substr($1,4,3) }'`
       rows=`ls -1 *_OUTPUT.TIF | awk '{ print substr($1,7,3) }'`
    for i in $(seq 1 `echo $rasters | awk '{ print NF}'`); do
         r=`echo $rasters | awk -v col=$i '{ print $col }'`
      path=`echo $paths | awk -v col=$i '{ print $col }'`
       row=`echo $rows | awk -v col=$i '{ print $col }'`
       # select by path and row from WRS shapefile
       rm -rf /tmp/221224_ls_grid.*
       ogr2ogr /tmp/221224_ls_grid.shp ~/PLJV/boundaries/landsat_grids/wrs1_asc_desc/wrs1_asc_desc.shp -sql "SELECT * FROM wrs1_asc_desc WHERE PATH = $path AND ROW = $row"
       gdalwarp -q -cutline /tmp/221224_ls_grid.shp $r $r"_cut.tif"
       mv $r"_cut.tif" $r
       echo -n "."
    done; echo " ";
# CONVERT TO POLYGONS
elif [[ $1 == "-c" ]]; then
    echo " -- converting slices to polygons and merging"
    for tif in `ls -1 *_OUTPUT.TIF`; do
      dst=`echo $tif | awk '{ print substr($1,1,index($1,".")-1)".shp" }'`
      if [ ! -r $dst ]; then
        gdal_polygonize.py $tif -f "ESRI Shapefile" "tmp_shp.shp"
        ogr2ogr $dst tmp_shp.shp -sql "SELECT * FROM tmp_shp WHERE DN = 1"
        #ogr2ogr $dst "tmp_shp.shp" -sql "SELECT * FROM tmp_shp WHERE ((OGR_GEOM_AREA)/1000) <= 1000" # surface area is less than ~10km2
        rm -rf tmp_shp.*
      else 
        echo " -- skipping : " $dst
      fi
    done
# MOSAIC RASTER SLICES
elif [[ $1 == "-m" ]]; then
    echo " -- merging raster slices to merged.tif"
    gdal_merge.py -o merged.tif -of GTiff -tap *.TIF
fi
