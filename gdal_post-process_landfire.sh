# calculate grass height (units are meters / 100)
# percent cover units are 0.N * 100
in=`echo "$1" | awk '{print tolower($0)}'`
if echo $in | grep -q "evh"; then
  rm -rf grass_height* shrub_height* tree_height*
  gdal_calc.py -A $1 --outfile=grass_height.tif --calc="(A<=103)*(A>=101)*(((2*(A-100)*0.25))-0.25)*100" #25,75,125
  gdal_calc.py -A $1 --outfile=shrub_height.tif --calc="(0*(A<104))+(0*(A>107))+(250*(A==104))+(750*(A==105))+(2000*(A==106))+(3250*(A==107))" --type="Int16"
  gdal_calc.py -A $1 --outfile=tree_height.tif --calc="0*(A<108)+(0*(A>111))+((A==108)*25)+((A==109)*750)+((A==110)*1750)+((A==111)*3750)" --type="Int16"
elif echo $in | grep -q "evh"; then
  rm -rf grass_perc_cover* shrub_perc_cover* tree_perc_cover*
  gdal_calc.py -A $1 --outfile=grass_perc_cover.tif --calc="(0*(A<121))+(0*(A>129))+(((A-120)*10)+5)" --type="Int16"
  gdal_calc.py -A $1 --outfile=shrub_perc_cover.tif --calc="(0*(A<111))+(0*(A>119))+(((A-110)*10)+5)" --type="Int16"
  gdal_calc.py -A $1 --outfile=tree_perc_cover.tif --calc="(0*(A<101))+(0*(A>109))+(((A-100)*10)+5)" --type="Int16"
fi
