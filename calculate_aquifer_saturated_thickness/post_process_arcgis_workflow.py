"""
Quick-and-dirty conversion of point shapefiles to topogrid using arcpy
Author: Kyle Taylor (kyle.taylor@pljv.org)
"""

import sys, os
import arcpy
from arcpy import env
from arcpy.sa import *

env.workspace = "."
arcpy.env.overwriteOutput = True

base_contours = "base_contour_lines.shp"  # base contour lines for TOPOGRID

if arcpy.CheckExtension("3D") == "Available":
    arcpy.CheckOutExtension("3D")
else:
    exit(" -- 3D Analyst extention is unavailable")


class ShapeFiles:
    def __init__(self, *args):
        """

        :rtype: vector
        """
        self.shapes = filter(lambda x: x.endswith('.shp'), args[0])
        if len(self.shapes) == 0:
            exit(" -- quitting: no shapefiles found in path")


#
# MAIN
#

if len(sys.argv) > 2:
    d = "."
    s = ShapeFiles(os.listdir("."))
else:
    d = sys.argv[1]
    env.workspace = d
    s = ShapeFiles(os.listdir(sys.argv[1]))

print "-- processing workflow"

for shape in s.shapes:
    # inPointElevations=TopoPointElevation([s, 'strtd_t'])
    inPointElevations = '%s strtd_t POINTELEVATION; %s ELEV CONTOUR' % (shape, base_contours)
    print inPointElevations
    r = d + "\\" + shape[0:4] + ".tif"
    o = arcpy.TopoToRaster_3d(inPointElevations, r)
    # crop and mask to the extent of the High Plains region shapefile.
