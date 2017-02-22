#!/usr/bin/python

import sys
import os
import gdal
import numpy

from scipy import ndimage

class Raster(object):
    '''
    raster
    '''
    def __init__(self, **kwargs):
        self._band = 1
        self._xmin = None
        self._x_size = None
        self._y_size = None
        self._ymax = None
        for i, arg in enumerate(kwargs):
            if arg == "band":
                self._band = kwargs[arg]

    def np_open(self, file_name=None, ndv=0):
        src_ds = gdal.Open(file_name, gdal.GA_ReadOnly)
        b = src_ds.GetRasterBand(self._band)
        b_ndv = b.GetNoDataValue()
        if b_ndv is not None:
            ndv = b_ndv
        return numpy.ma.masked_equal(b.ReadAsArray(), ndv)


def focal(img=None, *args):
    return ndimage.uniform_filter(img, (args[0], args[1]))


def array_to_raster(array):
    """Array > Raster
    Save a raster from a C order array.

    :param array: ndarray
    """
    dst_filename = '/a_file/name.tiff'

    # You need to get those values like you did.
    x_pixels = 16  # number of pixels in x
    y_pixels = 16  # number of pixels in y
    PIXEL_SIZE = 3  # size of the pixel...
    x_min = 553648
    y_max = 7784555  # x_min & y_max are like the "top left" corner.
    wkt_projection = 'a projection in wkt that you got from other file'

    driver = gdal.GetDriverByName('GTiff')

    dataset = driver.Create(
        dst_filename,
        x_pixels,
        y_pixels,
        1,
        gdal.GDT_Float32, )

    dataset.SetGeoTransform((
        x_min,  # 0
        PIXEL_SIZE,  # 1
        0,  # 2
        y_max,  # 3
        0,  # 4
        -PIXEL_SIZE))

    dataset.SetProjection(wkt_projection)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()  # Write to disk.
    return dataset, dataset.GetRasterBand(1)


if __name__ == "__main__":
    INPUT_RASTER = None
    WINDOW_DIMS = None
    TARGET_RECLASS_VALUE = None
    for i in range(0, len(sys.argv)):
        if sys.argv[i] == "-r":
            INPUT_RASTER = sys.argv[i + 1]
        elif sys.argv[i] == "-t":
             TARGET_RECLASS_VALUE = list(map(int,sys.argv[i + 1].split(',')))
        elif sys.argv[i] == "-mw":
            print(" -- calculating 1-km average raster")
            r_1km = focal(img=r, 1000, 1000)
            print(" -- calculating 5-km average raster")
            r_5km = focal(img=r, 5000, 5000)