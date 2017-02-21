#!/usr/bin/python

import sys
import os
import gdal
import numpy

from scipy import ndimage

def fn_getma(src_fn, bnum=1, ndv=0):
    src_ds = gdal.Open(src_fn, gdal.GA_ReadOnly)
    b = src_ds.GetRasterBand(bnum)
    b_ndv = b.GetNoDataValue()
    if (b_ndv is not None):
        ndv = b_ndv
    bma = numpy.ma.masked_equal(b.ReadAsArray(), ndv)
    return bma


def focal(img=None, nrow=1, ncol=1, *args):
    return ndimage.uniform_filter(img, (nrow, ncol))


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
    for i in range(0, len(sys.argv)):
        if (sys.argv[i] == "-r"):
            r = fn_getma(sys.argv[i + 1])
