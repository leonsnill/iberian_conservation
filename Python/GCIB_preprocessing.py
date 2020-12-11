# ======================================================================================================================
#
# GCIB Project PREPROCESSING
# Leon Nill, 2020-01-17
#
# ======================================================================================================================


# ----------------------------------------------------------------------------------------------------------------------
# IMPORT
# ----------------------------------------------------------------------------------------------------------------------
import os
import gdal
import numpy as np
from eopy import raster
from eopy import finder


# ----------------------------------------------------------------------------------------------------------------------
# RESAMPLING
# ----------------------------------------------------------------------------------------------------------------------
# reference image
reference_img = '/Users/leonnill/Documents/MSc_GCG/MSc_GCIB/Watermask/EUROPE_MASK_10km.tif'
roi = raster.extent([reference_img])

# LANDSAT STMs
f_landsat = finder('/Users/leonnill/Documents/MSc_GCG/MSc_GCIB/STMs/Europe', pattern='*.tif')

for f in f_landsat:
    inp = gdal.Open(f)
    #out = f.replace('5000m', '5km')
    #out = f.replace('EUROPE_5km', 'IBERIA_10km')
    #out = f.replace('europe', 'iberia')
    out = f.replace('5km', '10km')
    opt = gdal.TranslateOptions(projWin=roi, resampleAlg='average', noData=None, xRes=10000, yRes=10000)
    gdal.Translate(out, inp, options=opt)


# CHELSA BIOCLIM
f_bioclim = finder('/Users/leonnill/Documents/MSc_GCG/MSc_GCIB/CHELSA_BIOCLIM/EUROPE', pattern='*.tif')

for f in f_bioclim:
    inp = gdal.Open(f)
    out = f.replace('.tif', '_EUROPE_10km.tif')
    raster.reproject_raster(inp, 'average', out, ref=gdal.Open(reference_img))


# DEM
f_dem = '/Users/leonnill/Documents/MSc_GCG/MSc_GCIB/DEM/GTOPO30_1000m_EUROPE.tif'
out = f_dem.replace('1000m', '10km.tif')
opt = gdal.TranslateOptions(projWin=roi, resampleAlg='average',
                            noData=0, xRes=10000, yRes=10000, outputType=gdal.GDT_Int16)
gdal.Translate(out, gdal.Open(f_dem), options=opt)


# ----------------------------------------------------------------------------------------------------------------------
# LAND COVER TO FRACTION COVER
# ----------------------------------------------------------------------------------------------------------------------
#os.chdir('C:/Users/Leon/Google Drive/01_MSc_GCG/MSc6_GCIB/ursus_arctos')
# reference image
reference_img = r"Data\Mask\IBERIA_MASK_10km.tif"
xmin, ymax, xmax, ymin = raster.extent([gdal.Open(reference_img)])

warp = gdal.Warp(r"Data\Mask\roads-train-artificial_iberia_fraction_10km.tif", r"Data\Mask\roads-train-artificial_iberia_mask.tif", xRes=10000, yRes=10000,
          outputBounds=[xmin, ymin, xmax, ymax], resampleAlg='average', outputType=gdal.GDT_Int16, srcNodata=-9999)
warp = None



# read in lc-dataset and create binary mask from chosen lc-classes
lc_arr = gdal.Open(r"Data\Mask\artificial_mask.tif").ReadAsArray()
binary = np.where((lc_arr == 4) | (lc_arr == 5) | (lc_arr == 6), 1000, 0)
binary = binary.astype(np.int)
raster.array_to_geotiff(binary, outpath+'/binary.tif', inp_gdal=lc)

# open binary file and resample to target resolution
binary = gdal.Open(outpath+'/binary.tif')
opt = gdal.TranslateOptions(xRes=10000, yRes=10000, resampleAlg='average', outputType=gdal.GDT_Int16, noData=None)
gdal.Translate(outpath+'/iberia_lc_fraction_highveg_10km.tif', binary, options=opt)
