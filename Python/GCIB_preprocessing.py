# ======================================================================================================================
#
# GCIB Project PREPROCESSING
# Leon Nill, 2020-01-17
#
# ======================================================================================================================


# ----------------------------------------------------------------------------------------------------------------------
# IMPORT
# ----------------------------------------------------------------------------------------------------------------------
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
# reference image
reference_img = '/Users/leonnill/Documents/MSc_GCG/MSc_GCIB/Watermask/IBERIA_MASK_5km.tif'
roi = raster.extent([reference_img])

outpath = '/Users/leonnill/Documents/MSc_GCG/MSc_GCIB/LC_Fractions'
lc_file = '/Users/leonnill/Documents/LSNGIS/Data/Europe_LC_2015_Pflugmacher/europe_landcover_2015_RSE-Full3.tif'
lc = gdal.Open(lc_file)

# read in lc-dataset and create binary mask from chosen lc-classes
lc_arr = gdal.Open(lc_file).ReadAsArray()
binary = np.where((lc_arr == 4) | (lc_arr == 5) | (lc_arr == 6), 1000, 0)
binary = binary.astype(np.int)
raster.array_to_geotiff(binary, outpath+'/binary.tif', inp_gdal=lc)

# open binary file and resample to target resolution
binary = gdal.Open(outpath+'/binary.tif')
opt = gdal.TranslateOptions(xRes=10000, yRes=10000, resampleAlg='average', outputType=gdal.GDT_Int16, noData=None)
gdal.Translate(outpath+'/iberia_lc_fraction_highveg_10km.tif', binary, options=opt)



