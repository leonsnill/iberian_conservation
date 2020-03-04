#====================================================================================================#
#
# Title: Landsat and Sentinel-2 Image Compositing Tool
# Author: Leon Nill
# Last modified: 2019-06-20
#
#====================================================================================================#

'''
This tool allows for creating pixel-based Landsat image composites mostly based on the
approach of Griffiths et al. (2013): "A Pixel-Based Landsat Compositing Algorithm
for Large Area Land Cover Mapping".

The user can specify the calculation of either spectral-temporal metrics (STMs) (e.g. mean, min, ...)
or pixel-based composites based on scoring functions that determine the suitability of each pixel.

-- User Requirements --
SENSOR               [STRING] – Single sensors or combinations (S2_L1C, S2_L2A, LS, L5, L7, L8, SL)

TARGET_YEARS         [INT] – List of integer years.
SURR_YEARS           INT – 'surrounding years', i.e. should adjacent years be considered for compositing
MONTHLY              BOOLEAN – if True, a monthly iteration is used, if False, iteration is over chosen
                     day of years
SCORE                [STRING] – Score paramater used to create image composite in "qualityMosaic()"-function.
                     ('SCORE', 'NDVI') Selection is based on the maximum of the given parameter, e.g. max NDVI
TARGET_MONTHS_client [INT] – List of target months
STMs                 [ee.Reducer] STMs as ee.Reducer object(s), e.g. ee.Reducer.mean()

ROI                  [xMin, yMin, xMax, yMax] – List of corner coordinates, e.g. [22.26, -19.54, 22.94, -18.89]
ROI_NAME             STRING – Name of the study area which will be used for the output filenames
EPSG                 STRING - Coordinate System !Currently disabled and exports are in WGS84!
PIXEL_RESOLUTION     INT/FLOAT – Output pixelsize in meters

CLOUD_COVER          INT/FLOAT – Maximum cloud cover percentage of scenes considered in pre-selection
BANDS                [STRING] – List of string band-names for export image (B,G,R,NIR,SWIR1,SWIR2,NDVI,TCW, ...)

DOY_RANGE            INT – Offset in days to consider around target doy
REQ_DISTANCE_client  INT – Distance from clouds/ c. shadows in pixels at which optimal conditions are expected
MIN_DISTANCE_client  INT - Minimum distance from clouds/ c. shadows in pixels
W_DOYSCORE_client    FLOAT – Weight of day of year in scoring (0-1)
W_YEARSCORE_client   FLOAT – Weight of year in scoring (0-1)
W_CLOUDSCORE_client  FLOAT – Weight of cloud distance in scoring (0-1)

source activate nillpy
cd '/Users/leonnill/Google Drive/03_LSNRS/Code/01_python/learthengine/learthengine/composite'
python ImgComposite_OKV.py

'''

import ee
ee.Initialize()
import math
import numpy as np


# ====================================================================================================#
# INPUT
# ====================================================================================================#
SENSOR = 'LS'
CLOUD_COVER = 50
BANDS = ['TCB', 'TCG', 'TCW']  # 'B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2'  # 'TCB', 'TCG', 'TCW'
PIXEL_RESOLUTION = 5000

RESAMPLE = None #'bilinear'
REDUCE_RESOLUTION = None #ee.Reducer.mean().unweighted()
NATIVE_RESOLUTION = 30

#ROI = ee.Geometry.Rectangle([-9.3741, 43.8793, 3.5567, 40.8608])

ROI = ee.Geometry.Polygon([[-11,32.679],
                           [30.344,32.679],
                           [30.344,72.414],
                           [-11,72.414],
                           [-11,32.679]])
'''
ROI = ee.Geometry.Polygon([[-9.892,35.884],
                           [3.532,35.884],
                           [3.532,44.051],
                           [-9.892,44.051],
                           [-9.892,35.884]])
#-137.877345, 66.7003, -131.6107, 69.9396 MDR
#22.12, -20.17, 23.69, -18.84 OKV
#43.17, -21.96, 44.37, -21.2 MANGOKY
# 32.4894, 30.0307, 29.8332, 31.6516 NIL
# -9.355, 43.845, 3.399, 41.363 CANTABRIA
'''
ROI_NAME = 'EUROPE'
EPSG = 'EPSG:3035'


# --------------------------------------------------
# TIME
# --------------------------------------------------
TARGET_YEARS = [2016]
SURR_YEARS = 3

MONTHLY = False

if MONTHLY:
       TARGET_MONTHS_client = [1, 4]
       MONTHS_OFFSET = 1
       STMs = [ee.Reducer.median()]
else:
       TARGET_DOY_client = [197]
       # [16, 46, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350]
       DOY_RANGE = 40
       DOY_VS_YEAR = 20

       REQ_DISTANCE_client = 50
       MIN_DISTANCE_client = 10

       W_DOYSCORE_client = 0.7
       W_YEARSCORE_client = 0
       W_CLOUDSCORE_client = 0.3

       SCORE = 'STM_STD'
       BANDNAME = 'TC'

       STMs = [ee.Reducer.stdDev()]


# ====================================================================================================#
# FUNCTIONS
# ====================================================================================================#
# --------------------------------------------------
# RENAME BANDS
# --------------------------------------------------

'''
def fun_rename_bands_l57(img):
       bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
       vnirswir = img.select(bands).multiply(0.0001).rename(new_bands)
       qa = img.select(['pixel_qa'])
       return vnirswir.addBands(qa).copyProperties(img, ['system:time_start'])


def fun_rename_bands_l8(img):
       bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
       vnirswir = img.select(bands).multiply(0.0001).rename(new_bands)
       qa = img.select(['pixel_qa'])
       return vnirswir.addBands(qa).copyProperties(img, ['system:time_start'])

def fun_rename_bands_s2(img):
       bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'QA60', 'SCL']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2', 'QA60', 'SCL']
       return img.select(bands).rename(new_bands)
'''

def fun_rename_bands_l57(img):
       bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']
       return img.select(bands).rename(new_bands)


def fun_rename_bands_l8(img):
       bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'pixel_qa']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']
       return img.select(bands).rename(new_bands)


def fun_rename_bands_s2(img):
       bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'QA60', 'SCL']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2', 'QA60', 'SCL']
       return img.select(bands).rename(new_bands)

'''
def fun_rename_bands_l57(img):
       spectral_bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']
       spectral_bands_name = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
       spectral_img = img.select(spectral_bands).multiply(0.0001).rename(spectral_bands_name)
       add_bands = ['pixel_qa']
       add_img = img.select(add_bands).copyProperties(source=img).set('system:time_start', img.get('system:time_start'))
       return spectral_img.addBands(add_img)


def fun_rename_bands_l8(img):
       spectral_bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
       spectral_bands_name = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
       spectral_img = img.select(spectral_bands).multiply(0.0001).rename(spectral_bands_name)
       add_bands = ['pixel_qa']
       add_img = img.select(add_bands).copyProperties(source=img).set('system:time_start', img.get('system:time_start'))
       return spectral_img.addBands(add_img)


def fun_rename_bands_s2(img):
       spectral_bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
       spectral_bands_name = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
       spectral_img = img.select(spectral_bands).multiply(0.0001).rename(spectral_bands_name)
       add_bands = ['QA60', 'SCL']
       add_img = img.select(add_bands).copyProperties(source=img).set('system:time_start', img.get('system:time_start'))
       return spectral_img.addBands(add_img)
'''

# --------------------------------------------------
# CLOUD MASKING
# --------------------------------------------------


# Function to cloud mask Landsat 8 Surface Reflectance Products
def fun_mask_ls_sr(img):
       cloudShadowBitMask = ee.Number(2).pow(3).int()
       cloudsBitMask = ee.Number(2).pow(5).int()
       snowBitMask = ee.Number(2).pow(4).int()
       qa = img.select('pixel_qa')
       cloud = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(
              qa.bitwiseAnd(cloudsBitMask).eq(0)).rename('CLOUD')
       mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(
              qa.bitwiseAnd(cloudsBitMask).eq(0)).And(
              qa.bitwiseAnd(snowBitMask).eq(0))
       return img.addBands(cloud).updateMask(cloud).divide(10000).float()\
              .copyProperties(source=img).set('system:time_start', img.get('system:time_start'))  # mask --> cloud


# Function to cloud mask Landsat 5 & 7 Surface Reflectance Products (! deprecated !)
def fun_mask_l57_sr(img):
       cloudShadowBitMask = ee.Number(2).pow(3).int()
       cloudsBitMask = ee.Number(2).pow(5).int()
       qa = img.select('pixel_qa')
       mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(
              qa.bitwiseAnd(cloudsBitMask).eq(0)).rename('CLOUD')
       return img.addBands(mask).updateMask(mask).divide(10000).float()\
              .copyProperties(source=img).set('system:time_start', img.get('system:time_start'))


# Function to cloud mask Sentinel-2 L1C / L2A Products
def fun_mask_s2(img):
    qa = img.select('QA60')
    cloudBitMask = ee.Number(2).pow(10).int()
    cirrusBitMask = ee.Number(2).pow(11).int()
    mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0)).rename('CLOUD')
    return img.addBands(mask).updateMask(mask).divide(10000).float()\
              .copyProperties(source=img).set('system:time_start', img.get('system:time_start'))


def fun_mask_s2_scl(img):
    mask_kernel = ee.Kernel.square(radius=5)
    scl = img.select('SCL')
    mask = scl.neq(3).And(scl.neq(7)).And(
            scl.neq(8)).And(
            scl.neq(9)).And(
            scl.neq(10)).And(
            scl.neq(11))
    mask = mask.focal_mode(kernel=mask_kernel, iterations=1).rename('CLOUD')
    return img.addBands(mask).updateMask(mask).divide(10000).float()\
              .copyProperties(source=img).set('system:time_start', img.get('system:time_start'))



# --------------------------------------------------
# ADD BANDS
# --------------------------------------------------
def fun_add_doy_band(img):
       DOY_value = img.date().getRelative('day', 'year')
       DOY = ee.Image.constant(DOY_value).int().rename('DOY')
       DOY = DOY.updateMask(img.select('R').mask())
       return img.addBands(DOY)


def fun_doys(img):
       return ee.Feature(None, {'doy': img.date().getRelative('day', 'year')})


def fun_addyearband(img):
       YEAR_value = ee.Number.parse((img.date().format("YYYY")))
       YEAR = ee.Image.constant(YEAR_value).int().rename('YEAR')
       YEAR = YEAR.updateMask(img.select('R').mask())
       return img.addBands(YEAR)


def fun_addcloudband(img):
       CLOUD_MASK = img.mask().select('R')
       CLOUD_DISTANCE = CLOUD_MASK.Not()\
              .distance(ee.Kernel.euclidean(radius=REQ_DISTANCE, units='pixels'))\
              .rename('CLOUD_DISTANCE')
       CLIP_MAX = CLOUD_DISTANCE.lte(ee.Image.constant(REQ_DISTANCE))
       CLOUD_DISTANCE = CLOUD_DISTANCE.updateMask(CLIP_MAX)
       CLOUD_DISTANCE = CLOUD_DISTANCE.updateMask(CLOUD_MASK)
       return img.addBands(CLOUD_DISTANCE)


# --------------------------------------------------
# SCORING FUNCTIONS
# --------------------------------------------------
def fun_doyscore(img):

       DOYSCORE = img.expression(
              'exp(-0.5*pow((DOY-TARGET_DOY)/DOY_STD, 2))',
              {
                     'DOY': img.select('DOY'),
                     'DOY_STD': DOY_STD,
                     'TARGET_DOY': TARGET_DOY
              }
       ).rename('DOYSCORE')
       DOYSCORE = DOYSCORE.multiply(10000)
       return img.addBands(DOYSCORE)

def fun_doyscore_offset(DOY, TARGET_DOY, DOY_STD):
       return np.exp(-0.5*pow((DOY-TARGET_DOY)/DOY_STD, 2))

def fun_yearscore(img):
       YEAR = ee.Number.parse(img.date().format("YYYY"))
       YEAR_IMG = ee.Algorithms.If(YEAR.eq(TARGET_YEARS_OBJ),
              ee.Image.constant(1).multiply(10000).int().rename('YEARSCORE'),
              ee.Image.constant(DOYSCORE_OFFSET_OBJ).multiply(10000).int().rename('YEARSCORE'))
       return img.addBands(YEAR_IMG)

def fun_cloudscore(img):
       cloud_mask = img.mask().select('R')
       cloud_distance = img.select('CLOUD_DISTANCE')

       img_max = ee.Image.constant(REQ_DISTANCE)
       img_min = ee.Image.constant(MIN_DISTANCE)
       c = img_max.subtract(img_min).divide(ee.Image(2))
       b = cloud_distance.min(img_max)
       a = b.subtract(c).multiply(ee.Image(-0.2)).exp()
       e = ee.Image(1).add(a)

       cldDist = ee.Image(1).divide(e)
       masc_inv = cldDist.mask().Not()
       cldDist = cldDist.mask().where(1, cldDist)
       cldDist = cldDist.add(masc_inv)
       cldDist = cldDist.updateMask(cloud_mask).rename('CLOUDSCORE')
       cldDist = cldDist.multiply(10000)
       return img.addBands(cldDist)

def fun_score(img):
       SCORE = img.expression(
              'DOYSCORE*W_DOYSCORE + YEARSCORE*W_YEARSCORE + CLOUDSCORE*W_CLOUDSCORE',
              {
                     'DOYSCORE': img.select('DOYSCORE'),
                     'YEARSCORE': img.select('YEARSCORE'),
                     'CLOUDSCORE': img.select('CLOUDSCORE'),
                     'W_DOYSCORE': W_DOYSCORE,
                     'W_YEARSCORE': W_YEARSCORE,
                     'W_CLOUDSCORE': W_CLOUDSCORE
              }
       ).rename('SCORE')
       return img.addBands(SCORE)


# --------------------------------------------------
# INDICES
# --------------------------------------------------
def fun_ndvi(img):
       ndvi = img.normalizedDifference(['NIR', 'R']).rename('NDVI')
       #ndvi = ndvi.multiply(10000)
       return img.addBands(ndvi)


def fun_ndwi1(img):
       ndwi = img.normalizedDifference(['NIR', 'SWIR1']).rename('NDWI1')
       #ndwi = ndwi.multiply(10000)
       return img.addBands(ndwi)


def fun_ndwi2(img):
       ndwi = img.normalizedDifference(['G', 'NIR']).rename('NDWI2')
       #ndwi = ndwi.multiply(10000)
       return img.addBands(ndwi)

# Tasseled Cap Transformation (brightness, greenness, wetness) based on Christ 1985
def fun_tcg(img):
    tcg = img.expression(
                         'B*(-0.1603) + G*(-0.2819) + R*(-0.4934) + NIR*0.7940 + SWIR1*(-0.0002) + SWIR2*(-0.1446)',
                         {
                         'B': img.select(['B']),
                         'G': img.select(['G']),
                         'R': img.select(['R']),
                         'NIR': img.select(['NIR']),
                         'SWIR1': img.select(['SWIR1']),
                         'SWIR2': img.select(['SWIR2'])
                         }).rename('TCG')
    #tcg = tcg.multiply(10000)
    return img.addBands(tcg)

def fun_tcb(img):
    tcb = img.expression(
                         'B*0.2043 + G*0.4158 + R*0.5524 + NIR*0.5741 + SWIR1*0.3124 + SWIR2*0.2303',
                         {
                         'B': img.select(['B']),
                         'G': img.select(['G']),
                         'R': img.select(['R']),
                         'NIR': img.select(['NIR']),
                         'SWIR1': img.select(['SWIR1']),
                         'SWIR2': img.select(['SWIR2'])
                         }).rename('TCB')
    #tcb = tcb.multiply(10000)
    return img.addBands(tcb)

def fun_tcw(img):
       tcw = img.expression(
              'B*0.0315 + G*0.2021 + R*0.3102 + NIR*0.1594 + SWIR1*(-0.6806) + SWIR2*(-0.6109)',
              {
                     'B': img.select(['B']),
                     'G': img.select(['G']),
                     'R': img.select(['R']),
                     'NIR': img.select(['NIR']),
                     'SWIR1': img.select(['SWIR1']),
                     'SWIR2': img.select(['SWIR2'])
              }).rename('TCW')
       #tcw = tcw.multiply(10000)
       return img.addBands(tcw)


# ====================================================================================================#
# EXECUTE
# ====================================================================================================#
for year in TARGET_YEARS:
       if MONTHLY:
              year_min = year - SURR_YEARS
              year_max = year + SURR_YEARS

              for month in TARGET_MONTHS_client:

                     months_min = month - MONTHS_OFFSET
                     if months_min < 1:
                            months_min = 1
                     months_max = month + MONTHS_OFFSET
                     if months_max > 12:
                            months_max = 12

                     # --------------------------------------------------
                     # IMPORT ImageCollections
                     # --------------------------------------------------
                     imgCol_L5_SR = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(months_min, months_max, 'month')) \
                            .filter(ee.Filter.lt('CLOUD_COVER_LAND', CLOUD_COVER)) \
                            .map(fun_rename_bands_l57) \
                            .map(fun_mask_ls_sr)

                     imgCol_L7_SR = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(months_min, months_max, 'month')) \
                            .filter(ee.Filter.lt('CLOUD_COVER_LAND', CLOUD_COVER)) \
                            .map(fun_rename_bands_l57) \
                            .map(fun_mask_ls_sr)

                     imgCol_L8_SR = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(months_min, months_max, 'month')) \
                            .filter(ee.Filter.lt('CLOUD_COVER_LAND', CLOUD_COVER)) \
                            .map(fun_rename_bands_l8) \
                            .map(fun_mask_ls_sr)

                     imgCol_S2_L1C = ee.ImageCollection('COPERNICUS/S2') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(months_min, months_max, 'month')) \
                            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_COVER)) \
                            .map(fun_rename_bands_s2) \
                            .map(fun_mask_s2)

                     imgCol_S2_L2A = ee.ImageCollection('COPERNICUS/S2_SR') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(months_min, months_max, 'month')) \
                            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_COVER)) \
                            .map(fun_rename_bands_s2) \
                            .map(fun_mask_s2_scl)

                     # --------------------------------------------------
                     # MERGE imgCols
                     # --------------------------------------------------
                     if SENSOR == 'S2_L1C':
                            imgCol_SR = imgCol_S2_L1C
                     elif SENSOR == 'S2_L2A':
                            imgCol_SR = imgCol_S2_L2A
                     elif SENSOR == 'LS':
                            imgCol_SR = imgCol_L5_SR.merge(imgCol_L7_SR).merge(imgCol_L8_SR)
                            imgCol_SR = imgCol_SR.sort("system:time_start")
                     elif SENSOR == 'L8':
                            imgCol_SR = imgCol_L8_SR
                     elif SENSOR == 'L7':
                            imgCol_SR = imgCol_L7_SR
                     elif SENSOR == 'L5':
                            imgCol_SR = imgCol_L5_SR
                     elif SENSOR == 'SL8':
                            imgCol_SR = imgCol_L8_SR.merge(imgCol_S2_L2A)
                     elif SENSOR == 'SL':
                            imgCol_SR = imgCol_L5_SR.merge(imgCol_L7_SR).merge(imgCol_L8_SR).merge(imgCol_S2_L2A)
                     else:
                            imgCol_SR = None
                            print('No sensor specified!')

                     # --------------------------------------------------
                     # Calculate Indices
                     # --------------------------------------------------
                     imgCol_SR = imgCol_SR.map(fun_ndvi)
                     imgCol_SR = imgCol_SR.map(fun_ndwi1)
                     imgCol_SR = imgCol_SR.map(fun_ndwi2)
                     imgCol_SR = imgCol_SR.map(fun_tcg)
                     imgCol_SR = imgCol_SR.map(fun_tcb)
                     imgCol_SR = imgCol_SR.map(fun_tcw)

                     # --------------------------------------------------
                     # Add DOY, YEAR & CLOUD Bands to ImgCol
                     # --------------------------------------------------
                     imgCol_SR = imgCol_SR.map(fun_add_doy_band)
                     imgCol_SR = imgCol_SR.map(fun_addyearband)

                     # --------------------------------------------------
                     # SPECTRAL TEMPORAL METRICS
                     # --------------------------------------------------
                     for i in range(len(STMs)):
                            if i == 0:
                                   PAR = ee.Image(imgCol_SR.select(BANDS).reduce(STMs[i]))
                            else:
                                   PAR = PAR.addBands(ee.Image(imgCol_SR.select(BANDS).reduce(STMs[i])))

                     # Resample
                     if RESAMPLE:
                            PAR = PAR.resample(RESAMPLE).reproject(crs=EPSG, scale=PIXEL_RESOLUTION)
                     if REDUCE_RESOLUTION:
                            PAR = PAR.reduceResolution(reducer=REDUCE_RESOLUTION, maxPixels=1024) \
                                     .reproject(crs=EPSG, scale=PIXEL_RESOLUTION)
                     
                     #PAR = ee.Image(imgCol_SR.select(BANDS).reduce(ee.Reducer.median()));
                     PAR = PAR.multiply(10000)
                     PAR = PAR.int16()

                     if year_max-year_min == 0:
                            year_filename = str(year)
                     else:
                            year_filename = str(year_min)+'-'+str(year_max)

                     if months_max-months_min == 0:
                            out_file = SENSOR + '_STMs_' + ROI_NAME + '_' + str(PIXEL_RESOLUTION) + 'm_' + year_filename + '_' + str(month)
                     elif months_max-months_min == 11:
                            out_file = SENSOR + '_STMs_' + ROI_NAME + '_' + str(PIXEL_RESOLUTION) + 'm_' + year_filename
                     else:
                            out_file = SENSOR + '_STMs_' + ROI_NAME + '_' + str(PIXEL_RESOLUTION) + 'm_' + year_filename + '_' + str(months_min)+'-' \
                                    + str(months_max)

                     out = ee.batch.Export.image.toDrive(image=PAR, description=out_file,
                                                         scale=PIXEL_RESOLUTION,
                                                         maxPixels=1e13,
                                                         region=ROI['coordinates'][0],
                                                         crs=EPSG)
                     process = ee.batch.Task.start(out)
       else:
              for i in range(len(TARGET_DOY_client)):
                     # --------------------------------------------------
                     # Prepare variables
                     # --------------------------------------------------
                     # Define considered time intervals (min & max)
                     iter_target_doy = TARGET_DOY_client[i]
                     iter_target_doy_min = iter_target_doy - DOY_RANGE
                     if iter_target_doy_min < 1:
                            iter_target_doy_min = 1

                     iter_target_doy_max = iter_target_doy + DOY_RANGE
                     if iter_target_doy_max > 365:
                            iter_target_doy_max = 365

                     year_min = year - SURR_YEARS
                     year_max = year + SURR_YEARS

                     REQ_DISTANCE = ee.Number(REQ_DISTANCE_client)
                     MIN_DISTANCE = ee.Number(MIN_DISTANCE_client)

                     # --------------------------------------------------
                     # IMPORT ImageCollections
                     # --------------------------------------------------
                     imgCol_L5_SR = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(iter_target_doy_min, iter_target_doy_max, 'day_of_year')) \
                            .filter(ee.Filter.lt('CLOUD_COVER_LAND', CLOUD_COVER)) \
                            .map(fun_rename_bands_l57) \
                            .map(fun_mask_ls_sr)

                     imgCol_L7_SR = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(iter_target_doy_min, iter_target_doy_max, 'day_of_year')) \
                            .filter(ee.Filter.lt('CLOUD_COVER_LAND', CLOUD_COVER)) \
                            .map(fun_rename_bands_l57) \
                            .map(fun_mask_ls_sr)

                     imgCol_L8_SR = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(iter_target_doy_min, iter_target_doy_max, 'day_of_year')) \
                            .filter(ee.Filter.lt('CLOUD_COVER_LAND', CLOUD_COVER)) \
                            .map(fun_rename_bands_l8) \
                            .map(fun_mask_ls_sr)

                     imgCol_S2_L1C = ee.ImageCollection('COPERNICUS/S2') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(iter_target_doy_min, iter_target_doy_max, 'day_of_year')) \
                            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_COVER)) \
                            .map(fun_rename_bands_s2) \
                            .map(fun_mask_s2)

                     imgCol_S2_L2A = ee.ImageCollection('COPERNICUS/S2_SR') \
                            .filterBounds(ROI) \
                            .filter(ee.Filter.calendarRange(year_min, year_max, 'year')) \
                            .filter(ee.Filter.calendarRange(iter_target_doy_min, iter_target_doy_max, 'day_of_year')) \
                            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_COVER)) \
                            .map(fun_rename_bands_s2) \
                            .map(fun_mask_s2_scl)

                     # --------------------------------------------------
                     # MERGE imgCols
                     # --------------------------------------------------
                     if SENSOR == 'S2_L1C':
                            imgCol_SR = imgCol_S2_L1C
                     elif SENSOR == 'S2_L2A':
                            imgCol_SR = imgCol_S2_L2A
                     elif SENSOR == 'LS':
                            imgCol_SR = imgCol_L5_SR.merge(imgCol_L7_SR).merge(imgCol_L8_SR)
                            imgCol_SR = imgCol_SR.sort("system:time_start")
                     elif SENSOR == 'L8':
                            imgCol_SR = imgCol_L8_SR
                     elif SENSOR == 'L7':
                            imgCol_SR = imgCol_L7_SR
                     elif SENSOR == 'L5':
                            imgCol_SR = imgCol_L5_SR
                     elif SENSOR == 'SL8':
                            imgCol_SR = imgCol_L8_SR.merge(imgCol_S2_L2A)
                     elif SENSOR == 'SL':
                            imgCol_SR = imgCol_L5_SR.merge(imgCol_L7_SR).merge(imgCol_L8_SR).merge(imgCol_S2_L2A)
                     else:
                            imgCol_SR = None
                            print('No sensor specified!')

                     # --------------------------------------------------
                     # Calculate Indices
                     # --------------------------------------------------
                     imgCol_SR = imgCol_SR.map(fun_ndvi)
                     imgCol_SR = imgCol_SR.map(fun_ndwi1)
                     imgCol_SR = imgCol_SR.map(fun_ndwi2)
                     imgCol_SR = imgCol_SR.map(fun_tcg)
                     imgCol_SR = imgCol_SR.map(fun_tcb)
                     imgCol_SR = imgCol_SR.map(fun_tcw)

                     # --------------------------------------------------
                     # Add DOY, YEAR & CLOUD Bands to ImgCol
                     # --------------------------------------------------
                     imgCol_SR = imgCol_SR.map(fun_add_doy_band)
                     imgCol_SR = imgCol_SR.map(fun_addyearband)
                     imgCol_SR = imgCol_SR.map(fun_addcloudband)
                         
                     if SCORE == 'SCORE':
                            # --------------------------------------------------
                            # SCORING 1: DOY
                            # --------------------------------------------------
                            # add DOY-band to images in imgCol
                            DOYs = imgCol_SR.map(fun_doys).aggregate_array('doy').getInfo()

                            # retrieve target-DOY and DOY-std (client and server side)
                            # TARGET_DOY_client = int(np.median(DOYs))
                            TARGET_DOY = ee.Number(iter_target_doy)

                            DOY_STD_client = np.std(DOYs)
                            #DOY_STD_client = DOY_RANGE

                            DOY_STD = ee.Number(DOY_STD_client)

                            # add Band with final DOY score to every image in imgCol
                            imgCol_SR = imgCol_SR.map(fun_doyscore)

                            # --------------------------------------------------
                            # SCORING 2: YEAR
                            # --------------------------------------------------
                            # calculate DOY-score at maximum DOY vs Year threshold
                            DOYSCORE_OFFSET = fun_doyscore_offset(iter_target_doy - DOY_VS_YEAR,
                                                                  iter_target_doy, DOY_STD_client)
                            DOYSCORE_OFFSET_OBJ = ee.Number(DOYSCORE_OFFSET)
                            TARGET_YEARS_OBJ = ee.Number(year)

                            # add Band with final YEAR score to every image in imgCol
                            imgCol_SR = imgCol_SR.map(fun_yearscore)

                            # --------------------------------------------------
                            # SCORING 3: CLOUD DISTANCE
                            # --------------------------------------------------
                            imgCol_SR = imgCol_SR.map(fun_cloudscore)

                            # --------------------------------------------------
                            # FINAL SCORING
                            # --------------------------------------------------
                            W_DOYSCORE = ee.Number(W_DOYSCORE_client)
                            W_YEARSCORE = ee.Number(W_YEARSCORE_client)
                            W_CLOUDSCORE = ee.Number(W_CLOUDSCORE_client)

                            imgCol_SR = imgCol_SR.map(fun_score)

                            COMPOSITE = imgCol_SR.qualityMosaic(SCORE)
                            COMPOSITE = COMPOSITE.select(BANDS)
                            COMPOSITE = COMPOSITE.multiply(10000)
                            COMPOSITE = COMPOSITE.int16()
                                
                            if STMs is not None:
                                   for i in range(len(STMs)):
                                          COMPOSITE = COMPOSITE.addBands(ee.Image(imgCol_SR.select(BANDS).reduce(STMs[i])).int16())

                     elif SCORE == 'MAXNDVI':
                            COMPOSITE = imgCol_SR.qualityMosaic('NDVI')
                            COMPOSITE = COMPOSITE.select(BANDS)
                            COMPOSITE = COMPOSITE.multiply(10000)
                            COMPOSITE = COMPOSITE.int16()

                            if STMs is not None:
                                   for i in range(len(STMs)):
                                          COMPOSITE = COMPOSITE.addBands(ee.Image(imgCol_SR.select(BANDS).reduce(STMs[i])).int16())
                     
                     else:
                            if STMs is not None:
                                   for i in range(len(STMs)):
                                          if i == 0:
                                                 COMPOSITE = ee.Image(imgCol_SR.select(BANDS).reduce(STMs[i]))
                                          else:
                                                 COMPOSITE = COMPOSITE.addBands(ee.Image(imgCol_SR.select(BANDS).reduce(STMs[i])))
                                   # Resample
                                   if RESAMPLE:
                                          COMPOSITE = COMPOSITE.resample(RESAMPLE)
                                   if REDUCE_RESOLUTION:
                                          maxPixels_factor = math.ceil(PIXEL_RESOLUTION / NATIVE_RESOLUTION)
                                          COMPOSITE = COMPOSITE.reproject(crs=EPSG, scale=NATIVE_RESOLUTION)
                                          COMPOSITE = COMPOSITE.reduceResolution(reducer=REDUCE_RESOLUTION, bestEffort=False, maxPixels=maxPixels_factor*maxPixels_factor)

                                   
                                   COMPOSITE = COMPOSITE.multiply(10000)

                     if SURR_YEARS == 0:
                            year_filename = str(year)
                     else:
                            year_filename = str(year)+'-'+str(SURR_YEARS)

                     out_file = SENSOR + '_imgComposite_' + SCORE + '_' + BANDNAME + '_' + ROI_NAME + '_' + str(PIXEL_RESOLUTION) + 'm_' + year_filename + '_' + str(iter_target_doy)

                     out = ee.batch.Export.image.toDrive(image=COMPOSITE.toInt16(), description=out_file,
                                                         scale=PIXEL_RESOLUTION,
                                                         maxPixels=1e13,
                                                         region=ROI['coordinates'][0],
                                                         crs=EPSG)
                     process = ee.batch.Task.start(out)


# =====================================================================================================================#
# END
# =====================================================================================================================#
