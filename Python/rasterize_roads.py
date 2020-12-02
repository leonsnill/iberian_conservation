import gdal
import ogr, osr
import numpy as np
import eopy as eo

# rasterize function
def rasterize(layer, target_raster):

    xdim = target_raster.RasterXSize
    ydim = target_raster.RasterYSize
    gt = target_raster.GetGeoTransform()
    proj = target_raster.GetProjectionRef()

    # create target data source
    ds_target = gdal.GetDriverByName('MEM').Create(
        '', xdim, ydim, 1, gdal.GDT_Byte)

    # set transform and projection
    ds_target.SetGeoTransform(gt)
    ds_target.SetProjection(proj)

    img_crs = osr.SpatialReference()
    img_crs.ImportFromWkt(proj)  # raster crs

    gdal.RasterizeLayer(ds_target, [1], layer, burn_values=[1])

    return ds_target


# inp reference raster
file_ref = "Data/Mask/landcover_mask_empty.tif"

# input vector
dir_vector = r"Data\Vector\motorways_iberia.gpkg"

# different files
road_spain = ogr.Open(dir_vector+"/RT_VIARIA_CARRETERA/rt_tramo_vial.shp")
lyr_road_spain = road_spain.GetLayer()
road_port = ogr.Open(dir_vector+"/portugal-roads-shape/roads.shp")
lyr_road_port = road_port.GetLayer()

train_spain = ogr.Open(dir_vector+"/RT_FFCC/rt_tramofc_linea.shp")
lyr_train_spain = train_spain.GetLayer()
train_port = ogr.Open(dir_vector+"/portugal-railways-shape/railways.shp")
lyr_train_port = train_port.GetLayer()

# (1) railway raster
raster_train_spain = rasterize(lyr_train_spain, gdal.Open(file_ref))
raster_train_port = rasterize(lyr_train_port, gdal.Open(file_ref))

train_mask = np.logical_or(raster_train_spain.ReadAsArray(), raster_train_port.ReadAsArray())
eo.raster.array_to_geotiff(train_mask, "Data/Mask/train_mask.tif",
                           inp_img=gdal.Open(file_ref), nodata=0, compress=True)


attr = "ClaseD"
value = "Autopista"

# filter by attribute
lyr_road_spain.SetAttributeFilter("ClaseD = 'Autopista' or ClaseD = 'Autovia'")

# rasterize to layer
raster_road_spain = rasterize(lyr_road_spain, gdal.Open(file_ref))
raster_road_port = rasterize(lyr_road_port, gdal.Open(file_ref))

road_mask = np.logical_or(raster_road_spain.ReadAsArray(), raster_road_port.ReadAsArray())
eo.raster.array_to_geotiff(road_mask, "Data/Mask/road_mask.tif",
                           inp_img=gdal.Open(file_ref), nodata=0, compress=True)
