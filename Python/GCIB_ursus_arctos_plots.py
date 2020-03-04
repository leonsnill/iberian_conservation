import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import gdal
import numpy as np
import osr
import matplotlib as mpl

# ----------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------


def convertXY(xy_source, inp_proj, out_proj):
    # function to convert coordinates

    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size

    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inp_proj, out_proj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))

    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)

    return xx, yy


# get extent
def ext_min_max(ds):

    if isinstance(ds, str):
        ds = gdal.Open(ds)

    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()

    xres = gt[1]
    yres = gt[5]
    xmin = gt[0] #+ xres * 0.5
    xmax = gt[0] + (xres * ds.RasterXSize) #- xres * 0.5
    ymin = gt[3] + (yres * ds.RasterYSize) #+ yres * 0.5
    ymax = gt[3] #- yres * 0.5

    return (xmin, xmax, ymin, ymax)


def disc_color(cmap, n_distinct, col_list=True):

    cmap = plt.cm.get_cmap(cmap, n_distinct)
    newcolors = cmap(np.linspace(0, 1, n_distinct))

    if col_list:
        newcolors = mpl.colors.ListedColormap(newcolors)

    return newcolors


def cbar_ticks(a):

    a_unq = np.unique(a[~np.isnan(a)])
    a_n_unq = len(a_unq)
    a_min = np.nanmin(a)
    a_max = np.nanmax(a)
    loc_ticks = np.linspace(a_min, a_max, a_n_unq)

    pos_adj = np.abs(loc_ticks[0] - loc_ticks[1]) / 2
    vmin = a_min - pos_adj
    vmax = a_max + pos_adj

    return (a_n_unq, vmin, vmax, loc_ticks)


from matplotlib.patches import Polygon
def draw_screen_poly(ax, lats, lons, m, rec=True, col='red'):
    x, y = m(lons, lats)
    if rec:
        x = [x[0], x[0], x[1], x[1]]
        y = [y[0], y[1], y[1], y[0]]
    xy = zip(x,y)
    poly = Polygon(list(xy), fill=False, edgecolor=col, linewidth=3)
    ax.add_patch(poly)

def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax, extend='both')
    plt.sca(last_axes)
    return cbar

# ----------------------------------------------------------------------------------------------------------------------
# APPLY
# ----------------------------------------------------------------------------------------------------------------------

os.chdir('/Users/leonnill/Documents/MSc_GCG/MSc_GCIB')

# *******************************************************
# PREPARE DATA
# *******************************************************

ds = gdal.Open('SDMs/SDM_ensembles_mn_md_wmn_cav_sdp.tif')
mask = gdal.Open('Watermask/EUROPE_MASK_10km_to_SDMs.tif').ReadAsArray()

# get extent, resolution and projection
xmin, xmax, ymin, ymax = ext_min_max(ds)
xres, yres = ds.GetGeoTransform()[1], ds.GetGeoTransform()[5]
ds_proj = osr.SpatialReference(wkt=ds.GetProjection())

# create meshgrid of coordinates
msh_xy = np.mgrid[xmin:xmax:xres, ymax:ymin:yres]

# prepare data
arr = ds.ReadAsArray()
arr = arr[2]
arr = np.where(arr < -1, np.NaN, arr)
arr = np.where(mask == 0, np.NaN, arr)

# threshold occurence
th = np.where(arr >= 0.52, 1, np.NaN)

# occupied habitat
gbif = gdal.Open("/Users/leonnill/Documents/MSc_GCG/MSc_GCIB/GBIF/GBIF_ursusarctos_europe_1990-2020_presence_10km.tif")
gbif_a = gbif.ReadAsArray()
gbif_a = np.where(gbif_a < 1, np.NaN, gbif_a)
xmin, xmax, ymin, ymax = ext_min_max(gbif)
xres, yres = ds.GetGeoTransform()[1], gbif.GetGeoTransform()[5]
gbif_proj = osr.SpatialReference(wkt=gbif.GetProjection())

# create meshgrid of coordinates
gbif_xy = np.mgrid[xmin:xmax:xres, ymax:ymin:yres]

# *******************************************************
# MAP
# *******************************************************
# font size
mpl.rcParams.update({'font.size': 14})


bblons = [-9.5, 3.8]
bblats = [37.5, 45.7]

# initialise data
fig = plt.figure(figsize=(14, 8))
widths = [2.1, 2]
spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=widths)
ax1 = fig.add_subplot(spec[:, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax2.set_xticks([], [])
ax2.set_yticks([], [])
ax3 = fig.add_subplot(spec[1, 1])
ax3.set_xticks([], [])
ax3.set_yticks([], [])

fig.subplots_adjust(hspace=.1, wspace=.1)

m = Basemap(projection='laea', llcrnrlon=-10, llcrnrlat=34, urcrnrlon=45, urcrnrlat=70,
            resolution='l', epsg=3035, ax=ax1)
m2 = Basemap(projection='laea', llcrnrlon=bblons[0], llcrnrlat=bblats[0], urcrnrlon=bblons[1], urcrnrlat=bblats[1],
            resolution='h', epsg=3035, ax=ax2)
m3 = Basemap(projection='laea', llcrnrlon=bblons[0], llcrnrlat=bblats[0], urcrnrlon=bblons[1], urcrnrlat=bblats[1],
            resolution='h', epsg=3035, ax=ax3)


# Get the target projection from the basemap object
m_proj = osr.SpatialReference()
m_proj.ImportFromProj4(m.proj4string)

m2_proj = osr.SpatialReference()
m2_proj.ImportFromProj4(m2.proj4string)

# Convert from source projection to basemap projection
xx, yy = convertXY(msh_xy, ds_proj, m_proj)
xx2, yy2 = convertXY(msh_xy, ds_proj, m2_proj)
xx_gbif, yy_gbif = convertXY(gbif_xy, gbif_proj, m2_proj)

# color
n_unq, vmin, vmax, loc_ticks = cbar_ticks(arr)



colmap = disc_color('jet', n_unq-2, False)

pink = mpl.colors.to_rgba('lightpink')
grey = np.array([0.5, 0.5, 0.5, 1])
colmap = np.vstack([grey, colmap])
colmap = np.vstack([colmap, pink])
colmap = mpl.colors.ListedColormap(colmap)

red = mpl.colors.to_rgba('red')
colmap2 = np.vstack([grey, red])
colmap2 = mpl.colors.ListedColormap(colmap2)


# plot coastlines, meridians etc.
m.drawcoastlines(linewidth=1)
m.drawcountries(linewidth=1)
m.fillcontinents(color='lightgrey', zorder=0)
parallels = np.arange(-45, 90, 10)
meridians = np.arange(0, 360, 10)
m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=14)
m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=14)

parallels = np.arange(-45, 90, 5)
meridians = np.arange(0, 360, 5)
m2.drawcoastlines(linewidth=1)
m2.drawcountries(linewidth=1)
m2.fillcontinents(color='lightgrey', zorder=0)
m2.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=14)
m2.drawparallels(parallels, labels=[0, 1, 0, 0], fontsize=14)

m3.drawcoastlines(linewidth=1)
m3.drawcountries(linewidth=1)
m3.fillcontinents(color='lightgrey', zorder=0)
m3.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=14)
m3.drawparallels(parallels, labels=[0, 1, 0, 0], fontsize=14)

# plot the data (first layer)
#im1 = m.pcolormesh(xx, yy, arr.T, cmap=colmap, vmin=vmin, vmax=vmax)
im1 = m.pcolormesh(xx, yy, arr.T, cmap='jet', vmin=0.1, vmax=0.9)
im2 = m2.pcolormesh(xx2, yy2, arr.T, cmap='jet', vmin=0.1, vmax=0.9)
im3 = m3.pcolormesh(xx2, yy2, th.T, cmap=colmap, vmin=0, vmax=1)
im4 = m3.pcolormesh(xx_gbif, yy_gbif, gbif_a.T, cmap=colmap2, vmin=0, vmax=1)


cbar = colorbar(im1)

from matplotlib.font_manager import FontProperties
font0 = FontProperties()
font = font0.copy()
font1 = font0.copy()
font.set_style('normal')
font.set_weight('bold')
font.set_size(24)
font1.set_style('normal')
font1.set_weight('bold')
font1.set_size(14)

# Figure Sub Names
ax1.text(0.05, 0.95, 'A', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax1.transAxes, fontsize=20,
        horizontalalignment='center', verticalalignment='center', fontproperties=font)
ax2.text(0.05, 0.91, 'B', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax2.transAxes, fontsize=20,
        horizontalalignment='center', verticalalignment='center', fontproperties=font)
ax3.text(0.05, 0.91, 'C', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax3.transAxes, fontsize=20,
        horizontalalignment='center', verticalalignment='center', fontproperties=font)

ax1.text(0.5, 0.95, 'Habitat Suitability', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax1.transAxes, fontsize=14, horizontalalignment='center', fontproperties=font1)
ax1.text(0.5, 0.91, 'Ursus Arctos', style='italic',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax1.transAxes, fontsize=14, horizontalalignment='center')

#cbar = fig.colorbar(im1, ax=ax1)

#cbar.set_ticks(loc_ticks)
#cbar.ax.set_yticklabels(['Decolonized', 'Stable', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'Suitable'])

# map scale
#m.drawmapscale(1.4, 37, -3.25, 39.5, 500, barstyle='fancy', fontsize=14)
m.drawmapboundary(fill_color=None,zorder=0)
m2.drawmapboundary(fill_color=None,zorder=0)
m3.drawmapboundary(fill_color=None,zorder=0)

draw_screen_poly(ax1, bblats, bblons, m, rec=True)

labels = ['Occupied Habitat', 'Suitable Habitat [> 0.52]']
colors = ['red', 'lightpink']
patches = [
    mpl.patches.Patch(color=color, label=label)
    for label, color in zip(labels, colors)]
ax3.legend(patches, labels, loc='lower left', frameon=True, fontsize=14,
           framealpha=1, edgecolor='black', fancybox=False)


# plot / save figure
plt.savefig('/Users/leonnill/Desktop/Ursus_Arctos_SDM.png', dpi=300)
plt.close()
#plt.show()
#plt.close()





# ------------------------------
# DEM
# ------------------------------
dem = gdal.Open('DEM/GTOPO30_1000m_EUROPE_v2.tif')
#dem = gdal.Open('DEM/GTOPO30_10km.tif_EUROPE.tif')
dem_a = dem.ReadAsArray()
dem_a = np.where(dem_a <= 0, np.NaN, dem_a)

xmin, xmax, ymin, ymax = ext_min_max(dem)
xres, yres = dem.GetGeoTransform()[1], dem.GetGeoTransform()[5]
dem_proj = osr.SpatialReference(wkt=dem.GetProjection())

# create meshgrid of coordinates
dem_xy = np.mgrid[xmin:xmax:xres, ymax:ymin:yres]

# ------------------------------
# Presence - Absence
# ------------------------------
pre = gdal.Open('GBIF/UrsusArctos_EUROPE_10km_presence_thinned.tif')
pre_a = pre.ReadAsArray()
abs = gdal.Open('GBIF/UrsusArctos_EUROPE_10km_absence_thinned.tif')
abs_a = abs.ReadAsArray()

pre_a = np.where(pre_a < 1, np.NaN, pre_a)
abs_a = np.where(abs_a < 1, np.NaN, abs_a)

xmin, xmax, ymin, ymax = ext_min_max(pre)
xres, yres = pre.GetGeoTransform()[1], pre.GetGeoTransform()[5]
pre_proj = osr.SpatialReference(wkt=pre.GetProjection())

# create meshgrid of coordinates
pre_xy = np.mgrid[xmin:xmax:xres, ymax:ymin:yres]





# font size
mpl.rcParams.update({'font.size': 14})
mpl.rcParams['figure.dpi'] = 300

bblons = [-9.5, 3.8]
bblats = [35, 45.7]

# initialise data
fig = plt.figure(figsize=(17, 8))
widths = [1, 1]
spec = fig.add_gridspec(ncols=2, nrows=5, width_ratios=widths)
ax = fig.add_subplot(spec[:, 0])
ax2 = fig.add_subplot(spec[:4, 1])

ax3 = fig.add_subplot(spec[4, 1])
ax3.set_xticks([], [])
ax3.set_yticks([], [])
ax3.set_frame_on(False)

fig.subplots_adjust(hspace=.2, wspace=.05)



m = Basemap(projection='laea', llcrnrlon=-10, llcrnrlat=34, urcrnrlon=45, urcrnrlat=70,
            resolution='l', epsg=3035, ax=ax)
m2 = Basemap(projection='laea', llcrnrlon=bblons[0], llcrnrlat=bblats[0], urcrnrlon=bblons[1], urcrnrlat=bblats[1],
            resolution='h', epsg=3035, ax=ax2)
# Get the target projection from the basemap object
m_proj = osr.SpatialReference()
m_proj.ImportFromProj4(m.proj4string)

m2_proj = osr.SpatialReference()
m2_proj.ImportFromProj4(m2.proj4string)

# Convert from source projection to basemap projection
xx, yy = convertXY(dem_xy, dem_proj, m_proj)
xx2, yy2 = convertXY(pre_xy, pre_proj, m_proj)

xx_m2, yy_m2 = convertXY(dem_xy, dem_proj, m2_proj)
xx2_m2, yy2_m2 = convertXY(pre_xy, pre_proj, m2_proj)

# plot coastlines, meridians etc.
m.drawcoastlines(linewidth=1)
m.drawcountries(linewidth=1)
m.fillcontinents(color='lightgrey', zorder=0)
parallels = np.arange(-45, 90, 10)
meridians = np.arange(0, 360, 10)
m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=14)
m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=14)

parallels = np.arange(-45, 90, 5)
meridians = np.arange(0, 360, 5)
m2.drawcoastlines(linewidth=1)
m2.drawcountries(linewidth=1)
m2.fillcontinents(color='lightgrey', zorder=0)
m2.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=14)
m2.drawparallels(parallels, labels=[0, 1, 0, 0], fontsize=14)

im1 = m.pcolormesh(xx, yy, dem_a.T, cmap='terrain', vmin=0, vmax=2000)
im1_m2 = m2.pcolormesh(xx_m2, yy_m2, dem_a.T, cmap='terrain', vmin=0, vmax=2000)
cbar = colorbar(im1)
cbar.set_label('m a.s.l.')

red = mpl.colors.to_rgba('red')
grey = mpl.colors.to_rgba('black')
colmap = np.vstack([grey, red])
colmap = mpl.colors.ListedColormap(colmap)
im2 = m.pcolormesh(xx2, yy2, pre_a.T, cmap=colmap, vmin=0, vmax=1)
im2_m2 = m2.pcolormesh(xx2_m2, yy2_m2, pre_a.T, cmap=colmap, vmin=0, vmax=1)

im3 = m.pcolormesh(xx2, yy2, abs_a.T, cmap=colmap)
im3_m2 = m2.pcolormesh(xx2_m2, yy2_m2, abs_a.T, cmap=colmap)

m.drawmapboundary(fill_color=None,zorder=0)
m2.drawmapboundary(fill_color=None,zorder=0)

draw_screen_poly(ax, bblats, bblons, m, rec=True)

bblons2 = [-7.3, -3.6]
bblats2 = [41.8, 44.2]
draw_screen_poly(ax2, bblats2, bblons2, m2, rec=True, col='purple')


bblons3 = [-1.3, 2.5]
bblats3 = [41.7, 43.6]
draw_screen_poly(ax2, bblats3, bblons3, m2, rec=True, col='black')

m2.drawmapscale(1.3, 37.3, 4, 40, 500, barstyle='fancy', fontsize=14)

labels = ['Iberian Peninsula', 'Cantabrian Population', 'Pyrenean Population']
colors = ['red', 'purple', 'black']
patches = [
    mpl.patches.Patch(edgecolor=color, label=label, fill=False, lw=3)
    for label, color in zip(labels, colors)]
leg1 = ax3.legend(patches, labels, loc='upper right', frameon=True, fontsize=16,
           framealpha=0, edgecolor='black', fancybox=False)
labels = ['Presence', 'Absence']
colors = ['red', 'black']
patches = [
    mpl.patches.Patch(color=color, label=label, lw=3)
    for label, color in zip(labels, colors)]
ax3.legend(patches, labels, loc='upper left', frameon=True, fontsize=16,
           framealpha=0, edgecolor='black', fancybox=False)

ax3.add_artist(leg1)


from matplotlib.font_manager import FontProperties
font0 = FontProperties()
font = font0.copy()
font.set_style('normal')
font.set_weight('bold')
font.set_size(24)
# Figure Sub Names
ax.text(0.05, 0.95, 'A', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax.transAxes, fontsize=20,
        horizontalalignment='center', verticalalignment='center', fontproperties=font)
ax2.text(0.05, 0.95, 'B', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax2.transAxes, fontsize=20,
        horizontalalignment='center', verticalalignment='center', fontproperties=font)


#plt.show()
plt.savefig('/Users/leonnill/Desktop/Study_Area_Pre_Abs.png', dpi=fig.dpi, bbox_inches='tight')
plt.close()


# ----------------------------------------------------------------------------------------------------------------------
# DISPERSAL
# ----------------------------------------------------------------------------------------------------------------------


fart = gdal.Open('LC_Fractions/europe_lc_fraction_artificial_10km.tif')
xmin, xmax, ymin, ymax = ext_min_max(fart)
xres, yres = fart.GetGeoTransform()[1], fart.GetGeoTransform()[5]
fart_proj = osr.SpatialReference(wkt=fart.GetProjection())

# create meshgrid of coordinates
fart_xy = np.mgrid[xmin:xmax:xres, ymax:ymin:yres]
fart_a = fart.ReadAsArray()/1000
mwy = gdal.Open('Watermask/EUROPE_MASK_10km_MOTORWAYS.tif')
mwy_a = mwy.ReadAsArray()

disp_mask = np.where(fart_a > 0.05, 2, 0)
disp_mask = np.where(mwy_a == 2, disp_mask+1, disp_mask)
disp_mask = np.where(disp_mask == 0, np.NaN, disp_mask)


disp = gdal.Open('dispersal_ursus_arctos_10a_nobarrier.tif')
disp_fart = gdal.Open('dispersal_ursus_arctos_10a_onlylandbarrier.tif')
disp_mwyfart = gdal.Open('dispersal_ursus_arctos_10a_barrier.tif')

# get extent, resolution and projection
xmin, xmax, ymin, ymax = ext_min_max(disp)
xres, yres = disp.GetGeoTransform()[1], disp.GetGeoTransform()[5]
disp_proj = osr.SpatialReference(wkt=disp.GetProjection())

# create meshgrid of coordinates
disp_xy = np.mgrid[xmin:xmax:xres, ymax:ymin:yres]

disp_a = disp.ReadAsArray()
disp_a = np.where(disp_a<-1, np.NaN, disp_a)
disp_fart_a = disp_fart.ReadAsArray()
disp_fart_a = np.where(disp_fart_a<-1, np.NaN, disp_fart_a)
disp_mwyfart_a = disp_mwyfart.ReadAsArray()
disp_mwyfart_a = np.where(disp_mwyfart_a<-1, np.NaN, disp_mwyfart_a)




mpl.rcParams.update({'font.size': 12})
mpl.rcParams['figure.dpi'] = 300

bblons = [-9.2, 3.7]
bblats = [38.5, 46]

# initialise data
widths = [1, 1]
heights = [1, 1]
fig = plt.figure(figsize=(14, 8))
spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=widths, height_ratios=heights)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[1, 0])
ax4 = fig.add_subplot(spec[1, 1])

fig.subplots_adjust(hspace=.05, wspace=.05)

m1 = Basemap(projection='laea', llcrnrlon=bblons[0], llcrnrlat=bblats[0], urcrnrlon=bblons[1], urcrnrlat=bblats[1],
            resolution='l', epsg=3035, ax=ax1)
m2 = Basemap(projection='laea', llcrnrlon=bblons[0], llcrnrlat=bblats[0], urcrnrlon=bblons[1], urcrnrlat=bblats[1],
            resolution='l', epsg=3035, ax=ax2)
m3 = Basemap(projection='laea', llcrnrlon=bblons[0], llcrnrlat=bblats[0], urcrnrlon=bblons[1], urcrnrlat=bblats[1],
            resolution='l', epsg=3035, ax=ax3)
m4 = Basemap(projection='laea', llcrnrlon=bblons[0], llcrnrlat=bblats[0], urcrnrlon=bblons[1], urcrnrlat=bblats[1],
            resolution='l', epsg=3035, ax=ax4)


# Get the target projection from the basemap object
m_proj = osr.SpatialReference()
m_proj.ImportFromProj4(m1.proj4string)

# Convert from source projection to basemap projection
xx, yy = convertXY(disp_xy, disp_proj, m_proj)
xx2, yy2 = convertXY(fart_xy, fart_proj, m_proj)

# color
n_unq, vmin, vmax, loc_ticks = cbar_ticks(disp_a)
colmap = disc_color('jet', n_unq-2, False)
pink = mpl.colors.to_rgba('lightpink')
grey = np.array([0.5, 0.5, 0.5, 1])
colmap = np.vstack([grey, colmap])
colmap = np.vstack([colmap, pink])
colmap = mpl.colors.ListedColormap(colmap)

# plot coastlines, meridians etc.
parallels = np.arange(-45, 90, 4)
meridians = np.arange(0, 360, 4)

for m in [m1, m2, m3, m4]:
    m.drawcoastlines(linewidth=1)
    m.drawcountries(linewidth=1)
    m.fillcontinents(color='lightgrey', zorder=0)
    m.drawparallels(parallels)
    m.drawmeridians(meridians)

m2.drawparallels(parallels, fontsize=12, labels = [0, 1, 0, 0])
m4.drawparallels(parallels, fontsize=12, labels = [0, 1, 0, 0])
m3.drawmeridians(meridians, fontsize=12, labels = [0, 0, 0, 1])
m4.drawmeridians(meridians, fontsize=12, labels = [0, 0, 0, 1])



# colmap2
bl = mpl.colors.to_rgba('blue')
yl = mpl.colors.to_rgba('red')
og = mpl.colors.to_rgba('purple')
colmap2 = np.vstack([bl, yl, og])
colmap2 = mpl.colors.ListedColormap(colmap2)


# plot arrays
im1 = m1.pcolormesh(xx2, yy2, disp_mask.T, cmap=colmap2, vmin=1, vmax=3)
im2 = m2.pcolormesh(xx, yy, disp_a.T, cmap=colmap, vmin=vmin, vmax=vmax)
im3 = m3.pcolormesh(xx, yy, disp_fart_a.T, cmap=colmap, vmin=vmin, vmax=vmax)
im4 = m4.pcolormesh(xx, yy, disp_mwyfart_a.T, cmap=colmap, vmin=vmin, vmax=vmax)

# titles
from matplotlib.font_manager import FontProperties
font0 = FontProperties()
font = font0.copy()
font.set_style('normal')
font.set_weight('bold')
font.set_size(18)

# Figure Sub Names
ax1.text(0.025, 0.95, 'A   Barrier Mask', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax1.transAxes, fontsize=14,
        horizontalalignment='left', verticalalignment='top', fontproperties=font)
ax2.text(0.025, 0.95, 'B   No Barrier Scenario', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax2.transAxes, fontsize=14,
        horizontalalignment='left', verticalalignment='top', fontproperties=font)
ax3.text(0.025, 0.95, 'C   Artificial Surface Scenario', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax3.transAxes, fontsize=14,
        horizontalalignment='left', verticalalignment='top', fontproperties=font)
ax4.text(0.025, 0.95, 'D   Artificial + Motorway Scenario', style='normal',
        bbox={'facecolor': 'white', 'alpha': 0, 'pad': 10},
        transform=ax4.transAxes, fontsize=14,
        horizontalalignment='left', verticalalignment='top', fontproperties=font)

#m1.drawmapscale(-7.6, 40.8, 4, 40, 250, barstyle='fancy', fontsize=14)
m2.drawmapscale(-6.105, 39.6, 4, 40, 500, barstyle='fancy', fontsize=12)
m3.drawmapscale(-6.105, 39.6, 4, 40, 500, barstyle='fancy', fontsize=12)
m4.drawmapscale(-6.105, 39.6, 4, 40, 500, barstyle='fancy', fontsize=12)

labels = ['Motorways', 'Artificial', 'Both']
colors = ['blue', 'red', 'purple']
patches = [
    mpl.patches.Patch(color=color, label=label)
    for label, color in zip(labels, colors)]
ax1.legend(patches, labels, loc='upper right', frameon=True, fontsize=12,
           framealpha=1, edgecolor='black', fancybox=False)

p0 = ax3.get_position().get_points().flatten()
p1 = ax4.get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0], 0.05, p1[2]-p0[0], 0.025])
cbar = fig.colorbar(im2, cax=ax_cbar, orientation='horizontal')
cbar.set_ticks(loc_ticks)
cbar.ax.set_xticklabels(['Decolonized', 'Stable', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'Suitable'])

#plt.show()
plt.savefig('/Users/leonnill/Desktop/Dispersal.png', dpi=fig.dpi, bbox_inches='tight')
plt.close()


