import os
import cartopy
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from osgeo import gdal, osr
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib import colors
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
plt.rcParams.update({'font.size': 12})


# ----------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------
def _axes_to_lonlat(ax, coords):
    """(lon, lat) from axes coordinates."""
    display = ax.transAxes.transform(coords)
    data = ax.transData.inverted().transform(display)
    lonlat = ccrs.PlateCarree().transform_point(*data, ax.projection)

    return lonlat


def _upper_bound(start, direction, distance, dist_func):
    """A point farther than distance from start, in the given direction.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        direction  Nonzero (2, 1)-shaped array, a direction vector.
        distance:  Positive distance to go past.
        dist_func: A two-argument function which returns distance.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    if distance <= 0:
        raise ValueError(f"Minimum distance is not positive: {distance}")

    if np.linalg.norm(direction) == 0:
        raise ValueError("Direction vector must not be zero.")

    # Exponential search until the distance between start and end is
    # greater than the given limit.
    length = 0.1
    end = start + length * direction

    while dist_func(start, end) < distance:
        length *= 2
        end = start + length * direction

    return end


def _distance_along_line(start, end, distance, dist_func, tol):
    """Point at a distance from start on the segment  from start to end.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        end:       Outer bound on point's location.
        distance:  Positive distance to travel.
        dist_func: Two-argument function which returns distance.
        tol:       Relative error in distance to allow.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    initial_distance = dist_func(start, end)
    if initial_distance < distance:
        raise ValueError(f"End is closer to start ({initial_distance}) than "
                         f"given distance ({distance}).")

    if tol <= 0:
        raise ValueError(f"Tolerance is not positive: {tol}")

    # Binary search for a point at the given distance.
    left = start
    right = end

    while not np.isclose(dist_func(start, right), distance, rtol=tol):
        midpoint = (left + right) / 2

        # If midpoint is too close, search in second half.
        if dist_func(start, midpoint) < distance:
            left = midpoint
        # Otherwise the midpoint is too far, so search in first half.
        else:
            right = midpoint

    return right


def _point_along_line(ax, start, distance, angle=0, tol=0.01):
    """Point at a given distance from start at a given angle.

    Args:
        ax:       CartoPy axes.
        start:    Starting point for the line in axes coordinates.
        distance: Positive physical distance to travel.
        angle:    Anti-clockwise angle for the bar, in radians. Default: 0
        tol:      Relative error in distance to allow. Default: 0.01

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    # Direction vector of the line in axes coordinates.
    direction = np.array([np.cos(angle), np.sin(angle)])

    geodesic = cgeo.Geodesic()

    # Physical distance between points.
    def dist_func(a_axes, b_axes):
        a_phys = _axes_to_lonlat(ax, a_axes)
        b_phys = _axes_to_lonlat(ax, b_axes)

        # Geodesic().inverse returns a NumPy MemoryView like [[distance,
        # start azimuth, end azimuth]].
        return geodesic.inverse(a_phys, b_phys).base[0, 0]

    end = _upper_bound(start, direction, distance, dist_func)

    return _distance_along_line(start, end, distance, dist_func, tol)


def scale_bar(ax, location, length, metres_per_unit=1000, unit_name='km',
              tol=0.01, angle=0, color='black', linewidth=3, text_offset=0.005,
              ha='center', va='bottom', plot_kwargs=None, text_kwargs=None,
              **kwargs):
    """Add a scale bar to CartoPy axes.

    For angles between 0 and 90 the text and line may be plotted at
    slightly different angles for unknown reasons. To work around this,
    override the 'rotation' keyword argument with text_kwargs.

    Args:
        ax:              CartoPy axes.
        location:        Position of left-side of bar in axes coordinates.
        length:          Geodesic length of the scale bar.
        metres_per_unit: Number of metres in the given unit. Default: 1000
        unit_name:       Name of the given unit. Default: 'km'
        tol:             Allowed relative error in length of bar. Default: 0.01
        angle:           Anti-clockwise rotation of the bar.
        color:           Color of the bar and text. Default: 'black'
        linewidth:       Same argument as for plot.
        text_offset:     Perpendicular offset for text in axes coordinates.
                         Default: 0.005
        ha:              Horizontal alignment. Default: 'center'
        va:              Vertical alignment. Default: 'bottom'
        **plot_kwargs:   Keyword arguments for plot, overridden by **kwargs.
        **text_kwargs:   Keyword arguments for text, overridden by **kwargs.
        **kwargs:        Keyword arguments for both plot and text.
    """
    # Setup kwargs, update plot_kwargs and text_kwargs.
    if plot_kwargs is None:
        plot_kwargs = {}
    if text_kwargs is None:
        text_kwargs = {}

    plot_kwargs = {'linewidth': linewidth, 'color': color, **plot_kwargs,
                   **kwargs}
    text_kwargs = {'ha': ha, 'va': va, 'rotation': angle, 'color': color,
                   **text_kwargs, **kwargs}

    # Convert all units and types.
    location = np.asarray(location)  # For vector addition.
    length_metres = length * metres_per_unit
    angle_rad = angle * np.pi / 180

    # End-point of bar.
    end = _point_along_line(ax, location, length_metres, angle=angle_rad,
                            tol=tol)

    # Coordinates are currently in axes coordinates, so use transAxes to
    # put into data coordinates. *zip(a, b) produces a list of x-coords,
    # then a list of y-coords.
    ax.plot(*zip(location, end), transform=ax.transAxes, **plot_kwargs)

    # Push text away from bar in the perpendicular direction.
    midpoint = (location + end) / 2
    offset = text_offset * np.array([-np.sin(angle_rad), np.cos(angle_rad)])
    text_location = midpoint + offset

    # 'rotation' keyword argument is in text_kwargs.
    ax.text(*text_location, f"{length} {unit_name}", rotation_mode='anchor',
            transform=ax.transAxes, **text_kwargs)


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


# ======================================================================================================================
# FIGURE 1: Binary SDMs
# ======================================================================================================================
# input
ds1 = gdal.Open('Data/SDMs/ursusarctos_SDM_ensemble_mn_md_wmn_cav_sdp.tif')
img1 = ds1.ReadAsArray()
img1 = img1[2,:,:]
img1 = np.where(img1 > 0.46, 1, np.nan)  # lynx = 0.28 bear = 0.46 (weighted mean ensemble)

ds2 = gdal.Open('Data/SDMs/lynxpardinus_SDM_ensemble_mn_md_wmn_cav_sdp.tif')
img2 = ds2.ReadAsArray()
img2 = img2[2,:,:]
img2 = np.where(img2 > 0.28, 1, np.nan)  # lynx = 0.28 bear = 0.46 (weighted mean ensemble)

ds_mask = gdal.Open('Data/Mask/IBERIA_MASK_10km.tif')
mask = ds_mask.ReadAsArray()

# redlist
ds_redlist1 = gdal.Open(r"Data\RedList\ursusarctos_rasterized_10km.tif")
redlist1 = ds_redlist1.ReadAsArray()
redlist1 = np.where(redlist1 != 1, np.nan, 1)
ds_redlist2 = gdal.Open(r"Data\RedList\lynxpardinus_rasterized_10km.tif")
redlist2 = ds_redlist2.ReadAsArray()
redlist2 = np.where(redlist2 != 1, np.nan, 1)

# hillshade
ds_z = gdal.Open(r"C:\Users\Leon\Google Drive\01_MSc_GCG\MSc6_GCIB\Project\Data\DEM\SRTM_V3_90m_IBERIA_250m.tif")
z_data = ds_z.ReadAsArray()
#z_data = np.where(z_data == 0, np.nan, z_data)

from matplotlib.colors import LightSource

# Generate the hillshaded intensity
ls = LightSource()
intensity = ls.hillshade(z_data, vert_exag=1)
cmap = plt.cm.gist_gray
rgb = ls.shade(z_data, cmap=cmap, blend_mode='overlay',
                       vert_exag=1)


# get extent, resolution and projection
xmin, xmax, ymin, ymax = ext_min_max(ds_mask)
extent = (xmin, xmax, ymin, ymax)
xres, yres = ds_mask.GetGeoTransform()[1], ds_mask.GetGeoTransform()[5]
ds_proj = osr.SpatialReference(wkt=ds_mask.GetProjection())

# create meshgrid of coordinates
msh_xy = np.mgrid[xmin:xmax:xres, ymax:ymin:yres]

# extent and background
ext = (2590000.0, 3800000.0, 1520000.0, 2520000.0)

fig = plt.figure(figsize=(10, 6))
axes_class = (GeoAxes, dict(map_projection=ccrs.epsg(3035)))
grid = AxesGrid(fig, 111, axes_class=axes_class,
                      nrows_ncols = (1, 2),
                      axes_pad = 0.2,
                      label_mode = ''
                      )
fig.subplots_adjust(hspace=.1, wspace=.1)
ax1 = grid[0]
ax2 = grid[1]

ax1.set_extent(ext, crs=ccrs.epsg(3035))
ax2.set_extent(ext, crs=ccrs.epsg(3035))

ax1.imshow(rgb, extent=ext_min_max(ds_z), origin='upper', alpha=0.75)
ax2.imshow(rgb, extent=ext_min_max(ds_z), origin='upper', alpha=0.75)

ax1.add_feature(cartopy.feature.OCEAN, facecolor='white', zorder=1)
ax1.add_feature(cartopy.feature.COASTLINE, lw=1, linestyle='-', zorder=2)
ax1.add_feature(cartopy.feature.BORDERS, lw=1, linestyle='-', zorder=2)
ax2.add_feature(cartopy.feature.OCEAN, facecolor='white', zorder=1)
ax2.add_feature(cartopy.feature.COASTLINE, lw=1, linestyle='-', zorder=2)
ax2.add_feature(cartopy.feature.BORDERS, lw=1, linestyle='-', zorder=2)

gl1 = ax1.gridlines(draw_labels=True, linewidth=0.5, alpha=0.75, linestyle='-', color='grey', zorder=5)
gl2 = ax2.gridlines(draw_labels=True, linewidth=0.5, alpha=0.75, linestyle='-', color='grey', zorder=5)

gl1.bottom_labels = False
gl1.right_labels = False
gl1.left_labels = True
gl1.top_labels = True
gl2.bottom_labels = False
gl2.right_labels = False
gl2.left_labels = False
gl2.top_labels = True

gl1.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
gl1.ylocator = mticker.FixedLocator([36, 38, 40, 42, 44, 46])
gl1.xlabel_style = {'size': 10, 'color': 'black'}
gl1.ylabel_style = {'size': 10, 'color': 'black'}
gl2.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
gl2.ylocator = mticker.FixedLocator([36, 38, 40, 42, 44, 46])
gl2.xlabel_style = {'size': 10, 'color': 'black'}
gl2.ylabel_style = {'size': 10, 'color': 'black'}

# raster imgs
ax1.imshow(img1, cmap=mpl.colors.ListedColormap(mpl.colors.to_rgba('lightpink')),
          alpha=0.7, extent=ext_min_max(ds1), origin='upper', zorder=4)
ax1.imshow(redlist1, cmap=mpl.colors.ListedColormap(mpl.colors.to_rgba('red')),
          alpha=0.7, extent=ext_min_max(ds_redlist1), origin='upper', zorder=4)
ax2.imshow(img2, cmap=mpl.colors.ListedColormap(mpl.colors.to_rgba('lightpink')),
          alpha=0.7, extent=ext_min_max(ds2), origin='upper', zorder=4)
ax2.imshow(redlist2, cmap=mpl.colors.ListedColormap(mpl.colors.to_rgba('red')),
          alpha=0.7, extent=ext_min_max(ds_redlist2), origin='upper', zorder=4)

# legend
labels = ['Suitable', 'Occupied']
colors = ['lightpink', 'red']
patches = [
    mpl.patches.Patch(color=color, label=label)
    for label, color in zip(labels, colors)]
ax2.legend(patches, labels, loc='lower right', frameon=True, fontsize=13,
           framealpha=1, edgecolor='black', fancybox=False)

scale_bar(ax1, (0.05, 0.05), 250, zorder=6)
scale_bar(ax2, (0.05, 0.05), 250, zorder=6)

props = dict(boxstyle='square', facecolor='white', alpha=1)
# place a text box in upper left in axes coords
ax1.text(0.95, 0.95, "Ursus arctos", transform=ax1.transAxes, fontsize=12,
        verticalalignment='top', horizontalalignment='right', bbox=props, zorder=6)
ax2.text(0.95, 0.95, "Lynx pardinus", transform=ax2.transAxes, fontsize=12,
        verticalalignment='top', horizontalalignment='right', bbox=props, zorder=6)
ax1.text(0.025, 0.975, "A", transform=ax1.transAxes, fontsize=16,
                verticalalignment='top', horizontalalignment='left',
                fontweight='bold', zorder=6)
ax2.text(0.025, 0.975, "B", transform=ax2.transAxes, fontsize=16,
                verticalalignment='top', horizontalalignment='left',
                fontweight='bold', zorder=6)

plt.savefig("Plots/fig2_binary.pdf", dpi=300, bbox_inches='tight')
plt.close()


# ======================================================================================================================
# FIGURE 2: Continuous SDMs (Ensemble)
# ======================================================================================================================
# input
ds1 = gdal.Open('Data/SDMs/ursusarctos_SDM_ensemble_mn_md_wmn_cav_sdp.tif')
img1 = ds1.ReadAsArray()
img1 = img1[2,:,:]
img1 = np.where(img1 < 0, np.nan, img1)

ds2 = gdal.Open('Data/SDMs/lynxpardinus_SDM_ensemble_mn_md_wmn_cav_sdp.tif')
img2 = ds2.ReadAsArray()
img2 = img2[2,:,:]
img2 = np.where(img2 < 0, np.nan, img2)

ds_mask = gdal.Open('Data/Mask/IBERIA_MASK_10km.tif')
mask = ds_mask.ReadAsArray()

# extent and background
ext = (2590000.0, 3800000.0, 1520000.0, 2520000.0)

#fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.epsg(3035)}, figsize=(14, 10))
fig = plt.figure(figsize=(10, 6))
axes_class = (GeoAxes, dict(map_projection=ccrs.epsg(3035)))
grid = AxesGrid(fig, 111, axes_class=axes_class,
                      nrows_ncols = (1, 2),
                      axes_pad = 0.2,
                      label_mode = '',
                      cbar_location="right",
                      cbar_mode="single",
                      cbar_size="5%",
                      cbar_pad=0.2
                      )
fig.subplots_adjust(hspace=.1, wspace=.1)
ax1 = grid[0]
ax2 = grid[1]

#ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.epsg(3035))
#ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.epsg(3035))

ax1.set_extent(ext, crs=ccrs.epsg(3035))
ax2.set_extent(ext, crs=ccrs.epsg(3035))

ax1.add_feature(cartopy.feature.OCEAN, facecolor='white', zorder=1)
ax1.add_feature(cartopy.feature.COASTLINE, lw=1, linestyle='-', zorder=2)
ax1.add_feature(cartopy.feature.BORDERS, lw=1, linestyle='-', zorder=2)

ax2.add_feature(cartopy.feature.OCEAN, facecolor='white', zorder=1)
ax2.add_feature(cartopy.feature.COASTLINE, lw=1, linestyle='-', zorder=2)
ax2.add_feature(cartopy.feature.BORDERS, lw=1, linestyle='-', zorder=2)

gl1 = ax1.gridlines(draw_labels=True, linewidth=0.5, alpha=0.75, linestyle='-', color='grey', zorder=4)
gl2 = ax2.gridlines(draw_labels=True, linewidth=0.5, alpha=0.75, linestyle='-', color='grey', zorder=4)

gl1.bottom_labels = False
gl1.right_labels = False
gl1.left_labels = True
gl1.top_labels = True
gl2.bottom_labels = False
gl2.right_labels = False
gl2.left_labels = False
gl2.top_labels = True

gl1.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
gl1.ylocator = mticker.FixedLocator([36, 38, 40, 42, 44, 46])
gl1.xlabel_style = {'size': 10, 'color': 'black'}
gl1.ylabel_style = {'size': 10, 'color': 'black'}
gl2.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
gl2.ylocator = mticker.FixedLocator([36, 38, 40, 42, 44, 46])
gl2.xlabel_style = {'size': 10, 'color': 'black'}
gl2.ylabel_style = {'size': 10, 'color': 'black'}

# raster imgs
im1 = ax1.imshow(img1, cmap='jet', vmin=0, vmax=1, extent=ext_min_max(ds1), origin='upper', zorder=0)
im2 = ax2.imshow(img2, cmap='jet', vmin=0, vmax=1, extent=ext_min_max(ds2), origin='upper', zorder=0)

props = dict(boxstyle='square', facecolor='white', alpha=0.75)
# place a text box in upper left in axes coords
ax1.text(0.025, 0.975, "A", transform=ax1.transAxes, fontsize=16,
                verticalalignment='top', horizontalalignment='left',
                fontweight='bold')
ax2.text(0.025, 0.975, "B", transform=ax2.transAxes, fontsize=16,
                verticalalignment='top', horizontalalignment='left',
                fontweight='bold')

scale_bar(ax1, (0.05, 0.05), 250)
scale_bar(ax2, (0.05, 0.05), 250)

# legend
grid.cbar_axes[0].colorbar(im2)
#fig.subplots_adjust(right=0.9)
#cbar_ax = fig.add_axes([0.91, 0, 0.05, 1])  # [left, bottom, width, height]
#fig.colorbar(im2, ax=ax2)

plt.tight_layout(pad=0.05)

plt.savefig("Plots/fig_sdms.pdf", dpi=300, bbox_inches='tight')
plt.close()



# ======================================================================================================================
# FIGURE 2: Continuous SDMs (Single)
# ======================================================================================================================
# input
ds = gdal.Open('Data/SDMs/lynxpardinus_SDM_probability_models.tif')
img = ds.ReadAsArray()
algo_names = ['Generalized Linear Models', 'Generalized Additive Model', 'Bioclim', 'Boosted Regression Trees ',
              'Gaussian Process Regression']
subfig_names = ['A', 'B', 'C', 'D', 'E']
img = np.where(img < 0, np.nan, img)

ds_mask = gdal.Open('Data/Mask/IBERIA_MASK_10km.tif')
mask = ds_mask.ReadAsArray()

# extent and background
ext = (2590000.0, 3800000.0, 1520000.0, 2520000.0)

fig, axes = plt.subplots(3, 2, subplot_kw={'projection': ccrs.epsg(3035)}, figsize=(10, 12))
fig.tight_layout(pad=2.5)
fig.subplots_adjust(hspace=.1, wspace=.1)
for i, ax in enumerate(axes.flatten()):
    if i != 5:
        ax.set_extent(ext, crs=ccrs.epsg(3035))
        ax.add_feature(cartopy.feature.OCEAN, facecolor='white', zorder=1)
        ax.add_feature(cartopy.feature.COASTLINE, lw=1, linestyle='-', zorder=2)
        ax.add_feature(cartopy.feature.BORDERS, lw=1, linestyle='-', zorder=2)

        gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.8, linestyle='-', color='grey', zorder=5)
        gl.bottom_labels = False
        gl.right_labels = False
        #gl.left_labels = False
        gl.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
        gl.ylocator = mticker.FixedLocator([36, 38, 40, 42, 44, 46])
        gl.xlabel_style = {'size': 10, 'color': 'black'}
        gl.ylabel_style = {'size': 10, 'color': 'black'}

        # raster imgs
        im = ax.imshow(img[i], cmap='jet', vmin=0, vmax=1, extent=ext_min_max(ds), origin='upper', zorder=0)
        scale_bar(ax, (0.05, 0.05), 250)

        props = dict(boxstyle='square', facecolor='white', alpha=0.75)
        # place a text box in upper left in axes coords
        ax.text(0.95, 0.95, algo_names[i], transform=ax.transAxes, fontsize=12,
                verticalalignment='top', horizontalalignment='right', bbox=props)
        ax.text(0.025, 0.975, subfig_names[i], transform=ax.transAxes, fontsize=16,
                verticalalignment='top', horizontalalignment='left',
                fontweight='bold')

    if i == 5:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
        ax.set_frame_on(False)

plt.colorbar(im, ax=[ax_cb], location='left')

plt.savefig("Plots/fig_lynxpardinus_single_sdms.pdf", dpi=300)
plt.close()


# ======================================================================================================================
# FIGURE: Dispersal
# ======================================================================================================================
# input
disp1_bear = gdal.Open(r"Data\Dispersal\ursusarctos_th0.025\ursusarctos_th0.025.tif").ReadAsArray()
disp2_bear = gdal.Open(r"Data\Dispersal\ursusarctos_th0.05\ursusarctos_th0.05.tif").ReadAsArray()
disp3_bear = gdal.Open(r"Data\Dispersal\ursusarctos_th0.1\ursusarctos_th0.1.tif").ReadAsArray()
disp4_bear = gdal.Open(r"Data\Dispersal\ursusarctos_th0.25\ursusarctos_th0.25.tif").ReadAsArray()

disp1_lynx = gdal.Open(r"Data\Dispersal\lynxpardinus_th0.025\lynxpardinus_th0.025.tif").ReadAsArray()
disp2_lynx = gdal.Open(r"Data\Dispersal\lynxpardinus_th0.05\lynxpardinus_th0.05.tif").ReadAsArray()
disp3_lynx = gdal.Open(r"Data\Dispersal\lynxpardinus_th0.1\lynxpardinus_th0.1.tif").ReadAsArray()
disp4_lynx = gdal.Open(r"Data\Dispersal\lynxpardinus_th0.25\lynxpardinus_th0.25.tif").ReadAsArray()

#l_imgs = [disp1_bear, disp1_lynx, disp2_bear, disp2_lynx, disp3_bear, disp3_lynx, disp4_bear, disp4_lynx]
#l_imgs = [disp1_bear, disp2_bear, disp3_bear, disp4_bear]
l_imgs = [disp1_lynx, disp2_lynx, disp3_lynx, disp4_lynx]
l_imgs = [np.where(x < -1, np.NaN, x) for x in l_imgs]
l_imgs = [np.where(x == 20, 11, x) for x in l_imgs]

subfig_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
#threshold_names = ['> 2.5%', '> 2.5%', '> 5%', '> 5%', '> 10%', '> 10%', '> 25%', '> 25%']
threshold_names = ['> 2.5%', '> 5%', '> 10%', '> 25%']
species_names = ['Ursus arctos', 'Lynx pardinus', 'Ursus arctos', 'Lynx pardinus', 'Ursus arctos', 'Lynx pardinus',
                 'Ursus arctos', 'Lynx pardinus']

# color
def get_colmap(a):
    n_unq, vmin, vmax, loc_ticks = cbar_ticks(a)
    pink = mpl.colors.to_rgba('lightpink')
    decolonized = mpl.colors.to_rgba('blueviolet')
    if np.any(a == -1):
        colmap = disc_color('jet', n_unq - 2, False)
        colmap = np.vstack([decolonized, colmap])
    else:
        colmap = disc_color('jet', n_unq - 1, False)
    colmap = np.vstack([colmap, pink])
    colmap = mpl.colors.ListedColormap(colmap)
    return vmin, vmax, colmap, loc_ticks


ds_z = gdal.Open(r"E:\Meine Ablage\01_MSc_GCG\MSc6_GCIB\Project\Data\DEM\SRTM_V3_90m_IBERIA_250m.tif")
z_data = ds_z.ReadAsArray()
from matplotlib.colors import LightSource
# Generate the hillshaded intensity
ls = LightSource()
intensity = ls.hillshade(z_data, vert_exag=1)
cmap = plt.cm.gist_gray
rgb = ls.shade(z_data, cmap=cmap, blend_mode='overlay',
                       vert_exag=1)

ds_mask = gdal.Open('Data/Mask/IBERIA_MASK_10km.tif')
mask = ds_mask.ReadAsArray()

# extent and background
ext = (2590000.0, 3800000.0, 1520000.0, 2520000.0)

w = 7
p_left = 0.25
p_right = 0.25
w_cbar = 0.125
cbar_aspect = (w - p_left - p_right) / w_cbar


# original version with 4 scenarios
fig, axes = plt.subplots(2, 2, subplot_kw={'projection': ccrs.epsg(3035)}, figsize=(7, 6.5))  # 7, 13
fig.subplots_adjust(hspace=.05, wspace=.05)
for i, ax in enumerate(axes.flatten()):
    ax.set_extent(ext, crs=ccrs.epsg(3035))
    ax.add_feature(cartopy.feature.OCEAN, facecolor='white', zorder=1)
    ax.add_feature(cartopy.feature.COASTLINE, lw=1, linestyle='-', zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, lw=1, linestyle='-', zorder=2)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.8, linestyle='-', color='grey', zorder=5)
    gl.left_labels = False
    gl.bottom_labels = False
    gl.top_labels = False
    gl.right_labels = False
    if (i == 0) | (i == 1):
        gl.top_labels = True
    if (i == 0) | (i == 2) | (i == 4) | (i == 6):
        gl.left_labels = True
    gl.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
    gl.ylocator = mticker.FixedLocator([36, 38, 40, 42, 44, 46])
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}

    # raster imgs
    ax.imshow(rgb, extent=ext_min_max(ds_z), origin='upper', alpha=0.75)
    temp_img = l_imgs[i]
    vmin, vmax, temp_colmap, loc_ticks = get_colmap(temp_img)
    im = ax.imshow(temp_img, cmap=temp_colmap,
                   extent=ext_min_max(gdal.Open(r"Data\Dispersal\ursusarctos_th0.025\ursusarctos_th0.025.tif")),
                   origin='upper', zorder=1, alpha=0.75, vmin=vmin, vmax=vmax)
    if i == 0:
        vmin, vmax, temp_colmap, loc_ticks2 = get_colmap(temp_img)
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), orientation='horizontal',
                            fraction=0.08, pad=0.01, aspect=cbar_aspect)
    scale_bar(ax, (0.05, 0.05), 250, zorder=6)
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    #ax.text(0.95, 0.95, species_names[i], transform=ax.transAxes, fontsize=12,
    #    verticalalignment='top', horizontalalignment='right', bbox=props, zorder=6)
    ax.text(0.95, 0.05, threshold_names[i], transform=ax.transAxes, fontsize=12,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=props, zorder=6)
    ax.text(0.05, 0.95, subfig_names[i], transform=ax.transAxes, fontsize=16,
             verticalalignment='top', horizontalalignment='left',
             fontweight='bold', zorder=6)

'''
p0 = axes.flatten()[6].get_position().get_points().flatten()
p1 = axes.flatten()[7].get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0], 0.05, p1[2]-p0[0], 0.025])
ax_cb = ax_cbar.append_axes('bottom', size="2%", pad=0.5)
cbar = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
'''

#cbar = fig.colorbar(im, ax=axes.ravel().tolist(), orientation='horizontal',
#                    fraction=0.08, pad=0.01, aspect=cbar_aspect)

cbar.set_label('Year of Occupation', size=10)
cbar.set_ticks(loc_ticks2)
cbar.ax.set_xticklabels(['Decolonised', 'Stable', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'Uncolonised'])
cbar.ax.tick_params(labelsize=10)
#plt.colorbar(im, ax=[ax_cb], location='left')

plt.savefig("Plots/fig_dispersal_lynx.pdf", dpi=300, bbox_inches='tight')
plt.close()


# new version with 2 scenarios
#l_imgs = [disp1_bear, disp4_bear]
l_imgs = [disp1_lynx, disp4_lynx]
l_imgs = [np.where(x < -1, np.NaN, x) for x in l_imgs]
l_imgs = [np.where(x == 20, 11, x) for x in l_imgs]

subfig_names = ['A', 'B']
threshold_names = ['> 2.5%', '> 25%']

fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.epsg(3035)}, figsize=(8, 5))  # 7, 13
fig.subplots_adjust(hspace=.05, wspace=.05)
for i, ax in enumerate(axes.flatten()):
    ax.set_extent(ext, crs=ccrs.epsg(3035))
    ax.add_feature(cartopy.feature.OCEAN, facecolor='white', zorder=1)
    ax.add_feature(cartopy.feature.COASTLINE, lw=1, linestyle='-', zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, lw=1, linestyle='-', zorder=2)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.8, linestyle='-', color='grey', zorder=5)
    gl.left_labels = False
    gl.bottom_labels = False
    gl.top_labels = False
    gl.right_labels = False
    if (i == 0) | (i == 1):
        gl.top_labels = True
    if (i == 0) | (i == 2) | (i == 4) | (i == 6):
        gl.left_labels = True
    gl.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
    gl.ylocator = mticker.FixedLocator([36, 38, 40, 42, 44, 46])
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}

    # raster imgs
    ax.imshow(rgb, extent=ext_min_max(ds_z), origin='upper', alpha=0.75)
    temp_img = l_imgs[i]
    vmin, vmax, temp_colmap, loc_ticks = get_colmap(temp_img)
    im = ax.imshow(temp_img, cmap=temp_colmap,
                   extent=ext_min_max(gdal.Open(r"Data\Dispersal\ursusarctos_th0.025\ursusarctos_th0.025.tif")),
                   origin='upper', zorder=1, alpha=0.75, vmin=vmin, vmax=vmax)
    if i == 0:
        vmin, vmax, temp_colmap, loc_ticks2 = get_colmap(temp_img)
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), orientation='horizontal',
                            fraction=0.08, pad=0.01, aspect=cbar_aspect)
    #scale_bar(ax, (0.05, 0.05), 250, zorder=6)
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(0.95, 0.05, threshold_names[i], transform=ax.transAxes, fontsize=14,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=props, zorder=6)
    ax.text(0.05, 0.95, subfig_names[i], transform=ax.transAxes, fontsize=18,
             verticalalignment='top', horizontalalignment='left',
             fontweight='bold', zorder=6)

'''
p0 = axes.flatten()[6].get_position().get_points().flatten()
p1 = axes.flatten()[7].get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0], 0.05, p1[2]-p0[0], 0.025])
ax_cb = ax_cbar.append_axes('bottom', size="2%", pad=0.5)
cbar = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
'''

#cbar = fig.colorbar(im, ax=axes.ravel().tolist(), orientation='horizontal',
#                    fraction=0.08, pad=0.01, aspect=cbar_aspect)

cbar.set_label('Year of recolonisation', size=10)
cbar.set_ticks(loc_ticks2)
cbar.ax.set_xticklabels(['Decol.', 'Stable', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'Uncol.'])
cbar.ax.tick_params(labelsize=10)
#plt.colorbar(im, ax=[ax_cb], location='left')

plt.savefig("Plots/fig_dispersal_lynx_v2.pdf", dpi=300, bbox_inches='tight')
plt.close()




# -----------------
# Priorisation
# -----------------
vague_bear = gdal.Open(r"C:\Users\Leon\PycharmProjects\iberian_conservation\R\GAP_analysis\results\bear\dispersal_vague_noPA_repro.tif").ReadAsArray()
vague_lynx = gdal.Open(r"C:\Users\Leon\PycharmProjects\iberian_conservation\R\GAP_analysis\results\lynx\dispersal_vague_noPA_repro.tif").ReadAsArray()
low_bear = gdal.Open(r"C:\Users\Leon\PycharmProjects\iberian_conservation\R\GAP_analysis\results\bear\low_colonize_repro.tif").ReadAsArray()
low_lynx = gdal.Open(r"C:\Users\Leon\PycharmProjects\iberian_conservation\R\GAP_analysis\results\lynx\low_colonize_repro.tif").ReadAsArray()

l_imgs = [vague_bear, vague_lynx, low_bear, low_lynx]
l_imgs = [np.where(x == 0, np.NaN, x) for x in l_imgs]

subfig_names = ['A', 'B']
threshold_names = ['> 2.5%', '> 25%']

fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.epsg(3035)}, figsize=(8, 5))  # 7, 13
fig.subplots_adjust(hspace=.05, wspace=.05)
for i, ax in enumerate(axes.flatten()):
    ax.set_extent(ext, crs=ccrs.epsg(3035))
    ax.add_feature(cartopy.feature.OCEAN, facecolor='white', zorder=1)
    ax.add_feature(cartopy.feature.COASTLINE, lw=1, linestyle='-', zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, lw=1, linestyle='-', zorder=2)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.8, linestyle='-', color='grey', zorder=5)
    gl.left_labels = False
    gl.bottom_labels = False
    gl.top_labels = False
    gl.right_labels = False
    if (i == 0) | (i == 1):
        gl.top_labels = True
    if (i == 0) | (i == 2) | (i == 4) | (i == 6):
        gl.left_labels = True
    gl.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
    gl.ylocator = mticker.FixedLocator([36, 38, 40, 42, 44, 46])
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}

    # raster imgs
    ax.imshow(rgb, extent=ext_min_max(ds_z), origin='upper', alpha=0.75)

ax = axes.flatten()
im = ax[0].imshow(l_imgs[0], cmap=colors.ListedColormap(['red']),
               extent=ext_min_max(gdal.Open(r"Data\Dispersal\ursusarctos_th0.025\ursusarctos_th0.025.tif")),
               origin='upper', zorder=1, alpha=1)
im2 = ax[0].imshow(l_imgs[1], cmap=colors.ListedColormap(['yellow']),
               extent=ext_min_max(gdal.Open(r"Data\Dispersal\ursusarctos_th0.025\ursusarctos_th0.025.tif")),
               origin='upper', zorder=1, alpha=1)
im3 = ax[1].imshow(l_imgs[2], cmap=colors.ListedColormap(['red']),
               extent=ext_min_max(gdal.Open(r"Data\Dispersal\ursusarctos_th0.025\ursusarctos_th0.025.tif")),
               origin='upper', zorder=1, alpha=1)
im4 = ax[1].imshow(l_imgs[3], cmap=colors.ListedColormap(['yellow']),
               extent=ext_min_max(gdal.Open(r"Data\Dispersal\ursusarctos_th0.025\ursusarctos_th0.025.tif")),
               origin='upper', zorder=1, alpha=1)

import matplotlib.patches as mpatches
patch1 = mpatches.Patch(color='red', label='Ursus arctos')
patch2 = mpatches.Patch(color='yellow', label='Lynx pardinus')
all_handles = (patch1, patch2)
leg = ax[1].legend(handles=all_handles, loc=4, frameon=False, facecolor=None, fontsize=10)

plt.savefig("Plots/fig_prioritisation.pdf", dpi=300, bbox_inches='tight')
plt.close()




# *********************************************************************************************************************
# *********************************************************************************************************************






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

plt.show()
#plt.savefig('/Users/leonnill/Desktop/Dispersal.png', dpi=fig.dpi, bbox_inches='tight')
#plt.close()
