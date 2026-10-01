"""
demo_read_zarr_s3.py

s3://twc-graf-reforecast/2026-05-01.canonical/

This program will read selected GRAF reforecast files saved
in the S3 data store and plot.   It does not read upper-air
data with multiple vertical levels, but rather demo's
the capability with the single-variable data.

This version reads from the corrected "canonical" copy of the
GRAF reforecast (s3://twc-graf-reforecast/2026-05-01.canonical/),
which addresses issues found in the original release, including
temporally mis-ordered data, uncompressed storage, and missing
cases.  See the repository README for details.

You are likely to have to install some python libraries to
perform the reading, including ones to acess from s3.
There is also a "dateutils.py" file with some date processing
routines.  Less common libraries include s3fs, fsspec,
and boto3


Usage:

python demo_read_zarr_s3.py cyyyymmddhh clead cvar

where

cyyyymmddhh is YEAR/MONTH/DAY/HOUR of the
initial condition.

clead is the lead time in hours, e.g., 12, or
12.25 (for 15-min data) or 12.0833 (for 5-min data)

cvar is the variable name,
    t2m, u10, v10, mslp, dewpoint_2m, visibility
    ceiling_agl, cporain, cpoice, cposnow
    precipw, snowh, windgust10m, cape, ptype
    total_cloud_cover, rain_bucket, conv_bucket
    zrain_bucket, snow_bucket, apcp_bucket

Example:
    $python demo_read_zarr_s3.py 2023122012 12 apcp_bucket

Code by Tom Hamill, The Weather Company, 16 June 2025
Updated October 2026 to read from the corrected "canonical" copy
of the reforecast.

"""
from netCDF4 import Dataset
import os, sys
import numpy as np
import zarr
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
from mpl_toolkits.basemap import Basemap
from scipy.spatial import cKDTree
from datetime import datetime, timedelta
from dateutils import splitdate
from string import Formatter
import fsspec
import warnings
from numcodecs.abc import Codec
from numcodecs.compat import ensure_ndarray, ndarray_copy
from numcodecs.registry import register_codec

warnings.filterwarnings("ignore", category=SyntaxWarning)

# --------------------------------------------------------------------

class Normalize(Codec):
    """Simplified version of the scale-offset filter available in HDF5.
    Applies the transformation `(x - mean) * std` to all chunks. Results
    are rounded to the nearest integer but are not packed according to the
    minimum number of bits.

    Parameters
    ----------
    mean : float
        Value to subtract from data.
    std : float
        Value to multiply by data.
    """

    codec_id = 'normalize'

    def __init__(self, mean, std, dtype, astype=None):
        self.mean = mean
        self.std = std
        self.dtype = np.dtype(dtype)
        if astype is None:
            self.astype = self.dtype
        else:
            self.astype = np.dtype(astype)
        if self.dtype == np.dtype(object) or self.astype == np.dtype(object):
            raise ValueError('object arrays are not supported')

    def encode(self, buf):
        # normalise input
        arr = ensure_ndarray(buf).view(self.dtype)

        # flatten to simplify implementation
        arr = arr.reshape(-1, order='A')

        # compute std mean
        enc = (arr - self.mean) / self.std

        # convert dtype
        enc = enc.astype(self.astype, copy=False)

        return enc


    def decode(self, buf, out=None):
        # interpret buffer as numpy array
        enc = ensure_ndarray(buf).view(self.astype)

        # flatten to simplify implementation
        enc = enc.reshape(-1, order='A')

        # decode std mean
        dec = (enc * self.std) + self.mean

        # convert dtype
        dec = dec.astype(self.dtype, copy=False)

        # handle output
        return ndarray_copy(dec, out)


    def get_config(self):
        # override to handle encoding dtypes
        return dict(
            id=self.codec_id,
            std=self.std,
            mean=self.mean,
            dtype=self.dtype.str,
            astype=self.astype.str,
        )


    def __repr__(self):
        r = '%s(std=%s, mean=%s, dtype=%r' % (
            type(self).__name__,
            self.std,
            self.mean,
            self.dtype.str,
        )
        if self.astype != self.dtype:
            r += ', astype=%r' % self.astype.str
        r += ')'
        return r


register_codec(Normalize)

# --------------------------------------------------------------------

def plotting_information(cyyyymmddhh, clead):
    
    """
    define a dictionary for plotting information (plot titles, 
    contour levels, output file name, etc).
    """

    plot_info = {
        'title': {
            't2m': '2-meter temperature forecast, IC = '+cyyyymmddhh+', '+clead+' h',
            'u10': '10-meter u wind, IC = '+cyyyymmddhh+', '+clead+' h',
            'v10': '10-meter v wind, IC = '+cyyyymmddhh+', '+clead+' h',
            'mslp': 'Mean sea-level pressure, IC = '+cyyyymmddhh+', '+clead+' h',
            'dewpoint_2m': '2-meter dew point, IC = '+cyyyymmddhh+', '+clead+' h',
            'visibility': 'Surface visibility, IC = '+cyyyymmddhh+', '+clead+' h',
            'ceiling_agl': 'Cloud ceiling, IC = '+cyyyymmddhh+', '+clead+' h',
            'cporain': 'Conditional probability of rain, IC = '+cyyyymmddhh+', '+clead+' h',
            'cposnow': 'Conditional probability of snow, IC = '+cyyyymmddhh+', '+clead+' h',
            'cpoice': 'Conditional probability of ice, IC = '+cyyyymmddhh+', '+clead+' h',
            'precipw': 'Total column precipitable water, IC = '+cyyyymmddhh+', '+clead+' h',
            'snowh': 'Snow depth, IC = '+cyyyymmddhh+', '+clead+' h',
            'windgust10m': '10-m wind gust, IC = '+cyyyymmddhh+', '+clead+' h',
            'cape': 'CAPE, IC = '+cyyyymmddhh+', '+clead+' h',
            'ptype': 'Predominant precipitation type, IC = '+cyyyymmddhh+', '+clead+' h',
            'total_cloud_cover': 'Total cloud, IC = '+cyyyymmddhh+', '+clead+' h',
            'rain_bucket': '15-min total rainfall, IC = '+cyyyymmddhh+', '+clead+' h',
            'conv_bucket': '15-min convective precipitation, IC = '+cyyyymmddhh+', '+clead+' h',
            'zrain_bucket': '15-min freezing rain accumulation, IC = '+cyyyymmddhh+', '+clead+' h',
            'snow_bucket': '15-min snow liquid water, IC = '+cyyyymmddhh+', '+clead+' h',
            'apcp_bucket': '15-min total precipitation, IC = '+cyyyymmddhh+', '+clead+' h'},
        'plot_file': {
            't2m': 'GRAFR_t2m_'+cyyyymmddhh+'_'+clead+'h.png',
            'u10': 'GRAFR_u10m_'+cyyyymmddhh+'_'+clead+'h.png',
            'v10': 'GRAFR_uv10m_'+cyyyymmddhh+'_'+clead+'h.png',
            'mslp': 'GRAFR_MSLP_'+cyyyymmddhh+'_'+clead+'h.png',
            'dewpoint_2m': 'GRAFR_td2m_'+cyyyymmddhh+'_'+clead+'h.png',
            'visibility': 'GRAFR_visibility_'+cyyyymmddhh+'_'+clead+'h.png',
            'ceiling_agl': 'GRAFR_ceiling_AGL_'+cyyyymmddhh+'_'+clead+'h.png',
            'cporain': 'GRAFR_cporain_'+cyyyymmddhh+'_'+clead+'h.png',
            'cposnow': 'GRAFR_cposnow_'+cyyyymmddhh+'_'+clead+'h.png',
            'cpoice': 'GRAFR_cpoice_'+cyyyymmddhh+'_'+clead+'h.png',
            'precipw': 'GRAFR_precipw_'+cyyyymmddhh+'_'+clead+'h.png',
            'snowh': 'GRAFR_snowdepth_'+cyyyymmddhh+'_'+clead+'h.png',
            'windgust10m': 'GRAFR_10m_windgust_'+cyyyymmddhh+'_'+clead+'h.png',
            'cape': 'GRAFR_CAPE_'+cyyyymmddhh+'_'+clead+'h.png',
            'ptype': 'GRAFR_predom_precip_type_'+cyyyymmddhh+'_'+clead+'h.png',
            'total_cloud_cover': 'GRAFR_total_cloud_'+cyyyymmddhh+'_'+clead+'h.png',
            'rain_bucket': 'GRAFR_rain_bucket_'+cyyyymmddhh+'_'+clead+'h.png',
            'conv_bucket': 'GRAFR_conv_bucket_'+cyyyymmddhh+'_'+clead+'h.png',
            'zrain_bucket': 'GRAFR_zrain_bucket_'+cyyyymmddhh+'_'+clead+'h.png',
            'snow_bucket': 'GRAFR_snow_bucket_'+cyyyymmddhh+'_'+clead+'h.png',
            'apcp_bucket': 'GRAFR_apcp_bucket_'+cyyyymmddhh+'_'+clead+'h.png'},
        'clevs': {
            't2m':[-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45],
            'u10':[-25,-20,-15,-10,-5,0,5,10,15,20,25],
            'v10':[-25,-20,-15,-10,-5,0,5,10,15,20,25],
            'mslp':[980,984,988,992,996,1000,1004,1008,1012,1016,1020,1024,1028,1032],
            'dewpoint_2m': [-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45],
            'visibility': [0,100,500,1000,1500,2000,3000,5000,7500,10000,20000,30000],
            'ceiling_agl': [0,100,500,1000,1500,2000,3000,5000,7500,10000,20000,30000],
            'cporain': [0,5,10,20,30,40,50,60,70,80,90,95,97,100],
            'cposnow': [0,5,10,20,30,40,50,60,70,80,90,95,97,100],
            'cpoice': [0,5,10,20,30,40,50,60,70,80,90,95,97,100],
            'precipw': [0,1,2,3,5,7,10,12,15,20,25,30,40,60,80],
            'snowh': [0.0,0.01,0.02,0.03,0.05,0.07,0.1,0.15,0.2,0.3,0.4,0.5,1.0,2.0],
            'windgust10m': [0,1,2,3,5,6,8,10,12,15,20,25,30,40],
            'cape': [0,50,100,150,200,300,500,700,1000,1250,1500,2000,2500,3000,4000],
            'ptype': [0,1,2,3,4],
            'total_cloud_cover': [0,5,10,20,30,40,50,60,70,80,90,95,97,100],
            'rain_bucket': [0,0.1,0.2,0.3,0.5,0.7,1,2,3,5,7,10],
            'conv_bucket': [0,0.1,0.2,0.3,0.5,0.7,1,2,3,5,7,10],
            'zrain_bucket': [0,0.1,0.2,0.3,0.5,0.7,1,2,3,5,7,10],
            'snow_bucket': [0,0.1,0.2,0.3,0.5,0.7,1,2,3,5,7,10],
            'apcp_bucket': [0,0.1,0.2,0.3,0.5,0.7,1,2,3,5,7,10] },
        'scale_factor': {
            't2m':1.0,
            'u10':1.0,
            'v10':1.0,
            'mslp':0.01,
            'dewpoint_2m': 1.0,
            'visibility': 1.0,
            'ceiling_agl': 1.0,
            'cporain': 1.0,
            'cposnow': 1.0,
            'cpoice': 1.0,
            'precipw': 1.0,
            'snowh': 1.0,
            'windgust10m': 1.0,
            'cape': 1.0,
            'ptype': 1.0,
            'total_cloud_cover': 1.0,
            'rain_bucket': 1.0,
            'conv_bucket': 1.0,
            'zrain_bucket': 1.0,
            'snow_bucket': 1.0,
            'apcp_bucket': 1.0},
        'offset_factor': {
            't2m':273.15,
            'u10':0.0,
            'v10':0.0,
            'mslp':0.0,
            'dewpoint_2m': 273.15,
            'visibility': 0.0,
            'ceiling_agl': 0.0,
            'cporain': 0.0,
            'cposnow': 0.0,
            'cpoice': 0.0,
            'precipw': 0.0,
            'snowh': 0.0,
            'windgust10m': 0.0,
            'cape': 0.0,
            'ptype': 0.0,
            'total_cloud_cover': 0.0,
            'rain_bucket': 0.0,
            'conv_bucket': 0.0,
            'zrain_bucket': 0.0,
            'snow_bucket': 0.0,
            'apcp_bucket': 0.0},
        'clabel': {
            't2m':'2-meter temperature (deg C)',
            'u10':'10-m u-wind component (m/s)',
            'v10':'10-m v-wind component (m/s)',
            'mslp':'Mean sea-level pressure (hPa)',
            'dewpoint_2m': '2-meter dew point (deg C)',
            'visibility': 'Visibility (m)',
            'ceiling_agl': 'Ceiling (m)',
            'cporain': 'Conditional probability of rain (%)',
            'cposnow': 'Conditional probability of snow (%)',
            'cpoice': 'Conditional probability of ice (%)',
            'precipw': 'Total-column precipitable water (kg/m2)',
            'snowh': 'Snow depth (m)',
            'windgust10m': '10-m wind gust (m/s)',
            'cape': 'CAPE (J/kg)',
            'ptype': 'Precipitation type',
            'total_cloud_cover': 'Total cloud cover (%)',
            'rain_bucket': 'Accumulated non-convective rain (mm)',
            'conv_bucket': 'Accumulated convective rain (mm)',
            'zrain_bucket': 'Accumulated freezing rain (mm)',
            'snow_bucket': 'Accumulated snow liquid equivalent (mm)',
            'apcp_bucket': 'Accumulated total precipitation (mm)'}
        }

    for key in plot_info['title']:
        plot_info['title'][key] = 'GRAF reforecast ' + plot_info['title'][key]

    return plot_info

# --------------------------------------------------------------------

def build_or_load_conus_grid(m, lons, lats, cache_file='GRAF_CONUS_grid_cache.npz',
    nx=2000, ny=1500, margin_m=100000.):

    """
    The MPAS mesh underlying GRAF is global (nCells ~4.8M), but only a
    fraction of those points are needed to cover the CONUS plotting
    domain.  This precomputes, once, a nearest-neighbor mapping from a
    regular grid spanning the *entire* Basemap plot rectangle (m.xmin to
    m.xmax, m.ymin to m.ymax) back onto the mesh points that could serve
    as a nearest neighbor anywhere in it, and caches the mapping to disk
    since the mesh geometry and map projection never change from run to
    run.  Reading data for a new variable/time/case only requires
    gathering onto this fixed grid, not recomputing the mapping.

    Filtering candidate mesh points by a lon/lat box is unsafe here: the
    Lambert Conformal Conic domain is curved/trapezoidal in lon/lat space,
    so a simple box can miss real points needed near the map's corners
    and edges (Pacific coast, Canada, Caribbean), leaving those regions
    without a valid nearest neighbor.  Instead, project every point and
    keep those that land in projected (x, y) space within the map
    rectangle plus a flat margin.
    """

    if os.path.exists(cache_file):
        npz = np.load(cache_file)
        return (npz['inbox_idx'], npz['xg'], npz['yg'],
            npz['nn_idx'], npz['nn_dist'])

    print ('Building CONUS regridding cache (one-time cost): ', cache_file)

    x_all, y_all = m(lons, lats)
    x_all = np.asarray(x_all); y_all = np.asarray(y_all)
    inbox = (x_all >= m.xmin-margin_m) & (x_all <= m.xmax+margin_m) & \
        (y_all >= m.ymin-margin_m) & (y_all <= m.ymax+margin_m)
    inbox_idx = np.nonzero(inbox)[0]

    tree = cKDTree(np.column_stack([x_all[inbox_idx], y_all[inbox_idx]]))

    # ---- Grid spans the full Basemap plot rectangle, not just the
    #      bounding box of the candidate points, so pcolormesh fills the
    #      whole domain out to its true edges/corners.

    xg = np.linspace(m.xmin, m.xmax, nx)
    yg = np.linspace(m.ymin, m.ymax, ny)
    XG, YG = np.meshgrid(xg, yg)
    nn_dist, nn_idx = tree.query(np.column_stack([XG.ravel(), YG.ravel()]),
        k=1, workers=-1)

    np.savez(cache_file, inbox_idx=inbox_idx, xg=xg, yg=yg,
        nn_idx=nn_idx, nn_dist=nn_dist)
    return inbox_idx, xg, yg, nn_idx, nn_dist

# --------------------------------------------------------------------

def build_precip_cmap_norm(clevs):

    """
    Reproduce, via pcolormesh on the pre-gridded raster, the same color
    mapping that m.contourf(..., colors=colorst, extend='both') produced
    on the unstructured mesh: matplotlib actually uses the first
    (len(clevs)-1) colors of the oversized colorst palette as the
    interior colors, and re-uses the first/last of those as the
    "under"/"over" colors (rather than reserving two extra colors).
    """

    colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
        '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
        '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

    n_intervals = len(clevs) - 1
    colors_used = colorst[:n_intervals]
    return from_levels_and_colors(
        clevs, [colors_used[0]] + colors_used + [colors_used[-1]], extend='both')

# --------------------------------------------------------------------

def setup_animation_plot(m, XG, YG, clevs, cvar, clabel):

    """
    One-time setup for an animation: build the persistent figure/axes and
    apply the map projection/aspect/tick setup exactly once.  Because this
    axes object is created once and reused for every frame (rather than a
    fresh one per frame, as plot_data() does), it's immune to the
    axes-identity hash collision described in plot_data()'s comments --
    there's only ever one axes, so Basemap's internal per-axes cache
    always applies to it correctly.  It also means the static map
    background (coastlines/borders) and colorbar, which never change
    frame-to-frame, only ever get drawn once, which is faster too.

    Returns a dict of state that render_animation_frame() updates and
    saves each frame.
    """

    cmap, norm = build_precip_cmap_norm(clevs)

    fig = plt.figure(figsize=(9,9.))
    axloc = [0.02,0.08,0.96,0.84]
    ax1 = fig.add_axes(axloc)

    placeholder = np.zeros(XG.shape)
    quadmesh = ax1.pcolormesh(XG, YG, placeholder, cmap=cmap, norm=norm, shading='auto')

    ax1.set_xlim(m.xmin, m.xmax)
    ax1.set_ylim(m.ymin, m.ymax)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_autoscale_on(False)
    ax1.set_aspect('equal', anchor='C')

    m.drawcoastlines(linewidth=0.9,color='Gray',ax=ax1)
    m.drawcountries(linewidth=0.9,color='Gray',ax=ax1)
    m.drawstates(linewidth=0.7,color='Gray',ax=ax1)

    cax = fig.add_axes([0.02,0.055,0.96,0.02])
    cb = plt.colorbar(quadmesh,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=10)
    cb.set_label(clabel, fontsize=12)

    if cvar == 'ptype':
        cb.ax.set_xticklabels(['                     Rain',\
            '                        Ice/Graupel',\
            '                Sleet / Ice Pellets',\
            '                        Snow',' '])  # horizontal colorbar

    return dict(fig=fig, ax1=ax1, quadmesh=quadmesh)

# --------------------------------------------------------------------

def render_animation_frame(anim_state, data_array, title, plot_file, dpi=150):

    """ Update the persistent plot from setup_animation_plot() with one
    frame's data/title and save it. """

    anim_state['quadmesh'].set_array(np.asarray(data_array).ravel())
    anim_state['ax1'].set_title(title, fontsize=15, color='Black')
    anim_state['fig'].savefig(plot_file, dpi=dpi)
    print ('saving plot to file = ', plot_file)

# --------------------------------------------------------------------

def plot_data(data_array, XG, YG, m, title, plot_file, clevs, cvar, clabel, dpi=400, tight=True):

    #  --- plot hourly GRAF precipitation amount in CONUS/North American
    #      domain.  Save to png file.

    cmap, norm = build_precip_cmap_norm(clevs)

    fig = plt.figure(figsize=(9,9.))
    axloc = [0.02,0.08,0.96,0.84]
    ax1 = fig.add_axes(axloc)
    ax1.set_title(title, fontsize=15,color='Black')
    CS2 = ax1.pcolormesh(XG, YG, data_array, cmap=cmap, norm=norm, shading='auto')

    # ---- Force the view to the full Basemap plot rectangle so the data
    #      always fills the whole domain, regardless of autoscaling.

    ax1.set_xlim(m.xmin, m.xmax)
    ax1.set_ylim(m.ymin, m.ymax)

    # ---- The x/y coordinates here are map-projection meters, not
    #      anything meaningful to show a viewer.  Without this, matplotlib
    #      sometimes renders those tick labels (and a "1e6"-style
    #      scientific-notation offset label above the corner) if its
    #      formatter decides the label text is short enough to fit in the
    #      narrow figure margin -- appearing/disappearing between frames.

    ax1.set_xticks([])
    ax1.set_yticks([])

    # ---- Basemap's drawcoastlines/drawcountries/drawstates normally set
    #      the correct equal-aspect ratio for this axes themselves (via
    #      Basemap.set_axes_limits), but only the first time a given axes
    #      object is used -- it remembers which axes it's already set up
    #      in a set keyed by hash(ax), which for a plain object defaults
    #      to something derived from its memory address.  Across hundreds
    #      of frames, each with its own freshly created (and then closed)
    #      Figure/Axes, Python's allocator reuses freed addresses, so a
    #      brand-new axes can collide with a stale hash left behind by an
    #      long-since-closed one.  When that happens, Basemap thinks this
    #      axes was already set up and skips re-applying the aspect ratio,
    #      leaving that one frame vertically stretched/squashed.  Setting
    #      it explicitly here every time makes the plot immune to that.

    ax1.set_autoscale_on(False)
    ax1.set_aspect('equal', anchor='C')

    m.drawcoastlines(linewidth=0.9,color='Gray',ax=ax1)
    m.drawcountries(linewidth=0.9,color='Gray',ax=ax1)
    m.drawstates(linewidth=0.7,color='Gray',ax=ax1)
    #m.drawcounties(linewidth=0.1,color='LightGray')

    cax = fig.add_axes([0.02,0.055,0.96,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=10)
    cb.set_label(clabel, fontsize=12)

    if cvar == 'ptype':
        cb.ax.set_xticklabels(['                     Rain',\
            '                        Ice/Graupel',\
            '                Sleet / Ice Pellets',\
            '                        Snow',' '])  # horizontal colorbar

    # ---- set plot title

    # ---- bbox_inches='tight' crops to actual rendered content, whose
    #      extent can shift by a pixel or two between frames (e.g. tick
    #      label widths).  Animation frames need identical pixel
    #      dimensions to stack into a GIF, so skip the crop there.

    fig.savefig(plot_file, dpi=dpi, bbox_inches='tight' if tight else None)
    plt.close(fig)
    print ('saving plot to file = ',plot_file)
    print ('Done!')
    istat = 0
    return istat

# ------------------------------------------------------------------------

valid_variables = ['t2m', 'u10', 'v10', 'mslp', 'dewpoint_2m', 'visibility',\
    'ceiling_agl', 'cporain', 'cpoice', 'cposnow', 'precipw',\
    'snowh', 'windgust10m', 'cape', 'ptype', 'total_cloud_cover',\
    'rain_bucket','conv_bucket','zrain_bucket','snow_bucket','apcp_bucket']

# ---- variables actually stored in mpasout_05m.zarr; everything else is
#      in mpasout_15m.zarr.

bucket_05m_vars = ('rain_bucket', 'conv_bucket', 'zrain_bucket',
    'snow_bucket', 'apcp_bucket', 'total_cloud_cover', 'ptype')

# ---- precipitation-amount bucket variables: the raw 5-min buckets are
#      independent per-interval totals (they can drop back toward zero,
#      not a running accumulation -- confirmed by inspecting individual
#      mesh cells' raw values over time), so a single 5-min sample is a
#      noisy, hard-to-read snapshot.  read_timestep() sums three
#      consecutive 5-min buckets into a trailing 15-min accumulation for
#      these variables.  total_cloud_cover and ptype are also 5-min
#      bucket_05m_vars but aren't precipitation amounts, so they're left
#      as single-timestep reads.

precip_accum_vars = ('apcp_bucket', 'rain_bucket', 'conv_bucket',
    'zrain_bucket', 'snow_bucket')

ICfilename = 'GRAF_reforecast_initial_condition_dates_v2.txt'

# --------------------------------------------------------------------

def print_variable_help(cvar):
    print ('This variable ', cvar, ' is invalid.  Here are the valid ones ')
    print ('==================== variable list: ====================== ')
    print ('t2m: 2-meter temperature, valid every 15 min')
    print ('u10: 10-meter u-wind component, valid every 15 min')
    print ('v10: 10-meter v-wind component, valid every 15 min')
    print ('mslp: mean sea-level pressure, valid every 15 min')
    print ('dewpoint_2m: 2-meter dewpoint, valid every 15 min')
    print ('visibility: visibility, valid every 15 min')
    print ('ceiling_agl: ceiling above ground level, valid every 15 min')
    print ('cporain: conditional probability of rain, valid every 15 min')
    print ('cpoice: conditional probability of ice, valid every 15 min')
    print ('cposnow: conditional probability of snow, valid every 15 min')
    print ('precipw: total column precipitable water, valid every 15 min')
    print ('snowh: snow depth , valid every 15 min')
    print ('windgust10m: 10-m wind gust speed, valid every 15 min')
    print ('cape: CAPE, valid every 15 min')
    print ('ptype: predominant precipitation type, valid every 5 min')
    print ('total_cloud_cover: total cloud cover, valid every 5 min')
    print ('rain_bucket: total nonconvective rainfall, valid every 55 min')
    print ('conv_bucket: total convective rainvall, valid every 55 min')
    print ('zrain_bucket: total freezing rain , valid every 55 min')
    print ('snow_bucket: total liquid-equivalent snowfall, valid every 55 min')
    print ('apcp_bucket: total precipitation, valid every 55 min')
    print ('Exiting. ')

# --------------------------------------------------------------------

def lookup_case_directory(cyyyymmddhh, icfilename=ICfilename):

    """
    Look up cyyyymmddhh in the reforecast initial-condition list and
    return (directory_string, max_leadtime_hours).  Prints diagnostic
    info and returns None if the case day isn't found.
    """

    IClist = np.loadtxt(icfilename, usecols=0, dtype=int)
    print ('IClist[0:3] = ', IClist[0:3])
    leadtime_list = np.loadtxt(icfilename, usecols=1, dtype=int)
    ind = int(np.where(IClist == int(cyyyymmddhh))[0])

    if ind >= 0:
        directory_string = cyyyymmddhh+'_'+str(leadtime_list[ind])
        print (directory_string)
        return directory_string, int(leadtime_list[ind])
    else:
        cyyyymm = cyyyymmddhh[0:6]
        iyyyymm = int(cyyyymm)
        print ('That case day was not found.  ')
        print (' ')
        print ('Here are other GRAF casedays with the same year/month [yyyymm]: ',iyyyymm)
        print (' ')
        ndates = len(IClist)
        IClist_sub = []
        for idate in range(ndates):
            yyyymmddhh = IClist[idate]
            yyyymm_r = int(str(yyyymmddhh)[0:6])
            if yyyymm_r == iyyyymm:
                IClist_sub.append(yyyymmddhh)
        print (IClist_sub)
        return None, None

# --------------------------------------------------------------------

def load_static_latlon(latlonfile='rpm4km.static.latlon.nc'):

    """ Read lat/lon of the unstructured MPAS mesh grid points. """

    nc = Dataset(latlonfile, 'r')
    lats = nc.variables['latCell'][:]*180./3.141592656
    lons = nc.variables['lonCell'][:]*180./3.141592656
    print ('max, min lons = ', np.max(lons), np.min(lons))
    nc.close()
    return lons, lats

# --------------------------------------------------------------------

def open_case_store(directory_string, cvar):

    """
    Open the zarr store (05m- or 15m-cadence, depending on cvar) for one
    case, and decode + sort its Time coordinate.  Returns
    (zgroup, time_vals_sorted, sort_order) so that many timesteps can be
    read from the same open store without repeating this setup per frame.
    """

    if cvar in bucket_05m_vars:
        zarrdir = 's3://twc-graf-reforecast/2026-05-01.canonical/'+directory_string+'/mpasout_05m.zarr/'
    else:
        zarrdir = 's3://twc-graf-reforecast/2026-05-01.canonical/'+directory_string+'/mpasout_15m.zarr/'
    print (zarrdir)

    # --- Access the desired data directly through the zarr API rather than
    #     via xarray.open_zarr().  xarray builds a full Dataset (dask-backed
    #     variables + CF decoding) for every variable in the store, which
    #     costs several seconds of pure Python overhead regardless of how
    #     much data is actually requested.  Since only one variable is
    #     needed here, reading the zarr group directly skips that overhead.

    mapper = fsspec.get_mapper(zarrdir, anon=True)
    zgroup = zarr.open_consolidated(mapper, mode='r')

    # ---- Guard against any residual temporal mis-ordering: decode the
    #      Time coordinate and sort it.

    time_attrs = zgroup['Time'].attrs.asdict()
    time_ref = datetime.strptime(
        time_attrs['units'].split('since')[1].strip(), '%Y-%m-%d %H:%M:%S')
    time_vals = time_ref + np.array(
        [timedelta(minutes=int(minval)) for minval in zgroup['Time'][:]])
    sort_order = np.argsort(time_vals)
    time_vals_sorted = time_vals[sort_order]

    return zgroup, time_vals_sorted, sort_order

# --------------------------------------------------------------------

def read_timestep(zgroup, sort_order, time_index, cvar):

    """
    Read cvar from zgroup at one (already-time-sorted) timestep.  For the
    5-min precipitation bucket variables (precip_accum_vars), instead
    sums the current and two preceding 5-min buckets into a trailing
    15-min accumulation, clipped to whatever history is available at the
    start of the forecast (time_index 0 or 1).
    """

    if cvar in precip_accum_vars:
        start = max(0, time_index - 2)
        total = None
        for idx in range(start, time_index + 1):
            arr = zgroup[cvar][sort_order[idx], :]
            total = arr if total is None else total + arr
        return total
    else:
        actual_index = sort_order[time_index]
        return zgroup[cvar][actual_index, :]

# --------------------------------------------------------------------

def regrid_to_conus(data_array, inbox_idx, nn_idx, nn_dist, xg, yg, far_threshold_m=100000.):

    """
    Gather a full-mesh data_array onto the cached CONUS raster.  The MPAS
    mesh is global and variable-resolution: it's coarser (tens of km)
    outside the CONUS refinement region, which is normal, not a gap, so
    far_threshold_m only needs to catch genuine anomalies, not flag every
    coarse-mesh cell as missing.
    """

    gridded = data_array[inbox_idx][nn_idx].reshape(len(yg), len(xg))
    return np.ma.masked_where(nn_dist.reshape(len(yg), len(xg)) > far_threshold_m, gridded)

# --------------------------------------------------------------------

def make_basemap():
    return Basemap(rsphere=(6378137.00,6356752.3142),\
        resolution='l',area_thresh=1000.,projection='lcc',\
        lat_1=50.0,lat_2=20.0,lat_0=50.0,lon_0=255.-360.,\
        llcrnrlon=360-123.912,llcrnrlat=12.243,urcrnrlon=360.-39.897,\
        urcrnrlat=50.8549)  # GRAF CONUS domain.

# ------------------------------------------------------------------------

if __name__ == '__main__':

    cyyyymmddhh = sys.argv[1] # the desired initial cond. in year-month-day-hour format
    clead = sys.argv[2] # lead time in hours and fractions of hours, e.g., 12 or 12.25
    cvar = sys.argv[3]

    directory_string, max_leadtime = lookup_case_directory(cyyyymmddhh)
    if directory_string is None:
        sys.exit()

    try:
        idx = valid_variables.index(cvar)
    except ValueError:
        print_variable_help(cvar)
        sys.exit()

    plot_info = plotting_information(cyyyymmddhh, clead)
    title = plot_info['title'][cvar]
    plot_file = plot_info['plot_file'][cvar]
    clevs = plot_info['clevs'][cvar]
    scale_factor = plot_info['scale_factor'][cvar]
    offset_factor = plot_info['offset_factor'][cvar]
    clabel = plot_info['clabel'][cvar]

    try:
        lons, lats = load_static_latlon()

        if cvar in bucket_05m_vars:
            time_index = int(float(clead)*12) # was saved 12 times/hour, 5 min
        else:
            time_index = int(float(clead)*4) # was saved 4 times/hour, 15 min
        print ('time_index = ', time_index)

        zgroup, time_vals_sorted, sort_order = open_case_store(directory_string, cvar)

        yyyy, mm, dd, hh = splitdate(cyyyymmddhh)
        expected_time = datetime(yyyy, mm, dd, hh) + timedelta(hours=float(clead))
        selected_time = time_vals_sorted[time_index]
        print ('expected valid time = ', expected_time, ' ; selected valid time = ', selected_time)
        if abs((selected_time - expected_time).total_seconds()) > 60:
            print ('WARNING: selected timestep does not match expected valid time; '
                'data may be temporally mis-ordered.')

        data_array = read_timestep(zgroup, sort_order, time_index, cvar)
        print ('max, min data_array = ', np.max(data_array), \
            np.min(data_array))

        # --- Regrid onto the cached CONUS raster (see build_or_load_conus_grid)
        #     and plot.

        m = make_basemap()
        inbox_idx, xg, yg, nn_idx, nn_dist = build_or_load_conus_grid(m, lons, lats)
        XG, YG = np.meshgrid(xg, yg)

        gridded = regrid_to_conus(data_array, inbox_idx, nn_idx, nn_dist, xg, yg)

        istat_plot = plot_data((gridded-offset_factor)*scale_factor, \
            XG, YG, m, title, plot_file, clevs, cvar, clabel)

    except ValueError:
        print ('ValueError found.  Exiting.')
        sys.exit()



