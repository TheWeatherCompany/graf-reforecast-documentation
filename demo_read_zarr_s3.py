"""
demo_read_zarr_s3.py

This program will read selected GRAF reforecast files saved
in the S3 data store and plot.	 It does not read upper-air 
data with multiple vertical levels, but rather demo's 
the capability with the single-variable data.

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

"""
from netCDF4 import Dataset
import os, sys
import numpy as np
import zarr
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from dateutils import splitdate
from string import Formatter
import s3fs
import fsspec
import boto3 # AWS SDK for Python; may need to pip install boto3 on AWS
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

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
			'rain_bucket': '5-min cotal rainfall, IC = '+cyyyymmddhh+', '+clead+' h',
			'conv_bucket': '5-min convective precipitation, IC = '+cyyyymmddhh+', '+clead+' h',
			'zrain_bucket': '5-min freezing rain accumulation, IC = '+cyyyymmddhh+', '+clead+' h',
			'snow_bucket': '5-min snow liquid water, IC = '+cyyyymmddhh+', '+clead+' h',
			'apcp_bucket': '5-min total precipitation, IC = '+cyyyymmddhh+', '+clead+' h'},
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

	return plot_info

# --------------------------------------------------------------------

def plot_data(data_array, lons, lats, title, plot_file, clevs, cvar, clabel):

	#  --- plot hourly GRAF precipitation amount in CONUS/North American
	#	   domain.	Save to png file.

	colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
		'#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
		'#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
	m = Basemap(rsphere=(6378137.00,6356752.3142),\
		resolution='l',area_thresh=1000.,projection='lcc',\
		lat_1=50.0,lat_2=20.0,lat_0=50.0,lon_0=255.-360.,\
		llcrnrlon=360-123.912,llcrnrlat=12.243,urcrnrlon=360.-39.897,\
		urcrnrlat=50.8549)	# GRAF CONUS domain.
		
	x, y = m(lons, lats)
	
	fig = plt.figure(figsize=(9,9.))
	axloc = [0.02,0.08,0.96,0.84]
	ax1 = fig.add_axes(axloc)
	ax1.set_title(title, fontsize=15,color='Black')
	CS2 = m.contourf(x, y, data_array, clevs,\
		cmap=None, colors=colorst, extend='both', tri=True)
	#CS2 = m.contourf(x, y, data_array, clevs,\
	#	cmap=None, colors=colorst, extend='both')

	m.drawcoastlines(linewidth=0.9,color='Gray')
	m.drawcountries(linewidth=0.9,color='Gray')
	m.drawstates(linewidth=0.7,color='Gray')
	#m.drawcounties(linewidth=0.1,color='LightGray')

	cax = fig.add_axes([0.02,0.055,0.96,0.02])
	cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
		drawedges=True,ticks=clevs,format='%g')
	cb.ax.tick_params(labelsize=10)
	cb.set_label(clabel, fontsize=12)

	if cvar == 'ptype':
		cb.ax.set_xticklabels(['					 Rain',\
			'						 Ice/Graupel',\
			'				 Sleet / Ice Pellets',\
			'						 Snow',' '])  # horizontal colorbar

	# ---- set plot title

	fig.savefig(plot_file, dpi=400, bbox_inches='tight')
	print ('saving plot to file = ',plot_file)
	print ('Done!')
	istat = 0
	return istat

# ------------------------------------------------------------------------

cyyyymmddhh = sys.argv[1] # the desired initial cond. in year-month-day-hour format
clead = sys.argv[2] # lead time in hours and fractions of hours, e.g., 12 or 12.25
cvar = sys.argv[3]

# --- define a dictionary with plotting information.

valid_variables = ['t2m', 'u10', 'v10', 'mslp', 'dewpoint_2m', 'visibility',\
	'ceiling_agl', 'cporain', 'cpoice', 'cposnow', 'precipw',\
	'snowh', 'windgust10m', 'cape', 'ptype', 'total_cloud_cover',\
	'rain_bucket','conv_bucket','zrain_bucket','snow_bucket','apcp_bucket']

# ---- read the list of initial condition dates and lead times.

ICfilename = 'GRAF_reforecast_initial_condition_dates_v2.txt'
IClist = np.loadtxt(ICfilename,usecols=0,dtype=int)
print ('IClist[0:3] = ', IClist[0:3])
leadtime_list = np.loadtxt(ICfilename,usecols=1,dtype=int)
ind = int(np.where(IClist == int(cyyyymmddhh))[0])

if ind >= 0:
	directory_string = cyyyymmddhh+'_'+str(leadtime_list[ind])
	print (directory_string)
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
	
	
	sys.exit()

try:
	idx = valid_variables.index(cvar)
except ValueError:
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
	print ('ptype: predominant precipitation type, valid every 15 min')	 
	print ('total_cloud_cover: total cloud cover, valid every 5 min') 
	print ('rain_bucket: total nonconvective rainfall, valid every 55 min') 
	print ('conv_bucket: total convective rainvall, valid every 55 min') 
	print ('zrain_bucket: total freezing rain , valid every 55 min') 
	print ('snow_bucket: total liquid-equivalent snowfall, valid every 55 min') 
	print ('apcp_bucket: total precipitation, valid every 55 min') 
	print ('Exiting. ')
	sys.exit()
	
plot_info = plotting_information(cyyyymmddhh, clead)
title = plot_info['title'][cvar]
plot_file = plot_info['plot_file'][cvar]
clevs = plot_info['clevs'][cvar]
scale_factor = plot_info['scale_factor'][cvar]
offset_factor = plot_info['offset_factor'][cvar]
clabel = plot_info['clabel'][cvar]
	
try:
	
	# ---- Read from a netcdf file the lat/lon of the unstructured mesh
	#	   grid points.
	
	latlonfile = '/s3data/inference/corrdiff/graf/input/rpm4km.static.nc'
	nc = Dataset(latlonfile, 'r')
	lats = nc.variables['latCell'][:]*180./3.141592656
	lons = nc.variables['lonCell'][:]*180./3.141592656
	print ('max, min lons = ', np.max(lons), np.min(lons))
	npoints = len(lons)
	nc.close()
	
	# ---- See if cyyyymmddhh actually was an initial condition case in this.	
	#	   GRAF reforecast. If so, proceed with reading and generating plots.
	# 
	#cmd = 'aws s3 ls s3://twc-graf-reforecast/ | '+\
	#	'awk "{print $2}" | awk -F "/ {print $1}" > output.txt'
	#istat = os.system(cmd)
	#print (cmd)
	#infile = open("output.txt", "r")
	#line = infile.readline()
	#print (line)
	#infile.close()
	#i = line.index(cyyyymmddhh)
	# Connect to AWS S3 storage (anonymous login for public data)
	
	fs = s3fs.S3FileSystem(anon=True) #
	
	# ---- Different files are read depending on the variable, some of which 
	#	   are 5-min data, some are 15.
	
	if cvar == 'rain_bucket' or cvar == 'conv_bucket' or\
	cvar == 'zrain_bucket' or cvar == 'snow_bucket' or\
	cvar == 'apcp_bucket' or cvar == 'total_cloud_cover':
		#zarrdir = 's3://twc-graf-reforecast/'+line[i:i+13]+'/mpasout_05m.zarr/'
		zarrdir = 's3://twc-graf-reforecast/'+directory_string+'/mpasout_05m.zarr/'
	else:
		#zarrdir = 's3://twc-graf-reforecast/'+line[i:i+13]+'/mpasout_15m.zarr/'
		zarrdir = 's3://twc-graf-reforecast/'+directory_string+'/mpasout_15m.zarr/'
	print (zarrdir)
	
	# ---- Create a MutableMapping from the store URL
	
	mapper = fsspec.get_mapper(zarrdir, anon=True) #
	
	# ---- Determine the time index from the lead time in hours.fraction_of_hours
	
	if cvar == 'rain_bucket' or cvar == 'conv_bucket' or \
	cvar == 'zrain_bucket' or cvar == 'snow_bucket' or \
	cvar == 'total_cloud_cover' or cvar =='apcp_bucket':
		time_index = int(float(clead)*12) # was saved 12 times/hour, 5 min
	else:
		time_index = int(float(clead)*4) # was saved 4 times/hour, 15 min
	print ('time_index = ', time_index)
	
	# --- Access the desired data from 'data' object in zarr

	ds = xr.open_zarr(mapper, consolidated=True) 
	print ('ds.data_vars = ',ds.data_vars)
	print ('variable = ', cvar,' time index = ', time_index) 
	ds_slice = ds[cvar][time_index, :]
	print (ds_slice)
	data_array = ds_slice.data
	print ('max, min data_array = ', np.max(data_array), np.min(data_array))
	#print ('data_array[0::10000] = ', data_array[0::10000])
	#sys.exit()

	# --- Plot the data.
	
	istat_plot = plot_data((data_array-offset_factor)*scale_factor, \
		lons, lats, title, plot_file, clevs, cvar, clabel)
		
except ValueError: 
	print ('ValueError found.  Exiting.')
	sys.exit()



