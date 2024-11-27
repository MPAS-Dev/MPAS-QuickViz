#Program plots the mixed layer depth climatology

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Making pathway to folder with all data
directory 	= '../../Data/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	lon 		= fh.variables['lon'][:]	#Longitude
	lat 		= fh.variables['lat'][:]	#Latitude
	temp   		= fh.variables['TEMP'][:] 	#Temperature
	salt   		= fh.variables['SALT'][:] 	#Salinity
	mixed   	= fh.variables['MXL'][:] 	#Depth (m)
			
	fh.close()

	return lon, lat, temp, salt, mixed

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

files	= glob.glob(directory+'Data/Mixed_layer/E3SM_data_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  		= files[year_i][-7:-3]	
	year  		= int(date[0:4])
	time[year_i]	= year

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon, lat, temp, salt, mixed = ReadinData(files[0])

#-----------------------------------------------------------------------------------------

#Define empty array's
temp_all		= ma.masked_all((len(time)*12, len(lat), len(lon)))
salt_all		= ma.masked_all((len(time)*12, len(lat), len(lon)))
mixed_all		= ma.masked_all((len(time)*12, len(lat), len(lon)))

for file_i in range(len(files)):
	#Now determine for each month
	print(file_i)
	    
	lon, lat, temp, salt, mixed = ReadinData(files[file_i])
	
	for month_i in range(len(mixed)):
		#Add each month
		temp_all[file_i*12+month_i]	= temp[month_i]
		salt_all[file_i*12+month_i]	= salt[month_i]
		mixed_all[file_i*12+month_i]	= mixed[month_i]


temp_month	= ma.masked_all((12, len(lat), len(lon)))
salt_month	= ma.masked_all((12, len(lat), len(lon)))
mixed_month	= ma.masked_all((12, len(lat), len(lon)))

for month_i in range(12):
	#Loop over each month
	month_index		= np.arange(month_i, len(mixed_all), 12)
	temp_month[month_i]	= np.mean(temp_all[month_index], axis = 0)
	salt_month[month_i]	= np.mean(salt_all[month_index], axis = 0)
	mixed_month[month_i]	= np.mean(mixed_all[month_index], axis = 0)
	
#-----------------------------------------------------------------------------------------

mixed_crop			= 200
factor_mixed_crop		= 2
mixed_month[mixed_month > mixed_crop] 	= ((mixed_month[mixed_month > mixed_crop] - mixed_crop) / factor_mixed_crop) + mixed_crop

#-----------------------------------------------------------------------------------------

month	= ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

for month_i in range(12):

	#-----------------------------------------------------------------------------------------
	fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

	CS      = ax.contourf(lon, lat, temp_month[month_i] - temp_month[0], levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

	divider = make_axes_locatable(ax)
	ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
	fig.add_axes(ax_cb)

	cbar    = colorbar(CS, ticks = np.arange(-1, 1.01, 1), cax=ax_cb)
	cbar.set_label('Temperature difference ($^{\circ}$C)')

	gl = ax.gridlines(draw_labels=True)
	gl.top_labels = False
	gl.right_labels = False
	ax.set_extent([-70, 20, -71, 1], ccrs.PlateCarree())
	ax.coastlines('50m')
	ax.add_feature(cfeature.LAND, zorder=0)
	

	ax.set_title(month[month_i]+' minus January, E3SM Arctic')

	#-----------------------------------------------------------------------------------------
	
	fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

	CS      = ax.contourf(lon, lat, salt_month[month_i] - salt_month[0], levels = np.arange(-0.1, 0.101, 0.005), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

	divider = make_axes_locatable(ax)
	ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
	fig.add_axes(ax_cb)

	cbar    = colorbar(CS, ticks = np.arange(-0.1, 0.101, 0.1), cax=ax_cb)
	cbar.set_label('Salinity difference (g kg$^{-1}$)')

	gl = ax.gridlines(draw_labels=True)
	gl.top_labels = False
	gl.right_labels = False
	ax.set_extent([-70, 20, -71, 1], ccrs.PlateCarree())
	ax.coastlines('50m')
	ax.add_feature(cfeature.LAND, zorder=0)
	

	ax.set_title(month[month_i]+' minus January, E3SM Arctic')

	#-----------------------------------------------------------------------------------------
	
	fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

	CS      = ax.contourf(lon, lat, mixed_month[month_i], levels = np.arange(0, 400.1, 10), extend = 'max', cmap = 'Spectral_r', transform=ccrs.PlateCarree())

	divider = make_axes_locatable(ax)
	ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
	fig.add_axes(ax_cb)

	cbar    = colorbar(CS, ticks = [0, 100, 200, 300, 400], cax=ax_cb)
	cbar.ax.set_yticklabels([0, 100, 200, 400, 600])
	cbar.set_label('Mixed layer depth (m)')

	gl = ax.gridlines(draw_labels=True)
	gl.top_labels = False
	gl.right_labels = False
	ax.set_extent([-70, 20, -71, 1], ccrs.PlateCarree())
	ax.coastlines('50m')
	ax.add_feature(cfeature.LAND, zorder=0)
	

	ax.set_title(month[month_i]+', E3SM Arctic')
	show()

	#-----------------------------------------------------------------------------------------




