#Program plots the vertically averaged (upper 100 m) salinity in the Southern Oceaan

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from scipy.stats import genextreme
from matplotlib.colors import LogNorm
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable


#Making pathway to folder with all data
directory 	= '../../Data/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min	= 0
depth_max	= 100
year_start	= 1
year_end	= 5

files	= glob.glob(directory+'Data/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m/E3SM_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  = files[year_i][-7:-3]	
	year  = int(date[0:4])

	time[year_i] = year

time_start	= (np.abs(time - year_start)).argmin()
time_end	= (np.abs(time - year_end)).argmin() + 1
files		= files[time_start:time_end]

#-----------------------------------------------------------------------------------------

for file_i in range(len(files)):
	print(files[file_i])
	fh = netcdf.Dataset(files[file_i], 'r')

	lon	= fh.variables['lon'][:]	
	lat	= fh.variables['lat'][:]	
	salt	= fh.variables['SALT'][:]		#Salinity

	fh.close()

	if file_i == 0:
		salt_all	= ma.masked_all((len(files), len(lat), len(lon)))

	salt_all[file_i]	= salt

salt_all	= np.mean(salt_all, axis = 0)
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, salt_all, levels = np.arange(33, 37.1, 0.1), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(33, 37.1, 1), cax=ax_cb)
cbar.set_label('Salinity (g kg$^{-1}$)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)
ax.set_title('Salinity (0 - 100 m), LR-E3SM ('+str(year_start)+' - '+str(year_end)+')')

show()

