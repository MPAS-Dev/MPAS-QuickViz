#Program determines the MOV index

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
directory      = '/global/homes/r/rvwesten/E3SM/Data/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	lon 		= fh.variables['lon'][:]	#Longitude
	lat 		= fh.variables['lat'][:]	#Latitude
	temp   		= fh.variables['TEMP'][:] 	#Temperature
			
	fh.close()

	return lon, lat, temp

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

lat_min	= -30
lat_max	= -25
#-----------------------------------------------------------------------------------------

files	= glob.glob(directory+'Data/Mixed_layer/E3SM_data_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files)*12)

for year_i in range(len(files)):
	date  		= files[year_i][-7:-3]	
	year  		= int(date[0:4])
	
	for month_i in range(12):
		time[year_i*12+month_i]	= year + month_i / 12.0

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon, lat, temp = ReadinData(files[0])

lat_min_index	= (np.abs(lat - lat_min)).argmin()
lat_max_index	= (np.abs(lat - lat_max)).argmin()+1

#-----------------------------------------------------------------------------------------

#Define empty array's
temp_all		= ma.masked_all((len(time), len(lon)))

for file_i in range(len(files)):
	#Now determine for each month
	print(file_i)
	    
	lon, lat, temp = ReadinData(files[file_i])
	
	for month_i in range(len(temp)):
		#Add each month
		temp_all[file_i*12+month_i]	= np.mean(temp[month_i, lat_min_index:lat_max_index], axis = 0)

#Now remove the monthly mean
for month_i in range(12):
	time_index		= np.arange(month_i, len(time), 12)
	temp_mean		= np.mean(temp_all[time_index], axis = 0)
	temp_all[time_index]	= temp_all[time_index] - temp_mean
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= contourf(lon, time, temp_all, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'RdBu_r')
cbar	= colorbar(CS)

ax.set_ylim(500, 527)

show()



