#Program determines the ACC strength

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory 	= '../../Data/'

def ReadinData(filename, depth_min_index, depth_max_index):

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lat 		= fh.variables['lat'][:]					#Longitude
	depth   	= fh.variables['depth'][depth_min_index:depth_max_index] 	#Depth (m)
	layer		= fh.variables['layer'][depth_min_index:depth_max_index] 	#Layer thickness (m)
	grid_y		= fh.variables['DY'][:] 					#Meridional grid cell length (m)
	u_vel 		= fh.variables['UVEL'][depth_min_index:depth_max_index] 	#Zonal velocity (m/s)

	fh.close()

	return lat, depth, layer, grid_y, u_vel
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 6000

#-----------------------------------------------------------------------------------------

files	= glob.glob(directory+'Data/Drake_Passage/E3SM_data_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  		= files[year_i][-7:-3]	
	year  		= int(date[0:4])
	time[year_i]	= year

#-----------------------------------------------------------------------------------------

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

depth   	= fh.variables['depth'][:]	#Depth (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lat, depth, layer_field, grid_y, u_vel = ReadinData(files[0], depth_min_index, depth_max_index)

for lat_i in range(len(lat)):
	#Get all the layers which have a maximum depth below given range
	if np.sum(layer_field[:, lat_i]) > depth_max:
		#Adjust the last layer
		layer_field[-1, lat_i]	-= (np.sum(layer_field[:, lat_i]) - depth_max)

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_all		= ma.masked_all(len(time))

for time_i in range(len(time)):
	#Now determine for each month
	print(time_i)
	    
	lat, depth, layer_field_old, grid_y, u_vel = ReadinData(files[time_i], depth_min_index, depth_max_index)

	#Determine the meridional transport
	transport		= u_vel * layer_field * grid_y

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[time_i]	= np.sum(transport) / 1000000.0

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/ACC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'_m.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)

fh.variables['Transport'].long_name 	= 'Volume transport'

fh.variables['time'].units 		= 'Year'
fh.variables['Transport'].units 	= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time
fh.variables['Transport'][:]    	= transport_all

fh.close()
