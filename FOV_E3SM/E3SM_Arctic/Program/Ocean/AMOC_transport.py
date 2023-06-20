#Program determines the AMOC strength

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
	lon 		= fh.variables['lon'][:]					#Longitude
	depth   	= fh.variables['depth'][depth_min_index:depth_max_index] 	#Depth (m)
	layer		= fh.variables['layer'][depth_min_index:depth_max_index] 	#Layer thickness (m)
	grid_x		= fh.variables['DX'][:] 					#Zonal grid cell length (m)
	v_vel 		= fh.variables['VVEL'][depth_min_index:depth_max_index] 	#Meridional velocity (m/s)
	salt		= fh.variables['SALT'][depth_min_index:depth_max_index] 	#Salinity (g / kg)

	fh.close()

	return lon, depth, layer, grid_x, v_vel
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 1000

lat_FOV 	= 26
section_name	= 'FOV_section_26N'
#-----------------------------------------------------------------------------------------

files	= glob.glob(directory+'Data/'+section_name+'/E3SM_data_year_*.nc')
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
lon, depth, layer_field, grid_x, v_vel = ReadinData(files[0], depth_min_index, depth_max_index)

for lon_i in range(len(lon)):
	#Get all the layers which have a maximum depth below given range
	if np.sum(layer_field[:, lon_i]) > depth_max:
		#Adjust the last layer
		layer_field[-1, lon_i]	-= (np.sum(layer_field[:, lon_i]) - depth_max)

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_all		= ma.masked_all(len(time))

for time_i in range(len(time)):
	#Now determine for each month
	print(time_i)
	    
	lon, depth, layer_field_old, grid_x, v_vel = ReadinData(files[time_i], depth_min_index, depth_max_index)

	#Determine the meridional transport
	transport		= v_vel * layer_field * grid_x

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[time_i]	= np.sum(transport) / 1000000.0

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'_m.nc', 'w')

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
