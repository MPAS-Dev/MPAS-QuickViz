#Program determines the azonal (gyre) component at 34S

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

	lon 		= fh.variables['lon'][:]					                #Longitude
	depth   	= fh.variables['depth'][depth_min_index:depth_max_index] 	#Depth (m)
	layer		= fh.variables['layer'][depth_min_index:depth_max_index] 	#Layer thickness (m)
	grid_x		= fh.variables['DX'][:] 					                #Zonal grid cell length (m)
	v_vel 		= fh.variables['VVEL'][depth_min_index:depth_max_index] 	#Meridional velocity (m/s)
	salt		= fh.variables['SALT'][depth_min_index:depth_max_index] 	#Salinity (g / kg)

	fh.close()

	return lon, depth, layer, grid_x, v_vel, salt
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 6000

section_name	= 'FOV_section_34S'
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
lon, depth, layer_field, grid_x, v_vel, salt = ReadinData(files[0], depth_min_index, depth_max_index)

#Normalise layer field per layer
layer_field_area  = ma.masked_all(shape(layer_field))

for depth_i in range(len(depth)):
	#Determine the surface area
	layer_field_area[depth_i]  	= layer_field[depth_i] * grid_x

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_gyre_all		= ma.masked_all(len(time))

for file_i in range(len(files)):
	#Now determine for each month
	print(files[file_i])
	    
	lon, depth, layer_field, grid_x, v_vel, salt = ReadinData(files[file_i], depth_min_index, depth_max_index)

	#Determine the zonal means
	v_vel_zonal 	= np.mean(v_vel, axis = 1)
	salt_zonal	= np.mean(salt, axis = 1)
	
	v_vel_prime	= ma.masked_all(np.shape(v_vel))
	salt_prime	= ma.masked_all(np.shape(salt))
	
	for depth_i in range(len(depth)):
		#Determine the differences with respect to the zonal means
		v_vel_prime[depth_i]	= v_vel[depth_i] - v_vel_zonal[depth_i]
		salt_prime[depth_i]	= salt[depth_i] - salt_zonal[depth_i]

	#Now determine the azonal component (gyre, in Sv)
	transport_gyre_all[file_i]	= (-1.0 / 35.0) * np.sum(v_vel_prime * salt_prime * layer_field_area) / 10**6.0

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FOV_gyre_section_34S.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('F_gyre', float, ('time'), zlib=True)

fh.variables['F_gyre'].longname 	= 'Freshwater transport by gyre'

fh.variables['time'].units 		= 'Year'
fh.variables['F_gyre'].units 		= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time
fh.variables['F_gyre'][:] 		= transport_gyre_all

fh.close()
