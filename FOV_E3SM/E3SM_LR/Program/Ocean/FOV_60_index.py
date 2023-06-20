#Program determines the FOV index for 60N

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

section_name	= 'FOV_section_60N'
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
grid_x_norm	  = ma.masked_all((len(depth), len(lon)))

for depth_i in range(len(depth)):
	#Determine the surface area
	layer_field_area[depth_i]  	= layer_field[depth_i] * grid_x
    
	#Normalise the length
	grid_x_depth          		= ma.masked_array(grid_x, mask = v_vel[depth_i].mask)
	grid_x_norm[depth_i]  		= grid_x_depth / np.sum(grid_x_depth)

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_all		= ma.masked_all(len(time))
transport_salt_all	= ma.masked_all(len(time))

for file_i in range(len(files)):
	#Now determine for each month
	print(files[file_i])
	    
	lon, depth, layer_field, grid_x, v_vel, salt = ReadinData(files[file_i], depth_min_index, depth_max_index)

	#Determine the meridional transport
	transport	= v_vel * layer_field * grid_x

	#Determine the section averaged velocity (barotropic)
	vel_barotropic	= np.sum(transport) / np.sum(layer_field * grid_x)

	#Determine the overturning velocity (baroclinic)
	vel_baroclinic	 = v_vel - vel_barotropic

	#Determine the zonal means
	salt_zonal      = np.sum(salt * grid_x_norm, axis = 1)  - 35.0
	transport_clin	= np.sum(vel_baroclinic * layer_field * grid_x, axis = 1)

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[file_i]		= np.sum(transport) / 1000000.0
	        
	#Determine the total salinity transport
	transport_salt_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_clin * salt_zonal) / 1000000.0 

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FOV_index_'+section_name[4:]+'.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)
fh.createVariable('F_OV', float, ('time'), zlib=True)

fh.variables['Transport'].longname 	= 'Volume transport'
fh.variables['F_OV'].longname 		= 'Fresh water transport'

fh.variables['time'].units 		= 'Year'
fh.variables['Transport'].units 	= 'Sv'
fh.variables['F_OV'].units 		= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time
fh.variables['Transport'][:]    	= transport_all
fh.variables['F_OV'][:] 		= transport_salt_all

fh.close()
