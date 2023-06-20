#Program determines the FOV index for 34S and the difference components

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
grid_x_norm	  = ma.masked_all((len(depth), len(lon)))

for depth_i in range(len(depth)):
	#Determine the surface area
	layer_field_area[depth_i]  	= layer_field[depth_i] * grid_x
    
	#Normalise the length
	grid_x_depth          		= ma.masked_array(grid_x, mask = v_vel[depth_i].mask)
	grid_x_norm[depth_i]  		= grid_x_depth / np.sum(grid_x_depth)

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_all		     = ma.masked_all(len(time))
transport_salt_all	     = ma.masked_all(len(time))
transport_salt_ASW_all	 = ma.masked_all(len(time))
transport_salt_AIW_all	 = ma.masked_all(len(time))
transport_salt_NADW_all	 = ma.masked_all(len(time))
transport_salt_ABW_all	 = ma.masked_all(len(time))
salt_ASW_all		     = ma.masked_all(len(time))
salt_AIW_all		     = ma.masked_all(len(time))
salt_NADW_all		     = ma.masked_all(len(time))
salt_ABW_all		     = ma.masked_all(len(time))
vel_ASW_all		         = ma.masked_all(len(time))
vel_AIW_all		         = ma.masked_all(len(time))
vel_NADW_all		     = ma.masked_all(len(time))
vel_ABW_all		         = ma.masked_all(len(time))

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

	#-----------------------------------------------------------------------------------------
	#Get the water properties
	water_prop	= ma.masked_all((len(depth), len(lon)))

	#North Atlantic Deep Water (NADW) has negative meridional velocities
	depth_index_NADW = np.where((depth >= 700) & (transport_clin <= 0))[0][0]

	#Antarctic bottom water (ABW) is directly below the NADW, get the first index
	depth_index_ABW	= np.where((depth >= 3000) & (transport_clin >= 0))[0]

	if len(depth_index_ABW) == 0:
		#Assume below 4000m depth the ABW
		depth_index_ABW	= np.where(depth >= 4000)[0][0]
	else:
		depth_index_ABW	= depth_index_ABW[0]

	for depth_i in range(len(depth)):
			
		if depth_i < depth_index_NADW:
			#Surface water
			water_prop[depth_i]	= 1.0

		if depth[depth_i] >= 500 and depth_i < depth_index_NADW:
			#Antarctic Intermediate water
			water_prop[depth_i]	= 2.0
		
		if depth_i >= depth_index_NADW and depth_i < depth_index_ABW:
			#North Atlantic Deep Water (NADW)
			water_prop[depth_i]	= 3.0	

		if depth_i >= depth_index_ABW:
			#The ABW is defined below the NADW 
			water_prop[depth_i]	= 4.0	

	water_prop	= ma.masked_array(water_prop, mask = v_vel.mask)	

	#-----------------------------------------------------------------------------------------
	area_ASW	= ma.masked_where(water_prop != 1.0, layer_field_area)
	area_AIW	= ma.masked_where(water_prop != 2.0, layer_field_area)
	area_NADW	= ma.masked_where(water_prop != 3.0, layer_field_area)
	area_ABW	= ma.masked_where(water_prop != 4.0, layer_field_area)
	area_ASW	= area_ASW	/ np.sum(area_ASW)
	area_AIW	= area_AIW	/ np.sum(area_AIW)
	area_NADW	= area_NADW	/ np.sum(area_NADW)
	area_ABW	= area_ABW	/ np.sum(area_ABW)

	#Determine the spatial means
	vel_ASW_all[file_i]	= np.sum(vel_baroclinic * area_ASW)
	vel_AIW_all[file_i]	= np.sum(vel_baroclinic * area_AIW)
	vel_NADW_all[file_i]	= np.sum(vel_baroclinic * area_NADW)
	vel_ABW_all[file_i]	= np.sum(vel_baroclinic * area_ABW)
	salt_ASW_all[file_i]	= np.sum(salt * area_ASW)
	salt_AIW_all[file_i]	= np.sum(salt * area_AIW)
	salt_NADW_all[file_i]	= np.sum(salt * area_NADW)
	salt_ABW_all[file_i]	= np.sum(salt * area_ABW)

	#Determine the means over the water masses
	transport_ASW	= np.sum(ma.masked_where(water_prop != 1.0, vel_baroclinic * layer_field * grid_x), axis = 1)
	transport_AIW	= np.sum(ma.masked_where(water_prop != 2.0, vel_baroclinic * layer_field * grid_x), axis = 1)
	transport_NADW	= np.sum(ma.masked_where(water_prop != 3.0, vel_baroclinic * layer_field * grid_x), axis = 1)
	transport_ABW	= np.sum(ma.masked_where(water_prop != 4.0, vel_baroclinic * layer_field * grid_x), axis = 1)

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[file_i]		= np.sum(transport) / 1000000.0
	        
	#Determine the total salinity transport
	transport_salt_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_clin * salt_zonal) / 1000000.0 
	transport_salt_ASW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_ASW * salt_zonal) / 1000000.0 
	transport_salt_AIW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_AIW * salt_zonal) / 1000000.0 
	transport_salt_NADW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_NADW * salt_zonal) / 1000000.0 
	transport_salt_ABW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_ABW * salt_zonal) / 1000000.0 

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FOV_index_'+section_name[4:]+'.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)
fh.createVariable('F_OV', float, ('time'), zlib=True)
fh.createVariable('F_OV_ASW', float, ('time'), zlib=True)
fh.createVariable('F_OV_AIW', float, ('time'), zlib=True)
fh.createVariable('F_OV_NADW', float, ('time'), zlib=True)
fh.createVariable('F_OV_ABW', float, ('time'), zlib=True)
fh.createVariable('SALT_ASW', float, ('time'), zlib=True)
fh.createVariable('SALT_AIW', float, ('time'), zlib=True)
fh.createVariable('SALT_NADW', float, ('time'), zlib=True)
fh.createVariable('SALT_ABW', float, ('time'), zlib=True)
fh.createVariable('VVEL_ASW', float, ('time'), zlib=True)
fh.createVariable('VVEL_AIW', float, ('time'), zlib=True)
fh.createVariable('VVEL_NADW', float, ('time'), zlib=True)
fh.createVariable('VVEL_ABW', float, ('time'), zlib=True)

fh.variables['Transport'].longname 	= 'Volume transport'
fh.variables['F_OV'].longname 		= 'Fresh water transport'
fh.variables['F_OV_ASW'].longname 	= 'Fresh water transport (Atlantic Surface Water)'
fh.variables['F_OV_AIW'].longname 	= 'Fresh water transport (Antarctic Intermediate Water)'
fh.variables['F_OV_NADW'].longname 	= 'Fresh water transport (North Atlantic Deep Water)'
fh.variables['F_OV_ABW'].longname 	= 'Fresh water transport (Antarctic Bottom Water)'
fh.variables['SALT_ASW'].longname 	= 'Salinity (Atlantic Surface Water)'
fh.variables['SALT_AIW'].longname 	= 'Salinity (Antarctic Intermediate Water)'
fh.variables['SALT_NADW'].longname 	= 'Salinity (North Atlantic Deep Water)'
fh.variables['SALT_ABW'].longname 	= 'Salinity (Antarctic Bottom Water)'
fh.variables['VVEL_ASW'].longname 	= 'Meridional velocity (Atlantic Surface Water)'
fh.variables['VVEL_AIW'].longname 	= 'Meridional velocity (Antarctic Intermediate Water)'
fh.variables['VVEL_NADW'].longname 	= 'Meridional velocity (North Atlantic Deep Water)'
fh.variables['VVEL_ABW'].longname 	= 'Meridional velocity (Antarctic Bottom Water)'

fh.variables['time'].units 		= 'Year'
fh.variables['Transport'].units 	= 'Sv'
fh.variables['F_OV'].units 		= 'Sv'
fh.variables['F_OV_ASW'].units 		= 'Sv'
fh.variables['F_OV_AIW'].units 		= 'Sv'
fh.variables['F_OV_NADW'].units 	= 'Sv'
fh.variables['F_OV_ABW'].units 		= 'Sv'
fh.variables['SALT_ASW'].units 		= 'g/kg'
fh.variables['SALT_AIW'].units 		= 'g/kg'
fh.variables['SALT_NADW'].units 	= 'g/kg'
fh.variables['SALT_ABW'].units 		= 'g/kg'
fh.variables['VVEL_ASW'].units 		= 'cm/s'
fh.variables['VVEL_AIW'].units 		= 'cm/s'
fh.variables['VVEL_NADW'].units 	= 'cm/s'
fh.variables['VVEL_ABW'].units 		= 'cm/s'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time
fh.variables['Transport'][:]    	= transport_all
fh.variables['F_OV'][:] 		= transport_salt_all
fh.variables['F_OV_ASW'][:] 		= transport_salt_ASW_all
fh.variables['F_OV_AIW'][:] 		= transport_salt_AIW_all
fh.variables['F_OV_NADW'][:] 		= transport_salt_NADW_all
fh.variables['F_OV_ABW'][:] 		= transport_salt_ABW_all
fh.variables['SALT_ASW'][:] 		= salt_ASW_all
fh.variables['SALT_AIW'][:] 		= salt_AIW_all
fh.variables['SALT_NADW'][:] 		= salt_NADW_all
fh.variables['SALT_ABW'][:] 		= salt_ABW_all
fh.variables['VVEL_ASW'][:] 		= vel_ASW_all * 100.0
fh.variables['VVEL_AIW'][:] 		= vel_AIW_all * 100.0
fh.variables['VVEL_NADW'][:] 		= vel_NADW_all * 100.0
fh.variables['VVEL_ABW'][:] 		= vel_ABW_all * 100.0

fh.close()
