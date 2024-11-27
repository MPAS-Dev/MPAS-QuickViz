#Generates the FOV fields on the 0.5x0.5 rectangular grid
#Note that you also need the Layer grid file, first interpolate the variable timeMonthly_avg_layerThickness for one particular month
#Place this file in the corresponding directory
#You probably need to change the interpolated output directory path

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy.interpolate import griddata

directory_data  	= '/global/cfs/cdirs/m1199/e3sm-arrm-simulations/E3SMv2.1B60to10rA02/ocn/hist/'
directory 	        = '../../Data/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

files = glob.glob(directory_data+'E3SMv2.1B60to10rA02.mpaso.hist.am.timeSeriesStatsMonthly.*.nc')
files.sort()

files	= files[2340:]

print(files[0])

#-----------------------------------------------------------------------------------------

#Define empty array's
time            = np.zeros(len(files))

for year_i in range(len(files)):
        date  = files[year_i][-13:-3]
        year  = int(date[0:4])
        month = int(date[5:7])

        time[year_i] = year + (month-1) / 12.0
      
 
#-----------------------------------------------------------------------------------------

fh	= netcdf.Dataset(directory+'Data/Layer_grid.nc', 'r')

layer	= fh.variables['timeMonthly_avg_layerThickness'][0] 	#Layer thickness (m)

fh.close()

#Use general depth coordinate
depth	= np.zeros(len(layer)+1)

for depth_i in range(1, len(depth)):
	#Generate the depth boundaries
	depth[depth_i]	= depth[depth_i-1] + layer[depth_i-1, 87, 255]

#Take the mean to find general depth array
depth 	= 0.5 * (depth[1:] + depth[:-1])
		
#-----------------------------------------------------------------------------------------		
		
time_year   = ma.masked_all(int(len(time)/12))


for year_i in range(int(np.min(time)), int(np.min(time))+len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i - int(np.min(time))] = year_i
	files_month = glob.glob(directory_data+'E3SMv2.1B60to10rA02.mpaso.hist.am.timeSeriesStatsMonthly.'+str(year_i).zfill(4)+'-*.nc')
	files_month.sort()

	for file_i in range(len(files_month)):
		#Loop over each month
		os.system('ncremap -i '+files_month[file_i]+' -P mpas -m /global/cfs/cdirs/e3sm/diagnostics/mpas_analysis/maps/map_ARRM10to60E2r1_to_0.5x0.5degree_bilinear.nc -T . -O /global/homes/r/rvwesten/E3SM_Arctic/Program/Ocean -o Regrid_month.nc -v timeMonthly_avg_activeTracers_salinity,timeMonthly_avg_velocityMeridional')

		fh	= netcdf.Dataset('Regrid_month.nc', 'r')

		lon 	= fh.variables['lon'][:]
		lat	    = fh.variables['lat'][:]
		salt	= fh.variables['timeMonthly_avg_activeTracers_salinity'][0] 	#Salinity (g/kg)
		v_vel	= fh.variables['timeMonthly_avg_velocityMeridional'][0]		#Meridional velocity (m/s)

		fh.close()
		
		salt	= ma.masked_where(layer <= 0.0, salt)
		v_vel	= ma.masked_where(layer <= 0.0, v_vel)
		
		for lat_section in [-34, 26, 60]:
			#Get the lat index			
			lat_index	= (np.abs(lat - lat_section)).argmin()
			
			if lat_section == -34:
				#Section at 34S, start of Atlantic Sector
				lon_1, lon_2    = 250, 401
				section_name    = 'FOV_section_34S'
				
				if year_i == int(np.min(time)):
					#Get the layer for the section
					lon_34S		= lon[lon_1:lon_2]
					layer_34S	= layer[:, lat_index, lon_1:lon_2]
					dx_34S		= 6371000 * 2 * np.pi * np.cos(lat[lat_index] * np.pi / 180) * 0.5 / 360 + np.zeros(len(lon_34S))
					
			if lat_section == 26:
				#Section at 26N, RAPID array
				lon_1, lon_2    = 198, 335
				section_name    = 'FOV_section_26N'

				if year_i == int(np.min(time)):
					#Get the layer for the section
					lon_26N		= lon[lon_1:lon_2]
					layer_26N	= layer[:, lat_index, lon_1:lon_2]
					dx_26N		= 6371000 * 2 * np.pi * np.cos(lat[lat_index] * np.pi / 180) * 0.5 / 360 + np.zeros(len(lon_26N))
			if lat_section == 60:
				#Section at 60N, RAPID array
				lon_1, lon_2    = 230, 373
				section_name    = 'FOV_section_60N'
			
				if year_i == int(np.min(time)):
					#Get the layer for the section
					lon_60N		= lon[lon_1:lon_2]
					layer_60N	= layer[:, lat_index, lon_1:lon_2]
					dx_60N		= 6371000 * 2 * np.pi * np.cos(lat[lat_index] * np.pi / 180) * 0.5 / 360 + np.zeros(len(lon_60N))   
			if file_i == 0 and lat_section == -34:
				#Make empty arrays for the months
				v_vel_34S	= ma.masked_all((12, len(depth), lon_2 - lon_1))
				salt_34S	= ma.masked_all((12, len(depth), lon_2 - lon_1))

			if file_i == 0 and lat_section == 26:
				#Make empty arrays for the months
				v_vel_26N	= ma.masked_all((12, len(depth), lon_2 - lon_1))
				salt_26N	= ma.masked_all((12, len(depth), lon_2 - lon_1))

			if file_i == 0 and lat_section == 60:
				#Make empty arrays for the months
				v_vel_60N	= ma.masked_all((12, len(depth), lon_2 - lon_1))
				salt_60N	= ma.masked_all((12, len(depth), lon_2 - lon_1))
			
			if lat_section == -34:
				#Now save the data to the general array
				v_vel_34S[file_i]	= v_vel[:, lat_index, lon_1:lon_2]
				salt_34S[file_i]	= salt[:, lat_index, lon_1:lon_2]
				
			if lat_section == 26:
				#Now save the data to the general array
				v_vel_26N[file_i]	= v_vel[:, lat_index, lon_1:lon_2]
				salt_26N[file_i]	= salt[:, lat_index, lon_1:lon_2]

			if lat_section == 60:
				#Now save the data to the general array
				v_vel_60N[file_i]	= v_vel[:, lat_index, lon_1:lon_2]
				salt_60N[file_i]	= salt[:, lat_index, lon_1:lon_2]
		
		os.remove('Regrid_month.nc')												
	#------------------------------------------------------------------------------
	#Now convert to yearly averages
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])

	#Fill the array's with the same dimensions
	month_days_34S	= ma.masked_all((len(month_days), len(depth), len(lon_34S)))
	month_days_26N	= ma.masked_all((len(month_days), len(depth), len(lon_26N)))
	month_days_60N	= ma.masked_all((len(month_days), len(depth), len(lon_60N)))
		
	for month_i in range(len(month_days)):
		month_days_34S[month_i]		= month_days[month_i]
		month_days_26N[month_i]		= month_days[month_i]
		month_days_60N[month_i]		= month_days[month_i]
		
	#Now set mask
	month_days_34S	= ma.masked_array(month_days_34S, mask = salt_34S.mask)
	month_days_26N	= ma.masked_array(month_days_26N, mask = salt_26N.mask)
	month_days_60N	= ma.masked_array(month_days_60N, mask = salt_60N.mask)
		
	#Normalise the data
	month_days_34S	= month_days_34S / np.sum(month_days_34S, axis = 0)
	month_days_26N	= month_days_26N / np.sum(month_days_26N, axis = 0)				
	month_days_60N	= month_days_60N / np.sum(month_days_60N, axis = 0)
		
	#-----------------------------------------------------------------------------------------
	
	#Determine the time mean over the months of choice
	v_vel_34S	= np.sum(v_vel_34S * month_days_34S, axis = 0)
	salt_34S	= np.sum(salt_34S * month_days_34S, axis = 0)
	v_vel_26N	= np.sum(v_vel_26N * month_days_26N, axis = 0)
	salt_26N	= np.sum(salt_26N * month_days_26N, axis = 0)
	v_vel_60N	= np.sum(v_vel_60N * month_days_60N, axis = 0)
	salt_60N	= np.sum(salt_60N * month_days_60N, axis = 0)
	
	#-----------------------------------------------------------------------------------------

	filename = directory+'Data/FOV_section_34S/E3SM_data_year_'+str(year_i).zfill(4)+'.nc'

	fh = netcdf.Dataset(filename, 'w')

	fh.createDimension('depth', len(depth))
	fh.createDimension('lon', len(lon_34S))

	fh.createVariable('depth', float, ('depth'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('layer', float, ('depth', 'lon'), zlib=True)
	fh.createVariable('DX', float, ('lon'), zlib=True)
	fh.createVariable('VVEL', float, ('depth', 'lon'), zlib=True)
	fh.createVariable('SALT', float, ('depth', 'lon'), zlib=True)

	fh.variables['depth'].longname		= 'Depth from surface to midpoint of layer'
	fh.variables['layer'].longname		= 'Thickness of layer'
	fh.variables['lon'].longname		= 'Array of longtidues'
	fh.variables['DX'].longname	       = 'x-spacing'
	fh.variables['VVEL'].longname	       = 'Velocity in meridional direction'
	fh.variables['SALT'].longname	       = 'Salinity'

	fh.variables['depth'].units 		= 'm'
	fh.variables['layer'].units 		= 'm'
	fh.variables['lon'].units 		   = 'degrees E'
	fh.variables['DX'].units                = 'm'
	fh.variables['VVEL'].units               = 'm/s'
	fh.variables['SALT'].units               = 'g/kg'

	#Writing data to correct variable	
	fh.variables['depth'][:] 		= depth
	fh.variables['layer'][:] 		= layer_34S
	fh.variables['lon'][:]			= lon_34S
	fh.variables['DX'][:]            	= dx_34S
	fh.variables['VVEL'][:]              	= v_vel_34S
	fh.variables['SALT'][:]              	= salt_34S

	fh.close()

	#-----------------------------------------------------------------------------------------

	filename = directory+'Data/FOV_section_26N/E3SM_data_year_'+str(year_i).zfill(4)+'.nc'

	fh = netcdf.Dataset(filename, 'w')

	fh.createDimension('depth', len(depth))
	fh.createDimension('lon', len(lon_26N))

	fh.createVariable('depth', float, ('depth'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('layer', float, ('depth', 'lon'), zlib=True)
	fh.createVariable('DX', float, ('lon'), zlib=True)
	fh.createVariable('VVEL', float, ('depth', 'lon'), zlib=True)
	fh.createVariable('SALT', float, ('depth', 'lon'), zlib=True)

	fh.variables['depth'].longname		= 'Depth from surface to midpoint of layer'
	fh.variables['layer'].longname		= 'Thickness of layer'
	fh.variables['lon'].longname		= 'Array of longtidues'
	fh.variables['DX'].longname	       = 'x-spacing'
	fh.variables['VVEL'].longname	       = 'Velocity in meridional direction'
	fh.variables['SALT'].longname	       = 'Salinity'

	fh.variables['depth'].units 		= 'm'
	fh.variables['layer'].units 		= 'm'
	fh.variables['lon'].units 		   = 'degrees E'
	fh.variables['DX'].units                = 'm'
	fh.variables['VVEL'].units               = 'm/s'
	fh.variables['SALT'].units               = 'g/kg'

	#Writing data to correct variable	
	fh.variables['depth'][:] 		= depth
	fh.variables['layer'][:] 		= layer_26N
	fh.variables['lon'][:]			= lon_26N
	fh.variables['DX'][:]            	= dx_26N
	fh.variables['VVEL'][:]              	= v_vel_26N
	fh.variables['SALT'][:]              	= salt_26N

	fh.close()

	#-----------------------------------------------------------------------------------------

	filename = directory+'Data/FOV_section_60N/E3SM_data_year_'+str(year_i).zfill(4)+'.nc'

	fh = netcdf.Dataset(filename, 'w')

	fh.createDimension('depth', len(depth))
	fh.createDimension('lon', len(lon_60N))

	fh.createVariable('depth', float, ('depth'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('layer', float, ('depth', 'lon'), zlib=True)
	fh.createVariable('DX', float, ('lon'), zlib=True)
	fh.createVariable('VVEL', float, ('depth', 'lon'), zlib=True)
	fh.createVariable('SALT', float, ('depth', 'lon'), zlib=True)

	fh.variables['depth'].longname		= 'Depth from surface to midpoint of layer'
	fh.variables['layer'].longname		= 'Thickness of layer'
	fh.variables['lon'].longname		= 'Array of longtidues'
	fh.variables['DX'].longname	       = 'x-spacing'
	fh.variables['VVEL'].longname	       = 'Velocity in meridional direction'
	fh.variables['SALT'].longname	       = 'Salinity'

	fh.variables['depth'].units 		= 'm'
	fh.variables['layer'].units 		= 'm'
	fh.variables['lon'].units 		   = 'degrees E'
	fh.variables['DX'].units                = 'm'
	fh.variables['VVEL'].units               = 'm/s'
	fh.variables['SALT'].units               = 'g/kg'

	#Writing data to correct variable	
	fh.variables['depth'][:] 		= depth
	fh.variables['layer'][:] 		= layer_60N
	fh.variables['lon'][:]			= lon_60N
	fh.variables['DX'][:]            	= dx_60N
	fh.variables['VVEL'][:]              	= v_vel_60N
	fh.variables['SALT'][:]              	= salt_60N

	fh.close()	

