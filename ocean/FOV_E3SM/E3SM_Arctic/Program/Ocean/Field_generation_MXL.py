#Generates the mixed layer depth fields on the 0.5x0.5 rectangular grid
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

def Distance(lon1, lat1, lon2, lat2):
	"""Returns distance (m) of two points located at the globe
	coordinates need input in degrees"""

	lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2]) #Convert to radians

	#Haversine formula 
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	a = math.sin(dlat/2.0)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2.0)**2
	c = 2.0 * math.asin(sqrt(a)) 
	r = 6371000.0 # Radius of earth in meters
	
	return c * r #Distance between two points in meter

def GridCellComputer(longitude, latitude):
	"""Determines the area (m^2) per grid cell
	returns 2-D array (lat, lon) with the area per box"""

	#Define empty array for latitude per grid cell and the Area covered by the Ocean
	grid_x 	= np.zeros((len(latitude), len(longitude)))
	grid_y 	= np.zeros((len(latitude), len(longitude)))

	for lat_i in range(len(latitude)):

		#Determining zonal length (m), is latitude dependent, therefore, take middle of grid cell
		length_zonal_grid      	= Distance(0.0, latitude[lat_i], 0.08333206, latitude[lat_i]) 
		#Determining meriodinal length (m), is longitude independent
		length_meridional_grid	= Distance(0.0, 0.0, 0.0, 0.08333206)	

		grid_x[lat_i] 	= length_zonal_grid
		grid_y[lat_i] 	= length_meridional_grid

	return grid_x, grid_y
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

lat_min		= -71.25
lat_max		= 1.25
lon_min		= -70.25
lon_max		= 20.25
depth_level	= 500
	
files = glob.glob(directory_data+'E3SMv2.1B60to10rA02.mpaso.hist.am.timeSeriesStatsMonthly.*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time            = np.zeros(len(files))

for year_i in range(len(files)):
        date  = files[year_i][-13:-3]
        year  = int(date[0:4])
        month = int(date[5:7])

        time[year_i] = year + (month-1) / 12.0
      
 
#-----------------------------------------------------------------------------------------

fh		        = netcdf.Dataset(directory+'Data/Layer_grid.nc', 'r')

lon 		    = fh.variables['lon'][:]
lat		        = fh.variables['lat'][:]

lon_min_index	= (np.abs(lon - lon_min)).argmin()
lon_max_index	= (np.abs(lon - lon_max)).argmin()+1
lat_min_index	= (np.abs(lat - lat_min)).argmin()
lat_max_index	= (np.abs(lat - lat_max)).argmin()+1

lon		        = lon[lon_min_index:lon_max_index]
lat		        = lat[lat_min_index:lat_max_index]

layer		    = fh.variables['timeMonthly_avg_layerThickness'][0] 	#Layer thickness (m)

fh.close()

#Use general depth coordinate
depth		= np.zeros(len(layer)+1)

for depth_i in range(1, len(depth)):
	#Generate the depth boundaries
	depth[depth_i]	= depth[depth_i-1] + layer[depth_i-1, 87, 255]

#Take the mean to find general depth array
depth 	= 0.5 * (depth[1:] + depth[:-1])

depth_index	= (np.abs(depth-depth_level)).argmin()

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
		os.system('ncremap -i '+files_month[file_i]+' -P mpas -m /global/cfs/cdirs/e3sm/diagnostics/mpas_analysis/maps/map_ARRM10to60E2r1_to_0.5x0.5degree_bilinear.nc -T . -O /global/homes/r/rvwesten/E3SM_Arctic/Program/Ocean -o Regrid_month.nc -v timeMonthly_avg_dThreshMLD,timeMonthly_avg_activeTracers_salinity,timeMonthly_avg_activeTracers_temperature')

		fh	= netcdf.Dataset('Regrid_month.nc', 'r')

		mixed	= fh.variables['timeMonthly_avg_dThreshMLD'][0, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Temperature (m)
		salt	= fh.variables['timeMonthly_avg_activeTracers_salinity'][0, depth_index, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Salinity (g/kg)
		temp	= fh.variables['timeMonthly_avg_activeTracers_temperature'][0, depth_index, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Temperature (deg C)
				
		fh.close()
		
		if file_i == 0:
			#Empty array
			temp_month	= ma.masked_all((12, len(lat), len(lon)))
			salt_month	= ma.masked_all((12, len(lat), len(lon)))
			mixed_month	= ma.masked_all((12, len(lat), len(lon)))
		
		#Save the vertical velocity	
		temp_month[file_i]		= temp
		salt_month[file_i]		= salt
		mixed_month[file_i]		= mixed		
	#------------------------------------------------------------------------------

	
	filename 	= directory+'Data/Mixed_layer/E3SM_data_year_'+str(year_i).zfill(4)+'.nc'
	
	fh = netcdf.Dataset(filename, 'w')
  
	fh.createDimension('month', 12)
	fh.createDimension('lat', len(lat))
	fh.createDimension('lon', len(lon))
	  
	fh.createVariable('month', float, ('month'), zlib=True)
	fh.createVariable('lat', float, ('lat'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('TEMP', float, ('month', 'lat', 'lon'), zlib=True)
	fh.createVariable('SALT', float, ('month', 'lat', 'lon'), zlib=True)
	fh.createVariable('MXL', float, ('month', 'lat', 'lon'), zlib=True)

	fh.variables['month'].longname	 	= 'Month'
	fh.variables['lat'].longname	       = 'Array of latitudes'
	fh.variables['lon'].longname		= 'Array of longitudes'
	fh.variables['TEMP'].longname		= 'Potential temperature (500 m depth)'
	fh.variables['SALT'].longname		= 'Salinity (500 m depth)'
	fh.variables['MXL'].longname		= 'Mixed layer depth'
	  
	fh.variables['lat'].units 	        = 'degrees N'
	fh.variables['lon'].units 		     = 'degrees E'
	fh.variables['TEMP'].units        	 = 'degC'
	fh.variables['SALT'].units        	 = 'g/kg'
	fh.variables['MXL'].units        	 = 'm'
 
	#Writing data to correct variable	
	fh.variables['month'][:] 		    = np.arange(1,13)
	fh.variables['lat'][:] 	        	= lat
	fh.variables['lon'][:] 			    = lon
	fh.variables['TEMP'][:]            	= temp_month
	fh.variables['SALT'][:]            	= salt_month
	fh.variables['MXL'][:]            	= mixed_month

	fh.close()
	
