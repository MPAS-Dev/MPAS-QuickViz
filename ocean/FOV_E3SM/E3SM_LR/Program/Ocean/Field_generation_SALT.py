#Generates the vertically-integrated salinity fields on the 0.5x0.5 rectangular grid
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

directory_data  	= '/global/cfs/cdirs/m4259/E3SMv2_1/20220715.submeso.piControl.ne30pg2_EC30to60E2r2.chrysalis/ocn/hist/'
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
		length_zonal_grid      	= Distance(0.0, latitude[lat_i], np.mean(np.diff(longitude)), latitude[lat_i]) 
		#Determining meriodinal length (m), is longitude independent
		length_meridional_grid	= Distance(0.0, 0.0, 0.0, np.mean(np.diff(latitude)))	

		grid_x[lat_i] 	= length_zonal_grid
		grid_y[lat_i] 	= length_meridional_grid

	return grid_x, grid_y, grid_x * grid_y

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

lon_min		= -110
lon_max		= 143
lat_min		= -80
lat_max		= 25.5
depth_min	= 0
depth_max	= 100


files = glob.glob(directory_data+'20220715.submeso.piControl.ne30pg2_EC30to60E2r2.chrysalis.mpaso.hist.am.timeSeriesStatsMonthly.*.nc')
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

fh	           = netcdf.Dataset(directory+'Data/Layer_grid.nc', 'r')

lon 			= fh.variables['lon'][:]
lat			    = fh.variables['lat'][:]

grid_x, grid_y, area 	= GridCellComputer(lon, lat)

lon_min_index	= (np.abs(lon - lon_min)).argmin()
lon_max_index	= (np.abs(lon - lon_max)).argmin()+1
lat_min_index	= (np.abs(lat - lat_min)).argmin()
lat_max_index	= (np.abs(lat - lat_max)).argmin()+1

lon		        = lon[lon_min_index:lon_max_index]
lat		        = lat[lat_min_index:lat_max_index]
area		    = area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
layer	        = fh.variables['timeMonthly_avg_layerThickness'][0, :, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Layer thickness (m)

fh.close()

#Use general depth coordinate
depth		= np.zeros(len(layer)+1)

for depth_i in range(1, len(depth)):
	#Generate the depth boundaries
	depth[depth_i]	= depth[depth_i-1] + layer[depth_i-1, 87, 255]

#Take the mean to find general depth array
depth 	= 0.5 * (depth[1:] + depth[:-1])


depth_min_index	= (fabs(depth_min - depth)).argmin()
depth_max_index	= (fabs(depth_max - depth)).argmin() + 1
depth		= depth[depth_min_index:depth_max_index]
layer		= layer[depth_min_index:depth_max_index]

for lat_i in range(len(lat)):
	for lon_i in range(len(lon)):
		#Get all the layers which have a maximum depth below given range
		if np.sum(layer[:, lat_i, lon_i]) > depth_max:
			#Adjust the last layer
			layer[-1, lat_i, lon_i]  -= (np.sum(layer[:, lat_i, lon_i]) - depth_max)
       
#Get the total vertical extent for each layer
total_layer	= np.sum(layer, axis = 0)
volume		= total_layer * area
area		= ma.masked_array(area, mask = volume.mask)

for depth_i in range(len(depth)):
	#Normalise the field by its vertical extent
	layer[depth_i]	= layer[depth_i] / total_layer
				
#-----------------------------------------------------------------------------------------		
		
time_year   = ma.masked_all(int(len(time)/12))

for year_i in range(int(np.min(time)), int(np.min(time))+len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i - int(np.min(time))] = year_i
	files_month = glob.glob(directory_data+'20220715.submeso.piControl.ne30pg2_EC30to60E2r2.chrysalis.mpaso.hist.am.timeSeriesStatsMonthly.'+str(year_i).zfill(4)+'-*.nc')
	files_month.sort()

	for file_i in range(len(files_month)):
		#Loop over each month
		os.system('ncremap -i '+files_month[file_i]+' -P mpas -m /global/cfs/cdirs/e3sm/diagnostics/mpas_analysis/maps/map_EC30to60E2r2_to_0.5x0.5degree_bilinear.nc -T . -O /global/homes/r/rvwesten/E3SM_LR/Program/Ocean -o Regrid_month.nc -v timeMonthly_avg_activeTracers_salinity')

		fh	= netcdf.Dataset('Regrid_month.nc', 'r')

		salt	= fh.variables['timeMonthly_avg_activeTracers_salinity'][0, depth_min_index:depth_max_index, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Salinity (g/kg)

		fh.close()
		
		salt	= ma.masked_where(layer <= 0.0, salt)
		
		if file_i == 0:
			#Empty array
			salt_depth	= ma.masked_all((12, len(lat), len(lon)))
			
		#Get the vertical depth averaged salinity
		salt_depth[file_i]	= np.sum(salt * layer, axis = 0)
		
	#------------------------------------------------------------------------------
	#Now convert to yearly averages
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(lat), len(lon)))
		
	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]
		
	#Now set mask
	month_days_all	= ma.masked_array(month_days_all, mask = salt_depth.mask)
		
	#Normalise the data
	month_days_all	= month_days_all / np.sum(month_days_all, axis = 0)

	#Determine the time mean over the months of choice
	salt_depth	= np.sum(salt_depth * month_days_all, axis = 0)
			
	#-----------------------------------------------------------------------------------------
	
	filename 	= directory+'Data/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m/E3SM_data_year_'+str(year_i).zfill(4)+'.nc'
	
	fh = netcdf.Dataset(filename, 'w')

	fh.createDimension('lon', len(lon))
	fh.createDimension('lat', len(lat))

	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('lat', float, ('lat'), zlib=True)
	fh.createVariable('AREA', float, ('lat', 'lon'), zlib=True)
	fh.createVariable('VOLUME', float, ('lat', 'lon'), zlib=True)
	fh.createVariable('SALT', float, ('lat', 'lon'), zlib=True)

	fh.variables['lon'].longname	 	= 'Array of T-longtidues'
	fh.variables['lat'].longname	       = 'Array of T-latitudes'
	fh.variables['AREA'].longname	       = 'Area of T cells'
	fh.variables['VOLUME'].longname       = 'Volume of T cells'
	fh.variables['SALT'].longname	       = 'Depth-averaged salinity'

	fh.variables['lon'].units 		= 'degrees E'
	fh.variables['lat'].units 	        = 'degrees N'
	fh.variables['AREA'].units             = 'm^2'
	fh.variables['VOLUME'].units           = 'm^3'
	fh.variables['SALT'].units              = 'g/kg'

	#Writing data to correct variable	
	fh.variables['lon'][:] 			= lon
	fh.variables['lat'][:] 	     		= lat
	fh.variables['AREA'][:]           	= area
	fh.variables['VOLUME'][:]           	= volume
	fh.variables['SALT'][:]             	= salt_depth

	fh.close()


		

