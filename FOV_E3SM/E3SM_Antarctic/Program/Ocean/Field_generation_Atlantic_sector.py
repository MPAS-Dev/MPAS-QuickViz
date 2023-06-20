#Generates the Atlantic sector fields on the 0.5x0.5 rectangular grid
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

directory_data  	= '/pscratch/sd/a/abarthel/data/E3SMv2.1/20221116.CRYO1950.ne30pg2_SOwISC12to60E2r4.N2Dependent.submeso/archive/ocn/hist/'
directory 	        = '../../Data/'

def RHO_0(T, S):
    #Reference density which is not pressure dependent

	rho = (999.842594 + 6.793952 * 10**(-2.0) * T - 9.095290 * 10**(-3.0)*T**(2.0) + 1.001685 * 10**(-4.0) * T**(3.0) - 1.120083 * 10**(-6.0) * T**(4.0) + 6.536332 * 10**(-9.0) * T**(5.0)+ (8.25917 * 10**(-1.0) - 4.4490 * 10**(-3.0) * T + 1.0485 * 10**(-4.0) * T**(2.0) - 1.2580 * 10**(-6.0) * T**(3.0) + 3.315 * 10**(-9.0) * T**(4.0)) * S+ (- 6.33761 * 10**(-3.0) + 2.8441 * 10**(-4.0) * T - 1.6871 * 10**(-5.0) * T**(2.0) + 2.83258 * 10**(-7.0) * T**(3.0)) * S**(3.0/2.0)+ (5.4705 * 10**(-4.0) - 1.97975 * 10**(-5.0) * T + 1.6641 * 10**(-6.0) * T**(2.0) - 3.1203 * 10**(-8.0) * T**(3.0)) * S**(2.0) )

	return rho
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

lon_min		= -50
lon_max		= 20
lat_min		= -71
lat_max		= 1

files = glob.glob(directory_data+'20221116.CRYO1950.ne30pg2_SOwISC12to60E2r4.N2Dependent.submeso.chrysalis.mpaso.hist.am.timeSeriesStatsMonthly.*.nc')
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

layer	        = fh.variables['timeMonthly_avg_layerThickness'][0] 	#Layer thickness (m)

fh.close()

#Use general depth coordinate
depth		= np.zeros(len(layer)+1)

for depth_i in range(1, len(depth)):
	#Generate the depth boundaries
	depth[depth_i]	= depth[depth_i-1] + layer[depth_i-1, 87, 255]

#Take the mean to find general depth array
depth 	= 0.5 * (depth[1:] + depth[:-1])
layer	= layer[:, lat_min_index:lat_max_index, lon_min_index:lon_max_index]
				
#-----------------------------------------------------------------------------------------		
		
time_year   = ma.masked_all(int(len(time)/12))

for year_i in range(int(np.min(time)), int(np.min(time))+len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i - int(np.min(time))] = year_i
	files_month = glob.glob(directory_data+'20221116.CRYO1950.ne30pg2_SOwISC12to60E2r4.N2Dependent.submeso.chrysalis.mpaso.hist.am.timeSeriesStatsMonthly.'+str(year_i).zfill(4)+'-*.nc')
	files_month.sort()

	for file_i in range(len(files_month)):
		#Loop over each month
		os.system('ncremap -i '+files_month[file_i]+' -P mpas -m /global/cfs/cdirs/m4259/mapping_files/map_SOwISC12to60E2r4_to_0.5x0.5degree_bilinear.nc -T . -O /global/homes/r/rvwesten/E3SM/Program/Ocean -o Regrid_month.nc -v timeMonthly_avg_activeTracers_salinity,timeMonthly_avg_activeTracers_temperature,timeMonthly_avg_velocityZonal')

		fh	= netcdf.Dataset('Regrid_month.nc', 'r')

		salt	= fh.variables['timeMonthly_avg_activeTracers_salinity'][0, :, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Salinity (g/kg)
		temp	= fh.variables['timeMonthly_avg_activeTracers_temperature'][0, :, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Temperature (deg C)
		u_vel	= fh.variables['timeMonthly_avg_velocityZonal'][0, :, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Zonal velocity (m/s)
				
		fh.close()
		
		salt	= ma.masked_where(layer <= 0.0, salt)
		temp	= ma.masked_where(layer <= 0.0, temp)
		u_vel	= ma.masked_where(layer <= 0.0, u_vel)
		dens	= RHO_0(temp, salt)
		
		if file_i == 0:
			#Empty array
			salt_depth	= ma.masked_all((12, len(depth), len(lat)))
			temp_depth	= ma.masked_all((12, len(depth), len(lat)))
			u_vel_depth	= ma.masked_all((12, len(depth), len(lat)))
			dens_depth	= ma.masked_all((12, len(depth), len(lat)))
						
		#Get the zonal mean
		salt_depth[file_i]	= np.mean(salt, axis = 2)
		temp_depth[file_i]	= np.mean(temp, axis = 2)
		u_vel_depth[file_i]	= np.mean(u_vel, axis = 2)
		dens_depth[file_i]	= np.mean(dens, axis = 2)		
	#------------------------------------------------------------------------------
	#Now convert to yearly averages
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth), len(lat)))
		
	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]
		
	#Now set mask
	month_days_all	= ma.masked_array(month_days_all, mask = salt_depth.mask)
		
	#Normalise the data
	month_days_all	= month_days_all / np.sum(month_days_all, axis = 0)

	#Determine the time mean over the months of choice
	salt_depth	= np.sum(salt_depth * month_days_all, axis = 0)
	temp_depth	= np.sum(temp_depth * month_days_all, axis = 0)
	u_vel_depth	= np.sum(u_vel_depth * month_days_all, axis = 0)
	dens_depth	= np.sum(dens_depth * month_days_all, axis = 0)			
	#-----------------------------------------------------------------------------------------
	
	filename 	= directory+'Data/Atlantic_sector/E3SM_data_year_'+str(year_i).zfill(4)+'.nc'
	
	fh = netcdf.Dataset(filename, 'w')
  
	fh.createDimension('depth', len(depth))
	fh.createDimension('lat', len(lat))
	  
	fh.createVariable('depth', float, ('depth'), zlib=True)
	fh.createVariable('lat', float, ('lat'), zlib=True)
	fh.createVariable('SALT', float, ('depth', 'lat'), zlib=True)
	fh.createVariable('TEMP', float, ('depth', 'lat'), zlib=True)
	fh.createVariable('UVEL', float, ('depth', 'lat'), zlib=True)
	fh.createVariable('POT_DENS', float, ('depth', 'lat'), zlib=True)

	fh.variables['depth'].longname	 	= 'Mid-level depth'
	fh.variables['lat'].longname	       = 'Array of latitudes'
	fh.variables['SALT'].longname	       = 'Zonally-averaged salinity'
	fh.variables['TEMP'].longname	       = 'Zonally-averaged potential temperature'
	fh.variables['UVEL'].longname	       = 'Zonally-averaged zonal velocity'
	fh.variables['POT_DENS'].longname	   = 'Zonally-averaged potential density'
	  
	fh.variables['depth'].units 		= 'm'
	fh.variables['lat'].units 	        = 'degrees N'
	fh.variables['SALT'].units             = 'g/kg'
	fh.variables['TEMP'].units             = 'deg C'
	fh.variables['UVEL'].units             = 'm/s'
	fh.variables['POT_DENS'].units         = 'kg/m^3'

	  
	#Writing data to correct variable	
	fh.variables['depth'][:] 		= depth
	fh.variables['lat'][:] 	        = lat
	fh.variables['SALT'][:]        = salt_depth
	fh.variables['TEMP'][:]        = temp_depth
	fh.variables['UVEL'][:]        = u_vel_depth
	fh.variables['POT_DENS'][:]    = dens_depth


	fh.close()
		

