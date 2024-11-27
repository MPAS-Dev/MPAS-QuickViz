#Generates the 1027 kg/m^3 pycnocine depth fields on the 0.5x0.5 rectangular grid
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
lat_min		= -30
lat_max		= -25

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
		os.system('ncremap -i '+files_month[file_i]+' -P mpas -m /global/cfs/cdirs/m4259/mapping_files/map_SOwISC12to60E2r4_to_0.5x0.5degree_bilinear.nc -T . -O /global/homes/r/rvwesten/E3SM/Program/Ocean -o Regrid_month.nc -v timeMonthly_avg_activeTracers_salinity,timeMonthly_avg_activeTracers_temperature')

		fh	= netcdf.Dataset('Regrid_month.nc', 'r')

		salt	= fh.variables['timeMonthly_avg_activeTracers_salinity'][0, :, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Salinity (g/kg)
		temp	= fh.variables['timeMonthly_avg_activeTracers_temperature'][0, :, lat_min_index:lat_max_index, lon_min_index:lon_max_index] 	#Temperature (deg C)
				
		fh.close()
		
		salt	= ma.masked_where(layer <= 0.0, salt)
		temp	= ma.masked_where(layer <= 0.0, temp)
		dens	= RHO_0(temp, salt)

		if file_i == 0:
			#Empty array
			dens_depth	= ma.masked_all((12, len(lat), len(lon)))
			
		for lat_i in range(len(lat)):

			fig, ax = subplots()
			
			CS    = ax.contour(lon, depth, dens[:, lat_i], levels = [1027])

			close()

			for item in CS.collections:
				for i in item.get_paths():
					v = i.vertices
					x = v[:, 0]
					y = v[:, 1]


					for lon_i in range(len(lon)):
						#Check the location for each position
						if np.all(lon[lon_i] < x) or np.all(lon[lon_i] > x) or lon[lon_i] > 14.95:
							continue

						x_index         = (np.abs(lon[lon_i] - x)).argmin()

						if np.abs(lon[lon_i] - x[x_index]) > 0.1:
							#Too far apart
							continue


						dens_depth[file_i, lat_i, lon_i]   = y[x_index]
					
	#-----------------------------------------------------------------------------------------
	
	filename 	= directory+'Data/Pycnocline_depth/E3SM_data_year_'+str(year_i).zfill(4)+'.nc'

	fh = netcdf.Dataset(filename, 'w')

	fh.createDimension('month', 12)     
	fh.createDimension('lat', len(lat))
	fh.createDimension('lon', len(lon))

	fh.createVariable('month', float, ('month'), zlib=True)
	fh.createVariable('lat', float, ('lat'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('PD_depth', float, ('month', 'lat', 'lon'), zlib=True)

	fh.variables['lat'].longname	       = 'Array of latitudes'
	fh.variables['lat'].longname	       = 'Array of longitude'
	fh.variables['PD_depth'].longname   = 'Potetial density (1027) depth'

	fh.variables['lat'].units 	        = 'degrees N'
	fh.variables['lon'].units             = 'degrees E'
	fh.variables['PD_depth'].units         = 'm'


	#Writing data to correct variable	
	fh.variables['month'][:] 		= np.arange(12)+1
	fh.variables['lat'][:] 	        = lat
	fh.variables['lon'][:]        = lon
	fh.variables['PD_depth'][:]    = dens_depth


	fh.close()
	
	
	
