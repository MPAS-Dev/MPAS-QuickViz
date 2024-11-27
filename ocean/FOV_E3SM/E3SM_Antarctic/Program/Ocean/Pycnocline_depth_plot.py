#Program plots the pycnocline depth and Hovmoller diagram

from pylab import *
import netCDF4 as netcdf
import glob, os
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Making pathway to folder with all data
directory 	= '../../Data/'
    
#------------------------------------------------------------------------------
#--------------------------MAIN SCRIPT STARTS HERE-----------------------------
#------------------------------------------------------------------------------

#Only select for 27 years (similar to reanalysis)
time = np.zeros((526-500+1)*12)

for year_i in range(500, 527):
	#Loop over each year
	 
	filename = directory+'Data/Pycnocline_depth/E3SM_data_year_'+str(year_i).zfill(4)+'.nc'
	#------------------------------------------------------------------------------

	fh      = netcdf.Dataset(filename, 'r')
		    
	lon	    = fh.variables['lon'][:]     	
	lat	    = fh.variables['lat'][:]
	depth   = fh.variables['PD_depth'][:]

	fh.close()


	if year_i == 500:
		#Make empty array to save the data
		depth_all	= ma.masked_all((len(time), len(lat), len(lon)))

	for month_i in range(12):
		#Loop over each month
		time[(year_i-500)*12+month_i]      = year_i + month_i / 12.
		depth_all[(year_i-500)*12+month_i] = depth[month_i]
        
   
filename     = directory+'Data/Pycnocline_depth/E3SM_data_year_'+str(500).zfill(4)+'_January.nc'

fh           = netcdf.Dataset(filename, 'r')
		    
lon_plot	 = fh.variables['lon'][:]     	
lat_plot	 = fh.variables['lat'][:]
depth_plot   = fh.variables['PD_depth'][0]

fh.close()
	
#-----------------------------------------------------------------------------------------

#Take the meridional mean
depth_all   = np.mean(depth_all, axis = 1)

fig, ax     = subplots()

CS      = ax.contourf(lon, time, depth_all, levels = np.arange(0, 800.01, 25), extend = 'max', cmap = 'Spectral_r')
cbar    = colorbar(CS, ticks = np.arange(0, 800.01, 100))
cbar.set_label('Pycnocline depth (m)')

fig, ax     = subplots()
    
ax.plot(lon, np.mean(depth_all, axis = 0), '-k', linewidth = 2.0)
ax.set_xlim(-50, 20)
ax.set_ylim(1000, 0)
ax.set_ylabel('Depth (m)')
ax.grid()
ax.set_title('Isopycnal depth, E3SM (500 - 526)')

ax.set_xticks(np.arange(-50, 20.01, 10))
ax.set_xticklabels(['50$^{\circ}$W', '40$^{\circ}$W', '30$^{\circ}$W', '20$^{\circ}$W', '10$^{\circ}$W','0$^{\circ}$', '10$^{\circ}$E', '20$^{\circ}$E'])

ax2 	= fig.add_axes([0.15, 0.50, 0.45, 0.30], projection = ccrs.PlateCarree())

CS      = ax2.contourf(lon_plot, lat_plot, depth_plot, levels = np.arange(0, 800.1, 25), extend = 'max', cmap = 'Spectral_r', transform=ccrs.PlateCarree())
cbar    = colorbar(CS, ticks = np.arange(0, 800.01, 200), fraction=0.021, pad=0.04)
cbar.set_label('Isopycnal depth (m)', fontsize = 8)

ax2.plot([-50, 20], [-30, -30], '--k', transform=ccrs.PlateCarree(), zorder = 100)
ax2.plot([-50, 20], [-25, -25], '--k', transform=ccrs.PlateCarree(), zorder = 100)

gl = ax2.gridlines(draw_labels=False)
ax2.set_extent([-70, 30, -46, 1], ccrs.PlateCarree())
ax2.coastlines('50m')
ax2.add_feature(cfeature.LAND, zorder=0)
ax2.set_title('Isopycnal depth, E3SM (January 500)', fontsize = 10)

#-----------------------------------------------------------------------------------------


#Remove climatology
for month_i in range(12):
	time_index              = np.arange(month_i, len(time), 12)
	depth_mean              = np.mean(depth_all[time_index], axis = 0)
	depth_all[time_index]   = depth_all[time_index] - depth_mean

#-----------------------------------------------------------------------------------------
   
fig, ax     = subplots()

CS      = ax.contourf(lon, time, depth_all, levels = np.arange(-50, 50.01, 5), extend = 'both', cmap = 'RdBu_r')
cbar    = colorbar(CS, ticks = np.arange(-50, 50.01, 10))
cbar.set_label('Pycnocline depth anomaly (m)')

ax.set_xlim(-50, 20)
ax.set_ylim(500, 527-1/12.)
ax.set_ylabel('Model year')

ax.set_xticks(np.arange(-50, 20.01, 10))
ax.set_xticklabels(['50$^{\circ}$W', '40$^{\circ}$W', '30$^{\circ}$W', '20$^{\circ}$W', '10$^{\circ}$W','0$^{\circ}$', '10$^{\circ}$E', '20$^{\circ}$E'])

ax.plot([15, -47.6], [504.9, 512.8], '-k', linewidth = 3.0)
ax.set_title('Isopycnal depth (30$^{\circ}$S - 25$^{\circ}$S), E3SM')


show()
