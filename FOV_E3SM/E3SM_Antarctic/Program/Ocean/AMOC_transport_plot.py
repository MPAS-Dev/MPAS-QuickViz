#Program plots the AMOC strength

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory 	= '../../Data/'
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 1000


#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'_m.nc', 'r')
	
time		= fh.variables['time'][:]     	  	
transport	= fh.variables['Transport'][:]    	

fh.close()


fig, ax	= subplots()

ax.fill_between([-100, 2500], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')


ax.plot(time, transport, '-k', linewidth = 2.0)
ax.set_xlim(500, 600)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (sv)')
ax.set_xticks([500, 520, 540, 560, 580, 600])
ax.grid()

ax.set_title('AMOC strength, E3SM Antarctic')

show()
