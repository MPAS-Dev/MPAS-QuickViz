'''
sections.py
Compute zonal or meridional sections of variables on an mpas mesh
'''

#import numpy as np
#import matplotlib.pyplot as plt
#import time
#import xarray
#import os
#import subprocess

# User specified files and variables

# plot settings
axis_font = {'fontname':'Arial', 'size':'18'}    
title_font = {'fontname':'Arial', 'size':'32', 'color':'black', 'weight':'normal'}
matplotlib.rc('xtick', labelsize=28)
matplotlib.rc('ytick', labelsize=28)

#add nco stuff here
# Fields to process for each run, define baseline and experiment here and workdir and variables to extract

# What I actually did on the command line:
# mkdir regridded
# ncks -d Time,0,100,10 output/debugTracer.0001-01-01_00.00.00.nc regridded/temp.nc
# ncks -v latCell,lonCell -A init.nc regridded/temp.nc
# ncks -d Time,0,2 -v potentialDensity output/output.0001-01-01_00.00.00.nc regridded/tempDensity.nc
# ncks -v latCell,lonCell -A init.nc regridded/tempDensity.nc
# cd regridded
# ln -isf /usr/projects/climate/mpeterse/repos/APrime_Files/mapping/maps/* .
# ncremap -i temp.nc -o debugTracersLatLon.nc -P mpas -m map_oEC60to30v3_TO_0.5x0.5degree_blin.nc -R "--rgr lat_nm=latCell --rgr lon_nm=lonCell --rgr lat_nm_out=lat --rgr lon_nm_out=lon" -C
# ncremap -i tempDensity.nc -o densityLatLon.nc -P mpas -m map_oEC60to30v3_TO_0.5x0.5degree_blin.nc -R "--rgr lat_nm=latCell --rgr lon_nm=lonCell --rgr lat_nm_out=lat --rgr lon_nm_out=lon" -C

# I didn't get the density remapper to work.

# Mark easy plot section

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

lonRequest = [30.0]

ncfile1 = Dataset('debugLatLon.nc','r')
lat = ncfile1.variables['lat']
lon = ncfile1.variables['lon']
tracer1 = ncfile1.variables['tracer1']
tracer2 = ncfile1.variables['tracer2']
iLon = np.where(lon[:]>lonRequest)[0][0]


iLon = np.where(lon>30.0)[0]
plt.imshow(cellWidth)
plt.gca().invert_yaxis()
plt.colorbar
plt.xlabel('longitude')
plt.ylabel('latitude')
plt.title('cell width, km')
plt.savefig('cellWidth.png')

