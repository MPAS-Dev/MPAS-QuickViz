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
# axis_font = {'fontname':'Arial', 'size':'18'}    
# title_font = {'fontname':'Arial', 'size':'32', 'color':'black', 'weight':'normal'}
# matplotlib.rc('xtick', labelsize=28)
# matplotlib.rc('ytick', labelsize=28)

#add nco stuff here
# Fields to process for each run, define baseline and experiment here and workdir and variables to extract

'''
# What I actually did on the command line:
mkdir regridded
ncks -O output/debugTracer.0001-01-01_00.00.00.nc regridded/temp.nc
ncks -v latCell,lonCell -A init.nc regridded/temp.nc
cd regridded
ln -isf /usr/projects/climate/mpeterse/repos/APrime_Files/mapping/maps/* .
ncremap -i temp.nc -o debugTracersLatLon.nc -P mpas -m map_oEC60to30v3_TO_0.5x0.5degree_blin.nc -R "--rgr lat_nm=latCell --rgr lon_nm=lonCell --rgr lat_nm_out=lat --rgr lon_nm_out=lon" -C

ncks -d Time,0,2 -v potentialDensity output/output.0001-01-01_00.00.00.nc regridded/tempDensity.nc
ncks -v latCell,lonCell -A init.nc regridded/tempDensity.nc
ncremap -i tempDensity.nc -o densityLatLon.nc -P mpas -m map_oEC60to30v3_TO_0.5x0.5degree_blin.nc -R "--rgr lat_nm=latCell --rgr lon_nm=lonCell --rgr lat_nm_out=lat --rgr lon_nm_out=lon" -C
# I didn't get the density remapper to work.
'''

# Mark easy plot section

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

# Input arguments
varNames = ['potentialDensity','relativeSlopeTopOfCell','relativeSlopeTaperingCell','k33','tracer1','tracer2']
landValue = [1027, -0.1, -0.1, 0.0,-0.1,-0.1]

lonRequest = [-150.0, -120, 30]
latSpan = [-70,10]
iTime=4

ncfile1 = Dataset('debugTracersLatLon.nc','r')
lat = ncfile1.variables['lat']
lon = ncfile1.variables['lon']

fig = plt.gcf()
fig.set_size_inches(8.0,16.0)
nRow=len(varNames)
nCol=1
plt.clf()
iLon = np.where(lon[:]>lonRequest[0])[0][0]
iLat0 = np.where(lat[:]>latSpan[0])[0][0]
iLat1 = np.where(lat[:]>latSpan[1])[0][0]
print('lat',iLat0,iLat1)
for iRow in range(nRow):
    var = ncfile1.variables[varNames[iRow]]
    for iCol in range(nCol):
        #print('counter',2*iRow+iCol+1)
        plt.subplot(nRow, nCol, iRow*nCol+iCol+1)
        tmp = var[iTime,0:35,iLat0:iLat1,iLon]
        tmp2 = np.where(tmp>-1e20,tmp,landValue[iRow])
        ax = plt.imshow(tmp2,extent=[lat[iLat0],lat[iLat1],30,0])
        plt.set_cmap('gist_ncar')
        #plt.clim(-0.5,2.5)
        plt.colorbar()
        if iCol == 0:
            plt.ylabel(varNames[iRow])
        if iCol == nCol:
            plt.ylabel('level')
        if iRow == 0:
            plt.title('lon=' + str(lon[iLon]))
        if iRow == nRow-1:
            plt.xlabel('latitude')

#plt.colorbar()
plt.savefig('variables.png')

#plt.clf()
#ax = plt.imshow(tracer2[0,0,:,:]) 
#plt.gca().invert_yaxis()
#plt.clim(-0.5,2.5)
#plt.set_cmap('gist_ncar')
#plt.colorbar()
#plt.savefig('fig3.png')

