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
varNames = ['temperature','salinity','potentialDensity','relativeSlopeTopOfCell','relativeSlopeTaperingCell','tracer1','tracer2']
landValue = [-0.1,-0.1,1027, -0.1, -0.1, -0.1,-0.1]

lonRequest = [-150.0, -110]
latSpan = [-70,10]
layerSpan = [0,35]
iLayer = 22
iTime = 4

ncfile1 = Dataset('debugTracersLatLon.nc','r')
lat = ncfile1.variables['lat']
lon = ncfile1.variables['lon']

fig = plt.gcf()
fig.set_size_inches(16.0,16.0)
nRow=len(varNames)
nCol=3
plt.clf()
iLat0 = np.where(lat[:]>latSpan[0])[0][0]
iLat1 = np.where(lat[:]>latSpan[1])[0][0]
print('lat',iLat0,iLat1)
for iRow in range(nRow):
    var = ncfile1.variables[varNames[iRow]]
    for iCol in range(nCol):
        iLon = np.where(lon[:]>lonRequest[iCol-1])[0][0]

        #print('counter',2*iRow+iCol+1)
        plt.subplot(nRow, nCol, iRow*nCol+iCol+1)
        if iCol==0:
            tmp = var[iTime,iLayer,:,:]
            tmp2 = np.where(tmp>-1e20,tmp,landValue[iRow])
            ax = plt.imshow(tmp2,extent=[-180,180,-90,90],aspect=0.8)
        else:
            tmp = var[iTime,layerSpan[0]:layerSpan[1],iLat0:iLat1,iLon]
            tmp2 = np.where(tmp>-1e20,tmp,landValue[iRow])
            ax = plt.imshow(tmp2,extent=[lat[iLat0],lat[iLat1],layerSpan[1],layerSpan[0]])
        plt.set_cmap('gist_ncar')
        #plt.clim(-0.5,2.5)
        plt.colorbar()
        if iCol == 0:
            plt.title(varNames[iRow]+', layer '+str(iLayer))
            plt.gca().invert_yaxis()
        else:
            plt.title(varNames[iRow]+', lon='+str(lon[iLon]))

        if iRow == nRow-1:
            if iCol == 0:
                plt.xlabel('longitude')
            else:
                plt.xlabel('latitude')
        else:
            plt.gca().get_xaxis().set_visible(False)


#plt.colorbar()
plt.savefig('variables.png')

#plt.clf()
#ax = plt.imshow(tracer2[0,0,:,:]) 
#plt.gca().invert_yaxis()
#plt.clim(-0.5,2.5)
#plt.set_cmap('gist_ncar')
#plt.colorbar()
#plt.savefig('fig3.png')

