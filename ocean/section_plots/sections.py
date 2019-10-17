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
titleTxt = 'EC60to30, MPAS-Ocean stand alone, with Redi mixing on'
varNames = ['temperature','salinity','potentialDensity','relativeSlopeTopOfCell','relativeSlopeTaperingCell','tracer1','tracer2']
landValue = [-0.1,-0.1,1027, -0.1, -0.1, -0.1,-0.1]
meshFile = '../init.nc'

lonRequest = [-150.0, -110]
latSpan = [-70,10]
layerSpan = [0,35]
layer = [0, 22]
iTime = 4

# load mesh variables and time
ncfile1 = Dataset('debugTracersLatLon.nc','r')
xtime = ncfile1.variables['xtime']
lat = ncfile1.variables['lat']
lon = ncfile1.variables['lon']
ncMeshFile = Dataset(meshFile,'r')
refBottomDepth = ncMeshFile.variables['refBottomDepth']

fig = plt.gcf()
fig.set_size_inches(18.0,16.0)
nRows=len(varNames)
nCols=4
nLayerCols=2
iLat0 = np.where(lat[:]>latSpan[0])[0][0]
iLat1 = np.where(lat[:]>latSpan[1])[0][0]

plt.clf()
plt.subplot2grid((nRows+1, nCols), (0,0), colspan=3)
plt.text(.1,.6,titleTxt, fontsize=14, fontweight='bold', verticalalignment='bottom', horizontalalignment='left')
plt.text(.1,.4,'time='+str(xtime[2,:],'utf-8'), fontsize=12, fontweight='normal', verticalalignment='bottom', horizontalalignment='left')
plt.gca().axis('off')
for iRow in range(nRows):
    var = ncfile1.variables[varNames[iRow]]
    for iCol in range(nCols):
        plt.subplot2grid((nRows+1, nCols), (iRow+1,iCol) )
        plt.subplot(nRows+1, nCols, (iRow+1)*nCols+iCol+1)
        iLon = np.where(lon[:]>lonRequest[iCol-nLayerCols])[0][0]

        if iCol<nLayerCols:
            tmp = var[iTime,layer[iCol],:,:]
            tmp2 = np.where(tmp>-1e20,tmp,landValue[iRow])
            ax = plt.imshow(tmp2,extent=[-180,180,-90,90],aspect=1.0)
            plt.title(varNames[iRow]+', layer '+str(layer[iCol]))
            plt.gca().invert_yaxis()
            if varNames[iRow]=='salinity':
                plt.clim(33,37)
            if varNames[iRow]=='potentialDensity':
                plt.clim(1023, 1028)

            if iRow == nRows-1:
                plt.xlabel('longitude')
            else:
                plt.gca().get_xaxis().set_visible(False)

        else:
            tmp = var[iTime,layerSpan[0]:layerSpan[1],iLat0:iLat1,iLon]
            tmp2 = np.where(tmp>-1e20,tmp,landValue[iRow])
            ax = plt.imshow(tmp2,extent=[lat[iLat0],lat[iLat1],layerSpan[1],layerSpan[0]])
            plt.title(varNames[iRow]+', lon='+str(lon[iLon]))
            if iRow == nRows-1:
                plt.xlabel('latitude')

        plt.set_cmap('gist_ncar')
        #plt.clim(-0.5,2.5)
        plt.colorbar()


#plt.colorbar()
plt.tight_layout()
plt.savefig('variables.png')

#plt.clf()
#ax = plt.imshow(tracer2[0,0,:,:]) 
#plt.gca().invert_yaxis()
#plt.clim(-0.5,2.5)
#plt.set_cmap('gist_ncar')
#plt.colorbar()
#plt.savefig('fig3.png')

