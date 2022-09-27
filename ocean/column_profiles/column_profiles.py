#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s/'
#runDir = '/lustre/scratch5/turquoise/mpeterse/runs/220922_soma_find_noise/ocean/soma/32km/long/'
simName = 's02n'
iTime = 3
fileName = '/output.nc'
meshName = '/init.nc'
#domainName = 'soma'
domainName = 'EC60to30'

deg2rad = 3.14159/180.0
rad2deg = 180.0/3.14159
if domainName == 'EC60to30':
# global locations
# Atlantic
    lonMin = 360 -63
    lonMax = 360 -60
# Pacific
    #lonMin = 175
    #lonMax = 177
    latMin =  30
    latMax =  33
# Soma locations
elif domainName == 'soma':
    lonMin = 0.193*rad2deg
    lonMax = 0.210*rad2deg
    latMin = 0.477*rad2deg
    latMax = 0.490*rad2deg
else:
    print('incorrect domain')

#import os
#import glob
import matplotlib as mpl
from datetime import date
#mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as cols
#from matplotlib.pyplot import cm
#from matplotlib.colors import BoundaryNorm
#import cmocean
import xarray as xr
from netCDF4 import Dataset
#from mpas_analysis.shared.io.utility import decode_strings
#import gsw

mesh = xr.open_dataset(runDir+simName+meshName)
data = xr.open_dataset(runDir+simName+fileName)

z = mesh.refBottomDepth.values
nVertLevels = len(z)
nCells = data.dims['nCells']

latCell = mesh.variables['latCell']
lonCell = mesh.variables['lonCell']
maxLevelCell = mesh.variables['maxLevelCell']

cellList = np.where(np.logical_and(np.logical_and(np.logical_and(
    latCell>latMin*deg2rad, latCell<latMax*deg2rad), 
    lonCell>lonMin*deg2rad), lonCell<lonMax*deg2rad))[0]

fig = plt.figure(figsize=(20,12))
varNames = ['divergence','vertVelocityTop','vertTransportVelocityTop','temperature', 'salinity','density','pressure','zMid']
varNames = ['temperature', 'salinity']
varNames = ['divergence','vertVelocityTop','temperature', 'salinity']
for j in range(len(varNames)):
    plt.subplot(2,2,j+1)
    var = data.variables[varNames[j]]
    for i in range(len(cellList)):
        iCell = int(cellList[i])
        k = int(maxLevelCell[iCell])
        varData = var[iTime,iCell,0:k]
        #varData = np.where(varData>-1e20,varData,np.NAN)
        plt.plot(varData,np.arange(1,k+1),label='cell '+str(iCell))
    plt.gca().invert_yaxis()
    plt.title(varNames[j])
    plt.grid()
    plt.ylabel('vertical index, k')
    #plt.legend()

# add information to bottom of figure
today = date.today()
xtime = data.variables['xtime']
import codecs
strXTime = str(codecs.decode(xtime[iTime].values))[2:19]
lonMin = (lonMin+180.0)%360.0 - 180.0
plt.figtext(0.1,0.05,domainName+' '+simName+' '+strXTime+ 
    '  lon,lat: '+str(lonMin)+', '+str(latMin)+'       date: '+today.strftime("%d/%m/%Y"))

figfile = 'vert_profiles_' +simName+'_t'+str(iTime)+ '.png'
print(figfile)
plt.savefig(figfile, bbox_inches='tight')
plt.close()
