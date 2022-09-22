#!/usr/bin/env python
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s/'
simName = 's03g/'
fileName = 'output.nc'
meshName = 'init.nc'

lonMin = 0.193
lonMax = 0.210
latMin = 0.477
latMax = 0.490

#from __future__ import absolute_import, division, print_function, \
#    unicode_literals
import os
import glob
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cols
from matplotlib.pyplot import cm
from matplotlib.colors import BoundaryNorm
import cmocean
import xarray as xr
from netCDF4 import Dataset
from mpas_analysis.shared.io.utility import decode_strings
import gsw

print('reading data from:')
print(runDir+simName+meshName)
print(runDir+simName+fileName)
mesh = xr.open_dataset(runDir+simName+meshName)
data = xr.open_dataset(runDir+simName+fileName)

z = mesh.refBottomDepth.values
nVertLevels = len(z)
nCells = data.dims['nCells']

latCell = mesh.variables['latCell']
lonCell = mesh.variables['lonCell']
maxLevelCell = mesh.variables['maxLevelCell']

cellList = np.where(np.logical_and(np.logical_and(np.logical_and(latCell>latMin, latCell<latMax), lonCell>lonMin), lonCell<lonMax))[0]

iTime = 3
fig = plt.figure(figsize=(20,12))
varNames = ['divergence','vertVelocityTop','vertTransportVelocityTop','temperature', 'salinity','density','pressure','zMid']
varNames = ['divergence','vertVelocityTop','temperature', 'salinity']
for j in range(len(varNames)):
    plt.subplot(2,2,j+1)
    var = data.variables[varNames[j]]
    for i in range(len(cellList)):
        iCell = int(cellList[i])
        k = int(maxLevelCell[iCell])
        varData = var[iTime,iCell,0:k]
        #varData = np.where(varData>-1e20,varData,np.NAN)
        plt.plot(varData,np.arange(1,k+1))
    plt.gca().invert_yaxis()
    plt.title(varNames[j])
    plt.grid()

figfile = 'test.png'
plt.savefig(figfile, bbox_inches='tight')
plt.close()
