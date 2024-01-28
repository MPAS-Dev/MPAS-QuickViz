#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
simple scatter plot
Mark Petersen
Jan 2024
"""

############################## model files, run dirs
# running in /lustre/scratch5/mpeterse/runs/240126_wind_stress_test/ocean/global_ocean/QU240/WOA23/performance_test/forward
runDir = './'
meshFileName = 'mesh.nc'
outputFileName = 'output.nc'

import numpy as np
import xarray
import matplotlib as mpl
import matplotlib.pyplot as plt

deg2rad = 3.14159/180.0
rad2deg = 180.0/3.14159

print('read: '+runDir+meshFileName)
mesh = xarray.open_dataset(runDir+meshFileName)
nCells = mesh.sizes['nCells']
latCell = mesh.variables['latCell']
lonCell = mesh.variables['lonCell']

print('read: '+runDir+outputFileName)
output = xarray.open_dataset(runDir+outputFileName)

enable_plotting = True
if (enable_plotting):
    print('save plot')
    fig = plt.figure(figsize=(20,12))
    varNames =['velocityZonal','velocityMeridional']
    iLev = 0
    iTime= 0
    nCol = 1
    nRow = 2
    iVar = 0
    for iRow in range(nRow):
        for iCol in range(nCol):
            var = output.variables[varNames[iVar]]
            plt.subplot(nRow, nCol, iCol * nCol + iRow + 1)
            plt.scatter(lonCell*rad2deg, latCell*rad2deg ,c=var[iTime,:,iLev], cmap='bwr', s=1)
            plt.colorbar()
            plt.axis('off')
            plt.title(varNames[iVar])
            iVar += 1
    
    plt.savefig('Output.png')
