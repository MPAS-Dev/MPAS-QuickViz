#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/220927_EC30to60_zlevel_noPBC_noSmooth/ocean/global_ocean/EC30to60/PHC/init/initial_state/'
fileName = 'temperature.nc'
fileName = 'salinity.nc'
fileName = 'phc3.0_annual.nc'
fileNames = ['PotentialTemperature.01.filled.60levels.PHC.151106.nc','Salinity.01.filled.60levels.PHC.151106.nc']
newFileNames = ['PotentialTemperature.01.zonalAvg.60levels.PHC.220927.nc','Salinity.01.zonalAvg.60levels.PHC.220927.nc']

import matplotlib as mpl
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from netCDF4 import Dataset

#for fileName in fileNames:
varNames = ['TEMP','SALT']
for iVar in range(2):
   print(runDir+fileNames[iVar])
   data = xr.open_dataset(runDir+fileNames[iVar])
   var = data.variables[varNames[iVar]]
   for k in range(60):
       for i in range(180):
           var[k,i,:] = np.average(var[k,i,:])
   data.to_netcdf(path=runDir+newFileNames[iVar])
