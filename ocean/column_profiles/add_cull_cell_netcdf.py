#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s06k/ocean/global_ocean/QU240/mesh/base_mesh/'
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s06r/ocean/global_ocean/QU4/mesh/base_mesh/'
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s06o/ocean/global_ocean/QU30/mesh/base_mesh/'
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s06s/ocean/global_ocean/WC14/mesh/base_mesh/'
fileName = 'base_mesh.nc'
newFileName = 'base_mesh_with_cullCell.nc'

#import matplotlib as mpl
from datetime import date
import numpy as np
#import matplotlib.pyplot as plt
import xarray as xr
from netCDF4 import Dataset

deg2rad = 3.14159/180.0
rad2deg = 180.0/3.14159

# large version:
latMin = 8
latMax = 65
lonMin = -82 + 360
lonMax = -2 + 360

# small version:
latMin = 9
latMax = 43
lonMin = -98 + 360
lonMax = -6 + 360

print('read: '+runDir+fileName)
mesh = xr.open_dataset(runDir+fileName)
nCells = mesh.dims['nCells']
latCell = mesh.variables['latCell']
lonCell = mesh.variables['lonCell']

cullCell = np.ones(nCells, dtype=int)
print('before iCell loop')
#for iCell in range(nCells):
#    if np.logical_and(np.logical_and(np.logical_and(
#        latCell[iCell]>latMin*deg2rad, latCell[iCell]<=latMax*deg2rad), 
#        lonCell[iCell]>lonMin*deg2rad), lonCell[iCell]<=lonMax*deg2rad):
#        cullCell[iCell] = 0

cellList = np.where(np.logical_and(np.logical_and(np.logical_and(
    latCell>latMin*deg2rad, latCell<=latMax*deg2rad), 
    lonCell>lonMin*deg2rad), lonCell<=lonMax*deg2rad))[0]
cullCell[cellList] = 0

# remove triangle south of Central Americal
latMin = 0
latMax = 17
lonMin = -100 + 360
lonMax = -88.5 + 360
cellList = np.where(np.logical_and(np.logical_and(np.logical_and(
    latCell>latMin*deg2rad, latCell<=latMax*deg2rad), 
    lonCell>lonMin*deg2rad), lonCell<=lonMax*deg2rad))[0]
cullCell[cellList] = 1

latMin = 0
latMax = 15
lonMin = -89 + 360
lonMax = -83.8+ 360
cellList = np.where(np.logical_and(np.logical_and(np.logical_and(
    latCell>latMin*deg2rad, latCell<=latMax*deg2rad), 
    lonCell>lonMin*deg2rad), lonCell<=lonMax*deg2rad))[0]
cullCell[cellList] = 1

latMin = 0
latMax = 9
lonMin = -85 + 360
lonMax = -70 + 360
cellList = np.where(np.logical_and(np.logical_and(np.logical_and(
    latCell>latMin*deg2rad, latCell<=latMax*deg2rad), 
    lonCell>lonMin*deg2rad), lonCell<=lonMax*deg2rad))[0]
cullCell[cellList] = 1

# remove lake in Venezuela
latMin = 0
latMax = 10.8
lonMin = -74 + 360
lonMax = -70 + 360
cellList = np.where(np.logical_and(np.logical_and(np.logical_and(
    latCell>latMin*deg2rad, latCell<=latMax*deg2rad), 
    lonCell>lonMin*deg2rad), lonCell<=lonMax*deg2rad))[0]
cullCell[cellList] = 1

mesh["cullCell"]=(['nCells'],  cullCell)
mesh.to_netcdf(path=runDir+newFileName)
