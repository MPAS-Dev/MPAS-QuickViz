#!/usr/bin/env python
"""
Phillip J. Wolfram
12/03/2015
"""
import numpy as np
import xarray as xray
from mpas_xarray import preprocess_mpas, remove_repeated_time_index
import matplotlib as mpl
import socket
if socket.gethostname() == 'Shapiro' or socket.gethostname() == 'shapiro.lanl.gov':
  mpl.use('MacOSX')
else:
  mpl.use('Agg')
import matplotlib.pyplot as plt
from iotasks import timeit_context

with timeit_context('Load in datasets'):
  ds = xray.open_mfdataset('lagrPartTrack.0*nc', preprocess=preprocess_mpas)
  ds = remove_repeated_time_index(ds)
  mesh = xray.open_dataset('mesh.nc')

with timeit_context('Merge datasets'):
  dst = ds.merge(mesh)

with timeit_context('Set coordinates'):
  dst = dst.set_coords(['xCell','yCell'])

with timeit_context('Rechunk'):
  dst = dst.chunk({'Time':1000, 'nCells':2000})

with timeit_context('Compute mean depth averaging over time and x for each buoyancy surface'):
  meandepth = xray.concat([dst.buoyancySurfaceDepth[:,:,i].groupby('yCell').mean() \
      for i in np.arange(dst.buoyancySurfaceDepth.shape[2])], 'nBuoyancySurfaces')

with timeit_context('Compute mean zonal velocity averaging over time and x for each buoyancy surface'):
  meanU = xray.concat([dst.buoyancySurfaceVelocityZonal[:,:,i].groupby('yCell').mean() \
      for i in np.arange(dst.buoyancySurfaceVelocityZonal.shape[2])], 'nBuoyancySurfaces')

with timeit_context('Compute mean meridional velocity averaging over time and x for each buoyancy surface'):
  meanV = xray.concat([dst.buoyancySurfaceVelocityMeridional[:,:,i].groupby('yCell').mean() \
      for i in np.arange(dst.buoyancySurfaceVelocityMeridional.shape[2])], 'nBuoyancySurfaces')

with timeit_context('Perform the computations for mean depth'):
  meandepth.load()
  meandepth = meandepth.rename('buoyancySurfaceDepthMean')

with timeit_context('Perform the computations for zonal velocity'):
  meanU.load()
  meanU = meanU.rename('buoyancySurfaceVelocityZonalMean')

with timeit_context('Perform the computations for meridional velocity'):
  meanV.load()
  meanV = meanV.rename('buoyancySurfaceVelocityMeridionalMean')

with timeit_context('Compute std zonal velocity averaging over time and x for each buoyancy surface'):
  stdU = xray.concat([dst.buoyancySurfaceVelocityZonal[:,:,i].groupby('yCell').std() \
      for i in np.arange(dst.buoyancySurfaceVelocityZonal.shape[2])], 'nBuoyancySurfaces')

with timeit_context('Compute std zonal velocity averaging over time and x for each buoyancy surface'):
  stdV = xray.concat([dst.buoyancySurfaceVelocityMeridional[:,:,i].groupby('yCell').std() \
      for i in np.arange(dst.buoyancySurfaceVelocityMeridional.shape[2])], 'nBuoyancySurfaces')

with timeit_context('Perform the computations for zonal eddy velocity'):
  stdU.load()
  stdU = stdU.rename('buoyancySurfaceVelocityZonalStdDev')

with timeit_context('Perform the computations for meridional eddy velocity'):
  stdV.load()
  stdV = stdV.rename('buoyancySurfaceVelocityMeridionalStdDev')

with timeit_context('Append buoyancy surface values to the combined datasets'):
  comb = meandepth.to_dataset().merge(meanU.to_dataset())
  comb = comb.merge(meanV.to_dataset())
  comb = comb.merge(stdU.to_dataset())
  comb = comb.merge(stdV.to_dataset())
  comb = comb.merge(dst.buoyancySurfaceValues[0,:].to_dataset())

with timeit_context('Save the file'):
  comb.to_netcdf('buoyancySurfaceMeans.nc')

with timeit_context('Plot the file'):
  ds = xray.open_dataset('buoyancySurfaceMeans.nc')
  plt.figure()
  for i,a in enumerate(ds.buoyancySurfaceDepth):
        plt.plot(ds.yCell/1000,a,label=ds.buoyancySurfaceValues[i].values, lw=2)
        plt.xlabel('Meridional location (km)')
        plt.ylabel('Depth (m)')
        plt.legend(loc='best', prop={'size':10})
        plt.title('Mean buoyancy surface depths')
        plt.savefig('meandepths.png')

