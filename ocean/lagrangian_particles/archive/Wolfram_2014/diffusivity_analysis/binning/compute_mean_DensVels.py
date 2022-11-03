#!/usr/bin/env python
"""
Phillip J. Wolfram
12/04/2015
"""
import numpy as np
import netCDF4
import xray
from mpas_xray import preprocess_mpas, remove_repeated_time_index
import matplotlib as mpl
import socket
if socket.gethostname() == 'Shapiro' or socket.gethostname() == 'shapiro.lanl.gov':
  mpl.use('MacOSX')
else:
  mpl.use('Agg')
import matplotlib.pyplot as plt
from iotasks import timeit_context

with timeit_context('Load datasets'):
  ds = xray.open_mfdataset('output.01*nc', preprocess=preprocess_mpas)
  ds = remove_repeated_time_index(ds)

with timeit_context('Adjust coodrinates'):
  ds.yCell = ds.yCell.mean('Time')
  ds = ds.update({'yCell':ds.yCell})
  ds = ds.set_coords('yCell')

with timeit_context('Rechunk'):
  ds = ds.chunk({'Time':2000, 'nCells':2000})

with timeit_context('Lazy-compute time-mean, and x-averaged potentialDensity, zMid, zonal/meridional Vel'):
  potDens = xray.concat([ds.potentialDensity[:,:,i].groupby('yCell').mean() \
      for i in np.arange(ds.potentialDensity.shape[2])], 'nVertLevels')
  Uvel = xray.concat([ds.velocityZonal[:,:,i].groupby('yCell').mean() \
      for i in np.arange(ds.velocityZonal.shape[2])], 'nVertLevels')
  Vvel = xray.concat([ds.velocityMeridional[:,:,i].groupby('yCell').mean() \
      for i in np.arange(ds.velocityMeridional.shape[2])], 'nVertLevels')
  zMid = xray.concat([ds.zMid[:,:,i].groupby('yCell').mean() \
      for i in np.arange(ds.zMid.shape[2])], 'nVertLevels')

with timeit_context('Convert to dataset'):
  pds = potDens.to_dataset()
  pds = pds.merge(Uvel.to_dataset())
  pds = pds.merge(Vvel.to_dataset())
  pds = pds.merge(zMid.to_dataset())

with timeit_context('Load values'):
  pds.load()

with timeit_context('Save file'):
  pds.to_netcdf('potDens.nc')

with timeit_context('Plot'):
  data = netCDF4.Dataset('potDens.nc','r')
  y = data['yCell'][:]
  z = data['zMid'][:]
  den = data['potentialDensity'][:,:]
  uVel = data['velocityZonal'][:,:]
  vVel = data['velocityMeridional'][:,:]
  y = y*np.ones_like(z)/1000.
  den -= 1029
  plt.figure()
  plt.contourf(y,z,den, levels=np.arange(-1.,1.5,0.05))
  plt.colorbar()
  plt.contour(y,z,den,levels=np.arange(-1.,1.5,0.05),colors='k')
  plt.xlabel('Meridional location (km)')
  plt.ylabel('Depth (m)')
  plt.title('zonal-mean: potential density (kg/m3)')
  plt.savefig('potDensity.png')

  plt.figure()
  plt.contourf(y,z,uVel, levels=np.arange(-1.0,1.0,0.05),cmap='RdYlBu_r')
  plt.colorbar()
  plt.contour(y,z,uVel,levels=np.arange(-1.0,1.0,0.05),colors='k')
  plt.xlabel('Meridional location (km)')
  plt.ylabel('Depth (m)')
  plt.title('zonal-mean: zonal velocity(m/s)')
  plt.savefig('zonalVel.png')

  plt.figure()
  plt.contourf(y,z,vVel, levels=np.arange(-0.1,0.1,0.01),cmap='RdYlBu_r')
  plt.colorbar()
  plt.contour(y,z,vVel,levels=np.arange(-0.1,0.1,0.01),colors='k')
  plt.xlabel('Meridional location (km)')
  plt.ylabel('Depth (m)')
  plt.title('zonal-mean: meridional velocity(m/s)')
  plt.savefig('meridionalVel.png')
