#!/usr/bin/env python

import numpy as np
import xarray as xr
from iotasks import timeit_context

def compute_autocorrelation(lon0, lat0):

  # get the data files
  with timeit_context('open datafiles'):
    ds = xr.open_mfdataset('all_r*nc',concat_dim='Nr')
    ds = ds.chunk({'Np':2280})
    print ds.chunks

  # compute the coordinates (noting all particles start at same lonstart and latstart locations
  with timeit_context('set coordinates'):
    lonstart = ds.lon[0,0,0,:].drop(['Time','time','rlzn','Nr','Nb'])
    latstart = ds.lat[0,0,0,:].drop(['Time','time','rlzn','Nr','Nb'])

    # set the coordinates
    ds = ds.merge(xr.Dataset({'lonstart':lonstart,'latstart':latstart}))
    ds = ds.set_coords(['lonstart','latstart'])

  # remove periodic x direction and get the nearest points (only need latitude)
  with timeit_context('reduce dataset'):
    uniquey = np.unique(ds.latstart.values)
    dist = (uniquey - lat0)**2.0
    idx = np.where(dist == np.min(dist))[0]
    points = ds.sel(Np=(ds.latstart == uniquey[idx]))

  with timeit_context('load points'):
    points.load()

  # operate on points data
  def diff(dsT):
    return dsT.diff('Time',label='lower')

  with timeit_context('compute velocities'):
    dt = diff(points.Time)
    velx = diff(points.lon)/dt
    vely = diff(points.lat)/dt

    velx0 = velx.isel(Time=0)
    vely0 = vely.isel(Time=0)

  def avg(dsT):
    return dsT.mean(['Np','Nr'])

  with timeit_context('compute autocorrelation'):
    rhou = avg(velx0*velx)/avg(velx0*velx0)
    rhov = avg(vely0*vely)/avg(vely0*vely0)
    rho  = ( avg(velx0*velx) + avg(vely0*vely) + avg(vely0*velx) + avg(velx0*vely) )/( avg(velx0*velx0) + avg(vely0*vely0) )

  with timeit_context('compute autocorrelation'):
    TLu = (dt*rhou).sum('Time')
    TLv = (dt*rhov).sum('Time')
    TLc = (dt*rho ).sum('Time')

  with timeit_context('merge datasets'):
    rhos = xr.Dataset({'rhou': rhou, 'rhov': rhov, 'rho': rho, 'TLu':TLu, 'TLv':TLv, 'TLc':TLc})
    points = points.merge(rhos)

  with timeit_context('save to file'):
    points.to_netcdf('lat0_1500km.nc')

if __name__ == "__main__":
  compute_autocorrelation(lon0=500e3,lat0=1500e3)
