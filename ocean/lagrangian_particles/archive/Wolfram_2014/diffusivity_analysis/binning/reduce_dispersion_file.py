#!/usr/bin/env python

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import glob
from iotasks import timeit_context

# simple script to reduce files to only have vars in includevars
def reduce_dispersion_file(includevars=['meanu','eddyv','eddyspeed',\
    'K_yy','absK_yy','vp0vp','vp0vp0', \
    'xcluter','ycluster','npart','rlzn','yearoffset', 'time']):

  with timeit_context('building xarray dataset'):
    # get full xr dataset
    dslist = []
    nfiles = len(glob.glob('dispersion_calcs_rlzn0*layerrange_0000-0000.nc'))

    for i in np.arange(nfiles):
      ds = xr.open_mfdataset('dispersion_calcs_rlzn%04d_*nc'%(i))
      dslist.append(ds)
    dstotal = xr.concat(dslist,'Nr')

  with timeit_context('subsetting the data'):
    dstotal = dstotal.isel(Nt=slice(0,30))
    dstotal = dstotal[{'Nt-1':slice(0,29)}]

  with timeit_context('dropping variables'):
    allvars = dstotal.data_vars.keys()
    dropvars = np.setdiff1d(allvars, includevars)
    ds = dstotal.drop(dropvars)

  # potential bug here, may crash during serialization
  # may be related to xarray / dask bug
  with timeit_context('reshape'):
    ds = ds.transpose('Nt','Nt-1','Nr','Nb','Nc')

  with timeit_context('output to disk for 5yr'):
    ds.isel(Nr=slice(0,12*5)).mean('Nr').to_netcdf('diffusivity_results_5yr.nc')

  with timeit_context('output to disk for 10yr'):
    ds.mean('Nr').to_netcdf('diffusivity_results_10yr.nc')

  #with timeit_context('output to disk'):
  #  ds.to_netcdf('combined_dispersion_files.nc')


if __name__ == "__main__":
  reduce_dispersion_file()
