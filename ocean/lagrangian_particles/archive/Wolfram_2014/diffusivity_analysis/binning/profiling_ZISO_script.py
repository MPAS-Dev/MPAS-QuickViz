#!/usr/bin/env python

# initialize general libraries and setup notebook
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from multiprocessing.dummy import Pool as ThreadPool
import cluster_topology as CT
import GeneralTriangulation as GT
import effective_diffusivity as ED
from iotasks import savenpz

def split_array(array, npieces): #{{{
  nvalues = np.ceil(len(array)/npieces)
  split = []
  for i in np.arange(npieces):
    split.append(array[i*nvalues:((i+1)*nvalues)])
  return split #}}}

# target buoyancy surface
buoySurf = 1029.7

# target realization
rlzn = 0
nprocs = 20
procnum = 1

# gaussian properties
yloc = 1500e3
lgauss = 50e3

# mesh name
mname = 'mesh.nc'

# run intermediate tests
runtests = True

ds = np.load('all_rlzn_particle_data_rads.npz')
x = np.radians(ds['lonrad'])
y = np.radians(ds['latrad'])
mask = ds['notoutcropped']
dt = ds['dtdays']
buoyancySurface = ds['buoyancySurface']
print dt.shape, x.shape
tvec = np.zeros(dt.shape[1]+1)
tvec[1:] = np.cumsum(dt[rlzn,:]*86400.)

alayer = np.where(buoyancySurface == buoySurf)
print 'Buoyancy Surface of interest: ', buoyancySurface[alayer], 'among ', buoyancySurface

# total work
ytot = np.linspace(0,2000e3,201)

# note this is a shared-memory (threads) implementation.  This should be ok because we really only need one copy of x/y
def run_eff_diff(alayerID):
    # divide spatial work into chunks
    y0 = split_array(ytot, nprocs)[procnum]
    effdiff = ED.EffectiveDiffusivity(xt=np.squeeze(x[rlzn,:,alayerID,:]), yt=np.squeeze(y[rlzn,:,alayerID,:]), \
                                            mask=mask[rlzn,0,alayerID,:], \
                                            meshname=mname, calctype='ygauss', \
                                            y0=y0, concwidth=lgauss, \
                                            tparticles=tvec[:], advection=True, diffusion=True, \
                                            kappa=10., dt=60., dtout=1.*86400., \
                                            run=True)
    out = effdiff.output()
    savenpz('effective_diffusivity_rlzn_%03d_split_%04d_results_layer%02d.npz'%(rlzn, procnum, alayerID), **out)

    return effdiff

pool = ThreadPool(12)
results = pool.map(run_eff_diff, np.arange(11))
pool.close()
pool.join()
