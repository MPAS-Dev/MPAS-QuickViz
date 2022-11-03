#!/usr/bin/env python

# initialize general libraries and setup notebook
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import cluster_topology as CT
import GeneralTriangulation as GT
import effective_diffusivity as ED

# target buoyancy surface
buoySurf = 1026.85

# target realization
rlzn = 0

# gaussian properties
yloc = 1500e3
lgauss = 50e3

# mesh name
mname = 'mesh.nc'

# run intermediate tests
runtests = True

ds = np.load('SOMA.npz')
x = np.radians(ds['lon'])
y = np.radians(ds['lat'])
dt = ds['dtdays']
buoyancySurface = ds['buoyancySurface']
print dt.shape, x.shape
tvec = np.zeros(dt.shape[1]+1)
tvec[1:] = np.cumsum(dt[rlzn,:]*86400.)

alayer = np.where(buoyancySurface == buoySurf)
print 'Buoyancy Surface of interest: ', buoyancySurface[alayer], 'among ', buoyancySurface

effdiff = ED.EffectiveDiffusivity(xt=np.squeeze(x[rlzn,:,alayer,:]), yt=np.squeeze(y[rlzn,:,alayer,:]), \
                                        latlon=True, calctype='gauss', cleanupgrid=10, \
                                        tparticles=tvec, advection=True, diffusion=True, \
                                        kappa=10., dt=60., dtout=1.*86400., \
                                        x0=np.array([np.radians(-7.5)]), y0=np.array([np.radians(35.)]), \
                                        concwidth=100.e3, \
                                        run=True)
effdiff.plot_scalar_variance()
effdiff.plot_effective_diffusivity()
effdiff.plot_effective_lengths()
#plt.show()
