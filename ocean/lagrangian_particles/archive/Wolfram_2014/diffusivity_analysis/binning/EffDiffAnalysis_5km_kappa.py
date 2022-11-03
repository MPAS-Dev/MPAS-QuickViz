import sys
import socket
import numpy as np
import matplotlib as mpl
if socket.gethostname() == 'Shapiro' or socket.gethostname() == 'shapiro.lanl.gov':
    mpl.use('MacOSX')
else:
    mpl.use('Agg')
import xarray as xr
import GeneralTriangulation as GT
import effective_diffusivity as ed

# example call
# python EffDiffAnalysis_5km_long.py 1 '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/realizations/realization_24-01/' '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc'

# target buoyancy surface index
alayer = int(sys.argv[1])

# target realization
rlzn = 0

# mesh name
mname = sys.argv[3]

# open dataset and make assignments
ds = xr.open_dataset(sys.argv[2] + '/all_rlzn_particle_data_rads_rlzn0000.nc')
x = ds['lon'].values
y = ds['lat'].values
dt = ds['dtdays'].values
notoutcropped = ds['notoutcropped'].values
buoyancySurface = ds.Nb.values
tvec = np.hstack((np.array([0.0]), np.cumsum(ds['dtdays'][rlzn])))*86400.

print 'Computing on layer %s for potential density surface %s'%(alayer, buoyancySurface[alayer])
print '%d unique y values'%(len(np.unique(xr.open_dataset(mname).yCell)))
print ds

# compute effective diffusivity
effdiff = ed.EffectiveDiffusivity(np.squeeze(x[rlzn,:,alayer,:]),np.squeeze(y[rlzn,:,alayer,:]),
                    np.squeeze(notoutcropped[rlzn,0,alayer,:]),
                    tvec[:], meshname=mname, advection=True, diffusion=True, y0=1500e3,
                    kappa=100., calctype='ygrad', dt=12*60*60., dtout=1.*86400.,
                    saveloc=sys.argv[2] + '/double_kappa_eff_diff_conc_evolution_gradient_layer%04d/'%(alayer))
effdiff.compute_time_advancement()
effdiff.plot_scalar_variance()
effdiff.plot_effective_diffusivity()
effdiff.plot_effective_lengths()
