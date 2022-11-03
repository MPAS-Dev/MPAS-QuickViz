#!/usr/bin/env python
import netCDF4
import numpy as np

def build_particle_simple(f_grid, f_name, f_decomp, vertexReconst, horizTreatment,
        timeInt, verticalTreatment,buoySurf): #{{{

    nBuoySurf = buoySurf.shape[0]

    f_grid = netCDF4.Dataset(f_grid,'r')

    nparticles = len(f_grid.dimensions['nCells'])
    f_out = netCDF4.Dataset(f_name, 'w',format='NETCDF3_CLASSIC')
    f_out.createDimension('nParticles', nparticles*nBuoySurf)
    f_out.createDimension('nBuoyancySurfaces', nBuoySurf)
    f_out.createDimension('Time')

    f_out.createVariable('xParticle', 'f8', ('Time','nParticles'))
    f_out.createVariable('yParticle', 'f8', ('Time','nParticles'))
    f_out.createVariable('zParticle', 'f8', ('Time','nParticles'))
    f_out.createVariable('dtParticle', 'f8', ('Time','nParticles'))
    f_out.createVariable('buoyancyParticle', 'f8', ('Time','nParticles'))
    f_out.createVariable('currentBlock', 'i', ('Time', 'nParticles'))
    f_out.createVariable('currentCell', 'i', ('Time', 'nParticles'))

    f_out.createVariable('indexToParticleID', 'i', ('nParticles'))
    f_out.createVariable('vertexReconstMethod', 'i', ('Time','nParticles'))
    f_out.createVariable('timeIntegration', 'i', ('Time','nParticles'))
    f_out.createVariable('horizontalTreatment', 'i', ('Time','nParticles'))
    f_out.createVariable('verticalTreatment', 'i', ('Time','nParticles'))
    f_out.createVariable('indexLevel', 'i', ('Time','nParticles'))
    f_out.createVariable('buoyancySurfaceValues', 'f8', ('nBuoyancySurfaces'))

    # reset variables
    f_out.createVariable('resetTime', 'i', ('nParticles'))
    f_out.createVariable('currentBlockReset', 'i', ('nParticles'))
    f_out.createVariable('currentCellReset', 'i', ('nParticles'))
    f_out.createVariable('xParticleReset', 'f8', ('nParticles'))
    f_out.createVariable('yParticleReset', 'f8', ('nParticles'))
    f_out.createVariable('zParticleReset', 'f8', ('nParticles'))

    f_out.variables['xParticle'][0,:] = np.tile(f_grid.variables['xCell'][:],(nBuoySurf))
    f_out.variables['yParticle'][0,:] = np.tile(f_grid.variables['yCell'][:],(nBuoySurf))
    f_out.variables['zParticle'][0,:] = np.tile(f_grid.variables['zCell'][:],(nBuoySurf))
    f_out.variables['buoyancyParticle'][0,:] = (np.tile(buoySurf,(nparticles,1))).reshape(nparticles*nBuoySurf,order='F').copy()
    f_out.variables['buoyancySurfaceValues'][:] = buoySurf
    f_out.variables['dtParticle'][0,:] = 300.0 #5 min
    # assume single-processor mode for now
    f_out.variables['currentBlock'][:] = 0
    f_out.variables['resetTime'][:] = 1.0*24.0*60.0*60.0 # reset each day
    f_out.variables['vertexReconstMethod'][:] = vertexReconst
    f_out.variables['timeIntegration'][:] = timeInt
    f_out.variables['horizontalTreatment'][:] = horizTreatment
    f_out.variables['verticalTreatment'][:] = verticalTreatment
    f_out.variables['indexLevel'][:] = 1
    f_out.variables['indexToParticleID'][:] = np.arange(nparticles*nBuoySurf)

    # resets
    cellIndices = np.tile(np.arange(nparticles),(nBuoySurf))
    decomp = np.genfromtxt(f_decomp)
    f_out.variables['currentBlock'][0,:] = decomp[cellIndices]
    f_out.variables['currentBlockReset'][:] = decomp[cellIndices]
    # can't use these because the specific reset cell (local cell ID) is not known
    #f_out.variables['currentCell'][0,:] = cellIndices + 1
    #f_out.variables['currentCellReset'][:] = cellIndices + 1
    f_out.variables['currentCell'][0,:] = -1
    f_out.variables['currentCellReset'][:] = -1
    f_out.variables['xParticleReset'][:] = f_out.variables['xParticle'][0,:]
    f_out.variables['yParticleReset'][:] = f_out.variables['yParticle'][0,:]
    f_out.variables['zParticleReset'][:] = f_out.variables['zParticle'][0,:]

    f_out.close()

    return #}}}
build_particle_simple('mesh.nc', 'particles.nc', 'graph.info.part.2016', 1, 1, 2, 4, np.linspace(1028.74,1029.7,10))
