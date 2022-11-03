#!/usr/bin/env python
"""
    Build mask for outcropped particles and ascertain realizations.

    Phillip Wolfram
    LANL
    01/08/2016
"""
# import libraries / packages
import xray
import numpy as np
from scipy.spatial import cKDTree as KDTree
import pandas as pd
import calendar
import glob
import re
from mpas_xray import preprocess_mpas, remove_repeated_time_index
from iotasks import timeit_context, add_meta

def build_masks(particlefiles, Ntr, planar=False):  #{{{

    with timeit_context('Get file list'):
      files = sorted(glob.glob(particlefiles))
      print 'Using these files: ', files

    with timeit_context('Load in the xray datasets'):
      ds = xray.open_mfdataset(particlefiles, preprocess=preprocess_mpas, concat_dim='Time')
      mesh = xray.open_dataset('mesh.nc')

    with timeit_context('Compute mask for outcropped layers at cells'):
      # compute outcropped layers
      depth = ds.buoyancySurfaceDepth.values
      mask = np.ones_like(depth, dtype='bool')
      same = np.asarray(np.diff(depth,axis=-1), dtype='bool')
      mask[:,:,1:] *= same
      mask[:,:,:-1] *= same
      # mask is true for all cells on layers that are not outcropped ("good data")
      del depth, same
      mask = np.reshape(mask, (len(ds.Time), len(ds.nCells)*len(ds.nBuoyancySurfaces)))

    with timeit_context('Load particle data from xray dataset'):
      # make assignments and load data into memory
      x = ds.xParticle.values
      y = ds.yParticle.values
      if not planar:
        z = ds.zParticle.values

    # compute the particle mask for outcropped layers
    if planar: #{{{
        # mask all particles associated with an outcropped cell
        # get cell index corresponding to each particle
        tree = KDTree(zip(mesh.xCell.values.ravel(), mesh.yCell.values.ravel()))
        _, iCell = tree.query(np.concatenate((x[:,:,np.newaxis], y[:,:,np.newaxis]),axis=-1),k=1)
        del tree, x, y
      #}}}
    else: #{{{
      with timeit_context('Compute mask for spherical outcropped particles'):
        print 'WARNING: not implemented for spherical data!'
      #}}}

    with timeit_context('Compute mask for planar outcropped particles'):
      notoutcroppedNow = np.ones_like(iCell)
      # not the cleanest approach, probably could optimize the below
      for jj,cells in enumerate(iCell):
            notoutcroppedNow[jj,:] = mask[jj,cells]
      del mask, iCell

    with timeit_context('Saving data'):
      out = xray.Dataset({\
          'notoutcroppedNow':xray.DataArray(notoutcroppedNow, \
          coords=[(ds.Time.name,ds.Time.values), (ds.nParticles.name,ds.nParticles.values)]),\
          })
      out.attrs['files'] = str(files)
      out = add_meta(out)
      def date_str(fname):
        return re.findall('\d.*\d',fname)[0]
      if len(files) > 1:
        savename = 'valid_realization_masks_from_' + date_str(files[0]) + '_to_' + date_str(files[-1]) + '.nc'
      else:
        savename = 'valid_realization_mask.' + date_str(files[0]) + '.nc'
      print 'Saving to %s'%(savename)
      out.to_netcdf(savename)

    with timeit_context('Closing xray files'):
      ds.close()
      mesh.close()

    return  #}}}

if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="particlefile",
                      help="file/path for analysis",
                      metavar="FILE")
    parser.add_option("-n", "--nrealization", dest="nrealization",
                      help="number of time steps constituting a realization",
                      metavar="FILE")
    parser.add_option("-p", "--planar", dest="planar", help="if on a plane", metavar="BOOL")

    options, args = parser.parse_args()

    if not options.particlefile:
        assert False, 'Must specify a particle file!'
    print 'particle file name is ', options.particlefile
    if not options.nrealization:
        options.nrealization = 30
    else:
        options.nrealization = int(options.nrealization)
    print 'specified number of time steps for a realization is ', options.nrealization
    if not options.planar:
        print 'Warning: assuming data is not planar!'
        options.planar = False
    else:
        if options.planar == 'T':
            options.planar = True
            print 'Data specified to be planar.'
        elif options.planar == 'F':
            options.planar = False
            print 'Data specified to not be planar.'
        else:
            parser.error('Specify planar as "T" or "F"')

    build_masks(options.particlefile, options.nrealization, options.planar)
