#!/usr/bin/env python
"""
    Aggregate particle data in lat/lon (radians) coordinates
    in 'all_rlzn_particle_data_rads.npz' file.

    Phillip Wolfram
    LANL
    11/16/2015
"""
# import libraries / packages
import glob
import numpy as np
import pandas as pd
import calendar
import os
import netCDF4
from scipy.spatial import cKDTree as KDTree
from mpas_xarray import preprocess_mpas, remove_repeated_time_index
import xarray as xr
from latlon_coordinate_transforms import proj_lat_long_numexpr as proj_lat_long
from iotasks import timeit_context, savenpz

def aggregate_positions(particlefiles, Ntr, planar=False, rlznrange=None, rnumsavetemp=None, outdir='./', subset=None, resets=False):  #{{{
    """
    For all the particle file load the file and the time-based
    spread for the ensembles.

    Phillip Wolfram
    LANL
    12/01/2015
    """

    with timeit_context('Get file list'):
      files = sorted(glob.glob(particlefiles))
      print 'Using these files: ', files

    with timeit_context('Load in the xr datasets'):
      rlzns = xr.open_mfdataset(particlefiles, preprocess=preprocess_mpas, concat_dim='Time')
      mesh = xr.open_dataset('mesh.nc')

    with timeit_context('Remove duplicated times'):
      rlzns = remove_repeated_time_index(rlzns)

    with timeit_context('Load reset data from xr dataset'):
      nReset = rlzns.numTimesReset[:,0].values

    if resets:
      with timeit_context('Reduce dataset:'):
        # find all the "starting" points, assuming first one doesn't count
        # only need to do this for a single particle (other particles are on same reset freq)
        reset = np.diff(nReset,axis=0)
        reset = np.concatenate((np.zeros((1)), reset),axis=0)
        # exclude last entry (possibly incomplete)
        if np.max(reset):
          reset[np.where(reset)[0][-1]] = 0
          for i in np.where(reset)[0]:
            reset[i:(i+Ntr)] = 1
          rlzns = rlzns.isel(Time=np.where(reset > 0)[0])
        else:
          print 'No particle resets observed in the file'

      if subset is not None:
        # make sure zero is in subset (don't allow duplicates)
        subset = np.union1d(subset, np.array([0]))
        rlzns = rlzns.isel(Time=subset)
        Ntr = len(rlzns.coords['Time'])
      # must manually chunk because of a bug in the dask code (issue has been reported and fix is in the works)
      rlzns = rlzns.chunk({'Time':100})
      del nReset, reset

    with timeit_context('Get dimensions'):
      NtNtr = len(rlzns.coords['Time'])
      Nr = np.maximum(1,int(np.floor(NtNtr/Ntr)))
      NpNb = len(rlzns.coords['nParticles'])
      Nb = len(rlzns.coords['nBuoyancySurfaces'])
      Np = NpNb/Nb
      Nc = len(mesh.xCell)
      arraydims = {'Nr':Nr, 'Nb':Nb, 'Np':Np, 'Nc':Nc}
      print arraydims

    with timeit_context('Compute tree'):
      tree = KDTree(zip(mesh.xCell.values.ravel(), mesh.yCell.values.ravel()))

    with timeit_context('Compute buoyancySurfaces'):
      buoyancySurface = rlzns.buoyancySurfaceValues.isel(Time=0).values

    # compute realization range
    if rlznrange is None:
      rlznrange = np.arange(Nr)
    else:
      # sanitize to make sure only entries in range are computed
      rlznrange = np.intersect1d(rlznrange, np.arange(Nr))
    print 'Computing realizations for ', rlznrange

    # loop over each realization (and produce files)
    for rnum in rlznrange:
      with timeit_context('Computed realization %d'%(rnum)):

        with timeit_context('Compute time deltas'):
          # used to check time deltas in the case of uneven output with conversion to days
          dt = np.asarray(np.diff(np.reshape(rlzns.Time[rnum*Ntr:(rnum+1)*Ntr].values, \
              (1,Ntr)),axis=1),dtype='f8')/(1.0e9*86400.)
          if rlzns.config_calendar_type == u'gregorian_noleap':
            # compute days to shift
            subtractday = -1
            datetimes = pd.to_datetime(rlzns.Time[rnum*Ntr:(rnum+1)*Ntr].values)
            leapyear = np.asarray([calendar.isleap(at.year) for at in datetimes[:]])
            febsplit = np.concatenate((np.array([False]), \
                np.logical_and(datetimes[1:].month == 3, datetimes[:-1].month == 2)),axis=0)
            # reshape
            #datetimes = np.reshape(datetimes, (1,Ntr))
            leapyear = np.reshape(leapyear, (1,Ntr))
            febsplit = np.reshape(febsplit, (1,Ntr))
            # adjust for difference
            leapyear = leapyear[:,1:]
            febsplit = febsplit[:,1:]
            # apply shift
            dt += leapyear*febsplit*subtractday
          dtrlzn = dt
          del dt, leapyear, febsplit, datetimes

        with timeit_context('Compute mask for outcropped layers at cells'):
          # compute outcropped layers
          depth = rlzns.buoyancySurfaceDepth[rnum*Ntr:(rnum+1)*Ntr,:,:].values
          mask = np.ones_like(depth, dtype='bool')
          same = np.asarray(np.diff(depth,axis=-1), dtype='bool')
          mask[:,:,1:] *= same
          mask[:,:,:-1] *= same
          # mask is true for all cells on layers that are not outcropped ("good data")
          del depth, same

        with timeit_context('Perforing the projection'):
          if planar:
            lat = rlzns.yParticle[rnum*Ntr:(rnum+1)*Ntr,:].values
            lon = rlzns.xParticle[rnum*Ntr:(rnum+1)*Ntr,:].values
          #else:
          #  lat, lon = proj_lat_long(x[rnum*Ntr:(rnum+1)*Ntr,:].values.ravel(), y[rnum*Ntr:(rnum+1)*Ntr,:].values.ravel(), z[rnum*Ntr:(rnum+1)*Ntr,:].values.ravel())

        with timeit_context('Reshaping arrays'):
          lon = np.reshape(lon, (1, Ntr, Nb, Np))
          lat = np.reshape(lat, (1, Ntr, Nb, Np))
          mask = np.reshape(np.swapaxes(mask,1,2), (1, Ntr, Nb, Nc))

        # compute the particle mask for outcropped layers
        if planar: #{{{
          with timeit_context('Compute tree for planar mesh and time-varying particle IDs'):
            # mask all particles associated with an outcropped cell
            # get cell index corresponding to each particle
            # computed tree above
            # tree = KDTree(zip(mesh.xCell.values.ravel(), mesh.yCell.values.ravel()))
            _, iCell = tree.query(np.concatenate(\
                (lon[:,:,:,:,np.newaxis],\
                 lat[:,:,:,:,np.newaxis]),axis=-1),k=1)
            # del tree

          with timeit_context('Compute mask for planar outcropped particles'):
            assert Np == Nc, 'Number of particles is not equivalent to number of cells!  Computation of mask unclear.'
            iCell = np.reshape(iCell,(1*Ntr*Nb,Np))
            mask = np.reshape(mask,(1*Ntr*Nb,Nc))
            directmask = np.ones_like(iCell)
            # not the cleanest approach, probably could optimize the below
            for jj,cells in enumerate(iCell):
                  directmask[jj,:] = mask[jj,cells]
            del mask, iCell
            mask = np.reshape(directmask, (1,Ntr,Nb,Np))
            del directmask

          #}}}
        else: #{{{
          with timeit_context('Compute mask for spherical outcropped particles'):
            print 'WARNING: not implemented for spherical data!'
          #}}}

        with timeit_context('Getting time vector'):
          xtime = rlzns.Time[rnum*Ntr:(rnum+1)*Ntr].values
          yearoffset = rlzns.attrs['time_yearoffset']

        with timeit_context('Saving data'):
          if rnumsavetemp is None:
            rnumsave = rnum
          else:
            rnumsave = rnumsavetemp
          #savenpz('all_rlzn_particle_data_rads_rlzn%04d'%(rnum),lonrad=lon,latrad=lat,
          #    xtime=xtime, yearoffset=yearoffset,arraydims=arraydims, \
          #    buoyancySurface=buoyancySurface, dtdays=dtrlzn, notoutcropped=mask)
          xr.Dataset({'lon'           : (['Nr', 'Time', 'Nb', 'Np'], lon),  \
                      'lat'           : (['Nr', 'Time', 'Nb', 'Np'], lat),  \
                      'notoutcropped' : (['Nr', 'Time', 'Nb', 'Np'], mask), \
                      'dtdays'        : (['Nr', 'Nt-1'], dtrlzn),           \
                      },\
                      coords={'rlzn'       : (['Nr'], np.asarray([rnumsave])), \
                              'Nb'         : buoyancySurface,              \
                              'time'       : (['Time'],xtime),             \
                              'yearoffset' : yearoffset                    \
                              }\
                      )\
                      .to_netcdf(outdir + 'all_rlzn_particle_data_rads_rlzn%04d.nc'%(rnumsave))
          print 'Saved file to ' + 'all_rlzn_particle_data_rads_rlzn%04d.nc'%(rnumsave)

    with timeit_context('Closing xr files'):
      rlzns.close()
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
    parser.add_option("-r", "--rlznrange", dest="rlznrange", help="range of realizations to compute", metavar="str")
    parser.add_option("-s", dest="rnumsave", help="number of realizations to save", metavar="str")
    parser.add_option("-o", dest="outdir", help="location of output dir", metavar="str")
    parser.add_option('--subset', dest='subset', help='time indices of days to subset')

    options, args = parser.parse_args()

    if not options.particlefile:
        assert False, 'Must specify a particle file!'
    print 'particle file name is ', options.particlefile
    if not options.nrealization:
        options.nrealization = 12
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
    if not options.rnumsave:
        options.rnumsave = None
    else:
        options.rnumsave = int(options.rnumsave)
    if not options.rlznrange:
        options.rlznrange = None
    else:
        # e.g., via -r 'np.arange(0,10)'
        options.rlznrange = eval(options.rlznrange)
    if not options.outdir:
      options.outdir = './'
    else:
      print options.outdir
      #options.outdir += '/'
    if not options.subset:
      options.subset=None
    else:
      options.subset = eval(options.subset)

    aggregate_positions(options.particlefile, options.nrealization, options.planar, options.rlznrange, rnumsavetemp=options.rnumsave, outdir=options.outdir, subset=options.subset)
