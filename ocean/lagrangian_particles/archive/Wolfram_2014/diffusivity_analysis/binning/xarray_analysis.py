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
import calendar, datetime
import os
import netCDF4
from scipy.spatial import cKDTree as KDTree
from mpas_xarray import preprocess_mpas, remove_repeated_time_index
import xarray as xr
from latlon_coordinate_transforms import proj_lat_long_numexpr as proj_lat_long
from iotasks import timeit_context, savenpz
from latlon_coordinate_transforms import fix_periodic_timeseries, fix_periodicity_numexpr as fix_periodicity
from copy_netcdf_file_without_variables import copy_netcdf_file_without_variables
from cluster_topology import cluster_topology
from intersection_two_lines import first_line_intersection

def xarray_analysis(particlefiles, maskfiles, meshfile, radius, clusterid, layer=None):  #{{{
    """
    Compute analysis metrics for particle data.

    Phillip Wolfram
    LANL
    01/08/2016
    """

    with timeit_context('Get file lists'):
      pfiles = sorted(glob.glob(particlefiles))
      print 'Using these particles files for analysis: ', pfiles
      mfiles = sorted(glob.glob(maskfiles))
      print 'Using these mask files: ', mfiles

    with timeit_context('Load in the xr datasets'):
      particles = xr.open_mfdataset(pfiles, preprocess=preprocess_mpas)
      dsmask = xr.open_mfdataset(mfiles)
      mesh = xr.open_dataset(meshfile)
      # determine if on a plane
      if particles.on_a_sphere == "NO":
        planar = True
      if layer is not None:
        assert layer < len(particles.nBuoyancySurfaces), \
            'Specified layer is not within the range of the buoyancy surfaces'

    with timeit_context('Aggregate datasets'):
      ds = particles.merge(dsmask).merge(mesh)

    with timeit_context('Remove duplicated times'):
      ds = remove_repeated_time_index(ds)

    with timeit_context('Load reset data from xr dataset'):
      # note, assumes all particle resets are identical
      nReset = ds.numTimesReset.values[:,0]

    with timeit_context('Reduce dataset:'):
      # find all the "starting" points, assuming first one doesn't count
      # only need to do this for a single particle (other particles are on same reset freq)
      reset = np.zeros(nReset.shape[0])
      reset[1:] = np.diff(nReset)
      del nReset
      # exclude last entry
      reset[np.where(reset)[0][-1]] = 0
      # note: assumes that all reset frequencies are the same!
      Ntr = np.diff(np.where(reset))[0][0]
      assert len(np.unique(np.unique(np.diff(np.where(reset))[0])))==1, 'Reset frequencies may not be same length!'
      for i in np.where(reset)[0]:
        reset[i:(i+Ntr)] += 1
      # starting value is 2, values to use are 1s
      # selection technique for valid data
      ds = ds.isel(Time=np.where(reset > 0)[0])
      ## note "reset" variable is a mask from the full data set particles to ds

    with timeit_context('Get realization starting and elapsed times'):
      # make an array with the starting time indicies for each realization
      starttime = np.ones_like(ds.Time.values)
      st = np.nan
      for it,at in enumerate(reset[np.where(reset)]):
        # see above 'Reduce dataset' on why starting time is 2
        if at == 2:
          st = ds.Time.values[it]
        starttime[it] = st
      ds = ds.update({'starttime':('Time', starttime)})
      # get elapapsed time for the realizations
      elapsedtime = ds.Time - ds.starttime
      ds = ds.update({'elapsedtime':('Time', elapsedtime)})

    with timeit_context('Get dimensions'):
      NtNtr = len(ds.coords['Time'])
      Nr = int(np.floor(NtNtr/Ntr))
      NpNb = len(ds.coords['nParticles'])
      Nb = len(ds.coords['nBuoyancySurfaces'])
      Np = NpNb/Nb
      Nc = len(mesh.xCell)
      arraydims = {'Nr':Nr, 'Nb':Nb, 'Np':Np, 'Nc':Nc}

    with timeit_context('Compute time deltas'):
      # used to check time deltas in the case of uneven output with conversion to days
      dt = np.asarray(np.diff(np.reshape(ds.Time.values, (Nr,Ntr)),axis=1),dtype='f8')/(1.0e9*86400.)
      if ds.config_calendar_type == u'gregorian_noleap':
        # compute days to shift
        subtractday = -1
        datetimes = pd.to_datetime(ds.Time.values)
        leapyear = np.asarray([calendar.isleap(at.year) for at in datetimes[:]])
        febsplit = np.zeros_like(datetimes,dtype='bool')
        febsplit[1:] = np.logical_and(datetimes[1:].month == 3, datetimes[:-1].month == 2)
        # reshape
        datetimes = np.reshape(datetimes, (Nr,Ntr))
        leapyear = np.reshape(leapyear, (Nr,Ntr))
        febsplit = np.reshape(febsplit, (Nr,Ntr))
        # adjust for difference
        leapyear = leapyear[:,1:]
        febsplit = febsplit[:,1:]
        # apply shift
        dt += leapyear*febsplit*subtractday
      dt = np.reshape(np.hstack((np.nan*np.zeros((dt.shape[0]))[:,np.newaxis],dt)),Ntr*Nr)
      #dtarray = xr.DataArray(dt, coords=[('Nr',np.arange(Nr)), ('Ntr-1',np.arange(Ntr-1))])
      dtarray = xr.DataArray(dt, coords=[('Time',ds.Time.values)])
      ds= ds.merge(xr.Dataset({'deltaTimedays': dtarray}))
      del dt, leapyear, febsplit, datetimes

    with timeit_context('Referesh dimensions'):
      NtNtr = len(ds.coords['Time'])
      Nr = int(np.floor(NtNtr/Ntr))
      NpNb = len(ds.coords['nParticles'])
      Nb = len(ds.coords['nBuoyancySurfaces'])
      Np = NpNb/Nb
      Nc = len(mesh.xCell)
      arraydims = {'Nr':Nr, 'Nb':Nb, 'Np':Np, 'Nc':Nc}
      #reshape command example: np.reshape(var, (Nr, Ntr, Nb, Np))

    with timeit_context('Build up clusters for each mesh point and layer:'):

      # note, current initialization is consistent for the starting time over each
      # layer and has a single particle at each cell center
      for acell in clusterid: #np.arange(len(mesh.nCells)):
        print 'Operating on cell %d...'%(acell)
        if np.mod(acell,np.floor(len(mesh.nCells)/100))==0:
          print str(100*float(acell)/len(mesh.nCells)) + ' %'
        xc = mesh.xCell.values[acell]
        yc = mesh.yCell.values[acell]
        # transform all coordinates into correct coordinates if periodic
        x = mesh.xCell.values
        y = mesh.yCell.values
        if mesh.is_periodic == "YES":
          if mesh.x_period > 0:
            x = fix_periodicity(x,xc,mesh.x_period)
          if mesh.y_period > 0:
            y = fix_periodicity(y,yc,mesh.y_period)

        # using mesh, build up ring from cell neighbors to the particular point,
        # until the maximum distance is one cell distance larger than the radius
        # threshold
        basecells, meshedges, ringedges = cluster_topology(mesh, acell, x, y, radius)

        # assign cluster index and location
        clustercoord = [('cluster' , np.array([int(acell)]))]
        dscluster = xr.Dataset().assign(xCluster=xr.DataArray(np.array([x[acell]]), coords=clustercoord))\
                                .assign(yCluster=xr.DataArray(np.array([y[acell]]), coords=clustercoord))
        del x,y

        layerlist = []
        for layer in ds.nBuoyancySurfaces.values:
          print 'Operating on layer %d...'%(layer)
          dsmasks = xr.Dataset()

          # adust for variable layers
          cells = basecells + layer*len(mesh.nCells)

          #for maskcond in ['all','mask','fullmask','fullmaskstart']:
          for maskcond in ['fullmaskstart']:

            # compute mask conditions #{{{
            if maskcond == 'all':
              mask=True
              skipna=True
            elif maskcond == 'mask':
              # mask each particle and just don't use data that has outcropped
              # however it is used at a later time if it returns to the buoyancy surface
              mask = ds.notoutcroppedNow.isel(nParticles=cells).values.copy()
              skipna=True
            elif maskcond == 'fullmask':
              # excludes paths of particles for later times following a outcrop
              # but still compute statistics with remaining particles
              mask = 1*ds.notoutcroppedNow.isel(nParticles=cells).values.copy()
              startids = np.where(reset == 2)[0]
              for rlzn in np.arange(Nr):
                start = rlzn*Ntr
                end = start + Ntr
                mask[start:end,:] = mask[start:end,:].cumprod(axis=0)
                skipna=True
            elif maskcond == 'fullmaskstart':
              # excludes realizations that have outcropped particles
              mask = 1*ds.notoutcroppedNow.isel(nParticles=cells).values.copy()
              startids = np.where(reset == 2)[0]
              for rlzn in np.arange(Nr):
                start = rlzn*Ntr
                end = start + Ntr
                mask[start:end,:] = mask[start,:]
                skipna=True
            print 'Computing for ' + maskcond
            #}}}

            # get particle clusters (don't apply mask until periodicity is fixed)
            xParticle = ds.xParticle.isel(nParticles=cells)
            yParticle = ds.yParticle.isel(nParticles=cells)

            # correct periodicity  #{{{
            # (redundant computation in terms of whole dataset,
            # however, doesn't require that whole data set be processed which is a
            # more ammenable framework for parallel computation)
            # note this also forces a read into memory
            if mesh.x_period > 0:
              xParticle.values = fix_periodic_timeseries(xParticle.values, mesh.x_period)
            if mesh.y_period > 0:
              yParticle.values = fix_periodic_timeseries(yParticle.values, mesh.y_period)
            #}}}

            # mask outcropped particle positions (hopefully this is a lazy eval...)
            xParticle = xParticle.where(mask > 0)
            yParticle = yParticle.where(mask > 0)

            # fix chunks
            xParticle = xParticle.chunk({'Time':10000, 'nParticles':1000})
            yParticle = yParticle.chunk({'Time':10000, 'nParticles':1000})

            # absolute diffusivity #{{{
            xstart = xParticle.sel(Time=ds.starttime)
            xstart['Time'] = xParticle.Time.values
            dx = xParticle.chunk(xstart.chunks) - xstart

            ystart = yParticle.sel(Time=ds.starttime)
            ystart['Time'] = yParticle.Time.values
            dy = yParticle.chunk(ystart.chunks) - ystart

            def compute_mixing(dsA, dxA, dyA, name_append, pairs=False): #{{{
              """ factor=1 for absolute and relative dispersion, factor=1./4. for particle pairs """

              if pairs:
                factor=1./4.
              else:
                factor=1.

              dr = np.sqrt(dxA*dxA + dyA*dyA)

              # potentially use regular mean for each realization to
              # exclude outcropped data
              dxdx = factor*(dxA*dxA).mean('nParticles', skipna=skipna)
              dxdy = factor*(dxA*dyA).mean('nParticles', skipna=skipna)
              dydy = factor*(dyA*dyA).mean('nParticles', skipna=skipna)
              drdr = factor*(dr*dr).mean('nParticles', skipna=skipna)

              daystosec = 24*60*60.
              Kxx = 0.5*dxdx.chunk({'Time':len(dxdx.Time)}).diff('Time', label='upper')/(dsA.deltaTimedays*daystosec)
              Kxy = 0.5*dxdy.chunk({'Time':len(dxdy.Time)}).diff('Time', label='upper')/(dsA.deltaTimedays*daystosec)
              Kyy = 0.5*dydy.chunk({'Time':len(dydy.Time)}).diff('Time', label='upper')/(dsA.deltaTimedays*daystosec)
              Krr = 0.5*drdr.chunk({'Time':len(drdr.Time)}).diff('Time', label='upper')/(dsA.deltaTimedays*daystosec)

              # use nanmean over ensembles
              meandxdx = dxdx.groupby(dsA.elapsedtime.sel(Time=dsA.Time)).mean(skipna=True)
              meandxdy = dxdy.groupby(dsA.elapsedtime.sel(Time=dsA.Time)).mean(skipna=True)
              meandydy = dydy.groupby(dsA.elapsedtime.sel(Time=dsA.Time)).mean(skipna=True)
              meandrdr = drdr.groupby(dsA.elapsedtime.sel(Time=dsA.Time)).mean(skipna=True)

              meanKxx = Kxx.groupby(dsA.elapsedtime.sel(Time=dsA.Time[1:])).mean(skipna=True)
              meanKxy = Kxy.groupby(dsA.elapsedtime.sel(Time=dsA.Time[1:])).mean(skipna=True)
              meanKyy = Kyy.groupby(dsA.elapsedtime.sel(Time=dsA.Time[1:])).mean(skipna=True)
              meanKrr = Krr.groupby(dsA.elapsedtime.sel(Time=dsA.Time[1:])).mean(skipna=True)

              # group into a dataset
              dsO = xr.Dataset({'dxdx' + name_append : dxdx, \
                                'dxdy' + name_append : dxdy, \
                                'dydy' + name_append : dydy, \
                                'drdr' + name_append : drdr, \
                                'Kxx' + name_append : Kxx, \
                                'Kxy' + name_append : Kxy, \
                                'Kyy' + name_append : Kyy, \
                                'Krr' + name_append : Krr, \
                                'meandxdx' + name_append : meandxdx, \
                                'meandxdy' + name_append : meandxdy, \
                                'meandydy' + name_append : meandydy, \
                                'meandrdr' + name_append : meandrdr, \
                                'meanKxx' + name_append : meanKxx, \
                                'meanKxy' + name_append : meanKxy, \
                                'meanKyy' + name_append : meanKyy, \
                                'meanKrr' + name_append : meanKrr})

              return dsO #}}}

            dsabs = compute_mixing(ds, dx, dy, '_abs_' + maskcond)
            dsmasks = dsmasks.merge(dsabs)
            #}}}

            # full relative diffusivity #{{{
            cx = xParticle.mean('nParticles', skipna=True)
            cy = yParticle.mean('nParticles', skipna=True)

            dx = xParticle - cx
            dy = yParticle - cy

            dsrel = compute_mixing(ds, dx, dy, '_rel_' + maskcond)
            dsmasks = dsmasks.merge(dsrel)
            #}}}

            # compute relative = absolute diffusivity crossing point #{{{
            tens = np.asarray(ds.elapsedtime.groupby(ds.elapsedtime.sel(Time=ds.Time)).mean(), dtype='f8')/1e9
            def get_intersection(tensT, diffabsT, diffrelT, label, maskcondT): #{{{

              # first value is a nan for meanKrr etc
              tcross, Kcross = first_line_intersection(tensT[1:], \
                  diffabsT[label+'_abs_' + maskcondT].values[1:], diffrelT[label+'_rel_' + maskcondT].values[1:])

              clustercoord = [('cluster' , np.array([int(acell)]))]
              ds0 = xr.Dataset().assign(**{'crossingtime_' + label + '_' + maskcond : \
                                          xr.DataArray(tcross, coords=clustercoord)})\
                                .assign(**{'crossingK_' + label + '_' + maskcond : \
                                          xr.DataArray(Kcross, coords=clustercoord)})

              return ds0 #}}}

            for alabel in ['meanKrr']:
              print 'Computing for ', alabel
              dsmasks.merge(get_intersection(tens, dsabs, dsrel, alabel, maskcond))
            #}}}

            # velocity #{{{
            # must extend the full time dimension for this calculation, possible bottleneck or memory error location
            u = xParticle.chunk({'Time':len(xParticle.Time)}).diff('Time',label='lower')
            v = yParticle.chunk({'Time':len(yParticle.Time)}).diff('Time',label='lower')

            um = u.mean('nParticles', skipna=skipna)
            vm = v.mean('nParticles', skipna=skipna)

            up = np.sqrt((u - um).mean('nParticles', skipna=skipna)**2.0)
            vp = np.sqrt((v - vm).mean('nParticles', skipna=skipna)**2.0)

            # use nanmean over ensembles
            meanum = um.groupby(ds.elapsedtime.sel(Time=ds.Time[1:])).mean(skipna=True)
            meanvm = vm.groupby(ds.elapsedtime.sel(Time=ds.Time[1:])).mean(skipna=True)
            meanup = up.groupby(ds.elapsedtime.sel(Time=ds.Time[1:])).mean(skipna=True)
            meanvp = vp.groupby(ds.elapsedtime.sel(Time=ds.Time[1:])).mean(skipna=True)

            def compute_rho(velp, dsA): #{{{

              velp0 = velp.sel(Time=dsA.starttime[:-1])
              velp0['Time'] = velp.Time.values

              rho = (velp0*velp).groupby(dsA.elapsedtime.sel(Time=dsA.Time[:-1])).mean(skipna=True)\
                  /(velp0*velp0).groupby(dsA.elapsedtime.sel(Time=dsA.Time[:-1])).mean(skipna=True)

              # potential for last several values being bad because of bad vel estimates
              rho.values[-2:] = np.nan

              return rho #}}}

            def compute_TL(rho, dsA): #{{{

              dt = dsA.deltaTimedays.groupby(dsA.elapsedtime.sel('Time')).mean()
              dt.values[0] = dt.values[1]
              #daystosec = 24*60*60.
              #TL = np.cumsum((rho*dt).values*daystosec)
              # TL in days
              TL = np.cumsum((rho*dt).values)

              return TL #}}}

            rhou = compute_rho(up, ds)
            rhov = compute_rho(vp, ds)

            TLx = compute_TL(rhou, ds)
            TLy = compute_TL(rhov, ds)
            TL = 0.5*(TLx + TLy)

            dsmasks = dsmasks.merge(xr.Dataset({
                                              'meanum' + '_' + maskcond : meanum, \
                                              'meanvm' + '_' + maskcond : meanvm, \
                                              'meanup' + '_' + maskcond : meanup, \
                                              'meanvp' + '_' + maskcond : meanvp, \
                                              'rhou'   + '_' + maskcond : rhou, \
                                              'rhov'   + '_' + maskcond : rhov, \
                                              'TLx'    + '_' + maskcond : TLx, \
                                              'TLy'    + '_' + maskcond : TLy, \
                                              'TL'     + '_' + maskcond : TL, \
                                              }))
            #{{{
                                              #'u' + '_' + maskcond : u, \
                                              #'v' + '_' + maskcond : v, \
                                              #'um' + '_' + maskcond : um, \
                                              #'vm' + '_' + maskcond : vm, \
                                              #'up' + '_' + maskcond : up, \
                                              #'vp' + '_' + maskcond : vp, \
            #}}}
            #}}}

            # mesh relative diffusivity #{{{
            xp1 = xParticle.isel(nParticles=meshedges[:,0])
            xp1['nParticles'] = np.arange(meshedges.shape[0])

            yp1 = yParticle.isel(nParticles=meshedges[:,0])
            yp1['nParticles'] = np.arange(meshedges.shape[0])

            xp2 = xParticle.isel(nParticles=meshedges[:,1])
            xp2['nParticles'] = np.arange(meshedges.shape[0])

            yp2 = yParticle.isel(nParticles=meshedges[:,1])
            yp2['nParticles'] = np.arange(meshedges.shape[0])

            dx =  xp1 - xp2
            dy =  yp1 - yp2

            dsmasks = dsmasks.merge(compute_mixing(ds, dx, dy, '_mesh_' + maskcond, pairs=True))
            #}}}

            # ring relative diffusivity #{{{
            xp1 = xParticle.isel(nParticles=ringedges[:,0])
            xp1['nParticles'] = np.arange(ringedges.shape[0])

            yp1 = yParticle.isel(nParticles=ringedges[:,0])
            yp1['nParticles'] = np.arange(ringedges.shape[0])

            xp2 = xParticle.isel(nParticles=ringedges[:,1])
            xp2['nParticles'] = np.arange(ringedges.shape[0])

            yp2 = yParticle.isel(nParticles=ringedges[:,1])
            yp2['nParticles'] = np.arange(ringedges.shape[0])

            dx =  xp1 - xp2
            dy =  yp1 - yp2

            dsmasks = dsmasks.merge(compute_mixing(ds, dx, dy, '_ring_' + maskcond, pairs=True))
            #}}}

          layerlist.append(dsmasks)

        # aggregate all the different layers
        byID = pd.Index(ds.buoyancySurfaceValues.values[0,:],name='buoyancysurface')
        dscluster = dscluster.merge(xr.concat(layerlist, byID))
        outname = "mixing_analysis_cluster_{0:08d}.nc".format(acell)
        print 'Outputing data for acell=%d to file %s'%(acell,outname)
        dscluster.to_netcdf(outname)

        break

    with timeit_context('Closing xr files'):
      dscluster.close()
      for ads in layerlist:
        ads.close()
      ds.close()
      particles.close()
      dsmask.close()
      mesh.close()

    return  #}}}

if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-p", "--particles", dest="particlefile",
                      help="file/path for analysis",
                      metavar="FILE")
    parser.add_option("-m", "--masks", dest="maskfile",
                      help="maskfile computed with build_masks.py",
                      metavar="FILE")
    parser.add_option("-g", "--grid", dest="meshfile",
                      help="mesh.nc file",
                      metavar="FILE")
    parser.add_option("-r", "--radius", dest="radius",
                      help="radius to form clusters",
                      metavar="FLOAT")
    parser.add_option("-c", "--clusterid", dest="clusterid",
                      help="cluster id, e.g., 'np.arange(11,20)'",
                      metavar="INT")

    options, args = parser.parse_args()

    if not options.particlefile:
        assert False, 'Must specify a particle file!'
    print 'particle file name is ', options.particlefile
    if not options.maskfile:
        assert False, 'Must specify a mask file!'
    print 'mask file name is ', options.maskfile
    if not options.meshfile:
        assert False, 'Must specify a mesh file!'
    print 'mesh file name is ', options.meshfile
    if not options.radius:
      options.radius = 100.e3
    else:
      options.radius = float(options.radius)
    print 'Radius for clusters is %d km'%(options.radius/1e3)
    if not options.clusterid:
      options.clusterid = np.array([17250])
    else:
      option.radius = eval(option.clusterid)

    xarray_analysis(options.particlefile, options.maskfile, options.meshfile, options.radius, options.clusterid)
