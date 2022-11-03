#!/usr/bin/env python

# initialize general libraries and setup notebook
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from multiprocessing.dummy import Pool as ThreadPool
import cluster_topology as CT
import GeneralTriangulation as GT
import effective_diffusivity as ED
from iotasks import savenpz, timeit_context



def split_array(array, npieces): #{{{
    nvalues = np.floor(len(array)/npieces)
    split = []
    extra = 0
    for i in np.arange(npieces):
        if i == (npieces-1):
            extra = len(array) - nvalues*npieces
        split.append(array[i*nvalues:(((i+1)*nvalues+extra))])
    return split #}}}

def run(procnum, nprocs, rlzn):
  print 'Computing effective diffusivity for %d of %d processors for realization %d'%(procnum, nprocs, rlzn)

  # gaussian properties
  #yloc = 1500e3
  lgauss = 50e3

  # mesh name
  mname = 'mesh.nc'

  with timeit_context('Loading data'):
      ds = np.load('all_rlzn_particle_data_rads.npz')
      x = ds['lonrad']
      y = ds['latrad']
      mask = ds['notoutcropped']
      dt = ds['dtdays']
      #print dt.shape, x.shape
      tvec = np.zeros(dt.shape[1]+1)
      tvec[1:] = np.cumsum(dt[rlzn,:]*86400.)

  # {{{
  ## target buoyancy surface
  #buoySurf = 1029.7

  #buoyancySurface = ds['buoyancySurface']
  #alayer = np.where(buoyancySurface == buoySurf)
  #print 'Buoyancy Surface of interest: ', buoyancySurface[alayer], 'among ', buoyancySurface
  #}}}

  # total work
  ytot = np.linspace(0,2000e3,201)

  # note this is a shared-memory (threads) implementation.  This should be ok because we really only need one copy of x/y
  def run_eff_diff(alayerID): #{{{
    # divide spatial work into chunks
      y0 = split_array(ytot, nprocs)[procnum]
      print 'Computing for y0=%s, layer=%s, rlzn=%s'%(y0,alayerID,rlzn)
      # compute effective diffusivity for the spatial decomp, each layer, and the split spatial dimension
      effdiff = ED.EffectiveDiffusivity(\
          xt=np.squeeze(x[rlzn,:,alayerID,:]), yt=np.squeeze(y[rlzn,:,alayerID,:]), \
          mask=mask[rlzn,0,alayerID,:], \
          meshname=mname, calctype='ygauss', \
          y0=y0, concwidth=lgauss, \
          tparticles=tvec[:], advection=True, diffusion=True, \
          kappa=10., dt=60., dtout=1.*86400., \
          run=True)
      out = effdiff.output()
      # append driver metadata
      out['procnum'] = procnum
      out['nprocs'] = nprocs
      out['rlzn'] = rlzn

      savenpz('eff_diff_output/effective_diffusivity_rlzn_%03d_split_%04d_results_layer%02d.npz'%(rlzn, procnum, alayerID), **out)

      #return effdiff
      return #}}}

  ## run all the layers at once
  #with timeit_context('Running job from thread pool'):
  #    pool = ThreadPool(2)
  #    results = pool.map(run_eff_diff, np.arange(11))
  #    pool.close()
  #    pool.join()

  # super IO limited on mustang, must just run serially
  #for i in np.arange(11)[::-1][:2]:
  for i in np.arange(11):
      run_eff_diff(i)

if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-p", "--proc", dest="procnum",
            help="number of processor",
            metavar="INT")
    parser.add_option("-n", "--nprocs", dest="nprocs",
            help="number of total processors for spatial decomposition",
            metavar="INT")
    parser.add_option("-r", "--rlzn", dest="rlzn", help="rlzn number", metavar="INT")

    options, args = parser.parse_args()

    def ensure_input_ok(value, name):
        if not value:
            assert False, name + ' is required'
        value = int(value)
        print name + ' is ', value
        return value

    options.procnum = ensure_input_ok(options.procnum,'Processor number')
    options.nprocs  = ensure_input_ok(options.nprocs, 'Total number of processors')
    options.rlzn = ensure_input_ok(options.rlzn, 'Number of realizations')

    run(options.procnum, options.nprocs, options.rlzn)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
