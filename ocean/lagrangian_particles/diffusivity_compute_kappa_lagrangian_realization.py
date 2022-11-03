#!/usr/bin/env python
"""

    Compute single particle statistics for a realization.
    
    Phillip Wolfram
    LANL
    12/02/2014

"""
# import libraries / packages
import os
import multiprocessing as mp
from functools import partial
import numpy as np
import netCDF4

# function definitions
def get_ensemble_paths(rootdir, folder_prefix): #{{{
    fnames = []
    for root, dirs, files in os.walk(rootdir):
      if(root==rootdir):
        dirs.sort()
        for adir in dirs:
          if(adir.startswith(folder_prefix)):
            fnames.append(rootdir + '/' + adir)

    return fnames #}}}

def check_folder(folder): #{{{
  if not os.path.isdir(folder):
      os.mkdir(folder)
      return #}}}

def diffusivity(sigma, ts, te, tint):
  return 0.5*(sigma[te,:] - sigma[ts,:])/(float(te - ts)*float(tint)*24.*60.*60.)

def compute_diffusivity(root, prefix, ts, dt, tint):
  
  ensembles = get_ensemble_paths(root, prefix)
  ensembles.reverse()

  # loop over all files and aggregate data
  for aensemble in ensembles:
    print 'processing %s ' % (aensemble)
    data = np.load(aensemble+'/lagrangian_data.npz')
    print 'done'
    mux = data['mux']
    muy = data['muy']
    npart = data['Npart']
    sxx = data['dxdx_sum']/(npart-1)
    sxy = data['dxdy_sum']/(npart-1)
    syy = data['dydy_sum']/(npart-1)
    srr = data['drdr_sum']/(npart-1)


    for ats, adt, atint in zip(ts, dt, tint):
      print 'compute diffusivity'
      # ending time (for slice requires +1)
      ate = ats + adt + 1
      
      ## for moving cooordinate (could also slice from ats to ats + dt + 1) #{{{
      #mux_out = mux[ats:ate,:]
      #muy_out = muy[ats:ate,:] #}}}
      # starting coordinate is "base map"
      #mux_out = mux[0,:]
      #muy_out = muy[0,:]
      # starting point
      mux_out = mux[ats,:]
      muy_out = muy[ats,:]

      # just want starting value with slope...
      sxx_ts = sxx[ats,:]
      sxy_ts = sxy[ats,:]
      syy_ts = syy[ats,:]
      srr_ts = srr[ats,:]

      kappa_xx = diffusivity(sxx, ats, ats + adt, atint)
      kappa_xy = diffusivity(sxy, ats, ats + adt, atint)
      kappa_yy = diffusivity(syy, ats, ats + adt, atint)
      kappa_rr = diffusivity(srr, ats, ats + adt, atint)

      mux_ts = mux[ats,:]
      muy_ts = muy[ats,:]
      mux_te = mux[ats+adt,:]
      muy_te = muy[ats+adt,:]

      savename ='/realization_kappa-' + str(ats) + '-' + str(adt) + '-' + str(atint) 
      print 'saving %s in %s' % (savename, aensemble)
      np.savez(aensemble+savename + '.npz', mux_base=mux_out, muy_base=muy_out,
              mux_ts=mux_ts, muy_ts=muy_ts, mux_te=mux_te, muy_te=muy_te,
              mux_full=mux, muy_full=muy, sxx_ts=sxx_ts, sxy_ts=sxy_ts,
              syy_ts=syy_ts, srr_ts=srr_ts, kappa_xx=kappa_xx,
              kappa_xy=kappa_xy, kappa_yy=kappa_yy, kappa_rr=kappa_rr, 
              sxx_full = sxx, sxy_full=sxy, syy_full=syy, srr_full=srr)

  return 
 
if __name__ == "__main__":
  from optparse import OptionParser

  # Get command line parameters #{{{
  parser = OptionParser()
  parser.add_option("-r", "--root", dest="root",
      help="folder root holding analysis folders",
      metavar="FILE")
  parser.add_option("-p", "--prefix", dest="prefix",
      help="folder prefix for analysis",
      metavar="FILE")

  options, args = parser.parse_args()

  if not options.root:
    options.root = '.'
  if not options.prefix:
    options.prefix = "analyze_restart"
  #}}}

  ts = [5] #np.array([5, 5, 5, 5, 5, 15])
  dt = [1] #np.array([1, 2, 4, 8, 16, 15])
  tint = [2] #np.array([2, 2, 2, 2, 2, 2])

  compute_diffusivity(options.root, options.prefix, ts, dt, tint)
