#!/usr/bin/env python
"""

    Compute single particle statistics accross realizations.
    
    Phillip Wolfram
    LANL
    09/19/2014

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
    data = np.load(aensemble+'/ensemble.npz')
    print 'done'
    mux = data['mux']
    muy = data['muy']
    dxdx = data['dxdx']
    dxdy = data['dxdy']
    dydy = data['dydy']
    drdr = data['drdr']


    for ats, adt, atint in zip(ts, dt, tint):
      print 'compute diffusivity'
      # ending time (for slice requires +1)
      ate = ats + adt + 1
      
      ## for moving cooordinate (could also slice from ats to ats + dt + 1) #{{{
      #mux_out = mux[ats:ate,:]
      #muy_out = muy[ats:ate,:] #}}}
      # starting coordinate is "base map"
      mux_out = mux[0,:]
      muy_out = muy[0,:]

      ## could also take mean value because we are taking mean value for diffusivity #{{{
      #dxdx_out = np.mean(dxdx[ats:ate,:])
      #dxdy_out = np.mean(dxdy[ats:ate,:])
      #dydy_out = np.mean(dydy[ats:ate,:])
      #drdr_out = np.mean(drdr[ats:ate,:]) #}}}
      # just want starting value with slope...
      dxdx_ts = dxdx[ats,:]
      dxdy_ts = dxdy[ats,:]
      dydy_ts = dydy[ats,:]
      drdr_ts = drdr[ats,:]

      kappa_xx = diffusivity(dxdx, ats, ats + adt, atint)
      kappa_xy = diffusivity(dxdy, ats, ats + adt, atint)
      kappa_yy = diffusivity(dydy, ats, ats + adt, atint)
      kappa_rr = diffusivity(drdr, ats, ats + adt, atint)

      mux_ts = mux[ats,:]
      muy_ts = muy[ats,:]
      mux_te = mux[ats+adt,:]
      muy_te = muy[ats+adt,:]


      savename ='/kappa-' + str(ats) + '-' + str(adt) + '-' + str(atint) 
      print 'saving %s in %s' % (savename, aensemble)
      np.savez(aensemble+savename + '.npz', mux_base=mux_out, muy_base=muy_out, mux_ts=mux_ts, muy_ts=muy_ts, mux_te=mux_te, muy_te=muy_te, mux_full=mux, muy_full=muy, \
          dxdx_ts=dxdx_ts, dxdy_ts=dxdy_ts, dydy_ts=dydy_ts, drdr_ts=drdr_ts, kappa_xx=kappa_xx, kappa_xy=kappa_xy, kappa_yy=kappa_yy, kappa_rr=kappa_rr, \
          dxdx_full = dxdx, dxdy_full=dxdy, dydy_full=dydy, drdr_full=drdr)

      ## write ascii files for analysis in fortran
      #print 'writing ascii data'
      #with file(aensemble + savename + '_cluster_data.txt', 'w') as outfile:
      #  # want the total number of data points
      #  outfile.write('%d\n' % mux_out.shape[0])
      #  for n,(x,y,krr,kxx,kxy,kyy,srr,sxx,sxy,syy) in enumerate(zip(mux_out, muy_out, kappa_rr, kappa_xx, kappa_xy, kappa_yy, drdr_out, dxdx_out, dxdy_out, dydy_out)):
      #    outfile.write('%d %e %e %e %e %e %e %e %e %e %e\n' % (n,x,y,krr,kxx,kxy,kyy,srr,sxx,sxy,syy))
      #with file(aensemble + savename + '_interp_points.txt', 'w') as outfile:
      #  outfile.write('%d\n' % mux_out.shape[0])
      #  for x,y in zip(mux_out, muy_out):
      #    outfile.write('%e %e\n' % (x,y))
      #print 'done'

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
    parser.error("Root directory is a required input.")
  if not options.prefix:
    options.prefix = "analyze_output"
  #}}}

  ts = [5] #np.array([5, 5, 5, 5, 5, 15])
  dt = [1] #np.array([1, 2, 4, 8, 16, 15])
  tint = [2] #np.array([2, 2, 2, 2, 2, 2])

  compute_diffusivity(options.root, options.prefix, ts, dt, tint)
