#!/usr/bin/env python

# import libraries / packages
import os
import numpy as np
import netCDF4

# function definitions
def get_realization_paths(rootdir, folder_prefix): #{{{
    """
    Get paths for realization folders
    Phillip Wolfram, LANL
    09/19/2014
    """
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

def build_ensembles(root, prefix):
  """
  Build up ensembles from each realization
  Phillip Wolfram
  LANL
  09/29/2014
  """

  # get ready to build ensemble data sets
  check_folder(root +'/particle_ensembles/')
  
  file_names = get_realization_paths(root, prefix)

  num_rlzns = len(file_names)

  npart_sum = 0.0
  mux_sum  = 0.0
  muy_sum  = 0.0
  dxdx_sum = 0.0 
  dxdy_sum = 0.0 
  dydy_sum = 0.0 
  drdr_sum = 0.0 

  # loop over all files and aggregate data
  for num, folder in enumerate(file_names):
    print 'loading %s-th file in %s' % (str(num),folder)
    data = np.load(folder+'/lagrangian_data.npz')
    print 'done'
    npart = data['Npart']
    mux = data['mux']
    muy = data['muy']
    dxdx = data['dxdx_sum']
    dxdy = data['dxdy_sum']
    dydy = data['dydy_sum']
    drdr = data['drdr_sum']

    print 'aggregate data'
    npart_sum += npart
    mux_sum  += mux*npart
    muy_sum  += muy*npart
    dxdx_sum += dxdx*npart
    dxdy_sum += dxdy*npart
    dydy_sum += dydy*npart
    drdr_sum += drdr*npart

    print 'improved ensemble estimate'
    mux_ens = mux_sum/npart_sum
    muy_ens = muy_sum/npart_sum
    dxdx_ens = dxdx_sum/npart_sum
    dxdy_ens = dxdy_sum/npart_sum
    dydy_ens = dydy_sum/npart_sum
    drdr_ens = drdr_sum/npart_sum

    ensemble_folder = root +'/particle_ensembles'+'/ensemble'+"{0:03d}".format(num)+'/'
    print 'store ensemble values in %s' % (ensemble_folder)
    check_folder(ensemble_folder)
    np.savez(ensemble_folder+'ensemble.npz', mux=mux_ens, muy=muy_ens, dxdx=dxdx_ens, dxdy=dxdy_ens, dydy=dydy_ens, drdr=drdr_ens, npart=npart_sum)
 
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

  build_ensembles(options.root, options.prefix)
