#!/usr/bin/env python

# import libraries / packages
import os
import numpy as np
import netCDF4
from GeneralTriangulation import GeneralTriangulation

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

def error_estimates(root, prefix, case, nlayers):
  
  ensembles = get_ensemble_paths(root, prefix)
  nensembles = len(ensembles)

  # build up arrays for data analysis
  for aens in np.arange(nensembles):
    data = np.load(ensembles[aens] + '/' + case + '.npz')
    print 'on ensemble %d of %d' % (aens, nensembles)
    kappa_rr = data['kappa_rr']

    if aens == 0:
      x = data['mux_base']
      y = data['muy_base']
      
      npart_layers = x.shape[0]/nlayers

      xlayer = np.zeros((npart_layers, nlayers, nensembles))
      ylayer = np.zeros((npart_layers, nlayers, nensembles))
      kappalayer = np.zeros((npart_layers, nlayers, nensembles))
      std_to_mean = np.zeros((nlayers, nensembles))
      kappa_mean_change = np.zeros((npart_layers, nlayers, nensembles))
      kappa_std_change = np.zeros((npart_layers, nlayers, nensembles))
      kappa_mean = np.zeros((npart_layers, nlayers, nensembles))
      kappa_std = np.zeros((npart_layers, nlayers, nensembles))
      
    triang = []
    for alayer in np.arange(nlayers):
      xlayer[:,alayer, aens] = np.degrees(x[alayer*npart_layers:(alayer+1)*npart_layers])
      ylayer[:,alayer, aens] = np.degrees(y[alayer*npart_layers:(alayer+1)*npart_layers])
      triang.append(GeneralTriangulation(xlayer[:,alayer,aens], ylayer[:,alayer,aens]))

    for alayer in np.arange(nlayers):
      #kappalayer[:,alayer,aens] = kappa_rr[alayer*npart_layers:(alayer+1)*npart_layers]
      kappalayer[:,alayer,aens] = \
              triang[alayer].smooth_laplacian(kappa_rr[alayer*npart_layers:(alayer+1)*npart_layers], ntimes=10)
      kappa_std[:, alayer, aens] = np.std(kappalayer[:,alayer,:(aens+1)],axis=1)
      kappa_mean[:, alayer, aens] = np.mean(kappalayer[:,alayer,:(aens+1)],axis=1)
    
      basename = root + '/' + 'ensemble' + "%.4d" % (aens) + '_layer' + "%.4d" % (alayer)
      #plot_stats(basename, triang[alayer], kappa_std, kappa_mean)
      print 'done with layer %d for ensemble %d' % (alayer, aens)

  print 'saving...'
  np.savez('realization_error_estimates_data.npz',xlayer=xlayer, ylayer=ylayer, 
      kappa_std=kappa_std, kappa_mean=kappa_mean, kappalayer=kappalayer)
  print 'done'
 
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
  parser.add_option("-c", "--case", dest="case",
      help="case, prefix for *.npz file",
      metavar="STRING")
  parser.add_option("-l", "--numlayers", dest="nlayers",
      help="number of layers",
      metavar="INT")

  options, args = parser.parse_args()

  if not options.root:
    options.root = '.'
  if not options.prefix:
    options.prefix = "analyze_restart"
  if not options.case:
    options.case = "realization_kappa-5-1-2"
  if not options.nlayers:
    options.nlayers = 5
  
  #}}}

  error_estimates(options.root, options.prefix, options.case, options.nlayers)

