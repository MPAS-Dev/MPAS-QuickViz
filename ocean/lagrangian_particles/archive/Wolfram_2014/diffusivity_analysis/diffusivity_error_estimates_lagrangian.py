#!/usr/bin/env python

# import libraries / packages
import os
import numpy as np
import netCDF4

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
from GeneralTriangulation import GeneralTriangulation
from convert_ParaView_xml_to_matplotlib_colormap import convert_ParaView_xml_to_matplotlib_colormap

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

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

def plot_stats(basename, triang, std, mean):

  def basic_plot(var, varname):
    plt.figure() 
    triang.plot_scalar(var, cmap=convert_ParaView_xml_to_matplotlib_colormap('redblue.xml'))
    plt.title(varname)
    cmax = np.max(np.abs(var))
    plt.clim([-cmax, cmax])
    plt.colorbar()
    filename = basename + '_'+ varname + '.png'
    print 'saving %s' % filename
    plt.savefig(filename)

  basic_plot(mean, 'mean')
  basic_plot(std, 'std')
  basic_plot(std/mean, 'std_over_mean')

  plt.close('all')


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

      xlayer = np.zeros((npart_layers, nlayers))
      ylayer = np.zeros((npart_layers, nlayers))
      kappalayer = np.zeros((npart_layers, nlayers, nensembles))
      std_to_mean = np.zeros((nlayers, nensembles))
      kappa_mean_change = np.zeros((npart_layers, nlayers, nensembles))
      kappa_std_change = np.zeros((npart_layers, nlayers, nensembles))
      kappa_mean = np.zeros((npart_layers, nlayers, nensembles))
      kappa_std = np.zeros((npart_layers, nlayers, nensembles))
      
      triang = []
      for alayer in np.arange(nlayers):
        xlayer[:,alayer] = np.degrees(x[alayer*npart_layers:(alayer+1)*npart_layers])
        ylayer[:,alayer] = np.degrees(y[alayer*npart_layers:(alayer+1)*npart_layers])
        triang.append(GeneralTriangulation(xlayer[:,alayer], ylayer[:,alayer]))

    for alayer in np.arange(nlayers):
      kappalayer[:,alayer,aens] = triang[alayer].smooth_laplacian(kappa_rr[alayer*npart_layers:(alayer+1)*npart_layers], ntimes=10)
      #kappalayer[:,alayer,aens] = kappa_rr[alayer*npart_layers:(alayer+1)*npart_layers]
      kappa_std[:, alayer, aens] = np.std(kappalayer[:,alayer,:(aens+1)],axis=1)
      kappa_mean[:, alayer, aens] = np.mean(kappalayer[:,alayer,:(aens+1)],axis=1)
      if aens > 0:
        lastone = int(np.max([0,(aens-1)]))
        kappa_mean_change[:,alayer,aens] = (kappa_mean[:, alayer, aens] - kappa_mean[:, alayer,lastone]) / kappa_mean[:,alayer,lastone]
        kappa_std_change[:,alayer,aens] = (kappa_std[:, alayer, aens] - kappa_std[:,alayer,lastone]) / kappa_std[:,alayer,lastone]
      std_to_mean[alayer, aens] = np.nanmean(np.abs(kappa_std)) / np.nanmean(np.abs(kappa_mean))
    
      basename = root + '/' + 'ensemble' + "%.4d" % (aens) + '_layer' + "%.4d" % (alayer)
      #plot_stats(basename, triang[alayer], kappa_std, kappa_mean)
      print 'done with layer %d for ensemble %d' % (alayer, aens)

  np.savez('error_estimates_data.npz',xlayer=xlayer, ylayer=ylayer, std_to_mean=std_to_mean, kappa_mean_change=kappa_mean_change, kappa_std_change=kappa_std_change, \
      kappa_std=kappa_std, kappa_mean=kappa_mean, kappalayer=kappalayer)

  # plot overall stats now
  plt.figure()
  for alayer in np.arange(nlayers):
    plt.plot(std_to_mean[alayer,:], label='layer ' + str(alayer))
  #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.legend(loc=4)
  plt.ylabel('mean(|std_kappa|) / mean(|mean_kappa|)')
  plt.xlabel('Num ensembles')
  plt.savefig(root + '/ensemble_std_to_mean_kappa_comparision.png')
  plt.close()
  
  plt.figure()
  for alayer in np.arange(nlayers):
    plt.semilogy(np.nanmax(np.abs(kappa_mean_change[:,alayer,:])), label='layer ' + str(alayer))
  #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.legend()
  #plt.ylabel('mean | (kappa_mean_{i+1} - kappa_mean{i})/ kappa_mean_{i} | ')
  plt.ylabel('Max relative change in mean $\kappa$')
  plt.xlabel('Num ensembles')
  plt.savefig(root + '/ensemble_kappa_mean_change_comparision.png')
  plt.close()
  
  plt.figure()
  for alayer in np.arange(nlayers):
    plt.semilogy(np.nanmax(np.abs(kappa_std_change[:,alayer,:])), label='layer ' + str(alayer))
  #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.legend()
  #plt.ylabel('mean | (kappa_std_{i+1} - kappa_std{i})/ kappa_std_{i} | ')
  plt.ylabel('Max relative change in standard deviation of $\kappa$')
  plt.xlabel('Num ensembles')
  plt.savefig(root + '/ensemble_kappa_std_change_comparision.png')
  plt.close()

 
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
    parser.error("Root directory is a required input.")
  if not options.prefix:
    options.prefix = "analyze_restart"
  if not options.case:
    options.case = "kappa-5-1-2"
  if not options.nlayers:
    options.nlayers = 5
  
  #}}}

  error_estimates(options.root, options.prefix, options.case, options.nlayers)

