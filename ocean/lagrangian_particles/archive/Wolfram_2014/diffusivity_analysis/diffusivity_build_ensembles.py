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
            fnames.append(rootdir + '/' + adir + '/particles_statistics/')

    return fnames #}}}

def check_folder(folder):
  if not os.path.isdir(folder):
      os.mkdir(folder)
  
def output_data(num, xinterp, yinterp, clstr_num, clstr_meanx, clstr_meany, kappa_rr, kappa_xx, kappa_xy, kappa_yy, sig2_rr, sig2_xx, sig2_xy, sig2_yy, nlayers, ninterp_layers, npoints_layers): #{{{
  ensemble_folder = 'particle_statistics'+'/ensemble'+"{0:03d}".format(num)+'/'
  check_folder(ensemble_folder)
  for j in np.arange(nlayers):
    xil = xinterp[ninterp_layers*j:ninterp_layers*(j+1)] 
    yil = yinterp[ninterp_layers*j:ninterp_layers*(j+1)] 
    np.savetxt(ensemble_folder+'interp_point_layer'+"{0:03d}".format(j)+'.txt', np.vstack((xil,yil)).T, fmt=['%e','%e'],header=str(xil.shape[0]),comments='')
   
    nl = clstr_num[j,:(num+1)*npoints_layers].ravel()
    xl = clstr_meanx[j,:(num+1)*npoints_layers].ravel()
    yl = clstr_meany[j,:(num+1)*npoints_layers].ravel()
    krrl = kappa_rr[j,:(num+1)*npoints_layers].ravel()
    kxxl = kappa_xx[j,:(num+1)*npoints_layers].ravel()
    kxyl = kappa_xy[j,:(num+1)*npoints_layers].ravel()
    kyyl = kappa_yy[j,:(num+1)*npoints_layers].ravel()
    s2xx = sig2_xx[j,:(num+1)*npoints_layers].ravel()
    s2xy = sig2_xy[j,:(num+1)*npoints_layers].ravel()
    s2yy = sig2_yy[j,:(num+1)*npoints_layers].ravel()
    s2rr = sig2_rr[j,:(num+1)*npoints_layers].ravel()

    fmtstr = ['%e' for i in np.arange(11)]
    fmtstr[0] = '%d'
    tupledata = np.vstack((nl, xl, yl, krrl, kxxl, kxyl, kyyl, s2xx, s2xy, s2yy, s2rr)).T
    np.savetxt(ensemble_folder+'cluster_data_layer'+"{0:03d}".format(j)+'.txt', tupledata, fmt=fmtstr, header=str(nl.shape[0]), comments='')
  return num #}}}

def build_ensembles(root, prefix, nlayers, outrange):
  """
  Build up ensembles from each realization
  Phillip Wolfram
  LANL
  09/29/2014
  """

  file_names = get_realization_paths(root, prefix)

  num_rlzns = len(file_names)

  for num, folder in enumerate(file_names):
    print 'loading %s-th file in %s' % (str(num),folder)
    n = np.load(folder+'/n.npy')
    x = np.load(folder+'/x.npy')
    y = np.load(folder+'/y.npy')
    krr = np.load(folder+'/krr.npy')
    kxx = np.load(folder+'/kxx.npy')
    kxy = np.load(folder+'/kxy.npy')
    kyy = np.load(folder+'/kyy.npy')
    s2xx = np.load(folder+'/s2xx.npy')
    s2xy = np.load(folder+'/s2xy.npy')
    s2yy = np.load(folder+'/s2yy.npy')
    s2rr = np.load(folder+'/s2rr.npy')

    if (num == 0):
      xinterp = np.load(folder+'/xi.npy')
      yinterp = np.load(folder+'/yi.npy')
      
      Ninterp_points = xinterp.shape[0]
      ninterp_layers = Ninterp_points/nlayers
      
      Npoints = n.shape[0]
      if (Npoints % nlayers):
        print 'possible error with layer numbers'
      npoints_layers = Npoints/nlayers
      
      clstr_num   = np.zeros((nlayers,npoints_layers*num_rlzns))
      clstr_meanx = np.zeros((nlayers,npoints_layers*num_rlzns))
      clstr_meany = np.zeros((nlayers,npoints_layers*num_rlzns))
      kappa_rr    = np.zeros((nlayers,npoints_layers*num_rlzns))
      kappa_xx    = np.zeros((nlayers,npoints_layers*num_rlzns))
      kappa_xy    = np.zeros((nlayers,npoints_layers*num_rlzns))
      kappa_yy    = np.zeros((nlayers,npoints_layers*num_rlzns))
      sig2_xx     = np.zeros((nlayers,npoints_layers*num_rlzns))
      sig2_xy     = np.zeros((nlayers,npoints_layers*num_rlzns))
      sig2_yy     = np.zeros((nlayers,npoints_layers*num_rlzns))
      sig2_rr     = np.zeros((nlayers,npoints_layers*num_rlzns))

    for j in np.arange(nlayers):
      clstr_num[j,num*npoints_layers:(num+1)*npoints_layers] = n[j*npoints_layers:(j+1)*npoints_layers]
      clstr_meanx[j,num*npoints_layers:(num+1)*npoints_layers] = x[j*npoints_layers:(j+1)*npoints_layers]
      clstr_meany[j,num*npoints_layers:(num+1)*npoints_layers] = y[j*npoints_layers:(j+1)*npoints_layers]
      kappa_rr[j,num*npoints_layers:(num+1)*npoints_layers] = krr[j*npoints_layers:(j+1)*npoints_layers]
      kappa_xx[j,num*npoints_layers:(num+1)*npoints_layers] = kxx[j*npoints_layers:(j+1)*npoints_layers]
      kappa_xy[j,num*npoints_layers:(num+1)*npoints_layers] = kxy[j*npoints_layers:(j+1)*npoints_layers]
      kappa_yy[j,num*npoints_layers:(num+1)*npoints_layers] = kyy[j*npoints_layers:(j+1)*npoints_layers]
      sig2_rr[j,num*npoints_layers:(num+1)*npoints_layers] = s2rr[j*npoints_layers:(j+1)*npoints_layers]
      sig2_xx[j,num*npoints_layers:(num+1)*npoints_layers] = s2xx[j*npoints_layers:(j+1)*npoints_layers]
      sig2_xy[j,num*npoints_layers:(num+1)*npoints_layers] = s2xy[j*npoints_layers:(j+1)*npoints_layers]
      sig2_yy[j,num*npoints_layers:(num+1)*npoints_layers] = s2yy[j*npoints_layers:(j+1)*npoints_layers]

  # build ensemble data sets
  check_folder(root +'/particle_statistics/')
  for num in outrange:
    ensemble_folder = root +'/particle_statistics'+'/ensemble'+"{0:03d}".format(num)+'/'
    check_folder(ensemble_folder)
    for j in np.arange(nlayers):
      xil = xinterp[ninterp_layers*j:ninterp_layers*(j+1)] 
      yil = yinterp[ninterp_layers*j:ninterp_layers*(j+1)] 
      np.savetxt(ensemble_folder+'interp_point_layer'+"{0:03d}".format(j)+'.txt', np.vstack((xil,yil)).T, fmt=['%e','%e'],header=str(xil.shape[0]),comments='')
     
      nl = clstr_num[j,:(num+1)*npoints_layers].ravel()
      xl = clstr_meanx[j,:(num+1)*npoints_layers].ravel()
      yl = clstr_meany[j,:(num+1)*npoints_layers].ravel()
      krrl = kappa_rr[j,:(num+1)*npoints_layers].ravel()
      kxxl = kappa_xx[j,:(num+1)*npoints_layers].ravel()
      kxyl = kappa_xy[j,:(num+1)*npoints_layers].ravel()
      kyyl = kappa_yy[j,:(num+1)*npoints_layers].ravel()
      s2xx = sig2_xx[j,:(num+1)*npoints_layers].ravel()
      s2xy = sig2_xy[j,:(num+1)*npoints_layers].ravel()
      s2yy = sig2_yy[j,:(num+1)*npoints_layers].ravel()
      s2rr = sig2_rr[j,:(num+1)*npoints_layers].ravel()

      fmtstr = ['%e' for i in np.arange(11)]
      fmtstr[0] = '%d'
      tupledata = np.vstack((nl, xl, yl, krrl, kxxl, kxyl, kyyl, s2xx, s2xy, s2yy, s2rr)).T
      np.savetxt(ensemble_folder+'cluster_data_layer'+"{0:03d}".format(j)+'.txt', tupledata, fmt=fmtstr, header=str(nl.shape[0]), comments='')
 
  # attempt to parallelize output (limitation of individual machine)  #{{{
  ## choice of pool is based on max memory of machine
  #pool = mp.Pool(2)
  ## see http://timothyawiseman.wordpress.com/2012/12/21/a-really-simple-multiprocessing-python-example/ 
  ## for good tutorial on how to make this happen via partial for iteratable over single variable for map
  #output = partial(output_data, xinterp=xinterp, yinterp=yinterp, clstr_num=clstr_num, clstr_meanx=clstr_meanx, clstr_meany=clstr_meany, \
  #    kappa_rr=kappa_rr, kappa_xx=kappa_xx, kappa_xy=kappa_xy, kappa_yy=kappa_yy, \
  #    sig2_rr=sig2_rr, sig2_xx=sig2_xx, sig2_xy=sig2_xy, sig2_yy=sig2_yy, \
  #    nlayers=nlayers, ninterp_layers=ninterp_layers, npoints_layers=npoints_layers)
  ## just flip xrange to do hardest ones first (reversed xrange(num_rlzns))
  ##nums = pool.map(output, xrange(num_rlzns-1,-1,-1))
  #nums = pool.map(output, xrange(num_rlzns))
  #pool.close()
  #pool.join()
  #
  ##output_data(num_rlzns-1, xinterp=xinterp, yinterp=yinterp, clstr_num=clstr_num, clstr_meanx=clstr_meanx, clstr_meany=clstr_meany, \
  ##    kappa_rr=kappa_rr, kappa_xx=kappa_xx, kappa_xy=kappa_xy, kappa_yy=kappa_yy, \
  ##    sig2_rr=sig2_rr, sig2_xx=sig2_xx, sig2_xy=sig2_xy, sig2_yy=sig2_yy, \
  ##    nlayers=nlayers, ninterp_layers=ninterp_layers, npoints_layers=npoints_layers)

  ##for num in np.arange(num_rlzns):
  ##  p = mp.Process(target=output_data, args = (num,))
  ##  p.start()
  ##[p.join() for p in mp.active_children()]
  #}}}


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
  parser.add_option("-l", "--layers", dest="layers",
      help="number of particle layers",
      metavar="INTEGER")
  parser.add_option("-s", "--start", dest="start",
      help="starting number of ensemble",
      metavar="INTEGER")
  parser.add_option("-e", "--end", dest="end",
      help="ending number of ensemble",
      metavar="INTEGER")

  options, args = parser.parse_args()

  if not options.root:
    parser.error("Root directory is a required input.")
  if not options.prefix:
    options.prefix = "analyze_output"
  if options.layers:
    options.layers = int(options.layers)
  if options.start:
    options.start = int(options.start)
  if options.end:
    options.end = int(options.end)

  #}}}

  outrange = np.unique(np.array([options.start,options.end]))
  build_ensembles(options.root, options.prefix, options.layers, outrange)
