#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
from matplotlib.tri import Triangulation
#from matplotlib.tri import UniformTriRefiner
from numpy import isnan
from numpy import degrees
from lxml import etree
import os
from matplotlib.colors import LinearSegmentedColormap

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':22})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
from matplotlib import rc

def convert_ParaView_xml_to_matplotlib_colormap(fname):
    
    xmlfile = etree.parse(fname)
    root = xmlfile.getroot()

    x = []
    r = []
    g = []
    b = []
    for colormap in root.findall('ColorMap'):
        for point in colormap.findall('Point'):
            x.append(float(point.get('x')))
            r.append(float(point.get('r')))
            g.append(float(point.get('g')))
            b.append(float(point.get('b')))

    # now we have all the data in lists and we can build up the colormap
    cmdata = dict()
    cmdata['red'] = []
    cmdata['green'] = []
    cmdata['blue'] = []
    for xi,ri,gi,bi in zip(x,r,g,b):
        cmdata['red'].append((xi, ri, ri))
        cmdata['green'].append((xi, gi, gi))
        cmdata['blue'].append((xi, bi, bi))
    cmdata['red'] = tuple(cmdata['red'])
    cmdata['green'] = tuple(cmdata['green'])
    cmdata['blue'] = tuple(cmdata['blue'])
    cmapname = os.path.splitext(fname)[0]
    cm = LinearSegmentedColormap(cmapname, cmdata)

    return cm

def plot_save_fig(data, fileprefix, titlename, cs=None, nrm=None): #{{{
    # mask out bad values
    mask = np.where(isnan(data))
    data[mask] = 0
    #mask = np.where(data < 0)
    #data[mask] = 0
   
    # plot the data
    plt.figure()
    plt.gca().set_aspect('equal')
    plt.hold(True)
    #plt.triplot(triang, lw=0.5, color='white')
    #plt.tricontourf(triang, kappaIntopenmp)
    usecmap = plt.cm.spectral
    usecmap = plt.cm.gist_ncar
    usecmap = plt.cm.hsv
    usecmap = plt.cm.jet
    usecmap = convert_ParaView_xml_to_matplotlib_colormap('../../rainbow_desaturated.xml')
    plt.tripcolor(triang, data, shading='flat', cmap=usecmap, norm=nrm)
    if cs is not None:
      if nrm is not None:
        plt.clim([100.0,100000.0])
      else: 
        plt.clim([0.0,cs])
    plt.colorbar()
    plt.xlabel('$Longitude~(^\circ)$')
    plt.ylabel('$Latitude~(^\circ)$')
    plt.title(titlename)
    plt.savefig(fileprefix+layer+'.png')
    return #}}}

if __name__ == "__main__":
    from optparse import OptionParser

    # get command line parameters
    parser = OptionParser()
    parser.add_option("-c", "--colorscale", dest="colorscale",
            help="clim max")
    parser.add_option("-i", "--interp", dest="interp",
            help="interpolation file")
    parser.add_option("-k", "--kappa", dest="kappa",
            help="kappa file")
    parser.add_option("-l", "--layer", dest="layer",
            help="layer for selection from buoyancy surface")
            

    options, args  = parser.parse_args()

    if options.colorscale:
      options.colorscale = float(options.colorscale)
    if options.layer:
    layer = options.layer

    pInt = np.loadtxt(options.interp,skiprows=1)
    x = pInt[:,0]
    y = pInt[:,1]
    kappaIntopenmp_tot = np.loadtxt(options.kappa,skiprows=1)
    kappaIntopenmp = kappaIntopenmp_tot[:,0]
    kappaIntopenmp_xx = kappaIntopenmp_tot[:,1]
    kappaIntopenmp_xy = kappaIntopenmp_tot[:,2]
    kappaIntopenmp_yy = kappaIntopenmp_tot[:,3]
    sigmaInt_xx = kappaIntopenmp_tot[:,4]
    sigmaInt_xy = kappaIntopenmp_tot[:,5]
    sigmaInt_yy = kappaIntopenmp_tot[:,6]
    sigmaInt_rr = kappaIntopenmp_tot[:,7]

    # build triangulation
    triang = Triangulation(degrees(x),degrees(y))

    #norm=None
    norm=mpl.colors.LogNorm()
    buoySurf = np.array([1025.6, 1026.85, 1027.4, 1027.7, 1028.1])
    plot_save_fig(kappaIntopenmp, 'kappa', '$\kappa\;(m^2 s^{-1}):'+str(buoySurf[int(layer)])+'$',cs=options.colorscale, nrm=norm)
    #plot_save_fig(kappaIntopenmp_xx, 'kappa_xx', '$\kappa_{xx}\;(m^2 s^{-1})$',cs=options.colorscale, nrm=norm)
    #plot_save_fig(kappaIntopenmp_xy, 'kappa_xy', '$\kappa_{xy}\;(m^2 s^{-1})$',cs=options.colorscale, nrm=norm)
    #plot_save_fig(kappaIntopenmp_yy, 'kappa_yy', '$\kappa_{yy}\;(m^2 s^{-1})$',cs=options.colorscale, nrm=norm)
    #plot_save_fig(sigmaInt_xx, 'sigma_xx', '$\sigma_{xx}\;(m)$')
    #plot_save_fig(sigmaInt_xy, 'sigma_xy', '$\sigma_{xy}\;(m)$')
    #plot_save_fig(sigmaInt_yy, 'sigma_yy', '$\sigma_{yy}\;(m)$')
    #plot_save_fig(sigmaInt_rr, 'sigma_rr', '$\sigma_{rr}\;(m)$')

    #plt.show()


    # {{{
    #refiner = UniformTriRefiner(triang)
    #tri_refi, z_test_refi = refiner.refine_field(kappaIntopenmp,subdiv=3)

    #plt.figure()
    #plt.gca().set_aspect('equal')
    #plt.triplot(triang, lw=0.5, color='white')

    #plt.tricontourf(tri_refi, z_test_refi)
    #plt.colorbar()
    #plt.show()


    #plt.hist(kappaIntserial)
    #plt.show()

    #plt.scatter(x,y,c=kappaIntopenmp)
    #plt.colorbar()
    #plt.show()

    #}}}
    
    


