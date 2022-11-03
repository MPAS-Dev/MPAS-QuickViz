#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
from matplotlib.tri import Triangulation
#from matplotlib.tri import UniformTriRefiner
from numpy import isnan
from convert_ParaView_xml_to_matplotlib_colormap import convert_ParaView_xml_to_matplotlib_colormap      
from numpy import degrees, radians
import netCDF4
from subprocess import call
from scipy.spatial import ConvexHull

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':22})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
R = 6371229.0

def plot_cov_ellipse(cov, pos, volume=.5, ax=None, fc='none', ec=[0,0,0], a=1, lw=1, zorder=999):
    """
    Plots an ellipse enclosing *volume* based on the specified covariance
    matrix (*cov*) and location (*pos*). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    From: http://www.nhsilbert.net/source/2014/06/bivariate-normal-ellipse-plotting-in-python/

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        volume : The volume inside the ellipse; defaults to 0.5
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
    """

    import numpy as np
    from scipy.stats import chi2
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    kwrg = {'facecolor':fc, 'edgecolor':ec, 'alpha':a, 'linewidth':lw}

    # Width and height are "full" widths, not radius
    width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwrg)
    ellip.set_zorder(zorder)

    ax.add_artist(ellip)

if __name__ == "__main__":
    from optparse import OptionParser

    # get command line parameters
    parser = OptionParser()
    parser.add_option("-s", "--sf", dest="sf",
            help="scale factor")
    parser.add_option("-u", "--colorscalemax", dest="colorscale_max",
            help="clim max")
    parser.add_option("-m", "--colorscalemin", dest="colorscale_min",
            help="clim min")
    parser.add_option("-l", "--layer", dest="layer",
            help="layer #")
    parser.add_option("-d", "--declinations", dest="declinations",
            help="number of countours to plot")
            

    options, args  = parser.parse_args()

    if not options.sf:
        parser.error("scale factor is required")
    if options.colorscale_max:
      options.colorscale_max = float(options.colorscale_max)
    if options.colorscale_min:
      options.colorscale_min = float(options.colorscale_min)
    if options.declinations:
      options.declinations = float(options.declinations)
    if not options.layer:
      parser.error("layer required")
    layer = options.layer

    sf = float(options.sf)

    # load mean ssh
    depth_file = netCDF4.Dataset('../../buoyancySurface.nc','r')
    depth = depth_file.variables['buoyancySurfaceDepth'+layer]
    lonCell = depth_file.variables['lonCell']
    latCell = depth_file.variables['latCell']
    lon = degrees(np.mod(lonCell[:]+np.pi,2*np.pi)-np.pi)
    lat = degrees(latCell[:])
    triang = Triangulation(lon, lat)
    plt.figure()
    plt.gca().set_aspect('equal')
    plt.hold(True)
    usecmap = convert_ParaView_xml_to_matplotlib_colormap('../../rainbow_desaturated.xml')
    #plt.tripcolor(triang, depth[0,:], shading='flat', cmap=usecmap)
    #plt.clim([options.colorscale_min, options.colorscale_max])
    #plt.colorbar()
    c1 = plt.tricontour(triang, -depth[0,:], -np.linspace(options.colorscale_min, options.colorscale_max,options.declinations), colors='#7e7e7e', zorder=998)
    c2vals = -np.linspace(options.colorscale_min, options.colorscale_max,(options.declinations-1)*2+1)
    c2vals = c2vals[1::2]
    c2 = plt.tricontour(triang, -depth[0,:], c2vals, colors='#e0e0e0', zorder=997)
    plt.clabel(c1, fmt='%8.0f', colors='#7e7e7e', fontsize=10, zorder=998)
    #contour = plt.tricontour(triang, depth[0,:], np.linspace(-0.6,0.6,7), colors='#7e7e7e', zorder=997)
    #plt.clabel(contour, fmt = '%2.1f', colors = '#e0e0e0', fontsize=14,zorder=0)
    hull = ConvexHull(np.vstack((lon,lat)).T)
    for simplex in hull.simplices:
      plt.plot(lon[simplex],lat[simplex],'k-',lw=2,zorder=999)
    plt.xlabel('$Longitude~(^\circ)$')
    plt.ylabel('$Latitude~(^\circ)$')
    #plt.title(r"$\bf{\kappa}~("+str(sf)+r"\times 10^5~m^2s^{-1}) ~\textrm{and}~ \left< d \right>~(m)$")
    buoySurf = np.array([1025.6, 1026.85, 1027.4, 1027.7, 1028.1])
    plt.title(r"$\bf{\kappa}~("+str(sf)+r"\times 10^5~m^2s^{-1}):"+str(buoySurf[int(layer)])+"$")
    plt.hold(True)
    
    Ndim = 15
    x,y = plt.meshgrid(np.linspace(-15,15,Ndim),np.linspace(23,47,Ndim))
    dist = np.sqrt((13./16.*x)**2 + (y-35)**2.0)
    indomain = np.where(dist < 13.0)
    x = (x[indomain]).flatten()
    y = (y[indomain]).flatten()
    #plt.plot(x,y,'k.')
    # write to a file
    
    with file('ellipse_points.txt', 'w') as outfile:
      outfile.write('%d\n' % len(x))
      for xo,yo in zip(radians(x[:]),radians(y[:])):
        outfile.write('%e %e\n' % (xo,yo))
    # get interpolated points
    print call('../../process_diffusivity 1 buoyancySurfaceCluster'+str(layer)+'.txt ellipse_points.txt ellipse_result.txt',shell=True)

    # parameters affecting ploting of ellipses
    legcenter = np.array([15,22.5])
    textcenter = legcenter + np.array([-2,3])
    scaling = sf*1e5
    #plt.text(textcenter[0], textcenter[1], r"$\bf{\kappa}~("+str(sf)+r"~\times 10^5)$")

    cov = np.array([[4e5, 0],[0, 1e5]])*(1.0/R*180.0/np.pi)**2.0*scaling/sf
    plot_cov_ellipse(cov, legcenter, volume=0.5)
    plt.text(legcenter[0]+2.5, legcenter[1]-0.5, '4')
    plt.text(legcenter[0]-0.5, legcenter[1]+1.5, '1')
    
    kappaIntopenmp_tot = np.loadtxt('ellipse_result.txt' ,skiprows=1)
    #kappaIntopenmp = kappaIntopenmp_tot[:,0]
    kappaIntopenmp_xx = kappaIntopenmp_tot[:,1]
    kappaIntopenmp_xy = kappaIntopenmp_tot[:,2]
    kappaIntopenmp_yy = kappaIntopenmp_tot[:,3]
    #sigmaInt_xx = kappaIntopenmp_tot[:,4]
    #sigmaInt_xy = kappaIntopenmp_tot[:,5]
    #sigmaInt_yy = kappaIntopenmp_tot[:,6]

    Npoints = kappaIntopenmp_xx.shape[0]
    for i in np.arange(Npoints):
      kxx = kappaIntopenmp_xx[i]
      kxy = kappaIntopenmp_xy[i]
      kyy = kappaIntopenmp_yy[i]
      cov = np.array([[kxx,kxy],[kxy,kyy]])*(1.0/R*180/np.pi)**2.0*scaling
      #if(cov.min() >= 0): 
        #plt.plot(x[i],y[i],'.')
        #print cov
      plot_cov_ellipse(cov, np.array([x[i],y[i]]), volume=.5)

    plt.savefig('kappa_tensor'+str(layer)+'.png')
    plt.show()


    
    


