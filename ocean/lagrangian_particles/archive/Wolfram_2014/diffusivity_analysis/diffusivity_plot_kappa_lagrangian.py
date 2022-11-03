#!/usr/bin/env python

# import libraries / packages
from lxml import etree
import os
import numpy as np
from scipy.stats import chi2
import netCDF4

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
from matplotlib.tri import LinearTriInterpolator
from GeneralTriangulation import GeneralTriangulation as Triangulation
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import ConvexHull

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

Rearth = 6371229.0
# function definitions

def signed_log10(x): #{{{
  return np.sign(x)*np.log10(np.abs(x)) #}}}

def signed_sqrt(x): #{{{
  return np.sign(x)*np.sqrt(np.abs(x)) #}}}

def safe_log10(x): #{{{
  #return np.log10(np.maximum(x,np.zeros(x.shape))) #}}}
  return np.log10(x) #}}}

def plot_cov_ellipse(x, y, cov, *args, **kwargs): #{{{
  xp, yp = ellipse_points(x, y, cov, N=60)
  plt.plot(xp, yp, *args, **kwargs)
  return #}}}

def ellipse_points(xc, yc, cov, N=30, volume=None, nstd=1): #{{{
    """
    return collection of points for ellipse
    volume=0.118 for 1-sigma.
    Adapted from 
    http://www.nhsilbert.net/source/2014/06/bivariate-normal-ellipse-plotting-in-python/
    See also
    https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf

    """

    if volume is None:
        # set to be 1-sigma away
        rho = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
        r2  = 2*nstd**2.0/(1+rho)
        volume = chi2.cdf(r2,2)
   
    # compute the width, height, and rotation for the ellipse
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    vals, vecs = eigsorted(cov)

    theta = np.arctan2(*vecs[:,0][::-1])
    width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)

    # generate ellipse points
    t = np.linspace(0,2*np.pi,N)
    # generate points
    xp = width/2.0  * np.cos(t)
    yp = height/2.0 * np.sin(t)
    ## rotate
    pr = np.vstack((xp,yp))
    Rot = np.array([[np.cos(theta), - np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    pr = np.dot(Rot,pr)
    # translate
    xp = pr[0,:] + xc
    yp = pr[1,:] + yc

    return xp, yp #}}}

def test_ellipse_points(sigxx,sigyy, sigmaxy): #{{{
    rho = sigmaxy/np.sqrt(sigxx*sigyy)
    r2  = 2/(1+rho)
    vol = chi2.cdf(r2,2)
    x,y = ellipse_points(2.0,4.0,np.array([[sigxx,sigmaxy],[sigmaxy,sigyy]]),N=1000)
    width = max(x)-min(x)
    height = max(y)-min(y)
    plt.plot(x,y, 'k-')
    cov = np.cov(np.vstack((x,y)))
    print cov
    xt, yt = ellipse_points(x.mean(), y.mean(), cov)
    plt.plot(xt,yt,'o')
    # note, 'o''s dont have to be on line, but shape must be "parallel"

    plt.axis('equal')
    plt.show()
    return #}}}

def convert_ParaView_xml_to_matplotlib_colormap(fname): #{{{
    
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

    return cm #}}}

def get_ensemble_paths(rootdir, folder_prefix): #{{{
    fnames = []
    for root, dirs, files in os.walk(rootdir):
      if(root==rootdir):
        dirs.sort()
        for adir in dirs:
          if(adir.startswith(folder_prefix)):
            fnames.append(rootdir + '/' + adir)

    return fnames #}}}

def plot_save_fig(savename, triang, data, titlename, clim=None, \
    usecmap=convert_ParaView_xml_to_matplotlib_colormap('rainbow_desaturated.xml')): #{{{
   
    ## mask out bad values
    tri_dat = data[triang.triangles].min(axis=1)
    mask = np.where(np.isnan(tri_dat), 1, 0)
    #data[mask] = 0
    triang.set_mask(mask)
   
    # plot the data
    plt.figure()
    plt.gca().set_aspect('equal')
    plt.hold(True)
    #plt.triplot(triang, lw=0.5, color='white')
    #plt.tricontourf(triang, kappaIntopenmp)
    #usecmap = plt.cm.spectral
    #usecmap = plt.cm.gist_ncar
    #usecmap = plt.cm.hsv
    #usecmap = plt.cm.jet
    #usecmap = convert_ParaView_xml_to_matplotlib_colormap('rainbow_desaturated.xml')
    plt.tricontour(triang, data, colors='k', lw=1, levels=[0])
    plt.tripcolor(triang, data, shading='flat', cmap=usecmap)
    cb = plt.colorbar()
    if clim is not None:
        plt.clim(clim)
        crange = np.arange(clim[0],clim[1]+1)
        cb.set_ticks(crange)
    plt.xlabel('$Longitude~(^\circ)$')
    plt.ylabel('$Latitude~(^\circ)$')
    plt.title(titlename)
    plt.savefig(savename)
    print 'saved %s' % (savename)

    return #}}}

def plot_mean_layer_depth(layer, cs_min, cs_max, dec): #{{{
    # load mean ssh
    depth_file = netCDF4.Dataset('buoyancySurface.nc','r')
    depth = depth_file.variables['buoyancySurfaceDepth'+str(layer)]
    lonCell = depth_file.variables['lonCell']
    latCell = depth_file.variables['latCell']
    lon = np.degrees(np.mod(lonCell[:]+np.pi,2*np.pi)-np.pi)
    lat = np.degrees(latCell[:])
    triang = Triangulation(lon, lat)

    depth_filt = triang.smooth_laplacian(-depth[0,:], ntimes=10)
    
    #usecmap = convert_ParaView_xml_to_matplotlib_colormap('../../rainbow_desaturated.xml')
    #plt.tripcolor(triang, depth[0,:], shading='flat', cmap=usecmap)
    #plt.clim([options.colorscale_min, options.colorscale_max])
    #plt.colorbar()
    
    c1 = plt.tricontour(triang, depth_filt, -np.linspace(cs_min, cs_max, dec), colors='#7e7e7e', zorder=998)
    c2vals = -np.linspace(cs_min, cs_max ,(dec-1)*2+1)
    c2vals = c2vals[1::2]
    c2 = plt.tricontour(triang, depth_filt, c2vals, colors='#e0e0e0', zorder=997)
    plt.clabel(c1, fmt='%8.0f', colors='#7e7e7e', fontsize=10, zorder=998)
    
    #contour = plt.tricontour(triang, depth[0,:], np.linspace(-0.6,0.6,7), colors='#7e7e7e', zorder=997)
    #plt.clabel(contour, fmt = '%2.1f', colors = '#e0e0e0', fontsize=14,zorder=0)
    
    hull = ConvexHull(np.vstack((lon,lat)).T)
    for simplex in hull.simplices:
      plt.plot(lon[simplex],lat[simplex],'k-',lw=2,zorder=999)

    return  #}}}

def plot_ellipses(xd, yd, k_xxd, k_xyd, k_yyd, scaling, Ndim=15): #{{{
   
    # build the grid
    x,y = plt.meshgrid(np.linspace(-15,15,Ndim),np.linspace(23,47,Ndim))
    dist = np.sqrt((13./16.*x)**2 + (y-35)**2.0)
    indomain = np.where(dist < 13.0)
    x = (x[indomain]).flatten()
    y = (y[indomain]).flatten()

    # get new data on the grid
    triang = Triangulation(xd,yd)
    k_xxI = LinearTriInterpolator(triang, k_xxd)
    k_xyI = LinearTriInterpolator(triang, k_xyd)
    k_yyI = LinearTriInterpolator(triang, k_yyd)
    k_xx = k_xxI(x,y)
    k_xy = k_xyI(x,y)
    k_yy = k_yyI(x,y)

    # plot the ellipses for the data
    for i in np.arange(x.shape[0]):
      cov = scale_cov(np.array([[k_xx[i],k_xy[i]],[k_xy[i],k_yy[i]]]), scaling)
      plot_cov_ellipse(x[i], y[i], cov, 'k-', zorder=999)

    return #}}}

def scale_cov(cov, scaling): #{{{
  return (cov*(1.0/Rearth*180/np.pi))**2.0*scaling #}}}

def plot_ellipse_legend(sf): #{{{
    # parameters affecting ploting of ellipses
    legcenter = np.array([15, 22.5])
    textcenter = legcenter + np.array([-2,3])
    scaling = sf*10
    #plt.text(textcenter[0], textcenter[1], r"$\bf{\kappa}~("+str(sf)+r"~\times 10^5)$")

    cov = scale_cov(np.array([[5e4, 0],[0, 1e4]]), scaling)
    plot_cov_ellipse(legcenter[0], legcenter[1], cov, 'k-', lw=2)
    plt.text(legcenter[0]+2.3, legcenter[1]-0.4, '5')
    plt.text(legcenter[0]-0.2, legcenter[1]+0.7, '1')

    return scaling #}}}

def plot_save_fig_ellipse(savename, xint, yint, k_xx, k_xy, k_yy, layer, sf, cs_min, cs_max, dec): #{{{
    plt.figure()
    plt.gca().set_aspect('equal')
    plt.hold(True)

    plot_mean_layer_depth(layer, cs_min, cs_max, dec)

    scaling = plot_ellipse_legend(sf)

    plot_ellipses(xint, yint, k_xx, k_xy, k_yy, scaling)
    
    plt.xlabel('$Longitude~(^\circ)$')
    plt.ylabel('$Latitude~(^\circ)$')
    buoySurf = np.array([1025.6, 1026.85, 1027.4, 1027.7, 1028.1])
    plt.title(r"$\bf{\kappa}~("+str(sf)+r"\times 10^4~m^2s^{-1}):"+str(buoySurf[int(layer)])+"$")
    plt.savefig(savename)
    print 'saved %s' % (savename)
    
    return  #}}}

def plot_diffusivity(root, prefix, ts, dt, tint, nlayers): #{{{

  # multiple layer data #{{{
  sf = [1., 1., 1., 1., 1.]
  cs_min = [-30., -450., -600., -800., -1400.]
  cs_max = [0., -50., -50., -200., -900.]
  dec = [4, 9, 12, 13, 11]
  #}}}
  
  ensembles = get_ensemble_paths(root, prefix)
  ensembles.reverse()
  
  # loop over all files and aggregate data
  for aensemble in ensembles:
    print 'processing %s ' % (aensemble)

    for num, (ats, adt, atint) in enumerate(zip(ts, dt, tint)):
      #if num != 2:
      #  continue
      loadname ='/kappa-' + str(ats) + '-' + str(adt) + '-' + str(atint) + '.npz'
      print 'loading %s' % (loadname)
      data = np.load(aensemble+loadname)
      print 'done'

      mux = data['mux_ts']
      muy = data['muy_ts']
      dxdx = data['dxdx_ts']
      dxdy = data['dxdy_ts']
      dydy = data['dydy_ts']
      drdr = data['drdr_ts']
      kappa_xx = data['kappa_xx']
      kappa_xy = data['kappa_xy']
      kappa_yy = data['kappa_yy']
      kappa_rr = data['kappa_rr']
      
      for alayer in np.arange(nlayers):
        print 'plotting diffusivity for layer %d ' % (alayer)
        savename =  str(ats) + '-' + str(adt) + '-' + str(atint) + '-layer' + str(alayer) + '.png'
        Npartlayer = mux.shape[0]/nlayers
        
        x = np.degrees(mux[alayer*Npartlayer:(alayer+1)*Npartlayer])
        y = np.degrees(muy[alayer*Npartlayer:(alayer+1)*Npartlayer])
        triang = Triangulation(x,y)
        
        plt.close('all')
        
        k_xx = kappa_xx[alayer*Npartlayer:(alayer+1)*Npartlayer]
        k_xy = kappa_xy[alayer*Npartlayer:(alayer+1)*Npartlayer]
        k_yy = kappa_yy[alayer*Npartlayer:(alayer+1)*Npartlayer]
        k_rr = kappa_rr[alayer*Npartlayer:(alayer+1)*Npartlayer]
        
        
        plot_save_fig_ellipse(aensemble + '/' + 'kappa_tensor'+savename, x, y, k_xx, k_xy, k_yy, alayer, sf[alayer], cs_min[alayer], cs_max[alayer], dec[alayer])
        
        #plot_save_fig(aensemble + '/' + 'kappa_rr'+savename, triang, k_rr, '$\kappa_{rr}\;(m^2 s^{-1})$')
        #plot_save_fig(aensemble + '/' + 'kappa_xx'+savename, triang, k_xx, '$\kappa_{xx}\;(m^2 s^{-1})$')
        #plot_save_fig(aensemble + '/' + 'kappa_xy'+savename, triang, k_xy, '$\kappa_{xy}\;(m^2 s^{-1})$')
        #plot_save_fig(aensemble + '/' + 'kappa_yy'+savename, triang, k_yy, '$\kappa_{yy}\;(m^2 s^{-1})$')
       
        plot_save_fig(aensemble + '/' + 'kappa_rr_log'+savename, triang, signed_log10(k_rr),   \
            '$\\textrm{sgn}\left(\kappa_{rr}\\right)\log_{10} \left(|\kappa_{rr}|\\right) \;(m^2 s^{-1})$', clim=[-7,7],\
            usecmap=plt.get_cmap('gist_ncar'))
            #usecmap=convert_ParaView_xml_to_matplotlib_colormap('redgrayblue.xml'))
        plot_save_fig(aensemble + '/' + 'kappa_xx_log'+savename, triang, signed_log10(k_xx),   \
            '$\\textrm{sgn}\left(\kappa_{xx}\\right)\log_{10} \left(|\kappa_{xx}|\\right) \;(m^2 s^{-1})$', clim=[-7,7],\
            usecmap=plt.get_cmap('gist_ncar'))
            #usecmap=convert_ParaView_xml_to_matplotlib_colormap('redgrayblue.xml'))                                   
        plot_save_fig(aensemble + '/' + 'kappa_xy_log'+savename, triang, signed_log10(k_xy), \
            '$\\textrm{sgn}\left(\kappa_{xy}\\right)\log_{10} \left(|\kappa_{xy}|\\right) \;(m^2 s^{-1})$', clim=[-7,7],\
            usecmap=plt.get_cmap('gist_ncar'))
            #usecmap=convert_ParaView_xml_to_matplotlib_colormap('redgrayblue.xml'))                                  
        plot_save_fig(aensemble + '/' + 'kappa_yy_log'+savename, triang, signed_log10(k_yy), \
            '$\\textrm{sgn}\left(\kappa_{yy}\\right)\log_{10} \left(|\kappa_{yy}|\\right) \;(m^2 s^{-1})$', clim=[-7,7],\
            usecmap=plt.get_cmap('gist_ncar'))
            #usecmap=convert_ParaView_xml_to_matplotlib_colormap('redgrayblue.xml'))
        
        dxdxl = dxdx[alayer*Npartlayer:(alayer+1)*Npartlayer]
        dxdyl = dxdy[alayer*Npartlayer:(alayer+1)*Npartlayer]
        dydyl = dydy[alayer*Npartlayer:(alayer+1)*Npartlayer]
        drdrl = drdr[alayer*Npartlayer:(alayer+1)*Npartlayer]
        
        plot_save_fig(aensemble + '/' + 'std_rr_log'+savename, triang, signed_log10(signed_sqrt(drdrl)),   \
            '$\\textrm{sgn}\left(\sigma^2_{rr}\\right)\log_{10} \left(|\sigma^2_{rr}|\\right) \;(m^2 s^{-1})$', clim=[-7,7],\
            usecmap=plt.get_cmap('gist_ncar'))
            #usecmap=convert_ParaView_xml_to_matplotlib_colormap('redgrayblue.xml'))
        plot_save_fig(aensemble + '/' + 'std_xx_log'+savename, triang, signed_log10(signed_sqrt(dxdxl)),   \
            '$\\textrm{sgn}\left(\sigma^2_{xx}\\right)\log_{10} \left(|\sigma^2_{xx}|\\right) \;(m^2 s^{-1})$', clim=[-7,7],\
            usecmap=plt.get_cmap('gist_ncar'))
            #usecmap=convert_ParaView_xml_to_matplotlib_colormap('redgrayblue.xml'))
        plot_save_fig(aensemble + '/' + 'std_xy_log'+savename, triang, signed_log10(signed_sqrt(dxdyl)),   \
            '$\\textrm{sgn}\left(\sigma^2_{xy}\\right)\log_{10} \left(|\sigma^2_{xy}|\\right) \;(m^2 s^{-1})$', clim=[-7,7],\
            usecmap=plt.get_cmap('gist_ncar'))
            #usecmap=convert_ParaView_xml_to_matplotlib_colormap('redgrayblue.xml'))
        plot_save_fig(aensemble + '/' + 'std_yy_log'+savename, triang, signed_log10(signed_sqrt(dydyl)),   \
            '$\\textrm{sgn}\left(\sigma^2_{yy}\\right)\log_{10} \left(|\sigma^2_{yy}|\\right) \;(m^2 s^{-1})$', clim=[-7,7],\
            usecmap=plt.get_cmap('gist_ncar'))
            #usecmap=convert_ParaView_xml_to_matplotlib_colormap('redgrayblue.xml'))


  return  #}}}

 
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

  # cases #{{{
  ts = [5] #np.array([5, 5, 5, 5,  5, 15])
  dt = [1] #np.array([1, 2, 4, 8, 16, 15])
  tint = [2] #np.array([2, 2, 2, 2, 2, 2])
  nlayers = 5
  #}}}

  plot_diffusivity(options.root, options.prefix, ts, dt, tint, nlayers)
