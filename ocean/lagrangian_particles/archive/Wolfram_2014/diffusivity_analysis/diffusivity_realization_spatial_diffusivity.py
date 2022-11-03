#!/usr/bin/env python
"""

    Explore ways of computing a good estimate of diffusivity from single realization.
    
    Phillip Wolfram
    LANL
    02/27/2015

"""
# import libraries / packages
import os
import numpy as np
import numexpr as ne
import netCDF4
import cPickle as pickle
import scipy.spatial as spatial
from scipy.spatial import ConvexHull
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=5, metadata=dict(artist='Phillip J. Wolfram'))

# function definitions

def spherical_bearing(phi1, phi2, lam1, lam2): #{{{
  """
  compute the spherical bearing
      http://www.movable-type.co.uk/scripts/latlong.html
      
      where lambda is longitude and phi is latitude
  
      spherical_bearing returned is from N (-pi, pi)
  
      Phillip Wolfram
      LANL
      08/07/2014
  """
  return ne.evaluate("arctan2(sin(lam2-lam1)*cos(phi2), cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(lam2-lam1))") #}}}

def normalized_haversine_formula(phi1, phi2, lam1, lam2, r=1.0):  #{{{
    """
    compute the distance between two points via Haversine formula:
    http://www.movable-type.co.uk/scripts/latlong.html
    
    where lambda is longitude and phi is latitude

    c returned in non-dimensional units (radians)

    Phillip Wolfram
    LANL
    07/18/2014
    """
    a = ne.evaluate("sin((phi2-phi1)/2.0)**2 + cos(phi1) * cos(phi2) * sin((lam2-lam1)/2.0)**2")
    c = ne.evaluate("r * 2.0 * arctan2(sqrt(a), sqrt(1.0-a))")

    return c #}}}

def proj_lat_long(x, y, z):  #{{{
    """
    compute the latitude and longitude from
    the x,y,z points (follow's Doug's pnt.h)
    """
    plat = ne.evaluate("arcsin(z / sqrt(x ** 2 + y ** 2 + z ** 2))")
    plon = ne.evaluate("arctan2(y, x)")

    return plat, plon  #}}}

def proj_xyz(plat, plong, r):  #{{{
    """
    convert from latitude / longitude to xyz
    Phillip Wolfram
    LANL
    07/17/2014
    """
    x = ne.evaluate("r * cos(plong) * cos(plat)")
    y = ne.evaluate("r * sin(plong) * cos(plat)")
    z = ne.evaluate("r * sin(plat)")
    
    return x,y,z  #}}}


def compute_com_and_dispersion(plat, plon, r): #{{{
    """
    compute particle disperison (2nd moment) for cluster, using lat / lon for basis of 
    calculation to set coordinates for Kxx,Kxy,Kyy where x is zonal and y meridional

    dispersion units in m^2

    Phillip Wolfram
    LANL
    07/18/2014
    """
    PI = np.pi

    # compute center of mass for each time
    clat = np.ones(plat.shape)*np.mean(plat, axis=1)[:,np.newaxis,:]
    clon = np.ones(plon.shape)*np.mean(plon, axis=1)[:,np.newaxis,:]

    # compute distances in m from lat / long  (use COM for simplicity, but there will be some error because of coordinate transform)
    dx = normalized_haversine_formula(clat, clat, clon, plon, r)
    dy = normalized_haversine_formula(clat, plat, clon, clon, r)
    dr = normalized_haversine_formula(clat, plat, clon, plon, r)
    
    # fix orientation of points
    bearing = spherical_bearing(clat, plat, clon, plon)
    # because arctan2 returns results from -pi to pi for bearing, flip values to get right sign
    dx -= ne.evaluate("2*dx*(abs(bearing) > PI/2.0)")
    dy -= ne.evaluate("2*dy*(bearing < 0)")

    # store values at each time
    dxdx = np.sum(dx*dx, axis=1)
    dxdy = np.sum(dx*dy, axis=1)
    dydy = np.sum(dy*dy, axis=1)
    drdr = np.sum(dr*dr, axis=1)

    return np.squeeze(clon[:,0,:]), np.squeeze(clat[:,0,:]), dxdx, dxdy, dydy, drdr #}}}

def build_particle_file(fname_in, nlayers=5, slon=0, slat=0, layerrange=np.arange(0,5)):  #{{{
    
    # load the file database
    f_in = np.load(fname_in, 'r')
    plat = f_in['lat']
    plon = f_in['lon']

    Ntime = plon.shape[0]
    tvec = 2*np.arange(0,Ntime)
    Nparticles = plon.shape[1]
    Nensembles = plon.shape[2]
    rEarth = 6371220. # from netcdf file
    kclusters = np.array([10, 100, 500, 1000, 10000])
    Nclusters = kclusters.shape[0]
    
    # allocate memory
    mux      = np.zeros((Ntime, Nensembles, Nclusters))
    muy      = np.zeros((Ntime, Nensembles, Nclusters))
    dxdx_sum = np.zeros((Ntime, Nensembles, Nclusters))
    dxdy_sum = np.zeros((Ntime, Nensembles, Nclusters))
    dydy_sum = np.zeros((Ntime, Nensembles, Nclusters))
    drdr_sum = np.zeros((Ntime, Nensembles, Nclusters))
    Npart    = np.zeros((Ntime, Nensembles, Nclusters))
    
    # starting center points of cluster
    for nlayer in layerrange: 
        case = '_lon%f_lat%f_layer%d' % (slon, slat, nlayer)
        print case

        Npartlayers = Nparticles/nlayers
        # compute the clusters and compute dispersion in a single layer
        llon = plon[:,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:]
        llat = plat[:,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:]

        tree = spatial.cKDTree(np.vstack((llon[0,:,0],llat[0,:,0])).T)
        hull = ConvexHull(np.vstack((llon[0,:,0],llat[0,:,0])).T)
        
        for ac, k in enumerate(kclusters):
            _, particles = tree.query(np.vstack((slon,slat)).T,k=k) 
            particles = np.squeeze(particles)
            print 'k = %d ' % (k)
        
            # just do single ensemble for now
            lon = llon[:,particles,:]
            lat = llat[:,particles,:]

            # compute center of mass and relative dispersion to this center of mass
            mux[:,:,ac], muy[:,:,ac], \
            dxdx_sum[:,:,ac], dxdy_sum[:,:,ac], \
            dydy_sum[:,:,ac], drdr_sum[:,:,ac] = \
                    compute_com_and_dispersion(lat, lon, rEarth)

            # animate the cluster #{{{
            if True:
                # absolute
                movfig = plt.figure()
                ax = movfig.add_subplot(111, autoscale_on=False, xlim=(-20,20),ylim=(20,50))
                ax.hold(True)
                for simplex in hull.simplices:
                    ax.plot(llon[0,simplex,0],llat[0,simplex,0],'k-',lw=2,zorder=999)
                plt.xlabel('Lon')
                plt.ylabel('Lat')
                line , = ax.plot([],[], '.')
                time_text = ax.text(0.05,0.9, '', transform=ax.transAxes)
                time_template = 'time = %d days'
                def init():
                    line.set_data([],[])
                    time_text.set_text('')
                    return line
                def animate(t):
                    line.set_data(np.ravel(lon[t,:,:]),np.ravel(lat[t,:,:]))
                    time_text.set_text(time_template%(2*t))
                    return line
                im_ani = animation.FuncAnimation(movfig, animate, np.arange(0,Ntime), 
                        interval=50, init_func=init)
                im_ani.save('absolute_clusters%d' % (k) + case + '.mp4', writer=writer)
                
                # relative
                movfig = plt.figure()
                ax = movfig.add_subplot(111, autoscale_on=False, xlim=(-20,20),ylim=(20,50))
                ax.hold(True)
                for simplex in hull.simplices:
                    ax.plot(llon[0,simplex,0],llat[0,simplex,0],'k-',lw=2,zorder=999)
                plt.xlabel('Lon')
                plt.ylabel('Lat')
                line , = ax.plot([],[], '.')
                time_text = ax.text(0.05,0.9, '', transform=ax.transAxes)
                time_template = 'time = %d days'
                def init():
                    line.set_data([],[])
                    time_text.set_text('')
                    return line
                def animate(t):
                    line.set_data(np.ravel(lon[t,:,:]-mux[t,:,ac] + slon),np.ravel(lat[t,:,:]-muy[t,:,ac] + slat))
                    time_text.set_text(time_template%(2*t))
                    return line
                im_ani = animation.FuncAnimation(movfig, animate, np.arange(0,Ntime), 
                        interval=50, init_func=init)
                im_ani.save('relative_clusters%d' % (k) + case + '.mp4', writer=writer)

            #}}}

        # overview figure #{{{ 
        plt.figure(figsize= (32.  ,  16.95))
        plt.subplot(3,3,1)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            _, particles = tree.query(np.vstack((slon,slat)).T,k=k) 
            particles = np.squeeze(particles)
            plt.plot(llon[0,particles,0],llat[0,particles,0],'.', label="k=%1.0e" % (k), zorder=10-ac)
        plt.legend()

        plt.subplot(3,3,2)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            plt.plot(tvec, drdr_sum[:,0,ac]/k, label="k=%1.0e" % (k))
        plt.title('rr ens 0')
        plt.legend(loc='best')

        plt.subplot(3,3,3)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            plt.plot(tvec, np.mean(drdr_sum[:,:,ac],axis=1)/k, label="k=%1.0e" % (k))
        plt.title('rr')
        plt.legend(loc='best')
        
        plt.subplot(3,3,4)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            plt.plot(tvec, dxdx_sum[:,0,ac]/k, label="k=%1.0e" % (k))
        plt.title('xx ens 0')
        plt.legend(loc='best')
        
        plt.subplot(3,3,5)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            plt.plot(tvec, dydy_sum[:,0,ac]/k, label="k=%1.0e" % (k))
        plt.title('yy ens 0')
        plt.legend(loc='best')
        
        plt.subplot(3,3,6)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            plt.plot(tvec, dxdy_sum[:,0,ac]/k, label="k=%1.0e" % (k))
        plt.title('xy ens 0')
        plt.legend(loc='best')
        
        plt.subplot(3,3,7)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            plt.plot(tvec, np.mean(dxdx_sum[:,:,ac],axis=1)/k, label="k=%1.0e" % (k))
        plt.title('xx ')
        plt.legend(loc='best')
        
        plt.subplot(3,3,8)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            plt.plot(tvec, np.mean(dydy_sum[:,:,ac],axis=1)/k, label="k=%1.0e" % (k))
        plt.title('yy ')
        plt.legend(loc='best')
        
        plt.subplot(3,3,9)
        plt.hold(True)
        for ac, k in enumerate(kclusters):
            plt.plot(tvec, np.mean(dxdy_sum[:,:,ac],axis=1)/k, label="k=%1.0e" % (k))
        plt.title('xy ')
        plt.legend(loc='best')
       
        plt.savefig('overview' + case + '.png')
        #plt.show()
        #}}}

        plt.close('all')
    

    return  #}}}

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension",
                      metavar="FILE")
    parser.add_option("--slon", dest="slon",
            help="cluster starting longitude", metavar="FLOAT")
    parser.add_option("--slat", dest="slat",
            help="cluster starting latitude", metavar="FLOAT")
    parser.add_option("--nlayers", dest="nlayers",
            help="number of particle layers", metavar="INT")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")
    if not options.slon:
        parser.error("Need starting longitude.")
    if not options.slat:
        parser.error("Need starting latitude.")
    if not options.nlayers:
        options.nlayers = 5
    
    build_particle_file(options.inputfilename, int(options.nlayers), float(options.slon), float(options.slat))
