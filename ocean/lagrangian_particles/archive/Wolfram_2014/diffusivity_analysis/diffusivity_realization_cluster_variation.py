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
from GeneralTriangulation import GeneralTriangulation as Triangulation
import matplotlib 
import matplotlib.pyplot as plt
#from diffusivity_realization_data_lagrangian import normalized_haversine_formula 

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

def haversine_formula(phi1, phi2, lam1, lam2, r=1.0):  #{{{
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

    # scalar diffusivity
    #dr = r * normalized_haversine_formula(clat, plat, clon, plon)
    dr = haversine_formula(clat, plat, clon, plon, r)
    drdr_sum = np.sum(dr*dr, axis=1)
   
    # tensor diffusivity
    # compute distances in m from lat / long  (use COM for simplicity, but there will be some error because of coordinate transform)
    #dx = haversine_formula(clat, clat, clon, plon, r)
    #dy = haversine_formula(clat, plat, clon, clon, r)
    
    ## fix orientation of points
    #bearing = spherical_bearing(clat, plat, clon, plon)
    ## because arctan2 returns results from -pi to pi for bearing, flip values to get right sign
    #dx -= ne.evaluate("2*dx*(abs(bearing) > PI/2.0)")
    #dy -= ne.evaluate("2*dy*(bearing < 0)")

    ## store values at each time
    #dxdx_sum = np.sum(dx*dx, axis=1)
    #dxdy_sum = np.sum(dx*dy, axis=1)
    #dydy_sum = np.sum(dy*dy, axis=1)

    return np.squeeze(clon[:,0,:]), np.squeeze(clat[:,0,:]), drdr_sum #, dxdx_sum, dxdy_sum, dydy_sum  #}}}

def build_particle_file(fname_in, kclusters=100, nlayers=5, layerrange=np.arange(0,5), ti=5, deltat=2., Ndim=100):  #{{{
    
    # load the file database
    f_in = np.load(fname_in, 'r')
    plat = f_in['lat']
    plon = f_in['lon']
    
    Ntime = plon.shape[0]
    tvec = deltat*np.arange(0,Ntime)
    Nparticles = plon.shape[1]
    Npartlayers = Nparticles/nlayers
    Nensembles = plon.shape[2]
    rEarth = 6371220. # from netcdf file
   
    lonmin = np.min(plon[:])
    lonmax = np.max(plon[:])

    latmin = np.min(plat[:])
    latmax = np.max(plon[:])
    
    # cell clusters
    #allparticles = np.load('cellClusters.npy')
    #Nclusters = allparticles.shape[0]/5

    # build the grid
    x,y = np.meshgrid(np.linspace(-15,15,Ndim),np.linspace(23,47,Ndim))
    dist = np.sqrt((13./16.*x)**2 + (y-35)**2.0)
    indomain = np.where(dist < 13.0)
    x = (x[indomain]).flatten()
    y = (y[indomain]).flatten()
    Nclusters = x.ravel().shape[0]
    
    ## use existing particle locations
    #x = plon[0,:Nparticles/nlayers,0]
    #y = plat[0,:Nparticles/nlayers,0]
    #Nclusters = x.ravel().shape[0]


    print 'Nclusters = %d' % (Nclusters)
    
    # allocate memory #{{{
    mux      = np.zeros((Ntime, Nensembles, Nclusters))
    muy      = np.zeros((Ntime, Nensembles, Nclusters))
    dxdx_sum = np.zeros((Ntime, Nensembles, Nclusters))
    dxdy_sum = np.zeros((Ntime, Nensembles, Nclusters))
    dydy_sum = np.zeros((Ntime, Nensembles, Nclusters))
    drdr_sum = np.zeros((Ntime, Nensembles, Nclusters))
    #}}}

    # starting center points of cluster
    for nlayer in np.array([1]): #layerrange: 

        # compute the clusters and compute dispersion in a single layer
        llon = np.radians(plon[:,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:])
        llat = np.radians(plat[:,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:])

        hull = ConvexHull(np.vstack((llon[0,:,0],llat[0,:,0])).T)

        # custom defined clusters
        tree = spatial.cKDTree(np.vstack((llon[0,:,0],llat[0,:,0])).T)
        _, allparticles = tree.query(np.vstack((np.radians(x.ravel()),np.radians(y.ravel()))).T, k=kclusters)
       
        ## cell clusters (somewhat of a hack)
        #x = np.degrees(np.mean(llon[0,allparticles[:Nclusters,:],0],axis=1))
        #y = np.degrees(np.mean(llat[0,allparticles[:Nclusters,:],0],axis=1))

        for ac in np.arange(Nclusters):
            particles = allparticles[ac,:]
    
            # just do single ensemble for now
            lon = llon[:,particles,:]
            lat = llat[:,particles,:]

            # compute center of mass and relative dispersion to this center of mass
            mux[:,:,ac], muy[:,:,ac], \
                    drdr_sum[:,:,ac] = \
                    compute_com_and_dispersion(lat, lon, rEarth)
                    #, dxdx_sum[:,:,ac], dxdy_sum[:,:,ac], dydy_sum[:,:,ac] = \

        print 'finished %d layer' % (nlayer)
        #np.savez('var_sums', drdr_sum=drdr_sum, dxdx_sum=dxdx_sum, dxdy_sum=dxdy_sum, dydy_sum=dydy_sum, mux=mux, muy=muy)

    kappa_rr = np.log10(0.5*np.sum((drdr_sum[ti+1,:,:]-drdr_sum[ti,:,:])*kclusters, axis=0)/(kclusters*Nensembles-1)/(deltat*24.*60.*60.))

    triang = Triangulation(x,y, kappa_rr)
    #triang = Triangulation(x,y, np.log10(np.mean(drdr_sum[ti,:,:],axis=0)))
    #triang = Triangulation(x,y, np.log10(drdr_sum[ti,0,:]))

    plt.figure()
    triang.plot_scalar()
    #plt.clim(2,6)
    plt.colorbar()
    plt.xlim(-20,20)
    
    plt.figure()
    triang.plot_scalar(scalar=kappa_rr, nfilt=10)
    plt.clim(2,6)
    plt.colorbar()
    plt.xlim(-20,20)
    plt.show()
            
    return  #}}}

def original_results(nlayers=5): #{{{
    """ for debugging """ 

    layer = 1
    data = np.load('/Users/pwolfram/Documents/ReportsPresentations/SOMADiffusivityPaper/ProductionRuns/16km_export/particle_ensembles/ensemble000/ensemble.npz')
    drdr = data['drdr']
    x = data['mux']
    y = data['muy']

    Npl = drdr.shape[1]/nlayers

    plt.figure()
    triang = Triangulation(np.degrees(x[0,layer*Npl:(layer+1)*Npl]), np.degrees(y[0,layer*Npl:(layer+1)*Npl]), np.log10(drdr[5,layer*Npl:(layer+1)*Npl]))
    triang.plot_scalar()
    plt.colorbar()

    return #}}}

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension",
                      metavar="FILE")
    parser.add_option("-k", dest="kcluster",
            help="number of points for the cluster", metavar="FLOAT")
    parser.add_option("--nlayers", dest="nlayers",
            help="number of particle layers", metavar="INT")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")
    if not options.kcluster:
        parser.error("Need starting number for clusters.")
    if not options.nlayers:
        options.nlayers = 5
    
    #original_results()
    build_particle_file(options.inputfilename, int(options.kcluster), int(options.nlayers))
