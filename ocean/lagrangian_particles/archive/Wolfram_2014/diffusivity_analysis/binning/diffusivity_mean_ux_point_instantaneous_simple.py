#!/usr/bin/env python
"""

    Compute diffusivity via binning procedure
    
    Phillip Wolfram
    LANL
    03/03/2015

"""
# import libraries / packages
import os
import itertools
import numpy as np
import numexpr as ne
from scipy.spatial import ConvexHull, cKDTree as KDTree
from GeneralTriangulation import GeneralTriangulation as Triangulation
import matplotlib.pyplot as plt
from file_handling import check_folder
from latlon_coordinate_transforms import signed_distances_numexpr as signed_distances, haversine_formula_numexpr as haversine_formula
from cluster_formation import ParticleClusters

# function definitions

def test_get_unique_pairs():
    particles = np.array([0,1,2,3])
    print 'particle list = ', particles
    left, right = get_unique_pairs(particles)
    print np.vstack((left,right)).T
    print 'Should return 0,1  0,2  0,3  1,2  1,3  2,3'
    return 

def get_unique_pairs(particles):
    npairs = 0
    for a in itertools.combinations(particles,2):
        npairs+=1
    left = np.zeros((npairs),dtype='int')
    right = np.zeros((npairs),dtype='int')
    for i, a in enumerate(itertools.combinations(particles,2)):
        left[i] = a[0]
        right[i] = a[1]
    return left, right

def build_particle_file(fname_in, pointlon, pointlat, 
        ts=25, te=30, deltat=2., nlayers=5, layerrange=np.arange(0,5), degrees=False):  #{{{
    
    # load the file database
    f_in = np.load(fname_in, 'r')
    if degrees: 
        plat = np.radians(f_in['lat'])
        plon = np.radians(f_in['lon'])
    else:
        plat = f_in['latrad']
        plon = f_in['lonrad']
    plonVel = f_in['lonVel']
    platVel = f_in['lonVel']
    
    Ntime = plon.shape[0]
    tvec = deltat*np.arange(0,Ntime)
    Nparticles = plon.shape[1]
    Npartlayers = Nparticles/nlayers
    Nensembles = plon.shape[2]
    rEarth = 6371220. # from netcdf file
    
    case = 'CaseClusterRelDiffVel_lon=%f_lat=%f__ts=%d_te=%d_deltat=%f/'%(pointlon, pointlat, ts,te,deltat)
    print 'Case is ' + case
    check_folder(case)

    # convert points for radial clusters 
    x = np.radians(pointlon)
    y = np.radians(pointlat)

    # starting center points of cluster
    for nlayer in layerrange: 
        print 'on layer = ',nlayer

        # compute the clusters and compute dispersion in a single layer
        # need to rearrange data so that we have total realizations and cluster as index (with averaging over total realizations)
        # now it is of size (Np, Nt, Ne)
        llon = np.swapaxes(plon[ts:te+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llat = np.swapaxes(plat[ts:te+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llonVel = np.swapaxes(plonVel[ts:te+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llatVel = np.swapaxes(platVel[ts:te+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        
        # compute clusters sequentially to conserve on memory
        clusters = ParticleClusters(mode='cluster')

        # find cluster closest to point (only use first 5 points for simplicity)
        ac = np.argmin(np.nanmean(np.sqrt((plon[0,clusters.clusters[:Npartlayers,:5].astype('int'),0]-x)**2.0 + (plat[0,clusters.clusters[:Npartlayers,:5].astype('int'),0]-y)**2.0),axis=1)) 

        particles = clusters.get_cluster_index(ac)
        left, right = get_unique_pairs(particles)

        # get number of particles in the cluster (sum over Np and take extreems over Ne)
        npart = left.shape[0]

        # get pair distances
        sr = haversine_formula(llat[left,:,:],llat[right,:,:],llon[left,:,:],llon[right,:,:], rEarth)
        sx, sy = signed_distances(llat[left,:,:],llat[right,:,:],llon[left,:,:],llon[right,:,:], rEarth)

        u = llonVel[right,:,:] - llonVel[left,:,:]
        v = llatVel[right,:,:] - llatVel[left,:,:]
        #s = np.sqrt(u*u + v*v)
        s = ne.evaluate('sqrt(u*u + v*v)')
        
        K_rr = np.mean(np.mean(s*sr,axis=0),axis=1)
        K_xx = np.mean(np.mean(u*sx,axis=0),axis=1)
        K_xy = np.mean(np.mean(u*sy,axis=0),axis=1)
        K_yx = np.mean(np.mean(v*sx,axis=0),axis=1)
        K_yy = np.mean(np.mean(v*sy,axis=0),axis=1)

        np.savez(case+'layer%d'%(nlayer), npart=npart, K_rr=K_rr, K_xx=K_xx, K_xy=K_xy, K_yx=K_yx, K_yy=K_yy)

        def plot_K(K,name):
            plt.figure()
            plt.ylabel(name + ' (m^2/s)')
            plt.xlabel('time (days)')
            plt.plot(deltat*np.arange(K.shape[0]),K)
            savename = case + name +'layer%d.png'%(nlayer)
            print 'saving ' + savename
            plt.savefig(savename)
        plot_K(K_rr,'K_rr')
        plot_K(K_xx,'K_xx')
        plot_K(K_xy,'K_xy')
        plot_K(K_yx,'K_yx')
        plot_K(K_yy,'K_yy')

        print 'finished %d layer' % (nlayer)

            
    return  #}}}



if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser() #{{{
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension",
                      metavar="FILE")
    parser.add_option("-y", "--latitude", dest="lat",
            help="latitude in degress", metavar="FLOAT")
    parser.add_option("-x", "--longitude", dest="lon",
            help="longitude in degress", metavar="FLOAT")
    parser.add_option("--ts", "--startint", dest="ts",
            help="starting time index", metavar="INT")
    parser.add_option("--te", "--endtint", dest="te",
            help="ending time index", metavar="INT")
    parser.add_option("--dt", "--tsampfreq", dest="dt",
            help="time step (sampling interval) in days", metavar="FLOAT")
    parser.add_option("--nlayers", dest="nlayers",
            help="number of particle layers", metavar="INT")
    parser.add_option("-d", "--degrees", dest="degrees", help="Data in degrees? T or F", metavar="BOOL")
    
    options, args = parser.parse_args()

    if not options.inputfilename: 
        parser.error("Input filename is a required input.")
    if not options.lat:
        parser.error("Need latitude of point.")
    if not options.lon:
        parser.error("Need longitude of point.")
    if not options.ts:
        parser.error("Need starting time index (starting at 1).")
    if not options.te:
        parser.error("Need ending time index (starting at 1).")
    if not options.dt:
        parser.error("Need time sampling interval (days).")
    if not options.degrees:
        print 'Warning: assuming data is in radians!'
        options.degrees = False
    else:
        if options.degrees == 'T':
            options.degrees = True
        elif options.degrees == 'F':
            options.degrees = False
        else:
            parser.error('Specify degrees as "T" or "F"')
    if not options.nlayers:
        options.nlayers = 5
    #}}}
   
    print 'starting analysis'
    #original_results()
    build_particle_file(options.inputfilename, float(options.lon), float(options.lat), 
            int(options.ts)-1, int(options.te)-1, float(options.dt), int(options.nlayers), degrees=options.degrees)
