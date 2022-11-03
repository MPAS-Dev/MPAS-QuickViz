#!/usr/bin/env python
"""

    Compute diffusivity via binning procedure
    
    Phillip Wolfram
    LANL
    03/03/2015

"""
# plot on mustang
import matplotlib as mpl
mpl.use('Agg')
# import libraries / packages
import os
import numpy as np
import numexpr as ne
from scipy.spatial import ConvexHull, cKDTree as KDTree
from GeneralTriangulation import GeneralTriangulation as Triangulation
import matplotlib.pyplot as plt
from file_handling import check_folder
from latlon_coordinate_transforms import signed_distances_numpy as signed_distances, haversine_formula_numpy as haversine_formula
from cluster_formation import ParticleClusters
from compute_diffusivity import compute_diffusivity

def build_particle_file(fname_in, pointlon, pointlat, radius=100000, Nlen=5, ts=25, deltat=2., nlayers=5, layerrange=np.arange(0,5), Ndim=4, degrees=False):  #{{{
    
    # load the file database
    f_in = np.load(fname_in, 'r')
    if degrees: 
        plat = np.radians(f_in['lat'])
        plon = np.radians(f_in['lon'])
    else:
        plat = f_in['latrad']
        plon = f_in['lonrad']
    
    Ntime = plon.shape[0]
    tvec = deltat*np.arange(0,Ntime)
    Nparticles = plon.shape[1]
    Npartlayers = Nparticles/nlayers
    Nensembles = plon.shape[2]
    rEarth = 6371220. # from netcdf file
   
    latmin = np.min(plat[:])
    latmax = np.max(plon[:])

    # compute tree radius
    radiustree = radius / (rEarth * np.sin(np.maximum(np.abs(latmin),np.abs(latmax))))
    
    case = 'CaseBinningTensorPoint_lon=%f_lat=%f_r=%f_ts=%d_deltat=%f_Nlen=%d/'%(pointlon, pointlat, radius, ts, deltat, Nlen)

    print 'Case is ' + case
    check_folder(case)

    x = np.radians(pointlon)
    y = np.radians(pointlat)

    # starting center points of cluster
    for nlayer in layerrange: 

        # compute the clusters and compute dispersion in a single layer
        # need to rearrange data so that we have total realizations and cluster as index (with averaging over total realizations)
        # now it is of size (Np, Nt, Ne)
        llon = np.swapaxes(plon[ts:ts+Nlen+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llat = np.swapaxes(plat[ts:ts+Nlen+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)

        Npc=llon.shape[0]
        Ntc=llon.shape[1]
        Nec=llon.shape[2]
        
        # allocate memory #{{{
        meanK_rr       = np.zeros((Nlen))
        meanK_xx       = np.zeros((Nlen))
        meanK_xy       = np.zeros((Nlen))
        meanK_yy       = np.zeros((Nlen))
        
        K_rr = np.zeros((Nlen, Nec)) 
        K_xx = np.zeros((Nlen, Nec)) 
        K_xy = np.zeros((Nlen, Nec)) 
        K_yy = np.zeros((Nlen, Nec)) 
        Ncparticles   = np.zeros((Nlen, Nec)) 
        #}}}

        # compute clusters sequentially to conserve on memory
        clusters = ParticleClusters(llon[:,0,0], llat[:,0,0], mode='ball')

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        # compute clusters sequentially to conserve on memory
        localparticles = clusters.get_cluster_ball(x, y, radius=radiustree)
        dist = haversine_formula(y, llat[localparticles,0,0], x, llon[localparticles,0,0], rEarth) 
        particles = (dist < radius)

        for dtau in np.arange(Nlen):
            for ae in np.arange(Nec):
                # need two time steps to compute diffusivity
                lon1 = llon[particles,0 ,ae]
                lat1 = llat[particles,0 ,ae]
                lon2 = llon[particles,dtau,ae]
                lat2 = llat[particles,dtau,ae]
                # compute dispersion
                Ncparticles[dtau,ae], _, _, K_rr[dtau,ae], K_xx[dtau,ae], K_xy[dtau,ae], K_yy[dtau,ae] = \
                    compute_diffusivity(lon1, lat1, lon2, lat2, ((dtau+1)*deltat)*24.*60.*60.)
        # average diffusivity
        meanK_rr = np.nanmean(K_rr,axis=1)
        meanK_xx = np.nanmean(K_xx,axis=1)
        meanK_xy = np.nanmean(K_xy,axis=1)
        meanK_yy = np.nanmean(K_yy,axis=1)

        np.savez(case+'layer%d'%(nlayer), 
            meanK_rr=meanK_rr, meanK_xx=meanK_xx, meanK_xy=meanK_xy, meanK_yy=meanK_yy,
            K_rr=K_rr, K_xx=K_xx, K_xy=K_xy, K_yy=K_yy,
            Ncparticles=Ncparticles)

        def plot_K(K,name):
            plt.close('all')
            plt.figure()
            plt.plot(deltat*np.arange(Nlen),K)
            plt.xlabel('time')
            plt.ylabel(name)
            plt.savefig(case +'layer%d_'%(nlayer) + name+'.png')

        plot_K(meanK_rr,'meanK_rr')
        plot_K(meanK_xx,'meanK_xx')
        plot_K(meanK_xy,'meanK_xy')
        plot_K(meanK_yy,'meanK_yy')
        
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
    parser.add_option("-r", "--radius", dest="radius",
            help="radius of the cluster in km", metavar="FLOAT")
    parser.add_option("--ts", "--startint", dest="ts",
            help="starting time index", metavar="INT")
    parser.add_option("--dt", "--tsampfreq", dest="dt",
            help="time step (sampling interval) in days", metavar="FLOAT")
    parser.add_option("--nlayers", dest="nlayers",
            help="number of particle layers", metavar="INT")
    parser.add_option("-d", "--degrees", dest="degrees", help="Data in degrees? T or F", metavar="BOOL")
    parser.add_option("-l", "--layer", dest="layer", help="Layer number", metavar="INT")
    parser.add_option("-n", "--numlen", dest="nlen",
            help="number of dtau derivatives", metavar="FLOAT")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")
    if not options.radius:
        parser.error("Need starting number for radius.")
    if not options.nlen:
        parser.error("Need number of points for dtau derivatives.")
    if not options.lat:
        parser.error("Need latitude of point.")
    if not options.lon:
        parser.error("Need longitude of point.")
    if not options.ts:
        parser.error("Need starting time index (starting at 1).")
    if not options.dt:
        parser.error("Need time sampling interval (days).")
    if not options.nlayers:
        options.nlayers = 5
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
    if not options.layer:
      options.layer = np.arange(5)
    else:
      options.layer = np.array([int(options.layer)])

    #}}}
   
    print 'starting analysis'
    #original_results()
    build_particle_file(options.inputfilename, float(options.lon), float(options.lat), float(options.radius), float(options.nlen), 
            int(options.ts)-1, float(options.dt), int(options.nlayers), degrees=options.degrees,layerrange=options.layer)
