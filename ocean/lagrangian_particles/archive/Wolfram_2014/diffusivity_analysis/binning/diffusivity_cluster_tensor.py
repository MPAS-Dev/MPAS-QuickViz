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

# function definitions
def build_particle_file(fname_in, ts=25, te=30, deltat=2., nlayers=5, layerrange=np.arange(0,5), degrees=False):  #{{{

    
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
    
    case = 'CaseClusterTensor_ts=%d_te=%d_deltat=%f/'%(ts,te,deltat)
        
    # load cached clusters
    clusters = ParticleClusters(mode='cluster')
    Nclusters = clusters.clusters.shape[0]/nlayers

    print 'Nclusters = %d' % (Nclusters)
    print 'Case is ' + case
    check_folder(case)

    # starting center points of cluster
    for nlayer in layerrange: 

        # compute the clusters and compute dispersion in a single layer
        # need to rearrange data so that we have total realizations and cluster as index (with averaging over total realizations)
        # now it is of size (Np, Nt, Ne)
        llon = np.swapaxes(plon[ts:te+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llat = np.swapaxes(plat[ts:te+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)

        Npc=llon.shape[0]
        Ntc=llon.shape[1]
        Nec=llon.shape[2]
        
        # allocate memory #{{{
        meanK_rr       = np.zeros((Nclusters))
        meanK_xx       = np.zeros((Nclusters))
        meanK_xy       = np.zeros((Nclusters))
        meanK_yy       = np.zeros((Nclusters))
        meanx = np.zeros((Nclusters))
        meany = np.zeros((Nclusters))
        
        x = np.zeros((Nclusters, Ntc-1,Nec)) 
        y = np.zeros((Nclusters, Ntc-1,Nec)) 
        K_rr = np.zeros((Nclusters, Ntc-1,Nec)) 
        K_xx = np.zeros((Nclusters, Ntc-1,Nec)) 
        K_xy = np.zeros((Nclusters, Ntc-1,Nec)) 
        K_yy = np.zeros((Nclusters, Ntc-1,Nec)) 
        Ncparticles   = np.zeros((Nclusters, Ntc-1,Nec)) 
        #}}}

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        for ac in np.arange(Nclusters):
            particles = clusters.get_cluster_id(ac)
            # use starting positions

            for at in np.arange(Ntc-1):
                for ae in np.arange(Nec):
                    # need two time steps to compute diffusivity
                    lon1 = llon[particles,at  ,ae]
                    lon2 = llon[particles,at+1,ae]
                    lat1 = llat[particles,at  ,ae]
                    lat2 = llat[particles,at+1,ae]
                    # compute dispersion
                    Ncparticles[ac,at,ae], x[ac,at,ae], y[ac,at,ae], K_rr[ac,at,ae], K_xx[ac,at,ae], K_xy[ac,at,ae], K_yy[ac,at,ae] = \
                        compute_diffusivity(lon1, lat1, lon2, lat2, deltat*24.*60.*60.)
            # average diffusivity
            meanK_rr[ac] = np.nanmean(K_rr[ac,:,:].ravel())
            meanK_xx[ac] = np.nanmean(K_xx[ac,:,:].ravel())
            meanK_xy[ac] = np.nanmean(K_xy[ac,:,:].ravel())
            meanK_yy[ac] = np.nanmean(K_yy[ac,:,:].ravel())
            ## use mean cluster location at diffusivity computation time
            meanx[ac] = np.nanmean(x[ac,:,:].ravel())
            meany[ac] = np.nanmean(y[ac,:,:].ravel())

        np.savez(case+'layer%d'%(nlayer),
                meanx=np.degrees(meanx),meany=np.degrees(meany),
                x=np.degrees(x), y=np.degrees(y), meanK_rr=meanK_rr,
                meanK_xx=meanK_xx, meanK_xy=meanK_xy, meanK_yy=meanK_yy,
                K_rr=K_rr, K_xx=K_xx, K_xy=K_xy, K_yy=K_yy,
                Ncparticles=Ncparticles)
        
        print 'finished %d layer' % (nlayer)

    return  #}}}

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser() #{{{
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension",
                      metavar="FILE")
    parser.add_option("--ts", "--startint", dest="ts",
            help="starting time index", metavar="INT")
    parser.add_option("--te", "--endtint", dest="te",
            help="ending time index", metavar="INT")
    parser.add_option("--dt", "--tsampfreq", dest="dt",
            help="time step (sampling interval) in days", metavar="FLOAT")
    parser.add_option("--nlayers", dest="nlayers",
            help="number of particle layers", metavar="INT")
    parser.add_option("-d", "--degrees", dest="degrees", help="Data in degrees? T or F", metavar="BOOL")
    parser.add_option("-l", "--layer", dest="layer", help="Layer number", metavar="INT")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")
    if not options.ts:
        parser.error("Need starting time index (starting at 1).")
    if not options.te:
        parser.error("Need ending time index (starting at 1).")
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
    build_particle_file(options.inputfilename, int(options.ts)-1,
            int(options.te)-1, float(options.dt), int(options.nlayers),
            degrees=options.degrees,layerrange=options.layer)
