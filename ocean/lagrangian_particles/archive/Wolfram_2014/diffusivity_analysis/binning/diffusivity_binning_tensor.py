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
def compute_com_and_squaresums(plat, plon, r): #{{{
    """
    compute relative particle disperison (2nd moment) for cluster, using lat / lon for basis of 
    calculation to set coordinates for Kxx,Kxy,Kyy where x is zonal and y meridional

    dispersion units in m^2

    Phillip Wolfram
    LANL
    07/18/2014
    """

    # compute center of mass for each time
    clat = np.mean(plat, axis=0)
    clon = np.mean(plon, axis=0)

    # scalar diffusivity #{{{
    dr = haversine_formula(clat, plat, clon, plon, r)
    drdr_sum = np.sum(dr*dr, axis=0)
    #}}}
    
    # tensor diffusivity #{{{
    dx, dy = signed_distances(clat, plat, clon, plon, r)
    dxdx_sum = np.sum(dx*dx, axis=0)
    dxdy_sum = np.sum(dx*dy, axis=0)
    dydy_sum = np.sum(dy*dy, axis=0)
    #}}}

    return clon, clat, drdr_sum, dxdx_sum, dxdy_sum, dydy_sum #}}}

def build_particle_file(fname_in, radius=16000, ts=25, te=30, deltat=2., nlayers=5, layerrange=np.arange(0,5), Ndim=4, degrees=False):  #{{{
    
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
    
    # use existing particle locations
    x = plon[0,:Nparticles/nlayers,0]
    y = plat[0,:Nparticles/nlayers,0]
    Nclusters = x.ravel().shape[0]
    
    tri = Triangulation(x,y).coarsen()
    x = tri.x
    y = tri.y
    Nclusters = x.shape[0]
    case = 'CaseBinningTensor_r=%f_ts=%d_te=%d_deltat=%f_coarsenN=%d/'%(radius,ts,te,deltat, Nclusters)

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
        ensminpart = np.zeros((Nclusters))
        ensmaxpart = np.zeros((Nclusters))
        ensmeanpart = np.zeros((Nclusters))
        meanK_rr       = np.zeros((Nclusters))
        meanK_xx       = np.zeros((Nclusters))
        meanK_xy       = np.zeros((Nclusters))
        meanK_yy       = np.zeros((Nclusters))
        
        K_rr = np.zeros((Nclusters, Ntc-1,Nec)) 
        K_xx = np.zeros((Nclusters, Ntc-1,Nec)) 
        K_xy = np.zeros((Nclusters, Ntc-1,Nec)) 
        K_yy = np.zeros((Nclusters, Ntc-1,Nec)) 
        Ncparticles   = np.zeros((Nclusters, Ntc-1,Nec)) 
        #}}}

        # reshape into vectors
        llonv = np.reshape(llon,(Npc*Ntc*Nec))
        llatv = np.reshape(llat,(Npc*Ntc*Nec))

        # compute clusters sequentially to conserve on memory
        clusters = ParticleClusters(llonv, llatv, mode='ball')

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        for ac in np.arange(Nclusters):
            # compute clusters sequentially to conserve on memory
            localparticles = clusters.get_cluster_ball(x[ac], y[ac], radius=radiustree)
            dist = haversine_formula(y[ac], llatv[localparticles], x[ac], llonv[localparticles], rEarth) 
            particles = (dist < radius)

            # reshape subset index into prior form
            index = np.zeros((Npc*Ntc*Nec), dtype=bool)
            index[localparticles] = particles
            # get index chape back to match llon,llat
            index = index.reshape((Npc,Ntc,Nec))

            # get number of particles in the cluster (sum over Np and take extreems over Ne)
            ensminpart[ac] = np.min(np.sum(index, axis=0)[:])
            ensmaxpart[ac] = np.max(np.sum(index, axis=0)[:])
            ensmeanpart[ac] = np.mean(np.sum(index, axis=0)[:])

            for at in np.arange(Ntc-1):
                for ae in np.arange(Nec):
                    particles = index[:, at, ae]
                    # need two time steps to compute diffusivity
                    lon1 = llon[particles,at  ,ae]
                    lon2 = llon[particles,at+1,ae]
                    lat1 = llat[particles,at  ,ae]
                    lat2 = llat[particles,at+1,ae]
                    # compute dispersion
                    Ncparticles[ac,at,ae], _, _, K_rr[ac,at,ae], K_xx[ac,at,ae], K_xy[ac,at,ae], K_yy[ac,at,ae] = \
                        compute_diffusivity(lon1, lat1, lon2, lat2, deltat*24.*60.*60.)
            # average diffusivity
            meanK_rr[ac] = np.nanmean(K_rr[ac,:,:].ravel())
            meanK_xx[ac] = np.nanmean(K_xx[ac,:,:].ravel())
            meanK_xy[ac] = np.nanmean(K_xy[ac,:,:].ravel())
            meanK_yy[ac] = np.nanmean(K_yy[ac,:,:].ravel())

        np.savez(case+'layer%d'%(nlayer), x=np.degrees(x),y=np.degrees(y),
            meanK_rr=meanK_rr, meanK_xx=meanK_xx, meanK_xy=meanK_xy, meanK_yy=meanK_yy,
            K_rr=K_rr, K_xx=K_xx, K_xy=K_xy, K_yy=K_yy,
            ensmin=ensminpart, ensmax=ensmaxpart, ensmean=ensmeanpart,
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
    parser.add_option("-r", "--radius", dest="radius",
            help="radius of the cluster in km", metavar="FLOAT")
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
    if not options.radius:
        parser.error("Need starting number for radius.")
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
    build_particle_file(options.inputfilename, float(options.radius), 
            int(options.ts)-1, int(options.te)-1, float(options.dt), int(options.nlayers), degrees=options.degrees,layerrange=options.layer)
