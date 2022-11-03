#!/usr/bin/env python
"""

    Compute diffusivity via binning procedure
    
    Phillip Wolfram
    LANL
    03/03/2015

"""
# import libraries / packages
import os
import numpy as np
import numexpr as ne
from scipy.spatial import ConvexHull, cKDTree as KDTree
from GeneralTriangulation import GeneralTriangulation as Triangulation, test_coarsen
import matplotlib.pyplot as plt
from file_handling import check_folder
from latlon_coordinate_transforms import spherical_bearing, haversine_formula_numexpr as haversine_formula
from cluster_formation import ParticleClusters

# function definitions
def compute_rel_dispersion(plat, plon, r): #{{{
    """
    compute relative particle disperison (2nd moment) for cluster, using lat / lon for basis of 
    calculation to set coordinates for Kxx,Kxy,Kyy where x is zonal and y meridional

    dispersion units in m^2

    Phillip Wolfram
    LANL
    07/18/2014
    """

    # compute center of mass for each time
    clat = np.mean(plat, axis=0)[np.newaxis,:,:]
    clon = np.mean(plon, axis=0)[np.newaxis,:,:]

    # scalar diffusivity
    dr = haversine_formula(clat, plat, clon, plon, r)
    drdr_sum = np.nansum(dr*dr, axis=0)

    return drdr_sum #}}}

def compute_abs_dispersion(plat, plon, slat, slon, r): #{{{
    """
    compute absolute particle disperison 

    dispersion units in m^2

    Phillip Wolfram
    LANL
    03/05/2015
    """
    # scalar diffusivity
    dr = haversine_formula(slat, plat, slon, plon, r)
    
    drdr_sum = np.nansum(dr*dr, axis=0)

    return drdr_sum #}}}

def build_particle_file(fname_in, radius=16000, ts=25, te=30, deltat=2., nlayers=5, layerrange=np.arange(0,5), Ndim=20):  #{{{

    case = 'CaseAbsRel_r=%f_ts=%d_te=%d_deltat=%f_Ndim=%d/'%(radius,ts,te,deltat,Ndim)
    print 'Case is ' + case
    check_folder(case)
    
    # load the file database
    f_in = np.load(fname_in, 'r')
    plat = f_in['latrad']
    plon = f_in['lonrad']
    
    Ntime = plon.shape[0]
    tvec = deltat*np.arange(0,Ntime)
    Nparticles = plon.shape[1]
    Npartlayers = Nparticles/nlayers
    Nensembles = plon.shape[2]
    rEarth = 6371220. # from netcdf file
    
    # compute tree radius
    latmin = np.min(plat[:])
    latmax = np.max(plon[:])
    radiustree = radius / (rEarth * np.sin(np.maximum(np.abs(latmin),np.abs(latmax))))
    
    ## build the grid 
    #x,y = np.meshgrid(np.linspace(-15,15,Ndim),np.linspace(23,47,Ndim))
    #dist = np.sqrt((13./16.*x)**2 + (y-35)**2.0)
    #indomain = np.where(dist < 13.0)
    #x = np.radians((x[indomain]).flatten())
    #y = np.radians((y[indomain]).flatten())
    #Nclusters = x.ravel().shape[0]
    nlayer = 0
    x = np.swapaxes(plon[:,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
    y = np.swapaxes(plat[:,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)

    # coarsen grid twice (reduction of ~1/(6*6) points)
    tri = Triangulation(x[:,0,0],y[:,0,0]).coarsen().coarsen()
    x = tri.x
    y = tri.y
    Nclusters = x.shape[0]
    print 'Nclusters = %d' % (Nclusters)

    # allocate memory #{{{
    enspart = np.zeros((Nclusters))
    kappa_rr   = np.zeros((Nclusters))
    K_rr       = np.zeros((Nclusters))
    ratio      = np.zeros((Nclusters))
    #}}}
    
    # starting center points of cluster
    for nlayer in layerrange: 

        # compute the clusters and compute dispersion in a single layer
        # need to rearrange data so that we have total realizations and cluster as index (with averaging over total realizations)
        # now it is of size (Np, Nt, Ne)
        llon = np.swapaxes(plon[:,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llat = np.swapaxes(plat[:,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)

        #hull = ConvexHull(np.vstack((llon[:,0,0],llat[:,0,0])).T)

        # previous cluster assignments #{{{
        ## custom defined clusters (number of particles) 
        #tree = KDTree(np.vstack((llon[0,:,0],llat[0,:,0])).T)
        #_, allparticles = tree.query(np.vstack((x.ravel(),y.ravel())).T, k=kcluster)
        
        ## cell clusters (somewhat of a hack)
        #x = np.degrees(np.mean(llon[0,allparticles[:Nclusters,:],0],axis=1))
        #y = np.degrees(np.mean(llat[0,allparticles[:Nclusters,:],0],axis=1))
        #}}}
        
        # reshape into vectors
        llonv = llon[:,0,0]
        llatv = llat[:,0,0]

        ## compute clusters all at once
        clusters = ParticleClusters(llonv, llatv, x.ravel(), y.ravel(), radius=radiustree, mode='tree')

        # compute clusters sequentially to conserve on memory
        #clusters = ParticleClusters(llonv, llatv, mode='ball')

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        for ac in np.arange(Nclusters):
            # compute clusters all at once
            localparticles = clusters.clusters[ac]
            ## compute clusters sequentially to conserve on memory
            #localparticles = clusters.get_cluster(x[ac], y[ac], radius=radiustree)
            dist = haversine_formula(y[ac], llatv[localparticles], x[ac], llonv[localparticles], rEarth) 
            particles = (dist < radius)

            lon = llon[localparticles,:,:]
            lat = llat[localparticles,:,:]
            slon = llonv[localparticles]
            slat = llatv[localparticles]
            # compute dispersion
            drdr_sum_rel = compute_rel_dispersion(lat[particles,:,:], lon[particles,:,:], rEarth)
            drdr_sum_abs = compute_abs_dispersion(lat[particles,:,:], lon[particles,:,:], \
                    slat[particles,np.newaxis,np.newaxis], slon[particles,np.newaxis,np.newaxis], rEarth)
            # store values
            Ncparticles = np.sum(particles)
            enspart[ac] = Ncparticles
            dvarrel = (drdr_sum_rel[te,:]-drdr_sum_rel[ts,:])/(Ncparticles-1)
            dvarabs = (drdr_sum_abs[te,:]-drdr_sum_abs[ts,:])/(Ncparticles-1)
            # compute diffusivity
            def weight_func(data):
                return data
            def signed_log10(data):
                return np.sign(data)*np.log10(np.abs(data))
            kappa_rr[ac] = (0.5*np.mean(dvarabs)/(deltat*24.*60.*60.))
            K_rr[ac]     = (0.5*np.mean(dvarrel)/(deltat*24.*60.*60.))
            ratio[ac]    = K_rr[ac]/kappa_rr[ac]
            kappa_rr[ac] = signed_log10(kappa_rr[ac])
            K_rr[ac]     = signed_log10(K_rr[ac])

        np.savez('layer%d'%(nlayer), x=x,y=y, kappa_rr=kappa_rr, K_rr=K_rr, enspart=enspart)
        print 'finished %d layer' % (nlayer)

        triang = Triangulation(x,y)

        plt.figure()
        triang.plot_scalar(scalar=enspart)
        #plt.clim(2,6)
        plt.colorbar()
        plt.xlim(-20,20)
        plt.title('particles')
        plt.savefig(case+'nparticles' + '_layer%d'%(nlayer)+'.png')
        
        #plt.figure()
        #triang.plot_scalar(scalar=K_rr, cmap=plt.get_cmap('gist_ncar'))
        #plt.colorbar()
        #plt.clim(-7,7)
        #plt.xlim(-20,20)
        #plt.title('relative K_rr radius=' + str(radius/1000) + ' km')
        #plt.savefig(case+'kappa_rr_full.png')
        
        plt.figure()
        triang.plot_scalar(scalar=K_rr)
        plt.colorbar()
        plt.clim(1,6)
        plt.xlim(-20,20)
        plt.title('relative K_rr radius=' + str(radius/1000) + ' km')
        plt.savefig(case+'K_rr' + '_layer%d'%(nlayer)+'.png')
        
        plt.figure()
        triang.plot_scalar(scalar=kappa_rr)
        plt.colorbar()
        plt.clim(1,6)
        plt.xlim(-20,20)
        plt.title('absolute kappa_rr radius=' + str(radius/1000) + ' km')
        plt.savefig(case+'kappa_rr' + '_layer%d'%(nlayer)+'.png')

        fig = plt.figure()
        def plot_triang(c, name):
            triang.plot_scalar(scalar=c)
            plt.colorbar()
            plt.xlim(-20,20)
            plt.title(name + ' radius=' + str(radius/1000) + ' km')
        fig.add_subplot(1,3,1, adjustable='box', aspect='equal')
        plot_triang(K_rr,'relative K_rr')
        plt.clim(1,6)
        fig.add_subplot(1,3,2, adjustable='box', aspect='equal')
        plot_triang(kappa_rr,'absolute kappa_rr')
        plt.clim(1,6)
        fig.add_subplot(1,3,3, adjustable='box', aspect='equal')
        triang.plot_scalar(scalar=ratio)
        plt.title('K_rr/kappa_rr')
        plt.colorbar()
        plt.xlim(-20,20)
        plt.savefig(case+'KcompKappa' + '_layer%d'%(nlayer)+'.png')

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
    
    #original_results()
    build_particle_file(options.inputfilename, float(options.radius), 
            int(options.ts)-1, int(options.te)-1, float(options.dt), int(options.nlayers))
