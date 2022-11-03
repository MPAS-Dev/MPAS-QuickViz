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
from latlon_coordinate_transforms import spherical_bearing_numpy as spherical_bearing, haversine_formula_numpy as haversine_formula
from cluster_formation import ParticleClusters

# function definitions
def compute_com_and_dispersion(plat, plon, r): #{{{
    """
    compute relative particle disperison (2nd moment) for cluster, using lat / lon for basis of 
    calculation to set coordinates for Kxx,Kxy,Kyy where x is zonal and y meridional

    dispersion units in m^2

    Phillip Wolfram
    LANL
    07/18/2014
    """

    # compute center of mass for each time
    clat = np.nanmean(plat, axis=0)
    clon = np.nanmean(plon, axis=0)

    # scalar diffusivity
    dr = haversine_formula(clat, plat, clon, plon, r)
    drdr_sum = np.nansum(dr*dr, axis=0)

    return clon, clat, drdr_sum #}}}

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
   
    #lonmin = np.min(plon[:])
    #lonmax = np.max(plon[:])

    latmin = np.min(plat[:])
    latmax = np.max(plon[:])

    # compute tree radius
    radiustree = radius / (rEarth * np.sin(np.maximum(np.abs(latmin),np.abs(latmax))))
    
    ## build the grid
    #x,y = np.meshgrid(np.linspace(-15,15,Ndim),np.linspace(23,47,Ndim))
    #dist = np.sqrt((13./16.*x)**2 + (y-35)**2.0)
    #indomain = np.where(dist < 13.0)
    #x = np.radians((x[indomain]).flatten())
    #y = np.radians((y[indomain]).flatten())
    #case = 'Case_r=%f_ts=%d_te=%d_deltat=%f_Ndim=%d/'%(radius,ts,te,deltat,Ndim)
    #Nclusters = x.ravel().shape[0]
    
    # use existing particle locations
    x = plon[0,:Nparticles/nlayers,0]
    y = plat[0,:Nparticles/nlayers,0]
    Nclusters = x.ravel().shape[0]
    
    tri = Triangulation(x,y).coarsen()
    x = tri.x
    y = tri.y
    #plt.plot(x,y,'.')
    #plt.show()
    Nclusters = x.shape[0]
    case = 'Case_r=%f_ts=%d_te=%d_deltat=%f_coarsenN=%d/'%(radius,ts,te,deltat, Nclusters)

    print 'Nclusters = %d' % (Nclusters)
    print 'Case is ' + case
    check_folder(case)

    # allocate memory #{{{
    ensminpart = np.zeros((Nclusters))
    ensmaxpart = np.zeros((Nclusters))
    ensmeanpart = np.zeros((Nclusters))
    K_rr       = np.zeros((Nclusters))
    #}}}
    
    # starting center points of cluster
    #for nlayer in layerrange[-1::-1]: 
    for nlayer in layerrange: 

        # compute the clusters and compute dispersion in a single layer
        # need to rearrange data so that we have total realizations and cluster as index (with averaging over total realizations)
        # now it is of size (Np, Nt, Ne)
        llon = np.swapaxes(plon[ts:te+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llat = np.swapaxes(plat[ts:te+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)

        #hull = ConvexHull(np.vstack((llon[:,0,0],llat[:,0,0])).T)

        # previous cluster assignments #{{{
        ## custom defined clusters (number of particles) 
        #tree = KDTree(np.vstack((llon[0,:,0],llat[0,:,0])).T)
        #_, allparticles = tree.query(np.vstack((x.ravel(),y.ravel())).T, k=kcluster)
        
        ## cell clusters (somewhat of a hack)
        #x = np.degrees(np.mean(llon[0,allparticles[:Nclusters,:],0],axis=1))
        #y = np.degrees(np.mean(llat[0,allparticles[:Nclusters,:],0],axis=1))
        #}}}
        
        Npc=llon.shape[0]
        Ntc=llon.shape[1]
        Nec=llon.shape[2]

        # reshape into vectors
        llonv = np.reshape(llon,(Npc*Ntc*Nec))
        llatv = np.reshape(llat,(Npc*Ntc*Nec))

        ## compute clusters all at once
        #clusters = ParticleClusters(llonv, llatv, x.ravel(), y.ravel(), radius=radiustree, mode='tree')

        # compute clusters sequentially to conserve on memory
        clusters = ParticleClusters(llonv, llatv, mode='ball')

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        for ac in np.arange(Nclusters):
            ## compute clusters all at once
            #localparticles = clusters.clusters[ac]
            # compute clusters sequentially to conserve on memory
            localparticles = clusters.get_cluster(x[ac], y[ac], radius=radiustree)
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

            # allocate memory #{{{
            dvar        = np.zeros((Ntc-1,Nec)) 
            Ncparticles = np.zeros((Ntc-1,Nec)) 
            #}}}

            for at in np.arange(Ntc-1):
                for ae in np.arange(Nec):
                    particles = index[:, at, ae]
                    # need two time steps to compute diffusivity
                    lon = llon[particles,at:at+2,ae]
                    lat = llat[particles,at:at+2,ae]
                    # compute dispersion
                    _, _, drdr_sum = compute_com_and_dispersion(lat, lon, rEarth)
                    # store values
                    Ncparticles[at,ae] = np.sum(particles)
                    dvar[at,ae] = (drdr_sum[1]-drdr_sum[0])/(Ncparticles[at,ae]-1)
            # compute diffusivity
            def weight_func(data):
                return data
            def signed_log10(data):
                return np.sign(data)*np.log10(np.abs(data))
            K_rr[ac] = signed_log10(0.5*np.nansum((dvar[:])*weight_func(Ncparticles[:]))/(np.nansum(weight_func(Ncparticles[:]))) \
                    /(deltat*24.*60.*60.))

        np.savez(case+'layer%d'%(nlayer), x=np.degrees(x),y=np.degrees(y),
                    K_rr=K_rr, ensmin=ensminpart, ensmax=ensmaxpart,
                    ensmean=ensmeanpart, dvar=dvar, Ncparticles=Ncparticles)

        triang = Triangulation(np.degrees(x),np.degrees(y))

        plt.figure()
        triang.plot_scalar(scalar=ensminpart)
        #plt.clim(2,6)
        plt.colorbar()
        plt.xlim(-20,20)
        plt.title('min particles')
        plt.savefig(case + 'min' + '_layer%d'%(nlayer) + '.png')
        
        plt.figure()
        triang.plot_scalar(scalar=ensmaxpart)
        #plt.clim(2,6)
        plt.colorbar()
        plt.xlim(-20,20)
        plt.title('max particles')
        plt.savefig(case+'max' + '_layer%d'%(nlayer) + '.png')
        
        plt.figure()
        triang.plot_scalar(scalar=ensmeanpart)
        #plt.clim(2,6)
        plt.colorbar()
        plt.xlim(-20,20)
        plt.title('mean particles')
        plt.savefig(case+'mean' + '_layer%d'%(nlayer) + '.png')
        
        plt.figure()
        triang.plot_scalar(scalar=K_rr, cmap=plt.get_cmap('gist_ncar'))
        plt.colorbar()
        plt.clim(-7,7)
        plt.xlim(-20,20)
        plt.title('relative K_rr radius=' + str(radius/1000) + ' km')
        plt.savefig(case+'kappa_rr_full' + '_layer%d'%(nlayer) + '.png')
        
        plt.figure()
        triang.plot_scalar(scalar=K_rr)
        plt.colorbar()
        plt.clim(1,5)
        plt.xlim(-20,20)
        plt.title('relative K_rr radius=' + str(radius/1000) + ' km')
        plt.savefig(case+'kappa_rr' + '_layer%d'%(nlayer) + '.png')
        
        print 'finished %d layer' % (nlayer)

    #plt.show()
            
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
            int(options.ts)-1, int(options.te)-1, float(options.dt), int(options.nlayers), degrees=options.degrees, layerrange=options.layer)
