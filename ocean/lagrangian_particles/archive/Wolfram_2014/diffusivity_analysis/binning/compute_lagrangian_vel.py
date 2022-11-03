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
from file_handling import check_folder, file_exists
from latlon_coordinate_transforms import haversine_formula_numexpr as haversine_formula
from cluster_formation import ParticleClusters

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
   
    # cache the cluster information
    if file_exists('coarse_sampling.npz'):
        data = np.load('coarse_sampling.npz')
        x = data['x']
        y = data['y']
    else:
        # use existing particle locations
        x = plon[0,:Nparticles/nlayers,0]
        y = plat[0,:Nparticles/nlayers,0]

        tri = Triangulation(x,y).coarsen()
        x = tri.x
        y = tri.y
        np.savez('coarse_sampling.npz',x=x,y=y)

    #plt.plot(x,y,'.')
    #plt.show()
    Nclusters = x.shape[0]
    case = 'CaseEddySpeed_r=%f_ts=%d_te=%d_deltat=%f_coarsenN=%d/'%(radius,ts,te,deltat, Nclusters)

    print 'Nclusters = %d' % (Nclusters)
    print 'Case is ' + case
    check_folder(case)

    # compute speed
    speed = np.swapaxes(haversine_formula(plat[1:,:,:],plat[:-1,:,:],plon[1:,:,:],plon[:-1,:,:], rEarth),0,1)/(deltat*24.*60.*60.)

    # allocate memory #{{{
    ensminpart = np.zeros((Nclusters))
    ensmaxpart = np.zeros((Nclusters))
    ensmeanpart = np.zeros((Nclusters))
    eddyspeed   = np.zeros((Nclusters))
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
        print 'initializing clusters'
        clusters = ParticleClusters(llonv, llatv, mode='ball')
        #clusters = ParticleClusters(llonv, llatv, xloc=x, yloc=y, radius=radiustree, mode='tree')
        print 'done, computing values at each cluster'

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        for ac in np.arange(Nclusters):
            # compute clusters all at once
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

            eddyspeed[ac] = np.std(speed[index[:,:-1,:]].ravel())

        np.savez(case+'eddyspeed_layer%d'%(nlayer), x=np.degrees(x),y=np.degrees(y),
                    eddyspeed=eddyspeed, ensmin=ensminpart, ensmax=ensmaxpart, ensmean=ensmeanpart)

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
        triang.plot_scalar(scalar=eddyspeed)
        plt.colorbar()
        plt.clim(0,1.0)
        plt.xlim(-20,20)
        plt.title('eddy speed radius=' + str(radius/1000) + ' km')
        plt.savefig(case+'eddyspeed' + '_layer%d'%(nlayer) + '.png')
        
        
        print 'finished %d layer' % (nlayer)

    #plt.show()
            
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
            int(options.ts)-1, int(options.te)-1, float(options.dt), int(options.nlayers), degrees=options.degrees, layerrange=options.layer)
