#!/usr/bin/env python
"""

    Compute decorrelation length scale for a point.
    
    Phillip Wolfram
    LANL
    03/03/2015

"""
# import libraries / packages
import os
import numpy as np
import numexpr as ne
from scipy.spatial import ConvexHull, cKDTree as KDTree
from GeneralTriangulation import GeneralTriangulation as Triangulation
from matplotlib.tri import LinearTriInterpolator
import matplotlib.pyplot as plt
from file_handling import check_folder, file_exists
from latlon_coordinate_transforms import haversine_formula_numexpr as haversine_formula
from cluster_formation import ParticleClusters

# function definitions
def build_particle_file(fname_in, pointlon, pointlat, radius=150000., Nlen=20,
        ts=25, te=30, deltat=2., nlayers=5, layerrange=np.arange(0,5), Ndim=4, degrees=False):  #{{{
    
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
    
    lengths = np.linspace(30000,2*radius,Nlen)
    lengths = 2*(np.linspace(15000**0.5,radius**0.5,Nlen))**2.0
    Nradius = lengths.shape[0]

    # need to remove the mean velocity information on whole dataset
    # build up layers for interpolation of particle layers
    if file_exists('velocity_perturbation.npz'):
        data = np.load('velocity_perturbation.npz','r')
        plonVel = data['plonVel']
        platVel = data['platVel']
        print 'loaded velocity perturbation'
    else:    
        data = np.load('buoyancySurfaceStats.npz','r')
        triang = Triangulation(np.radians(data['lonCell']),np.radians(data['latCell']))
        print 'building up interpolants'
        for alayer in np.arange(nlayers):
            xint = plon[:,alayer*Npartlayers:(alayer+1)*Npartlayers,:].ravel()
            yint = plat[:,alayer*Npartlayers:(alayer+1)*Npartlayers,:].ravel()
            
            meanvel = LinearTriInterpolator(triang, data['velzonal_mean'][:,alayer])
            plonVel[:,alayer*Npartlayers:(alayer+1)*Npartlayers,:] = \
                    np.reshape(\
                    plonVel[:,alayer*Npartlayers:(alayer+1)*Npartlayers,:].ravel() - \
                    meanvel(xint,yint),(Ntime,Npartlayers,Nensembles))

            meanvel = LinearTriInterpolator(triang, data['velmerid_mean'][:,alayer])
            platVel[:,alayer*Npartlayers:(alayer+1)*Npartlayers,:] = \
                    np.reshape(\
                    platVel[:,alayer*Npartlayers:(alayer+1)*Npartlayers,:].ravel() - \
                    meanvel(xint,yint),(Ntime,Npartlayers,Nensembles))
            print 'layer %d finished'%(alayer)
       
        print 'saving velocity perturbations'
        np.savez('velocity_perturbation',plonVel=plonVel, platVel=platVel)
        print 'done saving velocity perturbations'
        # clean memory
        del meanvel, xint, yint, triang, data

    case = 'CaseDecorrelation_lon=%f_lat=%f_maxradius=%f_ts=%d_te=%d_deltat=%f_Nradius=%d_Nlen=%d/'%(pointlon, pointlat, radius,ts,te,deltat, Nradius, Nlen)
    print 'Case is ' + case
    check_folder(case)

    # allocate memory #{{{
    ensminpart  = np.zeros((Nradius))
    ensmaxpart  = np.zeros((Nradius))
    ensmeanpart = np.zeros((Nradius))
    Rlon        = np.zeros((Nradius))
    Rlat        = np.zeros((Nradius))
    eddyspeed = np.zeros((Nradius))
    #}}}

    # convert points for radial clusters 
    x = np.radians(pointlon)
    y = np.radians(pointlat)

    # starting center points of cluster
    for nlayer in layerrange: 

        # compute the clusters and compute dispersion in a single layer
        # need to rearrange data so that we have total realizations and cluster as index (with averaging over total realizations)
        # now it is of size (Np, Nt, Ne)
        llon = np.swapaxes(plon[[ts,te],nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llat = np.swapaxes(plat[[ts,te],nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llonVel = np.swapaxes(plonVel[[ts,te],nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)
        llatVel = np.swapaxes(platVel[[ts,te],nlayer*Npartlayers:(nlayer+1)*Npartlayers,:],0,1)

        Npc=llon.shape[0]
        Ntc=llon.shape[1]
        Nec=llon.shape[2]

        # reshape into vectors
        llonv = np.reshape(llon[:,0,:],(Npc*Nec))
        llatv = np.reshape(llat[:,0,:],(Npc*Nec))

        del llon, llat

        # compute clusters sequentially to conserve on memory
        clusters = ParticleClusters(llonv, llatv, mode='ball')

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        for al in np.arange(Nradius)[-1:0:-1]:
            # compute clusters sequentially to conserve on memory
            localparticles = clusters.get_cluster(x, y, radius=lengths[al]/2.0)
            dist = haversine_formula(y, llatv[localparticles], x, llonv[localparticles], rEarth) 
            particles = (dist < lengths[al]/2.0)

            # reshape subset index into prior form
            index = np.zeros((Npc*Nec), dtype=bool)
            index[localparticles] = particles
            # get index chape back to match llon,llat
            index = index.reshape((Npc,Nec))

            # get number of particles in the cluster (sum over Np and take extreems over Ne)
            ensminpart[al] = np.min(np.sum(index, axis=0)[:])
            ensmaxpart[al] = np.max(np.sum(index, axis=0)[:])
            ensmeanpart[al] = np.mean(np.sum(index, axis=0)[:])

            # allocate memory #{{{
            u0ulon          = np.nan*np.zeros((Npc,2,Nec)) 
            u0ulat          = np.nan*np.zeros((Npc,2,Nec)) 
            #}}}

            #plt.figure(1)
            #plt.hold(True)
            #plt.figure(2)
            #plt.hold(True)
            for ae in np.arange(Nec):
                particles = index[:, ae]
                def perturbation(vel):
                    return vel - np.mean(vel,axis=0)
                lonVel = (llonVel[particles,:,ae])
                latVel = (llatVel[particles,:,ae])

                u0ulon[particles,:,ae] = (lonVel[:,0])[:,np.newaxis]*lonVel
                u0ulat[particles,:,ae] = (latVel[:,0])[:,np.newaxis]*latVel
                #plt.figure(1)
                #plt.plot(llonVel[particles,0,ae],llonVel[particles,1,ae],'.')
                #plt.figure(2)
                #plt.plot(lonVel[:,0],lonVel[:,1],'.')

            #plt.show()
            #import pdb; pdb.set_trace()
            # compute autocorrelation
            Rlon[al] = np.nanmean(u0ulon[:,1,:].ravel())/np.nanmean(u0ulon[:,0,:].ravel())
            Rlat[al] = np.nanmean(u0ulat[:,1,:].ravel())/np.nanmean(u0ulat[:,0,:].ravel())
            eddyspeed[al] = np.std(np.sqrt(lonVel[:,0]*lonVel[:,0] + latVel[:,0]*latVel[:,0]))

        np.savez(case+'layer%d'%(nlayer), lenscale=lengths, Rlon=Rlon, Rlat=Rlat,eddyspeed=eddyspeed,
                ensmin=ensminpart, ensmax=ensmaxpart, ensmean=ensmeanpart)

        def plot_R(data, label):
            plt.figure()
            plt.plot(lengths/1000., data, 'o')
            plt.title('relative '+label+' at (%f,%f)'%(pointlon, pointlat))
            plt.ylabel('R')
            plt.ylim(0.0,1.0)
            plt.xlabel('L km')
            plt.savefig(case + label + '_layer%d'%(nlayer) + '.png')
        def plot_u(data, label):
            plt.figure()
            plt.plot(lengths/1000., data, 'o')
            plt.title('relative '+label+' at (%f,%f)'%(pointlon, pointlat))
            plt.ylabel('u_{rms} (m/s)')
            plt.ylim(0.0,1.0)
            plt.xlabel('L km')
            plt.savefig(case + label + '_layer%d'%(nlayer) + '.png')
        plot_R(Rlat, 'Rlat')
        plot_R(Rlon, 'Rlon')
        plot_u(eddyspeed, 'eddyspeed')
        
        plt.figure()
        plt.plot(lengths/1000., ensminpart, label='min')
        plt.plot(lengths/1000., ensmeanpart, label='mean')
        plt.plot(lengths/1000., ensmaxpart, label='max')
        plt.legend(loc='best')
        plt.ylabel('Number particles')
        plt.xlabel('L km')
        plt.title('cluster numbers at (%f,%f)'%(pointlon, pointlat))
        plt.savefig(case+'clusterstats' + '_layer%d'%(nlayer) + '.png')
        
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
    parser.add_option("-y", "--latitude", dest="lat",
            help="latitude in degress", metavar="FLOAT")
    parser.add_option("-x", "--longitude", dest="lon",
            help="longitude in degress", metavar="FLOAT")
    parser.add_option("-r", "--radius", dest="radius",
            help="radius of the cluster in km", metavar="FLOAT")
    parser.add_option("-l", "--numlen", dest="nlen",
            help="number of cluster lengths", metavar="FLOAT")
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
    if not options.radius:
        parser.error("Need starting number for radius.")
    if not options.nlen:
        parser.error("Need number of points for length discretization.")
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
   
    print 'starting analysis (should be in simluation root directory)'
    #original_results()
    build_particle_file(options.inputfilename, float(options.lon), float(options.lat), float(options.radius), float(options.nlen), 
            int(options.ts)-1, int(options.te)-1, float(options.dt), int(options.nlayers), degrees=options.degrees)
