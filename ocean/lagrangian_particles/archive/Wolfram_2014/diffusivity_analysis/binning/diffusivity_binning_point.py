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
from GeneralTriangulation import GeneralTriangulation as Triangulation
import matplotlib.pyplot as plt
from file_handling import check_folder
from latlon_coordinate_transforms import signed_distances_numexpr as signed_distances, haversine_formula_numexpr as haversine_formula
from cluster_formation import ParticleClusters

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
    
    Ntime = plon.shape[0]
    tvec = deltat*np.arange(0,Ntime)
    Nparticles = plon.shape[1]
    Npartlayers = Nparticles/nlayers
    Nensembles = plon.shape[2]
    rEarth = 6371220. # from netcdf file
    
    lengths = np.linspace(30000,2*radius,Nlen)
    lengths = 2*(np.linspace(15000**0.5,radius**0.5,Nlen))**2.0
    Nradius = lengths.shape[0]
   
    case = 'Case_lon=%f_lat=%f_maxradius=%f_ts=%d_te=%d_deltat=%f_Nradius=%d_Nlen=%d/'%(pointlon, pointlat, radius,ts,te,deltat, Nradius, Nlen)
    print 'Case is ' + case
    check_folder(case)

    # allocate memory #{{{
    ensminpart  = np.zeros((Nradius))
    ensmaxpart  = np.zeros((Nradius))
    ensmeanpart = np.zeros((Nradius))
    K_rr        = np.zeros((Nradius))
    K_xx        = np.zeros((Nradius))
    K_xy        = np.zeros((Nradius))
    K_yy        = np.zeros((Nradius))
    #}}}

    # convert points for radial clusters 
    x = np.radians(pointlon)
    y = np.radians(pointlat)

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

        # reshape into vectors
        llonv = np.reshape(llon,(Npc*Ntc*Nec))
        llatv = np.reshape(llat,(Npc*Ntc*Nec))

        # compute clusters sequentially to conserve on memory
        clusters = ParticleClusters(llonv, llatv, mode='ball')

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        for al in np.arange(Nradius):
            # compute clusters sequentially to conserve on memory
            localparticles = clusters.get_cluster(x, y, radius=lengths[al]/2.0)
            dist = haversine_formula(y, llatv[localparticles], x, llonv[localparticles], rEarth) 
            particles = (dist < lengths[al]/2.0)

            # reshape subset index into prior form
            index = np.zeros((Npc*Ntc*Nec), dtype=bool)
            index[localparticles] = particles
            # get index chape back to match llon,llat
            index = index.reshape((Npc,Ntc,Nec))

            # get number of particles in the cluster (sum over Np and take extreems over Ne)
            ensminpart[al] = np.min(np.sum(index, axis=0)[:])
            ensmaxpart[al] = np.max(np.sum(index, axis=0)[:])
            ensmeanpart[al] = np.mean(np.sum(index, axis=0)[:])

            # allocate memory #{{{
            dvar        = np.zeros((Ntc-1,Nec)) 
            dvarxx        = np.zeros((Ntc-1,Nec)) 
            dvarxy        = np.zeros((Ntc-1,Nec)) 
            dvaryy        = np.zeros((Ntc-1,Nec)) 
            Ncparticles = np.zeros((Ntc-1,Nec)) 
            #}}}

            #plt.figure(1) #{{{
            #plt.figure(2) #}}}
            for at in np.arange(Ntc-1):
                for ae in np.arange(Nec):
                    particles = index[:, at, ae]
                    # need two time steps to compute diffusivity
                    lon = llon[particles,at:at+2,ae]
                    lat = llat[particles,at:at+2,ae]
                    # compute dispersion
                    clon, clat, drdr_sum, dxdx_sum, dxdy_sum, dydy_sum = compute_com_and_squaresums(lat, lon, rEarth)
                    #plt.figure(1) #{{{
                    #plt.plot(lon[:,0],lat[:,0],'k.',alpha=0.5)
                    #plt.plot(lon[:,1],lat[:,1],'r.',alpha=0.5)
                    #plt.figure(2)
                    #plt.plot(lon[:,0]-clon[0],lat[:,0]-clat[0],'k.',alpha=0.5)
                    #plt.plot(lon[:,1]-clon[1],lat[:,1]-clat[1],'r.',alpha=0.5) #}}}
                    # store values
                    Ncparticles[at,ae] = np.sum(particles)
                    dvar[at,ae] = (drdr_sum[1]-drdr_sum[0])/(Ncparticles[at,ae]-1)
                    dvarxx[at,ae] = (dxdx_sum[1]-dxdx_sum[0])/(Ncparticles[at,ae]-1)
                    dvarxy[at,ae] = (dxdy_sum[1]-dxdy_sum[0])/(Ncparticles[at,ae]-1)
                    dvaryy[at,ae] = (dydy_sum[1]-dydy_sum[0])/(Ncparticles[at,ae]-1)
            # compute diffusivity
            def weight_func(data):
                return data
            def signed_log10(data):
                return np.sign(data)*np.log10(np.abs(data))
            K_rr[al] = signed_log10(0.5*np.nansum((dvar[:])*weight_func(Ncparticles[:]))/(np.nansum(weight_func(Ncparticles[:]))) \
                    /(deltat*24.*60.*60.))
            K_xx[al] = signed_log10(0.5*np.nansum((dvarxx[:])*weight_func(Ncparticles[:]))/(np.nansum(weight_func(Ncparticles[:]))) \
                    /(deltat*24.*60.*60.))
            K_xy[al] = signed_log10(0.5*np.nansum((dvarxy[:])*weight_func(Ncparticles[:]))/(np.nansum(weight_func(Ncparticles[:]))) \
                    /(deltat*24.*60.*60.))
            K_yy[al] = signed_log10(0.5*np.nansum((dvaryy[:])*weight_func(Ncparticles[:]))/(np.nansum(weight_func(Ncparticles[:]))) \
                    /(deltat*24.*60.*60.))
            #print K_rr[al] #{{{
            #plt.figure(1); plt.title('abs')
            #plt.figure(2); plt.title('rel')
            #plt.show() #}}}

        np.savez(case+'layer%d'%(nlayer), lenscale=lengths, K_rr=K_rr, K_xx=K_xx, K_xy=K_xy, K_yy=K_yy, 
                ensmin=ensminpart, ensmax=ensmaxpart, ensmean=ensmeanpart)

        def plot_K(data, label):
            plt.figure()
            plt.plot(lengths/1000., data, 'o')
            plt.title('relative '+label+' at (%f,%f)'%(pointlon, pointlat))
            plt.ylabel('log_10 ('+label + ' m^2/s)')
            #plt.ylim(1.0,5.0)
            plt.xlabel('L km')
            plt.savefig(case + label + '_layer%d'%(nlayer) + '.png')
        plot_K(K_rr, 'K_rr')
        plot_K(K_xx, 'K_xx')
        plot_K(K_xy, 'K_xy')
        plot_K(K_yy, 'K_yy')

        def plot_tmix(data, label):
            plt.figure()
            plt.plot(lengths/1000., (lengths**2.0/10**data)/(60*60*24), 'o')
            plt.ylabel(label + ' t_mix days')
            plt.xlabel('L km')
            plt.title('mixing times at (%f,%f)'%(pointlon, pointlat))
            plt.savefig(case+label+'_mixingtime' + '_layer%d'%(nlayer) + '.png')
        plot_tmix(K_rr, 'K_rr')
        plot_tmix(K_xx, 'K_xx')
        plot_tmix(K_xy, 'K_xy')
        plot_tmix(K_yy, 'K_yy')
       
        def plot_scaling(data, label):
            plt.figure()
            plt.plot(lengths/1000., ((10**data)/lengths**(4./3.)), 'o:')
            plt.ylabel("$K/L^{4/3}$")
            plt.xlabel('L km')
            plt.title('scaling at (%f,%f)'%(pointlon, pointlat))
            plt.savefig(case+label+'_scaling' + '_layer%d'%(nlayer) + '.png')
        plot_scaling(K_rr, 'K_rr')
        plot_scaling(K_xx, 'K_xx')
        plot_scaling(K_xy, 'K_xy')
        plot_scaling(K_yy, 'K_yy')

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
   
    print 'starting analysis'
    #original_results()
    build_particle_file(options.inputfilename, float(options.lon), float(options.lat), float(options.radius), float(options.nlen), 
            int(options.ts)-1, int(options.te)-1, float(options.dt), int(options.nlayers), degrees=options.degrees)
