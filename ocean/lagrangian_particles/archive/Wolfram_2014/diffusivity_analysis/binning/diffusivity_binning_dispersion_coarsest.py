#!/usr/bin/env python
"""

    Compute diffusivity via binning procedure
    
    Phillip Wolfram
    LANL
    03/03/2015

"""
## plot on mustang
#import matplotlib as mpl
#mpl.use('Agg')
# import libraries / packages
import os
import subprocess as sp
import numpy as np
import numexpr as ne
from scipy.spatial import ConvexHull, cKDTree as KDTree
from GeneralTriangulation import GeneralTriangulation as Triangulation
import matplotlib.pyplot as plt
import scipy.stats as stats
from file_handling import check_folder
from latlon_coordinate_transforms import signed_distances_numpy as signed_distances, haversine_formula_numexpr as haversine_formula
from cluster_formation import ParticleClusters
from compute_diffusivity import compute_dispersion_numexpr
from comp_funcs import SignedLog10
import seaborn as sns; plt.rcdefaults()

def build_particle_file(fname_in, radius=100000, Nlen=5, ts=25, deltat=2., nlayers=5, layerrange=np.arange(0,5), degrees=False):  #{{{
    
    # load the file database
    f_in = np.load(fname_in, 'r')
    if degrees: 
        lonname = 'lon'
        latname = 'lat'
        conv = np.radians
    else:
        lonname = 'lonrad'
        latname = 'latrad'
        def creturn(c): 
            return c
        conv = creturn
    
    Ntime = f_in[lonname].shape[0]
    tvec = deltat*np.arange(0,Ntime)
    Nparticles = f_in[lonname].shape[1]
    Npartlayers = Nparticles/nlayers
    Nensembles = f_in[lonname].shape[2]
    rEarth = 6371220. # from netcdf file
    
    data  = np.load('coarse_sampling.npz')
    #data  = np.load('coarser_sampling550.npz')
    #data  = np.load('tensor_points.npz')
    x = data['x']
    y = data['y']
    Nc = x.shape[0]

    case = 'CaseBinningEnsembledDispersion_r=%f_ts=%d_deltat=%f_Nlen=%d_Nc=%d/'%(radius, ts, deltat, Nlen, Nc)

    print 'Case is ' + case
    check_folder(case)


    # starting center points of cluster
    for nlayer in layerrange: 

        # compute the clusters and dispersion in a single layer
        # need to rearrange data so that we have total realizations and cluster as index (with averaging over total realizations)
        # now it is of size (Np, Nt, Ne)
        # could make sure this is the only data loaded from f_in['lon*'] etc to keep memory down
        llon = np.swapaxes(conv(f_in[lonname][ts:ts+Nlen+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:]),0,1)
        llat = np.swapaxes(conv(f_in[latname][ts:ts+Nlen+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:]),0,1)
        
        # allocate memory #{{{
        drdr  = np.zeros((Nc, Nlen)) 
        dxdx  = np.zeros((Nc, Nlen)) 
        dxdy  = np.zeros((Nc, Nlen)) 
        dydy  = np.zeros((Nc, Nlen)) 
        npart = np.zeros((Nc))
        #}}}

        particles = haversine_formula(y[:,np.newaxis],
                llat[np.newaxis,:,0,0], x[:,np.newaxis],
                llon[np.newaxis,:,0,0], rEarth) < radius

        # should check to make sure this is ok
        #assert np.min(haversine_formula(y, llat[particles,0,0], x, llon[particles,0,0], rEarth) < radius), 'Selection of radius is not correct'

        for acluster in np.arange(Nc):
            #print 'starting %d of %d clusters'%(acluster,Nc)
            # need two time steps to compute diffusivity
            lon = llon[particles[acluster,:], :, :]
            lat = llat[particles[acluster,:], :, :]
            # compute dispersion
            npart[acluster], drdr[acluster,:], dxdx[acluster,:], dxdy[acluster,:], dydy[acluster,:] = \
                compute_dispersion_numexpr(lat, lon, rEarth, ensmean=(0,2))

        np.savez(case+'layer%d'%(nlayer), ts=ts, te=ts+Nlen, deltat=deltat,
                drdr=drdr, dxdx=dxdx, dxdy=dxdy, dydy=dydy, npart=npart,
                x=x, y=y)

        data  = np.load(case+'layer%d.npz'%(nlayer))
        print data.files

        print 'finished %d layer with %d clusters' % (nlayer, Nc)

    return  case #}}}

def plot_stats(filename, layerrange): #{{{

    # color def #{{{
    darkgray = '#939393'
    medgray  = '#e9e9e9'
    lightgray = '#f6f6f6'
    #}}}


    for layer in layerrange:
        print 'on layer %d'%(layer)
        plt.close('all')
        saveprefix = filename+'/layer%d_'%(layer)
        def savefig(name):
            print 'saving %s'%(saveprefix + name)
            plt.savefig(saveprefix + name)

        ###################################################
        # load data
        ###################################################
        data = np.load(filename+'/layer%d.npz'%(layer))
        xs = data['xs']; ys = data['ys']
        deltat = data['deltat']
        drdr = data['drdr']; dxdx = data['dxdx']; dxdy = data['dxdy']; dydy = data['dydy']

        Nc = lon.shape[0]
        Nt = lon.shape[1]

        ###################################################
        # dispersion and diffusivity plots 
        ###################################################
        def plot_D(ax, D,name):
            Nt = D.shape[0]
            Ne = D.shape[1]
            t = np.arange(Nt)*deltat
            for i in np.arange(Ne):
                ax.plot(t, D[:,i], '-', color=darkgray)
            ax.plot(t,np.mean(D,axis=1),'k-',lw=3)
            ax.set_ylabel(name)
            ax.plot(t,0*t,'r-',lw=.5)
            #ax.grid('on')

        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True, figsize = (5,15))
        plot_D(ax1, drdr, r'$D^2_{rr}$')
        plot_D(ax2, dxdx, r'$D^2_{xx}$')
        plot_D(ax3, dxdy, r'$D^2_{xy}$')
        plot_D(ax4, dydy, r'$D^2_{yy}$')
        ymax = np.max([np.max(np.abs(ax.get_ylim())) for ax in [ax1,ax2,ax3,ax4]])
        for ax in [ax1,ax2,ax3,ax4]:
            ax.set_ylim(-ymax,ymax)
        plt.xlabel('time (days)')
        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        f.subplots_adjust(left=0.25)
        savefig('dispersion.png')

        def plot_K(ax, D, name):
            K = SignedLog10.convert(np.diff(D,axis=0)/(24.*60*60*deltat))
            Nt = K.shape[0]
            Ne = K.shape[1]
            t = np.arange(Nt)*deltat
            for i in np.arange(Ne):
                ax.plot(t, K[:,i], '-', color=darkgray)
            ax.plot(t,np.mean(K,axis=1),'k-',lw=3)
            ax.plot(t,0*t,'r-',lw=.5)
            ax.set_ylabel(SignedLog10.labelwrapper(name))

        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True, figsize = (5,15))
        plot_K(ax1, drdr, r'$K^2_{rr}$')
        plot_K(ax2, dxdx, r'$K^2_{xx}$')
        plot_K(ax3, dxdy, r'$K^2_{xy}$')
        plot_K(ax4, dydy, r'$K^2_{yy}$')
        ymax = np.max([np.max(np.abs(ax.get_ylim())) for ax in [ax1,ax2,ax3,ax4]])
        for ax in [ax1,ax2,ax3,ax4]:
            ax.set_ylim(-ymax,ymax)
        plt.xlabel('time (days)')
        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        f.subplots_adjust(left=0.15)
        savefig('diffusivity.png')

        return #}}}

if __name__ == "__main__":
    # analysis mode: time python  diffusivity_binning_stats.py -f all_rlzn_particle_data_rads.npz --ts 1 -n 31  --dt 2. -x 10.0 -y 32.5 -d F -r 100000. 
    # plotting mode: time python  diffusivity_binning_stats.py -f ./CaseBinningStats_lon\=10.000000_lat\=32.500000_r\=100000.000000_ts\=0_deltat\=2.000000_Nlen\=31/ -l 3 -p T

    from optparse import OptionParser
    parser = OptionParser() #{{{
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension or for analysis",
                      metavar="FILE")
    parser.add_option("-p", "--plot", dest="plot",
                     help="plot a file to analyze to build plots?", metavar="BOOL")
    parser.add_option("-r", "--radius", dest="radius",
            help="radius of the cluster in km", metavar="FLOAT")
    parser.add_option("--ts", "--startint", dest="ts",
            help="starting time index", metavar="INT")
    parser.add_option("--dt", "--tsampfreq", dest="dt",
            help="time step (sampling interval) in days", metavar="FLOAT")
    parser.add_option("--nlayers", dest="nlayers",
            help="number of particle layers", metavar="INT")
    parser.add_option('-c', "--compute", dest="compute", metavar="BOOL")
    parser.add_option("-d", "--degrees", dest="degrees", help="Data in degrees? T or F", metavar="BOOL")
    parser.add_option("-l", "--layer", dest="layer", help="Layer number", metavar="INT")
    parser.add_option("-n", "--numlen", dest="nlen",
            help="number of dtau derivatives", metavar="FLOAT")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename ('-f') is a required input.")
    if not options.layer:
      options.layer = np.arange(5)
    else:
      try:
        options.layer = eval(options.layer)
      except:
        options.layer = np.array([int(options.layer)])
    if not options.compute:
        options.compute == 'F'
    if not options.plot:
        options.plot == 'F'
    #}}}
        #}}}

    if options.compute == 'T':
        # parser options #{{{
        if not options.radius:
            parser.error("Need starting number for radius.")
        if not options.nlen:
            parser.error("Need number of points for dtau derivatives.")
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
        #}}}
        print 'starting analysis'
        #original_results()
        options.inputfilename = build_particle_file(options.inputfilename, float(options.radius), float(options.nlen), 
                int(options.ts)-1, float(options.dt), int(options.nlayers), degrees=options.degrees, layerrange=options.layer)
    if options.plot == 'T':
        print 'opening %s directory for plotting'%(options.inputfilename)
        plot_stats(options.inputfilename, options.layer)
        cmd = 'open ' + options.inputfilename + '/*png'
        sp.call(cmd,shell=True)

