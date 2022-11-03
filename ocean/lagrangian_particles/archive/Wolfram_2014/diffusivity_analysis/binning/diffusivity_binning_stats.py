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
from latlon_coordinate_transforms import signed_distances_numpy as signed_distances, haversine_formula_numpy as haversine_formula
from cluster_formation import ParticleClusters
from compute_diffusivity import compute_diffusivity, compute_com_dispersion
from comp_funcs import SignedLog10
import seaborn as sns; plt.rcdefaults()
import EllipsePlots
from PlotSettings import draw_arrow

def build_particle_file(fname_in, pointlon, pointlat, radius=100000, Nlen=5, ts=25, deltat=2., nlayers=5, layerrange=np.arange(0,5), degrees=False):  #{{{
    
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
    
    case = 'CaseBinningStats_lon=%f_lat=%f_r=%f_ts=%d_deltat=%f_Nlen=%d/'%(pointlon, pointlat, radius, ts, deltat, Nlen)

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
        drdr  = np.zeros((Ntc, Nec)) 
        dxdx  = np.zeros((Ntc, Nec)) 
        dxdy  = np.zeros((Ntc, Nec)) 
        dydy  = np.zeros((Ntc, Nec)) 
        clon = np.zeros((Ntc, Nec)) 
        clat = np.zeros((Ntc, Nec)) 
        Ncparticles   = np.zeros((Ntc, Nec)) 
        #}}}

        # compute clusters sequentially to conserve on memory
        clusters = ParticleClusters(llon[:,0,0], llat[:,0,0], mode='ball')

        # note loop is efficient because memory isn't unduely taxed by doing this in one call (also allows for larger cluster sizes)
        # compute clusters sequentially to conserve on memory
        localparticles = clusters.get_cluster_ball(x, y, radius=radiustree)
        dist = haversine_formula(y, llat[localparticles,0,0], x, llon[localparticles,0,0], rEarth) 
        particles = (dist < radius)
        particles = np.asarray(localparticles)[particles]

        assert np.min(haversine_formula(y, llat[particles,0,0], x, llon[particles,0,0], rEarth) < radius), 'Selection of radius is not correct'

        for dtau in np.arange(Ntc):
            for ae in np.arange(Nec):
                # need two time steps to compute diffusivity
                lon = llon[particles, dtau ,ae]
                lat = llat[particles, dtau ,ae]
                # compute dispersion
                Ncparticles[dtau,ae], clon[dtau,ae], clat[dtau,ae], drdr[dtau,ae], dxdx[dtau,ae], dxdy[dtau,ae], dydy[dtau,ae] = \
                    compute_com_dispersion(lon, lat)

        np.savez(case+'layer%d'%(nlayer), ts=ts, te=ts+Ntc, deltat=deltat,
                lon=np.degrees(llon[particles,:,:]), lat=np.degrees(llat[particles,:,:]), 
                clon=np.degrees(clon), clat=np.degrees(clat),
                drdr=drdr, dxdx=dxdx, dxdy=dxdy, dydy=dydy,
                Ncparticles=Ncparticles, xs=pointlon, ys=pointlat)
        data  = np.load(case+'layer%d.npz'%(nlayer))
        print data.files

        print 'finished %d layer with %d particles' % (nlayer, Ncparticles[0,0]) 

    return  case #}}}

def plot_stats(filename, layerrange, ext=7):

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
        clon=data['clon']; clat=data['clat']
        lon = data['lon']; lat = data['lat']
        xs = data['xs']; ys = data['ys']
        deltat = data['deltat']; ts = data['ts']
        drdr = data['drdr']; dxdx = data['dxdx']; dxdy = data['dxdy']; dydy = data['dydy']

        Np = lon.shape[0]
        Nt = lon.shape[1]
        Ne = lon.shape[2]

        ###################################################
        # spaghetti plots
        ###################################################
        plt.figure()
        # plot paths relative to COM at starting point
        nlon = lon - clon + xs
        nlat = lat - clat + ys
        for i in range(Np):
            for j in range(Ne):
                plt.plot(lon[i,:,j],lat[i,:,j],color=darkgray)
        for i in range(Np):
            for j in range(Ne):
                plt.plot(nlon[i,:,j],nlat[i,:,j],color=medgray)

        # plot center of cluster mass
        for i in range(Ne): plt.plot(clon[:,i],clat[:,i],'k-')
        com = plt.plot(np.mean(clon,axis=1),np.mean(clat,axis=1),'w-', lw=3)

        # plot circle for starting convex hull
        def plot_hull(x,y,*args,**kwargs):
            points = np.vstack((x,y)).T
            hull = ConvexHull(points)
            for simplex in hull.simplices:
                plt.plot(points[simplex,0],points[simplex,1], *args, **kwargs)
            return
        plot_hull(lon[:,0,0],lat[:,0,0],'b-',lw=1)
        plot_hull(nlon[:,-1,:].ravel(),nlat[:,-1,:].ravel(),'b--',lw=2)
        plot_hull(lon[:,-1,:].ravel(),lat[:,-1,:].ravel(),'b:',lw=2)

        savefig('spaghetti.png')

        ###################################################
        # dispersion and diffusivity plots 
        ###################################################
        def plot_D(ax, D,name):
            Nt = D.shape[0]
            Ne = D.shape[1]
            t = (ts + np.arange(Nt))*deltat
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

        ###################################################
        # ellipse plots
        ###################################################
        Kxx = np.mean(np.diff(dxdx,axis=0)/(24.*60*60*deltat), axis=1)
        Kxy = np.mean(np.diff(dxdy,axis=0)/(24.*60*60*deltat), axis=1)
        Kyy = np.mean(np.diff(dydy,axis=0)/(24.*60*60*deltat), axis=1)
        t = (ts + np.arange(Nt))*deltat
        
        plt.figure()
        Nx = 5
        Ny = 6
        degsize = 10
        escale = 1e-4

        plt.xlim(0,Nx*degsize)
        plt.ylim(0,Ny*degsize)
        N = Kxx.shape[0]
        for i in np.arange(Nx):
            for j in np.arange(Ny):
                k = i + Nx*j
                if k >= N:
                    break
                else:
                    cov = np.array([[Kxx[k], Kxy[k]],[Kxy[k], Kyy[k]]])*escale
                    ax = i*degsize + degsize/2.
                    ay = Ny*degsize - (j*degsize + degsize/2.)
                    EllipsePlots.plot_cov_ellipse(ax, ay, cov, 'k-', linewidth=2)
                    plt.text(i*degsize + degsize/10., Ny*degsize-(j*degsize + degsize/4.), 't = %.0f'%(t[k]))

        plt.grid('on')
        plt.gca().get_xaxis().set_ticklabels([])
        plt.gca().get_yaxis().set_ticklabels([])
        #plt.axis('off')

        # ellipse legend #{{{
        lb = 0.12; bb = 0.10; tb = 0.90; rb = 0.765; gap=0.02
        plt.subplots_adjust(left=lb, bottom=bb, right=rb, top=tb, wspace=gap, hspace=gap)
        # in center
        wp = ((rb-lb)-gap/2)/2.0
        hp =((tb-bb)-gap/2)/2.0
        sfx = 4.0
        sfy = 5.0
        iw = wp/sfx
        ih = hp/sfy
       
        xlim = plt.xlim()
        ylim = plt.ylim()
        xc = xlim[0] + 0.5*(xlim[-1]-xlim[0])
        yc = ylim[0] + 0.35*(ylim[-1]-ylim[0])
       
        ixi = 0.5*(rb+lb)-iw/2.0
        iyi = bb + 0.5*(tb-bb)-ih/2.0
        iyi = bb + 1.1*(tb-bb)-1.02*ih
        #iyi = bb + 0.0*(tb-bb)+0.025*ih
        cax = plt.gcf().add_axes([ixi, iyi, iw, ih])
        cax.set_xlim([xc - (xc - xlim[0])/sfx, xc + (xlim[-1]-xc)/sfx])
        cax.set_ylim([yc - (yc - ylim[0])/sfy, yc + (ylim[-1]-yc)/sfy])
        
        cov = np.array([[3.0, 0.0],[0.0, 1.0]])
        xp, yp = EllipsePlots.plot_cov_ellipse(xc, yc, cov, 'k-', linewidth=2)

        width = 0.3*np.sqrt(np.max(yp)-np.min(yp))
        draw_arrow((np.mean(xp), np.min(yp)), (np.mean(xp), np.max(yp)), width, line=True)
        draw_arrow((np.min(xp), np.mean(yp)), (np.max(xp), np.mean(yp)), width, line=True)

        cax.xaxis.set_major_locator(plt.NullLocator())
        cax.yaxis.set_major_locator(plt.NullLocator())
        cax.axis('off')
        #cax.patch.set_color('b')
        #cax.patch.set_alpha(0.3)
        # annotate ellipse
        cax.text(xc + 0.085*(xlim[-1]- xlim[0]), yc, r'$3 \times 10^{'+'%.0f'%(-np.log10(escale)) +r'}$'+ '$m^2s^{-1}$', ha='left', va='center')
        cax.text(xc , yc + 0.05*(ylim[-1] - ylim[0]), r'$1 \times 10^{'+'%.0f'%(-np.log10(escale)) +r'}$'+ '$m^2s^{-1}$', ha='center', va='bottom')
        #}}}
        savefig('tensorellipses.png')

        ####################################################
        ## scatter histograms
        ####################################################
        ##plt.xlim(xs-ext, xs+ext)
        ##plt.ylim(ys-ext, ys+ext)

        #def scatter_hist(x,y):
        #    def count(x,y):
        #        return x.shape[0]
        #    jp = sns.jointplot(x,y, marginal_kws={"fit": stats.norm})
        #    jp.annotate(count, template="{stat} = {val:d}", frameon=False)
        #    sns.despine()
        #
        #scatter_hist(lon[:,0,:],lat[:,0,:])
       
        ##rlzn = int(np.floor(np.random.rand(1)*30))
        ##scatter_hist(lon[:,-1,rlzn],lat[:,-1,rlzn])
        ##plt.suptitle('abs rlzn=%d'%(rlzn))
        ##scatter_hist(nlon[:,-1,rlzn],nlat[:,-1,rlzn])
        ##plt.suptitle('rel rlzn=%d'%(rlzn))
        #
        #scatter_hist(lon[:,-1,:],lat[:,-1,:])
        #plt.suptitle('abs') 
        #scatter_hist(nlon[:,-1,:],nlat[:,-1,:])
        #plt.suptitle('rel')
        
        #plt.show()


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
        #}}}
        print 'starting analysis'
        #original_results()
        options.inputfilename = build_particle_file(options.inputfilename, float(options.lon), float(options.lat), float(options.radius), float(options.nlen), 
                int(options.ts)-1, float(options.dt), int(options.nlayers), degrees=options.degrees, layerrange=options.layer)
    if options.plot == 'T':
        print 'opening %s directory for plotting'%(options.inputfilename)
        plot_stats(options.inputfilename, options.layer)
        cmd = 'open ' + options.inputfilename + '/*png'
        sp.call(cmd,shell=True)

