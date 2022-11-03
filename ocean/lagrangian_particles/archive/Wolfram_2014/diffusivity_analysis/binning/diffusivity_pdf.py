#!/usr/bin/env python
"""

    Compute particle pair statistic pdfs and decorrelation via cluster procedure

    Phillip Wolfram
    LANL
    04/01/2015

"""
## plot on mustang
#import matplotlib as mpl
#mpl.use('Agg')
# import libraries / packages
import numpy as np
from file_handling import check_folder, savenpz
from latlon_coordinate_transforms import signed_distances_numexpr as signed_distances, haversine_formula_numexpr as haversine_formula
from compute_diffusivity import compute_dispersion_numexpr, particle_pair_pdfs
# depreciated packages #{{{
#import os
#import subprocess as sp
#import numexpr as ne
#from scipy.spatial import ConvexHull, cKDTree as KDTree
#from cluster_formation import ParticleClusters
#from comp_funcs import SignedLog10
#import scipy.stats as stats
#from GeneralTriangulation import GeneralTriangulation as Triangulation
#import matplotlib.pyplot as plt
#import seaborn as sns; plt.rcdefaults()
#}}}

########################################
# helper functions 
########################################

# compute differences in error between each method
def perc_diff(A,B): return 100.*(A-B)/(0.5*(A+B))

def print_diff(name, acluster, Kclust, Kpair, tsd=5): #{{{
    Cmean = np.mean(Kclust, axis=-1)
    Pmean = np.mean(Kpair[acluster,:,:],axis=-1)
    percmean = perc_diff(Cmean, Pmean)
    print name + ' comparison for step %d between cluster and pairs:'%(tsd),
    print ' cluster = %.2e, pair = %.2e, diff=%.2e%%, maxdiff=%.2e%%'%\
            (Cmean[tsd], Pmean[tsd], percmean[tsd], \
            np.max(np.abs(perc_diff(Kclust,Kpair[acluster,:,:])[:])))
    return perc_diff(Kclust,Kpair[acluster,:,:]) #}}}

def check_pair_equality(name, acluster, lpdf, kpdf, k, tsd=5): #{{{
    pdfcomp = np.nansum(lpdf[acluster,:,:,:]*kpdf[acluster,:,:,:],axis=0)/np.nansum(lpdf[acluster,:,:,:],axis=0)
    pdfcompmean = np.mean(pdfcomp[tsd,:],axis=0)
    kmean = np.mean(k[acluster,tsd,:])
    print name + ' comparison for step %d between pdf and pairs:'%(tsd),
    print ' pairs=%.2e, pdf=%.2e, diff=%.2e%%, maxdiff=%.2e%%'%\
            (kmean, pdfcompmean, perc_diff(kmean,pdfcompmean), \
            np.max(np.abs(perc_diff(k[acluster,:,:],pdfcomp))))
    return  #}}}


def build_particle_file(particledata, samplepoints, lvec, radius=100000,
        Nlen=5, ts=25, deltat=2., nlayers=5, layerrange=np.arange(0,5),
        degrees=False, fraction=1.0):  #{{{

    # load the file database
    f_in = np.load(particledata, 'r')
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

    data  = np.load(samplepoints)
    x = data['x']
    y = data['y']
    Nc = x.shape[0]

    case = 'CaseClusterPDFsComplete_r=%0.2f_ts=%d_deltat=%0.2f_Nlen=%d_Nc=%d_frac=%0.2f/'%(radius, ts, deltat, Nlen, Nc, fraction)

    print 'Case is ' + case
    check_folder(case)

    # starting center points of cluster
    for nlayer in layerrange:

        llon = np.swapaxes(conv(f_in[lonname][ts:ts+Nlen+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:]),0,1)
        llat = np.swapaxes(conv(f_in[latname][ts:ts+Nlen+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:]),0,1)

        # compute velocity
        velx, vely = signed_distances(llat[:,1:,:], llat[:,:-1,:], llon[:,1:,:], llon[:,:-1,:], rEarth)
        velx /= (24.*60.*60.*deltat)
        vely /= (24.*60.*60.*deltat)

        # compute particle clusters
        particles = haversine_formula(y[:,np.newaxis],
                llat[np.newaxis,:,0,0], x[:,np.newaxis],
                llon[np.newaxis,:,0,0], rEarth) < radius

        ## should check to make sure this is ok if results strange #{{{
        #for acluster in np.arange(Nc):
        #    assert np.min(haversine_formula(y[acluster], \
        #        llat[particles[acluster,:],0,0], x[acluster], \
        #        llon[particles[acluster,:],0,0], rEarth) < radius), \
        #        'Selection of particles within radius is not correct'
        #}}}

        # allocate memory #{{{
        Nt = llon.shape[1]-1
        Ne = llon.shape[2]
        clat  = np.zeros((Nc, Nlen, Ne))
        clon  = np.zeros((Nc, Nlen, Ne))
        drdr  = np.zeros((Nc, Nlen, Ne))
        dxdx  = np.zeros((Nc, Nlen, Ne))
        dxdy  = np.zeros((Nc, Nlen, Ne))
        dydy  = np.zeros((Nc, Nlen, Ne))
        npart = np.zeros((Nc))
        meanu = np.zeros((Nc, Nt))
        eddyu = np.zeros((Nc, Nt))
        meanv = np.zeros((Nc, Nt))
        eddyv = np.zeros((Nc, Nt))
        rho   = np.zeros((Nc, Nt))
        rhou  = np.zeros((Nc, Nt))
        rhov  = np.zeros((Nc, Nt))
        meanspeed = np.zeros((Nc, Nt))
        eddyspeed = np.zeros((Nc, Nt))
        npart = np.zeros((Nc))

        # for construction of PDFs
        Nl = len(lvec)-1
        #print 'Nl=',Nl
        #print 'Nt=',Nt
        #print 'Nlen=',Nlen
        #print 'Ne=',Ne

        # n and n+1 terms
        # mean pair location
        xmPDF  = np.zeros((Nc,Nl,Nlen,Ne))
        ymPDF  = np.zeros((Nc,Nl,Nlen,Ne))
        # pair separation
        dxPDF  = np.zeros((Nc,Nl,Nlen,Ne))
        dyPDF  = np.zeros((Nc,Nl,Nlen,Ne))
        lPDF   = np.zeros((Nc,Nl,Nlen,Ne))
        thetaPDF=np.zeros((Nc,Nl,Nlen,Ne))
        # n+1/2 terms
        lmPDF  = np.zeros((Nc,Nl,Nt,Ne))
        # eddy velocity
        sePDF  = np.zeros((Nc,Nl,Nt,Ne))
        uePDF  = np.zeros((Nc,Nl,Nt,Ne))
        vePDF  = np.zeros((Nc,Nl,Nt,Ne))
        # mean velocity
        smPDF  = np.zeros((Nc,Nl,Nt,Ne))
        umPDF  = np.zeros((Nc,Nl,Nt,Ne))
        vmPDF  = np.zeros((Nc,Nl,Nt,Ne))
        # diffusivity
        KxxPDF = np.zeros((Nc,Nl,Nt,Ne))
        KxyPDF = np.zeros((Nc,Nl,Nt,Ne))
        KyyPDF = np.zeros((Nc,Nl,Nt,Ne))
        KrrPDF = np.zeros((Nc,Nl,Nt,Ne))
        KxxMean = np.zeros((Nc,Nt,Ne))
        KxyMean = np.zeros((Nc,Nt,Ne))
        KyyMean = np.zeros((Nc,Nt,Ne))
        KrrMean = np.zeros((Nc,Nt,Ne))
        npartpairs = np.zeros((Nc))
        npairs = np.zeros((Nc))

        # differences between particle pair and cluster methods
        percKrr = np.zeros((Nc,Nt,Ne))
        percKxx = np.zeros((Nc,Nt,Ne))
        percKxy = np.zeros((Nc,Nt,Ne))
        percKyy = np.zeros((Nc,Nt,Ne))
        #}}}

        for acluster in np.arange(Nc): #{{{
            # get particles in the cluster
            clusterpart = particles[acluster,:]

            # subsample particles
            clustids = np.where(clusterpart)[0]
            Ncorig = np.sum(clusterpart)
            # remove non-desired random fraction from cluster
            removeids = np.random.choice(clustids, np.floor(clustids.shape[0]*(1-fraction)), replace=False)
            clusterpart[removeids] = False
            # checks to make sure subsampling was correct
            Ncsample = np.sum(clusterpart)
            assert np.abs(Ncsample - np.round(fraction*Ncorig)) < 2, 'Error with subsampling'

            # compute the eddy speed and mean speed
            meanu[acluster,:] = np.mean(velx[clusterpart,:,:],axis=(0,2))
            meanv[acluster,:] = np.mean(vely[clusterpart,:,:],axis=(0,2))
            eddyu[acluster,:] = np.std(velx[clusterpart,:,:],axis=(0,2))
            eddyv[acluster,:] = np.std(vely[clusterpart,:,:],axis=(0,2))
            meanspeed[acluster,:] = np.sqrt(meanu[acluster,:]**2.0 + meanv[acluster,:]**2.0)
            eddyspeed[acluster,:] = np.sqrt(eddyu[acluster,:]**2.0 + eddyv[acluster,:]**2.0)

            # compute decorrelation coefficient rho (Nc, Nt).  Note (Np,Nt,Ne) is dim of velx
            up = velx[clusterpart,:,:] - meanu[acluster,:][np.newaxis,:,np.newaxis]
            vp = vely[clusterpart,:,:] - meanv[acluster,:][np.newaxis,:,np.newaxis]
            # could generalize to full tensor if desired
            rhou[acluster,:] = np.mean(up[:,0,:][:,np.newaxis,:]*up, axis=(0,2))/np.mean(up[:,0,:]*up[:,0,:])
            rhov[acluster,:] = np.mean(vp[:,0,:][:,np.newaxis,:]*vp, axis=(0,2))/np.mean(vp[:,0,:]*vp[:,0,:])
            rho[acluster,:] = np.sqrt((rhou[acluster,:]**2.0 + rhov[acluster,:]**2.0)/2.0)
            assert np.abs(rho[acluster,0] - 1) < 1e-14, 'rho not normalized correctly on cluster %d with %f!'%(acluster, rho[acluster,0]-1.0)

            # compute mean dispersion over particle clusters (result are size Ncluster, Ntime, Nensembles)
            npart[acluster], clat[acluster,:,:], clon[acluster,:,:], \
                    drdr[acluster,:,:], dxdx[acluster,:,:], \
                    dxdy[acluster,:,:], dydy[acluster,:,:] = \
                    compute_dispersion_numexpr(llat[clusterpart,:,:], \
                            llon[clusterpart,:,:], rEarth, ensmean=0)

            # compute particle pair PDFs
            npartpairs[acluster], npairs[acluster], xmPDF[acluster,:,:,:], ymPDF[acluster,:,:,:], \
                    dxPDF[acluster,:,:,:], dyPDF[acluster,:,:,:], \
                    lPDF[acluster,:,:,:], lmPDF[acluster,:,:,:],  \
                    sePDF[acluster,:,:,:], uePDF[acluster,:,:,:], vePDF[acluster,:,:,:], \
                    smPDF[acluster,:,:,:], umPDF[acluster,:,:,:], vmPDF[acluster,:,:,:], \
                    KxxPDF[acluster,:,:,:], KxyPDF[acluster,:,:,:], KyyPDF[acluster,:,:,:], KrrPDF[acluster,:,:,:], \
                    KxxMean[acluster,:,:], KxyMean[acluster,:,:], KyyMean[acluster,:,:], KrrMean[acluster,:,:], thetaPDF[acluster,:,:] = \
                    particle_pair_pdfs(llat[clusterpart,:,:], llon[clusterpart,:,:], lvec, deltat, rEarth)

            print 'finished %d cluster of %d clusters'%(acluster, Nc)

            print 'compute diffs'
            K_rr = 0.5*np.diff(drdr[acluster,:,:],axis=0)/(24.*60.*60.*deltat)
            K_xx = 0.5*np.diff(dxdx[acluster,:,:],axis=0)/(24.*60.*60.*deltat)
            K_xy = 0.5*np.diff(dxdy[acluster,:,:],axis=0)/(24.*60.*60.*deltat)
            K_yy = 0.5*np.diff(dydy[acluster,:,:],axis=0)/(24.*60.*60.*deltat)

            # print and store values
            percKrr[acluster,:,:] = print_diff('K_rr', acluster, K_rr, KrrMean)
            check_pair_equality('K_rr', acluster, lmPDF, KrrPDF, KrrMean)
            percKxx[acluster,:,:] = print_diff('K_xx', acluster, K_xx, KxxMean)
            check_pair_equality('K_xx', acluster, lmPDF, KxxPDF, KxxMean)
            percKxy[acluster,:,:] = print_diff('K_xy', acluster, K_xy, KxyMean)
            check_pair_equality('K_xy', acluster, lmPDF, KxyPDF, KxyMean)
            percKyy[acluster,:,:] = print_diff('K_yy', acluster, K_yy, KyyMean)
            check_pair_equality('K_yy', acluster, lmPDF, KyyPDF, KyyMean)
        #}}}

        # compute weight factors w_b for K_tot = sum w_b * K_b, where b are PDF bins
        weight = lmPDF/np.nansum(lmPDF,axis=1)[:,np.newaxis,:,:]
        meanweight = np.mean(weight,axis=-1)

        # compute diffusivity from dispersion 1/2 * d/dt (variance)
        K_rr = 0.5*np.diff(drdr,axis=1)/(24.*60.*60.*deltat)
        K_xx = 0.5*np.diff(dxdx,axis=1)/(24.*60.*60.*deltat)
        K_xy = 0.5*np.diff(dxdy,axis=1)/(24.*60.*60.*deltat)
        K_yy = 0.5*np.diff(dydy,axis=1)/(24.*60.*60.*deltat)

        # compute max global error
        def check_global_error(name, clust, pair):
            print name + ' global error maximum = %.2e%%'%(np.nanmax(np.abs(perc_diff(clust,pair))))

        check_global_error('K_rr', K_rr, KrrMean)
        check_global_error('K_xx', K_xx, KxxMean)
        check_global_error('K_xy', K_xy, KxyMean)
        check_global_error('K_yy', K_yy, KyyMean)


        # need to average over the Ne ensembles (first dimension is Nt) #{{{
        drdrmean = np.mean(drdr,axis=-1)
        dxdxmean = np.mean(dxdx,axis=-1)
        dxdymean = np.mean(dxdy,axis=-1)
        dydymean = np.mean(dydy,axis=-1)

        meanK_rr = np.mean(K_rr,axis=-1)
        meanK_xx = np.mean(K_xx,axis=-1)
        meanK_xy = np.mean(K_xy,axis=-1)
        meanK_yy = np.mean(K_yy,axis=-1)

        meanclat = np.mean(clat,axis=-1)
        meanclon = np.mean(clon,axis=-1)

        meanlPDF   = np.nanmean(lPDF, axis=-1)
        meanthetaPDF = np.nanmean(thetaPDF, axis=-1)
        meanlmPDF  = np.nanmean(lmPDF, axis=-1)
        meanxmPDF  = np.nanmean(xmPDF,axis=-1)
        meanymPDF  = np.nanmean(ymPDF,axis=-1)
        meandxPDF  = np.nanmean(dxPDF,axis=-1)
        meandyPDF  = np.nanmean(dyPDF,axis=-1)
        meansePDF  = np.nanmean(sePDF,axis=-1)
        meanuePDF  = np.nanmean(uePDF,axis=-1)
        meanvePDF  = np.nanmean(vePDF,axis=-1)
        meansmPDF  = np.nanmean(smPDF,axis=-1)
        meanumPDF  = np.nanmean(umPDF,axis=-1)
        meanvmPDF  = np.nanmean(vmPDF,axis=-1)
        meanKxxPDF = np.nanmean(KxxPDF,axis=-1)
        meanKxyPDF = np.nanmean(KxyPDF,axis=-1)
        meanKyyPDF = np.nanmean(KyyPDF,axis=-1)
        meanKrrPDF = np.nanmean(KrrPDF,axis=-1)
        meanKxxMean = np.nanmean(KxxMean,axis=-1)
        meanKxyMean = np.nanmean(KxyMean,axis=-1)
        meanKyyMean = np.nanmean(KyyMean,axis=-1)
        meanKrrMean = np.nanmean(KrrMean,axis=-1)

        meanpercKrr=np.mean(percKrr,axis=-1)
        meanpercKxx=np.mean(percKxx,axis=-1)
        meanpercKxy=np.mean(percKxy,axis=-1)
        meanpercKyy=np.mean(percKyy,axis=-1)
        #}}}

        # save data and make sure to document how files were produced
        savenpz(case+'diff_layer%d'%(nlayer), x=x, y=y, ts=ts, te=ts+Nlen, deltat=deltat,
                percKrr_full=percKrr, percKxx_full=percKxx, percKxy_full=percKxy, percKyy_full=percKyy,
                percKrr=meanpercKrr, percKxx=meanpercKxx, percKxy=meanpercKxy, percKyy=meanpercKyy,
                meta_info='percentage (%) difference between particle pair and cluster diffusivities')

        savenpz(case+'pairs_layer%d'%(nlayer), npartinpairs=npartpairs, npairs=npairs,
                x=x, y=y, ts=ts, te=ts+Nlen, deltat=deltat, bins=lvec,
                xm_full=xmPDF, ym_full=ymPDF, dx_full=dxPDF, dy_full=dyPDF,
                lm_full=lmPDF, l_full=lPDF, lm=meanlmPDF, l=meanlPDF,
                theta_full = thetaPDF, theta=meanthetaPDF, weight_full = weight, weight = meanweight,
                xm = meanxmPDF, ym=meanymPDF, dx=meandxPDF, dy=meandyPDF,
                se_full=sePDF, ue_full=uePDF, ve_full=vePDF,
                se=meansePDF, ue=meanuePDF, ve=meanvePDF,
                sm_full=smPDF, um_full=umPDF, vm_full=vmPDF,
                sm=meansmPDF, um=meanumPDF, vm=meanvmPDF,
                Kxx_full=KxxPDF, Kxy_full=KxyPDF, Kyy_full=KyyPDF, Krr_full=KrrPDF,
                Kxx=meanKxxPDF, Kxy=meanKxyPDF, Kyy=meanKyyPDF, Krr=meanKrrPDF,
                meanKxx=meanKxxMean, meanKxy=meanKxyMean, meanKyy=meanKyyMean, meanKrr=meanKrrMean,
                meanKxx_full=KxxMean, meanKxy_full=KxyMean, meanKyy_full=KyyMean, meanKrr_full=KrrMean,
                meta_info='values are for computations with particle pairs')

        savenpz(case+'layer%d'%(nlayer), ts=ts, te=ts+Nlen, deltat=deltat,
                clon=clon, clat=clat, meanclon=meanclon, meanclat=meanclat,
                drdr_full=drdr, dxdx_full=dxdx, dxdy_full=dxdy, dydy_full=dydy,
                drdr=drdrmean, dxdx=dxdxmean, dxdy=dxdymean, dydy=dydymean,
                K_rr=K_rr, K_xx=K_xx, K_xy=K_xy, K_yy=K_yy,
                meanK_rr=meanK_rr, meanK_xx=meanK_xx, meanK_xy=meanK_xy, meanK_yy=meanK_yy,
                rho=rho, rhou=rhou, rhov=rhov,
                eddyspeed=eddyspeed, meanspeed=meanspeed,
                meanu=meanu, meanv=meanv, eddyu=eddyu, eddyv=eddyv,
                npartcluster=npart, x=x, y=y, meta_info='values are for computations with clusters')

        data  = np.load(case+'layer%d.npz'%(nlayer))
        print 'clusters', data.files

        data  = np.load(case+'pairs_layer%d.npz'%(nlayer))
        print 'pairs', data.files

        print 'finished %d layer with %d clusters' % (nlayer, Nc)

    return  case #}}}

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser() #{{{
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension or for analysis",
                      metavar="FILE")
    parser.add_option("-s", "--samplepoints", dest="samplepoints", help="npz file with x and y sample points", metavar="FILE")
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
    parser.add_option("--frac", dest="fraction", help="Fraction of particles in radius to use", metavar="FLOAT")
    parser.add_option("--lvec", dest="lvec", help="Bin edges for construction of PDFs", metavar="FLOAT")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename ('-f') is a required input... e.g., -f all_rlzn_particle_data.npz")
    if not options.samplepoints:
        parser.error("Sample points input filename ('-s') is a required input... e.g., -s coarse_sampling.npz")
    if not options.layer:
      options.layer = np.arange(5)
    else:
      try:
        options.layer = eval(options.layer)
      except:
        options.layer = np.array([int(options.layer)])
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
    if not options.fraction:
        options.fraction = 1.0
    else:
        options.fraction = float(options.fraction)
    if not options.lvec:
        #options.lvec = np.linspace(0,3000,201)
        options.lvec = np.array([0.0, 1.0, 4.0, 16.0, 64.0, 256.0])*30
    else:
        options.lvec = eval(options.lvec)
    #}}}
    print 'starting analysis'
    options.inputfilename = build_particle_file(options.inputfilename, options.samplepoints, options.lvec, float(options.radius), float(options.nlen),
            int(options.ts)-1, float(options.dt), int(options.nlayers), degrees=options.degrees, layerrange=options.layer, fraction=options.fraction)

