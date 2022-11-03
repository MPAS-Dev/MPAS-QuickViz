#!/usr/bin/env python
"""

    Compute particle pair statistic pdfs and decorrelation via cluster procedure (SIMPLIFIED).
    This decomposes diffusivity based on the lenghts inherent in the pairs which are used to
    build the histogram.

    Phillip Wolfram
    LANL
    07/14/2015

"""
## plot on mustang
#import matplotlib as mpl
#mpl.use('Agg')
# import libraries / packages
import numpy as np
from file_handling import check_folder, savenpz
from latlon_coordinate_transforms import signed_distances_numexpr as signed_distances, haversine_formula_numexpr as haversine_formula
from compute_diffusivity import compute_dispersion_numexpr, particle_pair_pdfs_simple

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
        degrees=False, fraction=1.0, pairselection='all'):  #{{{

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

    print 'pairselection = %s'%(pairselection)
    if pairselection == 'hull':
        case = 'CaseClusterRingPDFsSimple_r=%0.2f_ts=%d_deltat=%0.2f_Nlen=%d_Nc=%d_frac=%0.2f/'%(radius, ts, deltat, Nlen, Nc, fraction)
    elif pairselection == 'nearest':
        case = 'CaseClusterNearestPDFsSimple_r=%0.2f_ts=%d_deltat=%0.2f_Nlen=%d_Nc=%d_frac=%0.2f/'%(radius, ts, deltat, Nlen, Nc, fraction)
    elif pairselection == 'all':
        case = 'CaseClusterPDFsSimple_r=%0.2f_ts=%d_deltat=%0.2f_Nlen=%d_Nc=%d_frac=%0.2f/'%(radius, ts, deltat, Nlen, Nc, fraction)
    else:
        assert False, 'pairselection= %s is not valid'%(pairselection)

    print 'Case is ' + case
    check_folder(case)

    # starting center points of cluster
    for nlayer in layerrange:

        llon = np.swapaxes(conv(f_in[lonname][ts:ts+Nlen+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:]),0,1)
        llat = np.swapaxes(conv(f_in[latname][ts:ts+Nlen+1,nlayer*Npartlayers:(nlayer+1)*Npartlayers,:]),0,1)

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
        Nb = len(lvec)-1

        # n+1/2 terms
        KrrPDF = np.zeros((Nc,Nb,Nt))
        npairs = np.zeros((Nc, Nb, Nt))

        #}}}

        for acluster in np.arange(Nc): #{{{
            # get particles in the cluster
            clusterpart = particles[acluster,:]

            # subsample particles
            clustids = np.where(clusterpart)[0]
            Ncorig = np.sum(clusterpart)
            # remove non-desired random fraction from cluster (requires numpy 1.7)
            removeids = np.random.choice(clustids, np.floor(clustids.shape[0]*(1-fraction)), replace=False)
            clusterpart[removeids] = False
            # checks to make sure subsampling was correct
            Ncsample = np.sum(clusterpart)
            assert np.abs(Ncsample - np.round(fraction*Ncorig)) < 2, 'Error with subsampling'

            # compute particle pair PDFs
            npairs[acluster,:,:], KrrPDF[acluster,:,:] = particle_pair_pdfs_simple(llat[clusterpart,:,:], llon[clusterpart,:,:], lvec, deltat, rEarth, pairselection)

            print 'finished %d cluster of %d clusters'%(acluster, Nc)
        #}}}

        # compute weight factors w_b for K_tot = sum w_b * K_b, where b are PDF bins
        weight = npairs/np.nansum(npairs,axis=1)[:,np.newaxis,:]

        # weighted, mean Krr
        Krr = np.nansum(KrrPDF*weight,axis=1)

        # save data and make sure to document how files were produced
        savenpz(case+'layer%d'%(nlayer), x=x, y=y, ts=ts, te=ts+Nlen, deltat=deltat, radius=radius, fraction=fraction,
            npairs=npairs, KrrPDF=KrrPDF, weight=weight, Krr=Krr, bins=lvec*1000., meta_info='cluster diffusivity PDF')

        data  = np.load(case+'layer%d.npz'%(nlayer))
        print 'clusters', data.files

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
    parser.add_option("--pairselection", dest="pairselection", help="only use hull pairs on boundary of cluster", metavar="BOOL")

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
    if not options.pairselection:
        options.pairselection = 'all'
    #}}}
    print 'starting analysis'
    options.inputfilename = build_particle_file(options.inputfilename, options.samplepoints, options.lvec, float(options.radius), float(options.nlen),
            int(options.ts)-1, float(options.dt), int(options.nlayers), degrees=options.degrees, layerrange=options.layer, fraction=options.fraction, pairselection=options.pairselection)

