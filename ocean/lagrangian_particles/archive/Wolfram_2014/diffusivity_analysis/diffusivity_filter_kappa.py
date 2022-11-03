#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from GeneralTriangulation import GeneralTriangulation
from scipy.stats import normaltest

def filter_kappa(afile, Nlayers):

    adata = np.load(afile)
    kappa_rr = adata['kappa_rr']
    kappa_xx = adata['kappa_xx']
    kappa_xy = adata['kappa_xy']
    kappa_yy = adata['kappa_yy']
    x = adata['mux_ts']
    y = adata['muy_ts']
    kappa_rr_filt = np.zeros(kappa_rr.shape)
    kappa_xx_filt = np.zeros(kappa_rr.shape)
    kappa_xy_filt = np.zeros(kappa_rr.shape)
    kappa_yy_filt = np.zeros(kappa_rr.shape)
    npart = x.shape[0]/Nlayers
    for alayer in np.arange(Nlayers):
        xl = x[(alayer)*npart:(alayer+1)*npart]
        yl = y[(alayer)*npart:(alayer+1)*npart]
        triang = GeneralTriangulation(xl, yl)

        kappa_rrl = kappa_rr[(alayer)*npart:(alayer+1)*npart]
        kappa_xxl = kappa_xx[(alayer)*npart:(alayer+1)*npart]
        kappa_xyl = kappa_xy[(alayer)*npart:(alayer+1)*npart]
        kappa_yyl = kappa_yy[(alayer)*npart:(alayer+1)*npart]

        kappa_rr_filt[(alayer)*npart:(alayer+1)*npart] = triang.smooth_laplacian(kappa_rrl, ntimes=10, mean=np.nanmean)
        kappa_xx_filt[(alayer)*npart:(alayer+1)*npart] = triang.smooth_laplacian(kappa_xxl, ntimes=10, mean=np.nanmean)
        kappa_xy_filt[(alayer)*npart:(alayer+1)*npart] = triang.smooth_laplacian(kappa_xyl, ntimes=10, mean=np.nanmean)
        kappa_yy_filt[(alayer)*npart:(alayer+1)*npart] = triang.smooth_laplacian(kappa_yyl, ntimes=10, mean=np.nanmean)
        print 'finished layer %d' % (alayer)
    np.savez(afile.replace('kappa','kappa_filtered'), kappa_rr=kappa_rr_filt, kappa_xx=kappa_xx_filt, kappa_xy=kappa_xy_filt, kappa_yy=kappa_yy_filt)

    return 

if __name__ == "__main__":
    from optparse import OptionParser
    
    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="afile",
            help="name of file to fileter",
            metavar="FILE")
    parser.add_option("-l", "--nlayers", dest="nlayers",
            help="number of layers",
            metavar="FILE")

    options, args = parser.parse_args()

    if not options.nlayers:
        options.nlayers = 5
    if not options.afile:
        options.afile = 'particle_ensembles/ensemble029/kappa-5-1-2.npz'

    filter_kappa(options.afile, options.nlayers)

