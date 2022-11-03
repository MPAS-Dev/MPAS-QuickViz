#!/usr/bin/env python

import matplotlib as mpl
import socket
if socket.gethostname() == 'Shapiro':
    mpl.use('MacOSX')
else:
    mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from GeneralTriangulation import GeneralTriangulation
from scipy.stats import normaltest
import PlotSettings
import netCDF4
from comp_funcs import SignedLog10
from matplotlib import rcParams
# necessary to get ipython notebook-like behavior
rcParams['figure.figsize'] = (6.0,4.0)
rcParams['figure.facecolor'] = 'white'
rcParams['figure.edgecolor'] = 'white'
rcParams['font.size'] = 10
rcParams['savefig.dpi'] = 72
rcParams['figure.subplot.bottom'] = .125
rcParams['savefig.dpi'] = 150
plt.close('all')
plotme= True

def save_fig(name,dpi=300, *args, **kwargs):
        print 'saving ', name,
        plt.savefig(name,dpi=dpi, bbox_inches='tight')
        print 'done.'

def make_plots(directory, res):
    """
    # inputs example
    directory = '/Users/pwolfram/Documents/ReportsPresentations/SOMADiffusivityPaper/ProductionRuns/Revisions/4km/CaseBinningEnsembledComplete_r=50000.00_ts=0_deltat=2.00_Nlen=16_Nc=5796_frac=1.00/'
    res = r"4\,km"
    """

    Nlayers = 11
    buoySurf = ['1028.5', '1028.65', '1028.8', '1028.95', '1029.1', '1029.25',
    '1029.4', '1029.55', '1029.7', '1029.85', '1030.0']
    nfilt = 0
    escale2 = 1e9
    es = 1e-4
    excludelayers = [1,2,4,5,6,8,9]
    excludelayers.reverse()
    xlim = np.linspace(0,1000000,5)
    ylim = np.linspace(0,2000000,5)

    # build the grid
    Ndim = 10
    x,y = mpl.pylab.meshgrid(np.linspace(0,1000000.,Ndim),np.linspace(0,2000000.,2*Ndim))
    trisimp = GeneralTriangulation(x.flatten(),y.flatten())

    reload(PlotSettings)
    data = np.load(directory + 'layer0.npz','r')
    t = (data['ts'] + np.arange(data['drdr'].shape[1]-1))*data['deltat']
    for at in np.arange(t.shape[0]):
        kappa_xx = []
        kappa_xy = []
        kappa_yy = []
        kappa_rr = []
        for alayer in np.arange(Nlayers):
            data = np.load(directory + 'layer%d.npz'%(alayer))
            #xl = np.degrees(data['x'])
            #yl = np.degrees(data['y'])
            xl = data['x']
            yl = data['y']
            grid = GeneralTriangulation(xl,yl)
            kappa_rr.append(GeneralTriangulation(xl, yl, scalar=grid.smooth_laplacian(data['meanK_rr'][:, at]),ntimes=nfilt))
            # also perform filtering as before
            kappa_xx.append(GeneralTriangulation(xl, yl, scalar=grid.smooth_laplacian(data['meanK_xx'][:, at], ntimes=nfilt)) >> trisimp)
            kappa_xy.append(GeneralTriangulation(xl, yl, scalar=grid.smooth_laplacian(data['meanK_xy'][:, at], ntimes=nfilt)) >> trisimp)
            kappa_yy.append(GeneralTriangulation(xl, yl, scalar=grid.smooth_laplacian(data['meanK_yy'][:, at], ntimes=nfilt)) >> trisimp)

        bsl = buoySurf[:]
        for aexcludelayer in excludelayers:
            kappa_xx.pop(aexcludelayer)
            kappa_xy.pop(aexcludelayer)
            kappa_yy.pop(aexcludelayer)
            kappa_rr.pop(aexcludelayer)
            bsl.pop(aexcludelayer)
        plt.close('all')
        PlotSettings.make_2x2_subplot_tensor(kappa_xx, kappa_xy, kappa_yy, bsl, 'X', 'Y',
                      r"$\kappa_{\mathrm{C}}\,(10^{"+"%d" %(-np.log10(es))+r"} \mathrm{m}^2 \mathrm{s}^{-1})\,\Delta x=\mathrm{"+res+"}\,t=" + "%.0f"%(t[at]) +r"\,\mathrm{days}$", escale=es*escale2,escale2=escale2, xlim=xlim,ylim=ylim, plothull=False)
        save_fig(directory + 'tensor_vertical_struct_t=%02d.png'%(t[at]))
        plt.close('all')
        PlotSettings.make_2x2_subplot(kappa_rr, bsl, 'X', 'Y', r"$\log_{10}\kappa_{\mathrm{C}}\,(\mathrm{m}^2 \mathrm{s}^{-1})\,\Delta x=\mathrm{"+res+r"}\,t=" + "%.0f"%(t[at]) +r"\,\mathrm{days}$", func=SignedLog10.convert, clim=np.linspace(2.0,6.0,5), contour=True,xlim=xlim,ylim=ylim)
        save_fig(directory + 'scalar_vertical_struct_t=%02d.png'%(t[at]))
        plt.close('all')
        PlotSettings.make_2x2_subplot(kappa_xx, bsl, 'X', 'Y', r"$\log_{10}\kappa_{\mathrm{xx}}\,(\mathrm{m}^2 \mathrm{s}^{-1})\,\Delta x=\mathrm{"+res+"}\,t=" + "%.0f"%(t[at]) +r"\,\mathrm{days}$", func=SignedLog10.convert, clim=np.linspace(2.0,6.0,5), contour=False,xlim=xlim,ylim=ylim)
        save_fig(directory + 'K_xx_vertical_struct_t=%02d.png'%(t[at]))
        PlotSettings.make_2x2_subplot(kappa_xy, bsl, 'X', 'Y', r"$\log_{10}\kappa_{\mathrm{xy}}\,(\mathrm{m}^2 \mathrm{s}^{-1})\,\Delta x=\mathrm{"+res+"}\,t=" + "%.0f"%(t[at]) +r"\,\mathrm{days}$", func=SignedLog10.convert, clim=np.linspace(-6.0,6.0,7), contour=False, cmapin='seismic',xlim=xlim,ylim=ylim)
        save_fig(directory + 'K_xy_vertical_struct_t=%02d.png'%(t[at]))
        PlotSettings.make_2x2_subplot(kappa_yy, bsl, 'X', 'Y', r"$\kappa_{\mathrm{yy}}\,(\mathrm{m}^2 \mathrm{s}^{-1})\,\Delta x=\mathrm{"+res+"}\,t=" + "%.0f"%(t[at]) +r"\,\mathrm{days}$", func=SignedLog10.convert, clim=np.linspace(2.0,6.0,5), contour=False,xlim=xlim,ylim=ylim)
        save_fig(directory + 'K_yy_vertical_struct_t=%02d.png'%(t[at]))
        plt.close('all')

    for at in np.arange(t.shape[0]):
        eddyspeed = []
        meanspeed = []
        rho = []
        for alayer in np.arange(Nlayers):
            data = np.load(directory + 'layer%d.npz'%(alayer))
            #xl = np.degrees(data['x'])
            #yl = np.degrees(data['y'])
            #xl = data['meanclon'][:,at]
            #yl = data['meanclat'][:,at]
            xl = data['x']
            yl = data['y']
            grid = GeneralTriangulation(xl,yl)
            eddyspeed.append(GeneralTriangulation(xl, yl, scalar=grid.smooth_laplacian(data['eddyspeed'][:,at]),ntimes=nfilt))
            meanspeed.append(GeneralTriangulation(xl, yl, scalar=grid.smooth_laplacian(data['meanspeed'][:,at]),ntimes=nfilt))
            rho.append(GeneralTriangulation(xl, yl, scalar=grid.smooth_laplacian(data['rho'][:,at]),ntimes=nfilt))


        bsl = buoySurf[:]
        for aexcludelayer in excludelayers:
            eddyspeed.pop(aexcludelayer)
            meanspeed.pop(aexcludelayer)
            bsl.pop(aexcludelayer)
        plt.close('all')
        PlotSettings.make_2x2_subplot(eddyspeed, bsl, 'X', 'Y', r"$u'\,(\mathrm{m} \mathrm{s}^{-1})\,\Delta x=\mathrm{"+res+r"}\,t=" + "%.0f"%(t[at]) +r"\,\mathrm{days}$", clim=np.linspace(0.0,0.75,4), contour=True, xlim=xlim, ylim=ylim)
        save_fig(directory + 'eddyspeed_vertical_struct_t=%02d.png'%(t[at]))
        plt.close('all')
        PlotSettings.make_2x2_subplot(meanspeed, bsl, 'X', 'Y', r"$\overline{u}\,(\mathrm{m} \mathrm{s}^{-1})\,\Delta x=\mathrm{"+res+r"}\,t=" + "%.0f"%(t[at]) +r"\,\mathrm{days}$", clim=np.linspace(0.0,1.0,5), contour=True, xlim=xlim, ylim=ylim)
        save_fig(directory + 'meanspeed_vertical_struct_t=%02d.png'%(t[at]))
        plt.close('all')
        PlotSettings.make_2x2_subplot(rho, bsl, 'X', 'Y', r"$\rho\,(-)\,\Delta x=\mathrm{"+res+r"}\,t=" + "%.0f"%(t[at]) +r"\,\mathrm{days}$", clim=np.linspace(0.0,1.0,5), contour=True, xlim=xlim, ylim=ylim)
        save_fig(directory + 'rho_vertical_struct_t=%02d.png'%(t[at]))
        plt.close('all')

    data = np.load(directory + 'layer0.npz','r')
    t = (data['ts'] + np.arange(data['eddyspeed'].shape[1]))*data['deltat']
    for alayer in np.arange(Nlayers):
        data = np.load(directory + 'layer%d.npz'%(alayer))
        rho = data['rho']
        plt.plot(t,np.mean(rho,axis=0),label='%s'%(buoySurf[alayer]))
    plt.ylim(0,1.0)
    plt.xlabel('Time (days)')
    plt.ylabel(r"$\rho$")
    plt.title('Decorrelation for each layer')
    #plt.legend(loc='best')
    save_fig(directory + 'rho_domain_mean_t=%02d.png'%(t[at]))

    # compute decorrelation time
    data = np.load(directory + 'layer0.npz','r')
    t = (data['ts'] + np.arange(data['eddyspeed'].shape[1]))*data['deltat']
    TL = []
    for alayer in np.arange(Nlayers):
        data = np.load(directory + 'layer%d.npz'%(alayer))
        xl = data['x']
        yl = data['y']
        grid = GeneralTriangulation(xl,yl)
        tau = 0.5*(np.sum(data['rhou'],axis=1) + np.sum(data['rhov'],axis=1))*data['deltat']
        TL.append(GeneralTriangulation(xl, yl, scalar=grid.smooth_laplacian(tau),ntimes=nfilt))

    bsl = buoySurf[:]
    for aexcludelayer in excludelayers:
        TL.pop(aexcludelayer)
        bsl.pop(aexcludelayer)
    plt.close('all')
    PlotSettings.make_2x2_subplot(TL, bsl, 'X', 'Y', r"$T_L$ $\Delta x=\mathrm{"+res+r"}$", clim=np.linspace(0.0,20.0,5), contour=True, xlim=xlim, ylim=ylim)
    save_fig(directory + 'TL.png')
    plt.close('all')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(usage="-d <directory> -r '4\,km'", description='Plot diffusivity results from diffusivity_binning_complete.py script')
    parser.add_argument('-d', '--directory', action='store', metavar='DIRECTORY', type=str, help='directory containing processed layer*.npz files')
    parser.add_argument('-r', '--res', action='store', metavar='STRING', type=str, help="string with name for case, e.g., '4\,km'")
    args = vars(parser.parse_args())

    make_plots(args['directory'], args['res'])


