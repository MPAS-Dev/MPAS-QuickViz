#!/usr/bin/env python
"""
Plots fitted effective diffusivity from cumulative concentrations.

Phillip J. Wolfram
10/13/2016
"""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import PiecewiseFit as PF
from PlotSettings import colormap
import glob

def compute_eff_diff(folder, Nb=51, Nk=21, nfiles=365, nstart=25, nend=-1):

    files = glob.glob(folder + 'conc_t*_days_tracer0000.npz')
    files.sort()

    if len(files) != nfiles:
        print 'Warning!!! There is not enough data to compute effective diffusivity in %s for folder %s'%(files, folder)
        return

    bins = np.linspace(0, 1, Nb)
    cbin = 0.5*(bins[1:] + bins[:-1])

    t = np.arange(len(files))

    Nt = np.load(files[0])['area'].shape[0]
    abins = np.linspace(0, 2, Nb)*1e12  #units of 10^6 km^2 or 10^12 m^2
    cfits = 0.5*(abins[1:] + abins[:-1])
    C = np.zeros((cfits.shape[0], len(files)))

    for ii, af in enumerate(files):
        df = np.load(af)
        conc = df['area']
        C[:,ii] = np.histogram(np.arange(conc.shape[0]), bins*Nt, weights=np.sort(conc))[0] /\
                  np.histogram(np.arange(conc.shape[0]), bins*Nt)[0]

    dcda, dcdt = np.gradient(C)
    partda = 2. * 1e6 * 1e3**2 / dcda.shape[0]
    partdt = np.gradient(t) * 86400.
    dcdt /= partdt
    dcda /= partda
    d2cda2 = np.gradient(dcda)[0] / partda

    fit = PF.PiecewisePolyFit(cfits, Nk)
    fit.set_knots(cfits[0], 0)
    fit.set_knots(cfits[-1], 0)
    K = fit.func_eval(cfits)
    dKda = fit.deriv_eval(cfits)
    #b = np.mean(dcdt[:,25:],axis=1)
    A = (dKda*np.mean(dcda[:,nstart:nend],axis=1)[:,np.newaxis] +
         K*(np.gradient(np.mean(dcda[:,nstart:nend],axis=1))/partda)[:,np.newaxis])
    b = np.mean(dcdt[:,nstart:nend],axis=1)
    fit.solve(A, b, bounds=(0, np.inf))

    Kefffitmean = fit(cfits)/ 1e6**2
    print Kefffitmean.max()
    Ax = cfits[:,np.newaxis]*np.ones_like(C)/1e12  #convert back to 10^12 m^2 or 10^6 km^2 for plotting simplicity
    np.savez(folder + '/fitted_eff_diff_after25day.npz', Ax=Ax[:,0], Kefffitmean=Kefffitmean)

def test_plot_eff_diff(folder):
    mpl.pylab.rcParams['figure.figsize'] = (18.0, 10.0) # Large figures
    axis_font = {'fontname':'Arial', 'size':'35'}
    title_font = {'fontname':'Arial', 'size':'50', 'color':'black', 'weight':'normal'}
    label_font = {'fontname':'Arial', 'size':'20'}
    mpl.rc('xtick', labelsize=30)
    mpl.rc('ytick', labelsize=30)

    # resolution for figures
    mpl.pylab.rcParams['savefig.dpi'] = 150

    # def load Kefffitmean
    df = np.load(folder + '/fitted_eff_diff_after25day.npz')
    Ax = df['Ax'][:]
    Kefffitmean = df['Kefffitmean'][:]

    plt.plot(Ax, Kefffitmean, 'k-', lw=10, label='$K_\mathrm{eff}^*$ time-mean fit',zorder=998)

    def plot_particle_stats(days=25):
        alayer = 8 # based on input sim
        dsf = np.load('/Users/pwolfram/Documents/ReportsPresentations/ZISOEffDiffCalc/full_dataset/FullDiffusivityYY_10yrs_%sdaysdata.npz'%(days))
        plt.plot(dsf['uniquey']/1e6, 10**dsf['vararray'][alayer,:], 'g-', lw=10,  label='$\kappa_{yy}$ at %s days'%days,zorder=999)

    plot_particle_stats(days=25)

    plt.legend(loc='upper left', ncol=2, fontsize=int(axis_font['size'])-10, frameon=False)
    plt.ylabel('Diffusivity (m$^2$s$^{-1}$)', **axis_font)
    plt.xlabel('$A$ ($10^6$ km$^2$)', **axis_font)
    plt.xlim(0,2)
    plt.ylim(0,1e4)
    plt.savefig(folder + '/method_comparison.png', bbox_inches='tight', dpi=300)


if __name__== '__main__':
    import sys
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    else:
        folder='/Users/pwolfram/Documents/ReportsPresentations/ZISOEffDiffCalc/' + \
               'areas_concgrad_tester_remap_betterconcgrad_rearrange/'
    print 'testing fitted effective diffusivity for folder %s'%folder

    compute_eff_diff(folder)
    #test_plot_eff_diff(folder)

