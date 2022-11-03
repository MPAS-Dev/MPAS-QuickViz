#!/usr/bin/env python
import numpy as np
import xray
import matplotlib as mpl
import socket
if socket.gethostname() == 'Shapiro' or socket.gethostname() == 'shapiro.lanl.gov':
  mpl.use('MacOSX')
else:
  mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy.spatial import cKDTree as KDTree
import glob
import re
daystos = 24.*60.*60.

def signed_log10(data):
  return np.sign(data)*np.log10(np.abs(data))

def main():
  # directories for comparisons
  fulldir = '/lustre/scratch1/turquoise/pwolfram/analysis/ZISO_analysis/full_10km/CaseBinningEnsembledComplete_r=100000.00_ts=0_deltat=1.00_Nlen=30_Nc=22800_frac=1.00/'
  highdir = '/lustre/scratch1/turquoise/pwolfram/analysis/ZISO_analysis/high_pass_10km/CaseBinningEnsembledComplete_r=100000.00_ts=0_deltat=1.00_Nlen=30_Nc=22800_frac=1.00/'
  lowdir  = '/lustre/scratch1/turquoise/pwolfram/analysis/ZISO_analysis/low_pass_10km/CaseBinningEnsembledComplete_r=100000.00_ts=0_deltat=1.00_Nlen=30_Nc=22800_frac=1.00/'

  ds = xray.open_dataset('buoyancySurfaceMeans.nc')
  buoyancySurface = ds.buoyancySurfaceValues.values

  files = glob.glob(fulldir + 'layer*.npz')
  layers = np.arange(len(files))
  print 'Calculating for layers : ', layers

  def aggregate_general_field(basedirlist, fielddata, mean=np.mean): #{{{
    """ example of fielddata = "data['mean_Krr']" """

    for alayer in layers:

      data = []
      for basedir in basedirlist:
        afile = basedir + 'layer'+str(alayer)+'.npz'
        data.append(np.load(afile))

      # assumes that all directories are on the same coordinate system
      y = data[0]['y']
      uniquey, ids = np.unique(y, return_index=True)
      t = np.arange(data[0]['ts'], data[0]['te']-1)*data[0]['deltat']

      field_all = np.nan*np.zeros((len(uniquey), len(t)))

      if alayer == layers[0]:
        fieldlist = np.nan*np.zeros((len(layers), len(uniquey), len(t)))

      field = eval(fielddata)

      for j, yids in enumerate(y[ids]):
        field_all[j,:] = mean(field[np.where(y==yids)[0],:],axis=0)

      fieldlist[alayer,:,:] = field_all

    return fieldlist #}}}

  # get TLs #{{{
  for alayer in layers:

    datafull = np.load(fulldir + 'layer'+str(alayer)+'.npz')
    datahigh = np.load(highdir + 'layer'+str(alayer)+'.npz')
    datalow  = np.load(lowdir  + 'layer'+str(alayer)+'.npz')

    # same grid for all of the directories
    y = datafull['y']
    uniquey, ids = np.unique(y, return_index=True)
    t = np.arange(datafull['ts'], datafull['te']-1)*datafull['deltat']

    # must compute separately because this is of a different dimension
    fullTL_all = np.nan*np.zeros((len(uniquey)))
    highTL_all = np.nan*np.zeros((len(uniquey)))
    lowTL_all  = np.nan*np.zeros((len(uniquey)))
    if alayer == layers[0]:
      fullTLlist = np.nan*np.zeros((len(layers), len(uniquey)))
      highTLlist = np.nan*np.zeros((len(layers), len(uniquey)))
      lowTLlist  = np.nan*np.zeros((len(layers), len(uniquey)))
    fullTL = 0.5*datafull['rhou'].shape[1]*(np.nanmean(datafull['rhou'],axis=1) + np.nanmean(datafull['rhov'],axis=1))*datafull['deltat']
    highTL = 0.5*datahigh['rhou'].shape[1]*(np.nanmean(datahigh['rhou'],axis=1) + np.nanmean(datahigh['rhov'],axis=1))*datahigh['deltat']
    lowTL  = 0.5* datalow['rhou'].shape[1]*(np.nanmean( datalow['rhou'],axis=1) + np.nanmean( datalow['rhov'],axis=1))* datalow['deltat']
    for j, yids in enumerate(y[ids]):
      fullTL_all[j] = np.mean(fullTL[np.where(y==yids)[0]],axis=0)
      highTL_all[j] = np.mean(highTL[np.where(y==yids)[0]],axis=0)
      lowTL_all[j]  = np.mean( lowTL[np.where(y==yids)[0]],axis=0)
    fullTLlist[alayer,:] = fullTL_all
    highTLlist[alayer,:] = highTL_all
    lowTLlist[alayer,:]  = lowTL_all
  #}}}

  # obviously non-pythonic code, but it does the job consisely (at cost of slowing code down)
  # 0 is full, 1 is high pass, 2 is low pass for
  basedirlist = [fulldir, highdir, lowdir]
  fullkrr        = "data[0]['meanK_rr']"
  fullkyy        = "data[0]['meanK_yy']"
  fullabskrr     = "data[0]['absmeanK_rr']"
  fullabskyy     = "data[0]['absmeanK_yy']"

  highkrr        = "data[1]['meanK_rr']"
  highkyy        = "data[1]['meanK_yy']"
  highabskrr     = "data[1]['absmeanK_rr']"
  highabskyy     = "data[1]['absmeanK_yy']"

  lowkrr        = "data[2]['meanK_rr']"
  lowkyy        = "data[2]['meanK_yy']"
  lowabskrr     = "data[2]['absmeanK_rr']"
  lowabskyy     = "data[2]['absmeanK_yy']"

  # interaction terms is full - high - low
  intkrr        = "(" + "data[0]['meanK_rr']    - data[1]['meanK_rr']     - data[2]['meanK_rr']"    + ")"
  intkyy        = "(" + "data[0]['meanK_yy']    - data[1]['meanK_yy']     - data[2]['meanK_yy']"    + ")"
  intabskrr     = "(" + "data[0]['absmeanK_rr'] - data[1]['absmeanK_rr']  - data[2]['absmeanK_rr']" + ")"
  intabskyy     = "(" + "data[0]['absmeanK_yy'] - data[1]['absmeanK_yy']  - data[2]['absmeanK_yy']" + ")"

  fullkrrlist        = aggregate_general_field(basedirlist, fullkrr)
  fullkyylist        = aggregate_general_field(basedirlist, fullkyy)
  fullabskrrlist     = aggregate_general_field(basedirlist, fullabskrr)
  fullabskyylist     = aggregate_general_field(basedirlist, fullabskyy)

  highkrrlist        = aggregate_general_field(basedirlist, highkrr)
  highkyylist        = aggregate_general_field(basedirlist, highkyy)
  highabskrrlist     = aggregate_general_field(basedirlist, highabskrr)
  highabskyylist     = aggregate_general_field(basedirlist, highabskyy)

  lowkrrlist        = aggregate_general_field(basedirlist, lowkrr)
  lowkyylist        = aggregate_general_field(basedirlist, lowkyy)
  lowabskrrlist     = aggregate_general_field(basedirlist, lowabskrr)
  lowabskyylist     = aggregate_general_field(basedirlist, lowabskyy)

  intkrrlist        = aggregate_general_field(basedirlist, intkrr)
  intkyylist        = aggregate_general_field(basedirlist, intkyy)
  intabskrrlist     = aggregate_general_field(basedirlist, intabskrr)
  intabskyylist     = aggregate_general_field(basedirlist, intabskyy)

  #vel = ds.buoyancySurfaceVelocityZonalMean.values
  depth = ds.buoyancySurfaceDepthMean.values
  y = uniquey*np.ones_like(depth)

  y=y.T/1000
  depth=depth.T

  def plot_time_varying_panels(vararray, contourrange, cbarticks, savename, cmapname='viridis'): #{{{
    fig = plt.figure(figsize=(5,15))
    for at,atime in enumerate(t):
      patches = []
      ax = plt.subplot(len(t),1,at+1)
      plt.contourf(y, depth, vararray[:,:,at].T, contourrange, cmap=cmapname)
      # plot mean buoyancy surface depths
      #for i,a in enumerate(ds.buoyancySurfaceDepthMean):
      #    plt.plot(ds.yCell/1000,a,'k-',alpha=0.15)
      if at > 0:
        ax.get_xaxis().set_visible(False)
      else:
        plt.gca().xaxis.set_ticks_position('top')
      if at+1 == len(t):
        cbar = plt.colorbar(orientation='horizontal', ticks=cbarticks)
      #plt.contour(y, depth, umlist[:,:,at].T, colors='k', levels=np.linspace(-1.0,1.0,2*12+1), lw=2)
      ax.text(0.05*np.max(y),0.2*np.min(depth), 'Time = %.1f'%(t[at]))
    plt.tight_layout(h_pad=0)
    plt.savefig(savename + '.png')
    return #}}}

  # plot time histories of diffusivities #{{{
  plot_time_varying_panels(signed_log10(fullkrrlist), np.linspace(0,6,2*12+1), [0,2,4,6],      'fullTimeVaryingKrr' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(fullkyylist), np.linspace(0,6,2*12+1), [0,2,4,6],      'fullTimeVaryingKyy' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(fullabskrrlist), np.linspace(0,6,2*12+1), [0,2,4,6],   'fullTimeVaryingAbsKrr' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(fullabskyylist), np.linspace(0,6,2*12+1), [0,2,4,6],   'fullTimeVaryingAbsKyy' , cmapname= 'magma' )

  plot_time_varying_panels(signed_log10(highkrrlist), np.linspace(0,6,2*12+1), [0,2,4,6],      'highTimeVaryingKrr' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(highkyylist), np.linspace(0,6,2*12+1), [0,2,4,6],      'highTimeVaryingKyy' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(highabskrrlist), np.linspace(0,6,2*12+1), [0,2,4,6],   'highTimeVaryingAbsKrr' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(highabskyylist), np.linspace(0,6,2*12+1), [0,2,4,6],   'highTimeVaryingAbsKyy' , cmapname= 'magma' )

  plot_time_varying_panels(signed_log10(lowkrrlist), np.linspace(0,6,2*12+1), [0,2,4,6],      'lowTimeVaryingKrr' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(lowkyylist), np.linspace(0,6,2*12+1), [0,2,4,6],      'lowTimeVaryingKyy' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(lowabskrrlist), np.linspace(0,6,2*12+1), [0,2,4,6],   'lowTimeVaryingAbsKrr' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(lowabskyylist), np.linspace(0,6,2*12+1), [0,2,4,6],   'lowTimeVaryingAbsKyy' , cmapname= 'magma' )

  plot_time_varying_panels(signed_log10(intkrrlist), np.linspace(0,6,2*12+1), [0,2,4,6],      'intTimeVaryingKrr' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(intkyylist), np.linspace(0,6,2*12+1), [0,2,4,6],      'intTimeVaryingKyy' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(intabskrrlist), np.linspace(0,6,2*12+1), [0,2,4,6],   'intTimeVaryingAbsKrr' , cmapname= 'magma' )
  plot_time_varying_panels(signed_log10(intabskyylist), np.linspace(0,6,2*12+1), [0,2,4,6],   'intTimeVaryingAbsKyy' , cmapname= 'magma' )
  #}}}

  plt.close('all')
  ##################################################
  ### specification of output time
  ##################################################
  for at in [10, 15, 20]:

  ##################################################
  ### YZ plots
  ##################################################
    def plot_yz_panel(vararray, contourrange, cbarticks, titlename, savename, cmapname='viridis'): #{{{
      plt.figure()
      plt.contourf(y, depth, vararray[:,:,at].T, contourrange, cmap=cmapname)
      plt.colorbar(ticks=cbarticks)
      cs = plt.contour(y, depth, vararray[:,:,at].T, colors='k', levels=contourrange, lw=2, alpha=0.3)
      plt.clabel(cs, inline=1, fontsize=6, alpha=0.3, fmt='%.2f')
      plt.xlabel('Y (km)')
      plt.ylabel('Z (m)')
      plt.xlim([0,2000])
      plt.ylim([-2500,0])
      plt.title(titlename + r" at t=%d days"%(int(t[at])))
      plt.savefig(savename + '_t=%d'%(int(t[at])) + '.png')
      return #}}}

    def plot_yz_panel_contour(vararray, vararray2, contourrange, contourrange2, cbarticks, titlename, savename, cmapname='viridis'): #{{{
      plt.figure()
      plt.contourf(y, depth, vararray[:,:,at].T, contourrange, cmap=cmapname)
      plt.colorbar(ticks=cbarticks)
      cs = plt.contour(y, depth, vararray[:,:,at].T, colors='k', levels=contourrange, lw=2, alpha=0.3)
      plt.clabel(cs, inline=1, fontsize=6, alpha=0.3, fmt='%.2f')
      cs = plt.contour(y, depth, vararray2[:,:,at].T, colors='k', levels=contourrange2, lw=2)
      plt.clabel(cs, inline=1, fontsize=6, fmt='%.2f')
      plt.xlabel('Y (km)')
      plt.ylabel('Z (m)')
      plt.xlim([0,2000])
      plt.ylim([-2500,0])
      plt.title(titlename + r" at t=%d days"%(int(t[at])))
      plt.savefig(savename + '_t=%d'%(int(t[at])) + '.png')
      return #}}}

  ##################################################
  ### YB plots
  ##################################################
    potentialDensityOffset = 1029.0
    buoysurf = buoyancySurface*np.ones_like(y) - potentialDensityOffset

    def plot_yb_panel(vararray, contourrange, cbarticks, titlename, savename, cmapname='viridis'): #{{{
      plt.figure()
      plt.contourf(y, buoysurf, vararray[:,:,at].T, contourrange, cmap=cmapname)
      plt.colorbar(ticks=cbarticks)
      cs = plt.contour(y, buoysurf, vararray[:,:,at].T, colors='k', levels=contourrange, lw=2, alpha=0.3)
      plt.clabel(cs, inline=1, fontsize=6, alpha=0.3, fmt='%.2f')
      plt.xlabel('Y (km)')
      plt.ylabel('Potential Density - 1029.0 (kg m$^{-3}$)')
      plt.gca().invert_yaxis()
      plt.xlim([0,2000])
      plt.title(titlename + r" at t=%d days"%(int(t[at])))
      plt.savefig(savename + '_t=%d'%(int(t[at])) + '.png')
      return #}}}

    def plot_yb_panel_contour(vararray, vararray2, contourrange, contourrange2, cbarticks, titlename, savename, cmapname='viridis'): #{{{
      plt.figure()
      plt.contourf(y, buoysurf, vararray[:,:,at].T, contourrange, cmap=cmapname)
      plt.colorbar(ticks=cbarticks)
      cs = plt.contour(y, buoysurf, vararray[:,:,at].T, colors='k', levels=contourrange, lw=2, alpha=0.3)
      plt.clabel(cs, inline=1, fontsize=6, alpha=0.3, fmt='%.2f')
      cs = plt.contour(y, buoysurf, vararray2[:,:,at].T, colors='k', levels=contourrange2, lw=2)
      plt.clabel(cs, inline=1, fontsize=6, fmt='%.2f')
      plt.xlabel('Y (km)')
      plt.ylabel('Potential Density - 1029.0 (kg m$^{-3}$)')
      plt.gca().invert_yaxis()
      plt.xlim([0,2000])
      plt.title(titlename + r" at t=%d days"%(int(t[at])))
      plt.savefig(savename + '_t=%d'%(int(t[at])) + '.png')
      return #}}}

  ##################################################
  ### Combination of YZ and YB plots plots
  ##################################################
    def make_panel_plots(vararray, contourrange, cbarticks, titlename, savename, cmapname='viridis'): #{{{
      plot_yz_panel(vararray, contourrange, cbarticks, titlename, savename + '_YZplot', cmapname)
      plot_yb_panel(vararray, contourrange, cbarticks, titlename, savename + '_YBplot', cmapname)
      plt.close('all')
      return #}}}

    def make_panel_contour_plots(vararray, vararray2, contourrange, contourrange2, cbarticks, titlename, savename, cmapname='viridis'): #{{{
      plot_yz_panel_contour(vararray, vararray2, contourrange, contourrange2, cbarticks, titlename, savename + '_YZplot', cmapname)
      plot_yb_panel_contour(vararray, vararray2, contourrange, contourrange2, cbarticks, titlename, savename + '_YBplot', cmapname)
      plt.close('all')
      return #}}}

    make_panel_plots(signed_log10(fullkrrlist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'fullDiffusivityRR', cmapname='magma')
    make_panel_plots(signed_log10(fullkyylist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{yy}$ (m$^2$ s$^{-1}$)", 'fullDiffusivityYY', cmapname='magma')

    make_panel_plots(signed_log10(fullabskrrlist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"Abs $\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'fullAbsDiffusivityRR', cmapname='magma')
    make_panel_plots(signed_log10(fullabskyylist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"Abs $\log_{10}\kappa_{yy}$ (m$^2$ s$^{-1}$)", 'fullAbsDiffusivityYY', cmapname='magma')

    make_panel_plots(signed_log10(highkrrlist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'highDiffusivityRR', cmapname='magma')
    make_panel_plots(signed_log10(highkyylist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{yy}$ (m$^2$ s$^{-1}$)", 'highDiffusivityYY', cmapname='magma')

    make_panel_plots(signed_log10(highabskrrlist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"Abs $\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'highAbsDiffusivityRR', cmapname='magma')
    make_panel_plots(signed_log10(highabskyylist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"Abs $\log_{10}\kappa_{yy}$ (m$^2$ s$^{-1}$)", 'highAbsDiffusivityYY', cmapname='magma')

    make_panel_plots(signed_log10(lowkrrlist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'lowDiffusivityRR', cmapname='magma')
    make_panel_plots(signed_log10(lowkyylist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{yy}$ (m$^2$ s$^{-1}$)", 'lowDiffusivityYY', cmapname='magma')

    make_panel_plots(signed_log10(lowabskrrlist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"Abs $\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'lowAbsDiffusivityRR', cmapname='magma')
    make_panel_plots(signed_log10(lowabskyylist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"Abs $\log_{10}\kappa_{yy}$ (m$^2$ s$^{-1}$)", 'lowAbsDiffusivityYY', cmapname='magma')

    make_panel_plots(signed_log10(intkrrlist), np.linspace(-6,6,4*6+1),     [-6,-4,-2,0,2,4,6], r"$\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'intDiffusivityRR', cmapname='Spectral')
    make_panel_plots(signed_log10(intkyylist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{yy}$ (m$^2$ s$^{-1}$)", 'intDiffusivityYY', cmapname='magma')

    make_panel_plots(signed_log10(intabskrrlist), np.linspace(-6,6,4*6+1),     [-6,-4,-2,0,2,4,6], r"Abs $\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'intAbsDiffusivityRR', cmapname='Spectral')
    make_panel_plots(signed_log10(intabskyylist), np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"Abs $\log_{10}\kappa_{yy}$ (m$^2$ s$^{-1}$)", 'intAbsDiffusivityYY', cmapname='magma')

if __name__ == "__main__":
  main()
