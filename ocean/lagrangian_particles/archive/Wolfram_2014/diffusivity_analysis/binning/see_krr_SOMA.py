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
  ds = xray.open_dataset('buoyancySurfaceMeans.nc')
  buoyancySurface = ds.buoyancySurfaceValues.values

  plt.figure()
  for i,a in enumerate(ds.buoyancySurfaceDepthMean):
    plt.plot(ds.yCell/1000,a,label=ds.buoyancySurfaceValues[i].values, lw=2)
    plt.xlabel('Meridional location (km)')
    plt.ylabel('Depth (m)')
    plt.legend(loc='best', prop={'size':10})
    plt.title('Mean buoyancy surface depths')
    plt.savefig('meandepths.png')

  files = glob.glob('layer*.npz')
  layers = sorted([int(re.findall('[0-9]+',afile)[0]) for afile in files])

  fig = plt.figure(figsize=(6,15))
  for alayer in layers:
    ax = plt.subplot(len(layers),1,alayer+1)

    afile = 'layer'+str(alayer)+'.npz'
    data = np.load(afile)
    y = data['y']
    uniquey, ids = np.unique(y, return_index=True)
    t = np.arange(data['ts'], data['te']-1)*data['deltat']

    krr_all =        np.nan*np.zeros((len(uniquey), len(t)))
    up_all =         np.nan*np.zeros((len(uniquey), len(t)))
    um_all =         np.nan*np.zeros((len(uniquey), len(t)))
    eddyU_all =      np.nan*np.zeros((len(uniquey), len(t)))
    eddyV_all =      np.nan*np.zeros((len(uniquey), len(t)))
    meanU_all =      np.nan*np.zeros((len(uniquey), len(t)))
    meanV_all =      np.nan*np.zeros((len(uniquey), len(t)))
    alpha_eddy_all = np.nan*np.zeros((len(uniquey), len(t)))
    alpha_mean_all = np.nan*np.zeros((len(uniquey), len(t)))
    rhomag_all =     np.nan*np.zeros((len(uniquey), len(t)))
    TLmag_all=       np.nan*np.zeros((len(uniquey), len(t)))
    sigmar_all =     np.nan*np.zeros((len(uniquey), len(t)))
    TL_all =         np.nan*np.zeros((len(uniquey)))
    AED_all   =      np.nan*np.zeros((len(uniquey), len(t)))
    AMD_all   =      np.nan*np.zeros((len(uniquey), len(t)))
    AMTL_all  =      np.nan*np.zeros((len(uniquey), len(t)))
    AETL_all  =      np.nan*np.zeros((len(uniquey), len(t)))
    ADTL_all  =      np.nan*np.zeros((len(uniquey), len(t)))
    DMDTL_all =      np.nan*np.zeros((len(uniquey), len(t)))
    DMMD_all  =      np.nan*np.zeros((len(uniquey), len(t)))
    DMED_all  =      np.nan*np.zeros((len(uniquey), len(t)))
    DMMTL_all =      np.nan*np.zeros((len(uniquey), len(t)))
    DMETL_all =      np.nan*np.zeros((len(uniquey), len(t)))
    LLM_all   =      np.nan*np.zeros((len(uniquey), len(t)))
    LLE_all   =      np.nan*np.zeros((len(uniquey), len(t)))
    LLKEV_all =      np.nan*np.zeros((len(uniquey), len(t)))
    LLKMV_all =      np.nan*np.zeros((len(uniquey), len(t)))
    LTKMV_all =      np.nan*np.zeros((len(uniquey), len(t)))
    LTKEV_all =      np.nan*np.zeros((len(uniquey), len(t)))
    LTDMV_all =      np.nan*np.zeros((len(uniquey), len(t)))
    LTDEV_all =      np.nan*np.zeros((len(uniquey), len(t)))
    LTDK_all  =      np.nan*np.zeros((len(uniquey), len(t)))

    if alayer == layers[0]:
      krrlist =        np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      eddyUlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      eddyVlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      meanUlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      meanVlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      umlist =         np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      uplist =         np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      alpha_eddylist = np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      alpha_meanlist = np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      rhomaglist =     np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      TLmaglist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      sigmarlist =     np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      TLlist =         np.nan*np.zeros((len(layers), len(uniquey)))
      AEDlist   =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      AMDlist   =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      AMTLlist  =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      AETLlist  =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      ADTLlist  =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      DMDTLlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      DMMDlist  =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      DMEDlist  =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      DMMTLlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      DMETLlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LLMlist   =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LLElist   =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LLKEVlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LLKMVlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LTKMVlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LTKEVlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LTDMVlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LTDEVlist =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))
      LTDKlist  =      np.nan*np.zeros((len(layers), len(uniquey), len(t)))

    krr = data['meanK_rr']
    um = data['meanspeed']
    up = data['eddyspeed']
    eddyU = data['eddyu']
    eddyV = data['eddyv']
    meanU = data['meanu']
    meanV = data['meanv']
    sigmar = np.sqrt(data['drdr'])[:,1:]
    TL = 0.5*data['rhou'].shape[1]*(np.nanmean(data['rhou'],axis=1) + np.nanmean(data['rhov'],axis=1))*data['deltat']
    rhomag = np.sqrt(data['rhou']**2 + data['rhov']**2)
    rhomag /= rhomag[:,0][:,np.newaxis]
    TLmag = np.cumsum(rhomag,axis=1)*data['deltat']
    alpha_eddy = krr/(up**2.0*TL[:,np.newaxis]*daystos)
    alpha_mean = krr/(um**2.0*TL[:,np.newaxis]*daystos)
    AED = krr/(up*sigmar)
    AMD = krr/(um*sigmar)
    AMTL = krr/(um**2.*TL[:,np.newaxis]*daystos)
    AETL = krr/(up**2.*TL[:,np.newaxis]*daystos)
    ADTL = krr/(sigmar**2./(TL[:,np.newaxis]*daystos))
    DMDTL = sigmar**2./(TL[:,np.newaxis]*daystos)
    DMMD = um*sigmar
    DMED = up*sigmar
    DMMTL = um**2.*TL[:,np.newaxis]*daystos
    DMETL = up**2.*TL[:,np.newaxis]*daystos
    LLM = um*TL[:,np.newaxis]*daystos
    LLE = up*TL[:,np.newaxis]*daystos
    LLKEV = krr/up
    LLKMV = krr/um
    LTKMV = (krr/um**2.)/(daystos)
    LTKEV = (krr/up**2.)/(daystos)
    LTDMV = (sigmar/um)/(daystos)
    LTDEV = (sigmar/up)/(daystos)
    LTDK = (sigmar**2/krr)/(daystos)

    for j, yids in enumerate(y[ids]):
      krr_all[j,:] =                np.mean(krr[np.where(y==yids)[0],:],axis=0)
      up_all[j,:] =                 np.mean(up[np.where(y==yids)[0],:],axis=0)
      um_all[j,:] =                 np.mean(um[np.where(y==yids)[0],:],axis=0)
      eddyU_all[j,:] =              np.mean(eddyU[np.where(y==yids)[0],:],axis=0)
      eddyV_all[j,:] =              np.mean(eddyV[np.where(y==yids)[0],:],axis=0)
      meanU_all[j,:] =              np.mean(meanU[np.where(y==yids)[0],:],axis=0)
      meanV_all[j,:] =              np.mean(meanV[np.where(y==yids)[0],:],axis=0)
      TL_all[j] =                   np.mean(TL[np.where(y==yids)[0]],axis=0)
      rhomag_all[j,:] =             np.nanmean(rhomag[np.where(y==yids)[0],:],axis=0)
      TLmag_all[j,:] =              np.nanmean(TLmag[np.where(y==yids)[0],:],axis=0)
      alpha_eddy_all[j,:] =         np.mean(alpha_eddy[np.where(y==yids)[0],:],axis=0)
      alpha_mean_all[j,:] =         np.mean(alpha_mean[np.where(y==yids)[0],:],axis=0)
      sigmar_all[j,:] =             np.mean(sigmar[np.where(y==yids)[0],:],axis=0)
      AED_all[j,:]   =              np.mean(AED[np.where(y==yids)[0],:],axis=0)
      AMD_all[j,:]   =              np.mean(AMD[np.where(y==yids)[0],:],axis=0)
      AMTL_all[j,:]  =              np.mean(AMTL[np.where(y==yids)[0],:],axis=0)
      AETL_all[j,:]  =              np.mean(AETL[np.where(y==yids)[0],:],axis=0)
      ADTL_all[j,:]  =              np.mean(ADTL[np.where(y==yids)[0],:],axis=0)
      DMDTL_all[j,:] =              np.mean(DMDTL[np.where(y==yids)[0],:],axis=0)
      DMMD_all[j,:]  =              np.mean(DMMD[np.where(y==yids)[0],:],axis=0)
      DMED_all[j,:]  =              np.mean(DMED[np.where(y==yids)[0],:],axis=0)
      DMMTL_all[j,:] =              np.mean(DMMTL[np.where(y==yids)[0],:],axis=0)
      DMETL_all[j,:] =              np.mean(DMETL[np.where(y==yids)[0],:],axis=0)
      LLM_all[j,:]   =              np.mean(LLM[np.where(y==yids)[0],:],axis=0)
      LLE_all[j,:]   =              np.mean(LLE[np.where(y==yids)[0],:],axis=0)
      LLKEV_all[j,:] =              np.mean(LLKEV[np.where(y==yids)[0],:],axis=0)
      LLKMV_all[j,:] =              np.mean(LLKMV[np.where(y==yids)[0],:],axis=0)
      LTKMV_all[j,:] =              np.mean(LTKMV[np.where(y==yids)[0],:],axis=0)
      LTKEV_all[j,:] =              np.mean(LTKEV[np.where(y==yids)[0],:],axis=0)
      LTDMV_all[j,:] =              np.mean(LTDMV[np.where(y==yids)[0],:],axis=0)
      LTDEV_all[j,:] =              np.mean(LTDEV[np.where(y==yids)[0],:],axis=0)
      LTDK_all[j,:]  =              np.mean(LTDK[np.where(y==yids)[0],:],axis=0)

    krrlist[alayer,:,:] =              krr_all
    uplist[alayer,:,:] =               up_all
    umlist[alayer,:,:] =               um_all
    eddyUlist[alayer,:,:] =            eddyU_all
    eddyVlist[alayer,:,:] =            eddyV_all
    meanUlist[alayer,:,:] =            meanU_all
    meanVlist[alayer,:,:] =            meanV_all
    rhomaglist[alayer,:,:] =           rhomag_all
    TLmaglist[alayer,:,:] =            TLmag_all
    alpha_meanlist[alayer,:,:] =       alpha_mean_all
    alpha_eddylist[alayer,:,:] =       alpha_eddy_all
    sigmarlist[alayer,:,:] =           sigmar_all
    TLlist[alayer,:] =                 TL_all
    AEDlist[alayer,:,:]   =            AED_all
    AMDlist[alayer,:,:]   =            AMD_all
    AMTLlist[alayer,:,:]  =            AMTL_all
    AETLlist[alayer,:,:]  =            AETL_all
    ADTLlist[alayer,:,:]  =            ADTL_all
    DMDTLlist[alayer,:,:] =            DMDTL_all
    DMMDlist[alayer,:,:]  =            DMMD_all
    DMEDlist[alayer,:,:]  =            DMED_all
    DMMTLlist[alayer,:,:] =            DMMTL_all
    DMETLlist[alayer,:,:] =            DMETL_all
    LLMlist[alayer,:,:]   =            LLM_all
    LLElist[alayer,:,:]   =            LLE_all
    LLKEVlist[alayer,:,:] =            LLKEV_all
    LLKMVlist[alayer,:,:] =            LLKMV_all
    LTKMVlist[alayer,:,:] =            LTKMV_all
    LTKEVlist[alayer,:,:] =            LTKEV_all
    LTDMVlist[alayer,:,:] =            LTDMV_all
    LTDEVlist[alayer,:,:] =            LTDEV_all
    LTDKlist[alayer,:,:]  =            LTDK_all

  # buoyancy surface plot #{{{
    T, Y = np.meshgrid(t,uniquey/1000.)
    plt.contourf(T, Y, signed_log10(krr_all), np.linspace(-6,6,2*12+1), cmap='RdYlBu_r')
    if alayer > 0:
      ax.get_xaxis().set_visible(False)
    else:
      plt.gca().xaxis.set_ticks_position('top')
    if alayer+1 == len(layers):
      plt.colorbar(orientation='horizontal')
    ax.text(0.05*len(t),0.8*np.max(uniquey/1000), buoyancySurface[alayer])
    ax.contour(T, Y, krr_all, [0],colors='w')

  plt.tight_layout(h_pad=0)
  plt.savefig('BuoyancySurfaceVaryingKrr.png')
  #}}}

  #vel = ds.buoyancySurfaceVelocityZonalMean.values
  depth = ds.buoyancySurfaceDepthMean.values
  y = uniquey*np.ones_like(depth)

  y=y.T/1000
  depth=depth.T
  #vel=vel.T

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

  plot_time_varying_panels(umlist, np.linspace(-0.4,0.4,4*4+1), [-0.4,-0.2,0,0.2,0.4], 'TimeVaryingZonalMeanVelocity', cmapname='Spectral')
  plot_time_varying_panels(signed_log10(krrlist), np.linspace(0,6,2*12+1), [0,2,4,6], 'TimeVaryingKrr', cmapname='magma')
  plot_time_varying_panels(rhomaglist, np.linspace(0,1.0,10+1), [0,0.2,0.4, 0.6, 0.8, 1.], 'TimeVaryingRho', cmapname='viridis')
  plot_time_varying_panels(TLmaglist, np.linspace(0,30,30+1), [0,5,10,15,20,25,30], 'TimeVaryingTL', cmapname='viridis')

  plt.close('all')
  ##################################################
  ### specification of output time
  ##################################################
  for at in [10, 15]:

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

    make_panel_plots(signed_log10(krrlist),                                         np.linspace(0,6,4*6+1),     [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'Diffusivity', cmapname='magma')
    make_panel_contour_plots(signed_log10(krrlist), uplist, np.linspace(0,6,4*6+1), np.linspace(0.0,0.4,2*4+1), [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'DiffusivityEddyVel', cmapname='magma')
    make_panel_contour_plots(signed_log10(krrlist), umlist, np.linspace(0,6,4*6+1), np.linspace(0.0,0.4,2*4+1), [0,1,2,3,4,5,6], r"$\log_{10}\kappa_{rr}$ (m$^2$ s$^{-1}$)", 'DiffusivityMeanVel', cmapname='magma')

    make_panel_plots(uplist, np.linspace(0,0.4,4*4+1), [0,0.1,0.2,0.3,0.4], r"$u'$ (m s$^{-1}$)", 'EddySpeed', cmapname='viridis')
    make_panel_plots(umlist, np.linspace(0,0.4,4*4+1), [0,0.1,0.2,0.3,0.4], r"$\overline{u}$ (m s$^{-1}$)", 'MeanSpeed', cmapname='viridis')

    make_panel_plots(eddyUlist, np.linspace(-0.4,0.4,4*4+1), [-0.4, -0.2, 0, 0.2, 0.4], r"$u' \cdot \hat{i}$ (m s$^{-1}$)", 'EddySpeedXcomp', cmapname='Spectral')
    make_panel_plots(eddyVlist, np.linspace(-0.4,0.4,4*4+1), [-0.4, -0.2, 0, 0.2, 0.4], r"$u' \cdot \hat{j}$ (m s$^{-1}$)", 'EddySpeedYcomp', cmapname='Spectral')

    make_panel_plots(meanUlist, np.linspace(-0.4,0.4,4*4+1), [-0.4, -0.2, 0, 0.2, 0.4], r"$\overline{u} \cdot \hat{i}$ (m s$^{-1}$)", 'MeanSpeedXcomp', cmapname='Spectral')
    make_panel_plots(meanVlist, np.linspace(-0.4,0.4,4*4+1), [-0.4, -0.2, 0, 0.2, 0.4], r"$\overline{u} \cdot \hat{j}$ (m s$^{-1}$)", 'MeanSpeedYcomp', cmapname='Spectral')

    make_panel_plots(rhomaglist, np.linspace(0,1.0,10+1), [0,0.2,0.4, 0.6, 0.8, 1.], r"$R(t)$ (-)", 'Rho', cmapname='viridis')

    make_panel_plots(np.ones_like(krrlist)*\
                             TLlist[:,:,np.newaxis],     np.linspace(0,20,20+1),  [0,5,10,15,20,40,60,80,100], r"$T_L$ (days)",                              'LagrangianTime',                  cmapname='viridis')
    make_panel_plots(LTKMVlist,                 np.linspace(0,100,10+1), [0,20,40,60,80,100], r"$\kappa_{rr}/\overline{u}$$^2$ (days)", 'LagrangianTimeKappaMeanVel',      cmapname='viridis')
    make_panel_plots(LTKEVlist,                 np.linspace(0,100,10+1), [0,20,40,60,80,100], r"$\kappa_{rr}/u'^2$ (days)",             'LagrangianTimeKappaEddyVel',      cmapname='viridis')
    make_panel_plots(LTDMVlist,                 np.linspace(0,100,10+1), [0,20,40,60,80,100], r"$\sigma_{rr}/\overline{u}$ (days)",     'LagrangianTimeDispersionMeanVel', cmapname='viridis')
    make_panel_plots(LTDEVlist,                 np.linspace(0,100,10+1), [0,20,40,60,80,100], r"$\sigma_{rr}/u'$ (days)",               'LagrangianTimeDispersionEddyVel', cmapname='viridis')
    make_panel_plots(signed_log10(LTDKlist),    np.linspace(0,6,4*6+1),  [0,1,2,3,4,5,6], r"$\log_{10}[\sigma_{rr}^2/\kappa_{rr}]$ (days)",    'LagrangianTimeDispersionKappa',   cmapname='viridis')
    #LTKMV = (krr/um**2.)/(daystos)
    #LTKEV = (krr/up**2.)/(daystos)
    #LTDMV = (sigmar/um)/(daystos)
    #LTDEV = (sigmar/up)/(daystos)
    #LTDK = (sigmar**2/krr)/(daystos)

    make_panel_plots(signed_log10(sigmarlist),       np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}(\sigma_{rr})$ (m)",              'LagrangianLengthDispersion', cmapname='viridis')
    make_panel_plots(signed_log10(LLKEVlist),        np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}(\kappa_{rr}/u')$ (m)",           'LagrangianLengthKappaEddyVel', cmapname='viridis')
    make_panel_plots(signed_log10(LLKMVlist),        np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}(\kappa_{rr}/\overline{u})$ (m)", 'LagrangianLengthKappaMeanVel', cmapname='viridis')
    make_panel_plots(signed_log10(LLElist),          np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}$[$u'T_L$] (m)",                  'LagrangianLengthEddy', cmapname='viridis')
    make_panel_plots(signed_log10(LLMlist),          np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}$[$\overline{u}T_L$] (m)",        'LagrangianLengthMean', cmapname='viridis')
    #LLKEV = krr/up
    #LLKMV = krr/um
    #LLE = up*TL[:,np.newaxis]*daystos
    #LLM = um*TL[:,np.newaxis]*daystos

    make_panel_plots(signed_log10(DMETLlist), np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}$[$u'^2T_L$] (m$^2$ s$^{-1}$)",             'DiffusivityModelEddyTL', cmapname='magma')
    make_panel_plots(signed_log10(DMMTLlist), np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}$[$\overline{u}$$^2T_L$] (m$^2$ s$^{-1}$)", 'DiffusivityModelMeanTL', cmapname='magma')
    make_panel_plots(signed_log10(DMEDlist),  np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}$[$\sigma_{rr} u'$] (m$^2$ s$^{-1}$)",              'DiffusivityModelEddyDispersion', cmapname='magma')
    make_panel_plots(signed_log10(DMMDlist),  np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}$[$\sigma_{rr} \overline{u}$] (m$^2$ s$^{-1}$)",    'DiffusivityModelMeanDispersion', cmapname='magma')
    make_panel_plots(signed_log10(DMDTLlist), np.linspace(0,6.0,4*6+1), [0,1,2,3,4,5,6], r"$\log_{10}$[$\sigma_{rr}^2 / T_L }$] (m$^2$ s$^{-1}$)",       'DiffusivityModelDispersionTL', cmapname='magma')
    #DMETL = up**2.*TL[:,np.newaxis]*daystos
    #DMMTL = um**2.*TL[:,np.newaxis]*daystos
    #DMED = up*sigmar
    #DMMD = um*sigmar
    #DMDTL = sigmar**2./(TL[:,np.newaxis]*daystos)

    make_panel_plots(signed_log10(alpha_eddylist), np.linspace(-3.0,3.0,2*3*2+1), [-3,-2,-1,0,1,2,3], r"$\log_{10}$[$\kappa_{rr}(u'^{2}T_L)^{-1}$] (m$^2$ s$^{-1}$)",           'alphaEddyTL', cmapname='Spectral')
    make_panel_plots(signed_log10(alpha_meanlist), np.linspace(-3.0,3.0,2*3*2+1), [-3,-2,-1,0,1,2,3], r"$\log_{10}$[$\kappa_{rr}(\overline{u}$$^2T_L)^{-1}$] (m$^2$ s$^{-1}$)", 'alphaMeanTL', cmapname='Spectral')
    make_panel_plots(signed_log10(ADTLlist), np.linspace(-3.0,3.0,2*3*2+1), [-3,-2,-1,0,1,2,3], r"$\log_{10}$[$\kappa_{rr}(\sigma_{rr}^2/T_L)^{-1}$]",           'alphaDispersionTL',   cmapname='Spectral')
    make_panel_plots(signed_log10(AETLlist),              np.linspace(-3.0,3.0,2*3*2+1), [-3,-2,-1,0,1,2,3], r"$\log_{10}$[$\kappa_{rr}(u'^{2}T_L)^{-1}$]",           'alphaEddyTL',         cmapname='Spectral')
    make_panel_plots(signed_log10(AMTLlist),              np.linspace(-3.0,3.0,2*3*2+1), [-3,-2,-1,0,1,2,3], r"$\log_{10}$[$\kappa_{rr}(\overline{u}$$^2T_L)^{-1}$]", 'alphaMeanTL',         cmapname='Spectral')
    make_panel_plots(signed_log10(AMDlist),       np.linspace(-3.0,3.0,2*3*2+1), [-3,-2,-1,0,1,2,3], r"$\log_{10}$[$\kappa_{rr}(\sigma_{rr} \overline{u}$$)^{-1}$]",  'alphaMeanDispersion', cmapname='Spectral')
    make_panel_plots(signed_log10(AEDlist),       np.linspace(-3.0,3.0,2*3*2+1), [-3,-2,-1,0,1,2,3], r"$\log_{10}$[$\kappa_{rr}(\sigma_{rr} u'$$)^{-1}$]",            'alphaEddyDispersion', cmapname='Spectral')
    #ADTL = krr/(sigmar**2./(TL[:,np.newaxis]*daystos))
    #AETL = krr/(up**2.*TL[:,np.newaxis]*daystos)
    #AMTL = krr/(um**2.*TL[:,np.newaxis]*daystos)
    #AMD = krr/(um*sigmar)
    #AED = krr/(up*sigmar)

if __name__ == "__main__":
  main()
