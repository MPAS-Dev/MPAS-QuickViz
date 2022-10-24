#!/usr/bin/env python
"""
Plot vertical sections (annual and seasonal climatologies),
given a certain transect mask
"""
############################## transects. Pick one:
transectNames = ['OSNAP section East']; include_colorbar = True
#transectNames = ['OSNAP section West']; include_colorbar = False

############################## months or seasons
seasonList = ['JFM', 'JAS', 'ANN']

############################## model files, run dirs, from cori for NARRM
### choose this one for the NARRM low resolution:
rr = '/global/cscratch1/sd/katsmith/archive/E3SMv2/'
simName = ['v2.LR.historical']
simShortName= ['v2.LR.historical']
meshfile = ['input_files/v2.LR.historical_0301.mpaso.rst.1980-01-01_00000.nc']
subdir = ['/0301/post/analysis/mpas_analysis/ts_1980-2014_climo_1980-2014/clim/mpas/avg/unmasked_EC30to60E2r2/']
maskfile = ['input_files/masks_v2.LR.historical.nc']

### choose this one for the NARRM historical:
rr = '/global/cscratch1/sd/katsmith/archive/E3SMv2/'
simName = ['v2.NARRM.historical']
simShortName= ['v2.NARRM.historical']
meshfile = ['input_files/v2.NARRM.historical_0301.mpaso.rst.1980-01-01_00000.nc']
subdir = ['/0301/post/analysis/mpas_analysis/ts_1980-2014_climo_1980-2014/clim/mpas/avg/unmasked_WC14to60E2r3/']
maskfile = ['input_files/masks_v2.NARRM.nc']

### choose this one for the NARRM historical:
rr='/global/cscratch1/sd/ethomas/MPAS_analysis/Q4_metric/'
simName = ['Interface_piControl_151-200']
simShortName= ['v2.NARRM.PI']
meshfile = ['input_files/v2.NARRM.historical_0301.mpaso.rst.1980-01-01_00000.nc']
subdir = ['/clim/mpas/avg/unmasked_WC14to60E2r3/']
maskfile = ['input_files/masks_v2.NARRM.nc']

############################## contours
sigma0contours = [24.0, 25.0, 26.0, 27.0, 27.2, 27.4, 27.6, 27.7, 27.75, 27.8, 27.82,  27.84, 27.86, 27.87, 27.88, 27.9, 27.95, 28.0, 28.05]

#from __future__ import absolute_import, division, print_function, \
#    unicode_literals
import os
import glob
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cols
from matplotlib.pyplot import cm
from matplotlib.colors import BoundaryNorm
import cmocean
import xarray as xr
from netCDF4 import Dataset
from mpas_analysis.shared.io.utility import decode_strings
import gsw

earthRadius = 6367.44

# Figure details
figdir = './verticalSections/'
if not os.path.isdir(figdir):
    os.makedirs(figdir)
if transectNames == 'OSNAP section East':
    figsize = (14, 6) # OSNAP East only
elif transectNames == 'OSNAP section West':
    figsize = (8, 6) # OSNAP West only w colorbar
    #figsize = (6, 6) # OSNAP West only no colorbar
else:
    figsize = (14, 6)

figdpi = 300
colorIndices0 = [0, 10, 28, 57, 85, 113, 125, 142, 155, 170, 198, 227, 242, 255]
#clevelsT = [-2.0, -1.8, -1.5, -1.0, -0.5, 0.0, 0.5, 2.0, 4.0, 6.0, 8.0, 10., 12.]
#clevelsS = [30.0, 31.0, 32.0, 33.0, 33.5, 34.0, 34.5, 34.8, 34.85, 34.9, 34.95, 35.0, 35.5]
# Better for OSNAP:
clevelsT = [-1.0, -0.5, 0.0, 0.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 8., 10., 12.]
#clevelsT = np.linspace(2.0, 6.0, 13)
clevelsS = [31.0, 33.0, 33.5, 33.8, 34.2, 34.6, 34.8, 34.85, 34.9, 34.95, 35.0, 35.2, 35.5]
#clevelsS = np.linspace(34.7, 35.2, 13)
#clevelsV = [-0.25, -0.2, -0.15, -0.1, -0.02, 0.0, 0.02, 0.1, 0.2, 0.3, 0.5]
clevelsV = np.linspace(-0.15,0.15, 31)
colormapT = plt.get_cmap('RdBu_r')
colormapS = cmocean.cm.haline
colormapV = plt.get_cmap('RdBu_r')
#colormapV = cmocean.cm.balance
#
underColor = colormapT(colorIndices0[0])
overColor = colormapT(colorIndices0[-1])
if len(clevelsT) + 1 == len(colorIndices0):
    # we have 2 extra values for the under/over so make the colormap
    # without these values
    colorIndices = colorIndices0[1:-1]
elif len(clevelsT) - 1 != len(colorIndices0):
    # indices list must be either one element shorter
    # or one element longer than colorbarLevels list
    raise ValueError('length mismatch between indices and '
                     'T colorbarLevels')
colormapT = cols.ListedColormap(colormapT(colorIndices))
colormapT.set_under(underColor)
colormapT.set_over(overColor)
underColor = colormapS(colorIndices0[0])
overColor = colormapS(colorIndices0[-1])
if len(clevelsS) + 1 == len(colorIndices0):
    # we have 2 extra values for the under/over so make the colormap
    # without these values
    colorIndices = colorIndices0[1:-1]
elif len(clevelsS) - 1 != len(colorIndices0):
    # indices list must be either one element shorter
    # or one element longer than colorbarLevels list
    raise ValueError('length mismatch between indices and '
                     'S colorbarLevels')
colormapS = cols.ListedColormap(colormapS(colorIndices))
colormapS.set_under(underColor)
colormapS.set_over(overColor)
#colormapV = cols.ListedColormap(colormapV(colorIndices))
#
cnormT = mpl.colors.BoundaryNorm(clevelsT, colormapT.N)
cnormS = mpl.colors.BoundaryNorm(clevelsS, colormapS.N)
cnormV = mpl.colors.BoundaryNorm(clevelsV, colormapV.N)

#sigma2contours = [35, 36, 36.5, 36.8, 37, 37.1, 37.2, 37.25, 37.44, 37.52, 37.6]
sigma2contours = None
#sigma0contours = np.arange(26.0, 28.0, 0.2) # Good for OSNAP, but not for all arcticSections
#sigma0contours = None

# Load in MPAS mesh and transect mask file
mesh = xr.open_dataset(meshfile[0])
mask = xr.open_dataset(maskfile[0])

allTransects = decode_strings(mask.transectNames)
if transectNames[0]=='all' or transectNames[0]=='StandardTransportSectionsRegionsGroup':
    transectNames = allTransects

# Get depth
z = mesh.refBottomDepth.values
nVertLevels = len(z)

refZMid = np.zeros(nVertLevels)
refZMid[0] =  0.5*z[0]
for k in range(1,nVertLevels):
    refZMid[k] =  0.5*(z[k-1] + z[k])

nTransects = len(transectNames)
maxCells = mask.dims['maxCellsInTransect']
# 3 is OSNAP East, 4 is OSNAP West
for iTransect in range(nTransects):
    # Identify transect
    transectName = transectNames[iTransect]
    transectIndex = allTransects.index(transectName)
    print('transect: ', transectName)

    # Choose mask for this particular transect
    transectmask = mask.isel(nTransects=transectIndex).squeeze()
    # Get a list of cells for this transect
    transectCells = transectmask.transectCellGlobalIDs.values
    transectCells = transectCells[np.where(transectCells > 0)[0]]
    ntransectCells = len(transectCells)

    # Get coordinates of each cell center and compute approximate spherical distance
    lonCells = mesh.lonCell.sel(nCells=transectCells-1).values
    latCells = mesh.latCell.sel(nCells=transectCells-1).values
    lonCells[lonCells>np.pi] = lonCells[lonCells>np.pi] - 2*np.pi

    bottomDepth = mesh.variables['bottomDepth'][transectCells-1]

    dist = [0]
    for iCell in range(1, ntransectCells):
        dx = (lonCells[iCell]-lonCells[0]) * np.cos(0.5*(latCells[iCell]+latCells[0]))
        dy = latCells[iCell]-latCells[0]
        norm = np.sqrt(dx**2 + dy**2)
        length = max(earthRadius * norm, dist[iCell-1])
        dist.append(length)
    [x, y] = np.meshgrid(dist, z)
    x = x.T

# compute weights for velocity perpendicular to the section:
    xWt = np.zeros(ntransectCells)
    yWt = np.zeros(ntransectCells)
    for iCell in range(1, ntransectCells-1):
        dx = (lonCells[iCell+1]-lonCells[iCell-1]) * np.cos(0.5*(latCells[iCell+1]+latCells[iCell-1]))
        dy = latCells[iCell+1]-latCells[iCell-1]
        norm = np.sqrt(dx**2 + dy**2)
        xWt[iCell] = dy/norm
        yWt[iCell] = dx/norm

    iCell = 0
    dx = (lonCells[iCell+1]-lonCells[iCell]) * np.cos(0.5*(latCells[iCell+1]+latCells[iCell]))
    dy = latCells[iCell+1]-latCells[iCell]
    norm = np.sqrt(dx**2 + dy**2)
    xWt[iCell] = dy/norm
    yWt[iCell] = dx/norm

    iCell = ntransectCells-1
    dx = (lonCells[iCell]-lonCells[iCell-1]) * np.cos(0.5*(latCells[iCell]+latCells[iCell-1]))
    dy = latCells[iCell]-latCells[iCell-1]
    norm = np.sqrt(dx**2 + dy**2)
    xWt[iCell] = dy/norm
    yWt[iCell] = dx/norm
    #print('xWt',xWt)
    #print('yWt',yWt)
    
    latmean = 180.0/np.pi*np.nanmean(latCells)
    lonmean = 180.0/np.pi*np.nanmean(lonCells)
    pressure = gsw.p_from_z(-z, latmean)
    pre = 'timeMonthly_avg_'

    # Load in T, S, and normalVelocity for each season, and plot them
    for season in seasonList:
        print('  season: ', season)
        for iSim in [0]: #range(len(simName)):
            fig = plt.figure(figsize=figsize, dpi=figdpi)
            print('    sim: ', simShortName[iSim])
            modeldir = rr + simName[iSim] + subdir[iSim]
#climofile = 'mpaso_ANN_198001_201412_climo.nc'; seasonName = 'ANN'
            print('modeldir, season',modeldir, season)

            modelfile = glob.glob('{}/mpaso_{}_*_climo.nc'.format(
                        modeldir, season))[0]
            meshSim = xr.open_dataset(meshfile[iSim])
            maxLevelCell = meshSim.maxLevelCell.sel(nCells=transectCells-1).values
            xr.Dataset.close(meshSim)
            # Initialize mask to True everywhere
            cellMask = np.ones((ntransectCells, nVertLevels), bool)
            for iCell in range(ntransectCells):
                # These become False if the second expression is negated (topography cells)
                cellMask[iCell, :] = np.logical_and(cellMask[iCell, :],
                                                     range(1, nVertLevels+1) <= maxLevelCell[iCell])
            print('modelfile',modelfile)
            ncid = Dataset(modelfile, 'r')

            #print('modelfile',modelfile)
            layerThickness = ncid.variables[pre+'layerThickness'][0, transectCells-1, :]
            zMid = np.zeros([ntransectCells,nVertLevels])
            for iCell in range(ntransectCells):
                zMid[iCell,:] = bottomDepth[iCell]
                zMid[iCell,0] =  0.5*layerThickness[iCell,0]
                for k in range(1,maxLevelCell[iCell]):
                   zMid[iCell,k] =  zMid[iCell,k-1] + 0.5*(layerThickness[iCell,k-1] + layerThickness[iCell,k])
            y = zMid


            # Load in T and S
            preT = pre + 'activeTracers_'
            temp = ncid.variables[preT + 'temperature'][0, transectCells-1, :]
            salt = ncid.variables[preT + 'salinity'][0, transectCells-1, :]
            velocityMeridional = ncid.variables[pre + 'velocityMeridional'][0, transectCells-1, :]
            velocityZonal = ncid.variables[pre + 'velocityZonal'][0, transectCells-1, :]
            velocityNormal = np.zeros([ntransectCells,nVertLevels])
            for iCell in range(ntransectCells):
                for k in range(1,maxLevelCell[iCell]):
                    # choose straight meridional or combo:
                    if transectNames == 'OSNAP section East':
                        velocityNormal[iCell,k] = velocityMeridional[iCell,k] # use meridional velocity on OSNAP East
                    elif transectNames == 'OSNAP section West':
                        # For OSNAP West, look at dx=1,dy=1, i.e. current directly east-southeastward
                        velocityNormal[iCell,k] = (-velocityZonal[iCell,k] + velocityMeridional[iCell,k])/np.sqrt(2)
                    else:
                        # use cell-by-cell weighting of velocity. This produced noisy plots.
                        velocityNormal[iCell,k] = xWt[iCell]*velocityZonal[iCell,k] + yWt[iCell]*velocityMeridional[iCell,k]
                    # nonlinear mapping to mimic https://www.nature.com/articles/s41467-021-23350-2/figures/1
                    if velocityNormal[iCell,k]>0.1:
                        velocityNormal[iCell,k] = 0.1 + 0.05/0.9*(velocityNormal[iCell,k] -0.1)
                    elif velocityNormal[iCell,k]<-0.1:
                        velocityNormal[iCell,k] = -0.1 + 0.05/0.9*(velocityNormal[iCell,k] +0.1)

            # Mask T,S values that fall on land and topography
            temp = np.ma.masked_array(temp, ~cellMask)
            salt = np.ma.masked_array(salt, ~cellMask)
            velocityMeridional = np.ma.masked_array(velocityMeridional, ~cellMask)
            velocityZonal = np.ma.masked_array(velocityZonal, ~cellMask)

            # Compute sigma's
            SA = gsw.SA_from_SP(salt, pressure[np.newaxis, :], lonmean, latmean)
            CT = gsw.CT_from_pt(SA, temp)
            sigma2 = gsw.density.sigma2(SA, CT)
            sigma0 = gsw.density.sigma0(SA, CT)

            #zmax = z[np.max(maxLevelCell)]
            zmax = np.max(bottomDepth)

            # Plot sections
            #  T first
            #figtitle = '{} Temperature, {}, {} years={}-{}'.format(
            #           simShortName[iSim], transectName, season, climoyearStart, climoyearEnd)
            #ax = plt.subplot(2,2,iSim*2+1)
            #ax.set_facecolor('darkgrey')
            #cf = ax.contourf(x, y, temp, cmap=colormapT, norm=cnormT, levels=clevelsT, extend='both')
            ##cf = ax.pcolormesh(x, y, temp, cmap=colormapT, norm=cnormT)
            #cax, kw = mpl.colorbar.make_axes(ax, location='right', pad=0.05, shrink=0.9)
            #cbar = plt.colorbar(cf, cax=cax, ticks=clevelsT, **kw)
            #cbar.ax.tick_params(labelsize=12, labelcolor='black')
            #cbar.set_label('C$^\circ$', fontsize=12, fontweight='bold')
            #if sigma2contours is not None:
            #    cs = ax.contour(x, y, sigma2, sigma2contours, colors='k', linewidths=1.5)
            #    cb = plt.clabel(cs, levels=sigma2contours, inline=True, inline_spacing=2, fmt='%2.1f', fontsize=9)
            #if sigma0contours is not None:
            #    cs = ax.contour(x, y, sigma0, sigma0contours, colors='k', linewidths=1.5)
            #    cb = plt.clabel(cs, levels=sigma0contours, inline=True, inline_spacing=2, fmt='%5.2f', fontsize=8)
            #!ax.set_ylim(0, zmax)
            #ax.set_ylim(0, 4000)
            #ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
            #ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
            #ax.set_title(figtitle, fontsize=12, fontweight='bold')
            #ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latCells[0]), xy=(0, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            #ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonCells[0]), xy=(0, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            #ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latCells[-1]), xy=(1, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            #ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonCells[-1]), xy=(1, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            #ax.invert_yaxis()

            #  then S or velocity
            figtitle = '{} {}, {}'.format(
                       simShortName[iSim], transectName, season)
            ax = plt.subplot(1,1,1)
            ax.set_facecolor('k')
            cf = ax.contourf(x, y, velocityNormal, cmap=colormapV, norm=cnormV, levels=clevelsV, extend='both')
            if include_colorbar:
                cax, kw = mpl.colorbar.make_axes(ax, location='right', pad=0.05, shrink=0.9)
                clevelVTicks = np.linspace(-0.14,0.14, 15)
                cbar = plt.colorbar(cf, cax=cax, ticks=clevelVTicks, **kw)
                yticklabelsV = np.round(np.linspace(-0.14,0.14, 15),2)
                yticklabelsV[0:2] = [-0.8,-0.4]
                yticklabelsV[13:15] = [0.4,0.8]
                cbar.set_ticklabels(yticklabelsV)
                cbar.ax.tick_params(labelsize=12, labelcolor='black')
            # add density contours:
            #if sigma2contours is not None:
            #    cs = ax.contour(x, y, sigma2, sigma2contours, colors='k', linewidths=1.5)
            #    cb = plt.clabel(cs, levels=sigma2contours, inline=True, inline_spacing=2, fmt='%2.1f', fontsize=9)
            #if sigma0contours is not None:
            #    cs = ax.contour(x, y, sigma0, sigma0contours, colors='k', linewidths=1.5)
            #    cb = plt.clabel(cs, levels=sigma0contours, inline=True, inline_spacing=2, fmt='%5.2f', fontsize=8)
            #ax.set_ylim(0, zmax)
            ax.set_ylim(0, 4000)
            if transectNames == 'OSNAP section East':
                ax.set_xlim(0, 2200) # OSNAP East only
            ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
            ax.set_title(figtitle, fontsize=12, fontweight='bold')
            ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latCells[0]), xy=(0, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonCells[0]), xy=(0, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latCells[-1]), xy=(1, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonCells[-1]), xy=(1, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            ax.invert_yaxis()

            ncid.close()
            # end for iSim in range(len(simName)):
            figfile = figdir + transectName.replace(' ', '') + '_v_' + simShortName[iSim] + '_' + season + '.png'
            plt.savefig(figfile, bbox_inches='tight')
            plt.close()
    # end for season in [1]: #months:
# end for iTransect in range(nTransects):
