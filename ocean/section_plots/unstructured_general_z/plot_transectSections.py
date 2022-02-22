#!/usr/bin/env python
"""
Plot vertical sections (annual and seasonal climatologies),
given a certain transect mask
"""
############################## transects
transectNames = ['all']
transectNames = ['Faroe Bank Ch N','Faroe Bank Ch','Faroe Shetland Ch']

############################## months or seasons
#seasonList = ['JFM', 'JAS', 'ANN']
seasonList = ['ANN']

############################## model files, run dirs
rr = '/lcrc/group/e3sm/ac.mpetersen/scratch/anvil/'
simName = ['20211006g.normal.WCYCL1850.ne30pg2_EC30to60E2r2.sigmaz.anvil/',
  '211006h.sigmaz-hMin1.WCYCL1850.ne30pg2_EC30to60E2r2.anvil/']
simShortName= ['z-level',
  'sigma-z']
meshfile = ['init_zlevel.nc','init_hMin1_14.nc']
# extra simName = '20211006f.init8.WCYCL1850.ne30pg2_EC30to60E2r2.sigmaz.anvil/'; simShortName='sigma-z_8'
subdir = 'yrs21-30/clim/mpas/avg/unmasked_EC30to60E2r2/'
#climofile = 'mpaso_JFM_002101_003003_climo.nc'; seasonName = 'JFM'
#climofile = 'mpaso_JAS_002107_003009_climo.nc'; seasonName = 'JAS'
pre = 'timeMonthly_avg_'
maskfile = 'mask.nc'
casename = ''; 'E3SM60to30' # no spaces
climoyearStart = 21
climoyearEnd = 30 

############################## variables
# not added yet

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
figdir = './verticalSections/{}'.format(casename)
if not os.path.isdir(figdir):
    os.makedirs(figdir)
figsize = (20, 12)
figdpi = 300
colorIndices0 = [0, 10, 28, 57, 85, 113, 125, 142, 155, 170, 198, 227, 242, 255]
#clevelsT = [-2.0, -1.8, -1.5, -1.0, -0.5, 0.0, 0.5, 2.0, 4.0, 6.0, 8.0, 10., 12.]
#clevelsS = [30.0, 31.0, 32.0, 33.0, 33.5, 34.0, 34.5, 34.8, 34.85, 34.9, 34.95, 35.0, 35.5]
# Better for OSNAP:
clevelsT = [-1.0, -0.5, 0.0, 0.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 8., 10., 12.]
#clevelsT = np.linspace(2.0, 6.0, 13)
clevelsS = [31.0, 33.0, 33.5, 33.8, 34.2, 34.6, 34.8, 34.85, 34.9, 34.95, 35.0, 35.2, 35.5]
#clevelsS = np.linspace(34.7, 35.2, 13)
clevelsV = [-0.25, -0.2, -0.15, -0.1, -0.02, 0.0, 0.02, 0.1, 0.2, 0.3, 0.5]
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
colormapV = cols.ListedColormap(colormapV(colorIndices))
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
mask = xr.open_dataset(maskfile)

allTransects = decode_strings(mask.transectNames)
if transectNames[0]=='all' or transectNames[0]=='StandardTransportSectionsRegionsGroup':
    transectNames = allTransects

# Get depth
z = mesh.refBottomDepth.values
nlevels = len(z)

nTransects = len(transectNames)
maxEdges = mask.dims['maxEdgesInTransect']
for iTransect in range(nTransects):
    # Identify transect
    transectName = transectNames[iTransect]
    transectIndex = allTransects.index(transectName)
    print('transect: ', transectName)

    # Choose mask for this particular transect
    transectmask = mask.isel(nTransects=transectIndex).squeeze()
    # Get a list of edges for this transect
    transectEdges = transectmask.transectEdgeGlobalIDs.values
    transectEdges = transectEdges[np.where(transectEdges > 0)[0]]
    ntransectEdges = len(transectEdges)

    # Get a list of cellsOnEdge pairs for each transect edge
    cellsOnEdge = mesh.cellsOnEdge.sel(nEdges=transectEdges-1).values

    # Create a land/topo mask for cellsOnEdge
    cellsOnEdge1 = cellsOnEdge[:, 0]
    cellsOnEdge2 = cellsOnEdge[:, 1]

    # Get edge signs for across-edge velocity direction
    edgeSigns = mask.transectEdgeMaskSigns.sel(nEdges=transectEdges-1, nTransects=transectIndex).values
    # Get coordinates of each edge center and compute approximate spherical distance
    lonEdges = mesh.lonEdge.sel(nEdges=transectEdges-1).values
    latEdges = mesh.latEdge.sel(nEdges=transectEdges-1).values
    lonEdges[lonEdges>np.pi] = lonEdges[lonEdges>np.pi] - 2*np.pi

    # load layer thickness

    BDOnCell1 = mesh.variables['bottomDepth'][cellsOnEdge1-1]
    BDOnCell2 = mesh.variables['bottomDepth'][cellsOnEdge2-1]
    bD = np.nanmean(np.array([BDOnCell1, BDOnCell2]), axis=0)

    dist = [0]
    for iEdge in range(1, ntransectEdges):
        dx = (lonEdges[iEdge]-lonEdges[iEdge-1]) * np.cos(0.5*(latEdges[iEdge]+latEdges[iEdge-1]))
        dy = latEdges[iEdge]-latEdges[iEdge-1]
        dist.append(earthRadius * np.sqrt(dx**2 + dy**2))
    dist = np.cumsum(dist)
    [x, y] = np.meshgrid(dist, z)
    x = x.T
    
    latmean = 180.0/np.pi*np.nanmean(latEdges)
    lonmean = 180.0/np.pi*np.nanmean(lonEdges)
    pressure = gsw.p_from_z(-z, latmean)

    # Load in T, S, and normalVelocity for each season, and plot them
    for season in seasonList:
        print('  season: ', season)
        fig = plt.figure(figsize=figsize, dpi=figdpi)
        for iSim in range(len(simName)):
            print('    sim: ', simShortName[iSim])
            modeldir = rr + simName[iSim] + subdir
            modelfile = glob.glob('{}/mpaso_{}_{:04d}??_{:04d}??_climo.nc'.format(
                        modeldir, season, climoyearStart, climoyearEnd))[0]
            ncid = Dataset(modelfile, 'r')
            #print('modelfile',modelfile)
            hOnCell1 = ncid.variables[pre + 'layerThickness'][0, cellsOnEdge1-1, :]
            hOnCell2 = ncid.variables[pre + 'layerThickness'][0, cellsOnEdge2-1, :]
            LTh = np.nanmean(np.array([hOnCell1, hOnCell2]), axis=0)
            zMid = np.zeros([ntransectEdges,nlevels])
            for iEdge in range(ntransectEdges):
                zMid[iEdge,0] =  0.5*LTh[iEdge,0]
                for k in range(1,nlevels):
                   zMid[iEdge,k] =  zMid[iEdge,k-1] + 0.5*(LTh[iEdge,k-1] + LTh[iEdge,k])
            y = zMid

            meshSim = xr.open_dataset(meshfile[iSim])
            maxLevelCell1 = meshSim.maxLevelCell.sel(nCells=cellsOnEdge1-1).values
            maxLevelCell2 = meshSim.maxLevelCell.sel(nCells=cellsOnEdge2-1).values
            xr.Dataset.close(meshSim)
            # Initialize mask to True everywhere
            cellMask1 = np.ones((ntransectEdges, nlevels), bool)
            cellMask2 = np.ones((ntransectEdges, nlevels), bool)
            for iEdge in range(ntransectEdges):
                # These become False if the second expression is negated (land cells)
                cellMask1[iEdge, :] = np.logical_and(cellMask1[iEdge, :],
                                                     cellsOnEdge1[iEdge, np.newaxis] > 0)
                cellMask2[iEdge, :] = np.logical_and(cellMask2[iEdge, :],
                                                     cellsOnEdge2[iEdge, np.newaxis] > 0)
                # These become False if the second expression is negated (topography cells)
                cellMask1[iEdge, :] = np.logical_and(cellMask1[iEdge, :],
                                                     range(1, nlevels+1) <= maxLevelCell1[iEdge])
                cellMask2[iEdge, :] = np.logical_and(cellMask2[iEdge, :],
                                                     range(1, nlevels+1) <= maxLevelCell2[iEdge])

            # Create a land/topo mask for transectEdges
            maxLevelEdge = []
            for iEdge in range(ntransectEdges):
                if cellsOnEdge1[iEdge]==0:
                    maxLevelEdge.append(maxLevelCell2[iEdge])
                elif cellsOnEdge2[iEdge]==0:
                    maxLevelEdge.append(maxLevelCell1[iEdge])
                else:
                    maxLevelEdge.append(np.min([maxLevelCell1[iEdge], maxLevelCell2[iEdge]]))
            # Initialize mask to True everywhere
            edgeMask = np.ones((ntransectEdges, nlevels), bool)
            for iEdge in range(ntransectEdges):
                # These become False if the second expression is negated (topography cells)
                edgeMask[iEdge, :] = np.logical_and(edgeMask[iEdge, :],
                                                    range(1, nlevels+1) <= maxLevelEdge[iEdge])
            ## Try loading in normalVelocity (on edge centers)
            #try:
            #    vel = ncid.variables['timeMonthly_avg_normalVelocity'][0, transectEdges-1, :]
            #    if 'timeMonthly_avg_normalGMBolusVelocity' in ncid.variables.keys():
            #        vel += ncid.variables['timeMonthly_avg_normalGMBolusVelocity'][0, transectEdges-1, :]
            #except:
            #    #print('*** normalVelocity variable not found: skipping it...')
            #    vel = None
            # Load in T and S (on cellsOnEdge centers)
            preT = pre + 'activeTracers_'
            tempOnCell1 = ncid.variables[preT + 'temperature'][0, cellsOnEdge1-1, :]
            tempOnCell2 = ncid.variables[preT + 'temperature'][0, cellsOnEdge2-1, :]
            saltOnCell1 = ncid.variables[preT + 'salinity'][0, cellsOnEdge1-1, :]
            saltOnCell2 = ncid.variables[preT + 'salinity'][0, cellsOnEdge2-1, :]

            # Mask T,S values that fall on land and topography
            tempOnCell1 = np.ma.masked_array(tempOnCell1, ~cellMask1)
            tempOnCell2 = np.ma.masked_array(tempOnCell2, ~cellMask2)
            saltOnCell1 = np.ma.masked_array(saltOnCell1, ~cellMask1)
            saltOnCell2 = np.ma.masked_array(saltOnCell2, ~cellMask2)
            # Interpolate T,S values onto edges
            temp = np.nanmean(np.ma.array([tempOnCell1, tempOnCell2]), axis=0)
            salt = np.nanmean(np.ma.array([saltOnCell1, saltOnCell2]), axis=0)

            # Compute sigma's
            SA = gsw.SA_from_SP(salt, pressure[np.newaxis, :], lonmean, latmean)
            CT = gsw.CT_from_pt(SA, temp)
            sigma2 = gsw.density.sigma2(SA, CT)
            sigma0 = gsw.density.sigma0(SA, CT)

            #zmax = z[np.max(maxLevelEdge)]
            zmax = np.max(bD)

            # Plot sections
            #  T first
            figtitle = '{} Temperature, {}, {} years={}-{}'.format(
                       simShortName[iSim], transectName, season, climoyearStart, climoyearEnd)
            ax = plt.subplot(2,2,iSim*2+1)
            ax.set_facecolor('darkgrey')
            cf = ax.contourf(x, y, temp, cmap=colormapT, norm=cnormT, levels=clevelsT, extend='both')
            #cf = ax.pcolormesh(x, y, temp, cmap=colormapT, norm=cnormT)
            cax, kw = mpl.colorbar.make_axes(ax, location='right', pad=0.05, shrink=0.9)
            cbar = plt.colorbar(cf, cax=cax, ticks=clevelsT, **kw)
            cbar.ax.tick_params(labelsize=12, labelcolor='black')
            cbar.set_label('C$^\circ$', fontsize=12, fontweight='bold')
            if sigma2contours is not None:
                cs = ax.contour(x, y, sigma2, sigma2contours, colors='k', linewidths=1.5)
                cb = plt.clabel(cs, levels=sigma2contours, inline=True, inline_spacing=2, fmt='%2.1f', fontsize=9)
            if sigma0contours is not None:
                cs = ax.contour(x, y, sigma0, sigma0contours, colors='k', linewidths=1.5)
                cb = plt.clabel(cs, levels=sigma0contours, inline=True, inline_spacing=2, fmt='%5.2f', fontsize=8)
            ax.set_ylim(0, zmax)
            ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
            ax.set_title(figtitle, fontsize=12, fontweight='bold')
            ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[0]), xy=(0, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[0]), xy=(0, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[-1]), xy=(1, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[-1]), xy=(1, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            ax.invert_yaxis()

            #  then S
            figtitle = '{} Salinity, {}, {} years={}-{}'.format(
                       simShortName[iSim], transectName, season, climoyearStart, climoyearEnd)
            ax = plt.subplot(2,2,iSim*2+2)
            ax.set_facecolor('darkgrey')
            # new mrp for colormap
            #if iSim<3: #==1:
            #    clevelsS = np.linspace(np.ma.min(salt), np.ma.max(salt), 13)
            #    colormapS = cols.ListedColormap(colormapS(colorIndices))
            #    colormapS.set_under(underColor)
            #    colormapS.set_over(overColor)
            #    cnormS = mpl.colors.BoundaryNorm(clevelsS, colormapS.N)
                # end new mrp
            cf = ax.contourf(x, y, salt, cmap=colormapS, norm=cnormS, levels=clevelsS, extend='both')
            cax, kw = mpl.colorbar.make_axes(ax, location='right', pad=0.05, shrink=0.9)
            cbar = plt.colorbar(cf, cax=cax, ticks=clevelsS, **kw)
            cbar.ax.tick_params(labelsize=12, labelcolor='black')
            cbar.set_label('psu', fontsize=12, fontweight='bold')
            cf = ax.plot(x, y)
            #if sigma2contours is not None:
            #    cs = ax.contour(x, y, sigma2, sigma2contours, colors='k', linewidths=1.5)
            #    cb = plt.clabel(cs, levels=sigma2contours, inline=True, inline_spacing=2, fmt='%2.1f', fontsize=9)
            #if sigma0contours is not None:
            #    cs = ax.contour(x, y, sigma0, sigma0contours, colors='k', linewidths=1.5)
            #    cb = plt.clabel(cs, levels=sigma0contours, inline=True, inline_spacing=2, fmt='%5.2f', fontsize=8)
            ax.set_ylim(0, zmax)
            ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
            ax.set_title(figtitle, fontsize=12, fontweight='bold')
            ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[0]), xy=(0, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[0]), xy=(0, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[-1]), xy=(1, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[-1]), xy=(1, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            ax.invert_yaxis()

            #  and finally normalVelocity (if vel is not None)
            #if vel is not None:
            #    # Mask velocity values that fall on land and topography
            #    vel = np.ma.masked_array(vel, ~edgeMask)
            #    # Get normalVelocity direction
            #    normalVel = vel*edgeSigns[:, np.newaxis]

            #    figtitle = 'Velocity ({}), {} ({}, years={}-{})'.format(
            #               transectName, season, casename, climoyearStart, climoyearEnd)
            #    figfile = '{}/Vel_{}_{}_{}_years{:04d}-{:04d}.png'.format(
            #               figdir, transectName.replace(' ', ''), casename, season, climoyearStart, climoyearEnd)
            #    fig = plt.figure(figsize=figsize, dpi=figdpi)
            #    ax = fig.subplot()
            #    ax.set_facecolor('darkgrey')
            #    cf = ax.contourf(x, y, normalVel, cmap=colormapV, norm=cnormV, levels=clevelsV)
            #    cax, kw = mpl.colorbar.make_axes(ax, location='right', pad=0.05, shrink=0.9)
            #    cbar = plt.colorbar(cf, cax=cax, ticks=clevelsV, **kw)
            #    cbar.ax.tick_params(labelsize=12, labelcolor='black')
            #    cbar.set_label('m/s', fontsize=12, fontweight='bold')
            #    ax.set_ylim(0, zmax)
            #    ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
            #    ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
            #    ax.set_title(figtitle, fontsize=12, fontweight='bold')
            #    ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[0]), xy=(0, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            #    ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[0]), xy=(0, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            #    ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[-1]), xy=(1, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            #    ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[-1]), xy=(1, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            #    ax.invert_yaxis()
            #    plt.savefig(figfile) #, bbox_inches='tight')
            #    plt.close()

            ncid.close()
        # end for iSim in range(len(simName)):
        figfile = figdir + transectName.replace(' ', '') + '_' + simShortName[0] + '_' + simShortName[1] + '_' + season + '.png'
        plt.savefig(figfile, bbox_inches='tight')
        plt.close()
    # end for season in [1]: #months:
# end for iTransect in range(nTransects):
