#!/usr/bin/env python
"""

Plot vertical sections (annual and seasonal climatologies),
given a certain transect mask

"""

# ensure plots are rendered on ICC
from __future__ import absolute_import, division, print_function, \
    unicode_literals
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

meshfile = 'mesh.nc'
maskfile = 'mask.nc'
casename = 'E3SM60to30' # no spaces
climoyearStart = 40
climoyearEnd = 59
modeldir = '/lustre/scratch4/turquoise/mpeterse/runs/211006_sigma-z_meshes/ocean/global_ocean/EC30to60/PHC/init/hMin1'
#seasons = ['JFM', 'JAS', 'ANN']
seasons = ['ANN']

## Options for transect names if maskfile=*_standardTransportSections.nc
# "Africa-Ant", "Agulhas", "Antilles Inflow", "Barents Sea Opening", "Bering Strait", "Davis Strait",
# "Drake Passage", "Florida-Bahamas", "Florida-Cuba", "Fram Strait", "Indonesian Throughflow",
# "Lancaster Sound", "Mona Passage", "Mozambique Channel", "Nares Strait", "Tasmania-Ant", "Windward Passage"
## Options for transect names if maskfile=*_arcticSections.nc
# "Barents Sea Opening", "Bering Strait", "Davis Strait", "Denmark Strait", "Fram Strait", 
# "Hudson Bay-Labrador Sea", "Iceland-Faroe-Scotland", "Lancaster Sound", "Nares Strait",
# "OSNAP section East", "OSNAP section West"
#transectNames = ['all']
#transectNames = ['Barents Sea Opening', 'Fram Strait']
#transectNames = ['Barents Sea Opening', 'Bering Strait', 'Davis Strait',
#                 'Denmark Strait', 'Fram Strait', 'Iceland-Faroe-Scotland']
transectNames = ['Atlantic Section OSNAP East']

# Figure details
figdir = './verticalSections/{}'.format(casename)
if not os.path.isdir(figdir):
    os.makedirs(figdir)
figsize = (10, 6)
figdpi = 300
colorIndices0 = [0, 10, 28, 57, 85, 113, 125, 142, 155, 170, 198, 227, 242, 255]
#clevelsT = [-2.0, -1.8, -1.5, -1.0, -0.5, 0.0, 0.5, 2.0, 4.0, 6.0, 8.0, 10., 12.]
#clevelsS = [30.0, 31.0, 32.0, 33.0, 33.5, 34.0, 34.5, 34.8, 34.85, 34.9, 34.95, 35.0, 35.5]
# Better for OSNAP:
clevelsT = [-1.0, -0.5, 0.0, 0.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 8., 10., 12.]
clevelsS = [31.0, 33.0, 33.5, 33.8, 34.2, 34.6, 34.8, 34.85, 34.9, 34.95, 35.0, 35.2, 35.5]
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
sigma0contours = [24.0, 25.0, 26.0, 27.0, 27.2, 27.4, 27.6, 27.8, 28.0]
#sigma0contours = None

# Load in MPAS mesh and transect mask file
mesh = xr.open_dataset(meshfile)
mask = xr.open_dataset(maskfile)

allTransects = decode_strings(mask.transectNames)
if transectNames[0]=='all' or transectNames[0]=='StandardTransportSectionsRegionsGroup':
    transectNames = allTransects

# Get depth
z = mesh.refBottomDepth.values
nlevels = len(z)

nTransects = len(transectNames)
maxEdges = mask.dims['maxEdgesInTransect']
for n in range(nTransects):
    # Identify transect
    transectName = transectNames[n]
    transectIndex = allTransects.index(transectName)
    print('Plotting sections for transect: ', transectName)

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
    maxLevelCell1 = mesh.maxLevelCell.sel(nCells=cellsOnEdge1-1).values
    maxLevelCell2 = mesh.maxLevelCell.sel(nCells=cellsOnEdge2-1).values
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

    # Get edge signs for across-edge velocity direction
    edgeSigns = mask.transectEdgeMaskSigns.sel(nEdges=transectEdges-1, nTransects=transectIndex).values
    # Get coordinates of each edge center and compute approximate spherical distance
    lonEdges = mesh.lonEdge.sel(nEdges=transectEdges-1).values
    latEdges = mesh.latEdge.sel(nEdges=transectEdges-1).values
    lonEdges[lonEdges>np.pi] = lonEdges[lonEdges>np.pi] - 2*np.pi
    dist = [0]
    for iEdge in range(1, ntransectEdges):
        dx = (lonEdges[iEdge]-lonEdges[iEdge-1]) * np.cos(0.5*(latEdges[iEdge]+latEdges[iEdge-1]))
        dy = latEdges[iEdge]-latEdges[iEdge-1]
        dist.append(earthRadius * np.sqrt(dx**2 + dy**2))
    dist = np.cumsum(dist)
    [x, y] = np.meshgrid(dist, z)
    x = x.T
    y = y.T
    # Check lon,lat of edges to make sure we have the right edges
    #print(180.0/np.pi*lonEdges)
    #print(180.0/np.pi*latEdges)
    # Check lon,lat of cells to make sure we have the right cellsOnEdge
    #print('lonCell0=', 180/np.pi*mesh.lonCell.sel(nCells=cellsOnEdge1-1).values)
    #print('latCell0=', 180/np.pi*mesh.latCell.sel(nCells=cellsOnEdge1-1).values)
    #print('lonCell1=', 180/np.pi*mesh.lonCell.sel(nCells=cellsOnEdge2-1).values)
    #print('latCell1=', 180/np.pi*mesh.latCell.sel(nCells=cellsOnEdge2-1).values)
    latmean = 180.0/np.pi*np.nanmean(latEdges)
    lonmean = 180.0/np.pi*np.nanmean(lonEdges)
    pressure = gsw.p_from_z(-z, latmean)

    # Load in T, S, and normalVelocity for each season, and plot them
    for s in [1]: #seasons:
        #print('   season: ', s)
        #modelfile = glob.glob('{}/mpaso_{}_{:04d}??_{:04d}??_climo.nc'.format(
        #            modeldir, s, climoyearStart, climoyearEnd))[0]
        modelfile = 'init_hMin1_14.nc'
        ncid = Dataset(modelfile, 'r')
        ## Try loading in normalVelocity (on edge centers)
        try:
            vel = ncid.variables['timeMonthly_avg_normalVelocity'][0, transectEdges-1, :]
            if 'timeMonthly_avg_normalGMBolusVelocity' in ncid.variables.keys():
                vel += ncid.variables['timeMonthly_avg_normalGMBolusVelocity'][0, transectEdges-1, :]
        except:
            print('*** normalVelocity variable not found: skipping it...')
            vel = None
        # Load in T and S (on cellsOnEdge centers)
        tempOnCell1 = ncid.variables['temperature'][0, cellsOnEdge1-1, :]
        tempOnCell2 = ncid.variables['temperature'][0, cellsOnEdge2-1, :]
        saltOnCell1 = ncid.variables['salinity'][0, cellsOnEdge1-1, :]
        saltOnCell2 = ncid.variables['salinity'][0, cellsOnEdge2-1, :]
        ncid.close()

        # Mask T,S values that fall on land and topography
        tempOnCell1 = np.ma.masked_array(tempOnCell1, ~cellMask1)
        tempOnCell2 = np.ma.masked_array(tempOnCell2, ~cellMask2)
        saltOnCell1 = np.ma.masked_array(saltOnCell1, ~cellMask1)
        saltOnCell2 = np.ma.masked_array(saltOnCell2, ~cellMask2)
        # Interpolate T,S values onto edges
        temp = np.nanmean(np.array([tempOnCell1, tempOnCell2]), axis=0)
        salt = np.nanmean(np.array([saltOnCell1, saltOnCell2]), axis=0)

        # Compute sigma's
        SA = gsw.SA_from_SP(salt, pressure[np.newaxis, :], lonmean, latmean)
        CT = gsw.CT_from_pt(SA, temp)
        sigma2 = gsw.density.sigma2(SA, CT)
        sigma0 = gsw.density.sigma0(SA, CT)

        zmax = z[np.max(maxLevelEdge)]

        # Plot sections
        #  T first
        figtitle = 'Temperature ({}), {} ({}, years={}-{})'.format(
                   transectName, s, casename, climoyearStart, climoyearEnd)
        figfile = '{}/Temp_{}_{}_{}_years{:04d}-{:04d}.png'.format(
                  figdir, transectName.replace(' ', ''), casename, s, climoyearStart, climoyearEnd)
        figfile = 'init_hMin1_14.png'
        fig = plt.figure(figsize=figsize, dpi=figdpi)
        ax = fig.add_subplot()
        ax.set_facecolor('darkgrey')
        cf = ax.contourf(x, y, temp, cmap=colormapT, norm=cnormT, levels=clevelsT, extend='max')
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
            cb = plt.clabel(cs, levels=sigma0contours, inline=True, inline_spacing=2, fmt='%2.1f', fontsize=9)
        ax.set_ylim(0, zmax)
        ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
        ax.set_title(figtitle, fontsize=12, fontweight='bold')
        ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[0]), xy=(0, -0.1), xycoords='axes fraction', ha='center', va='bottom')
        ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[0]), xy=(0, -0.15), xycoords='axes fraction', ha='center', va='bottom')
        ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[-1]), xy=(1, -0.1), xycoords='axes fraction', ha='center', va='bottom')
        ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[-1]), xy=(1, -0.15), xycoords='axes fraction', ha='center', va='bottom')
        ax.invert_yaxis()
        plt.savefig(figfile, bbox_inches='tight')
        plt.close()

        #  then S
        figtitle = 'Salinity ({}), {} ({}, years={}-{})'.format(
                   transectName, s, casename, climoyearStart, climoyearEnd)
        figfile = '{}/Salt_{}_{}_{}_years{:04d}-{:04d}.png'.format(
                  figdir, transectName.replace(' ', ''), casename, s, climoyearStart, climoyearEnd)
        fig = plt.figure(figsize=figsize, dpi=figdpi)
        ax = fig.add_subplot()
        ax.set_facecolor('darkgrey')
        cf = ax.contourf(x, y, salt, cmap=colormapS, norm=cnormS, levels=clevelsS, extend='max')
        cax, kw = mpl.colorbar.make_axes(ax, location='right', pad=0.05, shrink=0.9)
        cbar = plt.colorbar(cf, cax=cax, ticks=clevelsS, **kw)
        cbar.ax.tick_params(labelsize=12, labelcolor='black')
        cbar.set_label('psu', fontsize=12, fontweight='bold')
        if sigma2contours is not None:
            cs = ax.contour(x, y, sigma2, sigma2contours, colors='k', linewidths=1.5)
            cb = plt.clabel(cs, levels=sigma2contours, inline=True, inline_spacing=2, fmt='%2.1f', fontsize=9)
        if sigma0contours is not None:
            cs = ax.contour(x, y, sigma0, sigma0contours, colors='k', linewidths=1.5)
            cb = plt.clabel(cs, levels=sigma0contours, inline=True, inline_spacing=2, fmt='%2.1f', fontsize=9)
        ax.set_ylim(0, zmax)
        ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
        ax.set_title(figtitle, fontsize=12, fontweight='bold')
        ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[0]), xy=(0, -0.1), xycoords='axes fraction', ha='center', va='bottom')
        ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[0]), xy=(0, -0.15), xycoords='axes fraction', ha='center', va='bottom')
        ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[-1]), xy=(1, -0.1), xycoords='axes fraction', ha='center', va='bottom')
        ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[-1]), xy=(1, -0.15), xycoords='axes fraction', ha='center', va='bottom')
        ax.invert_yaxis()
        plt.savefig(figfile, bbox_inches='tight')
        plt.close()

        #  and finally normalVelocity (if vel is not None)
        if vel is not None:
            # Mask velocity values that fall on land and topography
            vel = np.ma.masked_array(vel, ~edgeMask)
            # Get normalVelocity direction
            normalVel = vel*edgeSigns[:, np.newaxis]

            figtitle = 'Velocity ({}), {} ({}, years={}-{})'.format(
                       transectName, s, casename, climoyearStart, climoyearEnd)
            figfile = '{}/Vel_{}_{}_{}_years{:04d}-{:04d}.png'.format(
                       figdir, transectName.replace(' ', ''), casename, s, climoyearStart, climoyearEnd)
            fig = plt.figure(figsize=figsize, dpi=figdpi)
            ax = fig.add_subplot()
            ax.set_facecolor('darkgrey')
            cf = ax.contourf(x, y, normalVel, cmap=colormapV, norm=cnormV, levels=clevelsV)
            cax, kw = mpl.colorbar.make_axes(ax, location='right', pad=0.05, shrink=0.9)
            cbar = plt.colorbar(cf, cax=cax, ticks=clevelsV, **kw)
            cbar.ax.tick_params(labelsize=12, labelcolor='black')
            cbar.set_label('m/s', fontsize=12, fontweight='bold')
            ax.set_ylim(0, zmax)
            ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
            ax.set_title(figtitle, fontsize=12, fontweight='bold')
            ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[0]), xy=(0, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[0]), xy=(0, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lat={:5.2f}'.format(180.0/np.pi*latEdges[-1]), xy=(1, -0.1), xycoords='axes fraction', ha='center', va='bottom')
            ax.annotate('lon={:5.2f}'.format(180.0/np.pi*lonEdges[-1]), xy=(1, -0.15), xycoords='axes fraction', ha='center', va='bottom')
            ax.invert_yaxis()
            plt.savefig(figfile, bbox_inches='tight')
            plt.close()
