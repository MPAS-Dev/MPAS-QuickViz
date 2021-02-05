from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cols
from matplotlib.pyplot import cm
from matplotlib.colors import from_levels_and_colors
from matplotlib.colors import BoundaryNorm
import cmocean

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf
from mpas_analysis.shared.io.utility import get_files_year_month, decode_strings

from geometric_features import FeatureCollection, read_feature_collection

from common_functions import hovmoeller_plot, add_inset

meshName = 'EC30to60E2r2'
restartFile = '/compyfs/inputdata/ocn/mpas-o/{}/ocean.EC30to60E2r2.200908.nc'.format(meshName)
regionMaskFile = '/compyfs/vene705/region_masks/{}_oceanOHCregions20201120.nc'.format(meshName)
featureFile = '/compyfs/vene705/region_masks/oceanOHCregions20201120.geojson'

modeldir = '/compyfs/zhen797/E3SM_simulations/20201108.alpha5_55_fallback.piControl.ne30pg2_r05_EC30to60E2r2-1900_ICG.compy/archive/ocn/hist'
runName = '20201108.alpha5_55_fallback.piControl.ne30pg2_r05_EC30to60E2r2-1900_ICG.compy'
runNameShort = 'alpha5_55_fallback'

outdir = './timeseries_data/{}'.format(runNameShort)
if not os.path.isdir(outdir):
    os.makedirs(outdir)
figdir = './timeseries/{}'.format(runNameShort)
if not os.path.isdir(figdir):
    os.makedirs(figdir)

if os.path.exists(restartFile):
    dsRestart = xarray.open_dataset(restartFile)
    dsRestart = dsRestart.isel(Time=0)
else:
    raise IOError('No MPAS restart/mesh file found')
areaCell = dsRestart.areaCell
if 'landIceMask' in dsRestart:
    # only the region outside of ice-shelf cavities
    openOceanMask = dsRestart.landIceMask == 0
else:
    openOceanMask = None
refBottomDepth = dsRestart.refBottomDepth
maxLevelCell = dsRestart.maxLevelCell
nVertLevels = dsRestart.sizes['nVertLevels']
vertIndex = xarray.DataArray.from_dict(
    {'dims': ('nVertLevels',), 'data': np.arange(nVertLevels)})
depthMask = (vertIndex < maxLevelCell).transpose('nCells', 'nVertLevels')

if os.path.exists(regionMaskFile):
    dsRegionMask = xarray.open_dataset(regionMaskFile)
    regionNames = decode_strings(dsRegionMask.regionNames)
    regionNames.append('Global')
    nRegions = np.size(regionNames)
else:
    raise IOError('No regional mask file found')

startYear = 1
endYear = 200
calendar = 'gregorian'

variables = [{'name': 'ohc',
              'title': 'OHC',
              'units': 'x10$^{22}$ J',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'colormap': cmocean.cm.balance,
              'clevels': [-2.4, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 2.4],
              'colorIndices': [0, 28, 57, 85, 113, 142, 170, 198, 227, 255],
              'fac': 1e-22*1026.0*3996.0}, # 1e-22*rho0*cp (where rho0=config_density0 and cp=config_specific_heat_sea_water)
             {'name': 'temperature',
              'title': 'Temperature',
              'units': '$^\circ$C',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'colormap': cmocean.cm.balance,
              'clevels': [-2.0, -1.0, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1.0, 2.0],
              'colorIndices': [0, 28, 57, 85, 113, 125, 130, 142, 170, 198, 227, 255],
              'fac': 1},
             {'name': 'salinity',
              'title': 'Salinity',
              'units': 'psu',
              'mpas': 'timeMonthly_avg_activeTracers_salinity',
              'colormap': cmocean.cm.balance,
              'clevels': [-1.0, -0.5, -0.2, -0.1, -0.02, 0, 0.02, 0.1, 0.2, 0.5, 1.0],
              'colorIndices': [0, 28, 57, 85, 113, 125, 130, 142, 170, 198, 227, 255],
              'fac': 1}]

startDate = '{:04d}-01-01_00:00:00'.format(startYear)
endDate = '{:04d}-12-31_23:59:59'.format(endYear)
years = range(startYear, endYear + 1)

variableList = [var['mpas'] for var in variables] + \
    ['timeMonthly_avg_layerThickness']

timeSeriesFile0 = '{}/OHC_T_S_trends_vsdepth'.format(outdir)

# Compute regional averages one year at a time
for year in years:

    timeSeriesFile = '{}_year{:04d}.nc'.format(timeSeriesFile0, year)

    if not os.path.exists(timeSeriesFile):
        print('\nComputing regional time series for year={}'.format(year))

        datasets = []
        for month in range(1, 13):
            inputFile = '{}/{}.mpaso.hist.am.timeSeriesStatsMonthly.{:04d}-{:02d}-01.nc'.format(
                modeldir, runName, year, month)
            if not os.path.exists(inputFile):
                raise IOError('Input file: {} not found'.format(inputFile))

            dsTimeSlice = open_mpas_dataset(fileName=inputFile,
                                            calendar=calendar,
                                            variableList=variableList,
                                            startDate=startDate,
                                            endDate=endDate)
            datasets.append(dsTimeSlice)
        # combine data sets into a single data set
        dsIn = xarray.concat(datasets, 'Time')

        # Global depth-masked layer thickness and layer volume
        layerThickness = dsIn.timeMonthly_avg_layerThickness
        layerThickness = layerThickness.where(depthMask, drop=False)
        layerVol = areaCell*layerThickness

        datasets = []
        regionIndices = []
        for regionName in regionNames:
            print('    region: {}'.format(regionName))

            # Compute region total area and, for regionName
            # other than Global, compute regional mask and
            # regionally masked layer volume
            if regionName=='Global':
                regionIndices.append(nRegions-1)

                totalArea = areaCell.sum()
                if year==years[0]:
                    print('      totalArea: {} mil. km^2'.format(1e-12*totalArea.values))
            else:
                regionIndex = regionNames.index(regionName)
                regionIndices.append(regionIndex)

                dsMask = dsRegionMask.isel(nRegions=regionIndex)
                cellMask = dsMask.regionCellMasks == 1
                if openOceanMask is not None:
                    cellMask = np.logical_and(cellMask, openOceanMask)

                localArea = areaCell.where(cellMask, drop=True)
                totalArea = localArea.sum()
                if year==years[0]:
                    print('      totalArea: {} mil. km^2'.format(1e-12*totalArea.values))
                localLayerVol = layerVol.where(cellMask, drop=True)

            # Temporary dsOut (xarray dataset) containing results for
            # all variables for one single region
            dsOut = xarray.Dataset()
            # Compute layer-volume weighted averages (or sums for OHC)
            for var in variables:
                outName = var['name']
                mpasVarName = var['mpas']
                units = var['units']
                factor = var['fac']
                description = var['title']

                timeSeries = dsIn[mpasVarName]
                timeSeries = timeSeries.where(depthMask, drop=False)
                if regionName=='Global':
                    timeSeries = (layerVol*timeSeries).sum(dim='nCells')
                    if outName!='ohc':
                        timeSeries = timeSeries / layerVol.sum(dim='nCells')
                else:
                    timeSeries = timeSeries.where(cellMask, drop=True)
                    timeSeries = (localLayerVol*timeSeries).sum(dim='nCells')
                    if outName!='ohc':
                        timeSeries = timeSeries / localLayerVol.sum(dim='nCells')

                dsOut[outName] = factor*timeSeries
                dsOut[outName].attrs['units'] = units
                dsOut[outName].attrs['description'] = description

            dsOut['totalArea'] = totalArea
            dsOut.totalArea.attrs['units'] = 'm^2'

            datasets.append(dsOut)

        # Combine data sets into a single data set for all regions
        dsOut = xarray.concat(datasets, 'nRegions')
        dsOut['refBottomDepth'] = refBottomDepth

        write_netcdf(dsOut, timeSeriesFile)
    else:
        print('Time series file already exists for year {}. Skipping it...'.format(year))

# Make plot
timeSeriesFiles = []
for year in years:
    timeSeriesFile = '{}_year{:04d}.nc'.format(timeSeriesFile0, year)
    timeSeriesFiles.append(timeSeriesFile)

if os.path.exists(featureFile):
    fcAll = read_feature_collection(featureFile)
else:
    raise IOError('No feature file found')

for regionIndex, regionName in enumerate(regionNames):
    print('    region: {}'.format(regionName))
    fc = FeatureCollection()
    for feature in fcAll.features:
        if feature['properties']['name'] == regionName:
            fc.add_feature(feature)
            break

    dsIn = xarray.open_mfdataset(timeSeriesFiles, combine='nested',
                                 concat_dim='Time', decode_times=False).isel(nRegions=regionIndex)

    movingAverageMonths = 12

    depths = dsIn.refBottomDepth.values[0]
    z = np.zeros(depths.shape)
    z[0] = -0.5 * depths[0]
    z[1:] = -0.5 * (depths[0:-1] + depths[1:])

    Time = dsIn.Time.values

    for var in variables:
        varName = var['name']

        clevels = var['clevels']
        colormap0 = var['colormap']
        colorIndices0 = var['colorIndices']
        underColor = colormap0(colorIndices0[0])
        overColor = colormap0(colorIndices0[-1])
        if len(clevels)+1 == len(colorIndices0):
            # we have 2 extra values for the under/over so make the colormap
            # without these values
            colorIndices = colorIndices0[1:-1]
        else:
            colorIndices = colorIndices0
        colormap = cols.ListedColormap(colormap0(colorIndices))
        colormap.set_under(underColor)
        colormap.set_over(overColor)
        cnorm = mpl.colors.BoundaryNorm(clevels, colormap.N)

        field = dsIn[varName]

        # Compute first-year average (note that this assumes monthly fields)
        fieldMean = field.isel(Time=range(12)).mean(dim='Time')

        # Compute moving average of the anomaly with respect to first-year average
        if movingAverageMonths != 1:
            N = movingAverageMonths
            movingAverageDepthSlices = []
            for nVertLevel in range(nVertLevels):
                depthSlice = field.isel(nVertLevels=nVertLevel) - fieldMean.isel(nVertLevels=nVertLevel)
                mean = pd.Series.rolling(depthSlice.to_series(), N,
                                         center=True).mean()
                mean = xarray.DataArray.from_series(mean)
                mean = mean[int(N / 2.0):-int(round(N / 2.0) - 1)]
                movingAverageDepthSlices.append(mean)
            field = xarray.DataArray(movingAverageDepthSlices)

        xLabel = 'Time (yr)'
        yLabel = '{} ({})'.format(var['title'], var['units'])
        title = '{} Anomaly, {}'.format(var['title'], regionName)
        figFileName = '{}/{}vsTimeDepth_{}.png'.format(figdir, varName,
                regionName[0].lower()+regionName[1:].replace(' ', ''))

        fig = hovmoeller_plot(Time[N-1:], z, field.values, colormap, cnorm, clevels,
                              title, xLabel, yLabel, calendar, colorbarLabel=var['units'],
                              titleFontSize=None, figsize=(15, 6), dpi=None)

        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        if regionName!='Global':
            add_inset(fig, fc, width=1.5, height=1.5, xbuffer=0.5, ybuffer=-1)

        plt.savefig(figFileName, dpi='figure', bbox_inches='tight',
                    pad_inches=0.1)
        plt.close()
