#!/usr/bin/env python
"""
    Name: calc_transformation_budget.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to calculate the SPNA transformation budget
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import xarray as xr
import pandas as pd
import yaml
import time
from datetime import datetime
import postprocesstools as pptools
import watermasstools as wmttools


def interpolate_to_edge(varCell, cellsOnEdge, subdomain):
    """Interpolate cell variable to edge
    """
    
    # Find cells on edge in subdomain
    index = list(np.array([np.where(np.isin(subdomain, pair))[0] for pair in cellsOnEdge]).T)

    # Interpolate to edge
    varEdge = (varCell[index[0], :] + varCell[index[1], :]) / 2
    
    return varEdge


def get_region_edges(regionName, regionMask, cellsOnEdge, lonEdge, latEdge):
    """Get open water edges of an MPAS-Ocean region.
    Code adapted from A. Barthel
    """
    
    # In/Out transect cutoffs (lon, lat)
    # Hardcoded for each region
    cutoffs = {
        'Labrador Sea'            : (-44.0, 60.4),
        'Irminger Iceland Rockall': ( -5.9, 60.3),
        'Nordic Seas'             : ( -3.5, 69.0),
        'Atlantic Basin'          : ( 30.0,  0.0),
    }
    
    # Find land edges over entire mesh
    landEdges = np.any(cellsOnEdge == -1, axis=1)
    
    # Find land edges over region
    getRidEdges = landEdges
    getRidEdges[~landEdges] = np.equal(*[regionMask[cellsOnEdge[~landEdges, col]] for col in (0, 1)])
    
    # Find open boundary edges and signs (positive INTO region)
    openBoundaryEdges, = np.where(~getRidEdges)
    openBoundarySigns = np.sign((regionMask[cellsOnEdge[~getRidEdges, 1]] - 0.5))
    
    # Determine which transects represent the poleward end of the region
    cutoff = cutoffs[regionName]
    polewardIndex = (lonEdge[openBoundaryEdges] > cutoff[0]) | (latEdge[openBoundaryEdges] > cutoff[1])
    
    return openBoundaryEdges, openBoundarySigns, polewardIndex


def get_transect_masks_from_regions(regionNames, regionMasks, cellsOnEdge, lonEdge, latEdge):
    """Hardcoded function to get transect masks from region masks
    """
    
    # Get transects
    transectMasks = {}
    transectPairs = [('Davis+Hudson', 'OSNAP West'), ('Iceland-Scotland', 'OSNAP East'), ('Fram+Barents+NS',)]
    for names, region, mask in zip(transectPairs, regionNames, regionMasks.T):
        edges, signs, ipoleward = get_region_edges(region, mask, cellsOnEdge, lonEdge, latEdge)
        for name, func, signchange in zip(names, ['array', 'invert'], [-1, 1]):
            index = getattr(np, func)(ipoleward)
            transectMasks[name] = {
                'index': edges[index],
                'sign': signchange * signs[index],
            }
    
    return transectMasks


def calc_volumetric_TS(temperature, salinity, volume, coords, bins):
    """Calculate volumetric TS for MPAS results. Uses `numpy.histogramdd` to
    create the volume-weighted 2D TS histogram.
    """

    # Loop through regions
    volumeTS = []
    for regionMask in coords['regionMasks'].T:
        S, T, V = [var[regionMask, :].ravel() for var in (salinity, temperature, volume)]
        vol, _ = np.histogramdd((T, S), weights=V, bins=bins)
        volumeTS.append(vol) # km3
    
    # Xarray output
    coordinates = {'regionNames': coords['regionNames'], 'temperatureBins': bins[0][:-1], 'salinityBins': bins[1][:-1]}
    variables = {'volumetricTS': (coordinates.keys(), np.array(volumeTS) * 1e-9)}
    ds_out = xr.Dataset(variables, coordinates)
    
    return ds_out


def calc_volume(volume, sigmaTheta, coords, regions, sigmarange=[27, 28.5], binsize=0.01):
    """Calculate time-averaged overturning transformation over specified sigma bins
    and return in sigma space as an `xr.Dataset`
    """
    
    # Build sigma variables
    sigmabins = wmttools.build_sigma_bins(sigmarange, binsize)
    
    # Initialize variables dict
    volumeSigma = []

    # Loop through sigmabins
    for sigmabin in sigmabins:
        
        sigmaMask = (sigmaTheta >= sigmabin) & (sigmaTheta <= sigmabin + binsize)
        
        volumeRegions = []
        for region in regions:
            index, = np.where(coords['regionNames'] == region)
            regionMask = coords['regionMasks'][:, index[0]][:, None] & sigmaMask
            volumeRegions.append(volume[regionMask].sum())
        volumeSigma.append(volumeRegions)

    # Xarray output
    coordinates = {'sigmaBins': sigmabins, 'regionNames': regions}
    variables = {'volume': (coordinates.keys(), np.array(volumeSigma) * 1e-9)}
    ds_out = xr.Dataset(variables, coordinates)
    
    return ds_out


def calc_overturning(transectVars, sigmarange=[27, 28.5], binsize=0.01):
    """Calculate time-averaged overturning transformation over specified sigma bins
    and return in sigma space as an `xr.Dataset`
    """
    
    # Build sigma variables
    sigmabins = wmttools.build_sigma_bins(sigmarange, binsize)
    
    # Initialize variables dict
    overturning = []

    # Loop through sigmabins
    for sigmabin in sigmabins:
        
        transformation = []
        for transect in transectVars:
            sigmaTheta, transport = [transectVars[transect][name] for name in ('sigmaTheta', 'transport')]
            sigmaMask = (sigmaTheta >= sigmabin) & (sigmaTheta <= sigmabin + binsize)
            transformation.append(transport[sigmaMask].sum())
        overturning.append(transformation)

    # Xarray output
    coordinates = {'sigmaBins': sigmabins, 'transectNames': list(transectVars.keys())}
    variables = {'overturningTransformation': (coordinates.keys(), np.array(overturning)[::-1, :].cumsum(axis=0)[::-1, :] * 1e-6)}
    ds_out = xr.Dataset(variables, coordinates)
    
    return ds_out


def load_coords(paths, bbox=[-70, 23, 44, 80]):
    """
    """

    # Load mask variables
    with xr.open_dataset(paths['maskfile']) as ds:
        regionNames = ds.regionNames.values.astype(str)
        regionMasks = ds.regionCellMasks.values

    # Load mesh variables
    with xr.open_dataset(paths['meshfile']) as ds:
        lonCell = np.rad2deg(ds.lonCell.values)
        latCell = np.rad2deg(ds.latCell.values)
        lonEdge = np.rad2deg(ds.lonEdge.values)
        latEdge = np.rad2deg(ds.latEdge.values)
        lonVertex = np.rad2deg(ds.lonVertex.values)
        latVertex = np.rad2deg(ds.latVertex.values)
        cellsOnEdge = ds.cellsOnEdge.values - 1
        verticesOnEdge = ds.verticesOnEdge.values - 1
        dvEdge = ds.dvEdge.values
        area = ds.areaCell.values
        depth = ds.refBottomDepth.values

    # Correct lons
    lonCell = np.where(lonCell > 180, lonCell - 360, lonCell)
    lonEdge = np.where(lonEdge > 180, lonEdge - 360, lonEdge)
    lonVertex = np.where(lonVertex > 180, lonVertex - 360, lonVertex)

    # Get transect masks
    transectMasks = get_transect_masks_from_regions(regionNames, regionMasks, cellsOnEdge, lonEdge, latEdge)

    # Build coords dict
    subdomain, = np.where((lonCell > bbox[0]) & (lonCell < bbox[1]) & (latCell > bbox[2]) & (latCell < bbox[3]))
    coords = {
        'regionNames': regionNames,
        'regionMasks': regionMasks[subdomain, :].astype(bool),
        'area': area[subdomain],
        'cellsOnEdge': cellsOnEdge,
        'verticesOnEdge': verticesOnEdge,
        'dvEdge': dvEdge,
        'lonEdge': lonEdge,
        'lonVertex': lonVertex,
        'depth': depth,
    }
    
    return coords, transectMasks, subdomain


def load_transect_vars(ds, sigmaTheta, layerThickness, transectMasks, coords, subdomain):
    """
    """
    
    # Get transport variables
    variables = {}
    for transect in transectMasks:

        # Get transect edge indices and signs
        index, sign = [transectMasks[transect][name] for name in ('index', 'sign')]

        # Get normal Velocity
        velocity = ds.timeMonthly_avg_normalVelocity[0, index, :].values

        # Get density and layer thickness and interpolate to edge
        cellsOnTransect = coords['cellsOnEdge'][index]
        sigmaThetaTransect = interpolate_to_edge(sigmaTheta, cellsOnTransect, subdomain)
        layerThicknessTransect = interpolate_to_edge(layerThickness, cellsOnTransect, subdomain)

        # Calculate overturning variables
        variables[transect] = {
            'sigmaTheta': sigmaThetaTransect,
            'transport': velocity * layerThicknessTransect * (coords['dvEdge'][index] * sign)[:, None],
        }
    
    return variables


def main_routine():
    """
    """
    
    # Save path
    savepath = '/pscratch/sd/b/bmoorema/results/aggregated/transformation/'

    # Loop though meshes
    for mesh in ['LR', 'HR']:
        
        # Get paths
        with open(f'../yaml/paths_{mesh}.yaml') as f:
            paths = yaml.safe_load(f)
        
        # Load coords
        coords, transectMasks, subdomain = load_coords(paths)

        # Loop through decades
        startyear = 1946
        for decade in [(1947, 1957), (1997, 2007)]:
            
            # Define results path
            decade_str = str(decade[0]) + '-' + str(decade[1])
            resultspath = paths['results'][decade_str] + '/' + paths['prefix']
            
            # Define times
            times = [datetime(year, month, 1) for year in range(*decade) for month in range(1, 13)]
            
            # Loop through times
            print(f'Loading {mesh}, decade {decade_str}...')
            wmt, n, starttime = [], len(times), time.time()
            for k, t in enumerate(times):

                # Load results
                filename = resultspath + f'.{t.year-startyear:04d}-{t.month:02d}-01.nc'
                with xr.open_dataset(filename) as ds:

                    # Get cell variables
                    prefix = 'timeMonthly_avg_'
                    sigmaTheta = ds[prefix + 'potentialDensity'][0, ...].values[subdomain, :] - 1000
                    layerThickness = ds[prefix + 'layerThickness'][0, ...].values[subdomain, :]
                    temperature = ds[prefix + 'activeTracers_temperature'][0, ...].values[subdomain, :]
                    salinity = ds[prefix + 'activeTracers_salinity'][0, ...].values[subdomain, :]

                    # Load transect variables
                    transectVars = load_transect_vars(ds, sigmaTheta, layerThickness, transectMasks, coords, subdomain)

                # Calculate volume
                volume = layerThickness * coords['area'][:, None]

                # Get state variables and buoyancy fluxes
                sigmaSurface, heatFactor, saltFactor = wmttools.calc_state_variables(salinity[:, 0], temperature[:, 0])
                fluxes = wmttools.build_combined_fluxes(ds, heatFactor, saltFactor, subdomain=subdomain)

                # Calculate 1D water mass transformation over regions
                ds_out = wmttools.calc_wmt(
                    fluxes, sigmaSurface, coords, regions=coords['regionNames'], sigmarange=[27, 28.5], binsize=0.01,
                )

                # Calculate overturning transformation over transects
                ds_out = xr.merge([ds_out, calc_overturning(transectVars)])

                # Calculate volume
                ds_out = xr.merge([ds_out, calc_volume(volume, sigmaTheta, coords, coords['regionNames'])])
                
                # Calculate volumetric TS
                bins = [np.arange(-3, 20.1, 0.1), np.arange(33, 37.01, 0.01)]
                ds_out = xr.merge([ds_out, calc_volumetric_TS(temperature, salinity, volume, coords, bins)])

                wmt.append(ds_out)
                
                # Update loop
                pptools.loopstatus(k, n, starttime)

            # Concatenate time
            wmt = xr.concat(wmt, pd.Index(times, name='time'))
            wmt.to_netcdf(savepath + f'transformationbudget_{mesh}_{decade_str}.nc')


if __name__ == "__main__":
    main_routine()
