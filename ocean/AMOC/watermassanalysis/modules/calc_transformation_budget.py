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
import os
from datetime import datetime
import postprocesstools as pptools
import watermasstools as wmttools


def get_time_params(years, startyear=1948):
    """Get time parameters based on requested simulation year range (e.g., [1, 11])
    """
    
    # Get decade string
    decade = '-'.join([str(startyear + year - 1) for year in years])

    # Get years and months arrays
    years, months = np.arange(*years), np.arange(1, 13)
    
    # Get time index
    timeindex = [datetime(startyear + year - 1, month, 1) for year in years for month in months]
    timeindex = pd.Index(timeindex, name='time')
    
    return years, months, decade, timeindex


def stack_arrays(arrays):
    """Stack list of unequal length arrays into 2D array, padding with NaNs
    """
    
    # Stack arrays
    arrays_stacked = np.stack(list(zip_longest(*arrays, fillvalue=np.nan)), axis=1)
    
    return arrays_stacked


def build_transect_coords(transectMasks, nVertLevels):
    """Get and stack transect coordinates and pad with NaNs for xarray output
    """
    
    # Convert transectMasks to list
    masks = list(transectMasks.values())
    
    # Loop through fields in a transect mask dict
    transectCoords = {}
    for name in masks[0].keys():
        
        # Stack transects for given field into list
        data = [mask[name] for mask in masks]
        
        # Preserve edges and signs stacked lists
        if name == 'edges':
            edges = data
        elif name == 'signs':
            signs = data
        
        # Pad the stack and convert to numpy array
        transectCoords[name] = stack_arrays(data)
    
    return transectCoords, edges, signs


def get_transects(array, edges, signs=None, cellsOnEdge=None, subdomain=None):
    """Get transects given a variable array and a list of edge ID arrays
    """
    
    # Edge velocity -> get signed values on transects
    if signs is not None:
        transects = [sign[:, None] * array[edge, :] for edge, sign in zip(edges, signs)]
    
    # Cell variable -> interpolate to transects
    else:
        transects = [pptools.interpolate_to_edge(array, cellsOnEdge[edge], subdomain) for edge in edges]
    
    return stack_arrays(transects)


def append_to_list(array, varName, data):
    """Append an array to a list called varName in the data dict
    """
    
    # Check if varName in data, then append array to data[varName]
    if varName not in data:
        data[varName] = []
    data[varName].append(array)
    
    return data


def load_variables(ds, data, varNames, edges, signs, cellsOnEdge, subdomain):
    """Load variables and return appended lists
    """
    
    # Loop through varTypes (cell2D, cell3D, edge3D)
    for varType in varNames:
    
        # Loop through varNames
        for varName in varNames[varType]:

            # Get varNameFull
            tag = 'activeTracers_' if varName in ['temperature', 'salinity'] else ''
            varNameFull = 'timeMonthly_avg_' + tag + varName

            # Skip variable if not in dataset (e.g. GM, salinityRestoring)
            if varNameFull not in ds:
                continue

            # Load into numpy
            array = ds[varNameFull][0, ...].values
            
            # 3D edge variables
            if varType == 'edge3D':
                
                # Get signed edge velocities on transects, stack, and append to list
                transects = get_transects(array, edges, signs=signs)
                data['transect'] = append_to_list(transects, varName, data['transect'])
            
            # Cell variables
            else:
                
                # Get array on subdomain
                array = [subdomain, ...]

                # 3D cell variables
                if varType == 'cell3D':

                    # Interpolate to transects, stack, and append to list
                    transects = get_transects(array, edges, cellsOnEdge=cellsOnEdge, subdomain=subdomain)
                    data['transect'] = append_to_list(transects, varName, data['transect'])

                    # Keep volume variables
                    if varName in ['temperature', 'salinity', 'layerThickness']:
                        data['3D'][varName] = array

                    # Get depth averages
                    array = np.array([array[:, :zindex].mean(axis=1) for zindex in zindexes])

                # Append to list
                data['2D'] = append_to_list(array, varName, data['2D'])
    
    return data


def bin_variables(
    data, fluxes, coords, transectCoords,
    temperatureBins=(-2, 13, 0.1),
    salinityBins=(33.5, 35.7, 0.01),
    densityBins=(26.5, 28.2, 0.01),
):
    """Bin area, volume, transport and WMT by density and TS
    """
    
    # Calculate area, volume, and transport
    area = coords['areaCell'] * 1e-6  # ----------------------------- km2
    volume = data['3D']['layerThickness'] * area[:, None] * 1e-3  # - km3
    transport = (
        data['transect']['normalVelocity'][-1] *
        data['transect']['layerThickness'][-1] *
        transectCoords['dvEdge'][..., None] * 1e-6  # --------------- Sv
    )

    # Get transport binned by density
    X = data['transect']['potentialDensity'][-1] - 1000
    _, binEdges = wmttools.get_bins_edges(*densityBins)
    h = [np.histogram(x.ravel(), bins=binEdges, weights=w.ravel())[0] for x, w in zip(X, transport)]
    data['binned']['densityTransport'].append(np.array(h))

    # Get WMT binned by density
    h = wmttools.calc_wmt(fluxes, coords, 'density', densityBins, regions=coords['regionNames'])
    data['WMT']['density'].append(h)

    # Define TS bins
    binNames, binArgs = ('temperature', 'salinity'), [temperatureBins, salinityBins]
    _, binEdges = zip(*[wmttools.get_bins_edges(*arg) for arg in binArgs])

    # Get area and volume binned by TS
    h_area, h_volume = [], []
    for regionMask in coords['regionCellMasks'].T:
        X = [data['3D'][name][regionMask, :] for name in binNames]
        x, w = [x[:, 0] for x in X], area[regionMask]
        h_area.append(np.histogram2d(*x, bins=binEdges, weights=w))
        x, w = [x.ravel() for x in X], volume[regionMask, :].ravel()
        h_volume.append(np.histogram2d(*x, bins=binEdges, weights=w))
    data['binned']['TSArea'].append(np.array(h_area))
    data['binned']['TSVolume'].append(np.array(h_volume))

    # Get transport binned by TS
    X = [data['transect'][name][-1] for name in names]
    h = [np.histogram2d(x.ravel(), y.ravel(), bins=binEdges, weights=w)[0] for x, y, w in zip(*X, transport)]
    data['binned']['TSTransport'].append(np.array(h))

    # Get WMT binned by TS
    h = wmttools.calc_wmt(fluxes, coords, names, args, regions=coords['regionNames'])
    data['WMT']['TS'].append(h)

    # Get 2D WMT binned by all
    for name in ('temperature', 'salinity', 'density'):  
        h = wmttools.calc_wmt(fluxes, coords, name, binArgs[name], remapvars=remapvars)
        data['2D'][name].append(h)
    
    return data


def aggregate_to_xarray(data):
    """
    """
    
    datasets = {}

    # ----- Concatenate 2D fields -------------
    ds = []
    for varName in data['2D']:
        da = []
        for array in data['2D'][varName]:
            if varName in varNames3D:
                array_remap = [pptools.remap(layer, coords['nCells'], **remapvars) for layer in array]
                array_remap = xr.concat(array_remap, depthindex)
            else:
                array_remap = pptools.remap(array, coords['nCells'], **remapvars)
            da.append(array_remap)
        da = xr.concat(da, timeindex)
        da.name = varName
        ds.append(da)
    datasets['2D'] = xr.merge(ds)

    # Merge binned WMT 2D
    datasets['2D'] = xr.merge([datasets['2D'], xr.merge(xr.concat(values, timeindex) for values in data['WMT2D'].values())])

    # ----- Concatenate transect fields -------
    dims = ['time', 'transectNames', 'nEdges', 'nVertLevels']
    coordinates = {dims[0]: timeindex, dims[1]: list(transectMasks.keys())}
    coordinates = coordinates | {key: (dims[1:3], values) for key, values in coordsTransect.items()}
    variables = {key: (dims, np.array(values)) for key, values in data['transect'].items()}
    datasets['transect'] = xr.Dataset(variables, coordinates)

    # ----- Concatenate binned quantites ------
    coordinates = {'time': timeindex, 'regionNames': coords['regionNames'], 'transectNames': list(transectMasks.keys())}
    coordinates = coordinates | {key + 'Bins': wmttools.get_bins_edges(*args)[0] for key, args in binArgs.items()}
    dims = [
        ['time', 'transectNames', 'densityBins'],                      # densityTransport
        ['time', 'transectNames', 'temperatureBins', 'salinityBins'],  # TSTransport
        ['time', 'regionNames', 'temperatureBins', 'salinityBins'],    # TSArea
        ['time', 'regionNames', 'temperatureBins', 'salinityBins'],    # TSVolume
    ]
    variables = {name: (dim, np.array(data['binned'][name])) for dim, name in zip(dims, data['binned'])}
    datasets['WMT'] = xr.Dataset(variables, coordinates)

    # Merge binned WMT
    datasets['WMT'] = xr.merge([datasets['WMT'], xr.merge(xr.concat(values, timeindex) for values in data['WMT'].values())])
    
    return datasets


def main(
    meshName='LR',
    years=[1, 11],
    depths=[0, 20, 100, 500, 1000]
    bbox=[-100, 40, 40, 85],
    yamlpath='../yaml/',
):
    """
    """
    
    # Get paths
    with open(os.path.join(yamlpath, f'paths_{meshName}.yaml'), 'r') as f:
        paths = yaml.safe_load(f)

    # Load variable definitions
    with open(os.path.join(yamlpath, 'variable_definitions.yaml'), 'r') as f:
        varNames = yaml.safe_load(f)

    # Time parameters
    years, months, decade, timeindex = get_time_params(years)

    # Results prefix
    prefix = os.path.join(paths['results'][decade], paths['prefix'])
    
    # Remapping variables
    remapvars = pptools.build_remapper(paths['meshfile'], bbox=bbox)

    # Load coords
    coords, transectMasks, subdomain = pptools.load_coords(meshName)
    cellsOnEdge, refBottomDepth = [coords[name] for name in ('cellsOnEdge', 'refBottomDepth')]

    # Get and stack transect coordinates and pad with NaNs for xarray output
    transectCoords, edges, signs = build_transect_coords(transectMasks, len(refBottomDepth))
    
    # Depth index and corresponding zindexes for depth averaging (still need weights!)
    depthindex = pd.Index(depths, name='depth')
    zindexes = [abs(refBottomDepth - depth).argmin() + 1 for depth in depths]

    # Loop through files
    for year in years:
        for month in months:
            f = f'{prefix}.{year:04d}-{month:02d}-01.nc'
            with xr.open_dataset(f) as ds:

                # Load variables
                data = load_variables(ds, data, varNames, edges, signs, cellsOnEdge, subdomain)
                
                # Calculate surface fluxes
                fluxes = wmttools.build_fluxes(ds, subdomain=subdomain)
                
                # Bin variables
                data = bin_variables(data, fluxes, coords, transectCoords)

    # Aggregate to xarray
    aggregate_to_xarray(data)
    
    # Save to file
    savepath = '/pscratch/sd/b/bmoorema/results/aggregated/new/new/'
    for tag in datasets:
        datasets[tag].to_netcdf(savepath + meshName + '_' + tag + '.nc')


if __name__ == "__main__":
    main()
