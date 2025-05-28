#!/usr/bin/env python
"""
    Name: load_variables.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Module to load variables
    for the ImPACTS Water Mass Analysis project.
"""

import os

import numpy as np
import xarray as xr
import yaml

import transects_MPAS
import transects_POP


def build_results_filename(params, datecode='0001-01', ftype='ocean'):
    """Build results filename from datecode
    """

    # Get params for building the filename structure
    keys = ['model_name', 'long_name', 'prefix', 'results_path']
    model, long_name, prefix, results_path = [params[key] for key in keys]

    # Parse multiple run directories
    if not isinstance(results_path, str):
        run_year = int(datecode[:4])
        for run_dir in results_path[1]:
            years = [int(year) for year in run_dir.split('_')[-2:]]
            if years[0] <= run_year <= years[1]:
                results_path = os.path.join(results_path[0], run_dir)
                break
        else:
            msg = (
                f'Date {datecode} not available in results record for simulation {long_name}'
            )
            raise ValueError(msg)

    # Build filename
    if model == 'MPAS':
        strings = [long_name, params[f'{ftype}_name'], prefix, f'{datecode}-01', 'nc']
    elif model == 'POP':
        strings = [long_name, prefix, datecode, 'nc']
    else:
        raise ValueError(f'Invalid model identifier: {model}')
    filename = os.path.join(results_path, '.'.join(strings))
    
    return filename


def load_variables(params, coords, datecode='0001-01', prefix='timeMonthly_avg_'):
    """Load variables and return appended lists
    """

    # Load YAML variable definitions
    with open('../yaml/variable_metadata.yaml', 'r') as f:
        definitions = yaml.safe_load(f)['variables']

    # Get model name
    run_vars = params['run']
    model = run_vars['model_name']
    if model == 'MPAS':
        space_ndim = 1
    elif model == 'POP':
        space_ndim = 2

    # Load results file
    filename = build_results_filename(run_vars, datecode)
    datasets = {'ocean': xr.open_dataset(filename)}

    # Parse simulation parameters
    data = _init_data(params)

    # Get depth-averaging parameters, get layerThickness edge fields if MPAS
    depthVars = _get_depth_avg_params(coords, params['depths'], datasets['ocean'])
    if model == 'MPAS' and 'masking' in params:
        data = _get_edge_fields(data, depthVars[0].T, coords, 'layerThickness')
    if 'binning' in params:
        data['3D']['volume'] = depthVars[0] * coords.areaCell.values[None, :] * 1e-9
    
    # Loop through variable names
    for varName in params['variables']:

        # Determine variable source file, and load if not already loaded
        defn = definitions[varName]
        ftype = defn['type'] if 'type' in defn else 'ocean'
        if ftype not in datasets:
            filename = build_results_filename(run_vars, datecode, ftype)
            datasets[ftype] = xr.open_dataset(filename)

        # Check that variable exists in results file
        defn = defn[model]
        varNameFull = prefix + defn['name']
        if varNameFull not in datasets[ftype]:
            raise KeyError(f'Variable {varNameFull} not found in results file type {ftype}')

        # Load into numpy, apply subdomain if indexed on nCells, apply conversion
        array = datasets[ftype][varNameFull][0, ...]
        dims = array.dims
        array = array.values
        if 'nCells' in dims:
            array = array[coords.nCells.values, ...]
        if 'conversion' in defn:
            array = array * defn['conversion']

        # Get AMOC and bypass remaining (hardcoded for Atlantic and refBottomDepth)
        if 'moc' in varName:
            if model == 'MPAS':
                array = array[0, ...].T
            elif model == 'POP':
                array = array[1, :, 1:, :].sum(axis=0).T
            data['meridional'][varName] = array
            continue
    
        # 3D cell variables
        if array.ndim > space_ndim:
    
            # Get edge fields
            if 'masking' in params:
                condition = model == 'MPAS' and ('Zon' in varName or 'Mer' in varName)
                if not condition:
                    data = _get_edge_fields(data, array, coords, varName)
    
            # MPAS: bypass remaining if edge3D, otherwise reorder to depth axis=0
            if model == 'MPAS':
                if 'Normal' in varName:
                    continue
                else:
                    array = array.T

            # Keep 3D variables for binning
            if varName in ['temperature', 'salinity', 'density']:
                if 'binning' in params:
                    data['3D'][varName] = array
    
            # Get depth-weighted averages
            array = _calc_depth_averages(array, *depthVars)
    
        # Append to list
        data['2D'][varName] = array
    
    return data


def add_surface_fluxes(data, params, discard_components=True):
    """Add surface fluxes
    """

    # Load YAML variable definitions
    with open('../yaml/variable_metadata.yaml', 'r') as f:
        definitions = yaml.safe_load(f)['variables']
    
    # Initialized totals
    ctgys = []
    for ctgy in ['Heat', 'Freshwater', 'Salt', 'Salinity']:
        key = f"total{ctgy}Flux"
        if key not in data:
            data[f"total{ctgy}Flux"] = 0
            ctgys.append(ctgy)
    
    # Iterate through variables
    for varName in params:
        defn = definitions[varName]
        if 'ctgy' in defn:
            ctgy = defn['ctgy'].capitalize()
            key = f'total{ctgy}Flux'
            if ctgy in ctgys and varName != key:
                data[key] = data[key] + data[varName]
                if discard_components:
                    _ = data.pop(varName)

    return data


def _init_data(params):
    """Initialize data dictionary
    """

    # Initialize data dictionary
    data = {'2D': {}}
    if 'binning' in params:
        data['3D'] = {}
    if any('moc' in var for var in params['variables']):
        data['meridional'] = {}
    if 'masking' in params:
        data['region'] = {}
        data['transect'] = {}

    return data


def _get_depth_avg_params(coords, depths, ds=None, prefix='timeMonthly_avg_'):
    """Get depth averaging parameters
    """

    # Get layerThickness
    model = coords.model_name
    if model == 'MPAS':
        layerThickness = ds[f'{prefix}layerThickness'][0, ...].values
        layerThickness = layerThickness[coords.nCells, :].T
    elif model == 'POP':
        layerThickness = coords.layerThickness.values[:, None, None]

    # Get bottom depth and vertical layer indexes
    bottomDepth = np.cumsum(layerThickness, axis=0)
    refBottomDepth = coords.refBottomDepth.values
    kindex = [abs(refBottomDepth - depth).argmin() for depth in depths]

    return layerThickness, bottomDepth, kindex


def _get_edge_fields(data, array, coords, varName):
    """Get edge fields
    """

    # Determine velocity boolean
    velocity = True if 'velocity' in varName else False

    # Get edge fields for region and transect
    model = coords.model_name
    for ctgy in ['region', 'transect']:
        if model == 'MPAS':
            data[ctgy][varName] = transects_MPAS.get_edge_fields(
                ctgy, array, coords, velocity=velocity,
            )
        elif model == 'POP':
            data[ctgy][varName] = transects_POP.get_edge_fields(
                ctgy, array, coords, velocity=velocity,
            )

    return data


def _calc_depth_averages(array, layerThickness, bottomDepth, kindex):
    """Calculate depth-weighted averages
    """

    # Calculate depth-weighted averages
    depthAvg = array * layerThickness
    depthAvg = [np.nansum(depthAvg[:k+1, ...], axis=0) / bottomDepth[k, :] for k in kindex]
    
    return np.array(depthAvg)
