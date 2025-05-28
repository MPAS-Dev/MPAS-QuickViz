#!/usr/bin/env python
"""
    Name: build_netcdf.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Build netCDF outputs using xarray.
"""

from datetime import datetime

import numpy as np
import xarray as xr
import yaml

from binning import get_bins_edges


def build_output_2D(datestamp, data, coords, params):
    """Build 2D output
    """

    # Load definitions
    with open('../yaml/variable_metadata.yaml', 'r') as f:
        definitions = yaml.safe_load(f)['variables']
    
    # Build coordinates dict
    coordinates = {'time': [datestamp]}
    if 'depths' in params:
        coordinates['depth'] = np.array(params['depths'], dtype='float')
    for coordName in ['lonCell', 'latCell']:
        coordinates[coordName] = coords[coordName]
    
    # Build variables dict
    variables = {'areaCell': coords.areaCell}
    for varName, values in data.items():
        dims = _get_dims_2D(params['run']['model_name'], values.ndim)
        variables[varName] = (dims, values[None, ...])
    
    # Build dataset
    ds = xr.Dataset(variables, coordinates)

    # Variable attributes
    for varName in data:
        defn = definitions[varName]
        if 'total' in varName and 'ctgy' in defn:
            attr = _get_flux_component_attr(varName, definitions)
            ds[varName].attrs['source'] = attr
        ds[varName].attrs['units'] = defn['units']

    # Global attributes
    ds = _add_global_attrs(ds, params['run'])

    return ds


def build_output_binned(datestamp, data, coords, params):
    """Build binned output
    """

    # Load definitions
    with open('../yaml/variable_metadata.yaml', 'r') as f:
        definitions = yaml.safe_load(f)['variables']
    
    # Build coordinates dict
    coordinates = {'time': [datestamp], 'regionNames': coords.regionNames}
    for name, args in params['binning'].items():
        if any(name in key.lower() for key in data):
            coordinates[f'{name}Bins'] = get_bins_edges(*args)[0]
    
    # Build variables dict
    variables = {}
    for varName, values in data.items():
        dims = _get_dims_binned(varName, params['binning'])
        variables[varName] = (dims, values[None, ...])
    
    # Build dataset
    ds = xr.Dataset(variables, coordinates)

    # Variable attributes
    for varName in data:
        if 'Transport' in varName:
            units = 'Sv'
        elif 'Volume' in varName:
            units = 'km3'
        else:
            msg = f'Variable name {varName} does not include Transport or Volume.'
            raise NameError(msg)
        ds[varName].attrs['units'] = units

    # Global attributes
    ds = _add_global_attrs(ds, params['run'])

    return ds


def build_output_transect(datestamp, data, coords, params):
    """Build binned output
    """

    # Load definitions
    with open('../yaml/variable_metadata.yaml', 'r') as f:
        definitions = yaml.safe_load(f)['variables']

    # Model name
    model = params['run']['model_name']
    
    # Initialize coordinates and variables dicts
    coordinates = {
        'time': [datestamp],
        'refBottomDepth': coords.refBottomDepth.values,
    }
    variables = {}

    # Build transect variables
    if 'transect' in data:

        # Add transectNames to coordinates
        coordinates['transectNames'] = coords.transectNames

        # Get relevant coordinate variables
        names = ['lon', 'lat', 'dv', 'dist', 'sign']
        if model == 'MPAS':
            names.append('edge')
        elif model == 'POP':
            names.extend(['j', 'i'])
        for name in names:
            key = f'{name}Transect'
            variables[key] = coords[key]

        # Get transect quanities from data
        dims = ['time', 'transectNames', 'nTransectEdges', 'refBottomDepth']
        for varName, values in data['transect'].items():
            variables[varName] = (dims, values[None, ...])

    # Build meridional variables
    if 'meridional' in data:
        coordinates['latMOC'] = coords.latMOC.values
        dims = ['time', 'latMOC', 'refBottomDepth']
        for varName, values in data['meridional'].items():
            variables[varName] = (dims, values[None, ...])
    
    # Build dataset
    ds = xr.Dataset(variables, coordinates)

    # Variable attributes
    for varName in data['transect'].keys() | data['meridional'].keys():
        defn = definitions[varName]
        ds[varName].attrs['units'] = defn['units']

    # Global attributes
    ds = _add_global_attrs(ds, params['run'])

    return ds


def _get_dims_2D(model, ndim):
    """Get dimensions for 2D variable
    """

    # Model-specific space dimensions
    if model == 'MPAS':
        space_dims, space_ndim = ['nCells'], 1
    elif model == 'POP':
        space_dims, space_ndim = ['nlat', 'nlon'], 2

    # Determine dims
    if ndim > space_ndim:
        dims = ['time', 'depth'] + space_dims
    else:
        dims = ['time'] + space_dims

    return dims


def _get_dims_binned(varName, params):
    """Get dimensions for binned variable
    """

    # Get dimensions
    for name in params:
        if name in varName.lower():
            dims = ['time', 'regionNames', f'{name}Bins']
            break
    else:
        if 'TS' in varName:
            dims = ['time', 'regionNames', 'temperatureBins', 'salinityBins']
        else:
            msg = f'Variable name {varName} does not conform to bin naming conventions.'
            raise NameError(msg)

    return dims


def _add_global_attrs(ds, params):
    """Add global attributes
    """
    
    # Global attributes
    for key in ['model', 'short', 'long', 'mesh']:
        ds.attrs[f'{key}_name'] = params[f'{key}_name']
    ds.attrs['date_created'] = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

    return ds


def _get_flux_component_attr(varName, definitions):
    """Get flux component attributes
    """

    # Search definitions for shared ctgy to varName
    attr = []
    for name, defn in definitions.items():
        conditions = (
            'ctgy' in defn and
            defn['ctgy'] == definitions[varName]['ctgy'] and
            name != varName
        )
        if conditions:
            attr.append(name)
    attr = ' + '.join(attr)

    return attr
