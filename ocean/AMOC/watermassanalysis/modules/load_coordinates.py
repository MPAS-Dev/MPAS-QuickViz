#!/usr/bin/env python
"""
    Name: load_coordinates.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to calculate the SPNA transformation budget
    for the ImPACTS Water Mass Analysis project.
"""

import yaml

import numpy as np
import xarray as xr

import transects_MPAS
import transects_POP
from load_variables import build_results_filename


def load_coords(params):
    """Load model coordinate variables
    """
    
    # Load YAML params
    run_vars = params['run']

    # Load coordinate variables
    coords = _load_coords_from_file(params)

    # Load regionCellMasks and transect masks
    if 'masking' in params:
        coords['regionCellMasks'] = _load_mask_from_file(params['masking'])
        coords = _build_transect_masks_from_regions(coords, params)

    # Trim to bbox
    if run_vars['model_name'] == 'MPAS' and 'bbox' in run_vars:
        coords = _trim_to_bbox(coords, run_vars['bbox'])

    # Convert regionCellMasks to bool and merge coords and transects
    coords['regionCellMasks'] = coords.regionCellMasks.astype(bool)
    
    # Global attributes
    coords.attrs['model_name'] = run_vars['model_name']
    coords.attrs['short_name'] = run_vars['short_name']
    coords.attrs['mesh_name'] = run_vars['mesh_name']
    
    return coords


def _load_coords_from_file(params):
    """Load coordinate variables from file
    """

    # Load YAML variable definitions
    with open('../yaml/variable_metadata.yaml', 'r') as f:
        definitions = yaml.safe_load(f)['coordinates']

    # Load mesh dataset and initialize coords dict
    run_vars = params['run']
    datasets = {'mesh': xr.open_dataset(run_vars['mesh_file'])}
    coords = {}

    # Iterate through coordinate names
    for name in params['coordinates']:

        # Determine coordinate source file, and load if not already loaded
        defn = definitions[name]
        ftype = defn['type'] if 'type' in defn else 'mesh'
        if ftype not in datasets:
            filename = build_results_filename(run_vars, ftype=ftype)
            datasets[ftype] = xr.open_dataset(filename)

        # Load and postprocess coordinate array
        defn = defn[run_vars['model_name']]
        coord = datasets[ftype][defn['name']].values
        if 'conversion' in defn:
            conversion = defn['conversion']
            if conversion == 'deg':
                coord = np.rad2deg(coord)
            else:
                coord = coord * conversion
        if 'lon' in name:
            coord = np.where(coord > 180, coord - 360, coord)
        if 'On' in name:
            coord = coord - 1
        coords[name] = (defn['dims'], coord)

    # Convert to xarray Dataset and add attributes
    coords = xr.Dataset(coords)
    for name in params['coordinates']:
        coords[name].attrs['units'] = definitions[name]['units']
    
    return coords


def _load_mask_from_file(params):
    """Load region cell mask from file
    """
    
    # Load regionCellMasks from file
    ds = xr.open_dataset(params['mask_file'])
    masks = ds.regionCellMasks
    if 'regionNames' not in masks.dims:
        masks.coords['regionNames'] = ds.regionNames.astype(str)
        masks = masks.swap_dims({'nRegions': 'regionNames'})
    if 'regions' in params:
        masks = masks.sel(regionNames=params['regions'])

    return masks


def _build_transect_masks_from_regions(coords, params):
    """Build transect masks from regions
    """

    # Get transects
    model = params['run']['model_name']
    filename = params['masking']['region_transects_file']
    if model == 'MPAS':
        coords = transects_MPAS.get_transect_masks_from_regions(
            coords, filename,
        )
    elif model == 'POP':
        coords = transects_POP.get_transect_masks_from_regions(
            coords, filename,
        )
    else:
        raise ValueError(f'Invalid model identifier: {model}')

    return coords


def _trim_to_bbox(coords, bbox):
    """Trim to bbox
    """

    # Trim to bbox
    lons, lats = coords.lonCell, coords.latCell
    subdomain, = np.where(
        (lons > bbox[0]) & (lons < bbox[1]) &
        (lats > bbox[2]) & (lats < bbox[3])
    )
    coords = coords.isel(nCells=subdomain)
    
    return coords
