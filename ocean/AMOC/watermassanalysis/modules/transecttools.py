#!/usr/bin/env python
"""
    Name: transecttools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for working with MPAS-Ocean transects
    for the ImPACTS Water Mass Analysis project.
"""

from itertools import zip_longest

import numpy as np
import xarray as xr


def get_transect_distance(lons, lats, d0=0, radius=6362):
    """Calculate the great circle (haversine) distance in kilometers
    between two points on the earth (specified in decimal degrees)
    """

    # Calculate cumulative haversine distance
    lons, lats = np.deg2rad(lons), np.deg2rad(lats)
    angle = np.arcsin(np.sqrt(
        np.sin(np.diff(lats)/2)**2
        + np.cos(lats[:-1])
        * np.cos(lats[1:])
        * np.sin(np.diff(lons)/2)**2
    ))
    distance = 2 * radius * angle
    distance = np.insert(np.cumsum(distance) + d0, 0, d0)

    return distance


def stack_arrays(arrays):
    """Stack list of unequal length arrays into 2D array, padding with NaNs
    """
    
    # Get fill value based on array shape
    shape = arrays[0].shape
    if len(shape) > 1:
        fillvalue = np.full(shape[1], np.nan)
    else:
        fillvalue = np.nan
    
    # Stack arrays
    arrays_stacked = np.stack(list(zip_longest(*arrays, fillvalue=fillvalue)), axis=1)
    
    return arrays_stacked


def build_transect_xarray(edgeMasks):
    """Stack transect coordinate variables and convert to xarray
    """

    # Loop through ctgys
    ds = {}
    for ctgy, masks in edgeMasks.items():

        # Build coordinates dict
        dims = [ctgy + 'Names', f'n{ctgy.capitalize()}Edges']
        coordinates = {dims[0]: list(masks.keys())}

        # Stack transects or region edges
        variables = {}
        for mask in masks.values():
            for key, values in mask.items():
                if key not in variables:
                    variables[key] = []
                variables[key].append(values)

        # Convert stacked transects to arrays and add xarray dims
        for key, values in variables.items():
            variables[key] = (dims, stack_arrays(values))

        # Create xr.Dataset and merge with coords
        ds[ctgy] = xr.Dataset(variables, coordinates)
    
    return ds
