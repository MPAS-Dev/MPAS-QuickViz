#!/usr/bin/env python
"""
    Name: binning.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Binning module
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np

from process_velocities import calc_transport


def get_bins_edges(start, stop, step):
    """Get bins and edges
    """
    
    # Get bins and edges
    bins = np.arange((stop - start) / step + 2) * step + start
    edges = bins - step / 2
    
    return bins[:-1], edges


def bin_by_region(data, coords, params):
    """Bin edge transport and cell volume by region
    """

    # Get boundary transport from normal velocities
    transports = calc_transport(data['region'], coords, 'region')

    # Get regionMasks and move region index to 0
    regionMasks = np.moveaxis(coords.regionCellMasks.values, -1, 0)

    # Bin by T,S,rho and TS
    data['binned'] = {}
    data = _bin_by_region_1D(data, transports, regionMasks, params)
    data = _bin_by_region_TS(data, transports, regionMasks, params)
    
    return data


def _bin_by_region_1D(data, transports, regionMasks, params):
    """Bin edge transport and cell volume by temperature, salinity, and density (1D)
    """

    # Get volume, and regionCellMasks
    volume = data['3D']['volume']
    
    # Iterate through bin dims
    for dimName in params:
    
        # Get arrays and edges
        normalize = 1000 if dimName == 'density' else 0
        arraysEdge = data['region'][dimName] - normalize
        arraysCell = data['3D'][dimName] - normalize
        _, binEdges = get_bins_edges(*params[dimName])
    
        # Bin by region
        H = {'Transport': [], 'Volume': []}
        iterables = zip(arraysEdge, transports, regionMasks)
        for arrayEdge, transport, regionMask in iterables:
            
            # Boundary transport
            h, _ = np.histogram(
                arrayEdge.ravel(),
                bins=binEdges, weights=transport.ravel(),
            )
            H['Transport'].append(h)
            
            # Volume
            h, _ = np.histogram(
                arraysCell[:, regionMask].ravel(),
                bins=binEdges, weights=volume[:, regionMask].ravel(),
            )
            H['Volume'].append(h)
    
        # Aggregate binned arrays
        dim = dimName.capitalize()
        for key, values in H.items():
            data['binned'][f'binned{key}{dim}'] = np.array(values)

    return data


def _bin_by_region_TS(data, transports, regionMasks, params):
    """Bin edge transport and cell volume by temperature-salinity (2D)
    """

    # Get transport, volume, and regionCellMasks
    volume = data['3D']['volume']
    
    # Get arrays and edges
    binEdges, arraysEdge, arraysCell = [], [], []
    for dimName in ['temperature', 'salinity']:
        arraysEdge.append(data['region'][dimName])
        arraysCell.append(data['3D'][dimName])
        binEdges.append(get_bins_edges(*params[dimName])[1])

    # Bin by region
    H = {'Transport': [], 'Volume': []}
    iterables = zip(*arraysEdge, transports, regionMasks)
    for arr1, arr2, transport, regionMask in iterables:

        # Boundary transport
        h, _, _ = np.histogram2d(
            *[arr.ravel() for arr in (arr1, arr2)],
            bins=binEdges, weights=transport.ravel(),
        )
        H['Transport'].append(h)

        # Volume
        h, _, _ = np.histogram2d(
            *[arr[:, regionMask].ravel() for arr in arraysCell],
            bins=binEdges, weights=volume[:, regionMask].ravel(),
        )
        H['Volume'].append(h)

    # Aggregate binned arrays
    for key, values in H.items():
        data['binned'][f'binned{key}TS'] = np.array(values)

    return data
