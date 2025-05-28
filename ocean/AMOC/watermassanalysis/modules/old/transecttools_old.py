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

import postprocesstools as pptools


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


def get_region_edges(regionCellMask, cellsOnEdge):
    """Get open water edges of an MPAS-Ocean masked region. Code adapted from Alice Barthel.
    Uses Python indexing, so subtract 1 from cellsOnEdge before calling this function.
    
    Inputs:
      regionCellMask: Region cell mask array (nCells,)
      cellsOnEdge: Cells on edge array, Python-indexed (nEdges, 2)
    
    Outputs:
      openBoundaryEdges: Edge indices (nBoundaryEdges,)
      openBoundarySigns: Edge signs (nBoundaryEdges,)
    """
    
    # Exclude land edges
    dropEdges = np.any(cellsOnEdge == -1, axis=1)
    
    # Exclude water edges that don't straddle the regionCellMask boundary (e.g. 0 neigboring 1)
    regionCellMaskNeighbors = [regionCellMask[cellsOnEdge[~dropEdges, col]] for col in (0, 1)]
    dropEdges[~dropEdges] = np.equal(*regionCellMaskNeighbors)
    
    # Find open boundary edges and signs (positive INTO region)
    openBoundaryEdges, = np.where(~dropEdges)
    openBoundarySigns = np.sign((regionCellMask[cellsOnEdge[~dropEdges, 1]] - 0.5))
    
    return openBoundaryEdges, openBoundarySigns


def get_region_edges_POP(regionCellMask):
    """Get edges of a POP masked region.
    """
    
    # Initialize lists
    j, i, vcomponents, signs = [], [], [], []
    
    # Velocity component corresponds to diff axis
    for vcomponent, axis in zip(['v', 'u'], [0, 1]):
        
        # Sign determines left/right or bottom/top region edge
        for sign in [1, -1]:
            
            # Find neigboring 0-1 pairs along axis
            ji = np.where(np.diff(regionCellMask, axis=axis, append=0) == sign)
            
            # Append indices, v-components and signs
            shape = ji[0].shape
            j.append(ji[0])
            i.append(ji[1])
            vcomponents.append(np.full(shape, vcomponent))
            signs.append(np.full(shape, sign))
    
    # Convert to arrays
    j, i, vcomponents, signs = [np.hstack(arrays) for arrays in (j, i, vcomponents, signs)]
    
    return j, i, vcomponents, signs


def get_neighbor_edges(edge, edges, verticesOnEdge, edgesOnVertex):
    """Return the neighbor edge indices of the given edge on the transect defined
    by the given index array. For MPAS meshes.
    """
    
    # Find all neighboring edges and compare to edges included in transect
    edgesAll = np.unique(edgesOnVertex[verticesOnEdge[edge, :], :])
    edgesNeighbor = list(set(edgesAll) & set(edges[edges != edge]))
    
    return edgesNeighbor


def get_neighbor_edges_POP(idx, j, i, vcomponents):
    """Return the neighbor edge indices of the given edge on the transect defined
    by the given j,i index arrays. For POP meshes. The velocity component at each
    j,i is needed to determine the applicable neighbor conditions.
    """

    # Order dims so that component dim is first (e.g., u -> i, v -> j)
    dims = [i, j] if vcomponents[idx] == 'u' else [j, i]
    
    # Get all possible neighbors adject to i and j
    neighbors = []
    for dim in dims:
        adjecent = [dim == dim[idx] + k for k in (-1, 0, 1)]
        neighbors.append(np.logical_or.reduce(adjecent))
    neighbors, = np.where(np.logical_and(*neighbors))
    neighbors = neighbors[neighbors != idx]
    
    # Get component logical and initialize conditions lists
    component_match = vcomponents[neighbors] == vcomponents[idx]
    conditions = [[component_match], [~component_match]]
    
    # Loop through dims, starting with component dim (e.g., u -> i, v -> j)
    for dim, func, k in zip(dims, ['equal', 'not_equal'], [-1, 1]):
    
        # If neighbor is same component:
        #   1. Neighbor component dim is the same
        #   2. Neighbor non-component dim is different
        conditions[0].append(getattr(np, func)(dim[neighbors], dim[idx]))
    
        # If neighbor is different component:
        #   1. Neighbor component dim is not behind
        #   2. Neighbor non-component dim is not ahead
        conditions[1].append(dim[neighbors] != dim[idx] + k)
    
    # Join conditions and select locs where conditions are true
    conditions = [np.logical_and.reduce(condition) for condition in conditions]
    conditions = np.logical_or(*conditions)
    neighbors = neighbors[conditions]

    return neighbors


def get_end_edge(edges, verticesOnEdge, edgesOnVertex):
    """Find the first end edge in edges. An end edge has only one neighbor.
    """

    # Search through edges to find the first terminal edge
    for edge in edges:
        edgeNeighbors = get_neighbor_edges(edge, edges, verticesOnEdge, edgesOnVertex)
        if len(edgeNeighbors) < 2:
            edgeEnd = edge
            break
    else:
        raise ValueError('No edge found!')
    
    return edgeEnd


def sort_edges(edges, verticesOnEdge, edgesOnVertex):
    """Sort edges surrounding a region and separate discontinuous
    sections with a fill value (-1). For MPAS meshes.
    """
    
    # Copy edges to running edgeList, initialize index_sorted,
    # and set first edge somewhere away from region
    edgeList, index_sorted, edge = np.copy(edges), [], -1
    
    while len(edgeList) > 0:
        
        # List of possible next edges
        edgesNeighbor = get_neighbor_edges(edge, edgeList, verticesOnEdge, edgesOnVertex)

        # Set edge to neighbor, or reset to new end if no edge neighbor found
        try:
            edge = edgesNeighbor[0]
        except IndexError:
            edge = get_end_edge(edgeList, verticesOnEdge, edgesOnVertex)
            index_sorted.append(-1)
        
        # Append edge and remove from running edgeList
        index_sorted.append(int(np.where(edges == edge)[0]))
        edgeList = edgeList[edgeList != edge]
    
    # Convert to numpy array
    index_sorted = np.array(index_sorted)
    
    return index_sorted


def sort_edges_POP(j, i, vcomponents, seed=0):
    """Sort edges surrounding a region. For POP meshes.
    """
    
    # Copy indices and v-components, initialize index_sorted,
    # and set first index to somewhere on the region
    j, i, vcomponents = [np.copy(array) for array in (j, i, vcomponents)]
    idx, nidx = seed, len(i)
    index_sorted = [idx]
    
    while len(index_sorted) < nidx:
        
        # Get edge neighbors (2 at start and then only 1 after)
        neighbors = get_neighbor_edges_POP(idx, j, i, vcomponents)

        # Remove idx from i, j and set next idx to neighbor
        i[idx] = -999
        j[idx] = -999
        idx = int(neighbors[0])
        index_sorted.append(idx)
    
    # Convert to numpy array
    index_sorted = np.array(index_sorted)
    
    return index_sorted


def get_edge_fields(ctgy, array, coords, interp=False):
    """Get edge fields given a variable array and a list of edge ID arrays
    """

    # Loop through transects
    ctgy = ctgy.capitalize()
    dims = {'dim': f'n{ctgy}Edges'}
    edgeFields = []
    iterables = [coords[name + ctgy] for name in ('edge', 'sign')]
    for edges, signs in zip(*iterables):
        
        # Trim NaNs
        edges, signs = [var.dropna(**dims).values.astype(int) for var in (edges, signs)]
        
        # Cell variable -> interpolate to transects
        if interp:
            edgeField = pptools.interpolate_to_edge(
                array, coords.cellsOnEdge[edges], coords.nCells,
            )
        
        # Edge velocity -> get signed values on transects
        else:
            edgeField = signs[:, None] * array[edges, :]
        
        # Append transect to list
        edgeFields.append(edgeField)
    
    return stack_arrays(edgeFields)


def get_edge_fields_POP(ctgy, array, coords, velocity=False):
    """
    """
    
    # Interpolate to edges given by vComponents
    array_interp = {}
    if velocity:
        axes, direction = [1, 2], 'down'
    else:
        axes, direction = [2, 1], 'up'
    for vComponent, axis in zip(['u', 'v'], axes):
        array_interp[vComponent] = pptools.interpolate_along_axis(array, axis, direction)
    
    # Build transects
    ctgy = ctgy.capitalize()
    dims = {'dim': f'n{ctgy}Edges'}
    edgeFields = []
    iterables = [coords[name + ctgy] for name in ('j', 'i', 'sign', 'vComponent')]
    for j, i, signs, vComponents in zip(*iterables):

        # Get transect coordinate arrays
        j, i, signs = [var.dropna(**dims).values.astype(int) for var in (j, i, signs)]
        
        # Get transect one point at a time
        edgeField = []
        for jj, ii, vComponent in zip(j, i, vComponents.values):
            edgeField.append(array_interp[vComponent][:, jj, ii])
        edgeField = np.array(edgeField)
        
        # If velocity
        if velocity:
            edgeField = signs[:, None] * edgeField
        
        # Append to list
        edgeFields.append(edgeField)
    
    return stack_arrays(edgeFields)


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


def get_transect_masks_from_regions(params, coords, meshName):
    """Get transect masks from region masks. For MPAS meshes.
    """
    
    # Get edge variables
    names = ['lonEdge', 'latEdge', 'dvEdge']
    lonEdge, latEdge, dvEdge = [coords[name].values for name in names]

    # Get pairing variables
    names = ['cellsOnEdge', 'verticesOnEdge', 'edgesOnVertex']
    cellsOnEdge, verticesOnEdge, edgesOnVertex = [coords[name].values for name in names]

    # Initialize edge masks dict
    edgeMasks = {'region': {}, 'transect': {}}

    # Loop through regions
    for regionName, regionDefs in params.items():

        # Unpack transect info in regionDefs
        transectNames, signChanges, orders, sections = list(regionDefs.values())
        orders, sections = orders[meshName], sections[meshName]

        # Get edge indices and signs bounding region, and sort edges
        regionCellMask = coords.regionCellMasks.sel(regionNames=regionName).values
        boundaryEdges, boundarySigns = get_region_edges(regionCellMask, cellsOnEdge)
        index = sort_edges(boundaryEdges, verticesOnEdge, edgesOnVertex)

        # Populate region edge mask
        idx = index[index != -1]
        edges, signs = boundaryEdges[idx], boundarySigns[idx]
        edgeMasks['region'][regionName] = {
            'edgeRegion': edges,
            'signRegion': signs,
            'dvRegion': dvEdge[edges],
        }

        # Split boundary into discrete transects; remove -1 separator and reverse order
        index = np.split(index, np.where(index == -1)[0][1:])
        index = [idx[idx >= 0][::order] for idx, order in zip(index, orders)]

        # For each transect, combine requested sections and apply sign change
        for section, transectName, signChange in zip(sections, transectNames, signChanges):
            
            # Get transect indices
            if hasattr(section, '__iter__'):
                idx = np.hstack([index[i] for i in section])
            else:
                idx = index[section]

            # Build transect coordinate variables
            edges, signs = boundaryEdges[idx], boundarySigns[idx] * signChange
            lons, lats, dv = lonEdge[edges], latEdge[edges], dvEdge[edges]
            distance = get_transect_distance(lons, lats, d0=dv[0] / 2 * 1e-3)
            edgeMasks['transect'][transectName] = {
                'edgeTransect': edges,
                'signTransect': signs,
                'lonTransect': lons,
                'latTransect': lats,
                'dvTransect': dv,
                'distTransect': distance,
            }
    
    # Stack transect coordinate variables and convert to xarray 
    ds = build_transect_xarray(edgeMasks)
    
    return ds


def get_transect_masks_from_regions_POP(params, coords, xshift=50):
    """Get transect masks from region masks. For POP meshes.
    """
    
    # Get lonEdge, latEdge and dvEdge as dictionaries of u,v components
    names = ['lon', 'lat', 'dx', 'dy']
    lons, lats, dx, dy = [coords[name + 'Edge'].values for name in names]
    lonEdge = {'u': lons, 'v': pptools.interpolate_along_axis(lons, 1, 'down')}
    latEdge = {'u': pptools.interpolate_along_axis(lats, 0, 'down'), 'v': lats}
    dvEdge = {'u': dy, 'v': dx}

    # Initialize edge masks dict
    coordNames = ['j', 'i', 'vComponent', 'sign', 'dv', 'lon', 'lat', 'dist']
    edgeMasks = {'region': {}, 'transect': {}}

    # Loop through regions
    for regionName, regionDefs in params.items():
        
        # Unpack transect info in regionDefs
        transectNames, signChanges, orders, indexes, seed = list(regionDefs.values())

        # Get edge indices and signs bounding region, and get sorting index
        regionCellMask = coords.regionCellMasks.sel(regionNames=regionName).values
        regionCellMask = np.roll(regionCellMask, xshift, axis=1)
        edgeCoords = list(get_region_edges_POP(regionCellMask))
        index = sort_edges_POP(*edgeCoords[:3], seed=seed)
        
        # Adjust i for xshift
        edgeCoords[1] = edgeCoords[1] - xshift
        invalid = edgeCoords[1] < 0
        edgeCoords[1][invalid] = edgeCoords[1][invalid] + regionCellMask.shape[1]
        
        # Sort indices and arrays
        edgeCoords = [var[index] for var in edgeCoords]
        
        # Get sorted lonEdge, latEdge and dvEdge
        edgeCoordsAux = [[], [], []]
        for j, i, vcomp in zip(*edgeCoords[:3]):
            edgeCoordsAux[0].append(dvEdge[vcomp][j, i])
            edgeCoordsAux[1].append(lonEdge[vcomp][j, i])
            edgeCoordsAux[2].append(latEdge[vcomp][j, i])
        edgeCoordsAux = [np.array(coord) for coord in edgeCoordsAux]

        # Add auxiliary coords to edgeCoords
        edgeCoords = edgeCoords + edgeCoordsAux

        # Populate region edge mask
        edgeMasks['region'][regionName] = {name + 'Region': var for name, var in zip(coordNames, edgeCoords)}
        
        # For each transect, combine requested sections and apply sign change
        for transectName, signChange, order, idx in zip(transectNames, signChanges, orders, indexes):
            
            # Build transect coordinate variables
            edgeCoordsTransect = [var[slice(*idx)][::order] for var in edgeCoords]
            edgeCoordsTransect[3] = edgeCoordsTransect[3] * signChange
            
            # Calculate distance
            dv, lons, lats = edgeCoordsTransect[4:]
            distance = get_transect_distance(lons, lats, d0=dv[0] / 2 * 1e-3)
            edgeCoordsTransect.append(distance)
            
            # Populate transect edge mask
            iterables = zip(coordNames, edgeCoordsTransect)
            edgeMasks['transect'][transectName] = {name + 'Transect': var for name, var in iterables}
    
    # Stack transect coordinate variables and convert to xarray 
    ds = build_transect_xarray(edgeMasks)
    
    return ds