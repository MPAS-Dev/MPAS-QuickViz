import numpy as np
import xarray as xr

from ..config.loader import load_config
import .utilities


def get_transect_masks_from_regions(coords, regionTransectFile, regionNames=None, xshift=50):
    """
    Get transect masks from region masks

    Parameters
    ----------
    coords : xarray.Dataset
        Dataset of coordinate variables

    regionTransectFile : str
        Name of region transect config file (see config dir)

    regionNames : list or tuple
        List of region names (default use regionNames in coords)

    xshift : int
        Number of indices to x-shift the regionCellMask (prevents broken regions)

    Returns
    -------
    coords : xarray.Dataset
        Dataset of coordinate variables, including transect coordinate variables
    """

    # Set regionList to coords if not provided
    if regionNames is None:
        regionList = coords.regionNames.values
    else:
        regionList = regionNames

    # Load region transect params
    regionTransects = load_config(regionTransectFile)

    # Unpack coordinate variables
    lonEdge = coords.lonEdge.values
    latEdge = coords.latEdge.values
    dx = coords.dxEdge.values
    dy = coords.dyEdge.values
    
    # Rebuild coordinate variables into u,v components
    lonEdge = {
        'u': lonEdge,
        'v': _interpolate_along_axis(lonEdge, 1, 'down'),
    }
    latEdge = {
        'u': _interpolate_along_axis(latEdge, 0, 'down'),
        'v': latEdge,
    }
    dvEdge = {
        'u': dy,
        'v': dx,
    }

    # Initialize edge masks dict
    edgeMasks = {'region': {}, 'transect': {}}

    # Loop through regions
    for regionName in regionList:

        # Get regionParams
        try:
            regionParams = regionTransects[regionName]
        except KeyError:
            msg = f'No entry for region {regionName} found in regionTransectFile'
            raise KeyError(msg)

        # Get edge indices and signs bounding region, and get sorting index
        regionCellMask = coords.regionCellMasks.sel(regionNames=regionName).values
        regionCellMask = np.roll(regionCellMask, xshift, axis=1)
        j, i, vcomps, signs = _get_region_edges(regionCellMask)
        index = _sort_edges(j, i, vcomps, seed=regionParams['seed'])
        
        # Adjust i for xshift
        i = i - xshift
        i[i < 0] = i[i < 0] + regionCellMask.shape[1]
        
        # Sort indices and arrays
        j, i, vcomps, signs = [var[index] for var in (j, i, vcomps, signs)]
        
        # Get sorted lonEdge, latEdge and dvEdge
        variables = [[], [], []]
        for jj, ii, vcomp in zip(j, i, vcomps):
            variables[0].append(lonEdge[vcomp][jj, ii])
            variables[1].append(latEdge[vcomp][jj, ii])
            variables[2].append(dvEdge[vcomp][jj, ii])
        lonEdgeRegion, latEdgeRegion, dvEdgeRegion = [np.array(var) for var in variables]

        # Populate region edge mask
        edgeMasks['region'][regionName] = {
            'jRegion': j,
            'iRegion': i,
            'vCompRegion': vcomps,
            'signRegion': signs,
            'dvRegion': dvEdgeRegion,
        }

        # For each transect, combine requested sections and apply sign change
        for transectName, transectParams in regionParams['transects'].items():

            # Get index slice and order
            idx = slice(*transectParams['segmentIndexes'])
            order = -1 if transectParams['reverseOrder'] else 1
            
            # Build transect coordinate variables and reverse order
            variables = [j, i, vcomps, signs, lonEdgeRegion, latEdgeRegion, dvEdgeRegion]
            jj, ii, vcomp, sign, lons, lats, dv = [var[idx][::order] for var in variables]

            # Change signs
            if transectParams['changeSign']:
                sign = -1 * sign

            # Get remaining coordinate variables
            distance = utilities.get_transect_distance(lons, lats, d0=dv[0] / 2 * 1e-3)
            
            # Populate transect edge mask
            edgeMasks['transect'][transectName] = {
                'jTransect': jj,
                'iTransect': ii,
                'vCompTransect': vcomp,
                'signTransect': sign,
                'lonTransect': lons,
                'latTransect': lats,
                'dvTransect': dv,
                'distTransect': distance,
            }
    
    # Stack transect coordinate variables and convert to xarray 
    ds = utilities.build_transect_xarray(edgeMasks)

    # Merge with coords and return
    coords = xr.merge([coords] + list(ds.values()))
    
    return coords


def get_edge_fields(ctgy, array, coords, velocity=False):
    """
    Get edge fields given a variable array and a list of edge ID arrays

    Parameters
    ----------
    ctgy : str
        Edges category (region or transect)

    array : numpy.ndarray
        Variable array (3D) from which to extract edge fields

    coords : xarray.Dataset
        Dataset of coordinate variables

    velocity : bool
        Use velocity conventions (POP B-grid)
        If False, interpolate cell quantities to edge

    Returns
    -------
    edgeFields : numpy.ndarray
        Array of the given variable on the set of region or transect edges
    """
    
    # Interpolate to edges given by v-components
    array_interp = {}
    if velocity:
        axes, direction = [1, 2], 'down'
    else:
        axes, direction = [2, 1], 'up'
    for vcomp, axis in zip(['u', 'v'], axes):
        array_interp[vcomp] = _interpolate_along_axis(array, axis, direction)
    
    # Build transects
    ctgy = ctgy.capitalize()
    dims = {'dim': f'n{ctgy}Edges'}
    edgeFields = []
    iterables = [coords[name + ctgy] for name in ('j', 'i', 'sign', 'vComp')]
    for j, i, signs, vcomps in zip(*iterables):

        # Get transect coordinate arrays
        j, i, signs = [var.dropna(**dims).values.astype(int) for var in (j, i, signs)]
        
        # Get transect one point at a time
        edgeField = []
        for jj, ii, vcomp in zip(j, i, vcomps.values):
            edgeField.append(array_interp[vcomp][:, jj, ii])
        edgeField = np.array(edgeField)
        
        # If velocity
        if velocity:
            edgeField = signs[:, None] * edgeField
        
        # Append to list
        edgeFields.append(edgeField)
    
    return utilities.stack_arrays(edgeFields)


def _interpolate_along_axis(array, axis, direction):
    """
    Interpolate array along axis.
    
    Direction `up` interpolates T points toward U points.
    Direction `down` interpolates U points toward T points.

    Parameters
    ----------
    array : numpy.ndarray
        Cell-centered array
    
    axis : int
        Axis along which to interpolate

    direction : str
        Direction of interpolation (up or down)

    Returns
    -------
    array_interp : numpy.ndarray
        Array interpolated to edge midpoints
    """
    
    # Get indexing based on interpolation direction
    n = array.shape[axis]
    if direction == 'up':
        idx, order = n - 1, -1
    else:
        idx, order = 0, 1
        
    # Expand array across the boundary of the interpolation dimension
    boundary = np.expand_dims(array.take(idx, axis), axis)
    array = [boundary, array][slice(None, None, order)]
    array = np.concatenate(array, axis)
    
    # Interpolate linearly
    array_interp = []
    for start in [0, 1]:
        array_interp.append(array.take(range(start, start + n), axis))
    array_interp = sum(array_interp) / 2
    
    return array_interp


def _get_region_edges(regionCellMask):
    """
    Get edges of a masked region.

    Parameters
    ----------
    regionCellMask : numpy.ndarray
        Region cell mask array (nlat, nlon)

    Returns
    -------
    j : numpy.ndarray
        Array of j indices that bound the region

    i : numpy.ndarray
        Array of i indices that bound the region

    vcomps : numpy.ndarray
        Array of v-components at each ji index pair (u or v)

    signs : numpy.ndarray
        Array of signs at each ji index pair (u or v). When multiplied,
        the signs transform the edge-normal velocity to positive-inward
    """
    
    # Initialize lists
    j, i, vcomps, signs = [], [], [], []
    
    # Velocity component corresponds to diff axis
    for vcomp, axis in zip(['v', 'u'], [0, 1]):
        
        # Sign determines left/right or bottom/top region edge
        for sign in [1, -1]:
            
            # Find neigboring 0-1 pairs along axis
            ji = np.where(np.diff(regionCellMask, axis=axis, append=0) == sign)
            
            # Append indices, v-components and signs
            shape = ji[0].shape
            j.append(ji[0])
            i.append(ji[1])
            vcomps.append(np.full(shape, vcomp))
            signs.append(np.full(shape, sign))
    
    # Convert to arrays
    j, i, vcomps, signs = [np.hstack(arrays) for arrays in (j, i, vcomps, signs)]
    
    return j, i, vcomps, signs


def _get_neighbor_edges(idx, j, i, vcomps):
    """
    Return the neighbor edge indices of the given edge on the
    transect defined by the given ji index arrays.
    
    The velocity component (u or v) at each ji is needed to
    determine the applicable neighbor conditions.

    Parameters
    ----------
    idx : int
        Index of the given edge
    
    j : numpy.ndarray
        Array of transect j indices

    i : numpy.ndarray
        Array of transect i indices

    vcomps : numpy.ndarray
        Array of v-components at each ji index pair (u or v)

    Returns
    -------
    neighbors : numpy.ndarray
        Array of neighbors to idx
    """

    # Order dims so that component dim is first (e.g., u -> i, v -> j)
    dims = [i, j] if vcomps[idx] == 'u' else [j, i]
    
    # Get all possible neighbors adject to i and j
    neighbors = []
    for dim in dims:
        adjecent = [dim == dim[idx] + k for k in (-1, 0, 1)]
        neighbors.append(np.logical_or.reduce(adjecent))
    neighbors, = np.where(np.logical_and(*neighbors))
    neighbors = neighbors[neighbors != idx]
    
    # Get component logical and initialize conditions lists
    component_match = vcomps[neighbors] == vcomps[idx]
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


def _sort_edges(j, i, vcomps, seed=0):
    """
    Sort edges surrounding a region.

    Parameters
    ----------
    j : numpy.ndarray
        Array of transect j indices

    i : numpy.ndarray
        Array of transect i indices

    vcomps : numpy.ndarray
        Array of v-components at each ji index pair (u or v)

    seed : int
        Seed index to begin sorting

    Returns
    -------
    index_sorted : numpy.ndarray
        Array of indices for sorting the edges ji arrays
    """
    
    # Copy indices and v-components, initialize index_sorted,
    # and set first index to somewhere on the region
    j, i, vcomps = [np.copy(array) for array in (j, i, vcomps)]
    idx, nidx = seed, len(i)
    index_sorted = [idx]
    
    while len(index_sorted) < nidx:
        
        # Get edge neighbors (2 at start and then only 1 after)
        neighbors = _get_neighbor_edges(idx, j, i, vcomps)

        # Remove idx from i, j and set next idx to neighbor
        i[idx] = -999
        j[idx] = -999
        idx = int(neighbors[0])
        index_sorted.append(idx)
    
    # Convert to numpy array
    index_sorted = np.array(index_sorted)
    
    return index_sorted
