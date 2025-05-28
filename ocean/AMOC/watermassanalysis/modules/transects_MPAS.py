import numpy as np
import xarray as xr

from ..config.loader import load_config
import .utilities


def get_transect_masks_from_regions(coords, regionTransectFile, regionNames=None):
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
    dvEdge = coords.dvEdge.values
    cellsOnEdge = coords.cellsOnEdge.values
    verticesOnEdge = coords.verticesOnEdge.values
    edgesOnVertex = coords.edgesOnVertex.values

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

        # Get edge indices and signs bounding region, and sort edges
        regionCellMask = coords.regionCellMasks.sel(regionNames=regionName).values
        boundaryEdges, boundarySigns = _get_region_edges(regionCellMask, cellsOnEdge)
        index = _sort_edges(boundaryEdges, verticesOnEdge, edgesOnVertex)

        # Get transect indices, edges, and signs
        idx = index[index != -1]
        edges, signs = boundaryEdges[idx], boundarySigns[idx]

        # Populate region edge mask
        edgeMasks['region'][regionName] = {
            'edgeRegion': edges,
            'signRegion': signs,
            'dvRegion': dvEdge[edges],
        }

        # Split boundary into discrete segments
        index = np.split(index, np.where(index == -1)[0][1:])

        # Get reverseOrder for each segment and convert from T/F to -1/1
        orders = -np.array(regionParams['reverseOrder']).astype(int)
        orders[orders >= 0] = 1

        # For each segment, remove -1 separator and reverse order 
        index = [idx[idx >= 0][::order] for idx, order in zip(index, orders)]

        # For each transect, combine requested sections and apply sign change
        for transectName, transectParams in regionParams['transects'].items():
            
            # Get transect indices, edges, and signs
            idx = np.hstack([index[i] for i in transectParams['segmentIndexes']])
            edges, signs = boundaryEdges[idx], boundarySigns[idx]

            # Change signs
            if transectParams['changeSign']:
                signs = -1 * signs

            # Get remaining coordinate variables
            lons, lats, dv = lonEdge[edges], latEdge[edges], dvEdge[edges]
            distance = utilities.get_transect_distance(lons, lats, d0=dv[0] / 2 * 1e-3)

            # Populate transect edge mask
            edgeMasks['transect'][transectName] = {
                'edgeTransect': edges,
                'signTransect': signs,
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
        Use velocity conventions (MPAS C-grid)
        If False, interpolate cell quantities to edge

    Returns
    -------
    edgeFields : numpy.ndarray
        Array of the given variable on the set of region or transect edges
    """

    # Loop through transects
    ctgy = ctgy.capitalize()
    dims = {'dim': f'n{ctgy}Edges'}
    edgeFields = []
    iterables = [coords[name + ctgy] for name in ('edge', 'sign')]
    for edges, signs in zip(*iterables):
        
        # Trim NaNs
        edges, signs = [var.dropna(**dims).values.astype(int) for var in (edges, signs)]

        # Edge velocity -> get signed values on transects
        if velocity:
            edgeField = signs[:, None] * array[edges, :]
        
        # Cell variable -> interpolate to transects
        else:
            edgeField = _interpolate_to_edge(
                array, coords.cellsOnEdge[edges], coords.nCells,
            )
        
        # Append transect to list
        edgeFields.append(edgeField)
    
    return utilites.stack_arrays(edgeFields)


def _interpolate_to_edge(varCell, cellsOnEdge, nCells):
    """
    Interpolate cell variable to edge

    Parameters
    ----------
    varCell : numpy.ndarray
        Array at cell centers

    cellsOnEdge : numpy.ndarray
        Coordinate variable defining cell IDs on each edge

    nCells : numpy.ndarray
        Cell IDs array

    Returns
    -------
    varCell : numpy.ndarray
        Array at cell edges
    """
    
    # Find cells on edge in subdomain
    index = list(np.array([np.where(np.isin(nCells, pair))[0] for pair in cellsOnEdge]).T)

    # Interpolate to edge
    varEdge = (varCell[index[0], :] + varCell[index[1], :]) / 2
    
    return varEdge


def _get_region_edges(regionCellMask, cellsOnEdge):
    """
    Get open water edges of an MPAS-Ocean masked region.

    Uses Python indexing, so subtract 1 from cellsOnEdge before
    calling this function. Code adapted from Alice Barthel.

    Parameters
    ----------
    regionCellMask : numpy.ndarray
        Region cell mask array (nCells,)

    cellsOnEdge : numpy.ndarray
        Cells on edge array, Python-indexed (nEdges, 2)

    Returns
    -------
    openBoundaryEdges : numpy.ndarray
        Edge indices (nBoundaryEdges,)

    openBoundarySigns : numpy.ndarray
        Edge signs (nBoundaryEdges,)
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


def _get_neighbor_edges(edge, edges, verticesOnEdge, edgesOnVertex):
    """
    Return the neighbor edge indices of the given edge on the
    transect defined by the given index array.

    Parameters
    ----------
    edge : int
        Given edge for which neighbor edges are requested

    edges : numpy.ndarray
        Edges array

    verticesOnEdge : numpy.ndarray
        Vertices on edge array, Python-indexed (nEdges, 2)

    edgesOnVertex : numpy.ndarray
        Edges on vertex array, Python-indexed (nVertices, 3)

    Returns
    -------
    edgesNeighbor : list
        List of neighbor edges
    """
    
    # Find all neighboring edges and compare to edges included in transect
    edgesAll = np.unique(edgesOnVertex[verticesOnEdge[edge, :], :])
    edgesNeighbor = list(set(edgesAll) & set(edges[edges != edge]))
    
    return edgesNeighbor


def _get_end_edge(edges, verticesOnEdge, edgesOnVertex):
    """
    Find the first end edge in edges. An end edge has only one neighbor.

    Parameters
    ----------
    edges : numpy.ndarray
        Edges array

    verticesOnEdge : numpy.ndarray
        Vertices on edge array, Python-indexed (nEdges, 2)

    edgesOnVertex : numpy.ndarray
        Edges on vertex array, Python-indexed (nVertices, 3)

    Returns
    -------
    edgeEnd : int
        Index of end edge
    """

    # Search through edges to find the first terminal edge
    for edge in edges:
        edgeNeighbors = _get_neighbor_edges(edge, edges, verticesOnEdge, edgesOnVertex)
        if len(edgeNeighbors) < 2:
            edgeEnd = edge
            break
    else:
        raise ValueError('No edge found!')
    
    return edgeEnd


def _sort_edges(edges, verticesOnEdge, edgesOnVertex):
    """
    Sort edges surrounding a region and separate discontinuous
    sections with a fill value (-1).

    Parameters
    ----------
    edges : numpy.ndarray
        Edges array

    verticesOnEdge : numpy.ndarray
        Vertices on edge array, Python-indexed (nEdges, 2)

    edgesOnVertex : numpy.ndarray
        Edges on vertex array, Python-indexed (nVertices, 3)

    Returns
    -------
    index_sorted : numpy.ndarray
        Array of indices for sorting the edges array
    """
    
    # Copy edges to running edgeList, initialize index_sorted,
    # and set first edge somewhere away from region
    edgeList = np.copy(edges)
    index_sorted = []
    edge = -1
    
    while len(edgeList) > 0:
        
        # List of possible next edges
        edgesNeighbor = _get_neighbor_edges(edge, edgeList, verticesOnEdge, edgesOnVertex)

        # Set edge to neighbor, or reset to new end if no edge neighbor found
        try:
            edge = edgesNeighbor[0]
        except IndexError:
            edge = _get_end_edge(edgeList, verticesOnEdge, edgesOnVertex)
            index_sorted.append(-1)
        
        # Append edge and remove from running edgeList
        index_sorted.append(int(np.where(edges == edge)[0]))
        edgeList = edgeList[edgeList != edge]
    
    # Convert to numpy array
    index_sorted = np.array(index_sorted)
    
    return index_sorted
