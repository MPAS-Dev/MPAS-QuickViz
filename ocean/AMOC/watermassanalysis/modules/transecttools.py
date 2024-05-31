#!/usr/bin/env python
"""
    Name: transecttools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for working with MPAS-Ocean transects
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import yaml


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


def get_neighbor_edges(edge, edges, verticesOnEdge, edgesOnVertex):
    """Return the neighbor edge indices of the given edge on the transect defined
    by the given index array.
    """
    
    # Find all neighboring edges and compare to edges included in transect
    edgesAll = np.unique(edgesOnVertex[verticesOnEdge[edge, :], :])
    edgesNeighbor = list(set(edgesAll) & set(edges[edges != edge]))
    
    return edgesNeighbor


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
    sections with a fill value (-1).
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


def get_transect_masks_from_regions(meshName, coords):
    """Get transect masks from region masks
    """
    
    # Get transect definitions for manually defining transects
    with open('../yaml/transect_definitions.yaml') as f:
        transectDefs = yaml.safe_load(f)
    
    # Loop through regions
    transectMasks = {}
    regionNames = ['Labrador Sea', 'Irminger-Iceland Basins', 'Nordic Seas']
    for regionName in regionNames:
        
        # Get regionCellMask from regionName
        col = np.where(coords['regionNames'] == regionName)[0][0]
        regionCellMask = coords['regionCellMasks'][:, col]
        
        # Get transect definitions for region
        names = ('transectNames', 'signChanges', 'orders', 'sections')
        transectNames, signChanges = [transectDefs[regionName][name] for name in names[:2]]
        orders, sections = [transectDefs[regionName][name][meshName] for name in names[2:]]
        
        # Get edge indices and signs bounding region
        edges, signs = get_region_edges(regionCellMask, coords['cellsOnEdge'])
        
        # Sort edges and split into discrete sections
        index = sort_edges(edges, coords['verticesOnEdge'], coords['edgesOnVertex'])
        index = np.split(index, np.where(index == -1)[0][1:])

        # Remove -1 separator and reverse order of each section
        index = [idx[idx > 0][::order] for idx, order in zip(index, orders)]
        
        # For each transect, combine requested sections and apply sign change
        for section, transectName, signChange in zip(sections, transectNames, signChanges):
            idx = np.hstack([index[i] for i in section])
            dvEdge = coords['dvEdge'][edges[idx]]
            transectMasks[transectName] = {
                'edges': edges[idx],
                'signs': signs[idx] * signChange,
                'distance': np.cumsum(dvEdge) - dvEdge / 2,
            }
    
    return transectMasks
