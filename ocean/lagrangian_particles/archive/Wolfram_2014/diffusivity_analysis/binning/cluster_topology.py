#!/usr/bin/env python
"""
Module computes cluster topologies from
MPAS grid data structure.

Phillip J Wolfram
01/15/2016
"""
import xarray as xr
import numpy as np
from itertools import combinations, chain
from scipy.misc import comb
import matplotlib.pyplot as plt
from latlon_coordinate_transforms import fix_periodicity_numexpr as fix_periodicity

def comb_index(n, k): #{{{
    count = comb(n, k, exact=True)
    index = np.fromiter(chain.from_iterable(combinations(range(n), k)), int, count=count*k)
    return index.reshape(-1, k) #}}}

def cluster_topology(mesh, acell, x, y, radius): #{{{
    """
    Computes cluster topology from mesh file, e.g., mesh.nc
    Returns:
      - cells, a list of cells making up the cluster
      - meshedges, a list of points within cells making up the nearest neighbor connections
      - ringedges, a list of points within cells making up the outer ring of the cluster
    """
    # now loop over all cell neighbors from point to accertain correct set of
    # particles
    cells = [acell]
    outercells = [acell]
    xc = x[acell]
    yc = y[acell]
    rr = 0
    while(rr < radius):
        tmp = []
        for ac in outercells:
            neighs = mesh.cellsOnCell[ac,:mesh.nEdgesOnCell[ac].values].values-1
            for cellneighs in np.unique(neighs[np.where(neighs >= 0)]):
                tmp.append(cellneighs)
            tmp = np.unique(tmp).tolist()
        # tmp has neighbors of all outer cells, but want it to store all known cells
        tmp += cells
        outercells = np.setxor1d(tmp,cells).tolist()
        cells = np.unique(tmp).tolist()
        # compute rr from maximum distance from center
        rr = np.min(np.sqrt((x[outercells]-xc)**2.0 + (y[outercells]-yc)**2.0))
    # remove particles that are too far away from the center (needed for nonuniform grid)
    inrad = np.sqrt((x[cells]-xc)**2.0 + (y[cells]-yc)**2.0) <= radius
    # final cell of cells in disk
    cells = np.asarray(cells)[inrad]

    # clean up hanging nodes from cluster
    connections = np.zeros((len(mesh.nCells)),dtype='i')
    cellsOnCell = mesh.cellsOnCell.values - 1
    for ac in cells:
        ns = mesh.nEdgesOnCell[ac].values
        connections[np.intersect1d(cells,cellsOnCell[ac,:ns])] += 1
    cells = cells[connections[cells] > 2]

    # get the mesh edges
    cellsOnEdge = mesh.cellsOnEdge.values - 1
    goodedges = np.logical_and(np.in1d(cellsOnEdge[:,0],cells), \
                               np.in1d(cellsOnEdge[:,1],cells))
    meshedges = cellsOnEdge[goodedges,:]
    cellid = np.arange(len(mesh.nCells))
    cellid[cells] = np.arange(len(cells))
    meshedges = cellid[meshedges]

    # get the ring edges
    ringedges = np.copy(goodedges)
    edgesOnCell = mesh.edgesOnCell.values - 1
    for ac in cells:
        ns = mesh.nEdgesOnCell[ac].values
        if np.prod(goodedges[edgesOnCell[ac,:ns]]):
            ringedges[edgesOnCell[ac,:ns]] = False
    ringedges = cellsOnEdge[ringedges,:]
    ringedges = cellid[ringedges]

    # all edge pairs
    #alledges = cells[comb_index(len(cells),2)]

    return cells, meshedges, ringedges #}}}

def test(mname='mesh.nc', radius = 1.0e5):
    mesh = xr.open_dataset(mname)
    clusters = np.arange(len(mesh.nCells))
    np.random.shuffle(clusters)
    for acell in clusters:
        xc = mesh.xCell.values[acell]
        yc = mesh.yCell.values[acell]
        # transform all coordinates into correct coordinates if periodic
        x = mesh.xCell.values
        y = mesh.yCell.values
        if mesh.is_periodic == "YES":
            if mesh.x_period > 0:
                x = fix_periodicity(x,xc,mesh.x_period)
            if mesh.y_period > 0:
                y = fix_periodicity(y,yc,mesh.y_period)

        cells, meshedges, ringedges = cluster_topology(mesh, acell, x, y, radius)

        plt.figure()
        plt.hold(True)
        plt.axis('equal')
        xc = x[cells]
        yc = y[cells]
        plt.plot(xc,yc,'.',zorder=1)
        for edge in meshedges:
            plt.plot(xc[edge],yc[edge],'k-',lw=2,zorder=0, alpha=0.3)
        for edge in ringedges:
            plt.plot(xc[edge],yc[edge],'r-',lw=1,zorder=2)
        plt.title('cell = %d'%acell)
        plt.savefig('cell%05d'%(acell)+'.png')
        plt.close()

if __name__ == "__main__":
    test()
