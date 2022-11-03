#!/usr/bin/env python

import netCDF4
import numpy as np

def build_halos(grid, graph):
    # load mesh data
    mdata = netCDF4.Dataset(grid, 'r')

    # load processor decomposition
    gdata = np.genfromtxt(graph, dtype='int')

    # build up the computational halos
    cellsOnCell = mdata.variables['cellsOnCell'][:]-1
    nEdgesOnCell = mdata.variables['nEdgesOnCell'][:]
    nprocs = gdata.max() + 1
    halos = np.zeros((nprocs,nprocs), dtype='int')
    for acell, nedges in enumerate(nEdgesOnCell):
        neighprocs = gdata[cellsOnCell[acell, :nedges]]
        halos[gdata[acell],neighprocs] = True
    # set diagonal to be false
    for arow in np.arange(nprocs):
        halos[arow,arow] = False
    np.savez('comp_halos.npz',halos=halos)
    return halos




if __name__ == "__main__":
    import sys
    halos = build_halos(sys.argv[1], sys.argv[2])
    import pdb; pdb.set_trace()
