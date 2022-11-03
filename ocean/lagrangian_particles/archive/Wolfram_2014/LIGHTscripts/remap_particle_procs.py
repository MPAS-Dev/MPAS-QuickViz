#!/usr/bin/env python
"""
    Routines to remap particles onto a new grid decomposition.
    Phillip Wolfram
    LANL
    08/19/2014
"""

import netCDF4
from scipy import spatial
import numpy as np

def remap_particles(fin, fpart, fdecomp):
    """
    load in particle positions, locations of grid cell centers, and decomposition 
    tree corresponding to fin.  

    The goal is to update particle field currentBlock to comply with the new grid
    as defined by fin and decomp.  FIN AND FDECOMP MUST BE COMPATIBLE!

    We assume that all particles will be within the domain such that a nearest
    neighbor search is sufficient to make the remap.

    Phillip Wolfram
    LANL
    08/19/2014
    """
    # load the files
    f_in = netCDF4.Dataset(fin, 'r')
    f_part = netCDF4.Dataset(fpart,'r+')

    # get the particle data
    xpart = f_part.variables['xParticle']
    ypart = f_part.variables['yParticle']
    zpart = f_part.variables['zParticle']
    currentBlock = f_part.variables['currentBlock']
    currentBlockReset = f_part.variables['currentBlockReset']
    try:
        currentCell = f_part.variables['currentCell']
    except:
        currentCell = f_part.createVariable('currentCell', 'i', ('nParticles'))
    try:
        currentCellReset = f_part.variables['currentCellReset']
    except:
        currentCellReset = f_part.createVariable('currentCellReset', 'i', ('nParticles'))

    # get the cell positions 
    xcell = f_in.variables['xCell']
    ycell = f_in.variables['yCell']
    zcell = f_in.variables['zCell']

    # build the spatial tree
    tree = spatial.cKDTree(np.vstack((xcell,ycell,zcell)).T)

    # get nearest cell for each particle
    dvEdge = f_in.variables['dvEdge']
    maxdist = 2.0*max(dvEdge[:])
    _, cellIndices = tree.query(np.vstack((xpart,ypart,zpart)).T,distance_upper_bound=maxdist,k=1)

    # load the decomposition
    decomp = np.genfromtxt(fdecomp)
    currentBlock[:] = decomp[cellIndices]
    currentBlockReset[:] = decomp[cellIndices]
    currentCell[:] = cellIndices + 1
    currentCellReset[:] = cellIndices + 1

    # close the files
    f_in.close()
    f_part.close()



if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="f_in",
                      help="file to open for appending \
                      particle data 'particle_' extension",
                      metavar="FILE")
    parser.add_option("-p", "--particles", dest="f_part",
                      help="file for particles",
                      metavar="FILE")
    parser.add_option("-d", "--decomp", dest="f_decomp",
                      help="file for decoposition, graph.info.part.# \
                      where # is the number of blocks for the \
                      decomposition",
                      metavar="FILE")

    options, args = parser.parse_args()

    if not options.f_in:
        parser.error("Filename is a required input.")
    if not options.f_part:
        parser.error("Particle filename is a required input.")
    if not options.f_decomp:
        parser.error("Decomposition filename is a required input.")


    remap_particles(options.f_in, options.f_part, options.f_decomp)


