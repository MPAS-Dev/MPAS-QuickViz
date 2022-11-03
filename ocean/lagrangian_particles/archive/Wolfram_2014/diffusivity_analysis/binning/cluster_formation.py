#!/usr/bin/env python
"""

    Compute particle clusters
    
    Phillip Wolfram
    LANL
    03/05/2015

"""
# import libraries / packages
import numpy as np
from scipy.spatial import cKDTree as KDTree
import matplotlib.pyplot as plt

class ParticleClusters(): #{{{
    """ParticleCluster class to perform clustering operations on particles
    Phillip J. Wolfram
    LANL
    03/05/2015

    Initializes particle cluster data structures based on particle locations.
    Should specify whether to use 'ball' (low-memory, high-comp) or 
    'tree' (high-memory, low-comp) structures.


    Note:
      Depending upon memory limitations, may choose to use this
      in either memory or computate efficiency modes.

    Args:
        partx(ndarray): Location of particles in x coordinate.
        party(ndarray): Location of particles in y coordinate.
        xloc(ndarray): Location of cluster centers in x coordinate.
        yloc(ndarray): Location of cluster centers in y coordinate.
        radius(float): Radius of clusters to be formed.
        mode(int): Use 'ball' or 'tree' structures.


    Attributes:
        radius(float): Radius of clusters to be formed.
        mode(int): Specify whether the class should operate in 'ball' or 'tree' mode.
        clusters(list): List of indicies of cluster particles
        allparticles(KDTree): KDTree of particles (if clusters are not directly computed).

    """
    def __init__(self, partx=None, party=None, xloc=None, yloc=None, radius=None, mode='ball'):
        self.mode = mode
        self.radius = radius
        
        allparticles = KDTree(np.vstack((partx, party)).T)

        if(mode is 'tree'):
            search = KDTree(np.vstack((xloc,yloc)).T)
            # can just use cluster data structure
            self.clusters = search.query_ball_tree(allparticles, radius)
            self.allparticles = None
        if(mode is 'cluster'):
            maxIndex = np.load('clusterDefinition/maxVertices.p')
            indexOnCluster = np.load('clusterDefinition/verticesOnCell.p')
            indexToParticle = np.load('clusterDefinition/verticesToParticle.p')
            clusters = np.nan*np.zeros((maxIndex.shape[0],maxIndex.max()))
            for aCluster, maxInd in enumerate(maxIndex):
                clusters[aCluster,0:maxInd] = indexToParticle[indexOnCluster[aCluster,0:maxInd]]
            self.clusters = clusters
        else:
            self.allparticles = allparticles
            self.clusters = None

    # multiple argument functions
    def get_cluster_ball(self, xloc, yloc, radius):
        """ Build cluster at (xloc,yloc) with radius and return list of particles in allparticles."""
        return self.allparticles.query_ball_point((xloc,yloc), radius)

    def get_cluster_index(self,ac):
        return self.clusters[ac,~np.isnan(self.clusters[ac,:])].astype('int')


#}}}

if __name__ == "__main__":
    print 'Library to build particle clusters'
