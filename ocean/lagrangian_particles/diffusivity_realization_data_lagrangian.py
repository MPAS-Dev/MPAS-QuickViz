#!/usr/bin/env python
"""

    Get data from each realization.  Script produces the following files
        com_time.p
        dl_time.p
    which cane be parsed by an ensembler to compute diffusivity.
    
    Phillip Wolfram
    LANL
    09/19/2014

"""
# import libraries / packages
import os
import numpy as np
import netCDF4
import _pickle as pickle

# function definitions

def spherical_bearing(phi1, phi2, lam1, lam2): #{{{
  """
  compute the spherical bearing
      http://www.movable-type.co.uk/scripts/latlong.html
      
      where lambda is longitude and phi is latitude
  
      spherical_bearing returned is from N (-pi, pi)
  
      Phillip Wolfram
      LANL
      08/07/2014
  """
  dphi = phi2 - phi1
  dlam = lam2 - lam1

  return np.arctan2(np.sin(dlam)*np.cos(phi2), np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*np.cos(dlam)) #}}}

def normalized_haversine_formula(phi1, phi2, lam1, lam2):  #{{{
    """
    compute the distance between two points via Haversine formula:
    http://www.movable-type.co.uk/scripts/latlong.html
    
    where lambda is longitude and phi is latitude

    c returned in non-dimensional units (radians)

    Phillip Wolfram
    LANL
    07/18/2014
    """
    dphi = phi2 - phi1 
    dlam = lam2 - lam1
    
    a = np.sin(dphi/2.0)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlam/2.0)**2
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0-a))

    return c # distance is c * Radius_of_earth #}}}

def proj_lat_long(x, y, z):  #{{{
    """
    compute the latitude and longitude from
    the x,y,z points (follow's Doug's pnt.h)
    """
    plat = np.arcsin(z / np.sqrt(x ** 2 + y ** 2 + z ** 2))
    plong = np.arctan2(y, x) 

    return plat, plong  #}}}

def proj_xyz(plat, plong, r):  #{{{
    """
    convert from latitude / longitude to xyz
    Phillip Wolfram
    LANL
    07/17/2014
    """
    x = r * np.cos(plong) * np.cos(plat)
    y = r * np.sin(plong) * np.cos(plat)
    z = r * np.sin(plat)
    
    return x,y,z  #}}}

def get_clusters(m_in, cluster_path): #{{{
    """
    Compute or load the cluster information to link particles
    to clusters to compute dispersion.

    Phillip Wolfram, LANL
    09/19/2014
    """
    # compute interior cells and masks for vertices 
    if cluster_path is None:
      print ('Computing particle clusters')
      verticesOnCell = m_in.variables['verticesOnCell']
      nEdgesOnCell   = m_in.variables['nEdgesOnCell']
      cellsOnVertex  = m_in.variables['cellsOnVertex'] 
      boundaryVertex = (cellsOnVertex[:,:] == 0)
      boundaryVertex = boundaryVertex.max(axis=1)
      
      nCells = len(m_in.dimensions['nCells'])
      inCells = np.ones(nCells+1)
      inCells[cellsOnVertex[boundaryVertex,:]-1] = 0
      boundaryCells = (inCells == 0)
      inCells = np.where(inCells > 0)
      verticesOnCell = verticesOnCell[inCells[0],:] - 1
      maxVertices = nEdgesOnCell[inCells[0]]
      verticesToParticle = np.cumsum(1-boundaryVertex) - 1

      # store cluster neighbor information for filter
      print ('saving verticesOnCell, maxVertices, verticesToParticle')
      pickle.dump(verticesToParticle, open( "verticesToParticle.p", "wb"))
      pickle.dump(verticesOnCell, open( "verticesOnCell.p", "wb"))
      pickle.dump(maxVertices, open( "maxVertices.p", "wb"))
      print ('done')
    else:
      print ('laoding cluster info')
      verticesToParticle = pickle.load(open(cluster_path+"/verticesToParticle.p","rb"))
      verticesOnCell = pickle.load(open(cluster_path+"/verticesOnCell.p","rb"))
      maxVertices = pickle.load(open(cluster_path+"/maxVertices.p","rb"))
      print ('done')

    return verticesOnCell, maxVertices, verticesToParticle #}}}

def compute_com_and_dispersion(plat, plon, r, indexOnCluster, maxIndex, indexToParticle): #{{{
    """
    compute particle disperison (2nd moment) for cluster, using lat / lon for basis of 
    calculation to set coordinates for Kxx,Kxy,Kyy where x is zonal and y meridional

    dispersion units in m^2

    Phillip Wolfram
    LANL
    07/18/2014
    """

    N = maxIndex.shape[0]
    clat = np.zeros(N)
    clon = np.zeros(N)
    dxdx = np.zeros(N)
    dxdy = np.zeros(N)
    dydy = np.zeros(N)
    drdr = np.zeros(N)
    Npart = np.zeros(N)
    
    print ('Computing COM and dispersion for ', N, ' clusters')
    for aCluster, maxInd in enumerate(maxIndex):
        # get points of cluster
        particles = indexToParticle[indexOnCluster[aCluster,0:maxInd]]
        pslat = plat[particles]
        pslon = plon[particles]
        
        # compute center of mass 
        clat[aCluster] = np.mean(pslat)
        clon[aCluster] = np.mean(pslon)

        # compute distances in m from lat / long  (use COM for simplicity, but there will be some error because of coordinate transform)
        dx = r * normalized_haversine_formula(clat[aCluster], clat[aCluster], clon[aCluster], pslon)
        dy = r * normalized_haversine_formula(clat[aCluster], pslat,          clon[aCluster], clon[aCluster])
        dr = r * normalized_haversine_formula(clat[aCluster], pslat,          clon[aCluster], pslon)
        
        # fix orientation of points
        bearing = spherical_bearing(clat[aCluster], pslat, clon[aCluster], pslon)
        # because arctan2 returns results from -pi to pi for bearing, flip values to get right sign
        dx -= 2*dx*(bearing < 0)
        dy -= 2*dy*(np.fabs(bearing) > np.pi/2.0)

        # store values
        Nparticles = len(particles)
        dxdx[aCluster] = sum(dx*dx)
        dxdy[aCluster] = sum(dx*dy)
        dydy[aCluster] = sum(dy*dy)
        drdr[aCluster] = sum(dr*dr)
        Npart[aCluster] = Nparticles

    return clon, clat, dxdx, dxdy, dydy, drdr, Npart #}}}

def build_particle_file(fname_in, mname_in, cluster_path):  #{{{
    """
    Take existing netCDF4 input f_in and produce
    vtu files for f_outbase with all the particle data

    Phillip Wolfram 
    LANL
    07/15/2014
    """
    # load the netCDF database
    f_in = netCDF4.Dataset(fname_in, 'r')
    m_in = netCDF4.Dataset(mname_in, 'r')

    # get relationship between a cluster and particles (group particles into clusters)
    indexOnCluster, maxIndex, indexToParticle = get_clusters(m_in, cluster_path)
    
    Ntime = len(f_in.dimensions['Time'])
    Nclusters = maxIndex.shape[0]
    rEarth = f_in.sphere_radius

    mux = np.zeros((Ntime,Nclusters))
    muy = np.zeros((Ntime,Nclusters))
    dxdx_sum = np.zeros((Ntime,Nclusters))
    dxdy_sum = np.zeros((Ntime,Nclusters))
    dydy_sum = np.zeros((Ntime,Nclusters))
    drdr_sum = np.zeros((Ntime,Nclusters))
    Npart = np.zeros((Ntime,Nclusters))

    for t in np.arange(Ntime):
        print ('Looping over step ', t, ' of ', Ntime)
      
        x = f_in.variables['xParticle'][t]
        y = f_in.variables['yParticle'][t]
        z = f_in.variables['zParticle'][t]
        
        plat, plon = proj_lat_long(x,y,z)

        # compute center of mass and relative dispersion to this center of mass
        mux[t,:], muy[t,:], dxdx_sum[t,:], dxdy_sum[t,:], dydy_sum[t,:], drdr_sum[t,:], Npart[t,:] = \
                compute_com_and_dispersion(plat, plon, rEarth, indexOnCluster, maxIndex, indexToParticle)
   
    print ('saving numpy arrays')
    np.savez('lagrangian_data.npz',mux=mux, muy=muy, dxdx_sum=dxdx_sum, dxdy_sum=dxdy_sum, dydy_sum=dydy_sum, drdr_sum=drdr_sum, Npart=Npart)
    print ('done')

    f_in.close()
    m_in.close()

    return  #}}}

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension",
                      metavar="FILE")
    

    parser.add_option("-m", "--mesh", dest="meshfilename",
                      help="file to open for appending \
                      mesh data 'mesh_' extension",
                      metavar="FILE")

    parser.add_option("-c", "--load_cluser_info_path", dest="clusterinfopath",
                      help="location of files with cluster info", 
                      metavar="DIR")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")
    if not options.meshfilename:
        parser.error("Mesh filename is a required input.")
    if not options.clusterinfopath:
        options.clusterinfopath = None

    build_particle_file(options.inputfilename,options.meshfilename, options.clusterinfopath)
