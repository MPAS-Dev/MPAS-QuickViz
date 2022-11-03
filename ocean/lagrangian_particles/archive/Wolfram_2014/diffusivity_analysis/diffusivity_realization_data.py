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
import numpy as np
import netCDF4
import cPickle as pickle

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
  sin   = np.sin
  cos   = np.cos
  atan2 = np.arctan2
  sqrt  = np.sqrt

  dphi = phi2 - phi1
  dlam = lam2 - lam1

  return atan2(sin(dlam)*cos(phi2), cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(dlam)) #}}}

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
    # aliases #{{{
    sin   = np.sin
    cos   = np.cos
    atan2 = np.arctan2
    sqrt  = np.sqrt
    #}}}

    dphi = phi2 - phi1 
    dlam = lam2 - lam1
    
    a = sin(dphi/2.0)**2 + cos(phi1) * cos(phi2) * sin(dlam/2.0)**2
    c = 2.0 * atan2(sqrt(a), sqrt(1.0-a))

    # distance is c * Radius_of_earth

    return c #}}}

def proj_lat_long(x, y, z):  #{{{
    """
    compute the latitude and longitude from
    the x,y,z points (follow's Doug's pnt.h)
    """
    plat = np.arcsin(z / np.sqrt(x ** 2 + y ** 2 + z ** 2))
    
    plong = np.arctan2(y, x) 
    # adjust for periodicity
    plong_fix = (plong < -np.pi) * (plong + 2 * np.pi) + \
        (plong > np.pi) * (plong - 2 * np.pi)
    plong = (plong_fix == 0) * plong + (plong_fix != 0) * plong_fix

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

def get_clusters(f_in, cluster_path): #{{{
    """
    Compute or load the cluster information to link particles
    to clusters to compute dispersion.

    Phillip Wolfram, LANL
    09/19/2014
    """
    # compute interior cells and masks for vertices 
    if cluster_path is None:
      print 'Computing particle clusters'
      verticesOnCell = f_in.variables['verticesOnCell']
      nEdgesOnCell   = f_in.variables['nEdgesOnCell']
      cellsOnVertex  = f_in.variables['cellsOnVertex'] 
      boundaryVertex = (cellsOnVertex[:,:] == 0)
      boundaryVertex = boundaryVertex.max(axis=1)
      
      nCells = len(f_in.dimensions['nCells'])
      inCells = np.ones(nCells+1)
      inCells[cellsOnVertex[boundaryVertex,:]-1] = 0
      boundaryCells = (inCells == 0)
      inCells = np.where(inCells > 0)
      verticesOnCell = verticesOnCell[inCells[0],:] - 1
      maxVertices = nEdgesOnCell[inCells[0]]
      verticesToParticle = np.cumsum(1-boundaryVertex) - 1

      # store cluster neighbor information for filter
      print 'saving verticesOnCell, maxVertices, verticesToParticle'
      pickle.dump(verticesToParticle, open( "verticesToParticle.p", "wb"))
      pickle.dump(verticesOnCell, open( "verticesOnCell.p", "wb"))
      pickle.dump(maxVertices, open( "maxVertices.p", "wb"))
      print 'done'
    else:
      print 'laoding cluster info'
      verticesToParticle = pickle.load(open(cluster_path+"/verticesToParticle.p","rb"))
      verticesOnCell = pickle.load(open(cluster_path+"/verticesOnCell.p","rb"))
      maxVertices = pickle.load(open(cluster_path+"/maxVertices.p","rb"))
      print 'done'

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

    dl_list = []
    N = indexOnCluster.shape[0]
    clat = np.zeros(N)
    clon = np.zeros(N)
    
    print 'Computing COM and dispersion for ', N, ' clusters'
    for aCluster, maxInd in enumerate(maxIndex):
        # get points of cluster
        particles = indexToParticle[indexOnCluster[aCluster,0:maxInd]]
        pslat = plat[particles]
        pslon = plon[particles]
        
        # compute center of mass 
        clat[aCluster] = np.mean(pslat)
        clon[aCluster] = np.mean(pslon)

        # compute distances in m from lat / long 
        dx = r * normalized_haversine_formula(             0,     0, clon[aCluster], pslon)
        dy = r * normalized_haversine_formula(clat[aCluster], pslat,              0,     0)
        # fix orientation of points
        bearing = spherical_bearing(clat[aCluster], pslat, clon[aCluster], pslon)
        dx -= 2*dx*(bearing < 0)
        dy -= 2*dy*(abs(bearing) > np.pi/2.0)

        dl_list.append(np.array([dx,dy]))

    clat = np.array(clat)
    clon = np.array(clon)

    return clon, clat, dl_list  #}}}

def build_particle_file(fname_in, fname_outbase, cluster_path):  #{{{
    """
    Take existing netCDF4 input f_in and produce
    vtu files for f_outbase with all the particle data

    Phillip Wolfram 
    LANL
    07/15/2014
    """
    # load the netCDF database
    f_in = netCDF4.Dataset(fname_in, 'r')
    
    Ntime = len(f_in.dimensions['Time'])
    rEarth         = f_in.sphere_radius

    # get relationship between a cluster and particles (group particles into clusters)
    indexOnCluster, maxIndex, indexToParticle = get_clusters(f_in, cluster_path)
   
    dl_time = []
    com_time = []
    for t in np.arange(Ntime):
        print 'Looping over step ', t, ' of ', Ntime
      
        x = f_in.variables['xParticle'][t]
        y = f_in.variables['yParticle'][t]
        z = f_in.variables['zParticle'][t]
        
        plat, plon = proj_lat_long(x,y,z)

        # compute center of mass and relative dispersion to this center of mass
        clon, clat, dl_list = compute_com_and_dispersion(plat, plon, rEarth, indexOnCluster, maxIndex, indexToParticle)

        # store data for each time in list 
        # dl is Ntime, nCluster, nParticlesInCluster where nParticlesInCluster is unstructured)
        # com is Ntime, nCluster, R2, all structured
        dl_time.append(dl_list)
        com_time.append(np.array([clon,clat]))

    # save data
    pickle.dump(dl_time, open( "dl_time.p", "wb"))
    pickle.dump(com_time, open( "com_time.p", "wb"))

    return  #}}}

if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension",
                      metavar="FILE")
    parser.add_option("-c", "--load_cluser_info_path", dest="clusterinfopath",
                      help="location of files with cluster info", 
                      metavar="DIR")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")
    if not options.clusterinfopath:
        options.clusterinfopath = None
    #}}}

    build_particle_file(options.inputfilename, 'particles', options.clusterinfopath) 

