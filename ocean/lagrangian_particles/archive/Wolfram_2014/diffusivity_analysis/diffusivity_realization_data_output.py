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
import cPickle as pickle

# function definitions

def compute_com_and_dispersion(plat, plon, r, indexOnCluster, maxIndex, indexToParticle): #{{{
    """
    compute particle disperison (2nd moment) for cluster, using lat / lon for basis of 
    calculation to set coordinates for Kxx,Kxy,Kyy where x is zonal and y meridional

    dispersion units in m^2

    Phillip Wolfram
    LANL
    07/18/2014
    """

    N = indexOnCluster.shape[0]
    clat = np.zeros(N)
    clon = np.zeros(N)
    sigma2xx = np.zeros(N)
    sigma2xy = np.zeros(N)
    sigma2yy = np.zeros(N)
    sigma2rr = np.zeros(N)
    
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

        Nparticles = len(particles)
        sigma2xx[aCluster] = sum(dx*dx)/(Nparticles-1)
        sigma2xy[aCluster] = sum(dx*dy)/(Nparticles-1)
        sigma2yy[aCluster] = sum(dy*dy)/(Nparticles-1)
        sigma2rr[aCluster] = sum(dx*dx + dy*dy)/(Nparticles-1)

    return clon, clat, sigma2xx, sigma2xy, sigma2yy, sigma2rr #}}}

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

def diff_interval_start_end(disp, tstart, tend, tint): #{{{
    """
    Compute the diffusivity assuming that time output frequency is in days so that tstart, tend, tint 
    are all in days.

    Phillip Wolfram
    LANL
    09/24/2014
    """
    return 0.5*(disp[tend] - disp[tstart])/((tend-tstart)*tint*86400.) #}}}

def compute_diffusivity(sigma2xx, sigma2xy, sigma2yy, sigma2rr, ts, te, tint): #{{{
    kappa_xx = diff_interval_start_end(sigma2xx, ts, te, tint)
    kappa_xy = diff_interval_start_end(sigma2xy, ts, te, tint)
    kappa_yy = diff_interval_start_end(sigma2yy, ts, te, tint)
    kappa_rr = diff_interval_start_end(sigma2rr, ts, te, tint)

    return kappa_rr, kappa_xx, kappa_xy, kappa_yy  #}}}

def output_particle_file(folder, mux, muy, kappa, kappa_xx, kappa_xy, kappa_yy, sigma2_xx, sigma2_xy, sigma2_yy, sigma2_rr, ts, te): #{{{
    """
    Write output to file for use with Fortran kernel processing.
    Phillip Wolfram
    LANL
    09/24/2014
    """
    if not os.path.isdir(folder):
        os.mkdir(folder)

    Nclusters = len(kappa)
    npoints = 0
    for acluster in np.arange(Nclusters):
      for atime in np.arange(ts[acluster],te[acluster]+1):
        npoints +=1
        
    # allocations #{{{
    kappa_tot = np.zeros(npoints)
    kappa_xx_tot = np.zeros(npoints)
    kappa_xy_tot = np.zeros(npoints)
    kappa_yy_tot = np.zeros(npoints)
    x_tot     = np.zeros(npoints)
    y_tot     = np.zeros(npoints)
    sigma2_xx_tot = np.zeros(npoints)
    sigma2_xy_tot = np.zeros(npoints)
    sigma2_yy_tot = np.zeros(npoints)
    sigma2_rr_tot = np.zeros(npoints)
    clusternum_tot = np.zeros(npoints)
    #}}}

    # build up data structure
    npoints = 0
    ncluster = 0
    for acluster in np.arange(Nclusters):
      ncluster += 1
      for atime in np.arange(ts[acluster],te[acluster]+1):
        kappa_tot[npoints] = kappa[acluster]
        kappa_xx_tot[npoints] = kappa_xx[acluster]
        kappa_xy_tot[npoints] = kappa_xy[acluster]
        kappa_yy_tot[npoints] = kappa_yy[acluster]
        x_tot[npoints]     = mux[atime, acluster]
        y_tot[npoints]     = muy[atime, acluster]
        sigma2_xx_tot[npoints] = sigma2_xx[atime, acluster]
        sigma2_xy_tot[npoints] = sigma2_xy[atime, acluster]
        sigma2_yy_tot[npoints] = sigma2_yy[atime, acluster]
        sigma2_rr_tot[npoints] = sigma2_rr[atime, acluster]
        clusternum_tot[npoints] = ncluster
        npoints+=1

    # write ascii files for analysis in fortran
    print 'writing ascii data to ', folder
    with file(folder + '/cluster_data.txt', 'w') as outfile:
      # want the total number of data points
      outfile.write('%d\n' % x_tot.shape[0])
      for n,x,y,k,kxx,kxy,kyy,sxx,sxy,syy, srr in zip(clusternum_tot, x_tot, y_tot, kappa_tot, kappa_xx_tot, kappa_xy_tot, kappa_yy_tot, sigma2_xx_tot, sigma2_xy_tot, sigma2_yy_tot, sigma2_rr_tot):
        outfile.write('%d %e %e %e %e %e %e %e %e %e %e\n' % (n,x,y,k,kxx,kxy,kyy,sxx,sxy,syy,srr))
    with file(folder + '/interp_points.txt', 'w') as outfile:
      outfile.write('%d\n' % len(mux[0,:]))
      for x,y in zip(mux[0,:],muy[0,:]):
        outfile.write('%e %e\n' % (x,y))
    print 'done'

    print 'saving numpy arrays'
    np.savez(folder+'/cluster_data.npz',n=clusternum_tot,x=x_tot,y=y_tot,\
        krr=kappa_tot,kxx=kappa_xx_tot,kxy=kappa_xy_tot,kyy=kappa_yy_tot,\
        s2xx=sigma2_xx_tot,s2xy=sigma2_xy_tot,s2yy=sigma2_yy_tot,s2rr=sigma2_rr_tot)
    np.savez(folder+'/interp_points.npz',x=mux[0,:],y=muy[0,:])
    np.save(folder+'/n.npy',clusternum_tot)
    np.save(folder+'/x.npy',x_tot)
    np.save(folder+'/y.npy',y_tot)
    np.save(folder+'/krr.npy',kappa_tot)
    np.save(folder+'/kxx.npy',kappa_xx_tot)
    np.save(folder+'/kxy.npy',kappa_xy_tot)
    np.save(folder+'/kyy.npy',kappa_yy_tot)
    np.save(folder+'/s2xx.npy',sigma2_xx_tot)
    np.save(folder+'/s2xy.npy',sigma2_xy_tot)
    np.save(folder+'/s2yy.npy',sigma2_yy_tot)
    np.save(folder+'/s2rr.npy',sigma2_rr_tot)
    np.save(folder+'/xi.npy',mux[0,:])
    np.save(folder+'/yi.npy',muy[0,:])
    print 'done'

    return  #}}}

def build_particle_file(fname_in, fname_outbase, cluster_path, outputfolder, ts, te, tint):  #{{{
    """
    Take existing netCDF4 input f_in and produce
    vtu files for f_outbase with all the particle data

    Phillip Wolfram 
    LANL
    07/15/2014
    """
    # load the netCDF database
    f_in = netCDF4.Dataset(fname_in, 'r')

    # get relationship between a cluster and particles (group particles into clusters)
    indexOnCluster, maxIndex, indexToParticle = get_clusters(f_in, cluster_path)
    
    Ntime = len(f_in.dimensions['Time'])
    #Nparticles = len(f_in.dimensions['nParticles'])
    Nclusters = indexOnCluster.shape[0]
    rEarth         = f_in.sphere_radius

    mux = np.zeros((Ntime,Nclusters))
    muy = np.zeros((Ntime,Nclusters))
    sigma2xx = np.zeros((Ntime,Nclusters))
    sigma2xy = np.zeros((Ntime,Nclusters))
    sigma2yy = np.zeros((Ntime,Nclusters))
    sigma2rr = np.zeros((Ntime,Nclusters))
    for t in np.arange(Ntime):
        print 'Looping over step ', t, ' of ', Ntime
      
        x = f_in.variables['xParticle'][t]
        y = f_in.variables['yParticle'][t]
        z = f_in.variables['zParticle'][t]
        
        plat, plon = proj_lat_long(x,y,z)

        # compute center of mass and relative dispersion to this center of mass
        mux[t,:], muy[t,:], sigma2xx[t,:], sigma2xy[t,:], sigma2yy[t,:], sigma2rr[t,:] = \
                compute_com_and_dispersion(plat, plon, rEarth, indexOnCluster, maxIndex, indexToParticle)
    
    kappa, kappa_xx, kappa_xy, kappa_yy = compute_diffusivity(sigma2xx, sigma2xy, sigma2yy, sigma2rr, ts, te, tint)

    e = np.ones(kappa.shape)
    output_particle_file(outputfolder, mux, muy, kappa, kappa_xx, kappa_xy, kappa_yy, sigma2xx, sigma2xy, sigma2yy, sigma2rr, ts*e, te*e)

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
    parser.add_option("-o", "--output_file", dest="outputfolder",
                      help="location of outputfolder",
                      metavar="FILE")
    parser.add_option("-e", "--end", dest="end",
                      help="step for ending time",
                      metavar="INTEGER")
    parser.add_option("-s", "--start", dest="start",
                      help="step for starting time",
                      metavar="INTEGER")
    parser.add_option("-i", "--interval", dest="interval",
                      help="time step interval",
                      metavar="INTEGER")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")
    if not options.clusterinfopath:
        options.clusterinfopath = None
    if not options.outputfolder:
        options.outputfolder = 'particles_statistics/'
    if not options.end:
        parser.error("Ending time step is a required input.")
    else:
        options.end= int(options.end)
    if not options.start:
        parser.error("Starting time step is a required input.")
    else:
        options.start = int(options.start)
    if not options.interval:
        parser.error("Time interval is a required input.")
    else:
        options.interval = int(options.interval)
    #}}}

    build_particle_file(options.inputfilename, 'particles', options.clusterinfopath, options.outputfolder, 
            options.start, options.end, options.interval) 

