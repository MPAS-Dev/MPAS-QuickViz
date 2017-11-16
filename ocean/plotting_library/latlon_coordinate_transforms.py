#!/usr/bin/env python
"""

    Coordinate transforms and distance calculations.
    Numpy and numexpr versions (using MKL) are available.
    Use numexpr for larger array sizes, numpy for smaller
    array sizes.

    Phillip Wolfram
    LANL
    03/13/2015

"""
# import libraries / packages
import numpy as np
import numexpr as ne
#from xray.ufuncs import fabs

rEarth = 6371220.  # from SOMA case

# function definitions

def test_signed_distances(): #{{{
    import matplotlib.pyplot as plt
    def test(x1,y1,x2,y2):
        sc = 0.05
        plt.figure()
        plt.hold(True)
        ax = plt.axes()
        ax.plot(x1,y1,'r.')
        ax.plot(x2,y2,'rx')
        def adjust_lims(sc): #{{{
            xr = np.abs(x2-x1)
            xmin, xmax = plt.xlim()
            plt.xlim(xmin-sc*xr,xmax+sc*xr)
            yr = np.abs(y2-y1)
            ymin, ymax = plt.ylim()
            plt.ylim(ymin-sc*yr,ymax+sc*yr)
            return sc*np.maximum(xr,yr)#}}}
        h  = adjust_lims(sc)
        ax.arrow(x1,y1,x2-x1,y2-y1, head_width=h, head_length=h)
        ax.set_aspect('equal')
        plt.xlabel('lon')
        plt.ylabel('lat')

    def test_distances(x1,y1,x2,y2):
        x1 = np.radians(x1)
        x2 = np.radians(x2)
        y1 = np.radians(y1)
        y2 = np.radians(y2)
        test(np.degrees(x1),np.degrees(y1),np.degrees(x2),np.degrees(y2))
        dxf, dyf = signed_distances_numexpr(y1,y2,x1,x2)
        dx, dy = signed_distances_numpy(y1,y2,x1,x2)
        print 'numpy   = ', dx,dy
        print 'numexpr = ', dxf, dyf
        print 'diff    = ', dx-dxf, dy-dyf
        plt.show()
    test_distances(0.0,0.0,90.,0.)
    test_distances(0.0,0.0,-90.,0.)
    test_distances(0.0,0.0,0.0,45.)
    test_distances(35.,30.,40.,0.)
    test_distances(0.,0.,-10.,10.)
    test_distances(0,0,5,10)
    test_distances(15.,30.,5.,10.)
    return #}}}

def signed_distances_numexpr(phi1, phi2, lam1, lam2, r=rEarth):  #{{{
    """
    Return signed dx and dy based on bearings where lambda is longitude and phi is latitude

    Phillip Wolfram
    LANL
    03/13/2015
    """
    dx = haversine_formula_numexpr(phi1, phi1, lam1, lam2, r)
    dy = haversine_formula_numexpr(phi1, phi2, lam1, lam1, r)

    # fix orientation of points
    bearing = spherical_bearing_numexpr(phi1, phi2, lam1, lam2)
    # because arctan2 returns results from -pi to pi for bearing, flip values to get right sign
    PI = np.pi
    dx -= ne.evaluate("2*dx*(bearing < 0)")
    dy -= ne.evaluate("2*dy*(abs(bearing) > PI/2.0)")

    return dx, dy #}}}

def signed_distances_numpy(phi1, phi2, lam1, lam2, r=rEarth):  #{{{
    """
    Return signed dx and dy based on bearings where lambda is longitude and phi is latitude

    Phillip Wolfram
    LANL
    03/13/2015
    """

    dx = haversine_formula_numpy(phi1, phi1, lam1, lam2, r)
    dy = haversine_formula_numpy(phi1, phi2, lam1, lam1, r)

    # fix orientation of points
    bearing = spherical_bearing_numpy(phi1, phi2, lam1, lam2)
    # because arctan2 returns results from -pi to pi for bearing, flip values to get right sign
    dx -= 2*dx*(bearing < 0)
    dy -= 2*dy*(abs(bearing) > np.pi/2.0)

    return dx, dy #}}}

def spherical_bearing_numexpr(phi1, phi2, lam1, lam2): #{{{
  """
  compute the spherical bearing
      http://www.movable-type.co.uk/scripts/latlong.html

      where lambda is longitude and phi is latitude

      spherical_bearing returned is from N (-pi, pi)

      Phillip Wolfram
      LANL
      08/07/2014
  """
  return ne.evaluate("arctan2(sin(lam2-lam1)*cos(phi2), cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(lam2-lam1))") #}}}

def spherical_bearing_numpy(phi1, phi2, lam1, lam2): #{{{
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

def haversine_formula_numexpr(phi1, phi2, lam1, lam2, r=rEarth):  #{{{
    """
    compute the distance between two points via Haversine formula:
    http://www.movable-type.co.uk/scripts/latlong.html

    where lambda is longitude and phi is latitude

    c returned in units of R (or radians if r=1.0)

    Phillip Wolfram
    LANL
    07/18/2014
    """
    a = ne.evaluate("sin((phi2-phi1)/2.0)**2 + cos(phi1) * cos(phi2) * sin((lam2-lam1)/2.0)**2")
    c = ne.evaluate("r * 2.0 * arctan2(sqrt(a), sqrt(1.0-a))")

    return c #}}}

def haversine_formula_numpy(phi1, phi2, lam1, lam2, r=rEarth):  #{{{
    """
    compute the distance between two points via Haversine formula:
    http://www.movable-type.co.uk/scripts/latlong.html

    where lambda is longitude and phi is latitude

    c returned in units of R (or radians if r=1.0)

    Phillip Wolfram
    LANL
    07/18/2014
    """
    dlam = lam2-lam1
    dphi = phi2-phi1

    a = np.sin(dphi/2.0)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlam/2.0)**2
    c = r * 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0-a))

    return c #}}}

def proj_lat_long_numexpr(x, y, z):  #{{{
    """
    compute the latitude and longitude from
    the x,y,z points (follow's Doug's pnt.h)
    """
    plat = ne.evaluate("arcsin(z / sqrt(x ** 2 + y ** 2 + z ** 2))")
    plon = ne.evaluate("arctan2(y, x)")

    return plat, plon  #}}}

def proj_lat_long_numpy(x, y, z):  #{{{
    """
    compute the latitude and longitude from
    the x,y,z points (follow's Doug's pnt.h)
    """
    plat = np.arcsin(z / np.sqrt(x ** 2 + y ** 2 + z ** 2))
    plong = np.arctan2(y, x)

    return plat, plong  #}}}

def proj_xyz_numexpr(plat, plong, r=rEarth):  #{{{
    """
    convert from latitude / longitude to xyz
    Phillip Wolfram
    LANL
    07/17/2014
    """
    x = ne.evaluate("r * cos(plong) * cos(plat)")
    y = ne.evaluate("r * sin(plong) * cos(plat)")
    z = ne.evaluate("r * sin(plat)")

    return x,y,z  #}}}

def proj_xyz_numpy(plat, plong, r=rEarth):  #{{{
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

def mid_point_numpy(lat, lon, r=rEarth): #{{{
    """
    compute the midpoint of of latitude and longitude for
    Np,Nsamples across Nsamples
    Phillip Wolfram
    LANL
    05/13/2015
    """
    x,y,z = proj_xyz_numpy(lat, lon, r)
    latm, lonm = proj_lat_long_numpy(np.mean(x,axis=-1), np.mean(y,axis=-1), np.mean(z,axis=-1))
    return latm, lonm #}}}

def mid_point_numexpr(lat, lon, r=rEarth, ensmean=(-1)): #{{{
    """
    compute the midpoint of of latitude and longitude for
    Np,Nsamples across Nsamples
    Phillip Wolfram
    LANL
    05/13/2015
    """
    x,y,z = proj_xyz_numexpr(lat, lon, r)
    latm, lonm = proj_lat_long_numexpr(np.mean(x,axis=ensmean), np.mean(y,axis=ensmean), np.mean(z,axis=ensmean))
    return latm, lonm #}}}

def test_midpoint(): #{{{
    # test the midpoint calculation

    def test(lat1,lon1,lat2,lon2, r=1):
        print 'point 1 =', lon1,lat1
        print 'point 2 =', lon2,lat2
        latm1, lonm1 = mid_point_numpy(np.vstack((lat1,lat2)).T,np.vstack((lon1,lon2)).T,r)
        latm2, lonm2 = mid_point_numexpr(np.vstack((lat1,lat2)).T,np.vstack((lon1,lon2)).T,r)
        print 'numpy   = ', lonm1, latm1
        print 'numexpr = ', lonm2, latm2
        print 'avg = ', 0.5*(lon1+lon2), 0.5*(lat1+lat2)
        print 'diff    = ', lonm1-lonm2, latm1-latm2

    test(np.radians(35.0),np.radians(101.),np.radians(40.0), np.radians(135.0))
    return #}}}

#def fix_periodicity_xray(px, xc, L): #{{{
#    """ fix periodicity similar to in mpas_vector_operations """
#    dist = px - xc
#    fix = fabs(dist) > L/2.0
#    pfix = px - fix*(dist/np.maximum(1e-15,np.abs(dist)))*L
#    return pfix #}}}

def fix_periodicity(px, xc, L): #{{{
    """ fix periodicity similar to in mpas_vector_operations """
    dist = px - xc
    fix = np.abs(dist) > L/2.0
    pfix = px - fix*(dist/np.abs(dist))*L
    idx = np.where(dist==0)
    pfix[idx] = (px*np.ones_like(pfix))[idx]
    return pfix #}}}

def fix_periodicity_numexpr(px, xc, L): #{{{
    """ fix periodicity similar to in mpas_vector_operations """
    pfix = ne.evaluate('px - (abs(px-xc) > L/2.0)*((px-xc)/abs(px-xc))*L')
    idx = np.where((px-xc) == 0)
    pfix[idx] = (px*np.ones_like(pfix))[idx]
    return pfix #}}}

def fix_periodic_timeseries(ts, L): #{{{
    """ fix periodic time series ts[times, values], assuming coherent dataset for times=0 (continguous data) """
    x0 = np.zeros((ts.shape[1]))
    for at in np.arange(ts.shape[0]-1)+1:
        ts[at,:] += x0
        dx = np.zeros((ts.shape[1]))
        dx += ((ts[at,:] - ts[at-1,:]) < -L/2.0)*L
        dx -= ((ts[at,:] - ts[at-1,:]) >  L/2.0)*L
        ts[at,:] += dx
        x0 += dx
    return ts #}}}

def fix_periodic_timeseries_dataarray(ts, L): #{{{
  """ Operates on a DataArray """
  var = ts.values
  var = fix_periodic_timeseries(var, L)
  ts[:] = var[:]
  return ts #}}}

def test_fix_periodicity(fix_periodicity_tester=fix_periodicity_numexpr): #{{{
    import matplotlib.pyplot as plt
    x = 2.0*np.pi*np.random.rand(100) - np.pi
    plt.figure()
    plt.subplot(1,2,1)
    plt.hist(x)
    x = fix_periodicity_tester(x, np.pi, 2.0*np.pi)
    plt.subplot(1,2,2)
    plt.hist(x)
    plt.show()
    return #}}}

def test_fix_periodicity_timeseries(N=100, L=2.0*np.pi): #{{{
    print "note that this test doesn't fully demonstrate robustness"
    import matplotlib.pyplot as plt
    x1 = np.linspace(0,3*L,N)
    x1p = np.mod(x1,L)
    x1c = fix_periodic_timeseries(np.vstack((x1p, x1p)).T, L)
    assert np.max(x1-x1c[:,0]) == 0 and np.min(x1c[:,0] == x1c[:,1]), 'Correction failed for periodic_timeseries'
    plt.figure()
    plt.plot(x1,np.ones((N)),'k.')
    plt.plot(x1p,2*np.ones((N)),'r.')
    plt.plot(x1c[:,0],3*np.ones((N)),'b.')
    plt.ylim([0,4])
    plt.show()
    x2 = np.pi*np.linspace(0,-3*L,N)
    x2p = np.mod(x2,L)
    x2c = fix_periodic_timeseries(np.vstack((x2p, x2p)).T, L)
    assert np.max(x2-x2c[:,0]) == 0 and np.min(x2c[:,0] == x2c[:,1]), 'Correction failed for periodic_timeseries'
    plt.figure()
    plt.plot(x2,np.ones((N)),'k.')
    plt.plot(x2p,2*np.ones((N)),'r.')
    plt.plot(x2c[:,0],3*np.ones((N)),'b.')
    plt.ylim([0,4])
    plt.show()
    return #}}}

if __name__ == "__main__":
    print 'Module to perform coordinate transforms, espeically for xyz to lat/lon and on lat/lon.'
    print 'Running tests'

    #test_signed_distances()
    #test_midpoint()
    #test_fix_periodicity()
    test_fix_periodicity_timeseries()

