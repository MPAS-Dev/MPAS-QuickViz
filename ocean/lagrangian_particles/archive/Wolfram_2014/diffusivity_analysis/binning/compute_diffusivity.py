#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from latlon_coordinate_transforms import signed_distances_numpy as signed_distances, signed_distances_numexpr, haversine_formula_numexpr, haversine_formula_numpy as haversine_formula, mid_point_numpy, mid_point_numexpr, spherical_bearing_numexpr as spherical_bearing
#import seaborn as sns
import scipy.stats as stats
rEarth = 6371220. # from netcdf file
# from http://stackoverflow.com/questions/16003217/n-d-version-of-itertools-combinations-in-numpy
from itertools import combinations, chain
from scipy.misc import comb
from fast_hist import fast_hist_weighted

def comb_index(n, k): #{{{
    count = comb(n, k, exact=True)
    index = np.fromiter(chain.from_iterable(combinations(range(n), k)),
                        int, count=count*k)
    return index.reshape(-1, k) #}}}

# function definitions

def hull_particle_pairs(plat, plon): #{{{
    """
    Compute particle pairs on convex hull of cluster of particles
    (plon,plat) for nparticles.  Return pairparticles (npairs,2)
    data structure which contains pairs of particles on the convex
    hull boundary.
    """
    from scipy.spatial import ConvexHull
    hull = ConvexHull(np.vstack((plon,plat)).T)
    return hull.simplices #}}}

def test_hull_particle_pairs(npart=100): #{{{
    print 'Testing ability to pick out particle pairs on boundary ... ',
    lon = np.random.rand(npart)
    lat = np.random.rand(npart)
    pairparticles = hull_particle_pairs(lat,lon)

    # plot for debugging purposes
    plt.figure()
    plt.plot(lon,lat,'.')
    for i in np.arange(pairparticles.shape[0]):
        plt.plot(lon[pairparticles[i,:]],lat[pairparticles[i,:]],'r-.')
        plt.text(np.mean(lon[pairparticles[i,:]]),np.mean(lat[pairparticles[i,:]]), str(i))
    plt.show()
    print 'finished.'
    return #}}}

def delaunay_particle_pairs(plat,plon): #{{{
    from matplotlib.tri import Triangulation
    tri = Triangulation(plon, plat)
    return tri.edges #}}}

def test_delaunay_particle_pairs(npart=100): #{{{
    print 'Testing ability to pick out particle pairs on boundary ... ',
    lon = np.random.rand(npart)
    lat = np.random.rand(npart)
    pairparticles = delaunay_particle_pairs(lat,lon)

    # plot for debugging purposes
    plt.figure()
    plt.plot(lon,lat,'r.')
    for i in np.arange(pairparticles.shape[0]):
        plt.plot(lon[pairparticles[i,:]],lat[pairparticles[i,:]],'k-', alpha=0.3)
    plt.show()
    print 'finished.'
    return #}}}

def particle_pair_pdfs_simple(plat, plon, lvec, deltat, r=rEarth, pairselection='all'): #{{{
    """
    Compute particle pair statistics for input of size (nparticles, Ntime, Nensembles)
    from (plat,plon) (nparticles,2) data corresponding to frequency (in days) of deltat
    bin on a histogram with bin-edges of lvec for particle pair separations.

    If pairselection designates which class of particle pairs to use for calculation.

    xm,ym is mean pair location
    dx,dy is coordinate distance between pairs
    l is distance between pairs


    Phillip Wolfram
    LANL
    07/14/2015
    """
    def align(A,B):
        return np.concatenate((A[:,:,:,np.newaxis],B[:,:,:,np.newaxis]),axis=3)

    def time_deriv(var):
        return np.diff(var,axis=1)/(24.*60.*60.*deltat)

    def hist(x,c):
        assert not np.any(x > lvec[-1]), 'max length = %f outside max bin edge of %f'%(np.max(x),lvec[-1])
        Np = x.shape[0]
        Nt = x.shape[1]
        Ne = x.shape[2]
        pdf = np.zeros((lvec.shape[0]-1, Nt))
        for i in np.arange(Nt):
          pdf[:,i] = fast_hist_weighted(x[:,i,:].ravel(), c[:,i,:].ravel(), lvec)
        return pdf

    if pairselection == 'hull':
        # only consider starting positions (same for starting time and all ensembles)
        pairparticles = hull_particle_pairs(plat[:,0,0],plon[:,0,0])
        # this has to use measurements on sphere
        lona = np.squeeze(plon[pairparticles[:,0],:,:])
        lata = np.squeeze(plat[pairparticles[:,0],:,:])
        lonb = np.squeeze(plon[pairparticles[:,1],:,:])
        latb = np.squeeze(plat[pairparticles[:,1],:,:])
        # plot the cluster boundaries and particles
        if False:
            # these are size (npairs, Ntime, Nensembles)
            # plot the rings
            def draw_ensemble(ens, color):
                fig = plt.figure()
                plt.ion()
                for i in np.arange(lona.shape[1]):
                    fig.clf()
                    plt.plot(plon[:,i,ens],plat[:,i,ens],color+'.')
                    x = np.vstack((lona[:,i,ens].T.ravel(),lonb[:,i,ens].T.ravel())).T
                    y = np.vstack((lata[:,i,ens].T.ravel(),latb[:,i,ens].T.ravel())).T
                    for j in np.arange(x.shape[0]):
                        plt.plot(x[j,:],y[j,:],color+'-',alpha=1.0-float(i)/float(x.shape[0]), lw = 1.0 + float(i)/float(x.shape[0]))
                    plt.title('ens = %d, i=%d'%(ens, i))
                    plt.draw()
                    plt.pause(1.0)
            for k in np.arange(plon.shape[2]):
                draw_ensemble(k,'k')
    elif pairselection == 'nearest':
        pairparticles = delaunay_particle_pairs(plat[:,0,0],plon[:,0,0])
        # this has to use measurements on sphere
        lona = np.squeeze(plon[pairparticles[:,0],:,:])
        lata = np.squeeze(plat[pairparticles[:,0],:,:])
        lonb = np.squeeze(plon[pairparticles[:,1],:,:])
        latb = np.squeeze(plat[pairparticles[:,1],:,:])
        if False:
            # these are size (npairs, Ntime, Nensembles)
            # plot the rings
            def draw_ensemble(ens, color):
                fig = plt.figure()
                plt.ion()
                for i in np.arange(lona.shape[1]):
                    fig.clf()
                    plt.plot(plon[:,i,ens],plat[:,i,ens],color+'.')
                    x = np.vstack((lona[:,i,ens].T.ravel(),lonb[:,i,ens].T.ravel())).T
                    y = np.vstack((lata[:,i,ens].T.ravel(),latb[:,i,ens].T.ravel())).T
                    for j in np.arange(x.shape[0]):
                        plt.plot(x[j,:],y[j,:],color+'-',alpha=0.1)
                    plt.title('ens = %d, i=%d'%(ens, i))
                    plt.draw()
                    plt.pause(1.0)
            for k in np.arange(plon.shape[2]):
                draw_ensemble(k,'k')

    elif pairselection == 'all':
        nparticles = plat.shape[0]
        pairs = comb_index(nparticles,2)
        npairs = pairs.shape[0]
        # this has to use measurements on sphere
        lona = np.squeeze(plon[pairs[:,0],:,:])
        lata = np.squeeze(plat[pairs[:,0],:,:])
        lonb = np.squeeze(plon[pairs[:,1],:,:])
        latb = np.squeeze(plat[pairs[:,1],:,:])
        # these are size (npairs, Ntime, Nensembles)
        if False:
            # these are size (npairs, Ntime, Nensembles)
            # plot the rings
            def draw_ensemble(ens, color):
                fig = plt.figure()
                plt.ion()
                for i in np.arange(lona.shape[1]):
                    fig.clf()
                    plt.plot(plon[:,i,ens],plat[:,i,ens],color+'.')
                    x = np.vstack((lona[:,i,ens].T.ravel(),lonb[:,i,ens].T.ravel())).T
                    y = np.vstack((lata[:,i,ens].T.ravel(),latb[:,i,ens].T.ravel())).T
                    for j in np.arange(x.shape[0]):
                        plt.plot(x[j,:],y[j,:],color+'-',alpha=0.1)
                    plt.title('ens = %d, i=%d'%(ens, i))
                    plt.draw()
                    plt.pause(1.0)
            for k in np.arange(plon.shape[2]):
                draw_ensemble(k,'k')

    else:
        assert False, 'pairselection = %s is not a known choice'%(pairselection)

    l = haversine_formula_numexpr(lata, latb, lona, lonb, r)
    # length in km
    lkm = 0.5*(l[:,1:,:]+l[:,0:-1,:])/1000.

    # diffusivity (note variance is half the cluster variance, e.g., LaCasce, 2008))
    # division by two accounts for full lengths used instead of half lengths with the cluster
    Krr = 0.5*time_deriv(0.5*l*l)

    # get number of pair lengths in bins
    lmPDF   = hist(lkm, np.ones_like(lkm))
    # compute mean values within each bin corresponding to particular scales
    KrrPDF  = hist(lkm, Krr)/lmPDF

    return lmPDF, KrrPDF #}}}

def particle_pair_pdfs(plat, plon, lvec, deltat, r=rEarth): #{{{
    """
    Compute particle pair statistics for input of size (nparticles, Ntime, Nensembles)
    from (plat,plon) (nparticles,2) data corresponding to frequency (in days) of deltat
    bin on a histogram with bin-edges of lvec for particle pair separations.

    xm,ym is mean pair location
    dx,dy is coordinate distance between pairs
    l is distance between pairs


    Phillip Wolfram
    LANL
    05/13/2015
    """
    nparticles = plat.shape[0]
    pairs = comb_index(nparticles,2)
    npairs = pairs.shape[0]
    # this has to use measurements on sphere
    lona = np.squeeze(plon[pairs[:,0],:,:])
    lata = np.squeeze(plat[pairs[:,0],:,:])
    lonb = np.squeeze(plon[pairs[:,1],:,:])
    latb = np.squeeze(plat[pairs[:,1],:,:])
    # these are size (npairs, Ntime, Nensembles)

    def align(A,B):
        return np.concatenate((A[:,:,:,np.newaxis],B[:,:,:,np.newaxis]),axis=3)

    xm, ym = mid_point_numexpr(align(lata,latb),align(lona,lonb), r)
    dx, dy = signed_distances_numexpr(lata, latb, lona, lonb, r)
    l     = haversine_formula_numexpr(lata, latb, lona, lonb, r)
    # length in km
    lkm = l/1000.

    # bearing with map to 0,pi from -pi,pi
    theta = spherical_bearing(lata, latb, lona, lonb)
    theta += (theta < 0)*np.pi

    def time_deriv(var):
        return np.diff(var,axis=1)/(24.*60.*60.*deltat)

    # these are size (npairs, Ntime-1, Nensembles)
    lkmm = 0.5*(lkm[:,1:,:]+lkm[:,0:-1,:])
    # mean velocity
    um = time_deriv(xm)
    vm = time_deriv(ym)
    sm = np.sqrt(um*um + vm*vm)
    # eddy velocity
    ue = time_deriv(dx)
    ve = time_deriv(dy)
    se = np.sqrt(ue*ue + ve*ve)
    # diffusivity (note variance is half the cluster variance, e.g., LaCasce, 2008))
    # division by two accounts for full lengths used instead of half lengths with the cluster
    Kxx = 0.5*time_deriv(0.5*dx*dx)
    Kxy = 0.5*time_deriv(0.5*dx*dy)
    Kyy = 0.5*time_deriv(0.5*dy*dy)
    Krr = 0.5*time_deriv(0.5*l*l)

    # average over pairs
    KxxMean = np.mean(Kxx,axis=0)
    KxyMean = np.mean(Kxy,axis=0)
    KyyMean = np.mean(Kyy,axis=0)
    KrrMean = np.mean(Krr,axis=0)

    # build histograms
    def hist(x,c):
        assert not np.any(x > lvec[-1]), 'max length = %f outside max bin edge of %f'%(np.max(x),lvec[-1])
        Np = x.shape[0]
        Nt = x.shape[1]
        Ne = x.shape[2]
        pdf = np.zeros((lvec.shape[0]-1, Nt, Ne))
        for i in np.arange(Nt):
            for j in np.arange(Ne):
                pdf[:,i,j],_ = np.histogram(x[:,i,j], weights=c[:,i,j], bins=lvec, density=False)
        return pdf

    # get number of pair lengths in bins
    lPDF    = hist(lkm,  np.ones_like(lkm))
    lmPDF   = hist(lkmm, np.ones_like(lkmm))
    # compute mean values within each bin corresponding to particular scales
    xmPDF   = hist(lkm,  xm )/lPDF
    ymPDF   = hist(lkm,  ym )/lPDF
    dxPDF   = hist(lkm,  dx )/lPDF
    dyPDF   = hist(lkm,  dy )/lPDF
    thetaPDF= hist(lkm,  theta)/lPDF

    sePDF   = hist(lkmm, se )/lmPDF
    uePDF   = hist(lkmm, ue )/lmPDF
    vePDF   = hist(lkmm, ve )/lmPDF

    smPDF   = hist(lkmm, sm )/lmPDF
    umPDF   = hist(lkmm, um )/lmPDF
    vmPDF   = hist(lkmm, vm )/lmPDF

    KxxPDF  = hist(lkmm, Kxx)/lmPDF
    KxyPDF  = hist(lkmm, Kxy)/lmPDF
    KyyPDF  = hist(lkmm, Kyy)/lmPDF
    KrrPDF  = hist(lkmm, Krr)/lmPDF

    return nparticles, npairs, xmPDF, ymPDF, dxPDF, dyPDF, lPDF, lmPDF, sePDF, uePDF, vePDF, smPDF, umPDF, vmPDF, KxxPDF, KxyPDF, KyyPDF, KrrPDF, KxxMean, KxyMean, KyyMean, KrrMean, thetaPDF #}}}


def compute_dispersion_planar_numexpr(plat, plon, ensmean=(0,2)): #{{{
    """
    Compute relative particle disperison (2nd moment) for cluster, using lat / lon in radians for basis of
    calculation to set coordinates for Dxx,Dxy,Dyy where x is zonal and y meridional.  Vectorized to
    work over all times.  lat and lon size are [Np,Nt,Ne].  ensmean setting averages over particles and realization.

    dispersion units in m^2

    Phillip Wolfram
    LANL
    03/27/2015
    """
    nparticles = plat.shape[0]

    # compute center of mass for each time and realization
    clon = np.mean(plon, axis=0)[np.newaxis,:,:]
    clat = np.mean(plat, axis=0)[np.newaxis,:,:]

    # scalar diffusivity #{{{
    dr = np.sqrt((clon-plon)**2.0 + (clat-plat)**2.0)
    drdr = np.mean(dr*dr, axis=ensmean)
    #}}}

    # tensor diffusivity #{{{
    dx = plon - clon
    dy = plat - clat
    dxdx = np.mean(dx*dx, axis=ensmean)
    dxdy = np.mean(dx*dy, axis=ensmean)
    dydy = np.mean(dy*dy, axis=ensmean)
    #}}}

    # absolute diffusivity #{{{
    dxabs = plon - plon[:,0,][:,np.newaxis,:]
    dyabs = plat - plat[:,0,][:,np.newaxis,:]
    drabs = np.sqrt(dxabs*dxabs + dyabs*dyabs)

    dxdxabs = np.mean(dxabs*dxabs, axis=ensmean)
    dxdyabs = np.mean(dxabs*dyabs, axis=ensmean)
    dydyabs = np.mean(dyabs*dyabs, axis=ensmean)
    drdrabs = np.mean(drabs*drabs, axis=ensmean)
    #}}}

    return nparticles, clat[0,:,:], clon[0,:,:], drdr, dxdx, dxdy, dydy, drdrabs, dxdxabs, dxdyabs, dydyabs #}}}

def compute_dispersion_numexpr(plat, plon, r=rEarth, ensmean=(0,2)): #{{{
    """
    Compute relative particle disperison (2nd moment) for cluster, using lat / lon in radians for basis of
    calculation to set coordinates for Dxx,Dxy,Dyy where x is zonal and y meridional.  Vectorized to
    work over all times.  lat and lon size are [Np,Nt,Ne].  ensmean setting averages over particles and realization.

    dispersion units in m^2

    Phillip Wolfram
    LANL
    03/27/2015
    """
    nparticles = plat.shape[0]

    # compute center of mass for each time and realization
    #clat = np.mean(plat, axis=0)[np.newaxis,:,:]
    #clon = np.mean(plon, axis=0)[np.newaxis,:,:]
    clat, clon = mid_point_numexpr(plat, plon, r, ensmean=(0))
    clat = clat[np.newaxis,:,:]
    clon = clon[np.newaxis,:,:]
    slat = plat[:,0,][:,np.newaxis,:]
    slon = plon[:,0,][:,np.newaxis,:]

    # scalar diffusivity #{{{
    dr = haversine_formula_numexpr(clat, plat, clon, plon, r)
    drdr = np.mean(dr*dr, axis=ensmean)
    #}}}

    # tensor diffusivity #{{{
    dx, dy = signed_distances_numexpr(clat, plat, clon, plon, r)
    dxdx = np.mean(dx*dx, axis=ensmean)
    dxdy = np.mean(dx*dy, axis=ensmean)
    dydy = np.mean(dy*dy, axis=ensmean)
    #}}}

    # absolute diffusivity #{{{
    absdr = haversine_formula_numexpr(slat, plat, slon, plon, r)
    absdrdr = np.mean(absdr*absdr, axis=ensmean)
    absdx, absdy = signed_distances_numexpr(slat, plat, slon, plon, r)
    absdxdx = np.mean(absdx*absdx, axis=ensmean)
    absdxdy = np.mean(absdx*absdy, axis=ensmean)
    absdydy = np.mean(absdy*absdy, axis=ensmean)
    #}}}

    return nparticles, clat[0,:,:], clon[0,:,:], drdr, dxdx, dxdy, dydy, absdrdr, absdxdx, absdxdy, absdydy #}}}

def compute_com_and_squaresums(plat, plon, r=rEarth): #{{{
    """
    Compute relative particle disperison (2nd moment) for cluster, using lat / lon in radians for basis of
    calculation to set coordinates for Kxx,Kxy,Kyy where x is zonal and y meridional

    dispersion units in m^2

    Phillip Wolfram
    LANL
    07/18/2014
    """

    # compute center of mass for each time
    clat = np.mean(plat, axis=0)
    clon = np.mean(plon, axis=0)

    # scalar diffusivity #{{{
    dr = haversine_formula(clat, plat, clon, plon, r)
    drdr_sum = np.sum(dr*dr, axis=0)
    #}}}

    # tensor diffusivity #{{{
    dx, dy = signed_distances(clat, plat, clon, plon, r)
    dxdx_sum = np.sum(dx*dx, axis=0)
    dxdy_sum = np.sum(dx*dy, axis=0)
    dydy_sum = np.sum(dy*dy, axis=0)
    #}}}

    return clon, clat, drdr_sum, dxdx_sum, dxdy_sum, dydy_sum #}}}

def compute_diffusivity(lon1, lat1, lon2, lat2, deltat, r=rEarth): #{{{
    """
    Compute diffusivity using particle cluster covariances between two time periods.

    Diffusivity in units of m^2/s assuming that deltat is in s.

    Phillip Wolfram
    LANL
    03/23/2015
    """

    nparticles = lon1.shape[0]
    clon1, clat1, drdr_sum1, dxdx_sum1, dxdy_sum1, dydy_sum1 = compute_com_and_squaresums(lat1, lon1, r)
    clon2, clat2, drdr_sum2, dxdx_sum2, dxdy_sum2, dydy_sum2 = compute_com_and_squaresums(lat2, lon2, r)

    K_rr = 0.5*(drdr_sum2-drdr_sum1)/(nparticles-1)/deltat
    K_xx = 0.5*(dxdx_sum2-dxdx_sum1)/(nparticles-1)/deltat
    K_xy = 0.5*(dxdy_sum2-dxdy_sum1)/(nparticles-1)/deltat
    K_yy = 0.5*(dydy_sum2-dydy_sum1)/(nparticles-1)/deltat

    clon = 0.5*(clon1+clon2)
    clat = 0.5*(clat1+clat2)

    return nparticles, clon, clat, K_rr, K_xx, K_xy, K_yy #}}}

def compute_com_dispersion(lon, lat, r=rEarth): #{{{
    """
    Compute dispersion using particle clusters.

    Dispersion is in units of m^2.

    Phillip Wolfram
    LANL
    03/25/2015
    """

    nparticles = lon.shape[0]
    clon, clat, drdr_sum, dxdx_sum, dxdy_sum, dydy_sum = compute_com_and_squaresums(lat, lon, r)

    drdr = drdr_sum / nparticles
    dxdx = dxdx_sum / nparticles
    dxdy = dxdy_sum / nparticles
    dydy = dydy_sum / nparticles

    return nparticles, clon, clat, drdr, dxdx, dxdy, dydy #}}}

# tests
def run_tests(): #{{{
    print '................................................................................'
    print '. Tests '
    print '................................................................................'
    test_delaunay_particle_pairs()
    test_hull_particle_pairs()
    zonal_one_dim()
    merid_one_dim()
    background_const_zonal_flow()
    zonal_one_dim_random_walk()
    merid_one_dim_random_walk()
    two_dim_random_walk_iso()
    two_dim_random_walk_aniso()
    two_dim_random_walk_iso_time_evolution()

    plt.show()
    return #}}}

def merid_one_dim(N=30): #{{{
    print '================================================================================'
    print ' merid_one_dim'
    print '================================================================================'

    lat1 = np.radians(np.random.randn(N))
    print 'variance = %e'%(np.var(lat1))
    lat2 = np.radians(2*np.random.randn(N))
    lon1 = 0*lat1
    lon2 = 0*lat1

    nparticles, _, _, K_rr, K_xx, K_xy, K_yy = compute_diffusivity(lon1,lat1,lon2,lat2, 1.0, r=1.0)
    print 'nparticles=%d K_rr=%e K_xx=%e K_xy=%e K_yy=%e'%(nparticles, K_rr, K_xx, K_xy, K_yy)

    #}}}

def zonal_one_dim(N=30): #{{{
    print '================================================================================'
    print ' zonal_one_dim'
    print '================================================================================'

    lon1 = np.radians(np.random.randn(N))
    print 'variance = %e'%(np.var(lon1))
    lon2 = np.radians(2*np.random.randn(N))
    lat1 = 0*lon1
    lat2 = 0*lon1

    nparticles, _, _, K_rr, K_xx, K_xy, K_yy = compute_diffusivity(lon1,lat1,lon2,lat2, 1.0, r=1.0)
    print 'nparticles=%d K_rr=%e K_xx=%e K_xy=%e K_yy=%e'%(nparticles, K_rr, K_xx, K_xy, K_yy)

    #}}}

def background_const_zonal_flow(N=30): #{{{
    print '================================================================================'
    print ' background_const_zonal_flow'
    print '================================================================================'
    r = np.random.rand(N)*3.0
    theta = np.random.rand(N)*2*np.pi
    jx  = np.random.rand(N)
    jy  = np.random.rand(N)
    lon1 = np.radians(np.cos(theta)*r + jx)
    lat1 = np.radians(np.sin(theta)*r + jy)

    dx = np.random.rand(1)*3.0

    lon2 = lon1 + np.radians(dx)
    lat2 = lat1

    nparticles, _, _, K_rr, K_xx, K_xy, K_yy = compute_diffusivity(lon1,lat1,lon2,lat2, 1.0, r=1.0)
    print 'nparticles=%d K_rr=%e K_xx=%e K_xy=%e K_yy=%e'%(nparticles, K_rr, K_xx, K_xy, K_yy)

    # could assert that these are sufficiently small enough
    #}}}

def zonal_one_dim_random_walk(N=100000): #{{{
    print '================================================================================'
    print ' zonal one dimensional random walk'
    print '================================================================================'

    r = np.random.rand(N)*3.0
    theta = np.random.rand(N)*2*np.pi
    jx  = np.random.rand(N)
    jy  = np.random.rand(N)
    lon1 = np.radians(np.cos(theta)*r + jx)
    lat1 = 0*lon1

    deltat = 2*24*60*60
    K = 1e5
    print 'Analytical K_xx = %e'%(K)
    r = 1/3. # e.g., Visser1997, Gross et al POD report
    dx = (2*np.random.rand(N)-1)*np.sqrt(2/r*deltat*K/(rEarth**2.0))

    lon2 = lon1 + dx
    lat2 = 0*lon1

    nparticles, _, _, K_rr, K_xx, K_xy, K_yy = compute_diffusivity(lon1,lat1,lon2,lat2, deltat)
    print 'nparticles=%d K_rr=%e K_xx=%e K_xy=%e K_yy=%e'%(nparticles, K_rr, K_xx, K_xy, K_yy)

    #}}}

def merid_one_dim_random_walk(N=100000): #{{{
    print '================================================================================'
    print ' merid one dimensional random walk'
    print '================================================================================'

    r = np.random.rand(N)*3.0
    theta = np.random.rand(N)*2*np.pi
    jx  = np.random.rand(N)
    jy  = np.random.rand(N)
    lat1 = np.radians(np.cos(theta)*r + jx)
    lon1 = 0*lat1

    deltat = 2*24*60*60
    K = 1e5
    print 'Analytical K_yy = %e'%(K)
    r = 1/3. # e.g., Visser1997, Gross et al POD report
    dy = (2*np.random.rand(N)-1)*np.sqrt(2/r*deltat*K/(rEarth**2.0))

    lon2 = 0*lat1
    lat2 = lat1 + dy

    nparticles, _, _, K_rr, K_xx, K_xy, K_yy = compute_diffusivity(lon1,lat1,lon2,lat2, deltat)
    print 'nparticles=%d K_rr=%e K_xx=%e K_xy=%e K_yy=%e'%(nparticles, K_rr, K_xx, K_xy, K_yy)

    #}}}

def two_dim_random_walk_iso(N=10000): #{{{
    print '================================================================================'
    print ' isotropic two dimensional random walk'
    print '================================================================================'

    lat0 = 35.;
    r = np.random.rand(N)*3.0
    theta = np.random.rand(N)*2*np.pi
    jx  = np.random.rand(N)
    jy  = np.random.rand(N)
    lon1 = np.radians(np.cos(theta)*r + jx)
    lat1 = np.radians(lat0 + np.sin(theta)*r + jy)

    deltat = 2*24*60*60
    K = 1e5
    print 'Analytical K_xx = %e'%(K)
    print 'Analytical K_yy = %e'%(K)
    r = 1/3. # e.g., Visser1997, Gross et al POD report, note K must be converted to degree^2/s units
    dx = (2*np.random.rand(N)-1)*np.sqrt(2/r*deltat*K/(rEarth*np.cos(np.abs(lat1)))**2.0)
    dy = (2*np.random.rand(N)-1)*np.sqrt(2/r*deltat*K/rEarth**2.0)

    lon2 = lon1 + dx
    lat2 = lat1 + dy

    nparticles, _, _, K_rr, K_xx, K_xy, K_yy = compute_diffusivity(lon1,lat1,lon2,lat2, deltat)
    print 'nparticles=%d K_rr=%e K_xx=%e K_xy=%e K_yy=%e'%(nparticles, K_rr, K_xx, K_xy, K_yy)

    #}}}

def two_dim_random_walk_aniso(N=10000): #{{{
    print '================================================================================'
    print ' aniosotropic two dimensional random walk'
    print '================================================================================'

    lat0 = 35.;
    r = np.random.rand(N)*3.0
    theta = np.random.rand(N)*2*np.pi
    jx  = np.random.rand(N)
    jy  = np.random.rand(N)
    lon1 = np.radians(np.cos(theta)*r + jx)
    lat1 = np.radians(lat0 + np.sin(theta)*r + jy)

    deltat = 2*24*60*60
    Kxx = 1e5
    Kyy = 1e4
    print 'Analytical K_xx = %e'%(Kxx)
    print 'Analytical K_yy = %e'%(Kyy)
    r = 1/3. # e.g., Visser1997, Gross et al POD report, note K must be converted to degree^2/s units
    dx = (2*np.random.rand(N)-1)*np.sqrt(2/r*deltat*Kxx/(rEarth*np.cos(np.abs(lat1)))**2.0)
    dy = (2*np.random.rand(N)-1)*np.sqrt(2/r*deltat*Kyy/rEarth**2.0)

    lon2 = lon1 + dx
    lat2 = lat1 + dy

    nparticles, _, _, K_rr, K_xx, K_xy, K_yy = compute_diffusivity(lon1,lat1,lon2,lat2, deltat)
    print 'nparticles=%d K_rr=%e K_xx=%e K_xy=%e K_yy=%e'%(nparticles, K_rr, K_xx, K_xy, K_yy)

    #}}}

def two_dim_random_walk_iso_time_evolution(N=200): #{{{
    print '================================================================================'
    print ' isotropic two dimensional random walk with time evolution'
    print '================================================================================'
    import seaborn as sns

    lat0 = 35.
    ddeg = 0.05
    rmax  = 100000
    lon1, lat1 = np.meshgrid(np.linspace(-ddeg, ddeg,N), np.linspace(lat0-ddeg, lat0+ddeg, N))
    dist = haversine_formula(lat0, lat1, 0.0, lon1, rEarth)
    particles  = dist < rmax
    lon1 = lon1[particles]
    lat1 = lat1[particles]
    N = lon1.shape[0]

    deltat = 10*24*60*60
    K = 1e3
    print 'Analytical K_xx = %e'%(K)
    print 'Analytical K_yy = %e'%(K)
    r = 1/3. # e.g., Visser1997, Gross et al POD report, note K must be converted to degree^2/s units
    lon2=lon1.copy()
    lat2=lat1.copy()
    nsteps = 100
    for i in np.arange(nsteps):
        dx = (2*np.random.rand(N)-1)*np.sqrt(2/r*deltat/nsteps*K/(rEarth*np.cos(np.abs(lat1)))**2.0)
        dy = (2*np.random.rand(N)-1)*np.sqrt(2/r*deltat/nsteps*K/rEarth**2.0)

        lon2 += dx
        lat2 += dy

    def scatter_hist(x,y):
        def count(x,y):
            return x.shape[0]
        x = np.degrees(x)
        y=np.degrees(y)
        jp = sns.jointplot(x,y, marginal_kws={"fit": stats.norm})
        jp.annotate(count, template="{stat} = {val:d}", frameon=False)
        sns.despine()


    nparticles, _, _, drdr, dxdx, dxdy, dydy = compute_com_dispersion(lon1, lat1)
    print '1st cluster: nparticles=%d D_rr=%e D_xx=%e K_xy=%e K_yy=%e'%(nparticles, drdr, dxdx, dxdy, dydy)
    scatter_hist(lon1,lat1)

    nparticles, _, _, drdr, dxdx, dxdy, dydy = compute_com_dispersion(lon2, lat2)
    print '1st cluster: nparticles=%d D_rr=%e D_xx=%e K_xy=%e K_yy=%e'%(nparticles, drdr, dxdx, dxdy, dydy)
    scatter_hist(lon2,lat2)

    nparticles, _, _, K_rr, K_xx, K_xy, K_yy = compute_diffusivity(lon1,lat1,lon2,lat2, deltat)
    print 'nparticles=%d K_rr=%e K_xx=%e K_xy=%e K_yy=%e'%(nparticles, K_rr, K_xx, K_xy, K_yy)

    return #}}}

if __name__=="__main__":
    print '********************************************************************************'
    print '* Library to compute diffusivity'
    print '********************************************************************************'
    run_tests()
