#!/usr/bin/env python
"""
Module to compute effective diffusivity

Phillip J Wolfram
01/21/2016
"""

import os
import numpy as np
import scipy.sparse as sp
import numexpr as ne
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime
import xarray as xr
from scipy.interpolate import interp1d
from GeneralTriangulation import GeneralTriangulation, edges_containing_points, \
        triangles_continaing_points, unstructured_interp, unstructured_remap
from latlon_coordinate_transforms import fix_periodic_timeseries, fix_periodicity, \
        haversine_formula_numexpr as haversine_formula
from cluster_topology import cluster_topology
from iotasks import agg_result_new_dim as agg_results, agg_cat_general, agg_result_new_dim
from comp_funcs import SignedLog10, ceilround

def isxr(atype):
    return isinstance(atype, xr.Dataset) or isinstance(atype, xr.DataArray)

def radially_constant_distribution_topology(tri, x0, y0, radius): #{{{
    """ returns clusters based on mesh topology for a given radius at (xloc,yloc) """


    tri.compute_com()
    tri.compute_le()
    x = tri.x
    y = tri.y

    if not tri.meshname:

        if tri.latlon:
            points = np.where(haversine_formula(y, y0, x, x0) < radius)[0]
        else:
            distsq = np.sqrt((x-x0)**2 + (y-y0)**2)
            points = np.where(distsq < radius)[0]

        clustertri = tri.triangles[triangles_continaing_points(points, tri.triangles),:]
        trim = GeneralTriangulation(x, y, tri=clustertri, latlon=tri.latlon)
        # remove hanging triangles
        clustertri = trim.triangles[np.where(np.sum(trim.neighbors != -1, axis=-1) > 1)[0],:]
        trim = GeneralTriangulation(x, y, tri=clustertri, latlon=tri.latlon)
        trim.compute_edgecellneighs()

        bndryedges = np.where(np.sum(trim.edgecellneighs == -1,axis=1))
        ringpoints = np.unique(trim.edges[bndryedges,:].ravel())

        meshpoints = np.unique(trim.triangles.ravel())

    else:

        assert tri.meshname is not None, 'Assumes that MPAS mesh was used to initialize particle grid, no mesh given'
        # find acell index closest to xloc, yloc (pick first of nearest)
        distsq = (x-x0)**2 + (y-y0)**2
        apoint = np.where(distsq == np.min(distsq))[0][0]

        points, meshpoints, ringpoints = cluster_topology(xr.open_dataset(tri.meshname), apoint, x, y, radius)

        meshpoints = points[meshpoints]

        ringpoints = points[ringpoints]

    # all edges that have vertices in ringpoints
    ringedges = edges_containing_points(ringpoints, tri.edges)

    # all edges that have vertices in meshpoints
    meshedges = edges_containing_points(meshpoints, tri.edges)

    centeredges = meshedges*np.logical_not(ringedges)

    # all cells that have vertices in points
    cellid = triangles_continaing_points(points, tri.triangles)

    return ringedges, centeredges, cellid #}}}

def meridionally_constant_distribution_topology(tri, yloc): #{{{
    """ computes the topology for a meridionally constant concentration """
    # adjust yloc
    minorshift = 0.01
    scy = 0.5*(np.min(tri.y[tri.triangles],axis=1) + np.max(tri.y[tri.triangles],axis=1))
    scyunq = np.unique(scy)
    if yloc < np.min(scy):
        yloc = scyunq[0] + minorshift*(scyunq[1] - scyunq[0])
    elif yloc > np.max(scy):
        yloc = scyunq[-1] - minorshift*(scyunq[-1] - scyunq[-2])

    assert not np.max(yloc == scy), 'Location and points are coincident, cannot form masks properly!'

    tri.compute_com()

    # get line of center cells
    dl = np.mean(np.diff(np.unique(scy)))
    cellid = np.logical_and(yloc - dl/2.0 < scy, yloc + dl/2.0 > scy)

    # get edge masks
    tri.compute_celledges()
    centeredges = tri.celledges[cellid,:].ravel()
    edges = np.zeros_like(tri.edges[:,0], dtype='bool')
    edges[centeredges] = True

    # center mask has non-zero cross product with (0,1)
    crossprod = tri.return_edge_normal()[0,centeredges]
    tol = np.mean(np.abs(crossprod))
    topbottommask = np.where(np.abs(crossprod) < tol)[0]
    boundaryedges = centeredges[topbottommask]

    centeredges = edges.copy()
    centeredges[boundaryedges] = False
    centeredges = np.where(centeredges)

    # handle case where the edges are on a boundary
    onbndry = np.min(tri.edgecellneighs[boundaryedges,:],axis=1) == -1
    boundaryedges = boundaryedges[np.logical_not(onbndry)]

    return boundaryedges, centeredges, cellid #}}}

def init_gaussian_conc(x, y, x0, y0, Lwidth, latlon=False): #{{{
    if latlon:
        r = haversine_formula(y, y0, x, x0)
    else:
        r = np.sqrt((x-x0)**2 + (y-y0)**2)
    return np.exp(-(r/Lwidth)**2) #}}}

def init_gaussian_line_conc(y, y0, Lwidth): #{{{
    return np.exp(-((y-y0)/Lwidth)**2) #}}}

def init_delta_line_conc(y, ids): #{{{
    delta = np.zeros_like(y)
    delta[ids] = 1
    return delta #}}}

def init_gradient_line_conc(y, concgrad): #{{{
    yunq = np.unique(y)
    maxy = np.max(yunq)
    miny = np.min(yunq)
    ymid = 0.5*(miny + maxy)
    Ly = maxy - miny
    conc = (y-miny)*concgrad
    minc = np.min(conc)
    maxc = np.max(conc)
    assert minc >= 0 and maxc <= 1.0, 'Concentration gradient is too strong to ensure that C in (0,1)'
    return conc #}}}

def simplify_plot(): #{{{
    plt.axis('off')
    plt.title('')
    # hard coded examples
    plt.ylim(0,2e6)
    plt.xlim(-2e6,2e6)
    plt.gca().set_aspect('equal', adjustable='box')
    return #}}}

class CummulativeIntegral: #{{{
    # class constructors #{{{
    def __init__(self, nbins=50, ntimes=1):#{{{
        """
        Parameters:

            nbins       number of bins for the histogram
            ntimes      number of time slices for scalar that this routine record

            # used later to store data
            yhist       histogrammed Y over time
            xhist       histogrammed X over time
            yint        scalar cummulative integral over time
            dycumsumdx  derivative in ycumsum by x, e.g.,

        """

        self.ntimes = ntimes

        # compute the bins, assumes that x coordinate normalized to (0,1)
        self.nbins = nbins
        self.dbin = 1.0/nbins
        self.bins = np.linspace(-self.dbin/2.0,1.0+self.dbin/2.0,self.nbins+1)

        # define the integrals
        self.y = np.nan*np.zeros((self.ntimes,self.nbins))
        self.yint = np.nan*np.zeros((self.ntimes,self.nbins))

        return #}}}
    #}}}

    # static methods #{{{
    @staticmethod
    def binned_mean(xin, yin, bins): #{{{
        """ Usage: CummulativeIntegral.binned_mean(xin, yin, bins) """

        sumdata = np.histogram(yin, bins=bins, weights=xin)[0]
        countdata = np.histogram(yin, bins=bins, weights=np.ones_like(xin))[0]

        mean = sumdata/countdata

        return mean #}}}
    #}}}

    # class methods #{{{
    def compute_integral_sum(self, xin, yin, timeindex): #{{{
        """
        Performs integral of
            int_{xin < bins} yin

            where yin = y*dA is assumed

            This is analagous to

            int_{c <= C(A,t)} (xi*dA)

            where the A discretization is set by c=xin, bins=C(A,t), and yi = (xi*dA)

            This is just a cummulative distribution of the histogram of yin over xin.

            timeindex is the time index for storage.
        """

        self.y[timeindex, :] = np.histogram(yin, bins=self.bins, weights=xin)[0]

        self.yint[timeindex, :] = np.cumsum(np.histogram(yin, bins=self.bins, weights=xin)[0])

        return #}}}

    def compute_xderivative(self, smoothing=None): #{{{
        """
        Computes
        d (\int_{xin < bins} yin) / dbins

        """
        import cv2

        if smoothing is not None:
            field = cv2.GaussianBlur(self.yint, smoothing, 0)
        else:
            field = self.yint

        xc = 0.5*(self.bins[:-1] + self.bins[1:])
        dyintdx = np.diff(field, axis=1)/np.diff(xc)

        return dyintdx #}}}

    def return_fractional_value(self, ignorelastbin, onderiv=True): #{{{
        """
        Returns the fractional value over the time axis (for each time slice)
        for y attribute.
        """

        yfrac = self.y.copy()

        if ignorelastbin:
            # don't consider all the the last bin (zeros if 0 <= C <= 1) as setup
            yfrac[:,0] = 0

        yfrac = yfrac / np.sum(yfrac,axis=1)[:,np.newaxis]

        if onderiv:
            # match to be same size as dyintdx from compute_xderivative
            # so that plotting can occur readily on same grid as for kappaeff
            # and Keff
            yfrac = 0.5*(yfrac[:,:-1] + yfrac[:,1:])

        return yfrac #}}}

    #}}}

# tests #{{{
def test_cummulative_integral(ntimes=200, ncells=1000): #{{{
    """
    Test considers an example problem for computing the effective diffusivity:

    c(x,t) = (1-x)*(1-t) + t/2 distribution in x in (0,1) and t in (0,1)

    |dc/dx|^2  = (t-1)^2

    for A = int_{c<C} dA

    where

    A(C,t) =       0                                   C <= t/2

                   (t+2(C-1))     t + 2(t/2 -1)
                 - ----------  +  -------------        t/2 < C < (1-t/2)
                    2(t-1)           2(t-1)

                   1                                   C >= (1-t/2)


    and

    dA/dc =         0               C <= t/2

                        1
                    - -----        t/2 < C < (1-t/2)
                      (t-1)

                    0               C >= (1-t/2)


    kappaeff = k * d/dc(int_{c<C} |dc/dx|^2 dA) * dA/dc

    where int_{c<C} |dc/dx|^2 dA = |dc/dx|^2 int_{c<C} dA for this case

    int_{c<C} |dc/dx|^2 dA  =  (t-1)^2 * A(c,t)

    Then,

    kappaeff = k * (t-1)^2 * (dA/dc)^2

         =              0               C <= t/2

                        1             t/2 < C < (1-t/2)

                        0               C >= (1-t/2)


    """
    kwargs = {'aspect':'auto', 'interpolation':'none', 'cmap':'viridis'}

    x = np.linspace(0,1,ncells)[np.newaxis,:]
    time = np.linspace(0,1,ntimes)[:,np.newaxis]

    area = np.ones((ntimes,ncells))/ncells
    dcdx2 = (time-1)**2.0
    gradc = dcdx2

    conc = (1-x)*(1-time) + time/2

    # plot concentration #{{{
    plt.imshow(conc, **kwargs)
    plt.title('c(x,t)')
    plt.xlabel('x')
    plt.ylabel('t')
    plt.colorbar()
    #}}}

    areaA = CummulativeIntegral(ntimes=ntimes)
    gradC = CummulativeIntegral(ntimes=ntimes)

    for tid, atime in enumerate(time):
        areaA.compute_integral_sum(area[tid,:], conc[tid,:], tid)
        gradC.compute_integral_sum(gradc[tid,:]*area[tid,:], conc[tid,:], tid)

    dareaAdc = areaA.compute_xderivative(smoothing=(5,3))
    dgradCdc = gradC.compute_xderivative(smoothing=(5,3))

    # get kappaeff
    kappaeff = dareaAdc*dgradCdc

    # print computed A #{{{
    plt.figure()
    plt.subplot(2,1,1)
    plt.imshow(areaA.yint, **kwargs)
    plt.title('A(c,t) numerical')
    plt.ylabel('t')
    plt.clim([0, 1.])
    plt.colorbar()

    # compare with exact A
    plt.subplot(2,1,2)
    ce = (0.5*(areaA.bins[:-1] + areaA.bins[1:])[:,np.newaxis]).T
    exact = np.nan*np.ones_like(ce*time)
    exact[np.where(ce <= (1-time/2))] = -((time + 2*(ce-1))/(2*(time-1)) - \
                                        (time + 2*(time/2-1))/(2*(time-1)))[np.where(ce < (1-time/2))]
    exact[np.where(ce < time/2)] = (0*time*ce)[np.where(ce < time/2)]
    exact[np.where(ce > 1-time/2)] = (np.ones_like(ce*time))[np.where(ce > 1-time/2)]
    plt.imshow(exact, **kwargs)
    plt.title('A(c,t) exact')
    plt.xlabel('c')
    plt.ylabel('t')
    plt.clim([0, 1.])
    plt.colorbar()
    #}}}

    # print computed dA/dc #{{{
    plt.figure()
    plt.subplot(2,1,1)
    plt.imshow(dareaAdc, **kwargs)
    plt.title('dA(c,t)/dc numerical')
    plt.ylabel('t')
    plt.clim([0, 50])
    plt.colorbar()

    # compare with exact dA/dc
    plt.subplot(2,1,2)
    dexactdc = 0*exact
    dexactdc[np.where(ce <= (1-time/2))] = -(np.ones_like(ce)/(time-1))[np.where(ce < (1-time/2))]
    dexactdc[np.where(ce < time/2)] = (0*time*ce)[np.where(ce < time/2)]
    dexactdc[np.where(ce > 1-time/2)] = (0*ce*time)[np.where(ce > 1-time/2)]

    plt.imshow(dexactdc, **kwargs)
    plt.title('dA(c,t)/dc exact')
    plt.xlabel('c')
    plt.ylabel('t')
    plt.clim([0, 50])
    plt.colorbar()
    #}}}

    # print computed kappaeff #{{{
    plt.figure()
    plt.subplot(2,1,1)
    plt.imshow(kappaeff, **kwargs)
    plt.title('kappaeff/kappa numerical')
    plt.ylabel('t')
    plt.clim([0,1.1])
    plt.colorbar()

    # compare with exact dA/dc
    plt.subplot(2,1,2)
    kappaeffexact = 0*exact
    kappaeffexact[np.where(ce <= (1-time/2))] = (np.ones_like(exact))[np.where(ce < (1-time/2))]
    kappaeffexact[np.where(ce < time/2)] = (0*time*ce)[np.where(ce < time/2)]
    kappaeffexact[np.where(ce > 1-time/2)] = (0*ce*time)[np.where(ce > 1-time/2)]

    plt.imshow(kappaeffexact, **kwargs)
    plt.title('kappaeff/kappa exact')
    plt.xlabel('c')
    plt.ylabel('t')
    plt.clim([0,1.1])
    plt.colorbar()



    #}}}

    return #}}}

#}}}

#}}}

class EffectiveDiffusivity: #{{{
    # class constructors #{{{
    def __init__(self, xt=None, yt=None, mask=None, tparticles=None, dt=60.0, dtout=1.*86400,
            calctype='ygauss', concwidth=50.e3, concgrad=1./2e6, kappa=10.0, x0=None, y0=None,
            meshname=None, cbins=50, run=False, latlon=False, cleanupgrid=False, dtreset=None,
            advection='lagrangian', advtype=None, diffusion=True, subdivide=None,
            simpleplot=True, saveloc='./'):  #{{{
        """
        Parameters:

            xt          particle timeseries in x in m

            yt          particle timeseries in y in m

            tparticles  time vector corresponding to particle positions xt,yt in days

            dt          time step for effective diffusivity calculation, e.g., 60 s

            dtout       output time, e.g., 1 day

            dtreset     reset time, (e.g., particles reset for 5 days + 1s in namelist
                                     is dtreset=5.*86400, assumes first particle entry in
                                     time is reset position)

            concwidth   standard deviation width used to compute initial Gaussian
                        concentration profile, e.g., 50km

            meshname    name of mesh for topology of particle control volumes

            kappa       background diffusivity, e.g., 10 m^2s^-1
                        Klocker, Andreas, et al. "Reconciling float-based and
                        tracer-based estimates of lateral diffusivities."
                        Journal of Marine Research 70.4 (2012): 569-602.

            cbins       number of bins for effective diffusivity calculation

            advection   string to specify if particle advection is to be used 'lagrangian', or
                        otherwise points to an eulerian dataset

            diffusion   boolean to specify if particle diffusion is to be used

            subdivide   number of times to subdivide input mesh / data

            calctype    'y'- compute effective diffusivity with initial distribution that
                        is constant in y, effective computes meridional component of
                        effective diffusivity

            x0,y0       location of starting point for effective diffusivity calculation

        """

        self.simpleplot = simpleplot
        self.saveloc = saveloc
        if not os.path.exists(saveloc):
            os.makedirs(saveloc)

        self.advection = advection
        self.advtype = advtype
        self.diffusion = diffusion

        self.kappa = float(kappa)

        self.cbins = cbins   # concentration bins

        # compute time properties and functions
        if tparticles is not None:
            self.tparticles = tparticles
            self.tend = tparticles[-1]
        else:
            print 'Warning: ending time was not specified via tparticles!'
            # assume there will only be an output of a single timestep
            self.tend = dtout

        assert dt < self.tend, 'Time step should be smaller than ending time'
        assert dt < dtout, 'Output time step must be larger than computational timestep'

        self.t = 0
        self.timeindex = 0
        self.dt = dt
        self.ntoutinterval = int(np.ceil(dtout/dt))
        self.dtout = self.ntoutinterval*dt
        self.dtreset = dtreset
        self.ntout = 1 + int(np.floor(self.tend/self.dtout))
        self.tvec = np.hstack((np.array([0.0]),
            np.cumsum(np.ones((int(np.ceil(self.tend/self.dt))))*self.dt)))

        if meshname is not None and (xt is not None or yt is not None):
            ds = xr.open_dataset(meshname)
            if ds.x_period:
                xt = fix_periodic_timeseries(xt, ds.x_period)
                self.x_period = ds.x_period
            else:
                self.x_period = None
            if ds.y_period:
                assert False, 'Implementation for y-periodicity is uncertain'
            ds.close()
        else:
            self.x_period = None
            self.y_period = None

        # perform subdivision to get enhanced resolution
        if subdivide is not None:
            pass


        # set up particle time integrators
        if not isxr(advection) and (advection == True or advection == 'lagrangian'):
            if advection and (xt is not None and yt is not None):
                # if resets are employed
                if self.dtreset is not None:

                    assert (float(self.dtreset)/float(self.dt)).is_integer(), \
                    'dtreset and dt need to be multiples of each other to perform a clean reset'

                    ttemp = np.sort(np.hstack(
                        (np.arange(self.dtreset+1,tparticles[-1], self.dtreset),
                        tparticles))
                        )
                    self.xp0= np.copy(xt[0,:])
                    self.yp0= np.copy(yt[0,:])
                    xtemp = self.xp0[np.newaxis,:]*np.ones((ttemp.shape[0],xt.shape[1]))
                    ytemp = self.yp0[np.newaxis,:]*np.ones((ttemp.shape[0],yt.shape[1]))
                    idx = np.where(np.in1d(ttemp, tparticles))
                    xtemp[idx,:] = xt[:,:]
                    ytemp[idx,:] = yt[:,:]
                    # need to fix periodicity here too...
                    sid = 0
                    newids = np.asarray(np.where(np.mod(ttemp,2))[0])
                    for eid in newids:
                        if self.x_period:
                            xtemp[sid:eid,:] = fix_periodic_timeseries(xtemp[sid:eid,:], self.x_period)
                        if ds.y_period:
                            assert False, 'Implementation for y-periodicity is uncertain'
                        sid = eid
                # no resets used
                else:
                    ttemp = tparticles
                    xtemp = xt
                    ytemp = yt
                self.xint = interp1d(ttemp, xtemp.T)
                self.yint = interp1d(ttemp, ytemp.T)
        else:
            if isinstance(advection, str) and not isxr(advection):
                vel = xr.open_dataset(advection)
            else:
                vel = advection
            if len(vel.nBuoyancySurfaces) == 3:
                u = vel.buoyancySurfaceVelocityZonal.isel(nBuoyancySurfaces=1).squeeze()
                v = vel.buoyancySurfaceVelocityMeridional.isel(nBuoyancySurfaces=1).squeeze()
            else:
                u = vel.buoyancySurfaceVelocityZonal.squeeze()
                v = vel.buoyancySurfaceVelocityMeridional.squeeze()

            self.uint = interp1d(tparticles, u.T)
            self.vint = interp1d(tparticles, v.T)
            if len(vel.nBuoyancySurfaces) == 3:
                buoyancySurfaceThicknessEst = 0.5*(vel.buoyancySurfaceDepth.isel(nBuoyancySurfaces=0) - 
                                                   vel.buoyancySurfaceDepth.isel(nBuoyancySurfaces=2))
                self.hint = interp1d(tparticles, buoyancySurfaceThicknessEst.squeeze().T)
            else:
                self.hint = interp1d(tparticles, np.ones_like(u.T))

        # set up mesh
        self.latlon = latlon
        if meshname:
            self.tri = GeneralTriangulation.from_MPAS_mesh(meshname, self.latlon, subdivide)
        else:
            self.tri = GeneralTriangulation(xt[0,:], yt[0,:], latlon=self.latlon)
            if cleanupgrid:
                badcells = np.where(np.sum(self.tri.neighbors == -1,axis=1))[0]
                self.tri.compute_com()
                # include their neighbors too (4 times)
                for i in np.arange(cleanupgrid):
                    badcells = np.unique(np.concatenate((badcells, np.unique(self.tri.neighbors[badcells,:]))))
                    # remove -1 not a number cells
                    if badcells[0] == -1:
                        badcells = badcells[1:]
                goodcells = np.ones(self.tri.triangles.shape[0], dtype='bool')
                goodcells[badcells] = False
                tricleaned = self.tri.triangles[goodcells,:]
                self.tri = GeneralTriangulation(x=xt[0,:], y=yt[0,:], tri=tricleaned, latlon=self.latlon)

        # store intial positions
        self.xinit = self.tri.x.copy()
        self.yinit = self.tri.y.copy()

        def expand_arrays(array, nsize): #{{{
            # if the value is a scalar or
            if type(array) == np.dtype('int') or type(array) == np.dtype('float'):
                newarray = array*np.ones(nsize)
            else:
                newarray = array
            return newarray #}}}

        if type(y0) in [np.float, np.int, np.float64, np.int64]:
            y0 = np.array([y0])
        #if type(x0) in [np.float, np.int, np.float64, np.int64]:
        #    x0 = np.array([x0])
        self.ntracers = len(y0)
        concwidth = expand_arrays(concwidth , self.ntracers)
        concgrad  = expand_arrays(concgrad,   self.ntracers)
        x0        = expand_arrays(x0        , self.ntracers)
        y0        = expand_arrays(y0        , self.ntracers)

        if mask is not None:
            trimask = self.tri.avg_from_tri_to_vertices(mask)
        else:
            trimask = 1

        # set up initial concentration corresponding to specified x0,y0
        for atracer in np.arange(self.ntracers):
            if calctype == 'ygauss':
                conc, qinit, centeredges, boundaryedges = \
                        self.setup_concentration(yloc=y0[atracer], concwidth=concwidth[atracer], disttype='ygauss')
            elif calctype == 'gauss':
                conc, qinit, centeredges, boundaryedges = \
                        self.setup_concentration(xloc=x0[atracer], yloc=y0[atracer], concwidth=concwidth[atracer], disttype='gauss')
            elif calctype == 'ydelta':
                conc, qinit, centeredges, boundaryedges = \
                        self.setup_concentration(yloc=y0[atracer], disttype='ydelta')
            elif calctype == 'ygrad':
                conc, qinit, centeredges, boundaryedges = \
                        self.setup_concentration(yloc=y0[atracer], concgrad=concgrad[atracer], disttype='ygrad')
            elif calctype == 'const':
                conc, qinit, centeredges, boundaryedges = self.setup_concentration(disttype='const')
            else:
                assert False, 'No other calctype of effective diffusivity has yet '\
                'been implemented for %s.'%(calctype)

            # initialize tracers
            if atracer == 0:
                self.tri.centeredges = []
                self.tri.boundaryedges = []
                self.conc  = np.nan*np.zeros((self.ntracers, conc.shape[0]))
                self.qinit = np.nan*np.zeros((self.ntracers))
                self.mask  = np.nan*np.zeros((self.ntracers))

            self.conc[atracer,:] = conc
            self.qinit[atracer]  = qinit
            self.mask[atracer]   = np.sum(conc*trimask)/np.sum(conc)

            self.tri.centeredges.append(centeredges)
            self.tri.boundaryedges.append(boundaryedges)

        self.calctype = calctype
        self.yloc = y0
        self.xloc = x0
        self.concwidth = concwidth

        if run:
            self.compute_time_advancement(plotconc=False)

        return #}}}
    #}}}

    # class methods #{{{
    def plot_concentration_triangles(self, trival, inkm=True, indays=True, alpha0=0.3, clims=[0,1.0], titleprefix='', simpleplot=False): #{{{

        if inkm:
            sflen = 1e3
        else:
            sflen = 1.0
        if indays:
            sftime = 24*60*60.
        else:
            sftime = 1.0

        triT = self.tri

        plt.figure()
        triT.plot_scalar(scalar=trival, cmap='viridis', alpha=alpha0)

        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/sflen))
        plt.gca().xaxis.set_major_formatter(ticks_x)

        ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/sflen))
        plt.gca().yaxis.set_major_formatter(ticks_y)

        plt.title(titleprefix + 't=%.1f days'%(self.t/sftime))
        if clims is not None:
            plt.clim(clims)

        if simpleplot:
            simplify_plot()
        else:
            plt.axis('equal')
            cb = plt.colorbar()
            cb.set_alpha(1)
            cb.draw_all()

        return #}}}

    def plot_concentration_vertices(self, trival=None, vertval=None, inkm=True, indays=True, #{{{
            alpha0=1.0, clims=[0.0, 1.0], cmap='viridis', titleprefix='', simpleplot=False, save=None):
        """ Plots the concentration on vertices via values on triangles, `trival`, or on vertices, `vertval`"""

        if inkm:
            sflen = 1e3
        else:
            sflen = 1.0
        if indays:
            sftime = 24*60*60.
        else:
            sftime = 1.0

        triT = self.tri

        plt.figure()

        if trival is not None:
            cP = triT.avg_to_vertices(trival)
        else:
            if vertval is not None:
                cP = vertval
            else:
                assert False, 'Must supply plot_concentration_vertices with value to plot'

        idx = np.argsort(cP)
        xp = triT.x[idx].copy()
        yp = triT.y[idx].copy()
        if self.x_period is not None:
            xp = np.mod(xp, self.x_period)

        plt.scatter(xp, yp, c=cP[idx], marker='.', \
                cmap=cmap, edgecolor='None', alpha=alpha0)

        if save is not None:
            np.savez(self.saveloc + save + 'conc_vertices_t%04d_days.npz'%(self.t/86400.), x=xp, y=yp, c=cP[idx])

        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/sflen))
        plt.gca().xaxis.set_major_formatter(ticks_x)

        ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/sflen))
        plt.gca().yaxis.set_major_formatter(ticks_y)

        plt.title(titleprefix + 't=%.1f days'%(self.t/sftime))
        if clims is not None:
            plt.clim(clims)
        else:
            plt.colorbar()

        if simpleplot:
            simplify_plot()
        else:
            plt.axis('equal')
            cb = plt.colorbar()
            cb.set_alpha(1)
            cb.draw_all()

        return #}}}

    def plot_effective_diffusivity(self, indays=True, alpha0=0.3, ignorelastbin=False, inverty=False): #{{{
        """ Plots effective diffusivity Hovmoller diagram """

        if indays:
            sftime = 24*60*60.
        else:
            sftime = 1.0

        for atracer in np.arange(self.ntracers):
            titleprefix = 'tracer %d: '%(atracer)
            plt.figure()
            #plt.hold(True)
            plt.contourf(self.caxis,self.outtime/sftime,SignedLog10.convert(self.Keff[atracer,:,:]),\
                    levels=np.linspace(2,5,20), cmap='magma', extend='both')
            plt.colorbar(ticks=[2,3,4,5])
            cs = plt.contour(self.caxis, self.outtime/sftime, self.area[atracer].return_fractional_value(ignorelastbin), \
                    colors='w', levels=np.linspace(0,1.0,10 + 1), lw=2, alpha=alpha0)
            plt.clabel(cs, inline=1, fontsize=6, alpha=alpha0, fmt='%.1f')
            plt.xlabel('C')
            plt.ylabel('t (days)')
            if inverty:
                plt.gca().invert_yaxis()
            plt.xlim([0,1.0])
            plt.title(r"$\log_{10}(K_\mathrm{eff})$ (m$^2$ s$^{-1}$)")
            plt.savefig(self.saveloc + '/keff_conc_tracer%0d.png'%(atracer))

            areacoord = 0.5*(self.area[atracer].yint[:,:-1] + self.area[atracer].yint[:,1:])/1e3/1e3/1e6
            timecoord = (self.outtime/sftime)[:,np.newaxis]*np.ones_like(areacoord)
            plt.figure()
            #plt.hold(True)
            plt.contourf(areacoord, timecoord, SignedLog10.convert(self.Keff[atracer,:,:]),\
                    levels=np.linspace(2,5,20), cmap='magma', extend='both')
            plt.colorbar(ticks=[2,3,4,5])
            cs = plt.contour(areacoord, timecoord, self.area[atracer].return_fractional_value(ignorelastbin), \
                    colors='w', levels=np.linspace(0,1.0,10 + 1), lw=2, alpha=alpha0)
            plt.clabel(cs, inline=1, fontsize=6, alpha=alpha0, fmt='%.1f')
            #plt.xlim([1.75e6,2e6])
            plt.ylabel('t (days)')
            if inverty:
                plt.gca().invert_yaxis()
            plt.xlabel('C')
            #plt.xlim([0,1.0])
            #plt.title(titleprefix + r"$\log_{10}(K_\mathrm{eff})$ (m$^2$ s$^{-1}$)")
            plt.xlabel('$A$ ($10^6$ km$^2$)')
            plt.clim([0,6])
            plt.title(r"$\log_{10}(K_\mathrm{eff})$ (m$^2$ s$^{-1}$)")
            plt.savefig(self.saveloc + '/keff_area_tracer%0d.png'%(atracer))

            np.savez(self.saveloc + '/keffdata_tracer%04d.npz'%(atracer),
                    caxis=self.caxis, t=self.outtime/sftime, keff=self.Keff[atracer,:,:],
                    area=self.area[atracer].return_fractional_value(ignorelastbin),
                    areacoord=areacoord, timecoord=timecoord, ybar=self.ybar[atracer, :, :])

        return #}}}

    def plot_conc_area(self, indays=True, ignorelastbin=True): #{{{
        """ Plots concentration area Hovmoller diagram, Ignores bin
        with zero by default to not artificially bias the results """

        if indays:
            sftime = 24*60*60.
        else:
            sftime = 1.0


        for atracer in np.arange(self.ntracers):
            titleprefix = 'tracer %d'%(atracer)

            xc = 0.5*(self.area[atracer].bins[:-1] + self.area[atracer].bins[1:])

            plt.figure()
            plt.contourf(xc, self.outtime/sftime, self.area[atracer].yint, cmap='viridis')
            plt.colorbar()
            plt.xlabel('C')
            plt.ylabel('t (days)')
            plt.xlim([0,1.0])
            plt.title(titleprefix + r"${A}$ (m$^2$)")

            plt.figure()
            fracarea = self.area[atracer].return_fractional_value(ignorelastbin, onderiv=False)
            maxfrac = ceilround(fracarea.max(), 1)
            plt.contourf(xc, self.outtime/sftime, fracarea, levels=np.linspace(0,maxfrac,2*int(maxfrac/0.1) + 1), cmap='viridis')
            #plt.clim([0, 1.])
            plt.colorbar()
            plt.xlabel('C')
            plt.ylabel('t (days)')
            plt.xlim([0,1.0])
            plt.title(titleprefix + r"$f_{A}$ (-)")

        return #}}}

    def plot_effective_lengths(self, indays=True, inkm=True): #{{{

        if indays:
            sftime = 24*60*60.
            timelabel = 't (days)'
        else:
            sftime = 1.0
            timelabel = 't (s)'
        if inkm:
            sflen = 1.e3
        else:
            sflen = 1.


        for atracer in np.arange(self.ntracers):
            titleprefix = 'tracer %d'%(atracer)

            plt.figure()
            #plt.hold(True)

            def local_plot(values, linetype, label, *args, **kwargs):
                plt.semilogy(self.outtime/sftime, values, linetype, label=label, *args, **kwargs)

            local_plot(self.Lmin[atracer]*np.ones_like(self.Leffinterp[atracer,:])/sflen, 'k-', 'Lmin', lw=3)
            local_plot(self.LeffStarconstgrad[atracer,:]/sflen,                  'r-', 'Leff* const grad', lw=3)
            local_plot(self.LeffStarvargrad[atracer,:]/sflen,                    'k-', 'Leff* var grad', lw=1)
            local_plot(self.Leffavg[atracer,:]/sflen,                            'g-', 'Leff average', lw=1)
            local_plot(self.Leffmax[atracer,:]/sflen,                            'b-', 'Leff max', lw=1)

            plt.legend(loc='best', fontsize=9)
            plt.xlabel(timelabel)
            if inkm:
                plt.ylabel('Effective lengths (km)')
            else:
                plt.ylabel('Effective lengths (m)')
            plt.ylim([np.round(np.log10(self.Lmin[atracer]/10)), plt.ylim()[-1]])
            plt.title(titleprefix)
            plt.savefig(self.saveloc + '/effective_lengths_tracer%0d.png'%(atracer))

        return #}}}

    def plot_scalar_variance(self): #{{{
        """ Plots global scalar variance time series """

        sftime = 24*60*60.

        for atracer in np.arange(self.ntracers):
            titleprefix = 'tracer %d: '%(atracer)
            plt.figure()
            plt.plot(self.outtime/sftime, self.concvar[atracer,:])
            # handle case where a blow-up has occurred
            if self.concvar[atracer,-1] > self.concvar[atracer,0]:
                plt.ylim([0,self.concvar[atracer,0]])
            plt.xlabel('t (days)')
            plt.ylabel('avg($C^2$)')
            plt.title(titleprefix + 'Scalar variance')
            plt.savefig(self.saveloc + '/scalar_variance_tracer%0d.png'%(atracer))

            plt.figure()
            plt.plot(self.outtime/sftime, self.concgrad[atracer,:])
            plt.xlabel('t (days)')
            #plt.ylabel(r"$\sum_i A_i |\nabla C |^2_i / \sum_i A_i$")
            plt.ylabel(r"avg$(|\nabla C |^2)$")
            plt.title(titleprefix + 'Concentration magnitude variance')
            plt.savefig(self.saveloc + '/concentration_magnitude_variance_tracer%0d.png'%(atracer))

            plt.figure()
            plt.plot(self.outtime/sftime, self.concgradmax[atracer,:])
            plt.xlabel('t (days)')
            plt.ylabel(r"max$(|\nabla C |^2)$")
            plt.title(titleprefix + 'Max concentration magnitude')
            plt.savefig(self.saveloc + '/concentration_magnitude_max_tracer%0d.png'%(atracer))

            plt.figure()
            # Klocker 2012 eq 21
            kappatot = -0.5*np.gradient(self.concvar[atracer,:], self.dtout)/self.concgrad[atracer,:]
            plt.plot(self.outtime/sftime, kappatot, '.-')
            plt.xlabel('t (days)')
            plt.ylabel(r"$\kappa_{tot}$")
            plt.title(titleprefix + 'Global diffusivity')
            plt.savefig(self.saveloc + '/global_diffusivity_tracer%0d.png'%(atracer))

        return #}}}

    def plot_initial_positions(self, sflen=1e3, sftime=24.*60*60): #{{{
        """ Plots initial positions on time-varying state of particle vertices """

        def local_plot(vals, label): #{{{
            self.plot_concentration_vertices(vertval=vals/sflen, clims=None, cmap='hsv')
            plt.suptitle(label)
            return #}}}

        local_plot(self.xinit, 'x(t=0)')
        local_plot(self.yinit, 'y(t=0)')
        local_plot(np.sqrt(self.xinit**2 + self.yinit**2), 'r(t=0)')


        # build triangulation based on final time
        tri = GeneralTriangulation(self.tri.x/sflen, self.tri.y/sflen)
        t1 = tri.triangles[:,0]
        t2 = tri.triangles[:,1]
        t3 = tri.triangles[:,2]

        # compute max initial separation on this triangulation
        xinit = self.xinit
        yinit = self.yinit
        d12 = np.sqrt((xinit[t1]-xinit[t2])**2 + (yinit[t1]-yinit[t2])**2)
        d23 = np.sqrt((xinit[t2]-xinit[t3])**2 + (yinit[t2]-yinit[t3])**2)
        d31 = np.sqrt((xinit[t1]-xinit[t3])**2 + (yinit[t1]-yinit[t3])**2)
        dmax = np.max(np.array([d12,d23,d31]),axis=0)/sflen

        # mask
        tri.compute_area()
        tri.set_mask(tri.area*sflen**2 > 5*self.tri.area.mean())

        plt.figure()
        tri.plot_scalar(dmax, cmap='viridis')
        plt.colorbar()
        plt.title("Adjacent particle's maximum distance apart at t=0. t=%.1f days "%(self.t/sftime))

        return #}}}

    def output(self): #{{{
        """
        A good way to use this is to do something like

        >> %time effdiff = ed.EffectiveDiffusivity(np.squeeze(x[rlzn,:,alyr,:]),np.squeeze(y[rlzn,:,alyr,:]), \
        >>                          mask=ds['notoutcropped'][rlzn,0,alyr,:], \
        >>                          tparticles=tvec, meshname=mname, y0=np.linspace(0,2e6,11), \
        >>                          advection=True, diffusion=True, kappa=10., \
        >>                          calctype='ygauss', concwidth=50e3, dt=60., dtout=1.*86400., run=True)
        >> out = effdiff.output()
        >> savenpz('effective_diffusivity_results_layer%02d.npz'%(alyr), **out)

        and then

        >> aggresult = ed.EffectiveDiffusivityAggregation.from_npz_files('effective_diffusivity_results_layer*npz')
        >> aggresult.plot_YB(buoyancySurface, tindex=15)
        """

        outs = dict()
        outs['yloc'] = self.yloc
        outs['mask'] = self.mask

        outs['concvar']     = self.concvar
        outs['concsum']     = self.concsum
        outs['concgrad']    = self.concgrad
        outs['concgradmax'] = self.concgradmax

        outs['Keff']       = self.Keff
        outs['Keffinterp'] = self.Keffinterp
        outs['Keffmax']    = self.Keffmax
        outs['Keffavg']    = self.Keffavg

        outs['kappaeff']       = self.kappaeff
        outs['kappaeffinterp'] = self.kappaeffinterp
        outs['kappaeffavg']    = self.kappaeffavg
        outs['kappaeffmax']    = self.kappaeffmax

        outs['Leffavg']           = self.Leffavg
        outs['Leffmax']           = self.Leffmax
        outs['Leffinterp']        = self.Leffinterp
        outs['LeffStarconstgrad'] = self.LeffStarconstgrad
        outs['LeffStarvargrad']   = self.LeffStarvargrad

        # only need one copy because they are axes or constants (duplicated)
        outs['outtime'] = self.outtime
        outs['caxis']   = self.caxis
        outs['Lmin']    = self.Lmin

        return outs #}}}

    def analysis_plots(self): #{{{
        print 'Ploting standard analysis'
        self.plot_scalar_variance()
        self.plot_effective_diffusivity()
        self.plot_effective_lengths()
        self.plot_conc_area()
        return #}}}

    def setup_concentration(self, xloc=None, yloc=None, concwidth=None, concgrad=None, disttype='ygauss'): #{{{
        """
        disttype corresponds to type of initial distribution profile:

            gauss    -  gaussian profile with concwidth specifying standard
                        deviation centered at xloc, yloc (isotropic)

            ygauss   -  gaussian profile with concwidth specifying standard
                        deviation centered at yloc (constant in x)

            ydelta  -   single strip of 1 concentration nearest to yloc, 0 otherwise

            ygrad   -   largescale gradient with concgrad specifying gradient in y
                        centered at yloc

            const   -   constant value

        """

        tri = self.tri

        if disttype[0] == 'y':
            boundaryedges, centeredges, cellid = \
                    meridionally_constant_distribution_topology(tri, yloc)

        elif disttype == 'gauss':
            boundaryedges, centeredges, cellid = \
                    radially_constant_distribution_topology(tri, xloc, yloc, concwidth)
        elif disttype == 'const':
            #Note, this strictly is not correct
            boundaryedges, centeredges, cellid = meridionally_constant_distribution_topology(tri, yloc)
        else:
            assert False, 'Must specify a known distribution type'

        tri.compute_com()
        tri.compute_le()
        tri.compute_edge_midpoint()

        # compute averaged concentration value on the edges from initial distribution
        # get initial distribution
        if disttype == 'const':
            initconcfunc = lambda x: np.ones_like(x)
            conc = initconcfunc(self.tri.cx)
            qinit = 1.0

        elif disttype == 'gauss':
            initconcfunc = init_gaussian_conc
            conc = initconcfunc(self.tri.cx, self.tri.cy, xloc, yloc, concwidth, self.latlon)
            qinit = np.sum(self.tri.le[boundaryedges]*\
                    initconcfunc(self.tri.ex[boundaryedges], self.tri.ey[boundaryedges], \
                                xloc, yloc, concwidth, self.latlon))\
                                /np.sum(self.tri.le[boundaryedges])

        elif disttype == 'ygauss':
            initconcfunc = init_gaussian_line_conc
            conc = initconcfunc(self.tri.cy, yloc, concwidth)
            qinit = np.sum(self.tri.le[boundaryedges]*\
                    initconcfunc(self.tri.ey[boundaryedges], yloc, concwidth))\
                                /np.sum(self.tri.le[boundaryedges])

        elif disttype == 'ydelta':
            initconcfunc = init_delta_line_conc
            conc = initconcfunc(self.tri.cy, cellid)
            # take the mid-point on each of the boundaryedges (0.5)
            qinit = 0.5

        elif disttype == 'ygrad':
            initconcfunc = init_gradient_line_conc
            conc = initconcfunc(self.tri.cy, concgrad)
            qinit = np.sum(self.tri.le[boundaryedges]*\
                                initconcfunc(self.tri.ey[boundaryedges], concgrad))\
                                /np.sum(self.tri.le[boundaryedges])

        else:
            assert False, 'Must specify known distribution type (disttype) for initial concentration'

        return conc, qinit, centeredges, boundaryedges #}}}

    def compute_time_advancement(self, plotconc=True, ignoreconcbounds=False): #{{{

        tri = self.tri

        if self.t == 0:
            self.timeindex = 0
            tri.compute_area()
            tri.compute_le()
            tri.compute_celledges()
            tri.compute_edgecellneighs()
            tri.compute_celledgedir()
            self.init_diagnostics()
            self.timestep_diagnostics(self.timeindex, plotconc)

        area     = tri.area[np.newaxis,:,np.newaxis]
        celledge = tri.celledges
        celledgedir  = tri.celledgedir[np.newaxis,:,:]
        kappa    = self.kappa
        normal   = tri.return_edge_normal()

        while self.t < self.tend:

            self.t += self.dt
            self.timeindex += 1

            # conc is fixed under advection but advection enhances mixing
            # increasing the gradients, hence it must go first
            if not isxr(self.advection) and (self.advection == True or self.advection == 'lagrangian'):
                tri.update_xy(self.xint(self.t), self.yint(self.t))
            elif self.advection:
                qedge = (np.mean(self.hint(self.t)[tri.edges]*self.uint(self.t)[tri.edges], axis=1)*normal[0,:] +
                         np.mean(self.hint(self.t)[tri.edges]*self.vint(self.t)[tri.edges], axis=1)*normal[1,:])[np.newaxis, :]
                vol = area*np.mean(self.hint(self.t)[tri.triangles], axis=1)[np.newaxis, :, np.newaxis]
                le = tri.le[celledge][np.newaxis,:,:]
                concedge = np.zeros((self.conc.shape[0], celledge.shape[0], celledge.shape[1]))
                edgecfl = qedge*tri.le/(1e-16 + np.mean(vol.squeeze()[tri.edgecellneighs], axis=1))[np.newaxis,:]*self.dt
                maxcfl = np.abs(edgecfl).max()
                assert maxcfl < 0.5, 'Max edge CFL is {}'.format(maxcfl)
                concedge = tri.return_edge_conc(self.conc, edgecfl, method=self.advtype)[:, celledge]
                qedge = qedge[:, celledge]
                div = ne.evaluate('sum(1.0/(vol+1e-16)*(le*qedge*concedge*celledgedir),axis=2)')
                self.conc -= self.dt*div
                # demonstrates rate of diffusion for peak and conservation of mass in domain
                #print self.conc.max(), (self.conc*area).sum()

            if self.diffusion:
                # explicit solve
                #dcdn = tri.return_edge_normal_gradient(self.conc, force=True)[:,celledge]
                #le = tri.le[celledge][np.newaxis,:,:]
                #div = ne.evaluate('sum(kappa/area*(le*dcdn*edgedir),axis=2)')
                #self.conc += self.dt*div

                # fully implicit solve
                _, F = tri.return_edge_normal_gradient_matrix(theta=1.0, kappa=kappa, dt=self.dt)
                #I = sp.eye(F.shape[0], format='csc')
                #solver = sp.linalg.factorized(I - kappa*self.dt*F)
                solver = sp.linalg.factorized(F)
                for ac in np.arange(self.conc.shape[0]):
                    self.conc[ac,:] = solver(self.conc[ac,:])
                    # not accurate enough for bounds preservation with tol=1e-5
                    #self.conc[ac,:], _ = sp.linalg.cg((I - kappa*self.dt*F), self.conc[ac,:])

                ## theta method
                #theta = 1
                #_, F = tri.return_edge_normal_gradient_matrix()
                #I = sp.eye(F.shape[0])
                #for ac in np.arange(self.conc.shape[0]):
                #    self.conc[ac,:] = sp.linalg.spsolve(\
                #            (I - (  theta)*kappa*self.dt*F), \
                #            (I + (1-theta)*kappa*self.dt*F).dot(self.conc[ac,:]))

            assert ignoreconcbounds or self.conc.ravel().max() <= 1.0 and self.conc.ravel().min() >= 0.0, \
                    'Concentration is not bounds-preserving! max=%f min=%f'\
                    %(self.conc.ravel().max(), self.conc.ravel().min())

            if np.mod(self.timeindex, self.ntoutinterval) == 0:
                self.timestep_diagnostics(int(np.floor(self.timeindex/self.ntoutinterval)), plotconc)

            if self.dtreset is not None:
                if np.mod(self.t, self.dtreset) == 0:
                    # note, presumes that self.dt will be some multiple of self.dtreset
                    # tested with an assert above
                    print '\t ...perform the particle reset via remapping approach (assuming x-periodicity)'

                    # get interpolation positions
                    triangles = tri.triangles
                    xtris = fix_periodicity(tri.x[triangles], tri.x[triangles[:,0]][:,np.newaxis], self.x_period)
                    xold = np.mod(np.mean(xtris, axis=1), self.x_period)
                    yold = np.mean(tri.y[triangles], axis=1)

                    xtris = fix_periodicity(self.xp0[triangles], self.xp0[triangles[:,0]][:,np.newaxis], self.x_period)
                    xnew = np.mod(np.mean(xtris, axis=1), self.x_period)
                    ynew = np.mean(self.yp0[triangles], axis=1)

                    # note, there may be some potential aliasing on the
                    # boundaries, but the concentration and gradient fields
                    # appear to be good so if this is occuring it is subtle and
                    # probably not an issue

                    # concs
                    for ac in np.arange(self.conc.shape[0]):
                        cold = self.conc[ac,:]
                        #cnew = unstructured_interp(xnew, ynew, xold, yold, cold, extrapolate=True)
                        cnew = unstructured_remap(xnew, ynew, xold, yold, cold)
                        self.conc[ac,:] = cnew

                    # update xy
                    tri.update_xy(self.xp0, self.yp0)

        # handle finalization
        print 'Concentration evolution is finished at %e days'%(self.t/86400.)
        self.finalize_diagnostics()

        return #}}}

    # diagnostic methods #{{{
    def init_diagnostics(self): #{{{

        # initialize structures to compute area-coordinate integrals
        self.area = []
        self.gradc2 = []
        self.ybar = np.zeros((self.ntracers, self.ntout, CummulativeIntegral().bins.shape[0]-1))
        for atracer in np.arange(self.ntracers):
            self.area.append(CummulativeIntegral(ntimes=self.ntout))
            self.gradc2.append(CummulativeIntegral(ntimes=self.ntout))

        self.LeffStarvargrad = np.zeros((self.ntracers, self.ntout))
        self.LeffStarconstgrad = np.zeros((self.ntracers, self.ntout))

        self.concvar     = np.zeros((self.ntracers, self.ntout))
        self.concgrad    = np.zeros((self.ntracers, self.ntout))
        self.concgradmax = np.zeros((self.ntracers, self.ntout))
        self.concsum     = np.zeros((self.ntracers, self.ntout))
        self.outtime     = np.zeros(( self.ntout))

        self.Lmin = np.zeros((self.ntracers))

        return #}}}

    def timestep_diagnostics(self, timeindex, plotconc): #{{{
        """
        Computes integral sums, e.g.,

        area =          int_{c<C} dA
        gradc2 =        int_{c<C} |grad c|^2 dA

        """
        print 'Performing diagnostics at time index = %d for t = %f days at %s'%\
                (self.timeindex, self.t/86400., datetime.datetime.now())
        tri = self.tri

        # check to make sure dt is reasonable
        if self.diffusion:
            dtmin = 2./9.*(tri.area.min()/tri.le.max())**2.0/self.kappa
            dtmean = 2./9.*(tri.area.mean()/tri.le.mean())**2.0/self.kappa
            dtokmean = self.dt < dtmean
            dtokmin = self.dt < dtmin
            if self.t == 0:
                if dtokmean:
                    print 'dt = %e s is too large compared to min of %e s for explicit solver, reduce by %d'\
                        %(self.dt, dtmin, int(np.ceil(self.dt/dtmin)))
            if not dtokmin:
                print 'dt = %e s may be too large compared to min of %e s, reduce by %d'\
                        %(self.dt, dtmin, int(np.ceil(self.dt/dtmin)))
                print 'dt = %e s may be too large compared to mean of %e s, reduce by %d'\
                        %(self.dt, dtmean, int(np.ceil(self.dt/dtmean)))

        self.outtime[timeindex] = self.t

        # note, also computes tri.dcdn on edges, which is needed for effective lengths
        #gradc2 = np.sum(tri.return_cell_gradient(self.conc,funcname='LSQ', force=True)**2.0,axis=0)
        gradc2 = tri.return_cell_gradient_magnitude(self.conc)

        for atracer in np.arange(self.ntracers):
            le = tri.le[tri.boundaryedges[atracer]]

            if timeindex == 0:
                # note, factor of two can occur because of flux out of top / bottom of cells for a strip (contour line)
                # because of finite width
                self.Lmin[atracer] = np.sum(le)

            # compute area-coordinate integral sums for effective diffusivity calc
            self.area[atracer].compute_integral_sum(tri.area, self.conc[atracer,:], timeindex)
            self.gradc2[atracer].compute_integral_sum(gradc2[atracer,:]*tri.area, self.conc[atracer,:], timeindex)
            self.ybar[atracer, timeindex, :] = CummulativeIntegral.binned_mean(np.mean(tri.y[tri.triangles],axis=1),
                                               self.conc[atracer,:], self.area[atracer].bins)

            # compute area-normalized conc variance (should decrease over time)
            weighted_avg = lambda val: np.sum(val*tri.area)/np.sum(tri.area)
            self.concvar[atracer, timeindex]     = weighted_avg(self.conc[atracer,:]**2.0)
            self.concsum[atracer, timeindex]     = weighted_avg(self.conc[atracer,:])
            self.concgrad[atracer, timeindex]    = weighted_avg(gradc2[atracer,:])
            self.concgradmax[atracer, timeindex] = np.max(gradc2[atracer,:])

            dcdn = tri.return_edge_normal_gradient(self.conc)
            absgrad = np.abs(dcdn[atracer,tri.boundaryedges[atracer]])
            # note, have excluded boundary edges from the computation
            if np.min(absgrad) == 0:
                idx = np.where(absgrad == 0)
                absgrad[idx] = 1
                print 'Warning! Zero concentration gradient obtained.  This will ' + \
                'contaminate LeffStarvargrad with nans! Replacing %d 0s with 1s.'%(len(idx[0]))

            # compute effective lengths using edge masks previously defined tri.boundaryedges
            self.LeffStarvargrad[atracer,timeindex]   = np.sqrt(np.sum(le/absgrad)*np.sum(le*absgrad))
            self.LeffStarconstgrad[atracer,timeindex] = np.sqrt(np.sum(le)**2)

            if plotconc:
                # note, FIX! save convention here is sketchy because each call to plot_concentration_vertices doesn't produce a unique save file
                self.plot_concentration_vertices(trival=self.conc[atracer,:], titleprefix='tracer %d: '%(atracer), \
                        simpleplot=self.simpleplot, save='conc_vert_tracer%04d_'%(atracer))
                np.savez(self.saveloc + 'conc_t%04d_days_tracer%04d.npz'%(self.t/86400., atracer), area=self.conc[atracer,:])
                plt.savefig(self.saveloc + 'conc_vertices_t%04d_days.png'%(self.t/86400.))
                # plot concentration gradient at vertices
                np.savez(self.saveloc + 'gradc2tris_t%04d_days_tracer%04d.npz'%(self.t/86400., atracer), area=gradc2[atracer,:])
                self.plot_concentration_vertices(trival=gradc2[atracer,:], titleprefix='tracer %d: '%(atracer), simpleplot=self.simpleplot,  clims=None, save='gradc2vert_tracer%04d_'%(atracer))
                plt.savefig(self.saveloc + 'gradc2_vertices_t%04d_days.png'%(self.t/86400.))

        if plotconc:
            # save computed areas
            triareas = tri.compute_area(force=True, output=True)
            np.savez(self.saveloc + 'areas_t%04d_days.npz'%(self.t/86400.), x=tri.x, y=tri.y, triangles=tri.triangles, area=triareas)
            self.plot_concentration_vertices(trival=triareas, titleprefix='', simpleplot=self.simpleplot, clims=None, save='area_vert_tracer%04d_'%(atracer))
            plt.savefig(self.saveloc + 'areas_vertices_t%04d_days.png'%(self.t/86400.))

        return #}}}

    def finalize_diagnostics(self): #{{{
        """
        Computes dA/dc and dS/dc to obtain kappaeff = kappa * dA/dc * dS/dc

            where S = int_{c < C} | grad c|^2 dA
        """

        Nc = self.area[0].y.shape[1]-1
        Nt = self.area[0].y.shape[0]

        self.kappaeff       = np.nan*np.zeros((self.ntracers, Nt, Nc))
        self.kappaeffinterp = np.nan*np.zeros((self.ntracers, Nt))
        self.kappaeffavg    = np.nan*np.zeros((self.ntracers, Nt))
        self.kappaeffmax    = np.nan*np.zeros((self.ntracers, Nt))

        smoothnc = int(self.area[0].nbins/10.)
        smoothnc += 1 - np.mod(smoothnc,2)
        # assumes that time output will be fairly coarse, so minimal or no smoothing in time
        smoothnt = 1
        smoothing = (smoothnc,smoothnt)

        for atracer in np.arange(self.ntracers):
            dareaAdc = self.area[atracer].compute_xderivative(smoothing=smoothing)
            dgradCdc = self.gradc2[atracer].compute_xderivative(smoothing=smoothing)

            self.kappaeff[atracer, :,:] = self.kappa*dareaAdc*dgradCdc

            # Thus, we can see that there are three ways to enhance mixing:
            # 1. background diffusivity
            # 2. concentration distribution (via diffusion) which results in a change in dareaAdc
            # 3. Stretching by advection that enhances gradients in the fluid in dgradCdc

            # get kappaeff at self.qinit to compute Leff = sqrt(kappaeff/kappa)
            area = 0.5*(self.area[atracer].yint[:,1:] + self.area[atracer].yint[:,:-1])
            xc = 0.5*(self.area[atracer].bins[:-1] + self.area[atracer].bins[1:])
            self.caxis = 0.5*(xc[:-1]+xc[1:])
            # handle out of range queries by snapping to the nearest points in histogram
            if self.qinit[atracer] > self.caxis[-1] and self.qinit[atracer] <= 1.0:
                self.kappaeffinterp[atracer,:] = self.kappaeff[atracer, :,-1]
            elif self.qinit[atracer] < self.caxis[0] and self.qinit[atracer] >= 0.0:
                self.kappaeffinterp[atracer, :] = self.kappaeff[atracer, :,0]
            else:
                # let nans be returned for queries outside range
                self.kappaeffinterp[atracer, :] = \
                        interp1d(self.caxis, self.kappaeff[atracer, :,:], axis=1, bounds_error=False)(self.qinit[atracer])

            self.kappaeffavg[atracer, :] = np.sum(self.kappaeff[atracer, :,:] * area, axis=1) / np.sum(area, axis=1)
            self.kappaeffmax[atracer, :] = np.max(self.kappaeff[atracer, :,:], axis=1)

        self.Keff       = self.kappaeff       / self.Lmin[:,np.newaxis,np.newaxis]**2.
        self.Keffavg    = self.kappaeffavg    / self.Lmin[:,np.newaxis]**2.
        self.Keffinterp = self.kappaeffinterp / self.Lmin[:,np.newaxis]**2.
        self.Keffmax    = self.kappaeffmax    / self.Lmin[:,np.newaxis]**2.

        self.Leffinterp = np.sqrt(self.kappaeffinterp / self.kappa)
        self.Leffavg    = np.sqrt(self.kappaeffavg    / self.kappa)
        self.Leffmax    = np.sqrt(self.kappaeffmax    / self.kappa)

        return #}}}
    #}}}

    #}}}

# class testing functions #{{{
def test_init_conc(effdiff): #{{{

    effdiff.plot_concentration_vertices(trival=effdiff.conc[0,:])
    plt.savefig('test_init_conc_vertices.png')

    effdiff.plot_concentration_triangles(trival=effdiff.conc[0,:])
    plt.savefig('test_init_conc_triangles.png')

    return #}}}

def test_edge_masks(effdiff): #{{{

    # test edge masks
    tri = effdiff.tri

    # get tri.ex, tri.ey
    tri.compute_edge_midpoint()
    tri.compute_le()

    plt.figure()
    plt.tripcolor(tri.x, tri.y, tri.triangles, facecolors=np.zeros(tri.triangles.shape[0]), edgecolors='k', alpha=0.3)
    plt.plot(tri.ex[tri.centeredges[0]],tri.ey[tri.centeredges[0]],'r.',alpha=0.5)
    plt.plot(tri.ex[tri.boundaryedges[0]],tri.ey[tri.boundaryedges[0]],'g.',alpha=0.5)
    dl = np.sqrt(3)/2*tri.le.mean()
    #plt.axis('equal')
    plt.ylim([effdiff.yloc-3*dl,effdiff.yloc+3*dl])
    plt.savefig('test_edge_masks.png')

    return #}}}

def run_tests(mname='mesh.nc', xloc=500.e3, yloc=1500.e3): #{{{

    print 'Testing Eulerian scalar advection...'
    mname = 'mesh_20km.nc'
    dt = 1.*60*60.
    ndays = 10.
    saveloc = './test_eulerian_advection/'
    dt = 0.1*60*60.
    msh = GeneralTriangulation.from_MPAS_mesh(meshname=mname)
    tvec = np.arange(0,ndays)*86400.
    msh.buoyancySurfaceVelocityZonal = 1*np.ones_like(msh.x)[np.newaxis,:]*np.ones_like(tvec[:,np.newaxis])
    msh.buoyancySurfaceVelocityMeridional = np.zeros_like(msh.x)[np.newaxis,:]*np.ones_like(tvec[:,np.newaxis])
    effdiff = EffectiveDiffusivity(meshname=mname, tparticles=tvec,
                                   kappa=0., calctype='gauss', advtype='laxwendroff',
                                   advection = msh, diffusion=False, y0=1500e3, x0=500e3,
                                   dt=dt, dtout=1.*86400., dtreset=None,
                                   saveloc=saveloc)
    #effdiff = EffectiveDiffusivity(meshname=mname, tparticles=tvec,
    #                               kappa=0., calctype='ygrad',
    #                               advection = msh, advtype='laxwendroff',
    #                               diffusion=False, y0=1500e3, x0=500e3,
    #                               dt=dt, dtout=1.*86400., dtreset=None,
    #                               saveloc=saveloc)
    effdiff.compute_time_advancement()
    #effdiff.plot_scalar_variance()
    #effdiff.plot_effective_diffusivity()
    #effdiff.plot_effective_lengths()
    #from cumulative_conc_eff_diff import compute_eff_diff, test_plot_eff_diff
    #compute_eff_diff(saveloc, nfiles=ndays, nstart=0)
    #test_plot_eff_diff(saveloc)
    #print 'Max Keff is %d'.format(np.abs(effdiff.Keff).max())
    print '... finished testing Eulerian scalar advection'
    import pdb; pdb.set_trace()

    # test functionality of effective diffusivity class to compute desired quantities
    print 'Testing isotropic gaussian...'
    effdiff = EffectiveDiffusivity(xt=None, yt=None, tparticles=None, \
            meshname=mname, x0=xloc, y0=yloc, calctype='gauss')
    test_init_conc(effdiff)
    test_edge_masks(effdiff)
    plt.show()

    print 'Testing gradient...'
    for y0 in [2e6*np.random.rand(1)]:
        effdiff = EffectiveDiffusivity(xt=None, yt=None, tparticles=None, meshname=mname, y0=y0,
                concgrad=0.999999/2e6, calctype='ygrad')
        test_init_conc(effdiff)
        test_edge_masks(effdiff)
    plt.show()

    print 'Testing ydelta...'
    for y0 in [0, 2e6, 50e3, 2e6*np.random.rand(1)]:
        effdiff = EffectiveDiffusivity(xt=None, yt=None, tparticles=None, meshname=mname, y0=y0, calctype='ydelta')
        test_init_conc(effdiff)
        test_edge_masks(effdiff)
    plt.show()

    print 'Testing gaussians...'
    for y0 in [0, 2e6, 50e3, 2e6*np.random.rand(1)]:
        effdiff = EffectiveDiffusivity(xt=None, yt=None, tparticles=None, meshname=mname, y0=y0, calctype='ygauss')
        test_init_conc(effdiff)
        test_edge_masks(effdiff)
    plt.show()

    # tests capability CummulativeIntegral class to compute effective diffusivity from data
    test_cummulative_integral()
    plt.show()

    return #}}}

#}}}
#}}}

class EffectiveDiffusivityAggregation: #{{{
    # constructors #{{{
    def __init__(self, npzdir, copyonly=['outtime','caxis']): #{{{
        import glob
        import re

        globs = npzdir + '/' + 'eff*npz'
        files = glob.glob(globs)

        # get range of rlzns, split, and layers
        rlzns = np.unique([int(re.findall( '_rlzn_[0-9]+_',afile)[0][6:9]) for afile in files])
        splits = np.unique([int(re.findall('_split_[0-9]+_',afile)[0][7:11]) for afile in files])
        layers = np.unique([int(re.findall('_layer[0-9]+',afile)[0][6:8]) for afile in files])

        rlzndat = []
        for arlzn in rlzns:
            # join along splits
            splitdat = []
            for asplit in splits:
                layerdat = []
                for alayer in layers:
                    afile = npzdir + '/' + 'effective_diffusivity_rlzn_%03d_split_%04d_results_layer%02d.npz'%(arlzn, asplit, alayer)
                    layerdat.append(np.load(afile))
                # join along layers
                splitdat.append(agg_cat_general(layerdat,axis=None, copyonly=copyonly + ['rlzn']))
            rlzndat.append(agg_cat_general(splitdat,axis=1, copyonly=copyonly + ['rlzn']))
        alldat = agg_cat_general(rlzndat,axis=None, copyonly=copyonly)

        self.aggresult = alldat

        return #}}}
    #}}}

    # methods #{{{
    def plot_YB(self, buoyancySurface, tindex=-1, indays=True, saveprefix='./'): #{{{

        # scale factors #{{{

        if indays:
            sftime = 24*60*60.
        else:
            sftime = 1.0

        #}}}

        masknan = np.asarray(self.aggresult['mask'] > (1-1e-4), dtype='f8')
        masknan[np.where(np.logical_not(masknan))] = np.nan


        # Keff #{{{
        plt.figure()

        palette = plt.cm.magma
        palette.set_bad('w',1.0)
        A = np.mean(np.ma.array(np.log10(masknan*self.aggresult['Keffavg'][:,:,:,tindex]), mask=np.isnan(masknan)),axis=0)

        contourrange = np.linspace(0,3,13) #np.linspace(0,6,4*6+1) #
        contourrange = np.linspace(0,6,4*6+1) #
        x = (self.aggresult['yloc'][0,:]/1e3)*np.ones_like(buoyancySurface[:,np.newaxis])
        y = (buoyancySurface -1029)[:,np.newaxis]*np.ones_like((self.aggresult['yloc'][0,:]/1e3))
        plt.contourf(x, y, A, contourrange, cmap=palette)
        #plt.colorbar(ticks=list(np.linspace(0,3,4)))
        plt.colorbar(ticks=list(np.linspace(0,6,7)))
        cs = plt.contour(x, y, A, colors='k', levels=contourrange, lw=2, alpha=0.3)
        plt.clabel(cs, inline=1, fontsize=6, alpha=0.3, fmt='%.2f')

        plt.xlabel('Y (km)')
        plt.ylabel('Potential density - 1029.0 (kg m$^{-3}$)')
        plt.gca().invert_yaxis()
        plt.axis('tight')
        plt.title('$\log_{10}{K_\mathrm{eff}}$ (m$^2$ s$^{-1}$) at t=%.1f days'%(self.aggresult['outtime'][tindex]/sftime))
        plt.savefig(saveprefix + 'Keff.png')
        #}}}


        # Leff #{{{
        plt.figure()

        palette = plt.cm.viridis
        palette.set_bad('w',1.0)
        A = np.mean(np.ma.array(np.log10(masknan*self.aggresult['Leffavg'][:,:,:,tindex]/1e3), mask=np.isnan(masknan)),axis=0)
        contourrange = np.linspace(3,5,11)
        x = (self.aggresult['yloc'][0,:]/1e3)*np.ones_like(buoyancySurface[:,np.newaxis])
        y = (buoyancySurface -1029)[:,np.newaxis]*np.ones_like((self.aggresult['yloc'][0,:]/1e3))
        plt.contourf(x, y, A, contourrange, cmap=palette)
        plt.colorbar(ticks=list(np.linspace(3,5,3)))
        cs = plt.contour(x, y, A, colors='k', levels=contourrange, lw=2, alpha=0.3)
        plt.clabel(cs, inline=1, fontsize=6, alpha=0.3, fmt='%.2f')

        plt.clim(3,5)
        plt.xlabel('Y (km)')
        plt.ylabel('Potential density - 1029.0 (kg m$^{-3}$)')
        plt.axis('tight')
        plt.gca().invert_yaxis()
        plt.title('$\log_{10}(L_\mathrm{eff})$ (km) at t=%.1f days'%(self.aggresult['outtime'][tindex]/sftime))
        plt.savefig(saveprefix + 'Leff.png')
        #}}}

        # Leff/Lmin #{{{
        plt.figure()
        palette = plt.cm.viridis
        palette.set_bad('w',1.0)
        A = np.mean(np.ma.array((masknan*self.aggresult['Leffavg'][:,:,:,tindex]/self.aggresult['Lmin'][:,:]), mask=np.isnan(masknan)),axis=0)

        contourrange = np.linspace(0,6,13)
        x = (self.aggresult['yloc'][0,:]/1e3)*np.ones_like(buoyancySurface[:,np.newaxis])
        y = (buoyancySurface -1029)[:,np.newaxis]*np.ones_like((self.aggresult['yloc'][0,:]/1e3))
        plt.contourf(x, y, A, contourrange, cmap=palette)
        plt.colorbar(ticks=list(np.linspace(0,6,7)))
        cs = plt.contour(x, y, A, colors='k', levels=contourrange, lw=2, alpha=0.3)
        plt.clabel(cs, inline=1, fontsize=6, alpha=0.3, fmt='%.2f')
        plt.clim(0,5)
        plt.xlabel('Y (km)')
        plt.ylabel('Potential density - 1029.0 (kg m$^{-3}$)')
        plt.axis('tight')
        plt.gca().invert_yaxis()
        plt.title('$L_\mathrm{eff} / L_\mathrm{min}$ at t=%.1f days'%(self.aggresult['outtime'][tindex]/sftime))
        plt.savefig(saveprefix + 'Leff_div_Lmin.png')
        #}}}

        return #}}}

    def plot_effective_lengths(self, indays=True, inkm=True, titleprefix='', simpleplot=False): #{{{
        """
        Plot results from EffectiveDiffusivity.output()
        Presumes on results from list of outputed dict of EffectiveDiffusivity.output(), e.g.,

        outs['yloc'] = self.yloc

        outs['concvar'] = self.concvar
        outs['concsum'] = self.concsum
        outs['concgrad'] = self.concgrad
        outs['concgradmax'] = self.concgradmax

        outs['Keff'] = self.Keff
        outs['Keffinterp'] = self.Keffinterp
        outs['Keffavg'] = self.Keffavg

        outs['Leffavg'] = self.Leffavg
        outs['Leffinterp'] = self.Leffinterp
        outs['LeffStarconstgrad'] = self.LeffStarconstgrad
        outs['LeffStarvargrad'] = self.LeffStarvargrad
        outs['Lmin'] = self.Lmin
        """

        # scale factors #{{{

        if indays:
            sftime = 24*60*60.
        else:
            sftime = 1.0

        if inkm:
            sflen = 1e3
        else:
            sflen = 1.0

        #}}}

        # get variables from results file #{{{

        # constant across list of results
        yloc = self.aggresult['yloc']
        outtime = self.aggresult['outtime']
        caxis = self.aggresult['caxis']
        Lmin = self.aggresult['Lmin']

        Leffinterp = self.aggresult['Leffinterp']
        Leffavg = self.aggresult['Leffavg']
        LeffStarvargrad = self.aggresult['LeffStarvargrad']
        LeffStarconstgrad = self.aggresult['LeffStarconstgrad']

        #}}}

        # LeffStarconstgrad (L_particle no conc) #{{{

        plt.figure()
        plt.contourf(yloc/sflen, outtime/sftime, LeffStarconstgrad.T/sflen, cmap='viridis')
        plt.colorbar()
        if inkm:
            plt.xlabel('y (km)')
            plt.title(titleprefix + r"$L_\mathrm{particle,\,no\,conc}$ km")
        else:
            plt.xlabel('y (m)')
            plt.title(titleprefix + r"$L_\mathrm{particle,\,no\,conc}$ m")
        if indays:
            plt.ylabel('t (days)')
        else:
            plt.ylabel('t (s)')

        #}}}

        if simpleplot:
            return

        # LeffStarvargrad L_particle conc grad #{{{

        plt.figure()
        plt.contourf(yloc/sflen, outtime/sftime, LeffStarvargrad.T/sflen, cmap='viridis')
        plt.colorbar()
        if inkm:
            plt.xlabel('y (km)')
            plt.title(titleprefix + r"$L_\mathrm{particle,\,conc\,grad}$ km")
        else:
            plt.xlabel('y (m)')
            plt.title(titleprefix + r"$L_\mathrm{particle,\,conc\,grad}$ m")
        if indays:
            plt.ylabel('t (days)')
        else:
            plt.ylabel('t (s)')

        #}}}

        # Leffavg L_eff, avg #{{{

        plt.figure()
        plt.contourf(yloc/sflen, outtime/sftime, Leffavg.T/sflen, cmap='viridis')
        plt.colorbar()
        if inkm:
            plt.xlabel('y (km)')
            plt.title(titleprefix + r"$L_\mathrm{eff,\,avg}$ km")
        else:
            plt.xlabel('y (m)')
            plt.title(titleprefix + r"$L_\mathrm{eff,\,avg}$ m")
        if indays:
            plt.ylabel('t (days)')
        else:
            plt.ylabel('t (s)')

        #}}}

        # Leffinterp L_interp #{{{

        plt.figure()
        plt.contourf(yloc/sflen, outtime/sftime, Leffinterp.T/sflen, cmap='viridis')
        plt.colorbar()
        if inkm:
            plt.xlabel('y (km)')
            plt.title(titleprefix + r"$L_\mathrm{interp}$ km")
        else:
            plt.xlabel('y (m)')
            plt.title(titleprefix + r"$L_\mathrm{interp}$ m")
        if indays:
            plt.ylabel('t (days)')
        else:
            plt.ylabel('t (s)')

        #}}}

        return #}}}

    def plot_effective_diffusivity(self, indays=True, inkm=True, titleprefix='', simpleplot=False): #{{{
        """
        Plot results from EffectiveDiffusivity.output()
        Presumes on results from list of outputed dict of EffectiveDiffusivity.output(), e.g.,

        outs['yloc'] = self.yloc

        outs['concvar'] = self.concvar
        outs['concsum'] = self.concsum

        outs['Keff'] = self.Keff
        outs['Keffinterp'] = self.Keffinterp
        outs['Keffavg'] = self.Keffavg

        outs['Leffavg'] = self.Leffavg
        outs['Leffinterp'] = self.Leffinterp
        outs['LeffStarconstgrad'] = self.LeffStarconstgrad
        outs['LeffStarvargrad'] = self.LeffStarvargrad
        outs['Lmin'] = self.Lmin
        """

        # scale factors #{{{

        if indays:
            sftime = 24*60*60.
        else:
            sftime = 1.0

        if inkm:
            sflen = 1e3
        else:
            sflen = 1.0

        #}}}

        # get variables from results file #{{{

        # constant across list of results
        yloc = self.aggresult['yloc']
        outtime = self.aggresult['outtime']
        caxis = self.aggresult['caxis']

        # variable across list of results
        concvar = self.aggresult['concvar']

        Keff = self.aggresult['Keff']
        Keffinterp = self.aggresult['Keffinterp']
        Keffavg = self.aggresult['Keffavg']
        Lmin = self.aggresult['Lmin']

        #}}}

        # Keffavg K_eff average #{{{

        plt.figure()
        plt.contourf(yloc/sflen, outtime/sftime, Keffavg.T, cmap='viridis')
        plt.colorbar()
        plt.xlabel('y (km)')
        plt.title(titleprefix + r"$K_\mathrm{eff}$ average m$^2$s$^{-1}$")
        if indays:
            plt.ylabel('t (days)')
        else:
            plt.ylabel('t (s)')

        #}}}

        if simpleplot:
            return


        # Keffinterp K_eff interpolated #{{{

        plt.figure()
        plt.contourf(yloc/sflen, outtime/sftime, Keffinterp.T, cmap='viridis')
        plt.colorbar()
        plt.xlabel('y (m)')
        plt.title(titleprefix + r"$K_\mathrm{eff}$ interpolated m$^2$ s$^{-1}$")
        if indays:
            plt.ylabel('t (days)')
        else:
            plt.ylabel('t (s)')

        #}}}

        return #}}}
    #}}}
#}}}

if __name__ == "__main__":

    run_tests()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
