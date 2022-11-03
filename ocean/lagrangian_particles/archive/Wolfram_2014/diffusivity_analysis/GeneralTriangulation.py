#!/usr/bin/env python
# Phillip J. Wolfram, 2014
import numpy as np
from numpy.linalg import lstsq
import numexpr as ne
import operator
import copy
import hashlib
import matplotlib as mplt
from matplotlib.tri import Triangulation, LinearTriInterpolator, CubicTriInterpolator
import matplotlib.pyplot as plt
from scipy import spatial
import scipy.sparse as sp
import pyamg
from latlon_coordinate_transforms import fix_periodicity_numexpr as fix_periodicity,  \
        proj_xyz_numexpr as proj_xyz, proj_lat_long_numexpr as proj_lat_long, rEarth, \
        mid_point_numexpr as mid_point, haversine_formula_numexpr as haversine_formula \

def edges_containing_points(points, edges): #{{{
    """ returns index to edges with indices in points """

    point1 = edges[:,0]
    point2 = edges[:,1]

    validedges = np.logical_and(np.in1d(point1,points),np.in1d(point2,points))

    return validedges #}}}

def triangles_continaing_points(points, triangles): #{{{
    """ returns index to triangles with indices in points """

    point1 = triangles[:,0]
    point2 = triangles[:,1]
    point3 = triangles[:,2]

    pintri1 = np.in1d(point1, points)
    pintri2 = np.in1d(point2, points)
    pintri3 = np.in1d(point3, points)

    validtri = np.logical_and(pintri1, np.logical_and(pintri2, pintri3))

    return validtri #}}}

def unstructured_remap(xnew, ynew, xold, yold, cold, method='linear'): #{{{
    """
    Remaps set of points to nearest neighbors, ensuring one-to-one behavior

    Phillip J. Wolfram
    07/22/2016
    """

    #tree = spatial.cKDTree(np.vstack((xold, yold)).T)

    #dist, idx = tree.query(np.vstack((xnew, ynew)).T, k=1, eps=0)
    #
    #dups = np.sort(idx , axis=None)
    #dups = dups[dups[1:] == dups[:-1]]
    #unique = np.setdiff1d(idx, dups)

    #uidx = np.in1d(idx, unique)


    #cnew = np.nan*np.zeros_like(cold)
    #cnew[uidx] = 1 #cold[idx[uidx]]

    #cnew = cold[idx]

    cstart = unstructured_interp(xnew, ynew, xold, yold, cold, method=method, extrapolate=True)
    idx = np.argsort(cstart)
    cnew = np.zeros_like(cold)*np.nan
    cnew[idx] = np.sort(cold)

    return cnew #idx, np.where(uidx)[0] #}}}

def unstructured_interp(xnew, ynew, xold, yold, cold, method='linear', extrapolate=False): #{{{
    """
    Performs an unstructured interpolate from unstructured 2D data
    to unstructured 2D data from cold data on (xold, yold) points
    to (xnew, ynew) points.

    Phillip J. Wolfram
    06/30/2016
    """
    if method == 'nearest':

        tree = spatial.cKDTree(np.vstack((xold,yold)).T)
        _, idx = tree.query(np.vstack((xnew, ynew)).T)
        cnew = cold[idx]

    elif method == 'linear':

        try:
          # build up triangular interpolator
          ti = LinearTriInterpolator(Triangulation(xold, yold), cold)
          cnew = ti(xnew, ynew).data
          if extrapolate:
              tree = spatial.cKDTree(np.vstack((xold,yold)).T)
              # handle case of out-of-bounds
              # via nearest neighbor interp
              nanvals = np.isnan(cnew)
              _, nanid = tree.query(np.vstack((xnew[nanvals], ynew[nanvals])).T)
              cnew[nanvals] = cold[nanid]
        except RuntimeError as err:
            print '\t\tHad an error with "%s" when using linear interpolation, trying nearest interpolation instead.'%(err)
            return unstructured_interp(xnew, ynew, xold, yold, cold, method='nearest')
    else:
        assert False, 'Method %s is not known'%(method)

    return cnew #}}}

def subdivide_triangulation(xin, yin, triin, subdivtype='edge'): #{{{
    """
    Subdivide a given triangulation specified by (xin, yin) points and
    triin trinagles.

    Phillip J. Wolfram
    07/05/2016
    """
    if subdivtype == 'center':
        xout, yout, triout = subdivide_triangulation_center(xin, yin, triin)
    elif subdivtype == 'edge':
        xout, yout, triout = subdivide_triangulation_edge(xin, yin, triin)

    return xout, yout, triout #}}}

def subdivide_triangulation_edge(xin, yin, triin): #{{{
    """
    Subdivide a given triangulation specified by (xin, yin) points and
    triin trinagles using edge centers as new points.

    Phillip J. Wolfram
    07/05/2016
    """

    # get triangulation for data structures
    tri = GeneralTriangulation(x=xin, y=yin, tri=triin)

    # cell centers that will be new set of points
    tri.compute_edge_midpoint()
    xe = tri.ex.copy()
    ye = tri.ey.copy()

    # new set of points
    npold = xin.shape[0]
    xout = np.hstack((xin, xe))
    yout = np.hstack((yin, ye))

    # new triangles are formed from cell-neighbors (points above ) to each edge and its points
    ntold = triin.shape[0]
    triout = np.zeros((4*ntold, 3))

    def compute_ordering(xT, yT, nT):
        v0 = np.array([(xT[nT[1]] - xT[nT[0]]), (yT[nT[1]] - yT[nT[0]])])
        v1 = np.array([(xT[nT[2]] - xT[nT[0]]), (yT[nT[2]] - yT[nT[0]])])
        if np.cross(v0,v1) > 0:
            # CCW orientation
            return np.array([nT[0], nT[1], nT[2]])
        else:
            return np.array([nT[0], nT[2], nT[1]])

    tri.compute_celledges()
    for atri in np.arange(ntold):
        # center
        triout[atri*4 + 0, :] = compute_ordering(xout, yout, tri.celledges[atri, :] + npold)

        # edges
        for trinum, (e1, e2) in enumerate(zip(tri.celledges[atri,:], np.roll(tri.celledges[atri,:],1))):
            # get common node
            n0 = np.intersect1d(tri.edges[e1,:], tri.edges[e2,:])
            # form triangle
            triout[atri*4 + trinum + 1, :] = compute_ordering(xout, yout, np.array([e1 + npold, e2 + npold, n0]))

    return xout, yout, triout #}}}

def subdivide_triangulation_center(xin, yin, triin): #{{{
    """
    Subdivide a given triangulation specified by (xin, yin) points and
    triin trinagles using cell center as new point.

    Phillip J. Wolfram
    07/05/2016
    """

    # get triangulation for data structures
    tri = GeneralTriangulation(x=xin, y=yin, tri=triin)

    # cell centers that will be new set of points
    tri.compute_com()
    xc = tri.cx.copy()
    yc = tri.cy.copy()

    # new set of points
    npold = xin.shape[0]
    xout = np.hstack((xin, xc))
    yout = np.hstack((yin, yc))

    # new triangles are formed from cell-neighbors (points above ) to each edge and its points
    ntold = triin.shape[0]
    triout = np.zeros((3*ntold, 3))

    tri.compute_celledges()
    tri.compute_edgecellneighs()

    atri = 0
    for theedge in np.arange(tri.edges.shape[0]):

        # cells and nodes for new triangles
        n0 = tri.edges[theedge,0]
        n1 = tri.edges[theedge,1]
        c0 = tri.edgecellneighs[theedge,0]
        c1 = tri.edgecellneighs[theedge,1]

        for center in [c0, c1]:
            if center > -1:
                # compute ordering
                v0 = np.array([(tri.x[n1] - tri.x[n0]), (tri.y[n1] - tri.y[n0])])
                v1 = np.array([(xc[center] - tri.x[n0]), (yc[center] - tri.y[n0])])
                if np.cross(v0,v1) > 0:
                    # CCW orientation
                    triout[atri, :] = np.array([n0, n1, center + npold])
                else:
                    triout[atri, :] = np.array([n1, n0, center + npold])
                atri += 1

    return xout, yout, triout#}}}

class GeneralTriangulation(Triangulation): #{{{

    def __init__(self, x=None, y=None, tri=None, scalar=None, latlon=False, *args, **kwargs): #{{{

        # extend
        self.meshname = None        # name of MPAS mesh that produced this triangulation
        self.area = None            # area of triangles
        self.weightarea = None      # apprortioned area of dual (non-exact)
        self.cx = None              # triangle center of mass y coord
        self.cy = None              # triangle center of mass x coord
        self.ex = None              # edge midpoint x coord
        self.ey = None              # edge midpoint y coord
        self.le = None              # length of edges
        self.celledges = None       # edges that make up a cell
        self.edgecellneighs = None  # cell neighbors to an edge
        self.nodeedgeneighs = None  # edge neighbors to a node
        self.celledgedir = None     # +1 if normal into cell on edge, -1 if normal out of cell
        self.interioredges = None
        self.areagradcoeff = None
        self.balls = None
        self.ballsradius = None
        self.interpmethod = LinearTriInterpolator
        self.interpolator = None
        self.interphash= None
        self.periodic = None
        self.Lx = None
        self.Ly = None
        self.latlon = latlon        # flag specifies whether the grid is a latlon grid
        if scalar is not None:
            self.scalar = scalar
        else:
            self.scalar = None
        self.mask = None

        # inherit if given a point set
        if x is not None and y is not None:
            super(GeneralTriangulation, self).__init__(x, y, triangles=tri)

        if self.latlon:
            self.convert_points_xyz()

        return #}}}

    @classmethod
    def from_MPAS_mesh(cls, meshname, varname=None, latlon=False, subdivide=None,
                       loadVars=None): #{{{
        import netCDF4

        ds = netCDF4.Dataset(meshname,'r')

        xv = ds.variables['xCell'][:]
        yv = ds.variables['yCell'][:]
        tri = ds.variables['cellsOnVertex'][:,:]-1
        # take care of edge vertexes (don't include partial triangles there)
        interior = np.where(np.sum(tri < 0, axis=1)  ==  0)[0]
        tri = tri[interior,:]
        if varname:
            varname = ds.variables[varname][interior]
        else:
            varname = None

        if subdivide is not None:
            for asubdiv in np.arange(subdivide):
                xv, yv, tri = subdivide_triangulation(xv, yv, tri)

        tri = cls(xv, yv, tri=tri, scalar=varname, latlon=latlon)

        if ds.is_periodic == "YES":
            tri.periodic = True
            tri.Lx = ds.x_period
            tri.Ly = ds.y_period

        tri.meshname = meshname

        if loadVars:
            for avar in loadVars:
                setattr(tri, avar, ds.variables[avar][:])

        ds.close()

        return tri #}}}

    # operator overloading #{{{

    # define scalar operations between different triangle types with scalars on the grid
    def __neg__(self): #{{{
        self.scalar = -self.scalar
        return self #}}}

    def __rshift__(self, other): #{{{
        new = GeneralTriangulation(other.x, other.y, other.scalar)
        new.scalar, new.mask = self.interp_triang(other, self.scalar, self.interpmethod)
        return new #}}}

    def __irshift__(self, other): #{{{
        other.scalar, other.mask = self.interp_triang(other, self.scalar, self.interpmethod)
        return other #}}}

    def __lshift__(self,other): #{{{
        new = GeneralTriangulation(self.x, self.y, self.scalar)
        new.scalar, new.mask = other.interp_triang(self, other.scalar, self.interpmethod)
        return new #}}}

    def __ilshift__(self,other): #{{{
        self.scalar, self.mask = other.interp_triang(self, other.scalar, self.interpmethod)
        return self #}}}

    # note, a different approach would be to create a superset of the points and interpolate
    # both values on to this triangulation, requiring more computation, however.
    # this is conceptually adding up the linear surfaces to create a new surface blending both
    def gen_op(self, other, op): #{{{
        new = GeneralTriangulation(self.x, self.y, tri=self.triangles, scalar=self.scalar)
        if isinstance(other, GeneralTriangulation):
            scalar, mask = other.interp_triang(new, other.scalar, other.interpmethod)
            # fix the mask
            new.mask = self.mask | mask
        else:
            scalar = other
            new.mask = self.mask
        new.scalar = op(new.scalar, scalar)
        return new #}}}

    def __add__(self, other):  #{{{
        return self.gen_op(other, operator.add) #}}}

    def __sub__(self, other): #{{{
        return self.gen_op(other, operator.sub) #}}}

    def __mul__(self, other): #{{{
        return self.gen_op(other, operator.mul) #}}}

    def __div__(self, other): #{{{
        return self.gen_op(other, operator.div) #}}}

    def __pow__(self, other): #{{{
        return self.gen_op(other, operator.pow) #}}}

    def igen_op(self, other, op): #{{{
        if isinstance(other, GeneralTriangulation):
            scalar, mask = other.interp_triang(self, other.scalar, self.interpmethod)
            # fix the mask
            self.mask = self.mask or other.mask
        else:
            scalar = other
            self.mask = self.mask
        self.scalar = op(self.scalar, scalar)
        return self #}}}

    def __iadd__(self, other):  #{{{
        return self.gen_op(other, operator.add) #}}}

    def __isub__(self, other): #{{{
        return self.gen_op(other, operator.sub) #}}}

    def __imul__(self, other): #{{{
        return self.gen_op(other, operator.mul) #}}}

    def __idiv__(self, other): #{{{
        return self.gen_op(other, operator.div) #}}}

    def __ipow__(self, other): #{{{
        return self.gen_op(other, operator.pow) #}}}
    #}}}

    # methods #{{{
    def update_xy(self, x, y): #{{{

        self.x = x
        self.y = y

        if self.latlon:
            self.convert_points_xyz()

        return #}}}

    def convert_points_xyz(self): #{{{

        # spherical case, y = lat, x = lon
        self.xs, self.ys, self.zs = proj_xyz(self.y, self.x)

        return #}}}

    def compute_area(self, force=False, output=False): #{{{
        if self.area is None or force:
            x = self.x[self.triangles]
            y = self.y[self.triangles]
            if self.latlon:
                # using http://mathworld.wolfram.com/LHuiliersTheorem.html

                arctan = np.arctan
                tan = np.tan

                a = haversine_formula(y[:,0], y[:,1], x[:,0], x[:,1]) / rEarth
                b = haversine_formula(y[:,1], y[:,2], x[:,1], x[:,2]) / rEarth
                c = haversine_formula(y[:,2], y[:,0], x[:,2], x[:,0]) / rEarth
                s = ne.evaluate("0.5*(a+b+c)")

                area = ne.evaluate("4.0*rEarth*rEarth*arctan(sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(.5*(s-b))*tan(.5*(s-c))))")

            else:
                if self.periodic:
                    if self.Lx:
                        x[:,1] = fix_periodicity(x[:,1],x[:,0],self.Lx)
                        x[:,2] = fix_periodicity(x[:,2],x[:,0],self.Lx)
                    if self.Ly:
                        y[:,1] = fix_periodicity(y[:,1],y[:,0],self.Ly)
                        y[:,2] = fix_periodicity(y[:,2],y[:,0],self.Ly)
                #area = 0.5*np.abs((Ax-Cx)*(By-Ay) - (Ax-Bx)*(Cy-Ay))
                area = 0.5*np.abs((x[:,0]-x[:,2])*(y[:,1]-y[:,0]) - (x[:,0]-x[:,1])*(y[:,2]-y[:,0]))

            if output:
              return area
            else:
              self.area = area
        return None #}}}

    def compute_weight_area(self): #{{{
        # need triangle areas
        self.compute_area()
        if self.weightarea is None:
            self.weightarea = np.zeros(self.scalar.shape)
            for atri, area in zip(self.triangles, self.area):
                # proportion 1/3 of each triangles area to
                # each of its precisely three nodes
                self.weightarea[atri] += 1.0/3.0*area
        return #}}}

    def nanmean(self):#{{{
        self.compute_weight_area()
        nanmean = np.nansum(self.scalar*self.weightarea, axis=0)/np.nansum(self.weightarea, axis=0)
        return nanmean #}}}

    def nanstddev(self):#{{{
        self.compute_weight_area()
        nanmean = self.nanmean()
        nanstd = np.sqrt(np.nansum(self.weightarea*(self.scalar - nanmean)**2.0, axis=0)/np.nansum(self.weightarea, axis=0))
        return nanstd #}}}

    def compute_com(self, force=False): #{{{
        # triangles are centered on the minimum vertex to compute the centroid
        if (self.cx is None or force) or (self.cy is None or force):

            if self.latlon:

                # assumes that x,y points have been converted to xs, ys, zs points
                self.cxs = np.mean(self.xs[self.triangles], axis=1)
                self.cys = np.mean(self.ys[self.triangles], axis=1)
                self.czs = np.mean(self.zs[self.triangles], axis=1)

                # y = lat, x = lon
                self.cy, self.cx = proj_lat_long(self.cxs, self.cys, self.czs)

            else:
                x = self.x[self.triangles]
                y = self.y[self.triangles]

                # planar case
                if self.periodic:
                    if self.Lx:
                        xc = np.min(x,axis=1)
                        x[:,0] = fix_periodicity(x[:,0],xc,self.Lx)
                        x[:,1] = fix_periodicity(x[:,1],xc,self.Lx)
                        x[:,2] = fix_periodicity(x[:,2],xc,self.Lx)

                    if self.Ly:
                        yc = np.min(y,axis=1)
                        y[:,0] = fix_periodicity(y[:,0],yc,self.Ly)
                        y[:,1] = fix_periodicity(y[:,1],yc,self.Ly)
                        y[:,2] = fix_periodicity(y[:,2],yc,self.Ly)

                self.cx = np.mean(x, axis=1)
                self.cy = np.mean(y, axis=1)

        return #}}}

    def compute_celledges(self, force=False): #{{{
        # computes the edges neighbors to each cell
        if self.celledges is None or force:
            ncells = self.triangles.shape[0]
            nedges = self.edges.shape[0]
            # initialize
            self.celledges = np.zeros((ncells,3), dtype='int64')
            # sorted edges
            esorted = np.sort(self.edges,axis=1)
            # use hash to do fast lookup
            hashtocell = {i + j*nedges : num for num,(i,j) in enumerate(esorted)}
            # sorted edges on triangles
            edges = np.zeros((nedges,3,2))
            e1s = np.sort(np.vstack((self.triangles[:,0], self.triangles[:,1])).T,axis=1)
            e2s = np.sort(np.vstack((self.triangles[:,1], self.triangles[:,2])).T,axis=1)
            e3s = np.sort(np.vstack((self.triangles[:,2], self.triangles[:,0])).T,axis=1)
            for acell in np.arange(ncells):
                self.celledges[acell,0] = hashtocell[e1s[acell,0] + nedges*e1s[acell,1]]
                self.celledges[acell,1] = hashtocell[e2s[acell,0] + nedges*e2s[acell,1]]
                self.celledges[acell,2] = hashtocell[e3s[acell,0] + nedges*e3s[acell,1]]
        return #}}}

    def compute_edgecellneighs(self, force=False): #{{{
        # computes the cell neighbors to each edge
        if self.edgecellneighs is None or force:
            # initialize to all -1 indicating no neighbor
            self.edgecellneighs = -np.ones(self.edges.shape, dtype='int64')
            self.compute_celledges()

            # expensive naive algorithm
            #for aedge in np.arange(self.edges.shape[0]):
            #    self.edgecellneighs[aedge,:] = np.where(self.celledges==aedge)[0]

            # the algorithm below only loops through the array once, not implicity N*N so it should go as N
            thecell = np.tile(np.arange(self.triangles.shape[0]),(1,3)).ravel()

            front = np.nan*np.zeros(self.edges.shape[0])
            back = np.nan*np.zeros(self.edges.shape[0])

            for acell, aedge in zip(thecell, self.celledges):
                front[aedge] = acell

            for acell, aedge in zip(np.flipud(thecell), np.flipud(self.celledges)):
                back[aedge] = acell

            self.edgecellneighs[:,0] = front
            self.edgecellneighs[:,1] = back

            self.edgecellneighs.sort(axis=-1)

            # clean up boundary edges
            self.edgecellneighs[np.where(self.edgecellneighs[:,0] == self.edgecellneighs[:,1]),1] = -1
        return #}}}

    def compute_celledgedir(self, force=False): #{{{
        if self.celledgedir is None or force:
            self.celledgedir = np.ones_like(self.celledges)
            for acell in np.arange(self.triangles.shape[0]):
                for aedge in np.arange(3):
                    self.celledgedir[acell, aedge] = -(2*(self.edgecellneighs[self.celledges[acell,aedge],0] == acell) - 1 - (self.edgecellneighs[self.celledges[acell,aedge],1] == -1))
        return #}}}

    def compute_nodeedgeneighs(self, force=False): #{{{
        if self.nodeedgeneighs is None or force:
            # computes the edge neighbors to each node
            self.nodeedgeneighs = [[] for nedges in np.unique(self.edges)]
            for ae, aedge in enumerate(self.edges):
                self.nodeedgeneighs[aedge[0]].append(ae)
                self.nodeedgeneighs[aedge[1]].append(ae)

            self.nodeedgeneighs = [np.asarray(edgeneighs, dtype='int') for edgeneighs in self.nodeedgeneighs]
        return #}}}

    def compute_edge_midpoint(self, force=False): #{{{
        if (self.ex is None or self.ey is None) or force:
            if self.latlon:
                self.ey, self.ex = mid_point(self.y[self.edges], self.x[self.edges], ensmean=-1)
            else:
                e1 = self.edges[:,0]
                e2 = self.edges[:,1]
                x1 = self.x[e1]
                x2 = self.x[e2]
                y1 = self.y[e1]
                y2 = self.y[e2]

                if self.periodic:
                    if self.Lx:
                        x2 = fix_periodicity(x2,x1,self.Lx)
                    if self.Ly:
                        y2 = fix_periodicity(y2,y1,self.Ly)
                self.ex = 0.5*(x1+x2)
                self.ey = 0.5*(y1+y2)

        return #}}}

    def compute_le(self, force=False): #{{{
        if self.le is None or force:

            e1 = self.edges[:,0]
            e2 = self.edges[:,1]
            x1 = self.x[e1]
            x2 = self.x[e2]
            y1 = self.y[e1]
            y2 = self.y[e2]

            if self.latlon:

                self.le = haversine_formula(y1, y2, x1, x2)
                # replaced by self.compute_edge_midpoint(), but very expensive unlike for plane
                #self.ey, self.ex = mid_point(self.y[self.edges], self.x[self.edges], ensmean=-1)

            else:

                if self.periodic:
                    if self.Lx:
                        x2 = fix_periodicity(x2,x1,self.Lx)
                    if self.Ly:
                        y2 = fix_periodicity(y2,y1,self.Ly)
                self.le = ne.evaluate('sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))')

        return #}}}

    def return_edge_normal(self): #{{{
        """ N.B. This is done in lat/lon space and is not strictly correct! """
        self.compute_le()
        e1 = self.edges[:,0]
        e2 = self.edges[:,1]
        x1 = self.x[e1]
        x2 = self.x[e2]
        y1 = self.y[e1]
        y2 = self.y[e2]
        if self.periodic:
            if self.Lx:
                x2 = fix_periodicity(x2,x1,self.Lx)
            if self.Ly:
                y2 = fix_periodicity(y2,y1,self.Ly)
        normal = np.array([-(y2-y1), (x2-x1)])
        normal /= np.sqrt(normal[0,:]**2.0 + normal[1,:]**2.0)
        self.compute_com()
        cx = self.cx
        cy = self.cy
        xm = 0.5*(x1 + x2)
        ym = 0.5*(y1 + y2)
        self.compute_edgecellneighs()
        vcx = cx[self.edgecellneighs[:,0]]
        vcy = cy[self.edgecellneighs[:,0]]
        if self.periodic:
            if self.Lx:
                vcx = fix_periodicity(vcx,xm,self.Lx)
            if self.Ly:
                vcy = fix_periodicity(vcy,ym,self.Ly)
        cnorm = np.array([(vcx - xm), (vcy - ym)])
        cnorm /= np.sqrt(cnorm[0,:]**2.0 + cnorm[1,:]**2.0)
        sign = np.sign(np.sum(cnorm * normal, axis=0))
        normal *= sign
        return normal #}}}

    def compute_area_gradient_coefficient(self, force=False): #{{{
        if (self.areagradcoeff is None and self.interioredges is None) or force:

            self.compute_area(force=force)
            self.compute_edgecellneighs()

            c1 = self.edgecellneighs[:,0]
            c2 = self.edgecellneighs[:,1]

            self.interioredges = np.where(c2 != -1)[0]

            A1 = self.area[c1[self.interioredges]]
            A2 = self.area[c2[self.interioredges]]

            self.areagradcoeff = 2.0/3.0 * (A1 + A2)

        return #}}}

    def return_edge_normal_gradient(self, conc, force=False): #{{{
        """ it is the responsibility of the caller to determine if compute_le needs updated """

        # assumes that le and areas have been precomputed!
        dcdn = np.zeros((conc.shape[0],self.edges.shape[0]))
        # the convention is that the positive gradient is
        # always into the first index of self.edgecellneighs

        # would need to update if we allow area to change (more general)
        self.compute_area_gradient_coefficient()
        self.compute_le(force=force)

        ecn = self.edgecellneighs
        ie = self.interioredges
        c1 = conc[:,ecn[ie,0]]
        c2 = conc[:,ecn[ie,1]]

        Ln = (self.areagradcoeff / self.le[self.interioredges])[np.newaxis,:]
        dcdn[:,ie] = ne.evaluate('(c1 - c2)/Ln')

        return dcdn #}}}

    def return_edge_conc(self, conc, cfledge, method='upwind'): #{{{

        cedge = np.zeros((conc.shape[0],self.edges.shape[0]))

        ecn = self.edgecellneighs
        ie = self.interioredges
        c1 = np.zeros_like(cfledge)
        c2 = np.zeros_like(cfledge)
        c1[:, ie] = conc[:,ecn[ie,0]]
        c2[:, ie] = conc[:,ecn[ie,1]]

        # always needed for generalized TVD, e.g., Casulli and Zanolli 2005, eq 12
        if method == 'upwind':
            cedge = ne.evaluate('(cfledge > 0)*c2 + (cfledge < 0)*c1')
        elif method == 'laxwendroff':
            phi = 1.0
            cedge = ne.evaluate('(cfledge > 0)*c2 + (cfledge < 0)*c1 + 0.5*phi*(c2-c1)*(2.0*(cfledge > 0) - 1.0)')
        #elif method == 'vanleer':
        #    r = ne.evaluate('1.0/(c2-c1)*(()*(cfledge < 0))')
        #    phi = ne.evaluate('(r + abs(r))/(1.0 + abs(r))')
        #    cedge = ne.evaluate('(cfledge > 0)*c2 + (cfledge < 0)*c1 + 0.5*phi*(c2-c1)*(2.0*(cfledge > 0) - 1.0)')
        elif method == None:
            raise ValueError, 'Must specify method for use with return_edge_conc!'
        else:
            raise ValueError, 'Method "{}" passed to return_edge_conc is unknown!'.format(method)

        return cedge #}}}
#}}}
    def return_edge_normal_gradient_matrix(self, theta=1.0, kappa=1.0, dt=1.0, force=False): #{{{
        """ it is the responsibility of the caller to determine if compute_le, update_xy need updated.
            matrix multiplication of matrix*conc should give same results as return_edge_normal_gradient"""

        # the convention is that the positive gradient is
        # always into the first index of self.edgecellneighs
        nt = self.triangles.shape[0]
        ne = self.edges.shape[0]
        matrix = sp.dok_matrix((ne, nt))
        cellmatrix = sp.dok_matrix((nt, nt))

        # would need to update if we allow area to change (more general)
        self.compute_area_gradient_coefficient()
        self.compute_le(force=force)

        ecn = self.edgecellneighs
        ie = self.interioredges
        ce = self.celledges
        le = self.le[self.interioredges]
        area = self.area

        Ln = (self.areagradcoeff / le)
        dcdncoeff = 1.0/Ln

        dcdncoefftot = np.zeros((ne))
        dcdncoefftot[self.interioredges] = dcdncoeff

        # loop over cells and build out the matrix
        matrix[ie, ecn[ie,0]] = dcdncoeff
        matrix[ie, ecn[ie,1]] = -dcdncoeff

        for aside in np.arange(3):
            ce = self.celledges[:, aside]
            ed = self.celledgedir[:, aside]

            cellmatrix[np.arange(nt), ecn[ce,0]] += theta*kappa*dt*dcdncoefftot[ce]*self.le[ce]*ed/area
            cellmatrix[np.arange(nt), ecn[ce,1]] -= theta*kappa*dt*dcdncoefftot[ce]*self.le[ce]*ed/area

        #matrix = matrix.tocsc()
        cellmatrix = sp.eye(cellmatrix.shape[0], format='csc') - cellmatrix.tocsc()
        return matrix, cellmatrix #}}}
#{{{
    def return_cell_gradient(self, conc, funcname='MAXnorm', force=False): #{{{
        """ it is the responsibility of the caller to determine if compute_le needs updated """
        # reconstruct from the normals normal gradients
        self.compute_le(force=force)
        self.compute_com(force=force)
        dcdn = self.return_edge_normal_gradient(conc, force=force)
        normal = self.return_edge_normal()

        # store gradient for further analysis
        self.dcdn = dcdn

        # is sharp and produces periodically sharp cell graidents in time
        # (not enough topological smoothing)
        def IED(self, normal, dcdn): #{{{
            self.compute_celledges()
            nx = normal[0,self.celledges]
            ny = normal[1,self.celledges]
            dcdn = dcdn[self.celledges]
            # analytical expressions for reconstructed gradient components
            def return_dcdx(nx, ny, dcdn): #{{{
                nx1 = nx[:,0]
                nx2 = nx[:,1]
                nx3 = nx[:,2]
                ny1 = ny[:,0]
                ny2 = ny[:,1]
                ny3 = ny[:,2]
                dcdn1 = dcdn[:,0]
                dcdn2 = dcdn[:,1]
                dcdn3 = dcdn[:,2]
                # return expression from GradientCalculation ipython notebook at https://gist.github.com/c061665cfbcd7ed0601b
                return ne.evaluate('1./3.*(dcdn1*ny2*(nx1*ny3 - nx3*ny1)*(nx2*ny3 - nx3*ny2) + dcdn1*ny3*(nx1*ny2 - nx2*ny1)*(nx2*ny3 - nx3*ny2) - dcdn2*ny1*(nx1*ny3 - nx3*ny1)*(nx2*ny3 - nx3*ny2) + dcdn2*ny3*(nx1*ny2 - nx2*ny1)*(nx1*ny3 - nx3*ny1) - dcdn3*ny1*(nx1*ny2 - nx2*ny1)*(nx2*ny3 - nx3*ny2) - dcdn3*ny2*(nx1*ny2 - nx2*ny1)*(nx1*ny3 - nx3*ny1))/((nx1*ny2 - nx2*ny1)*(nx1*ny3 - nx3*ny1)*(nx2*ny3 - nx3*ny2))') #}}}
            def return_dcdy(nx, ny, dcdn): #{{{
                nx1 = nx[:,0]
                nx2 = nx[:,1]
                nx3 = nx[:,2]
                ny1 = ny[:,0]
                ny2 = ny[:,1]
                ny3 = ny[:,2]
                dcdn1 = dcdn[:,0]
                dcdn2 = dcdn[:,1]
                dcdn3 = dcdn[:,2]
                # return expression from GradientCalculation ipython notebook at https://gist.github.com/c061665cfbcd7ed0601b
                return ne.evaluate('1./3.*((-dcdn1*nx2 + dcdn2*nx1)*(nx1*ny3 - nx3*ny1)*(nx2*ny3 - nx3*ny2) + (-dcdn1*nx3 + dcdn3*nx1)*(nx1*ny2 - nx2*ny1)*(nx2*ny3 - nx3*ny2) + (-dcdn2*nx3 + dcdn3*nx2)*(nx1*ny2 - nx2*ny1)*(nx1*ny3 - nx3*ny1))/((nx1*ny2 - nx2*ny1)*(nx1*ny3 - nx3*ny1)*(nx2*ny3 - nx3*ny2))') #}}}
            dcdx = return_dcdx(nx,ny,dcdn)
            dcdy = return_dcdy(nx,ny,dcdn)
            return dcdx, dcdy #}}}

        # smoothed IED approach (should have results in between IED and LSQ approaches)
        # still too noisy with not enough topological smoothing!!!
        def avgIED(self, normal, dcdn): #{{{
            def node_reconst(nx1, ny1, dcdn1, nx2, ny2, dcdn2): #{{{
                # return expression from GradientCalculation ipython notebook at https://gist.github.com/c061665cfbcd7ed0601b
                return np.array([(dcdn1*ny2 - dcdn2*ny1)/(nx1*ny2 - nx2*ny1)]),\
                        np.array([(-dcdn1*nx2 + dcdn2*nx1)/(nx1*ny2 - nx2*ny1)]) #}}}

            nx = normal[0,:]
            ny = normal[1,:]

            dcdxnode = np.zeros(len(self.x))
            dcdynode = np.zeros(len(self.x))
            nnodes = np.zeros(len(self.x))

            # compute weighted gradient for nodes
            # reconstruct values at each three vertices
            self.compute_celledges()
            e1 = self.celledges[:,0]
            e2 = self.celledges[:,1]
            e3 = self.celledges[:,2]

            # get nodes associated with edge pairs
            n1 = np.asarray([a[np.where(a[1:] == a[:-1])] \
                    for a in np.sort(np.hstack((self.edges[e1,:], self.edges[e2,:])), axis=1)]).T
            n2 = np.asarray([a[np.where(a[1:] == a[:-1])] \
                    for a in np.sort(np.hstack((self.edges[e2,:], self.edges[e3,:])), axis=1)]).T
            n3 = np.asarray([a[np.where(a[1:] == a[:-1])] \
                    for a in np.sort(np.hstack((self.edges[e3,:], self.edges[e1,:])), axis=1)]).T

            # accumulate computed gradients at nodes
            dcdx, dcdy = node_reconst(nx[e1],ny[e1],dcdn[e1],  nx[e2],ny[e2],dcdn[e2])
            for an, cx, cy in zip(n1, dcdx, dcdy):
                dcdxnode[an] += cx
                dcdynode[an] += cy
                nnodes[an] += 1

            dcdx, dcdy = node_reconst(nx[e2],ny[e2],dcdn[e2],  nx[e3],ny[e3],dcdn[e3])
            for an, cx, cy in zip(n2, dcdx, dcdy):
                dcdxnode[an] += cx
                dcdynode[an] += cy
                nnodes[an] += 1

            dcdx, dcdy = node_reconst(nx[e3],ny[e3],dcdn[e3],  nx[e1],ny[e1],dcdn[e1])
            for an, cx, cy in zip(n3, dcdx, dcdy):
                dcdxnode[an] += cx
                dcdynode[an] += cy
                nnodes[an] += 1

            # finish average at nodes
            dcdxnode /= nnodes
            dcdynode /= nnodes

            # average nodal values to get cell values
            dcdx = np.mean(dcdxnode[self.triangles],axis=1)
            dcdy = np.mean(dcdynode[self.triangles],axis=1)

            return dcdx, dcdy #}}}

        # maximum topological smoothing approach
        # this is the stablest approach and appears to give reasonable results
        def LSQ(self, normal, dcdn): #{{{
            ntracers = dcdn.shape[0]
            ncells = self.triangles.shape[0]
            dcdx = np.zeros((ntracers, ncells))
            dcdy = np.zeros((ntracers, ncells))
            # get the LSQ reconstructed gradient
            nx = normal[0,:]
            ny = normal[1,:]
            self.compute_nodeedgeneighs()

            for atracer in np.arange(ntracers):
                # compute least squares nodal values
                dcdxnode = np.zeros(len(self.nodeedgeneighs))
                dcdynode = np.zeros(len(self.nodeedgeneighs))
                for an, edgeneighs in enumerate(self.nodeedgeneighs):
                    # form the matrix for the inversion
                    nedges = len(edgeneighs)
                    assert nedges > 1, 'Must have two or more edges to do nodal inversion (nedges = %s)!'%(nedges)
                    a = np.zeros((nedges,2))
                    b = np.zeros((nedges))
                    for it, edge in enumerate(edgeneighs):
                        a[it,0] = nx[edge]
                        a[it,1] = ny[edge]
                        b[it] = dcdn[:,edge]

                    # perform least squares fit
                    x,_,_,_ = lstsq(a,b)
                    dcdxnode[an] = x[0]
                    dcdynode[an] = x[1]

                # average nodal values to get cell values
                dcdx[atracer,:] = np.mean(dcdxnode[self.triangles],axis=1)
                dcdy[atracer,:] = np.mean(dcdynode[self.triangles],axis=1)

            return dcdx, dcdy #}}}

        def MAXnorm(self, normal, dcdn, alpha=1.): #{{{
            self.compute_celledges()

            dcdn = alpha*np.max(np.abs(dcdn[:,self.celledges]),axis=2)

            return dcdn, 0.0*dcdn #}}}

        # evaluate using named reconstruction method
        reconst_method = eval(funcname)
        return np.stack((reconst_method(self, normal, dcdn))) #}}}

    def return_cell_gradient_magnitude(self, conc, force=False): #{{{
        """ returns | grad c |^2 """
        dcdn = self.return_edge_normal_gradient(conc, force=force)
        self.compute_area(force=force)
        self.compute_edgecellneighs()

        c1 = self.edgecellneighs[:,0]
        c2 = self.edgecellneighs[:,1]

        # fix ghost cells (mirrors area on boundary)
        ghostcells = np.where(c2 == -1)
        c2[ghostcells] = c1[ghostcells]

        A1 = self.area[c1]
        A2 = self.area[c2]

        edgearea = A1 + A2

        weights = edgearea[self.celledges]
        weights /= np.sum(weights, axis=1)[:,np.newaxis]

        gradc2 = 2.0*np.sum(dcdn[:,self.celledges]**2.0*weights, axis=-1)
        return gradc2 #}}}

    def compute_balls(self, radius): #{{{
        if (self.ballsradius != radius):
            if self.periodic:
                print 'NEED TO UPDATE FOR PERIODICITY!'
            tree = spatial.KDTree(zip(self.x, self.y))
            self.balls = tree.query_ball_point(zip(self.x, self.y), radius)
            self.ballsradius = radius
        return  #}}}

    def smooth_none(self, phi, ntimes, mean, radius): #{{{
        return phi #}}}

    def smooth_running(self, phi, radius=2.0, ntimes=1, mean=np.mean): #{{{
        """ select dimensionality for phi smoothing """
        dim = len(phi.shape)
        if dim == 1:
            phi = self.smooth_scalar_running(phi, radius, ntimes, mean)
        elif dim == 2:
            phi = self.smooth_vector_running(phi, radius, ntimes, mean)
        return phi #}}}

    def smooth_scalar_running(self, phi, radius, ntimes=1, mean=np.mean): #{{{
        self.compute_balls(radius)

        phinew = np.zeros(phi.shape)
        for i in np.arange(ntimes):
            for apoint, aball in enumerate(self.balls):
                phinew[apoint] = mean(phi[aball])
            phi = phinew.copy()

        return phi #}}}

    def smooth_vector_running(self, phi, radius, ntimes=1, mean=np.mean): #{{{
        self.compute_balls(radius)

        phinew = np.zeros(phi.shape)
        for i in np.arange(ntimes):
            for apoint, aball in enumerate(self.balls):
                phinew[apoint,:] = mean(phi[aball,:], axis=0)
            phi = phinew.copy()

        return phi #}}}

    def smooth_laplacian(self, phi, ntimes=1, mean=np.mean): #{{{
        """ select dimensionality for phi smoothing """
        dim = len(phi.shape)
        if dim == 1:
            phi = self.smooth_scalar_laplacian(phi, ntimes, mean)
        elif dim == 2:
            phi = self.smooth_vector_laplacian(phi,ntimes, mean)
        return phi #}}}

    def smooth_vector_laplacian(self, phi, ntimes=1, mean=np.mean): #{{{
        self.compute_area()

        for i in np.arange(ntimes):
            sumarea = np.zeros(self.x.shape)
            sumweights = np.zeros(phi.shape)
            for atri, area in zip(self.triangles, self.area):
                sumarea[atri] += area
                sumweights[atri,:] += mean(phi[atri,:],axis=0)*area
            phi = sumweights/sumarea[:,np.newaxis]

        return phi #}}}

    def smooth_scalar_laplacian(self, phi, ntimes=1, mean=np.mean): #{{{
        if ntimes > 0:
            self.compute_area()

        for i in np.arange(ntimes):
            sumarea = np.zeros(self.x.shape)
            sumweights = np.zeros(self.x.shape)
            for atri, area in zip(self.triangles, self.area):
                sumarea[atri] += area
                sumweights[atri] += mean(phi[atri])*area
            phi = sumweights/sumarea

        return phi #}}}

    def coarsen(self, test=False): #{{{
        """ built using pyamg via rootamg visualization example: https://github.com/pyamg/pyamg/wiki/Examples does not yet account for possible scalar interpolation"""
        self.compute_area()
        # build A
        Np = self.x.shape[0]
        A = sp.lil_matrix((Np,Np))
        for anum in np.arange(Np):
            adj = np.where(np.sum(anum == self.triangles,axis=1))
            for atri in np.sort(adj):
                A[anum,self.triangles[atri,:]] += self.area[atri,np.newaxis]
            A[anum,anum] *= -1
        #A = -A

        # Use Root-Node Solver
        mls = pyamg.rootnode_solver(A.tocsr(), max_levels=2, max_coarse=1, keep=True)

        # Grab the root-nodes (i.e., the C/F splitting)
        Cpts = mls.levels[0].Cpts

        # build coarse triangulation (subsample scalar too but could also perform interpolation
        # or averaging to get the same point.  Simple assignment approach is taken below.)
        if self.scalar is None:
            tri = GeneralTriangulation(self.x[Cpts], self.y[Cpts])
        else:
            tri = GeneralTriangulation(self.x[Cpts], self.y[Cpts], self.scalar[Cpts])

        if test:
            # build edge matrix
            E = np.vstack((A.tocoo().row,A.tocoo().col)).T  # edges of the matrix graph
            # AggOp[i,j] is 1 iff node i belongs to aggregate j
            AggOp = mls.levels[0].AggOp
            inneredges = AggOp.indices[E[:,0]] == AggOp.indices[E[:,1]]
            outeredges = -inneredges

            # return variables
            output = (tri, E, inneredges, outeredges)
        else:
            output = tri

        return output #}}}

    def plot_mesh(self): #{{{
        plt.triplot(self.x, self.y, self.triangles, 'go-')
        return #}}}

    def plot_scalar(self, scalar=None, nfilt=0, cmap='viridis', *args, **kwargs): #{{{
        """ function takes care of nans explicitly so that plot will always occur """
        if scalar is not None:
            scalar = self.smooth_laplacian(scalar,ntimes=nfilt, mean=np.nanmean)
        else:
            scalar = self.scalar
        if self.mask is None:
            mask = np.where(np.isnan(np.mean(scalar[self.triangles],axis=1)),1,0)
            # plot dots in masked triangles
            self.set_mask(1-mask)
            tris = self.get_masked_triangles()
            plt.plot(np.mean(self.x[tris],axis=1), np.mean(self.y[tris],axis=1), 'k.', markersize=0.1)
            self.set_mask(mask)
        plt.tripcolor(self,scalar,cmap=cmap, *args, **kwargs)
        return #}}}

    def get_interpolator(self, scalar, interptype): #{{{
        # hash grid and scalar (http://stackoverflow.com/questions/806151/how-to-hash-a-large-object-dataset-in-python/806342#806342)
        #ahash = hashlib.sha1(np.vstack((self.x, self.y, self.triangles.ravel(), scalar.copy()))).hexdigest()
        ahash = hashlib.sha1(np.vstack(scalar.copy())).hexdigest()
        if self.interphash != ahash:
            self.interphash = ahash
            self.interpolator = interptype(self, scalar)
        return #}}}

    def interp_triang(self, trinew, scalar, interptype=LinearTriInterpolator): #{{{
        """ interpolate scalar onto new triangulation trinew """
        if isinstance(trinew, GeneralTriangulation):
            if self.periodic:
                print 'NEED TO UPDATE FOR PERIODICITY!'
            self.get_interpolator(scalar, interptype)
            interpdata = self.interpolator(trinew.x, trinew.y)
            values = interpdata.data
            mask = np.where(np.sum(np.isnan(values[trinew.triangles]), axis=1),True,False)
        return values, mask #}}}

    def sample_region(self, cx, cy, radius=None, func=None, *kargs, **kwargs): #{{{
        """ a good choice for func is np.mean"""
        if self.periodic:
            print 'NEED TO UPDATE FOR PERIODICITY!'
        # find points within the region
        dist = (self.x - cx)**2.0 + (self.y - cy)**2.0
        if radius:
            points = np.where(dist < radius**2.0)
        else:
            points = np.where(dist == dist.min())

        funcvalue = None
        if self.scalar is not None:
            if len(self.scalar.shape) > 1:
                scalar = np.squeeze(self.scalar[points,:])
            else:
                scalar = np.squeeze(self.scalar[points])
            funcvalue = func(scalar, *kargs, **kwargs)

        return points, funcvalue  #}}}

    def avg_to_vertices(self, triscalar): #{{{
        """
        Averages a value defined on the triangle back to the vertices
        """
        self.compute_area()
        concvertices = np.zeros((self.x.shape[0],))
        areasum = np.zeros((self.x.shape[0],))
        # use np.add.at: http://stackoverflow.com/questions/24099404/numpy-array-iadd-and-repeated-indices/24100418#24100418
        np.add.at(concvertices, self.triangles, triscalar[:,np.newaxis]*self.area[:,np.newaxis])
        np.add.at(areasum, self.triangles, self.area[:,np.newaxis])
        concvertices /= areasum
        return concvertices #}}}

    def max_to_vertices(self, triscalar): #{{{
        """ Max from triangles to vertices """
        return self.extreme_to_vertices(triscalar, np.maximum) #}}}

    def min_to_vertices(self, triscalar): #{{{
        """ Min from triangles to vertices """
        return self.extreme_to_vertices(triscalar, np.minimum) #}}}

    def extreme_to_vertices(self, triscalar, extreme): #{{{
        """
        Max/min of a value defined on the triangle back to the vertices
        """
        concvertices = np.zeros((self.x.shape[0],))
        # use np.add.at: http://stackoverflow.com/questions/24099404/numpy-array-iadd-and-repeated-indices/24100418#24100418
        extreme.at(concvertices, self.triangles, triscalar[:,np.newaxis])
        return concvertices #}}}

    def avg_from_tri_to_vertices(self, vertscalar): #{{{
        """
        Averages a value defined on the vertices to triangles
        (not area-weighted)
        """
        concvertices = np.mean(vertscalar[self.triangles],axis=1)
        return concvertices #}}}
#}}}
    #}}}
#}}}
#{{{
def random_tri(N): #{{{
    x = np.random.random_sample(N)
    y = np.random.random_sample(N)
    values = np.random.randn(N)
    return GeneralTriangulation(x, y, scalar=values) #}}}

def test_shift_operators(N=50): #{{{

    tri1 = random_tri(N)
    tri2 = random_tri(N)
    tri1.scalar = tri1.x + 4*tri1.y
    tri2.scalar = tri2.y - tri2.x

    plt.figure()
    tri1.plot_scalar(shading='gouraud')
    plt.title('phi1')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    tri2.plot_scalar(shading='gouraud')
    plt.title('phi2')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    tri3 = tri2 >> tri1
    tri3.plot_scalar(shading='gouraud')
    plt.title('phi2 >> phi1')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    tri4 = tri2 << tri1
    tri4.plot_scalar(shading='gouraud')
    plt.title('tri2 << tri1')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    tri5 = copy.copy(tri2)
    tri5 <<= tri1
    tri5.plot_scalar(shading='gouraud')
    plt.title('tri2 <<= tri1')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    tri6 = copy.copy(tri2)
    tri6 >>= tri1
    tri6.plot_scalar(shading='gouraud')
    plt.title('tri2 >>= tri1')
    plt.clim([-2,2])
    plt.colorbar()

    return #}}}

def test_arithmetic_operators(N=50): #{{{

    tri1 = random_tri(N)
    tri2 = random_tri(N)
    tri1.scalar = tri1.x + tri1.y
    tri2.scalar = - 1.2*tri2.x + tri2.y
    triadd = copy.deepcopy(tri1)
    triadd += tri2
    triaddexact= copy.deepcopy(tri1)
    triaddexact.scalar = -0.2*triaddexact.x + 2*triaddexact.y

    plt.subplot(2,2,1)
    tri1.plot_scalar(shading='gouraud')
    plt.title('phi1')
    plt.clim([-2,2])
    plt.colorbar()

    plt.subplot(2,2,3)
    tri2.plot_scalar(shading='gouraud')
    plt.title('phi2')
    plt.clim([-2,2])
    plt.colorbar()

    plt.subplot(2,2,2)
    triaddexact.plot_scalar(shading='gouraud')
    plt.title('exact add')
    plt.clim([-2,2])
    plt.colorbar()

    plt.subplot(2,2,4)
    triadd.plot_scalar(shading='gouraud')
    plt.title('add')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()

    plt.subplot(3,1,1)
    triadd.plot_scalar(shading='gouraud')
    plt.title('2 * add')
    plt.clim([-2,2])
    plt.colorbar()

    plt.subplot(3,1,2)
    triadd *= 2.0
    triadd.plot_scalar(shading='gouraud')
    plt.title('add * 2.0')
    plt.clim([-4,4])
    plt.colorbar()

    plt.subplot(3,1,3)
    triadd = triadd / 2.0
    triadd.plot_scalar(shading='gouraud')
    plt.title('add/2.0')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    trimult = random_tri(N)
    trisquared3 = copy.deepcopy(trimult)
    trisquared = trimult * trimult
    trisquared2 = trimult ** 2.0
    trisquared3 *= trimult

    plt.subplot(3,1,1)
    trisquared.plot_scalar(shading='gouraud')
    plt.title('tri * tri')
    plt.clim([-2,2])
    plt.colorbar()
    plt.subplot(3,1,2)
    trisquared2.plot_scalar(shading='gouraud')
    plt.title('tri ** 2.0 ')
    plt.clim([-2,2])
    plt.colorbar()
    plt.subplot(3,1,3)
    trisquared3.plot_scalar(shading='gouraud')
    plt.title('tri *= tri')
    plt.clim([-2,2])
    plt.colorbar()

    return #}}}

def test_GeneralTriangulation(N=50, plotfunc = plt.tripcolor): #{{{

    tri = random_tri(N)
    # plot triangles and their area
    tri.compute_area()
    tri.compute_com()
    plt.hold(True)
    plt.plot(tri.x,tri.y,'.')
    for an, (ax, ay) in enumerate(zip(tri.x,tri.y)):
        plt.text(ax, ay, '%s'%(an))
    for an, aedge in enumerate(tri.edges):
        plt.plot(tri.x[aedge],tri.y[aedge],'k-')
        plt.text(np.mean(tri.x[aedge]), np.mean(tri.y[aedge]),'%d:%s'%(an, aedge))
    for an, atri in enumerate(tri.triangles):
        plt.text(np.mean(tri.x[atri]),np.mean(tri.y[atri]),'%d'%(an))
    #for cx,cy,area in zip(tri.cx, tri.cy, tri.area):
    #    plt.plot(cx,cy,'r.')
    #    plt.text(cx,cy,"%.2f" % (area))
    tri.compute_le()
    normal = tri.return_edge_normal()*tri.le/10.0
    xe = np.mean(tri.x[tri.edges],axis=1)
    ye = np.mean(tri.y[tri.edges],axis=1)
    for aedge in np.arange(tri.edges.shape[0]):
        #plt.plot(xe[aedge], ye[aedge],'r.')
        plt.arrow(xe[aedge], ye[aedge], normal[0,aedge], normal[1,aedge], head_width=0.01)

    plt.axis('equal')

    tri.compute_celledges()
    tri.compute_edgecellneighs()
    print 'edges neighbors to cell:'
    print tri.celledges
    print 'cell neighbors to edge'
    print tri.edgecellneighs

    # compute laplacian
    shtype = 'gouraud'
    shtype = 'flat'
    plt.figure()
    plotfunc(tri, tri.scalar, shading=shtype)
    plt.title('values')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    values_filt = tri.smooth_laplacian(tri.scalar)
    plotfunc(tri, values_filt,shading=shtype)
    plt.title('values smoothed laplacian 1')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    values_filt = tri.smooth_laplacian(values_filt)
    plotfunc(tri, values_filt,shading=shtype)
    plt.title('values smoothed laplacian 2 serial')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    values_filt = tri.smooth_laplacian(values_filt)
    plotfunc(tri, values_filt,shading=shtype)
    plt.title('values smoothed laplacian 3 serial')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    values_filt = tri.smooth_laplacian(tri.scalar, ntimes=3)
    plt.title('values smoothed laplacian 3')
    plotfunc(tri, values_filt,shading=shtype)
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    plotfunc(tri, tri.scalar, shading=shtype)
    plt.title('values')
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    values_filt = tri.smooth_running(tri.scalar, radius=0.5, ntimes=1)
    plt.title('values smoothed running 1')
    plotfunc(tri, values_filt,shading=shtype)
    plt.clim([-2,2])
    plt.colorbar()

    plt.figure()
    values_filt = tri.smooth_running(tri.scalar, radius=0.5, ntimes=2)
    plt.title('values smoothed running 2')
    plotfunc(tri, values_filt,shading=shtype)
    plt.clim([-2,2])
    plt.colorbar()

    return tri #}}}

def test_sample_region(): #{{{

    N = 50
    plotfunc = plt.tripcolor
    tri = random_tri(N)
    # plot triangles and their area
    plt.hold(True)
    plt.axis('equal')
    plt.plot(tri.x,tri.y,'.')

    xc = np.mean(tri.x)
    yc = np.mean(tri.y)
    radius = 0.3

    points, mean = tri.sample_region(xc, yc, radius, np.mean, axis=0)

    circ = plt.Circle((xc,yc),radius, color='k', fill=False)
    plt.gca().add_artist(circ)
    plt.plot(tri.x[points],tri.y[points],'rx', markersize=10)

    plt.figure()
    plt.axis('equal')
    plt.hold(True)
    tri.plot_scalar(shading='flat')

    plt.figure()
    plt.axis('equal')
    tri.scalar[points] = mean
    tri.plot_scalar(shading='flat')
    plt.colorbar()

    tri2 = GeneralTriangulation(np.random.random_sample(N), np.random.random_sample(N))

    points2, _ = tri2.sample_region(np.mean(tri2.x), np.mean(tri2.y),radius, None)
    plt.figure()
    plt.axis('equal')
    plt.hold(True)
    plt.plot(tri.x,tri.y,'.')
    circ = plt.Circle((xc,yc),radius, color='k', fill=False)
    plt.gca().add_artist(circ)
    plt.plot(tri.x[points],tri.y[points],'rx', markersize=10)

    return #}}}

def test_coarsen(x=np.random.rand(300), y=np.random.rand(300)): #{{{
    tri = GeneralTriangulation(x,y)
    # plot triangles and their area

    # coarsening step
    tric, E, inner, outer = tri.coarsen(test=True)

    plt.figure()
    plt.axis('equal')
    plt.hold(True)
    # from pyamg-examples/Rootnode/draw.py
    def lineplot(vertices, indices, linewidths=1):
        """Plot 2D line segments"""
        vertices = np.asarray(vertices)
        indices = np.asarray(indices)

        #3d tensor [segment index][vertex index][x/y value]
        lines = vertices[np.ravel(indices),:].reshape((indices.shape[0],2,2))

        col = mplt.collections.LineCollection(lines)
        col.set_color('k')
        col.set_linewidth(linewidths)

        sub = plt.gca()
        sub.add_collection(col,autolim=True)
        sub.autoscale_view()
    lineplot(np.vstack((tri.x,tri.y)).T, E[inner], linewidths=3.0)
    lineplot(np.vstack((tri.x,tri.y)).T, E[outer], linewidths=0.2)

    plt.scatter(tri.x, tri.y,   c='k', s=100.)
    plt.scatter(tric.x, tric.y, c='r', s=200.)

    return #}}}

def test_from_MPAS_mesh(fname): #{{{
    tri = GeneralTriangulation.from_MPAS_mesh(fname,'indexToVertexID')
    tri.plot_scalar(tri.scalar)
    plt.show()
    return #}}}

def test_subdivision(fname): #{{{

    # test to make sure triangles are subdivided
    tri = GeneralTriangulation.from_MPAS_mesh(fname)
    tris = GeneralTriangulation.from_MPAS_mesh(fname,subdivide=1)
    tris2 = GeneralTriangulation.from_MPAS_mesh(fname,subdivide=2)

    def cleanup():
        plt.axis('equal')
        plt.xlim(0,50000)
        plt.ylim(0,50000)

    plt.subplot(3,1,1)
    tri.plot_mesh()
    cleanup()

    plt.subplot(3,1,2)
    tris.plot_mesh()
    cleanup()

    plt.subplot(3,1,3)
    tris2.plot_mesh()
    cleanup()

    plt.show()

    # test to make sure winding is as expected for tris
    def plot_tri(triT, atri):
        plt.figure()
        plt.plot(triT.x[triT.triangles[atri,0]], triT.y[triT.triangles[atri,0]],'bx')
        plt.plot(triT.x[triT.triangles[atri,1]], triT.y[triT.triangles[atri,1]],'rs')
        plt.plot(triT.x[triT.triangles[atri,2]], triT.y[triT.triangles[atri,2]],'k.')

    plot_tri(tris,0)
    plot_tri(tris,1)

    plt.show()
    return #}}}
#}}}
def test_return_edge_normal_gradient_matrix(fname): #{{{
    tri = GeneralTriangulation.from_MPAS_mesh(fname)
    tri.compute_celledges()
    tri.compute_edgecellneighs()
    tri.compute_celledgedir()

    # compute a random field for the triangles
    conc = np.random.randn(tri.triangles.shape[0])

    dcdn = tri.return_edge_normal_gradient(conc[np.newaxis,:])
    matrix, cellmatrix = tri.return_edge_normal_gradient_matrix()
    dcdnmatrix = matrix.dot(conc)[np.newaxis,:]

    # optimization conversion
    cellmatrix -= sp.eye(cellmatrix.shape[0], format='csc')
    cellmatrix *= -1.0

    # make sure that dcdn - matrix*conc are approximately equal
    assert np.allclose(dcdn, dcdnmatrix), \
            'Matrix multiplication is not equal for return_edge_normal_gradient_matrix calculation'

    result = np.squeeze(np.sum((tri.le[tri.celledges]*dcdn[:,tri.celledges]*tri.celledgedir)/tri.area[np.newaxis, :, np.newaxis], axis=2))
    result2 = cellmatrix.dot(conc)
    assert np.allclose(result, result2), \
            'Matrix multiplication is not equal for full return_edge_normal_gradient_matrix calculation'
    print 'Matrix multiplication test passed'
    return #}}}

def test_return_gradient_magnitude(fname, concfunc): #{{{
    tri = GeneralTriangulation.from_MPAS_mesh(fname)
    tri.compute_celledges()
    tri.compute_edgecellneighs()
    tri.compute_celledgedir()

    # compute a random field for the triangles
    tri.compute_com()
    conc = concfunc(tri)[np.newaxis,:]
    oldgradc2 = np.sum(tri.return_cell_gradient(conc, funcname='LSQ', force=True)**2.0,axis=0)
    newgradc2 = tri.return_cell_gradient_magnitude(conc)

    def plot_scalar(tri, scalar, *args, **kwargs):
        plt.scatter(tri.cx, tri.cy, *args, s=1, c=scalar, edgecolors='face', **kwargs)

    plt.figure()
    plot_scalar(tri, conc[0,:])
    plt.colorbar()

    plt.figure()
    new = newgradc2[0,:]
    old = oldgradc2[0,:]
    scale = np.maximum(new.max(), old.max())
    print scale

    plt.subplot(3,1,1)
    plot_scalar(tri, old)
    plt.clim(0,scale)
    plt.colorbar()

    plt.subplot(3,1,2)
    plot_scalar(tri, new)
    plt.clim(0,scale)
    plt.colorbar()

    plt.subplot(3,1,3)
    plot_scalar(tri, np.log10(new/old), cmap='seismic')
    plt.clim(-1.0,1.0)
    plt.colorbar()

    plt.show()

    return #}}}

if __name__ == "__main__":

    #tri = test_GeneralTriangulation(plotfunc = plt.tricontourf)

    #test_shift_operators()

    #test_arithmetic_operators()

    #test_sample_region()

    #test_coarsen()

    #tri = test_GeneralTriangulation(N=10, plotfunc = plt.tricontourf)
    #plt.show()

    #tri = test_return_edge_normal_gradient_matrix('mesh.nc')

    tri = test_return_gradient_magnitude('mesh.nc', lambda x: (np.sin(2*np.pi*x.cy/1e6) + np.cos(2*np.pi*x.cx/1.e6)))

    #tri = test_from_MPAS_mesh('mesh.nc')

    #test_subdivision('mesh.nc')
