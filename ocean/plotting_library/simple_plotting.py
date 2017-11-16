#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from latlon_coordinate_transforms import fix_periodicity as fix_periodicity

def rad2deg(rad):
    return 180.0/np.pi*rad

def deg2rad(deg):
    return np.pi/180.0*deg

def convert_lonstr(lonstr):
    lon = float(lonstr[:-1])
    W = lonstr[-1] == 'W'
    if W:
        lon = 360.0 - lon
    return deg2rad(lon)

def convert_latstr(latstr):
    lat = float(latstr[:-1])
    S = latstr[-1] == 'S'
    if S:
        lat *= -1.0
    return deg2rad(lat)

def plot_poly(colors, cells, nvertices, x, y, xc=None, yc=None,
        xperiod=None, yperiod=None, colorbar=True, cmap='viridis'):
    """ from https://gist.github.com/pwolfram/2745eaccaf33f222fed6"""
    patches = []
    for ii,(ac,nv) in enumerate(zip(cells, nvertices)):
        ac = np.asarray(ac[:nv], dtype='int')
        xp = x[ac]
        yp = y[ac]
        if xc is not None:
            xp = fix_periodicity(xp, xc[ii], xperiod)
        if yc is not None:
            yp = fix_periodicity(yp, yc[ii], yperiod)
        patches.append(Polygon(zip(xp,yp)))

    pc = PatchCollection(patches, cmap=cmap, alpha=1.0)
    pc.set_array(colors)
    pc.set_edgecolor('face')
    pc.set_lw(0.1)

    plt.figure()
    ax = plt.gca()
    ax.add_collection(pc)
    ax.set_aspect('equal')
    ax.autoscale_view()
    if colorbar:
        plt.colorbar(pc)
    return pc

def plot_var(ds, variable, time=-1, maxdepth=-500.):
    zcoord = ds.refZMid.values
    nmax = np.minimum(ds.maxLevelCell.values.max(), np.where(zcoord < maxdepth)[0][0])
    if 'nCells' in variable.dims:
        for ii in np.arange(len(variable.nCells)):
            lon = rad2deg(ds.lonCell.values[ii])
            lat = rad2deg(ds.latCell.values[ii])
            plt.plot(variable[time,ii,:nmax],zcoord[:nmax], label='lon=%.4f, lat=%.4f'%(lon,lat))
    if 'nEdges' in variable.dims:
        for ii in np.arange(len(variable.nEdges)):
            lon = rad2deg(ds.lonEdge.values[ii])
            lat = rad2deg(ds.latEdge.values[ii])
            plt.plot(variable[time,ii,:nmax],zcoord[:nmax], label='lon=%.4f, lat=%.4f'%(lon,lat))
    if 'nVertices' in variable.dims:
        for ii in np.arange(len(variable.nVertices)):
            lon = rad2deg(ds.lonVertex.values[ii])
            lat = rad2deg(ds.latVertex.values[ii])
            plt.plot(variable[time,ii,:nmax],zcoord[:nmax], label='lon=%.4f, lat=%.4f'%(lon,lat))
    plt.ylabel('Depth (m)')
    plt.legend(loc='best',ncol=2, fontsize=12)
    attrs = variable.attrs
    #plt.xlabel(attrs['long_name'] + ' ' + attrs['units'])
    plt.xlabel(variable.name + ' t=%d'%(time))

def plot_horiz(ds, variable, atime=-1, maxLayers=50, layerDepth=None, lonlat=False, periodic=False):
    if layerDepth is not None:
        # get first layer underneath the layerDepth
      kkrange = [np.where(ds.refZMid < layerDepth)[0][0]]
    else:
        kkrange = np.arange(np.minimum(maxLayers,variable.shape[-1]))

    if 'nCells' in variable.dims:
        # plot cell values
        cells = ds.verticesOnCell.values - 1
        nvert = ds.nEdgesOnCell.values
        vertextype = 'Vertex'
    elif 'nVertices' in variable.dims:
        cells = ds.cellsOnVertex.values - 1
        nvert = np.asarray(3*np.ones((ds.cellsOnVertex.shape[0])), dtype='i')
        vertextype = 'Cell'
    elif 'nEdges' in variable.dims:
        def make_plot(var, kk=0):
            plt.figure()
            sc = plt.scatter(ds.xEdge, ds.yEdge,c=var.values, cmap='RdYlBu')
            #plt.title(var.attrs['long_name'] + ' k=%d (%f m)'%(kk, np.round(ds.refZMid[kk])))
            plt.title(var.name + ' t=%d k=%d (%.1f m)'%(atime, kk, np.round(ds.refZMid[kk])))
            plt.gca().set_facecolor('#939393')
            return sc
        if len(variable.shape) == 2:
            var = variable[atime,:]
            sc = make_plot(var)
        else:
            for kk in kkrange:
                var = variable[atime,:,kk]
                sc = make_plot(var,kk)
        return sc, None

    for kk in kkrange:
        plt.figure()
        if lonlat:
            x = rad2deg(ds['lon' + vertextype].values)
            y = rad2deg(ds['lat' + vertextype].values)
        else:
            x = ds['x' + vertextype].values
            y = ds['y' + vertextype].values

        if periodic:
            if 'nCells' in variable.dims:
                xc = ds.xCell.values
                yc = ds.yCell.values
            if 'nVertices' in variable.dims:
                xc = ds.xVertex.values
                yc = ds.yVertex.values
            xperiod = ds.x_period
            yperiod = ds.y_period
        else:
            xc = None
            yc = None
            xperiod = None
            yperiod = None
        if len(variable.shape) == 2:
            var = variable[atime,:]
        else:
            if 'Time' in variable.dims:
                var = variable[atime,:,kk]
            else:
                var = variable[:]
        pc = plot_poly(var.values, cells, nvert, x, y, xc, yc, xperiod, yperiod, 
                       colorbar=False)
        #plt.title(var.attrs['long_name'] + ' k=%d (%f m)'%(kk, np.round(ds.refZMid[kk])))
        try:
            plt.title(var.name + ' t=%d k=%d (%.1f m)'%(atime, kk, np.round(ds.refZMid[kk])))
        except:
            pass

        return pc, var.values

def test_plot_poly():
    """ see https://gist.github.com/pwolfram/2745eaccaf33f222fed6"""
    pass


if __name__ == "__main__":
    test_plot_poly()
