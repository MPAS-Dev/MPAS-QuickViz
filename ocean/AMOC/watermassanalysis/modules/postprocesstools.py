#!/usr/bin/env python
"""
    Name: postprocesstools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for model results postprocessing
    for the ImPACTS Water Mass Analysis project.
"""

import os
import time
from copy import deepcopy

import numpy as np
import pyremap
import xarray as xr
from scipy import signal


def loopstatus(k, n, starttime, interval=1):
    """Print loop progress at percentage intervals given an iteration k
    and a total number of iterations n
    """
    
    nscale = 100 / n
    percent = k * nscale
    if percent % interval < nscale:
        time_elapsed = time.time() - starttime
        msg = f'{int(percent)}% complete ... {time_elapsed:.0f} seconds elapsed'
        print(msg, flush=True)


def lowpass(data, cutoff=10, window_type='boxcar'):
    """Apply a Finite Impulse Response (FIR) lowpass filter according
    to the window_type and cutoff using a convolution algorithm
    """
    
    # Calculate low-pass signal
    window = signal.get_window(window_type, cutoff)
    filtered = np.convolve(data, window / sum(window), mode='same')
    
    return filtered


def downsample(array, widths=(5, 5)):
    """Downsample array using a groupby mean every w elements.
    Pads the array dimensions with nan so that `numpy.reshape` can be used
    in the groupby operation.
    """

    # Resample array
    pad = [(0, w-dim%w) for dim, w in zip(array.shape, widths)]
    array = np.pad(array, pad, constant_values=np.nan)
    args = [arg for dim, w in zip(array.shape, widths) for arg in (int(dim/w), w)]
    axis = tuple(i for i, e in enumerate(args) if e in widths)
    array = np.nanmean(array.reshape(*args), axis=axis)
    
    return array


def rotate_velocities(u, v, angle):
    """Rotate velocities
    """
    
    # Rotate velocity
    cosa, sina = np.cos(angle), np.sin(angle)
    u_rot = u * cosa - v * sina
    v_rot = u * sina + v * cosa
    
    return u_rot, v_rot


def interpolate_to_edge(varCell, cellsOnEdge, subdomain):
    """Interpolate cell variable to edge.
    """
    
    # Find cells on edge in subdomain
    index = list(np.array([np.where(np.isin(subdomain, pair))[0] for pair in cellsOnEdge]).T)

    # Interpolate to edge
    varEdge = (varCell[index[0], :] + varCell[index[1], :]) / 2
    
    return varEdge


def interpolate_along_axis(array, axis, direction):
    """Interpolate array along axis. Direction `up` interpolates
    T points toward U points. Direction `down` interpolates U
    points toward T points. For POP meshes.
    """
    
    # Get indexing based on interpolation direction
    n = array.shape[axis]
    if direction == 'up':
        idx, order = n - 1, -1
    else:
        idx, order = 0, 1
        
    # Expand array across the boundary of the interpolation dimension
    boundary = np.expand_dims(array.take(idx, axis), axis)
    array = [boundary, array][slice(None, None, order)]
    array = np.concatenate(array, axis)
    
    # Interpolate linearly
    array_interp = []
    for start in [0, 1]:
        array_interp.append(array.take(range(start, start + n), axis))
    array_interp = sum(array_interp) / 2
    
    return array_interp


def build_remapper(
    meshfile, grid=None, subdomain=None, bbox=None,
    mapping_path='/global/cfs/cdirs/m4259/mapping_files/',
):
    """Build a `pyremap.Remapper` object for the requested mesh. Hardcoded
    to use an existing mapping file to 0.5x0.5degree. Also constructs
    a full domain dataset with the requested varnames initialized to zeros.
    """
    
    # Build model-specific inDescriptors
    if 'mpas' in meshfile:
        meshName = meshfile.split('/')[-2]
        inDescriptor = pyremap.MpasMeshDescriptor(meshfile, meshName)
    elif 'pop' in meshfile:
        meshName = 'gx1v6' + grid
        inDescriptor = pyremap.LatLon2DGridDescriptor.read(
            meshfile,
            latVarName=grid + 'LAT',
            lonVarName=grid + 'LONG',
            regional=False,
        )
    else:
        raise ValueError('No string mpas or pop found in: {meshfile}')

    # Build outDescriptor and define mappingfile path
    outDescriptor = pyremap.get_lat_lon_descriptor(dLon=0.5, dLat=0.5)
    mappingfile = f'map_{meshName}_to_{outDescriptor.meshName}_bilinear.nc'
    mappingfile = os.path.join(mapping_path, mappingfile)
    remapper = pyremap.Remapper(inDescriptor, outDescriptor, mappingfile)

    # Define remap variables as dict
    remapvars = {'remapper': remapper}
    
    # Add bbox
    if bbox is not None:
        remapvars['bbox'] = bbox
    
    # Add subdomain
    if 'mpas' in meshfile and subdomain is not None:
        with xr.open_dataset(meshfile) as ds:
            da_full = xr.DataArray(np.zeros(ds.sizes['nCells']), dims='nCells')
        remapvars.update({'da_full': da_full, 'subdomain': subdomain})

    return remapvars


def remap(array_in, remapper, bbox=None, da_full=None, subdomain=None):
    """Remap `da_in` to lonlat using `remapper` object. `da_in` is on a
    subdomain so it must be populated into `da_full` before remapping
    """
    
    # Build input as xarray.Dataarray
    if subdomain is not None:
        da_in = deepcopy(da_full)
        da_in.loc[subdomain] = array_in
    else:
        da_in = xr.DataArray(array_in, dims=['nlat', 'nlon'])
    
    # Remap to lonlat
    da_out = remapper.remap(da_in)
    if bbox is not None:
        da_out = da_out.sel(lon=slice(*bbox[:2]), lat=slice(*bbox[2:]))
    
    return da_out