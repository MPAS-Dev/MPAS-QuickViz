#!/usr/bin/env python
"""
    Name: remapping.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for model results postprocessing
    for the ImPACTS Water Mass Analysis project.
"""

import os
from copy import deepcopy

import numpy as np
import pyremap
import xarray as xr


def build_remapper(
    params, resolution=0.5, grid=None, subdomain=None,
    mapping_path='/global/cfs/cdirs/m4259/mapping_files/',
):
    """Build a `pyremap.Remapper` object for the requested mesh. Hardcoded
    to use an existing mapping file to 0.5x0.5degree. Also constructs
    a full domain dataset with the requested varnames initialized to zeros.
    """

    # Unpack params
    names = ['model', 'meshName', 'meshFile']
    model, meshName, meshFile = [params[name] for name in names]
    
    # Build model-specific input descriptors
    if model == 'MPAS':
        inDescriptor = pyremap.MpasMeshDescriptor(meshFile, meshName)
    elif model == 'POP':
        meshName = meshName + grid
        inDescriptor = pyremap.LatLon2DGridDescriptor.read(
            meshFile,
            latVarName=grid + 'LAT',
            lonVarName=grid + 'LONG',
            regional=False,
        )
    else:
        raise ValueError(f"Invalid model identifier: {params['model']}")

    # Build output descriptor
    outDescriptor = pyremap.get_lat_lon_descriptor(
        dLon=resolution,
        dLat=resolution,
    )

    # Define mapping file path
    mappingfile = f'map_{meshName}_to_{outDescriptor.meshName}_bilinear.nc'
    mappingfile = os.path.join(mapping_path, mappingfile)

    # Build remapper object
    remapper = pyremap.Remapper(inDescriptor, outDescriptor, mappingfile)

    # Build mapping file
    if not os.path.exists(mappingfile):
        remapper.esmf_build_map(method='bilinear', mpi_tasks=1)

    # Define remap variables as dict
    remapvars = {'remapper': remapper}
    
    # Add bbox
    if 'bbox' in params:
        remapvars['bbox'] = bbox
    
    # Add subdomain
    if model == 'MPAS' and subdomain is not None:
        ds = xr.open_dataset(meshFile)
        da_full = xr.DataArray(np.zeros(ds.sizes['nCells']), dims='nCells')
        remapvars['da_full'] = da_full
        remapvars['subdomain'] = subdomain

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
