#!/usr/bin/env python
"""
    Name: postprocesstools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for model results postprocessing
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import xarray as xr
import pandas as pd
import yaml
import pyremap


def build_combined_variables(ds_in, varsfile='../yaml/variable_definitions.yaml'):
    """
    """
    
    # Open variable definitions
    with open(varsfile, 'r') as f:
        ctgys = yaml.safe_load(f)
    
    # Combine variables
    combined = {ctgy: sum([ds_in[name] for name in ctgys[ctgy]]) for ctgy in ctgys}
    combined['buoyancyHeatFlux'] = ds_in['heatFactor'] * combined['totalHeatFlux']
    combined['buoyancySaltFlux'] = ds_in['saltFactor'] * combined['totalSaltFlux']
    combined['buoyancyTotalFlux'] = combined['buoyancyHeatFlux'] + combined['buoyancySaltFlux']
    
    # Output combined variables as xarray dataset
    ds_out = xr.Dataset(combined)
    
    return ds_out


def build_remapper(
    mesh,
    mesh_path='/global/cfs/cdirs/e3sm/inputdata/ocn/mpas-o/',
    mapping_path='/global/cfs/cdirs/m4259/mapping_files/',
):
    """Build a `pyremap.Remapper` object for the requested mesh. Hardcoded
    to use an existing mapping file to 0.5x0.5degree. Also constructs
    a full domain dataset with the requested varnames initialized to nan.
    """
    
    # Mesh file name dictionary
    meshfiles = {
        'EC30to60E2r2': 'ocean.EC30to60E2r2.210210.nc',
        'oRRS18to6v3': 'oRRS18to6v3.171116.nc',
    }
    
    # Build remapper arguments
    meshfile = os.path.join(mesh_path, mesh, meshfiles[mesh])
    meshdescriptor = pyremap.MpasMeshDescriptor(meshfile, mesh)
    lonlatdescriptor = pyremap.get_lat_lon_descriptor(dLon=0.5, dLat=0.5)
    mappingfile = os.path.join(mapping_path, f'map_{mesh}_to_0.5x0.5degree_bilinear.nc')
    
    # Build nan array with len nCells
    with xr.open_dataset(meshfile) as ds:
        nan = np.empty(len(ds.nCells))
    nan[:] = np.nan
    
    # Define remap variables as dict
    remapvars = {
        'da_full': xr.DataArray(nan, dims='nCells'),
        'remapper': pyremap.Remapper(meshdescriptor, lonlatdescriptor, mappingfile),
    }
    
    return remapvars


def remap(da_in, da_full, remapper, bbox=[-100, 20, 0, 80]):
    """Remap `da_in` to lonlat using `remapper` object. `da_in` is on a
    subdomain so it must be populated into `da_full` before remapping
    """
    
    # Remap to lonlat
    da_full.loc[da_in.nCells] = da_in
    da_out = remapper.remap(da_full).sel(**bbox)
    
    return da_out


def build_sigma_bins(sigmarange, binsize):
    """Build sigma classes array given `sigmarange` and `binsize`
    Return as `pandas.Index` object
    """
    
    # Build sigma classes array
    nbins = int(abs(np.subtract(*sigmarange)) / binsize) + 1
    sigmabins = np.arange(nbins) * binsize + sigmarange[0]
    sigmabins = pd.Index(sigmabins, name='sigmaBins')
    
    return sigmabins


def calc_wmt(ds_in, remapvars=None, sigmarange=[21, 29], binsize=0.1):
    """Calculate time-averaged water mass transformation over specified sigma bins,
    remap to lon, lat and return in sigma-lonlat space as an `xr.Dataset`
    """
    
    # Build sigma variables
    sigma = ds_in.sigmaTheta
    sigmabins = build_sigma_bins(sigmarange, binsize)
    
    # Define unit scale (1D: Sv, 2D: 1e-6 Sv km-2)
    scale = 1e-6 if remapvars is None else 1e6
    
    # Initialize data storage
    data = {ctgy: [] for ctgy in ['heat', 'salt', 'total']}

    # Loop through sigmabins
    for sigmabin in sigmabins:

        # Create sigma mask
        mask = (sigma >= sigmabin) & (sigma <= sigmabin + binsize)

        # Loop through flux categories
        for ctgy in data:
            
            # Mask sigma and average over time, then assign nan to zero
            da = ds_in[f'buoyancy{ctgy.capitalize()}Flux'].where(mask)
            
            if remapvars is None:  # 1D output
                
                # Integrate over outcrop
                da = (da * ds_in.area).sum(dim='nCells')

            else:  # 2D output
                
                # Time average -> replace nan with zero -> remap to lonlat
                da = da.mean(dim='time').where(~np.isnan(da_in), 0)
                da = remap(da, **remapvars, bbox=bbox)
            
            # Append to list
            data[ctgy].append(da)

    # Concatenate DataArrays and calculate final quantities
    variables = {}
    for ctgy in datalists:
        transformation = xr.concat(data[ctgy], sigmabins) / binsize * scale
        variables[ctgy + 'Trans'] = transformation.isel(sigmaBins=slice(0, -1))
        variables[ctgy + 'Form'] = -transformation.diff(dim='sigmaBins', label='lower')
    
    return variables