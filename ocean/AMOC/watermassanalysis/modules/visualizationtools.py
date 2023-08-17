#!/usr/bin/env python
"""
    Name: visualizationtools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Visualization tools for the ImPACTS Water Mass Analysis project
"""

import numpy as np
import xarray as xr
import os
import pyremap
from matplotlib import pyplot as plt
from tqdm import tqdm


def parse_contour_levels(clims):
    """Parse `clims` params into arrays for setting contour levels and ticks.
    
    clims = [lower, upper, contour_inc, cbar_inc, residual_bound, contour_inc, cbar_inc]
    """
    
    # Return None if clims is None
    if clims is None:
        return None, None, None, None
    
    # Else parse contour levels and ticks
    else:
        
        # Increment to make upper bound inclusive
        i = 1e-6
    
        # Left panel
        l, u = clims[:2]
        levels1, ticks1 = [np.arange(l, u + i, d) for d in clims[2:4]]

        # Right panel
        u = clims[4]
        levels2, ticks2 = [np.arange(l, u + i, d) for l, d in zip([-u, 0], clims[5:])]
        ticks2 = np.unique(np.hstack([-ticks2, ticks2]))

    return levels1, levels2, ticks1, ticks2


def build_remapper(
    mesh, varnames,
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
    
    # Build remapper object
    meshfile = os.path.join(mesh_path, mesh, meshfiles[mesh])
    mappingfile = os.path.join(mapping_path, f'map_{mesh}_to_0.5x0.5degree_bilinear.nc')
    meshdescriptor = pyremap.MpasMeshDescriptor(meshfile, mesh)
    lonlatdescriptor = pyremap.get_lat_lon_descriptor(dLon=0.5, dLat=0.5)
    remapper = pyremap.Remapper(meshdescriptor, lonlatdescriptor, mappingfile)
    
    # Initialize full domain DataArray to nan
    with xr.open_dataset(meshfile) as ds:
        nan = np.empty(len(ds.nCells)) * np.nan
    da_full = xr.DataArray(nan, dims='nCells')
    
    return remapper, da_full


def build_variables_spatial(ds_in, varnames, timeranges, bbox=None, months=None):
    """Prebuild plotting variables to make plotting faster. Uses `pyremap` for
    remapping to lon lat.
    """
    
    # Define lon lat xarray slicing args if bbox is provided
    if bbox is not None:
        bbox_args = {'lon': slice(*bbox[:2]), 'lat': slice(*bbox[2:])}
    
    # Initialize output dict
    ds_lonlat = {mesh: {} for mesh in ds_in}
    
    # Loop through meshes
    for mesh in ds_in:
        
        # Build remapper objects
        remapper, da_full = build_remapper(mesh, varnames)
        
        # Loop through timeranges
        for timerange in timeranges:
            
            # Initialize lonlat dataset
            tstring = '-'.join(str(t.year) for t in timerange)
            ds_lonlat[mesh][tstring] = xr.Dataset()

            # Extract timerange and seasons
            ds_tslc = ds_in[mesh].sel(time=slice(*timerange))
            if months is not None:
                tindex = [month in months for month in ds_tslc.time.dt.month]
                ds_tslc = ds_tslc.sel(time=tindex)
            
            # Load variables on subdomain into full domain
            cellindex = ds_tslc.nCells.values
            for name in varnames:
                da_full.loc[{'nCells': cellindex}] = ds_tslc[name].mean(dim='time')
                ds_lonlat[mesh][tstring][name] = remapper.remap(da_full)
            
            # Slice according to bbox
            if bbox is not None:
                ds_lonlat[mesh][tstring] = ds_lonlat[mesh][tstring].sel(**bbox_args)
    
    return ds_lonlat


def plot_variable_spatial(ds, varname, units, scale=1, clims=None, cmap=None):
    """Plot specified variable over spatial region. Hard-coded for two timeranges,
    two meshes and the residual between meshes. These categories must be consistent
    with the structure of the `plotvars` dictionary.
    """
    
    # General definitions
    meshes, tstrs = list(ds), list(list(ds.values())[0])[::-1]
    levels1, levels2, ticks1, ticks2 = parse_contour_levels(clims)
    
    # Sigma contour level specs
    sigma_levels = {
        'levels'    : np.arange(25, 28, 0.5),
        'linestyles': ['-', '--', '-', '--', '-', '--'],
        'colors'    : ['k', 'k', 'gray', 'gray', 'lightgray', 'lightgray'],
    }
    
    # Make plot area
    fig, axs = plt.subplots(2, 3, figsize=(12, 8), gridspec_kw={'hspace': 0.05, 'wspace': 0.05})
    
    # Add titles
    for col, title in zip(axs.T, meshes + [f'{meshes[0]}-{meshes[1]}']):
        col[0].set_title(title)
    
    # Loop through time ranges
    for row, tstr in zip(axs, tstrs):
        
        # Print time range
        row[0].text(0.01, 0.8, tstr, transform=row[0].transAxes)
    
        # Loop through meshes
        for ax, mesh in zip(row, meshes):
            
            # Extract variables
            names = ('lon', 'lat', 'sigma', varname)
            lon, lat, sigma, variable = [ds[mesh][tstr][name] for name in names]

            # Plot variable
            c1 = ax.contourf(lon, lat, variable * scale, levels=levels1, cmap=cmap)

            # Plot sigma contours
            cs = ax.contour(lon, lat, sigma, **sigma_levels)
        
        # Plot residual
        residual = np.subtract(*[ds[mesh][tstr][varname] for mesh in meshes])
        c2 = row[2].contourf(lon, lat, residual * scale, levels=levels2, cmap='RdBu_r')
    
    # Remove ticks
    for ax in axs.ravel():
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

    # Add colorbars
    cax1 = fig.add_axes([0.13, 0.07, 0.5, 0.015])
    cax2 = fig.add_axes([0.67, 0.07, 0.21, 0.015])
    fig.colorbar(c1, cax=cax1, label=units, ticks=ticks1, orientation='horizontal')
    fig.colorbar(c2, cax=cax2, label=units, ticks=ticks2, orientation='horizontal')
    
    return fig, axs