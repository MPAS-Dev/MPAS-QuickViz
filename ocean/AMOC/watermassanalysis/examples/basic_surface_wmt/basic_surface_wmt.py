#!/usr/bin/env python
"""
    Name: basic_surface_wmt.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to calculate basic surface WMT
    from MPAS results
"""

import os
import sys
import yaml
from argparse import ArgumentParser
from datetime import datetime

import numpy as np
import pandas as pd
import pyremap
import xarray as xr

# Import modules from path
pkgroot = '/global/cfs/cdirs/m4259/bmoorema/MPAS-QuickViz/ocean/AMOC/watermassanalysis'
sys.path.append(os.path.join(pkgroot, 'modules'))
import watermasstools as wmttools


def define_args():
    """Define arguments for command-line use of this module
    """

    # Construct args object
    parser = ArgumentParser(description='Basic MPAS surface WMT calculator')
    parser.add_argument('meshname', help='Mesh name (LR or HR)')
    parser.add_argument('-c', '--coords', action='store_true', help='Build coordinates file')
    parser.add_argument('-f', '--filenumber', type=int, help='Process filenumber (zero start)')
    parser.add_argument('-p', '--path', type=str, default='./', help='Save path (default ./)')

    return parser.parse_args()


def build_coords_filename(params, meshName=None, fileType='mesh'):
    """Build mesh or mask filename
    """
    
    # Build coordinates filename
    filepath, filename = params[f'{fileType}Path'], params[f'{fileType}File']
    meshNameFull = params['meshName']
    if meshName is not None:
        filename, meshNameFull = filename[meshName], meshNameFull[meshName]
    if fileType == 'mesh':
        filepath = os.path.join(filepath, meshNameFull)
    filename = os.path.join(filepath, filename)
    
    return filename


def build_results_filename(filenumber, params, meshName, startyear=1948):
    """Build results filename and dates from filenumber
    """

    # Get params for building the filename structure
    keys = ['simName', 'prefix', 'resultsPath']
    simName, prefix, resultsPath = [params[key] for key in keys]

    # Apply meshName to params
    if meshName in ['LR', 'HR']:
        meshStr = '18to6v3' if meshName == 'HR' else '60to30E2r2'
        resultsPath = resultsPath[meshName]
        simName = f'{simName}_{meshStr}'
    elif meshName in ['v2', 'v2_1']:
        resultsPath = os.path.join(resultsPath, f'{meshName}Extension')
        simName = f'{meshName}.{simName}'
    
    # Parse dates from filenumber
    year, month = filenumber // 12 + 1, filenumber % 12 + 1
    datestamp = datetime(startyear + year - 1, month, 1)
    datecode = f'{year:04d}-{month:02d}-01'

    # Deal with multiple run directories
    if not isinstance(resultsPath, str):
        for runDir in resultsPath[1]:
            years = [int(year) for year in runDir.split('_')[-2:]]
            if years[0] <= year <= years[1]:
                resultsPath = os.path.join(resultsPath[0], runDir)
                break
        else:
            msg = (
                f'Date {datecode} corresponding to record number {filenumber} '
                f'not available in results record for simulation {simName}'
            )
            raise ValueError(msg)

    # Build filename
    filename = '.'.join([simName, prefix, datecode, 'nc'])
    filename = os.path.join(resultsPath, filename)
    
    return filename, datestamp, datecode


def build_remapper(
    meshFile, subdomain, bbox,
    mapping_path='/global/cfs/cdirs/m4259/mapping_files',
):
    """Build a `pyremap.Remapper` object for the requested mesh. Hardcoded
    to use an existing mapping file to 0.5x0.5degree. Also constructs
    a full domain dataset with the requested varnames initialized to zeros.
    """
    
    # Build descriptor and remapper objects
    meshName = meshFile.split('/')[-2]
    inDescriptor = pyremap.MpasMeshDescriptor(meshFile, meshName)
    outDescriptor = pyremap.get_lat_lon_descriptor(dLon=0.5, dLat=0.5)
    mappingfile = f'{mapping_path}/map_{meshName}_to_{outDescriptor.meshName}_bilinear.nc'
    remapper = pyremap.Remapper(inDescriptor, outDescriptor, mappingfile)
    
    # Build full xarray.DataArray of zeros
    with xr.open_dataset(meshFile) as ds:
        da_full = xr.DataArray(np.zeros(ds.sizes['nCells']), dims='nCells')

    # Define remap variables as dict
    remapvars = {
        'remapper': remapper,
        'bbox': bbox,
        'da_full': da_full,
        'subdomain': subdomain
    }

    return remapvars


def load_coords(meshName=None, savepath='./', bbox=[-100, 40, 0, 85]):
    """Get coordinates and mask variables need for WMT calculations
    as xarray.Dataset
    """
    
    # Load parameters.yaml
    with open('parameters-v2PiControl.yaml') as f:
        params = yaml.safe_load(f)

    # Get meshFile and maskFile names
    meshFile = build_coords_filename(params['run'], meshName)
    maskFile = build_coords_filename(params['run'], meshName, fileType='mask')
    
    coords = []
    
    # Load regionCellMask and reindex by regionNames
    with xr.open_dataset(maskFile) as ds:
        coord = ds.regionCellMasks
        coord.coords['regionNames'] = ds.regionNames.astype(str)
        coord = coord.swap_dims({'nRegions': 'regionNames'})
        coords.append(coord)

    # Load lonCell, latCell, areaCell and convert lonlat to degrees
    with xr.open_dataset(meshFile) as ds:
        for name in ['lonCell', 'latCell', 'areaCell']:
            coord = ds[name]
            if 'lon' in name or 'lat' in name:
                coord = np.rad2deg(coord)
                coord = coord.where(coord <= 180, coord - 360)
            coords.append(coord)

    # Merge coordinate arrays into xr.Dataset
    coords = xr.merge(coords)
    
    # Build and apply subdomain
    coords['nCells'] = coords.nCells
    subdomain = [
        coords.lonCell > bbox[0],
        coords.lonCell < bbox[1],
        coords.latCell > bbox[2],
        coords.latCell < bbox[3],
    ]
    subdomain = np.logical_and.reduce(subdomain)
    coords = coords.isel(nCells=subdomain)
    
    # Convert regionCellMasks to bool
    coords['regionCellMasks'] = coords.regionCellMasks.astype(bool)
    
    return coords
    
    # Save to netCDF
    coords.to_netcdf(os.path.join(savepath, f'{meshName}_coords.nc'))
    
    return coords


def process_monthly_file(filenumber, meshName, coords, savepath='./', bbox=[-100, 40, 0, 85]):
    """Calculate surface WMT for a monthly MPAS results file
    """
    
    # Sigmabin kwargs
    binArgs = (21.9, 28.11, 0.01)
    
    # Load parameters.yaml
    with open('parameters-v2PiControl.yaml') as f:
        params = yaml.safe_load(f)
        
    # Load coords file
    #coords = xr.open_dataset(os.path.join(savepath, f'{meshName}_coords.nc'))
    regionNames, subdomain = [coords[name].values for name in ('regionNames', 'nCells')]
    
    # Get remapvars
    #meshFile = build_coords_filename(params['run'], meshName)
    #remapvars = build_remapper(meshFile, subdomain, bbox)

    # Open results and load fluxes
    filename, datestamp, datecode = build_results_filename(filenumber, params['run'], meshName)
    with xr.open_dataset(filename) as ds:
        fluxes = wmttools.build_fluxes(ds, params['fluxes'], subdomain=subdomain)
    
    # Diagnostic output
    ds = wmttools.calc_wmt(fluxes, coords, 'density', binArgs, regionNames=regionNames)
    return ds

    # Calculate water mass transformation and save output to netCDF
    kwargsDict = {
        '1D': {'regionNames': regionNames},  # 1D kwargs (regional breakdown)
        '2D': {'remapvars': remapvars},      # 2D kwargs (use remapping)
    }
    for ctgy, kwargs in kwargsDict.items():
        ds = wmttools.calc_wmt(fluxes, coords, 'density', binArgs, **kwargs)
        ds = ds.expand_dims({'time': [datestamp]})
        filename = os.path.join(savepath, 'monthly_files', f'{meshName}_WMT{ctgy}_{datecode}.nc')
        ds.to_netcdf(filename, unlimited_dims='time')


if __name__ == "__main__":
    args = define_args()
    if args.coords:
        load_coords(args.meshname, args.path)
    else:
        process_monthly_file(args.filenumber, args.meshname, args.path)