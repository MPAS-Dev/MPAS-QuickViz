#!/usr/bin/env python
"""
    Name: build_mpas_2Dtimeavgs.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to build MPAS 2D time-averaged fields
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import xarray as xr
import pandas as pd
import pyremap
import yaml
import os
from argparse import ArgumentParser
from dateutil.parser import parse as dateparse
from tqdm import tqdm


def define_args():
    """Define arguments for command-line use of this module
    """
    
    # Optional arguments
    args = [
        ('-v', '--varnames'   , 'Variable names to include [var1,var2,...] or none', None),
        ('-t', '--timerange'  , 'Time range to average over [yyyymmdd,yyyymmdd]'   , None),
        ('-b', '--boundingbox', 'Bounding box limits [lon0,lon1,lat0,lat1]'        , None),
        ('-c', '--calctrans'  , 'Calculate 2D water mass transformation'           , 'store_true'),
    ]
    
    # Construct args object
    parser = ArgumentParser(description='Build MPAS 2D time-averaged fields')
    parser.add_argument('filename', help='Path to input file')
    for arg in args:
        parser.add_argument(*arg[:2], help=arg[2], action=arg[3])
    
    return parser.parse_args()


def parse_args(args, ds):
    """Parse `args` and assign defaults using the input `xr.Dataset`
    """
    
    # --- Variable names -----------------------
    if args.varnames is None:
        varnames = [name for name in ds if ds[name].dims == ('time', 'nCells')]
    else:
        varnames = args.varnames.split(',')
        if varnames[0].lower() == 'none':
            varnames = []
    
    # --- Time range ---------------------------
    if args.timerange is None:
        timerange = list(ds.time[[0, -1]].dt.date.values)
    else:
        timerange = [dateparse(date) for date in args.timerange.split(',')]
    
    # --- Bounding box -------------------------
    if args.boundingbox is None:
        bbox = [-100, 20, 0, 80]
    else:
        bbox = [float(lim) for lim in args.boundingbox.split(',')]
    bbox = {'lon': slice(*bbox[:2]), 'lat': slice(*bbox[2:])}
    
    # --- Path strings -------------------------
    path, prefix = os.path.split(args.filename)
    path = os.path.join(os.path.split(path)[0], 'lonlat')
    prefix = prefix.split('.')[0]
    mesh = prefix.split('_')[-1]
    ctgy = 'wmtf' if args.calctrans else 'vars'
    tstr = '_'.join(date.strftime('%Y%m%d') for date in timerange)
    outpath = os.path.join(path, f'{prefix}.mpas2Dtimeavg_{ctgy}_{tstr}.nc')
    
    return varnames, timerange, bbox, mesh, outpath


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
    
    # Build remapper object
    meshfile = os.path.join(mesh_path, mesh, meshfiles[mesh])
    mappingfile = os.path.join(mapping_path, f'map_{mesh}_to_0.5x0.5degree_bilinear.nc')
    meshdescriptor = pyremap.MpasMeshDescriptor(meshfile, mesh)
    lonlatdescriptor = pyremap.get_lat_lon_descriptor(dLon=0.5, dLat=0.5)
    remapper = pyremap.Remapper(meshdescriptor, lonlatdescriptor, mappingfile)
    
    # Initialize full domain DataArray to nan
    with xr.open_dataset(meshfile) as ds:
        nan = np.empty(len(ds.nCells))
    nan[:] = np.nan
    da_full = xr.DataArray(nan, dims='nCells')
    
    return remapper, da_full


def build_combined_variables(ds_in, varsfile='../yaml/variable_definitions.yaml'):
    """
    """
    
    # Open variable definitions
    with open(varsfile, 'r') as f:
        vardefs = yaml.safe_load(f)

    # Combined flux definitions
    ctgys = {
        'radiativeFlux'    : ['shortWaveHeatFlux', 'longWaveHeatFluxUp', 'longWaveHeatFluxDown'],
        'precipitationFlux': ['rainFlux', 'snowFlux'],
        'runoffFlux'       : ['riverRunoffFlux', 'iceRunoffFlux'],
        'totalHeatFlux'    : vardefs['heat'],
        'totalSaltFlux'    : vardefs['salt'],
    }
    
    # Combine variables
    combined = {ctgy: sum([ds_in[name] for name in ctgys[ctgy]]) for ctgy in ctgys}
    combined['buoyancyHeatFlux'] = ds_in['heatFactor'] * combined['totalHeatFlux']
    combined['buoyancySaltFlux'] = ds_in['saltFactor'] * combined['totalSaltFlux']
    combined['buoyancyTotalFlux'] = combined['buoyancyHeatFlux'] + combined['buoyancySaltFlux']
    
    # Output combined variables as xarray dataset
    ds_out = xr.Dataset(combined)
    
    return ds_out


def build_sigma_bins(sigmarange, binsize):
    """Build sigma classes array given `sigmarange` and `binsize`
    """
    
    # Build sigma classes array
    nbins = int(abs(np.subtract(*sigmarange)) / binsize) + 1
    sigmabins = np.arange(nbins) * binsize + sigmarange[0]
    
    return sigmabins


def mask_interval(array, lowerbound, binsize):
    """Return array mask for interval defined by lowerbound and binsize
    """
    
    # Build array mask
    mask = (array >= lowerbound) & (array <= lowerbound + binsize)
    
    return mask


def calc_transformation(
    ds_in, remapper, da_full, cellindex, bbox, sigmarange=[21, 29], binsize=0.1,
):
    """Calculate time-averaged water mass transformation over specified sigma bins,
    remap to lon, lat and return in sigma-lonlat space as an `xr.Dataset`
    """
    
    # Build sigma bins
    sigmabins = build_sigma_bins(sigmarange, binsize)
    sigmaindex = pd.Index(sigmabins, name='sigmaBins')
    
    # Initialize transformation storage lists
    datalists = {ctgy: [] for ctgy in ['heat', 'salt', 'total']}
    datavars = {}

    # Loop through sigmabins
    for sigmabin in sigmabins:

        # Create sigma mask
        mask = mask_interval(ds_in.sigmaTheta, sigmabin, binsize)

        # Loop through flux categories
        for ctgy in datalists:
            
            # Mask sigma and average over time, then assign nan to zero
            name = 'buoyancy' + ctgy.capitalize() + 'Flux'
            values = ds_in[name].where(mask).mean(dim='time')
            values = values.where(~np.isnan(values), 0)
            
            # Populate full DataArray -> remap to lonlat -> append to list
            da_full.loc[cellindex] = values
            da_lonlat = remapper.remap(da_full).sel(**bbox)
            datalists[ctgy].append(da_lonlat)

    # Concatenate DataArrays and calculate final quantities
    for ctgy in datalists:
        transformation = xr.concat(datalists[ctgy], sigmaindex) / binsize * 1e6 # 1e-12 Sv
        formation = -transformation.diff(dim='sigmaBins', label='lower')
        transformation = transformation.isel(sigmaBins=slice(0, -1))
        datavars.update({ctgy + 'Trans': transformation, ctgy + 'Form': formation})
    
    return datavars


def build_mpas_2D(args):
    """Run the build 2D variable fields routine. Uses `pyremap` for
    remapping to lon lat.
    """
    
    # Load aggregated dataset and build combined variables
    ds_in = xr.open_dataset(args.filename)
    ds_in = xr.merge([ds_in, build_combined_variables(ds_in)])
    
    # Parse args
    varnames, timerange, bbox, mesh, outpath = parse_args(args, ds_in)

    # Build remapper objects
    remapper, da_full = build_remapper(mesh)
    cellindex = {'nCells': ds_in.nCells.values}

    # Slice timerange
    ds_in = ds_in.sel(time=slice(*timerange))
    
    # Initialize storage dict
    dataarrays = {name: [] for name in varnames}
    if args.calctrans:
        for ctgy in ['Trans', 'Form']:
            dataarrays.update({name + ctgy: [] for name in ['heat', 'salt', 'total']})
    
    # Loop through months
    months = range(1, 13)
    for month in tqdm(months):
        
        # Isolate month
        ds = ds_in.sel(time=ds_in.time.dt.month == month)

        # Load variables on subdomain into full domain and remap to lonlat
        for name in varnames:
            da_full.loc[cellindex] = ds[name].mean(dim='time')
            dataarrays[name].append(remapper.remap(da_full).sel(**bbox))

        # Calculate spatial water mass transformation
        if args.calctrans:
            wmt = calc_transformation(ds, remapper, da_full, cellindex, bbox)
            for name in wmt:
                dataarrays[name].append(wmt[name])
    
    # Concatenate months and save to netCDF
    monthsindex = pd.Index(months, name='months')
    dataarrays = {name: xr.concat(dataarrays[name], monthsindex) for name in dataarrays}
    xr.Dataset(dataarrays).to_netcdf(outpath)


if __name__ == "__main__":
    build_mpas_2D(define_args())