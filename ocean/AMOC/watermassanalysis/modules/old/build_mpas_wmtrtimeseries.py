#!/usr/bin/env python
"""
    Name: build_mpas_wmtrtimeseries.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to build MPAS WMTR time series
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import xarray as xr
import pandas as pd
import os
from argparse import ArgumentParser
from dateutil.parser import parse as dateparse
from tqdm import tqdm
import postprocesstools as pptools


def define_args():
    """Define arguments for command-line use of this module
    """
    
    # Optional arguments
    args = [
        ('-v', '--varnames'   , 'Variable names to include [var1,var2,...] or none', None),
        ('-r', '--regions'    , 'Regions to include [region1,region2,...]'         , None),
        ('-c', '--calctrans'  , 'Calculate 2D water mass transformation'           , 'store_true'),
    ]
    
    # Construct args object
    parser = ArgumentParser(description='Build MPAS WMTR timeseries')
    parser.add_argument('filename', help='Path to input file')
    for arg in args:
        parser.add_argument(*arg[:2], help=arg[2], action=arg[3])
    
    return parser.parse_args()


def parse_args(args, ds):
    """Parse `args` and assign defaults using the input `xr.Dataset`
    """
    
    # --- Variable names -----------------------
    if args.varnames is None:
        dims = ['time', 'nCells']
        varnames = [name for name in ds if all(dim in ds[name].dims for dim in dims)]
    else:
        varnames = args.varnames.split(',')
        if varnames[0].lower() == 'none':
            varnames = []
    
    # --- Regions ------------------------------
    if args.regions is None:
        regions = list(ds.regionNames.values.astype(str))
    else:
        regions = args.varnames.split(',')
    
    # --- Path strings -------------------------
    ctgy = 'wmt' if args.calctrans else 'vars'
    desc = 'mpastimeseries_' + ctgy
    savepath, mesh = pptools.build_savepath_from_file(
        args.filename, desc, outdir='timeseries',
    )
    
    return varnames, regions, mesh, savepath


def build_mpas_timeseries(args):
    """Run the build MPAS time series routine.
    """
    
    # Load aggregated dataset and build combined variables
    ds_in = xr.open_dataset(args.filename)
    ds_in = xr.merge([ds_in, pptools.build_combined_variables(ds_in)])
    
    # Parse args
    varnames, regions, mesh, savepath = parse_args(args, ds_in)
    
    # Initialize storage dict
    dataarrays = {name: [] for name in varnames}
    if args.calctrans:
        ctgys, names = ['Trans', 'Form'], ['heat', 'salt', 'total']
        wmtnames = [name + ctgy for ctgy in ctgys for name in names]
        dataarrays.update({name: [] for name in wmtnames})
    
    # Loop through regions
    for region in tqdm(regions, desc='Loading variables by region'):
        
        # Mask area
        regionMask = ds_in.regionMasks.sel(regionNames=bytes(region, 'utf-8'))
        area_region = ds_in.area.where(regionMask)

        # Loop through variables
        for name in varnames:
            
            # Compute region average via area integral
            da = (ds_in[name] * area_region).sum(dim='nCells') / area_region.sum()
            dataarrays[name].append(da)    
    
    # Calculate water mass transformation/formation
    if args.calctrans:
        for date in tqdm(ds_in.time, desc='Loading wmt variables'):
            wmtvars = pptools.calc_wmt(ds_in.sel(time=date))
            for name in wmtvars:
                dataarrays[name].append(wmtvars[name])
    
    # Concatenate along region (variables) and time (wmtr vars)
    variables = {name: xr.concat(dataarrays[name], 'regionNames') for name in varnames}
    variables.update({name: xr.concat(dataarrays[name], 'time') for name in wmtnames})
    xr.Dataset(variables).to_netcdf(savepath)

if __name__ == "__main__":
    build_mpas_timeseries(define_args())