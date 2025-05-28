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
        dims = ['time', 'nCells']
        varnames = [name for name in ds if all(dim in ds[name].dims for dim in dims)]
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
    ctgy = 'wmt' if args.calctrans else 'vars'
    desc = 'mpas2Dtimeavg_' + ctgy
    savepath, mesh = pptools.build_savepath_from_file(
        args.filename, desc, timerange=timerange, outdir='lonlat',
    )
    
    return varnames, timerange, bbox, mesh, savepath


def build_mpas_2D(args):
    """Run the build 2D variable fields routine. Uses `pyremap` for
    remapping to lon lat.
    """
    
    # Load aggregated dataset and build combined variables
    ds_in = xr.open_dataset(args.filename)
    ds_in = xr.merge([ds_in, pptools.build_combined_variables(ds_in)])
    
    # Parse args
    varnames, timerange, bbox, mesh, savepath = parse_args(args, ds_in)

    # Build remapper objects
    remapvars = pptools.build_remapper(mesh, bbox=bbox)

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
            da = ds[name].mean(dim='time')
            if 'depths' in da.dims:
                da = [pptools.remap(da.sel(depths=depth), **remapvars) for depth in da.depths]
                da = xr.concat(da, 'depths') # **Check that 'depths' is sufficent for new index
            else:
                da = pptools.remap(da, **remapvars)
            dataarrays[name].append(da)

        # Calculate spatial water mass transformation
        if args.calctrans:
            wmtvars = pptools.calc_wmt(ds, remapvars=remapvars)
            for name in wmtvars:
                dataarrays[name].append(wmtvars[name])
    
    # Concatenate months and save to netCDF
    monthsindex = pd.Index(months, name='months')
    dataarrays = {name: xr.concat(dataarrays[name], monthsindex) for name in dataarrays}
    xr.Dataset(dataarrays).to_netcdf(savepath)


if __name__ == "__main__":
    build_mpas_2D(define_args())