#!/usr/bin/env python
"""
    Name: build_mpas_wmtrtimeseries.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to build MPAS water mass
    transformation and formation timeseries
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
        ('-v', '--varnames' , 'Variable names to include [var1,var2,...] or none', None),
        ('-t', '--timerange', 'Time range to average over [yyyymmdd,yyyymmdd]'   , None),
    ]
    
    # Construct args object
    parser = ArgumentParser(description='Build MPAS wmtr timeseries')
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
    
    # --- Path strings -------------------------
    path, prefix = os.path.split(args.filename)
    path = os.path.join(os.path.split(path)[0], 'wmtr')
    prefix = prefix.split('.')[0]
    mesh = prefix.split('_')[-1]
    tstr = '_'.join(date.strftime('%Y%m%d') for date in timerange)
    outpath = os.path.join(path, f'{prefix}.mpastimeseries_wmtr_{tstr}.nc')
    
    return varnames, timerange, mesh, outpath


def build_mpas_wmtr_timeseries(args):
    """Run the build 2D variable fields routine. Uses `pyremap` for
    remapping to lon lat.
    """
    
    # Load aggregated dataset and build combined variables
    ds_in = xr.open_dataset(args.filename)
    ds_in = xr.merge([ds_in, pptools.build_combined_variables(ds_in)])
    
    # Parse args
    varnames, timerange, mesh, outpath = parse_args(args, ds_in)
    
    # Initialize storage dict
    dataarrays = {name: [] for name in varnames}
    for ctgy in ['Trans', 'Form']:
        dataarrays.update({name + ctgy: [] for name in ['heat', 'salt', 'total']})
    
    # Loop through time
    for date in tqdm(ds_in.time):
        
        # Isolate time
        ds = ds_in.sel(time=date)

        # Load variables on subdomain into full domain and remap to lonlat
        for name in varnames:
            da = ds[name]
            dataarrays[name].append(da)

        # Calculate spatial water mass transformation
        if args.calctrans:
            wmtvars = pptools.calc_wmt(ds, remapvars=remapvars, bbox)
            for name in wmtvars:
                dataarrays[name].append(wmtvars[name])
    
    # Concatenate months and save to netCDF
    monthsindex = pd.Index(months, name='months')
    dataarrays = {name: xr.concat(dataarrays[name], monthsindex) for name in dataarrays}
    xr.Dataset(dataarrays).to_netcdf(outpath)


if __name__ == "__main__":
    build_mpas_2D(define_args())