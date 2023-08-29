#!/usr/bin/env python
"""
    Name: aggregate_mpas_2Dvariables.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to aggregate MPAS 2D variables
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import xarray as xr
import os
import sys
import yaml
from datetime import datetime
from fastjmd95 import rho, drhodt, drhods
from tqdm import tqdm


def load_paths_varnames(pathsfile, varsfile='../yaml/variable_definitions.yaml'):
    """Load paths and variable definitions for MPAS-Ocean 2D aggregation
    """
    
    # Load paths dict from yaml
    with open(pathsfile, 'r') as f:
        paths = yaml.safe_load(f)

    # Load variable definitions dict from yaml
    with open(varsfile, 'r') as f:
        vardefs = yaml.safe_load(f)
    
    # Flatten vardefs to get varnames list
    varnames = [name for names in vardefs.values() for name in names]

    return paths, varnames


def build_time_array(paths):
    """Build time array and daterange path defs from path dict names
    """
    
    # Build times and dateranges objects
    times, dateranges = [], {}
    for name in paths:
        years = [int(y) for y in name.split('-')]
        tlist = [datetime(y, m, 1) for y in range(*years) for m in range(1, 13)]
        times.append(np.array(tlist))
        dateranges[name] = times[-1][[0, -1]]
    times = np.hstack(times)
    
    return times, dateranges


def build_MPASO_filename(date, dateranges, paths, startyear=1946):
    """Build an MPAS-Ocean filename for a given date
    """
    
    # Determine which results path to use
    for name in dateranges:
        l, u = dateranges[name]
        if l <= date <= u:
            resultspath = paths['results'][name]
            break
    else:
        raise ValueError(f"Date {date.strftime('%Y-%m-%d')} outside of daterange")
    
    # Build filename
    suffix = f'.{date.year - startyear:04d}-{date.month:02d}-01.nc'
    filename = os.path.join(resultspath, paths['prefix'] + suffix)
    
    return filename


def load_MPASO_mesh(paths):
    """Load MPAS-Ocean mesh variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Define subdomain from maskfile
    with xr.open_dataset(paths['maskfile']) as ds:
        subdomain = ds.regionCellMasks.values[:, 0].astype(bool)

    # Load mesh variables on subdomain
    with xr.open_dataset(paths['meshfile']) as ds:
        coords = {
            'nCells': ds['nCells'].values[subdomain],
            'area': ds['areaCell'].values[subdomain],
            'lon': np.rad2deg(ds['lonCell'].values[subdomain]),
            'lat': np.rad2deg(ds['latCell'].values[subdomain]),
        }

    # Shift lon reference to -180, 180
    index = coords['lon'] > 180
    coords['lon'][index] = coords['lon'][index] - 360

    return coords, subdomain


def load_MPASO_results(resultsfile, varnames, subdomain, prefix='timeMonthly_avg_'):
    """Load MPAS-Ocean results variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Open results file and extract variables
    with xr.open_dataset(resultsfile) as ds:
        
        # Load 2-D variables
        results = {name: ds[prefix + name][0, :].values[subdomain] for name in varnames}
        
        # Load surface tracers
        for name in ['salinity', 'temperature']:
            results[name] = ds[prefix + 'activeTracers_' + name][0, :, 0].values[subdomain]

    return results


def calc_state_variables(S, T, P=0):
    """State Equation from Jackett and McDougal 1995
    Use package fastjmd95 by R. Abernathey and J. Busecke

         `pip install fastjmd95`
    """

    # Define constants
    rho0 = 1026.0     # Seawater density constant [kg m-3]
    cpsw = 3.996e3    # Heat capacity of seawater [J kg-1 K-1]

    # Calculate state variables
    statevars = {
        'sigmaTheta': rho(S, T, P) - 1000,
        'heatFactor': drhodt(S, T, P) / rho0 / cpsw,
        'saltFactor': -drhods(S, T, P) / rho0 * S,
    }

    return statevars


def write_output(variables, coords, paths):
    """Write results to netCDF. Filename conventions are taken from
    the results file prefix and the daterange
    """
    
    # Filename strings
    mesh = paths['meshfile'].split('/')[-2]
    prefix = '_'.join(paths['prefix'].split('_')[:3])
    datestr = '_'.join(time.strftime('%Y%m%d') for time in coords['time'][[0, -1]])
    
    # Build xarray dataset coordinate arrays
    coordinates = {
        'time'  : ('time'  , coords['time']),
        'nCells': ('nCells', coords['nCells']),
    }
    datavars = {
        'lon'   : ('nCells', coords['lon']),
        'lat'   : ('nCells', coords['lat']),
        'area'  : ('nCells', coords['area']),
    }
    datavars.update(
        {name: (['time', 'nCells'], np.vstack(variables[name])) for name in variables}
    )
    
    # Populate xarray dataset and save to netCDF
    filename = os.path.join(paths['out'], f'{prefix}_{mesh}.mpas2Daggregated_{datestr}.nc')
    xr.Dataset(datavars, coordinates).to_netcdf(filename)


def aggregate_mpas_2D(pathsfile):
    """Run the variable aggregation routine for years and paths dict specified.
    `paths` must include `results`, `prefix`, `meshfile` and `maskfile` fields
    """

    # Paths, varnames, mesh variables, time variables
    paths, varnames = load_paths_varnames(pathsfile)
    coords, subdomain = load_MPASO_mesh(paths)
    coords['time'], dateranges = build_time_array(paths['results'])

    # Initialize storage dictionary
    names = ['salinity', 'temperature', 'sigmaTheta', 'heatFactor', 'saltFactor']
    variables = {name: [] for name in names}
    variables.update({name: [] for name in varnames})

    # Loop through filenames
    for date in tqdm(coords['time']):
        
        # Load results
        filename = build_MPASO_filename(date, dateranges, paths)
        data = load_MPASO_results(filename, varnames, subdomain)
        
        # Calculate state variables and append to results
        statevars = calc_state_variables(data['salinity'], data['temperature'])
        data.update(statevars)
        
        # Append results to list
        for name in data:
            variables[name].append(data[name])
    
    # Save output
    write_output(variables, coords, paths)


if __name__ == "__main__":
    aggregate_mpas_2D(sys.argv[1])