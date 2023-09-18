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
import postprocesstools as pptools


def load_paths_vardefs(pathsfile, varsfile='../yaml/variable_definitions.yaml'):
    """Load paths and variable definitions for MPAS-Ocean 2D aggregation
    """
    
    # Load paths dict from yaml
    with open(pathsfile, 'r') as f:
        paths = yaml.safe_load(f)

    # Load variable definitions dict from yaml
    with open(varsfile, 'r') as f:
        vardefs = yaml.safe_load(f)

    return paths, vardefs


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


def get_depthint_params(ds, coords, subdomain):
    """Depth integral parameters for loading MPAS results 3D variables.
    Return indices and thickness referenced to `coords['depths']`.
    `coords['depths']` must start at zero, which is disregarded in the
    returned parameters.
    """
    
    # Define depth integration variables
    kindex = [abs(coords['depth'] - z).argmin() + 1 for z in coords['depths'][1:]]
    thickness = ds.timeMonthly_avg_layerThickness[0, :, :max(kindex)].values[subdomain, :]
    thicktotals = [np.sum(thickness[:, :k], axis=1) for k in kindex]
    
    return kindex, thickness, thicktotals


def load_MPASO_mesh(paths):
    """Load MPAS-Ocean mesh variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Define subdomain from maskfile
    with xr.open_dataset(paths['maskfile']) as ds:
        regionMasks = ds.regionCellMasks
        subdomain = regionMasks.sum(dim='nRegions', dtype=bool).values
        coords = {
            'regionMasks': regionMasks.values[subdomain, :],
            'regionNames': ds.regionNames.values,
        }

    # Load mesh variables on subdomain
    with xr.open_dataset(paths['meshfile']) as ds:
        coords.update({
            'depth': ds.refBottomDepth.values,
            'nCells': ds.nCells.values[subdomain],
            'area': ds.areaCell.values[subdomain],
            'lon': np.rad2deg(ds.lonCell.values[subdomain]),
            'lat': np.rad2deg(ds.latCell.values[subdomain]),
        })

    # Shift lon reference to -180, 180
    index = coords['lon'] > 180
    coords['lon'][index] = coords['lon'][index] - 360

    return coords, subdomain


def load_MPASO_results(
    resultsfile, vardefs, subdomain, coords=None, prefix='timeMonthly_avg_',
):
    """Load MPAS-Ocean results variables as defined by `vardefs`.
    If `coords` is specified, 3D variables will be averaged from
    the surface to the depths specified by `coords['depths']`.
    """
    
    # Open results file and extract variables
    results = {}
    with xr.open_dataset(resultsfile) as ds:
        
        # Depth integration parameters
        if coords is not None:
            kindex, thickness, thicktotals = get_depthint_params(ds, coords, subdomain)
        
        # Load 2-D variables
        for names in vardefs['2D'].values():
            for name in names:
                varname = prefix + name
                results[name] = ds[varname][0, :].values[subdomain]
        
        # Load 3-D variables
        for ctgy, names in vardefs['3D'].items():
            tag = ctgy + '_' if 'activeTracer' in ctgy else ''
            for name in names:
                varname = prefix + tag + name
                
                # Skip if no GM variables present
                if varname not in ds:
                    continue
                
                # Process variable
                else:

                    # Get the depth averages
                    if coords is not None:

                        # Load max depth into memory and then grab the surface field
                        variable = ds[varname][0, :, :max(kindex)].values[subdomain, :]
                        results[name] = [variable[:, 0]]

                        # Loop through depth floors
                        for k, thicktotal in zip(kindex, thicktotals):

                            # Average value via depth integral and append to list
                            varmean = np.sum(variable[:, :k] * thickness[:, :k], axis=1) / thicktotal
                            results[name].append(varmean)

                        # Concatenate depth averages (add zero dim for later concatenation)
                        results[name] = np.vstack(results[name]).T[None, ...]

                    # Just load the surface field
                    else:
                        results[name] = ds[varname][0, :, 0].values[subdomain]

    return results


def calc_state_variables(results, P=0):
    """State Equation from Jackett and McDougal 1995
    Use package fastjmd95 by R. Abernathey and J. Busecke

         `pip install fastjmd95`
    """

    # Define constants
    rho0 = 1026.0     # Seawater density constant [kg m-3]
    cpsw = 3.996e3    # Heat capacity of seawater [J kg-1 K-1]
    
    # Define surface salinity and temperature
    S, T = results['salinity'], results['temperature']
    S, T = [var[:, 0] if var.ndim == 2 else var for var in (S, T)]

    # Calculate state variables
    results.update({
        'sigmaTheta': rho(S, T, P) - 1000,
        'heatFactor': drhodt(S, T, P) / rho0 / cpsw,
        'saltFactor': -drhods(S, T, P) / rho0 * S,
    })

    return results


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
        'time'       : ('time'       , coords['time']),
        'nCells'     : ('nCells'     , coords['nCells']),
        'regionNames': ('regionNames', coords['regionNames']),
    }
    datavars = {
        'lon'        : ( 'nCells', coords['lon']),
        'lat'        : ( 'nCells', coords['lat']),
        'area'       : ( 'nCells', coords['area']),
        'regionMasks': (['nCells', 'regionNames'], coords['regionMasks']),
    }

    # Build xarray dataset variable arrays
    dims = ['time', 'nCells']
    for name in variables:
        variable = np.vstack(variables[name])
        
        # Dims for 3D vars with depth averages
        if variable.ndim == 3:
            datavars[name] = (dims + ['depths'], variable)
            if 'depths' not in coordinates:
                coordinates['depths'] = ('depths', coords['depths'])
        
        # Surface slices only
        else:
            datavars[name] = (dims, variable)
    
    # Populate xarray dataset and save to netCDF
    filename = os.path.join(paths['out'], f'{prefix}_{mesh}.mpas2Daggregated_{datestr}.nc')
    xr.Dataset(datavars, coordinates).to_netcdf(filename)


def aggregate_mpas_2D(pathsfile):
    """Run the variable aggregation routine for years and paths dict specified.
    `paths` must include `results`, `prefix`, `meshfile` and `maskfile` fields
    """

    # Paths, varnames, mesh variables, time variables
    paths, vardefs = load_paths_vardefs(pathsfile)
    coords, subdomain = load_MPASO_mesh(paths)
    coords['time'], dateranges = build_time_array(paths['results'])
    coords['depths'] = [0, 100, 500] # Averging depths (0 required)

    # Loop through filenames
    variables, ntime = {}, len(coords['time'])
    for t, date in enumerate(coords['time']):
        
        # Load results and calculate state variables
        filename = build_MPASO_filename(date, dateranges, paths)
        results = load_MPASO_results(filename, vardefs, subdomain, coords=coords)
        results = calc_state_variables(results)
        
        # Append results to list
        for name in results:
            if name not in variables:
                variables[name] = []
            variables[name].append(results[name])
        
        # Print status
        pptools.loopstatus(t, ntime)
    
    # Save output
    write_output(variables, coords, paths)


if __name__ == "__main__":
    aggregate_mpas_2D(sys.argv[1])