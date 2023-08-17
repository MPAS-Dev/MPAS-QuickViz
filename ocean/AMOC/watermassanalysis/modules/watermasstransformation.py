#!/usr/bin/env python
"""
    Name: watermasstransformation.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to perform water mass transformation calculations
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


def load_variable_definitions(
    vardefsfile='../yaml/variable_definitions.yaml',
):
    """Variable definitions for MPAS-Ocean
    """

    # Open yaml file and load fluxdefs dict
    with open(vardefsfile, 'r') as f:
        vardefs = yaml.safe_load(f)
    
    # Flatten vardefs to get unique varnames list
    varnames = [name for names in vardefs.values() for name in names]
    varnames = list(np.unique(np.hstack(varnames)))
    vardefs.pop('other')

    return varnames, vardefs


def load_MPASO_mesh(paths):
    """Load MPAS-Ocean mesh variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Load mask variables
    with xr.open_dataset(paths['maskfile']) as ds:
        coords = {
            'regionNames': ds['regionNames'].values.astype(str),
            'regionMasks': ds['regionCellMasks'].values.astype(bool),
        }
    
    # All points to be considered (for memory purposes)
    subdomain = np.sum(coords['regionMasks'], axis=1).astype(bool)
    coords['regionMasks'] = coords['regionMasks'][subdomain, :]
    
    # Load mesh variables
    with xr.open_dataset(paths['meshfile']) as ds:
        coords.update({
            'nCells': ds['nCells'].values[subdomain],
            'area': ds['areaCell'].values[subdomain],
            'lon': np.rad2deg(ds['lonCell'].values[subdomain]),
            'lat': np.rad2deg(ds['latCell'].values[subdomain]),
        })
    
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
        
        # Load flux variables
        results = {name: ds[prefix + name][0, :].values[subdomain] for name in varnames}
        
        # Load surface tracers
        for name in ['salinity', 'temperature']:
            results[name] = ds[prefix + 'activeTracers_' + name][0, :, 0].values[subdomain]

    return results


def build_sigma_bins(sigmarange, binsize):
    """Build sigma classes array
    """
    
    # Build sigma classes array
    nbins = int(abs(np.subtract(*sigmarange)) / binsize) + 1
    sigmabins = np.arange(nbins) * binsize + sigmarange[0]
    
    return sigmabins


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
        'sigma': rho(S, T, P) - 1000,
        'heat_factor': drhodt(S, T, P) / rho0 / cpsw,
        'salt_factor': -drhods(S, T, P) / rho0 * S,
    }

    return statevars


def calc_buoyancy_fluxes(results, fluxdefs):
    """Calculate buoyancy fluxes according to fluxdefs
    Also accounts for latent heat of fusion for melt fluxes.
    """
    
    # Latent heat of fusion [J kg-1]
    L_fusion = -3.337e5
    
    # Init buoyancy flux dict
    bfluxes = {}
    
    # Heat flux contributions
    for name in fluxdefs['heat']:
        bfluxes[name] = results['heat_factor'] * results[name]
    
    # Heat flux contributions due to melting/freezing
    for name in fluxdefs['melt']:
        name_appended = name.replace('Flux', 'HeatFlux')
        bfluxes[name_appended] = results['heat_factor'] * results[name] * L_fusion
    
    # Salt flux contributions
    for name in fluxdefs['salt']:
        bfluxes[name] = results['salt_factor'] * results[name]
    
    return bfluxes


def calc_transformation(bflux, area, mask, binsize):
    """Calculate water mass transformation per unit density
    using area integral
    """
    
    # Calculate transformation term F (in Sv)
    F = np.sum(bflux[mask] * area[mask]) / binsize * 1e-6
    
    return F


def write_output(mpasvars, wmtrvars, coords, paths):
    """Write wmtf results to netCDF. Filename conventions are taken from
    the results file prefix and the daterange
    """
    
    # Filename strings
    mesh = paths['meshfile'].split('/')[-2]
    prefix = '_'.join(paths['prefix'].split('_')[:3])
    datestr = '_'.join(time.strftime('%Y%m%d') for time in coords['time'][[0, -1]])
    
    # ---- Write MPAS variables to netCDF --------
    coordinates = {
        'time'       : ('time', coords['time']),
        'nCells'     : ('nCells', coords['nCells']),
        'regionNames': ('regionNames', coords['regionNames']),
    }
    datavars = {
        'lon'        : ('nCells', coords['lon']),
        'lat'        : ('nCells', coords['lat']),
        'area'       : ('nCells', coords['area']),
        'regionMasks': (['nCells', 'regionNames'], coords['regionMasks']),
    }
    datavars.update(
        {name: (['time', 'nCells'], np.vstack(mpasvars[name])) for name in mpasvars}
    )
    filename = os.path.join(paths['out'], f'{prefix}_{mesh}.mpasvars_{datestr}.nc')
    xr.Dataset(datavars, coordinates).to_netcdf(filename)
    # --------------------------------------------

    # ---- Write wmtr variables to netCDF --------
    coordinates = {
        'time'       : ('time', coords['time']),
        'regionNames': ('regionNames', coords['regionNames']),
        'sigmaBins'  : ('sigmaBins', coords['sigmaBins']),
    }
    dims = list(coordinates)
    shape = [len(coords[name]) for name in dims]
    datavars = (
        {name: (dims, np.array(wmtrvars[name]).reshape(shape)) for name in wmtrvars}
    )
    filename = os.path.join(paths['out'], f'{prefix}_{mesh}.wmtrvars_{datestr}.nc')
    xr.Dataset(datavars, coordinates).to_netcdf(filename)
    # --------------------------------------------


def run_wmtr(pathsfile, sigmarange=[19, 29], binsize=0.05):
    """Run the water mass transformation code for years and paths dict specified.
    `paths` must include `results`, `prefix`, `meshfile` and `maskfile` fields
    """

    # Open yaml file and load paths dict
    with open(pathsfile, 'r') as f:
        paths = yaml.safe_load(f)

    # Filenames, mesh variables, sigma bins, flux definitions
    varnames, fluxdefs = load_variable_definitions()
    coords, subdomain = load_MPASO_mesh(paths)
    coords['time'], dateranges = build_time_array(paths['results'])
    coords['sigmaBins'] = build_sigma_bins(sigmarange, binsize)

    # Initialize storage dictionaries
    mpasvars, wmtrvars = {}, {}

    # Loop through filenames
    for date in tqdm(coords['time']):
        
        # Load results
        filename = build_MPASO_filename(date, dateranges, paths)
        results = load_MPASO_results(filename, varnames, subdomain)
        
        # Calculate state variables and append to results
        statevars = calc_state_variables(results['salinity'], results['temperature'])
        results.update(statevars)
        
        # Append results to list
        for name in results:
            if name not in mpasvars:
                mpasvars[name] = []
            mpasvars[name].append(results[name])
        
        # Calculate buoyancy fluxes
        bfluxes = calc_buoyancy_fluxes(results, fluxdefs)

        # Loop through regions
        for rmask in coords['regionMasks'].T:

            # Apply region mask to flux variables
            bf_r = {name: bfluxes[name][rmask] for name in bfluxes}
            sigma_r, area_r = results['sigma'][rmask], coords['area'][rmask]

            # Loop through density bins
            for sigmabin in coords['sigmaBins']:

                # Create density mask
                dmask = (sigma_r >= sigmabin) & (sigma_r <= sigmabin + binsize)
                
                # Loop through buoyancy flux contributions
                for name in bf_r:
                    
                    # Calculate transformation term F (area integral over dmask)
                    F = calc_transformation(bf_r[name], area_r, dmask, binsize)
                    
                    # Append buoyancy flux contribution to list
                    if name not in wmtrvars:
                        wmtrvars[name] = []
                    wmtrvars[name].append(F)

    # Save output
    write_output(mpasvars, wmtrvars, coords, paths)


if __name__ == "__main__":
    run_wmtr(sys.argv[1])