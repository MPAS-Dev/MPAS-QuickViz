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
from dateutil.parser import parse
from fastjmd95 import rho, drhodt, drhods
from tqdm import tqdm


def load_flux_definitions():
    """Flux definitions for MPAS-Ocean
        TODO: include `seaIceSalinityFlux`
    """

    # Flux definitions (positive into the ocean)
    fluxdefs = {
        'heat': {
            'fluxes': [                 #### - Air-sea heat fluxes ----------- #######
                'shortWaveHeatFlux',    # Short wave radiation flux ---------- [W m-2]
                'longWaveHeatFluxUp',   # Upward long wave heat flux --------- [W m-2]
                'longWaveHeatFluxDown', # Downward long wave heat flux ------- [W m-2]
                'latentHeatFlux',       # Latent heat flux ------------------- [W m-2]
                'sensibleHeatFlux',     # Sensible heat flux ----------------- [W m-2]
            ],
            'meltfluxes': [             #### - Latent heat of melting fluxes - #######
                'snowFlux',             # Fresh water flux from snow --------- [kg m-2 s-1]
                'iceRunoffFlux',        # Fresh water flux from ice runoff --- [kg m-2 s-1]
            ],
        },
        'heat_ice': {
            'fluxes': [                 #### - Ice-ocean heat fluxes --------- #######
                'seaIceHeatFlux',       # Sea ice heat flux ------------------ [W m-2]
            ],
            'meltfluxes': [             #### - Latent heat of melting fluxes - #######
                'seaIceFreshWaterFlux', # Fresh water flux from sea ice ------ [kg m-2 s-1]
            ],
        },
        'salt': {
            'fluxes': [                 #### - Air-sea freshwater fluxes ----- #######
                'evaporationFlux',      # Evaporation flux ------------------- [kg m-2 s-1]
                'rainFlux',             # Fresh water flux from rain --------- [kg m-2 s-1]
                'snowFlux',             # Fresh water flux from snow --------- [kg m-2 s-1]
                'riverRunoffFlux',      # Fresh water flux from river runoff - [kg m-2 s-1]
                'iceRunoffFlux',        # Fresh water flux from ice runoff --- [kg m-2 s-1]
            ],
        },
        'salt_ice': {
            'fluxes': [                 #### - Ice-ocean freshwater fluxes --- #######
                'seaIceFreshWaterFlux', # Fresh water flux from sea ice ------ [kg m-2 s-1]
            ],
        },
    }
    
    # Flatten to list of unique flux variable names
    fluxnames = [component for flux in fluxdefs.values() for component in flux.values()]
    fluxnames = list(np.unique(np.hstack(fluxnames)))

    return fluxdefs, fluxnames


def build_sigma_bins(sigmarange, binsize):
    """Build sigma classes array
    """
    
    # Build sigma classes array
    nbins = int(abs(np.subtract(*sigmarange)) / binsize) + 1
    sigmabins = np.arange(nbins) * binsize + sigmarange[0]
    
    return sigmabins


def aggregate_fluxes(results, fluxdef):
    """Wrapper for summing flux variables according to those listed in fluxdef.
    Also accounts for latent heat of fusion for melt fluxes.
    """
    
    # Latent heat of fusion
    L_fusion = -3.337e5
    
    # Aggregate direct fluxes
    flux = sum([results[name] for name in fluxdef['fluxes']])
    
    # Aggregate indirect heat fluxes via melting
    if 'meltfluxes' in fluxdef:
        flux = flux + sum([results[name] for name in fluxdef['meltfluxes']]) * L_fusion
    
    return flux


def calc_buoyancy_fluxes(results, fluxdefs, P=0):
    """Calculate buoyancy fluxes following according to
    Jeong et al (2020), J Clim
    
    State Equation from Jackett and McDougal 1995
    Use package fastjmd95 by R. Abernathey and J. Busecke

         `pip install fastjmd95`
    """    
    
    # Define constants
    rho_0 = 1026    # --- MPAS value confirmed
    cpsw = 3.996e3  # --- MPAS value confirmed

    # Calculate state variables
    S, T = [results[name] for name in ('salinity', 'temperature')]
    sigma = rho(S, T, P) - 1000
    bfactors = {
        'salt': -drhods(S, T, P) / rho_0 * S,
        'heat': drhodt(S, T, P) / rho_0 / cpsw,
    }
    
    # Loop through flux categories
    bfluxes = {}
    for ctgy in fluxdefs:
        
        # Calculate buoyancy flux
        coeff = bfactors[ctgy.replace('_ice', '')]
        bfluxes[ctgy] = coeff * aggregate_fluxes(results, fluxdefs[ctgy])

    return bfluxes, sigma


def calc_tr(bflux, area, mask, binsize):
    """Calculate water mass transformation per unit density
    using area integral
    """
    
    # Calculate transformation term F (in Sv)
    F = np.sum(bflux[mask] * area[mask]) / binsize * 1e-6
    
    return F


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


def load_MPASO_results(resultsfile, fluxnames, subdomain, prefix='timeMonthly_avg_'):
    """Load MPAS-Ocean results variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Open results file and extract variables
    with xr.open_dataset(resultsfile) as ds:
        
        # Load flux variables
        results = {name: ds[prefix + name][0, :].values[subdomain] for name in fluxnames}
        
        # Load surface tracers
        for name in ['salinity', 'temperature']:
            results[name] = ds[prefix + 'activeTracers_' + name][0, :, 0].values[subdomain]

    return results


def write_output(mpasvars, wmtrvars, coords, paths):
    """Write wmtf results to netCDF. Filename conventions are taken from
    the results file prefix and the daterange
    """

    # Filename strings
    prefix = paths['prefix'].split('.')[0]
    datestr = '_'.join(time.strftime('%Y%m%d') for time in coords['time'][[0, -1]])
    
    # Write MPAS variables to netCDF
    filename = os.path.join(paths['out'], f'{prefix}.mpasvars_{datestr}.nc')
    coordinates = {name: (name, coords[name]) for name in ['time', 'nCells', 'regionNames']}
    datavars = {name: ('nCells', coords[name]) for name in ['lon', 'lat', 'area']}
    datavars['regionMasks'] = (['nCells', 'regionNames'], coords['regionMasks'])
    datavars.update({name: (['time', 'nCells'], np.vstack(mpasvars[name])) for name in mpasvars})
    ds = xr.Dataset(datavars, coordinates)
    ds.to_netcdf(filename)

    # Write wmtr variables to netCDF
    filename = os.path.join(paths['out'], f'{prefix}.wmtrvars_{datestr}.nc')
    dims = ['time', 'regionNames', 'sigmaBins']
    shape = [len(coords[name]) for name in dims]
    coordinates = {name: (name, coords[name]) for name in dims}
    datavars = {name: (dims, np.array(wmtrvars[name]).reshape(shape)) for name in wmtrvars}
    ds = xr.Dataset(datavars, coordinates)
    ds.to_netcdf(filename)


def run_wmtr(pathsfile, sigmarange=[19, 29], binsize=0.1):
    """Run the water mass transformation code for years and paths dict specified.
    `paths` must include `results`, `prefix`, `meshfile` and `maskfile` fields
    """

    # Open yaml file and load paths dict
    with open(pathsfile, 'r') as f:
        paths = yaml.safe_load(f)

    # Filenames, mesh variables, sigma bins, flux definitions
    fluxdefs, fluxnames = load_flux_definitions()
    coords, subdomain = load_MPASO_mesh(paths)
    coords['time'], dateranges = build_time_array(paths['results'])
    coords['sigmaBins'] = build_sigma_bins(sigmarange, binsize)

    # Initialize storage dictionaries
    mpasvars = {name: [] for name in fluxnames + ['salinity', 'temperature']}
    wmtrvars = {name: [] for name in fluxdefs}

    # Loop through filenames
    for date in tqdm(coords['time']):
        
        # Load results and append to list
        filename = build_MPASO_filename(date, dateranges, paths)
        results = load_MPASO_results(filename, fluxnames, subdomain)
        for name in results:
            mpasvars[name].append(results[name])
        
        # Calculate buoyancy fluxes
        bfluxes, sigma = calc_buoyancy_fluxes(results, fluxdefs)

        # Loop through regions
        for rmask in coords['regionMasks'].T:

            # Apply region mask to flux variables
            bf_r = {name: bfluxes[name][rmask] for name in bfluxes}
            sigma_r, area_r = sigma[rmask], coords['area'][rmask]

            # Loop through density bins
            for sigmabin in coords['sigmaBins']:

                # Create density mask and calculate transformation term F (area integral)
                dmask = (sigma_r >= sigmabin) & (sigma_r <= sigmabin + binsize)
                for name in fluxdefs:
                    wmtrvars[name].append(calc_tr(bf_r[name], area_r, dmask, binsize))

    # Save output
    write_output(mpasvars, wmtrvars, coords, paths)


if __name__ == "__main__":
    run_wmtr(sys.argv[1])