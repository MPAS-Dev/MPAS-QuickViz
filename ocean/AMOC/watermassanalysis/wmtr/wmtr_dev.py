# WMTR Development Module

import numpy as np
import xarray as xr
import os
import sys
import yaml
from datetime import datetime
from fastjmd95 import rho, drhodt, drhods
from tqdm import tqdm


def build_MPASO_filenames(years, paths, startyear=1946):
    """Build a list of sequential MPAS-Ocean filenames over a range of years
    """
    
    # Build times and filenames
    prefix = os.path.join(paths['results'], paths['prefix'])
    times, filenames = [], []
    for year in range(*years):
        for month in range(1, 13):
            times.append(datetime(year + startyear, month, 1))
            filenames.append(prefix + f'.{year:04d}-{month:02d}-01.nc')
    
    return np.array(times), filenames, len(filenames)


def build_sigma_bins(sigmarange, binsize):
    """Build sigma classes array
    """
    
    # Build sigma classes array
    nbins = int(abs(np.subtract(*sigmarange)) / binsize) + 1
    sigmabins = np.arange(nbins) * binsize + sigmarange[0]
    
    return sigmabins, nbins


def build_sigma_mask(sigma, bounds):
    """Build sigma mask
    """
    
    # Define bounds
    try:
        lower, upper = bounds
    except:
        lower, upper = bounds[0], 1e32

    # Build sigma mask
    mask = (sigma >= lower) & (sigma <= upper)
    
    return mask


def calc_buoyancy_fluxes(variables, P=0):
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
    S, T = [variables[name] for name in ('salinity', 'temperature')]
    sigma = rho(S, T, P) - 1000
    salt_factor = -drhods(S, T, P) / rho_0 * S
    heat_factor = drhodt(S, T, P) / rho_0 / cpsw
    
    # Calculate buoyancy fluxes
    bfluxes = {
        'heat': heat_factor * variables['heat'],
        'salt': salt_factor * variables['salt'],
        'heat_ice': heat_factor * variables['heat_ice'],
        'salt_ice': salt_factor * variables['salt_ice'],
    }

    return bfluxes, sigma


def calc_tr(bflux, area, mask, binsize):
    """Calculate water mass transformation per unit density
    using area integral
    """
    
    # Calculate transformation term F (in Sv)
    F = np.sum(bflux[mask] * area[mask]) / binsize * 1e-6
    
    return F


def load_flux_definitions():
    """Flux definitions for MPAS-Ocean
    """

    # Flux definitions
    fluxnames = {
        'heat': [
            'shortWaveHeatFlux',    # Short wave radiation flux ---------- [W m-2]
            'longWaveHeatFluxUp',   # Upward long wave heat flux --------- [W m-2]
            'longWaveHeatFluxDown', # Downward long wave heat flux ------- [W m-2]
            'latentHeatFlux',       # Latent heat flux ------------------- [W m-2]
            'sensibleHeatFlux',     # Sensible heat flux ----------------- [W m-2]
        ],
        'salt': [
            'evaporationFlux',      # Evaporation flux ------------------- [kg m-2 s-1]
            'rainFlux',             # Fresh water flux from rain --------- [kg m-2 s-1]
            'riverRunoffFlux',      # Fresh water flux from river runoff - [kg m-2 s-1]
            'snowFlux',             # Fresh water flux from snow --------- [kg m-2 s-1] (also affects heat)
        ],
        'heat_ice': [
            'seaIceHeatFlux',       # Sea ice heat flux ------------------ [W m-2]
        ],
        'salt_ice': [
            'seaIceSalinityFlux',   # Sea ice salinity flux -------------- [kg m-2 s-1] (PSU m s-1, check netcdf, pull request?)
            'seaIceFreshWaterFlux', # Fresh water flux from sea ice ------ [kg m-2 s-1]
            'iceRunoffFlux',        # Fresh water flux from ice runoff --- [kg m-2 s-1]
        ],
    }

    return fluxnames


def load_MPASO_mesh(paths):
    """Load MPAS-Ocean mesh variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Load coordinates
    with xr.open_dataset(paths['meshfile']) as ds:
        area = ds.areaCell.values
        #lon = np.rad2deg(ds.lonCell.values)
        #lat = np.rad2deg(ds.latCell.values)

    # Load region masks
    with xr.open_dataset(paths['maskfile']) as ds:
        regionnames = ds.regionNames.values.astype(str)
        regionmasks = ds.regionCellMasks.values.astype(bool).T
    
    # All points to be considered (for memory purposes)
    subdomain = np.sum(regionmasks, axis=0).astype(bool)
    area, regionmasks = area[subdomain], regionmasks[:, subdomain]

    return area, subdomain, regionnames, regionmasks, len(regionnames)


def load_MPASO_results(resultsfile, fluxnames, subdomain, prefix='timeMonthly_avg_'):
    """Load MPAS-Ocean results variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Open results file and extract variables
    variables = {}
    with xr.open_dataset(resultsfile) as ds:

        # Load surface tracers
        for name in ['salinity', 'temperature']:
            variables[name] = ds[prefix + 'activeTracers_' + name][0, :, 0].values[subdomain]

        # Load surface fluxes
        for name in fluxnames:
            variables[name] = sum([ds[prefix + nm][0, :].values[subdomain] for nm in fluxnames[name]])

    return variables


def parse_params(params_path):
    """Parse analysis parameters from the input yaml file
    """
    
    # Open yaml file to nested dict
    with open(params_path, 'r') as f:
        params = yaml.safe_load(f)
    
    # Parse analysis params from dict
    paths = params['paths']
    years = list(params['yearrange'].values())
    sigmarange = list(params['sigmaparams'].values())
    binsize = sigmarange.pop(2)
    
    return paths, years, sigmarange, binsize


def write_output(wmtr, coords, paths):
    """Write wmtf results to netCDF. Filename conventions are taken from
    the results file prefix and the daterange
    """

    # Date string for output filename
    datestr = '_'.join(time.strftime('%Y%m%d') for time in coords[0][[0, -1]])
    filename = paths['prefix'].split('.')[0] + f'.wmtr_{datestr}.nc'
    filename = os.path.join(paths['out'], filename)

    # Construct coordinates and variables dicts
    dims = ['time', 'regions', 'sigmabins']
    coordinates = {name: (name, coord) for name, coord in zip(dims, coords)}
    variables = {name: (dims, wmtr[name]) for name in wmtr}

    # Write to netCDF
    ds = xr.Dataset(variables, coordinates)
    ds.to_netcdf(filename)


def run_wmtr(params_path):
    """Run the water mass transformation code for years and paths dict specified.
    `paths` must include `results`, `prefix`, `meshfile` and `maskfile` fields
    """

    # Parse parameters
    paths, years, sigmarange, binsize = parse_params(params_path)

    # Filenames, mesh variables, sigma bins, flux definitions
    times, filenames, nfiles = build_MPASO_filenames(years, paths)
    area, subdomain, regionnames, regionmasks, nregions = load_MPASO_mesh(paths)
    sigmabins, nbins = build_sigma_bins(sigmarange, binsize)
    fluxnames = load_flux_definitions()        

    # Initialize transformation dict
    wmtr = {name: np.zeros((nfiles, nregions, nbins)) for name in fluxnames}

    # Loop through filenames
    for t, filename in zip(tqdm(range(nfiles)), filenames):

        # Load results and calculate buoyancy fluxes
        variables = load_MPASO_results(filename, fluxnames, subdomain)
        bfluxes, sigma = calc_buoyancy_fluxes(variables)

        # Loop through regions
        for r, rmask in zip(range(nregions), regionmasks):

            # Apply region mask to flux variables
            bf_r = {name: bfluxes[name][rmask] for name in bfluxes}
            sigma_r, area_r = sigma[rmask], area[rmask]

            # Loop through density bins
            for i in range(nbins):

                # Create density mask and calculate transformation term F (area integral)
                dmask = build_sigma_mask(sigma_r, sigmabins[i:i+2])
                for name in fluxnames:
                    wmtr[name][t, r, i] = calc_tr(bf_r[name], area_r, dmask, binsize)

    # Save output
    write_output(wmtr, [times, regionnames, sigmabins], paths)


if __name__ == "__main__":
    run_wmtr(sys.argv[1])