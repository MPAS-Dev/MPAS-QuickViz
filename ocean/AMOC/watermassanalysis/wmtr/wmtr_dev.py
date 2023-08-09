# WMTR Development Module

import numpy as np
import xarray as xr
import os
import sys
import yaml
import glob
from datetime import datetime
from dateutil.parser import parse
from fastjmd95 import rho, drhodt, drhods
from tqdm import tqdm


def build_MPASO_filenames(paths, startyear=1946):
    """Build a list of sequential MPAS-Ocean filenames over a range of years
    """
    
    # Build filenames
    filenames = []
    for path in paths['results']:
        filenames.extend(glob.glob(os.path.join(path, paths['prefix'] + '*')))
    filenames = np.array(sorted(filenames))
    
    # Find only unique filenames
    _, index = np.unique([os.path.split(f)[-1] for f in filenames], return_index=True)
    filenames = filenames[index]
    
    # Build times
    times = []
    for filename in filenames:
        time = parse(filename.split('.')[-2])
        times.append(time.replace(year=time.year + startyear))
    
    return filenames, np.array(times)


def build_sigma_bins(sigmarange, binsize):
    """Build sigma classes array
    """
    
    # Build sigma classes array
    nbins = int(abs(np.subtract(*sigmarange)) / binsize) + 1
    sigmabins = np.arange(nbins) * binsize + sigmarange[0]
    
    return sigmabins


def build_sigma_mask(sigma, sigmabin, binsize):
    """Build sigma mask
    """

    # Build sigma mask
    mask = (sigma >= sigmabin) & (sigma <= sigmabin + binsize)
    
    return mask


def calc_buoyancy_fluxes(results, P=0):
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
    salt_factor = -drhods(S, T, P) / rho_0 * S
    heat_factor = drhodt(S, T, P) / rho_0 / cpsw
    
    # Calculate buoyancy fluxes
    bfluxes = {
        'heat': heat_factor * results['heat'],
        'fresh': salt_factor * results['fresh'],
        'heat_ice': heat_factor * results['heat_ice'],
        'fresh_ice': salt_factor * results['fresh_ice'],
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
            'snowFlux',             # Fresh water flux from snow --------- [kg m-2 s-1]
            'iceRunoffFlux',        # Fresh water flux from ice runoff --- [kg m-2 s-1]
        ],
        'fresh': [
            'evaporationFlux',      # Evaporation flux ------------------- [kg m-2 s-1]
            'rainFlux',             # Fresh water flux from rain --------- [kg m-2 s-1]
            'riverRunoffFlux',      # Fresh water flux from river runoff - [kg m-2 s-1]
            'snowFlux',             # Fresh water flux from snow --------- [kg m-2 s-1]
            'iceRunoffFlux',        # Fresh water flux from ice runoff --- [kg m-2 s-1]
        ],
        'heat_ice': [
            'seaIceHeatFlux',       # Sea ice heat flux ------------------ [W m-2]
            'seaIceFreshWaterFlux', # Fresh water flux from sea ice ------ [kg m-2 s-1]
        ],
        'fresh_ice': [
            'seaIceFreshWaterFlux', # Fresh water flux from sea ice ------ [kg m-2 s-1]
        ],
        #'salt_ice': ['seaIceSalinityFlux'],   # Sea ice salinity flux -------------- [kg m-2 s-1]
    }

    return fluxnames


def load_MPASO_mesh(paths):
    """Load MPAS-Ocean mesh variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Load region masks
    with xr.open_dataset(paths['maskfile']) as ds:
        regionnames = ds.regionNames.values.astype(str)
        regionmasks = ds.regionCellMasks.values.astype(bool).T
    
    # All points to be considered (for memory purposes)
    subdomain = np.sum(regionmasks, axis=0).astype(bool)
    regionmasks = regionmasks[:, subdomain]
    
    # Load coordinates
    with xr.open_dataset(paths['meshfile']) as ds:
        lons = np.rad2deg(ds.lonCell.values)[subdomain]
        lats = np.rad2deg(ds.latCell.values)[subdomain]
        area = ds.areaCell.values[subdomain]
    
    # Correct lons
    lons[lons > 180] = lons[lons > 180] - 360

    return lons, lats, area, subdomain, regionnames, regionmasks


def load_MPASO_results(resultsfile, fluxnames, subdomain, prefix='timeMonthly_avg_'):
    """Load MPAS-Ocean results variables needed for calculating water mass
    transformation/formation. Use `'maxLevelCell` for landmask if needed.
    """

    # Latent heat of fusion (J kg-1)
    latent_heat_fusion = -3.337e5

    # Open results file and extract variables
    results = {}
    with xr.open_dataset(resultsfile) as ds:

        # Load surface tracers
        for name in ['salinity', 'temperature']:
            results[name] = ds[prefix + 'activeTracers_' + name][0, :, 0].values[subdomain]

        # Load surface fluxes
        meltfluxnames = ['snowFlux', 'iceRunoffFlux', 'seaIceFreshWaterFlux']
        for ctgy in fluxnames:

            # Loop through category
            fluxes = []
            for name in fluxnames[ctgy]:
                flux = ds[prefix + name][0, :].values[subdomain]

                # Apply melting to heat fluxes
                if ('heat' in ctgy) and (name in meltfluxnames):
                    flux = flux * latent_heat_fusion

                # Append to list
                fluxes.append(flux)

            # Sum fluxes in category
            results[ctgy] = sum(fluxes)

    return results


def parse_params(params_path):
    """Parse analysis parameters from the input yaml file
    """
    
    # Open yaml file to nested dict
    with open(params_path, 'r') as f:
        params = yaml.safe_load(f)
    
    # Parse analysis params from dict
    paths = params['paths']
    sigmarange = list(params['sigmaparams'].values())
    binsize = sigmarange.pop(2)
    
    return paths, sigmarange, binsize


def write_output(mpas_vars, wmtr_vars, coordinates, paths):
    """Write wmtf results to netCDF. Filename conventions are taken from
    the results file prefix and the daterange
    """

    # Date string for output filename
    datestr = '_'.join(time.strftime('%Y%m%d') for time in coordinates[0][[0, -1]])
    filename = paths['prefix'].split('.')[0] + f'.wmtr_{datestr}.nc'
    filename = os.path.join(paths['out'], filename)

    # Construct coordinates and variables dicts
    dims = ['Time', 'Regions', 'Sigmabins']
    coords = {name: (name, coord) for name, coord in zip(dims, coordinates)}

    # Reshape wmtr output and construct data_vars dict
    shape = [len(coord) for coord in coordinates]
    data_vars = {name + '_tr': (dims, np.array(wmtr_vars[name]).reshape(shape)) for name in wmtr_vars}
    
    # Add mpas_vars to data_vars dict
    coordnames = ['lon', 'lat', 'area']
    varnames = [name for name in list(mpas_vars) if name not in coordnames]
    data_vars.update({name: ('nCells', mpas_vars[name]) for name in coordnames})
    data_vars.update({name: (['Time', 'nCells'], np.vstack(mpas_vars[name])) for name in varnames})

    # Write to netCDF
    ds = xr.Dataset(data_vars, coords)
    ds.to_netcdf(filename)


def run_wmtr(params_path):
    """Run the water mass transformation code for years and paths dict specified.
    `paths` must include `results`, `prefix`, `meshfile` and `maskfile` fields
    """

    # Parse parameters
    paths, sigmarange, binsize = parse_params(params_path)

    # Filenames, mesh variables, sigma bins, flux definitions
    filenames, times = build_MPASO_filenames(paths)
    lons, lats, area, subdomain, regionnames, regionmasks = load_MPASO_mesh(paths)
    sigmabins = build_sigma_bins(sigmarange, binsize)
    fluxnames = load_flux_definitions()

    # Initialize storage dictionaries
    mpas_vars = {'lon': lons, 'lat': lats, 'area': area, 'salinity': [], 'temperature': []}
    mpas_vars.update({name: [] for name in fluxnames})
    wmtr_vars = {name: [] for name in fluxnames}

    # Loop through filenames
    for filename in tqdm(filenames):

        # Load results
        results = load_MPASO_results(filename, fluxnames, subdomain)
        for name in results:
            mpas_vars[name].append(results[name])
        
        # Calculate buoyancy fluxes
        bfluxes, sigma = calc_buoyancy_fluxes(results)

        # Loop through regions
        for rmask in regionmasks:

            # Apply region mask to flux variables
            bf_r = {name: bfluxes[name][rmask] for name in bfluxes}
            sigma_r, area_r = sigma[rmask], area[rmask]

            # Loop through density bins
            for sigmabin in sigmabins:

                # Create density mask and calculate transformation term F (area integral)
                dmask = build_sigma_mask(sigma_r, sigmabin, binsize)
                for name in fluxnames:
                    wmtr_vars[name].append(calc_tr(bf_r[name], area_r, dmask, binsize))

    # Save output
    write_output(mpas_vars, wmtr_vars, [times, regionnames, sigmabins], paths)


if __name__ == "__main__":
    run_wmtr(sys.argv[1])