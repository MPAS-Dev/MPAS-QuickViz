#!/usr/bin/env python
"""
    Name: watermasstools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for calculating water mass transformation
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import xarray as xr
import yaml
import fastjmd95 as jmd95
import postprocesstools as pptools


def build_sigma_bins(sigmarange, binsize):
    """Build sigma classes array given `sigmarange` and `binsize`
    Return as `pandas.Index` object
    """
    
    # Build sigma classes array
    nbins = int(abs(np.subtract(*sigmarange)) / binsize) + 1
    sigmabins = np.arange(nbins) * binsize + sigmarange[0]
    
    return sigmabins


def calc_state_variables(S, T, P=0):
    """State Equation from Jackett and McDougal 1995
    Use package fastjmd95 by R. Abernathey and J. Busecke
    """

    # Define constants
    rho0 = 1026.0     # Seawater density constant [kg m-3]
    cpsw = 3.996e3    # Heat capacity of seawater [J kg-1 K-1]

    # Calculate state variables
    sigmaTheta = jmd95.rho(S, T, P) - 1000
    heatFactor = jmd95.drhodt(S, T, P) / rho0 / cpsw
    saltFactor = -jmd95.drhods(S, T, P) / rho0

    return sigmaTheta, heatFactor, saltFactor


def build_combined_fluxes(
    ds_in, heatFactor, saltFactor, salinity, subdomain=None, prefix='timeMonthly_avg_',
    fluxdefs='../yaml/variable_combinations.yaml',
):
    """Build combined surface flux variables based on the combined definitions
    in the `fluxdefs` yaml file. Also build the buoyancy fluxes explicitly.
    """
    
    # Define constants
    rho0 = 1026.0     # Seawater density constant [kg m-3]
    rho_fw = 1e3      # Freshwater density constant [kg m3]
    
    # Define subdomain index
    index = subdomain if subdomain is not None else slice(None, None)
    
    # Open variable definitions
    with open(fluxdefs, 'r') as f:
        ctgys = yaml.safe_load(f)
    
    # Load fluxes from xarray, fill missing variables as zero
    fluxes = {}
    for ctgy, varNames in ctgys.items():
        fluxes[ctgy] = 0
        for varName in varNames:
            try:
                variable = ds_in[prefix + varName][0, :].values[index]
            except KeyError:
                variable = 0
            fluxes[ctgy] = fluxes[ctgy] + variable
    
    # Buoyancy flux due to heat
    fluxes['buoyancyHeatFlux'] = heatFactor * fluxes['totalHeatFlux']
    
    # Buoyancy flux due to salt
    fluxes['buoyancySaltFlux'] = saltFactor * (
        salinity * fluxes['totalFreshFlux'] -             # FW mass flux times salinity
        1e3 * rho_fw / rho0 * fluxes['totalSaltFlux'] -   # Salt mass flux times 1e3 * rho_fw / rho_sw
        rho_fw * fluxes['totalSalinityFlux']              # Salinity flux times rho_fw
    )
    
    # Total buoyancy flux
    fluxes['buoyancyTotalFlux'] = fluxes['buoyancyHeatFlux'] + fluxes['buoyancySaltFlux']
    
    return fluxes


def calc_wmt(
    fluxes, sigmaTheta, coords, regions=None,
    remapvars=None, sigmarange=[21, 29], binsize=0.1,
):
    """Calculate time-averaged water mass transformation over specified sigma bins,
    remap to lon, lat and return in sigma-lonlat space as an `xr.Dataset`
    """
    
    # Build sigma variables
    sigmabins = build_sigma_bins(sigmarange, binsize)
    
    # Define unit scale (1D: Sv, 2D: 1e-6 Sv km-2)
    scale = 1e-6 if remapvars is None else 1e6
    
    # Initialize variables dict
    ctgys = ['heat', 'salt', 'total']
    variables = {ctgy + 'Transformation': [] for ctgy in ctgys}
    
    # Initialize lonlat if 2D
    if remapvars is not None:
        lon, lat = None, None

    # ----------------------------------------------------------------------
    # Loop through sigmabins
    for sigmabin in sigmabins:

        # Create sigma mask
        sigmaMask = (sigmaTheta >= sigmabin) & (sigmaTheta <= sigmabin + binsize)

        # Loop through flux categories
        for ctgy in ctgys:
            
            # Extract flux and apply sigma mask
            flux = fluxes[f'buoyancy{ctgy.capitalize()}Flux'][sigmaMask]
            
            # --------------------------------------------------------------
            # Return spatial transformation map (2D output)
            if remapvars is not None:
                
                # Remap to lonlat, get coordinates and convert to numpy
                transformation = pptools.remap(flux, coords['nCells'][sigmaMask], **remapvars)
                if lon is None or lat is None:
                    lon, lat = [transformation[name].values for name in ('lon', 'lat')]
                transformation = transformation.values
            
            # --------------------------------------------------------------
            # Integrate over sigma outcrop (1D output)
            else:
                
                # Multiply by area
                flux = flux * coords['areaCell'][sigmaMask]
                
                # Calculate by region
                if regions is not None:
                    
                    # Loop over regions
                    transformation = []
                    for region in regions:
                        index, = np.where(coords['regionNames'] == region)
                        regionMask = coords['regionCellMasks'][sigmaMask, index[0]]
                        transformation.append(flux[regionMask].sum(axis=0))
                
                # Calculate global transformation only
                else:
                    transformation = flux.sum(axis=0)

            # --------------------------------------------------------------
            # Append to list
            variables[ctgy + 'Transformation'].append(transformation)

    # ----------------------------------------------------------------------
    # Prepare xarray output
    #
    # Define xr.Dataset coordinates
    coordinates = {'sigmaBins': sigmabins}
    if regions is not None:
        coordinates['regionNames'] = regions
    elif remapvars is not None:
        coordinates['lat'] = lat
        coordinates['lon'] = lon
    
    # Calculate final transformation and assign to tuple with dims
    dims = coordinates.keys()
    for name in variables:
        variables[name] = (dims, np.array(variables[name]) / binsize * scale)
    
    # Build xr.Dataset
    ds_out = xr.Dataset(variables, coordinates)
    
    return ds_out