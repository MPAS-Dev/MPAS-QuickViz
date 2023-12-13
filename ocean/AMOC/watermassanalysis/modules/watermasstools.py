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
from fastjmd95 import rho, drhodt, drhods
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


def build_combined_fluxes(
    ds_in, statevars, prefix='timeMonthly_avg_',
    fluxdefs='../yaml/variable_combinations.yaml',
):
    """Build combined surface flux variables based on the combined definitions
    in the `vardefs` yaml file. Also build the buoyancy fluxes explicitly.
    """
    
    # Open variable definitions
    with open(fluxdefs, 'r') as f:
        ctgys = yaml.safe_load(f)
    
    # Combine variables
    fluxes = {ctgy: sum([ds_in[prefix + name][0, ...] for name in ctgys[ctgy]]) for ctgy in ctgys}
    fluxes['buoyancyHeatFlux'] = statevars['heatFactor'] * fluxes['totalHeatFlux']
    fluxes['buoyancySaltFlux'] = statevars['saltFactor'] * fluxes['totalSaltFlux']
    fluxes['buoyancyTotalFlux'] = fluxes['buoyancyHeatFlux'] + fluxes['buoyancySaltFlux']
    
    return fluxes


def calc_wmt(fluxes, sigmaTheta, areaCell=None, remapvars=None, sigmarange=[21, 29], binsize=0.1):
    """Calculate time-averaged water mass transformation over specified sigma bins,
    remap to lon, lat and return in sigma-lonlat space as an `xr.Dataset`
    """
    
    # Build sigma variables
    sigmabins = build_sigma_bins(sigmarange, binsize)
    
    # Define unit scale (1D: Sv, 2D: 1e-6 Sv km-2)
    scale = 1e-6 if remapvars is None else 1e6
    
    # Initialize lists
    data = {ctgy: [] for ctgy in ['heat', 'salt', 'total']}

    # Loop through sigmabins
    for sigmabin in sigmabins:

        # Create sigma mask
        mask = (sigmaTheta >= sigmabin) & (sigmaTheta <= sigmabin + binsize)

        # Loop through flux categories
        for ctgy in data:
            
            # Mask sigma and average over time, then assign nan to zero
            da = fluxes[f'buoyancy{ctgy.capitalize()}Flux'].where(mask)
            
            if remapvars is None:  # 1D output
                
                # Integrate over outcrop
                da = (da * areaCell).sum(dim='nCells')

            else:  # 2D output
                
                # Time average -> replace nan with zero -> remap to lonlat
                da = da.mean(dim='time')
                da = pptools.remap(da.where(~np.isnan(da), 0), **remapvars)
            
            # Append to list
            data[ctgy].append(da)

    # Concatenate DataArrays and calculate final quantities
    wmt = {}
    for ctgy in data:
        transformation = xr.concat(data[ctgy], sigmabins) / binsize * scale
        wmt[ctgy + 'Trans'] = transformation #.isel(sigmaBins=slice(0, -1))
        #wmt[ctgy + 'Form'] = -transformation.diff(dim='sigmaBins', label='lower')
    
    return wmt