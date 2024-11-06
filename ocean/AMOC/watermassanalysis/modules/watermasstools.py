#!/usr/bin/env python
"""
    Name: watermasstools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for calculating water mass transformation
    for the ImPACTS Water Mass Analysis project.
"""

import yaml
from itertools import pairwise

import fastjmd95 as jmd95
import numpy as np
import pandas as pd
import xarray as xr

import postprocesstools as pptools


def define_constants():
    """Define constants for use in water mass transformation calculations
    """
    
    # Define constants
    cpsw = 3996.0     # Heat capacity of seawater [J kg-1 K-1]
    rho_sw = 1026.0   # Seawater density constant [kg m-3]
    rho_fw = 1000.0   # Freshwater density constant [kg m3]
    
    return cpsw, rho_sw, rho_fw


def calc_state_variables(S, T, P=0):
    """Calculate state variables following Jacket and McDougall 1995.
    """
    
    # Define constants
    cpsw, rho_sw, rho_fw = define_constants()
    
    # Calculate state variables
    density = jmd95.rho(S, T, P)
    alpha = -jmd95.drhodt(S, T, P) / rho_sw
    beta = jmd95.drhods(S, T, P) / rho_sw
    
    return density, alpha, beta


def get_bins_edges(start, stop, step):
    """Get bins and edges
    """
    
    # Get bins and edges
    bins = np.arange((stop - start) / step + 2) * step + start
    edges = bins - step / 2
    
    return bins[:-1], edges


def calc_bin_masks(array, edges):
    """Calculate array masks for each consecutive bin given by edges
    """
    
    # Calculate mask
    ndim = array.ndim
    masks = []
    for l, u in pairwise(edges):
        mask = (array >= l) & (array < u)
        if ndim == 1:
            mask, = np.where(mask)
        masks.append(mask)
    
    return masks


def parse_bin_args(binNames, binArgs, fluxNames):
    """Parse binning parameters and flux names from input args
    """

    # If multiple bin names assume TS space
    if isinstance(binNames, (list, tuple)):
        fluxNames = [name + 'Flux' for name in binNames]
        binSizes = [args[2] for args in binArgs]

    # Otherwise single binning dimension
    else:
        fluxNames = [name for name in fluxNames if binNames in name and 'Flux' in name]
        binSizes = [binArgs[2]] * len(fluxNames)
        binNames, binArgs = [binNames], [binArgs]

    return binNames, binArgs, binSizes, fluxNames


def build_fluxes(ds, fluxDefs, subdomain=None, zdim='nVertLevels', prefix='timeMonthly_avg_'):
    """Build surface flux variables based on the definitions in the
    `fluxNames` dict. Thermal and haline transformation terms are
    defined separately following Evans et al. 2014 JGR Oceans.
    """
    
    # Define constants
    cpsw, rho_sw, rho_fw = define_constants()
    
    # Load flux fields from dataset
    fluxes = {}
    for ctgy, fluxNames in fluxDefs.items():

        # Load surface salinity and temperature
        if ctgy in ['salinity', 'temperature']:
            fluxes[ctgy] = ds[prefix + fluxNames][0, ...].isel(**{zdim: 0}).values
            if subdomain is not None:
                fluxes[ctgy] = fluxes[ctgy][subdomain]

        # Load surface fluxes
        else:
            fluxes[ctgy] = 0
            for fluxName, varName in fluxNames.items():
                varNameFull = prefix + varName
                if varNameFull in ds:
                    flux = ds[varNameFull][0, :].values
                    if subdomain is not None:
                        flux = flux[subdomain]
                    if fluxName == 'total':
                        fluxes[ctgy] = flux
                    else:
                        fluxes[fluxName] = flux
                        if 'total' not in fluxNames:
                            fluxes[ctgy] = fluxes[ctgy] + flux
    
    # Calculate state variables
    density, alpha, beta = calc_state_variables(fluxes['salinity'], fluxes['temperature'])
    fluxes['density'] = density - 1000
    
    # Temperature mass flux [degC kg m-2 s-1]
    fluxes['temperatureFlux'] = fluxes['totalHeatFlux'] / cpsw
    
    # Salinity volume flux [PSU m s-1]
    fluxes['salinityFlux'] = (
        1e3 / rho_sw * fluxes['totalSaltFlux'] -                # g/kg SALT / rho_sw * kg SALT m-2 s-1
        fluxes['salinity'] / rho_fw * fluxes['totalFreshFlux']  # S / rho_fw * kg FW m-2 s-1
    )
    if 'totalSalinityFlux' in fluxes:
        fluxes['salinityFlux'] = fluxes['salinityFlux'] + fluxes['totalSalinityFlux']
    
    # Density fluxes [kg m-2 s-1]
    fluxes['densityHeatFlux'] = -alpha * fluxes['temperatureFlux']
    fluxes['densitySaltFlux'] = beta * rho_fw * fluxes['salinityFlux']
    fluxes['densityTotalFlux'] = fluxes['densityHeatFlux'] + fluxes['densitySaltFlux']
    
    # Convert temperature mass flux to volume flux [degC m s-1]
    fluxes['temperatureFlux'] = fluxes['temperatureFlux'] / density
    
    return fluxes


def calc_wmt(fluxes, coords, binNames, binArgs, remapvars=None, regionNames=None):
    """Calculate time-averaged water mass transformation over specified bins
    and return in sigma space as an `xr.Dataset`
    """
    
    # Parse binning parameters and flux names (assumes TS if binNames is iterable)
    binNames, binArgs, binSizes, fluxNames = parse_bin_args(binNames, binArgs, fluxes)
    
    # Build bins and edges
    bins, edges = zip(*[get_bins_edges(*args) for args in binArgs])
    
    # If 2D map, specify method and build bin masks manually
    if remapvars is not None:
        method = '2D'
        binArray = fluxes[binNames[0]]
        ndim = binArray.ndim
        binMasks = calc_bin_masks(binArray, edges[0])
        binIndex = pd.Index(bins[0], name=binNames[0] + 'Bins')
        remapvars_local = {key: value for key, value in remapvars.items() if key != 'subdomain'}
        subdomain = remapvars['subdomain'] if 'subdomain' in remapvars else None
    else:
        method = '1D'
        areaCell = coords.areaCell.values
    
    # Get region masks if regions requested, else default to open slice
    if regionNames is not None:
        regionMasks = []
        for regionName in regionNames:
            regionMask = coords.regionCellMasks.sel(regionNames=regionName).values.ravel()
            regionMasks.append(regionMask)
    else:
        regionMasks = [slice(None, None)]
    
    # Loop through fluxNames
    variables = {}
    for fluxName, binSize in zip(fluxNames, binSizes):
    
        # Init wmtName and get flux array from dict
        wmtName = fluxName.removesuffix('Flux') + 'Transformation'
        flux = fluxes[fluxName]
        
        # If 2D map, manually loop through binMasks and remap flux to lonlat
        if method == '2D':
            if 'Total' not in fluxName:
                wmt = []
                for binMask in binMasks:
                    if ndim > 1:
                        array = np.copy(flux)
                        array[np.isfinite(array) & ~binMask] = 0
                        subdomainMask = None
                    else:
                        array = flux[binMask]
                        subdomainMask = subdomain[binMask]
                    flux2D = pptools.remap(array, **remapvars_local, subdomain=subdomainMask)
                    wmt.append(flux2D)
                wmt = xr.concat(wmt, dim=binIndex) * 1e6   # 1e-6 Sv km-2
                wmt.name = wmtName
        
        # Otherwise area-integrate flux over regions and bin using `numpy` 1d or 2d histograms
        else:
            flux, wmt = (flux * areaCell).ravel(), []
            arrays = [fluxes[name].ravel() for name in binNames]
            for regionMask in regionMasks:
                x = [array[regionMask] for array in arrays]
                h, _ = np.histogramdd(x, bins=edges, weights=flux[regionMask])
                wmt.append(h)
            wmt = np.array(wmt).squeeze() * 1e-6   # Sv

        # Assign to variables dict
        variables[wmtName] = wmt / binSize

    # Build xarray output from remap
    if method == '2D':
        
        # Merge DataArrays from remap output
        ds = xr.merge(variables.values())
    
    # Build xarray output from concatenated histograms
    else:
        
        # Build coordinates dict
        coordinates = {}
        if regionNames is not None:
            coordinates['regionNames'] = regionNames
        for name, bn in zip(binNames, bins):
            coordinates[name + 'Bins'] = bn
    
        # Build variables dict and make DataSet
        dims = coordinates.keys()
        variables = {name: (dims, variables[name]) for name in variables}
        ds = xr.Dataset(variables, coordinates)
    
    return ds
