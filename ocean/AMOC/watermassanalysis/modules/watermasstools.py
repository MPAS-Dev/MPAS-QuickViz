#!/usr/bin/env python
"""
    Name: watermasstools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for calculating water mass transformation
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import xarray as xr
import pandas as pd
import yaml
import fastjmd95 as jmd95
import postprocesstools as pptools
from itertools import pairwise


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
    mask = [np.where((array >= l) & (array < u))[0] for l, u in pairwise(edges)]
    
    return mask


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


def build_fluxes(
    ds_in, subdomain=None, prefix='timeMonthly_avg_',
    fluxdefs='../yaml/variable_combinations.yaml',
):
    """Build surface flux variables based on the definitions in the
    `fluxdefs` yaml file. Thermal and haline transformation terms are
    defined separately following Evans et al. 2014 JGR Oceans.
    """
    
    # Define constants
    cpsw, rho_sw, rho_fw = define_constants()
    
    # Define subdomain index
    index = subdomain if subdomain is not None else slice(None, None)
    
    # Load salinity and temperature
    varNames = ['activeTracers_salinity', 'activeTracers_temperature']
    S, T = [ds_in[prefix + varName][0, :, 0].values[index] for varName in varNames]
    
    # Calculate state variables
    density, alpha, beta = calc_state_variables(S, T)
    
    # Open flux definitions yaml file
    with open(fluxdefs, 'r') as f:
        ctgys = yaml.safe_load(f)
    
    # Load fluxes, assign missing fluxes to zero
    fluxes = {}
    for ctgy, varNames in ctgys.items():
        fluxes[ctgy] = 0
        for varName in varNames:
            try:
                variable = ds_in[prefix + varName][0, :].values[index]
            except KeyError:
                variable = 0
            fluxes[ctgy] = fluxes[ctgy] + variable
    
    # Temperature mass flux [degC kg m-2 s-1]
    fluxes['temperatureFlux'] = fluxes['totalHeatFlux'] / cpsw
    
    # Salinity volume flux [PSU m s-1]
    fluxes['salinityFlux'] = (
        fluxes['totalSalinityFlux'] +              # PSU m s-1
        1e3 / rho_sw * fluxes['totalSaltFlux'] -   # g/kg SALT / rho_sw * kg SALT m-2 s-1
        S / rho_fw * fluxes['totalFreshFlux']      # S / rho_fw * kg FW m-2 s-1
    )
    
    # Density fluxes [kg m-2 s-1]
    fluxes['densityHeatFlux'] = -alpha * fluxes['temperatureFlux']
    fluxes['densitySaltFlux'] = beta * rho_fw * fluxes['salinityFlux']
    fluxes['densityTotalFlux'] = fluxes['densityHeatFlux'] + fluxes['densitySaltFlux']
    
    # Convert temperature mass flux to volume flux [degC m s-1]
    fluxes['temperatureFlux'] = fluxes['temperatureFlux'] / density
    
    # Add state variables to fluxes
    fluxes.update({'salinity': S, 'temperature': T, 'density': density - 1000})
    
    return fluxes


def calc_wmt(fluxes, coords, binNames, binArgs, remapvars=None, regions=None):
    """Calculate time-averaged water mass transformation over specified bins
    and return in sigma space as an `xr.Dataset`
    """
    
    # Unpack coordinate variables
    names = ['nCells', 'areaCell', 'regionNames', 'regionCellMasks']
    nCells, areaCell, regionNames, regionCellMasks = [coords[name] for name in names]
    
    # Parse binning parameters and flux names (assumes TS if binNames is iterable)
    binNames, binArgs, binSizes, fluxNames = parse_bin_args(binNames, binArgs, fluxes)
    
    # Build bins and edges
    bins, edges = zip(*[get_bins_edges(*args) for args in binArgs])
    
    # If 2D map, specify method and build bin masks manually
    if remapvars is not None:
        method = '2D'
        binMasks = calc_bin_masks(fluxes[binNames[0]], edges[0])
        binIndex = pd.Index(bins[0], name=binNames[0] + 'Bins')
    else:
        method = '1D'
    
    # Get region masks if regions requested, else default to open slice
    if regions is not None:
        regionMasks = [regionCellMasks[:, regionNames == region][:, 0] for region in regions]
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
            wmt = [pptools.remap(flux[mask], nCells[mask], **remapvars) for mask in binMasks]
            wmt = xr.concat(wmt, dim=binIndex) * 1e6   # 1e-6 Sv km-2
            wmt.name = wmtName + '2D'
        
        # Otherwise area-integrate flux over regions and bin using `numpy` 1d or 2d histograms
        else:
            flux, wmt = flux * areaCell, []
            for mask in regionMasks:
                x = [fluxes[name][mask] for name in binNames]
                h, _ = np.histogramdd(x, bins=edges, weights=flux[mask])
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
        if regions is not None:
            coordinates['regionNames'] = regions
        for name, bn in zip(binNames, bins):
            coordinates[name + 'Bins'] = bn
    
        # Build variables dict and make DataSet
        dims = coordinates.keys()
        variables = {name: (dims, variables[name]) for name in variables}
        ds = xr.Dataset(variables, coordinates)
    
    return ds
