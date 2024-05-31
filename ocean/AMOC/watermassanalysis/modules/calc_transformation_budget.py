#!/usr/bin/env python
"""
    Name: calc_transformation_budget.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to calculate the SPNA transformation budget
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np
import xarray as xr
import pandas as pd
import yaml
import time
from datetime import datetime
import postprocesstools as pptools
import watermasstools as wmttools


def calc_volumetric_TS(temperature, salinity, volume, coords, bins):
    """Calculate volumetric TS for MPAS results. Uses `numpy.histogramdd` to
    create the volume-weighted 2D TS histogram.
    """

    # Loop through regions
    volumeTS = []
    for regionMask in coords['regionMasks'].T:
        S, T, V = [var[regionMask, :].ravel() for var in (salinity, temperature, volume)]
        vol, _ = np.histogramdd((T, S), weights=V, bins=bins)
        volumeTS.append(vol) # km3
    
    # Xarray output
    coordinates = {'regionNames': coords['regionNames'], 'temperatureBins': bins[0][:-1], 'salinityBins': bins[1][:-1]}
    variables = {'volumetricTS': (coordinates.keys(), np.array(volumeTS) * 1e-9)}
    ds_out = xr.Dataset(variables, coordinates)
    
    return ds_out


def calc_volume(volume, sigmaTheta, coords, regions, sigmarange=[27, 28.5], binsize=0.01):
    """Calculate time-averaged overturning transformation over specified sigma bins
    and return in sigma space as an `xr.Dataset`
    """
    
    # Build sigma variables
    sigmabins = wmttools.build_sigma_bins(sigmarange, binsize)
    
    # Initialize variables dict
    volumeSigma = []

    # Loop through sigmabins
    for sigmabin in sigmabins:
        
        sigmaMask = (sigmaTheta >= sigmabin) & (sigmaTheta <= sigmabin + binsize)
        
        volumeRegions = []
        for region in regions:
            index, = np.where(coords['regionNames'] == region)
            regionMask = coords['regionMasks'][:, index[0]][:, None] & sigmaMask
            volumeRegions.append(volume[regionMask].sum())
        volumeSigma.append(volumeRegions)

    # Xarray output
    coordinates = {'sigmaBins': sigmabins, 'regionNames': regions}
    variables = {'volume': (coordinates.keys(), np.array(volumeSigma) * 1e-9)}
    ds_out = xr.Dataset(variables, coordinates)
    
    return ds_out


def calc_overturning(transectVars, sigmarange=[27, 28.5], binsize=0.01):
    """Calculate time-averaged overturning transformation over specified sigma bins
    and return in sigma space as an `xr.Dataset`
    """
    
    # Build sigma variables
    sigmabins = wmttools.build_sigma_bins(sigmarange, binsize)
    
    # Initialize variables dict
    overturning = []

    # Loop through sigmabins
    for sigmabin in sigmabins:
        
        transformation = []
        for transect in transectVars:
            sigmaTheta, transport = [transectVars[transect][name] for name in ('sigmaTheta', 'transport')]
            sigmaMask = (sigmaTheta >= sigmabin) & (sigmaTheta <= sigmabin + binsize)
            transformation.append(transport[sigmaMask].sum())
        overturning.append(transformation)

    # Xarray output
    coordinates = {'sigmaBins': sigmabins, 'transectNames': list(transectVars.keys())}
    variables = {'overturningTransformation': (coordinates.keys(), np.array(overturning)[::-1, :].cumsum(axis=0)[::-1, :] * 1e-6)}
    ds_out = xr.Dataset(variables, coordinates)
    
    return ds_out


def load_transect_vars(ds, sigmaTheta, layerThickness, transectMasks, coords, subdomain):
    """
    """
    
    # Get transport variables
    variables = {}
    for transect in transectMasks:

        # Get transect edge indices and signs
        index, sign = [transectMasks[transect][name] for name in ('index', 'sign')]

        # Get normal Velocity
        velocity = ds.timeMonthly_avg_normalVelocity[0, index, :].values

        # Get density and layer thickness and interpolate to edge
        cellsOnTransect = coords['cellsOnEdge'][index]
        sigmaThetaTransect = interpolate_to_edge(sigmaTheta, cellsOnTransect, subdomain)
        layerThicknessTransect = interpolate_to_edge(layerThickness, cellsOnTransect, subdomain)

        # Calculate overturning variables
        variables[transect] = {
            'sigmaTheta': sigmaThetaTransect,
            'transport': velocity * layerThicknessTransect * (coords['dvEdge'][index] * sign)[:, None],
        }
    
    return variables


def main_routine():
    """
    """
    
    # Save path
    savepath = '/pscratch/sd/b/bmoorema/results/aggregated/transformation/surface/'

    # Loop though meshes
    for mesh in ['LR', 'HR']:
        
        # Get paths
        with open(f'../yaml/paths_{mesh}.yaml') as f:
            paths = yaml.safe_load(f)
        
        # Load coords
        coords, transectMasks, subdomain = load_coords(paths)

        # Loop through decades
        startyear = 1946
        for decade in [(1947, 1957), (1997, 2007)]:
            
            # Define results path
            decade_str = str(decade[0]) + '-' + str(decade[1])
            resultspath = paths['results'][decade_str] + '/' + paths['prefix']
            
            # Define times
            times = [datetime(year, month, 1) for year in range(*decade) for month in range(1, 13)]
            
            # Loop through times
            print(f'Loading {mesh}, decade {decade_str}...')
            wmt, n, starttime = [], len(times), time.time()
            for k, t in enumerate(times):

                # Load results
                filename = resultspath + f'.{t.year-startyear:04d}-{t.month:02d}-01.nc'
                with xr.open_dataset(filename) as ds:

                    # Get cell variables
                    prefix = 'timeMonthly_avg_'
                    sigmaTheta = ds[prefix + 'potentialDensity'][0, ...].values[subdomain, :] - 1000
                    #layerThickness = ds[prefix + 'layerThickness'][0, ...].values[subdomain, :]
                    temperature = ds[prefix + 'activeTracers_temperature'][0, ...].values[subdomain, :]
                    salinity = ds[prefix + 'activeTracers_salinity'][0, ...].values[subdomain, :]

                    # Load transect variables
                    #transectVars = load_transect_vars(ds, sigmaTheta, layerThickness, transectMasks, coords, subdomain)

                # Calculate volume
                #volume = layerThickness * coords['area'][:, None]

                # Get state variables and buoyancy fluxes
                sigmaSurface, heatFactor, saltFactor, SSS = wmttools.calc_state_variables(salinity[:, 0], temperature[:, 0])
                fluxes = wmttools.build_combined_fluxes(ds, heatFactor, saltFactor, SSS, subdomain=subdomain)

                # Calculate 1D water mass transformation over regions
                ds_out = wmttools.calc_wmt(
                    fluxes, sigmaSurface, coords, regions=coords['regionNames'], sigmarange=[27, 28.5], binsize=0.01,
                )

                # Calculate overturning transformation over transects
                #ds_out = xr.merge([ds_out, calc_overturning(transectVars)])

                # Calculate volume
                #ds_out = xr.merge([ds_out, calc_volume(volume, sigmaTheta, coords, coords['regionNames'])])
                
                # Calculate volumetric TS
                #bins = [np.arange(-3, 20.1, 0.1), np.arange(33, 37.01, 0.01)]
                #ds_out = xr.merge([ds_out, calc_volumetric_TS(temperature, salinity, volume, coords, bins)])

                wmt.append(ds_out)
                
                # Update loop
                pptools.loopstatus(k, n, starttime)

            # Concatenate time
            wmt = xr.concat(wmt, pd.Index(times, name='time'))
            wmt.to_netcdf(savepath + f'transformationbudget_{mesh}_{decade_str}.nc')


if __name__ == "__main__":
    main_routine()
