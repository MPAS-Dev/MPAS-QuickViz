#!/usr/bin/env python
"""
    Name: calc_transformation_budget.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to calculate the SPNA transformation budget
    for the ImPACTS Water Mass Analysis project.
"""

import os
import sys
import yaml
from argparse import ArgumentParser
from datetime import datetime
from itertools import zip_longest
from tqdm import tqdm

import numpy as np
import xarray as xr
import pandas as pd

pkgroot = '/global/cfs/cdirs/m4259/bmoorema/MPAS-QuickViz/ocean/AMOC/watermassanalysis'
sys.path.append(os.path.join(pkgroot, 'modules'))

import postprocesstools as pptools
import transecttools as trtools
import watermasstools as wmttools


def define_args():
    """Define arguments for command-line use of this module
    """

    # Construct args object
    parser = ArgumentParser(description='Ocean model postprocessing suite')
    parser.add_argument('parameters', help='Parameters file path')
    parser.add_argument('filenumber', help='Process filenumber (zero start)')
    parser.add_argument('-p', '--path', type=str, default='./', help='Save path (default ./)')

    return parser.parse_args()


def build_results_filenames(filenumber, params, startyear=1948):
    """Build results filename and dates from filenumber
    """
    
    # Parse dates from filenumber
    year, month = filenumber // 12 + 1, filenumber % 12 + 1
    datestamp = datetime(startyear + year - 1, month, 1)
    datecode = f'{year:04d}-{month:02d}'

    # Get params for building the filename structure
    keys = ['model', 'simName', 'prefix', 'resultsPath']
    model, simName, prefix, resultsPath = [params[key] for key in keys]

    # Parse multiple run directories
    if not isinstance(resultsPath, str):
        for runDir in resultsPath[1]:
            years = [int(year) for year in runDir.split('_')[-2:]]
            if years[0] <= year <= years[1]:
                resultsPath = os.path.join(resultsPath[0], runDir)
                break
        else:
            msg = (
                f'Date {datecode} corresponding to record number {filenumber} '
                f'not available in results record for simulation {simName}'
            )
            raise ValueError(msg)

    # Build filenames
    if model == 'MPAS':
        filenames = {}
        for ctgy in ['ocean', 'seaice']:
            filename = '.'.join([simName, params[ctgy + 'Name'], prefix, f'{datecode}-01', 'nc'])
            filenames[ctgy] = os.path.join(resultsPath, filename)
    else:
        filenames = os.path.join(resultsPath, '.'.join([simName, datecode, 'nc']))
    
    return filenames, [datestamp, datecode]


def load_variables(ds, params, coords, depths, prefix='timeMonthly_avg_'):
    """Load variables and return appended lists
    """
    
    # Get kidx from depths
    kidx = [abs(coords.refBottomDepth.values - depth).argmin() for depth in depths]
    
    # Initialize data dict
    data = {ctgy: {} for ctgy in ['2D', '3D', 'region', 'transect', 'meridional']}
    
    # Get layer thickness from coords
    if coords.model == 'POP':
        layerThickness = coords.layerThickness.values[:, None, None]
        bottomDepth = np.cumsum(layerThickness, axis=0)
        data['3D']['volume'] = layerThickness * coords.areaCell.values[None, ...] * 1e-9
    
    # Loop through varTypes (cell2D, cell3D, edge3D, meridional, seaice)
    for varType in params:
    
        # Loop through varNames
        for varName, varDefs in params[varType].items():

            # Get varNameFull
            varNameFull = prefix + varDefs[0]
            conversion = varDefs[1]

            # Skip variable if not in dataset (e.g. GM, salinityRestoring)
            if varNameFull not in ds['ocean'] and varNameFull not in ds['seaice']:
                continue

            # Load into numpy
            key = 'seaice' if varType == 'seaice' else 'ocean'
            array = ds[key][varNameFull][0, ...].values * conversion
            
            # Get MOC and bypass remaining
            if varType == 'meridional':
                data['meridional'][varName] = array[0, ...].T
                continue
            
            # Apply subdomain in MPAS
            if coords.model == 'MPAS' and varType != 'edge3D':
                array = array[coords.nCells.values, ...]

            # 3D cell variables
            if '3D' in varType:

                # Get edge fields
                for ctgy in ['region', 'transect']:
                    if coords.model == 'MPAS':
                        interp = True if varType != 'edge3D' else False
                        data[ctgy][varName] = trtools.get_edge_fields(ctgy, array, coords, interp=interp)
                    elif coords.model == 'POP':
                        velocity = True if 'Velocity' in varName else False
                        data[ctgy][varName] = trtools.get_edge_fields_POP(ctgy, array, coords, velocity=velocity)
                
                # MPAS: bypass remaining if edge3D, otherwise reorder to depth axis=0
                if coords.model == 'MPAS':
                    if varType == 'edge3D':
                        continue
                    array = array.T

                # Process layerThickness and calculate volume [km3]
                if varName == 'layerThickness':
                    layerThickness = np.copy(array)
                    bottomDepth = np.cumsum(layerThickness, axis=0)
                    data['3D']['volume'] = \
                        layerThickness * coords.areaCell.values[None, :] * 1e-9
                    continue

                # Keep 3D variables for binning
                if varName in ['temperature', 'salinity', 'density']:
                    data['3D'][varName] = array

                # Get depth-weighted averages
                array = array * layerThickness
                array = np.array([np.nansum(array[:k+1, ...], axis=0) / bottomDepth[k, :] for k in kidx])

            # Append to list
            data['2D'][varName] = array
        
    # Rotate POP 2D velocities
    if coords.model == 'POP':
        angles = coords.angle.values[None, ...]
        for velType in ['Resolved', 'GM', 'Submeso']:
            uName, vName = [name + 'Velocity' + velType for name in ('u', 'v')]
            if uName in data['2D'] and vName in data['2D']:
                data['2D'][uName], data['2D'][vName] = pptools.rotate_velocities(
                    data['2D'][uName], data['2D'][vName], angles,
                )
        
    # Process edge velocities  
    for ctgy in ['region', 'transect']:
        Ctgy = ctgy.capitalize()
        
        # Aggregate velocity types
        velocityTotal = 0
        for velType in ['Resolved', 'GM', 'Submeso']:
            nName = 'velocityNormal' + velType
            uName, vName = [name + 'Velocity' + velType for name in ('u', 'v')]
            if uName in data[ctgy] and vName in data[ctgy]:
                data[ctgy][nName] = np.copy(data[ctgy]['uVelocity' + velType])
                idx = coords['vComponent' + Ctgy].values == 'v'
                data[ctgy][nName][idx, :] = data[ctgy]['vVelocity' + velType][idx, :]
            if nName in data[ctgy]:
                velocityTotal = velocityTotal + data[ctgy][nName]
        
        # Calculate transport
        dvEdge = coords['dv' + Ctgy].values[..., None]
        if coords.model == 'MPAS':
            layerThickness = data[ctgy]['layerThickness']
        elif coords.model == 'POP':
            layerThickness = coords.layerThickness.values[None, None, :]
        data[ctgy]['transport'] = velocityTotal * layerThickness * dvEdge * 1e-6
    
    return data


def bin_variables(data_in, binArgs, coords):
    """Bin variables by density and TS
    """
    
    # Binning parameters
    TSnames = ['temperature', 'salinity']
    params = {
        'Density': {
            'bins': [wmttools.get_bins_edges(*binArgs['density'])[1]],
        },
        'TS': {
            'bins': [wmttools.get_bins_edges(*binArgs[name])[1] for name in TSnames],
        },
    }
    
    # Initialize output dict
    data = {}
    
    # --- Transport along edges -------------
    # Loop through region edges and transects
    for ctgy in ['region', 'transect']:
        
        # Get binning arrays
        params['Density']['arrays'] = [data_in[ctgy]['density'] - 1000]
        params['TS']['arrays'] = [data_in[ctgy][name] for name in TSnames]
        
        # Loop through binning categories
        for dimName, lists in params.items():
            
            # Loop through transects and bin by transport
            H = []
            for transect in zip(*lists['arrays'], data_in[ctgy]['transport']):
                x = [x.ravel() for x in transect[:-1]]
                w = transect[-1].ravel()
                h, _ = np.histogramdd(x, bins=lists['bins'], weights=w)
                H.append(h)
            varName = f'binnedTransport{ctgy.capitalize()}{dimName}'
            data[varName] = np.array(H)
    
    # --- Area and volume over regions ------
    # Get areaCell and binning arrays
    areaCell = coords.areaCell.values * 1e-6
    params['Density']['arrays'] = [data_in['3D']['density'] - 1000]
    params['TS']['arrays'] = [data_in['3D'][name] for name in TSnames]
    
    # Loop through binning categories
    for dimName, lists in params.items():
    
        # Bin by area and volume over regions
        H = [[], []]
        for regionName in coords.regionNames.values:
            
            # Get region mask from region name
            regionMask = coords.regionCellMasks.sel(regionNames=regionName).values
            
            # Apply region mask to bin arrays
            arrays = [array[:, regionMask] for array in lists['arrays']]
        
            # Bin by area
            x = [x[0, :] for x in arrays]
            w = areaCell[regionMask]
            h, _ = np.histogramdd(x, bins=lists['bins'], weights=w)
            H[0].append(h)

            # Bin by volume
            x = [x.ravel() for x in arrays]
            w = data_in['3D']['volume'][:, regionMask].ravel()
            h, _ = np.histogramdd(x, bins=lists['bins'], weights=w)
            H[1].append(h)
    
        # Convert to numpy
        for name, h in zip(['Area', 'Volume'], H):
            varName = f'binned{name}Region{dimName}'
            data[varName] = np.array(h)
    
    return data


def calc_wmt(fluxes, coords, binArgs, remapvars):
    """Calculate 1D and 2D WMT for T, S and density
    """
    
    # Region names
    regionNames = coords.regionNames.values

    # Initialize data dict
    data = {'1D': {}, '2D': {}}

    # Get WMT binned by density
    data['1D']['densityTransformation'] = wmttools.calc_wmt(
        fluxes, coords, 'density', binArgs['density'],
        regionNames=regionNames,
    )

    # Get WMT binned by TS
    binNames = ['temperature', 'salinity']
    data['1D']['TSTransformation'] = wmttools.calc_wmt(
        fluxes, coords, binNames,
        [binArgs[name] for name in binNames],
        regionNames=regionNames,
    )
    
    # Merge 1D datasets
    data['1D'] = xr.merge(list(data['1D'].values()))

    # Get 2D WMT binned by all
    for binName in ('temperature', 'salinity', 'density'):
        data['2D'][binName + 'Transformation'] = wmttools.calc_wmt(
            fluxes, coords, binName, binArgs[binName],
            remapvars=remapvars,
        )

    return data


def build_output_datasets(data, coords, binArgs, depths, remapvars, datestrings, savepath='./'):
    """Build output data sets and save to netCDF
    """
    
    # Append subdirectory to savepath
    datestamp, datecode = datestrings
    savepath = os.path.join(savepath, 'monthly_files')

    # Save 2D fields
    ds = []
    for varName, values in data['2D'].items():
        if values.shape[0] == len(depths):
            values = [pptools.remap(layer, **remapvars) for layer in values]
            values = xr.concat(values, pd.Index(depths, name='depth'))
        else:
            values = pptools.remap(values, **remapvars)
        values.name = varName
        ds.append(values)
    ds = xr.merge(ds)
    filename = os.path.join(savepath, f'{coords.meshName}_2D_variables_{datecode}.nc')
    ds.expand_dims({'time': [datestamp]}).to_netcdf(filename, unlimited_dims='time')
    
    # Save 1D fields
    coordinates = coords[['regionNames', 'transectNames']]
    for key, values in binArgs.items():
        coordinates[key + 'Bins'] = wmttools.get_bins_edges(*values)[0]
    variables = {}
    binDims = {'Density': ['densityBins'], 'TS': ['temperatureBins', 'salinityBins']}
    for ctgy in ['region', 'transect']:
        for dimName, dimValues in binDims.items():
            dims = [ctgy + 'Names'] + dimValues
            for name in ['Area', 'Volume', 'Transport']:
                varName = f'binned{name}{ctgy.capitalize()}{dimName}'
                if varName in data['1D']:
                    variables[varName] = (dims, np.array(data['1D'][varName]))
    ds = xr.merge([xr.Dataset(variables, coordinates), data['wmt']['1D']])
    filename = os.path.join(savepath, f'{coords.meshName}_1D_{datecode}.nc')
    ds.expand_dims({'time': [datestamp]}).to_netcdf(filename, unlimited_dims='time')
    
    # Save transect fields
    names = ['latMOC', 'refBottomDepth'] + [coord for coord in coords if 'Transect' in coord]
    coordinates = coords[names]
    variables = {}
    dims = ['transectNames', 'nTransectEdges', 'refBottomDepth']
    for key, values in data['transect'].items():
        variables[key] = (dims, np.array(values))
    dims = ['latMOC', 'refBottomDepth']
    for key, values in data['meridional'].items():
        variables[key] = (dims, np.array(values))
    ds = xr.Dataset(variables, coordinates)
    filename = os.path.join(savepath, f'{coords.meshName}_transect_{datecode}.nc')
    ds.expand_dims({'time': [datestamp]}).to_netcdf(filename, unlimited_dims='time')
    
    # Save 2D WMT fields
    for ctgy, ds in data['wmt']['2D'].items():
        filename = os.path.join(savepath, f'{coords.meshName}_2D_{ctgy}_{datecode}.nc')
        ds.expand_dims({'time': [datestamp]}).to_netcdf(filename, unlimited_dims='time')

    
def process_monthly_file(meshName, filenumber, savepath='./'):
    """Process monthly file
    """
    
    # Get depths averages and binning params
    depths = [0, 100, 500, 1000, 1500]
    binArgs = {
        'temperature': (-2, 13, 0.1),
        'salinity': (33, 35.7, 0.01),
        'density': (26.5, 28.2, 0.01),
    }
    
    # Load coords file
    coords = xr.open_dataset(os.path.join(savepath, f'{meshName}_coords.nc'))
    
    # Load YAML params
    filename = os.path.join(pkgroot, 'yaml', f'parameters_{coords.model}.yaml')
    with open(filename, 'r') as f:
        params = yaml.safe_load(f)
    
    # Assign model specifics
    if coords.model == 'POP':
        grid, subdomain, varPrefix, zdim = 'T', None, '', 'z_t'
    elif coords.model == 'MPAS':
        grid, subdomain, varPrefix, zdim = None, coords.nCells.values, 'timeMonthly_avg_', 'nVertLevels'
    
    # Build remapping variables
    filename = build_coords_filename(params['run'], meshName=meshName)
    remapvars = pptools.build_remapper(
        filename, grid=grid, subdomain=subdomain, bbox=[-100, 45, 30, 85],
    )

    # Open results file, load and bin variables, get fluxes and calc wmt
    filenames, datestrings = build_results_filenames(filenumber, params['run'], meshName=meshName)
    ds = {name: xr.open_dataset(filename) for name, filename in filenames.items()}
    data = load_variables(ds, params['variables'], coords, depths, prefix=varPrefix)
    data['1D'] = bin_variables(data, binArgs, coords)
    fluxes = wmttools.build_fluxes(ds['ocean'], params['fluxes'], subdomain=subdomain, zdim=zdim, prefix=varPrefix)
    data['wmt'] = calc_wmt(fluxes, coords, binArgs, remapvars)
    exclude = [fluxes.pop(key) for key in ('temperature', 'salinity', 'density')]
    data['2D'].update(fluxes)
    
    # Build output datasets and save to netCDF
    build_output_datasets(data, coords, binArgs, depths, remapvars, datestrings, savepath=savepath)


if __name__ == "__main__":
    args = define_args()
    process_monthly_file(args.meshname, args.filenumber, args.path)
