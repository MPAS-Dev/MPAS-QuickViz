#!/usr/bin/env python
"""
    Name: calc_transformation_budget.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Executable module to calculate the SPNA transformation budget
    for the ImPACTS Water Mass Analysis project.
"""

from datetime import datetime

import numpy as np
import xarray as xr
import yaml

from load_coordinates import load_coords
from load_variables import load_variables, add_surface_fluxes
from process_velocities import process_velocities_POP
from binning import bin_by_region
import build_outputs


def process_results_file(paramsFile, filenumber, coordsFile=None):
    """Process results file
    """

    # Load YAML variable definitions
    with open('../yaml/variable_metadata.yaml', 'r') as f:
        definitions = yaml.safe_load(f)['variables']

    # Load YAML params
    with open(paramsFile, 'r') as f:
        params = yaml.safe_load(f)
    model = params['run']['model_name']

    # Load coords
    #filename = os.path.join(savepath, f"{params['run']['short_name']}_coords.nc")
    #if os.path.exists(filename):
    if coordsFile is not None:
        coords = xr.open_dataset(coordsFile)
    else:
        coords = load_coords(params)

    # Get datestamp and datecode
    datestamp, datecode = _get_datestamp(filenumber)

    # Load variables
    if model == 'MPAS':
        prefix = 'timeMonthly_avg_'
    elif model == 'POP':
        prefix = ''
    data = load_variables(params, coords, datecode, prefix=prefix)

    # Process velocity
    if model == 'POP' and any('velocity' in var for var in params['variables']):
        data = process_velocities_POP(data, coords)

    # Binning
    if 'binning' in params:
        data = bin_by_region(data, coords, params['binning'])

    # Add fluxes
    dfns, names = definitions.items(), params['variables']
    if any('ctgy' in dfn for name, dfn in dfns if name in names):
        data['2D'] = add_surface_fluxes(data['2D'], params['variables'])

    # Build outputs
    ds = {'2D': build_outputs.build_output_2D(
        datestamp, data['2D'], coords, params,
    )}
    if 'binning' in params:
        ds['binned'] = build_outputs.build_output_binned(
            datestamp, data['binned'], coords, params,
        )
    if 'masking' in params:
        ds['transect'] = build_outputs.build_output_transect(
            datestamp, data, coords, params,
        )

    return ds


def _get_datestamp(filenumber, startyear=1948):
    """Get datestamp from filenumber
    """

    # Parse dates from filenumber
    year = filenumber // 12 + 1
    month = filenumber % 12 + 1
    datestamp = datetime(startyear + year - 1, month, 1)
    datecode = f'{year:04d}-{month:02d}'

    return datestamp, datecode
        