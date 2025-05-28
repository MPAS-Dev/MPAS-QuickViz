#!/usr/bin/env python
"""
    Name: process_velocities.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Module to process velocities
    for the ImPACTS Water Mass Analysis project.
"""

import numpy as np


def calc_transport(data, coords, ctgy):
    """Calculate transport
    """
    
    # Get total velocity
    velocityTotal = [data[name] for name in data if 'velocityNormal' in name]
    velocityTotal = sum(velocityTotal)
    
    # Calculate transport
    model = coords.model_name
    dvEdge = coords[f'dv{ctgy.capitalize()}'].values[..., None]
    if model == 'MPAS':
        layerThickness = data['layerThickness']
    elif model == 'POP':
        layerThickness = coords.layerThickness.values[None, None, :]
    transport = velocityTotal * layerThickness * dvEdge * 1e-6

    return transport


def process_velocities_POP(data, coords):
    """Process velocities
    """

    # Rotate POP 2D velocities
    data['2D'] = _rotate_velocities(data['2D'], coords)

    # Process edge velocities
    for ctgy in ['region', 'transect']:
        if ctgy in data:
            data[ctgy] = _normal_from_uv(data[ctgy], coords, ctgy)

    return data


def _get_uv_names(data):
    """Get u,v name pairs from dict keys
    """

    # Get u,v name pairs
    comps = ['Zonal', 'Meridional']
    names = [[key for key in data if comp in key] for comp in comps]

    return names


def _rotate_velocities(data, coords):
    """Rotate velocities
    """

    # Get angles
    angles = coords.angle.values[None, ...]
    cosa, sina = np.cos(angles), np.sin(angles)

    # Rotate velocities
    for uname, vname in zip(*_get_uv_names(data)):
        u_rot = data[uname] * cosa - data[vname] * sina
        v_rot = data[uname] * sina + data[vname] * cosa
        data[uname], data[vname] = u_rot, v_rot

    return data


def _normal_from_uv(data, coords, ctgy):
    """Get normal velocities from u,v edge pairs
    """

    # Get normal velocities from u,v edge pairs
    vcomponents = coords[f'vComp{ctgy.capitalize()}'].values
    idx = vcomponents == 'v'
    for uname, vname in zip(*_get_uv_names(data)):
        nname = uname.replace('Zonal', 'Normal')
        data[nname] = np.copy(data[uname])
        data[nname][idx, :] = data[vname][idx, :]
        
    return data
