#!/usr/bin/env python
"""
    Name: utilities.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Shared tools for model results postprocessing
    for the ImPACTS Water Mass Analysis project.
"""

import os
import time
from copy import deepcopy

import numpy as np
import pyremap
import xarray as xr
from scipy import signal


def loopstatus(k, n, starttime, interval=1):
    """Print loop progress at percentage intervals given an iteration k
    and a total number of iterations n
    """
    
    nscale = 100 / n
    percent = k * nscale
    if percent % interval < nscale:
        time_elapsed = time.time() - starttime
        msg = f'{int(percent)}% complete ... {time_elapsed:.0f} seconds elapsed'
        print(msg, flush=True)


def lowpass(data, cutoff=10, window_type='boxcar'):
    """Apply a Finite Impulse Response (FIR) lowpass filter according
    to the window_type and cutoff using a convolution algorithm
    """
    
    # Calculate low-pass signal
    window = signal.get_window(window_type, cutoff)
    filtered = np.convolve(data, window / sum(window), mode='same')
    
    return filtered


def downsample(array, widths=(5, 5), func='mean'):
    """Downsample array using a groupby mean every w elements.
    Pads the array dimensions with nan so that `numpy.reshape` can be used
    in the groupby operation.
    """

    # Resample array
    pad = [(0, w-dim%w) for dim, w in zip(array.shape, widths)]
    array = np.pad(array, pad, constant_values=np.nan)
    args = [arg for dim, w in zip(array.shape, widths) for arg in (int(dim/w), w)]
    axis = tuple(i for i, e in enumerate(args) if e in widths)
    array = getattr(np, 'nan' + func)(array.reshape(*args), axis=axis)
    
    return array


def rotate_velocities(u, v, angle):
    """Rotate velocities
    """
    
    # Rotate velocity
    cosa, sina = np.cos(angle), np.sin(angle)
    u_rot = u * cosa - v * sina
    v_rot = u * sina + v * cosa
    
    return u_rot, v_rot
