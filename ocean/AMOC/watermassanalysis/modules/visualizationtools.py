#!/usr/bin/env python
"""
    Name: visualizationtools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Visualization tools for the ImPACTS Water Mass Analysis project
"""

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt, patches
from tqdm import tqdm


def parse_contour_levels(clims):
    """Parse `clims` params into arrays for setting contour levels and ticks.
    
    clims = [lower, upper, contour_inc, cbar_inc, residual_bound, contour_inc, cbar_inc]
    """
    
    # Return None if clims is None
    if clims is None:
        return None, None, None, None
    
    # Else parse contour levels and ticks
    else:
        
        # Increment to make upper bound inclusive
        i = 1e-6
    
        # Left panel
        l, u = clims[:2]
        levels1, ticks1 = [np.arange(l, u + i, d) for d in clims[2:4]]

        # Right panel
        u = clims[4]
        levels2, ticks2 = [np.arange(l, u + i, d) for l, d in zip([-u, 0], clims[5:])]
        ticks2 = np.unique(np.hstack([-ticks2, ticks2]))

    return levels1, levels2, ticks1, ticks2


def plot_variable_spatial(ds_in, varname, units, scale=1, month=None, clims=None, cmap=None):
    """Plot specified variable over spatial region. Hard-coded for two timeranges,
    two meshes and the residual between meshes. These categories must be consistent
    with the structure of the `plotvars` dictionary.
    """
    
    if month is None:
        month = slice(None, None)
    
    # General definitions
    meshes, tstrs = list(ds_in), list(list(ds_in.values())[0])[::-1]
    
    # Sigma contour level specs
    sigma_levels = {
        'levels'    : np.arange(25, 28, 0.5),
        'linestyles': ['-', '--', '-', '--', '-', '--'],
        'colors'    : ['k', 'k', 'gray', 'gray', 'lightgray', 'lightgray'],
    }
    
    # Make plot area
    fig, axs = plt.subplots(2, 3, figsize=(12, 8), gridspec_kw={'hspace': 0.05, 'wspace': 0.05})
    
    # Add titles
    for col, title in zip(axs.T, meshes + [f'{meshes[0]}-{meshes[1]}']):
        col[0].set_title(title)
    
    # Loop through time ranges
    for row, tstr in zip(axs, tstrs):
        
        # Print time range
        row[0].text(0.01, 0.8, tstr, transform=row[0].transAxes)
        
        # Loop through meshes
        residual = []
        for ax, mesh in zip(row, meshes):
            
            # Select months
            ds = ds_in[mesh][tstr].sel(months=month).mean(dim='months')
            
            # Extract variables
            names = ('lon', 'lat', 'sigmaTheta', varname)
            lon, lat, sigma, variable = [ds[name] for name in names]

            # Plot variable
            c1 = ax.pcolormesh(lon, lat, variable * scale, vmin=clims[0], vmax=clims[1], cmap=cmap)

            # Plot sigma contours
            cs = ax.contour(lon, lat, sigma, **sigma_levels)
            
            residual.append(variable)
        
        # Plot residual
        residual = np.subtract(*residual)
        c2 = row[2].pcolormesh(lon, lat, residual * scale, vmin=-clims[2], vmax=clims[2], cmap='RdBu_r')
    
    # Gray background and remove ticks
    for ax in axs.ravel():
        ax.set_xlim(-100, 20)
        ax.set_ylim(0, 80)
        ax.add_patch(patches.Rectangle([0, 0], 1, 1, fc='gray', transform=ax.transAxes, zorder=-10))
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

    # Add colorbars
    cax1 = fig.add_axes([0.14, 0.05, 0.48, 0.015])
    cax2 = fig.add_axes([0.66, 0.05, 0.23, 0.015])
    fig.colorbar(c1, cax=cax1, label=units, orientation='horizontal')
    fig.colorbar(c2, cax=cax2, label=units, orientation='horizontal')
    for x, mesh in zip([0.9, 0.1], meshes):
        cax2.text(x, 1.5, mesh, transform=cax2.transAxes, ha='center')
    
    # Add sigma contour legend
    for level, ls, color in zip(*list(sigma_levels.values())):
        axs[1, 0].plot(0, 0, ls, color=color, label=level)
    axs[1, 0].legend(ncols=7, loc=(0.08, -0.1))
    
    return fig, axs