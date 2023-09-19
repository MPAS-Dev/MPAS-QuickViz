#!/usr/bin/env python
"""
    Name: visualizationtools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Visualization tools for the ImPACTS Water Mass Analysis project
"""

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt, patches
from datetime import datetime
from tqdm import tqdm


def parse_plot_args(ds, depths, months, sigmabins, clims):
    """Parse arguments in the call to `plot_2Dvariable` into the
    needed variable types and assign defaults.
    """
    
    # Default values range
    range_def = [None, None]
    
    # Parse meshes and timestrings
    meshes = list(ds)
    timestrings = list(ds[meshes[0]])[::-1]
    
    # Parse months and sigmabins into slice args
    # name: (aggregation_function, values)
    sliceargs = {
        'depths'   : ('', 0 if depths is None else depths),
        'months'   : ('mean', slice(*range_def) if months is None else months),
        'sigmaBins': ('sum' , slice(*range_def) if sigmabins is None else sigmabins),
    }
    
    # Parse clims
    if clims is None:
        clims = [range_def, range_def]
    else:
        clims = [clims[:2], [-clims[2], clims[2]]]
    
    return meshes, timestrings, sliceargs, clims


def plot_2Dpanel(
    ax, ds, varname=None, scale=1, xlim=[-100, 20], ylim=[0, 80],
    clim=[None, None], cmap=None, plot_sigma=False,
):
    """Plot 2D variable given by `varname` and overlay sigma
    contours if requested. If no varname provided, `ds` is
    treated as an `xr.DataArray` and thus the plotting variable.
    """
    
    # Sigma contour styling
    sigma_levels = {
        'levels'    : np.arange(25, 28, 0.5),
        'linestyles': ['-', '--', '-', '--', '-', '--'],
        'colors'    : ['k', 'k', 'gray', 'gray', 'lightgray', 'lightgray'],
    }
    
    # Plot variable
    variable = ds[varname] if varname is not None else ds
    c = ax.pcolormesh(ds.lon, ds.lat, variable * scale, vmin=clim[0], vmax=clim[1], cmap=cmap)

    # Plot sigma contours and add dummy plots for legend
    if plot_sigma:
        ax.contour(ds.lon, ds.lat, ds.sigmaTheta, **sigma_levels)
        for level, ls, color in zip(*list(sigma_levels.values())):
            ax.plot(0, 0, ls, color=color, label=level)
    
    # Formatting
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.add_patch(patches.Rectangle([0, 0], 1, 1, fc='gray', transform=ax.transAxes, zorder=-10))
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    return c, variable


def plot_2Dvariable(
    frame, fig, ds_in, varname, units, scale=1, animdim=None,
    depths=None, months=None, sigmabins=None, clims=None, cmap=None,
):
    """Plot specified variable over spatial region. Hard-coded for two timeranges,
    two meshes and the residual between meshes. These categories must be consistent
    with the structure of the `plotvars` dictionary. If `animdim` is specified as
    either 'months' or 'sigmaBins', this function behaves as a frame generator for
    use with `matplotlib.animation.FuncAnimation`. Otherwise, the `frame` argument
    must still be provided, but has no effect.
    """
    
    # Parse plot arguments
    meshes, tstrs, sliceargs, clims = parse_plot_args(ds_in, depths, months, sigmabins, clims)
    
    # Animation frame labels
    if animdim == 'months':
        framelabel = datetime(1, frame+1, 1).strftime('%b')
    elif animdim == 'sigmaBins':
        value = float(ds_in[meshes[0]][tstrs[0]][animdim][frame])
        framelabel = f'$\sigma_{{\\theta}}$ = {value:0.1f} kg m$^{-3}$'
    
    # Make plot area and add titles
    fig.clf()
    xpos, ypos = [0, 0.33, 0.66], [0.56, 0.15]
    axs = np.vstack([[fig.add_subplot([x, y, 0.31, 0.39]) for x in xpos] for y in ypos])
    for col, title in zip(axs.T, meshes + [f'{meshes[0]}-{meshes[1]}']):
        col[0].set_title(title)
    
    # Loop through time ranges
    for row, tstr in zip(axs, tstrs):
        
        # Print time range
        row[0].text(0.01, 0.8, tstr, transform=row[0].transAxes)
        
        # --- Plot meshes ---------
        res = []
        for ax, mesh in zip(row, meshes):
            
            # Slice dataset
            ds = ds_in[mesh][tstr]
            if animdim is not None:
                ds = ds.isel({animdim: frame})
            for name, values in sliceargs.items():
                if name in ds.dims:
                    ds = ds.sel(**{name: values[1]})
                if name in ds.dims:
                    ds = getattr(ds, values[0])(dim=name, skipna=False)
            
            # Plot sliced dataset
            c1, variable = plot_2Dpanel(
                ax, ds, varname, scale=scale, clim=clims[0], cmap=cmap, plot_sigma=True,
            )
            res.append(variable)
        
        # --- Plot residual -------
        res = np.subtract(*res)
        c2, _ = plot_2Dpanel(row[2], res, scale=scale, clim=clims[1], cmap='RdBu_r')

    # Add colorbars
    cax1 = fig.add_axes([0.01, 0.1, 0.62, 0.015])
    cax2 = fig.add_axes([0.67, 0.1, 0.29, 0.015])
    fig.colorbar(c1, cax=cax1, label=units, orientation='horizontal')
    fig.colorbar(c2, cax=cax2, label=units, orientation='horizontal')
    for x, mesh in zip([0.9, 0.1], meshes):
        cax2.text(x, 1.5, mesh, transform=cax2.transAxes, ha='center')
    
    # Add sigma contour legend and frame label
    axs[1, 0].legend(ncols=7, loc=(0.08, -0.085))
    if animdim is not None:
        axs[0, 0].text(0.01, 0.95, framelabel, transform=axs[0, 0].transAxes)