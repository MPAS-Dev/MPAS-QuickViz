#!/usr/bin/env python
"""
    Name: visualizationtools.py
    Author: Ben Moore-Maley (bmoorema@lanl.gov)

    Visualization tools for the ImPACTS Water Mass Analysis project
"""

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt, patches, dates
from datetime import datetime


def slice_timeseries(da, region=None, sigmabins=None):
    """Slice time series `xr.DataArray` by either region or sigmabins
    """
    
    # Slice by region or sigmabins
    if 'regionNames' in da.dims:
        da = da.sel(regionNames=bytes(region, 'utf-8'))
        if 'depths' in da.dims: da = da.sel(depths=0)
    elif 'sigmaBins' in da.dims:
        da = da.sel(sigmaBins=slice(*sigmabins)).sum(dim='sigmaBins')
    
    return da


def build_lag_slices(lag):
    """Build slices to stagger arrays for a given `lag`
    """
    
    # Build slices
    slice2 = slice(abs(lag), -(abs(lag) + 1))
    slice1 = slice(slice2.start + lag, slice2.stop + lag)
    
    return [slice1, slice2]


def calc_lag_correlation(variables, lags, cutoff=None):
    """Calculate Pearson coefficient for each lag in `lags`
    and return as `np.array`
    """

    # Loop through lags and calculate Pearson coefficient
    r = []
    for lag in lags:
        slices = build_lag_slices(lag)
        x, y = [var[slc].values for var, slc in zip(variables, slices)]
        if cutoff is not None:
            index = y > cutoff
            x, y = x[index], y[index]
        r.append(np.corrcoef(x, y)[0, 1])
    r = np.array(r)

    return r


def plot_wmt_sigma(ds, varname, timeranges, months=None, ylims=None):
    """Plot water mass transformation and formation in sigma space
    time-averaged over each time range in `timeranges`.
    `varname` specifies the grouping (heat, salt or total) and
    `months` specifies which months to include in the time average.
    """
    
    # Default months
    if months is None:
        months = list(range(1, 13))
    if ylims is None:
        ylims = [(None, None), (None, None)]
    
    # Make figure
    fig, axs = plt.subplots(2, 1, figsize=(10, 5))
    palette = plt.get_cmap('tab10')

    # Loop through categories
    for ax, ctgy, title, ylim in zip(axs, ['Trans', 'Form'], ['Transformation', 'Formation'], ylims):

        # Loop through time ranges
        for timerange, ls in zip(timeranges, [':', '-']):
            tstr = '-'.join(str(t.year) for t in timerange)
            
            # Loop through meshes -> slice and time average -> plot
            for mesh, c in zip(ds, palette([0, 1])):
                da = ds[mesh][varname + ctgy].sel(time=slice(*timerange))
                da = da.sel(time=da.time.dt.month.isin(months)).mean(dim='time')
                ax.plot(da.sigmaBins, da, ls, color=c, label=f'{mesh} {tstr}')

        # Format panel
        ax.plot([21, 29], [0, 0], 'k--')
        ax.set_xlim([21.5, 28.3])
        ax.set_ylim(ylim)
        ax.set_ylabel('Sv')
        ax.set_title(title)

    # Final formatting
    axs[0].xaxis.set_ticklabels('')
    axs[1].set_xlabel('$\sigma_{\\theta}$ [kg m$^{-3}$]')
    axs[0].legend(loc=2)

    return fig, axs


def plot_wmt_timeseries(ds, varname, timeranges, sigmarange, ylim=[-100, 150]):
    """Plot formation time series for the formation component
    given by `varname` on the time ranges given by `timeranges`.
    `sigmarange` designates the range of sigma values to integrate
    over, determining the type of mode water in question.
    """

    # Make figure
    meshes = list(ds)
    fig, axs = plt.subplots(2, 2, figsize=(10, 5), gridspec_kw={'wspace': 0.05})
    palette = plt.get_cmap('BrBG').resampled(4)
    fcs, ecs = palette([1, 2]), palette([0, 3])

    # Loop through meshes
    for row, mesh in zip(axs, meshes):

        # Initial formatting
        row[0].set_ylabel('Sv')
        row[0].text(0.01, 0.03, mesh, transform=row[0].transAxes)
        row[1].yaxis.set_ticklabels('')

        # Loop through time ranges
        for ax, timerange in zip(row, timeranges):

            # Panel formatting
            ax.plot(timerange, [0, 0], 'k--')
            ax.set_xlim(timerange)
            ax.set_ylim(ylim)
            ax.xaxis.set_major_locator(dates.YearLocator(base=2))
            ax.xaxis.set_major_formatter(dates.DateFormatter('%Y'))
            if mesh == meshes[0]:
                ax.xaxis.set_ticklabels('')

            # Slice variable by time range and sigmarange
            da = ds[mesh][varname].sel(time=slice(*timerange)).sel(sigmaBins=slice(*sigmarange))

            # Plot positive (formation) and negative (destruction) components
            for label, func, fc, ec in zip(['Formation', 'Destruction'], ['greater', 'less'], fcs, ecs):
                da_plot = da.where(getattr(np, func)(da, 0)).sum(dim='sigmaBins')
                ax.fill_between(da_plot.time, da_plot, 0, fc=fc, ec=ec, alpha=0.6, label=label)

            # Plot total (net) formation
            da_plot = da.sum(dim='sigmaBins')
            ax.plot(da_plot.time, da_plot, color='r', label='Net')

    axs[1, 0].legend(ncols=3, loc=9)
    
    return fig, axs


def plot_wmt_hoffmueller(ds, varname, timeranges, xlim=[25, 28], clim=75):
    """Plot Hoffmueller time vs sigma figure for the WMT variable
    given by 'varname' over the time ranges in `timeranges`.
    """

    # Make figure and add xlabel and titles
    meshes = list(ds)
    fig, axs = plt.subplots(2, 3, figsize=(10, 10), gridspec_kw={'hspace': 0.05, 'wspace': 0.1})
    axs[1, 1].set_xlabel('$\sigma_{\\theta}$ [kg m$^{-3}$]')
    for col, title in zip(axs.T, meshes + [f'{meshes[0]}-{meshes[1]}']):
        col[0].set_title(title)

    # Loop through time ranges
    for row, timerange in zip(axs, timeranges[::-1]):

        # Plot variable on time range for each mesh
        res = []
        for ax, mesh in zip(row, meshes):
            da = ds[mesh][varname].sel(time=slice(*timerange))
            c1 = ax.pcolormesh(da.sigmaBins, da.time, da, cmap='BrBG_r', vmin=-clim, vmax=clim)
            res.append(da)

        # Plot residual between meshes
        res = np.subtract(*res)
        c2 = row[2].pcolormesh(res.sigmaBins, res.time, res, cmap='RdBu_r', vmin=-clim, vmax=clim)

        # Formatting
        for n, ax in enumerate(row):
            ax.set_xlim(xlim)
            if n !=0:
                ax.yaxis.set_ticklabels('')
            if timerange == timeranges[1]:
                ax.xaxis.set_ticklabels('')

    # Add colorbars
    cax1 = fig.add_axes([0.13, 0.04, 0.5, 0.01])
    cax2 = fig.add_axes([0.66, 0.04, 0.23, 0.01])
    fig.colorbar(c1, cax=cax1, label='Sv', orientation='horizontal')
    fig.colorbar(c2, cax=cax2, label='Sv', orientation='horizontal')
    for x, mesh in zip([0.9, 0.1], meshes):
        cax2.text(x, 1.5, mesh, transform=cax2.transAxes, ha='center')
    
    return fig, axs


def plot_lag_correlation(ax, ds, varnames, sigmabins, regions=None, cutoff=None, laglim=6):
    """Plot the correlation vs lag as defined by the Pearson coefficient for
    the variables given by `varnames`. `cutoff` defines the minimum value to
    included for the dependent variable.
    """
    
    # Define lags, meshes and default regions
    lags = list(range(-laglim, laglim+1))
    meshes = list(ds)
    if regions is None:
        regions = ds[meshes[0]].regionNames.values.astype(str)

    # Figure formatting
    palette = plt.get_cmap('tab10')
    xlim, ylim = [-laglim, laglim], [-0.8, 1]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.plot(xlim, [0, 0], 'k--')
    ax.plot([0, 0], ylim, 'k--')
    ax.set_xlabel('Lag [months]')
    ax.set_ylabel('Pearson coefficient')

    # Loop through meshes
    for mesh, c in zip(meshes, palette([0, 1])):
        
        # Slice variable 2
        var2 = slice_timeseries(ds[mesh][varnames[1]], sigmabins=sigmabins)

        # Loop through regions
        r = []
        for region in regions:

            # Slice variable 1 and append Pearson coefficient
            var1 = slice_timeseries(ds[mesh][varnames[0]], region=region)
            r.append(calc_lag_correlation([var1, var2], lags, cutoff=cutoff))

        # Calculate region stats
        r = np.vstack(r)
        r_mean = r.mean(axis=0)
        
        # Plot errorbar if multiple regions
        if len(regions) == 1:
            label, alpha = regions[0], 0.3
        else:
            r_err = [abs(getattr(r, func)(axis=0) - r_mean) for func in ('min', 'max')]
            ax.errorbar(lags, r_mean, yerr=r_err, color=c, capsize=2)
            label, alpha = mesh, 1
        
        # Plot mean
        ax.plot(lags, r_mean, '-o', color=c, label=label, alpha=alpha)


def plot_compare_variables(ax, ds, varnames, region=None, sigmabins=None, lims=None, lag=0):
    """Compare MPAS timeseries variables given by `varnames` on x, y.
    """
    
    # Define slices to account for lag, and add title
    slices = build_lag_slices(lag)

    # Loop through meshes
    for mesh in ds:
        
        # Slice variables
        variables = []
        for varname, slc in zip(varnames, slices):
            da = slice_timeseries(ds[mesh][varname], region=region, sigmabins=sigmabins)
            variables.append(da[slc])
        
        # Plot
        ax.plot(*variables, 'o', label=mesh)
    
    # Formatting
    lims = [None, None, None, None] if lims is None else lims
    ax.set_xlim(lims[:2])
    ax.set_ylim(lims[2:])
    ax.set_title(region)


def plot_correlation_multipanel(ds, varnames, sigmabins, regionlist, xlabel, lims=None, lag=0, cutoff=0):
    """Wrapper function around `plot_lag_correlation` and `plot_compare_variables`
    for a multipanel correlation summary between the variables given by `varnames`.
    """

    # Make plot
    fig = plt.figure(figsize=(10, 8))
    gs = plt.GridSpec(2, 3, hspace=0.3, wspace=0.08)
    props = {'fc': 'w', 'boxstyle': 'round', 'alpha': 0.7}

    # Top panel
    ax = fig.add_subplot(gs[0, :])
    for regions in regionlist:
        plot_lag_correlation(ax, ds, varnames, sigmabins, regions, cutoff)
    ax.legend()

    # Bottom panels
    panelregions = ['Labrador Sea', 'Irminger Sea', 'Greenland Sea']
    axs = [fig.add_subplot(gs[1, n]) for n in range(3)]
    for ax, region in zip(axs, panelregions):
        plot_compare_variables(ax, ds, varnames, region, sigmabins, lims, lag)
        if region != panelregions[0]: ax.yaxis.set_ticklabels('')
    axs[1].set_xlabel(xlabel)
    axs[0].set_ylabel('SPMW Formation [Sv]')
    axs[0].text(0.02, 0.03, f'Lag = {lag} months', transform=axs[0].transAxes, bbox=props)


def parse_2Dplot_args(ds, depths, months, sigmabins, clims):
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
    meshes, tstrs, sliceargs, clims = parse_2Dplot_args(ds_in, depths, months, sigmabins, clims)
    
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