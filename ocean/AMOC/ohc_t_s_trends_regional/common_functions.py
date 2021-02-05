"""
Functions for plotting
"""
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, others...

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import FuncFormatter, FixedLocator
import matplotlib.path
import cartopy
import cartopy.crs as ccrs
import shapely.geometry
from functools import partial
import xarray as xr
import pandas as pd
import numpy as np
import datetime
import netCDF4


def days_to_datetime(days, calendar='gregorian', referenceDate='0001-01-01'):
    """
    Convert days to ``datetime.datetime`` objects given a reference date and an
    MPAS calendar (either 'gregorian' or 'gregorian_noleap').

    Parameters
    ----------
    days : float or array-like of floats
        The number of days since the reference date.

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        A calendar to be used to convert days to a ``datetime.datetime``
        object.

    referenceDate : str, optional
        A reference date of the form::

            0001-01-01
            0001-01-01 00:00:00

    Returns
    -------
    datetime : `datetime.datetime` (or array-like of datetimes)
        The days since ``referenceDate`` on the given ``calendar``.

    Raises
    ------
    ValueError
        If an invalid ``days``, ``referenceDate`` or ``calendar`` is supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    # Convert from MPAS calendar to NetCDF4 calendar names.
    if calendar == 'gregorian_noleap':
        calendar = 'noleap'
    elif calendar != 'gregorian':
        raise ValueError('Unsupported calendar {}'.format(calendar))

    datetimes = netCDF4.num2date(days,
                                 'days since {}'.format(referenceDate),
                                 calendar=calendar)

    # convert to datetime.datetime
    if isinstance(datetimes, np.ndarray):
        newDateTimes = []
        for date in datetimes.flat:
            newDateTimes.append(_round_datetime(date))
        if len(newDateTimes) > 0:
            datetimes = np.reshape(np.array(newDateTimes),
                                      datetimes.shape)

    else:
        # Round datetimes to nearest second
        (year, month, day, hour, minute, second, microsecond) = \
            (datetimes.year, datetimes.month, datetimes.day, datetimes.hour, \
             datetimes.minute, datetimes.second, datetimes.microsecond)
        datetimes = datetime.datetime(year=year, month=month, day=day,
                                      hour=hour, minute=minute, second=second)
        add_seconds = int(1e-6 * microsecond + 0.5)
        datetimes = datetimes + datetime.timedelta(0, add_seconds)

    return datetimes


def date_to_days(year=1, month=1, day=1, hour=0, minute=0, second=0,
                 calendar='gregorian', referenceDate='0001-01-01'):
    """
    Convert a date to days since the reference date.

    Parameters
    ----------
    year, month, day, hour, minute, second : int, optional
        The date to be converted to days since ``referenceDate`` on the
        given ``calendar``.

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        A calendar to be used to convert days to a ``datetime.datetime``
        object.

    referenceDate : str, optional
        A reference date of the form::

            0001-01-01
            0001-01-01 00:00:00

    Returns
    -------
    days : float
        The days since ``referenceDate`` on the given ``calendar``.

    Raises
    ------
    ValueError
        If an invalid ``referenceDate`` or ``calendar`` is supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    date = datetime.datetime(year, month, day, hour, minute, second)

    # Convert from MPAS calendar to NetCDF4 calendar names.
    if calendar == 'gregorian_noleap':
        calendar = 'noleap'
    elif calendar != 'gregorian':
        raise ValueError('Unsupported calendar {}'.format(calendar))

    return netCDF4.date2num(date, 'days since {}'.format(referenceDate),
                            calendar=calendar)


def plot_xtick_format(calendar, minDays, maxDays, maxXTicks, yearStride=None):
    '''
    Formats tick labels and positions along the x-axis for time series plots

    Parameters
    ----------
    calendar : str
        the calendar to use for formatting the time axis

    minDays : float
        start time for labels

    maxDays : float
        end time for labels

    maxXTicks : int
        the maximum number of tick marks to display, used to sub-sample ticks
        if there are too many

    yearStride : int, optional
        the number of years to skip over between ticks
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def date_tick(days, pos, calendar='gregorian', includeMonth=True):
        days = np.maximum(days, 0.)
        date = days_to_datetime(days, calendar)
        if includeMonth:
            return '{:04d}-{:02d}'.format(date.year, date.month)
        else:
            return '{:04d}'.format(date.year)

    ax = plt.gca()

    start = days_to_datetime(np.amin(minDays), calendar=calendar)
    end = days_to_datetime(np.amax(maxDays), calendar=calendar)

    if yearStride is not None or end.year - start.year > maxXTicks / 2:
        if yearStride is None:
            yearStride = 1
        else:
            maxXTicks = None
        major = [date_to_days(year=year, calendar=calendar)
                 for year in np.arange(start.year, end.year + 1, yearStride)]
        formatterFun = partial(date_tick, calendar=calendar,
                               includeMonth=False)
    else:
        # add ticks for months
        major = []
        for year in range(start.year, end.year + 1):
            for month in range(1, 13):
                major.append(date_to_days(year=year, month=month,
                                          calendar=calendar))
        formatterFun = partial(date_tick, calendar=calendar,
                               includeMonth=True)

    ax.xaxis.set_major_locator(FixedLocator(major, maxXTicks))
    ax.xaxis.set_major_formatter(FuncFormatter(formatterFun))

    plt.setp(ax.get_xticklabels(), rotation=30)

    plt.autoscale(enable=True, axis='x', tight=True)


def subdivide_geom(geometry, geomtype, maxLength):
    '''
    Subdivide the line segments for a given set of geometry so plots are
    smoother
    '''
    # Authors
    # -------
    # Xylar Asay-Davis, Phillip J. Wolfram

    def subdivide_line_string(lineString, periodic=False):
        coords = list(lineString.coords)
        if periodic:
            # add periodic last entry
            coords.append(coords[0])

        outCoords = [coords[0]]
        for iVert in range(len(coords)-1):
            segment = shapely.geometry.LineString([coords[iVert],
                                                   coords[iVert+1]])
            if(segment.length < maxLength):
                outCoords.append(coords[iVert+1])
            else:
                # we need to subdivide this segment
                subsegment_count = int(np.ceil(segment.length/maxLength))
                for index in range(subsegment_count):
                    point = segment.interpolate(
                        float(index+1)/float(subsegment_count),
                        normalized=True)
                    outCoords.append(point.coords[0])

        if periodic:
            # remove the last entry
            outCoords.pop()
        return outCoords

    if geomtype == 'LineString':
        newGeometry = shapely.geometry.LineString(
            subdivide_line_string(geometry))
    elif geomtype == 'MultiLineString':
        outStrings = [subdivide_line_string(inLineString) for inLineString
                      in geometry]
        newGeometry = shapely.geometry.MultiLineString(outStrings)
    elif geomtype == 'Polygon':
        exterior = subdivide_line_string(geometry.exterior, periodic=True)
        interiors = [subdivide_line_string(inLineString, periodic=True)
                     for inLineString in geometry.interiors]
        newGeometry = shapely.geometry.Polygon(exterior, interiors)
    elif geomtype == 'MultiPolygon':
        polygons = []
        for polygon in geometry:
            exterior = subdivide_line_string(polygon.exterior, periodic=True)
            interiors = [subdivide_line_string(inLineString, periodic=True)
                         for inLineString in polygon.interiors]
            polygons.append((exterior, interiors))

        newGeometry = shapely.geometry.MultiPolygon(polygons)
    elif geomtype == 'Point':
        newGeometry = geometry
    else:
        print("Warning: subdividing geometry type {} is not supported.".format(
            geomtype))
        newGeometry = geometry

    return newGeometry


def timeseries_analysis_plot(dsvalues, N, title, xlabel, ylabel,
                             calendar, lineColors=None,
                             lineStyles=None, markers=None, lineWidths=None,
                             legendText=None, maxPoints=None,
                             titleFontSize=None, figsize=(15, 6), dpi=None,
                             firstYearXTicks=None, yearStrideXTicks=None,
                             maxXTicks=20, obsMean=None, obsUncertainty=None,
                             obsLegend=None, legendLocation='lower left'):
    """
    Plots the list of time series data sets.

    Parameters
    ----------
    dsvalues : list of xarray DataSets
        the data set(s) to be plotted

    N : int
        the numer of time points over which to perform a moving average

    title : str
        the title of the plot

    xlabel, ylabel : str
        axis labels

    calendar : str
        the calendar to use for formatting the time axis

    lineColors, lineStyles, markers, legendText : list of str, optional
        control line color, style, marker, and corresponding legend
        text.  Default is black, solid line with no marker, and no legend.

    lineWidths : list of float, optional
        control line width.  Default is 1.0.

    maxPoints : list of {None, int}, optional
        the approximate maximum number of time points to use in a time series.
        This can be helpful for reducing the number of symbols plotted if
        plotting with markers.  Otherwise the markers become indistinguishable
        from each other.

    titleFontSize : int, optional
        the size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    firstYearXTicks : int, optional
        The year of the first tick on the x axis.  By default, the first time
        entry is the first tick.

    yearStrideXTicks : int, optional
        The number of years between x ticks. By default, the stride is chosen
        automatically to have ``maxXTicks`` tick marks or fewer.

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the x axis.
        This may need to be adjusted depending on the figure size and aspect
        ratio.

    obsMean, obsUncertainty : list of float, optional
        Mean values and uncertainties for observations to be plotted as error
        bars. The two lists must have the same number of elements.

    obsLegend : list of str, optional
        The label in the legend for each element in ``obsMean`` (and
        ``obsUncertainty``)

    legendLocation : str, optional
        The location of the legend (see ``pyplot.legend()`` for details)

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The resulting figure
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Stephen Price

    if dpi is None:
        dpi = 200
    axis_font = {'size': 16}
    if titleFontSize is None:
        titleFontSize = 20
    title_font = {'size': titleFontSize,
                  'color': 'black',
                  'weight': 'normal'}

    fig = plt.figure(figsize=figsize, dpi=dpi)

    minDays = []
    maxDays = []
    labelCount = 0
    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        if N == 1 or N is None:
            mean = dsvalue
        else:
            mean = pd.Series.rolling(dsvalue.to_pandas(), N,
                                     center=True).mean()
            mean = xr.DataArray.from_series(mean)
        minDays.append(mean.Time.min())
        maxDays.append(mean.Time.max())

        if maxPoints is not None and maxPoints[dsIndex] is not None:
            nTime = mean.sizes['Time']
            if maxPoints[dsIndex] < nTime:
                stride = int(round(nTime / float(maxPoints[dsIndex])))
                mean = mean.isel(Time=slice(0, None, stride))

        if legendText is None:
            label = None
        else:
            label = legendText[dsIndex]
            labelCount += 1
        if lineColors is None:
            color = 'k'
        else:
            color = lineColors[dsIndex]
        if lineStyles is None:
            linestyle = '-'
        else:
            linestyle = lineStyles[dsIndex]
        if markers is None:
            marker = None
        else:
            marker = markers[dsIndex]
        if lineWidths is None:
            linewidth = 1.
        else:
            linewidth = lineWidths[dsIndex]

        plt.plot(mean['Time'].values, mean.values, color=color,
                 linestyle=linestyle, marker=marker, linewidth=linewidth,
                 label=label)

    if obsMean is not None:
        obsCount = len(obsMean)
        assert(len(obsUncertainty) == obsCount)

        # space the observations along the time line, leaving gaps at either
        # end
        start = np.amin(minDays)
        end = np.amax(maxDays)
        obsTimes = np.linspace(start, end, obsCount + 2)[1:-1]
        obsSymbols = ['o', '^', 's', 'D', '*']
        obsColors = ['b', 'g', 'c', 'm', 'r']
        for iObs in range(obsCount):
            if obsMean[iObs] is not None:
                symbol = obsSymbols[np.mod(iObs, len(obsSymbols))]
                color = obsColors[np.mod(iObs, len(obsColors))]
                fmt = '{}{}'.format(color, symbol)
                plt.errorbar(obsTimes[iObs],
                             obsMean[iObs],
                             yerr=obsUncertainty[iObs],
                             fmt=fmt,
                             ecolor=color,
                             capsize=0,
                             label=obsLegend[iObs])
                # plot a box around the error bar to make it more visible
                boxHalfWidth = 0.01 * (end - start)
                boxHalfHeight = obsUncertainty[iObs]
                boxX = obsTimes[iObs] + \
                    boxHalfWidth * np.array([-1, 1, 1, -1, -1])
                boxY = obsMean[iObs] + \
                    boxHalfHeight * np.array([-1, -1, 1, 1, -1])

                plt.plot(boxX, boxY, '{}-'.format(color), linewidth=3)
                labelCount += 1

    if labelCount > 1:
        plt.legend(loc=legendLocation)

    ax = plt.gca()

    if firstYearXTicks is not None:
        minDays = date_to_days(year=firstYearXTicks, calendar=calendar)

    plot_xtick_format(calendar, minDays, maxDays, maxXTicks,
                      yearStride=yearStrideXTicks)

    # Add a y=0 line if y ranges between positive and negative values
    yaxLimits = ax.get_ylim()
    if yaxLimits[0] * yaxLimits[1] < 0:
        x = ax.get_xlim()
        plt.plot(x, np.zeros(np.size(x)), 'k-', linewidth=1.2, zorder=1)

    if title is not None:
        plt.title(title, **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    return fig


def timeseries_analysis_plot_polar(dsvalues, N, title,
                                   lineColors=None, lineStyles=None,
                                   markers=None, lineWidths=None,
                                   legendText=None, titleFontSize=None,
                                   figsize=(15, 6), dpi=None):
    """
    Plots the list of time series data sets on a polar plot.

    Parameters
    ----------
    dsvalues : list of xarray DataSets
        the data set(s) to be plotted

    N : int
        the numer of time points over which to perform a moving average

    title : str
        the title of the plot

    lineColors, lineStyles, markers, legendText : list of str, optional
        control line color, style, marker, and corresponding legend
        text.  Default is black, solid line with no marker, and no legend.

    lineWidths : list of float, optional
        control line width.  Default is 1.0.

    titleFontSize : int, optional
        the size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The resulting figure
    """
    # Authors
    # -------
    # Adrian K. Turner, Xylar Asay-Davis

    daysInMonth = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    abrevMonthNames = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
                       "Sep", "Oct", "Nov", "Dec"]

    if dpi is None:
        dpi = 200
    if titleFontSize is None:
        titleFontSize = 20
    title_font = {'size': titleFontSize,
                  'color': black,
                  'weight': normal}

    fig = plt.figure(figsize=figsize, dpi=dpi)

    minDays = []
    maxDays = []
    labelCount = 0
    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        mean = pd.Series.rolling(dsvalue.to_pandas(), N, center=True).mean()
        mean = xr.DataArray.from_series(mean)
        minDays.append(mean.Time.min())
        maxDays.append(mean.Time.max())

        if legendText is None:
            label = None
        else:
            label = legendText[dsIndex]
            labelCount += 1
        if lineColors is None:
            color = 'k'
        else:
            color = lineColors[dsIndex]
        if lineStyles is None:
            linestyle = '-'
        else:
            linestyle = lineStyles[dsIndex]
        if markers is None:
            marker = None
        else:
            marker = markers[dsIndex]
        if lineWidths is None:
            linewidth = 1.
        else:
            linewidth = lineWidths[dsIndex]

        plt.polar((mean['Time'] / 365.0) * np.pi * 2.0, mean, color=color,
                  linestyle=linestyle, marker=marker, linewidth=linewidth,
                  label=label)

    if labelCount > 1:
        plt.legend(loc='lower left')

    ax = plt.gca()

    # set azimuthal axis formatting
    majorTickLocs = np.zeros(12)
    minorTickLocs = np.zeros(12)
    majorTickLocs[0] = 0.0
    minorTickLocs[0] = (daysInMonth[0] * np.pi) / 365.0
    for month in range(1, 12):
        majorTickLocs[month] = majorTickLocs[month - 1] + \
            ((daysInMonth[month - 1] * np.pi * 2.0) / 365.0)
        minorTickLocs[month] = minorTickLocs[month - 1] + \
            (((daysInMonth[month - 1] +
               daysInMonth[month]) * np.pi) / 365.0)

    ax.set_xticks(majorTickLocs)
    ax.set_xticklabels([])

    ax.set_xticks(minorTickLocs, minor=True)
    ax.set_xticklabels(abrevMonthNames, minor=True)

    if title is not None:
        plt.title(title, **title_font)

    return fig


def hovmoeller_plot(Time, z, field, colormap, cnorm, clevels,
                    title, xlabel, ylabel, calendar, colorbarLabel=None,
                    titleFontSize=None, figsize=(15, 6), dpi=None,
                    firstYearXTicks=None, yearStrideXTicks=None, maxXTicks=20):
    """
    Plots Hovmoeller graph

    Parameters
    ----------
    Time, z, field : numpy arrays

    title : str
        the title of the plot

    xlabel, ylabel : str
        axis labels

    calendar : str
        the calendar to use for formatting the time axis

    titleFontSize : int, optional
        the size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    firstYearXTicks : int, optional
        The year of the first tick on the x axis.  By default, the first time
        entry is the first tick.

    yearStrideXTicks : int, optional
        The number of years between x ticks. By default, the stride is chosen
        automatically to have ``maxXTicks`` tick marks or fewer.

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the x axis.
        This may need to be adjusted depending on the figure size and aspect
        ratio.

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The resulting figure
    """
    # Authors
    # -------
    # Milena Veneziani

    if dpi is None:
        dpi = 200
    axis_font = {'size': 16}
    if titleFontSize is None:
        titleFontSize = 20
    title_font = {'size': titleFontSize,
                  'color': 'black',
                  'weight': 'normal'}

    fig = plt.figure(figsize=figsize, dpi=dpi)

    [x ,y] = np.meshgrid(Time, z)

    cf = plt.contourf(x, y, field, cmap=colormap, norm=cnorm,
                      levels=clevels, extend='both')
    cbar = plt.colorbar(cf, orientation='vertical', spacing='uniform',
                        aspect=9, ticks=clevels)
#                        aspect=9, ticks=clevels, boundaries=clevels)
    if colorbarLabel is not None:
        cbar.set_label(colorbarLabel)

    minDays = np.min(Time)
    maxDays = np.max(Time)

    if firstYearXTicks is not None:
        minDays = date_to_days(year=firstYearXTicks, calendar=calendar)

    plot_xtick_format(calendar, minDays, maxDays, maxXTicks,
                      yearStride=yearStrideXTicks)

    if title is not None:
        plt.title(title, **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    return fig

def add_inset(fig, fc, latlonbuffer=45., polarbuffer=5., width=1.0,
              height=1.0, lowerleft=None, xbuffer=None, ybuffer=None,
              maxlength=1.):
    '''
    Plots an inset map showing the location of a transect or polygon.  Shapes
    are plotted on a polar grid if they are entirely poleward of +/-50 deg.
    latitude and with a lat/lon grid if not.

    Parameters
    ----------
    fig : ``matplotlib.figure.Figure``
        A matplotlib figure to add the inset to

    fc : ``geometric_features.FeatureCollection``
        A collection of regions, transects and/or points to plot in the inset

    latlonbuffer : float, optional
        The number of degrees lat/lon to use as a buffer around the shape(s)
        to plot if a lat/lon plot is used.

    polarbuffer : float, optional
        The number of degrees latitude to use as a buffer equatorward of the
        shape(s) in polar plots

    width, height : float, optional
        width and height in inches of the inset

    lowerleft : pair of floats, optional
        the location of the lower left corner of the axis in inches, default
        puts the inset in the upper right corner of ``fig``.

    xbuffer, ybuffer : float, optional
        right and top buffers from the top-right corner (in inches) if
        lowerleft is ``None``.

    maxlength : float or ``None``, optional
        Any segments longer than maxlength will be subdivided in the plot to
        ensure curvature.  If ``None``, no subdivision is performed.

    Returns
    -------
    inset : ``matplotlib.axes.Axes``
        The new inset axis
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def set_circular_boundary(ax):
        '''Set the boundary of the given axis to be circular (for a polar plot)'''

        # Compute a circle in axes coordinates, which we can use as a boundary
        # for the map. We can pan/zoom as much as we like - the boundary will be
        # permanently circular.
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = matplotlib.path.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)

    def get_bounds(fc):
        '''Compute the lon/lat bounding box for all transects and regions'''

        bounds = shapely.geometry.GeometryCollection()
        for feature in fc.features:
            shape = shapely.geometry.shape(feature['geometry'])
            shape_bounds = shapely.geometry.box(*shape.bounds)
            bounds = shapely.geometry.box(*(bounds.union(shape_bounds).bounds))
        return bounds.bounds

    minLon, minLat, maxLon, maxLat = get_bounds(fc)

    figsize = fig.get_size_inches()
    width /= figsize[0]
    height /= figsize[1]
    if lowerleft is None:
        if xbuffer is None:
            xbuffer = 0.1*width
        else:
            xbuffer /= figsize[0]
        if ybuffer is None:
            ybuffer = xbuffer*figsize[0]/figsize[1]
        else:
            ybuffer /= figsize[1]
        lowerleft = [1.0 - width - xbuffer, 1.0 - height - ybuffer]
    else:
        lowerleft = [lowerleft[0]/figsize[0], lowerleft[1]/figsize[1]]
    bounds = [lowerleft[0], lowerleft[1], width, height]

    if maxLat <= -50:
        # an Antarctic-focused map makes the most sense
        inset = fig.add_axes(bounds,
                             projection=ccrs.SouthPolarStereo())
        extent = [-180., 180., -90., max(-65., maxLat+polarbuffer)]
        set_circular_boundary(inset)
        xlocator = mticker.FixedLocator(np.linspace(-180., 180., 9))
        ylocator = mticker.FixedLocator(np.linspace(-90., -50., 9))
    elif minLat >= 50:
        # an Arctic-focused map makes the most sense
        inset = fig.add_axes(bounds,
                             projection=ccrs.NorthPolarStereo())
        extent = [-180, 180, min(65., minLat-polarbuffer), 90]
        set_circular_boundary(inset)
        xlocator = mticker.FixedLocator(np.linspace(-180., 180., 9))
        ylocator = mticker.FixedLocator(np.linspace(50., 90., 9))
    else:
        inset = fig.add_axes(bounds,
                             projection=ccrs.PlateCarree())
        extent = [max(-180., minLon-latlonbuffer),
                  min(180., maxLon+latlonbuffer),
                  max(-90., minLat-latlonbuffer),
                  min(90., maxLat+latlonbuffer)]
        xlocator = None
        ylocator = None

    # kind of like "top" justified -- graphics are toward the "north" end of
    # the subplot
    inset.set_anchor('N')

    inset.set_extent(extent,  ccrs.PlateCarree())
    inset.add_feature(cartopy.feature.LAND, zorder=1)
    inset.add_feature(cartopy.feature.OCEAN, zorder=0)

    gl = inset.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                         linewidth=0.5, color='gray', alpha=0.5,
                         linestyle='--')

    if xlocator is not None:
        gl.xlocator = xlocator

    if ylocator is not None:
        gl.ylocator = ylocator

    for feature in fc.features:
        geomtype = feature['geometry']['type']
        shape = shapely.geometry.shape(feature['geometry'])
        if maxlength is not None:
            shape = subdivide_geom(shape, shape.geom_type, maxlength)
        if geomtype in ['Polygon', 'MultiPolygon']:
            inset.add_geometries((shape,), crs=ccrs.PlateCarree(),
                                 edgecolor='blue', facecolor='blue', alpha=0.4,
                                 linewidth=1.)
        elif geomtype in ['Point', 'MultiPoint']:
            inset.add_geometries((shape,), crs=ccrs.PlateCarree(),
                                 edgecolor='none', facecolor='none', alpha=1.,
                                 markersize=3., markeredgecolor='k',
                                 markerfacecolor='k')
        else:
            inset.add_geometries((shape,), crs=ccrs.PlateCarree(),
                                 edgecolor='k', facecolor='none', alpha=1.,
                                 linewidth=1.)
            # put a red point at the beginning and a blue point at the end
            # of the transect to help show the orientation
            begin = shape.coords[0]
            end = shape.coords[-1]
            inset.plot(begin[0], begin[1], color='r', marker='o',
                       markersize=3., transform=ccrs.PlateCarree())
            inset.plot(end[0], end[1], color='g', marker='o',
                       markersize=3., transform=ccrs.PlateCarree())

    return inset
