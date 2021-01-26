from __future__ import absolute_import, division, print_function, \
    unicode_literals
import os
from netCDF4 import Dataset as netcdf_dataset
import numpy as np
import numpy.ma as ma
import xarray as xr
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import matplotlib.animation as animation
from matplotlib.pyplot import cm
from matplotlib.colors import from_levels_and_colors
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import cmocean

def _add_land_lakes_coastline(ax):
    land_50m = cfeature.NaturalEarthFeature(
            'physical', 'land', '50m', edgecolor='face',
            facecolor='lightgray', linewidth=0.5)
    lakes_50m = cfeature.NaturalEarthFeature(
            'physical', 'lakes', '50m', edgecolor='k',
            facecolor='aliceblue',
            linewidth=0.5)
    coast_50m = cfeature.NaturalEarthFeature(
            'physical', 'coastline', '50m', edgecolor='k',
            facecolor='None', linewidth=0.5)
    ax.add_feature(land_50m, zorder=2)
    ax.add_feature(lakes_50m, zorder=3)
    ax.add_feature(coast_50m, zorder=4)


meshfile = '/compyfs/inputdata/ocn/mpas-o/EC30to60E2r2/ocean.EC30to60E2r2.200908.nc'
#
runname = '20201030.alpha5_v1p-1_target.piControl.ne30pg2_r05_EC30to60E2r2-1900_ICG.compy'
modeldir = '/compyfs/malt823/E3SM_simulations/{}/archive/ocn/hist'.format(runname)
#runname = '20201022.alpha5_v1p-1_fallback.piControl.ne30pg2_r05_EC30to60E2r2-1900_ICG.compy'
#modeldir = '/compyfs/gola749/E3SM_simulations/{}/archive/ocn/hist'.format(runname)
#
#infiles = sorted(glob.glob('{}/{}.mpaso.hist.am.timeSeriesStatsMonthly.00*'.format(modeldir, runname)))[0]
#infiles = sorted(glob.glob('{}/{}.mpaso.hist.am.timeSeriesStatsMonthly.00*'.format(modeldir, runname)))[0:600]
infiles = sorted(glob.glob('{}/{}.mpaso.hist.am.timeSeriesStatsMonthlyMax.00*'.format(modeldir, runname)))[0:600]
#infiles = sorted(glob.glob('{}/{}.mpaso.hist.am.timeSeriesStatsMonthly.00*'.format(modeldir, runname)))
print('\ninfiles={}\n'.format(infiles))

#variable = 'temperatureForcingMLTend'
variable = 'mld'
variable = 'maxmld'
#variable = 'temperatureSurfaceFluxTendency'
#variable = 'temperatureShortWaveTendency'
#variable = 'temperatureHorizontalAdvectionTendency'
#variable = 'temperatureVerticalAdvectionTendency'
#variable = 'temperatureHorMixTendency'
#variable = 'temperatureVertMixTendency'
#variable = 'temperatureNonLocalTendency'
#variable = 'temperatureTotalAdvectionTendency' # derived variable
#variable = 'temperatureForcingTendency' # derived variable
#variable = 'temperatureSumTendencyTerms' # derived variable
#variable = 'temperatureTendency' # derived variable

figdir = './animations_xymaps'
if not os.path.isdir(figdir):
    os.makedirs(figdir)

figsize = [16, 12]
figdpi = 100

# z levels [m] (relevant for 3d variables)
#dlevels = [50.0, 100.0, 250.0, 500.0, 3000.0]
dlevels = [0.]

colorIndices0 = [0, 10, 28, 57, 85, 113, 142, 170, 198, 227, 242, 255]

variables = [{'name': 'temperatureSurfaceFluxTendency',
              'title': 'Surface flux tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': 'timeMonthly_avg_activeTracerSurfaceFluxTendency_temperatureSurfaceFluxTendency',
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureShortWaveTendency',
              'title': 'Penetrating shortwave flux tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': 'timeMonthly_avg_temperatureShortWaveTendency',
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureForcingTendency',
              'title': 'Total forcing tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': None,
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureHorizontalAdvectionTendency',
              'title': 'Horizontal advection tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': 'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_temperatureHorizontalAdvectionTendency',
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureVerticalAdvectionTendency',
              'title': 'Vertical advection tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': 'timeMonthly_avg_activeTracerVerticalAdvectionTendency_temperatureVerticalAdvectionTendency',
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureTotalAdvectionTendency',
              'title': 'Total advection tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': None,
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureHorMixTendency',
              'title': 'Horizontal mixing tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': 'timeMonthly_avg_activeTracerHorMixTendency_temperatureHorMixTendency',
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureVertMixTendency',
              'title': 'Vertical mixing tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': 'timeMonthly_avg_activeTracerVertMixTendency_temperatureVertMixTendency',
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureNonLocalTendency',
              'title': 'Non-local kpp flux tendency for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': 'timeMonthly_avg_activeTracerNonLocalTendency_temperatureNonLocalTendency',
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureSumTendencyTerms',
              'title': 'Sum of all tendency terms for temperature',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': None,
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureTendency',
              'title': 'Temperature tendency (derived)',
              'units': '$^\circ$C/s (x1e-6)',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'factor': 1e6,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 4],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'temperatureForcingMLTend',
              'title': 'Mixed Layer avg forcing tendency for temperature',
              'units': '$^\circ$C/s',
              'mpas': 'timeMonthly_avg_activeTracerForcingMLTend_temperatureForcingMLTend',
              'factor': 1,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4e-6, -3e-6, -2e-6, -1e-6, -0.5e-6, 0.0, 0.5e-6, 1e-6, 2e-6, 3e-6, 4e-6],
              'plot_anomalies': False,
              'is3d': False},
             {'name': 'temperatureHorAdvectionMLTend',
              'title': 'Mixed Layer avg Hadv tendency for temperature',
              'units': '$^\circ$C/s',
              'mpas': 'timeMonthly_avg_activeTracerHorAdvectionMLTend_temperatureHorAdvectionMLTend',
              'factor': 1,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4e-6, -3e-6, -2e-6, -1e-6, -0.5e-6, 0.0, 0.5e-6, 1e-6, 2e-6, 3e-6, 4e-6],
              'plot_anomalies': False,
              'is3d': False},
             {'name': 'temperatureVertAdvectionMLTend',
              'title': 'Mixed Layer avg Vadv tendency for temperature',
              'units': '$^\circ$C/s',
              'mpas': 'timeMonthly_avg_activeTracerVertAdvectionMLTend_temperatureVertAdvectionMLTend',
              'factor': 1,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4e-6, -3e-6, -2e-6, -1e-6, -0.5e-6, 0.0, 0.5e-6, 1e-6, 2e-6, 3e-6, 4e-6],
              'plot_anomalies': False,
              'is3d': False},
             {'name': 'temperatureHorMixMLTend',
              'title': 'Mixed Layer avg Hmix tendency for temperature',
              'units': '$^\circ$C/s',
              'mpas': 'timeMonthly_avg_activeTracerHorMixMLTend_temperatureHorMixMLTend',
              'factor': 1,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4e-7, -3e-7, -2e-7, -1e-7, -0.5e-7, 0.0, 0.5e-7, 1e-7, 2e-7, 3e-7, 4e-7],
              'plot_anomalies': False,
              'is3d': False},
             {'name': 'temperatureVertMixMLTend',
              'title': 'Mixed Layer avg Vmix tendency for temperature',
              'units': '$^\circ$C/s',
              'mpas': 'timeMonthly_avg_activeTracerVertMixMLTend_temperatureVertMixMLTend',
              'factor': 1,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4e-6, -3e-6, -2e-6, -1e-6, -0.5e-6, 0.0, 0.5e-6, 1e-6, 2e-6, 3e-6, 4e-6],
              'plot_anomalies': False,
              'is3d': False},
             {'name': 'temperatureNonLocalMLTend',
              'title': 'Mixed Layer avg Vmix-nonlocal tendency for temperature',
              'units': '$^\circ$C/s',
              'mpas': 'timeMonthly_avg_activeTracerNonLocalMLTend_temperatureNonLocalMLTend',
              'factor': 1,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4e-7, -3e-7, -2e-7, -1e-7, -0.5e-7, 0.0, 0.5e-7, 1e-7, 2e-7, 3e-7, 4e-7],
              'plot_anomalies': False,
              'is3d': False},
             {'name': 'temperatureTendML',
              'title': 'Mixed Layer avg temporal tendency of temperature',
              'units': '$^\circ$C/s',
              'mpas': 'timeMonthly_avg_activeTracersTendML_temperatureTendML',
              'factor': 1,
              'colormap': plt.get_cmap('RdBu_r'),
              'clevels': [-4e-6, -3e-6, -2e-6, -1e-6, -0.5e-6, 0.0, 0.5e-6, 1e-6, 2e-6, 3e-6, 4e-6],
              'plot_anomalies': False,
              'is3d': False},
             {'name': 'temperature',
              'title': 'Temperature',
              'units': '$^\circ$C',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'factor': 1,
              #'colormap': plt.get_cmap('RdBu_r'),
              #'clevels': [-1.8, -1.0, -0.5, 0.0, 0.5, 2.0, 4.0, 8.0, 12., 16., 22.],
              'colormap': cmocean.cm.balance,
              'clevels': [-2.0, -1.5, -1.0, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0, 1.5, 2.0],
              'plot_anomalies': True,
              'is3d': True},
             {'name': 'salinity',
              'title': 'Salinity',
              'units': 'PSU',
              'mpas': 'timeMonthly_avg_activeTracers_salinity',
              'factor': 1,
              'colormap': cmocean.cm.haline,
              'clevels': [27., 28., 29., 29.5, 30., 30.5, 31., 32., 33., 34., 35.],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'potentialDensity',
              'title': 'Potential Density',
              'units': 'kg m$^{-3}$',
              'mpas': 'timeMonthly_avg_potentialDensity',
              'factor': 1,
              'colormap': cmocean.cm.dense,
              'clevels': [24., 25.5, 25.9, 26.2, 26.5, 26.7, 26.8, 26.85, 26.9, 27.1, 27.75],
              'plot_anomalies': False,
              'is3d': True},
             {'name': 'maxmld',
              'title': 'Maximum Mixed Layer Depth',
              'units': 'm',
              'mpas': 'timeMonthlyMax_max_dThreshMLD',
              'factor': 1,
              #'colormap': plt.get_cmap('viridis'),
              'colormap': cmocean.cm.balance,
              'clevels': [10, 30, 50, 70, 100, 150, 200, 300, 500, 800, 2000],
              'plot_anomalies': False,
              'is3d': False},
             {'name': 'mld',
              'title': 'Mixed Layer Depth',
              'units': 'm',
              'mpas': 'timeMonthly_avg_dThreshMLD',
              'factor': 1,
              #'colormap': plt.get_cmap('viridis'),
              'colormap': cmocean.cm.balance,
              'clevels': [10, 30, 50, 70, 100, 150, 200, 300, 500, 800, 2000],
              'plot_anomalies': False,
              'is3d': False}]

# Identify dictionary for desired variable
vardict = next(item for item in variables if item['name'] == variable)

varname = vardict['name']
mpasvarname = vardict['mpas']
factor = vardict['factor']
plot_anomalies = vardict['plot_anomalies']
is3d = vardict['is3d']
vartitle = vardict['title']
varunits = vardict['units']
clevels = vardict['clevels']
colormap0 = vardict['colormap']
underColor = colormap0(colorIndices0[0])
overColor = colormap0(colorIndices0[-1])
if len(clevels) + 1 == len(colorIndices0):
    # we have 2 extra values for the under/over so make the colormap
    # without these values
    colorIndices = colorIndices0[1:-1]
elif len(clevels) - 1 != len(colorIndices0):
    # indices list must be either one element shorter
    # or one element longer than colorbarLevels list
    raise ValueError('length mismatch between indices and '
                     'colorbarLevels')
colormap = cols.ListedColormap(colormap0(colorIndices))
colormap.set_under(underColor)
colormap.set_over(overColor)
cnorm = cols.BoundaryNorm(clevels, colormap.N)

mesh = xr.open_dataset(meshfile)
lat = mesh.latCell.values
lon = mesh.lonCell.values
lat = np.rad2deg(lat)
lon = np.rad2deg(lon)
z = mesh.refBottomDepth.values
# Find model levels for each depth level
zlevels = np.zeros(np.shape(dlevels), dtype=np.int)
for id in range(len(dlevels)):
    dz = np.abs(z-dlevels[id])
    zlevels[id] = np.argmin(dz)
print('Model levels = ', z[zlevels])

ds = xr.open_mfdataset(infiles, combine='nested', concat_dim='Time')
ntime = ds.dims['Time']

if plot_anomalies:
    figtitle0 = 'Anomaly'
else:
    figtitle0 = ''

if is3d:
    for iz in range(len(dlevels)):
        figfile = '{}/{}{}_depth{:04d}_{}.mp4'.format(figdir, varname, figtitle0, int(dlevels[iz]), runname)
        figtitle0 = '{} {} (z={:5.1f} m)'.format(vartitle, figtitle0, z[zlevels[iz]])

        if varname=='temperatureTotalAdvectionTendency':
            mpasvarname1 = 'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_temperatureHorizontalAdvectionTendency'
            mpasvarname2 = 'timeMonthly_avg_activeTracerVerticalAdvectionTendency_temperatureVerticalAdvectionTendency'
            fld = np.add(ds[mpasvarname1][:, :, zlevels[iz]].values,
                         ds[mpasvarname2][:, :, zlevels[iz]].values)
        elif varname=='temperatureForcingTendency':
            mpasvarname1 = 'timeMonthly_avg_activeTracerSurfaceFluxTendency_temperatureSurfaceFluxTendency'
            mpasvarname2 = 'timeMonthly_avg_temperatureShortWaveTendency'
            fld = np.add(ds[mpasvarname1][:, :, zlevels[iz]].values,
                         ds[mpasvarname2][:, :, zlevels[iz]].values)
        elif varname=='temperatureTendency':
            temp = ds[mpasvarname][:, :, zlevels[iz]].values
            #fld = np.nan*np.ones(np.shape(temp))
            #fld[1:-1, :] = np.diff(temp, n=1, axis=0)
            fld = np.diff(temp, n=1, axis=0, prepend=np.nan*np.ones([1, np.shape(temp)[1]]))/86400. # assumes daily values
        elif varname=='temperatureSumTendencyTerms':
            mpasvarname1 = 'timeMonthly_avg_activeTracerSurfaceFluxTendency_temperatureSurfaceFluxTendency'
            mpasvarname2 = 'timeMonthly_avg_temperatureShortWaveTendency'
            mpasvarname3 = 'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_temperatureHorizontalAdvectionTendency'
            mpasvarname4 = 'timeMonthly_avg_activeTracerVerticalAdvectionTendency_temperatureVerticalAdvectionTendency'
            mpasvarname5 = 'timeMonthly_avg_activeTracerHorMixTendency_temperatureHorMixTendency'
            mpasvarname6 = 'timeMonthly_avg_activeTracerVertMixTendency_temperatureVertMixTendency'
            mpasvarname7 = 'timeMonthly_avg_activeTracerNonLocalTendency_temperatureNonLocalTendency'
            fld = np.add(ds[mpasvarname1][:, :, zlevels[iz]].values,
                         ds[mpasvarname2][:, :, zlevels[iz]].values)
            fld = np.add(fld, ds[mpasvarname3][:, :, zlevels[iz]].values)
            fld = np.add(fld, ds[mpasvarname4][:, :, zlevels[iz]].values)
            fld = np.add(fld, ds[mpasvarname5][:, :, zlevels[iz]].values)
            fld = np.add(fld, ds[mpasvarname6][:, :, zlevels[iz]].values)
            fld = np.add(fld, ds[mpasvarname7][:, :, zlevels[iz]].values)
        else:
            fld = ds[mpasvarname][:, :, zlevels[iz]].values
        fld = factor*fld
        if plot_anomalies:
            fld = fld - fld[0, :]
        #print(varname, int(dlevels[iz]), np.min(fld), np.max(fld))

        fig = plt.figure(figsize=figsize, dpi=figdpi)
        ax = plt.axes(projection=ccrs.Miller(central_longitude=0))
        _add_land_lakes_coastline(ax)
        data_crs = ccrs.PlateCarree()
        ax.set_extent([-180, 180, -90, 90], crs=data_crs)
        gl = ax.gridlines(crs=data_crs, color='k', linestyle=':', zorder=5)
        # This will work with cartopy 0.18:
        #gl.xlocator = mticker.FixedLocator(np.arange(-180., 181., 40.))
        #gl.ylocator = mticker.FixedLocator(np.arange(-80., 81., 20.))

        ax.cla()
        shading = ax.scatter(lon, lat, s=1.0, c=fld[0, :], cmap=colormap, norm=cnorm,
                             marker='o', transform=data_crs)
        cax, kw = mpl.colorbar.make_axes(ax, location='bottom', pad=0.03, shrink=0.9)
        cbar = plt.colorbar(shading, cax=cax, ticks=clevels, **kw)
        cbar.ax.tick_params(labelsize=14, labelcolor='black')
        cbar.set_label(varunits, fontsize=14)
        figtitle = '{} month={:d}'.format(figtitle0, 1)
        ax.set_title(figtitle, y=1.04, fontsize=18)
        #plt.savefig('tmp.png', bbox_inches='tight')

        def animate(i):
            ax.cla()
            shading = ax.scatter(lon, lat, s=1.0, c=fld[i, :], cmap=colormap, norm=cnorm,
                                 marker='o', transform=data_crs)
            figtitle = '{} month={:d}'.format(figtitle0, i+1)
            ax.set_title(figtitle, y=1.04, fontsize=16)

        interval = 100 #in seconds     
        #ani = animation.FuncAnimation(fig, animate, frames=range(5), interval=interval)
        ani = animation.FuncAnimation(fig, animate, frames=range(ntime), interval=interval)
        ani.save(figfile)
else:
    figfile = '{}/{}{}_{}.mp4'.format(figdir, varname, figtitle0, runname)
    figtitle0 = '{} {}'.format(vartitle, figtitle0)

    fld = ds[mpasvarname].values
    fld = factor*fld
    if plot_anomalies:
        fld = fld - fld[0, :]
    #print(varname, np.min(fld), np.max(fld))

    fig = plt.figure(figsize=figsize, dpi=figdpi)
    ax = plt.axes(projection=ccrs.Miller(central_longitude=0))
    _add_land_lakes_coastline(ax)
    data_crs = ccrs.PlateCarree()
    ax.set_extent([-180, 180, -90, 90], crs=data_crs)
    gl = ax.gridlines(crs=data_crs, color='k', linestyle=':', zorder=5)
    # This will work with cartopy 0.18:
    #gl.xlocator = mticker.FixedLocator(np.arange(-180., 181., 40.))
    #gl.ylocator = mticker.FixedLocator(np.arange(-80., 81., 20.))

    ax.cla()
    shading = ax.scatter(lon, lat, s=1.0, c=fld[0, :], cmap=colormap, norm=cnorm,
                         marker='o', transform=data_crs)
    cax, kw = mpl.colorbar.make_axes(ax, location='bottom', pad=0.03, shrink=0.9)
    cbar = plt.colorbar(shading, cax=cax, ticks=clevels, **kw)
    cbar.ax.tick_params(labelsize=14, labelcolor='black')
    cbar.set_label(varunits, fontsize=14)
    figtitle = '{} month={:d}'.format(figtitle0, 1)
    ax.set_title(figtitle, y=1.04, fontsize=18)
    plt.savefig('tmp.png', bbox_inches='tight')

    def animate(i):
        ax.cla()
        shading = ax.scatter(lon, lat, s=1.0, c=fld[i, :], cmap=colormap, norm=cnorm,
                             marker='o', transform=data_crs)
        figtitle = '{} month={:d}'.format(figtitle0, i+1)
        ax.set_title(figtitle, y=1.04, fontsize=18)

    interval = 100 #in seconds     
    #ani = animation.FuncAnimation(fig, animate, frames=range(5), interval=interval)
    ani = animation.FuncAnimation(fig, animate, frames=range(ntime), interval=interval)
    ani.save(figfile)
