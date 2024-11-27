#Program plots the resolution of the native grid

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
import matplotlib.tri as tri
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh	= netcdf.Dataset('/global/cfs/cdirs/m4259/mapping_files/map_SOwISC12to60E2r4_to_0.5x0.5degree_bilinear.nc', 'r')

lon	= fh.variables['xc_a'][:] * 180 / np.pi
lat	= fh.variables['yc_a'][:] * 180 / np.pi
area_a	= fh.variables['area_a'][:] * (180 / np.pi)**2.0

fh.close()

lon[lon > 180]	= lon[lon > 180] - 360.0

print(np.max(np.sqrt(area_a)))

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS 	= ax.tripcolor(lon, lat, np.sqrt(area_a), vmin=0, vmax=0.6, cmap='Spectral_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(0, 0.61, 0.2), cax=ax_cb)
cbar.set_label('Horizontal resolution ($^{\circ}$)')

ax.set_global()

ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
ax.coastlines()

ax.set_title('Grid resolution, E3SM Antarctic')
    
show()