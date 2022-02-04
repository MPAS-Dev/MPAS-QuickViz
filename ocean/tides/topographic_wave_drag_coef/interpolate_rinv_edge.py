# Author: Steven Brus
# Date: August, 2019
# Description: Interpolates CFSR atmospheric reanalysis data onto the MPAS-O mesh and 
#              creates an input file to support time varying atmospheric forcing in the model

# note, the original file jsl_lim24_inv_hrs.nc
# is available at https://drive.google.com/file/d/1CrMtAbiciJiozzGtZnG4X_VDZLiYW8Uo/view?usp=sharing

import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import glob
import pprint
import datetime
import os
import yaml
import subprocess
import argparse
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import interpolate
import write_forcing_file_rinv
plt.switch_backend('agg')
cartopy.config['pre_existing_data_dir'] = \
    os.getenv('CARTOPY_DIR', cartopy.config.get('pre_existing_data_dir'))

##################################################################################################
##################################################################################################

def interpolate_data_to_grid(grid_file,data_file,var):

  # Open files
  data_nc = netCDF4.Dataset(data_file,'r')
  grid_nc = netCDF4.Dataset(grid_file,'r')

  # Get grid from data file
  lon_data = data_nc.variables['Longitude'][:]
  lat_data = data_nc.variables['Latitude'][:]
  time = data_nc.variables['MT'][:]
  nsnaps = time.size
  nlon = lon_data.size
  nlat = lat_data.size
  print(np.amax(lon_data),'HYCOM max lon data',np.amin(lon_data),'min lon data')
  print(np.amax(lat_data),'HYCOM max lat data',np.amin(lat_data),'HYCOM min lat data')

  # Get grid from grid file 
  lon_grid = np.mod(grid_nc.variables['lonEdge'][:] + np.pi, 2.0 * np.pi) - np.pi
  lon_grid = lon_grid*180.0/np.pi
  lat_grid = grid_nc.variables['latEdge'][:]*180.0/np.pi
  print(np.amax(lon_grid),'max longitude in MPAS mesh',np.amin(lon_grid),'min lon mesh')
  print(np.amax(lat_grid),'max lat in MPAS mesh',np.amin(lat_grid),'min lat mesh')

  grid_points = np.column_stack((lon_grid,lat_grid))
  nEdges = lon_grid.size
  interp_data = np.zeros((nsnaps,nEdges))
  print(interp_data.shape,'interp')
  print(np.amin(lon_grid),np.amax(lon_grid))
  print(np.amin(lat_grid),np.amax(lat_grid))

  # Interpolate timesnaps
  print('Interpolating '+var)

  # Get data to interpolate
  data = data_nc.variables[var][:]

  # Interpolate data onto new grid
  interpolator = interpolate.RegularGridInterpolator((lon_data,lat_data),1/data[0,:,:].T,
                                                     bounds_error=False,fill_value=None,method='nearest')
  interp_data[0,:] = interpolator(grid_points)
   
 # iid3 = np.logical_and(inter_data>=-1, diff < 1)
 


  return lon_grid,lat_grid,interp_data,lon_data,lat_data,data

##################################################################################################
##################################################################################################

def plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,var_label,var_abrev):


  # Plot data
  fig = plt.figure()
  levels = np.linspace(np.amin(data),np.amax(data),100)
  ax0 = fig.add_subplot(2, 1, 1, projection=ccrs.PlateCarree())
  ds = 10                                # Downsample
  dsx = np.arange(0, lon_data.size, ds)  # data
  dsy = np.arange(0, lat_data.size, ds)  # to speed up
  dsxy = np.ix_(dsy, dsx)                # plotting
  cf = ax0.contourf(lon_data[dsx], lat_data[dsy], data[dsxy], levels=levels,
                    transform=ccrs.PlateCarree())
  ax0.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
  ax0.add_feature(cfeature.LAND, zorder=100)
  ax0.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax0.add_feature(cfeature.COASTLINE, zorder=101)
  ax0.set_title('data')
  cbar = fig.colorbar(cf,ax=ax0)
  cbar.set_label(var_label)

  # Plot interpolated data
  interp_data_plot = np.ma.masked_where(interp_data > 1e20, interp_data)
  ax1 = fig.add_subplot(2, 1, 2, projection=ccrs.PlateCarree())
  cf = ax1.scatter(lon_grid,lat_grid,s=0.5,c=interp_data_plot,marker='o',edgecolors='none',
                   vmin=np.amin(data),vmax=np.amax(data), transform=ccrs.PlateCarree())
  ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
  ax1.add_feature(cfeature.LAND, zorder=100)
  ax1.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax1.add_feature(cfeature.COASTLINE, zorder=101)
  ax1.set_title('interpolated data')
  cbar = fig.colorbar(cf,ax=ax1)
  cbar.set_label(var_label)

  # Save figure
  fig.savefig('test_rinv.png',dpi=500,bbox_inches='tight')
  plt.close()

##################################################################################################
##################################################################################################

if __name__ == '__main__':

  # Files to interpolate to/from
  grid_file = './initial_state.nc'
  rinv_file = './jsl_lim24_inv_hrs.nc'
  rinv_write_file = './topographic_wave_drag.nc'

  #Interpolation of rinv
  lon_grid,lat_grid,rinv_interp,lon_data,lat_data,rinv_data = interpolate_data_to_grid(grid_file,rinv_file,'rinv')

  # Plot rinv
  print(np.amin(rinv_data),'min data')
  print(np.amax(rinv_data),'max data')
  print(rinv_data.shape,'shape')
  #incorr_id = np.logical_or(rinv_data<24, rinv_data > 240)
  #rinv_data[incorr_id] = np.nan
  #print(np.amin(rinv_data),'min data')
  #print(np.amax(rinv_data),'max data')
  #print(rinv_data.shape,'shape')
  print('Plotting e-folding time: ')
  #plot_interp_data(lon_data,lat_data,rinv_data[0,:,:],lon_grid,lat_grid,rinv_interp[0,:],'e-folding time','rinv')

  # Write to NetCDF file
  subprocess.call(['rm',rinv_write_file])
  write_forcing_file_rinv.write_to_file(rinv_write_file,rinv_interp,'topographic_wave_drag')
