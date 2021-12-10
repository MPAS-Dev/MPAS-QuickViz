# Author : Nairita Pal
# Date: November, 2020
# Modified by Kristin Barton, June 2021
# Description: Interpolates CFSR atmospheric reanalysis data onto the MPAS-O mesh and 
#              creates an input file to support time varying atmospheric forcing in the model

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
import statistics
import math
import cmocean
import matplotlib
import sys
import getopt

#####
plt.switch_backend('agg')
cartopy.config['pre_existing_data_dir'] = \
    os.getenv('CARTOPY_DIR', cartopy.config.get('pre_existing_data_dir'))
cmap_reversed = matplotlib.cm.get_cmap('Spectral_r')
#####

def main(argv):
      
    # Setup for argument input
    resolution = ''
    saltype = ''
    
    # Get arguments -- optional (for specifying resolution and run type in title)
    try:
        opts, args = getopt.getopt(argv,"hr:s:c:",["resolution=","saltype==","constituent="])
    except getopt.GetoptError:
        print('plot_rmse_fixed.py -r <resolution> -s <saltype>')
    for opt, arg in opts:
        if opt in ('-h'):
            print('plot_rmse_fixed.py -r <resolution> -s <saltype>')
            sys.exit()
        elif opt in ("-r", "--resolution"):
            resolution = arg
            print('Resolution is ', resolution)
        elif opt in ("-s", "--saltype"):
            saltype = arg
            print('SAL type is ', saltype)

    # Initialize plotting variables
    LW   = 3                         # plot line width
    MS   = 1.5                       # plot symbol size 
    LF   = 30                        # label font size
    TW   = 2                         #Tick width
    TL   = 2                         #Tick length
    TF   = 8  
    i    = 0
    
    # Open data file
    data_file = 'analysis_members/harmonicAnalysis.nc'
    data_nc = netCDF4.Dataset(data_file,'r')
 
    lon_grid = np.mod(data_nc.variables['lonCell'][:] + np.pi, 2.0 * np.pi) - np.pi
    lon_grid = lon_grid*180.0/np.pi
    lat_grid = data_nc.variables['latCell'][:]*180.0/np.pi
    nCells = lon_grid.size
    data1 = np.zeros((nCells))
    data2 = np.zeros((nCells))
    depth = np.zeros((nCells))
    area = np.zeros((nCells))
    depth[:] = data_nc.variables['bottomDepth'][:]
    area[:] = data_nc.variables['areaCell'][:]
   
    data1_phase = np.zeros((nCells))
    data2_phase = np.zeros((nCells))
    cell = 0
    rmse_sum = 0
    count = 0
    shallow_rmse_amp = 0.0
    deep_rmse_amp = 0.0
    areaTotal = 0.0

    constituent_list = ['K1','M2','N2','O1','S2']
    constituent_num = 0
    subplot = 0

    # Use these to fix up the plots
    subplot_levels = [[np.linspace(0,0.65,16), np.linspace(0,0.65,16), np.linspace(0,0.13,16), np.linspace(0,0.13,16)], \
                      [np.linspace(0,1.4, 16),  np.linspace(0,1.4,16), np.linspace(0,0.22,16), np.linspace(0,0.5,16)], \
                      [np.linspace(0,0.22,16), np.linspace(0,0.22,16), np.linspace(0,0.05,16),np.linspace(0,0.025, 16)], \
                      [np.linspace(0,0.5, 16), np.linspace(0,0.5, 16), np.linspace(0,0.08,16),np.linspace(0,0.08,16)], \
                      [np.linspace(0,0.7, 16), np.linspace(0,0.7,16), np.linspace(0,0.5, 16), np.linspace(0,0.5,16)]]

    subplot_ticks = [[np.linspace(0,0.65, 10), np.linspace(0,0.65,10), np.linspace(0,0.13,10), np.linspace(0,0.13,10)], \
                      [np.linspace(0,1.4, 10),  np.linspace(0,1.4,10), np.linspace(0,0.22,10), np.linspace(0,0.25,10)], \
                      [np.linspace(0,0.22,10), np.linspace(0,0.22,10), np.linspace(0,0.05,10),np.linspace(0,0.05, 10)], \
                      [np.linspace(0,0.5, 10), np.linspace(0,0.5, 10), np.linspace(0,0.08,10),np.linspace(0,0.08,10)], \
                      [np.linspace(0,0.7, 10), np.linspace(0,0.7,10), np.linspace(0,0.5, 10), np.linspace(0,0.5,10)]]

    # Calculates in order: ['K1','M2','N2','O1','S2']
    for constituent_num in range(0,4):
        constituent = constituent_list[constituent_num]

        print(" ====== " + constituent + " Constituent ======")
        
        # Get data
        data1[:] = data_nc.variables[constituent+'Amplitude'][:]
        data1_phase[:] = data_nc.variables[constituent+'Phase'][:]*np.pi/180
        data2[:] = data_nc.variables[constituent+'AmplitudeTPXO8'][:]  
        data2_phase[:] = data_nc.variables[constituent+'PhaseTPXO8'][:]*np.pi/180

        # Calculate RMSE values
        rmse_amp = 0.5*(data1 - data2)**2
        rmse_com = 0.5*(data2**2 + data1**2) - data1*data2*np.cos(data2_phase - data1_phase)

        # Calculate mean (global) values
        rmse_amp_sum = 0
        rmse_com_sum = 0
        areaTotalAmp = 0
        areaTotalCom = 0
        for cell in range(0,nCells-1):
            if (depth[cell] > 20) and (rmse_com[cell] < 1000) and (rmse_amp[cell] < 1000):
                if (rmse_amp[cell] < 10000):
                    rmse_amp_sum += rmse_amp[cell]*area[cell]
                    areaTotalAmp += area[cell]
                if (rmse_com[cell] < 10000):
                    rmse_com_sum += rmse_com[cell]*area[cell]
                    areaTotalCom += area[cell]
        global_rmse_amp = np.sqrt(rmse_amp_sum / areaTotalAmp)
        global_rmse_com = np.sqrt(rmse_com_sum / areaTotalCom)
        print('Global RMSE (Amp) = ', global_rmse_amp)
        print('Global RMSE (Com) = ', global_rmse_com)

        # Calculate shallow RMSE (<=1000m)
        rmse_amp_sum = 0
        rmse_com_sum = 0
        areaTotalAmp = 0
        areaTotalCom = 0
        for cell in range(0,nCells-1):
            if (abs(lat_grid[cell]) < 66) and (depth[cell] < 1000) and (depth[cell] > 20):
                if (rmse_amp[cell] < 10000):
                    rmse_amp_sum += rmse_amp[cell]*area[cell]
                    areaTotalAmp += area[cell]
                if (rmse_com[cell] < 10000):
                    rmse_com_sum += rmse_com[cell]*area[cell]
                    areaTotalCom += area[cell]
        shallow_rmse_amp = np.sqrt(rmse_amp_sum / areaTotalAmp)
        shallow_rmse_com = np.sqrt(rmse_com_sum / areaTotalCom)
        print('Shallow RMSE (Amp) = ', shallow_rmse_amp)
        print('Shallow RMSE (Com) = ', shallow_rmse_com)

        # Calculate deep RMSE (>1000m)
        rmse_amp_sum = 0
        rmse_com_sum = 0
        areaTotal = 0
        for cell in range(0,nCells-1):
            if (abs(lat_grid[cell]) < 66) and (depth[cell] >= 1000):
                rmse_amp_sum += rmse_amp[cell]*area[cell]
                rmse_com_sum += rmse_com[cell]*area[cell]
                areaTotal += area[cell]
        deep_rmse_amp = np.sqrt(rmse_amp_sum / areaTotal)
        deep_rmse_com = np.sqrt(rmse_com_sum / areaTotal)
        print('Deep RMSE (Amp) = ', deep_rmse_amp)
        print('Deep RMSE (Com) = ', deep_rmse_com)

        rmse_amp = rmse_amp**0.5
        rmse_com = rmse_com**0.5

        # Plot data
        fig=plt.figure(figsize=(18,12))
        subplot_title = [constituent+' Amplitude (simulation) [m]', constituent+' Amplitude (TPXO8) [m]', \
                         constituent+' RMSE (Amplitude) [m]', constituent+' RMSE (Complex) [m]']

        # Setup the subplot
        for subplot in range(0,4) :
            ax = fig.add_subplot(2,2,subplot+1,projection = ccrs.PlateCarree())
            ax.set_title(subplot_title[subplot],fontsize=20)
            levels = subplot_ticks[constituent_num][subplot][:]
            if subplot == 0 :
                cf = ax.tricontourf(lon_grid,lat_grid,data1,levels=levels,
                    transform=ccrs.PlateCarree(),cmap=cmap_reversed)
                ax.tricontour(lon_grid,lat_grid,data1_phase,levels=10, linewidths=0.5, colors='k')
            elif subplot == 1 :
                iid = np.logical_and(data2_phase>=0, data2_phase < 360)
                cf = ax.tricontourf(lon_grid,lat_grid,data2,levels=levels,
                    transform=ccrs.PlateCarree(),cmap=cmap_reversed)
                ax.tricontour(lon_grid[iid],lat_grid[iid],data2_phase[iid],levels=10, linewidths=0.5, colors='k')
            elif subplot == 2 :
                cf = ax.tricontourf(lon_grid,lat_grid,rmse_amp,levels=levels,
                    transform=ccrs.PlateCarree(),cmap='OrRd')
            elif subplot == 3 :
                cf = ax.tricontourf(lon_grid,lat_grid,rmse_com,levels=levels,
                    transform=ccrs.PlateCarree(),cmap='OrRd')
            ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.LAND, zorder=100)
            ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
            ax.add_feature(cfeature.COASTLINE, zorder=101)
            ax.tick_params(axis='both', which='major', length=TL, width=TW, labelsize=TF)
            cbar = fig.colorbar(cf,ax=ax,ticks=levels.round(2),shrink=0.6)
            cbar.ax.tick_params(labelsize=16) 

        fig.tight_layout()
        fig.suptitle(resolution + ' with run type: ' + saltype + '\n' +
                        'Complex: Global Avg = ' + str(round(global_rmse_com*100,3)) + 'cm' + ';'\
                        ' Deep RMSE = ' + str(round(deep_rmse_com*100,3)) + 'cm' + \
                        '; Shallow RMSE = ' + str(round(shallow_rmse_com*100,3)) + 'cm', \
                        fontsize=20)
        plt.savefig(constituent+'_plot.png')
        plt.close()

if __name__ == '__main__':
    main(sys.argv[1:])
