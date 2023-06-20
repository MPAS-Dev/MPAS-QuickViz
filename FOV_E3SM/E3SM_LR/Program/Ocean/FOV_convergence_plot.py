#Program plots the freshwater convergence (34S and 60N)

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory 	= '../../Data/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	#Fresh water

	fh.close()

	return time, FOV

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time, FOV_34S	= ReadinData(directory+'Ocean/FOV_index_section_34S.nc')
time, FOV_60N	= ReadinData(directory+'Ocean/FOV_index_section_60N.nc')
#-----------------------------------------------------------------------------------------

FOV_34S_rean, FOV_60N_rean = -0.10138855319303171, -0.027075354933136512
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_rcp_34S		= ax.plot(time, FOV_34S, '-k', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$, E3SM')
graph_rcp_60N		= ax.plot(time, FOV_60N, '-b', linewidth = 1.5, label = '$F_{\mathrm{ovN}}$, E3SM')
graph_rcp_conver	= ax.plot(time, FOV_34S - FOV_60N, '-r', linewidth = 1.5, label = '$\Delta F_{\mathrm{ov}}$, E3SM')

graph_rean_34S		= ax.plot(time, np.zeros(len(time))+FOV_34S_rean, '--', color = 'gray', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$, Reanalysis')
graph_rean_60N		= ax.plot(time, np.zeros(len(time))+FOV_60N_rean, '--', color = 'cyan', linewidth = 1.5, label = '$F_{\mathrm{ovN}}$, Reanalysis')
graph_rean_conver	= ax.plot(time, np.zeros(len(time))+FOV_34S_rean - FOV_60N_rean, '--', color = 'firebrick', linewidth = 1.5, label = '$\Delta F_{\mathrm{ov}}$, Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 60)
ax.set_ylim(-0.5, 0.5)
ax.set_xticks([1, 10, 20, 30, 40, 50, 60])
ax.grid()

ax.fill_between([00, 100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      	= graph_rcp_34S + graph_rcp_60N + graph_rcp_conver
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)


graphs	      	= graph_rean_34S + graph_rean_60N + graph_rean_conver
legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax.legend(graphs, legend_labels, loc = 'upper right', ncol=1, framealpha = 1.0, numpoints = 1)
ax.add_artist(legend_1)


ax.set_title('f) Freshwater convergence, LR-E3SM')

show()
