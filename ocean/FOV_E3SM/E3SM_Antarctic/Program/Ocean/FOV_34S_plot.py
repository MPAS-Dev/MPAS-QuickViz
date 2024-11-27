#Program plots the F_ovS and the components

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats

#Making pathway to folder with all data
directory 	= '../../Data/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#MOC strength (Sv)
	FOV		    = fh.variables['F_OV'][:]	    #Fresh water
	FOV_ASW		= fh.variables['F_OV_ASW'][:]	#Fresh water
	FOV_AIW		= fh.variables['F_OV_AIW'][:]	#Fresh water
	FOV_NADW	= fh.variables['F_OV_NADW'][:]	#Fresh water
	FOV_ABW		= fh.variables['F_OV_ABW'][:]	#Fresh water
	salt_ASW	= fh.variables['SALT_ASW'][:]	#Salinity
	salt_AIW	= fh.variables['SALT_AIW'][:]	#Salinity
	salt_NADW	= fh.variables['SALT_NADW'][:]	#Salinity
	salt_ABW	= fh.variables['SALT_ABW'][:]	#Salininty
	vel_ASW		= fh.variables['VVEL_ASW'][:]	#Meridional velocity
	vel_AIW		= fh.variables['VVEL_AIW'][:]	#Meridional velocity
	vel_NADW	= fh.variables['VVEL_NADW'][:]	#Meridional velocity
	vel_ABW		= fh.variables['VVEL_ABW'][:]	#Meridional velocity

	fh.close()

	return time, transport, FOV, FOV_ASW, FOV_AIW, FOV_NADW, FOV_ABW, salt_ASW, salt_AIW, salt_NADW, salt_ABW, vel_ASW, vel_AIW, vel_NADW, vel_ABW


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

section_name	= 'section_34S'

#-----------------------------------------------------------------------------------------	

time, transport, FOV, FOV_ASW, FOV_AIW, FOV_NADW, FOV_ABW, salt_ASW, salt_AIW, salt_NADW, salt_ABW, vel_ASW, vel_AIW, vel_NADW, vel_ABW = ReadinData(directory+'Ocean/FOV_index_'+section_name+'.nc')

FOV_rean, FOV_ASW_rean, FOV_AIW_rean, FOV_NADW_rean, FOV_ABW_rean, FOV_rean_gyre = -0.10138855319303171, -0.12769111454122556, 0.12011490376119702, -0.10644935101861515, 0.012637008605611988, 0.2136790553107374

fh 	= netcdf.Dataset(directory+'Ocean/FOV_gyre_'+section_name+'.nc', 'r')

FOV_gyre	= fh.variables['F_gyre'][:]	#Fresh water

fh.close()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_FOV_all	= plot(time, FOV, '-k', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$, E3SM')
graph_FOV_gyre	= plot(time, FOV_gyre, '-r', linewidth = 1.5, label = '$F_{\mathrm{azS}}$, E3SM')
graph_rean_all	= plot(time, np.zeros(len(time))+FOV_rean, '--', color = 'gray', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$, Reanalysis')
graph_rean_gyre	= plot(time, np.zeros(len(time))+FOV_rean_gyre, '--', color = 'firebrick', linewidth = 1.5, label = '$F_{\mathrm{azS}}$, Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.5, 0.5)
ax.set_xlim(500, 600)
ax.grid()
ax.set_xticks([500, 520, 540, 560, 580, 600])

ax.fill_between([-100, 600], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      	= graph_FOV_all + graph_FOV_gyre
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)


graphs	      	= graph_rean_all + graph_rean_gyre
legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax.legend(graphs, legend_labels, loc = 'lower right', ncol=1, framealpha = 1.0, numpoints = 1)
ax.add_artist(legend_1)


ax.set_title('$F_{\mathrm{ovS}}$ and azonal (gyre) component ($F_{\mathrm{azS}}$), E3SM Antarctic')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_E3SM	= plot(time, FOV_ASW, '-k', linewidth = 1.5, label = 'E3SM')
graph_rean	= plot(time, np.zeros(len(time))+FOV_ASW_rean, '--', color = 'gray', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.5, 0.5)
ax.set_xlim(500, 600)
ax.grid()
ax.set_xticks([500, 520, 540, 560, 580, 600])

graphs	      = graph_E3SM + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('Atlantic Surface Water (ASW), E3SM Antarctic')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_E3SM	= plot(time, FOV_AIW, '-k', linewidth = 1.5, label = 'E3SM')
graph_rean	= plot(time, np.zeros(len(time))+FOV_AIW_rean, '--', color = 'gray', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.5, 0.5)
ax.set_xlim(500, 600)
ax.grid()
ax.set_xticks([500, 520, 540, 560, 580, 600])

graphs	      = graph_E3SM + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('Antarctic Intermediate Water (AIW), E3SM Antarctic')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_E3SM	= plot(time, FOV_NADW, '-k', linewidth = 1.5, label = 'E3SM')
graph_rean	= plot(time, np.zeros(len(time))+FOV_NADW_rean, '--', color = 'gray', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.5, 0.5)
ax.set_xlim(500, 600)
ax.grid()
ax.set_xticks([500, 520, 540, 560, 580, 600])

graphs	      = graph_E3SM + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('North Atlantic Deep Water (NADW), E3SM Antarctic')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_E3SM	= plot(time, FOV_ABW, '-k', linewidth = 1.5, label = 'E3SM')
graph_rean	= plot(time, np.zeros(len(time))+FOV_ABW_rean, '--', color = 'gray', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.5, 0.5)
ax.set_xlim(500, 600)
ax.grid()
ax.set_xticks([500, 520, 540, 560, 580, 600])

graphs	      = graph_E3SM + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('Antarctic Bottom Water (ABW), E3SM Antarctic')

show()
