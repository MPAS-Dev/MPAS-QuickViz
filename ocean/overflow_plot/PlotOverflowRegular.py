#!/usr/bin/env python
# usage:
# ./PlotOverflowRegular.py <field_name> <min> <max> <file_name>
# 
# Note that this requires zMid to be dumped in the streams file.
# 
# if <min> and <max> are both < 1.0E-6 then the native bounds of the field
# 		    will be used for the contour plot.

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import *

# input parameters
#run = 'cdgq'
#filename = 'output.overflow.' + run + '.nc'
filename = sys.argv[4]
xDepth = 16
field = sys.argv[1]
tmin = float(sys.argv[2])
tmax = float(sys.argv[3])

# load file
mpas = netcdf.netcdf_file(filename,'r')

x = mpas.variables['xCell'][:]
y = mpas.variables['yCell'][:]
z = mpas.variables['zMid'][:,:,:]
#zt = mpas.variables['zTop'][:,:,:]
t = mpas.variables[field][:,:,:]
nz = mpas.dimensions['nVertLevels']

ny = y.shape[0]/xDepth

yg = np.zeros(ny)
tg = np.zeros((nz,ny))

for i in np.arange(ny):
	j = 16*(i%ny)
	yg[i] = y[j]

nsteps = t.shape[0]
print 'num steps: ' + str(nsteps)

levs = np.linspace(tmin,tmax,101,endpoint=True)

# ensure layer thicknesses and bottom depth are correct

thick = mpas.variables['layerThickness'][0,:,:]
depth = mpas.variables['bottomDepth'][:]
maxLevelCell = mpas.variables['maxLevelCell'][:]
for c in np.arange(depth.shape[0]):
	h = depth[c]
	for l in np.arange(maxLevelCell[c]):
		h = h - thick[c,l]

	if abs(h) > 1.0e-8:
		print 'ssh error! ' + str(h)

for step in np.arange(nsteps):
	print 'plotting field for step ' + str(step)

	yy = np.zeros((nz,ny))
	zz = np.zeros((nz,ny))
	for i in np.arange(nz):
		yy[i][:] = yg[:]
		for j in np.arange(ny):
			zz[i][j] = z[step][xDepth*j][i]

	#for xInd in np.arange(xDepth):
	for xInd in np.arange(1):
		for i in np.arange(ny):
			for j in np.arange(nz):
				k = xDepth*i + xInd
				tg[j][i] = t[step][k][j]

		#outputfile = 'output/' + field + '_%.4d'%step + '_' + str(xInd) + '_' + run + '.png'
		outputfile = 'overflow/' + field + '_%.4d'%(step+40) + '.png'
		print 'output file: ' + outputfile

		if abs(tmin) > 1.0e-6 or abs(tmax) > 1.0e-6:
			plt.contourf(yy, zz, tg, levs)
		else:
			plt.contourf(yy, zz, tg, 100)

		#for lev in np.arange(zt.shape[2]):
		#	plt.plot(yg, zz[lev,:], '-', c='k')

		plt.ylim([-2000.0,10.0])
		#plt.colorbar(ticks=[10.0,12.0,14.0,16.0,18.0,20.0])
		plt.savefig(outputfile)
		plt.clf()
