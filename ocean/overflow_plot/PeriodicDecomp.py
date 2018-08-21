#!/usr/bin/env python

# This script generates a periodic decomposition for planar meshes
# with two processors in the x (periodic) dimension, and n_procs/2 
# in the y dimension.
# Allows for the use of periodic meshes where the gemoetry in the 
# periodic dimension (x) differs for the left and right processors.

import sys
import numpy as np
from scipy.io import *

n_procs = int(sys.argv[1])
infile = sys.argv[2]

if (n_procs%2 != 0):
	print 'ERROR: must use even number of processors'
	sys.exit()

mpas = netcdf.netcdf_file(infile,'r')

x = mpas.variables['xCell'][:]
xMin = np.amin(x)
dx = x[1] - x[0]
for i in np.arange(x.shape[0]) + 1:
	if x[i] < xMin + dx - 1.0e-4:
		nx = i
		break

ny = x.shape[0]/nx/(n_procs/2)

outfile  = 'graph.info.part.' + str(n_procs)

fp = open(outfile,'w')

p_ind = -2
for i in np.arange(x.shape[0]):
	iy = i/nx
	ix = i%nx

	if (iy%ny == 0 and ix == 0 and p_ind < n_procs - 2):
		p_ind = p_ind + 2

	if (ix < nx/2):
		p_out = p_ind
	else:
		p_out = p_ind + 1

	fp.write(str(p_out) + '\n')

fp.close()
