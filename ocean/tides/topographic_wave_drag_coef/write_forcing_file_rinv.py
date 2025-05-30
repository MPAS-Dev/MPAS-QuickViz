# Author: Steven Brus
# Date: April, 2020
# Description: This function writes time-varying forcing data to an input file for the model run.

# note, the original file jsl_lim24_inv_hrs.nc
# is available at https://drive.google.com/file/d/1CrMtAbiciJiozzGtZnG4X_VDZLiYW8Uo/view?usp=sharing

import os
import numpy as np
import netCDF4

##################################################################################################
##################################################################################################

def write_to_file(filename,data,var):

  if os.path.isfile(filename):
    data_nc = netCDF4.Dataset(filename,'a', format='NETCDF3_64BIT_OFFSET')
  else:
    data_nc = netCDF4.Dataset(filename,'w', format='NETCDF3_64BIT_OFFSET')

    # Find dimesions
    nEdges = data.shape[1]
    nsnaps = data.shape[0]

    # Declare dimensions
    data_nc.createDimension('nEdges',nEdges)
    data_nc.createDimension('StrLen',64)
    data_nc.createDimension('Time',None)

    # Create time variable
   # time = data_nc.createVariable('xtime','S1',('Time','StrLen'))
   # time[:,:] = netCDF4.stringtochar(xtime)

  # Set variables
  data_var = data_nc.createVariable(var,np.float64,('Time','nEdges'))
  data_var[:,:] = data[:,:]
  data_nc.close()

##################################################################################################
##################################################################################################
