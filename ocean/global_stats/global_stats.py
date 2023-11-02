'''
plot global stats
August 2016, Mark Petersen, LANL
'''

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

wd = '/lcrc/group/e3sm/ac.mpetersen/scratch/runs/'
outputDir='ab2_ECwISC'

#dirName1 = 'ocean_model_231030_0b367fb9_ch_gfortran_openmp_master_ec30to60_suite/ocean/global_ocean/EC30to60/WOA23/dynamic_adjustment/simulation/analysis_members'
dirName1 = 'ocean_model_231030_0b367fb9_ch_gfortran_openmp_master_ecwisc30to60_suite/ocean/global_ocean/ECwISC30to60/WOA23/dynamic_adjustment/simulation/analysis_members'
label1="split-explicit"

#dirName2 = 'ocean_model_231030_0b367fb9_ch_gfortran_openmp_AB2_ec30to60_suite_btrdt75/ocean/global_ocean/EC30to60/WOA23/dynamic_adjustment/simulation/analysis_members'
dirName2 = 'ocean_model_231030_0b367fb9_ch_gfortran_openmp_AB2_ecwisc30to60_suite_btrdt75/ocean/global_ocean/ECwISC30to60/WOA23/dynamic_adjustment/simulation/analysis_members'
label2="split-explicit-AB2"

fileName = 'globalStats.0001-01-31_00.00.00.nc'
titleText = 'MPAS-Ocean stand-alone spin-up'
varList = ['kineticEnergyCellMax', 'kineticEnergyCellAvg']
varListMinMax = ['temperature','salinity']
lt = '-:'

# open a the netCDF file for reading.
filePathName = wd + dirName1 + '/' + fileName
print('Reading data from: ' + filePathName)
ncfile1 = Dataset(filePathName,'r') 
# read the data in variable named 'data'.
t1 = ncfile1.variables['daysSinceStartOfSim'][:]

# open a the netCDF file for reading.
filePathName = wd + dirName2 + '/' + fileName
print('Reading data from: ' + filePathName)
ncfile2 = Dataset(filePathName,'r') 
# read the data in variable named 'data'.
t2 = ncfile2.variables['daysSinceStartOfSim'][:]
   
for var in varListMinMax:
   plt.clf()
   n=0
   dataMin = ncfile1.variables[var+'Min'][:]
   dataMax = ncfile1.variables[var+'Max'][:]
   dataAvg = ncfile1.variables[var+'Avg'][:]
   plt.plot(t1, dataMin, lt[n], label=label1)
   plt.plot(t1, dataAvg, lt[n], label=label1)
   plt.plot(t1, dataMax, lt[n], label=label1)
   n=1
   dataMin = ncfile2.variables[var+'Min'][:]
   dataMax = ncfile2.variables[var+'Max'][:]
   dataAvg = ncfile2.variables[var+'Avg'][:]
   plt.plot(t2, dataMin, '--', label=label2)
   plt.plot(t2, dataAvg, '--', label=label2)
   plt.plot(t2, dataMax, '--', label=label2)
   plt.xlabel('time, days')
   plt.ylabel(var)
   plt.title(titleText)
   plt.grid(True,which="both",ls="dotted")
   plt.xlim([0,350])
   plt.legend()
   plt.savefig('f/' + outputDir + '/' + var + '.png')
   print('Created plot: ' + wd + outputDir + '/' + var + '.png')

for var in varList:
   plt.clf()
   n=0
   data = ncfile1.variables[var][:]
   plt.plot(t1, data, lt[n], label=label1)
   n=1
   data = ncfile2.variables[var][:]
   plt.plot(t2, data, '--', label=label2)
   plt.xlabel('time, days')
   plt.ylabel(var)
   plt.title(titleText)
   plt.xlim([0,350])
   plt.legend()
   plt.grid(True,which="both",ls="dotted")
   plt.savefig('f/' + outputDir + '/' + var + '.png')
   print('Created plot: ' + wd + outputDir + '/' + var + '.png')
  
ncfile1.close()
ncfile2.close()

