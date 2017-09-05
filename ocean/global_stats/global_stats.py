'''
plot global stats
August 2016, Mark Petersen, LANL
'''

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

wd = '/lustre/scratch3/turquoise/mpeterse/runs/'
outputDir='c63n'
dirName1 = 'c62na/forward/analysis_members'
dirName2 = 'c63la/forward/analysis_members'
fileName = 'globalStats.0001-01-01_00.00.00.nc'
titleText = 'MPAS-Ocean stand-alone spin-up (c62n,c63l)'
varList = ['kineticEnergyCellMax', 'kineticEnergyCellAvg']
varListMinMax = ['temperature','salinity']
lt = '-:'

# open a the netCDF file for reading.
filePathName = wd + dirName1 + '/' + fileName
print 'Reading data from: ' + filePathName
ncfile1 = Dataset(filePathName,'r') 
# read the data in variable named 'data'.
t1 = ncfile1.variables['daysSinceStartOfSim']
label1="100 layer"

# open a the netCDF file for reading.
filePathName = wd + dirName2 + '/' + fileName
print 'Reading data from: ' + filePathName
ncfile2 = Dataset(filePathName,'r') 
# read the data in variable named 'data'.
t2 = ncfile2.variables['daysSinceStartOfSim']
label2="60 layer"
   
for var in varListMinMax:
   plt.clf()
   n=0
   dataMin = ncfile1.variables[var+'Min']
   dataMax = ncfile1.variables[var+'Max']
   dataAvg = ncfile1.variables[var+'Avg']
   plt.plot(t1, dataMin, lt[n], t1, dataAvg, lt[n], t1, dataMax, lt[n], label=label1)
   n=1
   dataMin = ncfile2.variables[var+'Min']
   dataMax = ncfile2.variables[var+'Max']
   dataAvg = ncfile2.variables[var+'Avg']
   plt.plot(t2, dataMin, '--', t2, dataAvg, '--', t2, dataMax, '--', label=label2)
   plt.xlabel('time, days')
   plt.ylabel(var)
   plt.title(titleText)
   plt.grid(True,which="both",ls="dotted")
   plt.xlim([0,350])
   plt.legend()
   plt.savefig('f/' + outputDir + '/' + var + '.png')
   print 'Created plot: ' + wd + outputDir + '/' + var + '.png'

for var in varList:
   plt.clf()
   n=0
   data = ncfile1.variables[var]
   plt.plot(t1, data, lt[n], label=label1)
   n=1
   data = ncfile2.variables[var]
   plt.plot(t2, data, '--', label=label2)
   plt.xlabel('time, days')
   plt.ylabel(var)
   plt.title(titleText)
   plt.xlim([0,350])
   plt.legend()
   plt.grid(True,which="both",ls="dotted")
   plt.savefig('f/' + outputDir + '/' + var + '.png')
   print 'Created plot: ' + wd + outputDir + '/' + var + '.png'
  
ncfile1.close()
ncfile2.close()

