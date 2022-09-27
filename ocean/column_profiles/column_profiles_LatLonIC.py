#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/220927_EC30to60_zlevel_noPBC_noSmooth/ocean/global_ocean/EC30to60/PHC/init/initial_state/'
simName = 's05a'
iTime = 0
fileName = 'temperature.nc'
fileName = 'salinity.nc'
domainName = 'EC60to30'

deg2rad = 3.14159/180.0
rad2deg = 180.0/3.14159
if domainName == 'EC60to30':
# global locations
# Atlantic
    lonMid = 360 -12
    latMid =  62
    lonWid = 1.4
    lonWid = 5
# Pacific
    #lonMid = 175
    #latMid =  30
    #lonWid =   2.0

    latWid = lonWid
    lonMin = lonMid-lonWid/2
    lonMax = lonMid+lonWid/2
    latMin = latMid-latWid/2
    latMax = latMid+latWid/2
else:
    print('incorrect domain')

import matplotlib as mpl
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from netCDF4 import Dataset

data = xr.open_dataset(runDir+fileName)

t_lon = data.variables['t_lon']
t_lat = data.variables['t_lat']
depth_t = data.variables['depth_t']
print(t_lon)
print(t_lat)
print(depth_t)

lonList = np.where(np.logical_and(t_lon>lonMin, t_lon<lonMax))[0]
latList = np.where(np.logical_and(t_lat>latMin, t_lat<latMax))[0]

print('lonList',lonList)
print('latList',latList)
fig = plt.figure(figsize=(20,12))
varNames = ['TEMP','SALT']
varNames = ['SALT']
fileNames = ['temperature.nc','salinity.nc']
for iVar in [0]: #range(len(varNames)):
    if iVar==0:
       var = data.variables[varNames[iVar]]
    elif iVar==1:
       dataS = xr.open_dataset(runDir+'salinity.nc')
       print(dataS)
       var = dataS.variables['SALT']
    plt.subplot(2,2,iVar+1)
    var = data.variables[varNames[iVar]]
    print(var)
    print(np.shape(var))
    print('variable:',varNames[iVar])
    for i in range(len(lonList)):
      for j in range(len(latList)):
        varData = var[:,latList[j],lonList[i]]
        print('lon,lat',t_lat[latList[j]],t_lon[lonList[i]])
        #var,Data = np.where(varData>-1e20,varData,np.NAN)
        print(varData)
        maxk=46
        plt.plot(varData[0:maxk],np.arange(maxk))
    plt.gca().invert_yaxis()
    plt.title(varNames[iVar])
    plt.grid()
    plt.ylabel('vertical index, k')
    #plt.legend()

# add information to bottom of figure
today = date.today()
try:
    xtime = data.variables['xtime']
    import codecs
    strXTime = str(codecs.decode(xtime[iTime].values))[2:19]
except:
    strXTime = 'init'
lonMid = (lonMid+180.0)%360.0 - 180.0
plt.figtext(0.1,0.92,domainName+' '+simName+' '+strXTime+ 
    '  lon,lat: '+str(lonMid)+', '+str(latMid)+'       date: '+today.strftime("%d/%m/%Y"))

figfile = 'vert_profiles_IC.png'
print(figfile)
plt.savefig(figfile, bbox_inches='tight')
plt.close()
