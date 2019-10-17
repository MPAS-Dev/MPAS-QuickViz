'''
sections.py
Compute zonal or meridional sections of variables on an mpas mesh

# Regrid file to lat/lon coordinates using unix commands:
mkdir regridded
ncks -O output/debugTracer.0001-01-01_00.00.00.nc regridded/temp.nc
ncks -v latCell,lonCell -A init.nc regridded/temp.nc
cd regridded
ln -isf /usr/projects/climate/mpeterse/repos/APrime_Files/mapping/maps/* .
ncremap -i temp.nc -o debugTracersLatLon.nc -P mpas -m map_oEC60to30v3_TO_0.5x0.5degree_blin.nc -R "--rgr lat_nm=latCell --rgr lon_nm=lonCell --rgr lat_nm_out=lat --rgr lon_nm_out=lon" -C
'''

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import datetime

# Input arguments
path = '/lustre/scratch4/turquoise/mpeterse/runs'
runName = 'redi26'
meshFile = 'init.nc'
dataFile = 'regridded/debugTracersLatLon.nc'
titleTxt = 'EC60to30, MPAS-Ocean stand alone, with Redi mixing on'
varNames = ['temperature','salinity','potentialDensity','relativeSlopeTopOfCell','relativeSlopeTaperingCell','tracer1','tracer2']
landValue = [-0.1,-0.1,1027, -0.1, -0.1, -0.1,-0.1]
layer = [0, 22] # layers for horizontal plots, columns 0 and 1.
layerSpan = [0,35] # vertical span of cross-section plots, columns 2 and 3.
lonRequest = [-150.0, -110]
latSpan = [-70,10]
iTime = 4

# load mesh variables and time
ncfile1 = Dataset(path+'/'+runName+'/'+dataFile,'r')
xtime = ncfile1.variables['xtime']
lat = ncfile1.variables['lat']
lon = ncfile1.variables['lon']
ncMeshFile = Dataset(path+'/'+runName+'/'+meshFile,'r')
refBottomDepth = ncMeshFile.variables['refBottomDepth']

# set up indices
nRows=len(varNames)
nCols=4
nLayerCols=2
iLat0 = np.where(lat[:]>latSpan[0])[0][0]
iLat1 = np.where(lat[:]>latSpan[1])[0][0]

# Set up figure, add title text
fig = plt.gcf()
fig.set_size_inches(18.0,16.0)
plt.clf()
plt.subplot2grid((nRows+1, nCols), (0,0), colspan=3)
plt.text(.1,.8,titleTxt, fontsize=14, fontweight='bold', verticalalignment='bottom', horizontalalignment='left')
plt.text(.1,.6,'time='+str(xtime[2,:],'utf-8'), fontsize=12, fontweight='normal', verticalalignment='bottom', horizontalalignment='left')
plt.text(.1,.4,'path: '+path+'/'+dataFile, fontsize=12, fontweight='normal', verticalalignment='bottom', horizontalalignment='left')
now = datetime.datetime.now()
plt.text(.1,.2,'created: '+now.strftime("%Y-%m-%d"), fontsize=12, fontweight='normal', verticalalignment='bottom', horizontalalignment='left')
plt.gca().axis('off')

# plot refBottomDepth
plt.subplot2grid((nRows+1, nCols), (0,3))
plt.plot(refBottomDepth[layerSpan[0]:layerSpan[1]+1],np.arange(layerSpan[0],layerSpan[1]+1))
plt.ylabel('vertical index')
plt.xlabel('depth, m')
plt.gca().invert_yaxis()
plt.gca().set_aspect(14)
plt.grid(True)

for iRow in range(nRows):
    var = ncfile1.variables[varNames[iRow]]
    for iCol in range(nCols):
        plt.subplot2grid((nRows+1, nCols), (iRow+1,iCol) )
        plt.subplot(nRows+1, nCols, (iRow+1)*nCols+iCol+1)
        iLon = np.where(lon[:]>lonRequest[iCol-nLayerCols])[0][0]

        # horizontal plots
        if iCol<nLayerCols:
            tmp = var[iTime,layer[iCol],:,:]
            tmp2 = np.where(tmp>-1e20,tmp,landValue[iRow])
            ax = plt.imshow(tmp2,extent=[-180,180,-90,90],aspect=1.0)
            plt.title(varNames[iRow]+', layer '+str(layer[iCol]))
            plt.gca().invert_yaxis()
            if varNames[iRow]=='salinity':
                plt.clim(33,37)
            if varNames[iRow]=='potentialDensity':
                plt.clim(1023, 1028)

            if iRow == nRows-1:
                plt.xlabel('longitude')
            else:
                plt.gca().get_xaxis().set_visible(False)

        # cross sections
        else:
            tmp = var[iTime,layerSpan[0]:layerSpan[1],iLat0:iLat1,iLon]
            tmp2 = np.where(tmp>-1e20,tmp,landValue[iRow])
            ax = plt.imshow(tmp2,extent=[lat[iLat0],lat[iLat1],layerSpan[1],layerSpan[0]])
            plt.title(varNames[iRow]+', lon='+str(lon[iLon]))
            if iRow == nRows-1:
                plt.xlabel('latitude')
        if iCol==2:
            plt.ylabel('vertical index')

        plt.set_cmap('gist_ncar')
        plt.colorbar()

plt.tight_layout()
plt.savefig('variables.png')
