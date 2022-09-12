import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.dates as mdates
import matplotlib as mpl
import xarray
import netCDF4 as nc
import glob

mpl.rcParams['figure.figsize'] = (20,12) # Large figures
mpl.rcParams['image.cmap'] = 'seismic'
#mpl.rcParams['image.cmap'] = 'Spectral'
dpi=200;

#wd = '/lustre/scratch5/turquoise/mpeterse/runs/s/'
#rundir = 'm01a/'
#subdir = 'hourlySnapshots/'
#mesh = xarray.open_dataset(wd+rundir+subdir+'init.nc')

#flist = sorted(glob.glob('/lcrc/group/e3sm/ac.abarthel/E3SMv2/20220517.v2.amoc.GMPASO-CORE.3Drest.T62_oRRS30to10v3.anvil/run/*mpaso.hist.am.mocStreamfunctionOutput.0013*')) #-- Snapshots
#flist2 = sorted(glob.glob('/lcrc/group/e3sm/ac.abarthel/E3SMv2/20220517.v2.amoc.GMPASO-CORE.3Drest.T62_oRRS30to10v3.anvil/run/*mpaso.hist.am.timeSeriesStatsMonthly.0013*')) # Monthly averages

flist = sorted(glob.glob('mpaso.hist.hourly.0001-02-01_00000.nc'))
nt=0


latMinDegArr = np.arange(5,86,10)
d2r=np.pi/180.
r2d=180./np.pi
width = 1.0
for lat in latMinDegArr:
    plt.clf()
    lonMinDeg = -170
    latMinDeg = lat
    lonMin = (lonMinDeg+360) * d2r
    latMin = latMinDeg * d2r
    lonMax = lonMin + width * d2r
    latMax = latMin + width * d2r
    index_provided = False
    if index_provided:
        cellList = 23883
        iCell = cellList
    else:
        mesh = xarray.open_dataset('../init.nc')
        latCell = mesh.variables['latCell']
        lonCell = mesh.variables['lonCell']
        cellList = np.where(np.logical_and(np.logical_and(np.logical_and(latCell>latMin, latCell<latMax), lonCell>lonMin), lonCell<lonMax))
        iCell = cellList[0][0]
        print('requested  lon,lat: ',lonMinDeg,latMinDeg)
        print('found cell lon,lat: ',float(lonCell[iCell])*r2d,float(latCell[iCell])*r2d)
        fCyclesPerDay = 2*np.sin(float(latCell[iCell]))
    
    f = nc.Dataset(flist[nt], 'r')
    f.variables['daysSinceStartOfSim'].set_auto_scale(False)
    daysSinceStartOfSim = f.variables['daysSinceStartOfSim'][:]
    
    print('create plot')
    titleTxt = [', initial', ', time 1', ', time 2']
    varNames =['vertVelocityTop','velocityZonal','velocityMeridional']
    nCol = 3
    nRow = 1
    levels = np.arange(5,60,5)
    colors = pl.cm.jet_r(np.linspace(0,1,len(levels)))
    for iRow in range(nRow):
        for iCol in range(nCol):
            var = f.variables[varNames[iCol]]
            plt.subplot(nCol, nRow, iRow * nRow + iCol + 1)
            iEnd=5*24
            for iLev in range(len(levels)):
                k = levels[iLev]
                plt.plot(daysSinceStartOfSim[0:iEnd],var[0:iEnd,iCell,k],label='k='+str(k),color=colors[iLev])
            plt.legend()
            plt.title(varNames[iCol]+ ' at lon,lat='+str(lonMinDeg) +', '+ str(latMinDeg)+' f='+str(round(fCyclesPerDay,3)))
            plt.grid()
            if iRow == nRow - 1:
                plt.ylabel('velocity, m/s')
            if iCol == nCol - 1:
                plt.xlabel('time [days]')
    
            # print('plotting: '+varNames[iRow] + titleTxt[iCol])
            # print(' min: ',np.min(var[timeArray[iCol], :]), ' max: ',np.max(var[timeArray[iCol], :]))
            # print(var[timeArray[iCol], 1:50])
            # print('np.shape(var)',np.shape(var))
            # varSliced = np.transpose(var[timeArray[iCol], :, :])
            # varMasked = ma.masked_where(varSliced < -1.0e33, varSliced)
            # ax = plt.imshow(varMasked)  # ,extent=[yMinkm,yMaxkm,zMin,zMax])
            # plt.axis('off')
    
    print('save plot')
    plt.savefig('Velocity_in_time_'+str(lonMinDeg)+ '_'+str(latMinDeg)+'.png')
