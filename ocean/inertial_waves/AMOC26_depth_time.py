import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
import xarray
import netCDF4 as nc
import glob
from datetime import date

mpl.rcParams['figure.figsize'] = (8,12) # Large figures
mpl.rcParams['image.cmap'] = 'seismic'
#mpl.rcParams['image.cmap'] = 'Spectral'
dpi=200;

#wd = '/lustre/scratch4/turquoise/mpeterse/runs/'
#rundir = '220621_EC30to60_moc01a/'
#subdir = 'ocean/global_ocean/EC30to60/PHC/dynamic_adjustment/simulation/'
#mesh = xarray.open_dataset(wd+rundir+subdir+'init.nc')

mesh = xarray.open_dataset('init.nc')
# compy
month = '0030-01-01'
rundir = '20220615.testG-MLE-bgrad0p75KPP.compy.TL319_EC30to60E2r2.compy'
shortName='Gcase EC'
filename = 'link/20220615.testG-MLE-bgrad0p75KPP.compy.TL319_EC30to60E2r2.compy.mpaso.hist.am.mocStreamfunctionOutput.'+month+'.nc'
filenameTime = 'link/20220615.testG-MLE-bgrad0p75KPP.compy.TL319_EC30to60E2r2.compy.mpaso.hist.am.globalStats.'+month+'.nc'
#flist2 = sorted(glob.glob('/lcrc/group/e3sm/ac.abarthel/E3SMv2/20220517.v2.amoc.GMPASO-CORE.3Drest.T62_oRRS30to10v3.anvil/run/*mpaso.hist.am.timeSeriesStatsMonthly.0013*')) # Monthly averages

#flist = sorted(glob.glob('analysis_members/mocStreamfunction.*.nc'))
#nt=0
ds = xarray.open_dataset(filename)

f = nc.Dataset(filenameTime, 'r')
f.variables['daysSinceStartOfSim'].set_auto_scale(False)
daysSinceStartOfSim = f.variables['daysSinceStartOfSim'][:]

day = 0 #Hovmuller for day 8
nt = len(ds.xtime.values)
t1 = 0#(day-1)*24
t2 = nt
cmax=20
levels = np.linspace(-cmax,cmax,22)
ax=plt.subplot(2,1,1)
ca=plt.contourf(np.arange(t1,t2)/24,-mesh.refBottomDepth,ds.mocStreamvalLatAndDepthRegion[t1:t2,0,:,127].T,
   levels=levels,extend='both')
#ax.set_xlim(40,44.5)
plt.xlabel('time [days]')
plt.ylabel('depth [m]')
plt.title('AMOC streamfunction [Sv] 26.5N '+shortName)
plt.colorbar(ticks=np.linspace(-cmax,cmax,9))

ax=plt.subplot(2,1,2)
k=44
ca=plt.plot(np.arange(t1,t2)/24,ds.mocStreamvalLatAndDepthRegion[t1:t2,0,:,127].max(axis=1))
plt.grid()
#ax.set_ylim(0,20)
#ax.set_xlim(40,44.5)
plt.xlabel('time [days]')
plt.ylabel('transport at 1000m depth, 26.5N [Sv] ')
plt.title('AMOC streamfunction [Sv] 26.5N'+shortName)

today = date.today()
plt.figtext(.03,.05,
   today.strftime("%b-%d-%Y") +'\n'+
   rundir +'\n'+
   'plot starts at: '+str(ds.xtime.values[0]))
plt.savefig('f/moc_26N.png')


#t1=200
#t2=224
#print( ds.mocStreamvalLatAndDepthRegion[0:1,0,:,112] )
##for ti in range(t1,t2):
#plt.plot(ds.mocStreamvalLatAndDepthRegion[t1:t2,0,:,112].T,-mesh.refBottomDepth)
#plt.savefig('f/moc_lines.png')

# for max over depth:
#plt.plot(ds.mocStreamvalLatAndDepthRegion[s1:s2,0,:,112].max(axis=1))

