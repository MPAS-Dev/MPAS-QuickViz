import netCDF4
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import numpy as np
import glob
import pprint
import datetime
import os
import yaml
import pprint
import subprocess
from mpl_toolkits.basemap import Basemap
from scipy import spatial 
plt.switch_backend('agg')
np.set_printoptions(threshold=np.nan)

pointstats_file = 'pointwiseStats.nc'
stations_file = 'NOAA-COOPS_stations/stations.txt'
obs_direc = 'NOAA-COOPS_stations/'
#stations_file = 'USGS_stations/stations.txt'
#obs_direc = 'USGS_stations/'
year = '2012'
min_date = '2012 10 20 00 00'
max_date = '2012 11 05 00 00'

################################################################################################
################################################################################################

def read_pointstats(pointstats_file):

  pointstats_nc = netCDF4.Dataset(pointstats_file,'r')
  data = {}
  data['date'] = pointstats_nc.variables['xtime'][:]
  data['datetime'] = []
  for date in data['date']:
    data['datetime'].append(datetime.datetime.strptime(''.join(date).strip(),'%Y-%m-%d_%H:%M:%S'))
  data['datetime'] = np.asarray(data['datetime'],dtype='O')
  data['lon'] = np.degrees(pointstats_nc.variables['lonCellPointStats'][:])
  data['lon'] = np.mod(data['lon']+180.0,360.0)-180.0
  data['lat'] = np.degrees(pointstats_nc.variables['latCellPointStats'][:])
  data['ssh'] = pointstats_nc.variables['sshPointStats'][:]

  return data  

################################################################################################
################################################################################################

def read_station_data(obs_file,min_date,max_date):

  frmt = '%Y %m %d %H %M'

  # Initialize variable for observation data
  obs_data = {}
  obs_data['ssh'] = []
  obs_data['datetime'] = []

  # Get data from observation file between min and max output times
  f = open(obs_file)
  obs = f.read().splitlines()
  for line in obs[1:]:
    if line.find('#') >= 0 or len(line.strip()) == 0 or not line[0].isdigit():
      continue
    try:                                                               # NOAA-COOPS format
      date = line[0:16]
      date_time = datetime.datetime.strptime(date,frmt)
      col = 5
      convert = 1.0
    except:                                                            # USGS station format
      date = line[0:19]
      date_time = datetime.datetime.strptime(date,'%m-%d-%Y %H:%M:%S')
      col = 2
      convert = 0.3048
    if date_time >= datetime.datetime.strptime(min_date,frmt) and \
       date_time <= datetime.datetime.strptime(max_date,frmt):
      obs_data['datetime'].append(date_time)
      obs_data['ssh'].append(line.split()[col])

  # Convert observation data and replace fill values with nan
  obs_data['ssh'] = np.asarray(obs_data['ssh'])
  obs_data['ssh'] = obs_data['ssh'].astype(np.float)*convert
  fill_val = 99.0
  obs_data['ssh'][obs_data['ssh'] >= fill_val] = np.nan

  obs_data['datetime'] = np.asarray(obs_data['datetime'],dtype='O')

  return obs_data

################################################################################################
################################################################################################

def read_station_file(station_file):

  # Read in stations names and location
  f = open(station_file)
  lines = f.read().splitlines()
  stations = {}
  stations['name'] = []
  stations['lon'] = []
  stations['lat'] = []
  for sta in lines:
    val = sta.split()
    stations['name'].append(val[2].strip("'"))
    stations['lon'].append(float(val[0]))
    stations['lat'].append(float(val[1]))
  nstations = len(stations['name'])
  stations['lon'] = np.asarray(stations['lon'])
  stations['lat'] = np.asarray(stations['lat'])

  return stations 

################################################################################################
################################################################################################

if __name__ == '__main__':

  pwd = os.getcwd()

  # Read in model point output data and create kd-tree 
  data = read_pointstats(pointstats_file)
  points = np.vstack((data['lon'],data['lat'])).T
  tree = spatial.KDTree(points)

  # Read in station file
  stations = read_station_file(stations_file)

  for i,sta in enumerate(stations['name']):
    print sta

    # Check if observation file exists
    obs_file = ""
    obs_file_check = obs_direc+sta+'_'+year+'.txt'
    if os.path.isfile(obs_file_check):
      obs_file = obs_file_check

    obs_file_check = obs_direc+sta+'.txt'
    if os.path.isfile(obs_file_check):
      obs_file = obs_file_check   

    # Skip to next iteration if not found
    if not obs_file:
      continue

    # Read in observed data
    obs_data = read_station_data(obs_file,min_date,max_date)

    # Find nearest model output point to station
    sta_lon = stations['lon'][i]
    sta_lat = stations['lat'][i]
    d,idx = tree.query(np.asarray([sta_lon,sta_lat]))

    # Create figure 
    fig = plt.figure(figsize=[6,4])
    gs = gridspec.GridSpec(nrows=2,ncols=2,figure=fig)

    # Plot station location
    ax = fig.add_subplot(gs[0,0])
    m = Basemap(projection='cyl',llcrnrlat=sta_lat-7.0 ,urcrnrlat=sta_lat+7.0,\
                                 llcrnrlon=sta_lon-10.0,urcrnrlon=sta_lon+10.0,resolution='l')
    m.fillcontinents(color='tan',lake_color='lightblue')
    m.drawcoastlines()
    ax.plot(sta_lon,sta_lat,'ro')
    ax.plot(data['lon'][idx],data['lat'][idx],'bo')

    # Plot local station location
    ax = fig.add_subplot(gs[0,1])
    m = Basemap(projection='cyl',llcrnrlat=sta_lat-1.75 ,urcrnrlat=sta_lat+1.75,\
                                 llcrnrlon=sta_lon-2.5,urcrnrlon=sta_lon+2.5,resolution='l')
    m.fillcontinents(color='tan',lake_color='lightblue')
    m.drawcoastlines()
    ax.plot(sta_lon,sta_lat,'ro')
    ax.plot(data['lon'][idx],data['lat'][idx],'bo')

    # Plot data
    ax = fig.add_subplot(gs[1,:])
    l1, = ax.plot(obs_data['datetime'],obs_data['ssh'],'r-')
    l2, = ax.plot(data['datetime'],data['ssh'][:,idx],'b-')
    ax.set_xlabel('time')
    ax.set_ylabel('ssh (m)')
    ax.set_xlim([datetime.datetime.strptime(min_date,'%Y %m %d %H %M'),datetime.datetime.strptime(max_date,'%Y %m %d %H %M')])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    lgd = plt.legend([l1,l2],['observed','modelled'],loc=9,bbox_to_anchor=(0.5,-0.5),ncol=2,fancybox=False,edgecolor='k')
    st = plt.suptitle('Station '+sta,y = 1.025,fontsize=16)
    fig.tight_layout()
    fig.savefig(sta+'.png',bbox_inches='tight',bbox_extra_artistis=(lgd,st,))
    plt.close()
