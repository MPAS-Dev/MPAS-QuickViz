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
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import spatial
plt.switch_backend('agg')
cartopy.config['pre_existing_data_dir'] = \
        os.getenv('CARTOPY_DIR', cartopy.config.get('pre_existing_data_dir'))


################################################################################################
################################################################################################

def read_pointstats(pointstats_file):

  pointstats_nc = netCDF4.Dataset(pointstats_file,'r')
  data = {}
  data['date'] = pointstats_nc.variables['xtime'][:]
  data['datetime'] = []
  for date in data['date']:
    d = b''.join(date).strip()
    data['datetime'].append(datetime.datetime.strptime(d.decode('ascii').strip('\x00'),'%Y-%m-%d_%H:%M:%S'))
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

def read_station_file(station_file,stations={}):

  # Initialize stations dictionary
  if len(stations) == 0:
    stations['name'] = []
    stations['lon'] = []
    stations['lat'] = []

  # Read in stations names and location
  f = open(station_file)
  lines = f.read().splitlines()
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

  # Read config file
  f = open(pwd+'/plot_ssh.config')
  cfg = yaml.load(f)
  pprint.pprint(cfg)

  # Read in model point output data and create kd-tree
  data = {}
  tree = {}

  for run in cfg['pointstats_file']:
    data[run] = read_pointstats(cfg['pointstats_file'][run])
    points = np.vstack((data[run]['lon'],data[run]['lat'])).T
    tree[run] = spatial.KDTree(points)

  # Read in station file
  stations = read_station_file(cfg['stations_file'])

  for i,sta in enumerate(stations['name']):
    print(sta)

    # Check if observation file exists
    obs_file = ""
    obs_file_check = cfg['obs_direc']+sta+'_'+cfg['year']+'.txt'
    if os.path.isfile(obs_file_check):
      obs_file = obs_file_check

    obs_file_check = cfg['obs_direc']+sta+'.txt'
    if os.path.isfile(obs_file_check):
      obs_file = obs_file_check

    # Skip to next iteration if not found
    if not obs_file:
      continue

    # Read in observed data and get coordinates
    obs_data = read_station_data(obs_file,cfg['min_date'],cfg['max_date'])
    sta_lon = stations['lon'][i]
    sta_lat = stations['lat'][i]

    # Create figure
    fig = plt.figure(figsize=[6,4])
    gs = gridspec.GridSpec(nrows=2,ncols=2,figure=fig)

    # Plot observation station location
    ax1 = fig.add_subplot(gs[0,0], projection=ccrs.PlateCarree())
    ax1.set_extent([sta_lon-10.0, sta_lon+10.00, sta_lat-7.0 , sta_lat+7.0], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.LAND, zorder=100)
    ax1.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
    ax1.coastlines('50m', zorder=101)
    ax1.plot(sta_lon,sta_lat,'C0o', zorder=102)

    # Plot local observation station location
    ax2 = fig.add_subplot(gs[0,1], projection=ccrs.PlateCarree())
    ax2.set_extent([sta_lon-2.5, sta_lon+2.5, sta_lat-1.75 , sta_lat+1.75], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND, zorder=100)
    ax2.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
    ax2.coastlines('50m', zorder=101)
    ax2.plot(sta_lon,sta_lat,'C0o', zorder=102)

    # Plot observed data
    ax3 = fig.add_subplot(gs[1,:])
    l1, = ax3.plot(obs_data['datetime'],obs_data['ssh'],'C0-')
    labels = ['Observed']
    lines = [l1]

    for i,run in enumerate(data):

      # Find closest output point to station location
      d,idx = tree[run].query(np.asarray([sta_lon,sta_lat]))

      # Plot output point location
      ax1.plot(data[run]['lon'][idx],data[run]['lat'][idx],'C'+str(i+1)+'o')
      ax2.plot(data[run]['lon'][idx],data[run]['lat'][idx],'C'+str(i+1)+'o')

      # Plot modelled data
      l2, = ax3.plot(data[run]['datetime'],data[run]['ssh'][:,idx],'C'+str(i+1)+'-')
      labels.append(run)
      lines.append(l2)

    # Set figure labels and axis properties and save
    ax3.set_xlabel('time')
    ax3.set_ylabel('ssh (m)')
    ax3.set_xlim([datetime.datetime.strptime(cfg['min_date'],'%Y %m %d %H %M'),datetime.datetime.strptime(cfg['max_date'],'%Y %m %d %H %M')])
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    lgd = plt.legend(lines,labels,loc=9,bbox_to_anchor=(0.5,-0.5),ncol=3,fancybox=False,edgecolor='k')
    st = plt.suptitle('Station '+sta,y = 1.025,fontsize=16)
    fig.tight_layout()
    fig.savefig(sta+'.png',bbox_inches='tight',bbox_extra_artists=(lgd,st,))
    plt.close()
