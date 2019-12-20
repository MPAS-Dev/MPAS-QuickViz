#!/usr/bin/env python
"""
Name: compute_transects.py
Author: Phillip J. Wolfram, Mark Petersen, Luke Van Roekel

Computes transport through sections.

Example call:
  ./compute_transects.py
  -k transect_masks.nc
  -m MPAS_mesh.nc
  -t 'RUN_PATH/analysis_members/timeSeriesStatsMonthly.*.nc'
  -n 'all'

To create the transect_masks.nc file, load e3sm-unified and:
   MpasMaskCreator.x MPAS_mesh.nc  transect_masks.nc -f transect_definitions.geojson
where the transect_definitions.geojson file includes a sequence of lat/lon points for each transect.
"""

# ensure plots are rendered on ICC
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from netCDF4 import Dataset
import glob
import platform

m3ps_to_Sv = 1e-6 # m^3/sec flux to Sverdrups

def get_mask_short_names(mask):
  shortnames = [str(aname.values)[:str(aname.values).find(',')].strip()
                for aname in mask.transectNames]
  mask['shortNames'] = xr.DataArray(shortnames, dims='nTransects')
  mask = mask.set_index(nTransects=['transectNames', 'shortNames'])
  return mask

def compute_transport(timeavg, mesh, mask, name='Drake Passage',output='transport.nc'):
  mesh = xr.open_dataset(mesh)
  mask = get_mask_short_names(xr.open_dataset(mask))

  if name.lower() == 'all':
    transectList = mask.shortNames[:].values
    condition = transectList != "Atlantic Transec"
    transectList = np.extract(condition, transectList)
  else:
    transectList = name.split(',')
    if platform.release()[0] == '3':
      for i in range(len(transectList)):
        transectList[i] = "b'" + transectList[i]

  print('Computing Transport for the following transects ',transectList)
  nTransects = len(transectList)
  maxEdges = mask.dims['maxEdgesInTransect']
# create empty t list for time
  t = []
# Compute refLayerThickness to avoid need for hist file
  refBottom = mesh.refBottomDepth.values
  nz = mesh.dims['nVertLevels']
  h = np.zeros(nz)
  h[0] = refBottom[0]
  for i in range(1,nz):
    h[i] = refBottom[i] - refBottom[i-1]

# Get a list of edges and total edges in each transect
  nEdgesInTransect = np.zeros(nTransects)
  edgeVals = np.zeros((nTransects,maxEdges))
  for i in range(nTransects):
    amask = mask.sel(shortNames=transectList[i]).squeeze()
    transectEdges = amask.transectEdgeGlobalIDs.values
    inds = np.where(transectEdges > 0)[0]
    nEdgesInTransect[i] = len(inds)
    transectEdges = transectEdges[inds]
    edgeVals[i,:len(inds)] = np.asarray(transectEdges-1, dtype='i')

  nEdgesInTransect = np.asarray(nEdgesInTransect, dtype='i')

# Create a list with the start and stop for transect bounds
  nTransectStartStop = np.zeros(nTransects+1)
  for j in range(1,nTransects+1):
    nTransectStartStop[j] = nTransectStartStop[j-1] + nEdgesInTransect[j-1]

  edgesToRead = edgeVals[0,:nEdgesInTransect[0]]
  for i in range(1,nTransects):
    edgesToRead = np.hstack([edgesToRead,edgeVals[i,:nEdgesInTransect[i]]])

  edgesToRead = np.asarray(edgesToRead, dtype='i')
  dvEdge = mesh.dvEdge.sel(nEdges=edgesToRead).values
  edgeSigns = np.zeros((nTransects,len(edgesToRead)))
  for i in range(nTransects):
    edgeSigns[i,:] = mask.sel(nEdges=edgesToRead, shortNames=transectList[i]).squeeze().transectEdgeMaskSigns.values

# Read time average files one at a time and slice
  fileList = sorted(glob.glob(timeavg))
  transport = np.zeros((len(fileList),nTransects))
  t = np.zeros(len(fileList))
  for i,fname in enumerate(fileList):
    ncid = Dataset(fname,'r')
    if 'timeMonthly_avg_normalTransportVelocity' in ncid.variables.keys():
      vel = ncid.variables['timeMonthly_avg_normalTransportVelocity'][0,edgesToRead,:]
    elif 'timeMonthly_avg_normalVelocity' in ncid.variables.keys():
      vel = ncid.variables['timeMonthly_avg_normalVelocity'][0,edgesToRead,:]
      if 'timeMonthly_avg_normalGMBolusVelocity' in ncid.variables.keys():
        vel += ncid.variables['timeMonthly_avg_normalGMBolusVelocity'][0,edgesToRead,:]
    else:
      raise KeyError('no appropriate normalVelocity variable found')
    t[i] = ncid.variables['timeMonthly_avg_daysSinceStartOfSim'][:] / 365.
    ncid.close()
#   Compute transport for each transect
    for j in range(nTransects):
      start = int(nTransectStartStop[j])
      stop = int(nTransectStartStop[j+1])
      transport[i,j] = (dvEdge[start:stop,np.newaxis]*h[np.newaxis,:]*vel[start:stop,:] \
          *edgeSigns[j,start:stop,np.newaxis]).sum()*m3ps_to_Sv

# Define some dictionaries for transect plotting
  obsDict = {'Drake Passage':[120,175],'Tasmania-Ant':[147,167],'Africa-Ant':None,'Antilles Inflow':[-23.1,-13.7], \
          'Mona Passage':[-3.8,-1.4],'Windward Passage':[-7.2,-6.8],'Florida-Cuba':[30,33],'Florida-Bahamas':[30,33], \
          'Indonesian Throughflow':[-21,-11],'Agulhas':[-90,-50],'Mozambique Channel':[-20,-8], \
          'Bering Strait':[0.17,1.49],'Lancaster Sound':None,'Fram Strait':None,'Robeson Channel':None,'Nares Strait':None}
  labelDict = {'Drake Passage':'drake','Tasmania-Ant':'tasmania','Africa-Ant':'africaAnt','Antilles Inflow':'Antilles', \
          'Mona Passage':'monaPassage','Windward Passage':'windwardPassage','Florida-Cuba':'floridaCuba',\
             'Florida-Bahamas':'floridaBahamas', \
          'Indonesian Throughflow':'indonesia','Agulhas':'agulhas','Mozambique Channel':'mozambique', \
          'Bering Strait':'beringstrait','Lancaster Sound':'lancaster','Fram Strait':'fram','Robeson Channel':'robeson','Nares Strait':'nares'}

  for i in range(nTransects):
    plt.figure()
    if platform.release()[0]=='3':
        searchString = transectList[i][2:]
    else:
        searchString = transectList[i]
    bounds = obsDict[searchString]
    title = 'Transport for '+searchString
    plt.plot(t,transport[:,i],'k',linewidth=2)
    if bounds is not None:
        plt.gca().fill_between(t, bounds[0]*np.ones_like(t), bounds[1]*np.ones_like(t), alpha=0.3, label='observations')
    plt.ylabel('Transport (Sv)',fontsize=32)
    plt.xlabel('Time (Years)',fontsize=32)
    plt.title(title,fontsize=32)
    plt.savefig('transport_'+labelDict[searchString]+'.png')

# Add calls to save transport and then can build up
  ncid=Dataset(output,mode='w',clobber=True, format='NETCDF3_CLASSIC')
  ncid.createDimension('Time',None)
  ncid.createDimension('nTransects',nTransects)
  ncid.createDimension('StrLen',64)
  transectNames=ncid.createVariable('TransectNames','c',('nTransects','StrLen'))
  times=ncid.createVariable('Time','f8','Time')
  transportOut=ncid.createVariable('Transport','f8',('Time','nTransects'))

  times[:] = t
  transportOut[:,:] = transport

  for i in range(nTransects):
    nLetters = len(transectList[i])
    transectNames[i,:nLetters] = transectList[i]
  ncid.close()


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("-o", "--output_file_pattern", dest="output_filename_pattern",
      help="MPAS Filename pattern for transport output.", metavar="NAME")
  parser.add_argument("-t", "--time_avg_file_pattern", dest="time_avg_filename_pattern",
      help="MPAS Filename pattern for time averaged AM output.", metavar="FILE",
      required=True)
  parser.add_argument("-m", "--mesh_file", dest="mesh_filename",
      help="MPAS Mesh filename.", required=True)
  parser.add_argument("-k", "--mask_file", dest="mask_filename",
      help="MPAS mask filename.", required=True)
  parser.add_argument("-n", "--name", dest="name",
      help="List of transect names for computation, or 'all'", metavar="NAME")
  args = parser.parse_args()

  compute_transport(timeavg=args.time_avg_filename_pattern,
      mesh=args.mesh_filename, mask=args.mask_filename, name=args.name,
      output=args.output_filename_pattern)
