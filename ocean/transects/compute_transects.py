#!/usr/bin/env python
"""
Name: compute_transects.py
Author: Phillip J. Wolfram, Mark Petersen

Computes transport through sections.

Example call:
  ./compute_transects.py
  -k /lustre/scratch3/turquoise/mpeterse/runs/c62n/ocean/global_ocean/EC_60to30km/spin_up/init_step2/EC60to30v3_transect_masks.nc
  -m /lustre/scratch2/turquoise/mpeterse/runs/c69z/init.nc
  -o '/lustre/scratch2/turquoise/mpeterse/runs/c69z/output/output.*.nc'
  -t '/lustre/scratch2/turquoise/mpeterse/runs/c69z/analysis_members/timeSeriesStatsMonthly.*.nc'
  -n 'Drake Passage'

"""

# ensure plots are rendered on ICC
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

m3ps_to_Sv = 1e-6 # m^3/sec flux to Sverdrups

def get_mask_short_names(mask):
  shortnames = [str(aname.values)[:str(aname.values).find(',')].strip()
                for aname in mask.transectNames]
  mask['shortNames'] = xr.DataArray(shortnames, dims='nTransects')
  mask = mask.set_index(nTransects=['transectNames', 'shortNames'])
  return mask

def compute_transport(output, timeavg, mesh, mask, name='Drake Passage'):
  output = xr.open_mfdataset(output, concat_dim='Time', decode_times=False)
  timeavg = xr.open_mfdataset(timeavg, concat_dim='Time', decode_times=False)
  mesh = xr.open_dataset(mesh)
  mask = get_mask_short_names(xr.open_dataset(mask))

  amask = mask.sel(shortNames=name).squeeze()
  amask = amask.where(amask.transectEdgeGlobalIDs > 0, drop=True)
  transectEdges = np.asarray(amask.transectEdgeGlobalIDs.values-1, dtype='i')

  edgesigns = mask.sel(nEdges=transectEdges, shortNames=name).squeeze().transectEdgeMaskSigns
  dvEdge = mesh.dvEdge.sel(nEdges=transectEdges)
  vel = timeavg.sel(nEdges=transectEdges).timeMonthly_avg_normalTransportVelocity
  h = output.refLayerThickness.sel(Time=0).squeeze()

  transport = (dvEdge*h*vel*edgesigns).sum(['nEdges','nVertLevels'])*m3ps_to_Sv
  t = timeavg.timeMonthly_avg_daysSinceStartOfSim.values/365.

  plt.figure()
  plt.plot(t, transport, 'k-', lw=3, label='monthly-averaged transport')
  plt.gca().fill_between(t, (134.-14.)*np.ones_like(t), (134.+14.)*np.ones_like(t), alpha=0.3, label='observations')
  plt.legend(loc='best', frameon=False)
  plt.ylabel('Transport (Sv)')
  plt.xlabel('yrs')
  plt.title('Transport for ' + name)
  plt.savefig('transport' + name.replace(' ', '') +'.png')


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("-o", "--output_file_pattern", dest="output_filename_pattern",
      help="MPAS Filename pattern for output.", metavar="FILE",
      required=True)
  parser.add_argument("-t", "--time_avg_file_pattern", dest="time_avg_filename_pattern",
      help="MPAS Filename pattern for time averaged AM output.", metavar="FILE",
      required=True)
  parser.add_argument("-m", "--mesh_file", dest="mesh_filename",
      help="MPAS Mesh filename.", required=True)
  parser.add_argument("-k", "--mask_file", dest="mask_filename",
      help="MPAS mask filename.", required=True)
  parser.add_argument("-n", "--name", dest="name",
      help="Transect name for computation", metavar="NAME")
  args = parser.parse_args()

  compute_transport(output=args.output_filename_pattern,
      timeavg=args.time_avg_filename_pattern,
      mesh=args.mesh_filename, mask=args.mask_filename, name=args.name)
