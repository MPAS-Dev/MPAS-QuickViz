#!/usr/bin/env python

import xarray as xr

def compute_diff(fullpath, highpath, lowpath, rnumsave, savefolder):

  print 'realization number %s'%(rnumsave)

  print 'loading files'
  full = xr.open_mfdataset(fullpath)
  high = xr.open_mfdataset(highpath)
  low  = xr.open_mfdataset(lowpath)

  print 'computing differential positions'
  # note, offset is necessary to make sure particle positions all start at same
  # place in time (it is ok to arbitrarily translate them)
  lon = full.lon - (high.lon + low.lon) + 2*full.sel(Time=0).lon
  lat = full.lat - (high.lat + low.lat) + 2*full.sel(Time=0).lat

  print 'saving dataset'
  # just use full values except for differenced values here
  xr.Dataset({'lon':lon, 'lat':lat, 'notoutcropped':full.notoutcropped, 'dtdays':full.dtdays}).to_netcdf(savefolder + 'all_rlzn_particle_data_rads_rlzn%04d.nc'%(rnumsave))

  print 'done'


if __name__ == "__main__":
    # Get command line parameters #{{{
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", "--full", dest="fullpath", help="path for advection with full velocity", metavar="FILE")
    parser.add_option("-e", "--high", dest="highpath", help="path for advection with high-passed velocity", metavar="FILE")
    parser.add_option("-m", "--low",  dest="lowpath",  help="path for advection with low-passed velocity", metavar="FILE")
    parser.add_option("-s", "--num", dest="rnumsave", help="number of realizations to save", metavar="STR")
    parser.add_option("-o", "--out", dest="savefolder", help="path for advection with inter velocity", metavar="STR")

    options, args = parser.parse_args()

    if not options.fullpath and not options.highpath and not options.lowpath:
      assert False, 'Must specify all files!'
    if not options.savefolder:
      options.savefolder = './'
    if not options.rnumsave:
      options.rnumsave = None
    else:
      options.rnumsave = int(options.rnumsave)
    #}}}
    compute_diff(options.fullpath, options.highpath, options.lowpath, options.rnumsave, options.savefolder)
