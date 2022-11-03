#!/usr/bin/env python
"""
    Ensemble realization data:
        dl_time.p
        com_time.p
    to produce ouput suitable for use with the Fortran Guassian Kernel smoother
    to produce a map of diffusivity
    
    Phillip Wolfram
    LANL
    09/19/2014

"""
# import libraries / packages
import numpy as np
import os 
import netCDF4

def num_ensembles(rootdir, folder_prefix): #{{{
    # get the number of realizations
    numRealizations = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for dir in dirs:
          # for each folder starting with the prefix
          if(dir.startswith(folder_prefix)):
            numRealizations += 1
    return numRealizations #}}}

def get_realization_paths(rootdir, folder_prefix): #{{{
    """
    Get paths for realization folders
    Phillip Wolfram, LANL
    09/19/2014
    """
    fnames = []
    for file in os.listdir(rootdir):
          if(file.startswith(folder_prefix)):
            fnames.append(rootdir + '/' + file)

    return fnames #}}}

def get_file(folder, file_prefix): #{{{
  for file in os.listdir(folder):
    if(file.startswith(file_prefix)):
      return folder + '/' + file #}}}

def check_folder(folder): #{{{
  if not os.path.isdir(folder):
      os.mkdir(folder) #}}}

def open_files(folder_names, file_prefix):
  files = []
  for folder in folder_names:
    files.append(netCDF4.Dataset(get_file(folder, file_prefix),'r'))
  return files

def close_files(files):
  for file in files:
    file.close()
  return

def aggregate_positions(rootdir, folder_prefix):  #{{{
    """
    For all analyze_* folders in rootdir, load 
    com_time pickles and aggregate the data

    Phillip Wolfram 
    LANL
    07/22/2014
    """

    folder_names = get_realization_paths(rootdir, folder_prefix)

    rlzns = open_files(folder_names, 'output')

    num_rlzns = len(folder_names)

    Nt = rlzns[0].variables['xParticle'].shape[0]
    Np = rlzns[0].variables['xParticle'].shape[1]
    lon = np.zeros((Nt,Np,num_rlzns))
    lat = np.zeros((Nt,Np,num_rlzns))

    for num,arlzn in enumerate(rlzns):
      print 'On realization %s of %s' % (num,num_rlzns)
      x = arlzn.variables['xParticle']
      y = arlzn.variables['yParticle']
      z = arlzn.variables['zParticle']

      latr, lonr = proj_lat_long(x[:,:], y[:,:], z[:,:])

      lon[:,:,num] = lonr[:,:]
      lat[:,:,num] = latr[:,:]

    lon = np.degrees(lon)
    lat = np.degrees(lat)

    #np.savez_compressed('all_rlzn_particle_data',lon=lon,lat=lat)
    np.savez('all_rlzn_particle_data',lon=lon,lat=lat)

    close_files(rlzns)

    return  #}}}

def proj_lat_long(x, y, z):  #{{{
    """
    compute the latitude and longitude from
    the x,y,z points (follow's Doug's pnt.h)
    """
    plat = np.arcsin(z / np.sqrt(x ** 2 + y ** 2 + z ** 2))
    plong = np.arctan2(y, x) 

    return plat, plong  #}}}

if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-r", "--root", dest="root",
                      help="folder root holding analysis folders",
                      metavar="FILE")
    parser.add_option("-p", "--prefix", dest="prefix",
                      help="folder prefix for analysis",
                      metavar="FILE")

    options, args = parser.parse_args()

    if not options.root:
      options.root = "."
    if not options.prefix:
      options.prefix = "analyze_"

    aggregate_positions(options.root, options.prefix)


