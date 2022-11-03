#!/usr/bin/env python
"""

    File handling module.
    
    Phillip Wolfram
    LANL
    03/03/2015

"""
# import libraries / packages
import numpy as np
import os
import sys
import socket
import netCDF4
import datetime

def calling_metadata(): #{{{
    """ returns the current directory, the command used to call script, and the host name of the machine  """
    return os.getcwd(), ' '.join(sys.argv), socket.gethostname() #}}}

def add_meta(ds): #{{{
    """ adds meta data to xray data set """
    ds.attrs['meta_cwd'] = os.getcwd()
    ds.attrs['meta_host'] = socket.gethostname()
    ds.attrs['meta_call'] = ' '.join(sys.argv)
    ds.attrs['meta_user'] = os.getenv('USER')
    ds.attrs['meta_time'] = datetime.datetime.now()
    return ds #}}}

def savenpz(filename, *args, **kwargs): #{{{
    """ wrapper for np.savez that stores meta data of how file was produced in terms of script call, current working directory, and host name """
    cwd, call, host = calling_metadata()
    np.savez(filename, meta_call=call, meta_cwd=cwd, meta_host=host, *args, **kwargs)
    return #}}}

def num_ensembles(rootdir, folder_prefix): #{{{
    """ 
    Compute number of realizations for realization folder prefix 'folder_prefix' in 'rootdir'
    Phillip Wolfram, LANL
    09/19/2014
    """
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
    Get paths for realization folders begining with 'folder_prefix' in 'rootdir'.
    Phillip Wolfram, LANL
    09/19/2014
    """
    fnames = []
    for file in os.listdir(rootdir):
          if(file.startswith(folder_prefix)):
            fnames.append(rootdir + '/' + file)
    return fnames #}}}

def get_file(folder, file_prefix): #{{{
    """
    Get file in 'folder' starting with 'file_prefix'.
    Phillip Wolfram, LANL
    09/19/2014
    """
    for file in os.listdir(folder):
        if(file.startswith(file_prefix)):
            return folder + '/' + file #}}}

def check_folder(folder): #{{{
    """
    Check to see if 'folder' exists.  If not, create it.
    Phillip Wolfram, LANL
    09/19/2014
    """
    if not os.path.isdir(folder):
        os.mkdir(folder) #}}}

def file_exists(fname): #{{{
    return os.path.isfile(fname) #}}}

def open_files(folder_names, file_prefix): #{{{
    """
    Open all files in 'folder_names' begining with
    'file_prefix' and return in list of opened netCDF4 databases.
    Phillip Wolfram, LANL
    09/19/2014
    """
    files = []
    for folder in folder_names:
        files.append(netCDF4.Dataset(get_file(folder, file_prefix),'r'))
    return files #}}}

def close_files(files): #{{{
    """
    Close all netCDF4 files in list 'files'.
    Phillip Wolfram, LANL
    09/19/2014
    """
    for file in files:
        file.close()
    return #}}}

if __name__ == "__main__":
    print 'Module with functions for file handling, particularly of realization data.'

