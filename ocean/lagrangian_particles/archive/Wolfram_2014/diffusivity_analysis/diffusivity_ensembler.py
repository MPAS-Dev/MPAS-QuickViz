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
rEarth = 6371220. 
# import libraries / packages
import numpy as np
import os 
import cPickle as pickle
from diffusivity_realization_data import normalized_haversine_formula

def aggregate_com(rootdir, folder_prefix):  #{{{
    """
    For all analyze_* folders in rootdir, load 
    com_time pickles and aggregate the data

    Phillip Wolfram 
    LANL
    07/22/2014
    """

    # make ensembles dir
    if not os.path.isdir(rootdir+'/ensembles'):
      os.mkdir(rootdir+'/ensembles')

    # get the number of realizations
    numRealizations = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for dir in dirs:
          # for each folder starting with the prefix
          if(dir.startswith(folder_prefix)):
            print dir
            numRealizations += 1
    print 'numRealizations = ', numRealizations

    ensembleNum = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for dir in dirs:
          # for each folder starting with the prefix
          if(dir.startswith(folder_prefix)):
            fname = rootdir + '/' + dir + '/com_time.p'
            print 'opening %s' % fname
            com = pickle.load(open(fname, "rb"))
            print 'done'
            if ensembleNum == 0:
              com_summed = com[:]
            else:
              com_summed = [a + b for a, b in zip(com_summed, com)]
            ensembleNum += 1

            com_avg = [a/ensembleNum for a in com_summed]

            #import pdb; pdb.set_trace()
            # save the data
            pickle.dump(com_summed, open(rootdir+'/ensembles/com_summed_'+str(ensembleNum)+'.p','wb'))
            pickle.dump(com_avg, open(rootdir+'/ensembles/com_avg_'+str(ensembleNum)+'.p','wb'))
            print 'ensembleNum = ', ensembleNum


    return  #}}}

def aggregate_dispersion(rootdir, folder_prefix):  #{{{
    """
    For all analyze_* folders in rootdir, load 
    com_time and dl_time pickles and aggregate the data

    Phillip Wolfram 
    LANL
    07/22/2014
    """

    # make ensembles dir
    if not os.path.isdir(rootdir+'/ensembles'):
      os.mkdir(rootdir+'/ensembles')

    # get the number of realizations
    numRealizations = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for dir in dirs:
          # for each folder starting with the prefix
          if(dir.startswith(folder_prefix)):
            print dir
            numRealizations += 1
    print 'numRealizations = ', numRealizations

    ensembleNum = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for dir in dirs:
          # for each folder starting with the prefix
          if(dir.startswith(folder_prefix)):
            fname = rootdir + '/' + dir + '/dl_time.p'
            print 'loading file % s' % fname
            # a ton of time is spent loading in the data
            dl = pickle.load(open(fname, "rb"))
            print 'done loading file'
            if ensembleNum == 0:
              drr_summed = np.zeros((len(dl), len(dl[0])))
              dxx_summed = np.zeros((len(dl), len(dl[0])))
              dxy_summed = np.zeros((len(dl), len(dl[0])))
              dyy_summed = np.zeros((len(dl), len(dl[0])))
              dl_count = np.zeros((len(dl), len(dl[0])))
            for i in np.arange(len(dl)):
              print 'time level %d of %d' %  (i, len(dl))
              for j in np.arange(len(dl[0])):
                #print 'cluster %d of %d' % (j, len(dl[0]))
                dl_count[i,j] +=  dl[i][j].shape[1]
                dx = dl[i][j][0]
                dy = dl[i][j][1]
                dr = np.sqrt(dx*dx + dy*dy)
                drr_summed[i,j] += sum(dr*dr)
                dxx_summed[i,j] += sum(dx*dx)
                dxy_summed[i,j] += sum(dx*dy)
                dyy_summed[i,j] += sum(dy*dy)
            ensembleNum += 1

            # get total dispersion using unbiased estimator
            # particle distance from center of mass
            dl_drr  = drr_summed/(dl_count-1)
            dl_dxx  = dxx_summed/(dl_count-1)
            dl_dxy  = dxy_summed/(dl_count-1)
            dl_dyy  = dyy_summed/(dl_count-1)

            # save the data
            print 'saving summed data'
            pickle.dump(drr_summed, open(rootdir+'/ensembles/dr2_summed_'+str(ensembleNum)+'.p','wb'))
            pickle.dump(dxx_summed, open(rootdir+'/ensembles/dxx_summed_'+str(ensembleNum)+'.p','wb'))
            pickle.dump(dxy_summed, open(rootdir+'/ensembles/dxy_summed_'+str(ensembleNum)+'.p','wb'))
            pickle.dump(dyy_summed, open(rootdir+'/ensembles/dyy_summed_'+str(ensembleNum)+'.p','wb'))
            print 'saving count data'
            pickle.dump(dl_count, open(rootdir+'/ensembles/dl_count_'+str(ensembleNum)+'.p','wb'))
            print 'saving dispersion data'
            pickle.dump(dl_drr, open(rootdir+'/ensembles/dl_drr_'+str(ensembleNum)+'.p','wb'))
            pickle.dump(dl_dxx, open(rootdir+'/ensembles/dl_dxx_'+str(ensembleNum)+'.p','wb'))
            pickle.dump(dl_dxy, open(rootdir+'/ensembles/dl_dxy_'+str(ensembleNum)+'.p','wb'))
            pickle.dump(dl_dyy, open(rootdir+'/ensembles/dl_dyy_'+str(ensembleNum)+'.p','wb'))
            print 'ensembleNum = ', ensembleNum

    return  #}}}

  
def diff_interval_start_end(disp, tstart, tend, tint): #{{{
    return 0.5*(disp[tend] - disp[tstart])/((tend-tstart)*tint*86400.) #}}}

def output_kappa(kappa, kappa_xx, kappa_xy, kappa_yy, ts, te, tint, com, sigma2_xx, sigma2_xy, sigma2_yy, interpLoc, folder): #{{{
  """
  Recognizing kappa is a moving timeseries (com) over time, use points
  between ts and te (indicies in time_vec) as data for a Gaussian weighted sum with dimensions given by
  sigma2's at interpLoc points

  Phillip Wolfram
  LANL
  07/24/2014
  """
  Nclusters = len(kappa)
  interval = 1
  npoints = 0
  for acluster in np.arange(Nclusters):
    for atime in np.arange(ts[acluster],te[acluster]+1,interval):
      npoints +=1
  # allocations #{{{
  kappa_tot = np.zeros(npoints)
  kappa_xx_tot = np.zeros(npoints)
  kappa_xy_tot = np.zeros(npoints)
  kappa_yy_tot = np.zeros(npoints)
  x_tot     = np.zeros(npoints)
  y_tot     = np.zeros(npoints)
  z_tot     = np.zeros(npoints)
  sigma2_xx_tot = np.zeros(npoints)
  sigma2_xy_tot = np.zeros(npoints)
  sigma2_yy_tot = np.zeros(npoints)
  clusternum_tot = np.zeros(npoints)
  #}}}
  npoints = 0
  ncluster = 0
  for acluster in np.arange(Nclusters):
    ncluster += 1
    for atime in np.arange(ts[acluster],te[acluster]+1,interval):
      kappa_tot[npoints:npoints+1] = kappa[acluster]
      kappa_xx_tot[npoints:npoints+1] = kappa_xx[acluster]
      kappa_xy_tot[npoints:npoints+1] = kappa_xy[acluster]
      kappa_yy_tot[npoints:npoints+1] = kappa_yy[acluster]
      x_tot[npoints:npoints+1]     = com[atime][0,acluster]
      y_tot[npoints:npoints+1]     = com[atime][1,acluster]
      sigma2_xx_tot[npoints:npoints+1] = sigma2_xx[atime][acluster]
      sigma2_xy_tot[npoints:npoints+1] = sigma2_xy[atime][acluster]
      sigma2_yy_tot[npoints:npoints+1] = sigma2_yy[atime][acluster]
      clusternum_tot[npoints:npoints+1] = ncluster
      npoints+=1

  # write ascii files for analysis in fortran
  print 'writing ascii data to ', folder
  with file(folder + '/cluster_data.txt', 'w') as outfile:
    # want the total number of data points
    outfile.write('%d\n' % x_tot.shape[0])
    for n,x,y,k,kxx,kxy,kyy,sxx,sxy,syy in zip(clusternum_tot, x_tot, y_tot, kappa_tot, kappa_xx_tot, kappa_xy_tot, kappa_yy_tot, sigma2_xx_tot, sigma2_xy_tot, sigma2_yy_tot):
      outfile.write('%d %e %e %e %e %e %e %e %e %e\n' % (n,x,y,k,kxx,kxy,kyy,sxx,sxy,syy))
  with file(folder + '/interp_points.txt', 'w') as outfile:
    outfile.write('%d\n' % len(interpLoc[0,:]))
    for x,y in zip(interpLoc[0,:],interpLoc[1,:]):
      outfile.write('%e %e\n' % (x,y))
  print 'done'

  return  #}}}

def write_output(ensembleDir, ensembleNum, folder, ts, te, tint):  #{{{
   """
   write output for use with Fortran guassian interpolation code

   Phillip Wolfram, LANL
   09/19/2014
   """
   if not os.path.isdir(ensembleDir+'/'+folder):
       os.mkdir(ensembleDir+'/'+folder)

   # load data for ensemble number
   print 'loading dispersion data'
   drr = pickle.load(open(ensembleDir+'/dl_drr_'+str(ensembleNum)+'.p','rb'))
   dxx = pickle.load(open(ensembleDir+'/dl_dxx_'+str(ensembleNum)+'.p','rb'))
   dxy = pickle.load(open(ensembleDir+'/dl_dxy_'+str(ensembleNum)+'.p','rb'))
   dyy = pickle.load(open(ensembleDir+'/dl_dyy_'+str(ensembleNum)+'.p','rb'))
   print 'loading COM data'
   com  = pickle.load(open(ensembleDir+'/com_avg_'+str(ensembleNum)+'.p','rb'))
   print 'loading num particles in ensemble'
   nparticles = pickle.load(open(ensembleDir+'/dl_count_'+str(ensembleNum)+'.p','rb'))
   print 'done'
   
   print 'computing kappa'
   kappa_simple =    diff_interval_start_end(drr, ts, te, tint)
   kappa_xx_simple = diff_interval_start_end(dxx, ts, te, tint)
   kappa_xy_simple = diff_interval_start_end(dxy, ts, te, tint)
   kappa_yy_simple = diff_interval_start_end(dyy, ts, te, tint)
   print 'done'

   # now output data for use with interpolation via guassian kernels in fortran code
   eones = np.ones(len(dxx[0]),dtype='int32')
   output_kappa(kappa_simple, kappa_xx_simple, kappa_xy_simple, kappa_yy_simple, ts*eones, te*eones, tint, com, dxx, dxy, dyy, com[0], ensembleDir + '/' + folder)
   print 'finished kappa calculation and outputed to file at ' + folder

   return #}}}

def dispersion_and_com_restart(rootdir, folder_prefix, newEnsembleNum): #{{{
    """
    dispersion_and_com_restart(rootdir, folder_prefix, newEnsembleNum):

    Incremental restart for dispersion calculation

    Phillip Wolfram 
    LANL
    07/22/2014
    """

    # make ensembles dir
    if not os.path.isdir(rootdir+'/ensembles'):
      os.mkdir(rootdir+'/ensembles')

    # get the number of realizations
    numRealizations = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for dir in dirs:
          # for each folder starting with the prefix
          if(dir.startswith(folder_prefix)):
            print dir
            numRealizations += 1
    print 'numRealizations = ', numRealizations

    ensembleNum = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for dir in dirs:
          # for each folder starting with the prefix
          if(dir.startswith(folder_prefix)):
            ensembleNum += 1
            if ensembleNum == newEnsembleNum :
              # this is the one that needs added
              fname = rootdir + '/' + dir + '/dl_time.p'
              print 'loading file % s' % fname
              dl = pickle.load(open(fname, "rb"))
              print 'done loading file'
              fname = rootdir + '/' + dir + '/com_time.p'
              print 'loading file % s' % fname
              com = pickle.load(open(fname, "rb"))
              print 'done loading file'
              print 'loading dl summed pickle'
              dr2_summed = pickle.load(open(rootdir+'/ensembles/dr2_summed_'+str(ensembleNum-1)+'.p','rb'))
              dxx_summed = pickle.load(open(rootdir+'/ensembles/dxx_summed_'+str(ensembleNum-1)+'.p','rb'))
              dxy_summed = pickle.load(open(rootdir+'/ensembles/dxy_summed_'+str(ensembleNum-1)+'.p','rb'))
              dyy_summed = pickle.load(open(rootdir+'/ensembles/dyy_summed_'+str(ensembleNum-1)+'.p','rb'))
              print 'loading dl count pickle'
              dl_count  = pickle.load(open(rootdir+'/ensembles/dl_count_'+str(ensembleNum-1)+'.p','rb'))
              print 'loading com summed pickle'
              com_summed = pickle.load(open(rootdir+'/ensembles/com_summed_'+str(ensembleNum-1)+'.p','rb'))

              for i in np.arange(len(dl)):
                print 'time level %d of %d' %  (i, len(dl))
                for j in np.arange(len(dl[0])):
                  #print 'cluster %d of %d' % (j, len(dl[0]))
                  dl_count[i,j] +=  dl[i][j].shape[1]
                  dx = dl[i][j][0]
                  dy = dl[i][j][1]
                  dr = np.sqrt(dx*dx + dy*dy)
                  dr2_summed[i,j] += sum(dr*dr)
                  dxx_summed[i,j] += sum(dx*dx)
                  dxy_summed[i,j] += sum(dx*dy)
                  dyy_summed[i,j] += sum(dy*dy)
              # particle distance from center of mass
              dl_disp = dr2_summed/(dl_count-1)
              dl_dxx  = dxx_summed/(dl_count-1)
              dl_dxy  = dxy_summed/(dl_count-1)
              dl_dyy  = dyy_summed/(dl_count-1)

              com_summed = [a + b for a, b in zip(com_summed, com)]
              com_avg = [a/ensembleNum for a in com_summed]

              # save the data
              print 'saving summed data'
              pickle.dump(dr2_summed, open(rootdir+'/ensembles/dr2_summed_'+str(ensembleNum)+'.p','wb'))
              pickle.dump(dxx_summed, open(rootdir+'/ensembles/dxx_summed_'+str(ensembleNum)+'.p','wb'))
              pickle.dump(dxy_summed, open(rootdir+'/ensembles/dxy_summed_'+str(ensembleNum)+'.p','wb'))
              pickle.dump(dyy_summed, open(rootdir+'/ensembles/dyy_summed_'+str(ensembleNum)+'.p','wb'))
              pickle.dump(com_summed, open(rootdir+'/ensembles/com_summed_'+str(ensembleNum)+'.p','wb'))
              print 'saving count data'
              pickle.dump(dl_count, open(rootdir+'/ensembles/dl_count_'+str(ensembleNum)+'.p','wb'))
              print 'saving dispersion data'
              pickle.dump(dl_disp, open(rootdir+'/ensembles/dl_disp_'+str(ensembleNum)+'.p','wb'))
              pickle.dump(dl_dxx, open(rootdir+'/ensembles/dl_dxx_'+str(ensembleNum)+'.p','wb'))
              pickle.dump(dl_dxy, open(rootdir+'/ensembles/dl_dxy_'+str(ensembleNum)+'.p','wb'))
              pickle.dump(dl_dyy, open(rootdir+'/ensembles/dl_dyy_'+str(ensembleNum)+'.p','wb'))
              print 'saving com data'
              pickle.dump(com_avg, open(rootdir+'/ensembles/com_avg_'+str(ensembleNum)+'.p','wb'))
              print 'ensembleNum = ', ensembleNum

    return  #}}}

def num_ensembles(rootdir, folder_prefix):
    # get the number of realizations
    numRealizations = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for dir in dirs:
          # for each folder starting with the prefix
          if(dir.startswith(folder_prefix)):
            numRealizations += 1
    return numRealizations

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
    parser.add_option("-e", "--end", dest="end",
                      help="step for ending time",
                      metavar="INTEGER")
    parser.add_option("-s", "--start", dest="start",
                      help="step for starting time",
                      metavar="INTEGER")
    parser.add_option("-i", "--interval", dest="interval",
                      help="time step interval",
                      metavar="INTEGER")

    options, args = parser.parse_args()

    if not options.root:
        parser.error("Root directory is a required input.")
    if not options.prefix:
        options.prefix = "analyze_output"
    if not options.end:
        parser.error("Ending time step is a required input.")
    else:
        options.end= int(options.end)
    if not options.start:
        parser.error("Starting time step is a required input.")
    else:
        options.start = int(options.start)
    if not options.interval:
        parser.error("Time interval is a required input.")
    else:
        options.interval = int(options.interval)

    #aggregate_com(options.root, options.prefix)
    #aggregate_dispersion(options.root, options.prefix)

    for i in (1 + np.arange(num_ensembles(options.root,options.prefix))):
        folder = '/diffusivityOutput'+ str(i) + '/'
        write_output('ensembles/', i, folder, options.start, options.end, options.interval)


