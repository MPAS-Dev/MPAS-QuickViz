#!/usr/bin/env python

import numpy as np
import xarray as xr
import glob
import rossby_radius_deformation

def compute_phase_speed():
  files = np.sort(glob.glob('lagrPartTrack.*nc'))
  files = [fn.replace('lagrPartTrack.','') for fn in files]
  outputfiles = ['../../output/output.' + fn for fn in files]

  ds = xr.open_mfdataset(outputfiles[:3], concat_dim='Time')
  yCell = ds.yCell.mean('Time')
  ds = ds.drop('yCell')
  ds = ds.merge({'yCell':yCell})
  ds = ds.set_coords('yCell')

  # need to compute depth mean flow because cw = U_{mean=zt} - beta LD^2
  # note that all are collapsed to yCell coordinates
  depthmeanflow = ds.velocityZonal.mean(['Time','nVertLevels'])\
                                  .groupby('yCell')\
                                  .mean()

  fmean = ds.fCell.groupby('yCell').mean()
  beta = fmean.diff('yCell')/ \
         ds.yCell.groupby('yCell').mean().diff('yCell')

  # note yet averaged to y
  ztop = ds.zTop.mean('Time')
  N = np.sqrt(ds.BruntVaisalaFreqTop.mean('Time'))
  omega = omega = 2*np.pi/86400.
  latEff = np.arcsin(fmean/(2*omega))*180.0/np.pi

  # need to ensure that N is not negative, if so, average top
  # layers until it is no longer negative to simulate the mixed layer




  import pdb; pdb.set_trace()


if __name__ == "__main__":
  compute_phase_speed()
