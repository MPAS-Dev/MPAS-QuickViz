#!/usr/bin/env bash
# Phillip J. Wolfram
# 01/27/2016

# compute the masks for the dataset
for i in lagrPartTrack.00*; do
  time python build_masks.py -f $i -p T
done

# compute the diffusivity realization for the dataset
time python xarray_analysis.py -p 'lagrPartTrack.0031-0*nc' -m 'valid_realization_mask_0031-0*nc' -g mesh.nc -r 1.e5


