#!/bin/bash

for i in 4km 4km_2dx 4km_8dx 8km 16km 32km 32km_8dx; do
  cp *.py *.sh *.xml *.npy *.npz $i
  cd $i
  #rm -rf Case*_r*; rm -rf slurm-* *.log
  msub run_binning_complete.sh
  cd ..
done
