#!/bin/bash

SCRIPT=run_binning_complete_pdf.sh

#for i in 4km 8km 16km 32km ; do
for i in 4km; do
  cp *.py *.sh *.xml *.npz $i
  cd $i
  #rm -rf Case*_r*; rm -rf slurm-* *.log
  if [ -z "$SLURM_JOBID" ]; then
    msub $SCRIPT
  else
    source $SCRIPT
  fi
  cd ..
done
