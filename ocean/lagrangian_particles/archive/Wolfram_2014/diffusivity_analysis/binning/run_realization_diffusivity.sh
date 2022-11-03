#!/usr/bin/env bash


# can do 12 at a time for 5km grid
# reasonable call:
# ./run_realization_diffusivity.sh 0 0 10
if [ $# -eq 4 ]; then
  rlzn=`printf "%04d" $1`
  startnum=$2
  endnum=$3

  for i in $(seq $startnum $endnum); do
    startlayer=$i
    stoplayer=$((i + 1))
    echo 'doing analysis for realization '$rlzn' on layer range '$startlayer' to '$stoplayer
    stdbuf -oL \
      time python diffusivity_binning_complete.py -f 'all_rlzn_particle_data_rads_folder_50day/all_rlzn_particle_data_rads_rlzn'$rlzn'.nc' \
      -s mesh10km.nc  -n $4 -l "np.arange($startlayer,$stoplayer)" -d F -r 100000.  --frac 1.0 -p T \
      &> out_rlzn_${rlzn}_layers${startlayer}_${stoplayer}.out &
  done

else
  echo 'wrong number of arguments'
fi
