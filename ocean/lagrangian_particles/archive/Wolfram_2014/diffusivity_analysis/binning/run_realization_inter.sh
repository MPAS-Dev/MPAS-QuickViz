#!/usr/bin/env bash

if [ $# -eq 1 ]; then
  rlzn=`printf "%04d" $1`

  mkdir all_rlzn_particle_data_rads_folder_50day/

  echo 'doing analysis to compute particle differences for realization '$rlzn
  stdbuf -oL \
    time python compute_position_diff.py -s $rlzn \
    --full '../full_analysis/all_rlzn_particle_data_rads_folder_50day/all_rlzn_particle_data_rads_rlzn'$rlzn'.nc' \
    --high '../high_analysis/all_rlzn_particle_data_rads_folder_50day/all_rlzn_particle_data_rads_rlzn'$rlzn'.nc' \
    --low '../low_analysis/all_rlzn_particle_data_rads_folder_50day/all_rlzn_particle_data_rads_rlzn'$rlzn'.nc' \
    --out './all_rlzn_particle_data_rads_folder_50day/' \
    &> out_difference_rlzn_${rlzn}.out &

else
  echo 'wrong number of arguments'
fi
