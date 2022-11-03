#!/usr/bin/env bash

setup_run() {
  jobdir='low_pass_realization_'$(printf "%02d" $1)'-'$(printf "%02d" $2)
  echo $jobdir
  mkdir $jobdir
  cd $jobdir
  mkdir restarts
  mv ../*nc restarts/
  cp -r ../template/* .
  timestamp=`ls restarts/restart.*.nc | grep -o [0-9].*[0-9]`
  mv restarts/particles.nc 'restarts/lagrangianParticleTrackingRestart.'${timestamp}'.nc'
  echo ' '`echo $timestamp | sed 's/\./:/g'` > Restart_timestamp
  msub *_realization_low_pass.sh
  cd ..
}

for i in {24..33}; do
  for j in {1..9}; do
    cp ../restarts/restart.00$i-0$j-01_00.00.00.nc .
    cp ../restarts/timeFiltersRestart.00$i-0$j-01_00.nc .
    setup_run $i $j
  done
  for j in {10..12}; do
    cp ../restarts/restart.00$i-$j-01_00.00.00.nc .
    cp ../restarts/timeFiltersRestart.00$i-$j-01_00.nc .
    setup_run $i $j
  done
done
