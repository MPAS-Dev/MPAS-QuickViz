#!/bin/tcsh
set usefulprocs = 6
set numEnsembles = 30
set count = 0
set totcount = 0
@ numEnsembles--
foreach i (analyze_*)
  cd $i
  pwd
  @ count++
  @ totcount++
  if ( $totcount > $numEnsembles) then 
    time python ../diffusivity_realization_data.py -f  buoyancySurface_particleData.nc -c ../clusterDefinition >& compute_realization_stats.log 
    sleep 60
    break
  endif
  if ( $count >= $usefulprocs ) then
    # wait on the (last) script to finish
    time python ../diffusivity_realization_data.py -f  buoyancySurface_particleData.nc -c ../clusterDefinition >& compute_realization_stats.log 
    sleep 60
    set count = 0
  else 
    # script should run in backgroud
    time python ../diffusivity_realization_data.py -f  buoyancySurface_particleData.nc -c ../clusterDefinition >& compute_realization_stats.log  &
  endif
  cd .. 
end

## then need to run (in serial for now)
#./aggregate_diffusivity.py -r '.' -p 'analyze_restart'
#
## now output the ensemble data
#./write_final_ensemble.py -e 30
#
## need to split up the cluster_data into multiple files
#
### take output and compute diffusivity
##module load intel
##make compile_intel
### copy in interpolation file interp_points_16km.txt
##cp ../repeated_buoyancy_16km/interp_points.txt  interp_points_16km.txt
### run the diffusivity calculation
##scaleFactor=1 make run
