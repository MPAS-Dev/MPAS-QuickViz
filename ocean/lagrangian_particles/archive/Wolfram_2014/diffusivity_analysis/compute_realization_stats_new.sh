#!/bin/tcsh
# Name of job script:
#MSUB -N compute_diff
# request resource: walltime
#MSUB -l walltime=02:00:00
#  request resource: number of nodes
#MSUB -l nodes=1:ppn=24
# define execution directory
#MSUB -d ./
# error and output in same file
#MSUB -j oe
#MSUB -m abe pwolfram@lanl.gov 
# define account associated with the job
# general climate
##MSUB -A s11_climate
# mustang job
#MSUB -A w14_mpaseddy
##MSUB -l depend=970675

module load python-epd
setenv PYTHONPATH /users/pwolfram/lib/python2.7/site-packages/

set totcount = 0
foreach i (analyze_*)
  if  ( ! (  -e $i/log.0000.out )) then
    @ totcount++
  endif
end 

set totcount = 12

set count = 0
foreach i (analyze_*)
  cd $i
  @ count++
  if  ( ! (  -e $i/log.0000.out )) then
    if ( $count >= $totcount ) then
      # wait on the (last) script to finish
      time ../diffusivity_realization_data_output.py -f output* -c ../clusterDefinition -s 15 -e 30 -i 2 >& compute_realization_stats.log
      sleep 60
      set count = 0 
    else 
      # script should run in backgroud
      time ../diffusivity_realization_data_output.py -f output* -c ../clusterDefinition -s 15 -e 30 -i 2 >& compute_realization_stats.log &
    endif
  endif
  cd .. 
end

## split into layers:
#/bin/bash -c 'for i in {0..4}; do echo 113056 > buoyancySurfaceCluster$i.txt; sed -n $((133056*i+2))',+133056p' cluster_data.txt >> buoyancySurfaceCluster$i.txt ; done'
#
## take output and compute diffusivity
#module load intel
#make compile_intel
## copy in interpolation file interp_points_16km.txt
#cp ../repeated_buoyancy_16km/interp_points.txt  interp_points_16km.txt
## run the diffusivity calculations
#/bin/bash -c 'module load intel; make compile_intel; for i in {0..4}; do time ./process_diffusivity 1 buoyancySurfaceCluster$i.txt interp_points_16km.txt kappaInterp_buoySurf$i.txt; done'
##time ./process_diffusivity 1 cluster_data.txt interp_points_16km.txt kappaInterp_result_openmp_sf1.txt
#
### build plots
#/bin/bash -c 'module load anconda; for i in {0..4}; do time ./plot_diffusivity_layers.py -l $0 ; done'
#
###!/bin/bash
## bash, but can have trouble because module doesn't get loaded
##source ~/module_source.sh
##module load python-epd
##export PYTHONPATH=/users/pwolfram/lib/python2.7/site-packages/
##for i in analyze_*; do
##  cd $i
##  ../compute_diffusivity.py -f output* &>compute_realization_stats.log &
##  cd .. 
##done
