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

module load anconda
setenv PYTHONPATH /users/pwolfram/lib/python2.7/site-packages/

set totcount = 0
foreach i (analyze_*)
  if  ( ! (  -e $i/log.0000.out )) then
    @ totcount++
  endif
end 

set totcount = 24

set count = 0
foreach i (analyze_*)
  cd $i
  if  ( ! (  -e log.0000.out )) then
  #if  ( ! (  -e lagrangian_data.npz )) then
    @ count++
    echo $i
    if ( $count >= $totcount ) then
      # wait on the (last) script to finish
      echo time ../diffusivity_realization_data_lagrangian.py -f output* -c ../clusterDefinition >& compute_realization_stats.log
      time ../diffusivity_realization_data_lagrangian.py -f output* -c ../clusterDefinition >& compute_realization_stats.log
      sleep 60
      set count = 0 
    else 
      # script should run in backgroud
      echo time ../diffusivity_realization_data_lagrangian.py -f output* -c ../clusterDefinition >& compute_realization_stats.log 
      time ../diffusivity_realization_data_lagrangian.py -f output* -c ../clusterDefinition >& compute_realization_stats.log &
    endif
  endif
  cd .. 
end
wait
