#! /bin/bash
# Name of job script:
#MSUB -N LPT_post
# request resource: walltime
#MSUB -l walltime=16:00:00
# define execution directory
#MSUB -d ./
# error and output in same file
#MSUB -j oe
#MSUB -m abe pwolfram@lanl.gov
# define account associated with the job
# general climate
##MSUB -A s11_climateacme
##MSUB -l nodes=1:ppn=16
## mustang job
#MSUB -A w14_mpaseddy
##  request resource: number of nodes
#MSUB -l nodes=1:ppn=24

echo 'loading modules'
module load python/2.7-anaconda-2.1.0
source activate LIGHT_analysis
echo 'done'

echo ${PWD}
echo `date`

# 30 GB mem buffer to prevent memory error issues
membuff=$((30*2**30))
mem=`cat /proc/meminfo | grep MemFree | awk '{print $2}'`

time python aggregate_particle_lonlat_rads.py -f 'lagrPartTrack.*nc' -p T -n 30

# compute diffusivity
for i in {0..1}; do
  time python diffusivity_binning_complete.py -f all_rlzn_particle_data_rads.npz -s mesh.nc --ts 1 -n 30  --dt 1. -d F -r 100000. -l "np.array([$i])" --frac 1.0 -p T &
  sleep 30s
  while [ $mem -gt $membuff ]; do
    sleep 30s
    mem=`cat /proc/meminfo | grep MemFree | awk '{print $2}'`
  done
done
wait

for i in {2..3}; do
  time python diffusivity_binning_complete.py -f all_rlzn_particle_data_rads.npz -s mesh.nc --ts 1 -n 30  --dt 1. -d F -r 100000. -l "np.array([$i])" --frac 1.0 -p T &
  sleep 30s
  while [ $mem -gt $membuff ]; do
    sleep 30s
    mem=`cat /proc/meminfo | grep MemFree | awk '{print $2}'`
  done
done
wait

# compute diffusivity
for i in {4..5}; do
  time python diffusivity_binning_complete.py -f all_rlzn_particle_data_rads.npz -s mesh.nc --ts 1 -n 30  --dt 1. -d F -r 100000. -l "np.array([$i])" --frac 1.0 -p T &
  sleep 30s
  while [ $mem -gt $membuff ]; do
    sleep 30s
    mem=`cat /proc/meminfo | grep MemFree | awk '{print $2}'`
  done
done
wait

# compute diffusivity
for i in {6..7}; do
  time python diffusivity_binning_complete.py -f all_rlzn_particle_data_rads.npz -s mesh.nc --ts 1 -n 30  --dt 1. -d F -r 100000. -l "np.array([$i])" --frac 1.0 -p T &
  sleep 30s
  while [ $mem -gt $membuff ]; do
    sleep 30s
    mem=`cat /proc/meminfo | grep MemFree | awk '{print $2}'`
  done
done
wait

for i in {8..9}; do
  time python diffusivity_binning_complete.py -f all_rlzn_particle_data_rads.npz -s mesh.nc --ts 1 -n 30  --dt 1. -d F -r 100000. -l "np.array([$i])" --frac 1.0 -p T &
  sleep 30s
  while [ $mem -gt $membuff ]; do
    sleep 30s
    mem=`cat /proc/meminfo | grep MemFree | awk '{print $2}'`
  done
done
wait

# compute diffusivity
for i in {10..10}; do
  time python diffusivity_binning_complete.py -f all_rlzn_particle_data_rads.npz -s mesh.nc --ts 1 -n 30  --dt 1. -d F -r 100000. -l "np.array([$i])" --frac 1.0 -p T &
  sleep 30s
  while [ $mem -gt $membuff ]; do
    sleep 30s
    mem=`cat /proc/meminfo | grep MemFree | awk '{print $2}'`
  done
done
wait


cd CaseBinningEnsembledComplete_r\=100000.00_ts\=0_deltat\=1.00_Nlen\=30_Nc\=22800_frac\=1.00/
time python ../meanDepths.py
time python ../see_krr.py
cd ../

time python diffusivity_binning_complete_plot.py -d CaseBinningEnsembledComplete_r\=100000.00_ts\=0_deltat\=1.00_Nlen\=30_Nc\=22800_frac\=1.00/ -r '10\,km'

echo `date`
