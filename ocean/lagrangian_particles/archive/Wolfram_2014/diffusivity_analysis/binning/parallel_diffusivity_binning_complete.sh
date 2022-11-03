#!/usr/bin/env bash
# define execution directory
#MSUB -d ./
#MSUB -j oe
##MSUB -A s11_climateacme
#MSUB -A w14_mpaseddy
#MSUB -l nodes=120:ppn=24
#MSUB -l walltime=03:00:00
#MSUB -l depend=244041

echo 'Started on '`date`

module load gcc openmpi

jobid=$SLURM_JOB_ID

# get files with environments for each node
rm -rf sadc_outfile
mpirun -npernode 1 -mca btl tcp,self --output-filename sadc_outfile /bin/bash -c "/bin/env "

# filter nodes to get hostnames
grep SLURMD_NODENAME sadc_outfile.1.* | sed 's/.*=//g' | sort > hostnames.txt
rm -rf sadc_outfile*

# this is the list of hostnames
#cat hostnames.txt

CPUPERCORE=`cat /proc/cpuinfo  | grep processor | wc | awk '{print $1}'`

i=0
for host in `cat hostnames.txt`; do 

  rlznnum=`printf "%d" $i`
  echo 'Launching job on node '$host' for iterate '$i' with realization='$rlznnum' output in file output'$i'-'$host'-'$jobid'.txt'
  #{{{
  # mpirun -np 1 -npernode 1 -host $host time python diffusivity_binning_complete.py -f all_rlzn_particle_data_rads.npz -s mesh.nc --ts 1 -n 30  --dt 1. -d F -r 100000. -l "np.array([$i])" --frac 1.0 -p T > output"$i"-"$host"-"$jobid".txt &
  #mpirun -np 1 -npernode 1 -host $host time python diffusivity_binning_complete.py -f 'all_rlzn_particle_data_rads_rlzn'${rlznnum}'.nc' -s mesh.nc  -n 30   -d F -r 100000.  --frac 1.0 -p T > output"$i"-"$host"-"$jobid".txt &
  #}}}
  # run each layer simultaneously on the node, but allow each script to have the full cores of each node
  # change run_realization_diffusivity.sh to accomodate serial (particle reset) or parallel realization modes via -n parameter
  mpirun -np 1 -npernode 1 -bind-to-core -cpus-per-proc ${CPUPERCORE} -host $host \
    time ./run_realization_diffusivity.sh ${rlznnum} 0 10 28\
    &> output"$i"-"$host"-"$jobid".txt &
  let i++
done
wait

rm hostnames.txt

echo 'Finished on '`date`
