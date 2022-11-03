#!/usr/bin/env bash
# define execution directory
#MSUB -d ./
#MSUB -j oe
##MSUB -A s11_climateacme
#MSUB -A w14_mpaseddy
#MSUB -l nodes=12:ppn=24
#MSUB -l walltime=16:00:00

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

Njobs=10
i=0
for host in `cat hostnames.txt`; do
  low=$((Njobs*i))
  high=$((Njobs*(i+1)))

  echo 'Launching job on node '$host' for iterate '$i' output in file output'$i'-'$host'-'$jobid'.txt from '\
    $low' to '$high' realizations'
  mpirun -np 1 -npernode 1 -host $host \
    time python aggregate_particle_lonlat_rads.py -f 'lagrPartTrack.*nc' -p T -n 30 -r "np.arange($low,$high)"  \
    &> output"$i"-"$host"-"$jobid".txt &

  let i++
done
wait

rm hostnames.txt

echo 'Finished on '`date`
