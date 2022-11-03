#!/usr/bin/env bash
#MSUB -j oe
#MSUB -A w14_mpaseddy
#MSUB -l nodes=11:ppn=24
#MSUB -l walltime=16:00:00
#MSUB -l depend=234584

casedir='/lustre/scratch3/turquoise/pwolfram/ZISO_5km/realizations/realization_33-07/'
currentdir=`pwd`
echo 'Started on '`date`' for '$casedir' starting in '$currentdir

module load gcc openmpi

cd $casedir

jobid=$SLURM_JOB_ID

# get files with environments for each node
rm -rf sadc_outfile
mpirun -npernode 1 -mca btl tcp,self --output-filename sadc_outfile /bin/bash -c "/bin/env "

# filter nodes to get hostnames
grep SLURMD_NODENAME sadc_outfile.1.* | sed 's/.*=//g' | sort > hostnames.txt
rm -rf sadc_outfile*

# this is the list of hostnames
cat hostnames.txt

nprocs=`wc hostnames.txt | awk '{print $1}'`

cd $currentdir

i=0
for host in `cat $casedir'hostnames.txt'`; do 

  echo 'Launching job on node '$host' for iterate '$i' out of '$nprocs' with output in file '$casedir'/'$outputdir'output'$i'-'$host'-'$jobid'.txt'
  mpirun -np 1 -npernode 1 -host $host time \
    python EffDiffAnalysis_5km_kappa.py $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' > $casedir/output"$i"-"$host"-"$jobid".txt &
    #python EffDiffAnalysis_5km_long.py $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' > $casedir/output"$i"-"$host"-"$jobid".txt &
    #python EffDiffAnalysis_5km_time.py $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' > $casedir/output"$i"-"$host"-"$jobid".txt &
    #python ZISO_driver.py -p $i -n $nprocs -r 22 > $outputdir/output"$i"-"$host"-"$jobid".txt &
  
  let i++
done
wait


echo 'Finished on '`date`
