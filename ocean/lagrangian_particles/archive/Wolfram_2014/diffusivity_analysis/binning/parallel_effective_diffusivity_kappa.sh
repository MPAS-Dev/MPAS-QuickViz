#!/usr/bin/env bash
#MSUB -j oe
#MSUB -A w14_mpaseddy
#MSUB -l nodes=11:ppn=24
#MSUB -l walltime=16:00:00
##MSUB -l depend=278928

hash=`date +%y%m%d%H%M%S%N`
casedir='/lustre/scratch3/turquoise/pwolfram/ZISO_5km/realizations/realization_24-01_15dayreset/'
scriptname='EffDiffAnalysis_5km_production_kappa_15dayreset.py'
kappa=250
currentdir=`pwd`
echo 'Started on '`date`' for '$casedir' starting in '$currentdir

module load gcc openmpi

cd $casedir

jobid=$SLURM_JOB_ID

# get files with environments for each node
rm -rf $hash'sadc_outfile'
mpirun -npernode 1 -mca btl tcp,self --output-filename $hash'sadc_outfile' /bin/bash -c "/bin/env "

# filter nodes to get hostnames
grep SLURMD_NODENAME $hash'sadc_outfile.1.'* | sed 's/.*=//g' | sort > $hash'hostnames.txt'
rm -rf $hash'sadc_outfile'*

# this is the list of hostnames
cat $hash'hostnames.txt'

nprocs=`wc hostnames.txt | awk '{print $1}'`

cd $currentdir

i=0
for host in `cat $casedir$hash'hostnames.txt'`; do 

  echo 'Launching job on node '$host' for iterate '$i' out of '$nprocs' with output in file '$casedir'/'$outputdir'output'$i'-'$host'-'$jobid'.txt'
  mpirun -np 1 -npernode 1 -host $host time \
    python $scriptname  $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' $kappa > $casedir/output"$i"-"$host"-"$jobid".txt &
    #python EffDiffAnalysis_5km_production_kappa_5dayreset.py  $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' $kappa > $casedir/output"$i"-"$host"-"$jobid".txt &
    #python EffDiffAnalysis_5km_kappa500.py $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' > $casedir/output"$i"-"$host"-"$jobid".txt &
    #python EffDiffAnalysis_5km_kappa25.py $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' > $casedir/output"$i"-"$host"-"$jobid".txt &
    #python EffDiffAnalysis_5km_kappa250.py $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' > $casedir/output"$i"-"$host"-"$jobid".txt &
    #python EffDiffAnalysis_5km_kappa.py $i $casedir '/lustre/scratch3/turquoise/pwolfram/ZISO_5km/mesh.nc' > $casedir/output"$i"-"$host"-"$jobid".txt &
  
  let i++
done
wait


echo 'Finished on '`date`
