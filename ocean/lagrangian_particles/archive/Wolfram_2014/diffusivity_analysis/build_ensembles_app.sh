#!/bin/bash
# build the mpi app to processes ensemble files in parallel
name=build_ensembles.app
submitname=submit_build_ensemble.tcsh
rm $name
range=$(($1/2))
for i in $(eval echo {0..$range}); do
  echo "-np 1 python diffusivity_build_ensembles.py -r /panfs/scratch3/vol5/pwolfram/SOMA_runs/repeated_buoyancy_32km_vertical_redo/ -p analyze_restart -l 5 -s $i -e $(($1-$i))" >> $name
done

cat $name

# build submit scipt

cat <<EOT > $submitname
#! /bin/tcsh
# Name of job script:
#MSUB -N build_ensm
# request resource: walltime
#MSUB -l walltime=0:30:00
#  request resource: number of nodes
#MSUB -l nodes=$(($range+1)):ppn=1
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

# for mustang
module purge
module load intel openmpi anconda

mpirun --app build_ensembles.app
EOT

cat $submitname
