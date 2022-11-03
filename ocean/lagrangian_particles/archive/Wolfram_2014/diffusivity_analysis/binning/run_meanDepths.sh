#!/usr/bin/env bash
# define execution directory
#MSUB -d ./
#MSUB -j oe
##MSUB -A s11_climateacme
#MSUB -A w14_mpaseddy
#MSUB -l nodes=01:ppn=24
#MSUB -l walltime=16:00:00

echo 'Started on '`date`

time python meanDepths.py

echo 'Finished on '`date`
