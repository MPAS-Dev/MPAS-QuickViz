#!/usr/bin/env bash
# define execution directory
#MSUB -d ./
#MSUB -j oe
##MSUB -A s11_climateacme
#MSUB -A w14_mpaseddy
#MSUB -l nodes=01:ppn=24
#MSUB -l walltime=16:00:00

echo 'Started on '`date`

echo time python reduce_dispersion_file.py
time python reduce_dispersion_file.py

echo 'Finished on '`date`
