#! /bin/bash
# Name of job script:
#MSUB -N LPT_post
# request resource: walltime
#MSUB -l walltime=16:00:00
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

echo 'loading modules'
module purge
module load anconda
echo 'done'

echo ${PWD}
echo `date`

time python diffusivity_binning_tensor.py -f  all_rlzn_particle_data.npz -r  100000 --ts 1 --te 16 --dt 2. -d T -l 0 &> kappa_tensor0.log &
time python diffusivity_binning_tensor.py -f  all_rlzn_particle_data.npz -r  100000 --ts 1 --te 16 --dt 2. -d T -l 1 &> kappa_tensor1.log &
time python diffusivity_binning_tensor.py -f  all_rlzn_particle_data.npz -r  100000 --ts 1 --te 16 --dt 2. -d T -l 2 &> kappa_tensor2.log &
time python diffusivity_binning_tensor.py -f  all_rlzn_particle_data.npz -r  100000 --ts 1 --te 16 --dt 2. -d T -l 3 &> kappa_tensor3.log &
time python diffusivity_binning_tensor.py -f  all_rlzn_particle_data.npz -r  100000 --ts 1 --te 16 --dt 2. -d T -l 4 &> kappa_tensor4.log &
wait

echo `date`
