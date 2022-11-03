#! /bin/bash
# Name of job script:
#MSUB -N LPT_pairs
# request resource: walltime
#MSUB -l walltime=3:00:00
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
module load python/2.7-anaconda-2.1.0
echo 'done'

echo ${PWD}
echo `date`

test=$(basename ${PWD}); 
resnum=${test%km*}
filtnum=${test#*km}
filtnum=${filtnum#_}
case=${resnum}'\,km\,'${filtnum}
echo ${case}

SCRIPTNAME=diffusivity_pdf.py
SCRIPTNAME=diffusivity_pdf_simple.py

DATAFILE=coarsest.npz
DATAFILE=coarse_sampling.npz

time python $SCRIPTNAME -f all_rlzn_particle_data.npz -s $DATAFILE --ts 1 -n 16  --dt 2. -d T -r 100000. -l 'np.array([0])' --frac 1.0 --hullpairs True &> kappa0.log &
time python $SCRIPTNAME -f all_rlzn_particle_data.npz -s $DATAFILE --ts 1 -n 16  --dt 2. -d T -r 100000. -l 'np.array([1])' --frac 1.0 --hullpairs True &> kappa1.log &
time python $SCRIPTNAME -f all_rlzn_particle_data.npz -s $DATAFILE --ts 1 -n 16  --dt 2. -d T -r 100000. -l 'np.array([2])' --frac 1.0 --hullpairs True &> kappa2.log &
time python $SCRIPTNAME -f all_rlzn_particle_data.npz -s $DATAFILE --ts 1 -n 16  --dt 2. -d T -r 100000. -l 'np.array([3])' --frac 1.0 --hullpairs True &> kappa3.log &
time python $SCRIPTNAME -f all_rlzn_particle_data.npz -s $DATAFILE --ts 1 -n 16  --dt 2. -d T -r 100000. -l 'np.array([4])' --frac 1.0 --hullpairs True &> kappa4.log &
wait
#time python convert_pdfs_vtk.py -d CaseClusterPDFsComplete_r\=100000.00_ts\=0_deltat\=2.00_Nlen\=16_Nc\=550_frac\=1.00/
#time python convert_pdfs_vtk.py -d CaseClusterPDFsSimple_r\=100000.00_ts\=0_deltat\=2.00_Nlen\=16_Nc\=550_frac\=1.00/
time python convert_pdfs_vtk.py -d CaseClusterRingPDFsSimple_r\=100000.00_ts\=0_deltat\=2.00_Nlen\=16_Nc\=550_frac\=1.00/

echo `date`
