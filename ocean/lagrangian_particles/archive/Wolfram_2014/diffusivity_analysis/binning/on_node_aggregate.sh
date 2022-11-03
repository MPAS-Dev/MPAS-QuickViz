#!/usr/bin/env bash
# define execution directory
#MSUB -d ./
#MSUB -j oe
##MSUB -A s11_climateacme
#MSUB -A w14_mpaseddy
#MSUB -l nodes=120:ppn=16
#MSUB -l walltime=00:15:00

echo 'Started on '`date`

#NFOLDERS=`ls -al -d realization_*/ | wc | awk '{print $1}'`
NFOLDERS=16

batch_of_runs() {
  STARTFOLDER=$1
  echo 'Starting realization number '$STARTFOLDER
  i=0
  while [ $i -lt $NFOLDERS ]; do

    rlznnum=$((i + STARTFOLDER))
    folder=`ls -al -d realization_*/ | awk '{print $9}' | awk "NR==$((rlznnum + 1))"`

    echo 'Launching job iterate '$rlznnum' output in file output'$rlznnum'-on_node.txt'
    time python aggregate_particle_lonlat_rads.py -f "$folder/analysis_members/lagrPartTrack.*nc" -p T -n 180 -r "np.arange(0,1)" -s $rlznnum  \
      &> output"$rlznnum"-on_node.txt &

    let i++
  done
  wait
}

for i in {0..7}; do
  batch_of_runs $((i*16))
done

echo 'Finished on '`date`
