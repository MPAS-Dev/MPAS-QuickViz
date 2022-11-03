#!/usr/bin/env bash
# define execution directory
#MSUB -d ./
#MSUB -j oe
##MSUB -A s11_climateacme
#MSUB -A w14_mpaseddy
#MSUB -l nodes=5:ppn=24
#MSUB -l walltime=08:00:00
##MSUB -l depend=232506

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

i=0
for host in `cat hostnames.txt`; do
  if [ $i -lt $SLURM_NNODES ]; then

    #folder=`ls -al -d lagrPartTrack_folder/lagrP*nc | awk '{print $9}' | awk "NR==$((i+1))"`
    #folder=`ls -l -d ../low_pass_realization_*/ | awk '{print $9}' | awk "NR==$((i+1))"`
    folder=`ls -l -d realization_24-01_*/ | awk '{print $9}' | awk "NR==$((i+1))"`
   
    # for effective diffusivity
    #folder=`ls -l -d realization_*-07/ | awk '{print $9}' | awk "NR==$((i+1))"`

    echo 'Launching job on node '$host' for iterate '$i' output in file output'$i'-'$host'-'$jobid'.txt for '\
      $folder
    mpirun -np 1 -npernode 1 -host $host \
      time python \
      aggregate_particle_lonlat_rads.py -f $folder'/analysis_members/lagrPartTrack*nc' -p T -n 365 -r "np.arange(0,1)" -s 0 -o $folder'/365days_' \
      &> output"$i"-"$host"-"$jobid".txt &
      
      # for 50 day calc
      #aggregate_particle_lonlat_rads.py -f "$folder/analysis_members/lagrPartTrack*nc" -p T -n 57 --subset "np.arange(40,57)" -r "np.arange(0,1)" -s $i  \
      
      # for effective diffusivity
      #aggregate_particle_lonlat_rads.py -f $folder'/analysis_members/lagrPartTrack*nc' -p T -n 365 -r "np.arange(0,1)" -s 0 -o $folder \
     
      # for diffusivity for simple set of locations
      #aggregate_particle_lonlat_rads.py -f "$folder" -p T -n 28 -r "np.arange(0,1)" -s $i  \

    let i++
  fi
done
wait

rm hostnames.txt
mkdir -p all_rlzn_particle_data_rads_folder_50day
mv all_*nc all_rlzn_particle_data_rads_folder_50day

echo 'Finished on '`date`
