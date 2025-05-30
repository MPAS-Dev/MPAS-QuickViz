#!/bin/tcsh 

foreach restartjob (`grep CAN analyze_restart.000*/slurm* | xargs -I {} dirname {}`)
  cd $restartjob
  echo 'restarting job in '$restartjob
  rm log.000*
  rm output*
  rm slurm*
  rm stats_*
  mv restart*_backup.nc restart*0.nc
  
  echo 'correcting time for submit script'
  sed -i 's/walltime=[0-9].*[0-9]/walltime=16:00:00/' 384_realization.sh
  
  echo 'submiting job in '$restartjob
  msub 384_realization.sh
  
  cd ../
end

