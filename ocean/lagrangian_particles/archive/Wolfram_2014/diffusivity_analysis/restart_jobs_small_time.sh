#!/bin/tcsh 

foreach restartjob (`mshow | grep pwolfram | grep Idle | grep "13:00:00" | awk '{print $1}'`)
  set restartdir = `checkjob $restartjob | grep SubmitDir | awk '{print $2}'`
  
  echo 'canceling job in '$restartdir
  canceljob $restartjob

  cd $restartdir
 
  echo 'fixing submit script'
  sed -i 's/walltime=[0-9].*[0-9]/walltime=16:00:00/' 384_realization.sh

  echo 'submiting job in '$restartdir
  msub 384_realization.sh
 
  cd ../
end

