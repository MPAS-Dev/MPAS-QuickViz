#!/bin/tcsh
./make_restart_dirs.sh
foreach i (analyze_*)
  cd $i
  cp ../template/* .
  msub 384_realization.sh
  cd .. 
end
