#!/bin/bash

#for i in surface_vertex_*backup.nc ; do 
#  j=${i#surface_vertex_}
#  file=${j%_backup.nc}.nc 
for i in restart*.nc ; do 
  file=${i}
  mkdir analyze_$file
  cd analyze_$file
  ln -s ../$i $file
  cd ../
done
