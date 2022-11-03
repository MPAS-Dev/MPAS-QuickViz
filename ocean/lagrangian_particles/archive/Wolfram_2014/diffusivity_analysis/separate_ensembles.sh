#!/bin/bash

cd ensembles/


for i in diffusivityOutput*; do
  cd $i

  # separate out the ensembles
  for i in {0..4}; do 
    echo 480816 > buoyancySurfaceCluster$i.txt; gsed -n $((480816*i+2))',+480815p' cluster_data.txt >> buoyancySurfaceCluster$i.txt ; 
  done

  # separate out the interp points too.
  for i in {0..4}; do 
    echo 30051 > interp$i.txt; gsed -n $((30051*i+2))',+30050p' interp_points.txt >> interp$i.txt ; 
  done

  cd ../

done
cd ../
