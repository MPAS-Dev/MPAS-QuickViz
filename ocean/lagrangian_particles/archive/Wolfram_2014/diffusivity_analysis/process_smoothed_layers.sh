#!/usr/bin/env bash

#split into layers
echo 'spliting layers'
../../split_layers.sh $1_cluster_data.txt $2

for alayer in $(seq 0 $2); do
  layer=`printf "%.6d" $alayer`
  echo 'processing ' $layer
  # Process with parallel kernel
  ../../mpi_process_diffusivity.sh $2 $1_cluster_data.txt_layer$layer $1_interp_points.txt_layer$layer $1_output.txt_layer$layer $1_layer$layer.app 1.0
  # Plot with
  echo 'plotting ' $layer
  ../../plot_diffusivity_layer.py -b $1 -l $layer
done
