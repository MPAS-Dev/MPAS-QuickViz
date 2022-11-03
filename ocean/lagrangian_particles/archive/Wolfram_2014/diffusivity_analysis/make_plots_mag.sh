#!/bin/bash

for i in {0..4}; do  ../../plot_diffusivity_layers.py -l $i -c 1e5; done
