#!/usr/bin/env bash

# make the links between the files and here for further analysis

for i in `ls -d ../../low_pass_realization*/analysis_members/lagrPartTrack.*nc`; do 
  ln -sf $i .
done

