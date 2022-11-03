#!/usr/bin/env bash

# used to provide input files to compute stratification, RRD, and cw

START=24
END=33

for i in `seq -s' ' ${START} ${END}`; do 
  echo ncra timeSeriesStats.00$i*.nc tavg.00$i.nc &
  ncra timeSeriesStats.00$i*.nc tavg.00$i.nc &
done
wait

echo ncra tavg.00[0-9][0-9].nc tavg.00${START}-00${END}.nc
ncra tavg.00[0-9][0-9].nc tavg.00$START-00$END.nc
