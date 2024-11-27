# E3SM_FOV

René M. van Westen

These directories contain Python (v3) scripts for plotting/analysing various E3SM model output.

Python scripts can be found in the directory 'Program'.
Model output can be found in the directory 'Data'.

The processed model output are stored as NETCDF files and using the relevant scripts one can regenerate all the figures.
I provided a selection of the model output (interpolated onto the 0.5x0.5 rectangular grid) and is only converted to yearly-averaged data (due to storage limitations). 
Some scripts (e.g., FOV_index.py and AMOC_transport.py) use the yearly-averaged model output, but you can use also the already available time series. 

The Field_generation_*.py scripts can only be used on perlmutter.nersc.gov machine, where you need the following conda environment: 
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
