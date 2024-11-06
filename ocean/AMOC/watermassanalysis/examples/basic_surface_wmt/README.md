# Basic surface water mass transformation for MPAS runs

The core surface water mass transformation (WMT) calculations are located in [`watermasstools.py`](https://github.com/MPAS-Dev/MPAS-QuickViz/blob/master/ocean/AMOC/watermassanalysis/modules/watermasstools.py). These calculations contain one additional dependency beyond the `e3sm_unified` environment

   - [`fastjmd95`](https://github.com/xgcm/fastjmd95) -- A numba-accelerated package for the Jackett and McDougall (1995) equation of state

### Serial usage

The command line executable module `basic_surface_wmt.py` is a postprocessing wrapper around the core WMT functions. This module is set up to accomplish two tasks

   1. Build a basic coordinates file `python basic_surface_wmt.py -c [MESHNAME]`
   2. Process a single monthly results file `python basic_surface_wmt.py -f [FILENUMBER] [MESHNAME]`

Here `MESHNAME` is either `LR` (`EC30to60E2r2`) or `HR` (`oRRS18to6v3`). Both options use the CORE-II E3SMv2 G-cases with the `20210421_sim7` tag. Different runs can be specified in `parameters.yaml` (small changes to `python basic_surface_wmt.py` will probably also be required).

### Parallel usage

The workflow is set up to process single monthly files such that each serial task can be distributed across multiple cpus on a single node. I use GNU Parallel for this. The general workflow is as follows

```
#!/bin/bash
#SBATCH --job-name=JOB_NAME
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --constraint=cpu
#SBATCH --exclusive
#SBATCH --output=JOB_NAME.o%j
#SBATCH --time=1:00:00
#SBATCH --account=ACCOUNT

source $HOME/.bashrc
module load parallel
conda activate ENV_NAME

meshName=LR
savePath="/path/to/${meshName}"
mkdir -p "${savePath}/monthly_files"
mkdir -p "${savePath}/concatenated"

# Calculate coords file (-c loads the coords file)
python basic_surface_wmt.py -p "${savePath}" -c "${meshName}"

# We use GNU Parallel to run the serial monthly processing across all cpus on our node
PARALLEL_OPTS="-N 1 --delay .2 -j $SLURM_NTASKS --joblog parallel-${SLURM_JOBID}.log"

# Here GNU Parallel distributes the individual SRUN tasks, so for SRUN
SRUN_OPTS="-N 1 -n 1 --exclusive"

# Process the monthly files, $(seq 0 119) does the first 10 years
parallel $PARALLEL_OPTS srun $SRUN_OPTS \
    python basic_surface_wmt.py -p "${savePath}" -f {} "${meshName}" ::: $(seq 0 119)

# Concatenate monthly files
ncrcat -h ${savePath}/monthly_files/${meshName}_WMT1D* \
    ${savePath}/concatenated/${meshName}_WMT1D_years1-10.nc
ncrcat -h ${savePath}/monthly_files/${meshName}_WMT2D* \
    ${savePath}/concatenated/${meshName}_WMT2D_years1-10.nc
```

### Visualization

A basic visualization example can be found in `basic_surface_wmt.ipynb`.