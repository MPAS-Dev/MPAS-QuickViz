#!/bin/tcsh 
module purge 
module load pgi

# get mean values for the buoyancy surface
/usr/projects/cesm/software/conejo/nco/bin/ncra -v buoyancySurfaceVelocityZonal,buoyancySurfaceVelocityMeridional,buoyancySurfaceDepth analyze_*/output*.nc buoyancySurface.nc 
/usr/projects/cesm/software/conejo/nco/bin/ncks -v lonCell,latCell,buoyancySurfaceValues analyze_restart.0005-01-01_00.00.00.nc/output.*.nc buoyancySurface.nc

# compute the buoyancy surface files
foreach i (analyze_*)
  cd $i
  /usr/projects/cesm/software/conejo/nco/bin/ncks -v xParticle,yParticle,zParticle,buoyancyParticle,zLevelParticle,lonVel,latVel,sumU,sumV,sumUU,sumUV,sumVV,buoyancySurfaceValues output.*.nc buoyancySurface_particleData.nc
  cd ..
end
