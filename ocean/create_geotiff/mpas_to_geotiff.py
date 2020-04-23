import yaml
import os
import pprint
import numpy as np
import gdal, osr
from netCDF4 import Dataset
from matplotlib.tri import Triangulation, LinearTriInterpolator

#######################################################################
#######################################################################

def write_to_geotiff(val,tri,nx,ny,bbox,fout):

  # Generate raster points to interpolate from mesh  
  xinterp = np.linspace(bbox[0],bbox[1],nx)
  yinterp = np.linspace(bbox[2],bbox[3],ny)
  Xinterp,Yinterp = np.meshgrid(xinterp,yinterp)

  # Calculate transformation factors 
  originX = xinterp[0]
  originY = yinterp[-1]
  xres = (xinterp[-1] - xinterp[0]) / float(nx)
  yres = (yinterp[-1] - yinterp[0]) / float(ny)

  # Intrepolate solution onto raster points
  interp = LinearTriInterpolator(tri,val)
  interp_vals = interp(Xinterp.ravel(),Yinterp.ravel()).data
  raster = np.reshape(interp_vals,Xinterp.shape)
  raster = raster[::-1]

  # Create and write geotiff
  driver = gdal.GetDriverByName('GTiff')
  outRaster = driver.Create(fout,nx,ny, 1, gdal.GDT_Float32)
  outRaster.SetGeoTransform((originX, xres, 0, originY, 0,-yres))
  outband = outRaster.GetRasterBand(1)
  outband.WriteArray(raster)
  outRasterSRS = osr.SpatialReference()
  outRasterSRS.ImportFromEPSG(4326)
  outRaster.SetProjection(outRasterSRS.ExportToWkt())
  outband.FlushCache()

#######################################################################
#######################################################################

if __name__ == "__main__":

  # Read in config file
  pwd = os.getcwd()
  f = open(pwd+'/mpas_to_geotiff.config')
  cfg = yaml.load(f,Loader=yaml.Loader)
  pprint.pprint(cfg)

  # Read data from output file
  ncfile = Dataset(cfg['output_file'],'r')
  var_dim = ncfile.variables[cfg['output_variable']].dimensions
  if var_dim == ('nCells',):
    var = ncfile.variables[cfg['output_variable']][:]
  elif var_dim == ('Time','nCells'): 
    var = ncfile.variables[cfg['output_variable']][cfg['time_index'],:]
  else:
    print('Incompatible output variable')
    raise SystemExit(0)

  if cfg['inundation']:
    bathy = ncfile.variables['bottomDepth'][:]
    var = var + bathy
    var[bathy > 0] = np.nan
    var[var < 0.01] = np.nan

  # Read data from mesh file
  ncmesh = Dataset(cfg['mesh_file'],'r')
  lon_mesh = np.rad2deg(np.mod(ncmesh.variables['lonCell'][:] + np.pi, 2.0*np.pi) - np.pi)
  lat_mesh = np.rad2deg(ncmesh.variables['latCell'][:])
  nEdgesOnCell = ncmesh.variables['nEdgesOnCell'][:]
  cellsOnCell = ncmesh.variables['cellsOnCell'][:,:]
 
  # Triangulate cells 
  triangles = Triangulation(lon_mesh,lat_mesh)

  # Compute triangulation mask (needs to be vectorized)
  mask = np.array(np.zeros((triangles.triangles.shape[0],)),dtype=bool)
  ntri = triangles.neighbors.shape[0]
  for i in range(ntri):
    k = 0
    for j in range(3):
      n = triangles.triangles[i,j]
      if nEdgesOnCell[n] != np.where(cellsOnCell[n,:] != 0)[0].size:  # Mask triangles
        k = k +1                                                      # containing
    if k == 3:                                                        # 3 boundary
      mask[i] = True                                                  # cells
  triangles.set_mask(mask)

  # Write out geotiff image
  output_name = cfg['output_variable']+'.tif'
  write_to_geotiff(var,triangles,cfg['nx'],cfg['ny'],cfg['bbox'],output_name)
