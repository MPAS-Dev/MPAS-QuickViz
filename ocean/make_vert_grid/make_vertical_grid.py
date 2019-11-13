def dz_z(z,Hmax,epsilon,dzmax):
    import numpy as np
    return dzmax*np.tanh(-z*np.pi/Hmax) + epsilon

from netCDF4 import Dataset
import numpy as np
"""
Write vertical grid to a netcdf file
"""
# the output array to write will be depth_t
dim_128_layer = 128
target_layers_arr = np.asarray([80])
# open a new netCDF file for writing.
ncfile = Dataset('MPAS-Ocean_vertical_grid.nc','w')
# create the depth_t dimension.
#ncfile.createDimension('dim_64_layer',target_layers_arr[0])
ncfile.createDimension('dim_80_layer',target_layers_arr[0])
#ncfile.createDimension('dim_128_layer',target_layers_arr[2])
# create the variable (4 byte integer in this case)
# first argument is name of variable, second is datatype, third is
# a tuple with the names of dimensions.
#refBottomDepth64 = ncfile.createVariable('refBottomDepth_64_layer',np.dtype('float64').char,('dim_64_layer'))
#refMidDepth64 = ncfile.createVariable('refMidDepth_64_layer',np.dtype('float64').char,('dim_64_layer'))
#refLayerThickness64 = ncfile.createVariable('refLayerThickness_64_layer',np.dtype('float64').char,('dim_64_layer'))

refBottomDepth80 = ncfile.createVariable('refBottomDepth_80_layer',np.dtype('float64').char,('dim_80_layer'))
refMidDepth80 = ncfile.createVariable('refMidDepth_80_layer',np.dtype('float64').char,('dim_80_layer'))
refLayerThickness80 = ncfile.createVariable('refLayerThickness_80_layer',np.dtype('float64').char,('dim_80_layer'))

#refBottomDepth128 = ncfile.createVariable('refBottomDepth_128_layer',np.dtype('float64').char,('dim_128_layer'))
#refMidDepth128 = ncfile.createVariable('refMidDepth_128_layer',np.dtype('float64').char,('dim_128_layer'))
#refLayerThickness128 = ncfile.createVariable('refLayerThickness_128_layer',np.dtype('float64').char,('dim_128_layer'))

epsilon = 1e-3
layer_min_thickness_arr = np.asarray([2.0])
layer_max_thickness_arr = np.asarray([223.0])
Hmax = 15000.0
Hmax2 = 6000.0
nLayers = 0

layerThickness = np.zeros((3,128))
dz = [epsilon]
z = [0]

for i,target_layers in enumerate(target_layers_arr):
    layer_min_thickness = layer_min_thickness_arr[i]
    layer_max_thickness = layer_max_thickness_arr[i]
    while nLayers != target_layers:
        zval = -epsilon
        dz = [epsilon]
        z = [0]
        nLayers = 0

        while zval > -Hmax2:
            difference = dz_z(zval,Hmax,epsilon,layer_max_thickness) - zval
            print difference,dz_z(zval,Hmax,epsilon,layer_max_thickness)
            while abs(difference) > 1e-2:
                zval -= epsilon
                difference = dz_z(zval,Hmax,epsilon,layer_max_thickness) -z[nLayers] + zval

            z.append(zval)
            dz.append(dz_z(zval,Hmax,epsilon,layer_max_thickness))
            nLayers += 1
            zval -= epsilon
            print 'nL = ',nLayers,layer_max_thickness

        dz_arr = np.asarray(dz)
        ind = abs(dz_arr - layer_min_thickness).argmin()

        dztemp = dz_arr[ind:]
        nLayers = len(dztemp)
        change = target_layers-nLayers
        layer_max_thickness -= float(change)
     
        print '*****len = ',len(dztemp)
    print np.sum(dztemp)
    layerThickness[i,:nLayers] = dztemp

# 64 layers:
#nVertLevels = 64
#botDepth = np.zeros(nVertLevels)
#midDepth = np.zeros(nVertLevels)
#botDepth[0] = layerThickness[0,0]
#midDepth[0] = 0.5*layerThickness[0,0]

#for i in range(1,nVertLevels):
#    botDepth[i] = botDepth[i-1] + layerThickness[0,i]
#    midDepth[i] = midDepth[i-1] + 0.5*(layerThickness[0,i] + layerThickness[0,i-1])
#
#refBottomDepth64[:] = botDepth
#refMidDepth64[:] = midDepth
#refLayerThickness64[:] = layerThickness[0,:nVertLevels]

# 80 layers:
nVertLevels = 80
layerThickness[1,:] = layerThickness[0,:]
botDepth = np.zeros(nVertLevels)
midDepth = np.zeros(nVertLevels)
botDepth[0] = layerThickness[1,0]
midDepth[0] = 0.5*layerThickness[1,0]

for i in range(1,nVertLevels):
    botDepth[i] = botDepth[i-1] + layerThickness[1,i]
    midDepth[i] = midDepth[i-1] + 0.5*(layerThickness[1,i] + layerThickness[1,i-1])

print botDepth
refBottomDepth80[:] = botDepth
refMidDepth80[:] = midDepth
refLayerThickness80[:] = layerThickness[1,:nVertLevels]

# 128 layers:
#nVertLevels = 128
#botDepth = np.zeros(nVertLevels)
#midDepth = np.zeros(nVertLevels)
#botDepth[0] = layerThickness[2,0]
#midDepth[0] = 0.5*layerThickness[2,0]

#for i in range(1,nVertLevels):
#    botDepth[i] = botDepth[i-1] + layerThickness[2,i]
#    midDepth[i] = midDepth[i-1] + 0.5*(layerThickness[2,i] + layerThickness[2,i-1])

#refBottomDepth128[:] = botDepth
#refMidDepth128[:] = midDepth
#refLayerThickness128[:] = layerThickness[2,:nVertLevels]
# close the file.
ncfile.close()
print '*** SUCCESS writing the vertical grid file!'
