from netCDF4 import Dataset
import numpy as np

F = 5 # Number of output files

# --- Open and read vars from NC file 
for nf in range(F):
    print('Processing ifile: output_'+str(nf+1)+'.nc')

    ncfile_init = Dataset('initial_state.nc','r')
    ncfile_out = Dataset('output_'+str(nf+1)+'.nc','r')
    
    xCell        = ncfile_init.variables['xCell']
    yCell        = ncfile_init.variables['yCell']
    xEdge        = ncfile_init.variables['xEdge']
    yEdge        = ncfile_init.variables['yEdge']
    areaCell     = ncfile_init.variables['areaCell']
    maxLevelCell = ncfile_init.variables['maxLevelCell']
    bottomDepth  = ncfile_init.variables['bottomDepth']

    xtime                 = ncfile_out.variables['xtime']
    hFull                 = ncfile_out.variables['layerThickness']
    densityFull           = ncfile_out.variables['density']
    kineticEnergyCellFull = ncfile_out.variables['kineticEnergyCell']
    vertTransportFull     = ncfile_out.variables['vertTransportVelocityTop']
    
    vertTransportFull = np.absolute(vertTransportFull)
    
    # --- Get array size
    itmp = hFull.shape
    K = itmp[2]
    T = itmp[0]
    nCells = itmp[1]
    
    ridgeDepth = 500.0
    bottomMax = np.max(bottomDepth)
    yMax = np.max(yEdge)
    yMin = np.min(yEdge)
    xMax = np.max(xEdge)

    gravity = 9.80616
    time = np.zeros(T)
    rpe = np.zeros(T)
    keMeanVolume = np.zeros(T)
    vertTransportVolume = np.zeros(T)
    vertTransportVolumeZ = np.zeros([K,T])
    
    for nt in range(T):
        print('Processing nt:',nt,'/',T-1)
    
        t = xtime[nt,:]
    
        h = hFull[nt,:,:]
        density = densityFull[nt,:,:]

        vol_1D = np.zeros([K*nCells])
        density_1D = np.zeros([K*nCells])

        i = 0
        for iCell in range(nCells):
            for k in range(maxLevelCell[iCell]):
                vol_1D[i] = h[iCell,k]*areaCell[iCell]
                density_1D[i] = density[iCell,k]
                i = i+1
        nCells_1D = i
        # --- Density sorting in ascending order
        sorted_ind = np.argsort(density_1D)
        density_sorted = np.zeros(nCells_1D)
        vol_sorted = np.zeros(nCells_1D)
    
        density_sorted = density_1D[sorted_ind]
        vol_sorted = vol_1D[sorted_ind]

        for j in range(len(sorted_ind)):
            density_sorted[j] = density_1D[sorted_ind[j]]
            vol_sorted[j] = vol_1D[sorted_ind[j]]
    
        rpe1 = np.zeros(nCells_1D)
        sillWidth = np.zeros(nCells_1D)
        yWidth = np.zeros(nCells_1D)
        zMid = np.zeros(nCells_1D)
        z = 0.0

        # --- RPE computation
        for i in range(nCells_1D):
            yWidth[i] = yMax - yMin
            area = yWidth[i]*xMax
            thickness = vol_sorted[i]/area
            zMid[i] = z-thickness/2.0
            z = z-thickness
            rpe1[i] = gravity*density_sorted[i]*(zMid[i]+bottomMax)*vol_sorted[i]

        rpe[nt] = np.sum(rpe1)/np.sum(areaCell)
        keMeanVolume[nt] = np.mean(kineticEnergyCellFull[nt,:,:])
        vertTransportVolume[nt] = np.mean(vertTransportFull[nt,:,:])
    
        print(rpe[nt])
    
    rpeNorm = (rpe-rpe[0])/rpe[0]
    
    print(rpeNorm)
    
    ncfile_init.close()
    ncfile_out.close()
    
    # --- Write in text
    file1 = open('rpe_'+str(nf+1)+'.txt','w')
    for nt in range(T):
        file1.write(str(rpeNorm[nt])+"\n")
    file1.close()
