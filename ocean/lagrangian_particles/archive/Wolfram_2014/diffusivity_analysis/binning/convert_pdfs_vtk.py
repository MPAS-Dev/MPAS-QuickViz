#!/usr/bin/env python

import numpy as np
from pyvtk import VtkData, UnstructuredGrid, Scalars, PointData
from matplotlib.tri import Triangulation

def convert_layers_vkt(datadir, nlayers):
    for alayer in np.arange(nlayers-1):
        print 'setting up vtk file for %d layer of %d'%(alayer, nlayers-1)
        datapairs = np.load(datadir + '/pairs_layer%d.npz'%(alayer))
        datacluster = np.load(datadir + '/pairs_layer%d.npz'%(alayer))
        datadiff    = np.load(datadir + '/diff_layer%d.npz'%(alayer))
        # data cursor
        data = datapairs
        x = np.degrees(data['x'])
        y = np.degrees(data['y'])
        npoints = data['x'].shape[0]
        bins = data['bins']
        bins = 0.5*(bins[:-1]+bins[1:])/100.
        nbins = bins.shape[0]
        Nt = data['te']-data['ts']
        
        triangt = Triangulation(x,y)
        eones = np.ones((nbins,npoints))
        xp = x[np.newaxis,:]*eones
        yp = y[np.newaxis,:]*eones
        zp = bins[:,np.newaxis]*eones
        points = [[ax,ay,ab] for ax,ay,ab in zip(xp.ravel(), yp.ravel(), zp.ravel())]
        wedges = []
        # 3D surfaces
        for zlayer in np.arange(nbins-1):
            for atrit in zip(triangt.triangles):
                wedges.append(np.squeeze(np.hstack((atrit + zlayer*npoints, atrit + (zlayer+1)*npoints))).tolist())

        exclusions = ['bins']
        for tout in np.arange(Nt):
            print str(tout) + ' of ' + str(Nt)

            # volume
            print 'building volume vtk'
            scalars=[]
            def sanitize_data_vol(data):
                return np.nan_to_num(data.ravel(order='F'))
            def sanitize_data_sur(data):
                return np.nan_to_num((data[:,np.newaxis]*np.ones((data.shape[0],nbins))).ravel(order='F'))

            scalars.append(Scalars(np.linalg.norm(np.asarray(points),axis=1), name='test',lookup_table='default'))
            for data in [datacluster, datadiff, datapairs]:
                for df in data.files:
                    if not '_full' in df and not df in exclusions:
                        dimension = len(data[df].shape)
                        dfdata = None
                        if __debug__:
                            print df, dimension, data[df].shape
                        if dimension == 1:
                            dfdata = sanitize_data_sur(data[df][:])
                        elif dimension == 2:
                            if tout < data[df].shape[1]:
                                dfdata = sanitize_data_sur(data[df][:,tout])
                        elif dimension == 3:
                            if tout < data[df].shape[2]:
                                dfdata = sanitize_data_vol(data[df][:,:,tout])
                        else: 
                            continue
                        if np.any(dfdata):
                            if __debug__:
                                print df, dfdata.shape
                            scalars.append(Scalars(dfdata, name=df, lookup_table='default'))
            # specialized calculations
            if tout < datapairs['Kxx'].shape[2]: 
                scalars.append(Scalars(sanitize_data_vol(datapairs['Kxx'][:,:,tout] +
                    datapairs['Kyy'][:,:,tout]), name='Krr', lookup_table='default'))
            if tout < datapairs['meanKxx'].shape[1]: 
                scalars.append(Scalars(sanitize_data_sur(datapairs['meanKxx'][:,tout] +
                    datapairs['meanKyy'][:,tout]), name='meanKrr', lookup_table='default'))

            vtk = VtkData(UnstructuredGrid(points, wedge=wedges), PointData(*scalars), 'pairs_layer%d'%(alayer))
            vtk.tofile(datadir+'/layer_vol_%d_%d'%(alayer,tout),'binary')



if __name__ == "__main__":

    from optparse import OptionParser

    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-d", "--datadir", dest="datadir",
                      help="datadir with npz files",
                      metavar="FILE")
    parser.add_option("-l", "--numlayers", dest="nlayers",
                      help="number of layers",
                      metavar="INT")

    options, args = parser.parse_args()

    if not options.datadir:
        parser.error("Must specify datadir for conversion.")
    if not options.nlayers:
        options.nlayers = 5
    else:
        options.nlayers = int(options.nlayers)
   
    convert_layers_vkt(options.datadir, options.nlayers)


