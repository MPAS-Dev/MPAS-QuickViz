#!/usr/bin/env python

import argparse

import numpy as np
import xarray as xr
from geometric_features import read_feature_collection
from geometric_features import GeometricFeatures
from mpas_tools.cime.constants import constants
from mpas_tools.logging import LoggingContext
from mpas_tools.io import write_netcdf
from mpas_tools.mesh.mask import compute_mpas_transect_masks
from mpas_tools.parallel import create_pool

from transect.vert import compute_transect


def combine_transect_datasets(ds_mesh, fc_transect, out_filename, pool,
                              logger):
    """
    Combine transects masks on cells and edges with a dataset for plotting
    that describes how the transect slices through cell and edge geometry.
    Add fields on edges and cells that define the (mean) distance along the
    transect for each cell or edge in the transect
    """

    earth_radius = constants['SHR_CONST_REARTH']

    ds_mask = compute_mpas_transect_masks(dsMesh=ds_mesh, fcMask=fc_transect,
                                          earthRadius=earth_radius,
                                          maskTypes=('cell', 'edge',),
                                          logger=logger,
                                          pool=pool, addEdgeSign=True)

    feature = fc_transect.features[0]
    geom_type = feature['geometry']['type']
    if geom_type == 'LineString':
        coordinates = [feature['geometry']['coordinates']]
    elif geom_type == 'MultiLineString':
        coordinates = feature['geometry']['coordinates']
    else:
        raise ValueError(
            f'Unexpected geometry type for the transect {geom_type}')

    lon = []
    lat = []
    for coords in coordinates:
        lon_local, lat_local = zip(*coords)
        lon.extend(lon_local)
        lat.extend(lat_local)
    lon = xr.DataArray(data=lon, dims='nTransectPoints')
    lat = xr.DataArray(data=lat, dims='nTransectPoints')

    layer_thickness = ds_mesh.layerThickness
    bottom_depth = ds_mesh.bottomDepth
    min_level_cell = ds_mesh.minLevelCell
    max_level_cell = ds_mesh.maxLevelCell

    ds_transect = compute_transect(lon, lat, ds_mesh, layer_thickness,
                                   bottom_depth, min_level_cell,
                                   max_level_cell, spherical=True)

    ds = ds_mask
    for var in ds_transect.data_vars:
        ds[var] = ds_transect[var]

    add_distance_field(ds, logger)

    write_netcdf(ds, out_filename)


def add_distance_field(ds, logger):
    """
    Add fields on edges and cells that define the (mean) distance along the
    transect for each cell or edge in the transect
    """

    dist_cell = np.zeros(ds.sizes['nCells'])
    count_cell = np.zeros(ds.sizes['nCells'], dtype=int)
    dist_edge = np.zeros(ds.sizes['nEdges'])
    count_edge = np.zeros(ds.sizes['nEdges'], dtype=int)

    logger.info('Adding transect distance fields on cells and edges...')

    for segment in range(ds.sizes['nSegments']):
        icell = ds.horizCellIndices.isel(nSegments=segment).values
        iedge = ds.horizEdgeIndices.isel(nSegments=segment).values
        # the distance for the midpoint of the segment is the mean
        # of the distances of the end points
        dist = 0.5 * (ds.dNode.isel(nHorizNodes=segment) +
                      ds.dNode.isel(nHorizNodes=segment + 1))
        dist_cell[icell] += dist
        count_cell[icell] += 1
        dist_edge[iedge] += dist
        count_edge[iedge] += 1

    mask = count_cell > 0
    dist_cell[mask] /= count_cell[mask]
    dist_cell[np.logical_not(mask)] = np.nan

    mask = count_edge > 0
    dist_edge[mask] /= count_edge[mask]
    dist_edge[np.logical_not(mask)] = np.nan

    ds['transectDistanceCell'] = ('nCells', dist_cell)
    ds['transectDistanceEdge'] = ('nEdges', dist_edge)
    logger.info('Done.')


def main():

    parser = argparse.ArgumentParser(description='''
        creates transect edge and cell masks along with edge sign and distance
        along the transect''')
    parser.add_argument('-m', dest='mesh_filename',
                        help='MPAS-Ocean horizontal and vertical filename',
                        required=True)
    parser.add_argument('-g', dest='geojson_filename',
                        help='Geojson filename with transect', required=False)
    parser.add_argument('-f', dest='feature_name',
                        help='Name of an ocean transect from '
                             'geometric_features',
                        required=False)
    parser.add_argument('-o', dest='out_filename',
                        help='Edge transect filename', required=True)
    parser.add_argument(
        "--process_count", required=False, dest="process_count", type=int,
        help="The number of processes to use to compute masks.  The "
             "default is to use all available cores")
    parser.add_argument(
        "--multiprocessing_method", dest="multiprocessing_method",
        default='forkserver',
        help="The multiprocessing method use for python mask creation "
             "('fork', 'spawn' or 'forkserver')")
    args = parser.parse_args()

    if args.geojson_filename is None and args.feature_name is None:
        raise ValueError('Must supply either a geojson file or a transect '
                         'name')

    if args.geojson_filename is not None:
        fc_transect = read_feature_collection(args.geojson_filename)
    else:
        gf = GeometricFeatures()
        fc_transect = gf.read(componentName='ocean', objectType='transect',
                              featureNames=[args.feature_name])

    ds_mesh = xr.open_dataset(args.mesh_filename)
    if 'Time' in ds_mesh.dims:
        ds_mesh = ds_mesh.isel(Time=0)

    pool = create_pool(process_count=args.process_count,
                       method=args.multiprocessing_method)

    with LoggingContext('create_transect_masks') as logger:

        combine_transect_datasets(ds_mesh, fc_transect, args.out_filename,
                                  pool, logger)


if __name__ == '__main__':
    main()
