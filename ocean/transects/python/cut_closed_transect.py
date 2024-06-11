#!/usr/bin/env python

import argparse

import numpy as np
from geometric_features import (
    GeometricFeatures,
    read_feature_collection
)
from shapely.geometry import (
    mapping,
    Polygon,
    shape
)


def cut_transect(fc_transect, lat, lon, size, out_filename):
    """
    Cut a square out of the given closed-loop transect to break the loop.
    """

    # find the closest point in the transect to the specificed lat/lon

    feature = fc_transect.features[0]
    coordinates = feature['geometry']['coordinates']
    feature_type = feature['geometry']['type']
    if feature_type == 'LineString':
        coordinates = [coordinates]
    elif feature_type != 'MultiLineString':
        raise ValueError(
            f'Unexpected geometry type for transect {feature_type}')

    min_dist = None
    center_lon = None
    center_lan = None
    for coords in coordinates:
        lon_local, lat_local = zip(*coords)
        dist = np.sqrt((np.array(lon_local) - lon)**2 +
                       (np.array(lat_local) - lat)**2)
        index_min = np.argmin(dist)
        if min_dist is None or dist[index_min] < min_dist:
            center_lon = lon_local[index_min]
            center_lan = lat_local[index_min]
            min_dist = dist[index_min]

    square = Polygon([(center_lon - 0.5 * size, center_lan - 0.5 * size),
                      (center_lon - 0.5 * size, center_lan + 0.5 * size),
                      (center_lon + 0.5 * size, center_lan + 0.5 * size),
                      (center_lon + 0.5 * size, center_lan - 0.5 * size),
                      (center_lon - 0.5 * size, center_lan - 0.5 * size)])

    feature = fc_transect.features[0]
    transect_shape = shape(feature['geometry'])
    transect_shape = transect_shape.difference(square)

    # now sort the coordinates so the start and end of the transect are at the
    # dividing point

    feature['geometry'] = mapping(transect_shape)

    feature_type = feature['geometry']['type']
    if feature_type == 'MultiLineString':
        coordinates = feature['geometry']['coordinates']

        # reorder the LineStrings so the first one starts right after the cut

        closest = None
        min_dist = None
        for index, coords in enumerate(coordinates):
            lon_first, lat_first = coords[0]
            dist = np.sqrt((lon_first - lon)**2 + (lat_first - lat)**2)
            if min_dist is None or dist < min_dist:
                closest = index
                min_dist = dist
        new_coords = list(coordinates[closest:])
        new_coords.extend(list(coordinates[:closest]))
        feature['geometry']['coordinates'] = tuple(new_coords)

    fc_transect.to_geojson(out_filename)


def main():

    parser = argparse.ArgumentParser(description='''
        cut the given transect loop as close as possible to the given
        latitude and longitude''')
    parser.add_argument('-g', dest='geojson_filename',
                        help='Geojson filename with transect', required=False)
    parser.add_argument('-f', dest='feature_name',
                        help='Name of an ocean transect from '
                             'geometric_features',
                        required=False)
    parser.add_argument('--lat', dest='lat', type=float,
                        help='The approx. latitude at which to cut the loop',
                        required=True)
    parser.add_argument('--lon', dest='lon', type=float,
                        help='The approx. longitude at which to cut the loop',
                        required=True)
    parser.add_argument('--size', dest='size', type=float,
                        help='The size in degrees of the square used to cut '
                             'the loop',
                        required=True)
    parser.add_argument('-o', dest='out_filename',
                        help='The geojson file with the cut transect to write '
                             'out',
                        required=True)
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

    cut_transect(fc_transect, args.lat, args.lon, args.size, args.out_filename)


if __name__ == '__main__':
    main()
