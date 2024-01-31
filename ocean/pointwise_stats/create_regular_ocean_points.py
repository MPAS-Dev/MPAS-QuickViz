'''
Create a file of regular lat/lon points for a points.geojson file
Mark Petersen
Jan 2024
'''
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

# set parameters:
lonStep = 2
latStep = lonStep
fileName = 'ocean_points_'+str(lonStep)+'x'+str(latStep)+'.geojson'

land_shp_fname = shpreader.natural_earth(resolution='50m',
                                       category='physical', name='land')

land_geom = unary_union(list(shpreader.Reader(land_shp_fname).geometries()))
land = prep(land_geom)

def is_land(x, y):
    return land.contains(sgeom.Point(x, y))

# This gives the expected results for two sample points:

file1 = open(fileName, "w")

file1.write('{\n')
file1.write('    "type": "FeatureCollection",\n')
file1.write('    "groupName": "unspecifiedGroupName",\n')
file1.write('    "features": [\n')

i=0
for lon in range(-180,180, lonStep):
  for lat in range(-90, 89, latStep):
    if not is_land(lon,lat):
      if i > 0:
        file1.write(',\n')
      file1.write('{ "type": "Feature", "properties": { "object": "point","name": "'+str(i)+'"}, "geometry": { "type": "Point", "coordinates": [ '+str(lon)+','+str(lat)+' ] } }')
      i += 1
file1.write('\n')

file1.write(']\n')
file1.write('}\n')

file1.close()

