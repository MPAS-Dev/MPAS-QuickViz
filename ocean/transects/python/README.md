# Python Transect Tools

## compute_transects.py

Computes transport through sections.

Example call:
```
  ./compute_transects.py
  -k transect_masks.nc
  -m MPAS_mesh.nc
  -t 'RUN_PATH/analysis_members/timeSeriesStatsMonthly.*.nc'
  -n 'all'
```
To create the `transect_masks.nc` file, load e3sm-unified and:
```
   MpasMaskCreator.x MPAS_mesh.nc  transect_masks.nc -f transect_definitions.geojson
```
where the `transect_definitions.geojson` file includes a sequence of lat/lon
points for each transect.

On LANL IC, example file is at
```
/usr/projects/climate/mpeterse/analysis_input_files/geojson_files/SingleRegionAtlanticWTransportTransects.geojson
```

## create_transect_masks.py

Requires a conda environment with:
* `python`
* `geometric_features`
* `matplotlib`
* `mpas_tools`
* `netcdf4`
* `numpy`
* `scipy`
* `shapely`
* `xarray`

The tools creates cell and edge masks, distance along the transect of cells
and edges in the mask, and the edge sign on edges.  It also includes
information (distance, cell and edge indices, interpolation weights, etc.)
along the transect itself to aid plotting.

The required inputs are an MPAS mesh file and a geojson file or the name of an
ocean transect from `geometric_features`.  The required output is a filename
with the masks and other information about the transect.

## cut_closed_transect.py

If a transect is a closed loop, the path of edges and edge signs don't work
correctly (the shortest path between the beginning and end of the transect is
trivial and involves a single edge). To avoid this, we provide a tool for
cutting a square (in lat/lon space) out of the transect to sever the loop.
The user provides a latitude and longitude (used to locate the closest point)
on the transect and the size of the square to cut out.

