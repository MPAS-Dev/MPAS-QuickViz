# Water Mass Analysis

Development space for code and visualization related to water mass analysis in MPAS-Ocean.

The primary method used in this analysis is the Walin air-sea flux framework ([Walin 1982, Tellus](https://doi.org/10.3402/tellusa.v34i2.10801); [Speer and Tziperman 1992](https://doi.org/10.1175/1520-0485(1992)022<0093:ROWMFI>2.0.CO;2)). Briefly, the mode water transformation $F$ and formation $M$ are defined in terms of the net air/sea density flux $f$ area-integrated over an isopycnal outcrop.

$$f = -\frac{\alpha\Phi_{\text{\heat}}}{C_p} + \beta S\Phi_{\text{\fresh}} \hspace{100pt} F(\rho) = \frac{\partial}{\partial\rho}\iint_{A}f\delta(\rho-\rho_0)dA \hspace{100pt} M(\rho) = - \frac{\partial F}{\partial\rho}$$

Below is a summary of the code that handles each step of the WMT calculations:

   * density calculations [`modules/aggregate_mpas_2Dvariables.py`](modules/aggregate_mpas_2Dvariables.py#L174)
   * surface flux definitions [`yaml/variable_combinations.yaml`](yaml/variable_combinations.yaml)
   * buoyancy fluxes [`modules/postprocesstools.py`](modules/postprocesstools.py#L116)
   * outcrop area-integration [`modules/postprocesstools.py`](modules/postprocesstools.py#L137)

### Postprocessing

This code space contains several functions and executable modules written to handle the MPAS-Ocean postprocessing steps for performing the WMT analysis:

[`modules/aggregate_mpas_2Dvariables.py`](modules/aggregate_mpas_2Dvariables.py)

Executable module to aggregate monthly MPAS-Ocean results files. Aggregation paths are stored in the [`yaml`](yaml) directory. To aggregate the `EC30to60E2r2` (LR) CORE-II results:

```
$ cd modules
$ python aggregate_maps_2Dvariables.py ../yaml/paths_LR.yaml
```

[`modules/build_mpas_wmtrtimeseries.py`](modules/build_mpas_wmtrtimeseries.py)

Executable module to spatially average the aggregated MPAS-Ocean files over regions, or over isopycnal outcrops for WMT quantities (if `-c` flag is used). `-v` and `-r` flags are also provided for requesting specific variables or regions, otherwise all are included.

```
$ cd modules
$ python build_mpas_wmtrtimeseries.py -c /path/to/aggregated_results_file.nc
```

[`modules/build_mpas_2Dtimeavgs.py`](modules/build_mpas_2Dtimeavgs.py)

Executable module to compute monthly average climatologies remapped to 2D lon-lat grids from the aggregated MPAS-Ocean files. 2D fields of WMT variables are also returned if the `-c` flag is used (following [Brambilla et al. 2008, JGR Oceans](https://doi.org/10.1029/2006JC004063); [Maze et al. 2009, JPO](https://doi.org/10.1175/2009JPO3985.1)). The `-t` flag is provided to specify the time range for the climatology, and `-v` and `-b` flags are provided for requesting specific variables or a custom bounding box. Defaults include all variables and times in the record, and the default bounding box returns the full North Atlantic region. To process the first decade of a simulation including WMT variables:

```
$ cd modules
$ python build_mpas_2Dtimeavgs.py -c -t 19470101,19561231 /path/to/aggregated_results_file.nc
```

### Visualization

Most useful visualization code eventually ends up in [`modules/visualizationtools.py`](modules/visualizationtools.py). These visualizations are summarized in two notebooks:

   * [`notebooks/mpasoLRHR_timeseries_19471956_19972006.ipynb`](notebooks/mpasoLRHR_timeseries_19471956_19972006.ipynb) Classic sigma-space WMT plots, time series plots, variable correlations
   * [`notebooks/mpasoLRHR_surfacefields_monthlyclimatology_19471956_19972006.ipynb`](https://nbviewer.org/github/MPAS-Dev/MPAS-QuickViz/blob/master/ocean/AMOC/watermassanalysis/notebooks/mpasoLRHR_surfacefields_monthlyclimatology_19471956_19972006.ipynb) Monthly climatologies of 2D fields (animated, [nbviewer](https://nbviewer.org/) required)