# Evaluation (Stage 03)

Project: SETx CMIP6 1 km daily downscaling, Southeast Texas

This directory holds the variable-evaluation code. The scripts measure how well
the 1 km dual-kriging (KrigR) downscaling of NEX-GDDP-CMIP6 reproduces observed
climate over Southeast Texas, variable by variable, and produce the manuscript
figures and metric tables. Each variable is compared against gridded references
(AORC, Daymet, ERA5 and the downscaled ERA5) and, where the variable has a
direct in-situ measurement, against GHCN-D station observations.

## The scripts

There is one self-contained script per variable. Each sets its variable, data
paths, units, physical bounds, and output directory in a configuration block
near the top, then runs the full evaluation for that variable.

| File | Variable | Figures / tables |
|------|----------|------------------|
| `tas_publication.py` | Daily mean temperature | 12 figures, 4 CSVs |
| `tasmax_publication.py` | Daily maximum temperature | 12 figures, 4 CSVs |
| `tasmin_publication.py` | Daily minimum temperature | 12 figures, 4 CSVs |
| `hurs_publication.py` | Relative humidity (derives a Daymet RH field via the Magnus formula) | 12 figures, 4 CSVs |
| `sfcwind_publication.py` | 10 m wind speed (note the lowercase filename) | 12 figures, 4 CSVs |
| `huss_publication.py` | Specific humidity (no station obs) | 9 figures, 2 CSVs |
| `rlds_publication.py` | Downwelling longwave radiation (no station obs) | 9 figures, 2 CSVs |
| `rsds_publication.py` | Downwelling shortwave radiation (no station obs) | 9 figures, 2 CSVs |

## Station observations and the figure count

Some variables have a direct GHCN-D station measurement and some do not. The
`HAS_STATION_OBS` flag in each script controls this. For `tas`, `tasmax`,
`tasmin`, `hurs`, and `sfcWind` the flag is true, and the script produces the
full set of 12 figures and 4 CSV tables, including the station-validation
figures. For `huss`, `rlds`, and `rsds` the flag is false: the script skips the
GHCN-D station figures and tables, falls back to AORC as the sole gridded
reference, and produces 9 figures and 2 CSV tables. The headers of those three
scripts state the count of 9 and mark the skipped items.

## How to run

Set the paths for your machine in the script's configuration block (see below),
then run one variable at a time.

```bash
python3 tas_publication.py --demo      # quick test on a small subset
python3 tas_publication.py             # full run: all figures and CSVs to OUT_DIR
```

Figures are written as paired PDF and PNG at 1200 DPI with serif fonts. The
scripts also run inside a Jupyter kernel, where each figure displays inline as it
is produced.

## Paths

The scripts do not carry real site paths. They use three descriptive
placeholders; set them for your machine in each script's configuration block:

- `/PATH/TO/ZARR_ROOT`: the root of the per-model 1 km downscaled Zarr stores,
  one subdirectory per variable.
- `/PATH/TO/EVAL_DATA`: the evaluation data tree, holding the ERA5 and
  downscaled-ERA5 files, the Daymet and NEX-native inputs, the GHCN-D station
  CSVs and metadata, the shapefiles, the WorldCover land cover, and the output
  directories.
- `/PATH/TO/AORC`: the AORC 1 km reference files.

## Environment

The scripts need Python 3.10 or newer with numpy, pandas, xarray, zarr,
netCDF4, scipy, matplotlib, and dask. Cartopy is optional (high-resolution
coastlines; a built-in Gulf Coast fallback covers the offline case), and GDAL or
rioxarray is optional (reading the WorldCover land-cover raster). CDO is used
outside Python. These dependencies are folded into the archive's top-level
`requirements.txt`.
