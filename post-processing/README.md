# Post Processing

The raw SETx netCDF files are output as individual timesteps per variable, per model, per scenario, resulting in an unwieldy collection of roughly 15 million files. To make the data usable, this post-processing step consolidates the individual netCDF files into chunked Zarr stores in four stages:

1. Build empty template Zarr stores with the correct metadata and coordinates.
2. Fill the templates with data from the raw netCDF/Zarr sources.
3. Scan the filled stores for corrupted or missing timesteps.
4. Recompute the bad timesteps and swap them into the store.

## 1. Building Template Zarr Stores

The raw data itself is metadata-poor, so rather than writing into it directly, we first generate empty target Zarr stores pre-filled with the correct coordinates, chunking, and CF-compliant metadata (`build_template_zarrs` in `data_cleaning_funcs.py`):

```python
build_template_zarrs("rsds", "/local1/SETx_UIFL_Zarr/")
```

This creates one Zarr store per model/scenario/time-range combination (`1950-2014` for the `hist` scenario, `2015-2100` for the SSP scenarios) for the given variable, with `time`/`lat`/`lon` coordinates and the correct calendar already set for each model.

## 2. Filling the Templates with Data

Once the templates exist, they need to be populated with data. `data_cleaning_funcs.py` provides two ways to do this, depending on the source:

- `fill_initial_dataset` scans a Zarr directory, locates the matching raw netCDF files in the SETx UIFL products/staging directories, and writes each one into its corresponding template store via `add_full_ds`.
- `write_ds` copies data from an existing, already-built Zarr store into the new template in chunks of 300 timesteps, for cases where the source data is itself already a consolidated Zarr store.

## 3. Identifying Corrupted or Missing Timesteps

Because so many individual file writes were required to build these stores, many of them were fully or partially corrupted. `scan_zarr` checks a Zarr store's variable for timesteps with an abnormal number of NaNs (relative to the most common NaN count across the whole store) and returns any anomalies it finds:

```python
orders = []
target_zarr_dir = "/local1/SETx_UIFL_Zarr/rsds"
for zarr_store in listdir(target_zarr_dir):
    orders.append(dask.delayed(scan_zarr)(f"{target_zarr_dir}/{zarr_store}", "rsds"))
    # orders.append(dask.delayed(scan_zarr)(f"{target_zarr_dir}/{zarr_store}", "rsds_SE"))

nan_report = dask.compute(*orders)
```

`nan_report` is a tuple of `(path, variable, anomalous_timesteps)` results, one entry per store scanned. `export_missing_slices` writes these results to a text file (one `model, variable, scenario, date` line per bad timestep), which is then used to drive re-computation of the affected timesteps.

## 4. Swapping In Recomputed Timesteps

The re-computed netCDF files come out in the same individual-timestep form as the raw data, but their filenames contain the date (after updates to the original Kriging pipeline). A separate suite of functions in `data_cleaning_funcs.py` is used to "swap" these corrected timesteps into the existing Zarr stores in place of the bad ones.

`get_slice_orders` scans a directory of recomputed netCDF files and matches each one to its target Zarr store, building a mapping of `zarr_path -> [(file_path, date, variable), ...]`. `swap_slice` then takes one such mapping entry and writes each listed file into the correct time slice of the target store. For example:

```python
from data_cleaning_funcs import *

target_zarr_dir = "/local1/SETx_UIFL_Zarr/rsds"
slice_orders = get_slice_orders("/local1/SETx_Swapping_temp/new_redown_rsds/1km_files", target_zarr_dir)
orders = []
for template_zarr_path in slice_orders:
    swap_args = slice_orders[template_zarr_path]
    orders.append(dask.delayed(swap_slice)(template_zarr_path, swap_args))
```

This can be parallelized with Dask for maximizing throughput:

```python
import dask
from dask.distributed import Client, LocalCluster
cluster = LocalCluster(n_workers=48, memory_limit="4GB", threads_per_worker=4, dashboard_address=':8002')
client = Client(cluster)

results = dask.compute(*orders)
```
