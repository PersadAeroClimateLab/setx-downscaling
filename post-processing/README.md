# Post Processing

The raw SETx netCDF files are output as individual timesteps per variable, per model, per scenario. The full raw collection is then an unweildly 15 million files, therefore post-processing was applied to concatenate these individual netCDF files into chunked zarr stores.

## Metadata Fixes

The raw data is metadata-poor. We generate target Zarr stores with pre-filled metadata:

```python
build_template_zarrs("rsds", "/local1/SETx_UIFL_Zarr/")
```

## Swapping Out Timesteps

Due to the high number of fileheads, many of the file writes were either fully or partially corrupted. An analysis of the zarr store using the `scan_zarr` was used to identify problematic timesteps for re-computation. 

```python
orders = []
target_zarr_dir = "/local1/SETx_UIFL_Zarr/rsds"
for zarr_store in listdir(target_zarr_dir):
    orders.append(dask.delayed(scan_zarr)(f"{target_zarr_dir}/{zarr_store}", "rsds"))
    # orders.append(dask.delayed(scan_zarr)(f"{target_zarr_dir}/{zarr_store}", "rsds_SE"))

nan_report = dask.compute(*orders)
```

The resulting `nan_report` dictionary is then used to recompute the timesteps. These re-computed netCDF files are still in the individual timestep form as the raw data, but the filename contains the corrected date after some updates to the original Kriging pipeline code. Therefore, a separate suite of functions in `data_cleaning_funcs.py` were used to "swap out" those timesteps. Example:

```python
from data_cleaning_funcs import *

target_zarr_dir = "/local1/SETx_UIFL_Zarr/rsds"
slice_orders = get_slice_orders("/local1/SETx_Swapping_temp/new_redown_rsds/1km_files", target_zarr_dir)
orders = []
for template_zarr_path in slice_orders:
    swap_args = slice_orders[template_zarr_path]
    orders.append(dask.delayed(swap_slice)(template_zarr_path, swap_args))
```

This can be parallelized with Dask:

```python
import dask
from dask.distributed import Client, LocalCluster
cluster = LocalCluster(n_workers=48, memory_limit="4GB", threads_per_worker=4, dashboard_address=':8002')
client = Client(cluster)

results = dask.compute(*orders)
```