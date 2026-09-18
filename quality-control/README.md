# Quality Control

Two independent approaches are used to audit the SETx-1km-CMIP6 downscaled Zarr stores in `/local1/SETx_UIFL_Zarr`: a packaged, automated checker (`SETx-UIFL-Downscaled-Quality-Control`) and a chain of standalone scripts for a deeper, manual audit (`full-scan-pipeline`).

## `SETx-UIFL-Downscaled-Quality-Control` (git submodule)

A self-contained, testable QC package. By default, most `git clone` commands will not initialize the submodule. Use `git submodule update --init` if it's not already populated. This package runs a fixed set of checks per store (catalog-driven with one `dask.compute` per store):

- **Missing data** via the store's own modal (most-common) footprint, not a fixed land mask.
- **Physical range** per-variable bounds, plus the `tasmin <= tas <= tasmax` consistency check.
- **Spatial artifacts**: banding, roughness, and repeated extreme values.

Results are written as a per-store JSON (status + scalar checks) and NetCDF (per-day time series) under `--out`, then summarized into a single `summary.csv`. It ships with a Dockerfile/`requirements.txt` (convertible to an Apptainer `.sif` for the target machine), a unit-tested `tests/` suite run against synthetic stores with injected faults, and `docs/llm/README.md` for use with an LLM coding agent before. See the submodule's own `README.md` for usage.

Use this when you want a repeatable, script-based pass/fail check with a stable spec.

## `full-scan-pipeline`

A chain of standalone scripts, run in sequence by `run_qc_pipeline.sh`, that scans every store far more exhaustively and ends in a curated list of bad dates plus per-date plots for a human to review and decide clip vs. re-downscale. Config (paths, worker count, thresholds) is overridden via environment variables (see the top of `run_qc_pipeline.sh`).

Expected outputs: `qc_zarr_out/SUMMARY.txt`, `qc_drilldown_out/DRILLDOWN_SUMMARY.txt`, `issues_by_variable/<var>_issues.csv`, `footprint_out/footprint_classification.csv` + `footprint_maps.pdf`, `bad_dates.csv`, and `bad_date_plots/<variable>/<model>_<scenario>_<date>.png`.

Use this for a comprehensive, exploratory audit. Its output is meant to be thoroughly reviewed by a person, not consumed as an automated pass/fail signal.
