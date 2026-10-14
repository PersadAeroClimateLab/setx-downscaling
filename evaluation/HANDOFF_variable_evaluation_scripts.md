# Handoff — SETx 1 km CMIP6 Variable-Evaluation Scripts

> Written in ASD-STE100 Simplified Technical English. Short sentences. Active voice. One instruction per line.
>
> **Purpose of this document.** This document hands off the variable-evaluation code for the SETx 1 km CMIP6 downscaling project. A new reader can start the work cold. The document maps every file, every rule, and every known trap.
>
> **Place in this archive.** In this archive the evaluation code sits in `03_evaluation/`, and the citation, license, and requirements files are shared at the archive root rather than repeated here. The site-specific absolute paths have been replaced with `/PATH/TO/...` placeholders; set them for your machine. `README.md` in this directory is the entry point; this handoff is the design reference behind it.

---

## 1. Scope

This handoff covers the **evaluation scripts only**. These scripts evaluate the downscaling skill for 8 climate variables. They are the `*_publication.py` suite, one self-contained script per variable.

This handoff does **not** cover:
- The Zarr QC pipeline (`run_qc_pipeline.sh` and its helpers).
- The two manuscript drafts (ESSD and Climate Dynamics).

The code state in this archive is the **v1.0 release**.

---

## 2. Project summary

- **Project.** Statistical downscaling evaluation. The method is 1 km dual-kriging (KrigR) applied to NEX-GDDP-CMIP6 model output over Southeast Texas.
- **Domain (SETx).** 28.5–33.75° N, 93.5–97.5° W.
- **Period.** Historical 1980–2014. Projections extend to 2100 under SSP scenarios.
- **Variables (8).** `tas`, `tasmax`, `tasmin`, `hurs`, `huss`, `sfcWind`, `rlds`, `rsds`.
- **Models (10).** CESM2, CMCC-ESM2, CNRM-ESM2-1, MRI-ESM2-0, NorESM2-MM, EC-Earth3, FGOALS-g3, GFDL-CM4, BCC-CSM2-MR, MPI-ESM1-2-HR.
- **Experiments (5).** historical, ssp126, ssp245, ssp370, ssp585.
- **Reference datasets.** AORC (1 km gridded reference), ERA5 native, ERA5 downscaled (ERA5-DS), Daymet, GHCN-D station observations, ESA WorldCover (land cover).
- **Group.** Persad Research Group, TACC.

---

## 3. File manifest

There is one self-contained script per variable. All eight sit in the
evaluation directory, together with this handoff and a local `README.md`.

| File | Role |
|---|---|
| `tas_publication.py` | Daily mean temperature. |
| `tasmax_publication.py` | Daily maximum temperature. |
| `tasmin_publication.py` | Daily minimum temperature. |
| `hurs_publication.py` | Relative humidity. Larger than the others (see note). |
| `sfcwind_publication.py` | 10 m wind speed. **Note the lowercase filename.** |
| `huss_publication.py` | Specific humidity. No station obs. |
| `rlds_publication.py` | Downward longwave radiation. No station obs. |
| `rsds_publication.py` | Downward shortwave radiation. No station obs. |

The citation, license, and requirements files are shared at the archive root,
not repeated here.

**Note on `hurs`.** The `hurs` script is about 70 lines longer than the others. This is correct. It derives a Daymet relative-humidity field at runtime via the Magnus formula, cached to a Zarr store. Do not treat the extra lines as a cleanup miss.

---

## 4. What differs between the scripts

The scripts share the same figure and metric logic. Each one sets, for its
variable:

- File paths (Zarr root, ERA5 files, AORC file, obs CSV, output directory).
- The Zarr glob pattern (the variable name inside the CMIP6 Zarr stores).
- Unit labels (for example °C, m/s, %, W/m², kg/kg).
- QC physical bounds (per variable).
- Scatter-plot axis overrides.
- Equivalence-band thresholds.
- Daymet handling. `hurs` derives a Magnus-formula RH field at runtime, cached
  to a Zarr store; the others read a Daymet file directly or skip Daymet.
- The `HAS_STATION_OBS` flag (see §7).

---

## 5. Checking the scripts

Confirm every script compiles before a run:

```bash
for f in *_publication.py; do python3 -c "import ast; ast.parse(open('$f').read()); print('$f OK')"; done
```

---

## 6. How to run one script

### 6.1 Configure the paths

Each script has a config block near the top. The block sets the variable, the data paths, and the output directory. Example (wind speed):

```python
LONGNAME = "Daily Average Wind Speed at 10m (m/s)"
ZARR_DIR        = Path("/PATH/TO/ZARR_ROOT/sfcWind")
ERA5_PATH       = Path(".../era5/sfcWind_era5.nc")
ERA5_DS_PATH    = Path(".../downs_eval/ERA5_sfcWind_1km_1980-2014.nc")
ERA5_DS_SE_PATH = Path(".../downs_eval/SE_ERA5_sfcWind_1km_1980-2014.nc")
AORC_PATH       = Path("/PATH/TO/AORC/AWND_10maboveground.nc")
OBS_CSV         = Path(".../GHCND_AWND_1980_2020_SETx.csv")
STN_META        = Path(".../SETX_stations.csv")
GHCND_INVENTORY = Path(".../ghcnd-stations.txt")
BOUNDARY_SHP    = Path(".../setx_boundary.shp")
LULC_PATH       = Path(".../worldcover/worldcover_setx.tif")
STATE_SHP       = Path(".../us_states.shp")   # fallback for the station map
OUT_DIR         = Path(".../outputs/sfcWind_publication")
```

> **All paths are site-specific.** They point to TACC and the local server. The script headers flag this. Change the paths for any other machine. Path hoisting to a config file is deferred (see §12).

### 6.2 Run a quick test

Use demo mode to test the script on a small subset. This is fast.

```bash
python3 tas_publication.py --demo
```

### 6.3 Run the full evaluation

Run the script without the flag to produce all figures and CSVs.

```bash
python3 tas_publication.py
```

Figures and CSVs go to `OUT_DIR`.

### 6.4 Notebook use

Each figure calls a `save_and_show()` helper. The helper runs `fig.savefig()`, then `plt.show()`, then `plt.close()`. This supports the VS Code / JupyterHub notebook workflow.

---

## 7. Figure catalog and the station-obs split

The figures follow **manuscript order**. Each script holds 28 figure functions. These map to 6 main figures, 3 optional subsection figures, and 7 SI (supplementary) figure groups.

### 7.1 The `HAS_STATION_OBS` flag

Some variables have direct GHCN-D station observations. Some do not. The `HAS_STATION_OBS` flag controls this.

| Group | Variables | `HAS_STATION_OBS` | Figures | CSVs |
|---|---|---|---|---|
| Station variables | `tas`, `tasmax`, `tasmin`, `hurs`, `sfcWind` | `True` | **12** | **4** |
| Non-station variables | `huss`, `rlds`, `rsds` | `False` | **9** | **2** |

For a non-station variable, the runtime branch skips all GHCN-D station figures and CSVs. The script header for those three states "9 figures" and marks the skipped items as `NOT PRODUCED - no GHCN-D obs`.

> **Accuracy trap (fixed).** An earlier header claimed "12 figures" on every script. That was false for the 3 non-station variables. Those three scripts now state the correct count (9). Keep this correct. A wrong count is a false metadata claim.

### 7.2 Known figure details

- **Fig 3** — Station-validation summary. A 1×3 row. Panel (a) label: `"Gridded Daily Mean Temperature vs GHCN-D Stations"`.
- **Fig 5** — Downscaling Added Value. A 2×2 layout. Panel (b) shows the KGE component decomposition (r, β, γ). Skill thresholds 0.41 and 0.75 appear on panel (c). Panel (d) has improvement-zone shading.
- **Taylor diagrams, KGE, CRMSE, bootstrap CI, per-pixel r p-value (`pix_r_pval`)** are standard metrics across the main and SI figures.
- **SI Fig A** uses an adaptive colorbar that works across variable magnitudes, from kg/kg to W/m².

---

## 8. GHCN-D station observation inputs

The station obs come from GHCN-D CSV files. The CSV columns are:

`station_id, station_name, latitude, longitude, elevation_m, state, date, value_raw, value_converted, unit, mflag, qflag, sflag`

- The `load_obs` function auto-detects the value column. It uses `value_converted`, which is already in the correct unit.
- Variable to CSV mapping:
  - `tas` ← TAVG CSV. For full 1980–2014 coverage, use (TMAX + TMIN)/2 (see §10).
  - `hurs` ← RHAV CSV.
  - `sfcWind` ← AWND CSV.

---

## 9. Data inputs and paths

| Item | Path |
|---|---|
| TACC project space | `/PATH/TO/PROJECT_SPACE/` |
| TACC scratch | `/PATH/TO/SCRATCH/` |
| Project data | `/PATH/TO/EVAL_DATA/` |
| Local server | thunderhead-4, `/PATH/TO/LOCAL/` |
| Primary Zarr root | `/PATH/TO/ZARR_ROOT` |
| sfcWind POST Zarr root | `/PATH/TO/ZARR_ROOT_POST` |
| Output root (pattern) | `.../outputs/<variable>_publication` |

**Data formats.** Zarr is the primary analysis format. NetCDF-4 is the source/reference format.

**Core libraries.** xarray, zarr, numpy, matplotlib, cartopy, pandas, scipy, dask. CDO is used outside Python.

---

## 10. Key technical learnings — read before you touch the code

Each item below is a real trap found during development. Keep each fix.

1. **cftime / noleap calendars.** CMIP6 uses noleap calendars. Open the data with `use_cftime=True`. Convert times element-wise with `pd.Timestamp(t.year, t.month, t.day)`. A Feb 28 → Mar 1 "gap" on a leap year is expected noleap behavior. It is not a data error.

2. **ERA5 radiation units.** ERA5-derived `rsds` and `rlds` store hourly accumulation in J/m²/hr. Divide by 3600 to get W/m². Apply this to ERA5 native, ERA5-DS, and ERA5-DS SE. Warning: the ERA5-DS metadata may wrongly claim the unit is already W m⁻². Trust the division, not the metadata.

3. **Cross-resolution comparison.** The ERA5 native grid (about 25×23) does not match the AORC grid (630×480). Do not compare them pixel by pixel. Use domain-averaged time series for any ERA5-native vs AORC comparison. A direct spatial compare raises a grid-mismatch error.

4. **MMM outlier handling.** The multi-model mean (MMM) must auto-detect and drop outlier models. Filter the model list (`MF`) **after** the load loop. If you keep a skipped model (for example CESM2) in the list, you get a `KeyError`.

5. **SE variable naming.** The per-model standard-error field in the Zarr store uses the pattern `{var}_SE` (for example `rsds_SE`, `hurs_SE`). Put `{var}_SE` as the **first** candidate in the SE name search. Case matters.

6. **K → °C conversion.** The temperature scripts (`tas`, `tasmax`, `tasmin`) convert Kelvin to Celsius. Guard this so it does not trigger falsely on data that is already in °C.

7. **GHCN-D coverage limits.** TAVG is reliable only from about 1991. RHAV is reliable only from about 2006 (ASOS / automated stations). For full 1980–2014 `tas` station coverage, use (TMAX + TMIN)/2, not TAVG.

8. **WorldCover land cover.** The 10 m WorldCover TIF is large. Downsample it at read time with GDAL `ReadAsArray(buf_xsize, buf_ysize)` (max about 2000 px per side) to avoid out-of-memory. Use `shading="nearest"` for the discrete class data.

9. **Coastline overlay.** The high-resolution coastline uses Cartopy Natural Earth 10 m data. A 92-waypoint Gulf Coast fallback is built in for offline or failed-download cases.

10. **Physical plausibility QC.** Clip gridded data to per-variable physical bounds. **Clip, do not NaN-mask.** NaN-masking creates holes in the fields.

11. **AST-based verification.** Do not judge "load-bearing" code by eye. Manual inspection is unreliable. Use Python's `ast` module to verify which variables are actually consumed downstream before you delete anything.

12. **Comment cleanup safety.** To strip comments, use Python's `tokenize` module, not a line-based `#` detector. A line-based detector corrupts lines that hold a `#` inside a string literal. After any such cleanup, prove the executable code is byte-identical and recompile all files.

---

## 11. Conventions

- **Communication.** Short, direct ("caveman speak") by default. Give detail only when asked.
- **Delivery.** Provide complete, copy-pasteable scripts. Do not provide partial diffs. Every script must compile clean and run as is.
- **Figures.** Output every figure as PDF and PNG at 1200 DPI. Use serif fonts (Times New Roman).
- **Figure reference.** Refer to a figure by its output number. The number maps to a specific function.
- **Consistency.** The scripts share the same figure and metric logic. When you change that shared logic, apply the same change to every script so the figures stay comparable across variables.

---

## 12. Release and open items

### Release state (v1.0)
- One self-contained script per variable.
- Code license: MIT. Figure license: CC-BY-4.0.
- Author, ORCID, and DOI fields are `PLACEHOLDER`. Fill them when confirmed.

### Open items to finish
1. **Verify `requirements.txt`.** The current version pins are invented, not measured. Run `pip freeze` on the real conda environment. Replace the pins with the true values before release.
2. **Fill the placeholders.** Replace every `PLACEHOLDER` (author, ORCID, DOI) once the details are confirmed.
3. **Path hoisting (deferred).** The TACC absolute paths are inline and site-specific. A move to a config file is deferred. Do it if the next project needs portability.

---

## 13. Quick-start checklist for the new project

1. Copy the evaluation files into the new project.
2. Edit the configuration block in each script for the new machine paths.
3. Compile-check all 8 scripts (see §5).
4. Run one script in demo mode (`--demo`) to confirm the data loads.
5. Read §10 before you change any shared logic.
6. Apply any shared-logic change to every script, so the figures stay comparable across variables.

---

*End of handoff.*
