# Downscaling pipeline

Project: SETx CMIP6 1 km daily downscaling, Southeast Texas

This directory holds the code that downscales NEX-GDDP-CMIP6 model output to
1 km daily fields over Southeast Texas. The method is universal kriging through
KrigR in R, run in two steps: 28 km to 9 km, then 9 km to 1 km. Every output
field is written with a standard-error (SE) field beside it.

Three pipelines live here, all in this one directory:

- Elevation family (tas, tasmin, tasmax, sfcWind), single covariate: the GMTED
  elevation DEM.
- VP family (huss, rsds, rlds), single covariate: huss uses vapor pressure,
  rsds and rlds use solar radiation (srad).
- hurs, the multi-covariate reference pipeline, three covariates: the DEM, the
  model tas field, and a monthly vapor-pressure climatology.

The elevation and VP families are separate copies of the same gated adaptive
engine. Each family has its own driver, R template, jobfile builder, and setup
script, with a name that carries the family label, so a change to one family
cannot touch the other and neither can be confused with the hurs pipeline. The
hurs pipeline is the reference the two gated families were derived from.

## Files in this directory

The scripts come in four roles. Each gated family has a full set; hurs has a
driver, template, and setup.

| File | Role |
|------|------|
| `custom_covariate_downscaling_elev_gated.py` | Elevation-family driver. Runs the two-step kriging per timestep, enforces the coverage gate, routes intermediates to `$TMPDIR`, and guards the R call with a timeout. |
| `custom_covariate_downscaling_vp_gated.py` | VP-family driver, a separate copy of the elevation driver. |
| `custom_covariate_hurs_downscaling.py` | hurs driver. Same two-step and gate logic, plus a result scan for physical bounds (0 to 100 %), for the three-covariate case. |
| `kriging_template_adaptive_elev_gated.R` | Elevation-family kriging template. Four-strategy cascade; emits the coverage tokens the driver reads. |
| `kriging_template_adaptive_vp_gated.R` | VP-family kriging template. |
| `kriging_template_multi_covariate.R` | hurs kriging template, `Data ~ DEM + layer.1 + layer.2`, with a nugget-floor fallback. |
| `create_launcher_jobfile_elev_gated.py` | Builds the TACC launcher jobfile for the elevation family; emits task lines that call the elevation driver. |
| `create_launcher_jobfile_vp_gated.py` | Builds the launcher jobfile for the VP family; emits task lines that call the VP driver. |
| `create_launcher_jobfile_adaptive_fixed.py` | The shared jobfile builder the hurs setup uses. Same per-task structure as the gated builders, but it emits the generic driver name `custom_covariate_downscaling.py`, which the hurs setup renames to the hurs driver. |
| `setup_elev_gated_tasks.sh` | Sets up one elevation variable: writes its gated config and calls the elevation builder. |
| `setup_vp_gated_tasks.sh` | Sets up huss, rsds, and rlds: writes each config (huss covariate = vp, rsds/rlds covariate = srad) and calls the VP builder. |
| `setup_hurs_multi_tasks.sh` | Sets up hurs: writes the hurs config, precomputes the tas monthly means (`cdo timmean`), and builds the jobfile and slurm submit script. |

The subdirectories are `configs/` (the per-variable JSON configs) and `utils/`
(the output-check helper). Both are described below.

A note on the hurs jobfile builder: `setup_hurs_multi_tasks.sh` writes its
launcher jobfile with the builder named in its `JOBFILE_BUILDER` variable,
`create_launcher_jobfile_adaptive_fixed.py` by default. That builder is the
shared, generic one in this directory; it only has to take
`<config.json> <output.job>` and emit one task line per timestep, so it can be
swapped for another that does the same. After the jobfile is written, the setup
normalizes the driver on every task line to `custom_covariate_hurs_downscaling.py`
with `sed`, rewriting the generic name `custom_covariate_downscaling.py` if the
builder emitted that. The bundled `create_launcher_jobfile_adaptive_fixed.py`
does emit the generic name, so that sed step is what points the hurs tasks at the
hurs driver. The gated elevation and VP setups need no such step, because their
builders write the family driver name into each task line directly.

## The coverage gate

The gate is the main thing the drivers add. After each kriging step (9 km, then
1 km) the driver measures coverage before it goes on. Coverage has to reach
`min_coverage` of the area, which is the same as saying NaN stays at or below
`1 - min_coverage`. If a step falls short, the driver deletes every companion
file for that timestep (9 km, SE_9 km, 1 km, and SE_1 km, under either scenario
spelling) and marks the timestep failed.

The default gate is 0.9988 on the domain denominator. The two denominators are:

- domain: the prediction cells that hold a finite covariate (the kriging holes
  only).
- grid: every cell in the target box.

Both the threshold and the denominator live in the config (`min_coverage` and
`coverage_denominator`), so changing either is a config edit, not a code change.

The standalone scanner in `utils/` (`utils/check_1km_coverage.sh`) uses a
different rule on purpose: the grid denominator at 89 %, run by hand after a run.
The production gate runs inside the driver at 99.88 % on the domain denominator,
after every downscale. One is a loose sweep of finished files; the other is a
strict per-step check.

One thing to watch: 0.9988 on the domain denominator was calibrated for hurs.
For the elevation family the GMTED covariate covers the whole grid, so the
domain is almost the same as the grid, and 99.88 % is a strict bar. If a variable
cannot fill 99.88 % of its land grid, lower `min_coverage` in that variable's
config. It is a one-line change.

## Method

The gated families use an adaptive four-strategy cascade: autoKrige, explicit
variogram models, autoKrige with a boosted `nmax`, and autoKrige with a reduced
`nmax` and a higher nugget. The template tries each in turn until one meets the
gate, and keeps the best field it produced.

hurs is the multi-covariate case. It krigs on three covariates at once with the
equation `Data ~ DEM + layer.1 + layer.2`, and the month of each timestep
selects the month-specific tas and vapor-pressure covariate. It is separate from
the adaptive single-covariate pipeline and the two should not be mixed. The hurs
nugget-floor cascade (Gaussian kept, a 2 %/5 % floor, lm as a last resort) was
not ported to the adaptive template, because it fixes a Gaussian zero-nugget
failure specific to the multi-covariate kriging.

## Scenario name matching

The output scenario token is matched against both spellings whenever a driver
looks for a prior file or deletes one: historical or hist, and ssp370 or ssp3.
A re-run therefore finds an earlier file however it was written, and a gate
failure clears the companion set under either spelling. New files are still
written with the short token. The check scripts behave the same way, and map a
truncated scenario back to the full name for the reprocessing list (ssp3 back to
ssp370, hist back to historical).

The drivers also route intermediate files to node-local `$TMPDIR` to cut the
metadata load on Lustre, apply a per-iteration wall-clock timeout to the R call
(`r_timeout_seconds`, default 3600), log the R return code and stderr on
failure, and copy a failing R script to `$SCRATCH` for inspection.

## configs/

The `configs/` directory holds one gated JSON config per gated-family variable,
seven in all:

- Elevation family: `cmip6_config_tas_gated.json`,
  `cmip6_config_tasmax_gated.json`, `cmip6_config_tasmin_gated.json`,
  `cmip6_config_sfcWind_gated.json`.
- VP family: `cmip6_config_huss_gated.json`, `cmip6_config_rsds_gated.json`,
  `cmip6_config_rlds_gated.json`.

Each config names its family's R template, its covariate, and the gate keys
(`min_coverage`, `coverage_denominator`), along with the working and source
directories, the Apptainer image, the domain extent, and the two target
resolutions. The elevation configs point at the GMTED covariate; the huss config
points at the vapor-pressure covariate; the rsds and rlds configs point at the
srad covariate.

The srad path in the rsds and rlds configs is the same as the vapor-pressure
path, in the same source directory and with the same naming:
`/PATH/TO/SOURCE_DATA/srad_daily_1980_2014_%d.nc` alongside
`/PATH/TO/SOURCE_DATA/vp_daily_1980_2014_%d.nc`.

There is no hurs config in this directory. `setup_hurs_multi_tasks.sh` writes
`cmip6_config_hurs_multi.json` at run time. One config still needs a check: the
sfcWind config was reconstructed (its config was missing on the HPC), so confirm
the sfcWind source data exists before a full run.

## utils/

The `utils/` directory holds two reprocessing helpers and its own README.
`check_missing_1km_files.sh` compares a master list of expected files against an
output folder and writes the entries that are still missing, treating an entry
as present only when both its variable file and its SE file are there.
`check_1km_coverage.sh` scans a directory of finished 1 km files, flags any whose
valid-pixel coverage falls below a threshold (grid denominator, 89 % by default),
and can delete the flagged timestep's companion set (`-D`) so it gets re-run.
Both write the four-column list (MODEL, VARIABLE, SCENARIO, DATE) that the setup
scripts read back, so the first finds files that were never produced and the
second finds files that came out too sparse. See `utils/README.md` for the
details.

## Environment

The kriging runs in R. On TACC it runs inside the rstudio-rocker Apptainer
image, loaded with `module load tacc-apptainer`, and the configs point at that
image through `/PATH/TO/rstudio-rocker.sif`. Do not load the CDO or NCO modules
at the same time as the container, since that mix breaks `ncdf4` loading; load
only `tacc-apptainer` for R.

If you run the pipeline anywhere that does not use that Apptainer container, you
need a local installation of R with the kriging packages: KrigR, sp, raster,
automap, and gstat. The surrounding steps also use CDO and NCO for file
operations, and the TACC launcher for the parallel task runs.

The working and source locations are written as placeholders, so set them for
your site before running:

- `/PATH/TO/WORKING_DIR`: holds the scripts, the GMTED `.RData` files, the
  generated configs and task files, and the `downscaled/` outputs.
- `/PATH/TO/SOURCE_DATA`: holds the CMIP6 inputs, the covariates, the shapefiles,
  and the `tas/archive` that hurs uses.
- `/PATH/TO/rstudio-rocker.sif`: the Apptainer/Singularity R image.

## Covariate paths

- huss: `/PATH/TO/SOURCE_DATA/vp_daily_1980_2014_%d.nc` (vapor pressure).
- rsds and rlds: `/PATH/TO/SOURCE_DATA/srad_daily_1980_2014_%d.nc` (srad). This
  is the same path as the vapor-pressure covariate, in the same source directory
  and with the same naming.
- tas, tasmin, tasmax, sfcWind: GMTED elevation (`gmted_km.RData` and the iter
  files).
- hurs: the GMTED DEM, the model tas archive under
  `/PATH/TO/SOURCE_DATA/tas/archive/`, and the monthly vapor-pressure file.

## Running the pipelines

Deploy the drivers, R templates, builders, setups, the configs you need, the
GMTED `.RData` files, and the missing/extraction lists to `/PATH/TO/WORKING_DIR/`.

Elevation family:

    cd /PATH/TO/WORKING_DIR
    module load launcher tacc-apptainer cdo nco ncl_ncarg
    ./setup_elev_gated_tasks.sh tas               # or tasmin/tasmax/sfcWind
    export LAUNCHER_WORKDIR=$PWD
    export LAUNCHER_SCHED=interleaved
    export LAUNCHER_JOB_FILE=tas_gated_tasks.job
    $LAUNCHER_DIR/paramrun

VP family:

    ./setup_vp_gated_tasks.sh                      # builds rlds, rsds, huss
    export LAUNCHER_JOB_FILE=huss_gated_tasks.job  # or rlds_/rsds_
    $LAUNCHER_DIR/paramrun

hurs:

    ./setup_hurs_multi_tasks.sh                    # uses a missing-file list
    export LAUNCHER_JOB_FILE=hurs_reprocessing_tasks.job
    $LAUNCHER_DIR/paramrun

A single-task test is the fastest way to confirm the gate before a full submit:

    python3 custom_covariate_downscaling_elev_gated.py temp_extractions/tas/config_000001.json
    # the log shows "Iteration 1 PASSED coverage gate" / "Iteration 2 PASSED coverage gate"

The extraction and missing-file lists use four columns, MODEL, VARIABLE,
SCENARIO, DATE, with the date written as YYYY-MM-DD. The four-column format
matters: a two-column list silently processes nothing.

## Before you run

- Run the single-task test on the HPC before a full submit.
- Confirm one file's coverage with `cdo info` before any bulk delete with the
  scanner's `-D` option.
- Confirm the sfcWind source data exists; its config is reconstructed.
- Decide the gate value for the elevation family (0.9988 on the domain
  denominator is strict when the domain is close to the grid).
