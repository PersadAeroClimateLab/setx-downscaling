# setx-downscaling

SETx Downscaling Pipeline Source.

This repository holds the SETx-UIFL Climate Theme pipeline. The pipeline
downscales and bias-corrects NEX-GDDP-CMIP6 climate model output to 1 km daily
fields over Southeast Texas. The 1 km fields feed the hydrological and
air-quality models that study how projected climate shifts affect the region
and beyond. For project background, see [setx-uifl.org](https://setx-uifl.org).

## What the pipeline does

The downscaling method is universal kriging through KrigR in R. Each field is
downscaled in two steps: 28 km to 9 km, then 9 km to 1 km. Every field carries
a standard-error (SE) field beside it. The raw output is one NetCDF file per
timestep, per variable, per model, per scenario. The later stages consolidate
these files into Zarr, check them, and repackage them for analysis and
distribution.

The work has four stages. Each stage is a top-level directory with its own
README.

| Stage | Directory | What it does |
|-------|-----------|--------------|
| 1. Downscaling | [`downscaling/`](downscaling/) | Krig NEX-GDDP-CMIP6 to 1 km intermediate daily NetCDF, with reprocessing helpers. |
| 2. Post-processing | [`post-processing/`](post-processing/) | Consolidate the per-timestep NetCDF files into Zarr, fix bad timesteps, and CMORize for distribution. |
| 3. Quality control | [`quality-control/`](quality-control/) | Audit the downscaled Zarr stores, and flag bad dates for a clip-or-redownscale decision. |
| 4. Evaluation | [`evaluation/`](evaluation/) | Measure downscaling skill per variable against reference datasets, and make the manuscript figures and tables. |

## Pipeline flow

```
        NEX-GDDP-CMIP6 (28 km daily)
                  |
                  v
        [1] downscaling/            universal kriging (KrigR), 28 km -> 9 km -> 1 km
                  |                 1 km intermediate daily NetCDF, one file per timestep (+ SE)
                  v
        [2] post-processing/        consolidate ~15 million NetCDF files into Zarr,
                  |                 scan and swap bad timesteps, CMORize to CMIP6 NetCDF
                  v
        [3] quality-control/        audit the Zarr stores (automated + full scan),
                  |                 write a reviewed list of bad dates
                  |  bad dates -> clip or re-downscale ---------> back to [1] / [2]
                  v
        [4] evaluation/             skill vs AORC, Daymet, ERA5, ERA5-DS, GHCN-D;
                                    figures and metric tables
```

Stages 2 and 3 both find defective timesteps. A bad timestep goes back to
downscaling for a re-run, and the recomputed slice is swapped into the Zarr
store. Evaluation reads the finished stores.

## Scope

**Variables (8).** The pipeline downscales and evaluates eight variables:

- Temperature: `tas`, `tasmax`, `tasmin`
- Wind: `sfcWind`
- Humidity: `huss` (specific), `hurs` (relative)
- Radiation: `rlds` (longwave), `rsds` (shortwave)

**Models (10).** The NEX-GDDP-CMIP6 subset used here is BCC-CSM2-MR, CESM2,
CMCC-ESM2, CNRM-ESM2-1, EC-Earth3, FGOALS-g3, GFDL-CM4, MPI-ESM1-2-HR,
MRI-ESM2-0, and NorESM2-MM. The scripts discover the models present at run
time, so they do not hold a fixed model list.

**Scenarios (5).** `historical`, `ssp126`, `ssp245`, `ssp370`, and `ssp585`.
Some files use short scenario identifiers (`hist`, `ssp3`). The scripts match
both spellings, and write new files with the full name.

**Domain.** Southeast Texas, about 28.5 to 33.8 degrees N and 93.5 to 97.5
degrees W.

**Period.** The consolidated stores hold 1950 to 2014 for the historical
scenario and 2015 to 2100 for the SSP scenarios. The evaluation focuses on the
1980 to 2014 historical window.

## The four stages

### 1. Downscaling (`downscaling/`)

Three pipelines share this directory:

- Elevation family (`tas`, `tasmin`, `tasmax`, `sfcWind`): one covariate, the
  GMTED elevation DEM.
- VP family (`huss`, `rsds`, `rlds`): one covariate. `huss` uses vapor
  pressure; `rsds` and `rlds` use solar radiation (srad).
- hurs (multi-covariate): three covariates, the DEM, the model `tas` field,
  and a monthly vapor-pressure climatology (`Data ~ DEM + layer.1 + layer.2`).

A coverage gate runs inside each driver. After each kriging step the driver
measures the valid-pixel coverage. If a step falls below the threshold, the
driver deletes the companion files for that timestep and marks it failed. The
default gate is 99.88 % on the domain denominator, set in the config.

See [`downscaling/README.md`](downscaling/README.md) for the file roles, the
method cascade, the configs, and the run commands.

### 2. Post-processing (`post-processing/`)

The raw output is roughly 15 million individual NetCDF files. This stage makes
the data usable in five steps: build empty template Zarr stores, fill them from
the raw sources, scan for corrupted or missing timesteps, swap in recomputed
timesteps, and CMORize the finished stores into CMIP6-style (CMOR/DRS) NetCDF.

See [`post-processing/README.md`](post-processing/README.md) for the function
reference and the CMORization detail.

### 3. Quality control (`quality-control/`)

Two independent approaches audit the Zarr stores:

- A packaged, automated checker, carried as a git submodule
  (`SETx-UIFL-Downscaled-Quality-Control`). It runs a fixed set of checks per
  store and writes a pass/fail summary.
- A `full-scan-pipeline` of standalone scripts, run in sequence. It scans every
  store far more exhaustively and ends in a reviewed list of bad dates plus
  per-date plots, for a person to decide clip against re-downscale.

The submodule is not populated by a normal clone. See [Getting started](#getting-started)
and [`quality-control/README.md`](quality-control/README.md).

### 4. Evaluation (`evaluation/`)

One self-contained script per variable measures how well the 1 km downscaling
reproduces observed climate. Each script compares the downscaled field against
gridded references (AORC, Daymet, ERA5, and downscaled ERA5) and, where the
variable has a direct station measurement, against GHCN-D station data. The
scripts produce the manuscript figures (paired PDF and PNG at 1200 DPI) and the
metric tables.

See [`evaluation/README.md`](evaluation/README.md) for the per-variable script
table and the run commands.

## Repository layout

```
setx-downscaling/
├── downscaling/              Stage 1: kriging to 1 km daily NetCDF
│   ├── configs/              7 per-variable gated JSON configs
│   ├── utils/                missing-file and coverage helpers (+ README)
│   ├── *.py  *.R  *.sh       drivers, kriging templates, jobfile builders, setups
│   └── README.md
├── post-processing/          Stage 2: consolidate to Zarr, fix, CMORize
│   ├── data_cleaning_funcs.py
│   ├── cmorization.py
│   └── README.md
├── quality-control/          Stage 3: audit the Zarr stores
│   ├── SETx-UIFL-Downscaled-Quality-Control/   git submodule (automated checker)
│   ├── full-scan-pipeline/                     standalone audit chain
│   └── README.md
├── evaluation/               Stage 4: skill vs references, manuscript figures
│   ├── *_publication.py      one script per variable (8)
│   └── README.md
├── LICENSE.md                BSD 3-Clause License (source code)
├── .gitmodules
└── .gitignore
```

## Getting started

Clone the repository with its submodule:

```bash
git clone --recurse-submodules https://github.com/PersadAeroClimateLab/setx-downscaling.git
```

If the repository is already cloned, populate the quality-control submodule:

```bash
git submodule update --init
```

Each stage sets its own paths and runs on its own. Read the stage README before
you run it.

### Environments

- **Downscaling (R).** The kriging runs in R. On TACC it runs inside the
  rstudio-rocker Apptainer image, with the TACC launcher for the parallel task
  runs, and CDO and NCO for file operations. If you do not run it at TACC as an
  Apptainer, you need a local R install with KrigR, sp, raster, automap, and
  gstat.
- **Post-processing, quality control, evaluation (Python).** These stages use
  Python 3.10 or newer with xarray, zarr, dask, numpy, pandas, netCDF4, scipy,
  and matplotlib. Evaluation also uses cartopy and rioxarray as optional extras.
  The quality-control submodule ships its own Dockerfile and requirements.

### Paths

The scripts use descriptive placeholders in place of site paths, such as
`/PATH/TO/WORKING_DIR`, `/PATH/TO/SOURCE_DATA`, `/PATH/TO/ZARR_ROOT`, and
`/PATH/TO/AORC`. Set them for your machine in each script or config before you
run. The stage READMEs list the placeholders each stage needs.

## Data products

- **1 km daily NetCDF** (stage 1): one file per timestep, per variable, per
  model, per scenario, each with an SE partner.
- **Consolidated Zarr stores** (stage 2): one store per model, scenario, and
  time range, chunked for analysis.
- **CMORized NetCDF** (stage 2): CMIP6-style (CMOR/DRS) output, split by year,
  for distribution.

## Reference datasets (evaluation)

Evaluation compares the downscaled fields against AORC, Daymet, ERA5, the
downscaled ERA5, and GHCN-D station observations.

## Script catalogue

[`SCRIPT_CATALOG.md`](SCRIPT_CATALOG.md) lists every script across the four
stages with its role.

## Contributing and feedback

Please report any bugs, ask questions, and make suggestions through the GitHub Issues form of this repository.

## Acknowledgements

This work was supported by the U.S. Department of Energy, Office of Science, Biological and 
Environmental Research Program's Southeast Texas Urban Integrated Field Laboratory under Award Number DE-SC0023216.

We thank the following people for their contributions to this project:

1. Dr. Geeta Persad, for her design input, guidance, and mentorship throughout the project.
2. The SETx-UIFL Climate Theme, for their input and expertise throughout the process.
3. The members of the Persad Aero-Climate Lab, for being a sounding board through the design, implementation, and production stages.

## License and citation

The source code is released under the BSD 3-Clause License (see [`LICENSE.md`](LICENSE.md)).
Downscaled data products are deposited separately. For citation, see
[setx-uifl.org](https://setx-uifl.org), and DOI:XXXXX.
