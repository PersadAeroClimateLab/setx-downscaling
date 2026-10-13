#!/bin/bash
#===============================================================================
# setup_hurs_multi_tasks.sh - Setup task files for MULTI-COVARIATE hurs
#
# hurs (relative humidity) is downscaled with THREE covariates:
#   - DEM (GMTED elevation)
#   - Temperature (pre-downscaled tas, MODEL- and SCENARIO-specific)
#   - Vapor pressure (monthly climatology)
# Kriging equation: Data ~ DEM + layer.1 + layer.2
#
# This is the multi-covariate pipeline (NOT the adaptive single-covariate one).
# It uses these files:
#   - create_launcher_jobfile_adaptive_fixed.py  (builds the launcher jobfile)
#   - kriging_template_multi_covariate.R          (multi-covariate kriging template)
#   - custom_covariate_hurs_downscaling.py             (per-task driver)
#
# It produces:
#   - cmip6_config_hurs_multi.json   (pipeline config)
#   - hurs_reprocessing_tasks.job    (launcher jobfile - run on idev OR via slurm)
#   - submit_hurs_launcher.slurm     (slurm batch submit script)
#
# It also PRECOMPUTES the tas monthly-mean covariate with cdo timmean (see the
# "tas monthly means" step). The mean files are written next to the source tas
# files and named tasmean_<same-name>. This removes the large per-timestep read
# and the long R runtime, because R then loads one field instead of a daily stack.
#
# Usage: ./setup_hurs_multi_tasks.sh [missing_file]
#   missing_file : date-based extraction list (default missing_files_hurs_v4.txt)
#                  Format: MODEL, VARIABLE, SCENARIO, DATE
#===============================================================================

WORKING_DIR="/PATH/TO/WORKING_DIR"
VAR="hurs"

# Pipeline files. Change JOBFILE_BUILDER if your builder has another name; it
# must accept "<config.json> <output.job>" and emit one
# "python3 custom_covariate_hurs_downscaling.py <config>" line per task.
JOBFILE_BUILDER="create_launcher_jobfile_adaptive_fixed.py"
R_TEMPLATE="kriging_template_multi_covariate.R"
DRIVER="custom_covariate_hurs_downscaling.py"

# Covariate locations
TAS_ARCHIVE_DIR="/PATH/TO/SOURCE_DATA/tas/archive"
# tas filename inside the config path (Python fills {scenario}/{model}; R fills %d
# and resolves the r* realization by glob). Keep in sync with the check below.
VP_PATTERN="/PATH/TO/SOURCE_DATA/vp_daily_1980_2014_%d.nc"

cd "$WORKING_DIR" || { echo "ERROR: cannot cd to $WORKING_DIR"; exit 1; }

MISSING_FILE="${1:-missing_files_${VAR}_v4.txt}"

echo "=============================================="
echo "Setting up $VAR task files (Multi-Covariate Kriging)"
echo "  Covariates: DEM + temperature (tas monthly mean) + vapor pressure"
echo "  Equation:   Data ~ DEM + layer.1 + layer.2"
echo "=============================================="
echo ""

#-------------------------------------------------------------------------------
# Check required pipeline files
#-------------------------------------------------------------------------------
for f in "$JOBFILE_BUILDER" "$DRIVER" "$R_TEMPLATE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file not found: $f"
        echo "       Copy it to $WORKING_DIR first."
        exit 1
    fi
done

if [[ ! -f "$MISSING_FILE" ]]; then
    echo "ERROR: extraction/missing file not found: $MISSING_FILE"
    echo "       Create it first (e.g. with utils/check_missing_1km_files.sh)."
    exit 1
fi

# GMTED covariate base files (required - kriging cannot run without them)
for g in gmted_28km_to_9km.RData gmted_9km_to_1km.RData; do
    if [[ ! -f "$g" ]]; then
        echo "ERROR: GMTED file not found: $g"
        exit 1
    fi
done
echo "Pipeline files and GMTED base found."
echo ""

#-------------------------------------------------------------------------------
# Count entries
#-------------------------------------------------------------------------------
ENTRY_COUNT=$(grep -cv "^#" "$MISSING_FILE" 2>/dev/null || echo "0")
echo "Found $ENTRY_COUNT entries in $MISSING_FILE"
echo ""

#-------------------------------------------------------------------------------
# Check covariate prerequisites
# Multi-covariate hurs FAILS on a timestep if the matching tas file is absent.
# tas is MODEL- and SCENARIO-specific, so check each pair in the list.
#-------------------------------------------------------------------------------
echo "Checking covariate prerequisites..."
MISSING_COV=0

# Unique MODEL + SCENARIO pairs in the extraction list
PAIRS=$(grep -v "^#" "$MISSING_FILE" | \
    awk -F',' 'NF>=4 {gsub(/ /,"",$1); gsub(/ /,"",$3); print $1","$3}' | sort -u)

if [[ -z "$PAIRS" ]]; then
    echo "  WARNING: no MODEL,SCENARIO pairs parsed from $MISSING_FILE"
else
    while IFS=',' read -r MODEL SCENARIO; do
        [[ -z "$MODEL" ]] && continue
        # tas filenames use the full CMIP6 scenario name (e.g. "historical").
        # The extraction list may use the short form ("hist"). Search both.
        case "$SCENARIO" in
            hist|historical) SCEN_CANDIDATES="historical hist" ;;
            *)               SCEN_CANDIDATES="$SCENARIO" ;;
        esac
        for m in $(seq 1 12); do
            found=0
            for SCEN in $SCEN_CANDIDATES; do
                GLOB="${TAS_ARCHIVE_DIR}/tas_NAM_SETx-UTAustin_KrigR_CMIP6_${SCEN}_${MODEL}_r*_day_v1_${m}.nc"
                if ls $GLOB >/dev/null 2>&1; then found=1; break; fi
            done
            if [[ "$found" -eq 0 ]]; then
                echo "  Missing tas: ${SCENARIO} ${MODEL} month ${m}"
                MISSING_COV=$((MISSING_COV + 1))
            fi
        done
    done <<< "$PAIRS"
fi

# Vapor pressure (same for all models)
for m in $(seq 1 12); do
    VP_FILE=$(printf "$VP_PATTERN" "$m")
    if [[ ! -f "$VP_FILE" ]]; then
        echo "  Missing vapor pressure: month ${m} ($VP_FILE)"
        MISSING_COV=$((MISSING_COV + 1))
    fi
done

if [[ "$MISSING_COV" -gt 0 ]]; then
    echo ""
    echo "  WARNING: $MISSING_COV covariate file(s) missing (listed above)."
    echo "  Those timesteps will fail until the covariates exist."
    echo "  tas must be downscaled per model+scenario before hurs can use it."
    echo "  Task files are still created so you can run the rest."
    echo ""
else
    echo "  All required tas and vapor-pressure covariates found."
    echo ""
fi

#-------------------------------------------------------------------------------
# Precompute tas monthly means (cdo timmean)
#
# hurs uses the temperature covariate as a monthly-mean field. Collapsing the
# daily tas stack to one field with cdo removes the large per-timestep read and
# the long R runtime. Each mean is written NEXT TO its source tas file and named
# tasmean_<source-name>. Existing (non-empty) means are skipped, so re-running
# this setup is cheap.
#
# NOTE: for a very large missing file this loops over many model+scenario+month
# files. cdo timmean is I/O-bound and low-memory, but if it is slow on the login
# node, run this setup inside an idev session.
#-------------------------------------------------------------------------------
echo "Precomputing tas monthly means (cdo timmean)..."

if ! command -v cdo >/dev/null 2>&1; then
    module load cdo 2>/dev/null || true
fi
if ! command -v cdo >/dev/null 2>&1; then
    echo "ERROR: cdo not found. Run 'module load cdo' then re-run this script."
    exit 1
fi

MEAN_MADE=0; MEAN_SKIP=0; MEAN_FAIL=0
if [[ -n "$PAIRS" ]]; then
    while IFS=',' read -r MODEL SCENARIO; do
        [[ -z "$MODEL" ]] && continue
        case "$SCENARIO" in
            hist|historical) SCEN_CANDIDATES="historical hist" ;;
            *)               SCEN_CANDIDATES="$SCENARIO" ;;
        esac
        for m in $(seq 1 12); do
            SRC=""
            for SCEN in $SCEN_CANDIDATES; do
                GLOB="${TAS_ARCHIVE_DIR}/tas_NAM_SETx-UTAustin_KrigR_CMIP6_${SCEN}_${MODEL}_r*_day_v1_${m}.nc"
                CAND=$(ls $GLOB 2>/dev/null | head -1)
                if [[ -n "$CAND" ]]; then SRC="$CAND"; break; fi
            done
            # No source tas for this model+scenario+month: already warned above.
            [[ -z "$SRC" ]] && continue

            BASE=$(basename "$SRC")
            MEAN="${TAS_ARCHIVE_DIR}/tasmean_${BASE#tas_}"

            if [[ -s "$MEAN" ]]; then
                MEAN_SKIP=$((MEAN_SKIP + 1))
                continue
            fi

            if cdo -s timmean "$SRC" "$MEAN" >/dev/null 2>&1; then
                MEAN_MADE=$((MEAN_MADE + 1))
            else
                echo "  WARNING: cdo timmean failed for $BASE"
                rm -f "$MEAN" 2>/dev/null
                MEAN_FAIL=$((MEAN_FAIL + 1))
            fi
        done
    done <<< "$PAIRS"
fi
echo "  tas means: made $MEAN_MADE, skipped(existing) $MEAN_SKIP, failed $MEAN_FAIL"
echo ""

#-------------------------------------------------------------------------------
# Clean old task files, create output/temp dirs
#-------------------------------------------------------------------------------
echo "Cleaning old task configs for $VAR..."
rm -f "temp_extractions/${VAR}/spec_"*.txt 2>/dev/null
rm -f "temp_extractions/${VAR}/config_"*.json 2>/dev/null
mkdir -p "temp_extractions/${VAR}"
mkdir -p "downscaled/${VAR}"
echo ""

#-------------------------------------------------------------------------------
# Create config
#
# The temperature covariate path points at the PRECOMPUTED monthly means
# (tasmean_...). Type stays "monthly_nc_stack"; the R template loads the single
# mean field (its mean() is a no-op on one layer).
#-------------------------------------------------------------------------------
CONFIG_FILE="cmip6_config_${VAR}_multi.json"
echo "Creating config: $CONFIG_FILE"

cat > "$CONFIG_FILE" << EOF
{
    "working_dir": "${WORKING_DIR}",
    "input_dir": "/PATH/TO/SOURCE_DATA",
    "output_dir": "${WORKING_DIR}/downscaled/${VAR}",
    "extraction_list": "${MISSING_FILE}",
    "variable": "${VAR}",
    "temp_dir": "${WORKING_DIR}/temp_extractions/${VAR}",
    "r_script_template": "${R_TEMPLATE}",
    "apptainer_image": "/PATH/TO/rstudio-rocker.sif",
    "required_modules": ["tacc-apptainer", "cdo", "launcher", "nco", "ncl_ncarg"],
    "shapefile_path": "/PATH/TO/SOURCE_DATA/${VAR}/ne_10m_admin_1_states_provinces.shp",
    "extent_lonmin": -97.5,
    "extent_lonmax": -93.5,
    "extent_latmin": 28.5,
    "extent_latmax": 33.8,
    "gmted_file": "${WORKING_DIR}/gmted_km.RData",
    "gmted_iter1": "${WORKING_DIR}/gmted_28km_to_9km.RData",
    "gmted_iter2": "${WORKING_DIR}/gmted_9km_to_1km.RData",
    "covariates": [
        {
            "name": "temperature",
            "type": "monthly_nc_stack",
            "path": "${TAS_ARCHIVE_DIR}/tasmean_NAM_SETx-UTAustin_KrigR_CMIP6_{scenario}_{model}_r*_day_v1_%d.nc",
            "description": "Precomputed monthly-MEAN temperature (cdo timmean) - MODEL- and SCENARIO-specific. r* matches the realization member.",
            "model_specific": true
        },
        {
            "name": "vapor_pressure",
            "type": "monthly_file",
            "path": "${VP_PATTERN}",
            "description": "Monthly vapor pressure climatology - same for all models",
            "model_specific": false
        }
    ],
    "kriging_equation": "Data ~ DEM + layer.1 + layer.2",
    "target_res_first": 0.08333,
    "target_res_second": 0.008334,
    "nmax": 75,
    "cores": 30,
    "min_coverage": 0.9988,
    "coverage_denominator": "domain",
    "r_timeout_seconds": 3600,
    "valid_min": 0,
    "valid_max": 100,
    "comment": "Multi-covariate two-step downscaling for hurs (DEM + tas monthly mean + vapor pressure). tas covariate must match model AND scenario. Coverage gate: each downscale must cover >= min_coverage of the area (NaN <= 1 - min_coverage), else both 9km and 1km (+SE) are deleted and the timestep fails. coverage_denominator: 'domain' = cells with complete covariates (recommended); 'grid' = whole bounding box."
}
EOF
echo ""

#-------------------------------------------------------------------------------
# Create launcher jobfile (one task per timestep)
#-------------------------------------------------------------------------------
JOB_FILE="${VAR}_reprocessing_tasks.job"
echo "Creating jobfile: $JOB_FILE"
python3 "$JOBFILE_BUILDER" "$CONFIG_FILE" "$JOB_FILE"

# The jobfile builder may still emit the old driver name
# (custom_covariate_downscaling.py). Point every task at this pipeline's driver.
if [[ -f "$JOB_FILE" ]]; then
    sed -i "s/custom_covariate_downscaling\.py/${DRIVER}/g" "$JOB_FILE"
fi
echo ""

#-------------------------------------------------------------------------------
# Create slurm submit script (batch alternative to idev)
#-------------------------------------------------------------------------------
SLURM_FILE="submit_${VAR}_launcher.slurm"
echo "Creating slurm submit script: $SLURM_FILE"

cat > "$SLURM_FILE" << EOF
#!/bin/bash
#SBATCH -J ${VAR}_downscale
#SBATCH -o ${VAR}_downscale.%j.out
#SBATCH -e ${VAR}_downscale.%j.err
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -A REPLACE_WITH_ALLOCATION
#
# EDIT BEFORE SUBMIT:
#   -A  : your TACC allocation
#   -N  : nodes,  -n : total parallel tasks,  -t : walltime
# Each task calls autoKrige with "cores" from the config (currently 30), so keep
# tasks-per-node low enough not to oversubscribe the node's cores.

module load launcher tacc-apptainer cdo nco ncl_ncarg

export LAUNCHER_WORKDIR=${WORKING_DIR}
export LAUNCHER_JOB_FILE=${WORKING_DIR}/${JOB_FILE}
export LAUNCHER_SCHED=interleaved   # workaround if the dynamic task service fails

\$LAUNCHER_DIR/paramrun
EOF
echo ""

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
echo "=============================================="
echo "SUMMARY"
echo "=============================================="
if [[ -f "$JOB_FILE" ]]; then
    TASKS=$(grep -c "python3 ${DRIVER}" "$JOB_FILE" 2>/dev/null || echo "0")
    echo "  Config file:  $CONFIG_FILE"
    echo "  Job file:     $JOB_FILE   ($TASKS tasks)"
    echo "  Slurm script: $SLURM_FILE"
    echo "  Temp dir:     temp_extractions/${VAR}"
    echo "  Output dir:   downscaled/${VAR}"
    echo "  tas means:    made $MEAN_MADE, skipped $MEAN_SKIP, failed $MEAN_FAIL"
    if [[ "$MISSING_COV" -gt 0 ]]; then
        echo "  Covariates:   $MISSING_COV missing (see warnings above)"
    fi
fi
echo ""

#-------------------------------------------------------------------------------
# Verification
#-------------------------------------------------------------------------------
echo "=============================================="
echo "VERIFICATION"
echo "=============================================="
SPEC_FILE="temp_extractions/${VAR}/spec_000001.txt"
if [[ -f "$SPEC_FILE" ]]; then
    echo "First spec file for $VAR:"
    cat "$SPEC_FILE"
    echo ""
fi

#-------------------------------------------------------------------------------
# Usage
#-------------------------------------------------------------------------------
echo "=============================================="
echo "USAGE"
echo "=============================================="
echo ""
echo "Option A - interactive (idev):"
echo "  idev -p normal -N 4 -n 8 -t 02:00:00"
echo "  cd $WORKING_DIR"
echo "  module load launcher tacc-apptainer cdo nco ncl_ncarg"
echo "  export LAUNCHER_WORKDIR=\$PWD"
echo "  export LAUNCHER_JOB_FILE=$JOB_FILE"
echo "  export LAUNCHER_SCHED=interleaved"
echo "  \$LAUNCHER_DIR/paramrun"
echo ""
echo "Option B - batch (slurm):"
echo "  # edit $SLURM_FILE first: set -A <allocation>, tune -N/-n/-t"
echo "  sbatch $SLURM_FILE"
echo ""
