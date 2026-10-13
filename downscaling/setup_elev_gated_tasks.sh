#!/bin/bash
#===============================================================================
# setup_elev_gated_tasks.sh - Set up task files for the elevation family
# Variables: tas, tasmin, tasmax, sfcWind  (covariate = GMTED elevation)
#
# Writes the config (kriging_template_adaptive_elev_gated.R + coverage-gate
# keys), runs the jobfile builder (create_launcher_jobfile_elev_gated.py), and
# emits <var>_gated_tasks.job that call custom_covariate_downscaling_elev_gated.py.
#
# Usage: ./setup_elev_gated_tasks.sh <variable> [missing_file]
#   ./setup_elev_gated_tasks.sh tas
#   ./setup_elev_gated_tasks.sh tasmin missing_files_tasmin_v4.txt
#   ./setup_elev_gated_tasks.sh sfcWind
#===============================================================================

WORKING_DIR="/PATH/TO/WORKING_DIR"
cd "$WORKING_DIR" || { echo "ERROR: cannot cd to $WORKING_DIR"; exit 1; }

JOBFILE_BUILDER="create_launcher_jobfile_elev_gated.py"
DRIVER="custom_covariate_downscaling_elev_gated.py"
R_TEMPLATE="kriging_template_adaptive_elev_gated.R"

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <variable> [missing_file]"
    echo "Supported variables: tas, tasmin, tasmax, sfcWind"
    exit 1
fi

VAR="$1"
MISSING_FILE="${2:-missing_files_${VAR}_v4.txt}"

case "$VAR" in
    tas|tasmin|tasmax|sfcWind)
        echo "Setting up GATED tasks for: $VAR"
        ;;
    *)
        echo "WARNING: Non-standard variable '$VAR'"
        echo "Standard variables: tas, tasmin, tasmax, sfcWind"
        read -p "Continue anyway? (y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
        ;;
esac

echo "=============================================="
echo "Setting up $VAR task files (Adaptive Kriging, GATED)"
echo "  R template: $R_TEMPLATE"
echo "  Driver:     $DRIVER"
echo "  Coverage gate: 0.9988 on the domain denominator (per config)"
echo "=============================================="
echo ""

# Check required files
for f in "$JOBFILE_BUILDER" "$DRIVER" "$R_TEMPLATE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file not found: $f"
        echo "       Copy it to $WORKING_DIR first."
        exit 1
    fi
done

for g in gmted_28km_to_9km.RData gmted_9km_to_1km.RData; do
    if [[ ! -f "$g" ]]; then
        echo "ERROR: GMTED file not found: $g"
        exit 1
    fi
done

if [[ ! -f "$MISSING_FILE" ]]; then
    echo "ERROR: Missing file not found: $MISSING_FILE"
    echo "       Create it first (e.g. with a coverage/missing-file check)."
    exit 1
fi

ENTRY_COUNT=$(grep -cv "^#" "$MISSING_FILE" 2>/dev/null || echo "0")
echo "Found $ENTRY_COUNT entries in $MISSING_FILE"
echo ""

# Clean old GATED task files for this variable
echo "Cleaning old gated task configs for $VAR..."
rm -f "temp_extractions/${VAR}/spec_"*.txt "temp_extractions/${VAR}/config_"*.json 2>/dev/null
mkdir -p "temp_extractions/${VAR}" "downscaled/${VAR}"
echo ""

# Create GATED config
CONFIG_FILE="cmip6_config_${VAR}_gated.json"
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
    "covariate_base_path": "${WORKING_DIR}/gmted_km.RData",
    "kriging_equation": "Data ~ layer",
    "target_res_first": 0.08333,
    "target_res_second": 0.008334,
    "nmax_base": 40,
    "cores": 30,
    "min_coverage": 0.9988,
    "coverage_denominator": "domain",
    "r_timeout_seconds": 3600,
    "comment": "Elevation-family GATED downscaling for ${VAR} (covariate = GMTED elevation). After EACH downscale the coverage must be >= min_coverage on the coverage_denominator, else both 9km and 1km (+SE) are deleted and the timestep fails. NOTE: for elevation vars domain ~= grid, so 0.9988 is strict; lower min_coverage here if needed.",
    "adaptive_parameters": {
        "description": "Automatically adjusts nugget and nmax based on correlation with elevation",
        "strong_correlation": {"threshold": ">=0.50", "nugget": 0.01, "nmax_multiplier": 1.0},
        "moderate_correlation": {"threshold": "0.30-0.49", "nugget": 0.02, "nmax_multiplier": 0.75},
        "weak_correlation": {"threshold": "<0.30", "nugget": 0.05, "nmax_multiplier": 0.50}
    }
}
EOF
echo ""

# Create GATED launcher jobfile
JOB_FILE="${VAR}_gated_tasks.job"
echo "Creating jobfile: $JOB_FILE"
python3 "$JOBFILE_BUILDER" "$CONFIG_FILE" "$JOB_FILE"

echo ""
echo "=============================================="
echo "SUMMARY"
echo "=============================================="
if [[ -f "$JOB_FILE" ]]; then
    TASKS=$(grep -c "python3 ${DRIVER}" "$JOB_FILE" 2>/dev/null || echo "0")
    echo "  Config file:  $CONFIG_FILE"
    echo "  Job file:     $JOB_FILE"
    echo "  Total tasks:  $TASKS"
    echo "  Temp dir:     temp_extractions/${VAR}"
    echo "  Output dir:   downscaled/${VAR}"
fi
echo ""

echo "=============================================="
echo "USAGE"
echo "=============================================="
echo "Run in idev session:"
echo "  cd $WORKING_DIR"
echo "  module load launcher tacc-apptainer cdo nco ncl_ncarg"
echo "  export LAUNCHER_WORKDIR=\$PWD"
echo "  export LAUNCHER_SCHED=interleaved"
echo "  export LAUNCHER_JOB_FILE=${JOB_FILE}"
echo "  \$LAUNCHER_DIR/paramrun"
echo ""
