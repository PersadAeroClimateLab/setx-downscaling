#!/bin/bash
#===============================================================================
# setup_vp_gated_tasks.sh - Set up task files for the VP family
# Variables: rlds, rsds, huss  (covariate: huss = vapor pressure, rsds/rlds = srad)
#
# Writes one config per variable (kriging_template_adaptive_vp_gated.R +
# coverage-gate keys), runs the jobfile builder
# (create_launcher_jobfile_vp_gated.py), and emits <var>_gated_tasks.job that
# call custom_covariate_downscaling_vp_gated.py.
#===============================================================================

WORKING_DIR="/PATH/TO/WORKING_DIR"
cd "$WORKING_DIR" || { echo "ERROR: cannot cd to $WORKING_DIR"; exit 1; }

JOBFILE_BUILDER="create_launcher_jobfile_vp_gated.py"
DRIVER="custom_covariate_downscaling_vp_gated.py"
R_TEMPLATE="kriging_template_adaptive_vp_gated.R"

# Covariate per variable within the VP family:
#   huss       -> vapor pressure (vp)
#   rsds, rlds -> solar radiation (srad)
# Same naming/structure as the vp covariate file; adjust the path if yours differs.
VP_COV="/PATH/TO/SOURCE_DATA/vp_daily_1980_2014_%d.nc"
SRAD_COV="/PATH/TO/SOURCE_DATA/srad_daily_1980_2014_%d.nc"
covariate_for() {
    case "$1" in
        huss)      echo "$VP_COV" ;;
        rsds|rlds) echo "$SRAD_COV" ;;
        *)         echo "$VP_COV" ;;
    esac
}

echo "=============================================="
echo "Setting up VP-family GATED task files (rlds, rsds, huss)"
echo "  R template: $R_TEMPLATE"
echo "  Driver:     $DRIVER"
echo "  Coverage gate: 0.9988 on the domain denominator (per config)"
echo "=============================================="
echo ""

# Check required pipeline files
for f in "$JOBFILE_BUILDER" "$DRIVER" "$R_TEMPLATE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file not found: $f"
        echo "       Copy it to $WORKING_DIR first."
        exit 1
    fi
done

# GMTED base files (needed even though the covariate here is VP: the two-step
# grids come from GMTED).
for g in gmted_28km_to_9km.RData gmted_9km_to_1km.RData; do
    if [[ ! -f "$g" ]]; then
        echo "ERROR: GMTED file not found: $g"
        exit 1
    fi
done

# Clean old GATED task files (leave any non-gated production files untouched)
echo "Cleaning old gated task configs..."
for VAR in rlds rsds huss; do
    rm -f "temp_extractions/${VAR}/spec_"*.txt "temp_extractions/${VAR}/config_"*.json 2>/dev/null
    mkdir -p "temp_extractions/${VAR}" "downscaled/${VAR}"
done
echo ""

# Emit one gated config per variable, then build its gated jobfile.
make_config() {
    local VAR="$1"
    local CONFIG_FILE="cmip6_config_${VAR}_gated.json"
    local COV
    COV="$(covariate_for "$VAR")"
    echo "Creating config: $CONFIG_FILE   (covariate: $COV)"
    cat > "$CONFIG_FILE" << EOF
{
    "working_dir": "${WORKING_DIR}",
    "input_dir": "/PATH/TO/SOURCE_DATA",
    "output_dir": "${WORKING_DIR}/downscaled/${VAR}",
    "extraction_list": "missing_files_${VAR}_v4.txt",
    "variable": "${VAR}",
    "temp_dir": "${WORKING_DIR}/temp_extractions/${VAR}",
    "r_script_template": "${R_TEMPLATE}",
    "apptainer_image": "/PATH/TO/rstudio-rocker.sif",
    "required_modules": ["tacc-apptainer", "cdo", "launcher", "nco", "ncl_ncarg"],
    "shapefile_path": "/PATH/TO/SOURCE_DATA/huss/ne_10m_admin_1_states_provinces.shp",
    "extent_lonmin": -97.5,
    "extent_lonmax": -93.5,
    "extent_latmin": 28.5,
    "extent_latmax": 33.8,
    "gmted_file": "${WORKING_DIR}/gmted_km.RData",
    "gmted_iter1": "${WORKING_DIR}/gmted_28km_to_9km.RData",
    "gmted_iter2": "${WORKING_DIR}/gmted_9km_to_1km.RData",
    "covariate_base_path": "${COV}",
    "kriging_equation": "Data ~ layer",
    "target_res_first": 0.08333,
    "target_res_second": 0.008334,
    "nmax_base": 40,
    "cores": 30,
    "min_coverage": 0.9988,
    "coverage_denominator": "domain",
    "r_timeout_seconds": 3600,
    "comment": "VP-family GATED downscaling for ${VAR} (covariate: ${COV}). huss uses vapor pressure; rsds/rlds use srad. After EACH downscale the coverage must be >= min_coverage on the coverage_denominator, else both 9km and 1km (+SE) are deleted and the timestep fails."
}
EOF
    python3 "$JOBFILE_BUILDER" "$CONFIG_FILE" "${VAR}_gated_tasks.job"
    echo ""
}

for VAR in rlds rsds huss; do
    echo "=============================================="
    echo "Creating ${VAR} gated tasks..."
    echo "=============================================="
    make_config "$VAR"
done

# Summary
echo "=============================================="
echo "TASK FILE SUMMARY"
echo "=============================================="
for VAR in rlds rsds huss; do
    JOB="${VAR}_gated_tasks.job"
    if [[ -f "$JOB" ]]; then
        TASKS=$(grep -c "python3 ${DRIVER}" "$JOB" 2>/dev/null || echo "0")
        echo "  $JOB : $TASKS tasks (temp_extractions/${VAR})"
    fi
done
echo ""

echo "=============================================="
echo "USAGE"
echo "=============================================="
echo "Run in idev session (one variable at a time):"
echo "  cd $WORKING_DIR"
echo "  module load launcher tacc-apptainer cdo nco ncl_ncarg"
echo "  export LAUNCHER_WORKDIR=\$PWD"
echo "  export LAUNCHER_SCHED=interleaved"
echo "  export LAUNCHER_JOB_FILE=rlds_gated_tasks.job   # or rsds_/huss_"
echo "  \$LAUNCHER_DIR/paramrun"
echo ""
