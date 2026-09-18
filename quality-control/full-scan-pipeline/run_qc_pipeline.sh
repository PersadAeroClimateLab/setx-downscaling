#!/usr/bin/env bash
# run_qc_pipeline.sh -- run the entire SETx-1km-CMIP6 Zarr QC end to end.
#
#   qc_zarr_all.py  ->  write_issue_csvs.py  ->  drilldown_qc.py
#     ->  classify_footprint.py  ->  make_bad_dates.sh  ->  plot_bad_dates.py
#
# footprint runs BEFORE the final bad_dates assembly so its inland-cutoff dates
# are included; it still gets a preliminary bad-values list for oob-cell location.
#
# Usage (long scan -> run under nohup):
#   nohup ./run_qc_pipeline.sh > qc_pipeline.log 2>&1 &
#   tail -f qc_pipeline.log
#
# Override any path/knob via env, e.g.:
#   ROOT=/data/zarr WORKERS=12 FRESH=1 ./run_qc_pipeline.sh
set -euo pipefail

# ---------------- config (edit or override via env) ----------------
ROOT="${ROOT:-/local1/SETx_UIFL_Zarr}"
SFCWIND_ROOT="${SFCWIND_ROOT:-/local1/SFCWIND_ZARR_POST}"
QC_OUT="${QC_OUT:-qc_zarr_out}"
ISSUES_DIR="${ISSUES_DIR:-issues_by_variable}"
DRILL_OUT="${DRILL_OUT:-qc_drilldown_out}"
FOOTPRINT_OUT="${FOOTPRINT_OUT:-footprint_out}"
PLOTS_OUT="${PLOTS_OUT:-bad_date_plots}"
BAD_DATES="${BAD_DATES:-bad_dates.csv}"
WORKERS="${WORKERS:-8}"
PY="${PY:-python3}"
FRESH="${FRESH:-1}"          # 1 = wipe the scan summary for a clean re-scan
# stores to skip in the scan: the pre-re-downscale primary CESM2 sfcWind (the good
# copy comes from SFCWIND_ROOT) and the pr stores (different Zarr format -> errors).
SKIP="${SKIP:-sfcWind:CESM2:*:primary,pr}"
# inland-cutoff detection (classify_footprint): a day missing >= INLAND_MISS_MIN
# reference-inland cells is a bad date. MAX_INLAND_DATES=0 emits ALL of them per
# store (comprehensive review); set >0 to emit only a representative sample.
INLAND_MISS_MIN="${INLAND_MISS_MIN:-300}"
MAX_INLAND_DATES="${MAX_INLAND_DATES:-0}"

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
step() { echo; echo "======================================================================"; \
         echo "  [$(date '+%H:%M:%S')]  $*"; \
         echo "======================================================================"; }
need() { [ -f "$DIR/$1" ] || { echo "ERROR: missing $DIR/$1"; exit 1; }; }
trap 'echo; echo "!! pipeline FAILED at line $LINENO (see above)"; exit 1' ERR

for s in qc_zarr_all.py write_issue_csvs.py drilldown_qc.py classify_footprint.py \
         make_bad_dates.sh plot_bad_dates.py; do need "$s"; done
[ -d "$ROOT" ] || { echo "ERROR: ROOT not found: $ROOT"; exit 1; }

# derived outputs are cheap -> always regenerate; the scan resumes unless FRESH=1
rm -rf "$ISSUES_DIR" "$DRILL_OUT" "$FOOTPRINT_OUT" "$PLOTS_OUT" \
       "$BAD_DATES" prelim_bad_dates.csv
[ "$FRESH" = "1" ] && rm -rf "$QC_OUT"
echo "config: ROOT=$ROOT  SFCWIND_ROOT=$SFCWIND_ROOT  WORKERS=$WORKERS  FRESH=$FRESH"
echo "        SKIP=$SKIP"

# ---------------- 1. scan every Zarr store ----------------
step "1/6  qc_zarr_all.py  (Zarr scan: missing/partial/empty, out-of-range, blowups)"
"$PY" "$DIR/qc_zarr_all.py" --root "$ROOT" --sfcwind-root "$SFCWIND_ROOT" \
      --out "$QC_OUT" --workers "$WORKERS" --skip "$SKIP"

# ---------------- 2. per-variable issue CSVs ----------------
step "2/6  write_issue_csvs.py  (one <var>_issues.csv per variable)"
"$PY" "$DIR/write_issue_csvs.py" --qc-dir "$QC_OUT" --out "$ISSUES_DIR"

# ---------------- 3. drill-down summary ----------------
step "3/6  drilldown_qc.py  (bucket stores: real_gap / bad_values / wobble)"
"$PY" "$DIR/drilldown_qc.py" --qc-dir "$QC_OUT" --out "$DRILL_OUT"

# ---------------- 4. footprint classification ----------------
# preliminary bad-values list so classify_footprint can locate oob cells
{ echo "variable,model,scenario,date,issue,category"; \
  awk -F, 'FNR>1{sub(/\r$/,""); if($NF=="bad_values") print}' "$ISSUES_DIR"/*_issues.csv \
      2>/dev/null || true; } > prelim_bad_dates.csv
step "4/6  classify_footprint.py  (edge-wobble vs inland-cutoff + inland dates)"
"$PY" "$DIR/classify_footprint.py" --root "$ROOT" --sfcwind-root "$SFCWIND_ROOT" \
      --bad-dates prelim_bad_dates.csv --out "$FOOTPRINT_OUT" \
      --inland-miss-min "$INLAND_MISS_MIN" --max-inland-dates "$MAX_INLAND_DATES"

# ---------------- 5. assemble the FINAL bad_dates (last) ----------------
step "5/6  make_bad_dates.sh  (union of ALL checks incl. inland cutoffs)"
bash "$DIR/make_bad_dates.sh" "$ISSUES_DIR" "$FOOTPRINT_OUT/inland_cutoff_dates.csv" "$BAD_DATES"

# ---------------- 6. plot every bad date for review ----------------
step "6/6  plot_bad_dates.py  (day-1/day/day+1 PNG per bad date)"
"$PY" "$DIR/plot_bad_dates.py" --bad-dates "$BAD_DATES" \
      --root "$ROOT" --sfcwind-root "$SFCWIND_ROOT" --out "$PLOTS_OUT"

# ---------------- done ----------------
step "PIPELINE COMPLETE"
echo "  scan verdict     : $QC_OUT/SUMMARY.txt"
echo "  drill-down       : $DRILL_OUT/DRILLDOWN_SUMMARY.txt"
echo "  per-variable     : $ISSUES_DIR/<var>_issues.csv"
echo "  footprint        : $FOOTPRINT_OUT/footprint_classification.csv (+ footprint_maps.pdf)"
echo "  bad dates        : $BAD_DATES"
echo "  review plots     : $PLOTS_OUT/<variable>/<model>_<scenario>_<date>.png"
echo
echo "Next: flip through $PLOTS_OUT/ by variable and decide clip vs re-downscale."
