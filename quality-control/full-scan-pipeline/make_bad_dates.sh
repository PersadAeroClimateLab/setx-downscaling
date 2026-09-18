#!/usr/bin/env bash
# make_bad_dates.sh -- assemble the FINAL bad_dates list from every QC output.
# Run this LAST, after qc_zarr_all -> write_issue_csvs -> drilldown -> classify_footprint,
# so it accounts for all checks. Then feed bad_dates.csv to plot_bad_dates.py.
#
# Sources:
#   issues_by_variable/*_issues.csv     every flagged date that is NOT benign wobble
#                                       (out-of-range, incoherent/blowup, constant,
#                                        empty, missing, large-partial, outlier)
#   footprint_out/inland_cutoff_dates.csv  days missing a chunk of the inland domain
#
# Usage:
#   ./make_bad_dates.sh [ISSUES_DIR] [INLAND_CSV] [OUT]
set -euo pipefail
ISSUES_DIR="${1:-issues_by_variable}"
INLAND_CSV="${2:-footprint_out/inland_cutoff_dates.csv}"
OUT="${3:-bad_dates.csv}"

{
  echo "variable,model,scenario,date,issue,category"
  # non-wobble flagged dates from the per-variable issue CSVs (strip trailing CR)
  awk -F, 'FNR>1{sub(/\r$/,""); if($NF!="wobble") print}' "$ISSUES_DIR"/*_issues.csv
  # inland-cutoff dates from classify_footprint (5 cols -> add category)
  if [ -f "$INLAND_CSV" ]; then
    awk -F, 'FNR>1{sub(/\r$/,""); print $0",inland_cutoff"}' "$INLAND_CSV"
  fi
} > "$OUT"

echo "wrote $OUT: $(( $(wc -l < "$OUT") - 1 )) bad-date rows"
echo "  by category:"
awk -F, 'FNR>1{sub(/\r$/,""); c[$NF]++} END{for(k in c) printf "    %-12s %d\n", k, c[k]}' "$OUT"
# To restrict to hard errors only, change the first awk filter to:  $NF=="bad_values"
