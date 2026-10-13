#!/bin/bash
#===============================================================================
# check_1km_coverage.sh
#
# Scan a directory for the downscaled 1km NetCDF files, report any whose
# valid-pixel coverage is BELOW a threshold (default 89%), and optionally
# DELETE the whole timestep set for each flagged file so it can be re-run.
#
#   - Includes files matching  *1km*.nc
#   - EXCLUDES the standard-error files whose name starts with  SE_
#   - Coverage = 100 * (gridsize - missing) / gridsize, read from `cdo info`
#     (missing = NaN / _FillValue cells). A file cdo cannot read is flagged.
#
# For each flagged 1km file  MODEL_SCEN_VAR_1km_YYYYMMDD.nc  the deletion set
# (same model, scenario, variable, date) is the four companion files:
#     MODEL_SCEN_VAR_1km_YYYYMMDD.nc        (the flagged 1km itself)
#     SE_MODEL_SCEN_VAR_1km_YYYYMMDD.nc
#     MODEL_SCEN_VAR_9km_YYYYMMDD.nc
#     SE_MODEL_SCEN_VAR_9km_YYYYMMDD.nc
# Exact filenames are used - no wildcard deletes. Deletion is DRY-RUN unless -D.
#
# Usage:
#   ./check_1km_coverage.sh [-t THRESHOLD] [-o OUTFILE] [-l LISTFILE] [-D] [-r] [DIR]
#     -t THRESHOLD  minimum acceptable coverage percent (default 89)
#     -o OUTFILE    write the flagged 1km file PATHS (one per line) to OUTFILE
#     -l LISTFILE   write the flagged files as a reprocessing extraction list
#                   in "MODEL, VARIABLE, SCENARIO, DATE" format
#     -D            DELETE the companion set for each flagged file (default:
#                   dry-run, which only prints what WOULD be deleted)
#     -e            exit non-zero if any file was flagged (for CI). DEFAULT is
#                   exit 0, so this is safe inside `cmd && cmd` chains.
#     -r            recurse into sub-directories
#     DIR           directory to scan (default: current directory)
#
# Filename scenario short -> full:  hist->historical, ssp1->ssp126,
#                                   ssp2->ssp245, ssp3->ssp370, ssp5->ssp585
#
# Examples:
#   ./check_1km_coverage.sh downscaled/hurs                       # scan only
#   ./check_1km_coverage.sh -l rerun.txt downscaled/hurs          # scan + re-run list
#   ./check_1km_coverage.sh -l rerun.txt downscaled/hurs          # (review the dry-run)
#   ./check_1km_coverage.sh -l rerun.txt -D downscaled/hurs       # scan + DELETE sets
#===============================================================================

set -u

THRESHOLD=89
OUTFILE=""
LISTFILE=""
RECURSE=0
DO_DELETE=0
EXIT_ON_FLAG=0     # default: always exit 0 so this is safe inside `cmd && cmd` chains

while getopts ":t:o:l:Derh" opt; do
  case "$opt" in
    t) THRESHOLD="$OPTARG" ;;
    o) OUTFILE="$OPTARG" ;;
    l) LISTFILE="$OPTARG" ;;
    D) DO_DELETE=1 ;;
    e) EXIT_ON_FLAG=1 ;;
    r) RECURSE=1 ;;
    h) sed -n '2,46p' "$0" | sed 's/^#\{1,\} \{0,1\}//'; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG (use -h for help)" >&2; exit 1 ;;
    :)  echo "Option -$OPTARG needs a value (use -h for help)" >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))
DIR="${1:-.}"

# cdo is required
if ! command -v cdo >/dev/null 2>&1; then
  module load cdo 2>/dev/null || true
fi
command -v cdo >/dev/null 2>&1 || { echo "ERROR: cdo not found. Run 'module load cdo' first." >&2; exit 1; }

[[ -d "$DIR" ]] || { echo "ERROR: not a directory: $DIR" >&2; exit 1; }

# Turn a 1km filename into a "MODEL, VARIABLE, SCENARIO, DATE" extraction line.
to_rerun_line() {
  local base front date8 var rest scen_short scen y m d model
  base="$(basename "$1")"; base="${base%.nc}"
  [[ "$base" == *_1km_* ]] || return 1
  front="${base%_1km_*}"; date8="${base##*_1km_}"
  [[ "$date8" =~ ^[0-9]{8}$ ]] || return 1
  var="${front##*_}"; rest="${front%_*}"; scen_short="${rest##*_}"; model="${rest%_*}"
  [[ -n "$model" && -n "$var" && -n "$scen_short" ]] || return 1
  case "$scen_short" in
    hist) scen="historical" ;;
    ssp1) scen="ssp126" ;;
    ssp2) scen="ssp245" ;;
    ssp3) scen="ssp370" ;;
    ssp5) scen="ssp585" ;;
    *)    scen="$scen_short"
          echo "  WARN: unknown scenario '$scen_short' in $base (passed through)" >&2 ;;
  esac
  y="${date8:0:4}"; m="${date8:4:2}"; d="${date8:6:2}"
  printf '%s, %s, %s, %s-%s-%s\n' "$model" "$var" "$scen" "$y" "$m" "$d"
}

# Delete (or, in dry-run, list) the companion set for one flagged 1km file.
# Uses EXACT filenames in the flagged file's own directory - never a wildcard.
handle_companions() {
  local f dir base stem date8 name path
  f="$1"; dir="$(dirname "$f")"; base="$(basename "$f")"; base="${base%.nc}"
  [[ "$base" == *_1km_* ]] || { echo "  WARN: cannot parse, skipped: $f" >&2; return; }
  stem="${base%_1km_*}"          # MODEL_SCEN_VAR
  date8="${base##*_1km_}"        # YYYYMMDD
  for name in \
      "${stem}_1km_${date8}.nc" \
      "SE_${stem}_1km_${date8}.nc" \
      "${stem}_9km_${date8}.nc" \
      "SE_${stem}_9km_${date8}.nc" ; do
    path="${dir}/${name}"
    if [[ -f "$path" ]]; then
      if [[ "$DO_DELETE" -eq 1 ]]; then
        if rm -f "$path"; then echo "    deleted:      $path"; deleted=$((deleted + 1))
        else                   echo "    FAILED delete: $path" >&2; fi
      else
        echo "    would delete: $path"; would=$((would + 1))
      fi
    fi
  done
}

# Collect candidate files: *1km*.nc, excluding basenames starting with SE_.
if [[ "$RECURSE" -eq 1 ]]; then
  mapfile -d '' FILES < <(find "$DIR" -type f -name '*1km*.nc' ! -name 'SE_*' -print0 | sort -z)
else
  mapfile -d '' FILES < <(find "$DIR" -maxdepth 1 -type f -name '*1km*.nc' ! -name 'SE_*' -print0 | sort -z)
fi

total_files=${#FILES[@]}
if [[ "$total_files" -eq 0 ]]; then
  echo "No *1km*.nc files (excluding SE_*) found in: $DIR"
  exit 0
fi

echo "Scanning $total_files file(s) in: $DIR"
echo "Flagging files with valid-pixel coverage < ${THRESHOLD}%"
if [[ "$DO_DELETE" -eq 1 ]]; then
  echo "DELETE MODE: companion sets (1km, SE_1km, 9km, SE_9km) WILL be removed."
else
  echo "DRY-RUN: no files will be deleted (add -D to delete)."
fi
echo "----------------------------------------------------------------------"

tmpout="$(mktemp)"; tmplist="$(mktemp)"
checked=0; flagged=0; errors=0; deleted=0; would=0

flag_file() {
  echo "$1" >> "$tmpout"
  if [[ -n "$LISTFILE" ]]; then
    to_rerun_line "$1" >> "$tmplist" || echo "  WARN: could not parse for re-run list: $1" >&2
  fi
  handle_companions "$1"
}

for f in "${FILES[@]}"; do
  read -r gridsize miss < <(cdo -s info "$f" 2>/dev/null \
      | awk '/^[[:space:]]*[0-9]+[[:space:]]*:/ {print $6, $7; exit}')

  if [[ -z "${gridsize:-}" || -z "${miss:-}" || "$gridsize" -eq 0 ]]; then
    echo "  ERROR  (unreadable, flagged): $f"
    flag_file "$f"
    errors=$((errors + 1)); flagged=$((flagged + 1)); continue
  fi

  checked=$((checked + 1))
  pct=$(awk -v g="$gridsize" -v m="$miss" 'BEGIN{ printf "%.2f", 100*(g-m)/g }')
  if awk -v p="$pct" -v t="$THRESHOLD" 'BEGIN{ exit !(p < t) }'; then
    printf "  %6.2f%%  %s\n" "$pct" "$f"
    flag_file "$f"
    flagged=$((flagged + 1))
  fi
done

echo "----------------------------------------------------------------------"
echo "Checked: $checked   Flagged (<${THRESHOLD}%): $flagged   Read errors: $errors"
if [[ "$DO_DELETE" -eq 1 ]]; then
  echo "Files deleted: $deleted"
else
  echo "Files that WOULD be deleted: $would   (re-run with -D to delete)"
fi

if [[ -n "$OUTFILE" ]]; then
  sort -u "$tmpout" > "$OUTFILE"
  echo "Flagged file paths written to: $OUTFILE"
fi
if [[ -n "$LISTFILE" ]]; then
  sort -u "$tmplist" > "$LISTFILE"
  echo "Reprocessing list (MODEL, VARIABLE, SCENARIO, DATE) written to: $LISTFILE ($(wc -l < "$LISTFILE") lines)"
fi
rm -f "$tmpout" "$tmplist"

# Exit 0 by default so this step is safe inside `cmd && cmd` chains (flagging is
# not an error). Use -e if you WANT a non-zero exit when files were flagged.
if [[ "$EXIT_ON_FLAG" -eq 1 ]]; then
  [[ "$flagged" -eq 0 ]]
else
  exit 0
fi
