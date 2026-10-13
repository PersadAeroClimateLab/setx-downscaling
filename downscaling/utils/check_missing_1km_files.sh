#!/bin/bash

################################################################################
# Script: check_missing_1km_files.sh
# Purpose: Check for missing 1km downscaled files by comparing expected files
#          from a master list against files actually present in a directory
#
# WRITE RULE:
#   An entry is written to the missing-file list if the variable file is
#   missing, OR the SE file is missing, OR both are missing. Only entries with
#   BOTH files present are treated as "found" and skipped.
#
# How it works:
#   - Checks BOTH (ssp1, ssp2, ssp3, ssp5) and (ssp126, ssp245, ssp370, ssp585) naming
#   - Checks BOTH date formats (YYYYMMDD and YYYY-MM-DD)
#   - Data lines stay clean 4-column (MODEL, VARIABLE, SCENARIO, DATE) so the
#     output feeds directly into setup_*.sh
#
# Usage: ./check_missing_1km_files_fixed.sh <master_list> <output_dir> [output_file]
################################################################################

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Function to display usage
usage() {
    echo "Usage: $0 <master_list> <output_dir> [output_file]"
    echo ""
    echo "Arguments:"
    echo "  master_list  - Input file with expected entries (e.g., tas_by_model.txt)"
    echo "  output_dir   - Directory containing actual *_1km_*.nc files"
    echo "  output_file  - (Optional) Output file for missing entries (default: missing_files.txt)"
    exit 1
}

# Check arguments
if [ $# -lt 2 ]; then
    echo -e "${RED}Error: Insufficient arguments${NC}"
    usage
fi

MASTER_LIST="$1"
OUTPUT_DIR="$2"
MISSING_FILE="${3:-missing_files.txt}"

# Validate inputs
if [ ! -f "$MASTER_LIST" ]; then
    echo -e "${RED}Error: Master list file not found: $MASTER_LIST${NC}"
    exit 1
fi

if [ ! -d "$OUTPUT_DIR" ]; then
    echo -e "${RED}Error: Output directory not found: $OUTPUT_DIR${NC}"
    exit 1
fi

# Print header
echo -e "${BLUE}╔═════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║      Checking for Missing 1km Downscaled Files (FIXED)          ║${NC}"
echo -e "${BLUE}╚═════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${YELLOW}Master list:${NC} $MASTER_LIST"
echo -e "${YELLOW}Output dir:${NC}  $OUTPUT_DIR"
echo -e "${YELLOW}Missing file:${NC} $MISSING_FILE"
echo ""

# Build a hash of existing files for fast lookup
echo -e "${BLUE}Building file index...${NC}"
declare -A existing_files

# Index all 1km files (excluding SE_ files)
while IFS= read -r filepath; do
    filename=$(basename "$filepath")
    # Skip SE_ files for main file lookup
    if [[ ! "$filename" =~ ^SE_ ]]; then
        existing_files["$filename"]=1
    fi
done < <(find "$OUTPUT_DIR" -maxdepth 1 -name "*_1km_*.nc" -type f)

# Index SE files separately
declare -A existing_se_files
while IFS= read -r filepath; do
    filename=$(basename "$filepath")
    existing_se_files["$filename"]=1
done < <(find "$OUTPUT_DIR" -maxdepth 1 -name "SE_*_1km_*.nc" -type f)

echo -e "  Found ${#existing_files[@]} variable files"
echo -e "  Found ${#existing_se_files[@]} SE files"
echo ""

# Initialize counters
total_expected=0
total_found=0
total_missing=0
missing_both=0
missing_var_only=0
missing_se_only=0
written_to_file=0

# Clear/create output file
echo "# Missing 1km Downscaled Files - Format: MODEL, VARIABLE, SCENARIO, DATE" > "$MISSING_FILE"
echo "# Rule: entry listed if variable file, SE file, or both are missing" >> "$MISSING_FILE"

echo -e "${BLUE}Processing entries...${NC}"

# Read the master list and check each entry
while IFS= read -r line; do
    # Skip empty lines and comments
    if [[ -z "$line" ]] || [[ "$line" =~ ^# ]]; then
        continue
    fi

    # Parse the line: MODEL, VARIABLE, SCENARIO, DATE
    line_clean=$(echo "$line" | sed 's/ *, */,/g')
    IFS=',' read -r model variable scenario date <<< "$line_clean"

    # Trim whitespace
    model=$(echo "$model" | xargs)
    variable=$(echo "$variable" | xargs)
    scenario=$(echo "$scenario" | xargs)
    date=$(echo "$date" | xargs)

    # Skip if any field is empty
    if [[ -z "$model" ]] || [[ -z "$variable" ]] || [[ -z "$scenario" ]] || [[ -z "$date" ]]; then
        continue
    fi

    ((total_expected++))

    # Convert date format: YYYY-MM-DD -> YYYYMMDD
    date_compact=$(echo "$date" | sed 's/-//g')

    # Generate ALL possible filename variations
    # Scenario variations: hist, ssp1/ssp126, ssp2/ssp245, ssp3/ssp370, ssp5/ssp585
    declare -a scenario_variants=()

    if [[ "$scenario" == "historical" ]]; then
        scenario_variants=("hist")
    elif [[ "$scenario" == "ssp126" ]]; then
        scenario_variants=("ssp126" "ssp1")
    elif [[ "$scenario" == "ssp245" ]]; then
        scenario_variants=("ssp245" "ssp2")
    elif [[ "$scenario" == "ssp370" ]]; then
        scenario_variants=("ssp370" "ssp3")
    elif [[ "$scenario" == "ssp585" ]]; then
        scenario_variants=("ssp585" "ssp5")
    else
        scenario_variants=("$scenario")
    fi

    # Check if file exists (any combination of scenario variant + date format)
    file_exists=false
    se_file_exists=false

    for scen in "${scenario_variants[@]}"; do
        # Check both date formats
        for dt in "$date_compact" "$date"; do
            var_file="${model}_${scen}_${variable}_1km_${dt}.nc"
            se_file="SE_${model}_${scen}_${variable}_1km_${dt}.nc"

            if [[ -n "${existing_files[$var_file]}" ]]; then
                file_exists=true
            fi
            if [[ -n "${existing_se_files[$se_file]}" ]]; then
                se_file_exists=true
            fi
        done
    done

    #--------------------------------------------------------------------------
    # WRITE RULE: found only if BOTH present; otherwise (var OR se OR both
    # missing) the entry is written to the missing-file list.
    #--------------------------------------------------------------------------
    if [ "$file_exists" = true ] && [ "$se_file_exists" = true ]; then
        ((total_found++))
    else
        ((total_missing++))

        # Write the clean original 4-column line for ANY missing case
        echo "$line" >> "$MISSING_FILE"
        ((written_to_file++))

        # Classify which file(s) were missing (for the summary only)
        if [ "$file_exists" = false ] && [ "$se_file_exists" = false ]; then
            ((missing_both++))
        elif [ "$file_exists" = false ]; then
            ((missing_var_only++))
        else
            ((missing_se_only++))
        fi
    fi

    # Progress indicator
    if (( total_expected % 500 == 0 )); then
        echo -ne "\r  Processed: ${total_expected} | Found: ${total_found} | Missing: ${total_missing}"
    fi

done < "$MASTER_LIST"

echo -ne "\r  Processed: ${total_expected} | Found: ${total_found} | Missing: ${total_missing}\n"
echo ""

# Calculate percentages
if [ $total_expected -gt 0 ]; then
    found_pct=$(awk "BEGIN {printf \"%.1f\", ($total_found/$total_expected)*100}")
    missing_pct=$(awk "BEGIN {printf \"%.1f\", ($total_missing/$total_expected)*100}")
else
    found_pct="0.0"
    missing_pct="0.0"
fi

# Print summary
echo -e "${BLUE}╔═════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                         SUMMARY                                 ║${NC}"
echo -e "${BLUE}╚═════════════════════════════════════════════════════════════════╝${NC}"
echo ""
printf "  ${YELLOW}Total expected:${NC}  %6d (100.0%%)\n" $total_expected
printf "  ${GREEN}Files found:${NC}     %6d (%5s%%) ${GREEN}[Both files present]${NC}\n" $total_found "$found_pct"
printf "  ${RED}Files missing:${NC}   %6d (%5s%%) ${RED}[One or both files missing - ALL written]${NC}\n" $total_missing "$missing_pct"
echo ""

if [ $total_missing -gt 0 ]; then
    echo -e "${CYAN}Breakdown of Missing Files (all written to list):${NC}"
    if [ $missing_both -gt 0 ]; then
        printf "  ${RED}Both files missing:${NC}      %6d ${GREEN}[WRITTEN]${NC}\n" $missing_both
    fi
    if [ $missing_var_only -gt 0 ]; then
        printf "  ${RED}Variable file only:${NC}      %6d ${GREEN}[WRITTEN]${NC}\n" $missing_var_only
    fi
    if [ $missing_se_only -gt 0 ]; then
        printf "  ${YELLOW}SE file only:${NC}            %6d ${GREEN}[WRITTEN]${NC}\n" $missing_se_only
    fi
    echo ""
fi

echo -e "${CYAN}Output File Summary:${NC}"
printf "  ${YELLOW}Written to file:${NC}         %6d entries\n" $written_to_file
echo ""

if [ $written_to_file -gt 0 ]; then
    echo -e "${YELLOW}Missing files written to:${NC} $MISSING_FILE"
    echo ""
    echo -e "${YELLOW}Preview (first 10):${NC}"
    grep -v "^#" "$MISSING_FILE" | head -10
else
    echo -e "${GREEN}✓ All entries complete (both variable and SE files present)!${NC}"
    rm -f "$MISSING_FILE"
fi

echo ""
echo -e "${BLUE}Done!${NC}"
