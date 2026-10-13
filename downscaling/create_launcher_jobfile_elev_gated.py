#!/usr/bin/env python3
"""
Build a TACC launcher job file for the elevation family (tas, tasmin, tasmax,
sfcWind). Reads a date-based extraction list (MODEL, VARIABLE, SCENARIO, DATE)
and emits one task line per timestep that calls the driver:

    python3 custom_covariate_downscaling_elev_gated.py <per-task config>

The covariate is GMTED elevation.

Usage:
    python3 create_launcher_jobfile_elev_gated.py <config.json> <output.job>
"""

import json
import sys
import os
from pathlib import Path
from collections import Counter

# Supported variables for this pipeline
ELEVATION_VARS = ['tas', 'tasmin', 'tasmax', 'sfcWind']

# The gated driver this pipeline runs. Change here if the driver is renamed.
DRIVER = "custom_covariate_downscaling_elev_gated.py"


def create_launcher_jobfile(config_file, output_jobfile):
    """Create launcher jobfile from date-based extraction list."""

    print(f"{'='*70}")
    print("ADAPTIVE KRIGING (GATED) - LAUNCHER JOBFILE CREATOR")
    print("Elevation-based Variables (tas, tasmin, tasmax, sfcWind)")
    print(f"{'='*70}")
    print(f"Config: {config_file}")
    print(f"Output: {output_jobfile}")
    print(f"Driver: {DRIVER}")
    print()

    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print(f"ERROR: Config file not found: {config_file}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"ERROR: Invalid JSON in config file: {e}")
        sys.exit(1)

    working_dir = Path(config.get('working_dir', '.'))
    extraction_list_name = config.get('extraction_list')
    r_template = config.get('r_script_template', 'kriging_template_adaptive_elev_gated.R')

    variable = config.get('variable')
    if not variable:
        for var in ELEVATION_VARS:
            if var in extraction_list_name:
                variable = var
                break
        if not variable:
            variable = 'tas'

    # Verify gated adaptive template
    if 'adaptive' not in r_template.lower():
        print(f"WARNING: Template name doesn't contain 'adaptive': {r_template}")
    if 'gated' not in r_template.lower():
        print(f"WARNING: Template name doesn't contain 'gated': {r_template}")
        print("         Make sure you are using kriging_template_adaptive_elev_gated.R")

    if not extraction_list_name:
        print("ERROR: 'extraction_list' not found in config")
        sys.exit(1)

    extraction_list = working_dir / extraction_list_name

    temp_dir = config.get('temp_dir')
    if not temp_dir:
        temp_dir = working_dir / f'temp_extractions/{variable}'
    else:
        temp_dir = Path(temp_dir)

    print(f"Working directory: {working_dir}")
    print(f"Variable: {variable}")
    print(f"Extraction list:   {extraction_list}")
    print(f"Temp directory:    {temp_dir}")
    print(f"R template:        {r_template}")
    print()

    if not extraction_list.exists():
        print(f"ERROR: Extraction list not found: {extraction_list}")
        sys.exit(1)

    specs = []
    print("Reading extraction list...")

    with open(extraction_list, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 4:
                model = parts[0]
                var = parts[1]
                scenario = parts[2]
                date_str = parts[3]
                try:
                    year, month, day = date_str.split('-')
                    month_number = int(month)
                    specs.append({
                        'model': model,
                        'variable': var,
                        'scenario': scenario,
                        'date': date_str,
                        'month_number': month_number,
                        'line': line
                    })
                    if len(specs) % 1000 == 0:
                        print(f"  Read {len(specs)} timesteps...")
                except ValueError as e:
                    print(f"  WARNING: Line {line_num}: Could not parse date '{date_str}': {e}")
                    continue
            else:
                print(f"  WARNING: Line {line_num}: Expected 4 columns, got {len(parts)}")

    if not specs:
        print("ERROR: No timesteps found in extraction list!")
        sys.exit(1)

    print(f"\nFound {len(specs)} timesteps to process")

    month_counts = Counter([s['month_number'] for s in specs])
    month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    print("\nDistribution by month:")
    for month in sorted(month_counts.keys()):
        print(f"  {month_names[month-1]:3s}: {month_counts[month]:4d} timesteps")

    dates = sorted([s['date'] for s in specs])
    print(f"\nDate range: {dates[0]} to {dates[-1]}")

    temp_dir.mkdir(exist_ok=True, parents=True)
    print(f"\nCreated directory: {temp_dir}")

    output_path = Path(output_jobfile)
    print(f"\nCreating launcher jobfile: {output_path}")

    with open(output_path, 'w') as f:
        f.write("# LAUNCHER JOBFILE - ADAPTIVE KRIGING (GATED) - ELEVATION FAMILY\n")
        f.write(f"# Variable: {variable}\n")
        f.write(f"# Driver: {DRIVER}\n")
        f.write("# Covariate: GMTED Elevation\n")
        f.write("# Production coverage gate enforced per iteration by the driver.\n")
        f.write(f"# Total tasks: {len(specs)}\n")
        f.write("#\n")

        for i, spec in enumerate(specs, 1):
            spec_file = temp_dir / f"spec_{i:06d}.txt"
            spec_config = temp_dir / f"config_{i:06d}.json"

            with open(spec_file, 'w') as sf:
                sf.write("# MODEL, VARIABLE, SCENARIO, DATE\n")
                sf.write(spec['line'] + "\n")

            task_config = config.copy()
            task_config['extraction_list'] = str(spec_file)
            with open(spec_config, 'w') as cf:
                json.dump(task_config, cf, indent=2)

            rel_config = os.path.relpath(spec_config, working_dir)
            cmd = f"cd {working_dir} && python3 {DRIVER} {rel_config}"
            f.write(f"{cmd}\n")

            if i % 1000 == 0:
                print(f"  Created task {i}/{len(specs)}...")

    print(f"\nCreated launcher jobfile: {output_jobfile}")
    print(f"Total tasks: {len(specs)}")
    print(f"Individual configs in: {temp_dir}")

    first_spec = temp_dir / "spec_000001.txt"
    if first_spec.exists():
        print("\nFirst spec file content:")
        with open(first_spec, 'r') as f:
            print(f.read().strip())

    print(f"\n{'='*70}")
    print("USAGE")
    print(f"{'='*70}")
    print("Run in idev session:")
    print(f"  cd {working_dir}")
    print("  module load launcher tacc-apptainer cdo nco ncl_ncarg")
    print("  export LAUNCHER_WORKDIR=$PWD")
    print(f"  export LAUNCHER_JOB_FILE={output_jobfile}")
    print("  export LAUNCHER_SCHED=interleaved")
    print("  $LAUNCHER_DIR/paramrun")
    print(f"{'='*70}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 create_launcher_jobfile_elev_gated.py <config.json> <output_jobfile>")
        print("\nExample:")
        print("  python3 create_launcher_jobfile_elev_gated.py \\")
        print("      cmip6_config_tas_gated.json \\")
        print("      tas_gated_tasks.job")
        sys.exit(1)

    create_launcher_jobfile(sys.argv[1], sys.argv[2])
