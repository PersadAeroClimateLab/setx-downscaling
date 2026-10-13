#!/usr/bin/env python3
"""
Build a TACC launcher job file for the VP family (huss, rsds, rlds). Reads a
date-based extraction list (MODEL, VARIABLE, SCENARIO, DATE) and emits one task
line per timestep that calls the driver:

    python3 custom_covariate_downscaling_vp_gated.py <per-task config>

Usage:
    python3 create_launcher_jobfile_vp_gated.py config.json output.job
"""

import json
import sys
import os

# The driver this pipeline runs.
DRIVER = "custom_covariate_downscaling_vp_gated.py"


def create_launcher_jobfile(config_path, output_job_path):
    """Create launcher jobfile from config."""

    with open(config_path, 'r') as f:
        config = json.load(f)

    working_dir = config.get('working_dir', os.getcwd())
    variable = config.get('variable', 'huss')
    extraction_list = config.get('extraction_list', 'extraction_list.txt')

    temp_dir = config.get('temp_dir', f'temp_extractions/{variable}')
    if not os.path.isabs(temp_dir):
        temp_dir = os.path.join(working_dir, temp_dir)

    print(f"Config: {config_path}")
    print(f"Variable: {variable}")
    print(f"Driver: {DRIVER}")
    print(f"Temp dir: {temp_dir}")

    os.makedirs(temp_dir, exist_ok=True)
    print(f"Created directory: {temp_dir}")

    if not os.path.isabs(extraction_list):
        extraction_list = os.path.join(working_dir, extraction_list)

    with open(extraction_list, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    # Parse timesteps - handle 4-column (MODEL, VARIABLE, SCENARIO, DATE)
    timesteps = []
    for line in lines:
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 4:
            timesteps.append({
                'model': parts[0],
                'variable': parts[1],
                'scenario': parts[2],
                'date': parts[3].split('#')[0].strip(),
                'original_line': line
            })
        elif len(parts) >= 2:
            # 2-column legacy: MODEL, DATE -> reconstruct full 4-column line
            timesteps.append({
                'model': parts[0],
                'date': parts[1].split('#')[0].strip(),
                'original_line': f"{parts[0]}, {variable}, historical, {parts[1].split('#')[0].strip()}"
            })

    print(f"Loaded {len(timesteps)} timesteps from {extraction_list}")
    if len(timesteps) > 0:
        print(f"  First entry: {timesteps[0]}")

    tasks = []
    for i, ts in enumerate(timesteps, 1):
        task_num = f"{i:06d}"

        # Spec file - PRESERVE the 4-column format
        spec_file = os.path.join(temp_dir, f"spec_{task_num}.txt")
        with open(spec_file, 'w') as f:
            f.write("# MODEL, VARIABLE, SCENARIO, DATE\n")
            f.write(ts['original_line'] + "\n")

        # Task-specific config
        task_config = config.copy()
        task_config['extraction_list'] = spec_file
        config_file = os.path.join(temp_dir, f"config_{task_num}.json")
        with open(config_file, 'w') as f:
            json.dump(task_config, f, indent=2)

        rel_config = os.path.relpath(config_file, working_dir)
        tasks.append(f"cd {working_dir} && python3 {DRIVER} {rel_config}")

        if i % 1000 == 0:
            print(f"  Created task {i}/{len(timesteps)}...")

    with open(output_job_path, 'w') as f:
        f.write("# LAUNCHER JOBFILE - ADAPTIVE KRIGING (GATED) - VP FAMILY\n")
        f.write(f"# Variable: {variable}\n")
        f.write(f"# Driver: {DRIVER}\n")
        f.write("# Covariate: monthly vapor pressure\n")
        f.write("# Production coverage gate enforced per iteration by the driver.\n")
        f.write(f"# Total tasks: {len(tasks)}\n")
        f.write("#\n")
        for task in tasks:
            f.write(task + "\n")

    print(f"Created launcher jobfile: {output_job_path}")
    print(f"Total tasks: {len(tasks)}")
    print(f"Individual configs in: {temp_dir}")

    if len(timesteps) > 0:
        first_spec = os.path.join(temp_dir, "spec_000001.txt")
        with open(first_spec, 'r') as f:
            print(f"First spec file content:\n{f.read().strip()}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 create_launcher_jobfile_vp_gated.py config.json output.job")
        sys.exit(1)
    create_launcher_jobfile(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
