#!/usr/bin/env python3
"""
Create launcher jobfile for adaptive kriging - FIXED VERSION

This script reads the config and uses the temp_dir setting to determine
where to put individual task configs and specs.

FIXED: Preserves original 4-column format: MODEL, VARIABLE, SCENARIO, DATE

Usage:
    python3 create_launcher_jobfile_adaptive_fixed.py config.json output.job
"""

import json
import sys
import os
from pathlib import Path

def create_launcher_jobfile(config_path, output_job_path):
    """Create launcher jobfile from config."""
    
    # Load config
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Get key settings from config
    working_dir = config.get('working_dir', os.getcwd())
    variable = config.get('variable', 'huss')
    extraction_list = config.get('extraction_list', 'extraction_list.txt')
    
    # CRITICAL FIX: Get temp_dir from config, with fallback
    temp_dir = config.get('temp_dir', f'temp_extractions/{variable}')
    
    # If temp_dir is absolute, use it; otherwise make it relative to working_dir
    if not os.path.isabs(temp_dir):
        temp_dir = os.path.join(working_dir, temp_dir)
    
    print(f"Config: {config_path}")
    print(f"Variable: {variable}")
    print(f"Temp dir: {temp_dir}")
    
    # Create temp directory
    os.makedirs(temp_dir, exist_ok=True)
    print(f"✓ Created directory: {temp_dir}")
    
    # Read extraction list
    if not os.path.isabs(extraction_list):
        extraction_list = os.path.join(working_dir, extraction_list)
    
    with open(extraction_list, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    # Parse timesteps - Handle 4-column format (MODEL, VARIABLE, SCENARIO, DATE)
    timesteps = []
    for line in lines:
        parts = [p.strip() for p in line.split(',')]
        
        # Handle different formats:
        # 4-column: MODEL, VARIABLE, SCENARIO, DATE
        # 2-column: MODEL, DATE (legacy)
        if len(parts) >= 4:
            # 4-column format: MODEL, VARIABLE, SCENARIO, DATE
            timesteps.append({
                'model': parts[0],
                'variable': parts[1],
                'scenario': parts[2],
                'date': parts[3].split('#')[0].strip(),  # Remove any trailing comments
                'original_line': line  # PRESERVE the original line!
            })
        elif len(parts) >= 2:
            # 2-column format: MODEL, DATE (legacy) - reconstruct full line
            timesteps.append({
                'model': parts[0],
                'date': parts[1].split('#')[0].strip(),
                'original_line': f"{parts[0]}, {variable}, historical, {parts[1].split('#')[0].strip()}"
            })
    
    print(f"Loaded {len(timesteps)} timesteps from {extraction_list}")
    
    if len(timesteps) > 0:
        print(f"  First entry: {timesteps[0]}")
    
    # Create individual config files for each task
    tasks = []
    for i, ts in enumerate(timesteps, 1):
        task_num = f"{i:06d}"
        
        # Create spec file - PRESERVE the original 4-column format
        spec_file = os.path.join(temp_dir, f"spec_{task_num}.txt")
        with open(spec_file, 'w') as f:
            f.write("# MODEL, VARIABLE, SCENARIO, DATE\n")
            f.write(ts['original_line'] + "\n")
        
        # Create task-specific config
        task_config = config.copy()
        task_config['extraction_list'] = spec_file
        
        config_file = os.path.join(temp_dir, f"config_{task_num}.json")
        with open(config_file, 'w') as f:
            json.dump(task_config, f, indent=2)
        
        # Create task command using RELATIVE path from working_dir
        rel_config = os.path.relpath(config_file, working_dir)
        tasks.append(f"cd {working_dir} && python3 custom_covariate_downscaling.py {rel_config}")
        
        if i % 1000 == 0:
            print(f"  Created task {i}/{len(timesteps)}...")
    
    # Write launcher jobfile with header
    with open(output_job_path, 'w') as f:
        f.write("# LAUNCHER JOBFILE - ADAPTIVE KRIGING\n")
        f.write(f"# Variable: {variable}\n")
        f.write("# Solution: Adaptive Regularization (Always uses covariate)\n")
        f.write(f"# Total tasks: {len(tasks)}\n")
        f.write("#\n")
        f.write("# Adaptive Parameters:\n")
        f.write("#   - Strong corr (>=0.50): nugget=0.01, nmax=40 (100%)\n")
        f.write("#   - Moderate (0.30-0.49): nugget=0.02, nmax=30 (75%)\n")
        f.write("#   - Weak (<0.30):         nugget=0.05, nmax=20 (50%)\n")
        f.write("#\n")
        f.write("# ALL tasks use Universal Kriging (Data ~ layer)\n")
        f.write("#\n")
        for task in tasks:
            f.write(task + "\n")
    
    print(f"✓ Created launcher jobfile: {output_job_path}")
    print(f"✓ Total tasks: {len(tasks)}")
    print(f"✓ Individual configs in: {temp_dir}")
    
    # Verify first spec file
    if len(timesteps) > 0:
        first_spec = os.path.join(temp_dir, "spec_000001.txt")
        with open(first_spec, 'r') as f:
            content = f.read().strip()
        print(f"✓ First spec file content:\n{content}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 create_launcher_jobfile_adaptive_fixed.py config.json output.job")
        sys.exit(1)
    
    config_path = sys.argv[1]
    output_job_path = sys.argv[2]
    
    create_launcher_jobfile(config_path, output_job_path)

if __name__ == "__main__":
    main()
