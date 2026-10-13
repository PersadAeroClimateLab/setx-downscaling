#!/usr/bin/env python3
"""
Custom Covariate Downscaling Wrapper
Handles two-step downscaling with month-specific covariates
"""

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
import logging
from datetime import datetime

class CustomCovariateDownscaler:
    """Handle two-step downscaling with custom covariates."""

    def __init__(self, config_file):
        """Initialize with configuration."""
        with open(config_file, 'r') as f:
            self.config = json.load(f)

        self.working_dir = Path(self.config['working_dir'])
        self.output_dir = Path(self.config['output_dir'])
        self.temp_dir = Path(self.config['temp_dir'])

        # Coverage gate: a downscale passes only if at least this fraction of the
        # area carries data (equivalently NaN <= 1 - min_coverage). Default 0.9988
        # (i.e. <= 0.12% holes). 'coverage_denominator' selects what "the area" is:
        #   "domain" = cells expected to have data (all covariates finite)
        #   "grid"   = every cell in the target bounding box
        self.min_coverage = float(self.config.get('min_coverage', 0.9988))
        self.coverage_denominator = str(
            self.config.get('coverage_denominator', 'domain')).lower()
        if self.coverage_denominator not in ('domain', 'grid'):
            self.coverage_denominator = 'domain'

        # Physical range used by the R result-scan (hurs = 0..100 %). A kriging
        # attempt whose values fall well outside this band is rejected and the
        # next method is tried; an accepted field is clamped to this range.
        self.valid_min = self.config.get('valid_min', 0)
        self.valid_max = self.config.get('valid_max', 100)

        # Per-iteration wall-clock limit for the R/kriging call. A hung task is
        # killed after this many seconds instead of blocking a launcher worker
        # forever. 0 or negative disables the guard. NOTE: 'timeout' can send
        # SIGKILL, but a process stuck in an uninterruptible filesystem wait
        # (Lustre stall) may still not die - this guards compute hangs, not a
        # hard $SCRATCH outage.
        self.r_timeout_seconds = int(self.config.get('r_timeout_seconds', 3600))

        # Node-local scratch for INTERMEDIATE files (extracted timestep + the
        # generated R script). Writing these to $TMPDIR (node-local disk) instead
        # of $SCRATCH removes ~2 small-file writes per timestep from Lustre and
        # cuts the metadata load that can stall a large launcher run. Final
        # outputs still go to output_dir on $SCRATCH.
        tmp_base = Path(os.environ.get('TMPDIR', '/tmp'))
        var_tag = str(self.config.get('variable', 'redown'))
        self.local_tmp = tmp_base / f"redown_{var_tag}_{os.getpid()}"

        # Create directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        try:
            self.local_tmp.mkdir(parents=True, exist_ok=True)
        except OSError:
            # Fall back to $SCRATCH temp_dir if node-local disk is unavailable.
            self.local_tmp = self.temp_dir

        # Setup logging
        log_file = self.working_dir / f"custom_downscaling_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Node-local temp dir: {self.local_tmp}")
        if self.r_timeout_seconds > 0:
            self.logger.info(f"R per-iteration timeout: {self.r_timeout_seconds}s")
        else:
            self.logger.info("R per-iteration timeout: disabled")

    def load_extraction_list(self):
        """Load extraction list with date-based format."""
        extraction_file = self.working_dir / self.config['extraction_list']
        specs = []

        self.logger.info(f"Loading extraction list: {extraction_file}")

        with open(extraction_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # Parse: MODEL, VARIABLE, SCENARIO, DATE
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 4:
                    model = parts[0]
                    variable = parts[1]
                    scenario = parts[2]
                    date_str = parts[3]

                    # Extract month from date (YYYY-MM-DD)
                    try:
                        year, month, day = date_str.split('-')
                        month_number = int(month)

                        specs.append({
                            'model': model,
                            'variable': variable,
                            'scenario': scenario,
                            'date': date_str,
                            'year': int(year),
                            'month_number': month_number,
                            'day': int(day)
                        })
                    except ValueError as e:
                        self.logger.warning(f"Line {line_num}: Could not parse date '{date_str}': {e}")
                        continue

        self.logger.info(f"Loaded {len(specs)} timesteps")
        return specs

    def create_r_script(self, spec, input_file, output_file, se_file,
                       target_res, iteration):
        """Create R script for one iteration."""

        template_file = self.working_dir / self.config['r_script_template']
        with open(template_file, 'r') as f:
            template = f.read()

        # Select GMTED file based on iteration
        if iteration == 1:
            gmted_file = self.config.get('gmted_iter1', self.config['gmted_file'])
        else:
            gmted_file = self.config.get('gmted_iter2', self.config['gmted_file'])

        # Handle covariate configurations with model substitution
        # Check if using single covariate (old format) or multiple covariates (new format)
        if 'covariates' in self.config:
            # Multi-covariate format - substitute {model}/{scenario} placeholders
            import json
            covariates = []
            for cov in self.config['covariates']:
                cov_copy = cov.copy()
                # Replace model and scenario placeholders with actual values.
                # {model}/{MODEL}  -> spec['model']    (e.g. MRI-ESM2-0)
                # {scenario}/{SCENARIO} -> spec['scenario'] (full word, e.g. historical)
                # NOTE: {variable} is intentionally NOT substituted here. A
                # covariate's own variable (e.g. tas) differs from the target
                # variable (e.g. hurs), so write the covariate variable literally
                # in the path. The r* realization wildcard is left in place and
                # resolved by glob inside the R template at load time.
                if 'path' in cov_copy:
                    original_path = cov_copy['path']
                    # tas covariate filenames may use either the full CMIP6
                    # scenario name ("historical") or the short form ("hist").
                    # For those, substitute the glob "hist*" so the R template's
                    # Sys.glob matches whichever form exists on disk. Other
                    # scenarios (ssp*) are substituted as-is.
                    if spec['scenario'] in ('hist', 'historical'):
                        scenario_for_cov = 'hist*'
                    else:
                        scenario_for_cov = spec['scenario']
                    cov_copy['path'] = (cov_copy['path']
                        .replace('{MODEL}', spec['model'])
                        .replace('{model}', spec['model'])
                        .replace('{SCENARIO}', scenario_for_cov)
                        .replace('{scenario}', scenario_for_cov))
                    if cov_copy['path'] != original_path:
                        self.logger.info(
                            f"Resolved covariate path: {cov['name']} for "
                            f"{spec['model']} ({spec['scenario']}) -> {cov_copy['path']}")
                covariates.append(cov_copy)
            covariate_configs = json.dumps(covariates)
        elif 'covariate_base_path' in self.config:
            # Single covariate format - convert to multi-covariate format
            import json
            single_cov = [{
                'name': 'covariate',
                'type': 'monthly_file',
                'path': self.config['covariate_base_path']
            }]
            covariate_configs = json.dumps(single_cov)
        else:
            self.logger.error("No covariate configuration found")
            return None

        # Replace placeholders
        script_content = template.replace('{{INPUT_FILE}}', str(input_file))
        script_content = script_content.replace('{{OUTPUT_FILE}}', str(output_file))
        script_content = script_content.replace('{{SE_FILE}}', str(se_file))
        script_content = script_content.replace('{{VARIABLE}}', spec['variable'])
        script_content = script_content.replace('{{DATE}}', spec['date'])
        script_content = script_content.replace('{{TARGET_RESOLUTION}}', str(target_res))
        script_content = script_content.replace('{{TIMESTEP_INDEX}}', 'NULL')
        # Accept either 'nmax' (multi-covariate configs) or 'nmax_base'
        # (adaptive single-covariate configs) so this driver works with both.
        nmax_value = self.config.get('nmax', self.config.get('nmax_base'))
        script_content = script_content.replace('{{NMAX}}', str(nmax_value))
        script_content = script_content.replace('{{CORES}}', str(self.config['cores']))
        script_content = script_content.replace('{{ITERATION}}', str(iteration))
        script_content = script_content.replace('{{MONTH_NUMBER}}', str(spec['month_number']))
        # {{COVARIATE_CONFIGS}} is inserted into an R double-quoted string:
        #   covariate_configs <- "{{COVARIATE_CONFIGS}}"
        # The JSON itself contains double quotes, which would terminate that
        # string early ("unexpected symbol" in R). Escape backslashes then
        # double quotes so the JSON survives as a valid R string literal that
        # fromJSON() can parse.
        covariate_configs_r = covariate_configs.replace('\\', '\\\\').replace('"', '\\"')
        script_content = script_content.replace('{{COVARIATE_CONFIGS}}', covariate_configs_r)
        script_content = script_content.replace('{{GMTED_FILE}}', gmted_file)
        script_content = script_content.replace('{{KRIGING_EQUATION}}', self.config['kriging_equation'])
        # Coverage gate parameters (R reports the numbers; the driver enforces).
        script_content = script_content.replace('{{MIN_COVERAGE}}', str(self.min_coverage))
        script_content = script_content.replace('{{COVERAGE_DENOMINATOR}}', self.coverage_denominator)
        # Result-scan physical bounds.
        script_content = script_content.replace('{{VALID_MIN}}', str(self.valid_min))
        script_content = script_content.replace('{{VALID_MAX}}', str(self.valid_max))

        # Spatial extent
        script_content = script_content.replace('{{EXTENT_LONMIN}}', str(self.config['extent_lonmin']))
        script_content = script_content.replace('{{EXTENT_LONMAX}}', str(self.config['extent_lonmax']))
        script_content = script_content.replace('{{EXTENT_LATMIN}}', str(self.config['extent_latmin']))
        script_content = script_content.replace('{{EXTENT_LATMAX}}', str(self.config['extent_latmax']))
        script_content = script_content.replace('{{SHAPEFILE_PATH}}', str(self.config.get('shapefile_path', 'NULL')))
        script_content = script_content.replace('{{WORKING_DIR}}', str(self.working_dir))

        # Write script to node-local temp (see local_tmp note in __init__).
        script_file = self.local_tmp / f"multi_cov_iter{iteration}_{spec['model']}_{spec['date'].replace('-', '')}.R"
        with open(script_file, 'w') as f:
            f.write(script_content)

        return script_file

    def run_r_script(self, script_file):
        """Execute an R script. Returns the CompletedProcess (never a bool).

        On failure, logs the return code plus BOTH stderr and the tail of
        stdout, so a killed or silent R process is never logged blank.
        """
        # Optional wall-clock guard: 'timeout -k 60 <N>' sends SIGTERM at N
        # seconds, then SIGKILL 60s later, so a compute hang cannot block a
        # worker forever. Return code 124 means the guard fired.
        r_cmd = f"apptainer exec {self.config['apptainer_image']} Rscript {script_file}"
        if self.r_timeout_seconds > 0:
            r_cmd = f"timeout -k 60 {self.r_timeout_seconds} {r_cmd}"

        cmd = [
            'bash', '-c',
            f"module load {' '.join(self.config['required_modules'])} && {r_cmd}"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 124:
            self.logger.error(
                f"R script TIMED OUT after {self.r_timeout_seconds}s (killed by guard); "
                f"treating as failure")
        if result.returncode != 0:
            self.logger.error(f"R script failed (returncode {result.returncode})")
            if result.stderr and result.stderr.strip():
                self.logger.error("  stderr:\n" + result.stderr.strip())
            if result.stdout and result.stdout.strip():
                tail = "\n".join(result.stdout.strip().splitlines()[-25:])
                self.logger.error("  stdout (tail):\n" + tail)

        return result

    def parse_coverage(self, stdout):
        """Parse the coverage tokens the R template prints.

        Returns (passed, cov_domain, cov_grid, cov_selected):
          passed       : bool from COVERAGE_PASS, or None if not present
          cov_domain   : fraction over the expected-data domain, or None
          cov_grid     : fraction over the whole grid, or None
          cov_selected : the fraction for the configured denominator, or None
        """
        if not stdout:
            return None, None, None, None

        def num(key):
            m = re.search(rf"{key}=([0-9]*\.?[0-9]+)", stdout)
            return float(m.group(1)) if m else None

        cov_domain = num("COVERAGE_DOMAIN")
        cov_grid = num("COVERAGE_GRID")
        cov_selected = cov_grid if self.coverage_denominator == 'grid' else cov_domain

        m = re.search(r"COVERAGE_PASS=(TRUE|FALSE)", stdout)
        passed = (m.group(1) == "TRUE") if m else None

        # Fallback: if R did not emit COVERAGE_PASS but did emit numbers,
        # decide here using the driver's own threshold.
        if passed is None and cov_selected is not None:
            passed = (cov_selected + 1e-9) >= self.min_coverage

        return passed, cov_domain, cov_grid, cov_selected

    def delete_outputs(self, files, reason):
        """Delete the given output files if they exist (used on gate failure)."""
        for f in files:
            try:
                p = Path(f)
                if p.exists():
                    p.unlink()
                    self.logger.info(f"  Deleted ({reason}): {p.name}")
            except OSError as exc:
                self.logger.warning(f"  Could not delete {f}: {exc}")

    def preserve_script(self, script_file, tag):
        """Copy a failing R script from node-local temp to $SCRATCH temp_dir.

        Node-local temp is wiped at process end, so on a failure we keep the
        exact generated script on $SCRATCH for inspection.
        """
        try:
            src = Path(script_file)
            if src.exists() and Path(self.temp_dir) != Path(self.local_tmp):
                dst = self.temp_dir / f"FAILED_{tag}_{src.name}"
                shutil.copy2(src, dst)
                self.logger.info(f"  Saved failing R script for inspection: {dst}")
        except OSError as exc:
            self.logger.warning(f"  Could not preserve R script {script_file}: {exc}")

    def cleanup_local_tmp(self):
        """Remove the node-local temp dir (intermediate files) at process end."""
        try:
            if Path(self.local_tmp) != Path(self.temp_dir) and Path(self.local_tmp).exists():
                shutil.rmtree(self.local_tmp, ignore_errors=True)
                self.logger.info(f"Cleaned node-local temp dir: {self.local_tmp}")
        except OSError as exc:
            self.logger.warning(f"Could not clean node-local temp dir: {exc}")

    def check_coverage_gate(self, result, iteration):
        """Given a finished R run, return True if this downscale passes the gate.

        A missing/failed R run, or coverage below the threshold, returns False.
        """
        if result is None or result.returncode != 0:
            self.logger.error(f"Iteration {iteration} failed (R error)")
            return False

        passed, cov_domain, cov_grid, cov_selected = self.parse_coverage(result.stdout)

        # Surface the cascade's own log (method chosen + per-attempt scans) even
        # on success, so the kriging method used is visible in the run log.
        if result.stdout:
            for line in result.stdout.splitlines():
                s = line.strip()
                if ("Method used:" in s) or ("scan[" in s) or ("Adaptive:" in s):
                    self.logger.info(f"Iteration {iteration} R: {s}")

        def pct(x):
            return "n/a" if x is None else f"{100 * x:.2f}%"

        self.logger.info(
            f"Iteration {iteration} coverage: "
            f"selected[{self.coverage_denominator}]={pct(cov_selected)} "
            f"(domain={pct(cov_domain)}, grid={pct(cov_grid)}); "
            f"threshold={100 * self.min_coverage:.2f}%")

        if passed is None:
            # No coverage numbers at all despite a zero exit code: cannot verify.
            self.logger.error(
                f"Iteration {iteration}: no coverage info in R output; "
                f"treating as FAILED coverage gate")
            return False

        if not passed:
            nan_txt = "n/a" if cov_selected is None else f"{100 * (1 - cov_selected):.2f}%"
            self.logger.error(
                f"Iteration {iteration} FAILED coverage gate "
                f"(coverage {pct(cov_selected)} < {100 * self.min_coverage:.2f}%, "
                f"NaN {nan_txt})")
            return False

        self.logger.info(f"Iteration {iteration} PASSED coverage gate")
        return True

    def process_spec(self, spec):
        """Process one timestep with two-step downscaling."""

        model = spec['model']
        variable = spec['variable']
        scenario = spec['scenario']
        date = spec['date']
        month = spec['month_number']

        self.logger.info("="*70)
        self.logger.info(f"Processing: {model} - {variable} - {date} (Month {month})")
        self.logger.info("="*70)

        # Construct input filename as a glob. The realization/grid member varies
        # by model, so use r*_g* rather than a hardcoded r1i1p1f1_gn. The list may
        # use the short scenario "hist" while filenames use "historical", so glob
        # hist* for those; other scenarios (ssp*) match exactly.
        if scenario in ('hist', 'historical'):
            scen_glob = 'hist*'
        else:
            scen_glob = scenario
        input_pattern = f"{variable}_day_{model}_{scen_glob}_r*_g*_all.nc"
        input_matches = sorted(Path(self.config['input_dir']).glob(input_pattern))

        if not input_matches:
            self.logger.error(
                f"Input file not found: {Path(self.config['input_dir']) / input_pattern}")
            return False
        if len(input_matches) > 1:
            self.logger.warning(
                f"Multiple input files match {input_pattern}; using first of {len(input_matches)}:")
            for m in input_matches:
                self.logger.warning(f"    - {m.name}")
        input_file = input_matches[0]
        self.logger.info(f"Input file: {input_file.name}")

        # Extract the specific timestep to node-local temp (not $SCRATCH).
        extracted_file = self.local_tmp / f"{model}_{scenario}_{variable}_{date.replace('-', '')}_extracted.nc"

        if not extracted_file.exists():
            self.logger.info(f"Extracting timestep for {date}...")
            if not self.extract_timestep(input_file, date, extracted_file):
                self.logger.error(f"Failed to extract timestep for {date}")
                return False
        else:
            self.logger.info(f"Using existing extracted file: {extracted_file.name}")

        scenario_short = 'hist' if scenario == 'historical' else scenario[:4]
        date_str = date.replace('-', '')

        # All four outputs for this timestep. On a coverage-gate failure at any
        # iteration, ALL of these that exist are deleted (delete both 9km and
        # 1km and their SE files), and the timestep is marked failed.
        target_res_1 = self.config['target_res_first']
        target_res_2 = self.config['target_res_second']
        output_9km = self.output_dir / f"{model}_{scenario_short}_{variable}_9km_{date_str}.nc"
        se_9km = self.output_dir / f"SE_{model}_{scenario_short}_{variable}_9km_{date_str}.nc"
        output_1km = self.output_dir / f"{model}_{scenario_short}_{variable}_1km_{date_str}.nc"
        se_1km = self.output_dir / f"SE_{model}_{scenario_short}_{variable}_1km_{date_str}.nc"
        all_outputs = [output_9km, se_9km, output_1km, se_1km]

        # ITERATION 1: coarse -> 9km
        self.logger.info(f"Iteration 1: Downscaling to 9km (using month {month} covariate)")

        if output_9km.exists():
            self.logger.info(
                f"9km output already exists: {output_9km.name} "
                f"(skipping downscale and coverage gate)")
        else:
            script_1 = self.create_r_script(spec, extracted_file, output_9km, se_9km,
                                           target_res_1, iteration=1)
            self.logger.info("Running iteration 1...")
            result_1 = self.run_r_script(script_1)

            # Coverage gate AFTER this downscale, BEFORE moving to iteration 2.
            if not self.check_coverage_gate(result_1, iteration=1):
                self.preserve_script(script_1, "iter1")
                self.delete_outputs(all_outputs, "iter1 gate/failure")
                self.logger.error(f"✗ Timestep {date} FAILED at iteration 1")
                return False
            self.logger.info(f"✓ Iteration 1 complete: {output_9km.name}")

        # ITERATION 2: 9km -> 1km
        self.logger.info(f"Iteration 2: Downscaling to 1km (using month {month} covariate)")

        if output_1km.exists():
            self.logger.info(
                f"1km output already exists: {output_1km.name} "
                f"(skipping downscale and coverage gate)")
        else:
            script_2 = self.create_r_script(spec, output_9km, output_1km, se_1km,
                                           target_res_2, iteration=2)
            self.logger.info("Running iteration 2...")
            result_2 = self.run_r_script(script_2)

            # Coverage gate AFTER this downscale. On failure, delete BOTH the
            # 9km and 1km outputs (and their SE files) for this timestep.
            if not self.check_coverage_gate(result_2, iteration=2):
                self.preserve_script(script_2, "iter2")
                self.delete_outputs(all_outputs, "iter2 gate/failure")
                self.logger.error(f"✗ Timestep {date} FAILED at iteration 2")
                return False
            self.logger.info(f"✓ Iteration 2 complete: {output_1km.name}")

        self.logger.info(f"✓ Successfully completed {date} (both downscales passed the coverage gate)")
        return True

    def extract_timestep(self, input_file, date, output_file):
        """Extract a single timestep from input file using CDO."""

        # Use CDO to extract by date
        cmd = [
            'bash', '-c',
            f"module load cdo && cdo seldate,{date} {input_file} {output_file}"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            self.logger.error(f"CDO extraction failed: {result.stderr}")
            return False

        if output_file.exists() and output_file.stat().st_size > 0:
            self.logger.info(f"✓ Extracted timestep to: {output_file.name}")
            return True
        else:
            self.logger.error(f"Extraction produced empty or missing file")
            return False

    def run(self):
        """Run downscaling for all specs."""

        specs = self.load_extraction_list()
        self.logger.info(f"Loaded {len(specs)} monthly batches to process\n")

        successful = 0
        failed = 0

        try:
            for i, spec in enumerate(specs, 1):
                self.logger.info(f"\nProcessing {i}/{len(specs)}")

                if self.process_spec(spec):
                    successful += 1
                else:
                    failed += 1
        finally:
            # Always wipe node-local intermediates, even on error/interrupt.
            self.cleanup_local_tmp()

        # Summary
        self.logger.info("\n" + "="*70)
        self.logger.info("BATCH PROCESSING SUMMARY")
        self.logger.info("="*70)
        self.logger.info(f"Total timesteps: {len(specs)}")
        self.logger.info(f"Successful: {successful}")
        self.logger.info(f"Failed: {failed}")

        # Show date range
        if specs:
            dates = sorted([s['date'] for s in specs])
            self.logger.info(f"Date range: {dates[0]} to {dates[-1]}")

            # Show distribution by month
            from collections import Counter
            month_counts = Counter([s['month_number'] for s in specs])
            self.logger.info("\nDistribution by month:")
            for month in sorted(month_counts.keys()):
                month_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'][month-1]
                self.logger.info(f"  {month_name}: {month_counts[month]} timesteps")

        self.logger.info("="*70)

        return failed == 0

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 custom_covariate_hurs_downscaling.py <config.json>")
        print("\nExample:")
        print("  python3 custom_covariate_hurs_downscaling.py cmip6_config_rsds_custom.json")
        sys.exit(1)

    downscaler = CustomCovariateDownscaler(sys.argv[1])
    success = downscaler.run()
    sys.exit(0 if success else 1)
