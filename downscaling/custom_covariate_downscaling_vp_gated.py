#!/usr/bin/env python3
"""
Two-step kriging driver (28km -> 9km -> 1km) for the VP family: huss (covariate
= vapor pressure), rsds and rlds (covariate = srad), monthly .nc. Pairs with
kriging_template_adaptive_vp_gated.R.

Behaviour:
  1. Coverage gate - after each downscale (9km, then 1km) coverage must be
     >= min_coverage of the area (NaN <= 1 - min_coverage). On failure all
     companion files of the timestep (9km, SE_9km, 1km, SE_1km, under either
     scenario spelling) are deleted and the timestep is marked failed.
     Default gate: 0.9988 on the "domain" denominator (finite-covariate cells).
  2. Intermediate files are written to node-local $TMPDIR.
  3. Per-iteration wall-clock timeout guard on the R call.
  4. run_r_script returns the CompletedProcess and logs return code + stderr +
     stdout tail on failure.
  5. A failing R script is copied to $SCRATCH for inspection.
  6. Output files are recognized and deleted under either scenario spelling
     (historical<->hist, ssp370<->ssp3).

Usage:
    python3 custom_covariate_downscaling_vp_gated.py <config.json>
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


class VpAdaptiveDownscalerGated:
    """Two-step adaptive-covariate downscaling with a production coverage gate."""

    def __init__(self, config_file):
        with open(config_file, 'r') as f:
            self.config = json.load(f)

        self.working_dir = Path(self.config['working_dir'])
        self.output_dir = Path(self.config['output_dir'])
        self.temp_dir = Path(self.config['temp_dir'])

        # Production coverage gate (config-driven). Default 0.9988 / domain.
        self.min_coverage = float(self.config.get('min_coverage', 0.9988))
        self.coverage_denominator = str(
            self.config.get('coverage_denominator', 'domain')).lower()
        if self.coverage_denominator not in ('domain', 'grid'):
            self.coverage_denominator = 'domain'

        # Per-iteration wall-clock guard (seconds). 0/negative disables it.
        self.r_timeout_seconds = int(self.config.get('r_timeout_seconds', 3600))

        # Node-local scratch for INTERMEDIATE files.
        tmp_base = Path(os.environ.get('TMPDIR', '/tmp'))
        var_tag = str(self.config.get('variable', 'redown'))
        self.local_tmp = tmp_base / f"redown_{var_tag}_{os.getpid()}"

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        try:
            self.local_tmp.mkdir(parents=True, exist_ok=True)
        except OSError:
            self.local_tmp = self.temp_dir

        log_file = self.working_dir / f"gated_downscaling_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Node-local temp dir: {self.local_tmp}")
        self.logger.info(
            f"Coverage gate: min_coverage={self.min_coverage} "
            f"denominator={self.coverage_denominator}")
        if self.r_timeout_seconds > 0:
            self.logger.info(f"R per-iteration timeout: {self.r_timeout_seconds}s")
        else:
            self.logger.info("R per-iteration timeout: disabled")

    @staticmethod
    def scenario_tokens(scenario):
        """Scenario spellings that may appear in an output filename.

        Globbing over BOTH spellings lets a re-run find/delete a prior file
        however it was written: historical<->hist, ssp370<->ssp3, etc.
        """
        if scenario in ('hist', 'historical'):
            return ['historical', 'hist']
        toks = [scenario]
        # full ssp name -> also its 4-char truncation (ssp370 -> ssp3)
        if scenario.startswith('ssp') and len(scenario) > 4 and scenario[:4] not in toks:
            toks.append(scenario[:4])
        return toks

    def load_extraction_list(self):
        """Load extraction list (MODEL, VARIABLE, SCENARIO, DATE)."""
        extraction_file = self.working_dir / self.config['extraction_list']
        specs = []
        self.logger.info(f"Loading extraction list: {extraction_file}")
        with open(extraction_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 4:
                    model, variable, scenario, date_str = parts[0], parts[1], parts[2], parts[3]
                    try:
                        year, month, day = date_str.split('-')
                        specs.append({
                            'model': model, 'variable': variable, 'scenario': scenario,
                            'date': date_str, 'year': int(year),
                            'month_number': int(month), 'day': int(day)
                        })
                    except ValueError as e:
                        self.logger.warning(f"Line {line_num}: Could not parse date '{date_str}': {e}")
                        continue
        self.logger.info(f"Loaded {len(specs)} timesteps")
        return specs

    def create_r_script(self, spec, input_file, output_file, se_file, target_res, iteration):
        """Create R script for one iteration."""
        template_file = self.working_dir / self.config['r_script_template']
        with open(template_file, 'r') as f:
            template = f.read()

        if iteration == 1:
            gmted_file = self.config.get('gmted_iter1', self.config['gmted_file'])
        else:
            gmted_file = self.config.get('gmted_iter2', self.config['gmted_file'])

        # Covariate config. Supports single (adaptive template: covariate_base_path)
        # and multi (covariates list) shapes; the template consumes whichever it has.
        covariate_configs = None
        if 'covariates' in self.config:
            covariates = []
            for cov in self.config['covariates']:
                cov_copy = cov.copy()
                if 'path' in cov_copy:
                    scen = 'hist*' if spec['scenario'] in ('hist', 'historical') else spec['scenario']
                    cov_copy['path'] = (cov_copy['path']
                        .replace('{MODEL}', spec['model']).replace('{model}', spec['model'])
                        .replace('{SCENARIO}', scen).replace('{scenario}', scen))
                covariates.append(cov_copy)
            covariate_configs = json.dumps(covariates)
        elif 'covariate_base_path' in self.config:
            single_cov = [{'name': 'covariate', 'type': 'monthly_file',
                           'path': self.config['covariate_base_path']}]
            covariate_configs = json.dumps(single_cov)
        else:
            self.logger.error("No covariate configuration found")
            return None

        s = template.replace('{{INPUT_FILE}}', str(input_file))
        s = s.replace('{{OUTPUT_FILE}}', str(output_file))
        s = s.replace('{{SE_FILE}}', str(se_file))
        s = s.replace('{{VARIABLE}}', spec['variable'])
        s = s.replace('{{DATE}}', spec['date'])
        s = s.replace('{{TARGET_RESOLUTION}}', str(target_res))
        s = s.replace('{{TIMESTEP_INDEX}}', 'NULL')
        nmax_value = self.config.get('nmax_base', self.config.get('nmax'))
        if nmax_value is None:
            self.logger.error("Config has neither 'nmax_base' nor 'nmax'")
            return None
        s = s.replace('{{NMAX}}', str(nmax_value))
        s = s.replace('{{CORES}}', str(self.config['cores']))
        s = s.replace('{{ITERATION}}', str(iteration))
        s = s.replace('{{MONTH_NUMBER}}', str(spec['month_number']))
        covariate_configs_r = covariate_configs.replace('\\', '\\\\').replace('"', '\\"')
        s = s.replace('{{COVARIATE_CONFIGS}}', covariate_configs_r)
        s = s.replace('{{COVARIATE_BASE_PATH}}', str(self.config.get('covariate_base_path', 'NULL')))
        s = s.replace('{{GMTED_FILE}}', gmted_file)
        s = s.replace('{{KRIGING_EQUATION}}', self.config['kriging_equation'])
        s = s.replace('{{MIN_COVERAGE}}', str(self.min_coverage))
        s = s.replace('{{COVERAGE_DENOMINATOR}}', self.coverage_denominator)
        s = s.replace('{{EXTENT_LONMIN}}', str(self.config['extent_lonmin']))
        s = s.replace('{{EXTENT_LONMAX}}', str(self.config['extent_lonmax']))
        s = s.replace('{{EXTENT_LATMIN}}', str(self.config['extent_latmin']))
        s = s.replace('{{EXTENT_LATMAX}}', str(self.config['extent_latmax']))
        s = s.replace('{{SHAPEFILE_PATH}}', str(self.config.get('shapefile_path', 'NULL')))
        s = s.replace('{{WORKING_DIR}}', str(self.working_dir))

        script_file = self.local_tmp / f"adaptive_iter{iteration}_{spec['model']}_{spec['date'].replace('-', '')}.R"
        with open(script_file, 'w') as f:
            f.write(s)
        return script_file

    def run_r_script(self, script_file):
        """Execute an R script. Returns the CompletedProcess (never a bool)."""
        r_cmd = f"apptainer exec {self.config['apptainer_image']} Rscript {script_file}"
        if self.r_timeout_seconds > 0:
            r_cmd = f"timeout -k 60 {self.r_timeout_seconds} {r_cmd}"
        cmd = ['bash', '-c',
               f"module load {' '.join(self.config['required_modules'])} && {r_cmd}"]
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
        """Parse the COVERAGE_* tokens the gated R template prints."""
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
        """Copy a failing R script from node-local temp to $SCRATCH temp_dir."""
        try:
            src = Path(script_file)
            if src.exists() and Path(self.temp_dir) != Path(self.local_tmp):
                dst = self.temp_dir / f"FAILED_{tag}_{src.name}"
                shutil.copy2(src, dst)
                self.logger.info(f"  Saved failing R script for inspection: {dst}")
        except OSError as exc:
            self.logger.warning(f"  Could not preserve R script {script_file}: {exc}")

    def cleanup_local_tmp(self):
        """Remove the node-local temp dir at process end."""
        try:
            if Path(self.local_tmp) != Path(self.temp_dir) and Path(self.local_tmp).exists():
                shutil.rmtree(self.local_tmp, ignore_errors=True)
                self.logger.info(f"Cleaned node-local temp dir: {self.local_tmp}")
        except OSError as exc:
            self.logger.warning(f"Could not clean node-local temp dir: {exc}")

    def check_coverage_gate(self, result, iteration):
        """Return True if this downscale passes the coverage gate, else False."""
        if result is None or result.returncode != 0:
            self.logger.error(f"Iteration {iteration} failed (R error)")
            return False
        passed, cov_domain, cov_grid, cov_selected = self.parse_coverage(result.stdout)
        if result.stdout:
            for line in result.stdout.splitlines():
                s = line.strip()
                if ("Strategy used:" in s) or ("STRATEGY" in s) or ("scan[" in s) \
                        or ("Adaptive parameters:" in s) or ("PASSED" in s) or ("FAILED" in s):
                    self.logger.info(f"Iteration {iteration} R: {s}")

        def pct(x):
            return "n/a" if x is None else f"{100 * x:.2f}%"

        self.logger.info(
            f"Iteration {iteration} coverage: "
            f"selected[{self.coverage_denominator}]={pct(cov_selected)} "
            f"(domain={pct(cov_domain)}, grid={pct(cov_grid)}); "
            f"threshold={100 * self.min_coverage:.2f}%")
        if passed is None:
            self.logger.error(
                f"Iteration {iteration}: no coverage info in R output; treating as FAILED")
            return False
        if not passed:
            nan_txt = "n/a" if cov_selected is None else f"{100 * (1 - cov_selected):.2f}%"
            self.logger.error(
                f"Iteration {iteration} FAILED coverage gate "
                f"(coverage {pct(cov_selected)} < {100 * self.min_coverage:.2f}%, NaN {nan_txt})")
            return False
        self.logger.info(f"Iteration {iteration} PASSED coverage gate")
        return True

    def process_spec(self, spec):
        """Process one timestep with two-step downscaling + coverage gate."""
        model = spec['model']
        variable = spec['variable']
        scenario = spec['scenario']
        date = spec['date']
        month = spec['month_number']

        self.logger.info("=" * 70)
        self.logger.info(f"Processing: {model} - {variable} - {date} (Month {month})")
        self.logger.info("=" * 70)

        # Locate the source .nc (flexible glob; hist* for historical).
        scen_glob = 'hist*' if scenario in ('hist', 'historical') else scenario
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

        extracted_file = self.local_tmp / f"{model}_{scenario}_{variable}_{date.replace('-', '')}_extracted.nc"
        if not extracted_file.exists():
            self.logger.info(f"Extracting timestep for {date}...")
            if not self.extract_timestep(input_file, date, extracted_file):
                self.logger.error(f"Failed to extract timestep for {date}")
                return False
        else:
            self.logger.info(f"Using existing extracted file: {extracted_file.name}")

        date_str = date.replace('-', '')
        # Canonical (short) scenario token used when WRITING new files.
        scenario_short = 'hist' if scenario == 'historical' else scenario[:4]
        # All scenario spellings recognized on disk (scenario globbing).
        scen_tokens = self.scenario_tokens(scenario)

        # Canonical write targets (short spelling).
        output_9km = self.output_dir / f"{model}_{scenario_short}_{variable}_9km_{date_str}.nc"
        se_9km = self.output_dir / f"SE_{model}_{scenario_short}_{variable}_9km_{date_str}.nc"
        output_1km = self.output_dir / f"{model}_{scenario_short}_{variable}_1km_{date_str}.nc"
        se_1km = self.output_dir / f"SE_{model}_{scenario_short}_{variable}_1km_{date_str}.nc"

        # Every spelling of every companion, for existence checks and deletion.
        def companions(res):
            out = []
            for tok in scen_tokens:
                out.append(self.output_dir / f"{model}_{tok}_{variable}_{res}_{date_str}.nc")
                out.append(self.output_dir / f"SE_{model}_{tok}_{variable}_{res}_{date_str}.nc")
            return out

        def existing_data(res):
            for p in companions(res):
                if not p.name.startswith('SE_') and p.exists():
                    return p
            return None

        all_outputs = companions('9km') + companions('1km')

        # ITERATION 1: coarse -> 9km
        self.logger.info(f"Iteration 1: Downscaling to 9km (using month {month} covariate)")
        existing_9 = existing_data('9km')
        if existing_9 is not None:
            self.logger.info(
                f"9km output already exists: {existing_9.name} (skipping downscale and gate)")
        else:
            script_1 = self.create_r_script(spec, extracted_file, output_9km, se_9km,
                                            self.config['target_res_first'], iteration=1)
            self.logger.info("Running iteration 1...")
            result_1 = self.run_r_script(script_1)
            if not self.check_coverage_gate(result_1, iteration=1):
                self.preserve_script(script_1, "iter1")
                self.delete_outputs(all_outputs, "iter1 gate/failure")
                self.logger.error(f"Timestep {date} FAILED at iteration 1")
                return False
            self.logger.info(f"Iteration 1 complete: {output_9km.name}")
            existing_9 = output_9km

        # ITERATION 2: 9km -> 1km (input is the 9km file found or just written)
        self.logger.info(f"Iteration 2: Downscaling to 1km (using month {month} covariate)")
        existing_1 = existing_data('1km')
        if existing_1 is not None:
            self.logger.info(
                f"1km output already exists: {existing_1.name} (skipping downscale and gate)")
        else:
            input_9km = existing_9 if existing_9 is not None else output_9km
            script_2 = self.create_r_script(spec, input_9km, output_1km, se_1km,
                                            self.config['target_res_second'], iteration=2)
            self.logger.info("Running iteration 2...")
            result_2 = self.run_r_script(script_2)
            if not self.check_coverage_gate(result_2, iteration=2):
                self.preserve_script(script_2, "iter2")
                self.delete_outputs(all_outputs, "iter2 gate/failure")
                self.logger.error(f"Timestep {date} FAILED at iteration 2")
                return False
            self.logger.info(f"Iteration 2 complete: {output_1km.name}")

        self.logger.info(f"Successfully completed {date} (both downscales passed the coverage gate)")
        return True

    def extract_timestep(self, input_file, date, output_file):
        """Extract a single timestep from input file using CDO."""
        cmd = ['bash', '-c', f"module load cdo && cdo seldate,{date} {input_file} {output_file}"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self.logger.error(f"CDO extraction failed: {result.stderr}")
            return False
        if output_file.exists() and output_file.stat().st_size > 0:
            self.logger.info(f"Extracted timestep to: {output_file.name}")
            return True
        self.logger.error("Extraction produced empty or missing file")
        return False

    def run(self):
        """Run downscaling for all specs."""
        specs = self.load_extraction_list()
        self.logger.info(f"Loaded {len(specs)} timesteps to process\n")
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
            self.cleanup_local_tmp()

        self.logger.info("\n" + "=" * 70)
        self.logger.info("BATCH PROCESSING SUMMARY")
        self.logger.info("=" * 70)
        self.logger.info(f"Total timesteps: {len(specs)}")
        self.logger.info(f"Successful: {successful}")
        self.logger.info(f"Failed: {failed}")
        if specs:
            dates = sorted([s['date'] for s in specs])
            self.logger.info(f"Date range: {dates[0]} to {dates[-1]}")
        self.logger.info("=" * 70)
        return failed == 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 custom_covariate_downscaling_vp_gated.py <config.json>")
        print("Example: python3 custom_covariate_downscaling_vp_gated.py cmip6_config_huss_gated.json")
        sys.exit(1)
    downscaler = VpAdaptiveDownscalerGated(sys.argv[1])
    success = downscaler.run()
    sys.exit(0 if success else 1)
