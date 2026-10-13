# utils

Two helper scripts for the downscaling stage. They find the 1 km outputs that
need (re)producing and feed the gaps back into reprocessing.

- `check_missing_1km_files.sh` finds timesteps that were never produced.
- `check_1km_coverage.sh` finds timesteps that were produced but hold too few
  valid pixels, and can delete them so they get re-run.

Both write the same four-column reprocessing list (MODEL, VARIABLE, SCENARIO,
DATE) that the setup scripts (`../setup_*.sh`) read back, so a check feeds
straight into a re-run.

## check_missing_1km_files.sh

Compares a master list of expected files against what is actually in an output
folder, and writes out the entries that still need to be produced. An entry
counts as present only when both of its files are there, the variable file and
its standard-error (SE) partner. If either one is missing, or both, the entry
goes on the missing list.

    ./check_missing_1km_files.sh <master_list> <output_dir> [output_file]

Arguments:

- `master_list`: the expected entries, four columns (MODEL, VARIABLE, SCENARIO,
  DATE), for example `tas_by_model.txt` from `metadata/`.
- `output_dir`: the folder that holds the actual `*_1km_*.nc` files.
- `output_file`: optional. Where to write the missing entries. Defaults to
  `missing_files.txt`.

What it handles:

- Both scenario spellings. It matches the short tokens (`ssp1`, `ssp2`, `ssp3`,
  `ssp5`) and the full names (`ssp126`, `ssp245`, `ssp370`, `ssp585`), and reads
  `historical` as `hist` in the filename.
- Both date forms in filenames, compact `YYYYMMDD` and dashed `YYYY-MM-DD`.
- The SE files. It indexes `SE_*_1km_*.nc` separately and needs both the
  variable file and its SE partner before it treats an entry as done.

Output:

- The missing list keeps the clean four-column format, so it feeds straight into
  the setup scripts as a reprocessing list.
- The script prints a summary: the totals found and missing, and a breakdown of
  how many entries were missing the variable file only, the SE file only, or
  both.
- If nothing is missing, it reports that every entry is complete and removes the
  empty output file.

Example:

    ./check_missing_1km_files.sh tas_by_model.txt /PATH/TO/WORKING_DIR/downscaled/tas missing_files_tas.txt

## check_1km_coverage.sh

Scans a directory of finished 1 km NetCDF files and flags any whose valid-pixel
coverage falls below a threshold (default 89 %). For each flagged file it can
delete the whole timestep set so the timestep can be re-run. It needs `cdo`.

    ./check_1km_coverage.sh [-t THRESHOLD] [-o OUTFILE] [-l LISTFILE] [-D] [-e] [-r] [DIR]

Options:

- `-t THRESHOLD`: minimum acceptable coverage percent (default 89).
- `-o OUTFILE`: write the flagged 1 km file paths, one per line.
- `-l LISTFILE`: write the flagged files as a four-column reprocessing list
  (MODEL, VARIABLE, SCENARIO, DATE), which the setup scripts read directly.
- `-D`: delete the companion set of each flagged file. The default is a dry run
  that only prints what would be deleted.
- `-e`: exit non-zero if anything was flagged (for CI). The default is exit 0,
  so the scan is safe inside `cmd && cmd` chains.
- `-r`: recurse into sub-directories.
- `DIR`: the directory to scan (default: the current directory).

How it works:

- It considers files matching `*1km*.nc` and skips the `SE_` standard-error
  files. Coverage is `100 * (gridsize - missing) / gridsize`, read from
  `cdo info` (missing = NaN or `_FillValue` cells). A file `cdo` cannot read is
  flagged.
- For a flagged file `MODEL_SCEN_VAR_1km_YYYYMMDD.nc`, the deletion set is the
  four companion files (`1km`, `SE_1km`, `9km`, `SE_9km`) for that same model,
  scenario, variable, and date. It uses exact filenames, never a wildcard, and
  deletes only with `-D`.
- When it writes the reprocessing list, it maps the short scenario token back to
  the full name (`hist` to `historical`, `ssp1` to `ssp126`, and so on).

Examples:

    ./check_1km_coverage.sh /PATH/TO/WORKING_DIR/downscaled/hurs                 # scan only
    ./check_1km_coverage.sh -l rerun.txt /PATH/TO/WORKING_DIR/downscaled/hurs    # scan + re-run list (review the dry run)
    ./check_1km_coverage.sh -l rerun.txt -D /PATH/TO/WORKING_DIR/downscaled/hurs # scan + delete the flagged sets

## The two together

`check_missing_1km_files.sh` catches what was never made; `check_1km_coverage.sh`
catches what was made but came out too sparse. Run both, concatenate their
four-column lists, and feed the result to the matching `setup_*.sh` to re-run the
gaps.
