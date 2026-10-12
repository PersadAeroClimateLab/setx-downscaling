#!/usr/bin/env python3
"""
compare_nc_zarr_diagnostics.py  --  NetCDF-vs-Zarr QC comparison for SETx-1km-CMIP6.

For every (variable, model, scenario) it pairs the published NetCDF year files
with the matching Zarr store and runs the SAME full-record diagnostics on both
sources, then compares them:

  A. INVENTORY   pairs found only in NC, only in Zarr, or in both.
  B. METADATA    calendar, units, grid shape, lat/lon values, lat/lon orientation,
                 SE variable present in each source.
  C. TIME AXIS   missing days (per source, against its own calendar), duplicate
                 dates, dates in NC but not in Zarr, dates in Zarr but not in NC,
                 NC file-name date range vs file content.
  D. CELL-BY-CELL (every common day, every cell -- no sampling)
                 NaN-footprint mismatch (NaN in one source, valid in the other),
                 value differences beyond tolerance, exact-equal cells, max and
                 mean absolute difference, and how many differing cells are
                 outside the physical bounds in each source (i.e. explained by a
                 QC clip). The same for the SE field.
  E. DIAGNOSTICS (computed separately on each source, then compared date by date)
                 EMPTY / PARTIAL / EXTRA-valid days (median NaN baseline),
                 NEG, BELOW_MIN, ABOVE_MAX, INF, FILL_LEAK (undeclared -3.4e38),
                 CONSTANT, BLOWUP (roughness >= 6x median), OUTLIER (|z| > 5),
                 SE_NEG, SE_ABOVE_CEIL, SE_COV_ANOM.
                 Each flagged date is labelled BOTH / NC_ONLY / ZARR_ONLY.

NC tree (walked recursively, parsed from the file name):
  <nc-root>/CMIP/<scenario>/<MODEL>/r1i1p1f1/day/<var>/v20260720/
      <var>_NAM_SETx-UTAustin_KrigR_CMIP6_<scenario>_<MODEL>_<member>_day_<YYYYMMDD>-<YYYYMMDD>.nc
  One file per year. Each file holds the field <var> AND its standard error SE_<var>.
  (The Zarr stores name the SE <var>_SE; both forms are paired automatically.)
  Physical bounds are the most recent set, from qc_zarr_all.py.
  When one file period exists in more than one version folder, the latest
  version is used (the count of superseded files is reported).

Zarr tree:
  <zarr-root>/<var>/<MODEL>_<hist|sspNNN>_<var>_<Y0>-<Y1>.zarr   (all 8 variables, incl. sfcWind)

Outputs (in --out):
  compare_summary.csv          one row per pair (resumable)
  <pair>_diff_dates.csv        every date with a time, footprint or value difference
  <pair>_flag_compare.csv      every date flagged by a diagnostic in either source
  <pair>_nc_files.csv          per-NC-file audit (name range vs content, duplicates)
  inventory.csv                all pairs and which side has them
  SUMMARY.txt                  roll-up

Run:
  python3 compare_nc_zarr_diagnostics.py --limit 1                     # smoke test
  nohup python3 -u compare_nc_zarr_diagnostics.py > cmp.log 2>&1 &     # full run
  python3 compare_nc_zarr_diagnostics.py --vars rlds --models BCC-CSM2-MR --scenarios ssp585
"""

import argparse
import csv
import os
import re
import time
import traceback
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore", category=RuntimeWarning)   # all-NaN reductions

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NC_ROOT_DEFAULT = "/Path/to/NetCDF/Files"
ZARR_ROOT_DEFAULT = "/Path/to/Zarr/Stores"

VARIABLES = ["hurs", "huss", "rlds", "rsds", "sfcWind", "tas", "tasmax", "tasmin"]
SCEN_ALIAS = {"hist": "historical", "historical": "historical"}

# Physical bounds (same as qc_zarr_all.py)
BOUNDS = {
    "tas": (180., 350.), "tasmax": (180., 350.), "tasmin": (180., 350.),
    "hurs": (0., 105.), "huss": (0., 0.06), "sfcWind": (0., 50.),
    "rlds": (40., 650.), "rsds": (0., 1500.),
}
# Calibrated SE ceilings (from se_distribution_diagnostic.py / qc_zarr_downscaled.py)
SE_CEILING = {
    "hurs": 3.51, "huss": 0.05, "rlds": 34.85, "rsds": 74.46,
    "sfcWind": 1.20, "tas": 4.65, "tasmax": 9.86, "tasmin": 4.68,
}

FILL_ABS = 1.0e30          # |x| >= this is a fill sentinel (-3.4e38, 1e20 etc. -> see note)
EMPTY_FRAC = 0.999         # >= this fraction NaN -> empty day
BLOWUP_FACTOR = 6.0        # roughness >= this x store median -> spatial blow-up
Z_OUTLIER = 5.0            # |z| of the domain mean -> outlier day
SE_HINTS = ("standard_error", "std_error", "stderror", "stderr", "se")

NC_RE = re.compile(
    r"^(?P<var>[A-Za-z]+)_NAM_SETx-UTAustin_KrigR_CMIP6_(?P<scen>[^_]+)_(?P<model>[^_]+)_"
    r"(?P<member>[^_]+)_day_(?P<start>\d{8})-(?P<end>\d{8})\.nc$")
ZARR_RE = re.compile(
    r"^(?P<model>.+?)_(?P<scen>[^_]+)_(?P<var>[^_]+)_(?P<y0>\d{4})-(?P<y1>\d{4})\.zarr$")
VERSION_RE = re.compile(r"^v\d{8}$")

CHECKS = ["EMPTY", "PARTIAL", "EXTRA", "NEG", "BELOW_MIN", "ABOVE_MAX", "INF",
          "FILL_LEAK", "CONSTANT", "BLOWUP", "OUTLIER",
          "SE_NEG", "SE_ABOVE_CEIL", "SE_INF", "SE_COV_ANOM"]

SUMMARY_FIELDS = [
    "pair", "var", "model", "scenario", "verdict", "reasons",
    # metadata
    "nc_files", "nc_superseded_files", "zarr_path",
    "calendar_nc", "calendar_zarr", "units_nc", "units_zarr",
    "grid_nc", "grid_zarr", "grid_match", "coord_max_diff", "lat_flipped",
    "se_var_nc", "se_var_zarr",
    # time
    "n_time_nc", "n_time_zarr", "n_common", "dates_nc_only", "dates_zarr_only",
    "missing_nc", "missing_zarr", "dup_nc", "dup_zarr", "nc_file_range_issues",
    # cell-by-cell (main)
    "cells_compared", "cells_exact", "cells_diff", "cells_diff_nc_oob", "cells_diff_zarr_oob",
    "days_with_diff", "max_abs_diff", "mean_abs_diff",
    "nan_nc_only_cells", "nan_zarr_only_cells", "days_nan_mismatch",
    "dmean_max_abs_diff", "dmean_corr",
    # cell-by-cell (SE)
    "se_cells_compared", "se_cells_diff", "se_days_with_diff", "se_max_abs_diff",
    "se_nan_nc_only_cells", "se_nan_zarr_only_cells",
    # per-source diagnostic summaries
    "normal_nan_nc", "normal_nan_zarr", "min_nc", "min_zarr", "max_nc", "max_zarr",
    "mean_nc", "mean_zarr",
] + [f"{c}_{s}" for c in CHECKS for s in ("nc", "zarr", "ncOnly", "zarrOnly")] + ["seconds"]


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------
def canon_scen(s):
    return SCEN_ALIAS.get(s.lower(), s.lower())


def discover_nc(root):
    """{(var, model, scen): {"files": [paths sorted by start], "superseded": n}}"""
    groups = defaultdict(dict)       # key -> {(start, end): (version, path)}
    superseded = defaultdict(int)
    if not os.path.isdir(root):
        return {}
    for dp, _dns, fns in os.walk(root):
        for fn in fns:
            m = NC_RE.match(fn)
            if not m:
                continue
            d = m.groupdict()
            key = (d["var"], d["model"], canon_scen(d["scen"]))
            vdir = os.path.basename(dp)
            ver = vdir if VERSION_RE.match(vdir) else "v00000000"
            per = (d["start"], d["end"])
            old = groups[key].get(per)
            if old is None:
                groups[key][per] = (ver, os.path.join(dp, fn))
            else:
                superseded[key] += 1
                if ver > old[0]:
                    groups[key][per] = (ver, os.path.join(dp, fn))
    out = {}
    for key, per_map in groups.items():
        files = [per_map[p][1] for p in sorted(per_map)]
        out[key] = {"files": files, "superseded": superseded[key]}
    return out


def _find_zarr(root):
    out = []
    if not os.path.isdir(root):
        return out
    for dp, dns, _ in os.walk(root):
        for d in [x for x in dns if x.endswith(".zarr")]:
            out.append(os.path.join(dp, d))
        dns[:] = [x for x in dns if not x.endswith(".zarr")]
    return sorted(out)


def discover_zarr(root):
    """{(var, model, scen): path} for every store under the one Zarr root."""
    out = {}
    for p in _find_zarr(root):
        m = ZARR_RE.match(os.path.basename(p))
        if m:
            out[(m["var"], m["model"], canon_scen(m["scen"]))] = p
    return out


# ---------------------------------------------------------------------------
# Opening
# ---------------------------------------------------------------------------
def _time_coder():
    return xr.coders.CFDatetimeCoder(use_cftime=True)


def _open_zarr_manual(path):
    """Fallback for stores the xarray backend cannot open (mixed v2/v3 metadata).
    Builds the Dataset from zarr-python arrays using v3 dimension_names or the v2
    _ARRAY_DIMENSIONS attribute, then CF-decodes it."""
    import zarr
    import dask.array as dska
    g = zarr.open_group(path, mode="r")
    variables = {}
    for name, arr in g.arrays():
        attrs = dict(arr.attrs)
        dims = getattr(getattr(arr, "metadata", None), "dimension_names", None)
        dims = dims or attrs.pop("_ARRAY_DIMENSIONS", None)
        attrs.pop("_ARRAY_DIMENSIONS", None)
        if not dims:
            dims = [f"{name}_dim{i}" for i in range(arr.ndim)]
        data = arr[...] if arr.ndim <= 1 else dska.from_array(arr, chunks=arr.chunks)
        variables[name] = xr.Variable(tuple(dims), data, attrs)
    ds = xr.Dataset(variables, attrs=dict(g.attrs))
    return xr.decode_cf(ds, decode_times=_time_coder())


def open_zarr(path):
    try:
        return xr.open_zarr(path, decode_times=_time_coder(), consolidated=None, chunks=None)
    except Exception:
        return _open_zarr_manual(path)


def open_nc(path, engine):
    return xr.open_dataset(path, engine=engine, decode_times=_time_coder(), chunks=None,
                           cache=False)


def _dims(da):
    t = next((d for d in ("time", "Time", "t") if d in da.dims), da.dims[0])
    la = next((d for d in ("lat", "latitude", "y") if d in da.dims), da.dims[-2])
    lo = next((d for d in ("lon", "longitude", "x") if d in da.dims), da.dims[-1])
    return t, la, lo


def _is_se_partner(name, main):
    n, m = name.lower(), main.lower()
    if n == m:
        return False
    for h in SE_HINTS:
        if n.endswith(h) and n[: -len(h)].rstrip("_-") == m:
            return True
        if n.startswith(h) and n[len(h):].lstrip("_-") == m:
            return True
    return False


def pick_vars(ds, hint):
    d3 = [v for v in ds.data_vars if ds[v].ndim == 3]
    main = hint if hint in d3 else next((v for v in d3 if not any(
        _is_se_partner(v, x) for x in d3)), None)
    se = next((v for v in d3 if main and _is_se_partner(v, main)), None)
    return main, se


def date_strings(tvals):
    out = []
    for t in tvals:
        if hasattr(t, "strftime"):
            out.append(t.strftime("%Y-%m-%d"))
        else:
            out.append(str(np.datetime_as_string(np.datetime64(t), unit="D")))
    return np.array(out)


def calendar_of(ds, tdim):
    enc = ds[tdim].encoding.get("calendar")
    if enc:
        return str(enc)
    try:
        return str(ds[tdim].dt.calendar)
    except Exception:
        return ds[tdim].attrs.get("calendar", "unknown")


def expected_dates(d0, d1, cal):
    """All YYYY-MM-DD strings from d0 to d1 (inclusive) in calendar `cal`."""
    try:
        rng = xr.date_range(d0, d1, freq="D", calendar=cal, use_cftime=True)
        return [f"{t:%Y-%m-%d}" for t in rng]
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Per-block numerics (all float64)
# ---------------------------------------------------------------------------
def clean_block(raw):
    """raw -> (clean float64 with NaN for fill/inf, per-day fill_leak, per-day inf)."""
    b = np.array(raw, dtype=np.float64, copy=True)
    n = b.shape[0]
    big = np.abs(b) >= FILL_ABS           # fill sentinels AND +/-inf (NaN compares False)
    inf = np.isinf(b)
    fill = big & ~inf
    b[big] = np.nan
    return b, fill.reshape(n, -1).sum(1), inf.reshape(n, -1).sum(1)


def _rough1(b, axis):
    """Mean |adjacent difference| along `axis`, per day, NaN-aware, in place-friendly."""
    d = np.abs(np.diff(b, axis=axis))
    ok = np.isfinite(d)
    d[~ok] = 0.0
    n = d.shape[0]
    cnt = ok.reshape(n, -1).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, d.reshape(n, -1).sum(1) / cnt, np.nan)


def day_diag(b, lo, hi):
    """Per-day diagnostic vectors for a cleaned block (n, lat, lon)."""
    n = b.shape[0]
    flat = b.reshape(n, -1)
    ok = ~np.isnan(flat)
    cnt = ok.sum(1)
    z = np.where(ok, flat, 0.0)
    s1 = z.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(cnt > 0, s1 / cnt, np.nan)
        z -= np.where(cnt > 0, mean, 0.0)[:, None]
        z[~ok] = 0.0
        var = np.where(cnt > 0, (z * z).sum(1) / cnt, np.nan)
    out = {"nan": flat.shape[1] - cnt, "mean": mean, "std": np.sqrt(var),
           "min": np.fmin.reduce(flat, axis=1), "max": np.fmax.reduce(flat, axis=1),
           "neg": (flat < 0).sum(1)}
    out["bmin"] = (flat < lo).sum(1) if lo is not None else np.zeros(n, np.int64)
    out["bmax"] = (flat > hi).sum(1) if hi is not None else np.zeros(n, np.int64)
    out["rough"] = 0.5 * (_rough1(b, 1) + _rough1(b, 2))
    return out


def cell_compare(a, z, lo, hi, rtol, atol):
    """Cell-by-cell comparison of cleaned blocks a (NC) and z (Zarr), per day."""
    n = a.shape[0]
    r = lambda m: m.reshape(n, -1).sum(1)
    fa, fz = np.isfinite(a), np.isfinite(z)
    both = fa & fz
    d = np.abs(a - z)
    d[~both] = 0.0
    with np.errstate(invalid="ignore"):
        diff = d > (atol + rtol * np.abs(a)) if rtol else d > atol   # NaN tol -> False
    n_both = r(both)
    out = {
        "compared": n_both,
        "exact": r(a == z),                       # NaN == NaN is False
        "diff": r(diff),
        "nan_nc_only": r(fz) - n_both,            # valid in Zarr, NaN in NC
        "nan_zarr_only": r(fa) - n_both,          # valid in NC, NaN in Zarr
        "maxd": np.where(n_both > 0, d.reshape(n, -1).max(1), np.nan),
        "sumd": d.reshape(n, -1).sum(1),
        "diff_nc_oob": np.zeros(n, np.int64),
        "diff_zarr_oob": np.zeros(n, np.int64),
    }
    if diff.any() and (lo is not None or hi is not None):
        lo_ = -np.inf if lo is None else lo
        hi_ = np.inf if hi is None else hi
        out["diff_nc_oob"] = r(diff & ((a < lo_) | (a > hi_)))
        out["diff_zarr_oob"] = r(diff & ((z < lo_) | (z > hi_)))
    return out


class SourceDiag:
    """Accumulates per-day diagnostic vectors for one source over its full record."""
    KEYS = ["nan", "mean", "min", "max", "std", "neg", "bmin", "bmax", "rough",
            "fill", "inf", "se_nan", "se_neg", "se_above", "se_inf"]

    def __init__(self, dates):
        self.dates = np.asarray(dates)
        T = len(dates)
        self.v = {k: np.full(T, np.nan) for k in self.KEYS}
        self.have_se = False

    def put(self, idx, diag, fill, inf, se=None):
        for k, a in diag.items():
            self.v[k][idx] = a
        self.v["fill"][idx] = fill
        self.v["inf"][idx] = inf
        if se is not None:
            self.have_se = True
            for k in ("se_nan", "se_neg", "se_above", "se_inf"):
                self.v[k][idx] = se[k]

    def flags(self, cells):
        """{check: bool array over this source's dates}."""
        v = self.v
        T = len(self.dates)
        nan = v["nan"]
        ok = np.isfinite(nan)
        normal = float(np.median(nan[ok])) if ok.any() else 0.0
        f = {}
        f["EMPTY"] = ok & (nan >= EMPTY_FRAC * cells)
        f["PARTIAL"] = ok & (nan > normal) & ~f["EMPTY"]
        f["EXTRA"] = ok & (nan < normal)
        f["NEG"] = v["neg"] > 0
        f["BELOW_MIN"] = v["bmin"] > 0
        f["ABOVE_MAX"] = v["bmax"] > 0
        f["INF"] = v["inf"] > 0
        f["FILL_LEAK"] = v["fill"] > 0
        f["CONSTANT"] = np.isfinite(v["min"]) & (v["min"] == v["max"])
        dm = v["mean"]
        mu, sd = np.nanmean(dm), np.nanstd(dm)
        f["OUTLIER"] = np.isfinite(dm) & (np.abs(dm - mu) / (sd + 1e-9) > Z_OUTLIER)
        rough = v["rough"]
        rfin = rough[np.isfinite(rough) & (rough > 0)]
        rmed = float(np.median(rfin)) if rfin.size else 0.0
        f["BLOWUP"] = np.isfinite(rough) & (rmed > 0) & (rough >= BLOWUP_FACTOR * rmed)
        if self.have_se:
            sn = v["se_nan"]
            sok = np.isfinite(sn)
            snorm = float(np.median(sn[sok])) if sok.any() else 0.0
            f["SE_NEG"] = v["se_neg"] > 0
            f["SE_ABOVE_CEIL"] = v["se_above"] > 0
            f["SE_INF"] = v["se_inf"] > 0
            f["SE_COV_ANOM"] = sok & (sn != snorm)
        else:
            for k in ("SE_NEG", "SE_ABOVE_CEIL", "SE_INF", "SE_COV_ANOM"):
                f[k] = np.zeros(T, bool)
        for k in f:
            f[k] = np.asarray(f[k], bool)
        stats = {"normal_nan": int(round(normal)),
                 "min": float(np.nanmin(v["min"])) if np.isfinite(v["min"]).any() else np.nan,
                 "max": float(np.nanmax(v["max"])) if np.isfinite(v["max"]).any() else np.nan,
                 "mean": float(np.nanmean(v["mean"])) if np.isfinite(v["mean"]).any() else np.nan}
        return f, stats


def se_day_stats(b, var):
    """Per-day SE counts for a cleaned SE block."""
    n = b.shape[0]
    flat = b.reshape(n, -1)
    ceil = SE_CEILING.get(var, np.inf)
    return {"se_nan": np.isnan(flat).sum(1), "se_neg": (flat < 0).sum(1),
            "se_above": (flat > ceil).sum(1)}


# ---------------------------------------------------------------------------
# Compare one pair
# ---------------------------------------------------------------------------
def _contig(idx):
    idx = np.asarray(idx)
    if idx.size and np.all(np.diff(idx) == 1):
        return slice(int(idx[0]), int(idx[-1]) + 1)
    return idx


def _grid(ds, main):
    t, la, lo = _dims(ds[main])
    return (np.asarray(ds[la].values, np.float64), np.asarray(ds[lo].values, np.float64),
            (t, la, lo))


def _load(da, dims, tsel, flip_lat):
    t, la, lo = dims
    x = da.transpose(t, la, lo).isel({t: tsel}).values
    return x[:, ::-1, :] if flip_lat else x


def compare_pair(key, nc_info, zarr_path, outdir, args):
    t_start = time.time()
    var, model, scen = key
    pair = f"{var}_{model}_{scen}"
    lo_b, hi_b = BOUNDS.get(var, (None, None))
    se_lo, se_hi = 0.0, SE_CEILING.get(var)
    row = {"pair": pair, "var": var, "model": model, "scenario": scen,
           "zarr_path": zarr_path,
           "nc_files": len(nc_info["files"]), "nc_superseded_files": nc_info["superseded"]}
    reasons = []

    # ---------------- Zarr side (header) ----------------
    zds = open_zarr(zarr_path)
    zmain, zse = pick_vars(zds, var)
    if zmain is None:
        raise RuntimeError(f"no 3-D data variable in {zarr_path}")
    zlat, zlon, zdims = _grid(zds, zmain)
    zdates = date_strings(zds[zdims[0]].values)
    row.update(calendar_zarr=calendar_of(zds, zdims[0]),
               units_zarr=zds[zmain].attrs.get("units", ""), se_var_zarr=zse or "",
               grid_zarr=f"{zlat.size}x{zlon.size}", n_time_zarr=len(zdates))
    m = ZARR_RE.match(os.path.basename(zarr_path.rstrip("/")))
    zexp = expected_dates(f"{m['y0']}-01-01", f"{m['y1']}-12-31", row["calendar_zarr"]) if m else None

    # ---------------- NC side (index pass: headers only) ----------------
    nc_entries = []           # (date, file_idx, local_idx)
    file_rows = []
    nc_cal = nc_units = nc_se = nc_main = None
    nlat = nlon = ndims = None
    for fi, f in enumerate(nc_info["files"]):
        with open_nc(f, args.nc_engine) as ds:
            mn, se = pick_vars(ds, var)
            if mn is None:
                file_rows.append([os.path.basename(f), "", "", 0, 0, "NO_DATA_VAR"])
                continue
            if nc_main is None:
                nc_main, nc_se = mn, se
                nlat, nlon, ndims = _grid(ds, mn)
                nc_cal = calendar_of(ds, ndims[0])
                nc_units = ds[mn].attrs.get("units", "")
            d = date_strings(ds[ndims[0]].values)
        fm = NC_RE.match(os.path.basename(f))
        s0 = f"{fm['start'][:4]}-{fm['start'][4:6]}-{fm['start'][6:]}"
        s1 = f"{fm['end'][:4]}-{fm['end'][4:6]}-{fm['end'][6:]}"
        exp = expected_dates(s0, s1, nc_cal or "standard")
        status = []
        if d.size == 0:
            status.append("EMPTY_FILE")
        else:
            if d[0] != s0 or d[-1] != s1:
                status.append("NAME_RANGE_MISMATCH")
            if exp is not None and len(d) != len(exp):
                status.append(f"N_DAYS {len(d)}!={len(exp)}")
            if np.any(d[1:] <= d[:-1]):
                status.append("NOT_MONOTONIC")
        file_rows.append([os.path.basename(f), d[0] if d.size else "", d[-1] if d.size else "",
                          len(d), len(exp) if exp else "", ";".join(status) or "OK"])
        for li, ds_ in enumerate(d):
            nc_entries.append((ds_, fi, li))
    if nc_main is None:
        raise RuntimeError("no readable NC data variable")

    with open(os.path.join(outdir, f"{pair}_nc_files.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "first_date", "last_date", "n_days", "n_expected", "status"])
        w.writerows(file_rows)
    row["nc_file_range_issues"] = sum(1 for r in file_rows if r[-1] != "OK")

    # dedupe NC dates (keep first occurrence in file order)
    nc_entries.sort(key=lambda e: (e[0], e[1], e[2]))
    seen, keep, dup_nc = set(), [], 0
    for e in nc_entries:
        if e[0] in seen:
            dup_nc += 1
            continue
        seen.add(e[0]); keep.append(e)
    ndates = np.array([e[0] for e in keep])
    row.update(calendar_nc=nc_cal, units_nc=nc_units, se_var_nc=nc_se or "",
               grid_nc=f"{nlat.size}x{nlon.size}", n_time_nc=len(ndates), dup_nc=dup_nc)

    # Zarr duplicates / index
    zidx_of = {}
    dup_z = 0
    for i, d in enumerate(zdates):
        if d in zidx_of:
            dup_z += 1
        else:
            zidx_of[d] = i
    row["dup_zarr"] = dup_z

    # missing days against own calendar
    nexp = None
    if keep:
        f0 = NC_RE.match(os.path.basename(nc_info["files"][0]))
        f1 = NC_RE.match(os.path.basename(nc_info["files"][-1]))
        nexp = expected_dates(f"{f0['start'][:4]}-{f0['start'][4:6]}-{f0['start'][6:]}",
                              f"{f1['end'][:4]}-{f1['end'][4:6]}-{f1['end'][6:]}", nc_cal)
    row["missing_nc"] = len(set(nexp) - set(ndates)) if nexp else ""
    row["missing_zarr"] = len(set(zexp) - set(zidx_of)) if zexp else ""

    nset, zset = set(ndates), set(zidx_of)
    common = sorted(nset & zset)
    nc_only_d = sorted(nset - zset)
    z_only_d = sorted(zset - nset)
    row.update(n_common=len(common), dates_nc_only=len(nc_only_d), dates_zarr_only=len(z_only_d))

    # ---------------- grid check ----------------
    flip = False
    grid_ok = (nlat.size == zlat.size and nlon.size == zlon.size)
    cmax = np.nan
    if grid_ok:
        if np.allclose(nlat, zlat, atol=args.coord_tol, rtol=0):
            pass
        elif np.allclose(nlat[::-1], zlat, atol=args.coord_tol, rtol=0):
            flip = True
        else:
            grid_ok = False
        if grid_ok:
            nl = nlat[::-1] if flip else nlat
            cmax = float(max(np.max(np.abs(nl - zlat)), np.max(np.abs(nlon - zlon))))
            grid_ok = np.allclose(nlon, zlon, atol=args.coord_tol, rtol=0)
    row.update(grid_match=grid_ok, coord_max_diff=cmax, lat_flipped=flip)
    cells = int(nlat.size * nlon.size)

    # ---------------- data pass ----------------
    ndiag, zdiag = SourceDiag(ndates), SourceDiag(zdates)
    nd_pos = {d: i for i, d in enumerate(ndates)}
    common_set = set(common)
    daily = {}               # date -> dict of comparison numbers (common days only)
    agg = defaultdict(float)
    agg["max_abs_diff"] = 0.0
    agg["se_max_abs_diff"] = 0.0
    do_se = bool(nc_se and zse)

    by_file = defaultdict(list)
    for d, fi, li in keep:
        by_file[fi].append((li, d))
    zda = zds[zmain]
    zse_da = zds[zse] if zse else None

    for fi in sorted(by_file):
        items = sorted(by_file[fi])
        with open_nc(nc_info["files"][fi], args.nc_engine) as ds:
            nda = ds[nc_main]
            nse_da = ds[nc_se] if nc_se and nc_se in ds else None
            for b0 in range(0, len(items), args.block_days):
                blk = items[b0:b0 + args.block_days]
                lidx = [li for li, _ in blk]
                bdates = [d for _, d in blk]
                a, a_fill, a_inf = clean_block(_load(nda, ndims, _contig(lidx), flip))
                pos = [nd_pos[d] for d in bdates]
                se_a = None
                if nse_da is not None:
                    sa, _, sa_inf = clean_block(_load(nse_da, ndims, _contig(lidx), flip))
                    se_a = se_day_stats(sa, var); se_a["se_inf"] = sa_inf
                ndiag.put(pos, day_diag(a, lo_b, hi_b), a_fill, a_inf, se_a)

                # common days -> read Zarr for the same dates
                cmask = np.array([d in common_set for d in bdates])
                if not cmask.any():
                    continue
                cd = [d for d, c in zip(bdates, cmask) if c]
                zi = [zidx_of[d] for d in cd]
                z, z_fill, z_inf = clean_block(_load(zda, zdims, _contig(zi), False))
                se_z = None
                if zse_da is not None:
                    sz, _, sz_inf = clean_block(_load(zse_da, zdims, _contig(zi), False))
                    se_z = se_day_stats(sz, var); se_z["se_inf"] = sz_inf
                zdiag.put(zi, day_diag(z, lo_b, hi_b), z_fill, z_inf, se_z)

                if not grid_ok:
                    continue
                ac = a[cmask]
                cc = cell_compare(ac, z, lo_b, hi_b, args.rtol, args.atol)
                scc = None
                if do_se and nse_da is not None and zse_da is not None:
                    scc = cell_compare(sa[cmask], sz, se_lo, se_hi, args.rtol, args.atol)
                for k, d in enumerate(cd):
                    rec = {k2: cc[k2][k] for k2 in cc}
                    if scc is not None:
                        rec.update({"se_" + k2: scc[k2][k] for k2 in scc})
                    daily[d] = rec
                for k2 in ("compared", "exact", "diff", "diff_nc_oob", "diff_zarr_oob",
                           "nan_nc_only", "nan_zarr_only", "sumd"):
                    agg[k2] += float(cc[k2].sum())
                agg["days_with_diff"] += int((cc["diff"] > 0).sum())
                agg["days_nan_mismatch"] += int(((cc["nan_nc_only"] + cc["nan_zarr_only"]) > 0).sum())
                if np.isfinite(cc["maxd"]).any():
                    agg["max_abs_diff"] = max(agg["max_abs_diff"], float(np.nanmax(cc["maxd"])))
                if scc is not None:
                    for k2 in ("compared", "diff", "nan_nc_only", "nan_zarr_only"):
                        agg["se_" + k2] += float(scc[k2].sum())
                    agg["se_days_with_diff"] += int((scc["diff"] > 0).sum())
                    if np.isfinite(scc["maxd"]).any():
                        agg["se_max_abs_diff"] = max(agg["se_max_abs_diff"],
                                                     float(np.nanmax(scc["maxd"])))

    # Zarr-only dates: diagnostics for the Zarr record still need them
    zonly_idx = sorted(zidx_of[d] for d in z_only_d)
    for b0 in range(0, len(zonly_idx), args.block_days):
        zi = zonly_idx[b0:b0 + args.block_days]
        z, z_fill, z_inf = clean_block(_load(zda, zdims, _contig(zi), False))
        se_z = None
        if zse_da is not None:
            sz, _, sz_inf = clean_block(_load(zse_da, zdims, _contig(zi), False))
            se_z = se_day_stats(sz, var); se_z["se_inf"] = sz_inf
        zdiag.put(zi, day_diag(z, lo_b, hi_b), z_fill, z_inf, se_z)
    zds.close()

    # ---------------- diagnostics per source, then compare ----------------
    nflags, nstats = ndiag.flags(cells)
    zcells = int(zlat.size * zlon.size)
    zflags, zstats = zdiag.flags(zcells)
    row.update(normal_nan_nc=nstats["normal_nan"], normal_nan_zarr=zstats["normal_nan"],
               min_nc=_r(nstats["min"]), min_zarr=_r(zstats["min"]),
               max_nc=_r(nstats["max"]), max_zarr=_r(zstats["max"]),
               mean_nc=_r(nstats["mean"]), mean_zarr=_r(zstats["mean"]))

    flag_rows = []
    for c in CHECKS:
        ns = set(ndates[nflags[c]].tolist())
        zs = set(np.asarray(zdates)[zflags[c]].tolist())
        row[f"{c}_nc"], row[f"{c}_zarr"] = len(ns), len(zs)
        row[f"{c}_ncOnly"], row[f"{c}_zarrOnly"] = len(ns - zs), len(zs - ns)
        for d in sorted(ns | zs):
            st = "BOTH" if (d in ns and d in zs) else ("NC_ONLY" if d in ns else "ZARR_ONLY")
            flag_rows.append([d, c, st])
        if ns != zs:
            reasons.append(f"{c}:ncOnly={len(ns - zs)},zarrOnly={len(zs - ns)}")
    with open(os.path.join(outdir, f"{pair}_flag_compare.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["date", "check", "status"]); w.writerows(flag_rows)

    # domain-mean series agreement on common days
    if common:
        ai = np.array([nd_pos[d] for d in common]); zi = np.array([zidx_of[d] for d in common])
        dn, dz = ndiag.v["mean"][ai], zdiag.v["mean"][zi]
        ok = np.isfinite(dn) & np.isfinite(dz)
        row["dmean_max_abs_diff"] = _r(float(np.max(np.abs(dn[ok] - dz[ok])))) if ok.any() else ""
        row["dmean_corr"] = (_r(float(np.corrcoef(dn[ok], dz[ok])[0, 1]))
                             if ok.sum() > 2 and np.std(dn[ok]) > 0 and np.std(dz[ok]) > 0 else "")

    # ---------------- per-date difference table ----------------
    with open(os.path.join(outdir, f"{pair}_diff_dates.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["date", "in_nc", "in_zarr", "nan_nc_only", "nan_zarr_only", "diff_cells",
                    "diff_nc_oob", "diff_zarr_oob", "max_abs_diff", "mean_nc", "mean_zarr",
                    "se_nan_nc_only", "se_nan_zarr_only", "se_diff_cells", "se_max_abs_diff"])
        for d in nc_only_d:
            w.writerow([d, 1, 0] + [""] * 12)
        for d in z_only_d:
            w.writerow([d, 0, 1] + [""] * 12)
        for d in common:
            r_ = daily.get(d)
            if r_ is None:
                continue
            bad = (r_["diff"] or r_["nan_nc_only"] or r_["nan_zarr_only"] or
                   r_.get("se_diff", 0) or r_.get("se_nan_nc_only", 0) or
                   r_.get("se_nan_zarr_only", 0))
            if not bad:
                continue
            w.writerow([d, 1, 1, int(r_["nan_nc_only"]), int(r_["nan_zarr_only"]),
                        int(r_["diff"]), int(r_["diff_nc_oob"]), int(r_["diff_zarr_oob"]),
                        _r(r_["maxd"]), _r(ndiag.v["mean"][nd_pos[d]]),
                        _r(zdiag.v["mean"][zidx_of[d]]),
                        int(r_.get("se_nan_nc_only", 0)), int(r_.get("se_nan_zarr_only", 0)),
                        int(r_.get("se_diff", 0)), _r(r_.get("se_maxd", np.nan))])

    # ---------------- summary row + verdict ----------------
    row.update(
        cells_compared=int(agg["compared"]), cells_exact=int(agg["exact"]),
        cells_diff=int(agg["diff"]), cells_diff_nc_oob=int(agg["diff_nc_oob"]),
        cells_diff_zarr_oob=int(agg["diff_zarr_oob"]), days_with_diff=int(agg["days_with_diff"]),
        max_abs_diff=_r(agg["max_abs_diff"]) if agg["compared"] else "",
        mean_abs_diff=_r(agg["sumd"] / agg["compared"]) if agg["compared"] else "",
        nan_nc_only_cells=int(agg["nan_nc_only"]), nan_zarr_only_cells=int(agg["nan_zarr_only"]),
        days_nan_mismatch=int(agg["days_nan_mismatch"]),
        se_cells_compared=int(agg["se_compared"]), se_cells_diff=int(agg["se_diff"]),
        se_days_with_diff=int(agg["se_days_with_diff"]),
        se_max_abs_diff=_r(agg["se_max_abs_diff"]) if agg["se_compared"] else "",
        se_nan_nc_only_cells=int(agg["se_nan_nc_only"]),
        se_nan_zarr_only_cells=int(agg["se_nan_zarr_only"]))

    if not grid_ok:
        reasons.insert(0, f"GRID {row['grid_nc']} vs {row['grid_zarr']} (cell compare skipped)")
    if str(row["units_nc"]) != str(row["units_zarr"]):
        reasons.append(f"units '{row['units_nc']}' vs '{row['units_zarr']}'")
    if bool(nc_se) != bool(zse):
        reasons.append(f"SE var nc='{nc_se or '-'}' zarr='{zse or '-'}'")
    if _cal_norm(nc_cal) != _cal_norm(row["calendar_zarr"]):
        reasons.append(f"calendar {nc_cal} vs {row['calendar_zarr']}")
    for k in ("dates_nc_only", "dates_zarr_only", "dup_nc", "dup_zarr", "missing_nc",
              "missing_zarr", "nc_file_range_issues", "nan_nc_only_cells",
              "nan_zarr_only_cells", "cells_diff", "se_cells_diff",
              "se_nan_nc_only_cells", "se_nan_zarr_only_cells"):
        v = row.get(k)
        if v not in ("", None, 0):
            reasons.append(f"{k}={v}")
    if not grid_ok:
        verdict = "GRID_MISMATCH"
    elif not reasons:
        verdict = "IDENTICAL"
    elif (row["cells_diff"] == 0 and row["nan_nc_only_cells"] == 0 and row["nan_zarr_only_cells"] == 0
          and row["se_cells_diff"] == 0 and row["dates_nc_only"] == 0 and row["dates_zarr_only"] == 0):
        verdict = "SAME_DATA_META_DIFF"      # data equal; metadata/flags differ
    else:
        verdict = "DIFFERENT"
    row["verdict"] = verdict
    row["reasons"] = "; ".join(reasons)
    row["seconds"] = round(time.time() - t_start, 1)
    return row


def _cal_norm(c):
    c = str(c).lower()
    return {"365_day": "noleap", "gregorian": "standard"}.get(c, c)


def _r(x, nd=6):
    try:
        x = float(x)
    except Exception:
        return ""
    return "" if not np.isfinite(x) else float(f"{x:.{nd}g}")


def _worker(key, nc_info, zpath, outdir, args):
    try:
        return compare_pair(key, nc_info, zpath, outdir, args)
    except Exception as e:
        var, model, scen = key
        return {"pair": f"{var}_{model}_{scen}", "var": var, "model": model, "scenario": scen,
                "zarr_path": zpath,
                "verdict": f"ERROR:{type(e).__name__}",
                "reasons": f"{e} | {traceback.format_exc(limit=3)!r}"[:2000]}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nc-root", default=NC_ROOT_DEFAULT)
    ap.add_argument("--zarr-root", default=ZARR_ROOT_DEFAULT)
    ap.add_argument("--out", default="cmp_nc_zarr_out")
    ap.add_argument("--vars", default=",".join(VARIABLES))
    ap.add_argument("--models", default="", help="comma list (default: all)")
    ap.add_argument("--scenarios", default="", help="comma list, e.g. historical,ssp585")
    ap.add_argument("--workers", type=int, default=4, help="pairs processed in parallel")
    ap.add_argument("--block-days", type=int, default=92,
                    help="days read per block (RAM ~ 6 x block x cells x 8 bytes per worker)")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=0.0)
    ap.add_argument("--coord-tol", type=float, default=1e-5, help="lat/lon tolerance (deg)")
    ap.add_argument("--nc-engine", default="netcdf4", choices=["netcdf4", "h5netcdf"])
    ap.add_argument("--limit", type=int, default=0, help="debug: cap number of pairs")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    only_vars = {v for v in args.vars.split(",") if v}
    only_models = {v for v in args.models.split(",") if v}
    only_scen = {canon_scen(v) for v in args.scenarios.split(",") if v}

    def wanted(k):
        return ((not only_vars or k[0] in only_vars) and
                (not only_models or k[1] in only_models) and
                (not only_scen or k[2] in only_scen))

    t0 = time.time()
    print(f"[discover] NC   under {args.nc_root}", flush=True)
    nc = {k: v for k, v in discover_nc(args.nc_root).items() if wanted(k)}
    print(f"[discover] Zarr under {args.zarr_root}", flush=True)
    zr = discover_zarr(args.zarr_root)
    zr = {k: v for k, v in zr.items() if wanted(k)}

    keys = sorted(set(nc) | set(zr))
    with open(os.path.join(args.out, "inventory.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["var", "model", "scenario", "in_nc", "in_zarr", "nc_files",
                    "nc_superseded", "zarr_path"])
        for k in keys:
            n = nc.get(k); z = zr.get(k)
            w.writerow([*k, int(n is not None), int(z is not None),
                        len(n["files"]) if n else 0, n["superseded"] if n else 0,
                        z or ""])
    both = [k for k in keys if k in nc and k in zr]
    nc_only = [k for k in keys if k in nc and k not in zr]
    z_only = [k for k in keys if k in zr and k not in nc]
    print(f"[inventory] pairs: both={len(both)}  NC-only={len(nc_only)}  Zarr-only={len(z_only)}",
          flush=True)

    jobs = [(k, nc[k], zr[k]) for k in both]
    if args.limit:
        jobs = jobs[:args.limit]

    summary_csv = os.path.join(args.out, "compare_summary.csv")
    done = set()
    if os.path.exists(summary_csv):
        with open(summary_csv) as f:
            done = {r["pair"] for r in csv.DictReader(f) if not r["verdict"].startswith("ERROR")}

    def pid(j):
        k = j[0]
        return f"{k[0]}_{k[1]}_{k[2]}"
    todo = [j for j in jobs if pid(j) not in done]
    print(f"[scan] {len(jobs)} pairs | {len(done)} already done | {len(todo)} to do", flush=True)

    new = not os.path.exists(summary_csv)
    fh = open(summary_csv, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
    if new:
        w.writeheader()
    # one-sided pairs get a row once (no heavy work)
    for k in nc_only + z_only:
        p = f"{k[0]}_{k[1]}_{k[2]}"
        if p in done:
            continue
        w.writerow({"pair": p, "var": k[0], "model": k[1], "scenario": k[2],
                    "verdict": "NC_ONLY" if k in nc else "ZARR_ONLY",
                    "zarr_path": zr.get(k, ""),
                    "nc_files": len(nc[k]["files"]) if k in nc else 0})
    fh.flush()

    n = 0
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = {ex.submit(_worker, k, ni, zp, args.out, args): k
                for (k, ni, zp) in todo}
        for fut in as_completed(futs):
            r = fut.result(); n += 1
            w.writerow({f: r.get(f, "") for f in SUMMARY_FIELDS}); fh.flush()
            v = r.get("verdict", "?")
            tag = {"IDENTICAL": "SAME ", "SAME_DATA_META_DIFF": "META "}.get(
                v, "ERR  " if v.startswith("ERROR") else "DIFF ")
            extra = "" if v == "IDENTICAL" else f"  {r.get('reasons', '')[:300]}"
            print(f"[{n}/{len(todo)}] {tag} {r['pair']}  ({r.get('seconds', '')}s){extra}",
                  flush=True)
    fh.close()
    write_report(summary_csv, args.out, time.time() - t0)


def write_report(summary_csv, outdir, secs):
    with open(summary_csv) as f:
        rows = list(csv.DictReader(f))
    # keep the latest row per pair (a rerun after ERROR appends a new row)
    last = {}
    for r in rows:
        last[r["pair"]] = r
    rows = list(last.values())

    def i(r, k):
        try:
            return int(float(r.get(k) or 0))
        except ValueError:
            return 0

    L = ["=" * 78, "  SETx-1km-CMIP6  NC vs Zarr  --  COMPARISON SUMMARY", "=" * 78]
    verd = defaultdict(int)
    for r in rows:
        verd[r["verdict"].split(":")[0]] += 1
    L.append("  verdicts: " + "  ".join(f"{k}={v}" for k, v in sorted(verd.items())))
    L.append("-" * 78)
    L.append(f"  {'var':8s} {'pairs':>5s} {'same':>5s} {'diff':>5s} {'1side':>5s} {'err':>4s}"
             f" {'dNC':>6s} {'dZarr':>6s} {'nanMis':>11s} {'diffCells':>12s} {'seDiff':>11s}")
    by_var = defaultdict(list)
    for r in rows:
        by_var[r["var"]].append(r)
    for v in sorted(by_var):
        rs = by_var[v]
        L.append(f"  {v:8s} {len(rs):5d} "
                 f"{sum(r['verdict'] == 'IDENTICAL' for r in rs):5d} "
                 f"{sum(r['verdict'] in ('DIFFERENT', 'GRID_MISMATCH', 'SAME_DATA_META_DIFF') for r in rs):5d} "
                 f"{sum(r['verdict'] in ('NC_ONLY', 'ZARR_ONLY') for r in rs):5d} "
                 f"{sum(r['verdict'].startswith('ERROR') for r in rs):4d} "
                 f"{sum(i(r, 'dates_nc_only') for r in rs):6d} "
                 f"{sum(i(r, 'dates_zarr_only') for r in rs):6d} "
                 f"{sum(i(r, 'nan_nc_only_cells') + i(r, 'nan_zarr_only_cells') for r in rs):11d} "
                 f"{sum(i(r, 'cells_diff') for r in rs):12d} "
                 f"{sum(i(r, 'se_cells_diff') for r in rs):11d}")
    L.append("-" * 78)
    L.append("  diagnostic flags that DISAGREE between sources (dates, all pairs):")
    for c in CHECKS:
        a = sum(i(r, f"{c}_ncOnly") for r in rows)
        b = sum(i(r, f"{c}_zarrOnly") for r in rows)
        if a or b:
            L.append(f"    {c:14s} NC-only={a:8d}   Zarr-only={b:8d}")
    L.append("-" * 78)
    bad = [r for r in rows if r["verdict"] not in ("IDENTICAL",)]
    for r in sorted(bad, key=lambda r: r["pair"]):
        L.append(f"  [{r['verdict']}] {r['pair']}")
        if r.get("reasons"):
            for item in r["reasons"].split("; ")[:12]:
                L.append(f"        - {item}")
    L.append("=" * 78)
    L.append("  columns: dNC/dZarr = dates only in NC / only in Zarr; nanMis = cells NaN in "
             "one source only; diffCells = cells beyond tolerance.")
    L.append("  drill down: <pair>_diff_dates.csv, <pair>_flag_compare.csv, <pair>_nc_files.csv")
    L.append(f"  runtime {secs / 60:.1f} min")
    report = "\n".join(L)
    print("\n" + report, flush=True)
    with open(os.path.join(outdir, "SUMMARY.txt"), "w") as f:
        f.write(report + "\n")


if __name__ == "__main__":
    main()
