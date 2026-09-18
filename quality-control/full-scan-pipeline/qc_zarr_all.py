#!/usr/bin/env python3
"""
qc_zarr_all.py  --  Zarr-only QC scanner for the SETx-1km-CMIP6 stores.

Author: Ifeanyichukwu C. Nduka (icnduka@jsg.utexas.edu)

Reads each Zarr store's CONTENT (no NetCDF comparison) and scans for:
  * MISSING days      -- calendar-complete date audit (store's own calendar)
  * EMPTY days        -- whole-field NaN
  * PARTIAL days      -- more NaN than the store's typical (median) footprint
  * EXTRA-valid days  -- fewer NaN than typical (spurious extended coverage,
                         e.g. the ssp585 2099-08-15 case)
  * NEGATIVE values   -- value < 0 (unphysical for wind/humidity/radiation)
  * OUT-OF-RANGE      -- below/above per-variable physical bounds
  * inf, constant slices, temporal outliers
  * SE variable (<var>_SE): negatives, inf, and its own coverage anomalies

Design lessons from prior QC iterations are baked in:
  - Coverage baseline is the MEDIAN per-day NaN count, NOT a union-of-finite
    land mask (union let one odd day flag every day partial).
  - One vectorised, memory-bounded pass over time (block reads).
  - Per-variable physical bounds.
  - Resumable: stores already in the summary CSV are skipped on restart.

Paths: ALL variables (including sfcWind) are scanned from --root; sfcWind is
ADDITIONALLY scanned from --sfcwind-root (the post-processed store). By default a
POST store SUPERSEDES the same (model,scenario) primary store, so the original
pre-re-downscale primary CESM2 sfcWind is dropped and won't confuse readers; the
nine other models' primary sfcWind are still scanned. Use --keep-superseded to
scan both (before/after; the POST copy is then tagged "_POST"). --vars restricts
to a variable list and --skip excludes stores by var[:model[:scenario[:source]]].

Run under nohup:
  nohup python3 -u qc_zarr_all.py > qc_zarr.log 2>&1 &
  tail -f qc_zarr.log
"""

import argparse, csv, os, sys, time, warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import xarray as xr

warnings.filterwarnings("ignore", category=RuntimeWarning)   # all-NaN-slice reductions

# ---- per-variable physical bounds (min, max) ----
BOUNDS = {
    "tas": (180., 350.), "tasmax": (180., 350.), "tasmin": (180., 350.),
    "hurs": (0., 105.), "huss": (0., 0.06), "sfcWind": (0., 50.),
    "rlds": (40., 650.), "rsds": (0., 1500.),
}
EMPTY_FRAC = 0.999          # >= this fraction NaN -> empty day
BLOCK = 365                 # timesteps per read
BLOWUP_FACTOR = 6.0         # day roughness >= this x the store median -> spatial blow-up
SUMMARY_FIELDS = [
    "store", "source", "var", "scenario", "calendar", "n_time", "n_expected", "n_missing",
    "n_duplicate", "time_monotonic", "normal_nan", "land_cells",
    "days_empty", "days_partial", "days_extra",
    "neg_cells", "below_min_cells", "above_max_cells", "inf_cells",
    "days_constant", "days_blowup", "rough_median", "days_outlier",
    "global_min", "global_max", "mean_val",
    "se_var", "se_neg", "se_inf", "se_anom_days", "verdict",
]


def find_zarr(root):
    out = []
    if not os.path.isdir(root):
        return out
    for dp, dns, _ in os.walk(root):
        for d in [x for x in dns if x.endswith(".zarr")]:
            out.append(os.path.join(dp, d))
        dns[:] = [x for x in dns if not x.endswith(".zarr")]
    return sorted(out)


def parse_name(path):
    p = os.path.basename(path)[:-5].split("_")
    return {"model": p[0] if p else "?",
            "scenario": p[1] if len(p) > 1 else "?",
            "var": p[2] if len(p) > 2 else "?",
            "range": p[3] if len(p) > 3 else ""}


def open_store(path):
    try:
        return xr.open_zarr(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True),
                            consolidated=False)
    except (AttributeError, TypeError):
        return xr.open_zarr(path, use_cftime=True, consolidated=False)


def pick_vars(ds, hint):
    d3 = [v for v in ds.data_vars if {"time", "lat", "lon"} <= set(ds[v].dims)]
    se = next((v for v in d3 if v.endswith("_SE") or "SE" in v), None)
    main = hint if hint in d3 else None
    if main is None:
        main = next((v for v in d3 if v != se), None)
    return main, se


def per_day_nan(da, T):
    """Vectorised per-day NaN count in one pass."""
    nan_ct = np.empty(T, np.int64)
    for t0 in range(0, T, BLOCK):
        b = np.asarray(da.isel(time=slice(t0, t0 + BLOCK)).values)
        nan_ct[t0:t0 + b.shape[0]] = np.isnan(b).reshape(b.shape[0], -1).sum(1)
    return nan_ct


def qc_store(path, outdir, store_id=None, source=""):
    name = store_id or os.path.basename(path)[:-5]
    meta = parse_name(path)
    ds = open_store(path)
    var, se = pick_vars(ds, meta["var"])
    da = ds[var]
    T, LA, LO = da.shape
    cells = LA * LO
    cal = getattr(getattr(ds.time, "dt", None), "calendar", None) or \
          ds.time.encoding.get("calendar", "noleap")
    dates = np.array([f"{t:%Y-%m-%d}" for t in ds.time.values])

    # ---- time completeness ----
    rng = meta["range"]
    if rng and "-" in rng and rng.split("-")[0].isdigit():
        y0, y1 = (int(x) for x in rng.split("-"))
    else:
        y0, y1 = int(dates[0][:4]), int(dates[-1][:4])
    try:
        exp = [f"{t:%Y-%m-%d}" for t in
               xr.date_range(f"{y0}-01-01", f"{y1}-12-31", freq="D", calendar=cal, use_cftime=True)]
    except Exception:
        exp = list(dates)
    have = set(dates)
    missing = [d for d in exp if d not in have]
    dup = sorted({d for d in dates if (dates == d).sum() > 1})
    times = ds.time.values
    monotonic = bool(np.all(times[1:] > times[:-1]))

    # ---- one vectorised pass ----
    vmin, vmax = BOUNDS.get(var, (None, None))
    nan_ct = np.empty(T, np.int64); dmean = np.empty(T); dmin = np.empty(T)
    dmax = np.empty(T); dstd = np.empty(T)
    n_neg = np.empty(T, np.int64); n_bmin = np.empty(T, np.int64)
    n_bmax = np.empty(T, np.int64); n_inf = np.empty(T, np.int64)
    rough = np.empty(T)                              # mean |adjacent-cell diff| (spatial roughness)
    for t0 in range(0, T, BLOCK):
        b = np.asarray(da.isel(time=slice(t0, t0 + BLOCK)).values, dtype=np.float64)
        s = slice(t0, t0 + b.shape[0]); flat = b.reshape(b.shape[0], -1)
        fin = np.isfinite(flat)
        nan_ct[s] = np.isnan(flat).sum(1)
        with np.errstate(all="ignore"):
            dmean[s] = np.nanmean(flat, 1); dmin[s] = np.nanmin(flat, 1)
            dmax[s] = np.nanmax(flat, 1); dstd[s] = np.nanstd(flat, 1)
            ry = np.nanmean(np.abs(b[:, 1:, :] - b[:, :-1, :]), axis=(1, 2))
            rx = np.nanmean(np.abs(b[:, :, 1:] - b[:, :, :-1]), axis=(1, 2))
            rough[s] = 0.5 * (ry + rx)
        n_neg[s] = (flat < 0).sum(1)
        n_bmin[s] = (flat < vmin).sum(1) if vmin is not None else 0
        n_bmax[s] = (flat > vmax).sum(1) if vmax is not None else 0
        n_inf[s] = np.isinf(flat).sum(1)

    # ---- coverage classification (median baseline) ----
    normal = int(round(np.median(nan_ct)))
    empty = nan_ct >= EMPTY_FRAC * cells
    partial = (nan_ct > normal) & ~empty
    extra = nan_ct < normal
    land_cells = cells - normal

    const = (dstd == 0) & np.isfinite(dstd)
    with np.errstate(invalid="ignore"):
        z = (dmean - np.nanmean(dmean)) / (np.nanstd(dmean) + 1e-9)
    outlier = np.abs(z) > 5
    # spatial-incoherence ("kriging blow-up"): roughness far above the store's own median.
    rfin = rough[np.isfinite(rough) & (rough > 0)]
    rmed = float(np.median(rfin)) if rfin.size else 0.0
    rough_ratio = np.where((rmed > 0) & np.isfinite(rough), rough / max(rmed, 1e-12), 0.0)
    blowup = np.isfinite(rough) & (rmed > 0) & (rough_ratio >= BLOWUP_FACTOR)
    fv = np.isfinite(dmin)
    gmin = float(np.nanmin(dmin[fv])) if fv.any() else float("nan")
    gmax = float(np.nanmax(dmax[fv])) if fv.any() else float("nan")

    # ---- SE ----
    se_neg = se_inf = 0; se_anom = []
    if se:
        se_nan = np.empty(T, np.int64)
        for t0 in range(0, T, BLOCK):
            sb = np.asarray(ds[se].isel(time=slice(t0, t0 + BLOCK)).values, dtype=np.float64)
            f = sb.reshape(sb.shape[0], -1)
            se_nan[t0:t0 + sb.shape[0]] = np.isnan(f).sum(1)
            se_neg += int((f < 0).sum()); se_inf += int(np.isinf(f).sum())
        se_norm = int(round(np.median(se_nan)))
        se_anom = [dates[i] for i in np.where(se_nan != se_norm)[0]]

    # ---- write per-store CSVs ----
    base = os.path.join(outdir, name)
    with open(base + "_missing_dates.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["date"]); w.writerows([[d] for d in missing])
    with open(base + "_incomplete_dates.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["date", "nan_count", "deviation", "class"])
        for i in range(T):
            if empty[i] or partial[i] or extra[i]:
                cls = "EMPTY" if empty[i] else ("PARTIAL" if partial[i] else "EXTRA")
                w.writerow([dates[i], int(nan_ct[i]), int(nan_ct[i] - normal), cls])
    with open(base + "_oob_dates.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["date", "reasons"])
        for i in range(T):
            r = []
            if n_neg[i]:  r.append(f"neg:{int(n_neg[i])}")
            if n_bmin[i]: r.append(f"below_min:{int(n_bmin[i])}")
            if n_bmax[i]: r.append(f"above_max:{int(n_bmax[i])}")
            if n_inf[i]:  r.append(f"inf:{int(n_inf[i])}")
            if const[i]:  r.append("constant")
            if blowup[i]: r.append(f"incoherent(rough={rough_ratio[i]:.0f}x)")
            if outlier[i]: r.append("outlier")
            if r: w.writerow([dates[i], ";".join(r)])

    row = {
        "store": name, "source": source, "var": var, "scenario": meta["scenario"], "calendar": cal,
        "n_time": T, "n_expected": len(exp), "n_missing": len(missing),
        "n_duplicate": len(dup), "time_monotonic": monotonic,
        "normal_nan": normal, "land_cells": land_cells,
        "days_empty": int(empty.sum()), "days_partial": int(partial.sum()),
        "days_extra": int(extra.sum()),
        "neg_cells": int(n_neg.sum()), "below_min_cells": int(n_bmin.sum()),
        "above_max_cells": int(n_bmax.sum()), "inf_cells": int(n_inf.sum()),
        "days_constant": int(const.sum()), "days_blowup": int(blowup.sum()),
        "rough_median": round(rmed, 4), "days_outlier": int(outlier.sum()),
        "global_min": round(gmin, 4), "global_max": round(gmax, 4),
        "mean_val": round(float(np.nanmean(dmean)), 4),
        "se_var": se or "", "se_neg": se_neg, "se_inf": se_inf,
        "se_anom_days": len(se_anom),
    }
    clean = (row["n_missing"] == 0 and row["n_duplicate"] == 0 and row["time_monotonic"] and
             row["days_empty"] == 0 and row["days_partial"] == 0 and row["days_extra"] == 0 and
             row["neg_cells"] == 0 and row["below_min_cells"] == 0 and row["above_max_cells"] == 0 and
             row["inf_cells"] == 0 and row["days_constant"] == 0 and
             row["days_blowup"] == 0 and se_neg == 0 and
             se_inf == 0 and len(se_anom) == 0)
    row["verdict"] = "CLEAN" if clean else "ISSUES"
    return row


def _worker(path, outdir, store_id, source):
    try:
        return qc_store(path, outdir, store_id, source)
    except Exception as e:
        return {"store": store_id, "source": source, "verdict": f"ERROR:{type(e).__name__}:{e}"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/local1/SETx_UIFL_Zarr")
    ap.add_argument("--sfcwind-root", default="/local1/SFCWIND_ZARR_POST")
    ap.add_argument("--out", default="qc_zarr_out")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0, help="debug: cap number of stores")
    ap.add_argument("--vars", default="",
                    help="only scan these variables (comma list), e.g. rlds,rsds")
    ap.add_argument("--skip", default="",
                    help="skip stores matching any spec (comma list). A spec is "
                         "var[:model[:scenario[:source]]], e.g. "
                         "sfcWind:CESM2:*:primary  or  pr  or  huss:MPI-ESM1-2-HR")
    ap.add_argument("--keep-superseded", action="store_true",
                    help="also scan primary stores that have a POST replacement "
                         "(default: POST supersedes primary, so the pre-re-downscale "
                         "primary CESM2 sfcWind is dropped)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    summary_csv = os.path.join(args.out, "qc_zarr_summary.csv")

    # ---- build filters ----
    only_vars = {v for v in args.vars.split(",") if v}
    skip_specs = [s.split(":") for s in args.skip.split(",") if s]

    def skipped(m, source):
        if only_vars and m["var"] not in only_vars:
            return True
        for spec in skip_specs:                       # spec = [var, model?, scen?, source?]
            fields = [m["var"], m["model"], m["scenario"], source]
            if all(sp in ("", "*") or sp == f for sp, f in zip(spec, fields)):
                return True
        return False

    # discover: ALL variables from --root (incl. sfcWind), PLUS sfcWind from
    # --sfcwind-root. By default a POST store SUPERSEDES the same-named primary
    # store (so the pre-re-downscale primary CESM2 sfcWind is dropped and doesn't
    # confuse the results); pass --keep-superseded to scan both. When both are
    # kept, the POST copy is tagged "_POST" to avoid a name collision.
    post_paths = find_zarr(args.sfcwind_root)
    post_keys = {(parse_name(p)["var"], parse_name(p)["model"], parse_name(p)["scenario"])
                 for p in post_paths}

    jobs = []                                          # (path, store_id, source)
    for p in find_zarr(args.root):
        m = parse_name(p); key = (m["var"], m["model"], m["scenario"])
        if skipped(m, "primary"):
            continue
        if key in post_keys and not args.keep_superseded:
            continue                                   # superseded by POST
        jobs.append((p, os.path.basename(p)[:-5], "primary"))
    primary_ids = {j[1] for j in jobs}
    for p in post_paths:
        m = parse_name(p)
        if skipped(m, "post"):
            continue
        sid = os.path.basename(p)[:-5]
        if sid in primary_ids:                         # collision (only in --keep-superseded)
            sid += "_POST"
        jobs.append((p, sid, "post"))
    seen = set()
    jobs = [j for j in jobs if not (j[1] in seen or seen.add(j[1]))]  # dedup by store_id
    if args.limit:
        jobs = jobs[:args.limit]

    done = set()
    if os.path.exists(summary_csv):
        with open(summary_csv) as f:
            done = {r["store"] for r in csv.DictReader(f)}
    todo = [j for j in jobs if j[1] not in done]
    print(f"[scan] {len(jobs)} stores found | {len(done)} already done | {len(todo)} to do",
          flush=True)

    new = not os.path.exists(summary_csv)
    fh = open(summary_csv, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
    if new:
        w.writeheader(); fh.flush()

    t0 = time.time(); n = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_worker, p, args.out, sid, src): sid for (p, sid, src) in todo}
        for fut in as_completed(futs):
            r = fut.result(); n += 1
            w.writerow({k: r.get(k, "") for k in SUMMARY_FIELDS}); fh.flush()
            v = r.get("verdict", "?")
            tag = "OK   " if v == "CLEAN" else ("ERR  " if v.startswith("ERROR") else "ISSUE")
            extra = "" if v in ("CLEAN",) else (
                f"  miss={r.get('n_missing')} empty={r.get('days_empty')} "
                f"part={r.get('days_partial')} extra={r.get('days_extra')} "
                f"neg={r.get('neg_cells')} oob={r.get('below_min_cells','')}+{r.get('above_max_cells','')}"
                f" blowup={r.get('days_blowup')}")
            print(f"[{n}/{len(todo)}] {tag} {r['store']}{extra}", flush=True)
    fh.close()

    # ---- overall summary ----
    rows = []
    with open(summary_csv) as f:
        rows = list(csv.DictReader(f))
    L = ["=" * 74, "  SETx-1km-CMIP6 Zarr QC — OVERALL SUMMARY", "=" * 74]
    issue = [r for r in rows if r.get("verdict") != "CLEAN"]
    for r in issue:
        L.append(f"  [{r['verdict']}] {r['store']}")
        if r.get("verdict", "").startswith("ERROR"):
            continue
        d = []
        for k, lab in [("n_missing", "missing days"), ("days_empty", "empty days"),
                       ("days_partial", "partial days"), ("days_extra", "extra-valid days"),
                       ("neg_cells", "negative cells"), ("below_min_cells", "below-min cells"),
                       ("above_max_cells", "above-max cells"), ("inf_cells", "inf cells"),
                       ("days_constant", "constant days"),
                       ("days_blowup", "spatial-blowup days"), ("se_neg", "SE negatives"),
                       ("se_anom_days", "SE anomaly days")]:
            v = r.get(k, "")
            if v not in ("", "0", 0):
                d.append(f"{v} {lab}")
        if r.get("n_duplicate") not in ("", "0"):
            d.append(f"{r['n_duplicate']} duplicate dates")
        for item in d:
            L.append(f"            - {item}")

    def tot(k):
        return sum(int(r.get(k) or 0) for r in rows if not r.get("verdict", "").startswith("ERROR"))
    # per-variable roll-up
    by_var = {}
    for r in rows:
        by_var.setdefault(r.get("var", "?"), []).append(r)
    L.append("-" * 74)
    L.append("  per-variable  (stores  |  issue-stores  |  missing  empty  partial  neg-cells)")
    for v in sorted(by_var):
        rs = by_var[v]
        iss = sum(1 for x in rs if x.get("verdict") != "CLEAN")
        L.append(f"    {v:8s}  {len(rs):3d}      {iss:3d}          "
                 f"{sum(int(x.get('n_missing') or 0) for x in rs):6d} "
                 f"{sum(int(x.get('days_empty') or 0) for x in rs):6d} "
                 f"{sum(int(x.get('days_partial') or 0) for x in rs):7d} "
                 f"{sum(int(x.get('neg_cells') or 0) for x in rs):9d}")
    L.append("-" * 74)
    L.append(f"  stores checked : {len(rows)}")
    L.append(f"  totals         : missing={tot('n_missing')}  empty={tot('days_empty')}  "
             f"partial={tot('days_partial')}  extra={tot('days_extra')}  "
             f"neg-cells={tot('neg_cells')}  oob={tot('below_min_cells')+tot('above_max_cells')}  "
             f"inf={tot('inf_cells')}")
    L.append("-" * 74)
    n_err = sum(1 for r in rows if r.get("verdict", "").startswith("ERROR"))
    if not issue:
        L.append(f"  RESULT: \u2713 NO ISSUES — all {len(rows)} stores passed every check.")
    else:
        L.append(f"  RESULT: \u2717 ISSUES in {len(issue)}/{len(rows)} stores"
                 + (f"  ({n_err} errored)" if n_err else "") + ".")
        L.append("          drill down via the per-store *_incomplete_dates.csv / *_oob_dates.csv")
    L.append("=" * 74)
    report = "\n".join(L)
    print("\n" + report, flush=True)
    with open(os.path.join(args.out, "SUMMARY.txt"), "w") as f:
        f.write(report + "\n")
    print(f"\n[done] {len(rows)} stores in {(time.time()-t0)/60:.1f} min -> {args.out}/", flush=True)


if __name__ == "__main__":
    main()
