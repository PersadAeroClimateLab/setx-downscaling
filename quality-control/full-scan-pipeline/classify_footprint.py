#!/usr/bin/env python3
"""
classify_footprint.py  --  Find days/stores that are MISSING a chunk of the inland
domain, measured against a complete reference variable (tas).

Why a reference: deriving the land mask from a store's own valid-fraction fails
when a region is absent on most days -- it then looks like ocean and is never
flagged (the bug that let the NE-corner-missing days slip through). tas is
spatially complete, so its finite footprint defines the true land domain; any
store missing part of that INLAND domain on a given day is a real coverage gap.

Per store (one pass) it measures, each day, how many cells of the reference
INLAND domain are NaN, and flags days missing >= --inland-miss-min. By default
ALL such dates are emitted (--max-inland-dates 0) for a comprehensive review;
set --max-inland-dates N to emit only a representative sample per store instead.
The full per-store count is always reported in footprint_classification.csv.

Coastal breathing right at the ocean mask is still reported but treated as benign
(edge_wobble). Optionally (--bad-dates) it also locates whether out-of-range
cells sit on the coast or form an inland block.

Run:
  python3 classify_footprint.py --root /local1/SETx_UIFL_Zarr \
      --sfcwind-root /local1/SFCWIND_ZARR_POST --bad-dates bad_dates.csv
"""

import argparse, csv, glob, os, warnings
import numpy as np
import xarray as xr
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

warnings.filterwarnings("ignore", category=RuntimeWarning)

BLOCK = 365
OCEAN_VF = 0.02          # valid on <2% of days -> permanent mask/ocean (for coastal calc)
BREATHE_HI = 0.99        # valid on <99% of days -> "breathing"
PERSIST_LO = 0.50        # ref inland cell present <50% of days -> systematic deficit (map)
EDGE_DIST = 3            # cells within this many px of the mask -> coastal
BLOCK_THRESH = 500       # a contiguous inland block this big -> real cutoff (oob location)
REF_VAR = "tas"          # complete reference variable that defines the land domain
BOUNDS = {"tas": (180., 350.), "tasmax": (180., 350.), "tasmin": (180., 350.),
          "hurs": (0., 105.), "huss": (0., 0.06), "sfcWind": (0., 50.),
          "rlds": (40., 650.), "rsds": (0., 1500.)}


def find_zarr(root):
    out = []
    if os.path.isdir(root):
        for dp, dns, _ in os.walk(root):
            for d in [x for x in dns if x.endswith(".zarr")]:
                out.append(os.path.join(dp, d))
            dns[:] = [x for x in dns if not x.endswith(".zarr")]
    return sorted(out)


def parse(path):
    p = os.path.basename(path)[:-5].split("_")
    return dict(model=p[0], scenario=p[1] if len(p) > 1 else "?",
                var=p[2] if len(p) > 2 else "?")


def open_store(path):
    try:
        return xr.open_zarr(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True),
                            consolidated=False)
    except (AttributeError, TypeError):
        return xr.open_zarr(path, use_cftime=True, consolidated=False)


def main_var(ds, hint):
    d3 = [v for v in ds.data_vars if {"time", "lat", "lon"} <= set(ds[v].dims)]
    se = next((v for v in d3 if v.endswith("_SE") or "SE" in v), None)
    return (hint if hint in d3 else next((v for v in d3 if v != se), None)), se


def classify_block(mask, dist2coast):
    lab, n = ndimage.label(mask)
    if n == 0:
        return 0, None, None, 0.0, np.zeros_like(mask, bool)
    sizes = ndimage.sum(np.ones_like(lab), lab, range(1, n + 1))
    k = int(np.argmax(sizes)) + 1
    bm = lab == k
    cy, cx = ndimage.center_of_mass(bm)
    return int(sizes[k - 1]), cy, cx, float(dist2coast[bm].mean()), bm


def build_ref_inland(root, n_sample=12):
    """Reference INLAND land domain from a complete tas store (union of finite)."""
    tas = [p for p in find_zarr(root) if parse(p)["var"] == REF_VAR]
    if not tas:
        return None, None, None
    ds = open_store(tas[0]); v, _ = main_var(ds, REF_VAR)
    T, LA, LO = ds[v].shape
    ref_land = np.zeros((LA, LO), bool)
    for i in np.unique(np.linspace(0, T - 1, n_sample).astype(int)):
        ref_land |= np.isfinite(np.asarray(ds[v].isel(time=int(i)).values))
    dist = ndimage.distance_transform_edt(ref_land)
    ref_inland = ref_land & (dist > EDGE_DIST)
    return ref_inland, ref_land, os.path.basename(tas[0])[:-5]


def analyse_store(path, ref_inland, max_dates, miss_min):
    m = parse(path)
    ds = open_store(path)
    var, _ = main_var(ds, m["var"])
    da = ds[var]; T, LA, LO = da.shape
    lat = np.asarray(ds.lat.values); lon = np.asarray(ds.lon.values)
    dates = np.array([f"{t:%Y-%m-%d}" for t in ds.time.values])
    use_ref = ref_inland is not None and ref_inland.shape == (LA, LO)

    finite = np.zeros((LA, LO), np.int64)
    miss = np.zeros(T, np.int64)                     # per-day inland cells missing vs reference
    for t0 in range(0, T, BLOCK):
        b = np.asarray(da.isel(time=slice(t0, t0 + BLOCK)).values)
        fin = np.isfinite(b)
        finite += fin.sum(0)
        if use_ref:
            miss[t0:t0 + b.shape[0]] = (~fin & ref_inland[None]).sum(axis=(1, 2))
    vf = finite / T

    land = vf >= OCEAN_VF
    dist2coast = ndimage.distance_transform_edt(land)
    coastal = int((land & (vf < BREATHE_HI) & (dist2coast <= EDGE_DIST)).sum())

    inland_day = miss >= miss_min
    n_inland_days = int(inland_day.sum())
    idx = np.where(inland_day)[0]
    if max_dates and len(idx) > max_dates:           # representative sample (avoid flooding)
        pick = np.unique(np.r_[np.linspace(0, len(idx) - 1, max_dates - 1).astype(int),
                               int(np.argmax(miss[idx]))])
        idx_emit = idx[pick]
    else:
        idx_emit = idx
    inland_dates = [dict(variable=var, model=m["model"], scenario=m["scenario"],
                         date=dates[i], issue=f"inland_missing({int(miss[i])}cells)")
                    for i in idx_emit]

    deficit = (ref_inland & (vf < PERSIST_LO)) if use_ref else np.zeros((LA, LO), bool)
    area, cy, cx, _, _ = classify_block(deficit, dist2coast)
    blk_lat = round(float(lat[int(round(cy))]), 2) if cy is not None else ""
    blk_lon = round(float(lon[int(round(cx))]), 2) if cx is not None else ""

    verdict = ("INLAND_CUTOFF" if n_inland_days > 0 else
               ("edge_wobble" if coastal else "clean"))
    row = dict(store=os.path.basename(path)[:-5], var=var, model=m["model"],
               scenario=m["scenario"], ref_used=use_ref, land_cells=int(land.sum()),
               inland_missing_days=n_inland_days,
               pct_days_missing=round(100 * n_inland_days / max(T, 1), 1),
               worst_missing_cells=int(miss.max()) if T else 0,
               deficit_block=area, block_lat=blk_lat, block_lon=blk_lon,
               coastal_cells=coastal, dates_emitted=len(inland_dates), verdict=verdict)
    return row, dict(vf=vf, lat=lat, lon=lon, deficit=deficit, land=land), inland_dates


def plot_vf(pdf, r, d):
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))
    im = a0.pcolormesh(d["lon"], d["lat"], np.where(d["land"], d["vf"], np.nan),
                       cmap="viridis", vmin=0, vmax=1, shading="auto")
    a0.set_aspect(1 / np.cos(np.deg2rad(float(d["lat"].mean()))))
    a0.set_title("valid fraction (share of days with data)")
    plt.colorbar(im, ax=a0, shrink=.8)
    a1.pcolormesh(d["lon"], d["lat"], np.where(d["land"], 0.15, np.nan),
                  cmap="Greys", vmin=0, vmax=1, shading="auto")
    a1.pcolormesh(d["lon"], d["lat"], np.where(d["deficit"], 1.0, np.nan),
                  cmap="autumn", vmin=0, vmax=1, shading="auto")
    a1.set_aspect(1 / np.cos(np.deg2rad(float(d["lat"].mean()))))
    a1.set_title(f"inland deficit vs tas (red)  block={r['deficit_block']} cells")
    fig.suptitle(f"[{r['verdict']}] {r['store']}   missing an inland chunk on "
                 f"{r['inland_missing_days']} days ({r['pct_days_missing']}%), "
                 f"worst {r['worst_missing_cells']} cells", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); pdf.savefig(fig); plt.close(fig)


def locate_bad_dates(path, dates, out_rows):
    m = parse(path); ds = open_store(path)
    var, _ = main_var(ds, m["var"]); vmin, vmax = BOUNDS.get(var, (None, None))
    tstr = np.array([f"{t:%Y-%m-%d}" for t in ds.time.values])
    samp = np.asarray(ds[var].isel(time=slice(0, min(60, ds.sizes["time"]))).values)
    land = np.isfinite(samp).any(0)
    dist2coast = ndimage.distance_transform_edt(land)
    for dt in dates:
        idx = np.where(tstr == dt)[0]
        if not idx.size:
            continue
        F = np.asarray(ds[var].isel(time=int(idx[0])).values)
        oob = np.isfinite(F)
        oob &= ((F < vmin) if vmin is not None else False) | ((F > vmax) if vmax is not None else False)
        if not oob.any():
            continue
        coastal = int((oob & (dist2coast <= EDGE_DIST)).sum())
        inland = int((oob & (dist2coast > EDGE_DIST)).sum())
        area, _, _, _, _ = classify_block(oob & (dist2coast > EDGE_DIST), dist2coast)
        out_rows.append(dict(variable=var, model=m["model"], scenario=m["scenario"],
                             date=dt, oob_cells=int(oob.sum()), coastal=coastal, inland=inland,
                             largest_inland_block=area,
                             verdict="inland_corruption" if area >= BLOCK_THRESH
                                     else ("edge_artifact" if coastal >= inland else "scattered_inland")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/local1/SETx_UIFL_Zarr")
    ap.add_argument("--sfcwind-root", default="/local1/SFCWIND_ZARR_POST")
    ap.add_argument("--vars", default="", help="comma list to restrict, e.g. hurs,rsds")
    ap.add_argument("--bad-dates", default="", help="bad_dates.csv to also locate oob cells")
    ap.add_argument("--inland-miss-min", type=int, default=200,
                    help="a day missing >= this many reference-inland cells is a bad date")
    ap.add_argument("--max-inland-dates", type=int, default=0,
                    help="max inland dates emitted PER STORE (0 = all, the default: "
                         "comprehensive but can be many); full count also reported")
    ap.add_argument("--out", default="footprint_out")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # reference inland domain from tas (spatially complete)
    ref_inland, ref_land, ref_name = build_ref_inland(args.root)
    if ref_inland is None:
        print("[ref] WARNING: no tas store found -> inland detection disabled", flush=True)
    else:
        print(f"[ref] inland domain from {ref_name}: "
              f"{int(ref_inland.sum())} inland cells (of {int(ref_land.sum())} land)", flush=True)

    stores = [p for p in find_zarr(args.root) if parse(p)["var"] not in ("sfcWind", "pr")]
    stores += find_zarr(args.sfcwind_root)
    if args.vars:
        keep = set(args.vars.split(","))
        stores = [p for p in stores if parse(p)["var"] in keep]
    print(f"[scan] {len(stores)} stores", flush=True)

    rows = []; inland_all = []
    pdf = PdfPages(os.path.join(args.out, "footprint_maps.pdf"))
    for i, p in enumerate(stores, 1):
        try:
            r, d, idates = analyse_store(p, ref_inland, args.max_inland_dates, args.inland_miss_min)
        except Exception as e:
            print(f"  ! {os.path.basename(p)}: {type(e).__name__}: {e}", flush=True); continue
        rows.append(r); inland_all.extend(idates)
        if r["verdict"] == "INLAND_CUTOFF":
            plot_vf(pdf, r, d)
        print(f"[{i}/{len(stores)}] {r['verdict']:14s} {r['store']}  "
              f"inland_missing_days={r['inland_missing_days']} ({r['pct_days_missing']}%) "
              f"emitted={r['dates_emitted']}", flush=True)
    pdf.close()

    with open(os.path.join(args.out, "footprint_classification.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    ipath = os.path.join(args.out, "inland_cutoff_dates.csv")
    with open(ipath, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["variable", "model", "scenario", "date", "issue"])
        for r in sorted(inland_all, key=lambda x: (x["variable"], x["model"], x["scenario"], x["date"])):
            w.writerow([r["variable"], r["model"], r["scenario"], r["date"], r["issue"]])
    print(f"[inland] {len(inland_all)} inland-cutoff dates (sampled) -> {ipath}", flush=True)

    if args.bad_dates and os.path.exists(args.bad_dates):
        want = {}
        with open(args.bad_dates) as f:
            for r in csv.DictReader(f):
                want.setdefault((r["variable"], r["model"], r["scenario"]), []).append(r["date"])
        loc = []
        pmap = {(parse(p)["var"], parse(p)["model"], parse(p)["scenario"]): p for p in stores}
        for key, dates in want.items():
            if key in pmap:
                try:
                    locate_bad_dates(pmap[key], dates, loc)
                except Exception as e:
                    print(f"  ! bad-date {key}: {type(e).__name__}", flush=True)
        if loc:
            with open(os.path.join(args.out, "bad_dates_located.csv"), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(loc[0].keys())); w.writeheader(); w.writerows(loc)
        print(f"[bad-dates] located {len(loc)} oob days -> bad_dates_located.csv", flush=True)

    ic = [r for r in rows if r["verdict"] == "INLAND_CUTOFF"]
    print(f"\n[done] {len(rows)} stores | INLAND_CUTOFF={len(ic)} "
          f"edge_wobble={sum(1 for r in rows if r['verdict']=='edge_wobble')} "
          f"clean={sum(1 for r in rows if r['verdict']=='clean')} -> {args.out}/", flush=True)
    for r in sorted(ic, key=lambda r: -r["inland_missing_days"])[:20]:
        print(f"    INLAND_CUTOFF {r['var']:6s} {r['model']:14s} {r['scenario']:11s} "
              f"missing {r['pct_days_missing']}% of days, worst {r['worst_missing_cells']} cells "
              f"@({r['block_lat']},{r['block_lon']})", flush=True)


if __name__ == "__main__":
    main()
