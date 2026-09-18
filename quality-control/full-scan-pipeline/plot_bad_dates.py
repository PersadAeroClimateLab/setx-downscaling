#!/usr/bin/env python3
"""
plot_bad_dates.py -- plot every date in bad_dates.csv (day-1 / day / day+1) into a
labeled directory for quick visual review, so you can decide clip vs re-downscale.

Author: Ifeanyichukwu C. Nduka (icnduka@jsg.utexas.edu)

For each bad date it saves  <out>/<variable>/<model>_<scenario>_<date>.png
with the three-panel neighbour view (NaN cells magenta, shared colour scale) and
the flagged issue(s) in the title. Stores are opened once and all their bad dates
plotted together, so this stays fast even for hundreds of dates.

Run:
  python3 plot_bad_dates.py --bad-dates bad_dates.csv \
      --root /local1/SETx_UIFL_Zarr --sfcwind-root /local1/SFCWIND_ZARR_POST
"""
import argparse, csv, os
from collections import defaultdict
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
    return (p[2] if len(p) > 2 else "?", p[0], p[1] if len(p) > 1 else "?")  # var, model, scen


def open_store(path):
    try:
        return xr.open_zarr(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True),
                            consolidated=False)
    except (AttributeError, TypeError):
        return xr.open_zarr(path, use_cftime=True, consolidated=False)


def main_var(ds, hint):
    d3 = [v for v in ds.data_vars if {"time", "lat", "lon"} <= set(ds[v].dims)]
    se = next((v for v in d3 if v.endswith("_SE") or "SE" in v), None)
    return hint if hint in d3 else next((v for v in d3 if v != se), None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bad-dates", default="bad_dates.csv")
    ap.add_argument("--root", default="/local1/SETx_UIFL_Zarr")
    ap.add_argument("--sfcwind-root", default="/local1/SFCWIND_ZARR_POST")
    ap.add_argument("--out", default="bad_date_plots")
    args = ap.parse_args()

    # (var, model, scenario) -> store path. Post overrides primary for sfcWind.
    paths = {}
    for p in find_zarr(args.root):
        paths[parse(p)] = p
    for p in find_zarr(args.sfcwind_root):
        paths[parse(p)] = p

    # group bad dates by store, merging issues per date
    by_store = defaultdict(lambda: defaultdict(list))
    with open(args.bad_dates) as f:
        for r in csv.DictReader(f):
            key = (r["variable"], r["model"], r["scenario"])
            by_store[key][r["date"]].append(r.get("issue", ""))
    n_dates = sum(len(v) for v in by_store.values())
    print(f"[plot] {n_dates} bad dates across {len(by_store)} stores -> {args.out}/", flush=True)

    made = miss = 0
    for key, datemap in sorted(by_store.items()):
        var, model, scen = key
        p = paths.get(key)
        if p is None:
            print(f"  ! no store for {key} ({len(datemap)} dates skipped)", flush=True)
            miss += len(datemap); continue
        ds = open_store(p); v = main_var(ds, var)
        alld = np.array([f"{t:%Y-%m-%d}" for t in ds.time.values])
        lat, lon = np.asarray(ds.lat.values), np.asarray(ds.lon.values)
        aspect = 1 / np.cos(np.deg2rad(float(lat.mean())))
        outdir = os.path.join(args.out, var); os.makedirs(outdir, exist_ok=True)
        cmap = plt.cm.viridis.copy(); cmap.set_bad("magenta")

        for date, issues in sorted(datemap.items()):
            w = np.where(alld == date)[0]
            if not w.size:
                miss += 1; continue
            c = int(w[0])
            idx = sorted(set(i for i in (c - 1, c, c + 1) if 0 <= i < len(alld)))
            fields = [np.asarray(ds[v].isel(time=i).values) for i in idx]
            fin = np.concatenate([F[np.isfinite(F)] for F in fields if np.isfinite(F).any()]) \
                if any(np.isfinite(F).any() for F in fields) else np.array([0.0, 1.0])
            vlo, vhi = np.percentile(fin, [2, 98])
            fig, axes = plt.subplots(1, len(idx), figsize=(4.3 * len(idx), 4.4), squeeze=False)
            im = None
            for ax, i, F in zip(axes[0], idx, fields):
                im = ax.pcolormesh(lon, lat, np.ma.masked_invalid(F), cmap=cmap,
                                   vmin=vlo, vmax=vhi, shading="auto")
                ax.set_aspect(aspect)
                ax.set_title(f"{alld[i]}" + ("  <-- flagged" if i == c else ""), fontsize=9)
            fig.colorbar(im, ax=axes[0].tolist(), shrink=.7,
                         label=f"{v}  (magenta = NaN)", location="right")
            fig.suptitle(f"{var}  {model}  {scen}  {date}   [{'; '.join(sorted(set(issues)))}]",
                         fontsize=11)
            fig.savefig(os.path.join(outdir, f"{model}_{scen}_{date}.png"),
                        dpi=110, bbox_inches="tight")
            plt.close(fig); made += 1
        print(f"  {var:7s} {model:14s} {scen:11s}  {len(datemap)} plots", flush=True)

    print(f"\n[done] {made} plots written, {miss} dates unresolved -> {args.out}/", flush=True)


if __name__ == "__main__":
    main()
