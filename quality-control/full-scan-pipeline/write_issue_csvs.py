#!/usr/bin/env python3
"""
write_issue_csvs.py  --  One issues CSV per variable from qc_zarr_all.py output.

Reads each store's *_incomplete_dates.csv and *_oob_dates.csv under --qc-dir and
emits, per variable, a flat list of every flagged (date, issue):

  columns: variable, model, scenario, date, issue, category

category lets you filter instantly:
  wobble      EXTRA-valid days and small PARTIAL days (benign footprint breathing)
  real_gap    EMPTY, MISSING, or PARTIAL losing > 2% of land
  bad_values  out-of-range cells (below_min / above_max / neg / inf) or constant slice
  outlier     domain-mean temporal outlier

'pr' is discarded.  Output: <out>/<variable>_issues.csv

Run:
  python3 write_issue_csvs.py --qc-dir qc_zarr_out --out issues_by_variable
"""

import argparse, csv, glob, os
from collections import defaultdict

TOTAL_CELLS = 630 * 480
PARTIAL_REAL_FRAC = 0.02


def parse_store(name):
    p = name.split("_")
    return (p[0] if p else "?",
            p[1] if len(p) > 1 else "?",
            p[2] if len(p) > 2 else "?")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-dir", default="qc_zarr_out")
    ap.add_argument("--out", default="issues_by_variable")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # land_cells per store (dedup) for the partial threshold
    land = {}
    scsv = os.path.join(args.qc_dir, "qc_zarr_summary.csv")
    if os.path.exists(scsv):
        with open(scsv) as f:
            for r in csv.DictReader(f):
                land.setdefault(r["store"], int(r.get("land_cells") or TOTAL_CELLS))

    by_var = defaultdict(list)
    stores = sorted(os.path.basename(p)[:-len("_incomplete_dates.csv")]
                    for p in glob.glob(os.path.join(args.qc_dir, "*_incomplete_dates.csv")))

    for st in stores:
        model, scen, var = parse_store(st)
        if var == "pr":
            continue
        lc = land.get(st, TOTAL_CELLS)

        # incomplete: EMPTY / PARTIAL / EXTRA
        ip = os.path.join(args.qc_dir, st + "_incomplete_dates.csv")
        if os.path.exists(ip):
            with open(ip) as f:
                for r in csv.DictReader(f):
                    cls = r.get("class", ""); dev = int(r.get("deviation") or 0)
                    if cls == "EMPTY":
                        issue, cat = "EMPTY", "real_gap"
                    elif cls == "PARTIAL":
                        issue = f"PARTIAL(missing={abs(dev)}cells)"
                        cat = "real_gap" if abs(dev) / max(lc, 1) > PARTIAL_REAL_FRAC else "wobble"
                    elif cls == "EXTRA":
                        issue, cat = f"EXTRA(extra={abs(dev)}cells)", "wobble"
                    else:
                        continue
                    by_var[var].append((var, model, scen, r["date"], issue, cat))

        # missing days
        mp = os.path.join(args.qc_dir, st + "_missing_dates.csv")
        if os.path.exists(mp):
            with open(mp) as f:
                for r in csv.DictReader(f):
                    by_var[var].append((var, model, scen, r["date"], "MISSING", "real_gap"))

        # out-of-range / constant / outlier
        op = os.path.join(args.qc_dir, st + "_oob_dates.csv")
        if os.path.exists(op):
            with open(op) as f:
                for r in csv.DictReader(f):
                    reasons = r.get("reasons") or ""
                    if not reasons:
                        continue
                    kinds = {t.split(":")[0] for t in reasons.split(";") if t}
                    if (kinds & {"below_min", "above_max", "neg", "inf", "constant"}
                            or any(k.startswith("incoherent") for k in kinds)):
                        cat = "bad_values"
                    elif "outlier" in kinds:
                        cat = "outlier"
                    else:
                        cat = "bad_values"
                    by_var[var].append((var, model, scen, r["date"], reasons, cat))

    # write per-variable CSVs (sorted model, scenario, date)
    print(f"{'variable':8s}  rows   bad_values  real_gap  outlier  wobble")
    for var in sorted(by_var):
        rows = sorted(by_var[var], key=lambda x: (x[1], x[2], x[3]))
        path = os.path.join(args.out, f"{var}_issues.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["variable", "model", "scenario", "date", "issue", "category"])
            w.writerows(rows)
        c = defaultdict(int)
        for x in rows:
            c[x[5]] += 1
        print(f"{var:8s}  {len(rows):6d}  {c['bad_values']:9d}  {c['real_gap']:8d}  "
              f"{c['outlier']:7d}  {c['wobble']:6d}   -> {var}_issues.csv")

    print(f"\n[done] -> {args.out}/  (one *_issues.csv per variable; filter on 'category')")


if __name__ == "__main__":
    main()
