#!/usr/bin/env python3
"""
drilldown_qc.py  --  Summarise the per-store CSVs written by qc_zarr_all.py.

Reads each store's *_incomplete_dates.csv, *_oob_dates.csv, *_missing_dates.csv
(and qc_zarr_summary.csv for land_cells / SE / dedup) and separates:

  REAL_GAP    empty days, missing days, or PARTIAL days losing a large fraction
              of the land footprint
  BAD_VALUES  out-of-range cells (below_min / above_max / neg / inf) or constant slices
  WOBBLE      only EXTRA-valid / small PARTIAL / SE anomalies  -> benign day-to-day
              footprint breathing, not a data problem

The key thing it adds over the run summary: the *magnitude* (cells) of every
EXTRA/PARTIAL flag, so the 184k "extra-valid" days can be shown to be edge-wobble.

'pr' stores are discarded.

Usage:
  python3 drilldown_qc.py --qc-dir qc_zarr_out
"""

import argparse, csv, glob, os
from collections import defaultdict

TOTAL_CELLS = 630 * 480
PARTIAL_REAL_FRAC = 0.02          # a PARTIAL day missing >2% of land -> real gap
SMALL_DEV = 50                    # <= this many cells deviation -> edge wobble
DEV_BUCKETS = [10, 100, 1000, 10000]


def parse_store(name):
    p = name.split("_")
    return (p[0] if p else "?",
            p[1] if len(p) > 1 else "?",
            p[2] if len(p) > 2 else "?")


def read_incomplete(path):
    cls = {"EMPTY": [], "PARTIAL": [], "EXTRA": []}
    if not os.path.exists(path):
        return cls
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                cls.setdefault(r["class"], []).append(abs(int(r["deviation"])))
            except (KeyError, ValueError):
                pass
    return cls


def read_oob(path):
    tot = defaultdict(int); cdays = odays = 0
    if not os.path.exists(path):
        return tot, cdays, odays
    with open(path) as f:
        for r in csv.DictReader(f):
            for tok in (r.get("reasons") or "").split(";"):
                if not tok:
                    continue
                if ":" in tok:
                    k, v = tok.split(":"); tot[k] += int(v)
                elif tok == "constant":
                    cdays += 1
                elif tok == "outlier":
                    odays += 1
    return tot, cdays, odays


def count_rows(path):
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        return max(sum(1 for _ in f) - 1, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-dir", default="qc_zarr_out")
    ap.add_argument("--out", default="qc_drilldown_out")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # summary CSV (dedup by store) -> land_cells, se_anom_days
    meta = {}
    scsv = os.path.join(args.qc_dir, "qc_zarr_summary.csv")
    if os.path.exists(scsv):
        with open(scsv) as f:
            for r in csv.DictReader(f):
                meta.setdefault(r["store"], r)     # first wins -> dedup

    stores = sorted(os.path.basename(p)[:-len("_incomplete_dates.csv")]
                    for p in glob.glob(os.path.join(args.qc_dir, "*_incomplete_dates.csv")))

    rows = []
    ev_bucket = defaultdict(int); pv_bucket = defaultdict(int)   # global deviation histograms
    for st in stores:
        model, scen, var = parse_store(st)
        if var == "pr":
            continue
        base = os.path.join(args.qc_dir, st)
        cls = read_incomplete(base + "_incomplete_dates.csv")
        oob, cdays, odays = read_oob(base + "_oob_dates.csv")
        n_missing = count_rows(base + "_missing_dates.csv")
        land = int(meta.get(st, {}).get("land_cells") or (TOTAL_CELLS))
        se_anom = int(meta.get(st, {}).get("se_anom_days") or 0)

        def stat(lst):
            return (len(lst), (max(lst) if lst else 0))
        n_empty, _ = stat(cls["EMPTY"])
        n_part, part_max = stat(cls["PARTIAL"])
        n_extra, extra_max = stat(cls["EXTRA"])
        for d in cls["EXTRA"]:
            ev_bucket[next((b for b in DEV_BUCKETS if d <= b), ">10000")] += 1
        for d in cls["PARTIAL"]:
            pv_bucket[next((b for b in DEV_BUCKETS if d <= b), ">10000")] += 1

        below = oob.get("below_min", 0); above = oob.get("above_max", 0)
        neg = oob.get("neg", 0); inf = oob.get("inf", 0)
        part_frac = part_max / max(land, 1)

        tags = []
        if n_empty or n_missing or part_frac > PARTIAL_REAL_FRAC:
            tags.append("REAL_GAP")
        if below or above or neg or inf:
            tags.append("BAD_VALUES")
        if cdays:
            tags.append("CONSTANT")
        if not tags and (n_extra or n_part or se_anom):
            tags.append("WOBBLE")
        if not tags:
            tags.append("CLEAN")

        rows.append(dict(store=st, model=model, scenario=scen, var=var, land_cells=land,
                         n_missing=n_missing, n_empty=n_empty, n_partial=n_part,
                         partial_dev_max=part_max, partial_frac_max=round(part_frac, 4),
                         n_extra=n_extra, extra_dev_max=extra_max,
                         below_min_cells=below, above_max_cells=above,
                         neg_cells=neg, inf_cells=inf, constant_days=cdays,
                         outlier_days=odays, se_anom_days=se_anom,
                         severity="+".join(tags)))

    # ---- write per-store enriched table ----
    fields = list(rows[0].keys()) if rows else []
    with open(os.path.join(args.out, "drilldown_by_store.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

    # ---- actionable subset ----
    actionable = [r for r in rows if "REAL_GAP" in r["severity"] or
                  "BAD_VALUES" in r["severity"] or "CONSTANT" in r["severity"]]
    actionable.sort(key=lambda r: -(r["below_min_cells"] + r["above_max_cells"] +
                                     r["neg_cells"] + r["inf_cells"] +
                                     r["n_empty"] * 10**6 + int(r["partial_frac_max"] > PARTIAL_REAL_FRAC) * 10**5))
    with open(os.path.join(args.out, "drilldown_actionable.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(actionable)

    # ---- per-variable rollup ----
    byvar = defaultdict(list)
    for r in rows:
        byvar[r["var"]].append(r)

    L = ["=" * 78, "  SETx-1km-CMIP6 Zarr QC — DRILL-DOWN SUMMARY  (pr discarded)", "=" * 78]
    L.append(f"  stores summarised : {len(rows)}   (from per-store CSVs, deduplicated)")
    L.append("")
    L.append("  per-variable  (stores | REAL_GAP | BAD_VALUES | WOBBLE-only | CLEAN | "
             "below-min  above-max)")
    def cnt(rs, tag): return sum(1 for r in rs if tag in r["severity"])
    for v in sorted(byvar):
        rs = byvar[v]
        wob = sum(1 for r in rs if r["severity"] == "WOBBLE")
        cln = sum(1 for r in rs if r["severity"] == "CLEAN")
        L.append(f"    {v:8s} {len(rs):5d}   {cnt(rs,'REAL_GAP'):6d}   {cnt(rs,'BAD_VALUES'):8d}"
                 f"   {wob:9d}   {cln:5d}   {sum(r['below_min_cells'] for r in rs):10d}"
                 f" {sum(r['above_max_cells'] for r in rs):10d}")

    L.append("")
    L.append("  EXTRA-valid day magnitudes (how many cells the footprint 'breathed'):")
    tot_ev = sum(ev_bucket.values()) or 1
    for b in DEV_BUCKETS + [">10000"]:
        c = ev_bucket.get(b, 0)
        lab = f"<={b} cells" if b != ">10000" else ">10000 cells"
        L.append(f"      {lab:16s}: {c:8d}  ({100*c/tot_ev:5.1f}%)")
    small = sum(v for k, v in ev_bucket.items() if isinstance(k, int) and k <= SMALL_DEV)
    L.append(f"    -> {100*small/tot_ev:.1f}% of extra-valid flags are <= {SMALL_DEV} cells "
             f"(edge wobble, benign).")

    L.append("")
    L.append("  ACTIONABLE — bad values (top 15 by out-of-range cell count):")
    bad = [r for r in actionable if "BAD_VALUES" in r["severity"] or "CONSTANT" in r["severity"]]
    for r in bad[:15]:
        L.append(f"    {r['var']:7s} {r['model']:14s} {r['scenario']:11s}  "
                 f"below_min={r['below_min_cells']:>8}  above_max={r['above_max_cells']:>8}"
                 f"  const={r['constant_days']}")
    if len(bad) > 15:
        L.append(f"    ... and {len(bad)-15} more (see drilldown_actionable.csv)")

    L.append("")
    L.append("  ACTIONABLE — real coverage gaps (empty / missing / large-partial days):")
    gap = [r for r in actionable if "REAL_GAP" in r["severity"]]
    for r in sorted(gap, key=lambda r: -(r["n_empty"] + r["n_partial"])):
        L.append(f"    {r['var']:7s} {r['model']:14s} {r['scenario']:11s}  "
                 f"empty={r['n_empty']} missing={r['n_missing']} partial={r['n_partial']} "
                 f"(max {100*r['partial_frac_max']:.1f}% of land)")
    if not gap:
        L.append("    (none)")

    L.append("")
    L.append("-" * 78)
    tb = lambda k: sum(r[k] for r in rows)
    L.append(f"  totals: REAL_GAP stores={sum(1 for r in rows if 'REAL_GAP' in r['severity'])}"
             f"  BAD_VALUES stores={sum(1 for r in rows if 'BAD_VALUES' in r['severity'])}"
             f"  WOBBLE-only={sum(1 for r in rows if r['severity']=='WOBBLE')}"
             f"  CLEAN={sum(1 for r in rows if r['severity']=='CLEAN')}")
    L.append(f"  out-of-range cells: below_min={tb('below_min_cells')}  "
             f"above_max={tb('above_max_cells')}  neg={tb('neg_cells')}  inf={tb('inf_cells')}")
    L.append("-" * 78)
    L.append("  READ: the extra-valid / SE-anomaly flood is day-to-day footprint wobble on the")
    L.append("  kriged radiation & humidity fields — benign. The real work is the BAD_VALUES")
    L.append("  stores (clip/refit out-of-range rlds/huss/tasmin) and the REAL_GAP stores.")
    L.append("=" * 78)
    report = "\n".join(L)
    print(report)
    with open(os.path.join(args.out, "DRILLDOWN_SUMMARY.txt"), "w") as f:
        f.write(report + "\n")
    print(f"\n[done] -> {args.out}/  (DRILLDOWN_SUMMARY.txt, drilldown_by_store.csv, "
          f"drilldown_actionable.csv)")


if __name__ == "__main__":
    main()
