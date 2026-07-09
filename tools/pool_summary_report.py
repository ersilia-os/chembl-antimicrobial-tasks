"""
Read-only cross-pathogen summary of the FINAL pools (after step 25).

For every pathogen/category, reports the pools we have, their chemical-space coverage (union of pool
molecules / category unique chem space), their activity ratios, and the cross-validated AUROC of each
pool's model. The final AUROC per pool is the grown-pool CV AUROC when its step-25 merge batch was
accepted, else the step-24 baseline (pool unchanged).

Inputs (output/stage4/<pathogen>/): 25_pool_summary.csv, 25_pools/<cat>/*.csv.gz, 23_chemspace_partition.csv
Outputs: docs/2026-07-03_final-pool-summary.csv (per pool) + docs/2026-07-03_final-pool-bycategory.csv
         and a printed per-pathogen/category table.

Usage: python tools/pool_summary_report.py
"""
import os
import glob

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
output_base = os.path.join(root, "..", "output")
docs_dir = os.path.join(root, "..", "docs")
DATE = "2026-07-03"
CATEGORIES = ["DR", "SP"]


def main():
    codes = sorted(os.path.basename(os.path.dirname(f))
                   for f in glob.glob(os.path.join(output_base, "stage4", "*", "25_pool_summary.csv")))
    perpool, bycat = [], []
    for p in codes:
        s = pd.read_csv(os.path.join(output_base, "stage4", p, "25_pool_summary.csv"))
        part = pd.read_csv(os.path.join(output_base, "stage4", p, "23_chemspace_partition.csv"))
        ntot = {r.category: int(r.n_total) for r in part.itertuples(index=False)}
        # final per-pool AUROC = grown if accepted else baseline
        s = s.copy()
        s["auroc"] = np.where(s.batch_accepted, s.grown_auroc, s.baseline_auroc)
        s["n_molecules"] = s.n_mol_after
        s["active_ratio"] = s.ar_after
        for cat in CATEGORIES:
            sc = s[s.category == cat]
            if sc.empty:
                continue
            union = set()
            for f in glob.glob(os.path.join(output_base, "stage4", p, "25_pools", cat, "*.csv.gz")):
                union |= set(pd.read_csv(f, usecols=["inchikey"]).inchikey)
            total = ntot.get(cat, 0)
            for r in sc.itertuples(index=False):
                perpool.append({"pathogen": p, "category": cat, "pool_id": r.pool_id,
                                "n_molecules": int(r.n_molecules), "active_ratio": round(r.active_ratio, 4),
                                "auroc": round(r.auroc, 4), "merge_accepted": bool(r.batch_accepted)})
            bycat.append({"pathogen": p, "category": cat, "n_pools": len(sc),
                          "coverage_pct": round(100 * len(union) / total, 1) if total else np.nan,
                          "n_mol_union": len(union), "chemspace": total,
                          "active_ratio_median": round(sc.active_ratio.median(), 3),
                          "active_ratio_min": round(sc.active_ratio.min(), 3),
                          "active_ratio_max": round(sc.active_ratio.max(), 3),
                          "auroc_median": round(sc.auroc.median(), 3),
                          "auroc_min": round(sc.auroc.min(), 3),
                          "auroc_max": round(sc.auroc.max(), 3)})
    # include kept step-28 low-data catch-all pools (a catch-all pools the whole category => 100% coverage)
    for p in sorted(os.path.basename(os.path.dirname(f))
                    for f in glob.glob(os.path.join(output_base, "stage4", "*", "26_cv_summary.csv"))):
        c = pd.read_csv(os.path.join(output_base, "stage4", p, "26_cv_summary.csv"))
        for r in c[c.kept].itertuples(index=False):
            perpool.append({"pathogen": p, "category": r.category, "pool_id": f"{r.category}_catchall",
                            "n_molecules": int(r.n_molecules), "active_ratio": round(r.active_ratio, 4),
                            "auroc": round(r.auroc, 4), "merge_accepted": True})
            bycat.append({"pathogen": p, "category": r.category, "n_pools": 1, "coverage_pct": 100.0,
                          "n_mol_union": int(r.n_molecules), "chemspace": int(r.n_molecules),
                          "active_ratio_median": round(r.active_ratio, 3), "active_ratio_min": round(r.active_ratio, 3),
                          "active_ratio_max": round(r.active_ratio, 3), "auroc_median": round(r.auroc, 3),
                          "auroc_min": round(r.auroc, 3), "auroc_max": round(r.auroc, 3)})

    pp = pd.DataFrame(perpool); bc = pd.DataFrame(bycat).sort_values(["pathogen", "category"])
    os.makedirs(docs_dir, exist_ok=True)
    pp.to_csv(os.path.join(docs_dir, f"{DATE}_final-pool-summary.csv"), index=False)
    bc.to_csv(os.path.join(docs_dir, f"{DATE}_final-pool-bycategory.csv"), index=False)

    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 60)
    print("=== FINAL POOLS by pathogen x category (after step 25) ===")
    print(bc[["pathogen", "category", "n_pools", "coverage_pct", "active_ratio_median",
              "active_ratio_min", "active_ratio_max", "auroc_median", "auroc_min", "auroc_max"]].to_string(index=False))
    print(f"\n=== TOTALS ===")
    print(f"pools: {len(pp)}  (DR {int((pp.category=='DR').sum())}, SP {int((pp.category=='SP').sum())})")
    for cat in CATEGORIES:
        c = pp[pp.category == cat]
        print(f"  {cat}: pools={len(c)}  AUROC median={c.auroc.median():.3f} "
              f"(>=0.7: {100*(c.auroc>=0.7).mean():.0f}%, >=0.8: {100*(c.auroc>=0.8).mean():.0f}%)  "
              f"active_ratio median={c.active_ratio.median():.3f}")
    bc2 = bc.copy()
    print(f"\nmean coverage: DR {bc2[bc2.category=='DR'].coverage_pct.mean():.1f}%  "
          f"SP {bc2[bc2.category=='SP'].coverage_pct.mean():.1f}%")
    print(f"\nwrote docs/{DATE}_final-pool-summary.csv (per pool) + docs/{DATE}_final-pool-bycategory.csv")


if __name__ == "__main__":
    main()
