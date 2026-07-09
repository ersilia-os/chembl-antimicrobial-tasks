"""
Read-only sweep to choose SINGLETON_MIN_CHEMSPACE_FRAC for step 24.

A step-24 singleton is promoted to a modelled pool only if it covers a meaningful fraction of the
pathogen/category unique chemical space. This report shows, per pathogen/category and per candidate
threshold (1% / 2% / 5%), how many singletons qualify and how much of the chem space they carry -
so the cutoff can be picked empirically. Reads only step-24 outputs (no modeling).

Inputs:  output/stage4/<p>/23_pool_summary.csv, output/stage4/<p>/23_chemspace_partition.csv
Outputs: docs/2026-07-02_singleton-inclusion.csv + stdout table

Usage:   python tools/singleton_inclusion_report.py
"""
import os
import glob

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
output_base = os.path.join(root, "..", "output")
docs_dir = os.path.join(root, "..", "docs")
THRESHOLDS = [0.01, 0.02, 0.05]
CATEGORIES = ["DR", "SP"]


def main():
    rows = []
    codes = sorted(os.path.basename(os.path.dirname(f))
                   for f in glob.glob(os.path.join(output_base, "stage4", "*", "23_pool_summary.csv")))
    for p in codes:
        summ = pd.read_csv(os.path.join(output_base, "stage4", p, "23_pool_summary.csv"))
        part = pd.read_csv(os.path.join(output_base, "stage4", p, "23_chemspace_partition.csv"))
        ntot = {r.category: int(r.n_total) for r in part.itertuples(index=False)}
        for cat in CATEGORIES:
            sub = summ[summ["category"] == cat]
            if sub.empty:
                continue
            total = ntot.get(cat, 0)
            multi = sub[~sub["is_singleton"]]
            sing = sub[sub["is_singleton"]].copy()
            sing["frac"] = sing["n_molecules"] / total if total else 0.0
            row = {"pathogen": p, "category": cat, "chemspace": total,
                   "n_multi_pools": len(multi), "n_singletons": len(sing)}
            for t in THRESHOLDS:
                q = sing[sing["frac"] >= t]
                row[f"n_sing_ge_{int(t*100)}pct"] = len(q)
                row[f"cov_ge_{int(t*100)}pct"] = round(q["n_molecules"].sum() / total, 4) if total else 0.0
            rows.append(row)
    df = pd.DataFrame(rows)
    os.makedirs(docs_dir, exist_ok=True)
    df.to_csv(os.path.join(docs_dir, "2026-07-02_singleton-inclusion.csv"), index=False)
    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 60)
    print(df.to_string(index=False))
    print("\nn_sing_ge_Xpct = singletons covering >= X% of chem space (would be modelled);"
          " cov_ge_Xpct = their summed molecules as fraction of chem space.")
    print(f"\nwrote {os.path.join(docs_dir, '2026-07-02_singleton-inclusion.csv')}")


if __name__ == "__main__":
    main()
