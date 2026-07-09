"""
Read-only reduction analysis for step-23 near-duplicate folding (transfer_pooling).

For every pathogen/category, quantify how much the pre-transfer near-duplicate fold reduces the
dataset count, and sweep the (overlap, discordance) cutoffs to inform the 90%/10% choice. Folding is
cheap (set operations, no modeling), and the expensive part (per-pair overlap + discordance) is
computed ONCE per category, then re-thresholded across the sweep.

Overlap = |A n B| / min(|A|,|B|) (containment). Edge iff overlap >= o AND shared-label discordance < d.
Groups = connected components. "Star" group = every member is >= o contained in the largest member
(benign: subsets of one library); a non-star group required transitive links (possible chaining).

Outputs (docs/):
    2026-07-02_neardup-reduction.csv     per pathogen/category x (o,d): datasets before/after, groups
    2026-07-02_neardup-reduction.png     figure (reduction at 0.90/0.10; sweep)
    prints the 0.90/0.10 table + sweep to stdout.

Usage:
    python tools/near_dup_report.py                 # all pathogens
    python tools/near_dup_report.py saureus mtuberculosis
"""
import os
import sys
import glob
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from pooling_utils import load_binarised_datasets, _neardup_edges

stylia.set_format("slide")
stylia.set_style("ersilia")

output_base = os.path.join(root, "..", "output")
docs_dir = os.path.join(root, "..", "docs")
DATE = "2026-07-02"
MAXPOST = 100
OVERLAPS = [0.80, 0.90, 0.95]
DISCORDS = [0.05, 0.10, 0.20]
LOCK = (0.90, 0.10)


def components(n, edges):
    if edges:
        r = [e[0] for e in edges] + [e[1] for e in edges]
        c = [e[1] for e in edges] + [e[0] for e in edges]
        adj = csr_matrix((np.ones(len(r)), (r, c)), shape=(n, n))
    else:
        adj = csr_matrix((n, n))
    _, comp = connected_components(adj, directed=False)
    groups = defaultdict(list)
    for i, ci in enumerate(comp):
        groups[ci].append(i)
    return list(groups.values())


def analyse(pathogen, category):
    ds, _ = load_binarised_datasets(os.path.join(output_base, "stage4", pathogen, "22_binarised", category))
    if not ds:
        return []
    n = len(ds)
    iksets = [d["ikset"] for d in ds]
    # all candidate pairs with overlap >= min(OVERLAPS); (i,j,overlap,discord)
    allpairs = _neardup_edges(ds, min(OVERLAPS), 1.01, MAXPOST)
    rows = []
    for o in OVERLAPS:
        for d in DISCORDS:
            edges = [(i, j) for (i, j, ov, dc) in allpairs if ov >= o and dc < d]
            groups = components(n, edges)
            sizes = sorted((len(g) for g in groups), reverse=True)
            n_after = len(groups)
            multi = [g for g in groups if len(g) > 1]
            # star test per multi group: all members >= o contained in the largest member
            n_nonstar = 0
            for g in multi:
                big = max(g, key=lambda i: len(iksets[i]))
                if not all(len(iksets[i] & iksets[big]) / min(len(iksets[i]), len(iksets[big])) >= o
                           for i in g if i != big):
                    n_nonstar += 1
            rows.append({"pathogen": pathogen, "category": category, "overlap": o, "discord": d,
                         "n_datasets": n, "n_after": n_after, "n_groups": len(multi),
                         "max_group": sizes[0] if sizes else 0, "n_nonstar_groups": n_nonstar})
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("pathogens", nargs="*")
    args = ap.parse_args()
    pathogens = args.pathogens or sorted(
        os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(output_base, "stage4", "*", "22_binarised")))

    all_rows = []
    for p in pathogens:
        for cat in ["DR", "SP"]:
            r = analyse(p, cat)
            all_rows += r
        print(f"  {p} done")
    df = pd.DataFrame(all_rows)
    os.makedirs(docs_dir, exist_ok=True)
    df.to_csv(os.path.join(docs_dir, f"{DATE}_neardup-reduction.csv"), index=False)

    lock = df[(df.overlap == LOCK[0]) & (df.discord == LOCK[1])].copy()
    lock["pct_reduction"] = (100 * (1 - lock.n_after / lock.n_datasets)).round(1)
    pd.set_option("display.width", 200); pd.set_option("display.max_rows", 60)
    print(f"\n=== Reduction at overlap={LOCK[0]}, discord={LOCK[1]} ===")
    print(lock[["pathogen", "category", "n_datasets", "n_after", "pct_reduction",
                "n_groups", "max_group", "n_nonstar_groups"]].to_string(index=False))
    tot = lock.groupby([]).agg(n=("n_datasets", "sum"), a=("n_after", "sum"))
    print(f"\nTOTAL datasets {int(lock.n_datasets.sum())} -> {int(lock.n_after.sum())} "
          f"({100*(1-lock.n_after.sum()/lock.n_datasets.sum()):.1f}% reduction); "
          f"non-star groups: {int(lock.n_nonstar_groups.sum())}")
    print("\n=== Sweep: total datasets-after across all pathogens (o rows x d cols) ===")
    piv = df.groupby(["overlap", "discord"]).n_after.sum().unstack("discord")
    print(piv.to_string())
    print(f"\nwrote {os.path.join(docs_dir, DATE + '_neardup-reduction.csv')}")
