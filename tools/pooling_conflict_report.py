"""
Read-only diagnostic report for step 23 pools.

For every pool built by step 23, this recomputes — WITHOUT modifying any pipeline
output — three things the user asked about:

  1. how many pools were generated (per pathogen / category, multi-dataset vs singleton);
  2. how many CONFLICTS exist "before handling", i.e. molecules that are tested in >= 2 of
     a pool's member datasets and whose per-dataset binary labels DISAGREE (both an active
     and an inactive call among the raw 22_binarised bins). Singleton pools cannot have
     conflicts (one dataset, one row per molecule);
  3. of those conflicted molecules, what % resolve to ACTIVE vs INACTIVE after step 23's
     value-resolution handling (average within (unit, direction) + re-binarise). The
     resolved label is read directly from the merged pool file, so the resolution matches
     exactly what step 23 wrote — this script does not re-implement the averaging.

Definition note: a conflict is measured on the STORED 22_binarised bins (each dataset's own
expert cutoff), which is the datasets' state "before handling". Concordant molecules (tested
multiple times but agreeing) are counted separately for context.

Inputs (output/stage4/<pathogen>/):
    23_pool_members.csv                              pool_id -> member dataset_id
    23_pool_summary.csv                              per-pool sizes / active ratio
    22_binarised/{DR,SP}/<dataset_id>.csv.gz         inchikey, bin  (raw per-dataset labels)
    23_pools/<cat>/<pool_id>.csv.gz                  inchikey, bin  (resolved labels)

Outputs (docs/):
    2026-07-01_step23-conflict-perpool.csv           one row per pool
    2026-07-01_step23-conflict-summary.csv           one row per pathogen x category
    (also prints the summary to stdout)

Usage:
    python tools/pooling_conflict_report.py                  # all pathogens with 25 outputs
    python tools/pooling_conflict_report.py saureus mtuberculosis
"""

import os
import sys
import glob
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

output_base = os.path.join(root, "..", "output")
docs_dir = os.path.join(root, "..", "docs")
os.makedirs(docs_dir, exist_ok=True)

CATEGORIES = ["DR", "SP"]
DATE = "2026-07-01"


def _read_bins(path):
    """inchikey -> bin (0/1) for a binarised or pool file; first row per ik wins."""
    df = pd.read_csv(path, usecols=["inchikey", "bin"]).drop_duplicates("inchikey")
    return dict(zip(df["inchikey"], df["bin"].astype(int)))


def analyse_pool(pathogen, category, pool_id, member_ids):
    """Per-molecule conflict accounting for one pool."""
    bdir = os.path.join(output_base, "stage4", pathogen, "22_binarised", category)

    # collect the raw per-dataset bins each molecule received across member datasets
    per_mol = defaultdict(list)               # inchikey -> [bin, bin, ...]
    for dsid in member_ids:
        f = os.path.join(bdir, dsid + ".csv.gz")
        if not os.path.exists(f):
            continue
        for ik, b in _read_bins(f).items():
            per_mol[ik].append(b)

    # resolved labels as written by step 23
    pool_file = os.path.join(output_base, "stage4", pathogen, "23_pools", category, pool_id + ".csv.gz")
    resolved = _read_bins(pool_file) if os.path.exists(pool_file) else {}

    n_multitested = n_concordant = n_conflict = 0
    n_conf_active = n_conf_inactive = n_conf_missing = 0
    for ik, bins in per_mol.items():
        if len(bins) < 2:
            continue
        n_multitested += 1
        if min(bins) == max(bins):
            n_concordant += 1
            continue
        n_conflict += 1                        # disagreeing raw calls -> a conflict
        rb = resolved.get(ik)
        if rb is None:
            n_conf_missing += 1
        elif rb == 1:
            n_conf_active += 1
        else:
            n_conf_inactive += 1

    n_mol = len(resolved) if resolved else len(per_mol)
    return {
        "pathogen": pathogen, "category": category, "pool_id": pool_id,
        "n_datasets": len(member_ids), "n_molecules": n_mol,
        "n_multitested": n_multitested, "n_concordant": n_concordant,
        "n_conflict": n_conflict,
        "n_conflict_active": n_conf_active, "n_conflict_inactive": n_conf_inactive,
        "n_conflict_unresolved": n_conf_missing,
        "pct_conflict_of_pool": round(100 * n_conflict / n_mol, 3) if n_mol else 0.0,
        "pct_conflict_active": round(100 * n_conf_active / n_conflict, 2) if n_conflict else np.nan,
        "pct_conflict_inactive": round(100 * n_conf_inactive / n_conflict, 2) if n_conflict else np.nan,
    }


def process(pathogen):
    mem_path = os.path.join(output_base, "stage4", pathogen, "23_pool_members.csv")
    if not os.path.exists(mem_path):
        return []
    members = pd.read_csv(mem_path)
    rows = []
    for (category, pool_id), grp in members.groupby(["category", "pool_id"]):
        rows.append(analyse_pool(pathogen, category, pool_id, grp["dataset_id"].tolist()))
    return rows


def summarise(perpool):
    """Aggregate per (pathogen, category)."""
    out = []
    for (p, c), g in perpool.groupby(["pathogen", "category"]):
        conf = int(g["n_conflict"].sum())
        act = int(g["n_conflict_active"].sum())
        ina = int(g["n_conflict_inactive"].sum())
        out.append({
            "pathogen": p, "category": c,
            "n_pools": len(g),
            "n_multi_dataset_pools": int((g["n_datasets"] > 1).sum()),
            "n_singleton_pools": int((g["n_datasets"] == 1).sum()),
            "total_molecules": int(g["n_molecules"].sum()),
            "total_multitested": int(g["n_multitested"].sum()),
            "total_conflicts": conf,
            "conflicts_active": act, "conflicts_inactive": ina,
            "conflicts_unresolved": int(g["n_conflict_unresolved"].sum()),
            "pct_conflict_active": round(100 * act / conf, 2) if conf else np.nan,
            "pct_conflict_inactive": round(100 * ina / conf, 2) if conf else np.nan,
        })
    return pd.DataFrame(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Step 24 pooling conflict report")
    ap.add_argument("pathogens", nargs="*", help="default: every pathogen with 23_pool_members.csv")
    args = ap.parse_args()

    pathogens = args.pathogens or sorted(
        os.path.basename(os.path.dirname(p))
        for p in glob.glob(os.path.join(output_base, "stage4", "*", "23_pool_members.csv")))

    all_rows = []
    for p in pathogens:
        r = process(p)
        if r:
            all_rows += r
            print(f"  {p}: {len(r)} pools analysed")
        else:
            print(f"  {p}: no 23_pool_members.csv - skipped")

    if not all_rows:
        print("No pools found. Run step 23 first.")
        sys.exit(0)

    perpool = pd.DataFrame(all_rows)
    summ = summarise(perpool)

    perpool_path = os.path.join(docs_dir, f"{DATE}_step23-conflict-perpool.csv")
    summ_path = os.path.join(docs_dir, f"{DATE}_step23-conflict-summary.csv")
    perpool.to_csv(perpool_path, index=False)
    summ.to_csv(summ_path, index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 30, "display.max_rows", 60):
        print("\n=== Per pathogen x category ===")
        print(summ.to_string(index=False))
        tot_conf = int(summ["total_conflicts"].sum())
        tot_act = int(summ["conflicts_active"].sum())
        tot_ina = int(summ["conflicts_inactive"].sum())
        print(f"\nTOTAL pools: {len(perpool)}   "
              f"multi-dataset: {int((perpool['n_datasets'] > 1).sum())}   "
              f"singletons: {int((perpool['n_datasets'] == 1).sum())}")
        print(f"TOTAL conflicts (before handling): {tot_conf}")
        if tot_conf:
            print(f"  -> resolved ACTIVE:   {tot_act} ({100*tot_act/tot_conf:.1f}%)")
            print(f"  -> resolved INACTIVE: {tot_ina} ({100*tot_ina/tot_conf:.1f}%)")
    print(f"\nwrote {perpool_path}\nwrote {summ_path}")
