import io
import random
import sys
import os
import zipfile

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import load_pathogen, load_expert_cutoffs
from dataset_utils import make_dataset_filename
from model_utils import load_ecfp_all, load_all_gz_csvs_from_zip, KFoldTrain, TrainRF

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print("Step 22: General organism-level modeling")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

decoy_ratio = 0.1

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------

print("Loading assay datasets metadata...")
cols = [
    "assay_id", "activity_type", "unit", "target_type_curated_extra",
    "dataset_type", "expert_cutoff",
]
assay_datasets = pd.read_csv(os.path.join(OUTPUT, "12_datasets.csv"))[cols]

# Keep only ORGANISM assays with quantitative or mixed datasets
assay_datasets = assay_datasets[
    (assay_datasets["target_type_curated_extra"] == "ORGANISM") &
    (assay_datasets["dataset_type"].isin(["quantitative", "mixed"]))
].copy()

print(f"  ORGANISM assay-cutoff rows: {len(assay_datasets)}")

print("Loading expert cutoffs...")
expert_cutoffs = load_expert_cutoffs(CONFIGPATH)

print("Loading ECFPs...")
ecfps = load_ecfp_all(os.path.join(DATAPATH, "chembl_processed", "06_chembl_ecfps.h5"))

print("Loading reference set and building decoy pool...")
reference_set = pd.read_csv(
    os.path.join(OUTPUT, "13_reference_set.csv.gz")
)["reference_compounds"].tolist()
x_ref = np.array([ecfps[cid] for cid in reference_set if cid in ecfps])
pathogen_compounds = set(
    pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"]
)
decoys_pool = set(i for i in ecfps if i not in pathogen_compounds)

print("Loading individual datasets from zip archives...")
dfs_qt = load_all_gz_csvs_from_zip(os.path.join(OUTPUT, "datasets", "datasets_qt.zip"))
dfs_mx = load_all_gz_csvs_from_zip(os.path.join(OUTPUT, "datasets", "datasets_mx.zip"))
print(f"  Quantitative: {len(dfs_qt)} | Mixed: {len(dfs_mx)}")

# ---------------------------------------------------------------------------
# Build one dataset per (activity_type, unit) pair at middle cutoff
# ---------------------------------------------------------------------------

activity_unit_pairs = (
    assay_datasets[["activity_type", "unit"]]
    .drop_duplicates()
    .reset_index(drop=True)
)
print(f"\nUnique (activity_type, unit) pairs in ORGANISM assays: {len(activity_unit_pairs)}")

results = []
out_datasets = {}

for _, pair in activity_unit_pairs.iterrows():
    activity_type = pair["activity_type"]
    unit = pair["unit"]
    unit_label = str(unit) if not (isinstance(unit, float) and np.isnan(unit)) else "None"

    # Get middle cutoff
    cutoff_key = (activity_type, unit, "ORGANISM", pathogen_code)
    cutoff_list = expert_cutoffs.get(cutoff_key)
    if not cutoff_list:
        print(f"  Skipping {activity_type} / {unit_label}: no expert cutoff found")
        continue
    mid_cutoff = cutoff_list[1] if len(cutoff_list) >= 2 else cutoff_list[0]

    # Select all assay rows for this pair at the middle cutoff
    unit_mask = (
        assay_datasets["unit"].eq(unit)
        if isinstance(unit, str)
        else assay_datasets["unit"].isna()
    )
    pair_rows = assay_datasets[
        (assay_datasets["activity_type"] == activity_type) &
        unit_mask &
        (assay_datasets["expert_cutoff"] == mid_cutoff)
    ]

    if len(pair_rows) == 0:
        print(f"  Skipping {activity_type} / {unit_label}: no assay rows at cutoff {mid_cutoff}")
        continue

    # Load datasets for all assays in this pair
    data_list = []
    for _, assay_row in pair_rows.iterrows():
        assay_id = assay_row["assay_id"]
        dataset_type = assay_row["dataset_type"]
        filename = make_dataset_filename(assay_id, activity_type, unit, dataset_type, mid_cutoff)
        dfs = dfs_qt if dataset_type == "quantitative" else dfs_mx
        if filename in dfs:
            data_list.append(dfs[filename][["compound_chembl_id", "smiles", "bin"]])

    if not data_list:
        print(f"  Skipping {activity_type} / {unit_label}: no dataset files found at cutoff {mid_cutoff}")
        continue

    data = pd.concat(data_list, ignore_index=True)

    # Deduplicate: active label wins; one row per compound
    data = (
        data.groupby("compound_chembl_id", as_index=False)
        .agg(bin=("bin", "max"), smiles=("smiles", "first"))
    )

    n_actives = int(data["bin"].sum())
    n_inactives = int((data["bin"] == 0).sum())
    n_compounds = len(data)
    n_assays = len(pair_rows)

    if n_actives == 0 or n_inactives == 0:
        print(f"  Skipping {activity_type} / {unit_label}: only one class present ({n_actives} actives, {n_inactives} inactives)")
        continue

    active_ratio = n_actives / n_compounds
    min_cpds = 1000 if active_ratio < 0.5 else 100
    if n_compounds < min_cpds:
        print(f"  Skipping {activity_type} / {unit_label}: too few compounds ({n_compounds} < {min_cpds})")
        continue
    if n_actives < 50:
        print(f"  Skipping {activity_type} / {unit_label}: too few actives ({n_actives})")
        continue

    print(f"  {activity_type} / {unit_label} | cutoff={mid_cutoff} | assays={n_assays} | cpds={n_compounds} | actives={n_actives} ({round(100*active_ratio,1)}%)")

    # Build X, Y for modeling (full data, no decoys yet)
    X = np.array(data["compound_chembl_id"].map(ecfps).tolist())
    Y = np.array(data["bin"].tolist())

    # Add decoys for modeling if active ratio > 0.5
    if active_ratio > 0.5:
        n_decoys = int(n_actives / decoy_ratio - n_inactives)
        print(f"    Adding {n_decoys} decoys for modeling")
        rng = random.Random(42)
        decoy_ids = rng.sample(list(decoys_pool), n_decoys)
        X_decoys = np.array([ecfps[i] for i in decoy_ids])
        X = np.vstack([X, X_decoys])
        Y = np.concatenate([Y, np.zeros(len(X_decoys), dtype=Y.dtype)])

    # Downsample negatives for modeling if active ratio < 5%
    active_ratio_model = n_actives / len(Y)
    if active_ratio_model < 0.05:
        n_neg_target = int(n_actives / 0.10) - n_actives
        neg_idx = np.where(Y == 0)[0]
        rng_ds = np.random.RandomState(42)
        sampled_neg = rng_ds.choice(neg_idx, size=min(n_neg_target, len(neg_idx)), replace=False)
        keep = np.sort(np.concatenate([np.where(Y == 1)[0], sampled_neg]))
        X, Y = X[keep], Y[keep]
        print(f"    Downsampled negatives for modeling: ratio {active_ratio_model:.3f} → {n_actives/len(Y):.3f}")

    avg_auroc, std_auroc = KFoldTrain(X, Y)
    print(f"    AUROC: {avg_auroc} ± {std_auroc}")

    results.append({
        "activity_type": activity_type,
        "unit": unit,
        "cutoff": mid_cutoff,
        "n_assays": n_assays,
        "n_compounds": n_compounds,
        "n_actives": n_actives,
        "n_inactives": n_inactives,
        "auroc": avg_auroc,
        "auroc_std": std_auroc,
    })

    unit_sanitized = str(unit).replace("/", "FwdS")
    dataset_name = f"ORG_{activity_type}_{unit_sanitized}_{mid_cutoff}"
    out_datasets[dataset_name] = data

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(OUTPUT, "22_general_model.csv"), index=False)
print(f"\nSaved {len(results_df)} modeled datasets to 22_general_model.csv")

zip_path = os.path.join(OUTPUT, "datasets", "22_general_datasets.zip")
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
    for name, df in out_datasets.items():
        filename = f"{name}.csv.gz"
        buf = io.BytesIO()
        df.to_csv(buf, index=False, compression="gzip")
        zf.writestr(filename, buf.getvalue())
print(f"Saved {len(out_datasets)} datasets to {zip_path}")

if len(results_df) > 0:
    print(f"\nAUROC summary: mean={results_df['auroc'].mean():.3f}, "
          f"min={results_df['auroc'].min():.3f}, max={results_df['auroc'].max():.3f}")
    print(f"Total compounds across all datasets: {results_df['n_compounds'].sum()}")
