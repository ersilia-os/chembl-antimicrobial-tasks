import io
import sys
import os
import zipfile

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH, DECOY_RATIO
from pathogen_utils import load_pathogen, load_expert_cutoffs
from dataset_utils import make_dataset_filename
from model_utils import load_ecfp_all, load_all_gz_csvs_from_zip, KFoldTrain, TrainRF, add_decoys, downsample_negatives

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print("Step 20: General organism-level modeling")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)


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

print("Loading decoy pool...")
pathogen_compounds = set(
    pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"]
)
decoys_pool = set(i for i in ecfps if i not in pathogen_compounds)

print("Loading individual datasets from zip archives...")
dfs_qt = load_all_gz_csvs_from_zip(os.path.join(OUTPUT, "12_datasets", "datasets_qt.zip"))
dfs_mx = load_all_gz_csvs_from_zip(os.path.join(OUTPUT, "12_datasets", "datasets_mx.zip"))
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

    # Deduplicate: active label wins; one row per compound, then per SMILES
    data = (
        data.groupby("compound_chembl_id", as_index=False)
        .agg(bin=("bin", "max"), smiles=("smiles", "first"))
    )
    data = (data
        .sort_values("bin", ascending=False)
        .drop_duplicates(subset="smiles", keep="first")
        .reset_index(drop=True))

    n_actives = int(data["bin"].sum())
    n_inactives = int((data["bin"] == 0).sum())
    n_compounds = len(data)
    n_assays = len(pair_rows)

    if n_actives == 0:
        print(f"  Skipping {activity_type} / {unit_label}: no actives found")
        continue

    active_ratio = n_actives / n_compounds
    print(f"  {activity_type} / {unit_label} | cutoff={mid_cutoff} | assays={n_assays} | cpds={n_compounds} | actives={n_actives} ({round(100*active_ratio,1)}%)")

    # Build X, Y for modeling (full data, no decoys yet)
    X = np.array(data["compound_chembl_id"].map(ecfps).tolist())
    Y = np.array(data["bin"].tolist())

    # Add decoys for modeling if no inactives or active ratio > 0.5
    if n_inactives == 0 or active_ratio > 0.5:
        X, Y, decoy_ids = add_decoys(X, Y, decoys_pool, ecfps, DECOY_RATIO)
        print(f"    Added {len(decoy_ids)} decoys for modeling")

    # Downsample negatives for modeling if active ratio < 5%
    active_ratio_model = n_actives / len(Y)
    X, Y = downsample_negatives(X, Y)
    if n_actives / len(Y) != active_ratio_model:
        print(f"    Downsampled negatives for modeling: ratio {active_ratio_model:.3f} → {n_actives/len(Y):.3f}")

    if n_actives > 10:
        avg_auroc, std_auroc = KFoldTrain(X, Y)
        print(f"    AUROC: {avg_auroc} ± {std_auroc}")
    else:
        avg_auroc, std_auroc = np.nan, np.nan
        print(f"    Skipping AUROC: too few actives ({n_actives}) for 4-fold CV")

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
# Build 20_all_positives.csv
# ---------------------------------------------------------------------------

all_actives = []
all_smiles_list = []
for _, row in assay_datasets.iterrows():
    filename = make_dataset_filename(row["assay_id"], row["activity_type"], row["unit"], row["dataset_type"], row["expert_cutoff"])
    dfs = dfs_qt if row["dataset_type"] == "quantitative" else dfs_mx
    if filename in dfs:
        df = dfs[filename]
        all_actives.append(df[df["bin"] == 1][["smiles"]])
        all_smiles_list.append(df[["smiles"]])

if all_actives:
    all_actives_df = pd.concat(all_actives, ignore_index=True)
    all_positives = (
        all_actives_df["smiles"]
        .value_counts()
        .rename_axis("smiles")
        .reset_index(name="occurrences")
    )
    all_positives.to_csv(os.path.join(OUTPUT, "20_all_positives.csv"), index=False)
    print(f"Saved {len(all_positives)} unique active SMILES to 20_all_positives.csv")
else:
    print("No actives found across ORGANISM assays; 20_all_positives.csv not written")

if all_smiles_list:
    all_smiles = (
        pd.concat(all_smiles_list, ignore_index=True)["smiles"]
        .dropna()
        .drop_duplicates()
        .reset_index(drop=True)
        .to_frame()
    )
    all_smiles.to_csv(os.path.join(OUTPUT, "20_all_smiles.csv"), index=False)
    print(f"Saved {len(all_smiles)} unique SMILES to 20_all_smiles.csv")
else:
    print("No SMILES found across ORGANISM assays; 20_all_smiles.csv not written")

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(OUTPUT, "20_general_datasets.csv"), index=False)
print(f"\nSaved {len(results_df)} modeled datasets to 20_general_datasets.csv")

zip_path = os.path.join(OUTPUT, "20_general_datasets.zip")
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
