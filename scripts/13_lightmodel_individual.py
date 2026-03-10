from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import random
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH
from pathogen_utils import load_pathogen
from dataset_utils import make_dataset_filename
from model_utils import (
    load_ecfp_all, load_data_from_zip, KFoldTrain, TrainRF,
)

# ---------------------------------------------------------------------------
# Modeling conditions
# ---------------------------------------------------------------------------

# Condition A: large, class-balanced datasets — no decoys needed
def condition_a(df):
    return (
        df["dataset_type"].isin(["quantitative", "mixed"])
        & (df["cpds_qt"] >= 1000)
        & (df["pos_qt"] >= 50)
        & df["ratio_qt"].between(0.001, 0.5, inclusive="both")
    )

# Condition B: active-enriched datasets — random ChEMBL decoys added to reach target ratio
def condition_b(df):
    return (
        df["dataset_type"].isin(["quantitative", "mixed"])
        & (df["pos_qt"] >= 100)
        & (df["ratio_qt"] >= 0.5)
    )

# Target active ratio when adding decoys for condition B
decoy_ratio = 0.1

conditions = {"A": condition_a, "B": condition_b}

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print("Step 13")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

# Load dataset summary from step 12
cols = [
    "assay_id", "activity_type", "unit", "target_type", "target_type_curated_extra",
    "cpds", "direction", "dataset_type", "expert_cutoff",
    "pos_qt", "ratio_qt", "cpds_qt", "pos_ql", "ratio_ql", "cpds_ql",
    "overlap_mx", "pos_mx", "ratio_mx", "cpds_mx",
]
assay_datasets = pd.read_csv(os.path.join(OUTPUT, "12_datasets.csv"))[cols]

pathogen_compounds = set(pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"])

# Load all Morgan fingerprints into memory
print("Loading ECFPs...")
ecfps = load_ecfp_all(os.path.join(DATAPATH, "chembl_processed", "06_chembl_ecfps.h5"))

# Decoy pool: all ChEMBL compounds not tested against this pathogen
decoys_pool = set(i for i in ecfps if i not in pathogen_compounds)

# Reference set: 10,000 randomly sampled decoys (compounds never screened against the pathogen).
# Using unseen compounds avoids data leakage — correlations between model predictions reflect
# shared biological signal rather than shared training data.
print(f"Sampling reference set from decoy pool ({len(decoys_pool)} compounds)...")
rng_ref = random.Random(42)
reference_set = rng_ref.sample(sorted(decoys_pool), min(10000, len(decoys_pool)))
pd.DataFrame(reference_set, columns=["reference_compounds"]).to_csv(
    os.path.join(OUTPUT, "13_reference_set.csv.gz"), index=False
)

# Reference feature matrix
x_ref = np.array([ecfps[cid] for cid in reference_set])

# Output directory for reference set predictions
PATH_TO_CORRELATIONS = os.path.join(OUTPUT, "correlations")
os.makedirs(PATH_TO_CORRELATIONS, exist_ok=True)

# ---------------------------------------------------------------------------
# Helper: model a single assay
# ---------------------------------------------------------------------------

def run_model(assay, label):
    """Load data, optionally add decoys, run cross-validation, train final model.

    Returns
    -------
    avg_auroc : float
    std_auroc : float
    rf : fitted RandomForestClassifier
    filename : str  (used to name the saved reference predictions)
    compound_ids : set  (compounds in this dataset, for coverage stats)
    """
    assay_id = assay.assay_id
    activity_type = assay.activity_type
    unit = assay.unit
    expert_cutoff = assay.expert_cutoff
    dataset_type = assay.dataset_type

    filename = make_dataset_filename(assay_id, activity_type, unit, dataset_type, expert_cutoff)
    zip_name = "datasets_qt.zip" if dataset_type == "quantitative" else "datasets_mx.zip"
    df = load_data_from_zip(os.path.join(OUTPUT, "datasets", zip_name), filename)

    compound_ids = set(df["compound_chembl_id"].astype(str))
    X = np.array(df["compound_chembl_id"].map(ecfps).tolist())
    Y = np.array(df["bin"].tolist())
    n_positives = int(Y.sum())

    print(f"  {assay_id} | {activity_type} | {unit} | cutoff={expert_cutoff}")
    print(f"    Compounds: {len(X)}, Positives: {n_positives} ({round(100 * n_positives / len(Y), 1)}%)")

    if label == "B":
        n_decoys = max(0, int(n_positives / decoy_ratio) - len(Y))
        print(f"    Adding {n_decoys} decoys from ChEMBL")
        rng = random.Random(42)
        decoy_ids = rng.sample(list(decoys_pool), n_decoys)
        X_decoys = np.array([ecfps[i] for i in decoy_ids])
        X = np.vstack([X, X_decoys])
        Y = np.concatenate([Y, np.zeros(len(X_decoys), dtype=Y.dtype)])
        print(f"    After decoys: {len(X)} compounds, {n_positives} positives ({round(100 * n_positives / len(Y), 1)}%)")

    avg_auroc, std_auroc = KFoldTrain(X, Y)
    print(f"    AUROC: {avg_auroc} ± {std_auroc}")

    rf = TrainRF(X, Y)
    return avg_auroc, std_auroc, rf, filename, compound_ids

# ---------------------------------------------------------------------------
# Main loop — one pass per condition
# ---------------------------------------------------------------------------

results = []
label_compounds = {label: set() for label in conditions}

for label, condition_fn in conditions.items():

    print(f"\nCondition {label}...")
    condition_datasets = assay_datasets[condition_fn(assay_datasets)].copy().reset_index(drop=True)
    condition_datasets["label"] = label

    os.makedirs(os.path.join(PATH_TO_CORRELATIONS, label), exist_ok=True)

    auroc_avg, auroc_std = [], []

    for _, assay in tqdm(condition_datasets.iterrows(), total=len(condition_datasets)):

        avg_auroc, std_auroc, rf, filename, compound_ids = run_model(assay, label)

        auroc_avg.append(avg_auroc)
        auroc_std.append(std_auroc)
        label_compounds[label].update(compound_ids)

        # Save reference set predictions
        y_prob_ref = rf.predict_proba(x_ref)[:, 1]
        ref_filename = filename.replace(".csv.gz", "_ref_probs.npz")
        np.savez_compressed(os.path.join(PATH_TO_CORRELATIONS, label, ref_filename), y_prob_ref=y_prob_ref)

    condition_datasets["avg"] = auroc_avg
    condition_datasets["std"] = auroc_std
    results.append(condition_datasets)

    n_datasets = len(condition_datasets)
    n_assays = len(set(map(tuple, condition_datasets[["assay_id", "activity_type", "unit"]].values)))
    coverage = round(100 * len(label_compounds[label]) / len(pathogen_compounds), 1)
    print(f"  Datasets: {n_datasets} | Assays: {n_assays} | Coverage: {coverage}%")

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

individual_lm = pd.concat(results, ignore_index=True)
individual_lm.to_csv(os.path.join(OUTPUT, "13_individual_LM.csv"), index=False)

all_compounds = set(cid for lab in label_compounds.values() for cid in lab)
print(f"\nOverall coverage (A+B): {round(100 * len(all_compounds) / len(pathogen_compounds), 1)}%")
print(f"Total datasets modeled: {len(individual_lm)}  {dict(Counter(individual_lm['label']))}")
