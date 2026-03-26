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

# Load expert cutoffs for hierarchical assignment
from pathogen_utils import load_expert_cutoffs
expert_cutoffs = load_expert_cutoffs(os.path.join(root, "..", "config"))

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
# Hierarchical assignment: determine which assays qualify for A vs B
# ---------------------------------------------------------------------------

def determine_hierarchical_assignment(assay_datasets, expert_cutoffs, pathogen_code):
    """Determine which assays qualify for A vs B based on middle cutoff evaluation."""
    qualified_for_a = set()
    qualified_for_b = set()
    qualified_for_neither = set()

    print("Hierarchical assignment analysis:")
    total_assays = 0
    middle_cutoff_available = 0
    fallback_used = 0

    for assay_key, group in assay_datasets.groupby(["assay_id", "activity_type", "unit"]):
        total_assays += 1
        target_type = group["target_type_curated_extra"].iloc[0]
        cutoff_key = (assay_key[1], assay_key[2], target_type, pathogen_code)
        cutoff_list = expert_cutoffs.get(cutoff_key)

        if cutoff_list and len(cutoff_list) >= 2:
            middle_cutoff_available += 1
            middle_cutoff = cutoff_list[1]  # Index 1 is middle cutoff
            middle_row = group[group["expert_cutoff"] == middle_cutoff]

            if not middle_row.empty:
                middle_data = pd.DataFrame([middle_row.iloc[0]])
                if condition_a(middle_data).iloc[0]:
                    qualified_for_a.add(assay_key)
                elif condition_b(middle_data).iloc[0]:  # ✅ Explicit check
                    qualified_for_b.add(assay_key)
                else:
                    qualified_for_neither.add(assay_key)  # ✅ Track neither case
        else:
            # Fallback for edge cases: use first available cutoff
            fallback_used += 1
            if not group.empty:
                fallback_data = pd.DataFrame([group.iloc[0]])
                if condition_a(fallback_data).iloc[0]:
                    qualified_for_a.add(assay_key)
                elif condition_b(fallback_data).iloc[0]:  # ✅ Explicit check
                    qualified_for_b.add(assay_key)
                else:
                    qualified_for_neither.add(assay_key)  # ✅ Track neither case

    print(f"  Total assays evaluated: {total_assays}")
    print(f"  Middle cutoff available: {middle_cutoff_available}")
    print(f"  Fallback used: {fallback_used}")
    print(f"  Qualified for A: {len(qualified_for_a)}")
    print(f"  Qualified for B: {len(qualified_for_b)}")
    print(f"  Qualified for neither: {len(qualified_for_neither)}")

    return qualified_for_a, qualified_for_b

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

# Show overall condition compliance before hierarchical assignment
print("Overall condition compliance (all datasets):")
all_condition_a = condition_a(assay_datasets)
all_condition_b = condition_b(assay_datasets)
print(f"  Total datasets: {len(assay_datasets)}")
print(f"  Meet condition A: {all_condition_a.sum()}")
print(f"  Meet condition B: {all_condition_b.sum()}")
print(f"  Meet both A and B: {(all_condition_a & all_condition_b).sum()}")
print(f"  Meet neither A nor B: {(~all_condition_a & ~all_condition_b).sum()}")

print("\n" + "="*60)
# Determine hierarchical assignment before modeling
qualified_for_a, qualified_for_b = determine_hierarchical_assignment(
    assay_datasets, expert_cutoffs, pathogen_code
)
print("="*60)

results = []
label_compounds = {label: set() for label in conditions}

for label, condition_fn in conditions.items():

    print(f"\nCondition {label} modeling...")

    if label == "A":
        # Only model assays pre-qualified for A
        mask = assay_datasets.apply(
            lambda x: (x["assay_id"], x["activity_type"], x["unit"]) in qualified_for_a,
            axis=1
        )
        pre_qualified = assay_datasets[mask]
        print(f"  Pre-qualified datasets: {len(pre_qualified)}")

        # Apply condition A requirements to pre-qualified datasets
        condition_datasets = pre_qualified[condition_fn(pre_qualified)]
        print(f"  After condition A filter: {len(condition_datasets)} (filtered out: {len(pre_qualified) - len(condition_datasets)})")

    elif label == "B":
        # Only model assays pre-qualified for B
        mask = assay_datasets.apply(
            lambda x: (x["assay_id"], x["activity_type"], x["unit"]) in qualified_for_b,
            axis=1
        )
        pre_qualified = assay_datasets[mask]
        print(f"  Pre-qualified datasets: {len(pre_qualified)}")

        # Apply condition B requirements to pre-qualified datasets
        condition_datasets = pre_qualified[condition_fn(pre_qualified)]
        print(f"  After condition B filter: {len(condition_datasets)} (filtered out: {len(pre_qualified) - len(condition_datasets)})")

    condition_datasets = condition_datasets.copy().reset_index(drop=True)
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
