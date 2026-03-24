from collections import defaultdict
import pandas as pd
import numpy as np
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH
from pathogen_utils import load_pathogen, load_expert_cutoffs

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print("Step 14: Selecting individual datasets")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

# Load ChEMBL data for pathogen
chembl_pathogen = pd.read_csv(os.path.join(OUTPUT, "08_chembl_cleaned_data.csv.gz"), low_memory=False)

# Build assay → compound set mapping, then free the full table
# Use None instead of NaN for consistent dictionary keys
assay_to_compounds = defaultdict(set)
for assay_id, activity_type, unit, compound_chembl_id in chembl_pathogen[["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"]].values:
    # Normalize NaN to None for consistent dictionary keys
    unit_key = None if pd.isna(unit) else unit
    assay_to_compounds[(assay_id, activity_type, unit_key)].add(compound_chembl_id)
del chembl_pathogen

# Load expert cutoffs and individual LM results
expert_cutoffs = load_expert_cutoffs(CONFIGPATH)
individual_lm = pd.read_csv(os.path.join(OUTPUT, "13_individual_LM.csv"))

# All pathogen compounds (for coverage stats)
pathogen_compounds = set(pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"])

# AUROC improvement threshold: use best cutoff if improvement > 10% AUROC
AUROC_IMPROVEMENT_THRESHOLD = 0.1

labels = ["A", "B"]
cols_to_keep = ["dataset_type", "pos_qt", "ratio_qt", "cpds_qt", "pos_ql", "ratio_ql", "cpds_ql", "overlap_mx", "pos_mx", "ratio_mx", "cpds_mx"]
keys = ["assay_id", "activity_type", "unit", "target_type_curated_extra"]

selected = []
original_compounds = {label: set() for label in labels}
selected_compounds = {label: set() for label in labels}

for label in labels:

    lm_label = individual_lm[individual_lm["label"] == label]
    assays_label = set(tuple(row) for row in lm_label[keys].values)

    already_selected = set()

    for assay_id, activity_type, unit, target_type in assays_label:

        # Normalize NaN to None for consistent dictionary lookup
        unit_key = None if pd.isna(unit) else unit
        key = (assay_id, activity_type, unit_key)
        original_compounds[label].update(assay_to_compounds[key])

        # Each assay selected at most once across cutoffs
        # Use original (assay_id, activity_type, unit) triplet for tracking selection
        selection_key = (assay_id, activity_type, unit)
        if selection_key in already_selected:
            continue

        # Mid cutoff is the reference cutoff (middle value in the expert list)
        cutoff_list = expert_cutoffs.get((activity_type, unit, target_type, pathogen_code))
        if not cutoff_list or len(cutoff_list) < 2:
            print(f"Warning: Missing cutoff for {activity_type}, {unit}, {target_type}")
            continue
        mid_cutoff = cutoff_list[1]

        # Get all cutoff results for this assay, sorted by AUROC
        if pd.isna(unit):
            unit_mask = lm_label["unit"].isna()
        else:
            unit_mask = lm_label["unit"].eq(unit)

        df = lm_label[
            (lm_label["assay_id"] == assay_id) &
            (lm_label["activity_type"] == activity_type) &
            unit_mask
        ].sort_values("avg", ascending=False).reset_index(drop=True)

        if len(df) == 0:
            print(f"Warning: No LM results for {assay_id}, {activity_type}, {unit}")
            continue

        best_auroc = df["avg"].iloc[0]
        best_cutoff = df["expert_cutoff"].iloc[0]

        mid_rows = df[df["expert_cutoff"] == mid_cutoff]
        mid_auroc = mid_rows["avg"].iloc[0] if len(mid_rows) > 0 else np.nan

        if len(mid_rows) == 0:
            print(f"Info: Mid cutoff {mid_cutoff} not found for {assay_id}, {activity_type}, {unit}, using best: {best_cutoff}")

        if best_auroc > 0.7:

            # Prefer the mid cutoff unless the best is substantially better (> threshold AUROC)
            if np.isnan(mid_auroc) or (best_auroc - mid_auroc) > AUROC_IMPROVEMENT_THRESHOLD:
                info = df[cols_to_keep].iloc[0].tolist()
                selected.append([label, assay_id, activity_type, unit, target_type, best_cutoff, best_auroc, False] + info)
            else:
                info = mid_rows[cols_to_keep].iloc[0].tolist()
                selected.append([label, assay_id, activity_type, unit, target_type, mid_cutoff, mid_auroc, True] + info)

            selected_compounds[label].update(assay_to_compounds[key])
            already_selected.add(selection_key)

selected_df = pd.DataFrame(selected, columns=[
    "label", "assay_id", "activity_type", "unit", "target_type", "cutoff", "auroc", "is_mid_cutoff",
    "dataset_type", "pos_qt", "ratio_qt", "cpds_qt", "pos_ql", "ratio_ql", "cpds_ql",
    "overlap_mx", "pos_mx", "ratio_mx", "cpds_mx",
])

selected_df.to_csv(os.path.join(OUTPUT, "14_individual_selected_LM.csv"), index=False)

# Sanity check: one dataset per assay triplet within each label
for label in labels:
    sub = selected_df[selected_df["label"] == label]
    assert len(set(tuple(row) for row in sub[["assay_id", "activity_type", "unit"]].values)) == len(sub), \
        f"Duplicate assay triplets in label {label}"

print("Chemical space coverage:")
for label in labels:
    before = round(100 * len(original_compounds[label]) / len(pathogen_compounds), 1)
    after = round(100 * len(selected_compounds[label]) / len(pathogen_compounds), 1)
    print(f"  {label}: {before}% → {after}%")
overall_before = round(100 * len(original_compounds["A"] | original_compounds["B"]) / len(pathogen_compounds), 1)
overall_after = round(100 * len(selected_compounds["A"] | selected_compounds["B"]) / len(pathogen_compounds), 1)
print(f"  Overall: {overall_before}% → {overall_after}%")
print(f"Selected datasets — A: {len(selected_df[selected_df['label'] == 'A'])}, B: {len(selected_df[selected_df['label'] == 'B'])}")
print(f"Datasets using mid cutoff: {selected_df['is_mid_cutoff'].sum()}/{len(selected_df)}")
