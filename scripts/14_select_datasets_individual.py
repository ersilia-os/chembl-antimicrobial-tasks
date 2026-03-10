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
assay_to_compounds = defaultdict(set)
for assay_id, activity_type, unit, compound_chembl_id in chembl_pathogen[["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"]].values:
    assay_to_compounds[(assay_id, activity_type, unit)].add(compound_chembl_id)
del chembl_pathogen

# Load expert cutoffs and individual LM results
expert_cutoffs = load_expert_cutoffs(CONFIGPATH)
individual_lm = pd.read_csv(os.path.join(OUTPUT, "13_individual_LM.csv"))

# All pathogen compounds (for coverage stats)
pathogen_compounds = set(pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"])

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

        key = (assay_id, activity_type, unit)
        original_compounds[label].update(assay_to_compounds[key])

        # Each assay selected at most once across cutoffs
        if key in already_selected:
            continue

        # Mid cutoff is the reference cutoff (middle value in the expert list)
        mid_cutoff = expert_cutoffs[(activity_type, unit, target_type, pathogen_code)][1]

        # Get all cutoff results for this assay, sorted by AUROC
        df = lm_label[
            (lm_label["assay_id"] == assay_id) &
            (lm_label["activity_type"] == activity_type) &
            (lm_label["unit"].eq(unit) if isinstance(unit, str) else lm_label["unit"].isna())
        ].sort_values("avg", ascending=False).reset_index(drop=True)

        best_auroc = df["avg"].iloc[0]
        best_cutoff = df["expert_cutoff"].iloc[0]

        mid_rows = df[df["expert_cutoff"] == mid_cutoff]
        mid_auroc = mid_rows["avg"].iloc[0] if len(mid_rows) > 0 else np.nan

        if best_auroc > 0.7:

            # Prefer the mid cutoff unless the best is substantially better (> 0.1 AUROC)
            if np.isnan(mid_auroc) or (best_auroc - mid_auroc) > 0.1:
                info = df[cols_to_keep].iloc[0].tolist()
                selected.append([label, assay_id, activity_type, unit, target_type, best_cutoff, best_auroc, False] + info)
            else:
                info = mid_rows[cols_to_keep].iloc[0].tolist()
                selected.append([label, assay_id, activity_type, unit, target_type, mid_cutoff, mid_auroc, True] + info)

            selected_compounds[label].update(assay_to_compounds[key])
            already_selected.add(key)

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
