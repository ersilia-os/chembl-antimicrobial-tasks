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

print("Step 16: Selecting merged datasets")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

# Load expert cutoffs
expert_cutoffs = load_expert_cutoffs(CONFIGPATH)

# Load merged LM results
merged_lm = pd.read_csv(os.path.join(OUTPUT, "15_merged_LM.csv"))

pathogen_compounds = set(pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"])

# ---------------------------------------------------------------------------
# Identify unique merged groups (base name + metadata)
# Each row in merged_lm is one (group, cutoff) combination.
# Base name is everything before the last underscore-separated cutoff suffix.
# ---------------------------------------------------------------------------

# Extract base name: name_ is "{base}_{cutoff}", e.g. "M_ORG0_1.0"
# First two underscore-parts are always the base ("M_ORG0" or "M_SP0").
merged_lm["base_name"] = merged_lm["name"].apply(lambda n: "_".join(n.split("_")[:2]))

groups = sorted(
    set(
        (base, act, unit, tt)
        for base, act, unit, tt in merged_lm[["base_name", "activity_type", "unit", "target_type_curated_extra"]].values
    )
)

cols_to_keep = [
    "name", "direction", "assay_type", "strain", "target_chembl_id",
    "n_assays", "n_cpds_union", "positives", "ratio", "avg", "std", "assay_keys",
]

selected = []

for base_name, activity_type, unit, target_type in groups:

    group_rows = merged_lm[merged_lm["base_name"] == base_name].reset_index(drop=True)

    # Mid cutoff: second value in the expert cutoffs list for this key
    cutoff_list = expert_cutoffs.get((activity_type, unit, target_type, pathogen_code), [])
    mid_cutoff = cutoff_list[1] if len(cutoff_list) >= 2 else np.nan

    df = group_rows.sort_values("avg", ascending=False).reset_index(drop=True)

    best_auroc = df["avg"].iloc[0]
    best_cutoff = df["expert_cutoff"].iloc[0]

    mid_rows = df[df["expert_cutoff"] == mid_cutoff]
    mid_auroc = mid_rows["avg"].iloc[0] if len(mid_rows) > 0 else np.nan

    if best_auroc <= 0.7:
        continue

    if np.isnan(mid_auroc) or (best_auroc - mid_auroc) > 0.1:
        info = df[cols_to_keep].iloc[0].tolist()
        selected.append([activity_type, unit, target_type, best_cutoff, best_auroc, False] + info)
    else:
        info = mid_rows[cols_to_keep].iloc[0].tolist()
        selected.append([activity_type, unit, target_type, mid_cutoff, mid_auroc, True] + info)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

selected_df = pd.DataFrame(selected, columns=[
    "activity_type", "unit", "target_type", "cutoff", "auroc", "is_mid_cutoff",
    "name", "direction", "assay_type", "strain", "target_chembl_id",
    "n_assays", "n_cpds_union", "positives", "ratio", "avg", "std", "assay_keys",
])
selected_df.to_csv(os.path.join(OUTPUT, "16_merged_selected_LM.csv"), index=False)

# Coverage from actual saved merged datasets (avoids assay-key lookup issues)
merged_dir = os.path.join(OUTPUT, "datasets", "M")
original_compounds = {"ORG": set(), "SP": set()}
selected_compounds = {"ORG": set(), "SP": set()}

if os.path.exists(merged_dir):
    for filename in os.listdir(merged_dir):
        if not filename.endswith(".csv.gz"):
            continue
        base = "_".join(filename.replace(".csv.gz", "").split("_")[:2])
        ty = "ORG" if "ORG" in base else "SP"
        df = pd.read_csv(os.path.join(merged_dir, filename))
        original_compounds[ty].update(df[df["smiles"] != "decoy"]["compound_chembl_id"].tolist())

for _, row in selected_df.iterrows():
    ty = "ORG" if "ORG" in row["name"] else "SP"
    filepath = os.path.join(merged_dir, f"{row['name']}.csv.gz")
    if os.path.exists(filepath):
        df = pd.read_csv(filepath)
        selected_compounds[ty].update(df[df["smiles"] != "decoy"]["compound_chembl_id"].tolist())

print("Chemical space coverage:")
for ty in ["ORG", "SP"]:
    before = round(100 * len(original_compounds[ty]) / len(pathogen_compounds), 1)
    after = round(100 * len(selected_compounds[ty]) / len(pathogen_compounds), 1)
    print(f"  {ty}: {before}% → {after}%")
overall_before = round(100 * len(original_compounds["ORG"] | original_compounds["SP"]) / len(pathogen_compounds), 1)
overall_after = round(100 * len(selected_compounds["ORG"] | selected_compounds["SP"]) / len(pathogen_compounds), 1)
print(f"  Overall: {overall_before}% → {overall_after}%")
print(f"Selected datasets — ORG: {len(selected_df[selected_df['target_type'] == 'ORGANISM'])}, SP: {len(selected_df[selected_df['target_type'] == 'SINGLE PROTEIN'])}")
print(f"Datasets using mid cutoff: {selected_df['is_mid_cutoff'].sum()}/{len(selected_df)}")
