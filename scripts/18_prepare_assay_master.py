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

print("Step 18: Preparing assay master table")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

keys = ["assay_id", "activity_type", "unit"]

expert_cutoffs = load_expert_cutoffs(CONFIGPATH)

# ---------------------------------------------------------------------------
# Load and merge assay metadata tables
# ---------------------------------------------------------------------------

columns_clusters = ["clusters_0.3", "clusters_0.6", "clusters_0.85"]
columns_parameters = [
    "organism_curated", "target_type_curated", "target_name_curated",
    "target_chembl_id_curated", "strain", "atcc_id", "mutations",
    "known_drug_resistances", "media",
]
columns_data_info = [
    "target_type_curated_extra", "dataset_type",
    "equal", "higher", "lower",
    "min_", "p1", "p25", "p50", "p75", "p99", "max_",
]

print("Loading data...")
assays_cleaned = pd.read_csv(os.path.join(OUTPUT, "08_assays_cleaned.csv"))
assays_clusters = pd.read_csv(os.path.join(OUTPUT, "10_assays_clusters.csv"))[keys + columns_clusters]
assays_parameters = pd.read_csv(os.path.join(OUTPUT, "09_assays_parameters.csv"))[keys + columns_parameters]
assay_data_info = pd.read_csv(os.path.join(OUTPUT, "12_assay_data_info.csv"))[keys + columns_data_info]

print("Merging tables...")
assays_master = assays_cleaned.merge(assays_clusters, on=keys, how="left", validate="1:1")
assays_master = assays_master.merge(assays_parameters, on=keys, how="left", validate="1:1")
assays_master = assays_master.merge(assay_data_info, on=keys, how="left", validate="1:1")

# Add evaluated cutoffs column
assays_master["evaluated_cutoffs"] = [
    ";".join(str(c) for c in expert_cutoffs[(a, u, tt, pathogen_code)])
    if (a, u, tt, pathogen_code) in expert_cutoffs else np.nan
    for a, u, tt in assays_master[["activity_type", "unit", "target_type_curated_extra"]].values
]

# ---------------------------------------------------------------------------
# Load modeling results to annotate each assay's pipeline status
# ---------------------------------------------------------------------------

def assay_keys_to_set(df, label=None):
    """Extract all (assay_id, activity_type, unit) tuples from an assay_keys column."""
    subset = df if label is None else df[df["label"] == label]
    return set(tuple(s.split("|")) for row in subset["assay_keys"] for s in row.split(";"))

individual_lm = pd.read_csv(os.path.join(OUTPUT, "13_individual_LM.csv"))
merged_lm = pd.read_csv(os.path.join(OUTPUT, "15_merged_LM.csv"))
individual_selected_lm = pd.read_csv(os.path.join(OUTPUT, "14_individual_selected_LM.csv"))
merged_selected_lm = pd.read_csv(os.path.join(OUTPUT, "16_merged_selected_LM.csv"))
final_datasets = pd.read_csv(os.path.join(OUTPUT, "17_final_datasets.csv"))

considered_a = set(tuple(r) for r in individual_lm[individual_lm["label"] == "A"][keys].values)
considered_b = set(tuple(r) for r in individual_lm[individual_lm["label"] == "B"][keys].values)
considered_m = assay_keys_to_set(merged_lm)

selected_a = set(tuple(r) for r in individual_selected_lm[individual_selected_lm["label"] == "A"][keys].values)
selected_b = set(tuple(r) for r in individual_selected_lm[individual_selected_lm["label"] == "B"][keys].values)
selected_m = assay_keys_to_set(merged_selected_lm)

final_considered_a = assay_keys_to_set(final_datasets, label="A")
final_considered_b = assay_keys_to_set(final_datasets, label="B")
final_considered_m = assay_keys_to_set(final_datasets, label="M")

final_only = final_datasets[final_datasets["selected"]].reset_index(drop=True)
final_selected_a = assay_keys_to_set(final_only, label="A")
final_selected_b = assay_keys_to_set(final_only, label="B")
final_selected_m = assay_keys_to_set(final_only, label="M")

# ---------------------------------------------------------------------------
# Annotate each assay with its pipeline status per label
# ---------------------------------------------------------------------------

def pipeline_comment(key, considered, selected, final_considered, final_selected):
    if key not in considered:
        return "Not considered"
    if key not in selected:
        return "Considered but not selected"
    if key not in final_considered:
        return "Considered and selected, but not in correlation analysis (non-ORGANISM)"
    if key not in final_selected:
        return "Considered and selected, but discarded due to high correlation"
    return "Considered, selected, and retained in final selection"

comments_a, comments_b, comments_m = [], [], []
for assay_id, activity_type, unit in assays_master[keys].values:
    key = (assay_id, activity_type, unit)
    comments_a.append(pipeline_comment(key, considered_a, selected_a, final_considered_a, final_selected_a))
    comments_b.append(pipeline_comment(key, considered_b, selected_b, final_considered_b, final_selected_b))
    comments_m.append(pipeline_comment(key, considered_m, selected_m, final_considered_m, final_selected_m))

assays_master["comment_A"] = comments_a
assays_master["comment_B"] = comments_b
assays_master["comment_M"] = comments_m

# ---------------------------------------------------------------------------
# Select and save final columns
# ---------------------------------------------------------------------------

all_cols = [
    "assay_id", "assay_type", "assay_organism", "target_organism",
    "organism_curated", "doc_chembl_id",
    "target_type", "target_type_curated", "target_type_curated_extra",
    "target_chembl_id", "target_chembl_id_curated", "target_name_curated",
    "bao_label", "source_label", "strain", "atcc_id", "mutations",
    "known_drug_resistances", "media",
    "activity_type", "unit", "activities", "nan_values", "cpds", "frac_cs",
    "direction", "act_flag", "inact_flag",
    "equal", "higher", "lower",
    "dataset_type", "evaluated_cutoffs",
    "min_", "p1", "p25", "p50", "p75", "p99", "max_",
    "comment_A", "comment_B", "comment_M",
]

assays_master[all_cols].to_csv(os.path.join(OUTPUT, "18_assays_master.csv"), index=False)
print("Assay master table saved.")
