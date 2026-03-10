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
assays_parameters = assays_parameters.drop_duplicates(subset=keys, keep="last").reset_index(drop=True) #TODO REMOVE ONCE LLM RERUN
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
# Add ratio columns from 12_datasets.csv
# ---------------------------------------------------------------------------

datasets_12 = pd.read_csv(os.path.join(OUTPUT, "12_datasets.csv")).sort_values(keys + ["expert_cutoff"])


def _agg_ratios(grp):
    qt_vals = grp["ratio_qt"].dropna()
    ql_vals = grp["ratio_ql"].dropna()
    mx_vals = grp["ratio_mx"].dropna()
    return pd.Series({
        "ratio_qt": ";".join(str(round(v, 5)) for v in qt_vals) if len(qt_vals) > 0 else np.nan,
        "ratio_ql": round(float(ql_vals.iloc[0]), 5) if len(ql_vals) > 0 else np.nan,
        "ratio_mx": round(float(mx_vals.iloc[0]), 5) if len(mx_vals) > 0 else np.nan,
    })


ratio_cols = datasets_12.groupby(keys, dropna=False).apply(_agg_ratios).reset_index()
assays_master = assays_master.merge(ratio_cols, on=keys, how="left")

# ---------------------------------------------------------------------------
# Load modeling results to annotate each assay's pipeline status
# ---------------------------------------------------------------------------

def _parse_assay_key(s):
    assay_id, activity_type, unit = s.split("|")
    return (assay_id, activity_type, np.nan if unit == "" else unit)


def assay_keys_to_set(df, label=None):
    """Extract all (assay_id, activity_type, unit) tuples from an assay_keys column."""
    subset = df if label is None else df[df["label"] == label]
    return set(_parse_assay_key(s) for row in subset["assay_keys"] for s in row.split(";"))

individual_lm = pd.read_csv(os.path.join(OUTPUT, "13_individual_LM.csv"))
merged_lm = pd.read_csv(os.path.join(OUTPUT, "15_merged_LM.csv"))
individual_selected_lm = pd.read_csv(os.path.join(OUTPUT, "14_individual_selected_LM.csv"))
merged_selected_lm = pd.read_csv(os.path.join(OUTPUT, "16_merged_selected_LM.csv"))
final_datasets = pd.read_csv(os.path.join(OUTPUT, "17_final_datasets.csv"))
correlations = pd.read_csv(os.path.join(OUTPUT, "17_dataset_correlations.csv"))  # fix 1

considered_a = set(tuple(r) for r in individual_lm[individual_lm["label"] == "A"][keys].values)
considered_b = set(tuple(r) for r in individual_lm[individual_lm["label"] == "B"][keys].values)
considered_m = assay_keys_to_set(merged_lm)

selected_a = set(tuple(r) for r in individual_selected_lm[individual_selected_lm["label"] == "A"][keys].values)
selected_b = set(tuple(r) for r in individual_selected_lm[individual_selected_lm["label"] == "B"][keys].values)
selected_m = assay_keys_to_set(merged_selected_lm)

cutoff_map = {tuple(r[:3]): r[3] for r in individual_selected_lm[keys + ["cutoff"]].values}
label_map = {tuple(r[:3]): r[0] for r in individual_selected_lm[["label"] + keys].values}

final_considered_a = assay_keys_to_set(final_datasets, label="A")
final_considered_b = assay_keys_to_set(final_datasets, label="B")
final_considered_m = assay_keys_to_set(final_datasets, label="M")

final_only = final_datasets[final_datasets["selected"]].reset_index(drop=True)
final_selected_a = assay_keys_to_set(final_only, label="A")
final_selected_b = assay_keys_to_set(final_only, label="B")
final_selected_m = assay_keys_to_set(final_only, label="M")

# ---------------------------------------------------------------------------
# Build lookup tables for enriched pipeline comments
# ---------------------------------------------------------------------------

def _norm_key(key):
    """Normalize NaN unit to None so the tuple is hashable and equality works."""
    a, at, u = key
    return (a, at, None if (not isinstance(u, str) and pd.isna(u)) else u)  # fix 2

# "Not modeled" reason: dataset_type and expert cutoff presence
dataset_type_map = {_norm_key(tuple(r[:3])): r[3] for r in assay_data_info[keys + ["dataset_type"]].values}

# Best AUROC per (label, assay_key) from step 13 — use _norm_key to avoid NaN!=NaN issue
best_auroc_map = {}
for lbl in ["A", "B"]:
    for (a, at, u), grp in individual_lm[individual_lm["label"] == lbl].groupby(keys, dropna=False):
        best_auroc_map[(lbl, _norm_key((a, at, u)))] = round(grp["avg"].max(), 3)

# Assays accepted in individual selection (for cascade comment in B and M)
accepted_individual = (
    set(_norm_key(tuple(r)) for r in individual_selected_lm[keys].values)
)

# Map (label, assay_key) -> dataset name in final_datasets (for correlation lookup)
key_to_final_name = {}
for lbl, assay_keys_str, name in final_datasets[["label", "assay_keys", "name"]].values:
    for s in assay_keys_str.split(";"):
        key_to_final_name[(lbl, _norm_key(_parse_assay_key(s)))] = name

# For each non-selected dataset, find the selected dataset it correlates most with (fix 3: no mutation)
selected_names = set(final_datasets[final_datasets["selected"]]["name"])
corr_cause = {}
if len(correlations) > 0:
    corr_scored = correlations.assign(
        score=(correlations["spearman"] + correlations["hit_overlap_1000"] + correlations["hit_overlap_100"]) / 3
    )
    for name, grp in corr_scored.groupby("name_1"):
        best = grp[grp["name_2"].isin(selected_names)].sort_values("score", ascending=False)
        if len(best) > 0:
            corr_cause[name] = best.iloc[0]["name_2"]

# All assays modeled in individual step (either label)
modeled_individual = (
    set(_norm_key(tuple(r)) for r in individual_lm[keys].values)
)

# ---------------------------------------------------------------------------
# Annotate each assay with its pipeline status per label
# ---------------------------------------------------------------------------

def pipeline_comment(key, target_type_extra, lbl,
                     considered, selected, final_considered, final_selected):
    nkey = _norm_key(key)
    assay_id, activity_type, unit = key

    if lbl == "M":
        if nkey not in considered:
            if nkey in accepted_individual:
                return "Not merged: accepted individually under condition A or B."
            if nkey not in modeled_individual:
                return "Not merged: not modeled individually (no expert cutoff or insufficient data)."
            return "Not merged: insufficient compound coverage or no compatible assays to group with."
        if nkey not in selected:
            return "Merged but not selected: AUROC below 0.70 threshold."
        if nkey not in final_considered:
            return "Merged and selected, not in correlation analysis (non-ORGANISM)."
        if nkey not in final_selected:
            name = key_to_final_name.get((lbl, nkey))
            cause = corr_cause.get(name)
            if cause:
                return f"Discarded: high correlation with dataset {cause}."
            return "Discarded: high correlation with another retained dataset."
        return "Retained in final selection."

    # Labels A and B
    if nkey not in considered:
        if lbl == "B" and nkey in accepted_individual:  # fix 6: cascade
            return "Not considered under condition B: already selected under condition A."
        if (activity_type, unit, target_type_extra, pathogen_code) not in expert_cutoffs:
            return "Not modeled: no expert cutoff defined for this activity/unit/target type."
        dtype = dataset_type_map.get(nkey, "none")
        if dtype == "none":
            return "Not modeled: no activity data available."
        if dtype == "qualitative":
            return "Not modeled: only qualitative data, no quantitative values to train on."
        return "Not modeled: too few positives (<50) at all evaluated cutoffs."

    if nkey not in selected:
        if lbl == "B" and nkey in accepted_individual:
            return "Not selected under condition B: already selected under condition A."
        auroc = best_auroc_map.get((lbl, nkey), np.nan)
        if not np.isnan(auroc):
            return f"Modeled but not selected: best AUROC {auroc} below 0.70 threshold."
        return "Modeled but not selected."

    if nkey not in final_considered:
        return "Selected but not in correlation analysis (non-ORGANISM target type)."

    if nkey not in final_selected:
        name = key_to_final_name.get((lbl, nkey))
        cause = corr_cause.get(name)
        if cause:
            return f"Discarded: high correlation with dataset {cause}."
        return "Discarded: high correlation with another retained dataset."

    return "Retained in final selection."


comments_a, comments_b, comments_m = [], [], []
for assay_id, activity_type, unit, target_type_extra in assays_master[keys + ["target_type_curated_extra"]].values:
    key = (assay_id, activity_type, unit)
    comments_a.append(pipeline_comment(key, target_type_extra, "A", considered_a, selected_a, final_considered_a, final_selected_a))
    comments_b.append(pipeline_comment(key, target_type_extra, "B", considered_b, selected_b, final_considered_b, final_selected_b))
    comments_m.append(pipeline_comment(key, target_type_extra, "M", considered_m, selected_m, final_considered_m, final_selected_m))

assays_master["comment_A"] = comments_a
assays_master["comment_B"] = comments_b
assays_master["comment_M"] = comments_m

assays_master["selected_cutoff"] = [cutoff_map.get(tuple(r), np.nan) for r in assays_master[keys].values]
assays_master["selected_label"] = [label_map.get(tuple(r), np.nan) for r in assays_master[keys].values]

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
    "dataset_type", "evaluated_cutoffs", "ratio_qt", "ratio_ql", "ratio_mx",
    "min_", "p1", "p25", "p50", "p75", "p99", "max_",
    "selected_cutoff", "selected_label",
    "comment_A", "comment_B", "comment_M",
]

assays_master[all_cols].to_csv(os.path.join(OUTPUT, "18_assays_master.csv"), index=False)
print("Assay master table saved.")
