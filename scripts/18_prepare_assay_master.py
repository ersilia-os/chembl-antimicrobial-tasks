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
    "assay_organism_curated", "target_type_curated", "target_name_curated",
    "target_chembl_id_curated", "assay_strain_curated", "atcc_id", "mutations",
    "known_drug_resistances", "culture_media",
]
columns_data_info = [
    "target_type_curated_extra", "dataset_type",
    "equal", "higher", "lower",
    "min_", "p1", "p25", "p50", "p75", "p99", "max_",
]

print("Loading data...")
assays_cleaned = pd.read_csv(os.path.join(OUTPUT, "08_assays_cleaned.csv"))
assays_clusters = pd.read_csv(os.path.join(OUTPUT, "10_assays_clusters.csv"))[keys + columns_clusters]
assays_parameters = pd.read_csv(os.path.join(OUTPUT, "09_assays_parameters_full.csv"))[keys + columns_parameters]
assay_data_info = pd.read_csv(os.path.join(OUTPUT, "12_assay_data_info.csv"))[keys + columns_data_info]

print("Merging tables...")
#assays_parameters = assays_parameters.drop_duplicates(subset=keys, keep="last").reset_index(drop=True) #TODO REMOVE ONCE LLM RERUN
assays_master = assays_cleaned.merge(assays_clusters, on=keys, how="left", validate="1:1")
assays_master = assays_master.merge(assays_parameters, on=keys, how="left", validate="1:1")
assays_master = assays_master.merge(assay_data_info, on=keys, how="left", validate="1:1")

# Add evaluated cutoffs column
assays_master["evaluated_cutoffs"] = [
    ";".join(str(c) for c in expert_cutoffs[(a, u, tt, pathogen_code)])
    if (a, u, tt, pathogen_code) in expert_cutoffs else np.nan
    for a, u, tt in assays_master[["activity_type", "unit", "target_type_curated_extra"]].values
]

# Build AUROC lookup from individual modeling results (loaded later in script)
# We need to load individual_lm early to build the AUROC lookup
individual_lm = pd.read_csv(os.path.join(OUTPUT, "13_individual_LM.csv"))

# Create AUROC lookup: (assay_id, activity_type, unit) -> {cutoff: auroc}
auroc_lookup = {}
for _, row in individual_lm.iterrows():
    key = (row["assay_id"], row["activity_type"], row["unit"])
    cutoff = row["expert_cutoff"]
    auroc = row["avg"]

    if key not in auroc_lookup:
        auroc_lookup[key] = {}
    auroc_lookup[key][cutoff] = auroc

# Add evaluated aurocs column - same format as cutoffs but with 2 decimal AUROC values
assays_master["evaluated_aurocs"] = [
    ";".join(f"{auroc_lookup.get((assay_id, a, u), {}).get(c, np.nan):.2f}"
             for c in expert_cutoffs[(a, u, tt, pathogen_code)])
    if (a, u, tt, pathogen_code) in expert_cutoffs and (assay_id, a, u) in auroc_lookup else np.nan
    for assay_id, a, u, tt in assays_master[["assay_id", "activity_type", "unit", "target_type_curated_extra"]].values
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
    if df is None or df.empty:
        return set()
    subset = df if label is None else df[df["label"] == label]
    return set(_parse_assay_key(s) for row in subset["assay_keys"] for s in row.split(";"))

merged_lm = pd.read_csv(os.path.join(OUTPUT, "15_merged_LM.csv"))
individual_selected_lm = pd.read_csv(os.path.join(OUTPUT, "14_individual_selected_LM.csv"))
merged_selected_lm = pd.read_csv(os.path.join(OUTPUT, "16_merged_selected_LM.csv"))
merging_analysis = pd.read_csv(os.path.join(OUTPUT, "15_merging_analysis.csv"))
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

# Assays accepted in individual selection (for hierarchical logic comment in B and M)
accepted_individual = (
    set(_norm_key(tuple(r)) for r in individual_selected_lm[keys].values)
)

# Map (label, assay_key) -> dataset name in final_datasets (for correlation lookup)
key_to_final_name = {}
for lbl, assay_keys_str, name in final_datasets[["label", "assay_keys", "name"]].values:
    for s in assay_keys_str.split(";"):
        key_to_final_name[(lbl, _norm_key(_parse_assay_key(s)))] = name

# Create mapping from assay key to merged group name for detailed M comments
assay_to_merged_group = {}
# From merged_lm (step 15) - all merged groups attempted
for _, row in merged_lm.iterrows():
    base_name = "_".join(row["name"].split("_")[:2])  # Extract M_ORG0 from M_ORG0_10.0
    for assay_key_str in row["assay_keys"].split(";"):
        assay_key = _norm_key(_parse_assay_key(assay_key_str))
        assay_to_merged_group[assay_key] = base_name

# For each non-selected dataset, find the selected dataset it correlates most with (fix 3: no mutation)
if len(final_datasets)>0:
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

# Create merging analysis lookup for detailed M comments
merging_failure_lookup = {}
if len(merging_analysis) > 0:
    for _, row in merging_analysis.iterrows():
        key = _norm_key((row["assay_id"], row["activity_type"], row["unit"]))
        merging_failure_lookup[key] = {
            "failure_reason": row["failure_reason"],
            "group_size": row.get("group_size", 0),
            "group_compounds": row.get("group_compounds", 0),
            "n_positives": row.get("n_positives", 0),
            "merged_group_name": row.get("merged_group_name", ""),
            "auroc": row.get("auroc", 0.0)
        }

# ---------------------------------------------------------------------------
# Annotate each assay with its pipeline status per label
# ---------------------------------------------------------------------------

def pipeline_comment(key, target_type_extra, lbl, considered, selected, final_considered, final_selected):
    """Generate standardized pipeline status comments for assays."""
    nkey = _norm_key(key)
    assay_id, activity_type, unit = key

    # Check expert cutoff availability
    has_cutoff = (activity_type, unit, target_type_extra, pathogen_code) in expert_cutoffs
    dtype = dataset_type_map.get(nkey, "none")

    if lbl == "M":
        return _generate_m_comment(nkey, dtype, has_cutoff, considered, selected, final_considered, final_selected)
    elif lbl == "B":
        return _generate_b_comment(nkey, dtype, has_cutoff, considered, selected, final_considered, final_selected)
    else:  # lbl == "A"
        return _generate_a_comment(nkey, dtype, has_cutoff, considered, selected, final_considered, final_selected)


def _generate_a_comment(nkey, dtype, has_cutoff, considered, selected, final_considered, final_selected):
    """Generate standardized comments for condition A (Individual Modeling - Large, Balanced Datasets)."""

    # A1. Pre-modeling Issues
    if nkey not in considered:
        if not has_cutoff:
            return "Not considered for A: no expert cutoff defined for this (activity_type, unit, target_type) combination"
        if dtype == "none":
            return "Not considered for A: no activity data available"
        if dtype == "qualitative":
            return "Not considered for A: only qualitative data available, requires quantitative values"
        # Default to insufficient compounds/positives for remaining cases
        return "Not considered for A: insufficient compounds (<1000 total)"

    # A2. Modeling Completed - Not Selected
    if nkey not in selected:
        auroc = best_auroc_map.get(("A", nkey), np.nan)
        if not np.isnan(auroc):
            return f"Modeled but not selected: best AUROC {auroc:.3f} below 0.70 threshold"
        return "Modeled but not selected: best AUROC below 0.70 threshold"

    # A3. Selected - Not in Final
    if nkey not in final_considered:
        return "Selected but excluded from correlation analysis (non-ORGANISM target type)"

    if nkey not in final_selected:
        name = key_to_final_name.get(("A", nkey))
        cause = corr_cause.get(name)
        if cause:
            return f"Discarded: high correlation with dataset {cause}"
        return "Discarded: high correlation with another retained dataset"

    # A4. Final Selection
    return "Retained in final selection"


def _generate_b_comment(nkey, dtype, has_cutoff, considered, selected, final_considered, final_selected):
    """Generate standardized comments for condition B (Individual Modeling - Active-Enriched Datasets)."""

    # B1. Hierarchical Pre-qualification
    if nkey not in considered:
        if nkey in accepted_individual:
            return "Not considered for B: already accepted individually under condition A"

        # B2. Pre-modeling Issues (same logic as A but different prefix)
        if not has_cutoff:
            return "Not considered for B: no expert cutoff defined for this (activity_type, unit, target_type) combination"
        if dtype == "none":
            return "Not considered for B: no activity data available"
        if dtype == "qualitative":
            return "Not considered for B: only qualitative data available, requires quantitative values"
        # Default to insufficient positives/ratio for remaining cases
        return "Not considered for B: insufficient positives (<100) at all evaluated cutoffs"

    # B3. Modeling Completed - Not Selected
    if nkey not in selected:
        auroc = best_auroc_map.get(("B", nkey), np.nan)
        if not np.isnan(auroc):
            return f"Modeled but not selected: best AUROC {auroc:.3f} below 0.70 threshold"
        return "Modeled but not selected: best AUROC below 0.70 threshold"

    # B4. Selected - Not in Final
    if nkey not in final_considered:
        return "Selected but excluded from correlation analysis (non-ORGANISM target type)"

    if nkey not in final_selected:
        name = key_to_final_name.get(("B", nkey))
        cause = corr_cause.get(name)
        if cause:
            return f"Discarded: high correlation with dataset {cause}"
        return "Discarded: high correlation with another retained dataset"

    # B5. Final Selection
    return "Retained in final selection"


def _generate_m_comment(nkey, dtype, has_cutoff, considered, selected, final_considered, final_selected):
    """Generate standardized comments for condition M (Merged Group Modeling)."""

    # Get the merged group name for this assay (if any)
    merged_group = assay_to_merged_group.get(nkey)

    # M1. Already accepted individually
    if nkey in accepted_individual:
        return "Not considered for M: already accepted individually under condition A or B"

    # M2. Not considered for merging - use detailed failure analysis
    if nkey not in considered:
        if dtype == "qualitative":
            return "Not considered for M: only qualitative data available, requires quantitative values"
        if not has_cutoff:
            return "Not considered for M: no expert cutoff defined for this (activity_type, unit, target_type) combination"
        if dtype == "none":
            return "Not considered for M: no activity data available"

        # Look up specific failure reason from merging analysis
        failure_info = merging_failure_lookup.get(nkey, {})
        failure_reason = failure_info.get("failure_reason", "insufficient_compatible_assays")

        if failure_reason == "insufficient_compatible_assays":
            group_size = failure_info.get("group_size", 0)
            return f"Not modeled: insufficient compatible assays for merging ({group_size} assay{'s' if group_size != 1 else ''} in group, need ≥2)"
        elif failure_reason == "insufficient_compounds_after_merging":
            group_compounds = failure_info.get("group_compounds", 0)
            return f"Not modeled: insufficient compounds after merging ({group_compounds} compounds, need >1000)"
        elif failure_reason == "insufficient_positives_after_merging":
            n_positives = failure_info.get("n_positives", 0)
            return f"Not modeled: insufficient positives after merging ({n_positives} positives, need >50)"
        else:
            return "Not modeled: insufficient compatible assays for merging"

    # M3. Modeled but not selected
    if nkey not in selected:
        group_info = f" in group {merged_group}" if merged_group else ""
        return f"Modeled but not selected{group_info}: AUROC below 0.70 threshold"

    # M4. Selected but excluded from correlation analysis
    if nkey not in final_considered:
        group_info = f" from group {merged_group}" if merged_group else ""
        return f"Selected{group_info} but excluded from correlation analysis (non-ORGANISM target type)"

    # M5. Selected but discarded due to correlation
    if nkey not in final_selected:
        name = key_to_final_name.get(("M", nkey))
        cause = corr_cause.get(name)
        group_info = f" from group {merged_group}" if merged_group else ""
        if cause:
            return f"Discarded{group_info}: high correlation with dataset {cause}"
        return f"Discarded{group_info}: high correlation with another retained dataset"

    # M6. Final selection
    group_info = f" from group {merged_group}" if merged_group else ""
    return f"Retained in final selection{group_info}"


comments_a, comments_b, comments_m = [], [], []
for assay_id, activity_type, unit, target_type_extra in assays_master[keys + ["target_type_curated_extra"]].values:
    key = (assay_id, activity_type, unit)
    comments_a.append(pipeline_comment(key, target_type_extra, "A", considered_a, selected_a, final_considered_a, final_selected_a))
    comments_b.append(pipeline_comment(key, target_type_extra, "B", considered_b, selected_b, final_considered_b, final_selected_b))
    comments_m.append(pipeline_comment(key, target_type_extra, "M", considered_m, selected_m, final_considered_m, final_selected_m))

assays_master["comment_A"] = comments_a
assays_master["comment_B"] = comments_b
assays_master["comment_M"] = comments_m

# Create evaluation context mappings (what was tested, not just what was selected)
eval_cutoff_map = {}
eval_label_map = {}

# Map individual datasets that were evaluated (A and B conditions from step 14)
for _, row in individual_selected_lm.iterrows():
    assay_key = _norm_key(tuple(row[keys]))
    eval_cutoff_map[assay_key] = row["cutoff"]
    eval_label_map[assay_key] = row["label"]

# Map merged datasets that were evaluated (M conditions from step 16)
for _, row in merged_selected_lm.iterrows():
    base_name = "_".join(row["name"].split("_")[:2])  # Extract M_ORG0 from M_ORG0_10.0
    for assay_key_str in row["assay_keys"].split(";"):
        assay_key = _norm_key(_parse_assay_key(assay_key_str))
        eval_cutoff_map[assay_key] = row["cutoff"]
        eval_label_map[assay_key] = base_name

# Create final selection set (for selected column only)
final_selected_assays = set()
for _, row in final_datasets[final_datasets["selected"]].iterrows():
    for assay_key_str in row["assay_keys"].split(";"):
        assay_key = _norm_key(_parse_assay_key(assay_key_str))
        final_selected_assays.add(assay_key)

# Set columns based on evaluation context and final selection
assays_master["selected"] = [
    _norm_key(tuple(r)) in final_selected_assays
    for r in assays_master[keys].values
]
assays_master["selected_cutoff"] = [
    eval_cutoff_map.get(_norm_key(tuple(r)), np.nan)
    for r in assays_master[keys].values
]
assays_master["selected_label"] = [
    eval_label_map.get(_norm_key(tuple(r)), np.nan)
    for r in assays_master[keys].values
]

# ---------------------------------------------------------------------------
# Add PubChem AID column
# ---------------------------------------------------------------------------

pubchem_file = os.path.join(CONFIGPATH, "pubchem_aids", f"chembl_assays_in_pubchem_{pathogen_code}.csv")
if os.path.exists(pubchem_file):
    pubchem_df = pd.read_csv(pubchem_file)[["AID", "Source ID"]].rename(columns={"Source ID": "assay_id"})
    pubchem_map = dict(zip(pubchem_df["assay_id"], pubchem_df["AID"]))
    assays_master["pubchem"] = assays_master["assay_id"].map(pubchem_map)
else:
    assays_master["pubchem"] = np.nan

# ---------------------------------------------------------------------------
# Add UniProt accession column
# ---------------------------------------------------------------------------

CHEMBL_ACTIVITIES = os.path.join(root, "..", "data", "chembl_activities")
target_dict = pd.read_csv(os.path.join(CHEMBL_ACTIVITIES, "target_dictionary.csv"))[["chembl_id", "tid"]]
target_comp = pd.read_csv(os.path.join(CHEMBL_ACTIVITIES, "target_components.csv"))[["tid", "component_id"]]
comp_seq    = pd.read_csv(os.path.join(CHEMBL_ACTIVITIES, "component_sequences.csv"))[["component_id", "accession"]]

target_acc = (
    target_dict
    .merge(target_comp, on="tid")
    .merge(comp_seq, on="component_id")
    .groupby("chembl_id")["accession"]
    .apply(lambda x: ";".join(x.dropna().unique()))
    .reset_index()
    .rename(columns={"accession": "uniprot_accession"})
)
target_acc = target_acc[target_acc["uniprot_accession"] != ""]

chembl_id_to_acc = dict(zip(target_acc["chembl_id"], target_acc["uniprot_accession"]))
assays_master["uniprot_accession"] = assays_master["target_chembl_id_curated"].map(chembl_id_to_acc)

# ---------------------------------------------------------------------------
# Select and save final columns
# ---------------------------------------------------------------------------

all_cols = [
    "assay_id", "assay_type", "assay_organism", "target_organism",
    "assay_organism_curated", "doc_chembl_id",
    "target_type", "target_type_curated", "target_type_curated_extra",
    "target_chembl_id", "target_chembl_id_curated", "uniprot_accession", "target_name_curated",
    "bao_label", "source_label", "assay_strain_curated", "atcc_id", "mutations",
    "known_drug_resistances", "culture_media",
    "activity_type", "unit", "activities", "nan_values", "cpds", "frac_cs",
    "direction", "act_flag", "inact_flag",
    "equal", "higher", "lower",
    "dataset_type", "evaluated_cutoffs", "evaluated_aurocs", "ratio_qt", "ratio_ql", "ratio_mx",
    "min_", "p1", "p25", "p50", "p75", "p99", "max_",
    "selected", "selected_cutoff", "selected_label",
    "comment_A", "comment_B", "comment_M",
    "pubchem",
]

assays_master[all_cols].to_csv(os.path.join(OUTPUT, "18_assays_master.csv"), index=False)
print("Assay master table saved.")
