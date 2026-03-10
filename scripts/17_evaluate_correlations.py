from collections import defaultdict, Counter
from scipy.stats import spearmanr
import numpy as np
import pandas as pd
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from pathogen_utils import load_pathogen
from dataset_utils import make_dataset_filename

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print("Step 17: Evaluating correlations among datasets")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)


def _parse_assay_key(s):
    assay_id, activity_type, unit = s.split("|")
    return (assay_id, activity_type, np.nan if unit == "" else unit)
path_to_correlations = os.path.join(OUTPUT, "correlations")
labels = ["A", "B", "M"]

# ---------------------------------------------------------------------------
# Load reference set predictions for all models
# ---------------------------------------------------------------------------

probs_ref = {}
for label in labels:
    probs_ref[label] = {}
    label_dir = os.path.join(path_to_correlations, label)
    if os.path.exists(label_dir):
        for fname in sorted(os.listdir(label_dir)):
            name = fname.replace("_ref_probs.npz", "")
            probs_ref[label][name] = np.load(os.path.join(label_dir, fname))["y_prob_ref"]

# ---------------------------------------------------------------------------
# Build assay → compound set mapping
# ---------------------------------------------------------------------------

chembl = pd.read_csv(os.path.join(OUTPUT, "08_chembl_cleaned_data.csv.gz"), low_memory=False)
assay_to_compounds = defaultdict(set)
for assay_id, activity_type, unit, compound_chembl_id in chembl[["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"]].values:
    assay_to_compounds[(assay_id, activity_type, unit)].add(compound_chembl_id)
del chembl

pathogen_compounds = set(pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"])

# ---------------------------------------------------------------------------
# Build unified final_datasets table (individual A/B + merged M)
# Only ORGANISM assays are included in the correlation analysis.
# ---------------------------------------------------------------------------

individual_selected = pd.read_csv(os.path.join(OUTPUT, "14_individual_selected_LM.csv"))
merged_selected = pd.read_csv(os.path.join(OUTPUT, "16_merged_selected_LM.csv"))
print(f"Total datasets: {len(individual_selected) + len(merged_selected)}")

individual_selected = individual_selected[individual_selected["target_type"] == "ORGANISM"].reset_index(drop=True)
merged_selected = merged_selected[merged_selected["target_type"] == "ORGANISM"].reset_index(drop=True)
print(f"ORGANISM datasets: {len(individual_selected) + len(merged_selected)}")

# Compound and positive counts: take the max across quantitative and mixed columns
individual_selected["cpds"] = individual_selected[["cpds_qt", "cpds_mx"]].max(axis=1)
individual_selected["positives"] = individual_selected[["pos_qt", "pos_mx"]].max(axis=1)

# Build canonical model names matching the filenames saved in step 13
individual_selected["name"] = [
    make_dataset_filename(assay_id, activity_type, unit, dataset_type, cutoff).replace(".csv.gz", "")
    for assay_id, activity_type, unit, dataset_type, cutoff
    in individual_selected[["assay_id", "activity_type", "unit", "dataset_type", "cutoff"]].values
]

# Normalise individual table to match merged schema
ind_tmp = individual_selected.drop(columns=[
    "assay_id", "is_mid_cutoff",
    "pos_qt", "ratio_qt", "cpds_qt",
    "pos_ql", "ratio_ql", "cpds_ql",
    "pos_mx", "ratio_mx", "cpds_mx", "overlap_mx",
]).copy()
ind_tmp["assay_keys"] = [
    "|".join(r) for r in individual_selected[["assay_id", "activity_type", "unit"]].values.astype(str)
]
ind_tmp["n_assays"] = 1

mrg_tmp = merged_selected.drop(columns=[
    "direction", "assay_type", "bao_label", "is_mid_cutoff",
    "ratio", "avg", "std", "strain", "target_chembl_id",
]).rename(columns={"n_cpds_union": "cpds"}).copy()
mrg_tmp["label"] = "M"
mrg_tmp["dataset_type"] = np.nan

final_datasets = pd.concat([ind_tmp, mrg_tmp], ignore_index=True)
final_datasets = final_datasets.sort_values(["label", "cpds"], ascending=[True, False]).reset_index(drop=True)

# ---------------------------------------------------------------------------
# Map each dataset name to its assay keys and compound set
# ---------------------------------------------------------------------------

name_to_assay_keys = {
    name: [_parse_assay_key(s) for s in assay_keys.split(";")]
    for name, assay_keys in zip(final_datasets["name"], final_datasets["assay_keys"])
}

name_to_compounds = {
    name: set(cpd for assay in assays for cpd in assay_to_compounds[assay])
    for name, assays in name_to_assay_keys.items()
}

# ---------------------------------------------------------------------------
# Pairwise correlation metrics
# ---------------------------------------------------------------------------

def hit_overlap_chance(probs1, probs2, top=100):
    """Normalised hit overlap above chance for the top-N predictions."""
    n = len(probs1)
    ind1 = set(np.argsort(probs1)[::-1][:top])
    ind2 = set(np.argsort(probs2)[::-1][:top])
    m = len(ind1 & ind2)
    expected = top * top / n
    return (m - expected) / (top - expected)

dataset_labels = final_datasets[["label", "name"]].values.tolist()

results = []
results_dict = {}
for label1, name1 in dataset_labels:
    for label2, name2 in dataset_labels:
        probs1 = probs_ref[label1][name1]
        probs2 = probs_ref[label2][name2]
        cpds1 = name_to_compounds[name1]
        cpds2 = name_to_compounds[name2]
        sp = round(spearmanr(probs1, probs2).statistic, 4)
        ho1000 = round(hit_overlap_chance(probs1, probs2, top=1000), 4)
        ho100 = round(hit_overlap_chance(probs1, probs2, top=100), 4)
        co = round(len(cpds1 & cpds2) / min(len(cpds1), len(cpds2)), 4)
        results.append([label1, name1, label2, name2, sp, ho1000, ho100, co])
        results_dict[(label1, name1, label2, name2)] = (sp, ho1000, ho100, co)

results_df = pd.DataFrame(results, columns=[
    "strategy_1", "name_1", "strategy_2", "name_2",
    "spearman", "hit_overlap_1000", "hit_overlap_100", "compound_overlap",
])
results_df.to_csv(os.path.join(OUTPUT, "17_dataset_correlations.csv"), index=False)

# ---------------------------------------------------------------------------
# Greedy dataset selection: discard datasets too correlated with already-selected ones
# ---------------------------------------------------------------------------

print("Prioritizing datasets...")
selected = []
for label, name in dataset_labels:
    keep = True
    for prev_label, prev_name in selected:
        sp, ho1000, ho100, co = results_dict[(label, name, prev_label, prev_name)]
        if (sp + ho1000 + ho100) / 3 > 0.5 and co > 0.5:
            keep = False
            break
    if keep:
        selected.append([label, name])

final_datasets["selected"] = [
    [label, name] in selected
    for label, name in final_datasets[["label", "name"]].values.tolist()
]
final_datasets.to_csv(os.path.join(OUTPUT, "17_final_datasets.csv"), index=False)

print(f"Selected datasets per label: {dict(Counter(final_datasets[final_datasets['selected']]['label']))}")

# ---------------------------------------------------------------------------
# Coverage summary
# ---------------------------------------------------------------------------

coverage = {label: set() for label in labels}
for label, assay_keys in final_datasets[["label", "assay_keys"]].values:
    for s in assay_keys.split(";"):
        coverage[label].update(assay_to_compounds[_parse_assay_key(s)])

for label in labels:
    pct = round(100 * len(coverage[label]) / len(pathogen_compounds), 1)
    print(f"  Coverage {label}: {pct}%")
all_covered = coverage["A"] | coverage["B"] | coverage["M"]
print(f"  Coverage ALL: {round(100 * len(all_covered) / len(pathogen_compounds), 1)}%")
