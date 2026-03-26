from collections import defaultdict, Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH
from pathogen_utils import load_pathogen, load_expert_cutoffs, extra_curation_target_type, add_target_type_curated
from dataset_utils import (
    get_assay_data, get_cut_value, count_relations,
    get_assay_data_quantitative, get_assay_data_qualitative,
    adjust_relation, disambiguate_compounds, binarize_with_expert_cutoff,
    set_variables_quantitative, set_variables_qualitative,
    get_activity_stats_quantitative, zip_and_remove, make_dataset_filename,
)

# Load pathogen info
pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print("Step 12")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

# Load and enrich assay table with curated target types
assays_cleaned = pd.read_csv(os.path.join(OUTPUT, "08_assays_cleaned.csv"))
parameters = pd.read_csv(os.path.join(OUTPUT, "09_assays_parameters_full.csv"))
assays_cleaned = add_target_type_curated(assays_cleaned, parameters)
assays_cleaned['target_type_curated_extra'] = [
    extra_curation_target_type(i, j)
    for i, j in zip(assays_cleaned['target_type'], assays_cleaned['target_type_curated'])
]

# Load activity records and expert cutoffs
os.makedirs(os.path.join(OUTPUT, 'datasets'), exist_ok=True)
print(f"Loading ChEMBL cleaned data for {pathogen_code}...")
chembl_pathogen = pd.read_csv(os.path.join(OUTPUT, "08_chembl_cleaned_data.csv.gz"), low_memory=False)
print(f"  Activities : {len(chembl_pathogen)}")
print(f"  Compounds  : {len(set(chembl_pathogen['compound_chembl_id']))}")
print(f"  Assays     : {len(assays_cleaned)}")

expert_cutoffs_map = load_expert_cutoffs(CONFIGPATH)


def _make_dataset_row(assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra,
                      activities, nan_values, cpds, direction, act_flag, inact_flag,
                      dataset_type, expert_cutoff,
                      pos_qt, ratio_qt, cpds_qt,
                      pos_ql, ratio_ql, cpds_ql,
                      overlap_mx, pos_mx, ratio_mx, cpds_mx):
    return [
        assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra,
        activities, nan_values, cpds, direction, act_flag, inact_flag,
        dataset_type, expert_cutoff,
        pos_qt, ratio_qt, cpds_qt,
        pos_ql, ratio_ql, cpds_ql,
        overlap_mx, pos_mx, ratio_mx, cpds_mx,
    ]


# Pre-index activity records by assay for fast lookup
assay_to_idx = defaultdict(list)
for i, assay_id in enumerate(chembl_pathogen["assay_chembl_id"].to_numpy()):
    assay_to_idx[assay_id].append(i)

datasets = []
assay_data_info = []

print("Preparing datasets for each assay...")
cols = ['compound_chembl_id', 'smiles', 'activity_type', 'value', 'relation', 'unit', 'text_flag']

for assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra, \
        activities, nan_values, cpds, direction, act_flag, inact_flag in tqdm(
        assays_cleaned[['assay_id', 'activity_type', 'unit', 'target_type',
                         'target_type_curated_extra', 'activities', 'nan_values', 'cpds',
                         'direction', 'act_flag', 'inact_flag']].values):

    # Extract records for this (assay_id, activity_type, unit) triplet
    tmp_df = chembl_pathogen.iloc[assay_to_idx[assay_chembl_id]]
    assay_data = get_assay_data(tmp_df, assay_chembl_id, activity_type, unit, cols)



    equal, lower, higher = count_relations(assay_data)

    # Default mixed stats
    positives_mixed, ratio_mixed, compounds_mixed, overlap_mixed = [np.nan] * 4

    # --- Qualitative dataset ---
    assay_data_qualitative = get_assay_data_qualitative(assay_data)
    if len(assay_data_qualitative) > 0:
        positives_qualitative, ratio_qualitative, compounds_qualitative = set_variables_qualitative(assay_data_qualitative)
    else:
        positives_qualitative, ratio_qualitative, compounds_qualitative = [np.nan] * 3

    # --- Quantitative dataset ---
    assay_data_quantitative_clean = get_assay_data_quantitative(assay_data)
    cutoffs = expert_cutoffs_map.get((activity_type, unit, target_type_curated_extra, pathogen_code), [np.nan])

    if len(assay_data_quantitative_clean) == 0:

        if np.isnan(compounds_qualitative):
            raise ValueError(f"Assay {assay_chembl_id} ({activity_type}, {unit}) has neither numeric values nor activity flags.")

        dataset_type = 'qualitative'
        positives_quantitative, ratio_quantitative, compounds_quantitative, activities_quantitative = [np.nan] * 4
        min_, p1, p25, p50, p75, p99, max_ = [np.nan] * 7
        expert_cutoff = np.nan

        assay_data_qualitative.to_csv(os.path.join(OUTPUT, 'datasets', make_dataset_filename(assay_chembl_id, activity_type, unit, 'qualitative')), index=False)
        datasets.append(_make_dataset_row(
            assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra,
            activities, nan_values, cpds, direction, act_flag, inact_flag,
            dataset_type, expert_cutoff,
            positives_quantitative, ratio_quantitative, compounds_quantitative,
            positives_qualitative, ratio_qualitative, compounds_qualitative,
            overlap_mixed, positives_mixed, ratio_mixed, compounds_mixed,
        ))

    else:

        for expert_cutoff in cutoffs:
            assay_data_quantitative = assay_data_quantitative_clean.copy()

            if np.isnan(expert_cutoff) or direction not in [-1, +1]:
                # Cannot binarize quantitatively
                compounds_quantitative = len(set(assay_data_quantitative['compound_chembl_id']))
                activities_quantitative = assay_data_quantitative['value'].tolist()
                min_, p1, p25, p50, p75, p99, max_ = get_activity_stats_quantitative(activities_quantitative)
                positives_quantitative = ratio_quantitative = np.nan

                if np.isnan(compounds_qualitative):
                    dataset_type = 'none'
                else:
                    dataset_type = 'qualitative'
                    assay_data_qualitative.to_csv(os.path.join(OUTPUT, 'datasets', make_dataset_filename(assay_chembl_id, activity_type, unit, 'qualitative')), index=False)

            else:
                # Full quantitative binarization
                cut = get_cut_value(assay_data_quantitative, direction)
                assay_data_quantitative = adjust_relation(assay_data_quantitative, direction, cut)
                assay_data_quantitative = disambiguate_compounds(assay_data_quantitative, direction)
                assay_data_quantitative = binarize_with_expert_cutoff(assay_data_quantitative, expert_cutoff, direction)
                positives_quantitative, ratio_quantitative, compounds_quantitative, activities_quantitative = set_variables_quantitative(assay_data_quantitative)
                min_, p1, p25, p50, p75, p99, max_ = get_activity_stats_quantitative(activities_quantitative)

                if np.isnan(compounds_qualitative):
                    dataset_type = 'quantitative'
                    assay_data_quantitative.to_csv(os.path.join(OUTPUT, 'datasets', make_dataset_filename(assay_chembl_id, activity_type, unit, 'quantitative', expert_cutoff)), index=False)

                else:
                    dataset_type = 'mixed'

                    overlap_mixed = round(
                        len(set(assay_data_quantitative['compound_chembl_id']) & set(assay_data_qualitative['compound_chembl_id'])) /
                        min(len(assay_data_quantitative), len(assay_data_qualitative)), 3
                    )

                    # Qualitative inactives not already in the quantitative set
                    qt_compounds = set(assay_data_quantitative['compound_chembl_id'])
                    assay_data_qualitative_mixed = assay_data_qualitative[
                        (~assay_data_qualitative['compound_chembl_id'].isin(qt_compounds)) &
                        (assay_data_qualitative['bin'] == 0)
                    ].reset_index(drop=True).copy()
                    assay_data_qualitative_mixed[['value', 'relation']] = np.nan

                    assay_data_quantitative_mixed = assay_data_quantitative.copy()
                    assay_data_quantitative_mixed[['text_flag', 'qualitative_label']] = np.nan

                    assay_data_mixed = pd.concat([assay_data_quantitative_mixed, assay_data_qualitative_mixed], axis=0).reset_index(drop=True)
                    positives_mixed, ratio_mixed, compounds_mixed, _ = set_variables_quantitative(assay_data_mixed)

                    assay_data_mixed.to_csv(os.path.join(OUTPUT, 'datasets', make_dataset_filename(assay_chembl_id, activity_type, unit, 'mixed', expert_cutoff)), index=False)

            datasets.append(_make_dataset_row(
                assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra,
                activities, nan_values, cpds, direction, act_flag, inact_flag,
                dataset_type, expert_cutoff,
                positives_quantitative, ratio_quantitative, compounds_quantitative,
                positives_qualitative, ratio_qualitative, compounds_qualitative,
                overlap_mixed, positives_mixed, ratio_mixed, compounds_mixed,
            ))

    assay_data_info.append([
        assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra,
        activities, cpds, dataset_type, equal, higher, lower,
        min_, p1, p25, p50, p75, p99, max_,
    ])


# Save summary tables
dataset_cols = [
    "assay_id", "activity_type", "unit", "target_type", "target_type_curated_extra",
    "activities", "nan_values", "cpds", "direction", "act_flag", "inact_flag",
    "dataset_type", "expert_cutoff",
    "pos_qt", "ratio_qt", "cpds_qt",
    "pos_ql", "ratio_ql", "cpds_ql",
    "overlap_mx", "pos_mx", "ratio_mx", "cpds_mx",
]
assay_info_cols = [
    "assay_id", "activity_type", "unit", "target_type", "target_type_curated_extra",
    "activities", "cpds", "dataset_type",
    "equal", "higher", "lower",
    "min_", "p1", "p25", "p50", "p75", "p99", "max_",
]

datasets_df = pd.DataFrame(datasets, columns=dataset_cols)
assay_data_info_df = pd.DataFrame(assay_data_info, columns=assay_info_cols)

datasets_df.to_csv(os.path.join(OUTPUT, '12_datasets.csv'), index=False)
assay_data_info_df.to_csv(os.path.join(OUTPUT, '12_assay_data_info.csv'), index=False)

# Archive individual dataset files
qt, ql, mx = zip_and_remove(os.path.join(OUTPUT, "datasets"))

print(f"Total assays     : {len(assay_data_info_df)}  {dict(Counter(assay_data_info_df['dataset_type']))}")
print(f"Total datasets   : {len(datasets_df)}  {dict(Counter(datasets_df['dataset_type']))}")
print(f"Archived files   : {qt} qt, {ql} ql, {mx} mx")
