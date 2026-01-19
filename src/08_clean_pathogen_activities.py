from collections import defaultdict
from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH

# Load pathogen info
pathogen_code = sys.argv[1]
df = pd.read_csv(os.path.join(root, "..", "config", 'pathogens.csv'))
row = df.loc[df["code"].eq(pathogen_code)]
if row.empty: 
    raise SystemExit(f"Unknown code: {pathogen_code}")
pathogen = row.iloc[0]["pathogen"]

print("Step 08")

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

def only_one(values, name):
    # Helper function - is there only a single value?
    if len(values) != 1:
        raise ValueError(f"Expected exactly one {name}, found {values}")
    return values[0]

def create_text_flag(ChEMBL_pathogen):

    cond_nan = (ChEMBL_pathogen['activity_comment'] == 0) & (ChEMBL_pathogen['standard_text'] == 0)
    cond_pos = (ChEMBL_pathogen['activity_comment'] == 1) | (ChEMBL_pathogen['standard_text'] == 1)
    cond_neg = (ChEMBL_pathogen['activity_comment'] == -1) | (ChEMBL_pathogen['standard_text'] == -1)

    # Detect row-level conflicts
    conflict = cond_pos & cond_neg
    if conflict.any():
        raise ValueError(
            "Conflicting labels (contains both 1 and -1):\n"
            + ChEMBL_pathogen.loc[conflict, ["compound_chembl_id", "activity_comment", "standard_text"]].head(20).to_string())

    # Assign row-level label
    ChEMBL_pathogen["text_flag"] = np.nan
    ChEMBL_pathogen.loc[cond_pos, "text_flag"] = 1
    ChEMBL_pathogen.loc[cond_neg, "text_flag"] = -1
    ChEMBL_pathogen.loc[cond_nan, "text_flag"] = 0

    # Remove original fields
    ChEMBL_pathogen = ChEMBL_pathogen.drop(columns=['activity_comment', 'standard_text'])

    return ChEMBL_pathogen


# Loading pathogen data
print(f"Loading ChEMBL preprocessed data for {pathogen_code}...")
ChEMBL_pathogen = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, f"{pathogen_code}_ChEMBL_raw_data.csv.gz"), low_memory=False)
ASSAYS_RAW = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, 'assays_raw.csv'))
print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
print(f"Original number of assays: {len(ASSAYS_RAW)} (unique: {len(set(ASSAYS_RAW['assay_id']))})")

# Remove non std compounds

# Discard activities with no value nor text_flag
ChEMBL_pathogen = ChEMBL_pathogen[(ChEMBL_pathogen['value'].isna() == False) | 
                                (ChEMBL_pathogen['text_flag'] != 0)].reset_index(drop=True)

print(f"Removing activities with no value nor text_flag...")
print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")

print(f"Keeping only those activities with consensus units")
CONSENSUS_UNITS = set(pd.read_csv(os.path.join(root, "..", "config", 'chembl_processed', "unit_conversion.csv"))['final_unit'])
ChEMBL_pathogen = ChEMBL_pathogen[(ChEMBL_pathogen['unit'].isin(CONSENSUS_UNITS) == True) |
                                    (ChEMBL_pathogen['unit'].isna() == True)].reset_index(drop=True)
print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")

# Get directions
DIRECTIONS = pd.read_csv(os.path.join(root, "..", "config", 'manual_curation', 'activity_std_units_curated_manual_curation.csv'))
DIRECTIONS = {(i,j): k for i,j,k in zip(DIRECTIONS['activity_type'], DIRECTIONS['unit'], DIRECTIONS['manual_curation']) if np.isnan(k) == False}
ChEMBL_pathogen['direction'] = [DIRECTIONS[(i,j)] if (i,j) in DIRECTIONS else np.nan 
                                for i,j in zip(ChEMBL_pathogen['activity_type'], ChEMBL_pathogen['unit'])]
count_directions = Counter(ChEMBL_pathogen['direction'].fillna('NaN'))
print(f"Directions assigned. Summary: {count_directions}")
print(f"Assigned directions [-1, 0, +1]: {round((count_directions[1] + count_directions[-1] + count_directions[0]) / len(ChEMBL_pathogen) * 100, 1)}%")
print(f"Assigned directions [-1, +1]: {round((count_directions[1] + count_directions[-1]) / len(ChEMBL_pathogen) * 100, 1)}%")

print(f"Keeping only activities with a direction [-1,+1] OR active/inactive text_flag")
ChEMBL_pathogen = ChEMBL_pathogen[(ChEMBL_pathogen['direction'].isin([1, -1]) == True) | 
                                    (ChEMBL_pathogen['text_flag'].isin([1, -1]))].reset_index(drop=True)
print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")


# Identify activity_type - unit pairs
print("Identifying activity_type - unit pairs...")
# Get pair counts
s = ChEMBL_pathogen[["activity_type", "unit"]]
out = (
s.value_counts(subset=["activity_type", "unit"], dropna=False)
    .reset_index(name="count")
    .sort_values("count", ascending=False, ignore_index=True))
total_count = out['count'].sum()
out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)

# Assign direction to activity_type_unit_pairs
out['direction'] = [DIRECTIONS[(i,j)] if (i,j) in DIRECTIONS else np.nan 
                                for i,j in zip(out['activity_type'], out['unit'])]

tmp = ChEMBL_pathogen[["activity_type", "unit", "text_flag"]].copy()
tmp["activity_type"] = tmp["activity_type"].fillna("")

# Define interesting index
tmp["has_unit"] = ~tmp["unit"].isna()
tmp["has_comment"] = ~(tmp['text_flag'] == 0)

# Count per activity_type x (has_unit, has_comment)
counts_long = (
    tmp.groupby(["activity_type", "has_unit", "has_comment"], dropna=False)
    .size()
    .reset_index(name="count"))

# Pivot to wide: one row per activity_type, 4 count columns
cols = ["unit_comment", "nounit_comment", "unit_nocomment", "nounit_no_comment"]
counts_wide = (
    counts_long
    .assign(
        bucket=np.select(
            [
                counts_long["has_unit"] & counts_long["has_comment"],
                ~counts_long["has_unit"] & counts_long["has_comment"],
                counts_long["has_unit"] & ~counts_long["has_comment"],
                ~counts_long["has_unit"] & ~counts_long["has_comment"],
            ],
            cols)))

counts_wide = counts_wide.pivot_table(index="activity_type", columns="bucket", values="count",
                fill_value=0, aggfunc="sum").reset_index()

# If any columns is not created (not appearences for that pathogen), create it manually with 0's
for c in cols:
    if c not in counts_wide.columns:
        counts_wide[c] = 0

# Sort by total counts
counts_wide["total_count"] = counts_wide[cols].sum(axis=1)
counts_wide = counts_wide.sort_values("total_count", ascending=False, ignore_index=True)

# Save pair summary
out.to_csv(os.path.join(root, "..", "output", pathogen_code, "activity_type_unit_pairs.csv"), index=False)

# Save pair summary
counts_wide.to_csv(os.path.join(root, "..", "output", pathogen_code, "activity_type_unit_comment.csv"), index=False)

# Save cleaned data
ChEMBL_pathogen.to_csv(os.path.join(root, "..", "output", pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), index=False)

# Get unique assays
assays = sorted(set(ChEMBL_pathogen['assay_chembl_id']))

ASSAYS_INFO = []
print("Collecting individual assay information...")
print(f"Number of unique assays: {len(assays)}")

# Get assay to index mapping
assay_to_idx = defaultdict(list)
for i, assay_id in enumerate(ChEMBL_pathogen["assay_chembl_id"].to_numpy()):
    assay_to_idx[assay_id].append(i)

# For each assay
for assay in tqdm(assays):


    # Get subset of assay data
    df_ = ChEMBL_pathogen.iloc[assay_to_idx[assay]]
    
    # Get values
    assay_type = list(set(df_['assay_type']))
    target_type = list(set(df_['target_type']))
    target_chembl_id = list(set(df_['target_chembl_id']))
    activity_types = list(set(df_['activity_type']))
    target_organism = list(set(df_['target_organism']))
    assay_organism = list(set(df_['assay_organism']))
    doc_chembl_id = list(set(df_['doc_chembl_id']))

    # Check coherence
    assay_type = only_one(assay_type, "assay_type")
    target_type = only_one(target_type, "target_type")
    target_chembl_id = only_one(target_chembl_id, "target_chembl_id")
    target_organism = only_one(target_organism, "target_organism")
    assay_organism = only_one(assay_organism, "assay_organism")
    doc_chembl_id = only_one(doc_chembl_id, "doc_chembl_id")

    # For each activity type
    for act_type in activity_types:

        df__ = df_[df_["activity_type"] == act_type]
        activity_type = list(set(df__['activity_type']))
        activity_type = only_one(activity_type, 'activity_type')
        units = list(set(df__['unit']))

        for u in units:


            # If unit is nan
            if pd.isna(u):
                df___ = df__[df__["unit"].isna()]
            else:
                df___ = df__[df__["unit"] == u]

            # Get metadata for that assay
            unit = list(set(df___['unit']))
            unit = only_one(unit, "unit")
            activities = len(df___)
            cpds = len(set(df___['compound_chembl_id']))
            nan_values = len(df___[df___['value'].isna()])
            direction = DIRECTIONS[(act_type, unit)] if (act_type, unit) in DIRECTIONS else np.nan
            text_flag = Counter(df___['text_flag'].tolist())
            text_flag = text_flag[-1] + text_flag[1]
            ASSAYS_INFO.append([assay, assay_type, assay_organism, doc_chembl_id, target_type, target_chembl_id, target_organism, activity_type, 
                                unit, activities, nan_values, cpds, direction, text_flag])
            

ASSAYS_INFO = pd.DataFrame(ASSAYS_INFO, columns=["assay_id", "assay_type", "assay_organism", "doc_chembl_id", "target_type", "target_chembl_id", "target_organism", "activity_type", 
                                                    "unit", "activities", 'nan_values', "cpds", "direction", "text_flag_counts"])
ASSAYS_INFO = ASSAYS_INFO.sort_values('cpds', ascending=False).reset_index(drop=True)

# Filter assays with too few compounds
ASSAYS_INFO = ASSAYS_INFO[ASSAYS_INFO['cpds'] > MIN_ASSAY_SIZE].reset_index(drop=True)

# Save assays info
ASSAYS_INFO.to_csv(os.path.join(root, "..", "output", pathogen_code, 'assays_cleaned.csv'), index=False)