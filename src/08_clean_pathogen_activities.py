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
df = pd.read_csv(os.path.join(CONFIGPATH, 'pathogens.csv'))
row = df.loc[df["code"].eq(pathogen_code)]
if row.empty: 
    raise SystemExit(f"Unknown code: {pathogen_code}")
pathogen = row.iloc[0]["pathogen"]

print("Step 08")

# Define output directory
OUTPUT = os.path.join(root, "..", "output")

def only_one(values, name):
    if len(values) != 1:
        raise ValueError(f"Expected exactly one {name}, found {values}")
    return values[0]

# Loading pathogen data
print(f"Loading ChEMBL preprocessed data for {pathogen_code}...")
ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_raw_data.csv.gz"), low_memory=False)
ASSAYS_RAW = pd.read_csv(os.path.join(OUTPUT, pathogen_code, 'assays_raw.csv'))
print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
print(f"Original number of assays: {len(ASSAYS_RAW)} (unique: {len(set(ASSAYS_RAW['assay_id']))})")

# Remove non std compounds
ChEMBL_pathogen = ChEMBL_pathogen[ChEMBL_pathogen['standardized_smiles'].isna() == False].reset_index(drop=True)
print("Removing activities having no standardized SMILES")
print(f"Number of activities (compounds): {len(ChEMBL_pathogen)} ({len(set(ChEMBL_pathogen['compound_chembl_id']))})")

# Discard activities with no value nor text_flag
ChEMBL_pathogen = ChEMBL_pathogen[(ChEMBL_pathogen['value'].isna() == False) | 
                                (ChEMBL_pathogen['text_flag'] != 0)].reset_index(drop=True)

print(f"Removing activities with no value nor text_flag...")
print(f"Number of activities (compounds): {len(ChEMBL_pathogen)} ({len(set(ChEMBL_pathogen['compound_chembl_id']))})")

# Discard non consensus units
CONSENSUS_UNITS = set(pd.read_csv(os.path.join(DATAPATH, 'chembl_processed', "unit_conversion.csv"))['final_unit'])
ChEMBL_pathogen = ChEMBL_pathogen[(ChEMBL_pathogen['unit'].isin(CONSENSUS_UNITS) == True) |
                                    (ChEMBL_pathogen['unit'].isna() == True)].reset_index(drop=True)
print(f"Keeping only those activities with consensus units")
print(f"Number of activities (compounds): {len(ChEMBL_pathogen)} ({len(set(ChEMBL_pathogen['compound_chembl_id']))})")

# Get directions
DIRECTIONS = pd.read_csv(os.path.join(CONFIGPATH, 'activity_std_units_curated_manual_curation.csv'))
DIRECTIONS = {(i,j): k for i,j,k in zip(DIRECTIONS['activity_type'], DIRECTIONS['unit'], DIRECTIONS['manual_curation_direction']) if np.isnan(k) == False}
ChEMBL_pathogen['direction'] = [DIRECTIONS[(i,j)] if (i,j) in DIRECTIONS else np.nan 
                                for i,j in zip(ChEMBL_pathogen['activity_type'], ChEMBL_pathogen['unit'])]
count_directions = Counter(ChEMBL_pathogen['direction'].fillna('NaN'))
print(f"Directions assigned. Summary: {count_directions}")
print(f"Assigned directions [-1, 0, +1]: {round((count_directions[1] + count_directions[-1] + count_directions[0]) / len(ChEMBL_pathogen) * 100, 1)}%")
print(f"Assigned directions [-1, +1]: {round((count_directions[1] + count_directions[-1]) / len(ChEMBL_pathogen) * 100, 1)}%")

# Keeping only directed activities
ChEMBL_pathogen = ChEMBL_pathogen[(ChEMBL_pathogen['direction'].isin([1, -1]) == True) | 
                                    (ChEMBL_pathogen['text_flag'].isin([1, -1]))].reset_index(drop=True)
print(f"Keeping only activities with a direction [-1,+1] OR active/inactive text_flag")
print(f"Number of activities (compounds): {len(ChEMBL_pathogen)} ({len(set(ChEMBL_pathogen['compound_chembl_id']))})")

# Saving cleaned pathogen activities
ChEMBL_pathogen.to_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), index=False)


print("Preparing activity - unit - text_comments report...")
cols = ["activity_type", "unit", "text_flag"]
df = ChEMBL_pathogen[cols]
flagged = df["text_flag"].isin([1, -1])

# Counting text flag
out = (df.assign(flagged=flagged.to_numpy())
      .groupby(["activity_type", "unit"], dropna=False)
      .agg(count=("flagged", "size"), comments=("flagged", "sum"))
      .reset_index()
      .sort_values("count", ascending=False, ignore_index=True))

# Cumulative proportion
total_count = out['count'].sum()
out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)

# Assign direction to activity_type_unit_pairs
out['direction'] = [DIRECTIONS[(i,j)] if (i,j) in DIRECTIONS else np.nan 
                                for i,j in zip(out['activity_type'], out['unit'])]

out.to_csv(os.path.join(OUTPUT, pathogen_code, "activity_type_unit_comments.csv"), index=False)

# Define chemical space
pathogen_chemical_space = set(pd.read_csv(os.path.join(OUTPUT, pathogen_code, "compound_counts.csv.gz"))['compound_chembl_id'])

# Get unique assays
assays = sorted(set(ChEMBL_pathogen['assay_chembl_id']))

# Get assay to source mapping
assay_data = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "assays.csv"), low_memory=False)[['chembl_id', 'src_id', 'bao_format']]
assay_to_src_id = {i: j for i,j in zip(assay_data['chembl_id'], assay_data['src_id'])}
assay_to_bao_format = {i: j for i,j in zip(assay_data['chembl_id'], assay_data['bao_format'])}

# Get source mappings
source = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "source.csv"))
src_id_to_src_short_name = {i: j for i,j in zip(source['src_id'], source['src_short_name'])}

# Get bioassay ontology mappings
BAO = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "bioassay_ontology.csv"))
bao_id_to_label = {i: j for i,j in zip(BAO['bao_id'], BAO['label'])}

# Get assay to index mapping
assay_to_idx = defaultdict(list)
for i, assay_id in enumerate(ChEMBL_pathogen["assay_chembl_id"].to_numpy()):
    assay_to_idx[assay_id].append(i)

ASSAYS_INFO = []
print("Cleaning individual assay information...")

# For each assay
for assay in tqdm(assays):

    # Get subset of assay data
    df_ = ChEMBL_pathogen.iloc[assay_to_idx[assay]]
    
    # Get values
    assay_type = list(set(df_['assay_type']))
    target_type = list(set(df_['target_type']))
    target_chembl_id = list(set(df_['target_chembl_id']))
    activity_types = list(set(df_['activity_type']))  # may be more than one
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

    # Get metadata
    bao_label = bao_id_to_label[assay_to_bao_format[assay]]
    source_label = src_id_to_src_short_name[assay_to_src_id[assay]]

    # For each activity type
    for act_type in activity_types:

        df__ = df_[df_["activity_type"] == act_type]
        activity_type = list(set(df__['activity_type']))
        activity_type = only_one(activity_type, 'activity_type')
        units = list(set(df__['unit']))  # may be more than one

        for u in units:

            # If unit is nan
            if pd.isna(u):
                df___ = df__[df__["unit"].isna()]
            else:
                df___ = df__[df__["unit"] == u]

            unit = list(set(df___['unit']))
            unit = only_one(unit, "unit")

            # Get metadata for that assay
            activities = len(df___)
            cpds = set(df___['compound_chembl_id'])
            nan_values = len(df___[df___['value'].isna()])
            text_flag = Counter(df___['text_flag'])
            act_flag = text_flag[1]
            inact_flag = text_flag[-1]
            inc_flag = text_flag[0]
            frac_cs = round(len(cpds.intersection(pathogen_chemical_space)) / len(pathogen_chemical_space), 5)
            direction = DIRECTIONS[(act_type, unit)] if (act_type, unit) in DIRECTIONS else np.nan

            
            ASSAYS_INFO.append([assay, assay_type, assay_organism, doc_chembl_id, target_type, target_chembl_id, target_organism, bao_label, source_label, 
                                activity_type, unit, activities, nan_values, len(cpds), act_flag, inact_flag, frac_cs, direction])


ASSAYS_INFO = pd.DataFrame(ASSAYS_INFO, columns=["assay_id", "assay_type", "assay_organism", "doc_chembl_id", "target_type", "target_chembl_id", "target_organism", "bao_label", 
                                                 "source_label", "activity_type", "unit", "activities", 'nan_values', "cpds", "act_flag", "inact_flag", "frac_cs", "direction"])

# Sort assays by compound count
ASSAYS_INFO = ASSAYS_INFO.sort_values('cpds', ascending=False).reset_index(drop=True)


# Save assays info
ASSAYS_INFO.to_csv(os.path.join(OUTPUT, pathogen_code, 'assays_cleaned.csv'), index=False)