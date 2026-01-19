from collections import defaultdict
from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
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

print("Step 07")

# Loading ChEMBL preprocessed data
print("Loading ChEMBL preprocessed data...")
ChEMBL = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "activities_preprocessed.csv"), low_memory=False)
print(f"Original size: {len(ChEMBL)}")

def get_pathogen_data(ChEMBL, pathogen, manual_assays):
    ChEMBL_pathogen = ChEMBL[ChEMBL['target_organism'].str.contains(pathogen, case=False, na=False) | 
                             ChEMBL['assay_organism'].str.contains(pathogen, case=False, na=False) |
                             ChEMBL['assay_chembl_id'].isin(manual_assays)].reset_index(drop=True)
    return ChEMBL_pathogen

def only_one(values, name):
    if len(values) != 1:
        raise ValueError(f"Expected exactly one {name}, found {values}")
    return values[0]

# Define variables
print(f"Filtering for pathogen: {pathogen}...")
PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code)
os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

# Read manual assays provided by the user
manual_assays = open(os.path.join(CONFIGPATH, 'assays', f"{pathogen_code}.csv"), "r").read()
manual_assays = set([j for i in manual_assays.split("\n") for j in i.split(",")])

# Get pathogen data
ChEMBL_pathogen = get_pathogen_data(ChEMBL, pathogen, manual_assays)
ChEMBL_pathogen.to_csv(os.path.join(PATH_TO_OUTPUT, f"{pathogen_code}_ChEMBL_raw_data.csv.gz"), index=False)
print(f"Number of activities: {len(ChEMBL_pathogen)}")
print(f"Number of unique compounds: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")

# Get organism metadata
df = dict(Counter(ChEMBL_pathogen['target_organism']))
df = pd.DataFrame([[i, df[i]] for i in sorted(df, key = lambda x: df[x], reverse=True)], columns=['organism', 'count'])
df.to_csv(os.path.join(PATH_TO_OUTPUT, "target_organism_counts.csv"), index=False)

# Get compound metadata and counts
compound_info = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "compound_info.csv"), low_memory=False)
ik_dict = dict(zip(compound_info["chembl_id"], compound_info["standard_inchi_key"]))
pair_counts = ChEMBL_pathogen[['compound_chembl_id', 'canonical_smiles']].value_counts().reset_index(name='count')
pair_counts["InChIKey"] = pair_counts["compound_chembl_id"].map(ik_dict)
cols = ["compound_chembl_id", 'InChIKey', 'canonical_smiles', 'count']
pair_counts[cols].to_csv(os.path.join(PATH_TO_OUTPUT, "compound_counts.csv.gz"), index=False)

# Define chemical space
pathogen_chemical_space = set(pair_counts['compound_chembl_id'])

# Get unique assays
assays = sorted(set(ChEMBL_pathogen['assay_chembl_id']))

# Get assay to index mapping
assay_to_idx = defaultdict(list)
for i, assay_id in enumerate(ChEMBL_pathogen["assay_chembl_id"].to_numpy()):
    assay_to_idx[assay_id].append(i)

ASSAYS_INFO = []
print("Collecting individual assay information...")

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


            ASSAYS_INFO.append([assay, assay_type, assay_organism, doc_chembl_id, target_type, target_chembl_id, target_organism, activity_type, unit, 
                                activities, nan_values, len(cpds), act_flag, inact_flag, frac_cs])


ASSAYS_INFO = pd.DataFrame(ASSAYS_INFO, columns=["assay_id", "assay_type", "assay_organism", "doc_chembl_id", "target_type", "target_chembl_id", "target_organism", 
                                                "activity_type", "unit", "activities", 'nan_values', "cpds", "act_flag", "inact_flag", "frac_cs"])

# Sort assays by compound count
ASSAYS_INFO = ASSAYS_INFO.sort_values('cpds', ascending=False).reset_index(drop=True)


# Save assays info
ASSAYS_INFO.to_csv(os.path.join(PATH_TO_OUTPUT, 'assays_raw.csv'), index=False)