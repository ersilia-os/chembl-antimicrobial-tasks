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
from default import CONFIGPATH, MIN_ASSAY_SIZE

# List of pathogens to process
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

# Helper function - is there only a single value?
def only_one(values, name):
    if len(values) != 1:
        raise ValueError(f"Expected exactly one {name}, found {values}")
    return values[0]

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

# For each pathogen
for pathogen in pathogens:

    print("\n\n\n")

    # Loading pathogen data
    pathogen_code = get_pathogen_code(pathogen)
    print(f"Loading ChEMBL preprocessed data for {pathogen_code}...")
    ChEMBL_pathogen = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, f"{pathogen_code}_ChEMBL_raw_data.csv.gz"), low_memory=False)
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
    ASSAYS_RAW = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, 'assays_raw.csv'))
    print(f"Original number of assays: {len(ASSAYS_RAW)}")

    # Converting activity types to their corresponding synonyms
    synonyms = pd.read_csv(os.path.join(root, "..", "config", "manual_curation", "synonyms.csv"))
    for activity, syns in zip(synonyms['activity'], synonyms['synonyms']):
        for syn in syns.split(";"):
            ChEMBL_pathogen.loc[ChEMBL_pathogen['activity_type'] == syn, 'activity_type'] = activity

    # Discard activities with no value nor act/inact flag in activity_comment not standard_text
    ChEMBL_pathogen = ChEMBL_pathogen[(ChEMBL_pathogen['value'].isna() == False) | 
                                    (ChEMBL_pathogen['activity_comment'] != 0) | 
                                    (ChEMBL_pathogen['standard_text'] != 0)].reset_index(drop=True)
    
    print(f"Removing activities with no value nor act/inact flag in activity_comment nor standard_test...")
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")

    # Identify canonical unit per activity type
    print("Identifying canonical unit per activity type...")
    # Get pair counts
    s = ChEMBL_pathogen[["activity_type", "unit"]]
    out = (
    s.value_counts(subset=["activity_type", "unit"], dropna=False)
        .reset_index(name="count")
        .sort_values("count", ascending=False, ignore_index=True))

    # Identify the most occurring pairs
    idx = out.groupby("activity_type")['count'].idxmax()
    out["canonical_unit"] = False
    out.loc[idx, "canonical_unit"] = True
    print(f"Number of unique activity type - unit pairs: {len(out)}")

    # Get canonical unit per activity type
    canonical = (
        out[out["canonical_unit"] == 1]
        .set_index("activity_type")[["unit"]])
    canonical_map = canonical["unit"].to_dict()
    ChEMBL_pathogen["canonical_unit"] = ChEMBL_pathogen["activity_type"].map(canonical_map)

    # Save pair summary
    out.to_csv(os.path.join(root, "..", "output", pathogen_code, "activity_type_unit_pairs.csv"), index=False)

    # Get directions
    DIRECTIONS = pd.read_csv(os.path.join(root, "..", "config", 'manual_curation', 'activity_std_units_curated_manual_curation.csv'))
    DIRECTIONS = {(i,j): k for i,j,k in zip(DIRECTIONS['activity_type'], DIRECTIONS['unit'], DIRECTIONS['manual_curation']) if np.isnan(k) == False}
    ChEMBL_pathogen['direction'] = [DIRECTIONS[(i,j)] if (i,j) in DIRECTIONS else np.nan 
                                    for i,j in zip(ChEMBL_pathogen['activity_type'], ChEMBL_pathogen['unit'])]
    print(f"Directions assigned. Summary: {Counter(ChEMBL_pathogen['direction'].fillna('NaN'))}")

    # Save cleaned data
    ChEMBL_pathogen.to_csv(os.path.join(root, "..", "output", pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), index=False)

    # Get unique assays
    assays = sorted(set(ChEMBL_pathogen['assay_chembl_id']))

    ASSAYS_INFO = []
    print("Collecting individual assay information...")
    print(f"Number of unique assays: {len(assays)}")

    # For each assay
    for assay in tqdm(assays):

        # Get subset of strain + assay data
        df_ = ChEMBL_pathogen[ChEMBL_pathogen["assay_chembl_id"] == assay]
        
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
                if type(u) != str:
                    df___ = df__[df__["unit"].isna()]
                else:
                    df___ = df__[df__["unit"] == u]
                unit = list(set(df___['unit']))
                unit = only_one(unit, "unit")
                activities = len(df___)
                cpds = len(set(df___['compound_chembl_id']))
                nan_values = len(df___[df___['value'].isna()])
                direction = DIRECTIONS[(act_type, unit)] if (act_type, unit) in DIRECTIONS else np.nan
                canonical_unit = canonical_map[act_type]
                ASSAYS_INFO.append([assay, assay_type, assay_organism, doc_chembl_id, target_type, target_chembl_id, 
                                    target_organism, activity_type, unit, canonical_unit, activities, nan_values, cpds, direction])

    ASSAYS_INFO = pd.DataFrame(ASSAYS_INFO, columns=["assay_id", "assay_type", "assay_organism", "doc_chembl_id", "target_type", "target_chembl_id", "target_organism", 
                                                        "activity_type", "unit", "canonical_unit", "activities", 'nan_values', "cpds", "direction"])
    ASSAYS_INFO = ASSAYS_INFO.sort_values('cpds', ascending=False).reset_index(drop=True)

    # Filter assays with too few compounds
    ASSAYS_INFO = ASSAYS_INFO[ASSAYS_INFO['cpds'] > MIN_ASSAY_SIZE].reset_index(drop=True)

    # Save assays info
    # ASSAYS_INFO.to_csv(os.path.join(root, "..", "output", pathogen_code, 'assays_cleaned.csv'), index=False)