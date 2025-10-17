from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH, MIN_ASSAY_SIZE

# Loading ChEMBL preprocessed data
print("Loading ChEMBL preprocessed data...")
ChEMBL = pd.read_csv(os.path.join(root, "..", "config", "chembl_processed", "activities_preprocessed.csv"), low_memory=False)
print(f"Original size: {len(ChEMBL)}")
print("Filtering out nan values...")
ChEMBL = ChEMBL[ChEMBL['value'].isna() == False].reset_index(drop=True)
print(f"Size after filtering nan values: {len(ChEMBL)}")

# Get pathogen data
pathogens = ["Mycobacterium tuberculosis"]

for pathogen in pathogens:

    print(f"Filtering for pathogen: {pathogen}...")
    pathogen_code = str(pathogen.split()[0][0] + pathogen.split()[1]).lower()
    PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code)
    os.makedirs(PATH_TO_OUTPUT, exist_ok=True)
    ChEMBL = ChEMBL[ChEMBL['target_organism'].str.contains(pathogen, case=False, na=False) | 
                    ChEMBL['assay_organism'].str.contains(pathogen, case=False, na=False)].reset_index(drop=True)

    print(f"Number of activities: {len(ChEMBL)}")
    print(f"Number of unique compounds: {len(set(ChEMBL['compound_chembl_id']))}") 
    df = dict(Counter(ChEMBL['target_organism']))
    df = pd.DataFrame([[i, df[i]] for i in sorted(df, key = lambda x: df[x], reverse=True)], columns=['organism', 'count'])
    df.to_csv(os.path.join(PATH_TO_OUTPUT, "target_organism_counts.csv"), index=False)

    # Helper function - is there only a single value?
    def only_one(values, name):
        if len(values) != 1:
            raise ValueError(f"Expected exactly one {name}, found {values}")
        return values[0]

    # Get unique assays
    assays = sorted(set(ChEMBL['assay_chembl_id']))

    ASSAYS_INFO = []
    print("Collecting individual assay information...")

    # For each assay
    for assay in tqdm(assays):

        # Get subset of strain + assay data
        df_ = ChEMBL[ChEMBL["assay_chembl_id"] == assay]
        
        # Get values
        assay_type = list(set(df_['assay_type']))
        target_type = list(set(df_['target_type']))
        target_chembl_id = list(set(df_['target_chembl_id']))
        activity_types = list(set(df_['activity_type']))
        target_organism = list(set(df_['target_organism']))
        assay_organism = list(set(df_['assay_organism']))

        # Check coherence
        assay_type = only_one(assay_type, "assay_type")
        target_type = only_one(target_type, "target_type")
        target_chembl_id = only_one(target_chembl_id, "target_chembl_id")
        target_organism = only_one(target_organism, "target_organism")
        assay_organism = only_one(assay_organism, "assay_organism")

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
                ASSAYS_INFO.append([assay, assay_type, assay_organism, target_type, target_chembl_id, target_organism, activity_type, unit, activities, cpds])

    ASSAYS_INFO = pd.DataFrame(ASSAYS_INFO, columns=["assay_id", "assay_type", "assay_organism", "target_type", "target_chembl_id", "target_organism", "activity_type", "unit", "activities", "cpds"])
    ASSAYS_INFO = ASSAYS_INFO.sort_values('cpds', ascending=False).reset_index(drop=True)

    # Filter assays with too few compounds
    ASSAYS_INFO = ASSAYS_INFO[ASSAYS_INFO['cpds'] > MIN_ASSAY_SIZE].reset_index(drop=True)

    # Save assays info
    ASSAYS_INFO.to_csv(os.path.join(PATH_TO_OUTPUT, 'assays.csv'), index=False)