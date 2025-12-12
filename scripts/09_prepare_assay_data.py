from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import sys
import os

def adjust_relation(ASSAY_DATA: pd.DataFrame, DIRECTION: int, CUT: float) -> pd.DataFrame:
    """
    Adjust relations in an assay DataFrame according to the biological direction.

    Parameters
    ----------
    ASSAY_DATA : pd.DataFrame
        Must contain the columns 'relation' and 'value'.
    DIRECTION : int
        +1 → higher = more active (e.g. % inhibition)
        -1 → lower = more active (e.g. IC50, MIC)
    CUT : float
        Extreme value used to replace censored measurements
        on the wrong side of the direction (min-1 or max+1)

    Returns
    -------
    pd.DataFrame
        Copy of ASSAY_DATA with adjusted relation and value.
    """

    df = ASSAY_DATA.copy()
    rel = df["relation"].astype(str)

    if DIRECTION == +1:

        # Higher = more active
        mask_gt = rel == ">"  # greater than
        mask_lt = rel == "<"  # lower than

        df.loc[mask_gt, "relation"] = "="
        df.loc[mask_lt, "relation"] = "="
        df.loc[mask_lt, "value"] = CUT

    elif DIRECTION == -1:

        # Lower = more active
        mask_lt = rel == "<"  # lower than
        mask_gt = rel == ">"  # greater than

        df.loc[mask_lt, "relation"] = "="
        df.loc[mask_gt, "relation"] = "="
        df.loc[mask_gt, "value"] = CUT

    else:

        raise ValueError(f"Invalid DIRECTION={DIRECTION}. Expected +1 or -1.")

    return df

def disambiguate_compounds(ASSAY_DATA: pd.DataFrame, DIRECTION: int) -> pd.DataFrame:

    """
    Select a single measurement per compound according to the biological direction.

    Parameters
    ----------
    ASSAY_DATA : pd.DataFrame
        Must contain the columns 'compound_chembl_id' and 'value'.
        Assumes all relations have already been adjusted.
    DIRECTION : int
        +1 → higher = more active (e.g. % inhibition)
        -1 → lower = more active (e.g. IC50, MIC)

    Returns
    -------
    pd.DataFrame
        A copy of ASSAY_DATA in which duplicated compounds 
        ('compound_chembl_id') are removed, keeping only the 
        most active measurement per compound (highest or lowest 
        depending on DIRECTION).
    """

    if DIRECTION not in [1, -1]:
        raise ValueError("DIRECTION must be +1 (higher = more active) or -1 (lower = more active).")
        
    df = ASSAY_DATA.copy()

    # Choose best measurement based on direction
    if DIRECTION == -1:
        # Lower = more active → keep minimum
        df_sorted = df.sort_values(by="value", ascending=True)
    elif DIRECTION == 1:
        # Higher = more active → keep maximum
        df_sorted = df.sort_values(by="value", ascending=False)

    # Keep the best row per compound_chembl_id
    df_best = df_sorted.drop_duplicates(subset="compound_chembl_id", keep="first")

    return df_best.reset_index(drop=True)

# Define root directory
# root = os.path.dirname(os.path.abspath(__file__))
root = "."

# List of pathogens to process
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"][8:9]

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

# For each pathogen
for pathogen in pathogens:

    print("\n\n\n")

    # Loading pathogen data
    pathogen_code = get_pathogen_code(pathogen)
    os.makedirs(os.path.join(OUTPUT, pathogen_code, 'datasets'), exist_ok=True)
    print(f"Loading ChEMBL cleaned data for {pathogen_code}...")
    ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
    ASSAYS_CLEANED = pd.read_csv(os.path.join(OUTPUT, pathogen_code, 'assays_cleaned.csv'))
    print(f"Cleaned number of assays: {len(ASSAYS_CLEANED)}")

    # Load expert cut-offs
    EXPERT_CUTOFFS = pd.read_csv(os.path.join(root, "..", "config", "manual_curation", "expert_cutoffs.csv"))
    EXPERT_CUTOFFS = {
        (a, b, c, d): e
        for a, b, c, d, e in EXPERT_CUTOFFS[
            ["activity_type", "unit", "target_type", "pathogen_code", "expert_cutoff"]
        ].values}

    # Define data ranges
    DATA_RANGES = []

    for assay_chembl_id, activity_type, unit, target_type, activities, cpds, direction in tqdm(ASSAYS_CLEANED[['assay_id', 'activity_type', 'unit',
                                                                                                            'target_type', 'activities', 'cpds', 'direction']].values[2:]):
            
        if direction not in [+1, -1]:
            if direction == 0:
                print(f"Direction is 0: ({activity_type}, {unit}). Omitting assay {assay_chembl_id}...")
            else:
                print(f"Direction not found for ({activity_type}, {unit}). Consider adding the direction manually in 'activity_std_units_curated_manual_curation.csv'")

        # Loading assay data
        cols = ['compound_chembl_id', 'canonical_smiles', 'activity_type', 'value', 'relation', 'unit']
        if type(unit) == str:
            ASSAY_DATA = ChEMBL_pathogen[(ChEMBL_pathogen['assay_chembl_id'] == assay_chembl_id) & 
                                        (ChEMBL_pathogen['activity_type'] == activity_type) & 
                                        (ChEMBL_pathogen['unit'] == unit)].reset_index(drop=True)[cols]
        else:
            ASSAY_DATA = ChEMBL_pathogen[(ChEMBL_pathogen['assay_chembl_id'] == assay_chembl_id) & 
                                        (ChEMBL_pathogen['activity_type'] == activity_type) & 
                                        (ChEMBL_pathogen['unit'].isna())].reset_index(drop=True)[cols]
        
        # Counter relations
        counter_relations = Counter(ASSAY_DATA['relation'].tolist())

        # Get value to adjust relations
        if direction == 1:
            CUT = min(ASSAY_DATA['value'])
        elif direction == -1:
            CUT = max(ASSAY_DATA['value'])
        else:
            CUT = None

        # Get expert cut-off if exists
        key = (activity_type, unit, target_type, pathogen_code)
        expert_cutoff = EXPERT_CUTOFFS[key] if key in EXPERT_CUTOFFS else np.nan

        if direction in [+1, -1]:

            # Adjust relation
            ASSAY_DATA = adjust_relation(ASSAY_DATA, direction, CUT)

            # Disambiguate duplicated compounds and returns 'sorted' data (depending on direction)
            ASSAY_DATA = disambiguate_compounds(ASSAY_DATA, direction)

            # Remove nans
            assay_activities = [i for i in ASSAY_DATA['value'].tolist() if np.isnan(i) == False]
            if len(assay_activities) == 0:
                assay_activities = [np.nan]

            # Binarization with expert cut-off
            if np.isnan(expert_cutoff) == False:
                if direction == +1:
                    ASSAY_DATA["bin"] = (ASSAY_DATA["value"] >= expert_cutoff).astype(int)
                else:
                    ASSAY_DATA["bin"] = (ASSAY_DATA["value"] <= expert_cutoff).astype(int)
                positives = Counter(ASSAY_DATA['bin'].tolist())[1]
                ratio = round(positives / len(ASSAY_DATA), 5)
            else:
                ASSAY_DATA['bin'] = [np.nan] * len(ASSAY_DATA)
                positives = np.nan
                ratio = np.nan

        else:

            # Take only equal relations
            assay_activities = [i for i,j in ASSAY_DATA[['value', 'relation']].values if np.isnan(i) == False and j == "="]
            if len(assay_activities) == 0:
                assay_activities = [np.nan]

            # Binarization with expert cut-off
            ASSAY_DATA['bin'] = [np.nan] * len(ASSAY_DATA)
            positives = np.nan
            ratio = np.nan
        
        # Calculate data
        min_ = round(np.min(assay_activities), 3)
        p1 = round(np.percentile(assay_activities, 1), 3)
        p25 = round(np.percentile(assay_activities, 25), 3)
        p50 = round(np.percentile(assay_activities, 50), 3)
        p75 = round(np.percentile(assay_activities, 75), 3)
        p99 = round(np.percentile(assay_activities, 99), 3)
        max_ = round(np.max(assay_activities), 3)

        # Relations
        equal = counter_relations["="]
        lower = counter_relations["<"]
        higher = counter_relations[">"]

        # Store data range
        DATA_RANGES.append([assay_chembl_id, activity_type, unit, target_type, activities, cpds, direction, equal, higher, 
                            lower, min_, p1, p25, p50, p75, p99, max_, expert_cutoff, positives, ratio])

        # Save data
        ASSAY_DATA.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{assay_chembl_id}_{activity_type}_{str(unit).replace('/', 'FwdS')}.csv.gz"), index=False)

DATA_RANGES = pd.DataFrame(DATA_RANGES, columns=["assay_chembl_id", "activity_type", "unit", "target_type", "activities", "cpds", "direction", "equal", "higher", 
                                            "lower", "min_", "p1", "p25", "p50", "p75", "p99", "max_", 'expert_cutoff', 'positives', 'ratio'])
DATA_RANGES.to_csv(os.path.join(OUTPUT, pathogen_code, 'assay_data_ranges.csv'), index=False)