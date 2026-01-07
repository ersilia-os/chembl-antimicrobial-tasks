from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import json
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
        on the wrong side of the direction (min or max)

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
root = os.path.dirname(os.path.abspath(__file__))

# List of pathogens to process
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]
pathogens = ["Acinetobacter baumannii", "Mycobacterium tuberculosis", "Klebsiella pneumoniae"]

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

# For each pathogen
for pathogen in pathogens:

    print("\n\n\n")
    pathogen_code = get_pathogen_code(pathogen)

    # Load cleaned assays
    ASSAYS_CLEANED = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, "assays_cleaned.csv"))

    # Define PATH to parameters
    PATH_TO_PARAMETERS = os.path.join(root, "..", "output", pathogen_code, 'parameters')

    TARGET_TYPE_CURATED = []

    # Iterating over assays
    for assay_id, activity_type, unit in ASSAYS_CLEANED[['assay_id', 'activity_type', 'unit']].values:

        # Prepare filename
        filename = "_".join([str(assay_id), str(activity_type), str(unit), 'parameters']) + ".json"
        
        # Read JSON file
        with open(os.path.join(PATH_TO_PARAMETERS, filename), "r") as file:
            par = json.load(file)

        # Store results
        TARGET_TYPE_CURATED.append(par['target_type_curated'])

    # Complete table
    ASSAYS_CLEANED['target_type_curated'] = TARGET_TYPE_CURATED

    # Loading pathogen data
    os.makedirs(os.path.join(OUTPUT, pathogen_code, 'datasets'), exist_ok=True)
    print(f"Loading ChEMBL cleaned data for {pathogen_code}...")
    ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
    print(f"Number of cleaned assays: {len(ASSAYS_CLEANED)}")

    # Load expert cut-offs
    EXPERT_CUTOFFS = pd.read_csv(os.path.join(root, "..", "config", "manual_curation", "expert_cutoffs.csv"))
    EXPERT_CUTOFFS = {
        (a, b, c, d): e
        for a, b, c, d, e in EXPERT_CUTOFFS[
            ["activity_type", "unit", "target_type", "pathogen_code", "expert_cutoff"]
        ].values}

    # Define data ranges
    DATA_RANGES = []

    for assay_chembl_id, activity_type, unit, target_type, target_type_curated, activities, cpds, direction in tqdm(ASSAYS_CLEANED[['assay_id', 'activity_type', 'unit', 'target_type',
                                                                                                            'target_type_curated', 'activities', 'cpds', 'direction']].values):

        # Loading assay data
        cols = ['compound_chembl_id', 'canonical_smiles', 'activity_type', 'value', 'relation', 'unit', 'activity_comment', 'standard_text']
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
            CUT = np.nan

        # Get expert cut-off if it exists
        key = (activity_type, unit, target_type_curated, pathogen_code)
        expert_cutoff = EXPERT_CUTOFFS[key] if key in EXPERT_CUTOFFS else np.nan

        if direction in [+1, -1]:

            # Adjust relation
            ASSAY_DATA = adjust_relation(ASSAY_DATA, direction, CUT)

            # Disambiguate duplicated compounds and returns 'sorted' data (depending on direction)
            # Numerical values are prioritized over nans
            ASSAY_DATA = disambiguate_compounds(ASSAY_DATA, direction)

            # Remove nans
            assay_activities = [i for i in ASSAY_DATA['value'].tolist() if np.isnan(i) == False]

            # Fully qualitative assay
            if len(assay_activities) == 0:
                assay_activities = [np.nan]
                dataset_type = 'qualitative'
            elif len(assay_activities) < len(ASSAY_DATA):
                dataset_type = 'mixed'
            else:
                dataset_type = 'quantitative'

            # Binarization with expert cut-off
            if np.isnan(expert_cutoff) == False:
                if direction == +1:
                    ASSAY_DATA["bin_quantitative"] = (ASSAY_DATA["value"] >= expert_cutoff).astype(int)
                else:
                    ASSAY_DATA["bin_quantitative"] = (ASSAY_DATA["value"] <= expert_cutoff).astype(int)
                positives = Counter(ASSAY_DATA['bin_quantitative'].tolist())[1]
                ratio = round(positives / len(ASSAY_DATA), 5)
            else:
                # Dataset could not be binarized using values due to missing expert cut-off
                ASSAY_DATA['bin_quantitative'] = [np.nan] * len(ASSAY_DATA)
                positives = np.nan
                ratio = np.nan

        else:

            # Qualitative assay
            dataset_type = 'qualitative'    

            # Not binarizing a null-direction assay
            ASSAY_DATA['bin_quantitative'] = [np.nan] * len(ASSAY_DATA)
            positives = np.nan
            ratio = np.nan
            assay_activities = [np.nan]


        # Getting qualitative bin
        ASSAY_DATA['bin_qualitative'] = [np.nan for i in range(len(ASSAY_DATA))]
        cond_nan = (ASSAY_DATA['activity_comment'] == 0) & (ASSAY_DATA['standard_text'] == 0)
        cond_pos = (ASSAY_DATA['activity_comment'] == 1) | (ASSAY_DATA['standard_text'] == 1)
        cond_neg = (ASSAY_DATA['activity_comment'] == -1) | (ASSAY_DATA['standard_text'] == -1)
        conflict = cond_pos & cond_neg
        if conflict.any():
            raise ValueError(
                "Conflicting labels (contains both 1 and -1):\n"
                + ASSAY_DATA.loc[conflict, ["activity_comment", "standard_text"]].head(20).to_string())
        ASSAY_DATA.loc[cond_pos, "bin_qualitative"] = 1
        ASSAY_DATA.loc[cond_neg, "bin_qualitative"] = 0

        # Positives and negatives qualitative
        positives_qual = Counter(ASSAY_DATA['bin_qualitative'].tolist())[1]
        negatives_qual = Counter(ASSAY_DATA['bin_qualitative'].tolist())[0]

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
        DATA_RANGES.append([assay_chembl_id, activity_type, unit, target_type, target_type_curated, activities, cpds, direction, equal, higher, 
                            lower, min_, p1, p25, p50, p75, p99, max_, expert_cutoff, positives, ratio, dataset_type, positives_qual, negatives_qual])

        # Save data
        if cpds >= 100:
            ASSAY_DATA.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{assay_chembl_id}_{activity_type}_{str(unit).replace('/', 'FwdS')}.csv.gz"), index=False)

    DATA_RANGES = pd.DataFrame(DATA_RANGES, columns=["assay_id", "activity_type", "unit", "target_type", "target_type_curated", "activities", "cpds", "direction", 
                                                    "equal", "higher", "lower", "min_", "p1", "p25", "p50", "p75", "p99", "max_", 'expert_cutoff', 'positives_quant', 'ratio_quant', 
                                                    'dataset_type', 'positives_qual', 'negatives_qual'])
    DATA_RANGES.to_csv(os.path.join(OUTPUT, pathogen_code, 'assays_data_ranges.csv'), index=False)
