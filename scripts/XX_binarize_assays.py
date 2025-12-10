import pandas as pd
import numpy as np
import sys
import os
# Define root directory
# root = os.path.dirname(os.path.abspath(__file__))
root = "."
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH

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

    # Loading pathogen data
    pathogen_code = get_pathogen_code(pathogen)
    print(f"Loading ChEMBL preprocessed data for {pathogen_code}...")
    ChEMBL = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, f"{pathogen_code}_ChEMBL_data.csv"), low_memory=False)
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL['compound_chembl_id']))}")
    ASSAYS = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, 'assays.csv'))
    print(f"Original number of assays: {len(ASSAYS)}")

    # Create output directory
    os.makedirs(os.path.join(OUTPUT, pathogen_code, "datasets"), exist_ok=True)

    break

# Get only assays with 100 or more compounds
ASSAYS = ASSAYS[ASSAYS['cpds'] >= 100].reset_index(drop=True)

# Get directions
DIRECTIONS = pd.read_csv(os.path.join(root, "..", "config", 'manual_curation', 'activity_std_units_curated_manual_curation.csv'))
DIRECTIONS = {(i,j): k for i,j,k in zip(DIRECTIONS['activity_type'], DIRECTIONS['unit'], DIRECTIONS['manual_curation']) if np.isnan(k) == False}

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

# Define some parameters
cols = ['compound_chembl_id', 'canonical_smiles', 'activity_type', 'value', 'relation', 'unit']
BINARIZATION_RESULTS = []
PERC = [5, 10]

for assay_chembl_id, activity_type, unit in zip(ASSAYS["assay_id"], ASSAYS["activity_type"], ASSAYS["unit"]):

    # Get direction
    if (activity_type, unit) in DIRECTIONS:
        direction = DIRECTIONS[(activity_type, unit)]
    else:
        # raise ValueError(f"Direction not found for ({activity_type}, {unit}). Consider adding the direction manually in 'activity_std_units_curated_manual_curation.csv'")
        print(f"Direction not found for ({activity_type}, {unit}). Consider adding the direction manually in 'activity_std_units_curated_manual_curation.csv'")

    if direction not in [+1, -1]:
        print(f"Direction is not +1/-1: ({activity_type}, {unit}). Omitting assay...")
        break

    # Loading assay data
    ASSAY_DATA = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_chembl_id) & (ChEMBL['activity_type'] == activity_type) & 
                        (ChEMBL['unit'] == unit)].reset_index(drop=True)[cols]


    # Get value to adjust relations
    if direction == 1:
        CUT = min(ASSAY_DATA['value']) - 1
    elif direction == -1:
        CUT = max(ASSAY_DATA['value']) + 1
    else:
        pass

    # Adjust relation
    ASSAY_DATA = adjust_relation(ASSAY_DATA, direction, CUT)

    # Disambiguate duplicated compounds
    ASSAY_DATA = disambiguate_compounds(ASSAY_DATA, direction)

    # For each percentile
    for perc in PERC:

        # Binarize
        N = int(len(ASSAY_DATA) * perc / 100)
        cut_off = ASSAY_DATA['value'][N]
        ASSAY_DATA['bin'] = [1] * N + [0] * (len(ASSAY_DATA) - N)
        actives, inactives = len(ASSAY_DATA[ASSAY_DATA['bin'] == 1]), len(ASSAY_DATA[ASSAY_DATA['bin'] == 0])

        # Save data
        ASSAY_DATA.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{assay_chembl_id}_{activity_type}_{unit}_{perc}.csv.gz"), index=False)

        # Save more data
        BINARIZATION_RESULTS.append([assay_chembl_id, activity_type, unit, perc, cut_off, len(ASSAY_DATA), actives, inactives, round(actives/len(ASSAY_DATA), 3), direction])

BINARIZATION_RESULTS = pd.DataFrame(BINARIZATION_RESULTS, columns=["assay_chembl_id", "activity_type", "unit", "perc", 
                                                                   "cut_off", "compounds", "actives", "inactives", "ratio", "direction"])
BINARIZATION_RESULTS.to_csv(os.path.join(OUTPUT, pathogen_code, "datasets_binarization.csv"), index=False)