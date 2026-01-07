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

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

def add_target_type_curated(ASSAYS_CLEANED, PATH_TO_PARAMETERS):
    """
    Add a `target_type_curated` column to ASSAYS_CLEANED by reading parameter JSON files.

    For each row in ASSAYS_CLEANED, a JSON file is opened using the pattern:
        "{assay_id}_{activity_type}_{unit}_parameters.json"

    The value stored under the key `"target_type_curated"` is extracted and appended
    to a list, which is then assigned as a new column in the dataframe.

    Parameters
    ----------
    ASSAYS_CLEANED : pandas.DataFrame
        DataFrame containing at least the columns: `assay_id`, `activity_type`, `unit`.
    PATH_TO_PARAMETERS : str
        Path to the directory containing the JSON parameter files.

    Returns
    -------
    pandas.DataFrame
        The same dataframe with an added `target_type_curated` column.
    """
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

    return ASSAYS_CLEANED

def load_expert_cutoffs(root):
    """
    Load expert cutoffs from the manual curation CSV and return them as a dictionary.

    The CSV is expected at:
        {root}/../config/manual_curation/expert_cutoffs.csv

    The returned dictionary maps:
        (activity_type, unit, target_type, pathogen_code) -> expert_cutoff

    Parameters
    ----------
    root : str
        Base path used to locate the config folder.

    Returns
    -------
    dict
        Dictionary of expert cutoffs keyed by
        (activity_type, unit, target_type, pathogen_code).
    """
    # Load expert cut-offs
    EXPERT_CUTOFFS = pd.read_csv(
        os.path.join(root, "..", "config", "manual_curation", "expert_cutoffs.csv")
    )

    EXPERT_CUTOFFS = {
        (a, b, c, d): e
        for a, b, c, d, e in EXPERT_CUTOFFS[
            ["activity_type", "unit", "target_type", "pathogen_code", "expert_cutoff"]
        ].values
    }

    return EXPERT_CUTOFFS

def get_assay_data(ChEMBL_pathogen, assay_chembl_id, activity_type, unit, cols):
    """
    Extract assay activity data for a given assay_chembl_id, activity_type, and unit.

    If `unit` is a string, the function filters rows where `unit` matches exactly.
    Otherwise, it filters rows where `unit` is missing (NaN).

    Parameters
    ----------
    ChEMBL_pathogen : pandas.DataFrame
        DataFrame containing ChEMBL pathogen activity records.
    assay_chembl_id : str
        Assay ChEMBL ID to filter on.
    activity_type : str
        Activity type to filter on (e.g., IC50, MIC).
    unit : str or None
        Unit to filter on; if not a string, NaN units are selected.
    cols : list
        List of columns to return.

    Returns
    -------
    pandas.DataFrame
        Filtered assay activity data with only the requested columns.
    """
    if type(unit) == str:
        ASSAY_DATA = ChEMBL_pathogen[
            (ChEMBL_pathogen['assay_chembl_id'] == assay_chembl_id) &
            (ChEMBL_pathogen['activity_type'] == activity_type) &
            (ChEMBL_pathogen['unit'] == unit)
        ].reset_index(drop=True)[cols]
    else:
        ASSAY_DATA = ChEMBL_pathogen[
            (ChEMBL_pathogen['assay_chembl_id'] == assay_chembl_id) &
            (ChEMBL_pathogen['activity_type'] == activity_type) &
            (ChEMBL_pathogen['unit'].isna())
        ].reset_index(drop=True)[cols]

    return ASSAY_DATA

def get_cut_value(ASSAY_DATA, direction):
    """
    Get a cutoff value from ASSAY_DATA to adjust relations based on direction.

    If direction == 1, returns the minimum value in ASSAY_DATA['value'].
    If direction == -1, returns the maximum value in ASSAY_DATA['value'].
    Otherwise, returns np.nan.

    Parameters
    ----------
    ASSAY_DATA : pandas.DataFrame
        DataFrame containing a 'value' column with numeric assay values.
    direction : int
        Direction indicator:
        - 1  -> use minimum value
        - -1 -> use maximum value
        - else -> np.nan

    Returns
    -------
    float
        Cutoff value computed from the 'value' column or np.nan.
    """
    if direction == 1:
        CUT = min(ASSAY_DATA['value'])
    elif direction == -1:
        CUT = max(ASSAY_DATA['value'])
    else:
        CUT = np.nan

    return CUT

def count_relations(ASSAY_DATA):
    """
    Count relation operators in ASSAY_DATA['relation'].

    Counts how many times each of the following appears:
        "=" , "<" , ">"

    Parameters
    ----------
    ASSAY_DATA : pandas.DataFrame
        DataFrame containing a 'relation' column.

    Returns
    -------
    tuple
        (equal, lower, higher) counts corresponding to "=", "<", ">".
    """
    counter_relations = Counter(ASSAY_DATA['relation'].tolist())
    equal = counter_relations["="]
    lower = counter_relations["<"]
    higher = counter_relations[">"]

    return equal, lower, higher

def get_assay_data_quantitative(ASSAY_DATA):
    """
    Return only rows in ASSAY_DATA with non-missing quantitative values.

    Filters ASSAY_DATA to keep rows where `value` is not NaN, and resets the index.

    Parameters
    ----------
    ASSAY_DATA : pandas.DataFrame
        DataFrame containing a 'value' column.

    Returns
    -------
    pandas.DataFrame
        Filtered dataframe containing only rows with non-null `value`.
    """
    ASSAY_DATA_QUANTITATIVE = ASSAY_DATA[ASSAY_DATA['value'].isna() == False].reset_index(drop=True)
    return ASSAY_DATA_QUANTITATIVE

def get_assay_data_qualitative(ASSAY_DATA):
    """
    Generate qualitative labels for assay data based on `activity_comment` and `standard_text`.

    This function:
    - Detects row-level conflicts where a row is simultaneously positive (1) and negative (-1).
    - Assigns a row-level qualitative label (1, -1, or 0).
    - Aggregates labels at the compound level to ensure each compound has a single consistent label.
    - Detects compound-level conflicts where the same compound is labeled both 1 and -1 across rows.
    - Assigns a final compound-level label to all rows.
    - Keeps only one row per compound.
    - Removes compounds with final label 0.

    Parameters
    ----------
    ASSAY_DATA : pandas.DataFrame
        DataFrame containing at least:
        - compound_chembl_id
        - activity_comment
        - standard_text

    Returns
    -------
    pandas.DataFrame
        Filtered dataframe with one row per compound and `qualitative_label` in {1, -1}.
    """
    # Qualitative view
    ASSAY_DATA_QUALITATIVE = ASSAY_DATA.copy()
    cond_nan = (ASSAY_DATA_QUALITATIVE['activity_comment'] == 0) & (ASSAY_DATA_QUALITATIVE['standard_text'] == 0)
    cond_pos = (ASSAY_DATA_QUALITATIVE['activity_comment'] == 1) | (ASSAY_DATA_QUALITATIVE['standard_text'] == 1)
    cond_neg = (ASSAY_DATA_QUALITATIVE['activity_comment'] == -1) | (ASSAY_DATA_QUALITATIVE['standard_text'] == -1)

    # Detect row-level conflicts
    conflict = cond_pos & cond_neg
    if conflict.any():
        raise ValueError(
            "Conflicting labels (contains both 1 and -1):\n"
            + ASSAY_DATA_QUALITATIVE.loc[conflict, ["compound_chembl_id", "activity_comment", "standard_text"]].head(20).to_string())

    # Assign row-level label
    ASSAY_DATA_QUALITATIVE["qualitative_label_row"] = np.nan
    ASSAY_DATA_QUALITATIVE.loc[cond_pos, "qualitative_label_row"] = 1
    ASSAY_DATA_QUALITATIVE.loc[cond_neg, "qualitative_label_row"] = -1
    ASSAY_DATA_QUALITATIVE.loc[cond_nan, "qualitative_label_row"] = 0

    # Aggregate to compound-level label
    compound_labels = ASSAY_DATA_QUALITATIVE.groupby("compound_chembl_id")["qualitative_label_row"].apply(set)

    # Detect compound-level conflicts (same compound has 1 and -1)
    compound_ids = compound_labels.index.tolist()
    label_sets = compound_labels.tolist()
    compound_conflict = [((1 in s) and (-1 in s)) for s in label_sets]
    if any(compound_conflict):
        bad = [(cid, s) for cid, s, c in zip(compound_ids, label_sets, compound_conflict) if c][:20]
        raise ValueError(
            "Conflicting compound labels (same compound has both 1 and -1 across rows):\n"
            + "\n".join([f"{cid}: {s}" for cid, s in bad]))

    # Final compound label: 1 > -1 > 0
    compound_final = [1 if 1 in s else (-1 if -1 in s else 0) for s in label_sets]
    compound_final = dict(zip(compound_ids, compound_final))

    # Assign back to all rows
    ASSAY_DATA_QUALITATIVE["qualitative_label"] = ASSAY_DATA_QUALITATIVE["compound_chembl_id"].map(compound_final)

    # Keep only one row per compound
    ASSAY_DATA_QUALITATIVE = ASSAY_DATA_QUALITATIVE.drop_duplicates(subset=["compound_chembl_id"]).reset_index(drop=True)

    # Remove compounds labeled as 0
    ASSAY_DATA_QUALITATIVE = ASSAY_DATA_QUALITATIVE[ASSAY_DATA_QUALITATIVE["qualitative_label"] != 0].reset_index(drop=True)

    # Binary label
    ASSAY_DATA_QUALITATIVE["bin"] = [0 if x == -1 else 1 for x in ASSAY_DATA_QUALITATIVE["qualitative_label"].tolist()]

    # Take only interesting columns
    cols = ["compound_chembl_id", "canonical_smiles", "activity_type", "unit", "activity_comment", "standard_text", "qualitative_label", 'bin'] 
    ASSAY_DATA_QUALITATIVE = ASSAY_DATA_QUALITATIVE[cols]

    return ASSAY_DATA_QUALITATIVE

def set_variables_quantitative(ASSAY_DATA_QUANTITATIVE):
    """
    Compute summary variables for quantitative assay data.

    Returns number of positives, positive ratio, and number of compounds (rows).
    If the dataframe is empty, returns np.nan for all values.

    Parameters
    ----------
    ASSAY_DATA_QUANTITATIVE : pandas.DataFrame
        DataFrame containing a 'bin' column (0/1).

    Returns
    -------
    tuple
        (positives_quantitative, ratio_quantitative, compounds_quantitative)
    """
    if len(ASSAY_DATA_QUANTITATIVE) > 0 and ASSAY_DATA_QUANTITATIVE['bin'].isna().any() == False:
            positives_quantitative = (ASSAY_DATA_QUANTITATIVE["bin"] == 1).sum()
            ratio_quantitative = round(positives_quantitative / len(ASSAY_DATA_QUANTITATIVE), 3)
            compounds_quantitative = len(ASSAY_DATA_QUANTITATIVE)
            activities_quantitative = ASSAY_DATA_QUANTITATIVE['value'].tolist()
    else:
        positives_quantitative = np.nan
        ratio_quantitative = np.nan
        compounds_quantitative = np.nan
        activities_quantitative = [np.nan]

    return positives_quantitative, ratio_quantitative, compounds_quantitative, activities_quantitative

def set_variables_qualitative(ASSAY_DATA_QUALITATIVE):
    """
    Compute summary variables for qualitative assay data.

    Returns number of positives, positive ratio, and number of compounds (rows).
    If the dataframe is empty, returns np.nan for all values.

    Parameters
    ----------
    ASSAY_DATA_QUALITATIVE : pandas.DataFrame
        DataFrame containing a 'bin' column (0/1).

    Returns
    -------
    tuple
        (positives_qualitative, ratio_qualitative, compounds_qualitative)
    """
    if len(ASSAY_DATA_QUALITATIVE) > 0:
        positives_qualitative = (ASSAY_DATA_QUALITATIVE["bin"] == 1).sum()
        ratio_qualitative = round(positives_qualitative / len(ASSAY_DATA_QUALITATIVE), 3)
        compounds_qualitative = len(ASSAY_DATA_QUALITATIVE)
    else:
        positives_qualitative = np.nan
        ratio_qualitative = np.nan
        compounds_qualitative = np.nan

    return positives_qualitative, ratio_qualitative, compounds_qualitative

def binarize_with_expert_cutoff(ASSAY_DATA_QUANTITATIVE, expert_cutoff, direction):
    """
    Binarize quantitative assay values using an expert cutoff.

    If `expert_cutoff` is not NaN:
      - direction == +1 -> bin = 1 if value >= expert_cutoff else 0
      - direction != +1 -> bin = 1 if value <= expert_cutoff else 0

    If `expert_cutoff` is NaN:
      - bin is set to NaN for all rows.

    Parameters
    ----------
    ASSAY_DATA_QUANTITATIVE : pandas.DataFrame
        DataFrame containing a numeric 'value' column.
    expert_cutoff : float
        Expert-defined cutoff value (may be NaN).
    direction : int
        Direction of binarization (+1 or -1).

    Returns
    -------
    pandas.DataFrame
        The same dataframe with a new column 'bin'.
    """
    if np.isnan(expert_cutoff) == False:
        if direction == +1:
            ASSAY_DATA_QUANTITATIVE["bin"] = (ASSAY_DATA_QUANTITATIVE["value"] >= expert_cutoff).astype(int)
        else:
            ASSAY_DATA_QUANTITATIVE["bin"] = (ASSAY_DATA_QUANTITATIVE["value"] <= expert_cutoff).astype(int)
    else:
        ASSAY_DATA_QUANTITATIVE["bin"] = [np.nan] * len(ASSAY_DATA_QUANTITATIVE)

    return ASSAY_DATA_QUANTITATIVE

def get_activity_stats_quantitative(activities_quantitative):
    """
    Compute summary statistics for quantitative activities.

    Calculates min, 1st percentile, 25th percentile, median (50th),
    75th percentile, 99th percentile, and max, rounded to 3 decimals.

    Parameters
    ----------
    activities_quantitative : array-like
        Iterable of numeric activity values.

    Returns
    -------
    tuple
        (min_, p1, p25, p50, p75, p99, max_)
    """
    min_ = round(np.min(activities_quantitative), 3)
    p1 = round(np.percentile(activities_quantitative, 1), 3)
    p25 = round(np.percentile(activities_quantitative, 25), 3)
    p50 = round(np.percentile(activities_quantitative, 50), 3)
    p75 = round(np.percentile(activities_quantitative, 75), 3)
    p99 = round(np.percentile(activities_quantitative, 99), 3)
    max_ = round(np.max(activities_quantitative), 3)

    return min_, p1, p25, p50, p75, p99, max_

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))

# List of pathogens to process
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]
pathogens = ["Acinetobacter baumannii", "Mycobacterium tuberculosis", "Klebsiella pneumoniae"]

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

# For each pathogen
for pathogen in pathogens[1:2]:

    # Get pathogen code
    pathogen_code = get_pathogen_code(pathogen)

    # Load cleaned assays
    ASSAYS_CLEANED = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, "assays_cleaned.csv"))

    # Define PATH to parameters
    PATH_TO_PARAMETERS = os.path.join(root, "..", "output", pathogen_code, 'parameters')
    # PATH_TO_PARAMETERS = os.path.join(root, "..", "output", pathogen_code, 'assay_parameters')

    # Get curated target type
    ASSAYS_CLEANED = add_target_type_curated(ASSAYS_CLEANED, PATH_TO_PARAMETERS)

    # Loading pathogen data
    os.makedirs(os.path.join(OUTPUT, pathogen_code, 'datasets'), exist_ok=True)
    print(f"Loading ChEMBL cleaned data for {pathogen_code}...")
    ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
    print(f"Number of cleaned assays: {len(ASSAYS_CLEANED)}")

    # Load expert cut-offs
    EXPERT_CUTOFFS = load_expert_cutoffs(root)

    # Define data ranges
    DATA_RANGES = []

    for assay_chembl_id, activity_type, unit, target_type, target_type_curated, activities, nan_values, cpds, direction, acc, stdtc in tqdm(ASSAYS_CLEANED[['assay_id', 
                                        'activity_type', 'unit', 'target_type','target_type_curated', 'activities', 'nan_values', 'cpds', 'direction', 
                                        'activity_comment_counts', 'standard_text_count']].values[:]):

        # Filtering [assay, activity_type, unit] data
        cols = ['compound_chembl_id', 'canonical_smiles', 'activity_type', 'value', 'relation', 'unit', 'activity_comment', 'standard_text']
        ASSAY_DATA = get_assay_data(ChEMBL_pathogen, assay_chembl_id, activity_type, unit, cols)
        
        # Count relations
        equal, lower, higher = count_relations(ASSAY_DATA)

        # Get expert cut-off if it exists
        key = (activity_type, unit, target_type_curated, pathogen_code)
        expert_cutoff = EXPERT_CUTOFFS[key] if key in EXPERT_CUTOFFS else np.nan

        # Quantitative view
        ASSAY_DATA_QUANTITATIVE = get_assay_data_quantitative(ASSAY_DATA)

        # Qualitative view
        ASSAY_DATA_QUALITATIVE = get_assay_data_qualitative(ASSAY_DATA)

        # Setting up some variables
        positives_qualitative, ratio_qualitative, compounds_qualitative = set_variables_qualitative(ASSAY_DATA_QUALITATIVE)

        # If direction is 1 or -1
        if direction in [+1, -1]:

            # Get value to adjust relations
            CUT = get_cut_value(ASSAY_DATA, direction)

            # Adjust relation
            ASSAY_DATA_QUANTITATIVE = adjust_relation(ASSAY_DATA_QUANTITATIVE, direction, CUT)

            # Disambiguate duplicated compounds and returns 'sorted' data (depending on direction)
            ASSAY_DATA_QUANTITATIVE = disambiguate_compounds(ASSAY_DATA_QUANTITATIVE, direction)

            # Binarization with expert cut-off
            ASSAY_DATA_QUANTITATIVE = binarize_with_expert_cutoff(ASSAY_DATA_QUANTITATIVE, expert_cutoff, direction)

            # Setting up some variables
            positives_quantitative, ratio_quantitative, compounds_quantitative, activities_quantitative = set_variables_quantitative(ASSAY_DATA_QUANTITATIVE)

            # Set the dataset type
            if np.isnan(expert_cutoff):
                if np.isnan(compounds_qualitative):
                    dataset_type = 'none'
                else:
                    dataset_type = 'qualitative'
            else:
                if np.isnan(compounds_qualitative):
                    dataset_type = 'quantitative'
                else:
                    dataset_type = 'mixed'
        else:

            # Qualitative assay
            dataset_type = 'qualitative'    

            # Setting up some variables
            positives_quantitative = np.nan
            ratio_quantitative = np.nan
            compounds_quantitative = np.nan
            activities_quantitative = np.nan

        # Get activity stats
        min_, p1, p25, p50, p75, p99, max_ = get_activity_stats_quantitative(activities_quantitative)

        # Store data range
        DATA_RANGES.append([assay_chembl_id, activity_type, unit, target_type, target_type_curated, activities, nan_values, cpds, direction, acc, stdtc, equal, higher, lower, dataset_type, 
                            expert_cutoff, positives_quantitative, ratio_quantitative, compounds_quantitative, min_, p1, p25, p50, p75, p99, max_, positives_qualitative, 
                            ratio_qualitative, compounds_qualitative])

        # Save data only if number of compounds is >= 100
        dataset_name = f"{assay_chembl_id}_{activity_type}_{str(unit).replace('/', 'FwdS')}"
        if compounds_quantitative >= 100 and np.isnan(expert_cutoff) == False:
            ASSAY_DATA_QUANTITATIVE.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{dataset_name}_qt.csv.gz"), index=False)
        if compounds_qualitative >= 100:
            ASSAY_DATA_QUALITATIVE.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{dataset_name}_ql.csv.gz"), index=False)


    DATA_RANGES = pd.DataFrame(DATA_RANGES, columns=["assay_id", "activity_type", "unit", "target_type", "target_type_curated", "activities", "nan_values", "cpds", "direction", 
                                                     'activity_comment_counts', 'standard_text_count', "equal", "higher", "lower", "dataset_type", "expert_cutoff", "pos_qt", 
                                                     "ratio_qt", "cpds_qt", "min_", "p1", "p25", "p50", "p75", "p99", "max_", "pos_ql", "ratio_ql", "cpds_ql"])
    DATA_RANGES.to_csv(os.path.join(OUTPUT, pathogen_code, 'assays_data.csv'), index=False)
    print("\n\n\n")
