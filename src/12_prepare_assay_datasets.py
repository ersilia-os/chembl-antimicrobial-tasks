from collections import defaultdict
from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import zipfile
import json
import sys
import os

# pd.set_option("display.max_columns", None)
# pd.set_option("display.max_rows", 50)
# pd.set_option("display.max_colwidth", None)
# pd.set_option("display.width", None)

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

def add_target_type_curated(ASSAYS_CLEANED, PARAMETERS):
    """
    Add a `target_type_curated` column to ASSAYS_CLEANED using curated assay parameters.

    For each row in ASSAYS_CLEANED, the function matches on:
        [`assay_id`, `activity_type`, `unit`]

    and pulls the corresponding value from the `target_type_curated` column in PARAMETERS.

    The function enforces that:
      - all keys in PARAMETERS exist in ASSAYS_CLEANED
      - all keys in ASSAYS_CLEANED exist in PARAMETERS

    Parameters
    ----------
    ASSAYS_CLEANED : pandas.DataFrame
        DataFrame containing at least the columns:
        `assay_id`, `activity_type`, `unit`.
    PARAMETERS : pandas.DataFrame
        DataFrame containing curated assay parameters, including
        `assay_id`, `activity_type`, `unit`, and `target_type_curated`.

    Returns
    -------
    pandas.DataFrame
        ASSAYS_CLEANED with an added `target_type_curated` column.
    """
    # Load only the columns we need from the parameters table
    PARAMETERS = PARAMETERS[["assay_id", "activity_type", "unit", "target_type_curated"]].copy()

    # Match writer behavior from Step 11: missing units are stored as empty strings
    ASSAYS_CLEANED = ASSAYS_CLEANED.copy()
    ASSAYS_CLEANED["unit"] = ASSAYS_CLEANED["unit"].fillna("")
    PARAMETERS["unit"] = PARAMETERS["unit"].fillna("")

    # Check that everything in PARAMETERS actually maps to ASSAYS_CLEANED
    if not PARAMETERS[["assay_id","activity_type","unit"]].merge(
        ASSAYS_CLEANED[["assay_id","activity_type","unit"]],
        on=["assay_id","activity_type","unit"],
        how="left",
        indicator=True
    )["_merge"].eq("both").all():
        raise ValueError("PARAMETERS contains keys not present in ASSAYS_CLEANED")
    
    # Check that everything in ASSAYS_CLEANED actually maps to PARAMETERS
    if not ASSAYS_CLEANED[["assay_id","activity_type","unit"]].merge(
        PARAMETERS[["assay_id","activity_type","unit"]],
        on=["assay_id","activity_type","unit"],
        how="left",
        indicator=True
    )["_merge"].eq("both").all():
        raise ValueError("ASSAYS_CLEANED contains keys not present in PARAMETERS")

    # Merge curated target type onto the cleaned assays table
    ASSAYS_CLEANED = ASSAYS_CLEANED.merge(PARAMETERS, on=["assay_id", "activity_type", "unit"], how="left", validate="1:1")

    # Unit "" back to nans
    ASSAYS_CLEANED["unit"] = ASSAYS_CLEANED["unit"].replace("", np.nan)

    return ASSAYS_CLEANED

def load_expert_cutoffs(CONFIGPATH):
    """
    Load expert cutoffs from the manual curation CSV and return them as a dictionary.

    The CSV is expected at:
        {CONFIGPATH}/manual_curation/expert_cutoffs.csv

    The returned dictionary maps:
        (activity_type, unit, target_type, pathogen_code) -> expert_cutoff

    Parameters
    ----------
    CONFIGPATH : str
        Path to the config folder.

    Returns
    -------
    dict
        Dictionary of expert cutoffs keyed by
        (activity_type, unit, target_type, pathogen_code).
    """
    # Load expert cut-offs
    EXPERT_CUTOFFS = pd.read_csv(os.path.join(CONFIGPATH, "expert_cutoffs.csv"))

    EXPERT_CUTOFFS = {
        (a, b, c, d): [float(k) for k in e.split(";")]
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
    ASSAY_DATA_QUANTITATIVE = ASSAY_DATA_QUANTITATIVE.drop(columns=['text_flag'])
    return ASSAY_DATA_QUANTITATIVE

def get_assay_data_qualitative(ASSAY_DATA):
    """
    Build a compound-level qualitative dataset from assay text flags.

    Aggregates `text_flag` values per compound, checks for conflicting
    qualitative labels, assigns a final compound-level label, and
    returns one row per compound with a binary activity label.

    Parameters
    ----------
    ASSAY_DATA : pandas.DataFrame
        DataFrame containing at least `compound_chembl_id` and `text_flag`.

    Returns
    -------
    pandas.DataFrame
        Qualitative dataset with one row per compound and a binary label.
    """

    ASSAY_DATA_QUALITATIVE = ASSAY_DATA.copy()

    # Aggregate to compound-level label
    compound_labels = ASSAY_DATA_QUALITATIVE.groupby("compound_chembl_id")["text_flag"].apply(set)

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
    cols = ["compound_chembl_id", "canonical_smiles", "activity_type", "unit", "text_flag", "qualitative_label", 'bin'] 
    ASSAY_DATA_QUALITATIVE = ASSAY_DATA_QUALITATIVE[cols]

    return ASSAY_DATA_QUALITATIVE
    
def set_variables_quantitative(ASSAY_DATA_QUANTITATIVE):
    """
    Compute basic statistics for a quantitative, binarized assay dataset.

    Parameters
    ----------
    ASSAY_DATA_QUANTITATIVE : pandas.DataFrame
        Quantitative assay data containing one row per compound and a
        binary `bin` column.

    Returns
    -------
    tuple
        (positives_quantitative, ratio_quantitative,
         compounds_quantitative, activities_quantitative)
    """
    positives_quantitative = (ASSAY_DATA_QUANTITATIVE["bin"] == 1).sum()
    ratio_quantitative = round(positives_quantitative / len(ASSAY_DATA_QUANTITATIVE), 5)
    compounds_quantitative = len(set(ASSAY_DATA_QUANTITATIVE['compound_chembl_id']))
    activities_quantitative = ASSAY_DATA_QUANTITATIVE['value'].tolist()
    assert compounds_quantitative == len(activities_quantitative)

    return positives_quantitative, ratio_quantitative, compounds_quantitative, activities_quantitative

def set_variables_qualitative(ASSAY_DATA_QUALITATIVE):
    """
    Compute basic statistics for a qualitative assay dataset.

    Parameters
    ----------
    ASSAY_DATA_QUALITATIVE : pandas.DataFrame
        Qualitative assay data containing one row per compound and a
        binary `bin` column.

    Returns
    -------
    tuple
        (positives_qualitative, ratio_qualitative, compounds_qualitative)
    """
    positives_qualitative = (ASSAY_DATA_QUALITATIVE["bin"] == 1).sum()
    ratio_qualitative = round(positives_qualitative / len(ASSAY_DATA_QUALITATIVE), 5)
    compounds_qualitative = len(set(ASSAY_DATA_QUALITATIVE['compound_chembl_id']))
    assert compounds_qualitative == len(ASSAY_DATA_QUALITATIVE)

    return positives_qualitative, ratio_qualitative, compounds_qualitative

def binarize_with_expert_cutoff(ASSAY_DATA_QUANTITATIVE, expert_cutoff, direction):
    """
    Binarize quantitative assay values using an expert-defined cutoff.

    Parameters
    ----------
    ASSAY_DATA_QUANTITATIVE : pandas.DataFrame
        Quantitative assay data containing a numeric `value` column.
    expert_cutoff : float
        Expert-defined activity threshold.
    direction : int
        +1 → higher values are more active
        -1 → lower values are more active

    Returns
    -------
    pandas.DataFrame
        Input dataframe with an added binary `bin` column.
    """
    if direction == +1:
        ASSAY_DATA_QUANTITATIVE["bin"] = (ASSAY_DATA_QUANTITATIVE["value"] >= expert_cutoff).astype(int)
    else:
        ASSAY_DATA_QUANTITATIVE["bin"] = (ASSAY_DATA_QUANTITATIVE["value"] <= expert_cutoff).astype(int)

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

def extra_curation_target_type(target_type, target_type_curated):
    """
    Post-process and enforce simple curation rules for ChEMBL assay target types.

    This function normalizes `target_type` and `target_type_curated` (strip + uppercase)
    and returns a constrained/standardized `target_type_curated` according to rules:

    - If target_type == "UNCHECKED": allow only {"ORGANISM", "SINGLE PROTEIN", "DISCARDED"};
      otherwise force "DISCARDED".
    - If target_type == "NOT-MOLECULAR": allow only {"ORGANISM", "DISCARDED"};
      otherwise force "DISCARDED".
    - If target_type is protein-related ({"SINGLE PROTEIN","PROTEIN COMPLEX","PROTEIN FAMILY"}):
      collapse to "SINGLE PROTEIN".
    - If target_type == "ORGANISM": return "ORGANISM".
    - All other target types are mapped to "DISCARDED".

    Parameters
    ----------
    target_type : str
        Raw ChEMBL "Target type" value.
    target_type_curated : str
        LLM- or human-proposed curated target type.

    Returns
    -------
    str
        The finalized curated target type: "SINGLE PROTEIN", "ORGANISM", or "DISCARDED".
    """

    if type(target_type_curated) != str:
        return 'DISCARDED'

    target_type = target_type.strip().upper()
    target_type_curated = target_type_curated.strip().upper()

    if target_type == 'UNCHECKED':
        if target_type_curated in ['ORGANISM', 'SINGLE PROTEIN', 'DISCARDED']:
            return target_type_curated
        else:
            return 'DISCARDED'
        
    elif target_type == 'NOT-MOLECULAR':
        if target_type_curated in ['ORGANISM', 'DISCARDED']:
            return target_type_curated
        else:
            return 'DISCARDED'
        
    elif target_type in ['SINGLE PROTEIN', 'PROTEIN COMPLEX', 'PROTEIN FAMILY']:
        return 'SINGLE PROTEIN'

    elif target_type in ['ORGANISM']:
        return 'ORGANISM'
    
    else:
        return 'DISCARDED'

def zip_and_remove(datasets_dir):
    """
    Create three zip archives in `datasets_dir`:
      - datasets_qt.zip containing all files ending with "_qt.csv.gz"
      - datasets_ql.zip containing all files ending with "_ql.csv.gz"
      - datasets_mx.zip containing all files ending with "_mx.csv.gz"

    After zipping, remove the original *_qt.csv.gz, *_ql.csv.gz, and *_mx.csv.gz files.

    Parameters
    ----------
    datasets_dir : str
        Directory containing the dataset files.

    Returns
    -------
    tuple[int, int, int]
        (n_qt_files, n_ql_files, n_mx_files)
    """
    qt_zip = os.path.join(datasets_dir, "datasets_qt.zip")
    ql_zip = os.path.join(datasets_dir, "datasets_ql.zip")
    mx_zip = os.path.join(datasets_dir, "datasets_mx.zip")

    if os.path.exists(qt_zip):
        os.remove(qt_zip)
    if os.path.exists(ql_zip):
        os.remove(ql_zip)
    if os.path.exists(mx_zip):
        os.remove(mx_zip)

    qt_files = [os.path.join(datasets_dir, i) for i in os.listdir(datasets_dir) if "_qt_" in i]  # cutoff specified in file name
    ql_files = [os.path.join(datasets_dir, i) for i in os.listdir(datasets_dir) if i.endswith("_ql.csv.gz")]
    mx_files = [os.path.join(datasets_dir, i) for i in os.listdir(datasets_dir) if "_mx_" in i]  # cutoff specified in file name

    with zipfile.ZipFile(qt_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as z:
        for fp in qt_files:
            z.write(fp, arcname=os.path.basename(fp))

    with zipfile.ZipFile(ql_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as z:
        for fp in ql_files:
            z.write(fp, arcname=os.path.basename(fp))

    with zipfile.ZipFile(mx_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as z:
        for fp in mx_files:
            z.write(fp, arcname=os.path.basename(fp))

    for fp in qt_files + ql_files + mx_files:
        os.remove(fp)

    return len(qt_files), len(ql_files), len(mx_files)


# Define root directory
# root = os.path.dirname(os.path.abspath(__file__))
root = '.'
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH

# Load pathogen info
# pathogen_code = sys.argv[1]
pathogen_code = 'mtuberculosis'
df = pd.read_csv(os.path.join(CONFIGPATH, 'pathogens.csv'))
row = df.loc[df["code"].eq(pathogen_code)]
if row.empty: 
    raise SystemExit(f"Unknown code: {pathogen_code}")
pathogen = row.iloc[0]["pathogen"]

print("Step 12")

# Define output directory
OUTPUT = os.path.join(root, "..", "output")

# Load cleaned assays
ASSAYS_CLEANED = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_cleaned.csv"))

# Define PATH to parameters
PARAMETERS = pd.read_csv(os.path.join(OUTPUT, pathogen_code, 'assays_parameters.csv'))

# Get curated target type
ASSAYS_CLEANED = add_target_type_curated(ASSAYS_CLEANED, PARAMETERS)

# Extra curation
ASSAYS_CLEANED['target_type_curated_extra'] = [extra_curation_target_type(i,j) for i,j in zip(ASSAYS_CLEANED['target_type'], ASSAYS_CLEANED['target_type_curated'])]

# Loading pathogen data
os.makedirs(os.path.join(OUTPUT, pathogen_code, 'datasets'), exist_ok=True)
print(f"Loading ChEMBL cleaned data for {pathogen_code}...")
ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
print(f"Number of cleaned assays: {len(ASSAYS_CLEANED)}")

# Load expert cut-offs
EXPERT_CUTOFFS = load_expert_cutoffs(CONFIGPATH)

# Get assay to index mapping
assay_to_idx = defaultdict(list)
for i, assay_id in enumerate(ChEMBL_pathogen["assay_chembl_id"].to_numpy()):
    assay_to_idx[assay_id].append(i)


# Define data ranges
DATASETS, ASSAY_DATA_INFO = [], []

print("Preparing datasets for each assay")

for assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra, activities, nan_values, cpds, direction, act_flag, inact_flag in tqdm(ASSAYS_CLEANED[['assay_id', 
                                    'activity_type', 'unit', 'target_type','target_type_curated_extra', 'activities', 'nan_values', 'cpds', 'direction', 
                                    'act_flag', 'inact_flag']].values[:]):

    # Filtering [assay, activity_type, unit] data
    cols = ['compound_chembl_id', 'canonical_smiles', 'activity_type', 'value', 'relation', 'unit', 'text_flag']
    tmp_df = ChEMBL_pathogen.iloc[assay_to_idx[assay_chembl_id]]    
    ASSAY_DATA = get_assay_data(tmp_df, assay_chembl_id, activity_type, unit, cols)

    # Get dataset name
    dataset_name = f"{assay_chembl_id}_{activity_type}_{str(unit).replace('/', 'FwdS')}"
    
    # Count relations
    equal, lower, higher = count_relations(ASSAY_DATA)

    # Mixed data, nan by default
    positives_mixed, ratio_mixed, compounds_mixed, overlap_mixed = [np.nan] * 4

    # Qualitative view
    ASSAY_DATA_QUALITATIVE = get_assay_data_qualitative(ASSAY_DATA)

    # Setting up some variables
    if len(ASSAY_DATA_QUALITATIVE) > 0:
        positives_qualitative, ratio_qualitative, compounds_qualitative = set_variables_qualitative(ASSAY_DATA_QUALITATIVE)
    else:
        positives_qualitative, ratio_qualitative, compounds_qualitative = [np.nan] * 3

    # Quantitative view
    ASSAY_DATA_QUANTITATIVE = get_assay_data_quantitative(ASSAY_DATA)

    # Get expert cut-offs if existing
    key = (activity_type, unit, target_type_curated_extra, pathogen_code)
    expert_cutoffs = EXPERT_CUTOFFS[key] if key in EXPERT_CUTOFFS else [np.nan]

    if len(ASSAY_DATA_QUANTITATIVE) == 0:

        if np.isnan(compounds_qualitative):
            raise ValueError("Dataset does not have numerical values nor activity flags. By definition, this is not possible at this stage. Please revise.")
        else:
            dataset_type = 'qualitative'

            # Quantitative binarization is not possible
            positives_quantitative, ratio_quantitative, compounds_quantitative, activities_quantitative = [np.nan] * 4
            min_, p1, p25, p50, p75, p99, max_ = [np.nan] * 7
            expert_cutoff = np.nan

            # Store data range
            DATASETS.append([assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra, activities, nan_values, cpds, direction, act_flag, 
                                inact_flag, dataset_type, expert_cutoff, positives_quantitative, ratio_quantitative, compounds_quantitative, positives_qualitative, 
                                ratio_qualitative, compounds_qualitative, overlap_mixed, positives_mixed, ratio_mixed, compounds_mixed])
            
            # Store dataset
            ASSAY_DATA_QUALITATIVE.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{dataset_name}_ql.csv.gz"), index=False)

    else:

        # For each expert_cutoff
        for expert_cutoff in expert_cutoffs:

            # If expert cutoff is nan
            if np.isnan(expert_cutoff) == True or direction not in [-1, +1]:

                # Quantitative binarization is not possible
                positives_quantitative = np.nan
                ratio_quantitative = np.nan
                compounds_quantitative = len(set(ASSAY_DATA_QUANTITATIVE['compound_chembl_id']))
                activities_quantitative = ASSAY_DATA_QUANTITATIVE['value'].tolist()
                min_, p1, p25, p50, p75, p99, max_ = get_activity_stats_quantitative(activities_quantitative)

                # Assess qualitative compounds
                if np.isnan(compounds_qualitative):
                    dataset_type = 'none'
                else:
                    dataset_type = 'qualitative'
                    
                    # Store dataset
                    ASSAY_DATA_QUALITATIVE.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{dataset_name}_ql.csv.gz"), index=False)

            else:
                
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

                # Get activity stats
                min_, p1, p25, p50, p75, p99, max_ = get_activity_stats_quantitative(activities_quantitative)

                if np.isnan(compounds_qualitative):
                    dataset_type = 'quantitative'

                    # Store dataset
                    ASSAY_DATA_QUANTITATIVE.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{dataset_name}_qt_{expert_cutoff}.csv.gz"), index=False)

                else:
                    dataset_type = 'mixed'

                    # Get overlap mixed
                    overlap_mixed = set(ASSAY_DATA_QUANTITATIVE['compound_chembl_id']).intersection(set(ASSAY_DATA_QUALITATIVE['compound_chembl_id']))
                    overlap_mixed = round(len(overlap_mixed) / min(len(ASSAY_DATA_QUANTITATIVE), len(ASSAY_DATA_QUALITATIVE)), 3)

                    # Get compounds in quantitative
                    qt_compounds = set(ASSAY_DATA_QUANTITATIVE['compound_chembl_id'])

                    # Prepare qualitative inactives
                    ASSAY_DATA_QUALITATIVE_MIXED = ASSAY_DATA_QUALITATIVE[(ASSAY_DATA_QUALITATIVE['compound_chembl_id'].isin(qt_compounds) == False) & 
                                                (ASSAY_DATA_QUALITATIVE['bin'] == 0)].reset_index(drop=True).copy()
                    ASSAY_DATA_QUALITATIVE_MIXED['value'] = np.nan
                    ASSAY_DATA_QUALITATIVE_MIXED['relation'] = np.nan

                    # Prepare quantitative compounds
                    ASSAY_DATA_QUANTITATIVE_MIXED = ASSAY_DATA_QUANTITATIVE.copy()
                    ASSAY_DATA_QUANTITATIVE_MIXED['text_flag'] = np.nan
                    ASSAY_DATA_QUANTITATIVE_MIXED['qualitative_label'] = np.nan

                    # Append inactive qualitatives to quantitative dataset [mixed dataset]
                    ASSAY_DATA_MIXED = pd.concat([ASSAY_DATA_QUANTITATIVE_MIXED, ASSAY_DATA_QUALITATIVE_MIXED], axis=0).reset_index(drop=True)

                    # Get metadata
                    positives_mixed, ratio_mixed, compounds_mixed, activities_mixed = set_variables_quantitative(ASSAY_DATA_MIXED)

                    # Store dataset
                    ASSAY_DATA_MIXED.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets', f"{dataset_name}_mx_{expert_cutoff}.csv.gz"), index=False)

            # Store data range
            DATASETS.append([assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra, activities, nan_values, cpds, direction, act_flag, 
                                inact_flag, dataset_type, expert_cutoff, positives_quantitative, ratio_quantitative, compounds_quantitative, positives_qualitative, 
                                ratio_qualitative, compounds_qualitative, overlap_mixed, positives_mixed, ratio_mixed, compounds_mixed])
            
    # Store assay data type and range
    ASSAY_DATA_INFO.append([assay_chembl_id, activity_type, unit, target_type, target_type_curated_extra, activities, cpds, dataset_type, equal, higher, lower, min_, p1, p25, p50, p75, p99, max_])
            
       
# Store data results
DATASETS = pd.DataFrame(DATASETS, columns=["assay_id", "activity_type", "unit", "target_type", "target_type_curated_extra", "activities", "nan_values", "cpds", "direction", 
                                                    'act_flag', 'inact_flag', "dataset_type", "expert_cutoff", "pos_qt", "ratio_qt", "cpds_qt", "pos_ql", "ratio_ql", "cpds_ql", 
                                                    "overlap_mx", "pos_mx", "ratio_mx", "cpds_mx"])

ASSAY_DATA_INFO = pd.DataFrame(ASSAY_DATA_INFO, columns=["assay_id", "activity_type", "unit", "target_type", "target_type_curated_extra", "activities", "cpds", "dataset_type", 
                                                         "equal", "higher", "lower", "min_", "p1", "p25", "p50", "p75", "p99", "max_"])

DATASETS.to_csv(os.path.join(OUTPUT, pathogen_code, 'datasets.csv'), index=False)
ASSAY_DATA_INFO.to_csv(os.path.join(OUTPUT, pathogen_code, 'assay_data_info.csv'), index=False)

# Zip and remove datasets
datasets_dir = os.path.join(OUTPUT, pathogen_code, "datasets")
qt, ql, mx = zip_and_remove(datasets_dir)

# Counting datasets and assays data types
counter_datasets = dict(Counter(DATASETS['dataset_type']))
counter_assays = dict(Counter(ASSAY_DATA_INFO['dataset_type']))
# assert len(ASSAYS_CLEANED) == int(counter['quantitative'] / 3 + counter['qualitative'] + counter['none'] + counter['mixed'] / 3)

print(f"Total number of assays: {len(ASSAY_DATA_INFO)}")
print(f"Types of assays: {counter_assays}")
print(f"Total number of datasets: {len(DATASETS)}")
print(f"Types of datasets: {counter_datasets}")