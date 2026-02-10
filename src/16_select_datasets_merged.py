from collections import defaultdict
import pandas as pd
import numpy as np
import zipfile
import gzip
import sys
import os

# pd.set_option("display.max_columns", None)
# pd.set_option("display.max_rows", 50)
# pd.set_option("display.max_colwidth", None)
# pd.set_option("display.width", None)

def load_expert_cutoffs(CONFIGPATH):
    """
    Load expert cutoffs from the manual curation CSV and return them as a dictionary.

    The CSV is expected at:
        {CONFIGPATH}/expert_cutoffs.csv

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

def get_filtered_data(individual_LM_LABEL, assay_id, activity_type, unit):
    if type(unit) == str:
        df = individual_LM_LABEL[(individual_LM_LABEL['assay_id'] == assay_id) & 
                                 (individual_LM_LABEL['activity_type'] == activity_type) & 
                                 (individual_LM_LABEL['unit'] == unit)].reset_index(drop=True)
    else:
        df = individual_LM_LABEL[(individual_LM_LABEL['assay_id'] == assay_id) & 
                                 (individual_LM_LABEL['activity_type'] == activity_type) & 
                                 (individual_LM_LABEL['unit'].isna())].reset_index(drop=True)
    return df

def load_data_from_zip(zip_path, filename):
    """Load a gzipped CSV file from a ZIP archive into a pandas DataFrame.

    Parameters
    ----------
    zip_path : str
        Path to the ZIP archive.
    filename : str
        Name of the gzipped CSV file inside the ZIP.

    Returns
    -------
    pandas.DataFrame
        Loaded data.
    """
    with zipfile.ZipFile(zip_path) as z:
        with z.open(filename) as raw:
            with gzip.open(raw, mode="rt") as f:
                df = pd.read_csv(f)
    return df

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

# Define root directory
# root = os.path.dirname(os.path.abspath(__file__))
root = "."
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH

# Load pathogen info
pathogen_code = sys.argv[1]
pathogen_code = "mtuberculosis"
df = pd.read_csv(os.path.join(CONFIGPATH, 'pathogens.csv'))
row = df.loc[df["code"].eq(pathogen_code)]
if row.empty: 
    raise SystemExit(f"Unknown code: {pathogen_code}")
pathogen = row.iloc[0]["pathogen"]

print(f"Step 16: Selecting merged datasets")

# Define output directory
OUTPUT = os.path.join(root, "..", "output")

# Load ChEMBL data for pathogen
ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)

# Get assay to index mapping
assay_to_idx = defaultdict(list)
for i, assay_id in enumerate(ChEMBL_pathogen["assay_chembl_id"].to_numpy()):
    assay_to_idx[assay_id].append(i)

# Load expert cut-offs
EXPERT_CUTOFFS = load_expert_cutoffs(CONFIGPATH)

# Load individual LM data
merged_LM = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "merged_LM.csv"))
MERGING_STRATEGIES = sorted(set(("_".join(i.split("_")[:2]), j, k, l) for i,j,k,l in zip(merged_LM['name'], merged_LM['activity_type'], merged_LM['unit'], merged_LM['target_type_curated_extra'])))

# Dict mapping assay_id, activity_type and unit to a set of compound ChEMBL IDs
ASSAY_TO_COMPOUNDS = defaultdict(set)
for assay_id, activity_type, unit, compound_chembl_id in ChEMBL_pathogen[["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"]].values:
    ASSAY_TO_COMPOUNDS[(assay_id, activity_type, unit)].add(compound_chembl_id)
del ChEMBL_pathogen

# Get all compounds for pathogen
compounds = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "compound_counts.csv.gz"))
compounds = set(compounds['compound_chembl_id'])

COLS_TO_KEEP = ["name", "direction", "assay_type", "bao_label", "strain", "target_chembl_id", "n_assays", "n_cpds_union", "positives", "ratio", "avg", "std", "assay_keys"]
SELECTED = []
ORIGINAL_COMPOUNDS = {i: set() for i in ['ORG', 'SP']}
SELECTED_COMPOUNDS = {i: set() for i in ['ORG', 'SP']}

for strategy in MERGING_STRATEGIES:

    activity_type = strategy[1]
    unit = strategy[2]
    target_type = strategy[3]
    strategy = strategy[0]

    ty = "".join([i for i in strategy.split("_")[1] if i not in "0123456789"])

    # Filter assays considered in merging strategy type
    merged_LM_strategy = merged_LM[merged_LM['name'].str.startswith(strategy)]
    assays_strategy = set([eval(assay) for row in merged_LM_strategy['assay_keys'].tolist() for assay in row.split(";")])

    # Get associated compounds
    for assay in assays_strategy:
        ORIGINAL_COMPOUNDS[ty].update(ASSAY_TO_COMPOUNDS[assay])

    # Define mid cutoff (available by definition)
    mid_cutoff = EXPERT_CUTOFFS[(activity_type, unit, target_type, pathogen_code)][1]

    # Sort by average AUROC
    df = merged_LM_strategy.sort_values("avg", ascending=False).reset_index(drop=True)

    # Get best auroc and best cutoff
    best_auroc = df["avg"].tolist()[0]
    best_cutoff = df["expert_cutoff"].tolist()[0]

    # Get mid auroc (if available)
    mid_auroc = df[df['expert_cutoff'] == mid_cutoff]["avg"].tolist()[0]

    # If the best dataset is modelable
    if best_auroc > 0.7:

        # If difference is quite high, keep best
        if np.isnan(mid_auroc) or (best_auroc - mid_auroc) > 0.1:
            INFO = df[COLS_TO_KEEP].values.tolist()[0]
            SELECTED.append([activity_type, unit, target_type, best_cutoff, best_auroc, False] + INFO)
            assay_keys = [eval(i) for i in SELECTED[-1][-1].split(";")]
            for assay in assay_keys:
                SELECTED_COMPOUNDS[ty].update(ASSAY_TO_COMPOUNDS[assay])

        # Otherwise, keep mid
        else:
            INFO = df[df['expert_cutoff'] == mid_cutoff][COLS_TO_KEEP].values.tolist()[0]
            SELECTED.append([activity_type, unit, target_type, mid_cutoff, mid_auroc, True] + INFO)
            assay_keys = [eval(i) for i in SELECTED[-1][-1].split(";")]
            for assay in assay_keys:
                SELECTED_COMPOUNDS[ty].update(ASSAY_TO_COMPOUNDS[assay])

# To pandas dataframe
SELECTED = pd.DataFrame(SELECTED, columns=['activity_type', 'unit', "target_type", 'cutoff', 'AUROC', "is_mid_cutoff", "name", "direction", 
                                           "assay_type", "bao_label", "strain", "target_chembl_id", "n_assays", "n_cpds_union", "positives", "ratio", "avg", "std", "assay_keys"])

# Save 
SELECTED.to_csv(os.path.join(OUTPUT, pathogen_code, "merged_selected_LM.csv"), index=False)

print("Chemical space coverage change:")
print(f"ORG: from {round(100 * len(ORIGINAL_COMPOUNDS['ORG']) / len(compounds), 1)}% to {round(100 * len(SELECTED_COMPOUNDS['ORG']) / len(compounds), 1)}%")
print(f"SP: from {round(100 * len(ORIGINAL_COMPOUNDS['SP']) / len(compounds), 1)}% to {round(100 * len(SELECTED_COMPOUNDS['SP']) / len(compounds), 1)}%")
print(f"Overall: from {round(100 * len(ORIGINAL_COMPOUNDS['ORG'].union(ORIGINAL_COMPOUNDS['SP'])) / len(compounds), 1)}% to {round(100 * len(SELECTED_COMPOUNDS['ORG'].union(SELECTED_COMPOUNDS['SP'])) / len(compounds), 1)}%")
print(f"Number of selected datasets (ORG): {len(SELECTED[SELECTED['target_type'] == 'ORG'])}")
print(f"Number of selected datasets (SP): {len(SELECTED[SELECTED['target_type'] == 'SP'])}")
print(f"Number of datasets with middle cut-off: {len(SELECTED[SELECTED['is_mid_cutoff'] == True])}/{len(SELECTED)}")