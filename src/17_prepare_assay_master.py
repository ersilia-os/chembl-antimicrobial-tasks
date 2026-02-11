from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import roc_auc_score
from IPython.display import display, HTML
from scipy.stats import spearmanr
from collections import Counter, defaultdict
import pandas as pd
from tqdm import tqdm
import numpy as np
import zipfile
import random
import gzip
import sys
import h5py
import os

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 50)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.width", None)

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

print("Step 17: towards the assay master table")

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

# Define output directory
OUTPUT = os.path.join(root, "..", "output")

# Shared columns
KEYS = ["assay_id", "activity_type", "unit"]

# Load expert cut-offs
EXPERT_CUTOFFS = load_expert_cutoffs(CONFIGPATH)

# Columns to take from different tables
COLUMNS_CLUSTERS = ["clusters_0.3", "clusters_0.6", "clusters_0.85"]
COLUMNS_PARAMETERS = ['organism_curated', 'target_type_curated', 'target_name_curated', "target_chembl_id_curated", "strain", "atcc_id", "mutations", "known_drug_resistances", "media"]
COLUMNS_DATA_INFO = ["target_type_curated_extra", "dataset_type", "equal", 'higher', 'lower', "min_", "p1", "p25", "p50", "p75", "p99", "max_"]

# Load assays info
print("Getting data")
ASSAYS_CLEANED = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_cleaned.csv"))
ASSAYS_CLUSTERS = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_clusters.csv"))[KEYS + COLUMNS_CLUSTERS]
ASSAYS_PARAMETERS = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_parameters.csv"))[KEYS + COLUMNS_PARAMETERS]
ASSAY_DATA_INFO = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assay_data_info.csv"))[KEYS + COLUMNS_DATA_INFO]

# Merge tables
print("Integrating data")
ASSAYS_MASTER = ASSAYS_CLEANED.merge(ASSAYS_CLUSTERS,on=KEYS, how="left", validate="1:1")
ASSAYS_MASTER = ASSAYS_MASTER.merge(ASSAYS_PARAMETERS, on=KEYS, how="left", validate="1:1")
ASSAYS_MASTER = ASSAYS_MASTER.merge(ASSAY_DATA_INFO, on=KEYS, how="left", validate="1:1")

# Add used cutoffs
CUTOFFS = []
for i in ASSAYS_MASTER[['activity_type', 'unit', 'target_type_curated_extra']].values:
    lab = tuple([i[0], i[1], i[2], pathogen_code])
    if lab in EXPERT_CUTOFFS:
        cutoffs = ";".join([str(c) for c in EXPERT_CUTOFFS[lab]])
    else:
        cutoffs = np.nan
    CUTOFFS.append(cutoffs)
ASSAYS_MASTER["evaluated_cutoffs"] = CUTOFFS

# Get considered in A, B or M
individual_LM = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "individual_LM.csv"))
merged_LM = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "merged_LM.csv"))
considered_A = set([tuple(i) for i in individual_LM[individual_LM['label'] == "A"][KEYS].values])
considered_B = set([tuple(i) for i in individual_LM[individual_LM['label'] == "B"][KEYS].values])
considered_M = set([eval(assay) for row in merged_LM['assay_keys'].tolist() for assay in row.split(";")])

# Selected in A, B or M
individual_selected_LM = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "individual_selected_LM.csv"))
merged_selected_LM = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "merged_selected_LM.csv"))
selected_A = set([tuple(i) for i in individual_selected_LM[individual_selected_LM['label'] == "A"][KEYS].values])
selected_B = set([tuple(i) for i in individual_selected_LM[individual_selected_LM['label'] == "B"][KEYS].values])
selected_M = set([eval(assay) for row in merged_selected_LM['assay_keys'].tolist() for assay in row.split(";")])

# Get final considered
final_datasets = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "final_datasets.csv"))
final_considered_A = set([eval(assay) for row in final_datasets[final_datasets['label'] == 'A']['assay_keys'].tolist() for assay in row.split(";")])
final_considered_B = set([eval(assay) for row in final_datasets[final_datasets['label'] == 'B']['assay_keys'].tolist() for assay in row.split(";")])
final_considered_M = set([eval(assay) for row in final_datasets[final_datasets['label'] == 'M']['assay_keys'].tolist() for assay in row.split(";")])

# Get final selected
final_datasets = final_datasets[final_datasets['selected']].reset_index(drop=True)
final_selected_A = set([eval(assay) for row in final_datasets[final_datasets['label'] == 'A']['assay_keys'].tolist() for assay in row.split(";")])
final_selected_B = set([eval(assay) for row in final_datasets[final_datasets['label'] == 'B']['assay_keys'].tolist() for assay in row.split(";")])
final_selected_M = set([eval(assay) for row in final_datasets[final_datasets['label'] == 'M']['assay_keys'].tolist() for assay in row.split(";")])

COMMENTS_A, COMMENTS_B, COMMENTS_M = [], [], []

for assay_id, activity_type, unit in ASSAYS_MASTER[KEYS].values:

    # Get assay identifier
    key = (assay_id, activity_type, unit)

    # A
    if key not in considered_A:
        COMMENTS_A.append("Not considered")
    elif key not in selected_A:
        COMMENTS_A.append("Considered but not selected")
    elif key not in final_considered_A:
        COMMENTS_A.append("Considered and selected, but is not ORGANISM assay")
    elif key not in final_selected_A:
        COMMENTS_A.append("Considered and selected, but discarded in the final selection due to correlations")
    elif key in final_selected_A:
        COMMENTS_A.append("Considered and selected, selected in the final selection")

    # B
    if key not in considered_B:
        COMMENTS_B.append("Not considered")
    elif key not in selected_B:
        COMMENTS_B.append("Considered but not selected")
    elif key not in final_considered_B:
        COMMENTS_B.append("Considered and selected, but is not ORGANISM assay")
    elif key not in final_selected_B:
        COMMENTS_B.append("Considered and selected, but discarded in the final selection due to correlations")
    elif key in final_selected_B:
        COMMENTS_B.append("Considered and selected, selected in the final selection")

    # M
    if key not in considered_M:
        COMMENTS_M.append("Not considered")
    elif key not in selected_M:
        COMMENTS_M.append("Considered but not selected")
    elif key not in final_considered_M:
        COMMENTS_M.append("Considered and selected, but is not ORGANISM assay")
    elif key not in final_selected_M:
        COMMENTS_M.append("Considered and selected, but discarded in the final selection due to correlations")
    elif key in final_selected_M:
        COMMENTS_M.append("Considered and selected, selected in the final selection")

ASSAYS_MASTER["comment_A"] = COMMENTS_A
ASSAYS_MASTER["comment_B"] = COMMENTS_B
ASSAYS_MASTER["comment_M"] = COMMENTS_M

# All columns to include
ALL_COLS = ["assay_id", "assay_type", "assay_organism", "target_organism", "organism_curated", "doc_chembl_id", "target_type", "target_type_curated", "target_type_curated_extra", 
          "target_chembl_id", "target_chembl_id_curated", "target_name_curated", "bao_label", "source_label", "strain", "atcc_id", "mutations", "known_drug_resistances", "media",
          "activity_type", "unit", "activities", "nan_values", "cpds", "frac_cs", "direction", "act_flag", 'inact_flag', "equal", "higher", "lower", "dataset_type", "evaluated_cutoffs", 
          "min_", "p1", "p25", "p50", "p75", "p99", "max_", 'comment_A', 'comment_B', 'comment_M']

# Select columns
ASSAYS_MASTER = ASSAYS_MASTER[ALL_COLS]

# Save results
ASSAYS_MASTER.to_csv(os.path.join(OUTPUT, pathogen_code, "assays_master.csv"), index=False)

print("Assay master table done!")