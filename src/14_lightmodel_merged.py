from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import roc_auc_score
from IPython.display import display, HTML
from scipy.stats import spearmanr
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import zipfile
import random
import gzip
import sys
import h5py
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

print("Step 14")

# Define output directory
OUTPUT = os.path.join(root, "..", "output")

# Shared columns
KEYS = ["assay_id", "activity_type", "unit"]

# Columns to take from datasets table
COLUMNS_DATASETS = ["equal", 'higher', 'lower', "target_type_curated_extra", "dataset_type", "cpds_qt", "min_", "p1", "p25", "p50", "p75", "p99", "max_", "pos_ql", "ratio_ql", "cpds_ql"]

def load_all_gz_csvs_from_zip(zip_path):
    """Read all ``*.csv.gz`` members from a ZIP into DataFrames.

    Parameters
    ----------
    zip_path : str | pathlib.Path
        Path to the ZIP archive.

    Returns
    -------
    dict[str, pandas.DataFrame]
        Mapping of ZIP member name -> loaded DataFrame.
    """
    dfs = {}
    with zipfile.ZipFile(zip_path, "r") as z:
        for name in z.namelist():
            if name.endswith(".csv.gz"):
                with z.open(name) as f:
                    dfs[name] = pd.read_csv(f, compression="gzip")
    return dfs

def get_all_results_from_individual_modeling(INDIVIDUAL_LM, LABELS=['A', 'B', 'C', 'D']):
    """Collect best AUROC (>0.7) per (assay_id, activity_type, unit) for each label.

    Returns
    -------
    RESULTS : dict[str, dict[tuple, list]]
        Per label: (assay_id, activity_type, unit) -> [expert_cutoff, best_auroc],
        considering only rows with AUROC > 0.7.
    CONSIDERED_ASSAYS : dict[str, set[tuple]]
        Per label: set of all (assay_id, activity_type, unit) keys encountered
        (no AUROC threshold applied).
    """
    RESULTS, CONSIDERED_ASSAYS = {}, {}
    for LABEL in LABELS:
        RESULTS[LABEL] = {}
        CONSIDERED_ASSAYS[LABEL] = set()
        rows = INDIVIDUAL_LM[INDIVIDUAL_LM[LABEL]][["assay_id", "activity_type", "unit", "expert_cutoff", f"{LABEL}_AVG"]].values
        for assay_id, activity_type, unit, expert_cutoff, auroc in rows:
            key = (assay_id, activity_type, unit)
            CONSIDERED_ASSAYS[LABEL].add(key)
            if auroc > 0.7:
                if key not in RESULTS[LABEL]:
                    RESULTS[LABEL][key] = [expert_cutoff, auroc]
                elif auroc > RESULTS[LABEL][key][1]:
                    RESULTS[LABEL][key] = [expert_cutoff, auroc]
    return RESULTS, CONSIDERED_ASSAYS

def where_considered(key, LABELS, CONSIDERED_ASSAYS):
    """Return labels (semicolon-separated) where `key` was considered; else NaN."""
    considered = []
    for LABEL in LABELS:
        if key in CONSIDERED_ASSAYS[LABEL]:
            considered.append(LABEL)
    if len(considered) > 0:
        return ";".join(considered)
    else:
        return np.nan
    
def where_accepted(key, LABELS, ACCEPTED_ASSAYS):
    """Return labels (semicolon-separated) where `key` was accepted; else NaN."""
    accepted = []
    for LABEL in LABELS:
        if key in ACCEPTED_ASSAYS[LABEL]:
            accepted.append(LABEL)
    if len(accepted) > 0:
        return ";".join(accepted)
    else:
        return np.nan

def get_filtered_assay_master_organism(assay_df, activity_type, unit, direction, assay_type, target_type_curated_extra, bao_label, strain):
    """Filter `assay_df` by metadata fields, treating non-string `unit` as missing (NaN)."""
    if type(unit) == str:
        if type(strain) == str:
            df = assay_df[(assay_df['activity_type'] == activity_type) & 
                        (assay_df['unit'] == unit) &
                        (assay_df['direction'] == direction) &
                        (assay_df['assay_type'] == assay_type) &
                        (assay_df['target_type_curated_extra'] == target_type_curated_extra) &
                        (assay_df['bao_label'] == bao_label) &
                        (assay_df['strain'] == strain)]
        else:
            df = assay_df[(assay_df['activity_type'] == activity_type) & 
                        (assay_df['unit'] == unit) &
                        (assay_df['direction'] == direction) &
                        (assay_df['assay_type'] == assay_type) &
                        (assay_df['target_type_curated_extra'] == target_type_curated_extra) &
                        (assay_df['bao_label'] == bao_label) &
                        (assay_df['strain'].isna())]
    else:
        if type(strain) == str:
            df = assay_df[(assay_df['activity_type'] == activity_type) & 
                        (assay_df['unit'].isna()) &
                        (assay_df['direction'] == direction) &
                        (assay_df['assay_type'] == assay_type) &
                        (assay_df['target_type_curated_extra'] == target_type_curated_extra) &
                        (assay_df['bao_label'] == bao_label) &
                        (assay_df['strain'] == strain)]
        else:
            df = assay_df[(assay_df['activity_type'] == activity_type) & 
                        (assay_df['unit'].isna()) &
                        (assay_df['direction'] == direction) &
                        (assay_df['assay_type'] == assay_type) &
                        (assay_df['target_type_curated_extra'] == target_type_curated_extra) &
                        (assay_df['bao_label'] == bao_label) &
                        (assay_df['strain'].isna())]
    return df

def get_filtered_assay_master_single_protein(assay_df, activity_type, unit, direction, assay_type, target_type_curated_extra, bao_label, strain, target_chembl_id):
    """Filter `assay_df` by metadata fields, treating non-string `unit` as missing (NaN)."""
    if type(unit) == str:
        if type(strain) == str:
            df = assay_df[(assay_df['activity_type'] == activity_type) & 
                        (assay_df['unit'] == unit) &
                        (assay_df['direction'] == direction) &
                        (assay_df['assay_type'] == assay_type) &
                        (assay_df['target_type_curated_extra'] == target_type_curated_extra) &
                        (assay_df['bao_label'] == bao_label) &
                        (assay_df['strain'] == strain) & 
                        (assay_df['target_chembl_id'] == target_chembl_id)]
        else:
            df = assay_df[(assay_df['activity_type'] == activity_type) & 
                        (assay_df['unit'] == unit) &
                        (assay_df['direction'] == direction) &
                        (assay_df['assay_type'] == assay_type) &
                        (assay_df['target_type_curated_extra'] == target_type_curated_extra) &
                        (assay_df['bao_label'] == bao_label) &
                        (assay_df['strain'].isna()) & 
                         (assay_df['target_chembl_id'] == target_chembl_id)]
    else:
        if type(strain) == str:
            df = assay_df[(assay_df['activity_type'] == activity_type) & 
                        (assay_df['unit'].isna()) &
                        (assay_df['direction'] == direction) &
                        (assay_df['assay_type'] == assay_type) &
                        (assay_df['target_type_curated_extra'] == target_type_curated_extra) &
                        (assay_df['bao_label'] == bao_label) &
                        (assay_df['strain'] == strain) & 
                        (assay_df['target_chembl_id'] == target_chembl_id)]
        else:
            df = assay_df[(assay_df['activity_type'] == activity_type) & 
                        (assay_df['unit'].isna()) &
                        (assay_df['direction'] == direction) &
                        (assay_df['assay_type'] == assay_type) &
                        (assay_df['target_type_curated_extra'] == target_type_curated_extra) &
                        (assay_df['bao_label'] == bao_label) &
                        (assay_df['strain'].isna()) & 
                        (assay_df['target_chembl_id'] == target_chembl_id)]
    return df

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

def load_ecfp_all(h5_path):
    """Load all ECFP (Morgan count) fingerprints.

    Parameters
    ----------
    h5_path : str
        Path to the HDF5 file containing datasets "SMILES" and "X_morgan".

    Returns
    -------
    dict[str, np.ndarray]
        Mapping {chembl_id: fingerprint (np.int8, shape (nBits,))}.
    """
    with h5py.File(h5_path, "r") as f:
        meta = f["SMILES"][:, 3].astype(str)
        fps  = f["X_morgan"][:]  # Load ALL

    return {cid: fp for cid, fp in zip(meta, fps)}

def KFoldTrain(X, Y, n_splits=4, n_estimators=100, random_state=42):
    """Stratified K-fold training/eval with RandomForest; returns mean AUROC and std.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix (n_samples, n_features).
    Y : np.ndarray
        Binary labels (n_samples,).
    n_splits : int
        Number of folds.
    n_estimators : int
        Number of trees in the random forest.
    random_state : int
        RNG seed (also used for fold shuffling).

    Returns
    -------
    tuple[float, float]
        (mean_auroc, std_auroc) rounded to 3 decimals.
    """
    def init_RF():
        return RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=None,
            min_samples_split=2,
            min_samples_leaf=1,
            max_features="sqrt",
            n_jobs=8,
            random_state=random_state,
        )

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    aurocs = []

    for train_idx, test_idx in skf.split(X, Y):
        X_train, X_test = X[train_idx], X[test_idx]
        Y_train, Y_test = Y[train_idx], Y[test_idx]
        rf = init_RF()
        rf.fit(X_train, Y_train)
        y_prob = rf.predict_proba(X_test)[:, 1]
        aurocs.append(roc_auc_score(Y_test, y_prob))

    return round(float(np.mean(aurocs)), 3), round(float(np.std(aurocs)), 3)

def TrainRF(X, Y, n_estimators=100):
    """Train a RandomForestClassifier on all provided data and return the fitted model.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix (n_samples, n_features).
    Y : np.ndarray
        Labels (n_samples,).

    Returns
    -------
    RandomForestClassifier
        Fitted classifier.
    """
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features="sqrt",
        n_jobs=8,
    )
    rf.fit(X, Y)
    return rf

def to_merge_unique_cpds(df, group_keys, assay_to_compounds):
    """
    Group assays by `group_keys` and compute:
      - n_assays: number of unique assays in the group
      - n_cpds_union: number of unique compounds across assays (set union)
      - assay_keys: ';'-separated tuple strings "(assay_id, activity_type, unit)" (last column)
    """

    def collect_assay_keys(block):
        """Return unique (assay_id, activity_type, unit) keys for this group."""
        keys = sorted({tuple(r) for r in block.values})
        return keys  # list[tuple]

    def union_size(keys):
        """Return size of union of compounds for the given assay keys."""
        u = set()
        for k in keys:
            u |= assay_to_compounds.get(k, set())
        return len(u)

    out = (df.groupby(group_keys, dropna=False)[["assay_id", "activity_type", "unit"]]
             .apply(collect_assay_keys)
             .reset_index(name="assay_keys"))

    out["n_assays"] = out["assay_keys"].apply(len)
    out["n_cpds_union"] = out["assay_keys"].apply(union_size)

    # store as ';'-separated tuple strings (easy round-trip via ast.literal_eval)
    out["assay_keys"] = out["assay_keys"].apply(lambda ks: ";".join(map(str, ks)))

    # make assay_keys the last column
    cols = [c for c in out.columns if c != "assay_keys"] + ["assay_keys"]
    out = out[cols]

    return out.sort_values("n_cpds_union", ascending=False).reset_index(drop=True)

def get_target_chembl_id(merging, target_type):
    if target_type == 'SINGLE PROTEIN':
        return merging.target_chembl_id
    else:
        return np.nan


RATIO = 0.1

# Set and create path to correlations
PATH_TO_CORRELATIONS = os.path.join(OUTPUT, pathogen_code, "correlations")
os.makedirs(os.path.join(PATH_TO_CORRELATIONS, "M"), exist_ok=True)

# Load assays info
print("Merging assay metadata")
ASSAYS_CLEANED = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_cleaned.csv"))
ASSAYS_PARAMETERS = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_parameters.csv"))
ASSAYS_DATASETS_ = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_datasets.csv"))
INDIVIDUAL_LM = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "individual_LM.csv"))

# Get assay to quantitative data info to collapse ASSAY_DATASETS_ (1 row per assay)
assay_to_qt_info = defaultdict(list)
for assay_id, activity_type, unit, expert_cutoff, ratio_qt in ASSAYS_DATASETS_[['assay_id', 'activity_type', 'unit', 'expert_cutoff', 'ratio_qt']].values:
    assay_to_qt_info[tuple([assay_id, activity_type, unit])].append([expert_cutoff, ratio_qt])

# Unique row per assay
ASSAYS_DATASETS = ASSAYS_DATASETS_[KEYS + COLUMNS_DATASETS].drop_duplicates().reset_index(drop=True)

# Get cutoffs and ratios
cutoffs = [";".join([str(j[0]) for j in assay_to_qt_info[tuple(i)]]) for i in ASSAYS_DATASETS[['assay_id', 'activity_type', 'unit']].values]
ratios = [";".join([str(j[1]) for j in assay_to_qt_info[tuple(i)]]) for i in ASSAYS_DATASETS[['assay_id', 'activity_type', 'unit']].values]
cutoffs = [i if i != 'nan' else np.nan for i in cutoffs]
ratios = [i if i != 'nan' else np.nan for i in ratios]

# Store results
ASSAYS_DATASETS.insert(8, 'cutoffs', cutoffs)
ASSAYS_DATASETS.insert(9, 'ratios', ratios)

# Merge everything
ASSAYS_MASTER = ASSAYS_CLEANED.merge(ASSAYS_PARAMETERS,on=KEYS, how="left", validate="1:1")
ASSAYS_MASTER = ASSAYS_MASTER.merge(ASSAYS_DATASETS,on=KEYS, how="left", validate="1:1")

# Dict mapping assay_id, activity_type and unit to a set of compound ChEMBL IDs
print("Mapping assays to compounds")
ChEMBL = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
ASSAY_TO_COMPOUNDS = defaultdict(set)
for assay_id, activity_type, unit, compound_chembl_id in ChEMBL[["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"]].values:
    ASSAY_TO_COMPOUNDS[(assay_id, activity_type, unit)].add(compound_chembl_id)
del ChEMBL

# Loading quantitative and qualitative datasets
print("Loading individual datasets")
qt_zip = os.path.join(OUTPUT, pathogen_code, "datasets", "datasets_qt.zip")
ql_zip = os.path.join(OUTPUT, pathogen_code, "datasets", "datasets_ql.zip")
dfs_qt = load_all_gz_csvs_from_zip(qt_zip)
dfs_ql = load_all_gz_csvs_from_zip(ql_zip)
print("Loaded quantitative:", len(dfs_qt), "datasets")
print("Loaded qualitative:", len(dfs_ql), "datasets")

# Get results from individual modeling ABCD
LABELS = ['A', 'B', 'C', 'D']
ACCEPTED_ASSAYS, CONSIDERED_ASSAYS = get_all_results_from_individual_modeling(INDIVIDUAL_LM, LABELS)

col_accepted, col_considered = [], []
for assay_id, activity_type, unit in ASSAYS_MASTER[["assay_id", "activity_type", "unit"]].values:
    # Get strategies in which this assay is considered and accepted
    key = tuple([assay_id, activity_type, unit])
    col_considered.append(where_considered(key, LABELS, CONSIDERED_ASSAYS))
    col_accepted.append(where_accepted(key, LABELS, ACCEPTED_ASSAYS))
ASSAYS_MASTER['Accepted'] = col_accepted
ASSAYS_MASTER['Considered'] = col_considered

# Reorder columns
ALL_COLS = ["assay_id", "assay_type", "assay_organism", "target_organism", "organism_curated", "doc_chembl_id", "target_type", "target_type_curated", "target_type_curated_extra", 
          "target_chembl_id", "target_chembl_id_curated", "target_name_curated", "bao_label", "source_label", "strain", "atcc_id", "mutations", "known_drug_resistances", "media",
          "activity_type", "unit", "activities", "nan_values", "cpds", "frac_cs", "direction", "act_flag", 'inact_flag', "equal", "higher", "lower", "dataset_type", "cutoffs", "ratios", 
          "cpds_qt", "pos_ql", "ratio_ql", "cpds_ql", "min_", "p1", "p25", "p50", "p75", "p99", "max_", 'Accepted', 'Considered']
ASSAYS_MASTER = ASSAYS_MASTER[ALL_COLS]

# Get accepted assays and accepted compounds in ABCD
accepted_assays = ASSAYS_MASTER[(ASSAYS_MASTER['Accepted'].isna() == False)][['assay_id', 'activity_type', 'unit']].values
accepted_compounds = set([j for i in accepted_assays for j in ASSAY_TO_COMPOUNDS[tuple(i)]])

# Loading Morgan fingerprints
print("Loading ECFPs...")
PATH_TO_ECFPs = os.path.join(DATAPATH, "chembl_processed", "ChEMBL_ECFPs.h5")
ecfps = load_ecfp_all(PATH_TO_ECFPs)

# Loading Reference set of compounds
print("Loading reference set of compounds")
REFERENCE_SET = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "reference_set.csv.gz"))['reference_compounds'].tolist()

# Prepare reference matrix of Morgan fingerprints
X_REF = np.array([ecfps[cid] for cid in REFERENCE_SET if cid in ecfps])

# Get all compounds from pathogen
compounds = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "compound_counts.csv.gz"))
compounds = set(compounds['compound_chembl_id'])

# Get ChEMBL compounds not tested against the pathogen
print("Defining decoys")
DECOYS_CHEMBL = set([i for i in ecfps if i not in compounds])

# Load expert cut-offs
EXPERT_CUTOFFS = load_expert_cutoffs(CONFIGPATH)

# Filtering assays
print("Identifying potential assays to merge")
print("Organisms...")
keys_organism = ["activity_type", "unit", "direction", "assay_type", "target_type_curated_extra", "bao_label", "strain"]
FILTERED_ASSAYS_ORGANISM = ASSAYS_MASTER[(ASSAYS_MASTER['Accepted'].isna()) & (ASSAYS_MASTER['target_type_curated_extra'] == 'ORGANISM')].copy()
TO_MERGE_ORGANISM = to_merge_unique_cpds(FILTERED_ASSAYS_ORGANISM, keys_organism, ASSAY_TO_COMPOUNDS)

print("Single proteins...")
keys_single_protein = ["activity_type", "unit", "direction", "assay_type", "target_type_curated_extra", "bao_label", "strain", 'target_chembl_id']
FILTERED_ASSAYS_SINGLE_PROTEIN = ASSAYS_MASTER[(ASSAYS_MASTER['Accepted'].isna()) & (ASSAYS_MASTER['target_type_curated_extra'] == 'SINGLE PROTEIN')].copy()
TO_MERGE_SINGLE_PROTEIN = to_merge_unique_cpds(FILTERED_ASSAYS_SINGLE_PROTEIN, keys_single_protein, ASSAY_TO_COMPOUNDS)

# Filtering only activity type - unit pairs relevant for merging
TO_MERGE_ORGANISM = TO_MERGE_ORGANISM[(TO_MERGE_ORGANISM['n_cpds_union'] > 1000) & (TO_MERGE_ORGANISM['n_assays'] > 1)].reset_index(drop=True)
TO_MERGE_SINGLE_PROTEIN = TO_MERGE_SINGLE_PROTEIN[(TO_MERGE_SINGLE_PROTEIN['n_cpds_union'] > 1000) & 
                                                  (TO_MERGE_SINGLE_PROTEIN['n_assays'] > 1) &
                                                  (TO_MERGE_SINGLE_PROTEIN['target_chembl_id'].isna() == False)].reset_index(drop=True)
TO_MERGE_ORGANISM['name'] = [f"M_O{r}" for r in range(len(TO_MERGE_ORGANISM))]
TO_MERGE_SINGLE_PROTEIN['name'] = [f"M_SP{r}" for r in range(len(TO_MERGE_SINGLE_PROTEIN))]

MERGED_COMPOUNDS = []
MERGED_LM = []
DATA = {"ORGANISM": TO_MERGE_ORGANISM, "SINGLE PROTEIN": TO_MERGE_SINGLE_PROTEIN}

for target_type in DATA:

    print(target_type)

    # Copy df
    data_target_type = DATA[target_type].copy()

    # Iterate over activity_type, unit
    for merging in data_target_type.itertuples():

        # Get data
        activity_type = merging.activity_type
        unit = merging.unit
        direction = float(merging.direction)
        assay_type = merging.assay_type
        target_type_curated_extra = merging.target_type_curated_extra
        bao_label = merging.bao_label
        strain = merging.strain
        target_chembl_id = get_target_chembl_id(merging, target_type)
        name = merging.name
        assay_keys = merging.assay_keys
        n_assays = merging.n_assays
        n_cpds_union = merging.n_cpds_union

        # Filter master table
        if target_type == 'ORGANISM':
            df = get_filtered_assay_master_organism(FILTERED_ASSAYS_ORGANISM, activity_type, unit, direction, assay_type, target_type_curated_extra, bao_label, strain)
        elif target_type == 'SINGLE PROTEIN':
            df = get_filtered_assay_master_single_protein(FILTERED_ASSAYS_SINGLE_PROTEIN, activity_type, unit, direction, assay_type, target_type_curated_extra, bao_label, strain, target_chembl_id)

        # Get quantitative and qualitative
        df_quant = df[(df['dataset_type'] == 'quantitative') | (df['dataset_type'] == 'mixed')].reset_index(drop=True)
        df_qual = df[(df['dataset_type'] == 'qualitative') | (df['dataset_type'] == 'mixed')].reset_index(drop=True)

        if len(df_quant) > 0:

            # QUANTITATIVE
            # For each expert cut-off
            for expert_cutoff in EXPERT_CUTOFFS[(activity_type, unit, target_type_curated_extra, pathogen_code)]:
                
                # Concatenate all files/data together
                assays = df_quant['assay_id'].tolist()
                files = [f"{i}_{activity_type}_{unit}_qt_{expert_cutoff}.csv.gz" for i in assays]
                data = [dfs_qt[f].assign(assay_id=a) for a, f in zip(assays, files)]
                data = pd.concat(data, ignore_index=True)
                if direction == -1:
                    data = data.sort_values("value", ascending=True).drop_duplicates("compound_chembl_id", keep="first").reset_index(drop=True)
                else:
                    data = data.sort_values("value", ascending=False).drop_duplicates("compound_chembl_id", keep="first").reset_index(drop=True)
                
                # Prepare matrices for training
                X = np.array(data['compound_chembl_id'].map(ecfps).to_list())
                Y = np.array(data['bin'].tolist())
                positives = sum(Y)

                if positives > 50:

                    print(f"Merging ... Activity type: {activity_type}, Unit: {unit}, Cutoff: {expert_cutoff}, Strain {strain}, Target ChEMBL ID ({target_chembl_id})")
                    print(f"\tTarget ChEMBL ID for ORGANISM assays is set to nan for simplicity")
                    print(f"\tCompounds: {len(X)}", f"Positives: {positives} ({round(100 * positives / len(Y), 1)}%)")

                    if positives / len(Y) > 0.5:

                        print(f"\tRatio too high: Adding random compounds from ChEMBL as decoys")
                        DECOYS = int(positives / RATIO - (len(Y) - 1))
                        print(f"\t{DECOYS} added decoys")
                        rng = random.Random(42)
                        DECOYS = rng.sample(list(DECOYS_CHEMBL), DECOYS)
                        X_decoys = np.array([ecfps[i] for i in DECOYS])
                        X = np.vstack([X, X_decoys])
                        Y = np.concatenate([Y, np.zeros(len(X_decoys), dtype=Y.dtype)])
                        positives = sum(Y)
                        print(f"\tCompounds: {len(X)}", f"Positives: {positives} ({round(100 * positives / len(Y),3)}%)")
                        decoy_df = pd.DataFrame({"compound_chembl_id": DECOYS, "bin": 0})
                        data = pd.concat([data, decoy_df], ignore_index=True)

                    # 4Fold Cros Validation
                    average_auroc, stds = KFoldTrain(X, Y, n_splits=4, n_estimators=100)
                    print(f"\tMean AUROC: {average_auroc} ± {stds}")
                    MERGED_LM.append([name, activity_type, unit, expert_cutoff, direction, assay_type, target_type_curated_extra, bao_label, strain, 
                                      target_chembl_id, n_assays, n_cpds_union, positives, round(positives/len(Y), 3), average_auroc, stds, assay_keys])
                    if average_auroc > 0.7:
                        MERGED_COMPOUNDS.extend(data['compound_chembl_id'].tolist())
                        RF = TrainRF(X, Y, n_estimators=100)
                        y_prob_ref = RF.predict_proba(X_REF)[:, 1]
                        filename = f'{name}_ref_probs.npz' 
                        np.savez_compressed(os.path.join(PATH_TO_CORRELATIONS, "M", filename), y_prob_ref=y_prob_ref)
                        # Save dataset
                        outdir = os.path.join(OUTPUT, pathogen_code, "datasets", "M")
                        os.makedirs(outdir, exist_ok=True)
                        data.to_csv(os.path.join(outdir, filename.replace("_ref_probs.npz", ".csv.gz")), index=False, compression="gzip")

                else:
                    print(f"Too few positive compounds for {activity_type}_{unit}_qt_{expert_cutoff}.. {strain} ... ({positives})")

        elif len(df_qual) > 0:

            # QUALITATIVE
            # Concatenate all files/data together
            assays = df_qual['assay_id'].tolist()
            files = [f"{i}_{activity_type}_{unit}_ql.csv.gz" for i in assays]
            data = [dfs_ql[f].assign(assay_id=a) for a, f in zip(assays, files)]
            data = pd.concat(data, ignore_index=True)
            if direction == -1:
                data = data.sort_values("value", ascending=True).drop_duplicates("compound_chembl_id", keep="first").reset_index(drop=True)
            else:
                data = data.sort_values("value", ascending=False).drop_duplicates("compound_chembl_id", keep="first").reset_index(drop=True)
            
            # Prepare matrices for training
            X = np.array(data['compound_chembl_id'].map(ecfps).to_list())
            Y = np.array(data['bin'].tolist())
            positives = sum(Y)

            if positives > 50:

                print(f"Merging ... Activity type: {activity_type}, Unit: {unit}, Cutoff: {expert_cutoff}, Strain {strain}, Target ChEMBL ID ({target_chembl_id})")
                print(f"\tTarget ChEMBL ID for ORGANISM assays is set to nan for simplicity")
                print(f"\tCompounds: {len(X)}", f"Positives: {positives} ({round(100 * positives / len(Y), 1)}%)")

                if positives / len(Y) > 0.5:

                    print(f"\tRatio too high: Adding random compounds from ChEMBL as decoys")
                    DECOYS = int(positives / RATIO - (len(Y) - 1))
                    print(f"\t{DECOYS} added decoys")
                    rng = random.Random(42)
                    DECOYS = rng.sample(list(DECOYS_CHEMBL), DECOYS)
                    X_decoys = np.array([ecfps[i] for i in DECOYS])
                    X = np.vstack([X, X_decoys])
                    Y = np.concatenate([Y, np.zeros(len(X_decoys), dtype=Y.dtype)])
                    positives = sum(Y)
                    print(f"\tCompounds: {len(X)}", f"Positives: {positives} ({round(100 * positives / len(Y),3)}%)")
                    decoy_df = pd.DataFrame({"compound_chembl_id": DECOYS, "bin": 0})
                    data = pd.concat([data, decoy_df], ignore_index=True)

                # 4Fold Cros Validation
                average_auroc, stds = KFoldTrain(X, Y, n_splits=4, n_estimators=100)
                print(f"\tMean AUROC: {average_auroc} ± {stds}")
                # In qualitative datasets, expert cutoff is nan
                MERGED_LM.append([name, activity_type, unit, np.nan, direction, assay_type, target_type_curated_extra, bao_label, strain, 
                                      target_chembl_id, n_assays, n_cpds_union, positives, round(positives/len(Y), 3), average_auroc, stds, assay_keys])
                if average_auroc > 0.7:
                    MERGED_COMPOUNDS.extend(data['compound_chembl_id'].tolist())
                    RF = TrainRF(X, Y, n_estimators=100)
                    y_prob_ref = RF.predict_proba(X_REF)[:, 1]
                    os.makedirs(os.path.join(PATH_TO_CORRELATIONS, "M"), exist_ok=True)
                    filename = f'{name}_ref_probs.npz' 
                    np.savez_compressed(os.path.join(PATH_TO_CORRELATIONS, "M", filename), y_prob_ref=y_prob_ref)
                    # Save dataset
                    outdir = os.path.join(OUTPUT, pathogen_code, "datasets", "M")
                    os.makedirs(outdir, exist_ok=True)
                    data.to_csv(os.path.join(outdir, filename.replace("_ref_probs.npz", ".csv.gz")), index=False, compression="gzip")

            else:
                print(f"Too few positive compounds for {activity_type}_{unit}_ql .. {strain} ... ({positives})")


        else:
            raise TypeError("Please revise. df_quant and df_qual are empty...")
        

MERGED_LM = pd.DataFrame(MERGED_LM, columns=["name", "activity_type", "unit", "expert_cutoff", "direction", "assay_type", "target_type_curated_extra", "bao_label", 
                                             "strain", "target_chembl_id", "n_assays", "n_cpds_union", "positives", "ratio", "average_auroc", "stds", "assay_keys"])
MERGED_LM.to_csv(os.path.join(OUTPUT, pathogen_code, "merged_LM.csv"), index=False)