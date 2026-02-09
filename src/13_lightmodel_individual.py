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

def get_reference_set_compounds(compounds):
    """
    Return a reference set of ChEMBL compound IDs from a compounds DataFrame.

    If the DataFrame has more than 10,000 rows, this returns a reduced reference
    set consisting of the first 5,000 and the last 5,000 `compound_chembl_id`
    values (as a list). Otherwise, it returns the full `compound_chembl_id`
    column.
    """
    if len(compounds) > 10000:
        return compounds['compound_chembl_id'][:5000].tolist() + compounds['compound_chembl_id'][-5000:].tolist()
    else:
        return compounds['compound_chembl_id']
    
def load_ecfp_subset_by_chembl_id(h5_path, chembl_id_set):
    """Load a subset of ECFP (Morgan count) fingerprints by ChEMBL ID.

    Parameters
    ----------
    h5_path : str
        Path to the HDF5 file containing datasets "SMILES" and "X_morgan".
    chembl_id_set : set[str] | iterable[str]
        ChEMBL IDs to keep.

    Returns
    -------
    dict[str, np.ndarray]
        Mapping {chembl_id: fingerprint (shape (nBits,))} for IDs present in the file.
        IDs in `chembl_id_set` that are not found are silently ignored.
    """
    chembl_id_set = set(chembl_id_set)
    with h5py.File(h5_path, "r") as f:
        ids = f["SMILES"][:, 3].astype(str)
        idx = np.flatnonzero(np.isin(ids, list(chembl_id_set)))
        fps = f["X_morgan"][idx]
    return {ids[i]: fp for i, fp in zip(idx, fps)}

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
        random_state=42
    )
    rf.fit(X, Y)
    return rf

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

print("Step 13")

# Define output directory
OUTPUT = os.path.join(root, "..", "output")

def condition_A(df):
    return (
        df["dataset_type"].isin(["quantitative", "mixed"])
        & (df["cpds_qt"] >= 1000)
        & (df["pos_qt"] >= 50)
        & (df["ratio_qt"].between(0.001, 0.5, inclusive="both")))

# def condition_B(df):
#     return (
#         df["dataset_type"].isin(["qualitative", "mixed"])
#         & (df["cpds_ql"] >= 1000)
#         & (df["pos_ql"] >= 50)
#         & (df["ratio_ql"].between(0.001, 0.5, inclusive="both")))

def condition_B(df):
    return (
        df["dataset_type"].isin(["quantitative", "mixed"])
        & (df["pos_qt"] >= 100)
        & (df["ratio_qt"] >= 0.5))

# def condition_D(df):
#     return (
#         df["dataset_type"].isin(["qualitative", "mixed"])
#         & (df["pos_ql"] >= 100)
#         & (df["ratio_ql"] >= 0.5))

RATIO = 0.1

# Create path to correlations
PATH_TO_CORRELATIONS = os.path.join(OUTPUT, pathogen_code, "correlations")
os.makedirs(PATH_TO_CORRELATIONS, exist_ok=True)

# Load assay datasets
COLS = ["assay_id", "activity_type", "unit", "target_type", "target_type_curated_extra", "cpds", "direction", "dataset_type", "expert_cutoff", 
        "pos_qt", "ratio_qt", "cpds_qt", "pos_ql", "ratio_ql", "cpds_ql", "overlap_mx", "pos_mx", "ratio_mx", "cpds_mx"]
print("Collecting datasets")
ASSAYS_DATASETS = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "datasets.csv"))[COLS]

# Create reference set of compounds per pathogen
compounds = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "compound_counts.csv.gz"))
print(f"Generating a reference set of compounds for {pathogen_code}")
REFERENCE_SET = get_reference_set_compounds(compounds)
pd.DataFrame(REFERENCE_SET, columns=['reference_compounds']).to_csv(os.path.join(OUTPUT, pathogen_code, "reference_set.csv.gz"), index=False)

# Get all compounds for pathogen
compounds = set(compounds['compound_chembl_id'])
print(f"Loading ECFPs...")

# Loading Morgan fingerprints
PATH_TO_ECFPs = os.path.join(DATAPATH, "chembl_processed", "ChEMBL_ECFPs.h5")
ecfps = load_ecfp_all(PATH_TO_ECFPs)

# Get ChEMBL compounds not tested against the pathogen
DECOYS_CHEMBL = set([i for i in ecfps if i not in compounds])

# Prepare reference matrix of Morgan fingerprints
X_REF = np.array([ecfps[cid] for cid in REFERENCE_SET if cid in ecfps])

CONDITIONS = {"A": condition_A, 
              "B": condition_B}

LABELS = sorted(CONDITIONS)
LABEL_COMPOUNDS = {lab: set() for lab in LABELS}
INDIVIDUAL_LM = []

# For each label
for LABEL in LABELS:

    print(f"Creating {LABEL} datasets...")

    # Identify condition
    CONDITION_DATASETS = ASSAYS_DATASETS[CONDITIONS[LABEL](ASSAYS_DATASETS)].copy().reset_index(drop=True)
    CONDITION_DATASETS['label'] = LABEL

    # Get AUROC INFO
    AUROC_AVG, AUROC_STD = [], []

    # Iterate over assays LABEL
    for c, assay in CONDITION_DATASETS.iterrows():

        # Load varibles
        assay_id, activity_type, unit, expert_cutoff = assay.assay_id, assay.activity_type, assay.unit, assay.expert_cutoff
        dataset_type = assay.dataset_type

        # Load data
        if dataset_type == 'quantitative':
            zip_path = os.path.join(OUTPUT, pathogen_code, "datasets", "datasets_qt.zip")
            filename = "_".join([str(assay_id), str(activity_type), str(unit), "qt", f"{expert_cutoff}.csv.gz"])
        elif dataset_type == 'mixed':
            zip_path = os.path.join(OUTPUT, pathogen_code, "datasets", "datasets_mx.zip")
            filename = "_".join([str(assay_id), str(activity_type), str(unit), "mx", f"{expert_cutoff}.csv.gz"])
        df = load_data_from_zip(zip_path, filename)

        # Add compounds
        cids = set(df["compound_chembl_id"].astype(str))
        LABEL_COMPOUNDS[LABEL].update(cids)

        # Prepare matrices
        X = np.array(df['compound_chembl_id'].map(ecfps).to_list())
        Y = np.array(df['bin'].tolist())
        positives = sum(Y)

        print(f"Assay ID: {assay_id}, Activity type: {activity_type}, Unit: {unit}, Cutoff: {expert_cutoff}")
        print(f"\tCompounds: {len(X)}", f"Positives: {positives} ({round(100 * positives / len(Y),3)}%)")

        if LABEL == 'B':  # with decoys

            print(f"\tAdding random compounds from ChEMBL as decoys")
            DECOYS = int(positives / RATIO - (len(Y) - 1))
            print(f"\t{DECOYS} added decoys")
            rng = random.Random(42)
            DECOYS = rng.sample(list(DECOYS_CHEMBL), DECOYS)
            X_decoys = np.array([ecfps[i] for i in DECOYS])
            X = np.vstack([X, X_decoys])
            Y = np.concatenate([Y, np.zeros(len(X_decoys), dtype=Y.dtype)])
            print(f"\tCompounds: {len(X)}", f"Positives: {positives} ({round(100 * positives / len(Y),3)}%)")

        # 4Fold Cros Validation
        average_auroc, stds = KFoldTrain(X, Y, n_splits=4, n_estimators=100)
        print(f"\tMean AUROC: {average_auroc} ± {stds}")
        AUROC_AVG.append(average_auroc)
        AUROC_STD.append(stds)

        # Train on full data and predict on reference set
        RF = TrainRF(X, Y, n_estimators=100)
        y_prob_ref = RF.predict_proba(X_REF)[:, 1]
        os.makedirs(os.path.join(PATH_TO_CORRELATIONS, LABEL), exist_ok=True)
        np.savez_compressed(os.path.join(PATH_TO_CORRELATIONS, LABEL, filename.replace(".csv.gz", "_ref_probs.npz")), y_prob_ref=y_prob_ref)

    CONDITION_DATASETS['avg'] = AUROC_AVG
    CONDITION_DATASETS['std'] = AUROC_STD
    INDIVIDUAL_LM.append(CONDITION_DATASETS.copy())

    considered_datasets = len(CONDITION_DATASETS)
    considered_assays = len(set([tuple(i) for i in CONDITION_DATASETS[['assay_id', 'activity_type', 'unit']].values]))

    print(f"Summary for {LABEL}...")
    print(f"Number of considered datasets: {considered_datasets}")
    print(f"Number of considered assays: {considered_assays}")
    print(f"Chemical space coverage: {round(100 * len(LABEL_COMPOUNDS[LABEL]) / len(compounds), 1)}%")

# Save results
INDIVIDUAL_LM = pd.concat(INDIVIDUAL_LM, ignore_index=True)
INDIVIDUAL_LM.to_csv(os.path.join(OUTPUT, pathogen_code, 'individual_LM.csv'), index=False)

all_cpds = set([cpd for lab in LABEL_COMPOUNDS for cpd in LABEL_COMPOUNDS[lab]])
print(f"Chemical space coverage (AB): {round(100 * len(all_cpds) / len(compounds), 1)}%")