"""
Modeling utilities for the light modeling pipeline (steps 13+).

Functions
---------
get_reference_set_compounds(compounds)
    Select up to 10,000 compound IDs as a pathogen reference set.

load_ecfp_all(h5_path)
    Load all Morgan fingerprints from an HDF5 file into memory.

load_ecfp_subset_by_chembl_id(h5_path, chembl_id_set)
    Load a subset of Morgan fingerprints by ChEMBL ID.

load_data_from_zip(zip_path, filename)
    Load a gzipped CSV from inside a zip archive.

load_all_gz_csvs_from_zip(zip_path)
    Load all *.csv.gz members from a ZIP archive into a dict of DataFrames.

KFoldTrain(X, Y, n_splits, n_estimators, random_state)
    Stratified K-fold cross-validation with Random Forest; returns mean AUROC ± std.

TrainRF(X, Y, n_estimators)
    Train a Random Forest on all data and return the fitted model.
"""

from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import numpy as np
import zipfile
import gzip
import h5py
import pandas as pd


def get_reference_set_compounds(compounds):
    """
    Return a reference set of ChEMBL compound IDs from a compounds DataFrame.

    If the DataFrame has more than 10,000 rows, returns the first 5,000 and
    last 5,000 compound_chembl_id values. Otherwise returns all.

    Parameters
    ----------
    compounds : pandas.DataFrame
        Must contain a 'compound_chembl_id' column.

    Returns
    -------
    list or pandas.Series
        Compound ChEMBL IDs for the reference set.
    """
    if len(compounds) > 10000:
        return compounds['compound_chembl_id'][:5000].tolist() + compounds['compound_chembl_id'][-5000:].tolist()
    return compounds['compound_chembl_id']


def load_ecfp_all(h5_path):
    """
    Load all Morgan fingerprints from an HDF5 file into memory.

    Parameters
    ----------
    h5_path : str
        Path to the HDF5 file containing datasets 'SMILES' and 'X_morgan'.

    Returns
    -------
    dict[str, np.ndarray]
        Mapping {chembl_id: fingerprint (shape (nBits,))}.
    """
    with h5py.File(h5_path, "r") as f:
        ids = f["SMILES"][:, 3].astype(str)
        fps = f["X_morgan"][:]
    return {cid: fp for cid, fp in zip(ids, fps)}


def load_ecfp_subset_by_chembl_id(h5_path, chembl_id_set):
    """
    Load a subset of Morgan fingerprints by ChEMBL ID.

    Parameters
    ----------
    h5_path : str
        Path to the HDF5 file containing datasets 'SMILES' and 'X_morgan'.
    chembl_id_set : iterable[str]
        ChEMBL IDs to load. IDs not found in the file are silently ignored.

    Returns
    -------
    dict[str, np.ndarray]
        Mapping {chembl_id: fingerprint (shape (nBits,))}.
    """
    chembl_id_set = set(chembl_id_set)
    with h5py.File(h5_path, "r") as f:
        ids = f["SMILES"][:, 3].astype(str)
        idx = np.flatnonzero(np.isin(ids, list(chembl_id_set)))
        fps = f["X_morgan"][idx]
    return {ids[i]: fp for i, fp in zip(idx, fps)}


def load_data_from_zip(zip_path, filename):
    """
    Load a gzipped CSV file from inside a ZIP archive into a DataFrame.

    Parameters
    ----------
    zip_path : str
        Path to the ZIP archive.
    filename : str
        Name of the gzipped CSV file inside the ZIP.

    Returns
    -------
    pandas.DataFrame
    """
    with zipfile.ZipFile(zip_path) as z:
        with z.open(filename) as raw:
            with gzip.open(raw, mode="rt") as f:
                return pd.read_csv(f)


def load_all_gz_csvs_from_zip(zip_path):
    """Load all *.csv.gz members from a ZIP archive into DataFrames.

    Parameters
    ----------
    zip_path : str
        Path to the ZIP archive.

    Returns
    -------
    dict[str, pandas.DataFrame]
        Mapping of member filename -> loaded DataFrame.
    """
    dfs = {}
    with zipfile.ZipFile(zip_path, "r") as z:
        for name in z.namelist():
            if name.endswith(".csv.gz"):
                with z.open(name) as f:
                    dfs[name] = pd.read_csv(f, compression="gzip")
    return dfs


def KFoldTrain(X, Y, n_splits=4, n_estimators=100, random_state=42):
    """
    Stratified K-fold cross-validation with Random Forest; returns mean AUROC ± std.

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
        RNG seed for both fold shuffling and the classifier.

    Returns
    -------
    tuple[float, float]
        (mean_auroc, std_auroc) rounded to 3 decimals.
    """
    def _make_rf():
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
        rf = _make_rf()
        rf.fit(X[train_idx], Y[train_idx])
        y_prob = rf.predict_proba(X[test_idx])[:, 1]
        aurocs.append(roc_auc_score(Y[test_idx], y_prob))

    return round(float(np.mean(aurocs)), 3), round(float(np.std(aurocs)), 3)


def TrainRF(X, Y, n_estimators=100):
    """
    Train a RandomForestClassifier on all provided data and return the fitted model.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix (n_samples, n_features).
    Y : np.ndarray
        Binary labels (n_samples,).
    n_estimators : int
        Number of trees.

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
        random_state=42,
    )
    rf.fit(X, Y)
    return rf
