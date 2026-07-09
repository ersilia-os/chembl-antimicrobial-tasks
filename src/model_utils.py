"""
Modeling utilities for the light modeling pipeline (steps 13+).

Functions
---------
load_ecfp_all(h5_path)
    Load all Morgan fingerprints from an HDF5 file into memory.

load_data_from_zip(zip_path, filename)
    Load a gzipped CSV from inside a zip archive.

load_all_gz_csvs_from_zip(zip_path)
    Load all *.csv.gz members from a ZIP archive into a dict of DataFrames.

add_decoys(X, Y, decoys_pool, ecfps, decoy_ratio, random_state)
    Add random decoy compounds to reach a target active fraction.

downsample_negatives(X, Y, target_ratio, random_state)
    Downsample negatives (when actives < 5%) to reach a target active fraction.

KFoldTrain(X, Y, n_splits, n_estimators, random_state)
    Stratified K-fold cross-validation with Random Forest; returns mean AUROC ± std.

TrainRF(X, Y, n_estimators)
    Train a Random Forest on all data and return the fitted model.

cv_oof_predict(X, Y, n_splits, n_estimators, random_state)
    Stratified K-fold out-of-fold probability predictions (one per sample).

youden_threshold(y_true, y_prob)
    Youden-index-optimal probability cutoff from a single ROC.

balanced_agreement(y_true, y_pred)
    Mean recall over the classes present in y_true (single-class safe).
"""

import os
import sys

from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
import numpy as np
import zipfile
import random
import gzip
import h5py
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, root)
from default import RANDOM_SEED


def load_ecfp_all(h5_path):
    """
    Load all Morgan fingerprints from an HDF5 file into memory.

    Parameters
    ----------
    h5_path : str
        Path to the HDF5 file containing datasets 'smiles' and 'X_morgan'.

    Returns
    -------
    dict[str, np.ndarray]
        Mapping {chembl_id: fingerprint (shape (nBits,))}.
    """
    with h5py.File(h5_path, "r") as f:
        ids = f["smiles"][:, 3].astype(str)
        fps = f["X_morgan"][:]
    return {cid: fp for cid, fp in zip(ids, fps)}


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


def add_decoys(X, Y, decoys_pool, ecfps, decoy_ratio, random_state=RANDOM_SEED):
    """Add random decoy compounds to reduce the active ratio to `decoy_ratio`.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix.
    Y : np.ndarray
        Binary labels.
    decoys_pool : iterable
        Pool of compound IDs to sample decoys from.
    ecfps : dict
        Mapping {compound_id: fingerprint}.
    decoy_ratio : float
        Target active fraction after adding decoys (e.g. 0.1 means 10% actives).
    random_state : int
        RNG seed for reproducibility.

    Returns
    -------
    X : np.ndarray
    Y : np.ndarray
    decoy_ids : list
        IDs of the sampled decoy compounds.
    """
    n_positives = int(Y.sum())
    n_decoys = max(0, int(n_positives / decoy_ratio) - len(Y))
    rng = random.Random(random_state)
    decoy_ids = rng.sample(list(decoys_pool), n_decoys)
    X_decoys = np.array([ecfps[i] for i in decoy_ids])
    X = np.vstack([X, X_decoys])
    Y = np.concatenate([Y, np.zeros(len(X_decoys), dtype=Y.dtype)])
    return X, Y, decoy_ids


def downsample_negatives(X, Y, target_ratio=0.10, random_state=RANDOM_SEED):
    """Downsample negatives when the active ratio is below 5% to reach `target_ratio`.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix.
    Y : np.ndarray
        Binary labels.
    target_ratio : float
        Desired active fraction after downsampling.
    random_state : int
        RNG seed for reproducibility.

    Returns
    -------
    X : np.ndarray
    Y : np.ndarray
    """
    n_positives = int(Y.sum())
    active_ratio = n_positives / len(Y)
    if active_ratio >= 0.05:
        return X, Y
    n_neg_target = int(n_positives / target_ratio) - n_positives
    neg_idx = np.where(Y == 0)[0]
    rng_ds = np.random.RandomState(random_state)
    sampled_neg = rng_ds.choice(neg_idx, size=min(n_neg_target, len(neg_idx)), replace=False)
    keep = np.sort(np.concatenate([np.where(Y == 1)[0], sampled_neg]))
    return X[keep], Y[keep]


def KFoldTrain(X, Y, n_splits=4, n_estimators=100, random_state=RANDOM_SEED):
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


def cv_oof_predict(X, Y, n_splits=5, n_estimators=100, random_state=RANDOM_SEED):
    """Stratified K-fold out-of-fold (OOF) probabilities with Random Forest.

    Each sample is predicted exactly once, as a test example, by a forest trained on the other
    folds (5-fold => 80/20 train/test per fold). Returns an array of P(active) aligned to `Y`,
    suitable for a single pooled ROC / Youden cutoff and for actives-vs-inactives distribution plots.

    Parameters
    ----------
    X : np.ndarray (n_samples, n_features)
    Y : np.ndarray (n_samples,)  binary labels
    n_splits, n_estimators, random_state : CV / RF settings (same RF config as KFoldTrain).

    Returns
    -------
    np.ndarray (n_samples,)  out-of-fold P(active) per sample.
    """
    Y = np.asarray(Y).astype(int)
    oof = np.full(len(Y), np.nan, dtype=float)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    for train_idx, test_idx in skf.split(X, Y):
        rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=None, min_samples_split=2,
                                    min_samples_leaf=1, max_features="sqrt", n_jobs=8,
                                    random_state=random_state)
        rf.fit(X[train_idx], Y[train_idx])
        oof[test_idx] = rf.predict_proba(X[test_idx])[:, 1]
    return oof


def youden_threshold(y_true, y_prob):
    """Youden-index-optimal probability cutoff from a single ROC.

    J = sensitivity + specificity - 1 = tpr - fpr, maximised over ROC thresholds. Returns the
    probability cutoff and the sensitivity/specificity/J achieved there.

    Returns
    -------
    tuple[float, float, float, float]  (cutoff, sensitivity, specificity, J)
    """
    fpr, tpr, thr = roc_curve(y_true, y_prob)
    j = tpr - fpr
    # roc_curve prepends a sentinel threshold at index 0 (np.inf / max(score)+1, tpr=fpr=0,
    # J=0); skip it so a poor model can't return an out-of-range cutoff. Clamp to [0,1].
    k = 1 + int(np.argmax(j[1:])) if len(thr) > 1 else 0
    cutoff = float(min(max(thr[k], 0.0), 1.0))
    return cutoff, float(tpr[k]), float(1.0 - fpr[k]), float(j[k])


def balanced_agreement(y_true, y_pred):
    """Mean recall over the classes PRESENT in y_true (single-class safe).

    Both classes present -> balanced accuracy; all-inactive set -> specificity; all-active -> sensitivity.
    Used to score how well a pool model (thresholded at its cutoff) calls a leftover dataset's actives
    active AND its inactives inactive, without relying on AUROC (undefined for single-class sets).
    """
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    recalls = [(y_pred[y_true == c] == c).mean() for c in (0, 1) if (y_true == c).any()]
    return float(np.mean(recalls)) if recalls else float("nan")


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
        random_state=RANDOM_SEED,
    )
    rf.fit(X, Y)
    return rf
