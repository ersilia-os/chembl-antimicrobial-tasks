"""
Dataset preparation utilities for step 12.

Functions
---------
get_assay_data(chembl_pathogen, assay_chembl_id, activity_type, unit, cols)
    Extract rows for a single (assay_id, activity_type, unit) triplet.

get_cut_value(assay_data, direction)
    Return the extreme value used to replace censored measurements.

count_relations(assay_data)
    Count "=", "<", ">" relation operators.

get_assay_data_quantitative(assay_data)
    Filter to rows with a non-null numeric value.

get_assay_data_qualitative(assay_data)
    Build a compound-level binary dataset from text flags.

adjust_relation(assay_data, direction, cut)
    Correct censored values on the wrong side of the biological direction.

disambiguate_compounds(assay_data, direction)
    Keep the single most active measurement per compound.

binarize_with_expert_cutoff(assay_data_quantitative, expert_cutoff, direction)
    Apply an expert cutoff to produce a binary label column.

set_variables_quantitative(assay_data_quantitative)
    Compute summary statistics for a binarized quantitative dataset.

set_variables_qualitative(assay_data_qualitative)
    Compute summary statistics for a qualitative dataset.

get_activity_stats_quantitative(activities_quantitative)
    Compute percentile statistics for a list of activity values.

make_dataset_filename(assay_id, activity_type, unit, dataset_type, expert_cutoff)
    Build the canonical filename for a dataset file produced by step 12.

zip_and_remove(datasets_dir)
    Archive _qt_, _ql, and _mx_ dataset files into zip archives and remove originals.
"""

from collections import Counter
import pandas as pd
import numpy as np
import zipfile
import os


def get_assay_data(chembl_pathogen, assay_chembl_id, activity_type, unit, cols):
    """
    Extract assay activity data for a given assay_chembl_id, activity_type, and unit.

    If `unit` is a string, the function filters rows where `unit` matches exactly.
    Otherwise, it filters rows where `unit` is missing (NaN).

    Parameters
    ----------
    chembl_pathogen : pandas.DataFrame
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
    if isinstance(unit, str):
        assay_data = chembl_pathogen[
            (chembl_pathogen['assay_chembl_id'] == assay_chembl_id) &
            (chembl_pathogen['activity_type'] == activity_type) &
            (chembl_pathogen['unit'] == unit)
        ].reset_index(drop=True)[cols]
    else:
        assay_data = chembl_pathogen[
            (chembl_pathogen['assay_chembl_id'] == assay_chembl_id) &
            (chembl_pathogen['activity_type'] == activity_type) &
            (chembl_pathogen['unit'].isna())
        ].reset_index(drop=True)[cols]

    return assay_data


def get_cut_value(assay_data, direction):
    """
    Get a cutoff value from assay_data to adjust relations based on direction.

    If direction == 1, returns the minimum value in assay_data['value'].
    If direction == -1, returns the maximum value in assay_data['value'].
    Otherwise, returns np.nan.

    Parameters
    ----------
    assay_data : pandas.DataFrame
        DataFrame containing a 'value' column with numeric assay values.
    direction : int
        Direction indicator: 1 → use minimum value, -1 → use maximum value.

    Returns
    -------
    float
        Cutoff value computed from the 'value' column or np.nan.
    """
    if direction == 1:
        return min(assay_data['value'])
    elif direction == -1:
        return max(assay_data['value'])
    else:
        return np.nan


def count_relations(assay_data):
    """
    Count relation operators in assay_data['relation'].

    Parameters
    ----------
    assay_data : pandas.DataFrame
        DataFrame containing a 'relation' column.

    Returns
    -------
    tuple
        (equal, lower, higher) counts corresponding to "=", "<", ">".
    """
    counter_relations = Counter(assay_data['relation'].tolist())
    return counter_relations["="], counter_relations["<"], counter_relations[">"]


def get_assay_data_quantitative(assay_data):
    """
    Return only rows in assay_data with non-missing quantitative values.

    Parameters
    ----------
    assay_data : pandas.DataFrame
        DataFrame containing a 'value' column.

    Returns
    -------
    pandas.DataFrame
        Filtered dataframe containing only rows with non-null `value`,
        with the `text_flag` column dropped.
    """
    assay_data_quantitative = assay_data[assay_data['value'].notna()].reset_index(drop=True)
    return assay_data_quantitative.drop(columns=['text_flag'])


def get_assay_data_qualitative(assay_data):
    """
    Build a compound-level qualitative dataset from assay text flags.

    Aggregates `text_flag` values per compound, checks for conflicting
    qualitative labels, assigns a final compound-level label, and
    returns one row per compound with a binary activity label.

    Parameters
    ----------
    assay_data : pandas.DataFrame
        DataFrame containing at least `compound_chembl_id` and `text_flag`.

    Returns
    -------
    pandas.DataFrame
        Qualitative dataset with one row per compound and a binary label.

    Raises
    ------
    ValueError
        If the same compound has both active (1) and inactive (-1) labels.
    """
    assay_data_qualitative = assay_data.copy()

    # Aggregate to compound-level label
    compound_labels = assay_data_qualitative.groupby("compound_chembl_id")["text_flag"].apply(set)
    compound_ids = compound_labels.index.tolist()
    label_sets = compound_labels.tolist()

    # Detect compound-level conflicts — drop and warn rather than crash
    compound_conflict = [((1 in s) and (-1 in s)) for s in label_sets]
    if any(compound_conflict):
        bad = [cid for cid, c in zip(compound_ids, compound_conflict) if c]
        print(f"  Warning: {len(bad)} compound(s) with conflicting active/inactive labels dropped: {bad[:5]}{'...' if len(bad) > 5 else ''}")
        conflict_set = set(bad)
        assay_data_qualitative = assay_data_qualitative[~assay_data_qualitative["compound_chembl_id"].isin(conflict_set)].reset_index(drop=True)
        compound_labels = assay_data_qualitative.groupby("compound_chembl_id")["text_flag"].apply(set)
        compound_ids = compound_labels.index.tolist()
        label_sets = compound_labels.tolist()

    # Final compound label: 1 > -1 > 0
    compound_final = dict(zip(compound_ids, [1 if 1 in s else (-1 if -1 in s else 0) for s in label_sets]))
    assay_data_qualitative["qualitative_label"] = assay_data_qualitative["compound_chembl_id"].map(compound_final)

    # One row per compound, remove neutral (0)
    assay_data_qualitative = assay_data_qualitative.drop_duplicates(subset=["compound_chembl_id"]).reset_index(drop=True)
    assay_data_qualitative = assay_data_qualitative[assay_data_qualitative["qualitative_label"] != 0].reset_index(drop=True)

    # Binary label
    assay_data_qualitative["bin"] = [0 if x == -1 else 1 for x in assay_data_qualitative["qualitative_label"].tolist()]

    cols = ["compound_chembl_id", "smiles", "activity_type", "unit", "text_flag", "qualitative_label", "bin"]
    return assay_data_qualitative[cols]


def adjust_relation(assay_data, direction, cut):
    """
    Adjust relations in an assay DataFrame according to the biological direction.

    Censored measurements on the wrong side of the direction are replaced with
    the extreme value `cut` and their relation set to "=".

    Parameters
    ----------
    assay_data : pd.DataFrame
        Must contain the columns 'relation' and 'value'.
    direction : int
        +1 → higher = more active; -1 → lower = more active.
    cut : float
        Extreme value (min or max) used to replace wrong-side censored measurements.

    Returns
    -------
    pd.DataFrame
        Copy of assay_data with adjusted relation and value.
    """
    df = assay_data.copy()
    rel = df["relation"].astype(str)

    if direction == +1:
        # Higher = more active: ">" is fine, "<" is wrong-side → replace
        df.loc[rel == ">", "relation"] = "="
        df.loc[rel == "<", "relation"] = "="
        df.loc[rel == "<", "value"] = cut

    elif direction == -1:
        # Lower = more active: "<" is fine, ">" is wrong-side → replace
        df.loc[rel == "<", "relation"] = "="
        df.loc[rel == ">", "relation"] = "="
        df.loc[rel == ">", "value"] = cut

    else:
        raise ValueError(f"Invalid direction={direction}. Expected +1 or -1.")

    return df


def disambiguate_compounds(assay_data, direction):
    """
    Select a single measurement per compound according to the biological direction.

    Parameters
    ----------
    assay_data : pd.DataFrame
        Must contain 'compound_chembl_id' and 'value'.
        Assumes all relations have already been adjusted.
    direction : int
        +1 → keep maximum value; -1 → keep minimum value.

    Returns
    -------
    pd.DataFrame
        One row per compound, keeping the most active measurement.
    """
    if direction not in [1, -1]:
        raise ValueError("direction must be +1 (higher = more active) or -1 (lower = more active).")

    ascending = direction == -1  # min for -1, max for +1
    return (assay_data
            .sort_values("value", ascending=ascending)
            .drop_duplicates(subset="compound_chembl_id", keep="first")
            .reset_index(drop=True))


def binarize_with_expert_cutoff(assay_data_quantitative, expert_cutoff, direction):
    """
    Binarize quantitative assay values using an expert-defined cutoff.

    Parameters
    ----------
    assay_data_quantitative : pandas.DataFrame
        Must contain a numeric `value` column.
    expert_cutoff : float
        Expert-defined activity threshold.
    direction : int
        +1 → values >= cutoff are active; -1 → values <= cutoff are active.

    Returns
    -------
    pandas.DataFrame
        Input dataframe with an added binary `bin` column.
    """
    if direction == +1:
        assay_data_quantitative["bin"] = (assay_data_quantitative["value"] >= expert_cutoff).astype(int)
    else:
        assay_data_quantitative["bin"] = (assay_data_quantitative["value"] <= expert_cutoff).astype(int)

    return assay_data_quantitative


def binarize_conservative(assay_data_quantitative, expert_cutoff, direction):
    """Qualifier-aware, conservative binarisation ("prove-the-hit" rule).

    Unlike `binarize_with_expert_cutoff`, this respects the censoring qualifier
    in the `relation` column and only labels a compound ACTIVE (1) when the
    measurement *proves* it lies on the active side of the cutoff. Ambiguous or
    wrong-side censored records (whose true value could sit either side of the
    cutoff) are labelled INACTIVE (0) rather than guessed. Rows are kept, never
    dropped. This strictly shrinks the active set relative to a bare threshold.

    Value-less rows (`value` is NaN) are the recovered qualitative-inactive
    molecules step 22 carries into the dataset: they have no measured value to
    prove a hit, so they are always labelled INACTIVE (0). This is enforced
    explicitly via `value.notna()` rather than relying on NaN comparison
    semantics.

    Only '=', '>' and '<' relations occur in the curated data.

    Parameters
    ----------
    assay_data_quantitative : pandas.DataFrame
        Must contain numeric `value` and a `relation` column ('=', '>' or '<').
    expert_cutoff : float
        Activity threshold.
    direction : int
        +1 -> higher = more active (active iff proven value >= cutoff);
        -1 -> lower = more active (active iff proven value <= cutoff).

    Returns
    -------
    pandas.DataFrame
        Input dataframe with an added binary `bin` column.
    """
    value = assay_data_quantitative["value"]
    relation = assay_data_quantitative["relation"].astype(str)
    has_value = value.notna()  # value-less qualitative-inactive rows are never active
    if direction == +1:
        # proven active only if '=' or '>' AND value already at/above the cutoff
        active = has_value & relation.isin(["=", ">"]) & (value >= expert_cutoff)
    else:
        # proven active only if '=' or '<' AND value already at/below the cutoff
        active = has_value & relation.isin(["=", "<"]) & (value <= expert_cutoff)
    assay_data_quantitative["bin"] = active.astype(int)
    return assay_data_quantitative


def set_variables_quantitative(assay_data_quantitative):
    """
    Compute basic statistics for a quantitative, binarized assay dataset.

    Parameters
    ----------
    assay_data_quantitative : pandas.DataFrame
        One row per compound with a binary `bin` column and a `value` column.

    Returns
    -------
    tuple
        (positives, ratio, n_compounds, activities_list)
    """
    positives = (assay_data_quantitative["bin"] == 1).sum()
    n_compounds = len(set(assay_data_quantitative['compound_chembl_id']))
    activities = assay_data_quantitative['value'].tolist()
    assert n_compounds == len(activities)
    return positives, round(positives / len(assay_data_quantitative), 5), n_compounds, activities


def set_variables_qualitative(assay_data_qualitative):
    """
    Compute basic statistics for a qualitative assay dataset.

    Parameters
    ----------
    assay_data_qualitative : pandas.DataFrame
        One row per compound with a binary `bin` column.

    Returns
    -------
    tuple
        (positives, ratio, n_compounds)
    """
    positives = (assay_data_qualitative["bin"] == 1).sum()
    n_compounds = len(set(assay_data_qualitative['compound_chembl_id']))
    assert n_compounds == len(assay_data_qualitative)
    return positives, round(positives / len(assay_data_qualitative), 5), n_compounds


def get_activity_stats_quantitative(activities_quantitative):
    """
    Compute summary statistics for a list of quantitative activity values.

    Returns min, p1, p25, p50, p75, p99, max, rounded to 3 decimals.

    Parameters
    ----------
    activities_quantitative : array-like

    Returns
    -------
    tuple
        (min_, p1, p25, p50, p75, p99, max_)
    """
    return (
        round(np.min(activities_quantitative), 3),
        round(np.percentile(activities_quantitative, 1), 3),
        round(np.percentile(activities_quantitative, 25), 3),
        round(np.percentile(activities_quantitative, 50), 3),
        round(np.percentile(activities_quantitative, 75), 3),
        round(np.percentile(activities_quantitative, 99), 3),
        round(np.max(activities_quantitative), 3),
    )


def make_dataset_filename(assay_id, activity_type, unit, dataset_type, expert_cutoff=None):
    """
    Build the filename for a dataset file produced by step 12.

    Parameters
    ----------
    assay_id : str
    activity_type : str
    unit : str or None
        Accepted for call-site compatibility but not included in the filename.
    dataset_type : str
        One of 'quantitative', 'qualitative', 'mixed'.
    expert_cutoff : float or None
        Required for quantitative and mixed types.

    Returns
    -------
    str
        e.g. 'CHEMBL12345_IC50_qt_1.0.csv.gz'
    """
    base = f"{assay_id}_{activity_type}"
    if dataset_type == "quantitative":
        return f"{base}_qt_{expert_cutoff}.csv.gz"
    elif dataset_type == "qualitative":
        return f"{base}_ql.csv.gz"
    elif dataset_type == "mixed":
        return f"{base}_mx_{expert_cutoff}.csv.gz"
    else:
        raise ValueError(f"Unknown dataset_type: {dataset_type!r}")


def zip_and_remove(datasets_dir):
    """
    Archive dataset files into three zip archives and remove the originals.

    Creates:
      - datasets_qt.zip  (files containing "_qt_")
      - datasets_ql.zip  (files ending with "_ql.csv.gz")
      - datasets_mx.zip  (files containing "_mx_")

    Parameters
    ----------
    datasets_dir : str
        Directory containing the dataset files.

    Returns
    -------
    tuple[int, int, int]
        (n_qt_files, n_ql_files, n_mx_files)
    """
    qt_files = [os.path.join(datasets_dir, f) for f in os.listdir(datasets_dir) if "_qt_" in f]
    ql_files = [os.path.join(datasets_dir, f) for f in os.listdir(datasets_dir) if f.endswith("_ql.csv.gz")]
    mx_files = [os.path.join(datasets_dir, f) for f in os.listdir(datasets_dir) if "_mx_" in f]

    for zip_path, files in [
        (os.path.join(datasets_dir, "datasets_qt.zip"), qt_files),
        (os.path.join(datasets_dir, "datasets_ql.zip"), ql_files),
        (os.path.join(datasets_dir, "datasets_mx.zip"), mx_files),
    ]:
        if os.path.exists(zip_path):
            os.remove(zip_path)
        with zipfile.ZipFile(zip_path, mode="w", compression=zipfile.ZIP_DEFLATED) as z:
            for fp in files:
                z.write(fp, arcname=os.path.basename(fp))

    for fp in qt_files + ql_files + mx_files:
        os.remove(fp)

    return len(qt_files), len(ql_files), len(mx_files)


# ---------------------------------------------------------------------------
# Curation helpers (steps 22–24) — merged from the former curation_v2_utils.
# The v2 unit of work is the (assay_id, activity_type) pair; molecules are
# compared by InChIKey. Activity categories are assigned from an explicit,
# user-reviewed inclusion list keyed by (activity_type, unit):
#   DR = dose-response / concentration endpoints (e.g. MIC, IC50; umol.L-1);
#   SP = single-point readouts (% and zone-of-inhibition mm);
#   anything not in config/activity_type_categories.csv -> UNCATEGORIZED.
# ---------------------------------------------------------------------------

# Filename-safe slugs for units (binarisation output, step 23).
UNIT_SLUG = {"umol.L-1": "uM", "%": "pct", "mm": "mm"}


def load_activity_categories(config_path):
    """Load the (activity_type, unit) -> category inclusion list (step 22).

    Reads config/activity_type_categories.csv. Returns a dict keyed by
    (activity_type, unit) -> 'DR' or 'SP'. Combinations absent from the file are
    treated as UNCATEGORIZED by `assign_activity_category`.
    """
    df = pd.read_csv(os.path.join(config_path, "activity_type_categories.csv"))
    return {(str(a), str(u)): str(c)
            for a, u, c in zip(df["activity_type"], df["unit"], df["category"])}


def assign_activity_category(activity_type, unit, categories):
    """Return the activity category ('DR', 'SP', or 'UNCATEGORIZED') for an
    (activity_type, unit) pair, using the inclusion map from
    `load_activity_categories`."""
    return categories.get((str(activity_type), str(unit)), "UNCATEGORIZED")


def per_pair_molecule_flags(acts, key):
    """Per (pair, inchikey) presence flags used by the step-21 curation.

    `acts` must carry the `key` columns plus 'inchikey', 'value' and 'text_flag'.
    Rows with a null InChIKey are ignored. Returns a DataFrame with columns
    `key + ['inchikey', 'has_num', 'has_active', 'has_inactive']`:
        has_num       the molecule has >=1 numeric record in the pair
        has_active    it has >=1 text active flag (text_flag == 1)
        has_inactive  it has >=1 text inactive flag (text_flag == -1)
    Shared by `build()` (recovered-inactive recovery) and the attrition path so
    both derive the same per-pair molecule membership from one implementation.
    """
    a = acts[acts["inchikey"].notna()].copy()
    a["_num"] = a["value"].notna()
    a["_act"] = a["text_flag"] == 1
    a["_ina"] = a["text_flag"] == -1
    return (a.groupby(key + ["inchikey"])
            .agg(has_num=("_num", "any"), has_active=("_act", "any"),
                 has_inactive=("_ina", "any"))
            .reset_index())


def recovered_inactive_mask(molflags):
    """Boolean mask over `per_pair_molecule_flags` rows selecting recovered
    qualitative-inactive molecules: no numeric record, a purely-inactive text
    flag, and no active flag (curator rule, 2026-07-01)."""
    return ((~molflags["has_num"]) & molflags["has_inactive"]
            & (~molflags["has_active"]))


def dataset_inchikey_sets(acts, key):
    """Per-pair set of the InChIKeys that actually land in the written dataset =
    numeric molecules UNION recovered qualitative inactives.

    Returns a dict {pair_tuple: set(inchikey)}. This is the single canonical
    "dataset membership" used for every molecule count/union in the attrition
    figures, so plots match the datasets step 21 writes (not the broader set of
    all tested molecules, nor the numeric-only subset).
    """
    mf = per_pair_molecule_flags(acts, key)
    keep = mf[mf["has_num"] | (mf["has_inactive"] & ~mf["has_active"])]
    return {k: set(g["inchikey"]) for k, g in keep.groupby(key)}


def dedup_most_active_by_inchikey(df, direction, preserve_qualifier=False):
    """Keep one row per InChIKey, the single most active measurement.

    Mirrors the step-12 intra-assay disambiguation but matches molecules on
    `inchikey` instead of `compound_chembl_id`. Censored measurements on the wrong
    side of the biological direction are first reset to the group extreme via
    `adjust_relation` (needs only the direction, not any cutoff).

    Parameters
    ----------
    df : pandas.DataFrame
        Numeric records for one (assay_id, activity_type) pair. Must contain
        'inchikey', 'value' and 'relation'.
    direction : int
        +1 -> higher = more active (keep maximum); -1 -> lower = more active
        (keep minimum). Any other value: keep the first record per InChIKey.
    preserve_qualifier : bool, default False
        If False (legacy behaviour), the kept row carries the `adjust_relation`
        output: wrong-side censored values are replaced by the group extreme and
        every relation becomes "=". If True, selection of the most-active record
        is unchanged, but the kept row's ORIGINAL `value` and `relation` are
        restored so downstream binarisation can reason about the qualifier
        (a '>' survives only when a molecule has no better record, i.e. a pure
        non-inhibitor). Used by step 22 to enable conservative binarisation.
    """
    df = df.copy()
    if preserve_qualifier:
        df["_orig_value"] = df["value"]
        df["_orig_relation"] = df["relation"].astype(str)

    if direction in (1, -1):
        cut = get_cut_value(df, direction)
        df = adjust_relation(df, direction, cut)
        ascending = direction == -1  # minimum for -1, maximum for +1
        # mergesort is stable -> deterministic tie-breaking (reproducible)
        df = df.sort_values("value", ascending=ascending, kind="mergesort")
    # else: undefined direction (should not occur for numeric data); keep first.

    kept = df.drop_duplicates(subset="inchikey", keep="first").reset_index(drop=True)
    if preserve_qualifier:
        kept["value"] = kept["_orig_value"]
        kept["relation"] = kept["_orig_relation"]
        kept = kept.drop(columns=["_orig_value", "_orig_relation"])
    return kept


def load_category_cutoffs(config_path):
    """Load the single fixed binarisation cutoff per unit (step 23).

    Reads config/category_cutoffs.csv (columns: category, unit, cutoff) and returns
    a `{unit: cutoff}` dict. The direction (which side is active) is taken from each
    dataset's own `direction` column, not from here.
    """
    df = pd.read_csv(os.path.join(config_path, "category_cutoffs.csv"))
    return {str(u): float(c) for u, c in zip(df["unit"], df["cutoff"])}
