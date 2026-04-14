"""
Shared utilities for the per-pathogen pipeline (steps 07+).

Functions
---------
load_pathogen(pathogen_code)
    Resolve a pathogen code to its full name from config/pathogens.csv.

load_manual_assays(pathogen_code)
    Load the list of manually curated assay ChEMBL IDs for a pathogen.

load_assay_metadata()
    Load lookup dicts for assay source labels and BAO ontology labels.

build_assays_info(df, pathogen_chemical_space, assay_to_src_id,
                  assay_to_bao_format, src_id_to_src_short_name,
                  bao_id_to_label, directions=None)
    Build the per-(assay_id, activity_type, unit) summary table used
    by both assays_raw.csv (step 07) and assays_cleaned.csv (step 08).

load_expert_cutoffs(CONFIGPATH)
    Load expert binarization cutoffs from config/expert_cutoffs.csv.

extra_curation_target_type(target_type, target_type_curated)
    Enforce simplified target type rules (SINGLE PROTEIN / ORGANISM / DISCARDED).

add_target_type_curated(assays_cleaned, parameters)
    Merge curated target type from step 09 parameters onto cleaned assays.
"""

from collections import defaultdict, Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import re
import os

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, root)
from default import DATAPATH, CONFIGPATH


def harmonize(x):
    """Normalize a standard_type string to a canonical activity type.

    Strips underscores, spaces, dots, slashes and uppercases the result,
    collapsing variants such as 'IC 50', 'ic_50', 'IC/50' → 'IC50'.
    Used in steps 01 and 05.
    """
    return re.sub(r"[_\s./\\]", "", str(x).upper().strip())


def load_pathogen(pathogen_code):
    """Return the full pathogen name for a given code.

    Raises SystemExit with a clear message if the code is not found.
    """
    df = pd.read_csv(os.path.join(CONFIGPATH, "pathogens.csv"))
    row = df.loc[df["code"].eq(pathogen_code)]
    if row.empty:
        raise SystemExit(f"Unknown pathogen code: '{pathogen_code}'. "
                         f"Valid codes: {sorted(df['code'].tolist())}")
    return row.iloc[0]["pathogen"]


def load_manual_assays(pathogen_code):
    """Return the set of manually curated assay ChEMBL IDs for a pathogen.

    Reads config/assays/<pathogen_code>.csv, which may contain assay IDs
    as comma- or newline-separated values.
    Returns an empty set if the file does not exist.
    """
    path = os.path.join(CONFIGPATH, "assays", f"{pathogen_code}.csv")
    if not os.path.exists(path):
        return set()
    manual_assays = open(os.path.join(path), "r").read()
    ids = set([j for i in manual_assays.split("\n") for j in i.split(",")])
    return ids


def load_assay_metadata():
    """Load lookup dicts for assay source labels and BAO ontology labels.

    Returns
    -------
    assay_to_src_id : dict  {assay_chembl_id -> src_id}
    assay_to_bao_format : dict  {assay_chembl_id -> bao_format}
    src_id_to_src_short_name : dict  {src_id -> short_name}
    bao_id_to_label : dict  {bao_id -> label}
    """
    assay_data = pd.read_csv(
        os.path.join(DATAPATH, "chembl_activities", "assays.csv"), low_memory=False
    )[["chembl_id", "src_id", "bao_format"]]
    assay_to_src_id = dict(zip(assay_data["chembl_id"], assay_data["src_id"]))
    assay_to_bao_format = dict(zip(assay_data["chembl_id"], assay_data["bao_format"]))

    source = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "source.csv"))
    src_id_to_src_short_name = dict(zip(source["src_id"], source["src_short_name"]))

    bao = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "bioassay_ontology.csv"))
    bao_id_to_label = dict(zip(bao["bao_id"], bao["label"]))

    return assay_to_src_id, assay_to_bao_format, src_id_to_src_short_name, bao_id_to_label


def _only_one(values, name):
    """Assert exactly one unique value exists and return it."""
    if len(values) != 1:
        raise ValueError(f"Expected exactly one {name}, found: {values}")
    return values[0]


def build_assays_info(df, pathogen_chemical_space,
                      assay_to_src_id, assay_to_bao_format,
                      src_id_to_src_short_name, bao_id_to_label,
                      directions=None):
    """Build a per-(assay_id, activity_type, unit) summary DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Activity records for a single pathogen (raw or cleaned).
    pathogen_chemical_space : set
        Set of compound_chembl_ids representing the pathogen chemical space.
    assay_to_src_id, assay_to_bao_format, src_id_to_src_short_name,
    bao_id_to_label : dict
        Metadata lookup dicts from load_assay_metadata().
    directions : dict or None
        If provided, {(activity_type, unit): direction} mapping used to add
        a 'direction' column to the output (step 08). Omitted in step 07.

    Returns
    -------
    pd.DataFrame sorted by compound count descending.
    """
    assays = sorted(set(df["assay_chembl_id"]))

    assay_to_idx = defaultdict(list)
    for i, assay_id in enumerate(df["assay_chembl_id"].to_numpy()):
        assay_to_idx[assay_id].append(i)

    columns = [
        "assay_id", "assay_type", "assay_organism", "assay_tax_id", "assay_strain", "assay_confidence_score",
        "doc_chembl_id",
        "target_type", "target_organism", "target_tax_id", "target_chembl_id", "tid",
        "bao_label", "source_label",
        "activity_type", "unit", "activities", "nan_values",
        "cpds", "act_flag", "inact_flag", "frac_cs",
    ]
    if directions is not None:
        columns.append("direction")

    rows = []
    for assay in tqdm(assays):
        df_ = df.iloc[assay_to_idx[assay]]

        assay_type = _only_one(list(set(df_["assay_type"])), "assay_type")
        target_type = _only_one(list(set(df_["target_type"])), "target_type")
        target_chembl_id = _only_one(list(set(df_["target_chembl_id"])), "target_chembl_id")
        target_organism = _only_one(list(set(df_["target_organism"])), "target_organism")
        assay_organism = _only_one(list(set(df_["assay_organism"])), "assay_organism")
        doc_chembl_id = _only_one(list(set(df_["doc_chembl_id"])), "doc_chembl_id")
        tid           = _only_one(df_["tid"].unique().tolist(),                    "tid")
        assay_tax_id  = _only_one(df_["assay_tax_id"].unique().tolist(),           "assay_tax_id")
        assay_strain  = _only_one(df_["assay_strain"].unique().tolist(),           "assay_strain")
        assay_conf    = _only_one(df_["assay_confidence_score"].unique().tolist(), "assay_confidence_score")
        target_tax_id = _only_one(df_["target_tax_id"].unique().tolist(),          "target_tax_id")

        bao_label = bao_id_to_label.get(assay_to_bao_format.get(assay), np.nan)
        source_label = src_id_to_src_short_name.get(assay_to_src_id.get(assay), np.nan)

        for act_type in sorted(set(df_["activity_type"])):
            df__ = df_[df_["activity_type"] == act_type]
            activity_type = _only_one(list(set(df__["activity_type"])), "activity_type")

            for u in sorted(set(df__["unit"]), key=lambda x: (pd.isna(x), x)):
                df___ = df__[df__["unit"].isna()] if pd.isna(u) else df__[df__["unit"] == u]
                unit = _only_one(list(set(df___["unit"])), "unit")

                activities = len(df___)
                cpds = set(df___["compound_chembl_id"])
                nan_values = int(df___["value"].isna().sum())
                text_flag = Counter(df___["text_flag"])
                act_flag = text_flag[1]
                inact_flag = text_flag[-1]
                frac_cs = round(len(cpds & pathogen_chemical_space) / len(pathogen_chemical_space), 5)

                row = [
                    assay, assay_type, assay_organism, assay_tax_id, assay_strain, assay_conf,
                    doc_chembl_id,
                    target_type, target_organism, target_tax_id, target_chembl_id, tid,
                    bao_label, source_label,
                    activity_type, unit, activities, nan_values,
                    len(cpds), act_flag, inact_flag, frac_cs,
                ]
                if directions is not None:
                    row.append(directions.get((act_type, unit), np.nan))

                rows.append(row)

    result = pd.DataFrame(rows, columns=columns)
    return result.sort_values("cpds", ascending=False).reset_index(drop=True)


def load_expert_cutoffs(CONFIGPATH):
    """
    Load expert cutoffs from config/expert_cutoffs.csv.

    Returns a dictionary mapping
        (activity_type, unit, target_type, pathogen_code) -> list of float cutoffs.

    Parameters
    ----------
    CONFIGPATH : str
        Path to the config folder.
    """
    cutoffs_df = pd.read_csv(os.path.join(CONFIGPATH, "expert_cutoffs.csv"))
    return {
        (a, b, c, d): [float(k) for k in e.split(";")]
        for a, b, c, d, e in cutoffs_df[
            ["activity_type", "unit", "target_type", "pathogen_code", "expert_cutoff"]
        ].values
    }


def extra_curation_target_type(target_type, target_type_curated):
    """
    Enforce simplified curation rules for ChEMBL assay target types.

    Collapses the LLM-curated target type to one of:
    SINGLE PROTEIN, ORGANISM, or DISCARDED.

    Rules
    -----
    - UNCHECKED     → allow ORGANISM / SINGLE PROTEIN / DISCARDED; else DISCARDED
    - NON-MOLECULAR → allow ORGANISM / DISCARDED; else DISCARDED
    - SINGLE PROTEIN / PROTEIN COMPLEX / PROTEIN FAMILY → SINGLE PROTEIN
    - ORGANISM → ORGANISM
    - anything else → DISCARDED

    Parameters
    ----------
    target_type : str
        Raw ChEMBL target type.
    target_type_curated : str
        LLM- or human-proposed curated target type (from step 09).

    Returns
    -------
    str
        "SINGLE PROTEIN", "ORGANISM", or "DISCARDED".
    """
    if not isinstance(target_type_curated, str):
        return 'DISCARDED'

    target_type = target_type.strip().upper()
    target_type_curated = target_type_curated.strip().upper()

    if target_type == 'UNCHECKED':
        return target_type_curated if target_type_curated in {'ORGANISM', 'SINGLE PROTEIN', 'DISCARDED'} else 'DISCARDED'

    elif target_type == 'NON-MOLECULAR':
        return target_type_curated if target_type_curated in {'ORGANISM', 'DISCARDED'} else 'DISCARDED'

    elif target_type in {'SINGLE PROTEIN', 'PROTEIN COMPLEX', 'PROTEIN FAMILY'}:
        return 'SINGLE PROTEIN'

    elif target_type == 'ORGANISM':
        return 'ORGANISM'

    else:
        return 'DISCARDED'


def add_target_type_curated(assays_cleaned, parameters):
    """
    Merge `target_type_curated` from the step 09 parameters file onto assays_cleaned.

    Matches on (assay_id, activity_type, unit) and enforces a bijective join —
    raises ValueError if either table has keys not present in the other.

    Parameters
    ----------
    assays_cleaned : pandas.DataFrame
        Output of step 08; must contain assay_id, activity_type, unit.
    parameters : pandas.DataFrame
        Output of step 09; must contain assay_id, activity_type, unit,
        target_type_curated.

    Returns
    -------
    pandas.DataFrame
        assays_cleaned with an added `target_type_curated` column.
    """
    parameters = parameters[["assay_id", "activity_type", "unit", "target_type_curated"]].copy()

    # Step 09 writes missing units as empty strings
    assays_cleaned = assays_cleaned.copy()
    assays_cleaned["unit"] = assays_cleaned["unit"].fillna("")
    parameters["unit"] = parameters["unit"].fillna("")

    keys = ["assay_id", "activity_type", "unit"]

    # Step 09 may write duplicate rows (e.g. from resume runs); keep the last
    n_before = len(parameters)
    parameters = parameters.drop_duplicates(subset=keys, keep="last").reset_index(drop=True)
    if len(parameters) < n_before:
        print(f"  Warning: removed {n_before - len(parameters)} duplicate rows from parameters (step 09 resume artefact)")

    left_only = parameters[keys].merge(assays_cleaned[keys], on=keys, how="left", indicator=True)
    if not left_only["_merge"].eq("both").all():
        raise ValueError("parameters contains keys not present in assays_cleaned")

    right_only = assays_cleaned[keys].merge(parameters[keys], on=keys, how="left", indicator=True)
    if not right_only["_merge"].eq("both").all():
        raise ValueError("assays_cleaned contains keys not present in parameters")

    assays_cleaned = assays_cleaned.merge(parameters, on=keys, how="left", validate="1:1")
    assays_cleaned["unit"] = assays_cleaned["unit"].replace("", np.nan)
    return assays_cleaned
