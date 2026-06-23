"""
Shared utilities for the per-pathogen pipeline (steps 07+).

Functions
---------
load_pathogen(pathogen_code)
    Resolve a pathogen code to its full name from config/pathogens.csv.

load_manual_assays(pathogen_code)
    Load the list of manually curated assay ChEMBL IDs for a pathogen.

load_excluded_assays(pathogen_code)
    Load the set of assay ChEMBL IDs to exclude for a pathogen.

load_pubchem_assays(pathogen_code)
    Return the set of ChEMBL assay IDs that have a PubChem counterpart.

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
from default import DATAPATH, CONFIGPATH, STRAIN_CATALOG_PREFIXES


def harmonize(x):
    """Normalize a standard_type string to a canonical activity type.

    Strips underscores, spaces, dots, slashes and uppercases the result,
    collapsing variants such as 'IC 50', 'ic_50', 'IC/50' → 'IC50'.
    Used in steps 01 and 05.
    """
    return re.sub(r"[_\s./\\]", "", str(x).upper().strip())


# Catalog acronyms whose designator is numeric (the number alone identifies the
# strain) vs. WHO/FDA whose designator may be a letter ("WHO X"). They are matched
# differently so a leading acronym is only stripped when a real designator follows
# (avoids mangling a strain name that merely starts with the same letters).
_STRAIN_PREFIXES_NUMERIC = sorted(
    [p for p in STRAIN_CATALOG_PREFIXES if p not in ("WHO", "FDA")], key=len, reverse=True
)
_STRAIN_PREFIXES_ALNUM = sorted(
    [p for p in STRAIN_CATALOG_PREFIXES if p in ("WHO", "FDA")], key=len, reverse=True
)
_STRAIN_SEP = r"[\s\-_./\\]"
_STRAIN_PREFIX_RE = re.compile(
    r"\b(?:"
    # numeric-catalog acronyms: strip only when a number follows, optionally after a
    # short letter sub-series (e.g. "ATCC BAA-1605" -> "BAA-1605", "ATCC 25922" -> "25922")
    r"(?:" + "|".join(_STRAIN_PREFIXES_NUMERIC) + r")(?=" + _STRAIN_SEP + r"*(?:[A-Z]+" + _STRAIN_SEP + r"*)?\d)"
    # WHO/FDA: strip when any alphanumeric designator follows ("WHO X" -> "X")
    r"|(?:" + "|".join(_STRAIN_PREFIXES_ALNUM) + r")(?=" + _STRAIN_SEP + r"*[A-Z0-9])"
    r")" + _STRAIN_SEP + r"*",
    re.IGNORECASE,
)


def normalize_strain(x):
    """Canonical merge key for a curated strain name (step 15).

    Uppercases, strips leading culture-collection catalog acronyms
    (``STRAIN_CATALOG_PREFIXES`` — ATCC, DSM, NCTC, ... plus WHO/FDA), and removes
    spaces, hyphens, dots, underscores and slashes, so formatting and catalog-prefix
    variants collapse to a single key::

        'BAA-1605', 'BAA 1605', 'ATCC BAA1605'  -> 'BAA1605'
        'ATCC 25922', 'NCTC 25922', '25922'      -> '25922'
        'WHO X', 'WHO-X'                         -> 'X'

    Numeric-catalog acronyms are only stripped when a digit follows, and WHO/FDA
    only when an alphanumeric designator follows, so a strain name that merely
    begins with those letters is left intact. Returns '' for empty/NaN input.
    """
    if x is None:
        return ""
    s = str(x).upper().strip()
    if s == "" or s == "NAN":
        return ""
    stripped = _STRAIN_PREFIX_RE.sub("", s)
    key = re.sub(_STRAIN_SEP, "", stripped)
    if key == "":  # value was nothing but a catalog acronym — keep it as-is
        key = re.sub(_STRAIN_SEP, "", s)
    return key


def prefix_key(x):
    """Prefix-PRESERVING key: uppercase, strip, remove separators but KEEP the
    leading collection acronym. Unlike `normalize_strain`, this does not fold the
    catalog prefix, so it can distinguish the same number across collections —
    `LMG 17978` → `LMG17978` vs `ATCC 17978` → `ATCC17978`. Used only by the curated
    `prefixed` split rules, to pull a genuinely-different same-number strain out of
    its number group before the prefix is stripped. Returns '' for empty/NaN.
    """
    if x is None:
        return ""
    s = str(x).upper().strip()
    if s == "" or s == "NAN":
        return ""
    return re.sub(_STRAIN_SEP, "", s)


def load_strain_equivalences(pathogen_code, config_path=CONFIGPATH):
    """Load the curated strain alias map for a pathogen, used by `strain_norm_key`.

    Reads `config/strains/strain_equivalences_<pathogen>.csv` (built by
    `tools/build_strain_equivalences.py` from BacDive + a manual seed). Returns a
    `(exact, contains, prefixed)` tuple of rules, keyed by `member_key`:
      - `exact`    : normalized base key → canonical label (the common case).
      - `contains` : list of (substring, label); a base key *containing* the
                     substring maps to the label (e.g. any "...BCG..." → "BCG").
      - `prefixed` : prefix-preserving key → label; matches the original value
                     *with* its collection acronym, to SPLIT same-number strains
                     that differ by collection (e.g. `LMG 17978` ≠ `ATCC 17978`).
    Returns None when no file exists (clean no-op for uncurated pathogens).
    """
    path = os.path.join(config_path, "strains", f"strain_equivalences_{pathogen_code}.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    exact, contains, prefixed = {}, [], {}
    for _, r in df.iterrows():
        label = r["canonical_label"]
        mt = str(r.get("match_type", "exact")).strip().lower()
        if mt == "contains":
            contains.append((str(r["member_key"]), label))
        elif mt == "prefixed":
            prefixed[str(r["member_key"])] = label
        else:
            exact[str(r["member_key"])] = label
    return exact, contains, prefixed


def strain_norm_key(strain_curated, atcc_id, equivalences=None):
    """Vectorized strain merge key used in steps 15/18/21.

    Layered resolution, in order:
      1. `normalize_strain(assay_strain_curated)` — formatting/catalog-prefix variants.
      2. fall back to `normalize_strain(atcc_id)` when the curated name is empty
         (closes the "split-column" gap where identity lives only in the ATCC field).
      3. if `equivalences` is given, resolve to a canonical strain label:
         `prefixed` SPLIT rules first (matched on the original value WITH its
         collection acronym, so e.g. `LMG 17978` is pulled out before its prefix is
         stripped), then `exact` (cross-collection / name↔number, e.g. ATCC 25922 =
         DSM 1103; H37Rv = ATCC 27294), then `contains` (substring, e.g. BCG).

    Accepts pandas Series, returns a Series. With `equivalences=None` it behaves
    exactly as before (layers 1–2 only).
    """
    cur = strain_curated.fillna("").astype(str)
    atc = atcc_id.fillna("").astype(str)
    norm = cur.map(normalize_strain)
    base = norm.where(norm != "", atc.map(normalize_strain))
    if not equivalences:
        return base
    exact, contains, prefixed = equivalences
    eff = cur.where(norm != "", atc)  # the original value that produced the base key

    def resolve(eff_val, base_key):
        if prefixed:
            pk = prefix_key(eff_val)
            if pk in prefixed:
                return prefixed[pk]
        if base_key in exact:
            return exact[base_key]
        for sub, label in contains:
            if sub and sub in base_key:
                return label
        return base_key

    return pd.Series([resolve(e, b) for e, b in zip(eff, base)], index=base.index)


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

    Reads config/assays_to_include/<pathogen_code>.csv, which may contain assay IDs
    as comma- or newline-separated values.
    Returns an empty set if the file does not exist.
    """
    path = os.path.join(CONFIGPATH, "assays_to_include", f"{pathogen_code}.csv")
    if not os.path.exists(path):
        return set()
    manual_assays = open(os.path.join(path), "r").read()
    ids = set([j for i in manual_assays.split("\n") for j in i.split(",")])
    return ids


def load_excluded_assays(pathogen_code):
    """Return the set of assay ChEMBL IDs to exclude for a pathogen.

    Reads config/assays_to_exclude/<pathogen_code>.csv (one ID per line).
    Returns an empty set if the file does not exist.
    """
    path = os.path.join(CONFIGPATH, "assays_to_exclude", f"{pathogen_code}.csv")
    if not os.path.exists(path):
        return set()
    return set(line.strip() for line in open(path) if line.strip())


def load_pubchem_assays(pathogen_code):
    """Return the set of ChEMBL assay IDs that have a PubChem counterpart.

    Reads config/pubchem_aids/chembl_assays_in_pubchem_<pathogen_code>.csv
    and returns the values in the 'Source ID' column.
    Returns an empty set if the file does not exist or has no data rows.
    """
    path = os.path.join(CONFIGPATH, "pubchem_aids", f"chembl_assays_in_pubchem_{pathogen_code}.csv")
    if not os.path.exists(path):
        return set()
    df = pd.read_csv(path)
    if "Source ID" not in df.columns or df.empty:
        return set()
    return set(df["Source ID"].dropna())


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
        (activity_type, unit, target_type, pathogen_code) -> list of float cutoffs,
    ordered so that index 0 is the most lenient threshold and index -1 is the most
    stringent. For direction == -1 (e.g. MIC, IC50), the ascending numeric list is
    reversed; for direction == 1 (e.g. % inhibition), it is kept as-is.

    Parameters
    ----------
    CONFIGPATH : str
        Path to the config folder.
    """
    cutoffs_df = pd.read_csv(os.path.join(CONFIGPATH, "expert_cutoffs.csv"))
    result = {}
    for _, row in cutoffs_df.iterrows():
        key = (row["activity_type"], row["unit"], row["target_type"], row["pathogen_code"])
        cutoffs = [float(k) for k in str(row["expert_cutoff"]).split(";")]
        if row["direction"] == -1.0:
            cutoffs = cutoffs[::-1]
        result[key] = cutoffs
    return result


def reference_cutoff(cutoff_list, unit):
    """Return the preferred (reference) cutoff used by dataset selection.

    Percent endpoints (unit == "%") prefer the most lenient cutoff (index 0, e.g. 50):
    50% effect remains the default "active" call even though it is now the lowest of the
    three. All other endpoints use the positional middle cutoff (index 1). Returns None
    for an empty list.
    """
    if not cutoff_list:
        return None
    idx = 0 if unit == "%" else 1
    return cutoff_list[idx] if len(cutoff_list) > idx else cutoff_list[-1]


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
