import pandas as pd
import numpy as np
import sys
import os
import re

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH


# --- Helper functions ---

def calculate_pchembl(uM):
    """Convert a concentration in uM to a pChEMBL value (clipped to [1, 9])."""
    try:
        value = uM * 1e-6
        if value == 0:
            value = 1e-10
        return np.clip(-np.log10(value), 1, 9)
    except:
        return np.nan


def create_text_flag(df):
    """Combine activity_comment and standard_text binary columns into a single text_flag.

    Returns 1 (active), -1 (inactive), or 0 (unknown/no signal).
    Raises if a row has conflicting signals (both 1 and -1).
    Drops the source columns after merging.
    """
    cond_pos = (df["activity_comment"] == 1) | (df["standard_text"] == 1)
    cond_neg = (df["activity_comment"] == -1) | (df["standard_text"] == -1)
    cond_nan = (df["activity_comment"] == 0) & (df["standard_text"] == 0)

    conflict = cond_pos & cond_neg
    if conflict.any():
        raise ValueError(
            "Conflicting labels (contains both 1 and -1):\n"
            + df.loc[conflict, ["compound_chembl_id", "activity_comment", "standard_text"]].head(20).to_string()
        )

    df["text_flag"] = pd.Series(pd.NA, index=df.index, dtype="Int64")
    df.loc[cond_pos, "text_flag"] = 1
    df.loc[cond_neg, "text_flag"] = -1
    df.loc[cond_nan, "text_flag"] = 0

    return df.drop(columns=["activity_comment", "standard_text"])


# --- Main script ---

print("Step 05")

# 1. Load data and drop rows with no SMILES
print("Loading data...")
activities_all_raw = pd.read_csv(
    os.path.join(DATAPATH, "chembl_processed", "04_activities_all_raw.csv"),
    low_memory=False
)
nans = activities_all_raw["smiles"].isna().sum()
print(f"Removing {nans} activities with no associated SMILES...")
activities_all_raw = activities_all_raw[activities_all_raw["smiles"].notna()].reset_index(drop=True)

# 2. Flag activity comments (+1 active, -1 inactive, 0 unknown)
print("Flagging activity comments...")
activity_comments = pd.read_csv(os.path.join(CONFIGPATH, "activity_comments_manual_curation.csv"), low_memory=False)
comments_act = set(activity_comments[activity_comments["manual_curation_activity"] == 1]["activity_comment"])
comments_inact = set(activity_comments[activity_comments["manual_curation_activity"] == -1]["activity_comment"])

def flag_comment(val):
    val = str(val).strip().lower()
    if val == "nan":
        return 0
    if val in comments_act:
        return 1
    if val in comments_inact:
        return -1
    return 0

activities_all_raw["activity_comment_bin"] = activities_all_raw["activity_comment"].map(flag_comment)

# 3. Flag standard text values (+1 active, -1 inactive, 0 unknown)
print("Flagging standard text values...")
standard_text = pd.read_csv(os.path.join(CONFIGPATH, "standard_text_manual_curation.csv"), low_memory=False)
text_act = set(standard_text[standard_text["manual_curation_activity"] == 1]["standard_text_value"])
text_inact = set(standard_text[standard_text["manual_curation_activity"] == -1]["standard_text_value"])

def flag_text(val):
    val = str(val)
    if val == "nan":
        return 0
    if val in text_act:
        return 1
    if val in text_inact:
        return -1
    return 0

activities_all_raw["standard_text_bin"] = activities_all_raw["standard_text_value"].map(flag_text)

# 4. Convert units and values using ucum_manual.csv
print("Standardizing units and converting values...")
ucum = pd.read_csv(os.path.join(CONFIGPATH, "ucum_manual.csv"))
unit_to_final = dict(zip(ucum["units"], ucum["final_unit"]))
unit_to_formula = dict(zip(ucum["units"], ucum["transformer"]))

new_values, new_units = [], []
for _, mw, std_value, std_unit in activities_all_raw[["activity_id", "mw", "standard_value", "standard_units"]].itertuples(index=False):
    final_unit = unit_to_final.get(std_unit, np.nan)
    new_units.append(final_unit)

    if pd.notna(std_value):
        formula = unit_to_formula.get(std_unit, np.nan)
        if pd.notna(formula):
            new_values.append(eval(formula, {"standard_value": std_value, "molecular_weight": mw}))
        else:
            new_values.append(std_value)
    else:
        new_values.append(np.nan)

activities_all_raw["converted_values"] = new_values
activities_all_raw["converted_units"] = new_units

# 5. Standardize relations (collapse >=, >>, ~ etc. to =, >, <)
print("Standardizing relations...")
RELATIONS = {"=": "=", ">": ">", "<": "<", ">>": ">", ">=": ">", "<<": "<", "<=": "<", np.nan: "=", "~": "="}
activities_all_raw["relation"] = activities_all_raw["standard_relation"].map(RELATIONS)
print(f"  Old: {activities_all_raw['standard_relation'].value_counts(dropna=False).to_dict()}")
print(f"  New: {activities_all_raw['relation'].value_counts(dropna=False).to_dict()}")

# 6. Calculate pChEMBL from converted uM values
print("Calculating pChEMBL values...")
activities_all_raw["pchembl_calculated"] = activities_all_raw.apply(
    lambda r: calculate_pchembl(r["converted_values"])
    if pd.notna(r["converted_values"]) and r["converted_units"] == "umol.L-1"
    else np.nan,
    axis=1
)

# 7. Convert doc_id to doc_chembl_id
print("Converting doc IDs...")
docs = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "docs.csv"), low_memory=False)
doc_id_to_chembl_id = dict(zip(docs["doc_id"], docs["chembl_id"]))
print(f"  {len(docs)} docs, {len(doc_id_to_chembl_id)} mappings")
del docs
activities_all_raw["doc_chembl_id"] = activities_all_raw["doc_id"].map(doc_id_to_chembl_id)

# 8. Harmonize standard_type -> harmonized_type
def harmonize(x):
    return re.sub(r"[_\s./\\]", "", str(x).upper().strip())

activities_all_raw["harmonized_type"] = activities_all_raw["standard_type"].map(harmonize)

# 9. Rename and drop raw columns
activities_all_raw = activities_all_raw.drop(columns=[
    "standard_relation", "standard_value", "standard_units", "standard_type",
    "activity_comment", "standard_text_value", "doc_id"
])
activities_all_raw = activities_all_raw.rename(columns={
    "activity_comment_bin": "activity_comment",
    "standard_text_bin": "standard_text",
    "converted_values": "value",
    "converted_units": "unit",
    "harmonized_type": "activity_type",
    "pchembl_value": "pchembl",
})

# 10. Map activity type synonyms to canonical names
print("Mapping activity type synonyms...")
synonyms = pd.read_csv(os.path.join(CONFIGPATH, "synonyms.csv"))
for canonical, syns in zip(synonyms["activity_type"], synonyms["synonyms"]):
    for syn in syns.split(";"):
        activities_all_raw.loc[activities_all_raw["activity_type"] == syn.strip(), "activity_type"] = canonical

# 10. Merge comment and text flags into a single text_flag column
activities_all_raw = create_text_flag(activities_all_raw)

# 11. Save diagnostic summary of text flags per (activity_type, unit)
flagged = activities_all_raw["text_flag"].isin([1, -1])
text_flag_summary = (
    activities_all_raw[["activity_type", "unit", "text_flag"]]
    .assign(flagged=flagged.to_numpy())
    .groupby(["activity_type", "unit"], dropna=False)
    .agg(count=("flagged", "size"), comments=("flagged", "sum"))
    .reset_index()
    .sort_values("count", ascending=False, ignore_index=True)
)
text_flag_summary.to_csv(
    os.path.join(DATAPATH, "chembl_processed", "05_activity_std_units_curated_comments.csv"), index=False
)

# 12. Save final preprocessed activities
print("Saving results...")
activities_all_raw.to_csv(
    os.path.join(DATAPATH, "chembl_processed", "05_activities_preprocessed.csv"), index=False
)
print("Done.")
