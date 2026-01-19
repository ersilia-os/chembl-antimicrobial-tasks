from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import os
import re

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH

# Load data
print("Step 05")
print("Loading data...")
activities_all_raw = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "activities_all_raw.csv"), low_memory=False)

# 2. Flagging activity comments
activity_comments_bin = pd.read_csv(os.path.join(CONFIGPATH, "activity_comments_manual_curation.csv"), low_memory=False)
activity_comments_act = set(activity_comments_bin[activity_comments_bin['manual_curation_activity'] == 1]['activity_comment'])
activity_comments_inact = set(activity_comments_bin[activity_comments_bin['manual_curation_activity'] == -1]['activity_comment'])

# 3. Flagging standard text
standard_text_bin = pd.read_csv(os.path.join(CONFIGPATH, "standard_text_manual_curation.csv"), low_memory=False)
standard_text_act = set(standard_text_bin[standard_text_bin['manual_curation_activity'] == 1]['standard_text_value'])
standard_text_inact = set(standard_text_bin[standard_text_bin['manual_curation_activity'] == -1]['standard_text_value'])

# 4. Unit conversion
unit_conversion = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "unit_conversion.csv"))
standard_unit_to_final_unit = {i: j for i,j in zip(unit_conversion['standard_units'], unit_conversion['final_unit'])}
standard_unit_to_conversion_formula = {i: j for i,j in zip(unit_conversion['standard_units'], unit_conversion['conversion_formula'])}

# 6. Relations py dict
RELATIONS = {"=": "=",
             ">": ">",
             "<": "<",
             ">>": ">",
             ">=": ">",
             "<<": "<",
             "<=": "<",
             np.nan: "=",
             "~": "="}

def convert_relation(i, RELATIONS):
    return RELATIONS[i]

# 7. pChEMBL calculation
def calculate_pchembl(uM):
    try:
        value = uM * 1e-6
        pchembl_value = np.clip(-np.log10(value), 1, 9)
        return pchembl_value
    except:
        return np.nan
    
# 8. Converting doc_id to doc_chembl_id
docs = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "docs.csv"), low_memory=False)
doc_id_to_doc_chembl_id = {i: j for i, j in zip(docs['doc_id'], docs["chembl_id"])}

# 11. Creating text_flag

def create_text_flag(df):

    cond_nan = (df['activity_comment'] == 0) & (df['standard_text'] == 0)
    cond_pos = (df['activity_comment'] == 1) | (df['standard_text'] == 1)
    cond_neg = (df['activity_comment'] == -1) | (df['standard_text'] == -1)

    # Detect row-level conflicts
    conflict = cond_pos & cond_neg
    if conflict.any():
        raise ValueError(
            "Conflicting labels (contains both 1 and -1):\n"
            + df.loc[conflict, ["compound_chembl_id", "activity_comment", "standard_text"]].head(20).to_string())

    # Assign row-level label
    df["text_flag"] = pd.Series(pd.NA, index=df.index, dtype="Int64")
    df.loc[cond_pos, "text_flag"] = 1
    df.loc[cond_neg, "text_flag"] = -1
    df.loc[cond_nan, "text_flag"] = 0

    # Remove original fields
    df = df.drop(columns=['activity_comment', 'standard_text'])
    
    return df

# 1. Removing those activities having no canonical smiles
nans = len(activities_all_raw[activities_all_raw['canonical_smiles'].isna()])
print(f"Removing activities having no associated canonical smiles ({nans}) ...")
activities_all_raw = activities_all_raw[activities_all_raw['canonical_smiles'].isna() == False].reset_index(drop=True)

# 2. Cleaning activity comments
print("Cleaning activity comments...")
NEW_ACTIVITY_COMMENT = []
for act_comment in tqdm(activities_all_raw['activity_comment']):
    act_comment = str(act_comment).strip().lower()
    if act_comment == 'nan':
        NEW_ACTIVITY_COMMENT.append(0)
    elif act_comment in activity_comments_act:
        NEW_ACTIVITY_COMMENT.append(1)
    elif act_comment in activity_comments_inact:
        NEW_ACTIVITY_COMMENT.append(-1)
    else:
        NEW_ACTIVITY_COMMENT.append(0)

activities_all_raw['new_activity_comment'] = NEW_ACTIVITY_COMMENT
print(f"New activity comments: {dict(Counter(activities_all_raw['new_activity_comment']))}")


# 3. Cleaning standard text
print("Cleaning standard text...")
NEW_STANDARD_TEXT = []
for std_text_value in tqdm(activities_all_raw['standard_text_value']):
    std_text_value = str(std_text_value)
    if std_text_value == 'nan':
        NEW_STANDARD_TEXT.append(0)
    elif std_text_value in standard_text_act:
        NEW_STANDARD_TEXT.append(1)
    elif std_text_value in standard_text_inact:
        NEW_STANDARD_TEXT.append(-1)
    else:
        NEW_STANDARD_TEXT.append(0)

activities_all_raw['new_standard_text_value'] = NEW_STANDARD_TEXT
print(f"New standard text: {dict(Counter(activities_all_raw['new_standard_text_value']))}")


# 4. Standardizing units and values
NEW_VALUES, NEW_UNITS = [], []
print("Standardizing units and converting values")
for mw, std_value, std_unit in tqdm(activities_all_raw[['standardized_MW', 'standard_value', 'standard_units']].values):

    # Get conversion formula
    if std_unit in standard_unit_to_conversion_formula:
        conversion_formula = standard_unit_to_conversion_formula[std_unit]
    else:
        # Only when std_unit is nan
        conversion_formula = np.nan

    # Get final_unit
    if std_unit in standard_unit_to_final_unit:
        final_unit = standard_unit_to_final_unit[std_unit]
    else:
        # Only when std_unit is nan
        final_unit = np.nan
    NEW_UNITS.append(final_unit)

    # Get new value
    if str(std_value) != 'nan':
        if str(conversion_formula) != 'nan':
            data = {'standard_value': std_value, 'molecular_weight': mw}
            new_value = eval(conversion_formula, data)
            NEW_VALUES.append(new_value)
        else:
            NEW_VALUES.append(std_value)
    else:
        NEW_VALUES.append(np.nan)

# Save in df
activities_all_raw['converted_values'] = NEW_VALUES
activities_all_raw['converted_units'] = NEW_UNITS

# Saving CSV file with converted units and values
d = dict(Counter(activities_all_raw['converted_units']))
df = [[unit, d[unit]] for unit in sorted(d, key=lambda x: d[x])[::-1]]
df = pd.DataFrame(df, columns=['unit', 'count'])
total_count = np.sum(df['count'])
df['cumulative_prop'] = (df['count'].cumsum() / total_count).round(3)
df.to_csv(os.path.join(DATAPATH, "chembl_processed", "converted_units.csv"), index=False)

# Dict mapping old units with new units
new_unit_to_old_units = {i: set() for i in set(NEW_UNITS)}
for i,j in zip(activities_all_raw['converted_units'], activities_all_raw['standard_units']):
    new_unit_to_old_units[i].add(j)

df = [[unit, len(new_unit_to_old_units[unit]), " ; ".join([str(k) for k in new_unit_to_old_units[unit]])] 
      for unit in sorted(new_unit_to_old_units, key=lambda x: len(new_unit_to_old_units[x]))[::-1]]
df = pd.DataFrame(df, columns=['unit', 'count', 'old_units'])
df.to_csv(os.path.join(DATAPATH, "chembl_processed", "converted_units_map.csv"), index=False)


# 5. Harmonizing activity types
print("Harmonizing activity types")
def harmonize_act_type(act_type):
    # _, spaces, / and \
    return re.sub(r"[_\s./\\]", "", str(act_type).upper().strip())

HARMONIZED_TYPES = [harmonize_act_type(i) for i in tqdm(activities_all_raw['standard_type'])]
activities_all_raw['harmonized_type'] = HARMONIZED_TYPES

harmonized_types_to_types = {i: set() for i in set(HARMONIZED_TYPES)}
for ty, harm_ty in zip(activities_all_raw['standard_type'], activities_all_raw['harmonized_type']):
    harmonized_types_to_types[harm_ty].add(ty)

df = [[ty, len(harmonized_types_to_types[ty]), " ; ".join([str(k) for k in harmonized_types_to_types[ty]])] 
      for ty in sorted(harmonized_types_to_types, key=lambda x: len(harmonized_types_to_types[x]))[::-1]]
df = pd.DataFrame(df, columns=['type', 'count', 'old_types'])
df.to_csv(os.path.join(DATAPATH, "chembl_processed", "harmonized_types_map.csv"), index=False)

# 6. Clean relations
print("Cleaning relations...")
activities_all_raw["new_relation"] = [convert_relation(i, RELATIONS) for i in tqdm(activities_all_raw["standard_relation"])]
print(f"Old relations: {dict(Counter(activities_all_raw['standard_relation']))}")
print(f"New relations: {dict(Counter(activities_all_raw['new_relation']))}")


# 7. Calculating pChEMBL
print("Calculating pChEMBL values...")
calculated_pChEMBLs = []
for unit, value, pch in tqdm(activities_all_raw[['converted_units', 'converted_values', 'pchembl_value']].values):
    if str(value) != 'nan' and unit == "umol.L-1":
        value = calculate_pchembl(value)
        calculated_pChEMBLs.append(value)
    else:
        calculated_pChEMBLs.append(np.nan)
activities_all_raw['calculated_pChEMBLs'] = calculated_pChEMBLs

# 8. Converting doc_id to doc_chembl_id
docs = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "docs.csv"), low_memory=False)
doc_id_to_doc_chembl_id = {i: j for i, j in zip(docs['doc_id'], docs["chembl_id"])}
len_docs = len(docs)
del docs
print(f"Converting Doc IDs...")
print(f"Number of docs: {len_docs}; Number of mappings: {len(doc_id_to_doc_chembl_id)}")
activities_all_raw['doc_id'] = [doc_id_to_doc_chembl_id[i] for i in activities_all_raw['doc_id']]


# 9. Column renaming and cleanup
del activities_all_raw['standard_relation']
del activities_all_raw['standard_value']
del activities_all_raw['standard_units']
del activities_all_raw['standard_type']
del activities_all_raw['activity_comment']
del activities_all_raw['standard_text_value']

activities_all_raw = activities_all_raw.rename(columns={
        "new_relation": "relation",
        "new_activity_comment": "activity_comment",
        "new_standard_text_value": "standard_text",
        "converted_values": "value",
        "converted_units": "unit",
        "harmonized_type": "activity_type",
        "pchembl_value": "pchembl",
        "calculated_pChEMBLs": "pchembl_calculated",
        "doc_id": "doc_chembl_id"
        })

# 10. Converting activity types to their corresponding synonyms
print("Mapping activity types to their synonyms")
synonyms = pd.read_csv(os.path.join(CONFIGPATH, "synonyms.csv"))
for activity, syns in zip(synonyms['activity_type'], synonyms['synonyms']):
    for syn in syns.split(";"):
        activities_all_raw.loc[activities_all_raw['activity_type'] == syn, 'activity_type'] = activity


# 11. Creating a curated activity type - standard units file for manual curation
s = activities_all_raw[["activity_type", "unit"]].astype("string").fillna("")
out = (
s.value_counts(subset=["activity_type", "unit"], dropna=False)
    .reset_index(name="count")
    .sort_values("count", ascending=False, ignore_index=True)
)
total_count = out['count'].sum()
out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)
out['manual_curation_direction'] = np.nan
out.to_csv(os.path.join(DATAPATH, "chembl_processed", "activity_std_units_curated.csv"), index=False)
print(f"Total number of unique activity type - standard unit pairs: {len(out)}")

# More on step 11.
# Creating text flag
activities_all_raw = create_text_flag(activities_all_raw)
cols = ["activity_type", "unit", "text_flag"]
df = activities_all_raw[cols]
flagged = df["text_flag"].isin([1, -1])

# Counting text flag
out = (
    df.assign(flagged=flagged.to_numpy())
      .groupby(["activity_type", "unit"], dropna=False)
      .agg(count=("flagged", "size"), comments=("flagged", "sum"))
      .reset_index()
      .sort_values("count", ascending=False, ignore_index=True)
)

out.to_csv(os.path.join(DATAPATH, "chembl_processed", "activity_std_units_curated_comments.csv"), index=False)

print("Saving results...")
activities_all_raw.to_csv(os.path.join(DATAPATH, 'chembl_processed', 'activities_preprocessed.csv'), index=False)

# print("Removing intermediate activity data")
# os.remove(os.path.join(DATAPATH, "chembl_processed", "activities_all_raw.csv"))