from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
import re

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))

# Load data
activities_all_raw = pd.read_csv(os.path.join(root, "..", "config", "chembl_processed", "activities_all_raw.csv"), low_memory=False)

# 1. Flagging activity comments
activity_comments_bin = pd.read_csv(os.path.join(root, "..", "config", "manual_curation", "activity_comments_manual_curation.csv"), low_memory=False)
activity_comments_act = set(activity_comments_bin[activity_comments_bin['manual_curation'] == 1]['activity_comment'])
activity_comments_inact = set(activity_comments_bin[activity_comments_bin['manual_curation'] == -1]['activity_comment'])

# 2. Flagging standard text
standard_text_bin = pd.read_csv(os.path.join(root, "..", "config", "manual_curation", "standard_text_manual_curation.csv"), low_memory=False)
standard_text_act = set(standard_text_bin[standard_text_bin['manual_curation'] == 1]['standard_text_value'])
standard_text_inact = set(standard_text_bin[standard_text_bin['manual_curation'] == -1]['standard_text_value'])

# 3. Unit conversion
unit_conversion = pd.read_csv(os.path.join(root, "..", "config", "chembl_processed", "unit_conversion.csv"))
standard_unit_to_final_unit = {i: j for i,j in zip(unit_conversion['standard_units'], unit_conversion['final_unit'])}
standard_unit_to_conversion_formula = {i: j for i,j in zip(unit_conversion['standard_units'], unit_conversion['conversion_formula'])}

# 5. Relations py dict
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

# 6. pChEMBL calculation
def calculate_pchembl(uM):
    try:
        value = uM * 1e-6
        pchembl_value = np.clip(-np.log10(value), 1, 9)
        return pchembl_value
    except:
        return np.nan



# 1. Cleaning activity comments
print("Cleaning activity comments...")
NEW_ACTIVITY_COMMENT = []
for act_comment in tqdm(activities_all_raw['activity_comment']):
    if str(act_comment) == 'nan':
        NEW_ACTIVITY_COMMENT.append(0)
    elif act_comment in activity_comments_act:
        NEW_ACTIVITY_COMMENT.append(1)
    elif act_comment in activity_comments_inact:
        NEW_ACTIVITY_COMMENT.append(-1)
    else:
        NEW_ACTIVITY_COMMENT.append(0)

activities_all_raw['new_activity_comment'] = NEW_ACTIVITY_COMMENT
print(f"New activity comments: {dict(Counter(activities_all_raw['new_activity_comment']))}")


# 2. Cleaning standard text
print("Cleaning standard text...")
NEW_STANDARD_TEXT = []
for std_text_value in tqdm(activities_all_raw['standard_text_value']):
    if str(std_text_value) == 'nan':
        NEW_STANDARD_TEXT.append(0)
    elif std_text_value in standard_text_act:
        NEW_STANDARD_TEXT.append(1)
    elif std_text_value in standard_text_inact:
        NEW_STANDARD_TEXT.append(-1)
    else:
        NEW_STANDARD_TEXT.append(0)

activities_all_raw['new_standard_text_value'] = NEW_STANDARD_TEXT
print(f"New standard text: {dict(Counter(activities_all_raw['new_standard_text_value']))}")


# 3. Harmonizing units and values
NEW_VALUES, NEW_UNITS = [], []
print("Harmonizing units and converting values")
for mw, std_value, std_unit in tqdm(activities_all_raw[['MW', 'standard_value', 'standard_units']].values):

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
df.to_csv(os.path.join(root, "..", "config", "chembl_processed", "converted_units.csv"), index=False)

# Dict mapping old units with new units
new_unit_to_old_units = {i: set() for i in set(NEW_UNITS)}
for i,j in zip(activities_all_raw['converted_units'], activities_all_raw['standard_units']):
    new_unit_to_old_units[i].add(j)

df = [[unit, len(new_unit_to_old_units[unit]), " ; ".join([str(k) for k in new_unit_to_old_units[unit]])] 
      for unit in sorted(new_unit_to_old_units, key=lambda x: len(new_unit_to_old_units[x]))[::-1]]
df = pd.DataFrame(df, columns=['unit', 'count', 'old_units'])
df.to_csv(os.path.join(root, "..", "config", "chembl_processed", "converted_units_map.csv"), index=False)


# 4. Harmonizing activity types
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
df.to_csv(os.path.join(root, "..", "config", "chembl_processed", "harmonized_types_map.csv"), index=False)

# 5. Clean relations
print("Cleaning relations...")
activities_all_raw["new_relation"] = [convert_relation(i, RELATIONS) for i in tqdm(activities_all_raw["standard_relation"])]
print(f"Old relations: {dict(Counter(activities_all_raw['standard_relation']))}")
print(f"New relations: {dict(Counter(activities_all_raw['new_relation']))}")


# 6. Calculating pChEMBL
print("Calculating pChEMBL values...")
calculated_pChEMBLs = []
for unit, value, pch in tqdm(activities_all_raw[['converted_units', 'converted_values', 'pchembl_value']].values):
    if str(value) != 'nan' and unit == "umol.L-1":
        value = calculate_pchembl(value)
        calculated_pChEMBLs.append(value)
    else:
        calculated_pChEMBLs.append(np.nan)

activities_all_raw['calculated_pChEMBLs'] = calculated_pChEMBLs

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
        "calculated_pChEMBLs": "pchembl_calculated"
        })

activities_all_raw.to_csv(os.path.join(root, "..", "config", 'chembl_processed', 'activities_preprocessed.csv'), index=False)