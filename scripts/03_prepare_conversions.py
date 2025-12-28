import pandas as pd
import numpy as np
import sys
import os

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH

# Load data
all_units = pd.read_csv(os.path.join(CONFIGPATH, "chembl_processed", "UnitStringValidations.csv"), skiprows=1)  # file standard_units.csv to https://ucum.nlm.nih.gov/
del all_units['cumulative_prop']
del all_units['Validation Result']
del all_units['Notes']
all_units = all_units[all_units['standard_units'].isna() == False].reset_index(drop=True)

# Load curated data
ucum_GT = pd.read_csv(os.path.join(CONFIGPATH, "manual_curation", "ucum_GT.csv"))

# Valid units
unit_to_valid_unit = {}
for unit, valid_unit in ucum_GT[['units', "val_unit"]].values:
    unit_to_valid_unit[unit] = valid_unit

# Final units
unit_to_final_unit = {}
for unit, final_unit in ucum_GT[['units', "final_unit"]].values:
    unit_to_final_unit[unit] = final_unit

# Conversion formula GT
unit_to_conv_formula = {}
for unit, conv_formula in ucum_GT[['units', "transformer"]].values:
    unit_to_conv_formula[unit] = conv_formula


# mg/dl	
unit_to_final_unit['mg/dl'] = "umol.L-1"
unit_to_conv_formula['mg/dl'] = "standard_value * 10000 /molecular_weight"

# uM
unit_to_final_unit['uM'] = "umol.L-1"

# uG
unit_to_conv_formula['ug'] = "standard_value/1000"

# umol/ml
unit_to_conv_formula['umol/ml'] = "standard_value*1000"

# M l-1
unit_to_conv_formula['M l-1'] = "standard_value*1000000"

# /uL
unit_to_conv_formula['/uL'] = "standard_value*1000000"

# Getting valid units
all_units['valid_unit'] = [unit_to_valid_unit[i] if i in unit_to_valid_unit else np.nan for i in all_units['standard_units']]

# Getting final units
all_units['final_unit'] = [unit_to_final_unit[i] if i in unit_to_final_unit else np.nan for i in all_units['standard_units']]

# Getting conversion formula
all_units['conversion_formula'] = [unit_to_conv_formula[i] if i in unit_to_conv_formula else np.nan for i in all_units['standard_units']]

# Get proportion
total_count = all_units['count'].sum()
all_units['cumulative_prop'] = (all_units['count'].cumsum() / total_count).round(3)

# Fill nans as standard values
final_units, conversion_formula = [], []

for std_unit, valid_unit, final_unit, conv_formula in all_units[["standard_units", "valid_unit", "final_unit", "conversion_formula"]].values:

    # Adjust final units
    if pd.isna(final_unit):
        final_units.append(std_unit)
    else:
        final_units.append(final_unit)

    # Adjust conversion formula
    if pd.isna(conv_formula):
        conversion_formula.append("standard_value")
    else:
        conversion_formula.append(conv_formula)

# Fill nans as standard values
all_units['final_unit'] = final_units
all_units['conversion_formula'] = conversion_formula

# Save results
all_units.to_csv(os.path.join(CONFIGPATH, "chembl_processed", "unit_conversion.csv"), index=False)