from rdkit.Chem import Descriptors
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import os


root = os.path.dirname(os.path.abspath(__file__))

print("Loading data...")

# Load ChEMBL activities
activities = pd.read_csv(os.path.join(root, "..", "config", "chembl_activities", "activities.csv"), low_memory=False)

# Filter columns
columns = ['activity_id', 'assay_id', 'molregno','standard_relation', 'standard_value', 'standard_units', 'standard_type', 'activity_comment',
           'data_validity_comment', 'pchembl_value','standard_upper_value','standard_text_value', 'action_type']
activities = activities[columns]

# Load assays
assays = pd.read_csv(os.path.join(root, "..", "config", "chembl_activities", "assays.csv"), low_memory=False)

# Load targets
targets = pd.read_csv(os.path.join(root, "..", "config", "chembl_activities", "target_dictionary.csv"), low_memory=True)

# Load compounds
compounds = pd.read_csv(os.path.join(root, "..", "config", "chembl_processed", "compound_info.csv"), low_memory=True)

print("Merging data...")

# Merge with assays
NEW_ACTIVITIES = activities.merge(assays[['assay_id', 'chembl_id', 'assay_type', 'confidence_score', 'tid', 'assay_organism']].
                                  rename(columns={'chembl_id': 'assay_chembl_id', 'confidence_score': 'assay_confidence_score'}), on='assay_id', how='left')

# Merge with targets
NEW_ACTIVITIES = NEW_ACTIVITIES.merge(targets[['tid', 'target_type', "organism", "chembl_id", 'tax_id']].
                                      rename(columns={'chembl_id': 'target_chembl_id', 'organism': 'target_organism', 'tax_id': 'target_tax_id'}), on='tid', how='left')

# Merge with compounds
NEW_ACTIVITIES = NEW_ACTIVITIES.merge(compounds[["molregno", 'chembl_id', 'MW', "canonical_smiles"]].rename(columns={'chembl_id': 'compound_chembl_id'}), on='molregno', how='left')

# Specify final columns
FINAL_COLUMNS = [
    'activity_id',
    'assay_id',
    'assay_chembl_id',
    'assay_type',
    'assay_confidence_score',
    'assay_organism',
    'tid',
    'target_type',
    'target_organism',
    'target_chembl_id',
    'target_tax_id',
    "compound_chembl_id",
    "canonical_smiles",
    "MW",
    'standard_relation',
    'standard_value', 
    'standard_units', 
    'standard_type', 
    'activity_comment',
    'pchembl_value', 
    'standard_text_value'
    ]

# Select final columns
NEW_ACTIVITIES = NEW_ACTIVITIES[FINAL_COLUMNS]

print("Saving results...")

# Create file
NEW_ACTIVITIES.to_csv(os.path.join(root, "..", "config", "chembl_processed", "activities_all_raw.csv"), index=False)