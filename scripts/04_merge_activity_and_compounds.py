import pandas as pd
import sys
import os


root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH

print("Step 04")
print("Loading data...")

# Load ChEMBL activities
activities = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "activities.csv"), low_memory=False)

# Filter columns
columns = ['activity_id', 'assay_id', 'molregno','standard_relation', 'standard_value', 'standard_units', 'standard_type', 'activity_comment',
           'data_validity_comment', 'pchembl_value','standard_upper_value','standard_text_value', 'action_type', 'bao_endpoint']
activities = activities[columns]

# Load assays
assays = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "assays.csv"), low_memory=False)

# Load targets
targets = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "target_dictionary.csv"), low_memory=True)

# Load compounds
compounds = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "02_compound_info.csv"), low_memory=True)
compounds_standardized = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "03_compound_info_standardized.csv"), low_memory=True)
compounds['standardized_smiles'] = compounds_standardized['standardized_smiles']
compounds['standardized_MW'] = compounds_standardized['standardized_MW']

print("Merging data...")

# Merge with assays
new_activities = activities.merge(assays[['assay_id', 'chembl_id', 'assay_type', 'confidence_score', 'tid', 'assay_organism', 'assay_tax_id', 'assay_strain', 'doc_id']].
                                  rename(columns={'chembl_id': 'assay_chembl_id', 'confidence_score': 'assay_confidence_score'}), on='assay_id', how='left')

# Merge with targets
new_activities = new_activities.merge(targets[['tid', 'target_type', "organism", "chembl_id", 'tax_id']].
                                      rename(columns={'chembl_id': 'target_chembl_id', 'organism': 'target_organism', 'tax_id': 'target_tax_id'}), on='tid', how='left')

# Merge with compounds
new_activities = new_activities.merge(compounds[["molregno", 'chembl_id', "standardized_MW", "standardized_smiles"]].rename(columns={'chembl_id': 'compound_chembl_id'}), on='molregno', how='left')

# Specify final columns
final_columns = [
    'activity_id',
    'assay_id',
    'assay_chembl_id',
    'assay_type',
    'assay_confidence_score',
    'assay_organism',
    'assay_tax_id',
    'assay_strain',
    'doc_id',
    'tid',
    'target_type',
    'target_organism',
    'target_chembl_id',
    'target_tax_id',
    "compound_chembl_id",
    "standardized_smiles",
    "standardized_MW",
    'standard_relation',
    'standard_value', 
    'standard_units', 
    'standard_type', 
    'activity_comment',
    'pchembl_value', 
    'standard_text_value',
    'bao_endpoint'
    ]

# Select final columns
new_activities = new_activities[final_columns]

new_activities = new_activities.rename(columns={'standardized_smiles': 'smiles', "standardized_MW": "mw"})

print("Saving results...")

# Create file
new_activities.to_csv(os.path.join(DATAPATH, "chembl_processed", "04_activities_all_raw.csv"), index=False)