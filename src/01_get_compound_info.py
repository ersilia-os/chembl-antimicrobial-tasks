from rdkit.Chem import Descriptors
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import os


root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH

# Define output path
print("Step 01")
output_dir = os.path.join(DATAPATH, "chembl_processed")
os.makedirs(output_dir, exist_ok=True)

# Load tables
print("Loading compound tables")
df1 = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "compound_structures.csv"), low_memory=False)
df2 = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "molecule_dictionary.csv"), low_memory=False)

# Merge tables
print("Merging compound tables")
df_merged = df1.merge(df2[['molregno', 'chembl_id']], on='molregno', how='left')

# Calculate Molecular Weight
MW = []
print("Calculating Molecular Weight for all compounds")
for smiles in tqdm(df_merged['canonical_smiles'].tolist()):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        MW.append(round(mw, 3))
    except:
        MW.append(np.nan)

# Add new column
df_merged['MW'] = MW

# Save file
df_merged = df_merged.sort_values('molregno').reset_index(drop=True)
df_merged.to_csv(os.path.join(output_dir, "compound_info.csv"), index=False)


