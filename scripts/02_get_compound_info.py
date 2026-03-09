"""
Step 02 — Build compound info table with molecular weights.

Merges compound_structures.csv and molecule_dictionary.csv (from step 00) on
molregno to produce a unified compound table that includes chembl_id, InChI,
InChI key, canonical SMILES, and RDKit-calculated molecular weight (MW).

Input:  data/chembl_activities/compound_structures.csv
        data/chembl_activities/molecule_dictionary.csv
Output: data/chembl_processed/02_compound_info.csv
"""

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

print("Step 02")
output_dir = os.path.join(DATAPATH, "chembl_processed")
os.makedirs(output_dir, exist_ok=True)
print("Loading compound tables")
df1 = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "compound_structures.csv"), low_memory=False)
df2 = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "molecule_dictionary.csv"), low_memory=False)
print("Merging compound tables")
df_merged = df1.merge(df2[['molregno', 'chembl_id']], on='molregno', how='left')
MW = []
print("Calculating Molecular Weight for all compounds")
for smiles in tqdm(df_merged['canonical_smiles'].tolist()):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)  # raises if mol is None (invalid SMILES)
        MW.append(round(mw, 3))
    except:
        MW.append(np.nan)  # NaN for missing or unparseable SMILES
df_merged['MW'] = MW

# Save file
df_merged = df_merged.sort_values('molregno').reset_index(drop=True)
df_merged.to_csv(os.path.join(output_dir, "02_compound_info.csv"), index=False)


