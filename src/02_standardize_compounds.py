from chembl_structure_pipeline import standardizer
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Descriptors
from tqdm import tqdm
import datamol as dm
import os
import pandas as pd
import sys
RDLogger.DisableLog('rdApp.*')  # Disable RDKit warnings

def get_canonical_smiles_datamol(smiles):
    try:
        mol = dm.to_mol(smiles)
        mol = dm.fix_mol(mol)
        mol = dm.sanitize_mol(mol)
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return mol, smiles
    except:
        return None, ""

def get_canonical_smiles_rdkit(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return mol, smiles
    except:
        return None, ""

def get_canonical_smiles(smiles):
    smiles = str(smiles).strip()
    try:
        mol, canonical_smiles = get_canonical_smiles_datamol(smiles)
        if mol is not None:
            return mol, canonical_smiles
    except Exception:
        pass
    try:
        mol, canonical_smiles = get_canonical_smiles_rdkit(smiles)
        return mol, canonical_smiles
    except Exception:
        return None, ""
    
def get_standardized_smiles(mol):
    try:
        mol, _ = standardizer.get_parent_mol(mol)
        mol = standardizer.standardize_mol(mol)
        standardized_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return mol, standardized_smiles
    except:
        return None, ""

def calculate_mw(mol):
    try:
        mw = Descriptors.MolWt(mol)
        return str(round(mw, 3))
    except:
        return None

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH

print("Step 02")
print("Loading compound SMILES")
compounds = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "compound_info.csv"))
SMILES = compounds['canonical_smiles']

OUTPUT = []

print("Standardizing compounds and recalculating Molecular Weight")
for smiles in tqdm(SMILES):

    # Get canonical SMILES
    mol, canonical_smiles = get_canonical_smiles(smiles)

    # Get standardized SMILES
    mol, standardized_smiles = get_standardized_smiles(mol)

    # Calculate mw
    mw = calculate_mw(mol)

    # Store results
    OUTPUT.append([standardized_smiles, mw])

OUTPUT = pd.DataFrame(OUTPUT, columns=["standardized_smiles", 'standardized_MW'])
OUTPUT.to_csv(os.path.join(DATAPATH, "chembl_processed", "compound_info_standardized.csv"), index=False)