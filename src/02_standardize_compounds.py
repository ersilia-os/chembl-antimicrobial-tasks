from chembl_structure_pipeline import standardizer
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
    mol = dm.to_mol(smiles)
    if not mol:
        return None
    mol = dm.fix_mol(mol)
    mol = dm.sanitize_mol(mol)
    if not mol:
        return None
    smi = dm.to_smiles(mol)
    if not smi:
        return None
    m2 = Chem.MolFromSmiles(smi)
    if not m2:
        return None
    return Chem.MolToSmiles(m2, canonical=True, isomericSmiles=True)

def get_canonical_smiles_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

def get_canonical_smiles(smiles):
    smiles = str(smiles).strip()
    try:
        canonical_smiles = get_canonical_smiles_datamol(smiles)
        if canonical_smiles:
            return canonical_smiles
    except Exception:
        pass
    try:
        return get_canonical_smiles_rdkit(smiles)
    except Exception:
        return None
    
def get_standardized_smiles(smiles):
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return None
    m = standardizer.standardize_mol(m)
    m, _ = standardizer.get_parent_mol(m)
    if not m:
        return None
    m = standardizer.standardize_mol(m)
    if not m:
        return None
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)

def deal_with_nones(string):
    if string is None:
        return ""
    return string

def calculate_mw(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    mw = Descriptors.MolWt(mol)
    return str(round(mw, 3))


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
    canonical_smiles = get_canonical_smiles(smiles)

    # Get standardized SMILES
    standardized_smiles = get_standardized_smiles(canonical_smiles)

    # Calculate mw
    mw = calculate_mw(standardized_smiles)

    # Deal with None's
    canonical_smiles = deal_with_nones(canonical_smiles)
    standardized_smiles = deal_with_nones(standardized_smiles)
    mw = deal_with_nones(mw)

    # Store results
    OUTPUT.append([standardized_smiles, mw])

OUTPUT = pd.DataFrame(OUTPUT, columns=["standardized_smiles", 'standardized_MW'])
OUTPUT.to_csv(os.path.join(DATAPATH, "chembl_processed", "compound_info_standardized.csv"))