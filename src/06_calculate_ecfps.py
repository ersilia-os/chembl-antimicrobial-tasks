from rdkit.Chem import rdFingerprintGenerator
from tqdm import tqdm
from rdkit import Chem
import pandas as pd
import numpy as np
import h5py
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH

def clip_sparse(vect, nBits=2048):
    MAX_I8 = 127
    arr = np.zeros(nBits, dtype=np.int8)
    for i, v in vect.GetNonzeroElements().items():
        arr[i] = min(v, MAX_I8)
    return arr

def smiles_to_ecfp(SMILES, radius=3, nBits=2048):
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)
    OUTPUT_SMILES, X = [], []
    for info in tqdm(SMILES):
        info = [str(i) for i in info]
        smi = info[4]
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                ecfp = mfpgen.GetCountFingerprint(mol)
                ecfp = clip_sparse(ecfp, nBits=nBits)
                X.append(ecfp)
                OUTPUT_SMILES.append(info)
        except Exception:
            continue
    assert len(OUTPUT_SMILES) == len(X), "Row mismatch between X and OUTPUT_SMILES"   
    return OUTPUT_SMILES, np.array(X, dtype=np.int8)

print("Step 06")

# Read compound_info
print("Loading data")
compound_info = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "compound_info.csv"))
compounds_standardized = pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "compound_info_standardized.csv"), low_memory=True)
compound_info['standardized_smiles'] = compounds_standardized['standardized_smiles']
compound_info['standardized_MW'] = compounds_standardized['standardized_MW']
compound_info["molregno"] = compound_info["molregno"].astype(str)
compound_info = compound_info[compound_info['standardized_smiles'].isna() == False].reset_index(drop=True)

# Get only useful data
SMILES = compound_info[["molregno", "standard_inchi", "standard_inchi_key", "chembl_id", "standardized_smiles"]].values.tolist()

# Calculate Morgan
print("Calculating Morgan Fingerprints...")
SMILES, X_Morgan = smiles_to_ecfp(SMILES)

print(f"Original number of compounds: {len(compound_info)}")
print(f"Final number of compounds: {len(SMILES), X_Morgan.shape}")

print("Saving results to H5 file...")
with h5py.File(os.path.join(DATAPATH, "chembl_processed", "ChEMBL_ECFPs.h5"), "w") as f:
    smiles_dt = h5py.string_dtype(encoding="utf-8")
    f.create_dataset("SMILES", data=np.asarray(SMILES, dtype=object), dtype=smiles_dt, compression="gzip", chunks=True)
    f.create_dataset("X_morgan", data=X_Morgan.astype(np.int8), compression="gzip", chunks=True)