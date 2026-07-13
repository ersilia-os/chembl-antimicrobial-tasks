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
from default import STAGE1PATH, FP_RADIUS, FP_BITS

def clip_sparse(vect, nBits=FP_BITS):
    MAX_I8 = 127
    arr = np.zeros(nBits, dtype=np.int8)
    for i, v in vect.GetNonzeroElements().items():
        arr[i] = min(v, MAX_I8)
    return arr

def smiles_to_ecfp(smiles, radius=FP_RADIUS, nBits=FP_BITS):
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)
    output_smiles, X = [], []
    for info in tqdm(smiles):
        info = [str(i) for i in info]
        smi = info[4]
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                ecfp = mfpgen.GetCountFingerprint(mol)
                ecfp = clip_sparse(ecfp, nBits=nBits)
                X.append(ecfp)
                output_smiles.append(info)
        except Exception:
            continue
    assert len(output_smiles) == len(X), "Row mismatch between X and output_smiles"   
    return output_smiles, np.array(X, dtype=np.int8)

print("Step 06")

# Read compound_info
print("Loading data")
compound_info = pd.read_csv(os.path.join(STAGE1PATH, "02_compound_info.csv"))
compounds_standardized = pd.read_csv(os.path.join(STAGE1PATH, "03_compound_info_standardized.csv"), low_memory=True)
# 03 writes standardized structures row-aligned to 02 with no join key, so re-attach by position.
# Guard: a future filter/sort/dedup in step 03 would silently misalign every molecule's structure;
# assert the row counts still match so that breakage is loud, not silent. (A keyed merge on
# molregno is the proper fix, deferred to the next stage-1 re-run.)
assert len(compounds_standardized) == len(compound_info), (
    f"03_compound_info_standardized rows ({len(compounds_standardized)}) != 02_compound_info rows "
    f"({len(compound_info)}) - positional re-attach would misalign structures; re-run step 03")
compound_info['standardized_smiles'] = compounds_standardized['standardized_smiles']
compound_info['standardized_MW'] = compounds_standardized['standardized_MW']
compound_info["molregno"] = compound_info["molregno"].astype(str)
compound_info = compound_info[compound_info['standardized_smiles'].isna() == False].reset_index(drop=True)

# Get only useful data
smiles = compound_info[["molregno", "standard_inchi", "standard_inchi_key", "chembl_id", "standardized_smiles"]].values.tolist()

# Calculate Morgan
print("Calculating Morgan Fingerprints...")
smiles, X_morgan = smiles_to_ecfp(smiles)

print(f"Original number of compounds: {len(compound_info)}")
print(f"Final number of compounds: {len(smiles), X_morgan.shape}")

print("Saving results to H5 file...")
with h5py.File(os.path.join(STAGE1PATH, "06_chembl_ecfps.h5"), "w") as f:
    smiles_dt = h5py.string_dtype(encoding="utf-8")
    f.create_dataset("smiles", data=np.asarray(smiles, dtype=object), dtype=smiles_dt, compression="gzip", chunks=True)
    f.create_dataset("X_morgan", data=X_morgan.astype(np.int8), compression="gzip", chunks=True)