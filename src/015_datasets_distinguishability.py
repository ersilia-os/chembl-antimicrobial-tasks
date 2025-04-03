from rdkit import Chem
from rdkit.Chem import AllChem
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import random
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import argparse

random_seed = 54

np.random.seed(random_seed)
random.seed(random_seed)

parser = argparse.ArgumentParser("Estimate distinguishability of those datasets labeled as non-modelable")
parser.add_argument("--pathogen_code", type=str)
parser.add_argument("--output_dir", type=str)
args = parser.parse_args()

pathogen_code = args.pathogen_code
data_dir = args.output_dir
tasks_dir = os.path.join(data_dir, pathogen_code, "013_raw_tasks")

def get_binary_fingerprints_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024), dtype=int)
    return fp

def load_fingerprints(data_dir):
    X = np.load(os.path.join(data_dir, pathogen_code, "014_fingerprints.npy"))
    with open(os.path.join(data_dir, pathogen_code, "014_fingerprints_inchikeys.txt"), "r") as f:
        keys = f.read().splitlines()
    return X, keys

def load_random_compounds_from_chembl(chembl_data, N):
    np.random.seed(42)
    chembl_data = pd.read_csv(chembl_data, sep='\t', low_memory=False)
    chembl_data = chembl_data[(~chembl_data['Smiles'].isna()) & (chembl_data['Type'] == 'Small molecule')]['Smiles'].tolist()
    chembl_data = np.random.choice(chembl_data, N, replace=False).tolist()
    return chembl_data

# Fingerprints are generated in 014_datasets_modelability.py, we just load them here
X, inchikeys = load_fingerprints(data_dir)

# Load random compounds from ChEMBL and generate fingerprints
chembl_random = load_random_compounds_from_chembl("../data/chembl_35_smallmolecules.tsv", 10)