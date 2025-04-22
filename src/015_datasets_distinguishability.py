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
N = 200000
print(f"Sampling {N} compounds from ChEMBL to perform distinguishability studies")
chembl_random = sorted(load_random_compounds_from_chembl("../data/chembl_35_smallmolecules.tsv", N))
chembl_random_fps = np.array([get_binary_fingerprints_from_smiles(i) for i in chembl_random])

# Get list of tasks to distinguish
modelability = pd.read_csv(os.path.join(data_dir, pathogen_code, "014_modelability.csv"))
tasks_dist = sorted(modelability['task'])
# tasks_dist = sorted(modelability[modelability['auroc_avg'] > 0.7]['task'])
print(len(tasks_dist))


def distinguishability(df, X, inchikeys, chembl_random_fps):
    inchikeys_ = list(df['inchikey'])
    columns = list(df.columns)
    assert len(columns) == 3, "The dataframe must have 3 columns"
    y = np.array(df[columns[-1]], dtype=int)
    indices = {}
    for i, ik in enumerate(inchikeys):
        indices[ik] = i
    idxs = [indices[ik] for ik in inchikeys_]
    X = X[idxs]
    pos_idxs = np.where(y == 1)[0]
    print("Loaded dataset with {0} samples".format(X.shape[0]))
    print("{0} are positive".format(len(pos_idxs)))
    if len(chembl_random_fps) >= 4 * len(pos_idxs):
        np.random.seed(42)
        random_idxs = np.random.choice(chembl_random_fps.shape[0], 4 * len(pos_idxs), replace=False).tolist()
        chembl_random_fps = chembl_random_fps[random_idxs]
    print(f"{len(chembl_random_fps)} randomly sampled negatives")
    X = np.concatenate((X[pos_idxs], chembl_random_fps), axis=0)
    y = np.concatenate((y[pos_idxs], np.array([0] * len(chembl_random_fps))))
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aurocs = []
    for train, test in tqdm(skf.split(X, y)):
        clf = RandomForestClassifier(n_estimators=100, n_jobs=8, random_state=42)
        print("Fitting model")
        clf.fit(X[train], y[train])
        try:
            aurocs += [roc_auc_score(y[test], clf.predict_proba(X[test])[:, 1])]
        except:
            print("Caution. AUROC calculation failed. Probably no positives or negatives are found.")
            aurocs += [np.nan]
        print("AUROC", aurocs[-1])
    results = {"auroc_avg": round(np.mean(aurocs), 4),
               "auroc_std": round(np.std(aurocs), 4),
               "num_samples": X.shape[0],
               "num_pos_samples": np.sum(y),
               "pos:neg": round(np.sum(y) / (X.shape[0] - np.sum(y)), 4)}
    return results
    


R = []
for task in tasks_dist:
    print("Distinguishing task", task)
    df = pd.read_csv(os.path.join(tasks_dir, task + ".csv"))
    results = distinguishability(df, X, inchikeys, chembl_random_fps)
    R += [(task, results["auroc_avg"], results["auroc_std"], results["num_samples"], results["num_pos_samples"], results["pos:neg"])]

pd.DataFrame(R, columns=["task", "auroc_avg", "auroc_std", "num_samples", "num_pos_samples", "pos:neg"]).to_csv(os.path.join(data_dir, pathogen_code, "015_distinguishability.csv"), index=False)
