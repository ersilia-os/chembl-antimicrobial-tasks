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
import joblib

random_seed = 54

np.random.seed(random_seed)
random.seed(random_seed)

parser = argparse.ArgumentParser("Estimate modelability of datasets")
parser.add_argument("--pathogen_code", type=str)
parser.add_argument("--output_dir", type=str)
parser.add_argument("--organism", action="store_true", help="Flag for organism task")
parser.add_argument("--protein", action="store_true", help="Flag for protein task")
args = parser.parse_args()

pathogen_code = args.pathogen_code
data_dir = args.output_dir

if args.organism == True:
    tasks_dir = os.path.join(data_dir, "013a_raw_tasks_MOD")
elif args.protein == True:
    tasks_dir = os.path.join(data_dir, "013b_raw_tasks_MOD")

def get_binary_fingerprints_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024), dtype=int)
    return fp

ik_smi_pairs = []
for f in os.listdir(tasks_dir):
    print("Reading task", f)
    df = pd.read_csv(os.path.join(tasks_dir, f))
    for i, row in df.iterrows():
        ik_smi_pairs.append((row['inchikey'], row['smiles']))
ik_smi_pairs = list(set(ik_smi_pairs))

X = np.zeros((len(ik_smi_pairs), 1024), dtype=int)
keys = []
for i, (ik, smi) in tqdm(enumerate(ik_smi_pairs)):
    X[i] = get_binary_fingerprints_from_smiles(smi)
    keys.append(ik)
np.save(os.path.join(data_dir, "014_fingerprints.npy"), X)
with open(os.path.join(data_dir, "014_fingerprints_inchikeys.txt"), "w") as f:
    for k in keys:
        f.write(k + "\n")

def load_fingerprints(data_dir):
    X = np.load(os.path.join(data_dir, "014_fingerprints.npy"))
    with open(os.path.join(data_dir, "014_fingerprints_inchikeys.txt"), "r") as f:
        keys = f.read().splitlines()
    return X, keys

X, inchikeys = load_fingerprints(data_dir)

def modelability(df, X, inchikeys):
    inchikeys_ = list(df['inchikey'])
    columns = list(df.columns)
    assert len(columns) == 3, "The dataframe must have 3 columns"
    y = np.array(df[columns[-1]], dtype=int)
    indices = {}
    for i, ik in enumerate(inchikeys):
        indices[ik] = i
    idxs = [indices[ik] for ik in inchikeys_]
    X = X[idxs]
    print("Ready to model dataset with {0} samples".format(X.shape[0]))
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

def save_model(df, X, inchikeys):
    inchikeys_ = list(df['inchikey'])
    columns = list(df.columns)
    assert len(columns) == 3, "The dataframe must have 3 columns"
    y = np.array(df[columns[-1]], dtype=int)
    indices = {}
    for i, ik in enumerate(inchikeys):
        indices[ik] = i
    idxs = [indices[ik] for ik in inchikeys_]
    X = X[idxs]
    print("Ready to save full model for dataset with {0} samples".format(X.shape[0]))
    # skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    # aurocs = []
    # for train, test in tqdm(skf.split(X, y)):
    clf = RandomForestClassifier(n_estimators=100, n_jobs=8, random_state=42)
    print("Fitting model")
    clf.fit(X, y)
    try:
        auroc = roc_auc_score(y, clf.predict_proba(X)[:, 1])
    except:
        print("Caution. AUROC calculation failed. Probably no positives or negatives are found.")
        auroc = np.nan
    print("AUROC", auroc)
    results = {"auroc": round(auroc, 4),
               "num_samples": X.shape[0],
               "num_pos_samples": np.sum(y),
               "pos:neg": round(np.sum(y) / (X.shape[0] - np.sum(y)), 4)}
    return clf, results

# Directory to save the models
models_dir = os.path.join(data_dir, "014_models_MOD")
os.makedirs(models_dir, exist_ok=True)

R = []
R_models = []
for l in sorted(os.listdir(tasks_dir)):
    print("Modeling task", l)
    df = pd.read_csv(os.path.join(tasks_dir, l))
    results = modelability(df, X, inchikeys)
    fname = l[:-4]
    R += [(fname, results["auroc_avg"], results["auroc_std"], results["num_samples"], results["num_pos_samples"], results["pos:neg"])]
    print("Saving full model for", l)
    clf, results = save_model(df, X, inchikeys)
    R_models += [(fname, results["auroc"], results["num_samples"], results["num_pos_samples"], results["pos:neg"])]
    joblib.dump(clf, os.path.join(models_dir, fname + ".joblib"), compress=9)

pd.DataFrame(R, columns=["task", "auroc_avg", "auroc_std", "num_samples", "num_pos_samples", "pos:neg"]).to_csv(os.path.join(data_dir, "014_modelability.csv"), index=False)
pd.DataFrame(R_models, columns=["task", "auroc", "num_samples", "num_pos_samples", "pos:neg"]).to_csv(os.path.join(data_dir, "014_models_MOD.csv"), index=False)

