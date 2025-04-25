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

def load_fingerprints(data_dir):
    X = np.load(os.path.join(data_dir, "014_fingerprints.npy"))
    with open(os.path.join(data_dir, "014_fingerprints_inchikeys.txt"), "r") as f:
        keys = f.read().splitlines()
    return X, keys

def load_random_compounds_from_chembl(chembl_data, N):
    np.random.seed(42)
    chembl_data = pd.read_csv(chembl_data, sep='\t', low_memory=False)
    chembl_data = chembl_data[(~chembl_data['Smiles'].isna()) & (chembl_data['Type'] == 'Small molecule')][['Smiles', 'Inchi Key']].values
    indices = np.random.choice(chembl_data.shape[0], size=N, replace=False)
    chembl_data = chembl_data[indices]
    return chembl_data[:,0], chembl_data[:,1]

# Fingerprints are generated in 014_datasets_modelability.py, we just load them here
X, inchikeys = load_fingerprints(data_dir)

# Load random compounds from ChEMBL and generate fingerprints
N = 250000  # set to e.g. 250k
print(f"Sampling {N} compounds from ChEMBL to perform distinguishability studies")
chembl_random_SMILES, chembl_random_IK = load_random_compounds_from_chembl("../data/chembl_35_smallmolecules.tsv", N)
chembl_random_fps = np.array([get_binary_fingerprints_from_smiles(i) for i in chembl_random_SMILES])


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
    
def save_model_and_task(df, X, inchikeys, chembl_random_fps, chembl_random_IK, chembl_random_SMILES, task, dist_tasks_dir):
    """
    df: task dataset - ['inchikey', 'smiles', <activity>]
    X: mtx of Morgan Fps as Bit Vect including all compounds associated to a pathogen i.e. in all tasks
    inchikeys: IKs all compounds associated to a pathogen i.e. in all tasks
    chembl_random_fps: For a random subsample of ChEMBL compounds, mtx of their Morgan Fps as Bit Vect
    chembl_random_IK: For a random subsample of ChEMBL compounds, IKs
    chembl_random_SMILES: For a random subsample of ChEMBL compounds, SMILES
    task: task we're modelling/distinguishing
    dist_tasks_dir: directory to store the twin DIS version of the task MOD
    """
    inchikeys_ = list(df['inchikey'])
    smiles_ = list(df['smiles'])
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
    X = X[pos_idxs]
    y = y[pos_idxs]
    smiles_ = np.array(smiles_)[pos_idxs]
    inchikeys_ = np.array(inchikeys_)[pos_idxs]
    print("Now the loaded dataset has {0} samples".format(X.shape[0]))
    if len(chembl_random_fps) >= 4 * len(pos_idxs):
        np.random.seed(42)
        random_idxs = np.random.choice(chembl_random_fps.shape[0], 4 * len(pos_idxs), replace=False).tolist()
        chembl_random_fps_subset = chembl_random_fps[random_idxs]
        chembl_random_SMILES_subset = chembl_random_SMILES[random_idxs]
        chembl_random_IK_subset = chembl_random_IK[random_idxs]
    else:
        chembl_random_fps_subset = chembl_random_fps
        chembl_random_SMILES_subset = chembl_random_SMILES
        chembl_random_IK_subset = chembl_random_IK
    print(f"{len(chembl_random_fps)} randomly sampled negatives")
    # Merge task positives with ChEMBL negatives
    X = np.concatenate((X, chembl_random_fps_subset), axis=0)
    y = np.concatenate((y, np.array([0] * len(chembl_random_fps_subset))))
    smiles_ = list(smiles_) + list(chembl_random_SMILES_subset)
    inchikeys_ = list(inchikeys_) + list(chembl_random_IK_subset)
    df_DIS = pd.DataFrame({columns[0]: inchikeys_, columns[1]: smiles_, columns[2]: y})
    df_DIS.to_csv(os.path.join(dist_tasks_dir, task + ".csv"), index=False)
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
models_dir = os.path.join(data_dir, "015_models_DIS")
os.makedirs(models_dir, exist_ok=True)

# Directory to save the tasks
dist_tasks_dir = os.path.join(data_dir, "015_raw_tasks_DIS")
os.makedirs(dist_tasks_dir, exist_ok=True)


R = []
R_models = []
for task in sorted(os.listdir(tasks_dir)[:10]):
    task = task.replace(".csv", "")
    print("Distinguishing task", task)
    df = pd.read_csv(os.path.join(tasks_dir, task + ".csv"))
    results = distinguishability(df, X, inchikeys, chembl_random_fps)
    R += [(task, results["auroc_avg"], results["auroc_std"], results["num_samples"], results["num_pos_samples"], results["pos:neg"])]
    print("Saving full model for", task)
    clf, results = save_model_and_task(df, X, inchikeys, chembl_random_fps, chembl_random_IK, chembl_random_SMILES, task, dist_tasks_dir)
    R_models += [(task, results["auroc"], results["num_samples"], results["num_pos_samples"], results["pos:neg"])]
    joblib.dump(clf, os.path.join(models_dir, task + ".joblib"), compress=9)

pd.DataFrame(R, columns=["task", "auroc_avg", "auroc_std", "num_samples", "num_pos_samples", "pos:neg"]).to_csv(os.path.join(data_dir, "015_distinguishability.csv"), index=False)
pd.DataFrame(R_models, columns=["task", "auroc", "num_samples", "num_pos_samples", "pos:neg"]).to_csv(os.path.join(data_dir, "015_models_DIS.csv"), index=False)
