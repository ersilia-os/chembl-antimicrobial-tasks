from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from lazyqsar.agnostic import LazyBinaryClassifier
import pandas as pd
import numpy as np
import lazyqsar
import pickle
import h5py
import tqdm
import sys
import os

alpha = int(sys.argv[1])

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))

# Load pickle
pathogen_code, file, column = pickle.load(open(os.path.join(root, "..", "..", "tmp", "models_to_train.pkl"), "rb"))[alpha]

# Load data
sys.stderr.write(str(pathogen_code) + "\n" + str(file) + " -- " + str(column))
sys.stderr.write("\n\n")
sys.stderr.flush()
data = pd.read_csv(os.path.join(root, "..", "..", "output", pathogen_code, "datasets", file))
X = data['canonical_smiles'].astype(str).tolist()
Y = data['bin'].astype(np.int8).to_numpy()

# Define stratified 5 fold CV
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Load descriptors
sys.stderr.write("Loading descriptors...")
sys.stderr.write("\n\n")
sys.stderr.flush()
with h5py.File(os.path.join(root, "..", "..", "output", pathogen_code, "descriptors.h5"), "r") as f:
    SMILES = f['SMILES'][:]
    X_Morgan = f['X_Morgan'][:]
    ### add here rdkit descriptors

# Define dict mapping smiles to morgan fingerprints
SMILES_TO_MORGAN = {
    smiles.decode("utf-8"): fp
    for (chembl_id, smiles), fp in zip(SMILES, X_Morgan)}

# For each split
AUROCS = []
for c, (train, test) in enumerate(skf.split(X, Y)):
    
    sys.stderr.write(f"Training split {c}...")
    sys.stderr.write("\n")
    sys.stderr.flush()

    # Get train/test indices
    X_train, Y_train = [X[i] for i in train], Y[train]
    X_test, Y_test = [X[i] for i in test], Y[test]

    # Get Morgan Fingerprints and train model
    X_train = np.array([SMILES_TO_MORGAN[smi] for smi in X_train])
    X_test  = np.array([SMILES_TO_MORGAN[smi] for smi in X_test])
    model = LazyBinaryClassifier(mode="default")
    model.fit(X=X_train, y=Y_train)

    # Store results
    auroc = round(roc_auc_score(Y_test, model.predict_proba(X_test)[:, 1]), 3)
    AUROCS.append(auroc)
    sys.stderr.write(f"AUROC: {auroc}")
    sys.stderr.write("\n\n")
    sys.stderr.flush()

# Train full model
X = np.array([SMILES_TO_MORGAN[smi] for smi in X])
model = LazyBinaryClassifier(mode="default")
model.fit(X=X, y=Y)

# Save model
model.save(os.path.join(root, "..", "..", "output", pathogen_code, "models", file.replace(".csv.gz", ".zip")))

# Save CV results
with open(os.path.join(root, "..", "..", "output", pathogen_code, "models", file.replace('.gz', '')), "w") as outfile:
    outfile.write(",".join([str(i) for i in AUROCS]))