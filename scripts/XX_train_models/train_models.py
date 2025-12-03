from sklearn.model_selection import StratifiedKFold
from lazyqsar.qsar import LazyBinaryQSAR
from lazyqsar.agnostic import LazyBinaryClassifier
import pandas as pd
import numpy as np
import lazyqsar
import h5py
import tqdm
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))

# Load data
data = pd.read_csv(os.path.join(root, "..", "..", "output", "mtuberculosis", "datasets", "CHEMBL4649948_IC50_umol.L-1_perc_1.csv.gz"))
X = data['canonical_smiles'].tolist()
Y = data['bin'].tolist()
X = np.array(X, dtype='str')
Y = np.array(Y, dtype=np.int8)

# Define stratified 5 fold CV
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Load descriptors
with h5py.File(os.path.join(root, "..", "..", "output", "mtuberculosis", "descriptors.h5"), "r") as f:
    SMILES = f['SMILES'][:]
    X_Morgan = f['X_Morgan'][:]

# Define dict mapping smiles to morgan fingerprints
SMILES_TO_MORGAN = {i[1].decode('utf-8'): j for i,j in zip(SMILES, X_Morgan)}

# For each split
for train, test in skf.split(X, Y):
    print(len(train), len(test))

    # Get train/test indices
    X_train, Y_train = X[train], Y[train]
    X_test, Y_test = X[test], Y[test]

    # Get Morgan Fingerprints and train model
    X_train = np.array([SMILES_TO_MORGAN[i] for i in X_train])
    X_test = np.array([SMILES_TO_MORGAN[i] for i in X_test])
    model = LazyBinaryClassifier(mode="fast")
    model.fit(X=X_train, y=Y_train)

    break