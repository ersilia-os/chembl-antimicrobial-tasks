# Use conda env with bitbirch installed
from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import bblean
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import load_pathogen

# Load pathogen info
pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

# Thresholds - Tanimoto Coefficient
thrs = [0.3, 0.6, 0.85]
OUTPATH = os.path.join(root, "..", "output", pathogen_code)

print("Step 10")

print(f"Loading cleaned chembl preprocessed data for {pathogen_code}...")
chembl = pd.read_csv(os.path.join(OUTPATH, f"08_chembl_cleaned_data.csv.gz"), low_memory=False)
print(f"Calculating clusters for pathogen: {pathogen_code}...")
assays_info = pd.read_csv(os.path.join(OUTPATH, '08_assays_cleaned.csv'))
assays_info = assays_info[['assay_id', 'activity_type', 'unit', 'activities', 'cpds']].copy()
print(f"Original number of assays: {len(assays_info)}")
clusters = []

# For each assay
for assay_id, activity_type, unit in tqdm(assays_info[['assay_id', 'activity_type', 'unit']].values):

    try:

        # If unit is nan
        if pd.isna(unit):
            compounds_assay = chembl[(chembl['assay_chembl_id'] == assay_id) & (chembl['activity_type'] == activity_type) & 
                                    (chembl['unit'].isna() == True)]['smiles'].tolist()
        else:
            compounds_assay = chembl[(chembl['assay_chembl_id'] == assay_id) & (chembl['activity_type'] == activity_type) & 
                                    (chembl['unit'] == unit)]['smiles'].tolist()

        # Get list of unique compounds
        compounds_assay = list(set(compounds_assay))

        # Calculate ECFP4s
        fps = bblean.fps_from_smiles(compounds_assay, kind="ecfp4", pack=True, n_features=2048)

        # Clustering using BitBirch
        thr_to_clusters = {}
        for c, thr in enumerate(thrs):
            try:
                bb_tree = bblean.BitBirch(branching_factor=128, threshold=thr, merge_criterion="diameter")
                bb_tree.fit(fps)
                clusters = len(bb_tree.get_cluster_mol_ids())
                thr_to_clusters[thr] = clusters
            except:
                thr_to_clusters[thr] = np.nan
        clusters.append([thr_to_clusters[thr] for thr in thrs])

    except:

        clusters.append([np.nan] * len(thrs))

for c, thr in enumerate(thrs):
    assays_info[f'clusters_{thr}'] = np.array(clusters)[:, c]

# Save cluster results
assays_info.to_csv(os.path.join(root, OUTPATH, '10_assays_clusters.csv'), index=False)
