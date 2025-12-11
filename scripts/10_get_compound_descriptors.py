from lazyqsar.descriptors.chemeleon import ChemeleonDescriptor
from lazyqsar.descriptors.morgan import MorganFingerprint
from lazyqsar.descriptors.rdkit_descriptors import RDKitDescriptor
import pandas as pd
import numpy as np
import h5py
import tqdm
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))

# List of pathogens to process
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"][8:9]

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

for pathogen in pathogens:

    # Loading pathogen data
    pathogen_code = get_pathogen_code(pathogen)
    print(f"Loading ChEMBL cleaned data for {pathogen_code}...")
    ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")


    # ChEMBL ID to SMILES
    ChEMBL_id_to_SMILES = {i:j for i,j in zip(ChEMBL_pathogen['compound_chembl_id'], ChEMBL_pathogen['canonical_smiles'])}
    ids = sorted(ChEMBL_id_to_SMILES)
    SMILES = [[i, ChEMBL_id_to_SMILES[i]] for i in ids]

    # Calculate Morgan
    print("Calculating Morgan Fingerprints...")
    X_morgan = MorganFingerprint().transform([i[1] for i in SMILES])

    # Calculate rdkit
    print("Calculating RDKit descriptors...")
    X_rdkit = RDKitDescriptor().transform([i[1] for i in SMILES])

    # # Calculate Chemeleon
    # print("Calculating Chemeleon embeddings...")
    # X_chemeleon = ChemeleonDescriptor().transform([i[1] for i in SMILES])

    print("Saving results to H5 file...")
    with h5py.File(os.path.join(OUTPUT, pathogen_code, "descriptors.h5"), "w") as f:
        dt = h5py.string_dtype(encoding='utf-8')
        f.create_dataset("SMILES", data=SMILES, dtype=dt, compression="gzip")
        f.create_dataset("X_morgan", data=X_morgan.astype(np.int8), compression="gzip")
        f.create_dataset("X_rdkit", data=X_rdkit.astype(np.float32), compression="gzip")
        # f.create_dataset("X_chemeleon", data=X_chemeleon)
