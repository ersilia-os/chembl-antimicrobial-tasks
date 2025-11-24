from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

# List of pathogens to process
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

ASSAYS = []

# For each pathogen
for pathogen in tqdm(pathogens):

    # Get pathogen code
    pathogen_code = get_pathogen_code(pathogen)

    # Get assay results
    assays = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, "assays.csv"))
    ASSAYS.append(assays)

    # Create descriptions dir
    os.makedirs(os.path.join(root, "..", "output", pathogen_code, "descriptions"), exist_ok=True)

# Concat assays
ASSAYS = pd.concat(ASSAYS, ignore_index=True).values
ASSAYS = np.array_split(ASSAYS, 50000)[:5]

print(f"Number of assays: {len([0 for i in ASSAYS for j in i])}")
print(f"Number of assay tasks: {len(ASSAYS)}")

# Save pickle file
pickle.dump(ASSAYS, open(os.path.join(root, "..", "tmp", "assays.pkl"), "wb"))