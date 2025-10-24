# Use a conda env with ollama installed
# Run this code on a GPU machine
from collections import Counter
from zipfile import ZipFile, ZIP_DEFLATED
from tqdm import tqdm
import pandas as pd
import numpy as np
import random
import ollama
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH
random.seed(42)

# List of pathogens to process
pathogens = ["Mycobacterium tuberculosis"]

# Load assays and docs information
assays = pd.read_csv(os.path.join(CONFIGPATH, "chembl_activities", "assays.csv"), low_memory=False)
docs = pd.read_csv(os.path.join(CONFIGPATH, "chembl_activities", "docs.csv"), low_memory=False)
assay_type_map = {"F": "Functional", "B": "Binding", "T": "Toxicity", "A": "ADME", "P": "Physicochemical", "U": "Uncategorized"}

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

SYSTEM = f"""
You are a ChEMBL biodata curator. Your task is to write a complete, accurate and standardized description of a given biological assay.

Formatting instructions:
- The description must be structured into three paragraphs (enumerated below), each of 80-120 words.
- Each paragraph must begin with a bold markdown title in the exact format:
  **1. Assay description** \newline
  **2. Outcome interpretation** \newline
  **3. Results and insights** \newline
- The assay description paragraph must explain the objective of the assay, the experimental system, and methodology. Specify the biological target and target type, \
    the pathogen under study (and the corresponding strain, if available), the assay format (e.g., cell-based, binding), the detection method, and any relevant experimental \
    conditions (e.g., temperature, compound concentration).
- The outcome interpretation paragraph must describe how assay outputs are measured and interpreted. Specify how results relate to biological activity or target modulation, \
    controls, reference compounds, signal thresholds and normalization steps. Identify the direction of the biological activity: (-1) if lower values lead to higher activity \
    e.g., IC50 or percentage of survival; (+1) if higher values result in higher activity e.g., percentage of growth inhibition or percentage of effect; (0) if it’s \
    inconclusive or not trivial, e.g., clearance or solubility). Notice that artefacts may eventually appear e.g., data with percentages of growth inhibition < 0% or > 100%. \
    Consider those artefacts as noise.  
- The results and insights paragraph must summarize typical activity ranges, notable behaviors (e.g., agonists, inhibitors), data quality and curation notes that support \
    integration and reproducibility. Highlight meaningful observations from the distribution of activity data. It must be coherent with the outcome interpretation paragraph \
    (i.e. the direction of the biological activity). Additionally, elaborate on the chemical diversity of the assay compounds, basing your description on the number of observed \
    clusters at different ECFP4 Tanimoto similarity cut-offs (e.g., if 10 compounds lead to 10 clusters at a 0.3 ECFP4 Tanimoto similarity cut-off, the assay is \
    chemically diverse. On the contrary, if 10 compounds lead to a single cluster at a ECFP4 0.85 Tanimoto similarity cut-off, the set is probably a chemical series). \
    Finally, you will be provided 20 sampled smiles maximum (top10 and 10 randomly selected) from the assay with their corresponding activity values: try to identify \
    common scaffolds or functional groups related with antimicrobial activity (e.g., quinolones or sulfonamides), particularly among active compounds (notice \
    that actives will have lower values if direction is (-1) or higher if direction is (+1)).  
- Separate paragraphs with a single blank line.
- Use only standard ASCII spacing for all numbers, units, and symbols.
- Insert commas in numbers when necessary (e.g., 1,000; 100,000).
- Do not insert non-breaking spaces, narrow spaces, or special typographic characters.
- Do not use tables or hidden formatting.
- Use scientific and formal language.
- Avoid speculation, informal expressions or fabricated information.
- If any relevant data is missing (reported as ‘nan’), state “not reported” rather than inventing details.
- Do NOT write Q&A, lists, bullets, or add any heading besides the three above.
- Think as much as needed.

"""

# For each pathogen
for pathogen in pathogens:
    
    # Get assays info
    pathogen_code = get_pathogen_code(pathogen)
    ASSAYS_INFO = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, 'assays.csv'))

    # Get clusters info
    CLUSTERS_INFO = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, 'assays_clusters.csv'))

    # Create output directory
    PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code)
    os.makedirs(os.path.join(PATH_TO_OUTPUT, "descriptions"), exist_ok=True)

    # Loading ChEMBL data for that pathogen
    print(f"Loading ChEMBL preprocessed data for {pathogen_code}...")
    ChEMBL = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, f"{pathogen_code}_ChEMBL_data.csv"), low_memory=False)
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL['compound_chembl_id']))}")

    for assay_data in ASSAYS_INFO[['assay_type', 'assay_organism', 'target_type', 'target_organism', 'activity_type', 'unit', 'activities', 'cpds', 'assay_id']].values[:]:

        # Getting assay ID and doc information
        assay_id = assay_data[8]
        doc_id = assays[assays['chembl_id'] == assay_id]['doc_id'].tolist()[0]

        # Getting ChEMBL bioactivities, compounds and clusters
        if type(assay_data[5]) == str:
            assay_activities = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == assay_data[4]) & (ChEMBL['unit'] == assay_data[5])]["value"].astype(float).tolist()
            compounds = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == assay_data[4]) & (ChEMBL['unit'] == assay_data[5])]['canonical_smiles'].tolist()
            relations = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == assay_data[4]) & (ChEMBL['unit'] == assay_data[5])]['relation'].tolist()
            clusters = CLUSTERS_INFO[(CLUSTERS_INFO['assay_id'] == assay_id) & (CLUSTERS_INFO['activity_type'] == assay_data[4]) & (CLUSTERS_INFO['unit'] == assay_data[5])]
        else:
            assay_activities = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == assay_data[4]) & (ChEMBL['unit'].isna())]["value"].astype(float).tolist()
            compounds = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == assay_data[4]) & (ChEMBL['unit'].isna())]['canonical_smiles'].tolist()
            relations = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == assay_data[4]) & (ChEMBL['unit'].isna())]['relation'].tolist()
            clusters = CLUSTERS_INFO[(CLUSTERS_INFO['assay_id'] == assay_id) & (CLUSTERS_INFO['activity_type'] == assay_data[4]) & (CLUSTERS_INFO['unit'].isna())]

        # Get cluster data
        clusters = clusters[["clusters_0.3", "clusters_0.6", "clusters_0.85"]]
        assert len(clusters) == 1
        clusters = clusters.values[0]

        # Get compounds
        compounds_notnans = [f"{i} --> {assay_data[4]} {j} {k} {assay_data[5]}" for i,j,k in zip(compounds, relations, assay_activities) if np.isnan(k) == False]
        compounds_notnans = sorted(compounds_notnans, key=lambda x: float(x.split()[-2]))[::-1]
        compounds_nans = [f"{i} --> {assay_data[4]} {j} {k} {assay_data[5]}" for i,j,k in zip(compounds, relations, assay_activities) if np.isnan(k) == True]
        COMPOUNDS = []
        if len(compounds_notnans) >= 20:
            COMPOUNDS.extend(compounds_notnans[:10] + random.sample(compounds_notnans[10:], 10))
        else:
            COMPOUNDS.extend(compounds_notnans + compounds_nans)  

        # Getting activities that are nans
        assay_activities_nans = [i for i in assay_activities if np.isnan(i)]  # Caution if the number of non-nans activities is 0, the percentile calculation will fail
        assay_activities = [i for i in assay_activities if np.isnan(i) == False]

        result = {
            "Assay ChEMBL ID": assay_id,
            "Assay type": assay_type_map[assay_data[0]],
            "Assay organism": assay_data[1],
            "Assay description": assays[assays['chembl_id'] == assay_id]['description'].tolist()[0],
            "Assay strain": assays[assays['chembl_id'] == assay_id]['assay_strain'].tolist()[0],
            "Assay category": assays[assays['chembl_id'] == assay_id]['assay_category'].tolist()[0],
            "Assay test type": assays[assays['chembl_id'] == assay_id]['assay_test_type'].tolist()[0],
            "Assay cell type": assays[assays['chembl_id'] == assay_id]['assay_cell_type'].tolist()[0],
            "Document title": docs[docs['doc_id'] == doc_id]['title'].tolist()[0],
            "Document abstract": docs[docs['doc_id'] == doc_id]['abstract'].tolist()[0],
            "Document journal": docs[docs['doc_id'] == doc_id]['journal'].tolist()[0],
            "Document PubMed ID": docs[docs['doc_id'] == doc_id]['pubmed_id'].tolist()[0],
            "Document DOI": docs[docs['doc_id'] == doc_id]['doi'].tolist()[0],
            "Target type": assay_data[2],
            "Target organism": assay_data[3],
            "Activity type": assay_data[4],
            "Unit": assay_data[5],
            "Number of activities": len(assay_activities),
            "Number of activities with nan value": len(assay_activities_nans),
            "Number of compounds": assay_data[7],
            "Activity stats": {
                "Percentile 1": round(np.percentile(assay_activities, 1), 3),
                "Percentile 25": round(np.percentile(assay_activities, 25), 3),
                "Mean": round(np.mean(assay_activities), 3),
                "Median": round(np.percentile(assay_activities, 50), 3),
                "Percentile 75": round(np.percentile(assay_activities, 75), 3),
                "Percentile 99": round(np.percentile(assay_activities, 99), 3)},
            "Relation stats": dict(Counter(relations)),
            "Number of compound clusters at a ECFP4 Tanimoto similarity cut-off of 0.3": clusters[0],
            "Number of compound clusters at a ECFP4 Tanimoto similarity cut-off of 0.6": clusters[1],
            "Number of compound clusters at a ECFP4 Tanimoto similarity cut-off of 0.85:": clusters[2],
            "Example smiles": "\n" + "\n".join(COMPOUNDS)
        }

        result = "\n".join([i + ": " + str(result[i]) for i in result])
        USER = f"""Below you will find enumerated annotations from the assay under study.\n\n{result}\n\nUsing the information provided, return a standardized description for the assay."""

        if type(assay_data[5]) == str:
            assay_data[5] = assay_data[5].replace('/', 'FwdS')

        # Print data
        with open(os.path.join(PATH_TO_OUTPUT, "descriptions", f"{assay_id}_{assay_data[4]}_{assay_data[5]}_input.txt"), "w") as f:
            f.write(USER)
        
        # Non streaming call
        response = ollama.generate(model='gpt-oss:20b', prompt=SYSTEM + USER, stream=False, think=True)

        # Print response
        with open(os.path.join(PATH_TO_OUTPUT, "descriptions", f"{assay_id}_{assay_data[4]}_{assay_data[5]}_output.txt"), "w") as f:
            f.write(response.response)

        # Create a zip that bundles both generated files for this assay
        base = f"{assay_id}_{assay_data[4]}_{assay_data[5]}"
        in_path = os.path.join(PATH_TO_OUTPUT, "descriptions", f"{base}_input.txt")
        out_path = os.path.join(PATH_TO_OUTPUT, "descriptions", f"{base}_output.txt")
        zip_path = os.path.join(PATH_TO_OUTPUT, "descriptions", f"{base}.zip")

        with ZipFile(zip_path, "w", compression=ZIP_DEFLATED, compresslevel=9) as zf:
            zf.write(in_path, arcname=f"{base}_input.txt")
            zf.write(out_path, arcname=f"{base}_output.txt")

        # Remove text files
        os.remove(in_path)
        os.remove(out_path)

        print(f"✓ Completed {assay_id}")