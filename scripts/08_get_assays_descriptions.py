# Use a conda env with bitbirch installed
from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import ollama
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH

# Loading ChEMBL preprocessed data
print("Loading ChEMBL preprocessed data...")
ChEMBL = pd.read_csv(os.path.join(root, "..", "config", "chembl_processed", "activities_preprocessed.csv"), low_memory=False)
print(f"Original size: {len(ChEMBL)}")
print(f"Filtering out nan values...")
ChEMBL = ChEMBL[ChEMBL['value'].isna() == False].reset_index(drop=True)
print(f"Size after filtering nan values: {len(ChEMBL)}")

# List of pathogens to process
pathogens = ["Mycobacterium tuberculosis"]

# Load assays and docs information
assays = pd.read_csv(os.path.join(CONFIGPATH, "chembl_activities", "assays.csv"), low_memory=False)
docs = pd.read_csv(os.path.join(CONFIGPATH, "chembl_activities", "docs.csv"), low_memory=False)
assay_type_map = {"F": "Functional", "B": "Binding", "T": "Toxicity", "A": "ADME", "P": "Physicochemical", "U": "Uncategorized"}

SYSTEM = """
You are a ChEMBL biodata curator. Your task is to write a complete, accurate and standardized description of a given biological assay.

Formatting instructions:
- The description must be structured into three paragraphs (enumerated below), each of 80-120 words.
- Each paragraph must begin with a bold markdown title in the exact format:
  **1. Assay description** \newline
  **2. Outcome interpretation** \newline
  **3. Results and insights** \newline
- The assay description paragraph must explain the objective of the assay, the experimental system, and methodology. Specify the biological target, assay format (e.g., cell-based, binding), detection method, and any relevant experimental conditions (e.g., temperature, compound concentration).
- The outcome interpretation paragraph must describe how assay outputs are measured and interpreted. Specify how results relate to biological activity or target modulation, the direction of the biological activity (-1 if lower values lead to higher activity e.g., IC50; +1 if higher values result in higher activity e.g., percent. inhibition or effect; 0 if it’s inconclusive, e.g., clearance or solubility), controls, reference compounds, signal thresholds and normalization steps.
- The results and insights paragraph must summarize typical activity ranges, notable behaviors (e.g., agonists, inhibitors), data quality and curation notes that support integration and reproducibility. Highlight meaningful observations from the distribution of activity data. It must be coherent with the outcome interpretation paragraph.
- Separate paragraphs with a single blank line.
- Use only standard ASCII spacing for all numbers, units, and symbols.
- Insert commas in numbers when necessary (e.g., 1,000; 100,000).
- Do not insert non-breaking spaces, narrow spaces, or special typographic characters.
- Do not use tables or hidden formatting.
- Use scientific and formal language.
- Avoid speculation, informal expressions or fabricated information.
- If any relevant data is missing (reported as ‘nan’), state “not reported” rather than inventing details.
- Do NOT write Q&A, lists, bullets, or add any heading besides the three above.
- Do not use thinking mode.

"""

# For each pathogen
for pathogen in pathogens:
    
    # Get assays info
    pathogen_code = str(pathogen.split()[0][0] + pathogen.split()[1]).lower()
    ASSAYS_INFO = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, 'assays.csv'))

    # Create output directory
    PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code)
    os.makedirs(os.path.join(PATH_TO_OUTPUT, "descriptions"), exist_ok=True)

    for i in ASSAYS_INFO[['assay_type', 'assay_organism', 'target_type', 'target_organism', 'activity_type', 'unit', 'activities', 'cpds', 'assay_id']].values:

        assay_id = i[8]
        doc_id = assays[assays['chembl_id'] == assay_id]['doc_id'].tolist()[0]
        if type(i[5]) == str:
            assay_activities = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == i[4]) & (ChEMBL['unit'] == i[5])]["value"].astype(float).tolist()
        else:
            assay_activities = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == i[4]) & (ChEMBL['unit'].isna())]["value"].astype(float).tolist()

        result = {
            "Assay ChEMBL ID": assay_id,
            "Assay type": assay_type_map[i[0]],
            "Assay organism": i[1],
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
            "Target type": i[2],
            "Target organism": i[3],
            "Activity type": i[4],
            "Unit": i[5],
            "Number of activities": i[6],
            "Number of compounds": i[7],
            "Stats": {
                "Percentile 1": round(np.percentile(assay_activities, 1), 3),
                "Percentile 25": round(np.percentile(assay_activities, 25), 3),
                "Mean": round(np.mean(assay_activities), 3),
                "Median": round(np.percentile(assay_activities, 50), 3),
                "Percentile 75": round(np.percentile(assay_activities, 75), 3),
                "Percentile 99": round(np.percentile(assay_activities, 99), 3)
            }
        }

        result = "\n".join([i + ": " + str(result[i]) for i in result])
        USER = f"""Below you will find enumerated annotations from the assay under study.\n\n{result}\n\nUsing the information provided, return a standardized description for the assay."""

        # Print data
        with open(os.path.join(PATH_TO_OUTPUT, "descriptions", f"{assay_id}_{i[4]}_{i[5]}_input.txt"), "w") as f:
            f.write(USER)
        
        # Non streaming call
        response = ollama.generate(model='gpt-oss:20b', prompt=SYSTEM + USER, stream=False, think=True)

        # Print response
        with open(os.path.join(PATH_TO_OUTPUT, "descriptions", f"{assay_id}_{i[4]}_{i[5]}_output.txt"), "w") as f:
            f.write(response.response)

        print(f"✓ Completed {assay_id}")


