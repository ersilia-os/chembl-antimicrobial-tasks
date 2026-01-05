from collections import Counter
from pydantic import BaseModel
from zipfile import ZipFile, ZIP_DEFLATED
from tqdm import tqdm
import pandas as pd
import numpy as np
import ollama
import json
import sys
import os

import subprocess
print("=== Python GPU Check ===\n\n")
try:
    result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
    print(result.stdout)
except Exception as e:
    print(f"nvidia-smi failed: {e}\n\n")

# Define root directory
# root = os.path.dirname(os.path.abspath(__file__))
root = "."
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH

# Load assays and docs information
assays = pd.read_csv(os.path.join(CONFIGPATH, "chembl_activities", "assays.csv"), low_memory=False)
docs = pd.read_csv(os.path.join(CONFIGPATH, "chembl_activities", "docs.csv"), low_memory=False)
assay_type_map = {"F": "Functional", "B": "Binding", "T": "Toxicity", "A": "ADME", "P": "Physicochemical", "U": "Uncategorized"}

# List of pathogens
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]
pathogens = ["Acinetobacter baumannii", "Mycobacterium tuberculosis", "Klebsiella pneumoniae"]

class Parameters(BaseModel):
    organism: str
    target_type: str
    strain: str
    atcc_id: str
    mutations: list[str]
    known_drug_resistances: list[str]
    media: str

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

for pathogen in pathogens:

    # Creating output directory
    print(f"Processing pathogen: {pathogen}...")
    pathogen_code = get_pathogen_code(pathogen)
    PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code, "parameters")
    os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

    # Loading assay data
    ASSAYS = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, "assays_cleaned.csv"))

    # For each assay
    for idx, ASSAY in tqdm(ASSAYS.iterrows()):

        # Getting doc_id
        doc_id = assays[assays['chembl_id'] == ASSAY.assay_id]['doc_id'].tolist()[0]

        result = {
            "Assay ChEMBL ID": ASSAY.assay_id,
            "Assay type": assay_type_map[ASSAY.assay_type],
            "Assay organism": ASSAY.assay_organism,
            "Assay description": assays[assays['chembl_id'] == ASSAY.assay_id]['description'].tolist()[0],
            "Assay strain": assays[assays['chembl_id'] == ASSAY.assay_id]['assay_strain'].tolist()[0],
            "Assay category": assays[assays['chembl_id'] == ASSAY.assay_id]['assay_category'].tolist()[0],
            "Assay test type": assays[assays['chembl_id'] == ASSAY.assay_id]['assay_test_type'].tolist()[0],
            "Assay cell type": assays[assays['chembl_id'] == ASSAY.assay_id]['assay_cell_type'].tolist()[0],
            "Document title": docs[docs['doc_id'] == doc_id]['title'].tolist()[0],
            "Document abstract": docs[docs['doc_id'] == doc_id]['abstract'].tolist()[0],
            "Document journal": docs[docs['doc_id'] == doc_id]['journal'].tolist()[0],
            "Document PubMed ID": docs[docs['doc_id'] == doc_id]['pubmed_id'].tolist()[0],
            "Document DOI": docs[docs['doc_id'] == doc_id]['doi'].tolist()[0],
            "Target type": ASSAY.target_type,
            "Target organism": ASSAY.target_organism,
            "Activity type": ASSAY.activity_type,
            "Unit": ASSAY.unit,
            "Number of activities": ASSAY.activities,
            "Number of activities with nan value": ASSAY.nan_values,
            "Number of compounds": ASSAY.cpds,
            "Direction of biological activity:": ASSAY.direction, 
            "Number of compounds being active/inactive according to activity_comment": ASSAY.activity_comment_counts,
            "Number of compounds being active/inactive according to standard_text": ASSAY.standard_text_count,
        }

        result = "\n\t".join(["- " + i + ": " + str(result[i]) for i in result])

        PROMPT = f"""
        You are an information extraction assistant specialized in analyzing biochemical data from ChEMBL. Below, you will find a set of assay annotations from a single ChEMBL assay under study. 

        Your job is to return ONLY a JSON object with these keys:
        - organism (string)
        - target_type (string)
        - strain (string)
        - atcc_id (string)
        - mutations (array of strings)
        - known_drug_resistances (array of strings)
        - media (string)

        Return EXACTLY this JSON structure (fill values, keep keys unchanged):

        {{
        "organism": "",
        "target_type": "",
        "strain": "",
        "atcc_id": "",
        "mutations": [],
        "known_drug_resistances": [],
        "media": ""
        }}

            Rules:

        - "organism" refers to the particular organism under study in the assay. If already specified and coherent with the rest of the data, leave it as is.
        - “organism” should be the species/cell line name (e.g., Mycobacterium tuberculosis, Homo sapiens), NOT the strain identifier.
        - "target_type" should only be modified if its current value is UNCHECKED and the assay annotations clearly indicate that it should be one of: SINGLE PROTEIN, CELL-LINE, or ORGANISM.
        - "strain" refers to the particular strain under study in the assay.
        - "strain" refers only to biological strain names (e.g., H37Rv, K12, PAO1). Do NOT include culture collection/catalog identifiers (e.g, ATCC, DSM or NCTC related identifiers or catalog numbers).
        - "atcc_id" refers to the specific ATCC (American Type Culture Collection) identifier, if provided. Otherwise, leave it empty. 
        - "mutations" should include specific genetic variants or engineered changes if mentioned; otherwise [].
        - "known_drug_resistances" should list drug resistances of the strain used in the assay; if only general mentions exist, use [].
        - "media" refers to the growth or culture medium (e.g., Middlebrook 7H9 broth, Lowenstein–Jensen, etc.).
        - If a field is missing or not stated: use "" for strings and [] for arrays.
        - Do not include any other keys or any extra text before or after the JSON.
        - Do not use markdown fences.
        - Must be valid JSON
        - Use double quotes
        - Always include all keys
        - Only extract information explicitly stated in the provided assay annotations. Do not infer or use external background knowledge.

            Assay annotations:

        {result}"""
        
        schema = Parameters.model_json_schema()
        response = ollama.chat(
            messages=[
            {'role': 'user','content': PROMPT}],
        model='gpt-oss:20b',options={'temperature': 0, 'num_ctx': 12288}, keep_alive="1h",
        format=schema)

        # Parse JSON safely
        js = json.loads(response.message['content'])

        # Some validation
        expected = {"organism", "target_type", "strain", "atcc_id", "mutations", "known_drug_resistances", "media"}
        assert set(js.keys()) == expected, f"Unexpected keys: {set(js.keys())}"

        # Add metadata
        js["assay_id"] = ASSAY.assay_id
        js["assay_type"] = ASSAY.assay_type
        js["activity_type"] = ASSAY.activity_type
        js["unit"] = ASSAY.unit

        # Write to a JSON file
        out_path = os.path.join(PATH_TO_OUTPUT, "_".join([ASSAY.assay_id, str(ASSAY.activity_type), str(ASSAY.unit)]) + "_parameters.json")
        with open(out_path, "w") as outfile:
            json.dump(js, outfile, indent=2)