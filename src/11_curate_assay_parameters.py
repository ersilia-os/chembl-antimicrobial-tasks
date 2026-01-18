# This script needs to be run with a GPU machine available
from collections import Counter
from pydantic import BaseModel
from zipfile import ZipFile, ZIP_DEFLATED
from tqdm import tqdm
import pandas as pd
import numpy as np
import shutil
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
    target_type_curated: str
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
    PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code, "assay_parameters")
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
        - target_type_curated (string)
        - strain (string)
        - atcc_id (string)
        - mutations (array of strings)
        - known_drug_resistances (array of strings)
        - media (string)

        Return EXACTLY this JSON structure (fill values, keep keys unchanged):

        {{
        "organism": "",
        "target_type_curated": "",
        "strain": "",
        "atcc_id": "",
        "mutations": [],
        "known_drug_resistances": [],
        "media": ""
        }}

        Rules:

        1. Only extract information explicitly stated in the provided assay annotations. Do NOT infer or use external background knowledge.
        2. If a field is missing or not stated, use "" for strings and [] for arrays. Do not use null.
        3. Use double quotes. Output must be valid JSON. Do not include any other keys or any extra text before or after the JSON. Do not use markdown fences.

        Field definitions:

        - "organism": refers to the species/cell line name (e.g., "Mycobacterium tuberculosis", "Homo sapiens"). Do NOT use strain identifiers here.
        - "target_type_curated": refers to target type and should only be modified if target_type is UNCHECKED and the assay annotations clearly indicate that it should be one of: SINGLE PROTEIN or ORGANISM. If there is not enough information/evidence to set it to SINGLE PROTEIN or ORGANISM, leave it UNCHECKED. If the original target_type IS NOT UNCHECKED, target_type_curated should equate to target_type.
        - "strain": refers only to biological strain names (e.g., H37Rv, K12, PAO1). Do NOT include culture collection/catalog identifiers (e.g, ATCC, DSM or NCTC related identifiers or catalog numbers).
        - "atcc_id": refers to the specific ATCC (American Type Culture Collection) identifier, if explicitly stated. Otherwise, leave it empty. Explicitly include "ATCC " before the numerical ID. 
        - "mutations": include ONLY explicit mutations or variants stated in the text. Mutation format MUST be: one-letter amino acid + position (integer) + one-letter amino acid (e.g., "S450L"). If the text uses a longer form (e.g., Ser450Leu) and the conversion is explicit/unambiguous, convert it to the one-letter format.
        - "known_drug_resistances": list drugs for which resistance is explicitly stated (e.g., ["rifampicin", "isoniazid"]). Do NOT infer resistance from mutations. Only include resistances explicitly stated.
        - "media": refers to the growth or culture medium explicitly stated (e.g., Middlebrook 7H9 broth, Lowenstein–Jensen, etc.).

        Important considerations:

        - If assay_type is BINDING and target_type is UNCHECKED, target_type_curated can only be set to UNCHECKED or SINGLE PROTEIN, not ORGANISM.

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
        expected = {"organism", "target_type_curated", "strain", "atcc_id", "mutations", "known_drug_resistances", "media"}
        assert set(js.keys()) == expected, f"Unexpected keys: {set(js.keys())}"

        # Write to a JSON file
        out_path = os.path.join(PATH_TO_OUTPUT, "_".join([ASSAY.assay_id, str(ASSAY.activity_type), str(ASSAY.unit)]) + "_parameters.json")
        with open(out_path, "w") as outfile:
            json.dump(js, outfile, indent=2)

    # Compress all JSON files in a ZIP file
    zip_path = os.path.join(root, "..", "output", pathogen_code, "assay_parameters.zip")
    with ZipFile(zip_path, "w", compression=ZIP_DEFLATED) as zipf:
        for fname in os.listdir(PATH_TO_OUTPUT):
            if fname.endswith(".json"):
                full_path = os.path.join(PATH_TO_OUTPUT, fname)
                # store inside zip without the full directory path
                zipf.write(full_path, arcname=fname)

    # Remove the whole directory after zipping
    shutil.rmtree(PATH_TO_OUTPUT)