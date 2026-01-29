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
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH

# Load pathogen info
pathogen_code = sys.argv[1]
df = pd.read_csv(os.path.join(CONFIGPATH, 'pathogens.csv'))
row = df.loc[df["code"].eq(pathogen_code)]
if row.empty: 
    raise SystemExit(f"Unknown code: {pathogen_code}")
pathogen = row.iloc[0]["pathogen"]

print("Step 11 - Needs to run on a GPU-enabled machine.")

# Set path to output
PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code)

# Load assays and docs information
assays = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "assays.csv"), low_memory=False)
docs = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "docs.csv"), low_memory=False)
assay_type_map = {"F": "Functional", "B": "Binding", "T": "Toxicity", "A": "ADME", "P": "Physicochemical", "U": "Uncategorized"}

class Parameters(BaseModel):
    organism_curated: str
    target_type_curated: str
    target_name_curated: str
    target_chembl_id_curated: str
    strain: str
    atcc_id: str
    mutations: list[str]
    known_drug_resistances: list[str]
    media: str

COLS = ["assay_id", "activity_type", "unit", "organism_curated", "target_type_curated", "target_name_curated", "target_chembl_id_curated", 
        "strain", "atcc_id", "mutations", "known_drug_resistances", "media"]

# Creating output directory
print(f"Processing pathogen: {pathogen}...")

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
        "Number of compounds being active/inactive according to text comments": ASSAY.act_flag + ASSAY.inact_flag,
    }

    result = "\n\t".join(["- " + i + ": " + str(result[i]) for i in result])

    PROMPT = f"""
    You are an information extraction assistant specialized in analyzing biochemical data from ChEMBL associated to pathogens of global health concern. 
    Below, you will find a set of assay annotations from a single ChEMBL assay under study. 

    Your job is to return ONLY a JSON object with these keys:
    - organism_curated (string)
    - target_type_curated (string)
    - target_name_curated (string)
    - target_chembl_id_curated (string)
    - strain (string)
    - atcc_id (string)
    - mutations (array of strings)
    - known_drug_resistances (array of strings)
    - media (string)

    Return EXACTLY this JSON structure (fill values, keep keys unchanged):

    {{
    "organism_curated": "",
    "target_type_curated": "",
    "target_name_curated": "",
    "target_chembl_id_curated": "",
    "strain": "",
    "atcc_id": "",
    "mutations": [],
    "known_drug_resistances": [],
    "media": ""
    }}

    Rules:

    1. Only extract information explicitly stated in the provided assay annotations. For all fields except target_type_curated, only extract what is explicitly 
        stated. For target_type_curated, you MAY classify using the mapping rules below based on the provided annotations.
    2. If a field is missing or not stated, use "" for strings and [] for arrays. Do not use null.
    3. Use double quotes. Output must be valid JSON. Do not include any other keys or any extra text before or after the JSON. Do not use markdown fences.
    4. Do not output placeholders such as "unknown", "not stated", "not specified", or "N/A". Use "" or [] as specified.
    5. Do not include any comma characters (",") inside string values. If the source text contains commas, replace them with semicolons (";").

    Field definitions:

    - "organism_curated": refers to the biological species explicitly stated (e.g., "Mycobacterium tuberculosis", "Homo sapiens"). 
    Do NOT include strain identifiers.
    - "target_type_curated": refers to the type of target. Target types you might find: SINGLE PROTEIN, ORGANISM, CELL-LINE, PROTEIN COMPLEX, PROTEIN-PROTEIN INTERACTION, 
    PROTEIN FAMILY, TISSUE, SELECTIVITY GROUP, NUCLEIC-ACID, PROTEIN COMPLEX GROUP, SMALL MOLECULE, CHIMERIC PROTEIN, OLIGOSACCHARIDE, UNKNOWN, SUBCELLULAR, MACROMOLECULE, 
    PROTEIN NUCLEIC-ACID COMPLEX, LIPID, METAL, 3D CELL CULTURE, PHENOTYPE, NON-MOLECULAR, ADMET, UNCHECKED, NO TARGET. Your job is to assign unlabeled assays with their
    corresponding target type, following the mapping rules stated below. 

        if target_type == UNCHECKED --> target_type_curated: SINGLE PROTEIN, ORGANISM or DISCARDED
        if target_type == NON-MOLECULAR --> target_type_curated: ORGANISM or DISCARDED
        if target_type != UNCHECKED AND target_type != NON-MOLECULAR --> target_type_curated: [same as target_type] OR DISCARDED
    
    
    So, "target_type_curated" MUST be one of: DISCARDED, SINGLE PROTEIN, ORGANISM, or exactly the provided "Target type" string (verbatim).
    You should set "target_type_curated" to "DISCARDED" if (a) the annotations explicitly indicate the assay is not a bioactivity assay (e.g., transcriptomics), or (b) 
    the annotations do not explicitly provide enough information to assign any assayed biological entity, meaning both of the following are missing:
        --> an explicit target name or target identifier (e.g., protein name, target ChEMBL ID)
        --> an explicit assayed organism or cell-line entity.

    - "target_name_curated": refers to the target name explicitly stated in the annotations (e.g., "Mycobacterium tuberculosis" or "Lysine--tRNA ligase"). 
    If none is explicitly stated, leave it empty. If the assay target is a protein/complex, do NOT set target_name_curated to the organism name.
    - "target_chembl_id_curated": if provided in assay annotations, the target ChEMBL ID (e.g., CHEMBL360 for Mycobacterium tuberculosis or 
    CHEMBL3301561 for Lysine--tRNA ligase from Plasmodium falciparum). Otherwise leave it empty.
    - "strain": refers only to biological strain names (e.g., H37Rv, K12, PAO1). Do NOT include culture collection/catalog identifiers (e.g, ATCC, DSM or NCTC 
    related identifiers or catalog numbers).
    - "atcc_id": refers to the specific ATCC (American Type Culture Collection) identifier, if explicitly stated. Otherwise, leave it empty. 
    Extract only the numeric part and format exactly as "ATCC <number>". 
    - "mutations": include ONLY explicit mutations or variants stated in the text. Mutation format MUST be: one-letter amino acid + position (integer) + 
    one-letter amino acid (e.g., "S450L"). If the text uses a longer form (e.g., Ser450Leu) and the conversion is explicit/unambiguous, convert it to the one-letter format.
    - "known_drug_resistances": list drugs for which resistance is explicitly stated (e.g., ["rifampicin", "isoniazid"]). Do NOT infer resistance from mutations. 
    Only include resistances explicitly stated.
    - "media": refers to the growth or culture medium explicitly stated (e.g., Middlebrook 7H9 broth, Lowenstein–Jensen, etc.).

    Important considerations:

    - If assay_type is Binding and target_type is UNCHECKED, then target_type_curated can only be set to DISCARDED or SINGLE PROTEIN, not ORGANISM.

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
    COLS = ["organism_curated", "target_type_curated", "target_name_curated", "target_chembl_id_curated",
                "strain", "atcc_id", "mutations", "known_drug_resistances", "media"]
    assert set(js.keys()) == set(COLS), f"Unexpected keys: {set(js.keys())}"

    if idx == 0:
        with open(os.path.join(PATH_TO_OUTPUT, 'assays_parameters.csv'), "w") as file:
            file.write(",".join(["assay_id", "activity_type", "unit"] + COLS))
            file.write("\n")

    with open(os.path.join(PATH_TO_OUTPUT, 'assays_parameters.csv'), "a") as file:
            if type(ASSAY.unit) == str:
                unit = ASSAY.unit
            else:
                unit = ""
            file.write(",".join([ASSAY.assay_id, ASSAY.activity_type, unit] + 
                                [js[i] if type(js[i]) == str else ";".join(js[i]) for i in COLS]))
            file.write("\n")