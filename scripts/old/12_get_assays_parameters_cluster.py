# Use a conda env with ollama installed
# Run this code on a GPU machine
from collections import Counter
from zipfile import ZipFile, ZIP_DEFLATED
from pydantic import BaseModel
from tqdm import tqdm
import pandas as pd
import numpy as np
import zipfile
import random
import ollama
import json
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH

# List of pathogens to process
pathogens = ["Mycobacterium tuberculosis"]

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

for pathogen in pathogens:

    print(f"Processing pathogen: {pathogen}")
    pathogen_code = get_pathogen_code(pathogen)

    # Get assay data
    ASSAYS = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, "assays.csv"), low_memory=False)[:]

    # Create output directory
    PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code, "parameters")
    os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

    # For each assay
    for assay_id, assay_type, act_type, unit in zip(ASSAYS["assay_id"], ASSAYS["assay_type"], ASSAYS["activity_type"], ASSAYS["unit"]):

        if type(unit) != str:
            unit = 'nan'
        else:
            unit = unit.replace('/', 'FwdS').replace(" ", "__")

        print(f'Assay: {"_".join([assay_id, assay_type, act_type, unit])}')

        # Reading input data file from previous step
        with zipfile.ZipFile(os.path.join(root, "..", "output", pathogen_code, "descriptions", "_".join([assay_id, act_type, unit]) + ".zip"), 'r') as zip_ref:
            with zip_ref.open("_".join([assay_id, act_type, unit]) + "_input.txt") as f:
                input_data = f.read().decode('utf-8')

        PROMPT = f"""
        You are an information extraction assistant specialized in analyzing biochemical data.

        Return ONLY a JSON object with these keys:
        - organism (string)
        - target_type (string)
        - strain (string)
        - mutations (array of strings)
        - known_drug_resistances (array of strings)
        - media (string)

         Rules:

        - If a field is missing or not stated: use "" for strings and [] for arrays.
        - Do not include any other keys or any extra text before or after the JSON.
        - Do not use markdown fences.
        - "Mutations" should include specific genetic variants or engineered changes if mentioned; otherwise [].
        - "Known drug resistances" should list drug resistances of the strain used in the assay; if only general mentions exist, use [].
        - "Media" refers to the growth or culture medium (e.g., Middlebrook 7H9 broth, Lowenstein–Jensen, etc.).

        All available assay annotations are enumerated below:
        {input_data}

        """
        
        class Parameters(BaseModel):
            organism: str
            strain: str
            mutations: list[str]
            known_drug_resistances: list[str]
            media: str
        schema = Parameters.model_json_schema()

        response = ollama.chat(
         messages=[
            {'role': 'user','content': PROMPT}],
        model='gpt-oss:20b',options={'temperature': 0, 'num_ctx': 12288}, keep_alive="1h",
        format=schema)

        # Parse JSON safely
        js = json.loads(response.message['content'])

        # Some validation
        expected = {"organism", "strain", "mutations", "known_drug_resistances", "media"}
        assert set(js.keys()) == expected, f"Unexpected keys: {set(js.keys())}"

        # Add metadata
        js["assay_id"] = assay_id
        js["assay_type"] = assay_type
        js["activity_type"] = act_type
        js["unit"] = unit

        # Write to a JSON file
        out_path = os.path.join(PATH_TO_OUTPUT, "_".join([assay_id, act_type, unit]) + "_parameters.json")
        with open(out_path, "w") as outfile:
            json.dump(js, outfile, indent=2)