# Use a conda env with ollama installed
# Run this code on a GPU machine
from collections import Counter
from zipfile import ZipFile, ZIP_DEFLATED
from tqdm import tqdm
import pandas as pd
import numpy as np
import zipfile
import random
import ollama
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
    ASSAYS = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, "assays.csv"), low_memory=False)[10:]

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
        Read the assay annotations and return a single CSV line with the following 5 columns, in this exact order and separated by pipes (|):
        - Organism
        - Strain
        - Mutations
        - Known drug resistances
        - Media

        All available assay annotations are enumerated below:
        {input_data}

        Rules for the output:
        - If any field is missing or not stated, leave it blank.
        - "Mutations" should include specific genetic variants or engineered changes if mentioned; otherwise leave blank.
        - "Known drug resistances": list drug resistances of the strain used in the assay; if only general mentions exist, leave blank.
        - "Media" refers to the growth or culture medium (e.g., Middlebrook 7H9 broth, Lowenstein–Jensen, etc.).
        - Output exactly one line, no header, no extra text, no quotes, no trailing pipes (|) and, specially, no tabs nor pipes (|) within individual columns.
        - The final output must have exactly 5 columns in the specified order and, therefore, EXACTLY 4 pipes (|).
        - Triple check that the output format is correct before returning it. An output with less or more than 4 pipes (|) is not valid.

        Examples of VALID outputs:
        Mycobacterium tuberculosis|H37Rv|||Middlebrook 7H9 broth
        Mycobacterium tuberculosis||||

        Examples of INVALID outputs:
        Mycobacterium tuberculosis|H37Rv||||  ← too many pipes (INVALID)
        Mycobacterium tuberculosis|H37Rv| ← too few pipes (INVALID)

        """

        # Non streaming call
        response = ollama.generate(model='gpt-oss:20b', prompt=PROMPT, stream=False, think=False)
        result = response.response.strip().split("|")

        # Check number of columns in response
        if len(result) != 5:
            print(f"Error: Expected 5 columns but got {len(result)}. Response was: {response.response}")
            break

        # Add extra info
        result = [assay_id, assay_type, act_type, unit] + result

        # Save result
        with open(os.path.join(PATH_TO_OUTPUT, "_".join([assay_id, act_type, unit])) + "_parameters.csv", "w") as outfile:
            outfile.write("|".join([str(i) for i in result]))