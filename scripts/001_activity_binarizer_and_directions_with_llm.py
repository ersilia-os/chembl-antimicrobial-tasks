import os
import sys
import pandas as pd
import collections
import random
import subprocess
import requests
import time
from tqdm import tqdm
# from llama_index.core.llms import ChatMessage, MessageRole
# from llama_index.llms.llamafile import Llamafile

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH

# from default import CONFIGPATH, LLM_BIN_FILENAME
# llm_bin_path = os.path.abspath(os.path.join(root, "..", "bin", LLM_BIN_FILENAME))
# print(f"LLM bin path: {llm_bin_path}")
# llm_model = Llamafile()

# --- REPLACE the llamafile import + instantiation with this block ----------------
# from llama_index.llms.llamafile import Llamafile   # ← remove this

OLLAMA_URL   = os.getenv("OLLAMA_URL", "http://localhost:11434/api/generate")
OLLAMA_MODEL = os.getenv("OLLAMA_MODEL", "gemma3:4b")  # OLLAMA_DEBUG=1 ollama serve

def classify_activity_comments(input_file, file_path, rewrite=False):
    """
    Classify activity comments as active (1), inactive (-1) or inconclusive (0).
    """
    df = pd.read_csv(input_file, dtype=str)[:55]
    columns = list(df.columns)
    if "activity_classified" not in columns:
        df["activity_classified"] = ""
        df.to_csv(file_path, index=False)
    if rewrite:
        df["activity_classified"] = ""
        df.to_csv(file_path, index=False)
    texts = df["activity_comment"].tolist()
    classifications_ = df["activity_classified"].tolist()
    classifications = []
    for classification in classifications_:
        if str(classification) == "nan":
            classifications.append("")
        else:
            classifications.append(str(classification))
    i = 0
    current_classifications = [x for x in classifications]
    print(len(texts), len(classifications))
    for comment, classification in tqdm(zip(texts, classifications)):
        comment = str(comment)
        if comment == "nan":
            current_classifications[i] = "0"
            i += 1
            continue
        if comment.isnumeric():
            current_classifications[i] = "0"
            i += 1
            continue
        if comment.startswith("Median N="):
            current_classifications[i] = "0"
            i += 1
            continue
        if "no significant biotransformation" in comment.lower():
            current_classifications[i] = "-1"
            i += 1
            continue
        if ". significant biotransformation" in comment.lower():
            current_classifications[i] = "1"
            i += 1
            continue
        if comment.startswith("inhibitor ["):
            current_classifications[i] = "1"
            i += 1
            continue
        if comment.startswith("substrate ["):
            current_classifications[i] = "0"
            i += 1
            continue
        if "putative" in comment.lower() and "is biotransformed by" in comment.lower():
            current_classifications[i] = "1"
            i += 1
            continue
        if "biotransformation" in comment.lower() and "could not be" in comment.lower():
            current_classifications[i] = "0"
            i += 1
            continue
        if "biotransformation" in comment.lower() and "was proven to be" in comment.lower():
            current_classifications[i] = "1"
            i += 1
            continue
        if classification != "":
            i += 1
            continue

        
        prompt_text = (
            f"You are a biomedical data curator from ChEMBL. Your task is to classify comments about the bioactivity of compounds.\n"
            f"Activity comments should be classified as 1. This includes terms like 'active', 'significant activity', 'inhibition', 'binding', 'toxic' or any type of 'toxicity', 'agonist', 'antagonist', 'weakly active', 'positive', 'substantial inhibition' or any evidence that the compound has an observable biological effect etc.\n"
            f"Inactivity or no activity comments should be classified as -1. This includes terms like 'inactive', 'no activity', 'not active', 'no inhibition', 'non-toxic', 'negative', absence of an observable biological effect etc.\n"
            f"Inconclusive or irrelevant comments should be classified as 0. This includes empty text, scores or numbers missing units, terms like 'inconclusive', 'not measured', etc., terms that simply explain the type of measurement (e.g. '% of inhibition'), terms that mention statistics from the experiment (e.g. mean, median, summarised, no. measurements) and terms that seem unrelated to the measured activity.\n"
            f"Simply answer 1, 0 or -1. Do not explain anything. Just give the numerical value.\n"
            f"This is the comment you have to curate:"
            f"- {comment}\n"
        )

        payload = {
        "model": OLLAMA_MODEL,
        "prompt": prompt_text,
        "stream": False,
        "options": {"temperature": 0}}

        r = requests.post(OLLAMA_URL, json=payload, timeout=120)
        r.raise_for_status()
        prediction = (r.json().get("response") or "").strip()
        prediction = prediction.replace("<|eot_id|>", "").replace("<end_of_turn>", "").strip()

        print(f"Prediction: {prediction} / Comment: {comment}")
        if prediction == "1":
            current_classifications[i] = "1"
        elif prediction == "-1":
            current_classifications[i] = "-1"
        elif prediction == "0":
            current_classifications[i] = "0"
        else:
            current_classifications[i] = "To resolve"
        i += 1
        df["activity_classified"] = current_classifications
        df.to_csv(file_path, index=False)
    df["activity_classified"] = current_classifications
    df.to_csv(file_path, index=False)

def classify_activity_standards_with_direction(input_file, file_path, rewrite=False):
    """
    Classify activity standards with the right direction (-1: the lower the more active, 1: the higher the more active, 0: inconclusive).
    """
    df = pd.read_csv(input_file, dtype=str)
    columns = list(df.columns)
    if "activity_direction" not in columns:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    if rewrite:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    current_directions_ = df["activity_direction"].tolist()
    current_directions = []
    for direction in current_directions_:
        if str(direction) == "nan":
            current_directions += [""]
        else:
            current_directions += [str(direction)]
    for i, v in enumerate(df.values):
        std_type = v[1]
        definition = v[2]
        std_unit = v[3]
        if current_directions[i] != "":
            continue
        prompt_text = (   ### REFINE PROMPT, PARTICULARLY THE FIRST SENTENCE
            f"You are a biomedical data curator from ChEMBL. Your task is to classify comments types of activity measures of bioactivity assays.\n"
            f"When the lower the value, the higher the activity (including toxicity), classify as -1. For example, in an IC50 measure, smaller concentrations indicate more potency. Same for smaller concentrations leading to higher toxicity and growth inhibition, inhibition constants (Ki) etc.\n"
            f"When the higher the value, the higher the activity (including toxicity), classify as 1. For example, in a percentage inhibition, larger values indicate more potency. Likewise, high zone of inhibition, etc.\n"
            f"When the direction is inconclusive, classify as 0. This includes assays related to pharmacokinetics in the human body (e.g. plasma protein binding, bioavailability, fraction unbound in plasma, etc.), or physicochemical measurements not related to antipathogen activity data (e.g., any kind of solubility).\n"
            f"Simply answer 1, 0 or -1. Do not explain anything. Just give the numerical value.\n"
            f"This is the activity type and units you have to curate:"
            f"- {std_type} : {definition}. Units: {std_unit}\n"
        )
        payload = {
        "model": OLLAMA_MODEL,
        "prompt": prompt_text,
        "stream": False,
        "options": {"temperature": 0}}

        r = requests.post(OLLAMA_URL, json=payload, timeout=120)
        r.raise_for_status()
        prediction = (r.json().get("response") or "").strip()
        prediction = prediction.replace("<|eot_id|>", "").replace("<end_of_turn>", "").strip()

        print(f"Prediction: {prediction} / {std_type} : {definition}. Units: {std_unit}")
        if prediction == "1":
            current_directions[i] = "1"
        elif prediction == "-1":
            current_directions[i] = "-1"
        elif prediction == "0":
            current_directions[i] = "0"
        else:
            current_directions[i] = "0"
        df["activity_direction"] = current_directions
        df.to_csv(file_path, index=False)
    df["activity_direction"] = current_directions
    df.to_csv(file_path, index=False)

def get_sample_assay_descriptions_for_standard_units():
    print("Getting assay descriptions...")
    df_ass = pd.read_csv(os.path.join(CONFIGPATH,"chembl_activities", "assay_descriptions.csv"), low_memory=False)
    assay_descs = {}
    for v in df_ass[["assay_id", "description"]].values:
        assay_descs[v[0]] = v[1]
    print("Getting activities")
    df_act = pd.read_csv(os.path.join(CONFIGPATH, "chembl_activities","activities.csv"), low_memory=False)
    print("Done reading activities")
    std_unit_assays = collections.defaultdict(list)
    for v in tqdm(df_act[["assay_id", "standard_type", "standard_units"]].values):
        std_unit_assays[(v[1], v[2])] += [v[0]]
    R = []
    print(std_unit_assays[('Potency', "nM")])
    for k,v in tqdm(std_unit_assays.items()):
        random.shuffle(v)
        v = [assay_descs[x] for x in v[:10]]
        v = [x for x in v if str(v) != "" and str(v) != "nan"]
        v = v[:3]
        if len(v) == 2:
            v = v + [""]
        elif len(v) == 1:
            v = v + ["", ""]
        else:
            pass
        R += [[k[0], k[1], v[0], v[1], v[2]]]
    df = pd.DataFrame(R, columns = ["standard_type", "standard_unit", "assay_description_1", "assay_description_2", "assay_description_3"])
    df.to_csv(os.path.join(CONFIGPATH,"llm_processed", "activity_std_units_with_3_assay_descriptions.csv"), index=False)

def classify_all_activity_standards_with_direction(input_file, file_path, rewrite=False):
    """
    Classify all activity standard units with the right direction (-1: the lower the more active, 1: the higher the more active, 0: inconclusive).
    """
    df = pd.read_csv(input_file, dtype=str)[:1]
    da = pd.read_csv(os.path.join(CONFIGPATH,"llm_processed", "activity_std_units_with_3_assay_descriptions.csv"), dtype=str)
    std_unit_descriptions = {}
    for v in da.values:
        std_unit_descriptions[(v[0], v[1])] = [v[2], v[3], v[4]]
    columns = list(df.columns)
    if "activity_direction" not in columns:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    if rewrite:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    current_directions_ = df["activity_direction"].tolist()
    current_directions = []
    for direction in current_directions_:
        if str(direction) == "nan":
            current_directions += [""]
        else:
            current_directions += [str(direction)]
    for i, v in tqdm(enumerate(df.values)):
        std_type = v[0]
        std_unit = v[1]
        descriptions = std_unit_descriptions[(std_type, std_unit)]
        if current_directions[i] != "":
            continue
        prompt_text = (
            f"You are a biomedical data curator from ChEMBL. Your task is to classify types of bioactivity measures\n"
            f"You will be provided with an activity type and units, and three assay descriptions where those units are used. You can use these assay descriptions to get hints about the activity type and units. *The most important determinant for your decision should be the activity type and unit.*\n"
            f"Assign -1 when:\n"
            f"- The lower the value, the higher the bioactivity. That is, smaller values indicate more potency.\n"
            f"- Examples: IC50 (uM, nM, etc.), MIC, MIC90, LD50, CC50 etc.\n"
            f"Assign 1 when:\n"
            f"- The higher the value, the higher the activity. That is, larger values indicate more potency.\n"
            f"- Examples: % inhibition, activity (percentage or proportion (no unit)), IZ (zone of inhibition, in mm), etc.\n"
            f"Assign 0 when:\n"
            f"- The direction is inconclusive. That is, it is not clear if higher or lower values are preferrable\n"
            f"- Examples: logP, other physicochemical measurements, some pharmacokinetic properties, etc.\n"
            f"Attention: Check the units. Sometimes the same activity type can have different directions depending on the units. For example, -log (minus log, as in pIC50) will reverse the IC50 activity trend from -1 to 1.\n"
            f"Here are some assay descriptions you can use to know more about the activity type unit:\n"
            f"1. {descriptions[0]}\n"
            f"2. {descriptions[1]}\n"
            f"3. {descriptions[2]}\n"
            f"This is the activity type and units you have to classify:"
            f"- Activity type: {std_type}. Unit: {std_unit}\n"
            f"Simply answer 1, 0 or -1. Do not explain anything. Do not give more than one value. Just give the numerical value.\n"
        )

        print(prompt_text)
        payload = {
        "model": OLLAMA_MODEL,
        "prompt": prompt_text,
        "stream": False,
        "options": {"temperature": 0}}

        r = requests.post(OLLAMA_URL, json=payload, timeout=120)
        r.raise_for_status()
        prediction = (r.json().get("response") or "").strip()
        prediction = prediction.replace("<|eot_id|>", "").replace("<end_of_turn>", "").strip()


        print(f"Prediction: {prediction} / {std_type} Units: {std_unit} {descriptions[0]}.")
        if prediction == "1":
            current_directions[i] = "1"
        elif prediction == "-1":
            current_directions[i] = "-1"
        elif prediction == "0":
            current_directions[i] = "0"
        else:
            print("To revise")
            print(prediction)
            current_directions[i] = "To revise"
        df["activity_direction"] = current_directions
        df.to_csv(file_path, index=False)
    df["activity_direction"] = current_directions
    df.to_csv(file_path, index=False)

def process_standard_text(input_file, file_path):
    df = pd.read_csv(input_file, dtype=str)
    manual_dict = {"Active".lower():1,  
                "Compound NOT metabolized".lower():-1,
                "Compound metabolized".lower():1,
                "Could not be calculated".lower():0,
                "IC50 could not be calculated".lower():0,
                "Inactive".lower():-1,
                "Linear Fit".lower():0,
                "NOT MEASURED".lower():0,
                "Nil".lower():-1,
                "No binding".lower():-1,
                "Not Active".lower():-1,
                "Not Determined".lower():0
                }
    st_text = df["activity_comment"].tolist()
    final_dict = {}
    for text in st_text:
        if str(text).lower() in manual_dict:  
            final_dict[text] = manual_dict[str(text)]
        else:
            print("This comment is not processed: ", text)
    df_ = pd.DataFrame(final_dict.items(), columns=["standard_text_value", "standard_text_classification"])
    df_.to_csv(file_path, index=False)

def assess_overlap_comments(input_file, file_path):  #rewrite this function

    df1 = pd.read_csv(input_file)
    df2 = pd.read_csv(file_path)
    df1 = df1[:len(df2)]
    assert df1["activity_comment"].equals(df2["activity_comment"])
    matches = [i == j for i,j in zip(df1['activity_classified'], df2['manual_curation'])]
    print(f"{sum(matches) / len(matches) * 100} % of matches with manually curated annotations")

def assess_overlap_act_type_units(input_file, file_path):  # as well as this one

    df1 = pd.read_csv(input_file)
    df2 = pd.read_csv(file_path)
    df1 = df1[:len(df2)]
    assert df1["activity_comment"].equals(df2["activity_comment"])
    matches = [i == j for i,j in zip(df1['activity_classified'], df2['manual_curation'])]
    print(f"{sum(matches) / len(matches) * 100} % of matches with manually curated annotations")


if __name__ == "__main__":

    # print("Starting with activity comments classification...")
    # input_file = os.path.join(CONFIGPATH, "chembl_activities", "activity_comments.csv")
    # os.makedirs(os.path.join(CONFIGPATH,"llm_processed"), exist_ok=True)
    # file_path= os.path.join(CONFIGPATH,"llm_processed", "activity_comments.csv")
    # manual_curation_file = os.path.join(CONFIGPATH, "manual_curation", "activity_comments_manual_curation.csv")
    # classify_activity_comments(input_file, file_path, rewrite=False)
    # print("Activity comments classification done.")
    # assess_overlap_comments(file_path, manual_curation_file)

    # print("Starting with activity standards classification (direction)...")
    # input_file = os.path.join(CONFIGPATH, "chembl_activities", "activity_stds_lookup.csv")
    # file_path = os.path.join(CONFIGPATH,"llm_processed", "activity_stds_lookup.csv")
    # classify_activity_standards_with_direction(input_file,file_path, rewrite=False)
    # print("Activity standards classification done.")

    # print("Getting sample assay descriptions for standard units...")
    # get_sample_assay_descriptions_for_standard_units()
    # print("Sample assay descriptions for standard units done.")

    print("Starting with all activity standards classification (direction)...")
    input_file = os.path.join(CONFIGPATH, "chembl_activities","activity_std_units.csv")
    file_path = os.path.join(CONFIGPATH,"llm_processed", "activity_std_units.csv")
    classify_all_activity_standards_with_direction(input_file, file_path, rewrite=False)
    manual_curation_file = os.path.join(CONFIGPATH, "manual_curation", "activity_std_units_manual_curation.csv")
    assess_overlap_comments(file_path, manual_curation_file)
    # input_file = os.path.join(CONFIGPATH, "chembl_activities","standard_text.csv")
    # file_path = os.path.join(CONFIGPATH,"llm_processed", "standard_text.csv")
    # # process_standard_text()