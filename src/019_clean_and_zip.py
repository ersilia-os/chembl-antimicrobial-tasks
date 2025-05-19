import os
import shutil
import pandas as pd
import argparse
from zipfile import ZipFile

parser = argparse.ArgumentParser(description='Wrap up tasks and clean output folder')
parser.add_argument("--pathogen_code", type=str)
parser.add_argument("--output_dir", type=str)
parser.add_argument("--organism", action="store_true", help="Flag for organism task")
parser.add_argument("--protein", action="store_true", help="Flag for protein task")
args = parser.parse_args()

pathogen_code = args.pathogen_code
data_dir = args.output_dir

if args.organism == True:
    tasks_dir = os.path.join(data_dir, "013a_raw_tasks_MOD")
    label = "organism"
elif args.protein == True:
    tasks_dir = os.path.join(data_dir, "013b_raw_tasks_MOD")
    label = "protein"

# Load selected tasks
selected_tasks = pd.read_csv(os.path.join(data_dir, "018_selected_tasks_FINAL.csv"))
output_zip = os.path.join(data_dir, f"{pathogen_code}_{label}_tasks.zip")
tasks_dir_dis = os.path.join(data_dir, "015_raw_tasks_DIS")

# Create a zip file of the selected tasks
with ZipFile(output_zip, 'w') as zipf:
    for task, sel in zip(selected_tasks["task"], selected_tasks["SELECTED"]):
        print(task)
        if sel == 1:
            filename = f"{task}_{sel}.csv"
            filepath = os.path.join(tasks_dir, task + ".csv")
            zipf.write(filepath, arcname=filename)
        else:
            filename = f"{task}_{sel}.csv"
            filepath = os.path.join(tasks_dir_dis, task + ".csv")
            zipf.write(filepath, arcname=filename)

# # Clean up the output directory
# shutil.rmtree(os.path.join(data_dir, "013a_raw_tasks_MOD"))
# shutil.rmtree(os.path.join(data_dir, "014_models_MOD"))
# shutil.rmtree(os.path.join(data_dir, "015_models_DIS"))
# shutil.rmtree(os.path.join(data_dir, "015_raw_tasks_DIS"))

# os.remove(os.path.join(data_dir, "014_fingerprints.npy"))
# os.remove(os.path.join(data_dir, "014_models_MOD.csv"))
# os.remove(os.path.join(data_dir, "015_models_DIS.csv"))
# os.remove(os.path.join(data_dir, "018_aggregation_scores.csv"))

