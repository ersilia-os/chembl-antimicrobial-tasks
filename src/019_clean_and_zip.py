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

print(data_dir)
print(output_zip)
print(tasks_dir_dis)

# Create a zip file of the selected tasks
with ZipFile(output_zip, 'w') as zipf:
    for task, sel in zip(selected_tasks["task"], selected_tasks["SELECTED"]):
        if sel == 1:
            filename = f"{task}_{sel}.csv"
            filepath = os.path.join(tasks_dir, filename)
            zipf.write(filepath, arcname=filename)
        else:
            filename = f"{task}_{sel}.csv"
            filepath = os.path.join(tasks_dir_dis, filename)
            zipf.write(filepath, arcname=filename)

# if args.flush:
#     shutil.rmtree(os.path.join(data_dir, pathogen_code, "013_raw_tasks"))
#     os.remove(os.path.join(data_dir, pathogen_code, "016_selected_tasks.csv"))
#     os.remove(os.path.join(data_dir, pathogen_code, f"012_{pathogen_code}_cleaned.csv"))
#     os.remove(os.path.join(data_dir, pathogen_code, "014_modelability.csv"))