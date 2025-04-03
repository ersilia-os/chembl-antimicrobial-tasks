import os
import shutil
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Wrap up tasks and clean output folder')
parser.add_argument('--pathogen_code', type=str, help='Pathogen code')
parser.add_argument('--output_dir', type=str, help='Output folder')
parser.add_argument('--flush', action='store_true', help='Flush output folder')
args = parser.parse_args()

pathogen_code = args.pathogen_code
data_dir = args.output_dir

### SUMMARY FOR ALL TASKS ###

ds = pd.read_csv(os.path.join(data_dir, pathogen_code, "016_selected_tasks.csv"))

if os.path.exists(os.path.join(data_dir, pathogen_code, "017_tasks")):
    shutil.rmtree(os.path.join(data_dir, pathogen_code, "017_tasks"))
os.makedirs(os.path.join(data_dir, pathogen_code, "017_tasks"))

for task in ds["task"].tolist():
    shutil.copy(os.path.join(data_dir, pathogen_code, "013_raw_tasks", task+".csv"), os.path.join(data_dir, pathogen_code, "017_tasks", task+".csv"))

ds.to_csv(os.path.join(data_dir, pathogen_code, "016_tasks_summary.csv"), index=False)


### SUMMARY FOR THE ORGANISM ###

mols_original = pd.read_csv(os.path.join(data_dir, pathogen_code, f"011_{pathogen_code}_original.csv"), low_memory=False)
mols_cleaned = pd.read_csv(os.path.join(data_dir, pathogen_code, f"012_{pathogen_code}_cleaned.csv"), low_memory=False)
num_tasks = pd.read_csv(os.path.join(data_dir, pathogen_code, "014_modelability.csv"), low_memory=False)
selected_tasks = pd.read_csv(os.path.join(data_dir, pathogen_code, "016_selected_tasks.csv"), low_memory=False)
df = pd.DataFrame([[pathogen_code, len(mols_original), len(mols_cleaned), len(num_tasks), len(selected_tasks)]], 
                  columns=['pathogen code', 'num mols original', 'num mols cleaned', 'num tasks', 'num modelable tasks'])

df.to_csv(os.path.join(data_dir, pathogen_code, f"017_{pathogen_code}_summary.csv"), index=False)

if args.flush:
    shutil.rmtree(os.path.join(data_dir, pathogen_code, "013_raw_tasks"))
    os.remove(os.path.join(data_dir, pathogen_code, "016_selected_tasks.csv"))
    os.remove(os.path.join(data_dir, pathogen_code, f"012_{pathogen_code}_cleaned.csv"))
    os.remove(os.path.join(data_dir, pathogen_code, "014_modelability.csv"))