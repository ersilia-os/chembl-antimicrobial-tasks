from scipy.stats import spearmanr, kendalltau, pearsonr
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import joblib
import argparse

parser = argparse.ArgumentParser(f"Select non redundant tasks")
parser.add_argument("--pathogen_code", type=str)
parser.add_argument("--output_dir", type=str)
parser.add_argument("--organism", action="store_true", help="Flag for organism task")
parser.add_argument("--protein", action="store_true", help="Flag for protein task")
args = parser.parse_args()


pathogen_code = args.pathogen_code
data_dir = args.output_dir

if args.organism == True:
    tasks_dir = os.path.join(data_dir, "013a_raw_tasks_MOD")
elif args.protein == True:
    tasks_dir = os.path.join(data_dir, "013b_raw_tasks_MOD")

selected_tasks = pd.read_csv(os.path.join(data_dir, "017_selected_tasks_final.csv"))
tasks = selected_tasks['task'].tolist()
sel_reason = selected_tasks['SELECTED'].tolist()


def get_AGG_SCORE(tasks_dir, tasks):

    INACTIVES, AGG_SCORE = set(), {}

    for task in tasks:

        # Load data
        data = pd.read_csv(os.path.join(tasks_dir, task + ".csv"))
        activity = data.columns[-1]

        # Store actives in dict
        actives = sorted((set(data[data[activity] == 1]['inchikey'])))
        for act in actives:
            if act not in AGG_SCORE:
                AGG_SCORE[act] = []
            AGG_SCORE[act].append(task)

        # Store inactives in a set
        inactives = set(data[data[activity] == 0]['inchikey'])
        INACTIVES = INACTIVES.union(inactives)

    # Add inactives in AGG_SCORE - only if not already there
    for inact in INACTIVES:
        if inact not in AGG_SCORE:
            AGG_SCORE[inact] = []

    return AGG_SCORE

def calculate_reference_set(df, N=10000):
    df = df.copy()[['inchikey', 'num_tasks']]
    if len(df) < N:
        return df
    else:
        positives = df[df['num_tasks'] != 0][: int(N/2)]
        negatives = df[df['num_tasks'] == 0]
        if len(negatives) > int(N/2):
            negatives = negatives.sample(n=int(N/2), random_state=42)
        return pd.concat([positives, negatives])
    
def load_fingerprints(data_dir):
    X = np.load(os.path.join(data_dir, "014_fingerprints.npy"))
    with open(os.path.join(data_dir, "014_fingerprints_inchikeys.txt"), "r") as f:
        keys = f.read().splitlines()
    return X, keys

print("Calculating AGG Score")

# Get AGG Scores
AGG_SCORE = get_AGG_SCORE(tasks_dir, tasks)

# Sort by AGG Score
compounds_sorted = sorted(AGG_SCORE, key=lambda x: len(AGG_SCORE[x]), reverse=True)

# Store in a df
df = []
for cpd in compounds_sorted:
    df.append([cpd, len(AGG_SCORE[cpd]), ";".join(AGG_SCORE[cpd])]) 

df = pd.DataFrame(df, columns=['inchikey', 'num_tasks', 'tasks'])
df.to_csv(os.path.join(data_dir, '018_aggregation_scores.csv'), index=False)

print("Generating reference set")

# Generate reference set
ref_set = calculate_reference_set(df, 10000)

# Storing reference set
ref_set.to_csv(os.path.join(data_dir, '018_ref_set.csv'), index=False)  # Maybe not needed

# Getting reference IKs
ref_IK = ref_set["inchikey"].tolist()

# Loading fingerprints
X, inchikeys = load_fingerprints(data_dir)

# Preparing submatrix
indices = {}
for i, ik in enumerate(inchikeys):
    indices[ik] = i
idxs = [indices[ik] for ik in ref_IK]
X = X[idxs]

# Check that models are overfitted
auroc_cutoff = 0.99
models_MOD = pd.read_csv(os.path.join(data_dir, "014_models_MOD.csv"))['auroc'].tolist()
models_DIS = pd.read_csv(os.path.join(data_dir, "015_models_DIS.csv"))['auroc'].tolist()
print(f"{len([i for i in models_MOD if i > auroc_cutoff]) / len(models_MOD) * 100}% of the full models_MOD have AUROC > {auroc_cutoff}")
print(f"{len([i for i in models_DIS if i > auroc_cutoff]) / len(models_DIS) * 100}% of the full models_DIS have AUROC > {auroc_cutoff}")


print("Making predictions on the reference set of compounds")

PREDICTIONS = {}

# For each task
for task, reason in tqdm(zip(tasks, sel_reason)):

    # Load overfitted models
    if reason == 1:
        model = joblib.load(os.path.join(data_dir, "014_models_MOD", task + ".joblib"))
    elif reason == 2 or reason == 3:
        model = joblib.load(os.path.join(data_dir, "015_models_DIS", task + ".joblib"))

    # Make predictions and store them
    PREDICTIONS[task] = model.predict_proba(X)[:,1]


CORRELATIONS = []

# Calculate correlations between each pair of tasks
for c, task1 in enumerate(tasks):
    for task2 in tasks[c:]:
        p1 = PREDICTIONS[task1]
        p2 = PREDICTIONS[task2]
        s = spearmanr(p1, p2)
        p = pearsonr(p1, p2)
        k = kendalltau(p1, p2)
        CORRELATIONS.append([task1, task2, s.statistic, s.pvalue, p.statistic, p.pvalue, k.statistic, k.pvalue])

# Store results
CORRELATIONS = pd.DataFrame(CORRELATIONS, columns=['task1', 'task2', 'Spearman', 'Spearman pvalue',
                                                   'Pearson', 'Pearson pvalue', 'Kendall', 'Kendall pvalue'])
CORRELATIONS.to_csv(os.path.join(data_dir, "018_correlations.csv"), index=False)

print("Correlations between tasks have been calculated!")

