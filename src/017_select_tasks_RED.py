import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import collections
import argparse
from default_parameters import MAX_TASKS_PER_PATHOGEN


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

print("Selecting tasks by REDUNDANCY")

dm = pd.read_csv(os.path.join(data_dir, "016_selected_tasks_MOD_DIS.csv"))
dm = dm[dm['SELECTED'] != 0]
ratios = []
priorities = []
for r in dm.iterrows():
    r = r[1]
    ratios += [r["num_pos_samples"] / r["num_samples_MOD"]]
    priorities += [int(r["task"][0])]
dm["ratio"] = ratios
dm["priority"] = priorities

# Get the selected ones
selected_tasks = set(dm[dm['SELECTED'] != 0]['task'])

positive_sets = {}
for task in os.listdir(tasks_dir):
    fname = task[:-4]
    if fname not in selected_tasks:
        continue
    df = pd.read_csv(os.path.join(tasks_dir, task))
    columns = list(df.columns)
    c = columns[-1]
    inchikeys = df[df[c] == 1]["inchikey"].tolist()
    positive_sets[fname] = set(inchikeys)


def positive_overlaps(positive_sets):
    tasks = sorted(list(positive_sets.keys()))
    R = []
    for task1 in tasks:
        for task2 in tasks:
            if task1 >= task2:
                continue
            o = len(positive_sets[task1].intersection(positive_sets[task2]))
            n1 = len(positive_sets[task1])
            n2 = len(positive_sets[task2])
            oi = o / min(n1, n2)
            t = o / (len(positive_sets[task1].union(positive_sets[task2])))
            p1 = int(task1[0])
            p2 = int(task2[0])
            dm_ = dm[dm["task"] == task1]
            for v in dm_.values:
                auroc1 = v[1]
                n_total1 = v[3]
                break
            dm_ = dm[dm["task"] == task2]
            for v in dm_.values:
                auroc2 = v[1]
                n_total2 = v[3]
                break
            r = [task1, task2, n1, n2, o, oi, t, p1, p2, auroc1, auroc2, n_total1, n_total2]
            R += [r]
    return pd.DataFrame(R, columns=["task1", "task2", "n1", "n2", "overlap", "overlap_index", "jaccard_index", "priority1", "priority2", "auroc1", "auroc2", "n_total1", "n_total2"])

dp = positive_overlaps(positive_sets)
dp = dp.sort_values("jaccard_index", ascending=False)
dp.to_csv(os.path.join(data_dir, "017_task_overlap.csv"), sep=',', index=False)

to_remove_1 = set()
for r in dp[dp["jaccard_index"] > 0.8].iterrows():
    r = r[1]
    if r["priority1"] < r["priority2"]:
        to_remove_1.add(r["task2"])
    elif r["priority1"] > r["priority2"]:
        to_remove_1.add(r["task1"])
    else:
        if r["n_total1"] > r["n_total2"]*1.25:
            to_remove_1.add(r["task2"])
        elif r["n_total2"] > r["n_total1"]*1.25:
            to_remove_1.add(r["task1"])
        else:
            if r["auroc1"] < r["auroc2"]:
                to_remove_1.add(r["task2"])
            else:
                to_remove_1.add(r["task1"])

to_remove_1 = list(to_remove_1)

valid_tasks = selected_tasks - set(to_remove_1)
valid_tasks = sorted(list(valid_tasks))

lb = np.percentile(dm["num_pos_samples"], 10)
ub = np.percentile(dm["num_pos_samples"], 90)

dm = dm[dm["task"].isin(valid_tasks)]

dm = dm[dm["num_pos_samples"] >= lb]
dm = dm[dm["num_pos_samples"] <= ub]

to_remove_2 = []

percentile_names = collections.defaultdict(list)
for fname in dm["task"].tolist():
    if "_percentile_" in fname:
        agg_name = fname.split("_percentile_")[0]
        value = int(fname.split("_percentile_")[1].split("_")[0])
        percentile_names[agg_name] += [(fname, value)]

for k,v in percentile_names.items():
    if len(v) > 1:
        for x in v:
            if x[1] == 50:
                to_remove_2 += [x[0]]

to_remove_2 = set(to_remove_2)
dm = dm[~dm["task"].isin(to_remove_2)]

to_remove_3 = []

percentage_activity_names = collections.defaultdict(list)
for fname in dm["task"].tolist():
    if "_percentile_" in fname:
        continue
    if "_percentage_activity_" in fname:
        agg_name = fname.split("_percentage_activity_")[0]
        value = int(fname.split("_percentage_activity_")[1].split("_")[0])
        percentage_activity_names[agg_name] += [(fname, value)]

for k,v in percentage_activity_names.items():
    if len(v) > 1:
        for x in v:
            if x[1] == 50:
                to_remove_3 += [x[0]]

# From lists to sets
to_remove_1 = set(to_remove_1)
to_remove_2 = set(to_remove_2)
to_remove_3 = set(to_remove_3)

# Store results
dm = pd.read_csv(os.path.join(data_dir, "016_selected_tasks_MOD_DIS.csv"))
dm = dm[dm['SELECTED'] != 0]
info = []
for r in dm['task'].tolist():
    if r in to_remove_1:
        info.append('redundancy')
    elif r in to_remove_2:
        info.append('percentile')
    elif r in to_remove_3:
        info.append('percentage')
    else:
        info.append(1)

dm['RED'] = info
dm['priority'] = [int(i[0]) for i in dm['task']]


dm = dm.sort_values(by = ["priority", "auroc_avg_MOD"], ascending=[True, False]).head(MAX_TASKS_PER_PATHOGEN)
dm.to_csv(os.path.join(data_dir, "017_selected_tasks_RED.csv"), index=False)

print(f"Number of tasks being selected by RED: {len(dm[dm['RED'] == 1])}")

dm = dm[dm["RED"] == 1]
dm.to_csv(os.path.join(data_dir, "017_selected_tasks_RED_filtered.csv"), index=False)
