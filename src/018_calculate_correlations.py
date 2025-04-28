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