"""
Step 12 (v2 pipeline) - Calculate Ersilia model eos1klk for all pathogen molecules.

Builds a SINGLE cross-pathogen, deduplicated SMILES list from every pathogen's
`07_all_smiles.csv` (the deduplicated SMILES space of each pathogen, step 07) and
runs the Ersilia Model Hub model `eos1klk` (lazychemvis-reference-library) on it.
Molecules overlap across pathogens, so the union is computed once and predicted once
rather than per pathogen (556,004 unique vs 897,186 summed at time of writing).
Purely additive: nothing in steps 00-20 is modified.

Why 07_all_smiles.csv: it is the per-pathogen deduplicated SMILES space (step 07),
so the union is every molecule ever associated with any supported pathogen.

The model is fetched from DockerHub if it is not already available locally. Prediction
runs in resumable chunks of BATCH_SIZE molecules: each chunk is written to its own
`.h5` file, and chunks whose output already exists are skipped on re-run. All chunk
outputs are then merged into a single consolidated `.h5`. HDF5 is used (rather than CSV)
because the per-molecule output of a reference/embedding model is wide and CSV would be
far larger on disk.

Inputs (output/stage2/<pathogen_code>/ for every code in config/pathogens.csv):
    07_all_smiles.csv     deduplicated SMILES for that pathogen (single `smiles` column)

Outputs (output/stage2/general/):
    12_all_smiles_unique.csv           the cross-pathogen unique SMILES list (one column)
    12_eos1klk_predictions.h5          consolidated predictions for all molecules

The resumable intermediates (predictions/output_NNNNN.h5, batch_inputs/input_NNNNN.csv,
12_run_manifest.csv and the run log) are written during the run and then deleted once the
consolidated .h5 is complete. The consolidated .h5 is only written on success, so if a run
fails partway the per-chunk outputs are kept and a re-run resumes from where it stopped.

Run in the `ersilia` conda env (it has the ersilia CLI plus pandas/h5py/numpy):
    conda run -n ersilia python scripts/12_calculate_eos1klk.py
"""

import os
import sys
import glob
import subprocess

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH

# --- constants --------------------------------------------------------------
ERSILIA_MODEL = "eos1klk"          # lazychemvis-reference-library
BATCH_SIZE = 1000                  # molecules per resumable chunk

# --- paths ------------------------------------------------------------------
OUTPUT_ROOT = os.path.join(root, "..", "output")
OUTPUT = os.path.join(OUTPUT_ROOT, "stage2", "general")
PRED_DIR = os.path.join(OUTPUT, "predictions")
INPUT_DIR = os.path.join(OUTPUT, "batch_inputs")
os.makedirs(OUTPUT, exist_ok=True)
os.makedirs(PRED_DIR, exist_ok=True)
os.makedirs(INPUT_DIR, exist_ok=True)

UNIQUE_CSV = os.path.join(OUTPUT, "12_all_smiles_unique.csv")
MERGED_H5 = os.path.join(OUTPUT, f"12_{ERSILIA_MODEL}_predictions.h5")
MANIFEST = os.path.join(OUTPUT, "12_run_manifest.csv")


# --- ersilia CLI helpers ----------------------------------------------------
def ersilia(args, capture=False):
    """Invoke the ersilia CLI. Streams output unless capture=True."""
    cmd = ["ersilia"] + args
    if capture:
        return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    subprocess.run(cmd, check=True)
    return None


def model_is_local():
    """True if ERSILIA_MODEL is already fetched (appears in the local catalog)."""
    try:
        out = ersilia(["catalog", "--local"], capture=True)
    except subprocess.CalledProcessError:
        return False
    return ERSILIA_MODEL in out


# --- step 1: build the cross-pathogen unique SMILES list --------------------
def build_unique_smiles():
    pathogens = pd.read_csv(os.path.join(CONFIGPATH, "pathogens.csv"))
    seen = set()
    print("Step 12 - building cross-pathogen unique SMILES list")
    for code in pathogens["code"]:
        path = os.path.join(OUTPUT_ROOT, "stage2", code, "07_all_smiles.csv")
        if not os.path.exists(path):
            print(f"  {code}: 07_all_smiles.csv MISSING - skipped")
            continue
        df = pd.read_csv(path)
        col = "smiles" if "smiles" in df.columns else df.columns[0]
        vals = set(df[col].dropna().astype(str))
        seen |= vals
        print(f"  {code}: {len(vals)} smiles")
    # sorted for a deterministic, reproducible ordering (no stochastic step here)
    unique = sorted(seen)
    pd.DataFrame({"smiles": unique}).to_csv(UNIQUE_CSV, index=False)
    print(f"  union: {len(unique)} unique smiles -> {UNIQUE_CSV}")
    return unique


# --- step 2: predict in resumable chunks ------------------------------------
def run_chunks(smiles):
    n = len(smiles)
    n_chunks = (n + BATCH_SIZE - 1) // BATCH_SIZE
    print(f"\nPredicting {ERSILIA_MODEL} over {n} molecules in {n_chunks} "
          f"chunk(s) of {BATCH_SIZE}")

    manifest, done = [], 0
    for i in range(n_chunks):
        lo, hi = i * BATCH_SIZE, min((i + 1) * BATCH_SIZE, n)
        in_csv = os.path.join(INPUT_DIR, f"input_{i:05d}.csv")
        out_h5 = os.path.join(PRED_DIR, f"output_{i:05d}.h5")
        status = "done"
        if os.path.exists(out_h5):
            status = "skipped (exists)"
        else:
            pd.DataFrame({"smiles": smiles[lo:hi]}).to_csv(in_csv, index=False)
            # internal batch_size == chunk size: one chunk is one ersilia run
            ersilia(["run", "-i", in_csv, "-o", out_h5, "-b", str(BATCH_SIZE)])
            done += 1
        manifest.append({"chunk": i, "row_start": lo, "row_end": hi,
                         "n": hi - lo, "output": os.path.basename(out_h5),
                         "status": status})
        if (i + 1) % 25 == 0 or i + 1 == n_chunks:
            print(f"  {i + 1}/{n_chunks} chunks ({done} run this session)")

    pd.DataFrame(manifest).to_csv(MANIFEST, index=False)
    print(f"  manifest -> {MANIFEST}")
    return [os.path.join(PRED_DIR, m["output"]) for m in manifest]


# --- step 3: merge per-chunk h5 into one consolidated h5 --------------------
def merge_h5(batch_files):
    """Concatenate row-aligned datasets (Keys/Inputs/Values) across chunk files;
    copy constant datasets (e.g. Features/column names) once from the first file."""
    import h5py

    batch_files = [f for f in batch_files if os.path.exists(f)]
    if not batch_files:
        print("  no chunk outputs to merge")
        return

    # discover which datasets are row-aligned by comparing first-axis length to
    # this file's row count (taken from the first dataset that looks like a vector)
    def nrows(h5):
        # row count = max first-axis length present in the file
        return max(int(h5[k].shape[0]) for k in h5.keys())

    with h5py.File(batch_files[0], "r") as f0:
        names = list(f0.keys())
        r0 = nrows(f0)
        row_aligned = {k: (f0[k].shape[0] == r0) for k in names}
        constants = {k: f0[k][()] for k, aligned in row_aligned.items() if not aligned}

    collected = {k: [] for k, a in row_aligned.items() if a}
    total = 0
    for path in batch_files:
        with h5py.File(path, "r") as f:
            for k in collected:
                collected[k].append(f[k][()])
        total += len(collected[next(iter(collected))][-1])

    with h5py.File(MERGED_H5, "w") as out:
        for k, parts in collected.items():
            out.create_dataset(k, data=np.concatenate(parts, axis=0))
        for k, v in constants.items():
            out.create_dataset(k, data=v)
    print(f"  merged {len(batch_files)} chunk(s), {total} rows -> {MERGED_H5}")


# --- cleanup: drop resumable intermediates once the merged .h5 is complete --
def cleanup_intermediates():
    """Remove the resumable intermediates now that the consolidated .h5 is written:
    the per-chunk `predictions/` and `batch_inputs/` directories, the run manifest,
    and the run log. The unique SMILES list and the merged predictions .h5 are kept."""
    import shutil

    removed = []
    for d in (PRED_DIR, INPUT_DIR):
        if os.path.isdir(d):
            shutil.rmtree(d)
            removed.append(os.path.basename(d) + "/")
    for f in [MANIFEST] + glob.glob(os.path.join(OUTPUT, "*run*.log")):
        if os.path.exists(f):
            os.remove(f)
            removed.append(os.path.basename(f))
    print(f"  cleaned intermediates: {', '.join(removed) if removed else '(none)'}")


# --- main -------------------------------------------------------------------
def main():
    smiles = build_unique_smiles()

    if not model_is_local():
        print(f"\nFetching {ERSILIA_MODEL} from DockerHub...")
        ersilia(["fetch", ERSILIA_MODEL, "--from_dockerhub"])
    else:
        print(f"\n{ERSILIA_MODEL} already available locally")

    print(f"Serving {ERSILIA_MODEL}...")
    ersilia(["serve", ERSILIA_MODEL])
    try:
        batch_files = run_chunks(smiles)
        print("\nMerging chunk outputs...")
        merge_h5(batch_files)
    finally:
        print("Closing model server...")
        ersilia(["close"])

    # Clean up intermediates only after a successful run (the consolidated .h5 is
    # written); if the run failed earlier, the per-chunk outputs stay for a resume.
    if os.path.exists(MERGED_H5):
        print("\nCleaning up intermediate batches...")
        cleanup_intermediates()
    print("Done.")


if __name__ == "__main__":
    main()
