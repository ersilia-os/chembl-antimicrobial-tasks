"""Manually append step-09 rows for potency assays recovered by the 2026-07-15 config edits
(XC50/EC90/MBC/MBC90/MFC given direction -1 + DR category), WITHOUT re-running the GPU/Ollama
step 09.

Rationale: step 09 only adds LLM-*curated* target fields on top of the raw ChEMBL assay metadata
already present in 07_assays_raw.csv. For the recovered assays ~98.5% are raw target_type=ORGANISM,
for which extra_curation_target_type() returns ORGANISM deterministically regardless of the LLM value.
So we synthesise the step-09 rows with the approved rule:
    curated field = raw field   (organism / strain / target_type / target_chembl_id;
                                 target_name_curated = target_organism)
    LLM-only fields (atcc_id, mutations, known_drug_resistances, culture_media) = blank
Only assays NOT already present in 09 are appended (assays that also had a surviving type — MIC etc. —
are already enriched; step 21 enriches per assay_id). The 18 UNCHECKED/NON-MOLECULAR assays fall out
correctly (11 are target-directed histidine-kinase/AcrB; 7 are genuine but n=3 < MIN_MOLECULES).

Idempotent: re-running skips assays now present in 09. Original 09 files are backed up once to
09_assays_parameters_full.csv.bak_preMBC. A manifest of appended (pathogen, assay_id) rows is written
to 09_manual_appended_assays.csv per pathogen for provenance (these rows are overwritten if step 09
is ever genuinely re-run).
"""
import os, sys, shutil
import pandas as pd
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))
S2 = os.path.join(root, "..", "output", "stage2")

RECOVERED_TYPES = {"XC50", "EC90", "MBC", "MBC90", "MFC"}
RECOVERED_UNIT = "umol.L-1"

# curated <- raw mapping (approved rule)
CURATED_FROM_RAW = {
    "assay_organism_curated": "assay_organism",
    "assay_strain_curated": "assay_strain",
    "target_type_curated": "target_type",
    "target_chembl_id_curated": "target_chembl_id",
    "target_name_curated": "target_organism",
}
BLANK_LLM_COLS = ["atcc_id", "mutations", "known_drug_resistances", "culture_media"]


def process(pathogen):
    d = os.path.join(S2, pathogen)
    f09 = os.path.join(d, "09_assays_parameters_full.csv")
    f07 = os.path.join(d, "07_assays_raw.csv")
    if not (os.path.exists(f09) and os.path.exists(f07)):
        return None

    p09 = pd.read_csv(f09, low_memory=False)
    cols09 = list(p09.columns)
    have = set(p09["assay_id"].unique())

    ar = pd.read_csv(f07, low_memory=False)
    rec = ar[(ar["activity_type"].astype(str).str.upper().isin(RECOVERED_TYPES)) &
             (ar["unit"] == RECOVERED_UNIT)].copy()
    new = rec[~rec["assay_id"].isin(have)].copy()
    if not len(new):
        return dict(pathogen=pathogen, added_rows=0, added_assays=0)

    # build 09-format rows: raw cols carry over (same names present in 07_assays_raw);
    # curated cols mirror raw; LLM-only cols blank; direction from config (-1 for these potencies).
    for cur, raw in CURATED_FROM_RAW.items():
        new[cur] = new[raw]
    for c in BLANK_LLM_COLS:
        new[c] = np.nan
    new["direction"] = -1  # step 21 recomputes live from config anyway

    # any 09 column not produced above (should be none) -> blank, then order exactly like 09
    for c in cols09:
        if c not in new.columns:
            new[c] = np.nan
    new = new[cols09]

    # backup original once, then append
    bak = f09 + ".bak_preMBC"
    if not os.path.exists(bak):
        shutil.copy2(f09, bak)
    out = pd.concat([p09, new], ignore_index=True)
    out.to_csv(f09, index=False)

    # provenance manifest
    man = new[["assay_id", "activity_type", "unit", "target_type", "assay_organism"]].drop_duplicates()
    man.to_csv(os.path.join(d, "09_manual_appended_assays.csv"), index=False)

    return dict(pathogen=pathogen, added_rows=len(new),
                added_assays=new["assay_id"].nunique(),
                new_09_total=len(out))


if __name__ == "__main__":
    codes = sys.argv[1:] or sorted(
        p for p in os.listdir(S2) if os.path.isdir(os.path.join(S2, p)) and p != "general")
    print(f"{'pathogen':14} {'added_rows':>10} {'added_assays':>12} {'09_total_after':>15}")
    tot_r = tot_a = 0
    for p in codes:
        r = process(p)
        if r is None:
            continue
        tot_r += r["added_rows"]; tot_a += r["added_assays"]
        if r["added_rows"]:
            print(f"{p:14} {r['added_rows']:10d} {r['added_assays']:12d} {r.get('new_09_total',''):>15}")
    print(f"\nTOTAL appended: {tot_r} rows / {tot_a} assays")
