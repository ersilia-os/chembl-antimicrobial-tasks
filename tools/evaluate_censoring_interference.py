"""
evaluate_censoring_interference.py - how often does step-22 qualifier handling
change a binary label near the activity cutoff?

Motivation
----------
Step 22 (`adjust_relation`, inside `dedup_most_active_by_inchikey`) rewrites
censored measurements ('>' / '<') to an exact value:
  - wrong-side censored (e.g. '>X' for a lower-is-better endpoint) has its VALUE
    replaced with the dataset extreme, so the numeric value genuinely changes;
  - informative-side censored (e.g. '<X') keeps its value X but drops the qualifier.
For a regression target this is lossy. For our BINARY use it only matters if the
rewrite lands on the wrong side of the cutoff. This tool quantifies that, per
activity cutoff level, using the SAME low/middle/high cutoffs the pipeline uses
(config/expert_cutoffs.csv). Those cutoffs vary by target_type (ORGANISM vs
SINGLE PROTEIN) and are uniform across pathogens; step 21 keeps only ORGANISM
datasets, so we always use the ORGANISM triple (e.g. MIC/IC50 = 5;10;20 uM).

Two complementary metrics, per censored molecule that SURVIVES dedup (i.e. whose
final label rests on a censored record):

  1. Ambiguity straddle - the cutoff C falls inside the true-value half-line the
     qualifier defines, so the true label is genuinely unknown and step 21 forces
     a guess:
         '>X'  ambiguous at C  <=>  X <  C     (true value in (X, inf))
         '<X'  ambiguous at C  <=>  X >  C     (true value in (-inf, X))
     (Direction-independent: it is pure geometry of the half-line vs C.)

  2. Face-value label change - the binarised label from the TRANSFORMED value
     differs from the label you would assign to the RAW stated value at face
     value (qualifier ignored). Only wrong-side records (value replaced) can
     trigger this; it is the literal "the number moved across the cutoff" count.

Worked example (lower-is-better, '>19'):
    C=10 -> not ambiguous (19>10, true value >19>10 certainly inactive), no face change.
    C=20 -> ambiguous (true value could be in (19,20] active or >20 inactive);
            face change too (raw 19<=20 active, transformed extreme > 20 inactive).

This tool ONLY reads and reports. It changes nothing and applies no decision.

Run:
    python tools/evaluate_censoring_interference.py
"""

import os
import sys
import glob
import datetime

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import CONFIGPATH, KEPT_TARGET_TYPE
from dataset_utils import get_cut_value, adjust_relation, UNIT_SLUG
from pathogen_utils import load_expert_cutoffs

OUTPUT_ROOT = os.path.join(root, "..", "output")
DOCS_DIR = os.path.join(root, "..", "docs")
os.makedirs(DOCS_DIR, exist_ok=True)

LEVELS = ["low", "middle", "high"]
KEY = ["assay_chembl_id", "activity_type"]


# --- cutoff levels (ORGANISM triple from config/expert_cutoffs.csv) -----------

def levels_from_list(cutoff_list):
    """Map an ordered expert-cutoff list to {low, middle, high} exactly as step 19."""
    n = len(cutoff_list)
    if n == 0:
        return None
    if n == 1:
        return {"low": None, "middle": cutoff_list[0], "high": None}
    if n == 2:
        return {"low": cutoff_list[0], "middle": None, "high": cutoff_list[1]}
    return {"low": cutoff_list[0], "middle": cutoff_list[1], "high": cutoff_list[-1]}


# --- provenance reconstruction ----------------------------------------------

def kept_with_provenance(g_src, direction):
    """Reproduce step-22 dedup while preserving each surviving row's ORIGINAL
    value and relation. Returns one row per InChIKey with columns:
    inchikey, value (transformed/final), orig_value, orig_relation.
    Mirrors dataset_utils.dedup_most_active_by_inchikey (verified equivalent by
    tools/validate_datasets.py check C)."""
    g = g_src[["inchikey", "value", "relation"]].copy()
    g["orig_value"] = g["value"]
    g["orig_relation"] = g["relation"].astype(str)
    if direction in (1, -1):
        cut = get_cut_value(g, direction)          # min for +1, max for -1 (original values)
        g = adjust_relation(g, direction, cut)     # rewrites value/relation for '>' / '<'
        g = g.sort_values("value", ascending=(direction == -1), kind="mergesort")
    return g.drop_duplicates(subset="inchikey", keep="first").reset_index(drop=True)


def classify(qualifier, x, v_final, cutoff, direction):
    """Return (ambiguous, face_change) for one censored survivor at one cutoff."""
    if qualifier == ">":
        ambiguous = x < cutoff
    elif qualifier == "<":
        ambiguous = x > cutoff
    else:
        ambiguous = False
    if direction == -1:                 # active iff value <= cutoff
        lbl_face = x <= cutoff
        lbl_trans = v_final <= cutoff
    else:                               # active iff value >= cutoff
        lbl_face = x >= cutoff
        lbl_trans = v_final >= cutoff
    return bool(ambiguous), bool(lbl_face != lbl_trans)


# --- source loading ---------------------------------------------------------

def load_source_numeric(pathogen):
    """Numeric activity records as step 21 sees them: value present, inchikey mapped."""
    out = os.path.join(OUTPUT_ROOT, "stage2", pathogen)
    acts = pd.read_csv(
        os.path.join(out, "08_chembl_cleaned_data.csv.gz"), low_memory=False,
        usecols=["compound_chembl_id", "assay_chembl_id", "activity_type",
                 "value", "relation"])
    counts = pd.read_csv(os.path.join(out, "07_compound_counts.csv.gz"),
                         usecols=["compound_chembl_id", "inchikey"])
    cid2ik = dict(zip(counts["compound_chembl_id"], counts["inchikey"]))
    acts["inchikey"] = acts["compound_chembl_id"].map(cid2ik)
    return acts[acts["value"].notna()].copy()


# --- per-pathogen analysis --------------------------------------------------

def analyse_pathogen(pathogen, expert_cutoffs):
    out = os.path.join(OUTPUT_ROOT, "stage4", pathogen)
    summary = pd.read_csv(os.path.join(out, "21_curation_summary.csv"), low_memory=False)
    retained = summary[(~summary["discarded"]) &
                       (summary["category"].isin(["DR", "SP"]))].reset_index(drop=True)

    # group source numeric records by pair once
    num = load_source_numeric(pathogen)
    groups = {k: v for k, v in num.groupby(KEY)}

    rows = []                 # per-dataset detail
    examples = []             # (text) worked examples: label flips across levels
    n_skipped = 0             # datasets without an expert-cutoff triple

    for r in retained.itertuples():
        key = (r.assay_id, r.activity_type)
        direction = int(r.direction)
        cutoff_list = expert_cutoffs.get(
            (r.activity_type, r.unit, KEPT_TARGET_TYPE, pathogen))
        levels = levels_from_list(cutoff_list) if cutoff_list else None
        if levels is None or key not in groups:
            n_skipped += 1
            continue

        kept = kept_with_provenance(groups[key], direction)
        n_mol = len(kept)
        cens = kept[kept["orig_relation"].isin([">", "<"])]
        n_cens = len(cens)
        # wrong-side (value actually changed): '>' for -1, '<' for +1
        wrong = ">" if direction == -1 else "<"
        n_valuechanged = int((cens["orig_relation"] == wrong).sum())

        row = dict(pathogen=pathogen, assay_id=r.assay_id,
                   activity_type=r.activity_type, unit=r.unit, category=r.category,
                   direction=direction, n_mol=n_mol, n_censored=n_cens,
                   n_value_changed=n_valuechanged)

        for level in LEVELS:
            C = levels[level]
            if C is None:
                row[f"cutoff_{level}"] = None
                row[f"ambig_{level}"] = None
                row[f"facechg_{level}"] = None
                continue
            amb = fc = 0
            for c in cens.itertuples():
                a, f = classify(c.orig_relation, c.orig_value, c.value, C, direction)
                amb += a
                fc += f
            row[f"cutoff_{level}"] = C
            row[f"ambig_{level}"] = amb
            row[f"facechg_{level}"] = fc

        rows.append(row)

        # collect a worked example: a molecule ambiguous at one level but not another
        if len(examples) < 2 and n_cens:
            for c in cens.itertuples():
                flags = {lvl: classify(c.orig_relation, c.orig_value, c.value,
                                       levels[lvl], direction)[0]
                         for lvl in LEVELS if levels[lvl] is not None}
                if len(set(flags.values())) > 1:      # differs across levels
                    cuts = {lvl: levels[lvl] for lvl in LEVELS if levels[lvl] is not None}
                    examples.append(
                        f"pair {r.assay_id} / {r.activity_type} ({r.category}, "
                        f"unit={r.unit}, direction={direction})\n"
                        f"molecule {c.inchikey}: raw '{c.orig_relation}{c.orig_value}' "
                        f"-> transformed value {c.value}\n"
                        f"cutoffs {cuts}\n"
                        f"ambiguous per level: {flags}")
                    break

    return rows, examples, n_skipped, len(retained)


# --- reporting --------------------------------------------------------------

def summarise(detail):
    """Aggregate per (pathogen, category, level)."""
    agg = []
    for level in LEVELS:
        for (pathogen, cat), sub in detail.groupby(["pathogen", "category"]):
            has = sub[f"cutoff_{level}"].notna()
            sub = sub[has]
            if not len(sub):
                continue
            n_mol = int(sub["n_mol"].sum())
            n_cens = int(sub["n_censored"].sum())
            n_vc = int(sub["n_value_changed"].sum())
            amb = int(sub[f"ambig_{level}"].sum())
            fc = int(sub[f"facechg_{level}"].sum())
            agg.append(dict(
                pathogen=pathogen, category=cat, level=level,
                n_datasets=len(sub), n_molecules=n_mol, n_censored=n_cens,
                n_value_changed=n_vc, n_ambiguous=amb, n_face_change=fc,
                ambig_pct_of_mol=round(100 * amb / n_mol, 3) if n_mol else 0.0,
                facechg_pct_of_mol=round(100 * fc / n_mol, 3) if n_mol else 0.0))
    return pd.DataFrame(agg)


def verify_current_outputs(pathogen):
    """Post-fix verification: read the ACTUAL step-22 outputs and count actives
    whose qualifier does not prove a hit (the anti-conservative residual).
    After the conservative-binarisation change this should be 0. Returns
    (n_active, n_unproven_active) or None if outputs predate the `relation` column."""
    out = os.path.join(OUTPUT_ROOT, "stage4", pathogen)
    bs = pd.read_csv(os.path.join(out, "22_binarisation_summary.csv"))
    bin_dir = os.path.join(out, "22_binarised")
    n_active = n_unproven = 0
    for r in bs.itertuples():
        slug = UNIT_SLUG.get(r.unit, str(r.unit))
        fp = os.path.join(bin_dir, r.category,
                          f"{r.assay_id}_{str(r.activity_type).replace('/', '-')}_{slug}.csv.gz")
        if not os.path.exists(fp):
            continue
        b = pd.read_csv(fp)
        if "relation" not in b.columns:
            return None  # legacy output without qualifiers; nothing to verify
        rel = b["relation"].astype(str)
        C, d = r.cutoff, int(r.direction)
        proven = (rel.isin(["=", "<"]) & (b["value"] <= C)) if d == -1 \
            else (rel.isin(["=", ">"]) & (b["value"] >= C))
        active = b["bin"] == 1
        n_active += int(active.sum())
        n_unproven += int((active & ~proven).sum())
    return n_active, n_unproven


def main():
    pathogens = sorted(
        os.path.basename(os.path.dirname(p))
        for p in glob.glob(os.path.join(OUTPUT_ROOT, "stage4", "*", "22_binarisation_summary.csv")))
    print(f"Pathogens with step-22 output: {pathogens}\n")
    expert_cutoffs = load_expert_cutoffs(CONFIGPATH)

    all_detail, all_examples = [], []
    for pathogen in pathogens:
        rows, examples, n_skipped, n_ret = analyse_pathogen(pathogen, expert_cutoffs)
        all_detail.extend(rows)
        all_examples.extend((pathogen, e) for e in examples)
        print(f"{pathogen}: {len(rows)}/{n_ret} DR/SP datasets analysed "
              f"({n_skipped} skipped - no expert cutoff for that activity_type)")

    detail = pd.DataFrame(all_detail)
    agg = summarise(detail)

    # post-fix verification: residual anti-conservative actives in actual outputs
    print("\nVerification (actual step-22 outputs, prove-the-hit rule):")
    verify_rows = []
    for pathogen in pathogens:
        res = verify_current_outputs(pathogen)
        if res is None:
            print(f"  {pathogen}: legacy output (no `relation` column) - re-run step 22 to verify")
            verify_rows.append((pathogen, None, None))
            continue
        n_active, n_unproven = res
        flag = "OK" if n_unproven == 0 else "!! RESIDUAL"
        print(f"  {pathogen}: actives={n_active}  unproven_actives={n_unproven}  [{flag}]")
        verify_rows.append((pathogen, n_active, n_unproven))

    today = datetime.date.today().isoformat()
    csv_path = os.path.join(DOCS_DIR, f"{today}_censoring-interference.csv")
    md_path = os.path.join(DOCS_DIR, f"{today}_censoring-interference.md")
    detail.to_csv(csv_path, index=False)

    # console + markdown summary
    print("\n" + "=" * 78)
    print("Interference summary (per pathogen x category x cutoff level)")
    print("=" * 78)
    print(agg.to_string(index=False))

    lines = [f"# Censoring interference near activity cutoffs (steps 22/23)", "",
             f"Generated: {today}  ", f"Tool: `tools/evaluate_censoring_interference.py`  ",
             f"Cutoffs: low/middle/high from `config/expert_cutoffs.csv` "
             f"(ORGANISM triple per activity_type/unit, as step 19; "
             f"uniform across pathogens).", "",
             "Metrics are over molecules whose FINAL binary label rests on a censored",
             "('>'/'<') record that survived step-22 dedup. **Ambiguous** = the cutoff",
             "falls inside the qualifier's true-value range (label genuinely unknown);",
             "this is an inherent property of the data, independent of labelling policy.",
             "**Face change** = a bare-threshold rewrite would give a different label",
             "than the raw stated value (the interference the OLD binarisation caused).",
             "The pipeline now labels conservatively (prove-the-hit), so ambiguous",
             "records are all inactive - see the verification section below.", "",
             "## Summary (per pathogen x category x level)", "",
             "| Pathogen | Cat | Level | Datasets | Molecules | Censored | Value-changed "
             "| Ambiguous | Ambig % | Face-change | Face % |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for a in agg.itertuples():
        lines.append(f"| {a.pathogen} | {a.category} | {a.level} | {a.n_datasets} "
                     f"| {a.n_molecules} | {a.n_censored} | {a.n_value_changed} "
                     f"| {a.n_ambiguous} | {a.ambig_pct_of_mol} | {a.n_face_change} "
                     f"| {a.facechg_pct_of_mol} |")
    lines += ["", "## Verification: residual unproven actives in actual outputs", "",
              "Actives whose qualifier does not prove a hit under the conservative rule "
              "(should be 0 post-fix; `n/a` = legacy output without `relation`).", "",
              "| Pathogen | Actives | Unproven actives |", "|---|---|---|"]
    for pathogen, na, nu in verify_rows:
        lines.append(f"| {pathogen} | {'n/a' if na is None else na} "
                     f"| {'n/a' if nu is None else nu} |")
    lines += ["", "## Worked examples (label flips across cutoff levels)", ""]
    for pathogen, e in all_examples[:12]:
        lines += [f"**{pathogen}**", "", "```", e, "```", ""]
    lines += [f"", f"Per-dataset detail: `{os.path.basename(csv_path)}`"]
    with open(md_path, "w") as fh:
        fh.write("\n".join(lines))

    if all_examples:
        print("\nWorked examples (ambiguous at some levels, not others):")
        for pathogen, e in all_examples[:6]:
            print(f"\n  [{pathogen}]\n" + "\n".join("    " + ln for ln in e.splitlines()))

    print(f"\nDetail CSV : {csv_path}")
    print(f"Report     : {md_path}")


if __name__ == "__main__":
    main()
