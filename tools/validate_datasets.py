"""
validate_datasets.py - integrity & correctness validation for pipeline steps 21 & 22.

QC tool (NOT a pipeline step). Confirms that:
  - step 21 (21_build_curated_datasets.py) loses no rows, does not
    cross-contaminate assays, deduplicates by InChIKey correctly (keeps the
    most-active record), and keeps its summary table in sync with the datasets
    written to disk; and
  - step 22 (22_binarise_datasets.py) is purely additive (drops nothing),
    preserves molecule<->value identity, and binarises with the exact
    direction/cutoff rule it claims to.

Design: every check runs to completion and COLLECTS violations (nothing crashes
on the first failure). The console prints per-pathogen worked examples so the
behaviour is eyeball-verifiable, and a PASS/FAIL markdown report is written to
docs/. Exit code is non-zero if any check FAILs, so this can later gate CI.

This tool ONLY reads and validates. It never edits data/output and never applies
any threshold or drop decision (per CLAUDE.md human-sign-off rules). A failing
check is surfaced, not silently "fixed".

Run:
    python tools/validate_datasets.py
"""

import os
import sys
import glob
import datetime

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (
    CONFIGPATH, MIN_MOLECULES, COADD_SOURCE_LABEL, KEPT_TARGET_TYPE, RANDOM_SEED,
)
from dataset_utils import (
    dedup_most_active_by_inchikey, binarize_conservative,
    load_category_cutoffs, UNIT_SLUG,
)
from pathogen_utils import load_pubchem_assays

# --- config -----------------------------------------------------------------
OUTPUT_ROOT = os.path.join(root, "..", "output")

# Validate every pathogen that has step-21 output; pass pathogen codes as CLI
# args to restrict (e.g. `python tools/validate_datasets.py mtuberculosis hpylori`).
def discover_pathogens():
    return sorted(os.path.basename(os.path.dirname(p))
                  for p in glob.glob(os.path.join(OUTPUT_ROOT, "stage4", "*", "21_curation_summary.csv")))

PATHOGENS = sys.argv[1:] if len(sys.argv) > 1 else discover_pathogens()

# Reconstruction sample size for the (expensive) dedup re-run check (C).
DEDUP_SAMPLE = 60

DOCS_DIR = os.path.join(root, "..", "docs")
os.makedirs(DOCS_DIR, exist_ok=True)

RNG = np.random.default_rng(RANDOM_SEED)


# --- result accumulator -----------------------------------------------------

class Report:
    """Collects per-(pathogen, check) results and renders console + markdown."""

    def __init__(self):
        self.rows = []          # list of dicts: pathogen, check, title, passed, n_violations, note
        self.examples = []      # list of (pathogen, title, text) worked examples

    def add(self, pathogen, check, title, passed, n_violations=0, note=""):
        self.rows.append(dict(pathogen=pathogen, check=check, title=title,
                              passed=bool(passed), n_violations=int(n_violations),
                              note=str(note)))
        flag = "PASS" if passed else "FAIL"
        extra = "" if passed else f"  ({n_violations} violations)"
        print(f"    [{flag}] {check} {title}{extra}"
              + (f"  - {note}" if note else ""))

    def add_example(self, pathogen, title, text):
        self.examples.append((pathogen, title, text))
        print(f"\n    --- worked example: {title} ---\n"
              + "\n".join("      " + ln for ln in text.splitlines()) + "\n")

    def all_passed(self):
        return all(r["passed"] for r in self.rows)

    def write_markdown(self, path):
        today = datetime.date.today().isoformat()
        lines = [f"# Validation report - pipeline steps 21 & 22", "",
                 f"Generated: {today}  ",
                 f"Tool: `tools/validate_datasets.py`  ",
                 f"Pathogens: {', '.join(sorted({r['pathogen'] for r in self.rows}))}", ""]

        overall = "ALL CHECKS PASS" if self.all_passed() else "SOME CHECKS FAILED"
        lines += [f"**Overall: {overall}**", ""]

        # summary table per pathogen
        for pathogen in sorted({r["pathogen"] for r in self.rows}):
            prows = [r for r in self.rows if r["pathogen"] == pathogen]
            lines += [f"## {pathogen}", "",
                      "| Check | Title | Result | Violations | Note |",
                      "|---|---|---|---|---|"]
            for r in prows:
                res = "PASS" if r["passed"] else "**FAIL**"
                lines.append(f"| {r['check']} | {r['title']} | {res} "
                             f"| {r['n_violations']} | {r['note']} |")
            lines.append("")

        # worked examples
        if self.examples:
            lines += ["## Worked examples", ""]
            for pathogen, title, text in self.examples:
                lines += [f"### {pathogen} - {title}", "", "```", text, "```", ""]

        with open(path, "w") as fh:
            fh.write("\n".join(lines))


# --- helpers ----------------------------------------------------------------

def dataset_fname(assay_id, activity_type):
    """Reproduce the step-22/23 filename sanitisation of activity_type."""
    return f"{assay_id}_{str(activity_type).replace('/', '-')}"


def load_source_records(pathogen):
    """Rebuild the source activity records as step 21 sees them: map
    compound_chembl_id -> inchikey via 07_compound_counts. Returns ALL records
    (with a `text_flag` column); the numeric subset is `value.notna()`. Callers
    use the full set for contamination (recovered qualitative inactives are
    non-numeric) and the numeric subset for dedup checks.
    """
    out = os.path.join(OUTPUT_ROOT, "stage2", pathogen)
    acts = pd.read_csv(
        os.path.join(out, "08_chembl_cleaned_data.csv.gz"), low_memory=False,
        usecols=["compound_chembl_id", "assay_chembl_id", "activity_type",
                 "value", "relation", "text_flag"])
    counts = pd.read_csv(os.path.join(out, "07_compound_counts.csv.gz"),
                         usecols=["compound_chembl_id", "inchikey"])
    cid2ik = dict(zip(counts["compound_chembl_id"], counts["inchikey"]))
    acts["inchikey"] = acts["compound_chembl_id"].map(cid2ik)
    return acts


# --- STEP 21 CHECKS ---------------------------------------------------------

def check_step21(pathogen, rep):
    out = os.path.join(OUTPUT_ROOT, "stage4", pathogen)
    summary = pd.read_csv(os.path.join(out, "21_curation_summary.csv"), low_memory=False)
    datasets_dir = os.path.join(out, "21_datasets")

    retained = summary[~summary["discarded"]].reset_index(drop=True)
    discarded = summary[summary["discarded"]]

    # A1 accounting -----------------------------------------------------------
    xor_ok = (summary["discarded"] == (summary["discard_reason"].fillna("") != "")).all()
    reason_sum = int((summary["discard_reason"].fillna("") != "").sum())
    a1_ok = xor_ok and (reason_sum == len(discarded)) and \
        (len(summary) == len(retained) + len(discarded))
    rep.add(pathogen, "A1", "accounting (retained XOR reason; counts sum)",
            a1_ok, 0 if a1_ok else 1,
            f"{len(summary)} pairs = {len(retained)} retained + {len(discarded)} discarded")

    # A2 coverage (summary <-> files) ----------------------------------------
    expected = {dataset_fname(r.assay_id, r.activity_type) for r in retained.itertuples()}
    actual = {os.path.basename(p)[:-len(".csv.gz")]
              for p in glob.glob(os.path.join(datasets_dir, "*.csv.gz"))}
    missing = expected - actual          # retained but no file
    orphan = actual - expected           # file but not retained
    a2_ok = not missing and not orphan
    rep.add(pathogen, "A2", "coverage summary<->files", a2_ok,
            len(missing) + len(orphan),
            f"missing files: {len(missing)}, orphan files: {len(orphan)}")

    # A3 attrition monotonic --------------------------------------------------
    stats = pd.read_csv(os.path.join(out, "21_curation_stats.csv"))
    a3_viol = 0
    for col in ["datasets", "chemical_space"]:
        seq = stats[col].tolist()
        a3_viol += sum(1 for i in range(1, len(seq)) if seq[i] > seq[i - 1])
    rep.add(pathogen, "A3", "attrition monotonic non-increasing",
            a3_viol == 0, a3_viol)

    # Row-level checks over each retained dataset: B1, B2, B3, B4, B5, B6 -----
    # Source records: ALL records (numeric + qualitative) for contamination, and
    # the numeric subset for dedup. Recovered qualitative-inactive molecules are
    # non-numeric, so contamination must be checked against the full record set.
    recs = load_source_records(pathogen)
    num = recs[recs["value"].notna()].copy()
    src_sets = (recs.groupby(["assay_chembl_id", "activity_type"])["inchikey"]
                .agg(lambda s: set(s.dropna())).to_dict())
    has_dataset_col = "n_mol_dataset" in retained.columns

    b1 = b2 = b3 = b5 = b6 = 0
    b1_ex, b2_ex, b3_ex, b5_ex, b6_ex = [], [], [], [], []
    ret_lookup = {(r.assay_id, r.activity_type): r for r in retained.itertuples()}

    for (assay_id, activity_type), r in ret_lookup.items():
        fp = os.path.join(datasets_dir, dataset_fname(assay_id, activity_type) + ".csv.gz")
        if not os.path.exists(fp):
            continue  # already flagged by A2
        df = pd.read_csv(fp)

        # B1 InChIKey uniqueness
        if df["inchikey"].duplicated().any():
            b1 += 1
            if len(b1_ex) < 5:
                b1_ex.append(f"{assay_id}/{activity_type}: "
                             f"{int(df['inchikey'].duplicated().sum())} dup inchikeys")
        # B2 row count == n_mol_dataset (numeric molecules + recovered inactives)
        expected_n = int(r.n_mol_dataset) if has_dataset_col else int(r.n_mol_numeric)
        if len(df) != expected_n:
            b2 += 1
            if len(b2_ex) < 5:
                b2_ex.append(f"{assay_id}/{activity_type}: file={len(df)} "
                             f"vs n_mol_dataset={expected_n}")
        # B3 subset of source (no cross-assay contamination) - full record set
        src = src_sets.get((assay_id, activity_type), set())
        extra = set(df["inchikey"].dropna()) - src
        if extra:
            b3 += 1
            if len(b3_ex) < 5:
                b3_ex.append(f"{assay_id}/{activity_type}: {len(extra)} inchikeys not in source")
        # B5 no nulls in identity columns. value/relation MAY be NaN for recovered
        # qualitative-inactive molecules (labelled inactive downstream); inchikey
        # and smiles must always be present. `direction` is only required for
        # DR/SP pairs (the ones binarised); UNCATEGORIZED pairs legitimately lack it.
        nullcols = [c for c in ["inchikey", "smiles"] if df[c].isna().any()]
        if r.category in ("DR", "SP") and df["direction"].isna().any():
            nullcols.append("direction")
        if nullcols:
            b5 += 1
            if len(b5_ex) < 5:
                b5_ex.append(f"{assay_id}/{activity_type}: nulls in {nullcols}")
        # B6 qualifiers valid: relation must be '=', '>', '<' or NaN (recovered
        # inactive). Any other operator is invalid.
        if "relation" in df.columns:
            bad = df.loc[df["relation"].notna() &
                         ~df["relation"].astype(str).isin(["=", ">", "<"]), "relation"]
            if len(bad):
                b6 += 1
                if len(b6_ex) < 5:
                    b6_ex.append(f"{assay_id}/{activity_type}: "
                                 f"{bad.value_counts(dropna=False).to_dict()}")

    rep.add(pathogen, "B1", "InChIKey unique per dataset", b1 == 0, b1,
            "; ".join(b1_ex))
    rep.add(pathogen, "B2", "row count == n_mol_dataset", b2 == 0, b2,
            "; ".join(b2_ex))
    rep.add(pathogen, "B3", "no cross-assay contamination", b3 == 0, b3,
            "; ".join(b3_ex))
    rep.add(pathogen, "B5", "no nulls in inchikey/smiles (+direction for DR/SP)",
            b5 == 0, b5, "; ".join(b5_ex))
    rep.add(pathogen, "B6", "qualifiers valid ({=,>,<} or NaN)",
            b6 == 0, b6, "; ".join(b6_ex))

    # B4 scalar unit/direction on retained pairs ------------------------------
    # only DR/SP pairs are binarised and thus require a scalar unit & direction;
    # UNCATEGORIZED retained pairs (no cutoff) may legitimately have NaN unit/direction.
    drsp_ret = retained[retained["category"].isin(["DR", "SP"])]
    b4_bad = drsp_ret[drsp_ret["unit"].isna() | drsp_ret["direction"].isna()]
    rep.add(pathogen, "B4", "retained DR/SP pairs have scalar unit & direction",
            len(b4_bad) == 0, len(b4_bad))

    # D filter correctness on retained pairs ----------------------------------
    pubchem_ids = load_pubchem_assays(pathogen)
    d_viol = 0
    d_ex = []
    for r in retained.itertuples():
        reasons = []
        if r.target_type_curated_extra != KEPT_TARGET_TYPE:
            reasons.append(f"not {KEPT_TARGET_TYPE}")
        if not (r.n_numeric > 0):
            reasons.append("no numeric")
        if r.n_numeric_units > 1:
            reasons.append("multi-unit")
        if r.source_label == COADD_SOURCE_LABEL:
            reasons.append("CoADD")
        if r.assay_id in pubchem_ids:
            reasons.append("PubChem")
        size = int(r.n_mol_dataset) if has_dataset_col else int(r.n_mol_numeric)
        if not (size > MIN_MOLECULES):
            reasons.append(f"<= {MIN_MOLECULES} mols")
        if reasons:
            d_viol += 1
            if len(d_ex) < 5:
                d_ex.append(f"{r.assay_id}/{r.activity_type}: {reasons}")
    rep.add(pathogen, "D", "no retained pair violates a curation filter",
            d_viol == 0, d_viol, "; ".join(d_ex))

    # C dedup most-active: reconstruct on a sample and compare ----------------
    ret_keys = list(ret_lookup.keys())
    if ret_keys:
        idx = RNG.choice(len(ret_keys), size=min(DEDUP_SAMPLE, len(ret_keys)),
                         replace=False)
        c_viol = 0
        c_ex = []
        for i in idx:
            assay_id, activity_type = ret_keys[i]
            r = ret_lookup[(assay_id, activity_type)]
            if pd.isna(r.direction):
                continue  # UNCATEGORIZED pair (no direction); not dedup-sorted / not binarised
            fp = os.path.join(datasets_dir, dataset_fname(assay_id, activity_type) + ".csv.gz")
            if not os.path.exists(fp):
                continue
            got = pd.read_csv(fp)
            # compare only the numeric molecules (exclude recovered inactives,
            # which carry NaN value/relation and are not part of the dedup).
            got = got[got["value"].notna()]
            g = num[(num["assay_chembl_id"] == assay_id) &
                    (num["activity_type"] == activity_type)]
            recon = dedup_most_active_by_inchikey(
                g[["inchikey", "value", "relation"]].copy(), int(r.direction),
                preserve_qualifier=True)
            # compare value per inchikey
            a = got.set_index("inchikey")["value"].sort_index()
            b = recon.set_index("inchikey")["value"].sort_index()
            if not (a.index.equals(b.index) and np.allclose(a.values, b.values)):
                c_viol += 1
                if len(c_ex) < 5:
                    c_ex.append(f"{assay_id}/{activity_type}: dedup mismatch")
        rep.add(pathogen, "C", f"dedup keeps most-active (sample of {len(idx)})",
                c_viol == 0, c_viol, "; ".join(c_ex))

        # worked example: a molecule that had >1 raw record in some retained pair
        _dedup_worked_example(pathogen, num, ret_lookup, datasets_dir, rep)
        # worked example: a censored (>/<) source record resolved by adjust_relation
        _qualifier_worked_example(pathogen, num, ret_lookup, datasets_dir, rep)

    return summary, retained


def _qualifier_worked_example(pathogen, num, ret_lookup, datasets_dir, rep):
    """Find a retained pair with a wrong-side censored singleton and show how
    adjust_relation resolved it: raw '>'/'<' value -> group extreme, relation '='."""
    for (assay_id, activity_type), r in ret_lookup.items():
        if pd.isna(r.direction):
            continue
        direction = int(r.direction)
        wrong = ">" if direction == -1 else "<"   # wrong-side qualifier for this direction
        g = num[(num["assay_chembl_id"] == assay_id) &
                (num["activity_type"] == activity_type)]
        vc = g["inchikey"].value_counts()
        singles = set(vc[vc == 1].index)
        cand = g[(g["relation"].astype(str) == wrong) & (g["inchikey"].isin(singles))]
        if len(cand) == 0:
            continue
        row = cand.iloc[0]
        ik = row["inchikey"]
        chosen = pd.read_csv(os.path.join(
            datasets_dir, dataset_fname(assay_id, activity_type) + ".csv.gz"))
        out = chosen[chosen["inchikey"] == ik][["value", "relation"]].iloc[0]
        text = (f"pair: {assay_id} / {activity_type}   direction={direction}\n"
                f"molecule: {ik}\n\n"
                f"raw source record   : value={row['value']}  relation='{row['relation']}'  "
                f"(wrong-side censored non-inhibitor)\n"
                f"preserved in output : value={out['value']}  relation='{out['relation']}'  "
                f"(qualifier kept; conservative binarisation labels it inactive "
                f"unless it proves active)")
        rep.add_example(pathogen, "step 21 preserves censoring qualifier", text)
        return


def _dedup_worked_example(pathogen, num, ret_lookup, datasets_dir, rep):
    """Find one retained pair where some InChIKey has multiple raw records and
    show raw values/relations -> chosen value, so 'keep most active' is visible."""
    for (assay_id, activity_type), r in ret_lookup.items():
        if pd.isna(r.direction):
            continue
        g = num[(num["assay_chembl_id"] == assay_id) &
                (num["activity_type"] == activity_type)]
        vc = g["inchikey"].value_counts()
        multi = vc[vc > 1]
        if len(multi) == 0:
            continue
        ik = multi.index[0]
        raw = g[g["inchikey"] == ik][["value", "relation"]].reset_index(drop=True)
        fp = os.path.join(datasets_dir, dataset_fname(assay_id, activity_type) + ".csv.gz")
        chosen = pd.read_csv(fp)
        chosen = chosen[chosen["inchikey"] == ik]["value"]
        direction = int(r.direction)
        rule = ("lower = more active (keep min)" if direction == -1
                else "higher = more active (keep max)")
        text = (f"pair: {assay_id} / {activity_type}   direction={direction} ({rule})\n"
                f"molecule: {ik}\n\n"
                f"raw records:\n{raw.to_string(index=False)}\n\n"
                f"value kept in 21_datasets: {chosen.iloc[0]}")
        rep.add_example(pathogen, "step 21 dedup keeps most-active record", text)
        return


# --- STEP 22 CHECKS ---------------------------------------------------------

def check_step22(pathogen, summary22, retained22, rep):
    out = os.path.join(OUTPUT_ROOT, "stage4", pathogen)
    bin_summary = pd.read_csv(os.path.join(out, "22_binarisation_summary.csv"))
    cutoffs = load_category_cutoffs(CONFIGPATH)
    datasets22 = os.path.join(out, "21_datasets")
    bin_dir = os.path.join(out, "22_binarised")

    # E3 coverage: retained DR/SP pairs <-> binarised files -------------------
    ret_drsp = retained22[retained22["category"].isin(["DR", "SP"])]
    expected = set()
    for r in ret_drsp.itertuples():
        slug = UNIT_SLUG.get(r.unit, str(r.unit))
        expected.add((r.category, dataset_fname(r.assay_id, r.activity_type) + f"_{slug}"))
    actual = set()
    for cat in ["DR", "SP"]:
        for p in glob.glob(os.path.join(bin_dir, cat, "*.csv.gz")):
            actual.add((cat, os.path.basename(p)[:-len(".csv.gz")]))
    missing = expected - actual
    orphan = actual - expected
    rep.add(pathogen, "E3", "coverage retained DR/SP <-> binarised files",
            not missing and not orphan, len(missing) + len(orphan),
            f"missing: {len(missing)}, orphan: {len(orphan)}")

    # Per-binarised-file checks: E1, E2, F1, F4 ------------------------------
    e1 = e2 = f1 = f4 = f5 = n_legacy = 0
    e1_ex, e2_ex, f1_ex, f5_ex = [], [], [], []
    dr_example = sp_example = None

    for r in bin_summary.itertuples():
        slug = UNIT_SLUG.get(r.unit, str(r.unit))
        fp_bin = os.path.join(bin_dir, r.category,
                              dataset_fname(r.assay_id, r.activity_type) + f"_{slug}.csv.gz")
        fp_src = os.path.join(datasets22,
                              dataset_fname(r.assay_id, r.activity_type) + ".csv.gz")
        if not (os.path.exists(fp_bin) and os.path.exists(fp_src)):
            continue
        b = pd.read_csv(fp_bin)
        s = pd.read_csv(fp_src)

        # E1 row conservation
        if len(b) != len(s):
            e1 += 1
            if len(e1_ex) < 5:
                e1_ex.append(f"{r.assay_id}/{r.activity_type}: bin={len(b)} vs src={len(s)}")

        # E2 molecule<->value identity preserved (NaN values for recovered
        # inactives must match on both sides -> equal_nan=True)
        mb = b.set_index("inchikey")["value"].sort_index()
        ms = s.set_index("inchikey")["value"].sort_index()
        if not (mb.index.equals(ms.index)
                and np.allclose(mb.values, ms.values, equal_nan=True)):
            e2 += 1
            if len(e2_ex) < 5:
                e2_ex.append(f"{r.assay_id}/{r.activity_type}: inchikey/value mismatch")

        # qualifier-dependent checks (F1/F5/G) need the `relation` column; legacy
        # outputs (pre conservative-binarisation) lack it -> skip gracefully.
        if "relation" not in b.columns or "relation" not in s.columns:
            n_legacy += 1
            continue

        # F1 recompute bin independently (conservative rule) and compare
        cutoff = cutoffs.get(r.unit)
        direction = int(s["direction"].iloc[0])
        recomputed = binarize_conservative(
            s[["value", "relation"]].copy(), cutoff, direction)["bin"].values
        got = b.set_index("inchikey").reindex(s["inchikey"])["bin"].values \
            if b["inchikey"].is_unique and s["inchikey"].is_unique else None
        if got is None or not np.array_equal(got, recomputed):
            f1 += 1
            if len(f1_ex) < 5:
                f1_ex.append(f"{r.assay_id}/{r.activity_type}: bin recompute mismatch")

        # F4 label domain
        if not b["bin"].isin([0, 1]).all():
            f4 += 1

        # F5 no unproven actives: every bin==1 must be PROVEN active by its
        # qualifier (independent of the binarisation function).
        rel = b["relation"].astype(str)
        if direction == -1:
            proven = rel.isin(["=", "<"]) & (b["value"] <= cutoff)
        else:
            proven = rel.isin(["=", ">"]) & (b["value"] >= cutoff)
        unproven_active = int(((b["bin"] == 1) & ~proven).sum())
        if unproven_active:
            f5 += 1
            if len(f5_ex) < 5:
                f5_ex.append(f"{r.assay_id}/{r.activity_type}: {unproven_active} unproven actives")

        # collect one DR and one SP non-degenerate example for direction sanity
        if dr_example is None and r.category == "DR" and not r.degenerate:
            dr_example = (r, b, direction, cutoff)
        if sp_example is None and r.category == "SP" and not r.degenerate:
            sp_example = (r, b, direction, cutoff)

    rep.add(pathogen, "E1", "binarisation drops no rows", e1 == 0, e1, "; ".join(e1_ex))
    rep.add(pathogen, "E2", "molecule<->value identity preserved", e2 == 0, e2,
            "; ".join(e2_ex))
    rep.add(pathogen, "F1", "bin recomputed (conservative) matches stored",
            f1 == 0, f1, "; ".join(f1_ex))
    rep.add(pathogen, "F5", "no unproven actives (prove-the-hit)", f5 == 0, f5,
            "; ".join(f5_ex))
    if n_legacy:
        rep.add(pathogen, "F0", "step-23 outputs carry `relation` (conservative)",
                False, n_legacy,
                f"{n_legacy} legacy datasets without `relation` - re-run step 22")
    rep.add(pathogen, "F4", "bin in {0,1}", f4 == 0, f4)

    # F2 cutoff provenance ----------------------------------------------------
    f2_viol = 0
    for r in bin_summary.itertuples():
        if not np.isclose(r.cutoff, cutoffs.get(r.unit, np.nan)):
            f2_viol += 1
    rep.add(pathogen, "F2", "summary cutoff == config cutoff", f2_viol == 0, f2_viol)

    # F3 summary self-consistency --------------------------------------------
    f3_viol = 0
    f3_ex = []
    for r in bin_summary.itertuples():
        ok = (r.n_active + r.n_inactive == r.n_molecules)
        ratio_ok = (r.n_molecules == 0) or np.isclose(
            r.active_ratio, round(r.n_active / r.n_molecules, 5))
        degen_ok = (bool(r.degenerate) == (r.n_active == 0 or r.n_inactive == 0))
        if not (ok and ratio_ok and degen_ok):
            f3_viol += 1
            if len(f3_ex) < 5:
                f3_ex.append(f"{r.assay_id}/{r.activity_type}")
    rep.add(pathogen, "F3", "summary counts/ratio/degenerate consistent",
            f3_viol == 0, f3_viol, "; ".join(f3_ex))

    # G direction sanity worked examples (qualifier-aware / conservative) -----
    for label, ex in [("DR", dr_example), ("SP", sp_example)]:
        if ex is None:
            continue
        r, b, direction, cutoff = ex
        rel = b["relation"].astype(str)
        if direction == -1:
            rule = f"active iff relation in {{=,<}} and value <= {cutoff}"
            proven = rel.isin(["=", "<"]) & (b["value"] <= cutoff)
        else:
            rule = f"active iff relation in {{=,>}} and value >= {cutoff}"
            proven = rel.isin(["=", ">"]) & (b["value"] >= cutoff)
        # every active must be proven; conservative inactives may sit on the
        # value side of the cutoff when their qualifier does not prove a hit.
        ok = bool(((b["bin"] == 1) == proven).all())
        preview = b.head(8)[["inchikey", "value", "relation", "bin"]].copy()

        def _reason(v, rl, bn):
            side = "<=" if direction == -1 else ">="
            if bn == 1:
                return f"'{rl}{v}' proves value {side} {cutoff} -> active"
            beyond = (direction == -1 and v > cutoff) or (direction == +1 and v < cutoff)
            if beyond:
                return f"value {v} beyond cutoff -> inactive"
            return f"'{rl}{v}' does not prove a hit -> inactive (conservative)"

        preview["reason"] = [_reason(v, rl, bn) for v, rl, bn in
                             zip(preview["value"], preview["relation"], preview["bin"])]
        n_conservative = int(((b["bin"] == 0) & ~(
            (b["value"] > cutoff) if direction == -1 else (b["value"] < cutoff))).sum())
        text = (f"pair: {r.assay_id} / {r.activity_type}   category={label}   "
                f"unit={r.unit}   direction={direction}   cutoff={cutoff}\n"
                f"rule: {rule}\n"
                f"inactives whose value is on the active side but qualifier "
                f"does not prove a hit (kept inactive, conservative): {n_conservative}\n\n"
                f"{preview.to_string(index=False)}")
        rep.add_example(pathogen, f"step 22 conservative labelling ({label})", text)
        rep.add(pathogen, "G", f"conservative labelling ({label}): actives all proven",
                ok, 0 if ok else 1)


# --- main -------------------------------------------------------------------

def main():
    rep = Report()
    for pathogen in PATHOGENS:
        out = os.path.join(OUTPUT_ROOT, "stage4", pathogen)
        if not os.path.exists(os.path.join(out, "21_curation_summary.csv")):
            print(f"\n== {pathogen}: SKIPPED (no step-21 output) ==")
            rep.add(pathogen, "--", "step 21/22 outputs present", False, 1,
                    "no 21_curation_summary.csv")
            continue
        print(f"\n{'=' * 70}\n== {pathogen} ==\n{'=' * 70}")
        print("\n  Step 22 checks:")
        summary22, retained22 = check_step21(pathogen, rep)
        if os.path.exists(os.path.join(out, "22_binarisation_summary.csv")):
            print("\n  Step 23 checks:")
            check_step22(pathogen, summary22, retained22, rep)
        else:
            print(f"\n  Step 23: SKIPPED (no 22_binarisation_summary.csv)")
            rep.add(pathogen, "--", "step 22 output present", False, 1,
                    "no 22_binarisation_summary.csv")

    report_path = os.path.join(DOCS_DIR,
                               f"{datetime.date.today().isoformat()}_step21-22-validation.md")
    rep.write_markdown(report_path)

    print(f"\n{'=' * 70}")
    n_fail = sum(1 for r in rep.rows if not r["passed"])
    print(f"Report written to: {report_path}")
    print(f"Overall: {'ALL PASS' if n_fail == 0 else f'{n_fail} checks FAILED'}")
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
