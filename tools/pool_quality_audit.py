"""
pool_quality_audit.py — audit Stage-4 final pools that fall below the project AUROC threshold.

For every pathogen with Stage-4 output, this reports, per (pathogen, category) and aggregated per
pathogen, how many FINAL pools do not reach the AUROC threshold (AUROC_MIN_THRESHOLD in src/default.py,
0.70) and how much chemical space those sub-threshold pools cover. Two coverage views are given:

  - GROSS      : the union of molecules in the sub-threshold pools, regardless of whether the same
                 molecules are also modelled by a >= threshold pool. "What fraction of chem-space passes
                 through a weak model."
  - EXCLUSIVE  : the subset of those molecules that appear in NO pool that reaches the threshold anywhere
                 in the pathogen (any category). "What fraction of chem-space would lose all modelled
                 coverage if the weak pools were dropped."

Both are expressed against two denominators: the pathogen-CATEGORY chemical space and the TOTAL-pathogen
chemical space (DR u SP unique molecules).

Final pools = steps 25 (leftover-merged modelled pools) and 26 (low-data catch-alls).
  - Step 25 shipped AUROC = grown_auroc if the batch was accepted, else baseline_auroc.
  - Step 26 catch-alls are gated at 0.70: kept ones are >= threshold; those dropped with reason
    "low_auroc" are counted as sub-threshold (their molecules = the category's whole binarised union).
    Catch-alls dropped for "too_small"/"not_cv_able" were never CV-scored (no AUROC) and are reported
    separately as unmodelled, not folded into the sub-threshold coverage math.
  - Fallback: if a pathogen has step-24 modelled pools but no step-25 output (step 25 did not run /
    produced nothing), step-24 pools are used as the modelled set and flagged in the report.

Chemical space per category = union of InChIKeys across output/stage4/<code>/22_binarised/<category>/*.csv.gz
(this equals the step-23 chemspace_partition n_total; computed here directly so the tool also works for
pathogens with no partition file). Total pathogen chem-space = union across all categories.

Outputs (overwritten in place on every run):
  docs/pool_quality_audit.md                                  human-readable report
  output/stage4/general_plots/pool_quality_audit.csv          machine-readable per-row table

Run:  python tools/pool_quality_audit.py
"""

import os
import sys
import glob
import datetime
from collections import defaultdict

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import AUROC_MIN_THRESHOLD  # noqa: E402

REPO = os.path.join(root, "..")
STAGE4 = os.path.join(REPO, "output", "stage4")
DOCS = os.path.join(REPO, "docs")
CSV_OUT = os.path.join(STAGE4, "general_plots", "pool_quality_audit.csv")
REPORT_OUT = os.path.join(DOCS, "pool_quality_audit.md")

THR = AUROC_MIN_THRESHOLD  # "reaches threshold" == auroc >= THR; "below" == auroc < THR

# folders that are not pathogens
NON_PATHOGEN = {"general_plots", "general"}


# --- small IO helpers --------------------------------------------------------
def read_iks(path):
    """Set of inchikeys from a gzipped pool csv (empty set if missing)."""
    if not os.path.exists(path):
        return set()
    return set(pd.read_csv(path, usecols=["inchikey"])["inchikey"])


def category_iks(pdir):
    """{category: set(inchikeys)} from 22_binarised/<cat>/*.csv.gz — the modelable chem-space."""
    out = {}
    bdir = os.path.join(pdir, "22_binarised")
    if not os.path.isdir(bdir):
        return out
    for cat in sorted(os.listdir(bdir)):
        cdir = os.path.join(bdir, cat)
        if not os.path.isdir(cdir):
            continue
        iks = set()
        for f in glob.glob(os.path.join(cdir, "*.csv.gz")):
            iks |= read_iks(f)
        out[cat] = iks
    return out


def pct(n, d):
    return round(100.0 * n / d, 3) if d else 0.0


# --- per-pathogen audit ------------------------------------------------------
def audit_pathogen(code):
    """Return (rows, below_pools, unmodelled, meta) for one pathogen, or None if no chem-space."""
    pdir = os.path.join(STAGE4, code)
    cat_iks = category_iks(pdir)
    if not cat_iks:
        return None
    cat_total = {c: len(s) for c, s in cat_iks.items()}
    pathogen_iks = set().union(*cat_iks.values())
    pathogen_total = len(pathogen_iks)

    good_iks = set()                 # union of molecules in >= threshold final pools (all categories)
    bad_by_cat = defaultdict(set)    # union of molecules in sub-threshold final pools, per category
    below_pools = []                 # one record per sub-threshold pool
    unmodelled = []                  # catch-all categories never CV-scored (no AUROC)

    # --- modelled pools: prefer step 25, else fall back to step 24 ---
    f25 = os.path.join(pdir, "25_pool_summary.csv")
    f24 = os.path.join(pdir, "24_cv_summary.csv")
    if os.path.exists(f25):
        modelled_source = "25"
        df = pd.read_csv(f25)
        for r in df.itertuples(index=False):
            shipped = r.grown_auroc if r.batch_accepted else r.baseline_auroc
            iks = read_iks(os.path.join(pdir, "25_pools", r.category, f"{r.pool_id}.csv.gz"))
            if pd.notna(shipped) and shipped >= THR:
                good_iks |= iks
            else:
                bad_by_cat[r.category] |= iks
                below_pools.append({"category": r.category, "pool_id": r.pool_id, "step": 25,
                                    "auroc": round(float(shipped), 4) if pd.notna(shipped) else None,
                                    "n_molecules": len(iks)})
    elif os.path.exists(f24):
        modelled_source = "24 (step 25 not run)"
        df = pd.read_csv(f24)
        for r in df[df["modelled"]].itertuples(index=False):
            iks = read_iks(os.path.join(pdir, "23_pools", r.category, f"{r.pool_id}.csv.gz"))
            if pd.notna(r.auroc) and r.auroc >= THR:
                good_iks |= iks
            else:
                bad_by_cat[r.category] |= iks
                below_pools.append({"category": r.category, "pool_id": r.pool_id, "step": 24,
                                    "auroc": round(float(r.auroc), 4) if pd.notna(r.auroc) else None,
                                    "n_molecules": len(iks)})
    else:
        modelled_source = "none (catch-alls only)"

    # --- step 26 catch-alls ---
    f26 = os.path.join(pdir, "26_cv_summary.csv")
    if os.path.exists(f26):
        for r in pd.read_csv(f26).itertuples(index=False):
            if bool(r.kept):
                good_iks |= read_iks(os.path.join(pdir, "26_pools", r.category,
                                                  f"{r.category}_catchall.csv.gz"))
            elif str(r.reason) == "low_auroc":
                iks = cat_iks.get(r.category, set())   # catch-all merges the whole category
                bad_by_cat[r.category] |= iks
                below_pools.append({"category": r.category, "pool_id": f"{r.category}_catchall",
                                    "step": 26, "auroc": round(float(r.auroc), 4) if pd.notna(r.auroc) else None,
                                    "n_molecules": len(iks)})
            else:
                unmodelled.append({"category": r.category, "reason": str(r.reason),
                                   "n_molecules": int(r.n_molecules)})

    # --- per-category coverage rows ---
    rows = []
    for cat in sorted(bad_by_cat):
        gross = bad_by_cat[cat]
        excl = gross - good_iks
        n_below = sum(1 for p in below_pools if p["category"] == cat)
        rows.append({
            "pathogen": code, "category": cat, "modelled_source": modelled_source,
            "n_pools_below": n_below,
            "gross_mol": len(gross), "cat_total": cat_total.get(cat, 0), "pathogen_total": pathogen_total,
            "gross_pct_cat": pct(len(gross), cat_total.get(cat, 0)),
            "gross_pct_pathogen": pct(len(gross), pathogen_total),
            "excl_mol": len(excl),
            "excl_pct_cat": pct(len(excl), cat_total.get(cat, 0)),
            "excl_pct_pathogen": pct(len(excl), pathogen_total),
        })

    # --- pathogen aggregate row ---
    gross_all = set().union(*bad_by_cat.values()) if bad_by_cat else set()
    excl_all = gross_all - good_iks
    agg = {
        "pathogen": code, "category": "ALL", "modelled_source": modelled_source,
        "n_pools_below": len(below_pools),
        "gross_mol": len(gross_all), "cat_total": "", "pathogen_total": pathogen_total,
        "gross_pct_cat": "", "gross_pct_pathogen": pct(len(gross_all), pathogen_total),
        "excl_mol": len(excl_all), "excl_pct_cat": "",
        "excl_pct_pathogen": pct(len(excl_all), pathogen_total),
    }
    meta = {"modelled_source": modelled_source, "pathogen_total": pathogen_total,
            "n_good_covered": len(good_iks)}
    return rows, agg, below_pools, unmodelled, meta


# --- report rendering --------------------------------------------------------
def md_table(headers, rows):
    line = lambda cells: "| " + " | ".join(str(c) for c in cells) + " |"
    out = [line(headers), line(["---"] * len(headers))]
    out += [line(r) for r in rows]
    return "\n".join(out)


def render_report(all_rows, all_agg, all_below, all_unmodelled):
    today = datetime.date.today().isoformat()
    n_path = len(all_agg)
    n_affected = sum(1 for a in all_agg if a["n_pools_below"] > 0)
    total_below = sum(a["n_pools_below"] for a in all_agg)

    L = []
    L.append("# Stage-4 pool quality audit — sub-threshold final pools and their chem-space coverage")
    L.append("")
    L.append(f"**Generated:** {today} by `tools/pool_quality_audit.py` (regenerate any time).")
    L.append(f"**AUROC threshold:** `AUROC_MIN_THRESHOLD = {THR}` (from `src/default.py`). "
             f"A pool *reaches* the threshold when its shipped CV AUROC ≥ {THR}; it is *sub-threshold* below that.")
    L.append("")
    L.append(f"Across **{n_path} pathogens**, **{n_affected}** have at least one sub-threshold final pool "
             f"({total_below} such pools in total).")
    L.append("")
    L.append("## What this measures")
    L.append("")
    L.append("Final pools are the shipped datasets from **step 25** (leftover-merged modelled pools) and "
             "**step 26** (low-data catch-alls). For each pathogen, sub-threshold pools are counted and "
             "their molecule coverage reported in two views:")
    L.append("")
    L.append("- **Gross** — union of molecules in the sub-threshold pools, *ignoring* whether the same "
             "molecules are also covered by a ≥-threshold pool. *(What passes through a weak model.)*")
    L.append("- **Exclusive** — the subset in **no** ≥-threshold pool anywhere in the pathogen. *(What "
             "would lose all modelled coverage if the weak pools were dropped.)*")
    L.append("")
    L.append("Each is given against two denominators: the pathogen-**category** chemical space and the "
             "**total-pathogen** chemical space (DR ∪ SP unique InChIKeys, computed from "
             "`22_binarised/`). Step-25 shipped AUROC = `grown_auroc` if the batch was accepted else "
             "`baseline_auroc`. Step-26 catch-alls dropped with reason `low_auroc` count as sub-threshold "
             "(their molecules = the category's whole binarised union); catch-alls dropped as "
             "`too_small`/`not_cv_able` were never scored and are listed separately, not in the coverage math.")
    L.append("")

    # ---- Table 1: per-pathogen summary ----
    L.append("## Per-pathogen summary")
    L.append("")
    L.append("`excl %` / `gross %` are of the **total pathogen** chem-space. `modelled source` flags "
             "pathogens where step 25 did not run.")
    L.append("")
    rows = []
    for a in sorted(all_agg, key=lambda x: (-x["n_pools_below"], x["pathogen"])):
        rows.append([a["pathogen"], a["n_pools_below"], a["gross_mol"], f'{a["gross_pct_pathogen"]}%',
                     a["excl_mol"], f'{a["excl_pct_pathogen"]}%', a["pathogen_total"],
                     a["modelled_source"]])
    L.append(md_table(["pathogen", "#pools <thr", "gross mol", "gross % path", "excl mol",
                       "excl % path", "pathogen chem-space", "modelled source"], rows))
    L.append("")

    # ---- Table 2: per-category detail (affected only) ----
    detail = [r for r in all_rows]
    if detail:
        L.append("## Per-category detail (pathogens with sub-threshold pools)")
        L.append("")
        rows = []
        for r in sorted(detail, key=lambda x: (x["pathogen"], x["category"])):
            rows.append([r["pathogen"], r["category"], r["n_pools_below"],
                         r["gross_mol"], f'{r["gross_pct_cat"]}%', f'{r["gross_pct_pathogen"]}%',
                         r["excl_mol"], f'{r["excl_pct_cat"]}%', f'{r["excl_pct_pathogen"]}%'])
        L.append(md_table(["pathogen", "category", "#pools", "gross mol", "gross % cat",
                           "gross % path", "excl mol", "excl % cat", "excl % path"], rows))
        L.append("")

    # ---- Table 3: the specific sub-threshold pools ----
    flat = [(code, p) for code, plist in all_below for p in plist]
    if flat:
        L.append("## Sub-threshold pools")
        L.append("")
        rows = []
        for code, p in sorted(flat, key=lambda x: (x[0], x[1]["category"], x[1]["pool_id"])):
            rows.append([code, p["category"], p["pool_id"], f'step {p["step"]}',
                         p["auroc"] if p["auroc"] is not None else "—", p["n_molecules"]])
        L.append(md_table(["pathogen", "category", "pool", "source", "AUROC", "molecules"], rows))
        L.append("")

    # ---- Table 4: unmodelled catch-all categories (context) ----
    flatu = [(code, u) for code, ulist in all_unmodelled for u in ulist]
    if flatu:
        L.append("## Catch-all categories never modelled (no AUROC — context, not in coverage math)")
        L.append("")
        L.append("These categories had a catch-all attempt dropped before CV for insufficient data, so "
                 "their molecules have no model at all — a different gap from low-quality models.")
        L.append("")
        rows = [[code, u["category"], u["reason"], u["n_molecules"]]
                for code, u in sorted(flatu, key=lambda x: (x[0], x[1]["category"]))]
        L.append(md_table(["pathogen", "category", "reason", "molecules"], rows))
        L.append("")

    L.append("---")
    L.append("")
    L.append("*Sources: `output/stage4/<code>/{25_pool_summary,24_cv_summary,26_cv_summary}.csv`, the "
             "pool molecule files under `25_pools/`, `23_pools/`, `26_pools/`, and the category "
             "chem-space from `22_binarised/`. Threshold from `src/default.py`. Machine-readable table: "
             "`output/stage4/general_plots/pool_quality_audit.csv`.*")
    return "\n".join(L) + "\n"


# --- main --------------------------------------------------------------------
def main():
    codes = sorted(d for d in os.listdir(STAGE4)
                   if os.path.isdir(os.path.join(STAGE4, d)) and d not in NON_PATHOGEN)

    all_rows, all_agg, all_below, all_unmodelled = [], [], [], []
    for code in codes:
        res = audit_pathogen(code)
        if res is None:
            print(f"  {code}: no Stage-4 chem-space (skipped)")
            continue
        rows, agg, below, unmodelled, meta = res
        all_rows.extend(rows)
        all_agg.append(agg)
        all_below.append((code, below))
        all_unmodelled.append((code, unmodelled))
        print(f"  {code}: {agg['n_pools_below']} sub-threshold pools | "
              f"gross {agg['gross_mol']} ({agg['gross_pct_pathogen']}% of pathogen), "
              f"exclusive {agg['excl_mol']} ({agg['excl_pct_pathogen']}%) | source {meta['modelled_source']}")

    # CSV (per-category rows + aggregate rows)
    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    cols = ["pathogen", "category", "modelled_source", "n_pools_below", "gross_mol", "cat_total",
            "pathogen_total", "gross_pct_cat", "gross_pct_pathogen", "excl_mol", "excl_pct_cat",
            "excl_pct_pathogen"]
    pd.DataFrame(all_rows + all_agg, columns=cols).to_csv(CSV_OUT, index=False)

    # Markdown report
    os.makedirs(DOCS, exist_ok=True)
    with open(REPORT_OUT, "w") as fh:
        fh.write(render_report(all_rows, all_agg, all_below, all_unmodelled))

    print(f"\nReport : {os.path.relpath(REPORT_OUT, REPO)}")
    print(f"Table  : {os.path.relpath(CSV_OUT, REPO)}")


if __name__ == "__main__":
    main()
