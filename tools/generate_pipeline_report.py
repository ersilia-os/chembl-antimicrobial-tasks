"""
Cross-pathogen pipeline status report.

Produces one CSV row per pathogen summarising how molecules and assays move through
the ChEMBL antimicrobial-tasks pipeline: how many are assessed, how they are classified
(SINGLE PROTEIN / ORGANISM / discarded), and how many survive to the final selected
A/B/M datasets. Intended as a before/after snapshot for pipeline-improvement tracking.

This is a reporting utility (not part of the numbered per-pathogen pipeline / run_all.sh).
It is read-only over output/ and config/.

Usage
-----
    python tools/generate_pipeline_report.py [output_csv]

Default output: docs/260622_report.csv

Key definitions (confirmed with the user)
------------------------------------------
- Assay counting unit: each (assay_id, activity_type, unit) row in 18_assays_master.csv.
- Molecules: cleaned activities (08), deduplicated by compound_chembl_id.
- ORG actives/inactives: binarized datasets (12_datasets bin column). The ORG totals are
  taken at the *reference* cutoff (reference_cutoff in src/pathogen_utils.py, the same
  convention used by steps 20/21). "Kept" = molecules present in a selected dataset.
- "Kept" datasets: selected == True in 17_final_datasets.csv (post correlation-dedup).
  M datasets pool several assays, counted via their assay_keys.
"""

import os
import sys
import io
import csv
import gzip
import zipfile

import pandas as pd
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from pathogen_utils import load_pathogen, load_expert_cutoffs, reference_cutoff  # noqa: E402

CONFIGPATH = os.path.join(root, "..", "config")
OUTPUT_ROOT = os.path.join(root, "..", "output")
DOCS_DIR = os.path.join(root, "..", "docs")
os.makedirs(DOCS_DIR, exist_ok=True)

DEFAULT_OUT = os.path.join(DOCS_DIR, "260622_report.csv")

# Sentinel used so that NaN units/activity types match consistently across tables.
NA = "<NA>"

COLUMNS = [
    "pathogen_code",
    "n_total_mols",
    "n_unique_mols",
    "n_total_mols_org",
    "n_unique_mols_org",
    "n_assays_org",
    "n_assays_sp",
    "n_assays_discarded_class",
    "n_unique_actives_org",
    "n_unique_inactives_org",
    "n_assays_a",
    "n_assays_b",
    "n_assays_m",                # number of selected M (merged) datasets
    "n_assays_underlying_m",     # number of underlying assays pooled into selected M datasets
    "n_assays_not_kept",
    "frac_discarded",
    "frac_unique_mols_discarded",
    "frac_unique_actives_discarded",
    "n_strains",
    # ---- extra metrics ----
    "n_assays_raw",
    "n_assays_cleaned",
    "frac_unique_inactives_discarded",
    "mean_auroc_selected",
    "pct_total_chem_space_covered",   # kept ORG mols / full 07 compound universe (step-21 base)
    "pct_org_chem_space_covered",     # kept ORG mols / all ORG-assessed unique mols (incl. ql, none)
    "pct_org_qtmx_chem_space_covered",  # kept ORG mols / quantitative+mixed ORG mols (coverable universe)
    # ---- small ORG assays: < 20 compounds (deduplicated, ORG-only) ----
    "n_assays_org_lt20", "pct_assays_org_lt20",
    "mols_org_lt20", "pct_chemspace_org_lt20",
    "mols_org_excl_lt20", "pct_chemspace_org_excl_lt20",
    # ---- qualitative-only ORG assays (excluded from A/B and from merging) ----
    "n_assays_org_ql", "pct_assays_org_ql",
    "mols_org_ql", "pct_chemspace_org_ql",
    "mols_org_excl_ql", "pct_chemspace_org_excl_ql",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _norm(series):
    """Normalise a key column: strip, fill NaN with a sentinel string."""
    return series.astype("string").fillna(NA).str.strip()


def _parse_bins(data_bytes):
    """Parse gzipped-CSV bytes into a list of (compound_chembl_id, is_active) pairs.

    Uses the csv module (fast + quote-safe) rather than pandas, because we read tens of
    thousands of tiny dataset files. The `bin` column is located by header name (it is not
    always the last column). Returns [] on any malformed/empty input.
    """
    try:
        text = gzip.decompress(data_bytes).decode("utf-8")
    except (OSError, EOFError, UnicodeDecodeError):
        return []
    reader = csv.reader(io.StringIO(text))
    try:
        header = next(reader)
    except StopIteration:
        return []
    try:
        ci = header.index("compound_chembl_id")
        bi = header.index("bin")
    except ValueError:
        return []
    pairs = []
    hi = max(ci, bi)
    for row in reader:
        if len(row) <= hi:
            continue
        pairs.append((row[ci], row[bi] == "1"))
    return pairs


def _bins_from_zip(zf, names, member):
    """Read (compound_chembl_id, is_active) pairs from a member of an open ZipFile."""
    if member not in names:
        return []
    return _parse_bins(zf.read(member))


def _bins_from_file(path):
    """Read (compound_chembl_id, is_active) pairs from a .csv.gz file on disk."""
    if not os.path.exists(path):
        return []
    with open(path, "rb") as fh:
        return _parse_bins(fh.read())


def _cutoff_str(cut):
    """Format a numeric cutoff the way step-12 writes it into filenames (e.g. 10.0)."""
    return str(float(cut))


def _union(sets):
    """Union of an iterable of sets (empty set if none)."""
    out = set()
    for s in sets:
        out |= s
    return out


def _fill_bucket(out, tag, n_assays_bucket, union_set, excl_set, n_assays_org, n_unique_org):
    """Populate the six per-bucket columns (assay count + chem-space share) for a tag."""
    out[f"n_assays_org_{tag}"] = int(n_assays_bucket)
    out[f"pct_assays_org_{tag}"] = round(100 * n_assays_bucket / n_assays_org, 1) if n_assays_org else 0
    out[f"mols_org_{tag}"] = len(union_set)
    out[f"pct_chemspace_org_{tag}"] = round(100 * len(union_set) / n_unique_org, 1) if n_unique_org else 0
    out[f"mols_org_excl_{tag}"] = len(excl_set)
    out[f"pct_chemspace_org_excl_{tag}"] = round(100 * len(excl_set) / n_unique_org, 1) if n_unique_org else 0


# ---------------------------------------------------------------------------
# Per-pathogen metric computation
# ---------------------------------------------------------------------------


def compute_pathogen(code, expert):
    out = {c: 0 for c in COLUMNS}
    out["pathogen_code"] = code
    for fcol in ["frac_discarded", "frac_unique_mols_discarded", "frac_unique_actives_discarded",
                 "frac_unique_inactives_discarded", "mean_auroc_selected",
                 "pct_total_chem_space_covered", "pct_org_chem_space_covered",
                 "pct_org_qtmx_chem_space_covered"]:
        out[fcol] = np.nan

    pdir = os.path.join(OUTPUT_ROOT, code)
    if not os.path.isdir(pdir):
        print(f"  [{code}] output folder missing — emitting zeros")
        return out

    master_path = os.path.join(pdir, "18_assays_master.csv")
    cleaned_path = os.path.join(pdir, "08_chembl_cleaned_data.csv.gz")
    datasets_path = os.path.join(pdir, "12_datasets.csv")
    final_path = os.path.join(pdir, "17_final_datasets.csv")

    if not os.path.exists(master_path):
        print(f"  [{code}] 18_assays_master.csv missing — emitting zeros")
        return out

    master = pd.read_csv(master_path, low_memory=False)
    master["unit_k"] = _norm(master["unit"])
    master["activity_k"] = _norm(master["activity_type"])
    master["assay_k"] = _norm(master["assay_id"])

    # ---- Classification counts (per assay x activity x unit row) ----
    tt = master["target_type_curated_extra"]
    out["n_assays_org"] = int((tt == "ORGANISM").sum())
    out["n_assays_sp"] = int((tt == "SINGLE PROTEIN").sum())
    out["n_assays_discarded_class"] = int((tt == "DISCARDED").sum() + tt.isna().sum())

    # ORG (assay, activity, unit) keys
    org_mask = tt == "ORGANISM"
    org_keys = set(
        map(tuple, master.loc[org_mask, ["assay_k", "activity_k", "unit_k"]].drop_duplicates().values)
    )

    # ---- Strains (ORG assays only) ----
    strain_col = "assay_strain_norm" if "assay_strain_norm" in master.columns else "assay_strain_curated"
    if strain_col in master.columns:
        org_strains = master.loc[org_mask, strain_col].dropna()
        org_strains = org_strains[org_strains.astype(str).str.strip() != ""]
        out["n_strains"] = int(org_strains.nunique())

    # ---- Funnel context: raw / cleaned assay counts (per assay x activity x unit) ----
    raw_assays_path = os.path.join(pdir, "07_assays_raw.csv")
    cleaned_assays_path = os.path.join(pdir, "08_assays_cleaned.csv")
    if os.path.exists(raw_assays_path):
        out["n_assays_raw"] = int(len(pd.read_csv(raw_assays_path, usecols=["assay_id"])))
    if os.path.exists(cleaned_assays_path):
        out["n_assays_cleaned"] = int(len(pd.read_csv(cleaned_assays_path, usecols=["assay_id"])))

    # ---- Molecule counts from cleaned activities ----
    org_assessed = set()
    qtmx_u = set()   # unique molecules in quantitative+mixed ORG assays (coverable universe)
    if os.path.exists(cleaned_path):
        cleaned = pd.read_csv(
            cleaned_path,
            usecols=["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"],
            low_memory=False,
        )
        out["n_total_mols"] = int(len(cleaned))
        out["n_unique_mols"] = int(cleaned["compound_chembl_id"].nunique())

        # 08 cleaned data keys assays by the numeric assay_id; the master / dataset tables
        # use the ChEMBL string id (assay_chembl_id), so join on that.
        cleaned["unit_k"] = _norm(cleaned["unit"])
        cleaned["activity_k"] = _norm(cleaned["activity_type"])
        cleaned["assay_k"] = _norm(cleaned["assay_chembl_id"])
        keytups = list(map(tuple, cleaned[["assay_k", "activity_k", "unit_k"]].values))
        is_org = pd.Series([k in org_keys for k in keytups], index=cleaned.index)
        org_rows = cleaned[is_org].copy()
        out["n_total_mols_org"] = int(len(org_rows))
        out["n_unique_mols_org"] = int(org_rows["compound_chembl_id"].nunique())
        org_assessed = set(org_rows["compound_chembl_id"].unique())

        # ---- ORG assay-size breakdown + qualitative-only chem space ----
        # compound set per ORG (assay, activity, unit)
        org_rows["_k"] = list(zip(org_rows["assay_k"], org_rows["activity_k"], org_rows["unit_k"]))
        cset = org_rows.groupby("_k")["compound_chembl_id"].apply(set)
        cset_d = cset.to_dict()
        sizes = cset.map(len)
        n_org_assays = out["n_assays_org"]
        n_org_unique = out["n_unique_mols_org"]

        small_u = _union(cset[sizes < 20].values)
        large_u = _union(cset[sizes >= 20].values)   # kept internally to find small-exclusive mols
        _fill_bucket(out, "lt20", int((sizes < 20).sum()), small_u, small_u - large_u,
                     n_org_assays, n_org_unique)

        # qualitative-only ORG assays: compounds reachable only via qualitative data
        # (A/B and merging both consume quantitative/mixed only).
        dt = master.loc[org_mask, "dataset_type"]
        org_master_keys = list(zip(master.loc[org_mask, "assay_k"],
                                   master.loc[org_mask, "activity_k"],
                                   master.loc[org_mask, "unit_k"]))
        ql_keys = {k for k, t in zip(org_master_keys, dt) if t == "qualitative"}
        qtmx_keys = {k for k, t in zip(org_master_keys, dt) if t in ("quantitative", "mixed")}
        ql_u = _union(cset_d[k] for k in ql_keys if k in cset_d)
        qtmx_u = _union(cset_d[k] for k in qtmx_keys if k in cset_d)
        _fill_bucket(out, "ql", len(ql_keys), ql_u, ql_u - qtmx_u, n_org_assays, n_org_unique)
    else:
        print(f"  [{code}] 08 cleaned data missing — molecule counts = 0")

    # ---- ORG actives / inactives from binarized datasets (reference cutoff) ----
    active_org = set()
    binarized_org = set()
    if os.path.exists(datasets_path):
        dsets = pd.read_csv(datasets_path, low_memory=False)
        dsets = dsets[dsets["target_type_curated_extra"] == "ORGANISM"]
        suffix_map = {"quantitative": "qt", "mixed": "mx", "qualitative": "ql"}
        groups = (
            dsets[dsets["dataset_type"].isin(suffix_map)]
            [["assay_id", "activity_type", "unit", "dataset_type"]]
            .drop_duplicates()
        )
        # Cutoffs actually present per (assay, activity, unit), precomputed once
        # (avoids an O(groups x rows) scan inside the loop).
        avail_cuts = (
            dsets.dropna(subset=["expert_cutoff"])
            .groupby(["assay_id", "activity_type", "unit"])["expert_cutoff"]
            .apply(lambda s: sorted(set(s)))
            .to_dict()
        )

        # Build the reference-cutoff filename for each ORG group, grouped by zip.
        targets = {"qt": [], "mx": [], "ql": []}
        for assay, act, unit, dtype in groups.itertuples(index=False):
            suffix = suffix_map[dtype]
            if suffix == "ql":
                targets["ql"].append(f"{assay}_{act}_ql.csv.gz")
                continue
            cut_list = expert.get((act, unit, "ORGANISM", code))
            if cut_list:
                cut = reference_cutoff(cut_list, unit)
            else:  # fallback: middle of cutoffs present for this group
                avail = avail_cuts.get((assay, act, unit), [])
                cut = avail[len(avail) // 2] if avail else None
            if cut is None:
                continue
            targets[suffix].append(f"{assay}_{act}_{suffix}_{_cutoff_str(cut)}.csv.gz")

        for suffix in ["qt", "mx", "ql"]:
            if not targets[suffix]:
                continue
            zpath = os.path.join(pdir, "12_datasets", f"datasets_{suffix}.zip")
            if not os.path.exists(zpath):
                continue
            with zipfile.ZipFile(zpath) as zf:
                names = set(zf.namelist())
                for fname in targets[suffix]:
                    for cid, is_active in _bins_from_zip(zf, names, fname):
                        # Restrict to molecules actually assessed against the pathogen
                        # (drops any decoys / non-assessed compounds).
                        if cid not in org_assessed:
                            continue
                        binarized_org.add(cid)
                        if is_active:
                            active_org.add(cid)

        out["n_unique_actives_org"] = int(len(active_org))
        out["n_unique_inactives_org"] = int(len(binarized_org - active_org))
    else:
        print(f"  [{code}] 12_datasets.csv missing — actives/inactives = 0")

    # ---- Kept (selected) datasets: A / B / M counts and kept molecule sets ----
    kept_keys = set()       # (assay, activity, unit) triples retained in A/B/M
    kept_mols = set()       # ORG-assessed molecules present in any selected dataset
    kept_actives = set()    # ORG-active molecules present in any selected dataset
    if os.path.exists(final_path):
        final = pd.read_csv(final_path)
        if "selected" in final.columns and len(final) > 0:
            final["selected"] = final["selected"].astype(bool)
            sel = final[final["selected"]]
            out["n_assays_a"] = int((sel["label"] == "A").sum())
            out["n_assays_b"] = int((sel["label"] == "B").sum())
            out["n_assays_m"] = int((sel["label"] == "M").sum())   # number of selected M datasets
            if "auroc" in sel.columns and len(sel) > 0:
                aurocs = sel["auroc"].dropna()
                if len(aurocs) > 0:
                    out["mean_auroc_selected"] = round(float(aurocs.mean()), 3)

            m_keys = set()
            for _, r in sel.iterrows():
                for k in str(r["assay_keys"]).split(";"):
                    parts = k.split("|")
                    if len(parts) == 3:
                        a, act, unit = (p.strip() for p in parts)
                        unit = unit if unit and unit.lower() != "nan" else NA
                        trip = (a, act, unit)
                        kept_keys.add(trip)
                        if r["label"] == "M":
                            m_keys.add(trip)
            # underlying assays pooled across all selected M datasets
            out["n_assays_underlying_m"] = int(len(m_keys))

            # Read each selected dataset file to collect kept molecules / actives.
            # Group names by location so each zip is opened only once.
            by_zip = {"qt": [], "mx": [], "ql": []}
            m_names = []
            for name in sel["name"].astype(str):
                if name.startswith("M_"):
                    m_names.append(name)
                elif "_qt_" in name:
                    by_zip["qt"].append(name)
                elif "_mx_" in name:
                    by_zip["mx"].append(name)
                else:
                    by_zip["ql"].append(name)

            def _absorb(pairs):
                for cid, is_active in pairs:
                    if cid not in org_assessed:
                        continue
                    kept_mols.add(cid)
                    if is_active:
                        kept_actives.add(cid)

            for name in m_names:
                _absorb(_bins_from_file(os.path.join(pdir, "12_datasets", "M", name + ".csv.gz")))
            for suffix, names_list in by_zip.items():
                if not names_list:
                    continue
                zpath = os.path.join(pdir, "12_datasets", f"datasets_{suffix}.zip")
                if not os.path.exists(zpath):
                    continue
                with zipfile.ZipFile(zpath) as zf:
                    znames = set(zf.namelist())
                    for name in names_list:
                        _absorb(_bins_from_zip(zf, znames, name + ".csv.gz"))
    else:
        print(f"  [{code}] 17_final_datasets.csv missing — A/B/M counts = 0")

    # ---- Discarded / kept fractions ----
    out["n_assays_not_kept"] = int(out["n_assays_org"] - len(kept_keys))
    if out["n_assays_org"] > 0:
        out["frac_discarded"] = round(out["n_assays_not_kept"] / out["n_assays_org"], 4)

    kept_mols_org = kept_mols & org_assessed
    if out["n_unique_mols_org"] > 0:
        out["frac_unique_mols_discarded"] = round(
            (out["n_unique_mols_org"] - len(kept_mols_org)) / out["n_unique_mols_org"], 4
        )
    if out["n_unique_actives_org"] > 0:
        kept_a = len(kept_actives & active_org)  # guarantee subset of total actives
        out["frac_unique_actives_discarded"] = round(
            (out["n_unique_actives_org"] - kept_a) / out["n_unique_actives_org"], 4
        )
    if out["n_unique_inactives_org"] > 0:
        inactive_org = binarized_org - active_org
        kept_i = len(kept_mols_org & inactive_org)  # pure inactives retained in a selected dataset
        out["frac_unique_inactives_discarded"] = round(
            (out["n_unique_inactives_org"] - kept_i) / out["n_unique_inactives_org"], 4
        )

    # ---- Chemical-space coverage (same numerator, two denominators) ----
    # pct_total: kept ORG mols over the full 07_compound_counts universe (the set step 21's
    #   funnel uses) — % of ALL pathogen compounds covered by kept datasets.
    # pct_org: kept ORG mols over ORG-assessed unique mols — exact complement of
    #   frac_unique_mols_discarded.
    counts_path = os.path.join(pdir, "07_compound_counts.csv.gz")
    if os.path.exists(counts_path):
        universe = pd.read_csv(counts_path, usecols=["compound_chembl_id"])["compound_chembl_id"].nunique()
        if universe > 0:
            out["pct_total_chem_space_covered"] = round(100 * len(kept_mols_org) / universe, 2)
    if out["n_unique_mols_org"] > 0:
        out["pct_org_chem_space_covered"] = round(100 * len(kept_mols_org) / out["n_unique_mols_org"], 2)
    # coverage restricted to the coverable universe (quantitative+mixed ORG molecules)
    if len(qtmx_u) > 0:
        out["pct_org_qtmx_chem_space_covered"] = round(100 * len(kept_mols_org & qtmx_u) / len(qtmx_u), 2)

    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    out_csv = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT
    expert = load_expert_cutoffs(CONFIGPATH)
    pathogens = pd.read_csv(os.path.join(CONFIGPATH, "pathogens.csv"))

    rows = []
    for code in pathogens["code"]:
        name = load_pathogen(code)
        print(f"Processing {code} ({name}) ...")
        rows.append(compute_pathogen(code, expert))

    report = pd.DataFrame(rows, columns=COLUMNS)
    report.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv} ({len(report)} pathogens).")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(report.to_string(index=False))


if __name__ == "__main__":
    main()
