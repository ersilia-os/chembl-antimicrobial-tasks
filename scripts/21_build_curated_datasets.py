"""
Step 21 - Build curated per-(assay_id, activity_type) datasets, then plot.

Two stages, run back-to-back per pathogen:

1. BUILD (`build()`): reads the existing per-pathogen outputs (steps 07/08/09) and
   produces the curated, deduplicated, quantitative-only datasets that feed the
   binarisation/pooling stages (first step of the Stage-4 signal branch).
2. PLOT (`process()`): per-pathogen curation attrition & characterization figures.
   Cross-pathogen comparison figures now live in scripts/27_general_plots.py.

Inputs (output/stage2/<pathogen>/):
    08_chembl_cleaned_data.csv.gz   activity records (incl. value, relation, unit,
                                    direction, text_flag, doc_chembl_id)
    09_assays_parameters_full.csv   assay-level curated fields (source_label, target/strain
                                    curation, atcc_id, mutations, ...) — Stage 2, so Stage 4
                                    does not depend on the Stage-3 assays_master table
    07_compound_counts.csv.gz       compound_chembl_id -> inchikey map
    config/pubchem_aids/...         ChEMBL assays that have a PubChem counterpart

Config-derived fields are recomputed live here rather than trusted from the cached
upstream files, so a curator's config edit propagates without re-running steps 08/19:
    - direction                 re-mapped from config/activity_std_units_manual_curation.csv
                                by (activity_type, unit) (mirrors step 08).
    - assay_strain_norm         recomputed via strain_norm_key() + config/strains/ (step 19).
    - target_type_curated_extra recomputed via extra_curation_target_type() (step 13).

Curation (one unit of work = one (assay_id, activity_type) pair):
    - InChIKey dedup, keeping the single most active record per molecule.
    - Step 1: keep only ORGANISM (target_type_curated_extra). SP/other kept in the
      summary, flagged discarded.
    - Step 2: keep only pairs with numeric (quantitative) values. If a pair has both
      numeric and text records, the numeric ones are used; text-flag-only *pairs* are
      discarded.
    - Step 3: discard pairs whose numeric records span >1 unit (incompatible, not pooled).
    - Step 4: discard UNCATEGORIZED pairs — (activity_type, unit) not in the DR/SP
      inclusion list (config/activity_type_categories.csv); no defined binarisation cutoff.
    - Step 5: discard pairs with an undefined biological direction (not exactly one of
      {-1, +1}); a valid single direction is required for the dedup and step-22 binarisation.
    - Step 6: discard CO-ADD (source_label).
    - Step 7: discard PubChem-superseded assays.
    - Step 8: discard pairs with <= MIN_MOLECULES unique molecules, counting the
      *dataset* size n_mol_dataset (measured molecules + recovered inactives, below).

Recovered qualitative inactives (curator decision, 2026-07-01): within a retained
pair, a molecule (InChIKey) with NO numeric record is carried into the dataset as an
INACTIVE iff its qualitative records are purely inactive (text_flag == -1, no active
flag). A molecule with ANY active flag but no value is discarded (can't be placed
without a measurement); flagless value-less records carry no signal and stay dropped.
Recovered inactives count toward the MIN_MOLECULES size threshold and are written with
value/relation = NaN; step 23 labels them inactive via binarize_conservative.

Outputs (output/stage4/<pathogen>/):
    21_datasets/<assay_id>_<activity_type>.csv.gz   one file per retained pair
    21_curation_summary.csv                         one row per pair (all pairs,
                                                    incl. discarded + reason)
    21_curation_overview.png                        per-pathogen overview figure
    21_curation_attrition_waterfall.png             per-pathogen attrition waterfall
    21_curation_stats.csv                           per-pathogen attrition stats
    21_curation_distributions.csv                   molecule-reuse & datasets-per-document
                                                    distributions (the two figure panels
                                                    not reconstructable from the CSVs above)
    21_curation_categories.csv                      retained DR/SP split by measure
                                                    (datasets, molecules, chem space) — H panel

Attrition uses the ordered curation stages of the build. Because the build records,
per (assay_id, activity_type) pair ("assay" in the v2 nomenclature), the *first*
failing step (`discard_step`; 0 = retained), a pair "survives stage k" iff
`discard_step == 0 or discard_step > k`. Stages are nested, so the metrics
decrease monotonically.

All plotting goes through stylia in the publication house style (print format, article
palette), applied via src/plot_base.py; colours come from the shared semantic map in
src/plot_colors.py so DR/SP and active/inactive stay consistent across every figure.

Runs on ONE pathogen per invocation (like steps 07-20); loop over pathogens with
scripts/run_step.sh. Cross-pathogen comparison figures are a separate step (23).

Usage:
    python scripts/21_build_curated_datasets.py <pathogen_code>
"""

import os
import sys
import glob
from collections import Counter

import numpy as np
import pandas as pd
from tqdm import tqdm

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from pathogen_utils import (
    load_pathogen, load_pubchem_assays,
    strain_norm_key, load_strain_equivalences, extra_curation_target_type,
)
from default import CONFIGPATH, MIN_MOLECULES, COADD_SOURCE_LABEL, KEPT_TARGET_TYPE, STAGE4PATH
from dataset_utils import (
    assign_activity_category, load_activity_categories, dedup_most_active_by_inchikey,
)
from plot_utils import (
    curation_stats, cat_label,
    CATEGORY_ORDER, OVERVIEW_DROP,
)
import plot_colors as pc
from plot_base import BasePlot, render_figure, figure_legend

# --- build module-level constants -------------------------------------------
KEY = ["assay_chembl_id", "activity_type"]

# Assay-level curated fields carried over from step 09 (Stage 2; constant per assay).
ENRICH_COLS = [
    "target_type_curated_extra", "source_label", "doc_chembl_id",
    "assay_organism", "assay_organism_curated", "target_organism",
    "target_type", "target_type_curated", "target_chembl_id_curated",
    "target_name_curated", "bao_label", "assay_strain_curated",
    "assay_strain_norm", "atcc_id", "mutations",
]
# The 2 config-derived cols (target_type_curated_extra, assay_strain_norm) are recomputed live
# below, so they are NOT pulled from the source table — the other 13 come from step 09 (Stage 2).
ENRICH_SOURCE_COLS = [c for c in ENRICH_COLS
                      if c not in ("target_type_curated_extra", "assay_strain_norm")]


# --- BUILD ------------------------------------------------------------------

def build(pathogen_code):
    pathogen = load_pathogen(pathogen_code)

    OUTPUT = os.path.join(STAGE4PATH, pathogen_code)
    stage2_dir = os.path.join(root, "..", "output", "stage2", pathogen_code)
    DATASETS_DIR = os.path.join(OUTPUT, "21_datasets")
    os.makedirs(DATASETS_DIR, exist_ok=True)
    # Clear prior dataset files so the directory always matches the current run.
    for f in glob.glob(os.path.join(DATASETS_DIR, "*.csv.gz")):
        os.remove(f)

    print(f"Step 21 (build) - {pathogen_code} ({pathogen})")

    # --- load ---------------------------------------------------------------
    acts = pd.read_csv(os.path.join(stage2_dir, "08_chembl_cleaned_data.csv.gz"), low_memory=False)
    params = pd.read_csv(os.path.join(stage2_dir, "09_assays_parameters_full.csv"), low_memory=False)
    counts = pd.read_csv(os.path.join(stage2_dir, "07_compound_counts.csv.gz"),
                         usecols=["compound_chembl_id", "inchikey"])
    pubchem_ids = load_pubchem_assays(pathogen_code)

    cid2ik = dict(zip(counts["compound_chembl_id"], counts["inchikey"]))
    acts["inchikey"] = acts["compound_chembl_id"].map(cid2ik)
    acts["is_numeric"] = acts["value"].notna()
    acts["is_text"] = acts["text_flag"].isin([1, -1])

    # LIVE recompute of direction from config (do NOT trust cached 08 value).
    # Mirrors scripts/08_clean_pathogen_activities.py:46-48; config overrides, else keep
    # the cached 08 direction so we never drop a direction step 08 had (NaN-safe fallback).
    _dir_cfg = pd.read_csv(os.path.join(CONFIGPATH, "activity_std_units_manual_curation.csv"))
    _dir_map = {(i, j): k for i, j, k in zip(
        _dir_cfg["activity_type"], _dir_cfg["unit"], _dir_cfg["manual_curation_direction"])
        if np.isnan(k) == False}
    _remapped = pd.Series(
        [_dir_map.get((i, j), np.nan) for i, j in zip(acts["activity_type"], acts["unit"])],
        index=acts.index)
    acts["direction"] = _remapped.where(_remapped.notna(), acts["direction"])

    n_missing_ik = int(acts["inchikey"].isna().sum())
    print(f"  activities: {len(acts)} | assays: {acts['assay_chembl_id'].nunique()} "
          f"| compounds w/o InChIKey: {n_missing_ik}")
    print(f"  PubChem-superseded assays listed: {len(pubchem_ids)}")

    # --- per-pair aggregation -----------------------------------------------
    num = acts[acts["is_numeric"]]

    # Recover qualitative-inactive molecules (curator decision, 2026-07-01):
    # within a pair, a molecule (InChIKey) with NO numeric record is kept as an
    # inactive iff its qualitative records are purely inactive (text_flag == -1,
    # no active flag). Any active flag -> discarded (can't place it without a
    # value). Flagless value-less records carry no signal and stay dropped.
    acts_ik = acts[acts["inchikey"].notna()].copy()
    acts_ik["_active_flag"] = acts_ik["text_flag"] == 1
    acts_ik["_inactive_flag"] = acts_ik["text_flag"] == -1
    molflags = (acts_ik.groupby(KEY + ["inchikey"])
                .agg(has_num=("is_numeric", "any"),
                     has_active=("_active_flag", "any"),
                     has_inactive=("_inactive_flag", "any"))
                .reset_index())
    recov = molflags[(~molflags["has_num"]) & molflags["has_inactive"] & (~molflags["has_active"])]
    disc_active = molflags[(~molflags["has_num"]) & molflags["has_active"]]
    n_recov = recov.groupby(KEY).size().rename("n_recovered_inactive")
    n_disc_active = disc_active.groupby(KEY).size().rename("n_text_only_active_discarded")
    # one representative record (compound_chembl_id, smiles) per recovered InChIKey
    rep = acts_ik.drop_duplicates(KEY + ["inchikey"])[KEY + ["inchikey", "compound_chembl_id", "smiles"]]
    recov_rows = recov[KEY + ["inchikey"]].merge(rep, on=KEY + ["inchikey"], how="left")
    recov_by_key = {k: g for k, g in recov_rows.groupby(KEY)}

    pairs = (acts.groupby(KEY)
             .agg(n_activities=("activity_id", "size"),
                  n_numeric=("is_numeric", "sum"),
                  n_text=("is_text", "sum"),
                  n_mol_all=("inchikey", "nunique"))
             .reset_index())

    # numeric-only stats per pair
    num_units = (num.groupby(KEY)["unit"]
                 .agg(lambda s: sorted(set(s.dropna().astype(str)))))
    num_dir = (num.groupby(KEY)["direction"]
               .agg(lambda s: sorted(set(s.dropna()))))
    num_nmol = num.groupby(KEY)["inchikey"].nunique()

    pairs = pairs.set_index(KEY)
    pairs["numeric_units"] = num_units
    pairs["numeric_directions"] = num_dir
    pairs["n_mol_numeric"] = num_nmol
    pairs = pairs.reset_index()

    pairs["numeric_units"] = pairs["numeric_units"].apply(lambda x: x if isinstance(x, list) else [])
    pairs["numeric_directions"] = pairs["numeric_directions"].apply(lambda x: x if isinstance(x, list) else [])
    pairs["n_mol_numeric"] = pairs["n_mol_numeric"].fillna(0).astype(int)
    pairs["n_numeric_units"] = pairs["numeric_units"].apply(len)

    # recovered qualitative-inactive counts, and the dataset size that now drives
    # the MIN_MOLECULES filter (measured molecules + recovered inactives).
    pairs = pairs.set_index(KEY)
    pairs["n_recovered_inactive"] = n_recov
    pairs["n_text_only_active_discarded"] = n_disc_active
    pairs = pairs.reset_index()
    pairs["n_recovered_inactive"] = pairs["n_recovered_inactive"].fillna(0).astype(int)
    pairs["n_text_only_active_discarded"] = pairs["n_text_only_active_discarded"].fillna(0).astype(int)
    pairs["n_mol_dataset"] = pairs["n_mol_numeric"] + pairs["n_recovered_inactive"]

    # single unit / direction (only meaningful when exactly one is present)
    pairs["unit"] = pairs["numeric_units"].apply(lambda u: u[0] if len(u) == 1 else np.nan)
    pairs["direction"] = pairs["numeric_directions"].apply(lambda d: d[0] if len(d) == 1 else np.nan)
    # category from the (activity_type, unit) inclusion list (config-driven, user-reviewed)
    activity_categories = load_activity_categories(CONFIGPATH)
    pairs["category"] = [assign_activity_category(a, u, activity_categories)
                         for a, u in zip(pairs["activity_type"], pairs["unit"])]

    # flags surfaced for review (not in themselves the discard reason, except multi-unit)
    pairs["flag_multi_unit"] = pairs["n_numeric_units"] > 1
    pairs["flag_uncategorized"] = (pairs["n_numeric"] > 0) & (pairs["category"] == "UNCATEGORIZED")

    # --- enrichment from step-09 assay parameters (assay-level) -------------
    ma = params[["assay_id"] + ENRICH_SOURCE_COLS].drop_duplicates("assay_id")
    pairs = pairs.merge(ma, left_on="assay_chembl_id", right_on="assay_id", how="left")
    pairs = pairs.drop(columns=["assay_id"])

    # LIVE recompute of the 2 config-derived assay-level fields (not carried from 09).
    pairs["assay_strain_norm"] = strain_norm_key(
        pairs["assay_strain_curated"], pairs["atcc_id"],
        load_strain_equivalences(pathogen_code))            # mirrors step 19:47-49
    pairs["target_type_curated_extra"] = [
        extra_curation_target_type(tt, ttc)
        for tt, ttc in zip(pairs["target_type"], pairs["target_type_curated"])
    ]                                                       # mirrors step 13:33-35

    # --- ordered curation filters -------------------------------------------
    # discard_step: 0 = retained; otherwise the first failing step (reason recorded).
    # Steps 4 (uncategorized) and 5 (undefined direction) are endpoint-quality gates:
    # they run after single-unit (unit-dependent) and before the source/size gates, so
    # every retained pair is DR/SP with a defined single direction (guards the write-step
    # dedup and step 22's binarisation).
    def decide(r):
        if r["target_type_curated_extra"] != KEPT_TARGET_TYPE:
            return 1, f"not ORGANISM ({r['target_type_curated_extra']})"
        if r["n_numeric"] == 0:
            return 2, "no quantitative data (qualitative-only)"
        if r["n_numeric_units"] > 1:
            return 3, "multiple units (incompatible, not pooled)"
        if r["category"] == "UNCATEGORIZED":
            return 4, "uncategorized ((activity_type, unit) not in DR/SP)"
        if pd.isna(r["direction"]) or int(r["direction"]) not in (-1, 1):
            return 5, f"undefined direction ({r['direction']})"
        if r["source_label"] == COADD_SOURCE_LABEL:
            return 6, "Co-ADD"
        if r["assay_chembl_id"] in pubchem_ids:
            return 7, "superseded by PubChem"
        if r["n_mol_dataset"] <= MIN_MOLECULES:
            return 8, f"<= {MIN_MOLECULES} molecules"
        return 0, ""

    decided = pairs.apply(decide, axis=1, result_type="expand")
    pairs["discard_step"] = decided[0]
    pairs["discard_reason"] = decided[1]
    pairs["discarded"] = pairs["discard_step"] != 0

    retained = pairs[~pairs["discarded"]].copy()
    print(f"  pairs total: {len(pairs)} | retained: {len(retained)} "
          f"| discarded: {int(pairs['discarded'].sum())}")
    print("  discard reasons:", dict(Counter(pairs.loc[pairs['discarded'], 'discard_reason'])))
    print(f"  flags -> multi_unit: {int(pairs['flag_multi_unit'].sum())} "
          f"| uncategorized (retained-eligible): {int(pairs['flag_uncategorized'].sum())}")
    print(f"  qual-inactive recovered (retained pairs): {int(retained['n_recovered_inactive'].sum())} "
          f"| qual-active text-only discarded: {int(pairs['n_text_only_active_discarded'].sum())}")

    # --- write retained datasets --------------------------------------------
    keep_cols = ["inchikey", "compound_chembl_id", "smiles", "value", "relation",
                 "unit", "activity_type", "direction"]
    retained_keys = set(map(tuple, retained[KEY].values))
    dir_lookup = {tuple(k): d for k, d in zip(retained[KEY].values, retained["direction"].values)}
    cat_lookup = {tuple(k): c for k, c in zip(retained[KEY].values, retained["category"].values)}
    unit_lookup = {tuple(k): u for k, u in zip(retained[KEY].values, retained["unit"].values)}

    written = 0
    for key, g in tqdm(num.groupby(KEY), desc="  writing datasets"):
        if key not in retained_keys:
            continue
        direction = dir_lookup[key]
        ded = dedup_most_active_by_inchikey(g[keep_cols], direction, preserve_qualifier=True)
        ded["category"] = cat_lookup[key]
        assay_chembl_id, activity_type = key
        # append recovered qualitative-inactive molecules (no numeric value); they
        # are labelled inactive downstream by binarize_conservative (step 23).
        rr = recov_by_key.get(key)
        if rr is not None and len(rr):
            add = pd.DataFrame({
                "inchikey": rr["inchikey"].values,
                "compound_chembl_id": rr["compound_chembl_id"].values,
                "smiles": rr["smiles"].values,
                "value": np.nan,
                "relation": np.nan,
                "unit": unit_lookup[key],
                "activity_type": activity_type,
                "direction": direction,
                "category": cat_lookup[key],
            })
            ded = pd.concat([ded, add], ignore_index=True)
        fname = f"{assay_chembl_id}_{str(activity_type).replace('/', '-')}.csv.gz"
        ded.to_csv(os.path.join(DATASETS_DIR, fname), index=False)
        written += 1
    print(f"  datasets written: {written}")
    assert written == len(retained), \
        f"wrote {written} dataset files but {len(retained)} pairs were retained"

    # --- write summary ------------------------------------------------------
    pairs = pairs.rename(columns={"assay_chembl_id": "assay_id"})
    pairs["numeric_units"] = pairs["numeric_units"].apply(lambda x: ";".join(x))
    pairs["numeric_directions"] = pairs["numeric_directions"].apply(
        lambda x: ";".join(str(v) for v in x))

    summary_cols = (
        ["assay_id", "activity_type", "unit", "category", "direction"]
        + ENRICH_COLS
        + ["n_activities", "n_numeric", "n_text", "n_mol_all", "n_mol_numeric",
           "n_recovered_inactive", "n_text_only_active_discarded", "n_mol_dataset",
           "n_numeric_units", "numeric_units", "numeric_directions",
           "flag_multi_unit", "flag_uncategorized",
           "discarded", "discard_step", "discard_reason"]
    )
    pairs = pairs[summary_cols].sort_values(
        ["discarded", "n_mol_numeric"], ascending=[True, False]).reset_index(drop=True)
    pairs.to_csv(os.path.join(OUTPUT, "21_curation_summary.csv"), index=False)
    print(f"  summary -> {os.path.join(OUTPUT, '21_curation_summary.csv')}")

    # reconciliation guard: every pair is either retained or has exactly one reason
    assert (pairs["discarded"] == (pairs["discard_reason"] != "")).all(), \
        "discarded flag and reason are inconsistent"
    print("  OK - reconciliation passed")


# --- individual-plot panels (each draws into a provided ax) ------------------
# Colours come from the shared semantic palette (plot_colors); article style is
# applied by plot_base at import.

class AttritionBarPlot(BasePlot):
    """Per-stage attrition bar for one metric (datasets / chemical space / ...)."""

    def __init__(self, ax, labels, values, color_key, ylabel, title, abc=None, pct_of=None):
        super().__init__(ax, abc)
        x = np.arange(len(labels))
        ax.bar(x, values, color=pc.color(color_key))
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        if pct_of:
            ax.set_ylim(0, max(values) * 1.18)  # headroom so % labels stay inside
            for xi, v in zip(x, values):
                ax.text(xi, v, f"{100 * v / pct_of:.0f}%", ha="center", va="bottom")
        self._label(ylabel=ylabel, title=title)
        self.data = pd.DataFrame({"stage": labels, "value": values})


class HistPlot(BasePlot):
    """Distribution of a per-dataset / per-molecule count (optionally log-x)."""

    def __init__(self, ax, data, color_key, xlabel, title, abc=None, log_x=False, bins=30):
        super().__init__(ax, abc)
        data = np.asarray([d for d in data if d is not None and not pd.isna(d)], dtype=float)
        if len(data) == 0:
            self._empty(title=title, xlabel=xlabel, ylabel="count")
            return
        if log_x:
            data = data[data > 0]
            bins = np.logspace(np.log10(data.min()), np.log10(max(data.max(), 2)), bins)
            ax.set_xscale("log")
        ax.hist(data, bins=bins, color=pc.color(color_key))
        self._label(xlabel=xlabel, ylabel="count", title=title)
        self.data = pd.DataFrame({"value": data})


class CategoryBarPlot(BasePlot):
    """Retained-dataset counts by activity category (DR / SP / UNCATEGORIZED)."""

    def __init__(self, ax, counts, title, abc=None):
        super().__init__(ax, abc)
        keys = [c for c in CATEGORY_ORDER if c in counts] + \
               [c for c in counts if c not in CATEGORY_ORDER]
        if not keys:
            self._empty(title=title, ylabel="datasets")
            return
        values = [counts[k] for k in keys]
        # keep DR/SP on their semantic colours (cobalt/tangerine) so they match every
        # other figure; fall back to the categorical palette for anything else.
        fallback = pc.category_palette(len(keys))
        colors = [pc.color(k) if k in ("DR", "SP") else
                  (pc.color("neutral") if k == "UNCATEGORIZED" else fallback[i])
                  for i, k in enumerate(keys)]
        x = np.arange(len(keys))
        ax.bar(x, values, color=colors)
        ax.set_xticks(x)
        ax.set_xticklabels([cat_label(k) for k in keys], rotation=45, ha="right")
        self._label(ylabel="datasets", title=title)
        self.data = pd.DataFrame({"category": [cat_label(k) for k in keys], "datasets": values})


class BinnedBarPlot(BasePlot):
    """Counts of a skewed per-item integer value grouped into a few explicit bins —
    clearer than a log-x histogram dominated by the '1' bar. `bins` is a list of
    (label, lo, hi) with hi=None meaning open-ended (value >= lo)."""

    def __init__(self, ax, data, bins, color_key, xlabel, title, abc=None):
        super().__init__(ax, abc)
        vals = np.asarray([v for v in data if v is not None and not pd.isna(v)], dtype=float)
        counts = [int(((vals >= lo) & (vals <= hi)).sum()) if hi is not None
                  else int((vals >= lo).sum()) for _, lo, hi in bins]
        x = np.arange(len(bins))
        ax.bar(x, counts, color=pc.color(color_key))
        ax.set_xticks(x)
        ax.set_xticklabels([b[0] for b in bins])
        self._label(xlabel=xlabel, ylabel="count", title=title)
        self.data = pd.DataFrame({"bin": [b[0] for b in bins], "count": counts})


class GroupedSharePlot(BasePlot):
    """DR-vs-SP split of the retained data across three measures (datasets, molecules,
    chemical space), each stacked to 100% within its measure. Molecules use the same
    'sum of per-dataset uniques' definition as the Total-molecules bar; chemical space
    is unique InChIKeys per category (DR/SP can overlap via shared molecules, so that
    group's split is of the summed per-category counts, not the retained union)."""

    def __init__(self, ax, datasets, mols, chem, abc=None):
        super().__init__(ax, abc)
        measures = [("Datasets", datasets), ("Molecules", mols), ("Chem space", chem)]
        x = np.arange(len(measures))
        dr_ct = [d.get("DR", 0) for _, d in measures]
        sp_ct = [d.get("SP", 0) for _, d in measures]
        tot = [(a + b) or 1 for a, b in zip(dr_ct, sp_ct)]
        dr_frac = [100 * a / t for a, t in zip(dr_ct, tot)]
        sp_frac = [100 * b / t for b, t in zip(sp_ct, tot)]
        ax.bar(x, dr_frac, color=pc.color("DR"), label="Dose Response")
        ax.bar(x, sp_frac, bottom=dr_frac, color=pc.color("SP"), label="Single Point")
        ax.set_ylim(0, 100)
        ax.set_xticks(x)
        ax.set_xticklabels([m[0] for m in measures], rotation=45, ha="right")
        ax.legend(loc="lower right")  # DR/SP key inside the panel, over the DR area
        self._label(ylabel="% of retained (DR+SP)", title="Retained DR/SP split")
        self.data = pd.DataFrame({"measure": [m[0] for m in measures],
                                  "DR": dr_ct, "SP": sp_ct})


class WaterfallPlot(BasePlot):
    """Attrition as a descending staircase: each nested stage's survivors as a
    grounded bar, with the amount lost since the previous stage drawn as a faint
    cap on top reaching the prior level. % labels are relative to stage 0."""

    def __init__(self, ax, labels, values, ylabel, title, abc=None):
        super().__init__(ax, abc)
        values = [float(v) for v in values]
        x = np.arange(len(values))
        surv = pc.color("survivor")
        lost = pc.color("lost")
        for i, v in enumerate(values):
            ax.bar(x[i], v, color=surv)
            if i > 0:
                drop = values[i - 1] - v
                if drop > 0:
                    ax.bar(x[i], drop, bottom=v, color=lost)
        base = values[0] or 1.0
        ax.set_ylim(0, (max(values) or 1.0) * 1.15)
        for i, v in enumerate(values):
            ax.text(x[i], v, f"{100 * v / base:.0f}%", ha="center", va="bottom")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        self._label(ylabel=ylabel, title=title)
        self.data = pd.DataFrame({"stage": labels, "value": values})


# --- per-pathogen processing ------------------------------------------------

def process(pathogen_code):
    pathogen = load_pathogen(pathogen_code)
    out = os.path.join(STAGE4PATH, pathogen_code)
    print(f"Step 21 (plot) - {pathogen_code} ({pathogen})")

    stats = curation_stats(pathogen_code)

    # per-stage attrition series computed once in curation_stats (single source of
    # truth); consumed here for the per-pathogen figures / stats CSV.
    stage_labels = stats["stage_labels"]
    stage_assays = stats["stage_assays"]
    stage_cs = stats["stage_cs"]
    stage_mols = stats["stage_mols"]
    stage_strains = stats["stage_strains"]
    total_cs = stats["total_cs"]

    sizes = stats["sizes"]
    reuse = stats["reuse"]
    # per-category (DR/SP) retained metrics for the H panel (datasets, molecules, chem space)
    cat_datasets = stats["category_datasets"]
    cat_mols = stats["category_mols"]
    cat_chem = stats["category_chem"]

    # datasets per document (per-pathogen figure only; not carried in the results dict)
    summ = pd.read_csv(os.path.join(out, "21_curation_summary.csv"), low_memory=False)
    retained = summ[~summ["discarded"]].copy()
    docs_per_assay = retained.groupby("doc_chembl_id").size().tolist()

    # explicit bins for the two skewed count distributions (F, G): compact bars beat a
    # log-x histogram dominated by the '1' bar.
    REUSE_BINS = [("1", 1, 1), ("2", 2, 2), ("3–5", 3, 5), (">5", 6, None)]
    DOC_BINS = [("1", 1, 1), ("2–3", 2, 3), ("4–10", 4, 10), (">10", 11, None)]

    # overview attrition bars use a reduced stage set (drop the near-flat Single-unit and
    # DR/SP stages); the waterfall below keeps the full cascade.
    ov = [i for i, l in enumerate(stage_labels) if l not in OVERVIEW_DROP]
    ov_labels = [stage_labels[i] for i in ov]
    ov_assays = [stage_assays[i] for i in ov]
    ov_cs = [stage_cs[i] for i in ov]
    ov_mols = [stage_mols[i] for i in ov]
    ov_strains = [stage_strains[i] for i in ov]

    # overview figure (2x4): 4 per-stage attrition bars + 3 distributions + the DR/SP
    # split panel (H). H carries its own bottom-right DR/SP legend.
    render_figure(out, "21_curation_overview", 2, 4, [
        lambda ax: AttritionBarPlot(ax, ov_labels, ov_assays, "datasets",
                                    "datasets", "Datasets per stage", abc="A"),
        lambda ax: AttritionBarPlot(ax, ov_labels, ov_cs, "chem_space",
                                    "unique InChIKeys", "Chemical space", abc="B",
                                    pct_of=total_cs or None),
        lambda ax: AttritionBarPlot(ax, ov_labels, ov_mols, "molecules",
                                    "molecules (sum)", "Total molecules", abc="C"),
        lambda ax: AttritionBarPlot(ax, ov_labels, ov_strains, "strains",
                                    "unique strains", "Strains (normalised)", abc="D"),
        lambda ax: HistPlot(ax, sizes, "sizes",
                            "molecules / dataset", "Retained dataset sizes", abc="E", log_x=True),
        lambda ax: BinnedBarPlot(ax, reuse, REUSE_BINS, "reuse",
                                 "datasets / molecule", "Molecule reuse", abc="F"),
        lambda ax: BinnedBarPlot(ax, docs_per_assay, DOC_BINS, "docs",
                                 "datasets / document", "Datasets per document", abc="G"),
        lambda ax: GroupedSharePlot(ax, cat_datasets, cat_mols, cat_chem, abc="H"),
    ])

    # attrition waterfall (datasets + chemical space), nested-stage staircase, with a
    # single shared legend placed below the panels (never over the bars).
    def _wf_legend(fig, panels):
        figure_legend(fig, [("retained", pc.color("survivor")),
                            ("lost at stage", pc.color("lost"))],
                      loc="lower center", anchor=(0.5, -0.02), ncol=2)
    render_figure(out, "21_curation_attrition_waterfall", 1, 2, [
        lambda ax: WaterfallPlot(ax, stage_labels, stage_assays,
                                 "Datasets", "Dataset attrition", abc="A"),
        lambda ax: WaterfallPlot(ax, stage_labels, stage_cs,
                                 "Unique InChIKeys", "Chemical-space attrition", abc="B"),
    ], legend=_wf_legend)

    # companion "ready dataset" for the two distributions not carried in any existing
    # CSV (molecule reuse, datasets per document); the attrition bars/waterfall are
    # already in 21_curation_stats.csv and retained sizes in 21_curation_summary.
    dist_rows = ([{"distribution": "molecule_reuse", "value": v} for v in reuse]
                 + [{"distribution": "datasets_per_document", "value": v} for v in docs_per_assay])
    pd.DataFrame(dist_rows, columns=["distribution", "value"]).to_csv(
        os.path.join(out, "21_curation_distributions.csv"), index=False)

    # stats CSV
    org_cs = stage_cs[1] if len(stage_cs) > 1 else 0
    stats_df = pd.DataFrame({"stage": stage_labels, "datasets": stage_assays,
                             "chemical_space": stage_cs, "total_molecules": stage_mols,
                             "strains": stage_strains})
    stats_df["pct_total_cs"] = (100 * stats_df["chemical_space"] / total_cs).round(2) if total_cs else 0
    stats_df["pct_org_cs"] = (100 * stats_df["chemical_space"] / org_cs).round(2) if org_cs else 0
    stats_df.to_csv(os.path.join(out, "21_curation_stats.csv"), index=False)

    # companion source data for the H panel (retained DR/SP split by measure). Molecules =
    # sum of per-dataset uniques (same definition as the "Total molecules" bar); chem_space
    # = unique InChIKeys per category (DR/SP overlap on shared molecules).
    pd.DataFrame([{"category": c,
                   "datasets": cat_datasets.get(c, 0),
                   "molecules": cat_mols.get(c, 0),
                   "chem_space": cat_chem.get(c, 0)} for c in ["DR", "SP"]]).to_csv(
        os.path.join(out, "21_curation_categories.csv"), index=False)

    return stats


if __name__ == "__main__":
    pathogen_code = sys.argv[1]
    # Publication style (also applied at plot_base import; re-assert here so the import
    # order relative to plot_utils, which sets slide+ersilia, can never win).
    stylia.set_format("print")
    stylia.set_style("article")
    build(pathogen_code)
    process(pathogen_code)
