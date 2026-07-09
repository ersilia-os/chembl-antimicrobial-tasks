"""
Step 27 - Cross-pathogen comparison figures (ALL Stage-4 general/cross-pathogen plotting).

The single terminal general-plots step: every cross-pathogen figure lives here, drawn from
the per-pathogen summary files the earlier steps (21-26) wrote. Nothing here re-runs a model;
it only aggregates on-disk results. Per-pathogen figures stay with their own steps; step 28
(trace/flow audit) is a separate report.

Consolidates three previously-scattered sources: the curation/binarisation comparisons (this
step's originals), the pool-summary figures (moved from the transfer step), and the low-data
catch-all overview (moved from the low-data step).

Aggregators curation_stats/binarisation_stats live in src/plot_utils.py and return the same
results dicts steps 21/22 build, so figures never diverge from the computed statistics.

Outputs (output/stage4/general_plots/):
    27_curation_frac_chembl.png       ChEMBL chemical space covered by the pathogens
    27_curation_overlap.png           pairwise molecule overlap (cleaned + retained)
    27_curation_target_type.png       target-type composition (assays + unique molecules)
    27_curation_discard.png           discard cascade (assays + molecule instances)
    27_curation_retained.png          retained-data characterization (chem space, sizes,
                                      category, reuse)
    27_curation_wholecell_sizes.png   Whole-cell dataset-size distribution (pre-selection)
    27_chembl_space.json              cached full-ChEMBL bioactive compound count
    27_binarisation_ratio.png         active-ratio box (DR/SP) + class composition
    27_binarisation_agreement.png     hit reproducibility / contradictions / DR↔SP agreement
    27_binarisation_active_overlap.png  active-molecule Jaccard heatmaps (DR, SP)
    27_cutoff_sensitivity.png         active ratio at each low/middle/high candidate cutoff,
                                      per pathogen and unit (DR µM / SP % / SP mm); current
                                      global default starred. Diagnostic for per-pathogen
                                      cutoff overrides - changes no label.
    27_cutoff_sensitivity.csv         tidy (pathogen, category, unit, tier, cutoff, counts,
                                      active_ratio, is_default) table behind that figure
    27_eos1klk_all_pathogens.png      eos1klk embedding of retained molecules, coloured by pathogen
    27_eos1klk_tmap.png               t-MAP small multiples (retained molecules), one panel per pathogen
    27_pool_partition.png             step-23 molecule-coverage split (multi/singleton/not-considered)
    27_pool_active_ratios.png         step-23 active ratio of considered pools (multi + singleton) across pathogens
    27_pool_overview.png              condensed per-pathogen pools (3x2, DR/SP): step-23 pool size,
                                      step-24 CV AUROC, step-24 CV outcome counts
    27_merge_overview.png             step-25 merge per pathogen (3x2, DR/SP): baseline→grown AUROC,
                                      leftover outcome mix, molecules absorbed
    27_final_data_overlap.png         final-pool overlap per pathogen (2x2, DR/SP): pairwise Jaccard
                                      overlap % + molecule-level cross-pool label conflict
                                      (conflicting / total unique molecules)
    27_final_data_overlap.csv         tidy per-pool-pair table (n_shared, jaccard, discrepancy)
    27_activity_ratio.png             active ratio through the pipeline (2x2): stage distribution
                                      22→26 (DR/SP) + per-pathogen seed→final pool ratio
    27_activity_ratio_per_pathogen.png  small multiples: 22→26 active-ratio flow per pathogen (DR/SP)
    27_lowdata_overview.png           catch-all AUROC per deferred category (from step 26)
    27_master_table.csv               pooling backbone: one row per binarised dataset traced
                                      fold→pool(23)→CV(24)→merge(25)/catch-all(26) with final status
    27_pipeline_summary.png           cross-step summary (2x3): consolidation funnel, coverage,
                                      final-pool AUROC + active ratio, per-pathogen fate + pool counts

Usage:
    python scripts/27_general_plots.py                          # all processed pathogens
    python scripts/27_general_plots.py saureus mtuberculosis    # an explicit subset
"""

import os
import sys
import json
import glob
import itertools
from collections import defaultdict

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from plot_utils import (
    curation_stats, binarisation_stats,
    load_coords, embed_figure,
    pathogen_colors,
    cat_label, size_bins,
    RANDOM_SEED, CHEMBL_ACTIVITIES_05, CAP_POINTS,
    TARGET_TYPES, TARGET_TYPE_LABELS, TT_MOL_REGIONS,
    STEP_LABELS, STEP_ORDER,
    SIZE_LABELS,
    CAT_DISPLAY, COMPOSITION, HIT_LABELS,
)
import plot_colors as pc
from plot_base import figure_legend, inset_legend, compact_yaxis
from default import STAGE4PATH, STAGE2PATH, DATABASE_NAME, CONFIGPATH
from pathogen_utils import load_pathogen
from dataset_utils import (
    UNIT_SLUG, CUTOFF_TIERS,
    load_category_cutoffs, load_category_cutoff_candidates,
)


# ============================================================================
# CURATION cross-pathogen helpers (per-pathogen figures are in step 21)
# ============================================================================

def plot_overlap_heatmap(ax, codes, sets, title="Shared retained molecules", abc=None):
    """Pairwise Jaccard of molecule sets (cleaned, retained, or active); diagonal
    masked. The colour scale is capped at the 98th percentile of off-diagonal values
    so the (mostly small) shared fractions stay legible."""
    n = len(codes)
    M = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a, b = sets[i], sets[j]
            u = len(a | b)
            M[i, j] = (len(a & b) / u) if u else 0.0
    off = M[~np.isnan(M)]
    vmax = max(float(np.nanpercentile(off, 98)) if off.size else 1.0, 1e-3)
    cm = pc.fading("total").cmap.copy()
    cm.set_bad("#eeeeee")
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cm, vmin=0, vmax=vmax, aspect="equal")
    ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xticks(np.arange(n)); ax.set_xticklabels(codes, rotation=90)
    ax.set_yticks(np.arange(n)); ax.set_yticklabels(codes)
    stylia.label(ax, xlabel="", ylabel="",
                 title=f"{title} (Jaccard; scale ≤ {vmax:.2f})", abc=abc)


# --- story figure: chemical-space overview (Panels A-D) ---------------------

def chembl_bioactive_count(cache_path):
    """Unique compounds with any bioactivity record (the Panel A denominator).
    Computed once from 05_activities_preprocessed.csv (~20 s over 24M rows) and
    cached as JSON stamped with the ChEMBL DB version (`DATABASE_NAME`, e.g.
    'chembl_36'). The cache is reused only while that stamp matches the current
    DATABASE_NAME, so a Stage-1 refresh against a new ChEMBL release forces a
    recompute rather than silently reusing a stale count. Returns None if 05 is
    absent. (No explicit ChEMBL version string exists in the pipeline; the DB name
    is the only version signal available.)"""
    if os.path.exists(cache_path):
        with open(cache_path) as fh:
            cached = json.load(fh)
        if cached.get("chembl_db") == DATABASE_NAME:
            return cached["bioactive_compounds"]
        print(f"  [story] cache stamped {cached.get('chembl_db')!r} != current "
              f"{DATABASE_NAME!r}; recomputing")
    if not os.path.exists(CHEMBL_ACTIVITIES_05):
        print(f"  [story] {CHEMBL_ACTIVITIES_05} not found - Panel A will be blank")
        return None
    print("  [story] computing full-ChEMBL bioactive chemical space from 05 ...")
    ids = set()
    for chunk in pd.read_csv(CHEMBL_ACTIVITIES_05, usecols=["compound_chembl_id"],
                             dtype=str, chunksize=2_000_000):
        ids.update(chunk["compound_chembl_id"].dropna().unique())
    n = len(ids)
    with open(cache_path, "w") as fh:
        json.dump({"bioactive_compounds": n, "chembl_db": DATABASE_NAME,
                   "source": "output/stage1/05_activities_preprocessed.csv"}, fh)
    print(f"  [story] bioactive ChEMBL compounds: {n:,} ({DATABASE_NAME}, cached -> {cache_path})")
    return n


def plot_chembl_coverage(ax, related, total, n_pathogens, abc):
    """Donut: fraction of the bioactive ChEMBL chemical space covered by the union
    of all pathogens' cleaned molecules."""
    if not total:
        stylia.label(ax, xlabel="", ylabel="", title="ChEMBL coverage (n/a)", abc=abc)
        return
    rest = max(total - related, 0)
    ax.pie([related, rest], colors=[pc.color("related"), pc.lighten("neutral", 0.6)],
           startangle=90, counterclock=False, wedgeprops={"width": 0.42})
    ax.text(0, 0, f"{100 * related / total:.1f}%", ha="center", va="center")
    stylia.label(ax, xlabel="", ylabel="",
                 title=f"Chem-space related to the {n_pathogens} pathogens\n"
                       f"({related:,} of {total:,} bioactive ChEMBL compounds)", abc=abc)


def molecule_region_colors():
    """Region -> colour (from the shared semantic palette): Whole-cell turquoise,
    shared orchid (between the two pure types), Protein target amber, Other silver."""
    return dict(zip(TT_MOL_REGIONS, pc.region_colors()))


def plot_targettype_molecules(ax, codes, results, abc, normalize=True, show_y=True):
    """Horizontal stacked bars (pathogens on y): cleaned chemical space split into
    Whole-cell, the Whole-cell+Protein-target overlap, Protein target, and Other.
    normalize=True -> fractions summing to 1; False -> absolute molecule counts."""
    cmap = molecule_region_colors()
    y = np.arange(len(codes))
    lefts = np.zeros(len(codes))
    for lab in TT_MOL_REGIONS:
        raw = np.array([r["tt_combo_counts"].get(lab, 0) for r in results], dtype=float)
        denom = np.array([max(sum(r["tt_combo_counts"].values()), 1) for r in results])
        vals = raw / denom if normalize else raw
        ax.barh(y, vals, left=lefts, color=cmap[lab], label=lab)
        lefts += vals
    ax.set_yticks(y); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    xlabel = "Fraction of chem-space" if normalize else "Unique molecules"
    title = "Chem-space by target type" if normalize else "Total unique molecules"
    stylia.label(ax, xlabel=xlabel, ylabel="", title=title, abc=abc)


def plot_targettype_assays(ax, codes, results, colors, abc, normalize=True, show_y=True):
    """Horizontal stacked bars (pathogens on y): assays per target type (disjoint).
    normalize=True -> fractions summing to 1; False -> absolute assay counts."""
    y = np.arange(len(codes))
    lefts = np.zeros(len(codes))
    for i, T in enumerate(TARGET_TYPES):
        raw = np.array([r["tt_assay_counts"].get(T, 0) for r in results], dtype=float)
        denom = np.array([max(sum(r["tt_assay_counts"].values()), 1) for r in results])
        vals = raw / denom if normalize else raw
        ax.barh(y, vals, left=lefts, color=colors[i], label=TARGET_TYPE_LABELS[T])
        lefts += vals
    ax.set_yticks(y); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    xlabel = "Fraction of assays" if normalize else "Number of assays"
    title = "Assays by target type" if normalize else "Total assays"
    stylia.label(ax, xlabel=xlabel, ylabel="", title=title, abc=abc)


# --- attrition figure helpers (horizontal, pathogens on y) ------------------

def plot_discard(ax, codes, results, field, colors, abc, normalize, entity, show_y=True):
    """Horizontal stacked bars of the curation outcome per pathogen, including the
    'retained' segment. field = 'step_counts' (by assay) or 'step_mols' (by
    molecule-occurrence). normalize -> fractions; else absolute counts."""
    y = np.arange(len(codes))
    lefts = np.zeros(len(codes))
    for color, st in zip(colors, STEP_ORDER):
        raw = np.array([r[field].get(st, 0) for r in results], dtype=float)
        denom = np.array([max(sum(r[field].values()), 1) for r in results])
        vals = raw / denom if normalize else raw
        ax.barh(y, vals, left=lefts, color=color, label=STEP_LABELS[st])
        lefts += vals
    ax.set_yticks(y); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    pre = "Fraction of" if normalize else "Number of"
    # entity (assays / molecule instances) is already on the x-axis, so keep the title short
    stylia.label(ax, xlabel=f"{pre} {entity}", ylabel="", title="Discard reason", abc=abc)


def plot_box_h(ax, codes, data_lists, color, xlabel, title, abc, rng,
               logx=False, show_y=True):
    """Horizontal empty boxplots with jittered points (subsampled to CAP_POINTS)."""
    box_color = pc.color("ink")
    pos = np.arange(len(codes))
    ax.boxplot(data_lists, positions=pos, vert=False, showfliers=False, widths=0.6,
               patch_artist=False, boxprops={"color": box_color},
               whiskerprops={"color": box_color}, capprops={"color": box_color},
               medianprops={"color": box_color})
    for i, data in enumerate(data_lists):
        d = np.asarray([v for v in data if v is not None and not pd.isna(v)], dtype=float)
        if logx:
            d = d[d > 0]
        if len(d) == 0:
            continue
        if len(d) > CAP_POINTS:
            # random sample for the bulk, but always keep the extreme values so outlier
            # datasets (e.g. pfalciparum / mtb large sets) are never sampled away
            srt = np.sort(d)
            d = np.concatenate([rng.choice(d, size=CAP_POINTS, replace=False),
                                srt[:10], srt[-10:]])
        yj = pos[i] + rng.uniform(-0.18, 0.18, size=len(d))
        ax.scatter(d, yj, s=6, alpha=0.45, color=color, edgecolors="none", zorder=1)
    if logx:
        ax.set_xscale("log")
    ax.set_yticks(pos); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel=xlabel, ylabel="", title=title, abc=abc)


def plot_pct_retained(ax, codes, results, abc, show_y=True):
    """Grouped horizontal bars: % of chemical space retained, measured against the
    pathogen's total cleaned space and against its Whole-cell (ORGANISM) space."""
    y = np.arange(len(codes))
    h = 0.4
    ax.barh(y - h / 2, [r["pct_total_cs"] for r in results], height=h,
            color=pc.color("total"), label="of total")
    ax.barh(y + h / 2, [r["pct_org_cs"] for r in results], height=h,
            color=pc.color("wholecell"), label="of Whole-cell")
    ax.set_yticks(y); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel="Retained (%)", ylabel="",
                 title="Chemical space retained", abc=abc)


def plot_reuse(ax, codes, results, abc, show_y=True):
    """Horizontal stacked bars: fraction of retained molecules appearing in 1, 2, or
    >=3 retained datasets. A larger reused (2 / >=3) share = more redundancy across
    datasets (more poolable); a large '1 dataset' share = mostly dataset-private."""
    rc = pc.reuse_colors()
    bins = [("1 dataset", lambda v: v == 1, rc[0]),
            ("2 datasets", lambda v: v == 2, rc[1]),
            ("≥3 datasets", lambda v: v >= 3, rc[2])]
    y = np.arange(len(codes))
    lefts = np.zeros(len(codes))
    for lab, pred, color in bins:
        fracs = np.array([sum(1 for v in r["reuse"] if pred(v)) / (len(r["reuse"]) or 1)
                          for r in results])
        ax.barh(y, fracs, left=lefts, color=color, label=lab)
        lefts += fracs
    ax.set_yticks(y); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel="Fraction of unique molecules", ylabel="",
                 title="Molecule reuse", abc=abc)


def plot_category_h(ax, codes, results, abc, show_y=True):
    """Horizontal stacked bars: DR vs SP split of retained datasets. (UNCATEGORIZED is
    discarded in step 21, so retained datasets are only DR/SP.)"""
    cats = ["DR", "SP"]
    colors = {"DR": pc.color("DR"), "SP": pc.color("SP")}
    y = np.arange(len(codes))
    lefts = np.zeros(len(codes))
    for cat in cats:
        raw = np.array([r["category_counts"].get(cat, 0) for r in results], dtype=float)
        denom = np.array([max(sum(r["category_counts"].values()), 1) for r in results])
        ax.barh(y, raw / denom, left=lefts, color=colors[cat], label=cat_label(cat))
        lefts += raw / denom
    ax.set_yticks(y); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel="Fraction of datasets", ylabel="",
                 title="Dataset category", abc=abc)


# --- Whole-cell dataset-size figure helpers ---------------------------------

def plot_size_box(ax, codes, size_lists, bin_colors, abc, rng, show_y=True):
    """Horizontal boxplot of Whole-cell dataset sizes (log x) with jittered scatter
    points coloured by their size bin."""
    box_color = pc.color("ink")
    pos = np.arange(len(codes))
    ax.boxplot(size_lists, positions=pos, vert=False, showfliers=False, widths=0.6,
               patch_artist=False, boxprops={"color": box_color},
               whiskerprops={"color": box_color}, capprops={"color": box_color},
               medianprops={"color": box_color})
    for i, sizes in enumerate(size_lists):
        d = np.asarray([v for v in sizes if v and not pd.isna(v)], dtype=float)
        d = d[d > 0]
        if len(d) == 0:
            continue
        if len(d) > CAP_POINTS:
            d = rng.choice(d, size=CAP_POINTS, replace=False)
        cols = [bin_colors[k] for k in size_bins(d)]
        yj = pos[i] + rng.uniform(-0.18, 0.18, size=len(d))
        ax.scatter(d, yj, s=4, alpha=0.5, color=cols, edgecolors="none", zorder=1)
    ax.set_xscale("log")
    ax.set_yticks(pos); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel="Molecules / dataset", ylabel="",
                 title="Whole-cell dataset sizes", abc=abc)


def plot_size_stack(ax, codes, size_lists, bin_colors, abc, weight, normalize, show_y=True):
    """Horizontal stacked bars by size bin. weight='assays' counts datasets per bin;
    weight='molecules' sums dataset sizes per bin. normalize -> fractions summing to 1."""
    nb = len(SIZE_LABELS)
    mat = np.zeros((len(codes), nb))
    for i, sizes in enumerate(size_lists):
        d = np.asarray([v for v in sizes if v and not pd.isna(v)], dtype=float)
        d = d[d > 0]
        idx = size_bins(d)
        for k in range(nb):
            sel = d[idx == k]
            mat[i, k] = len(sel) if weight == "assays" else sel.sum()
    if normalize:
        tot = mat.sum(axis=1, keepdims=True)
        tot[tot == 0] = 1
        mat = mat / tot
    y = np.arange(len(codes))
    lefts = np.zeros(len(codes))
    for k, (lab, color) in enumerate(zip(SIZE_LABELS, bin_colors)):
        ax.barh(y, mat[:, k], left=lefts, color=color, label=lab)
        lefts += mat[:, k]
    ax.set_yticks(y); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    if weight == "assays":
        xlabel, title = "Fraction of assays", "Assays by dataset size"
    elif normalize:
        xlabel, title = "Fraction of molecule instances (non-unique)", "Molecule instances by dataset size"
    else:
        xlabel, title = "Molecule instances (non-unique)", "Molecule instances by dataset size"
    stylia.label(ax, xlabel=xlabel, ylabel="", title=title, abc=abc)


def make_curation_comparison(results, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    results = sorted(results, key=lambda r: r["n_retained"], reverse=True)
    codes = [r["code"] for r in results]
    rng = np.random.default_rng(RANDOM_SEED)

    # 1. ChEMBL coverage donut (square, standalone).
    sets08 = [r["all_iks"] for r in results]
    retsets = [r["retained_iks"] for r in results]
    related = len(set().union(*sets08)) if sets08 else 0
    total_chembl = chembl_bioactive_count(os.path.join(out_dir, "27_chembl_space.json"))
    fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
    plot_chembl_coverage(axs.next(), related, total_chembl, len(results), abc=None)
    stylia.save_figure(os.path.join(out_dir, "27_curation_frac_chembl.png"))

    # 2. Chemical-space overlap heatmaps (square): cleaned vs retained molecules.
    fig, axs = stylia.create_figure(1, 2, width=1.0, height=0.55)
    plot_overlap_heatmap(axs.next(), codes, sets08, title="Overlap (cleaned)", abc="A")
    plot_overlap_heatmap(axs.next(), codes, retsets, title="Overlap (retained)", abc="B")
    stylia.save_figure(os.path.join(out_dir, "27_curation_overlap.png"))

    # 3. ATTRITION — split into three compact thematic figures (each 1x4, pathogens on
    #    y, labels only on the leftmost panel, single legend row directly below). No
    #    empty legend-strip rows, so panels stay tight.
    tt_colors = pc.targettype_colors()                   # Whole-cell, Protein target, Other
    step_colors = pc.discard_step_colors()               # STEP_ORDER order (9)
    region_cmap = molecule_region_colors()               # 4 molecule regions

    # (a) target-type composition: fraction + count, by assay and by unique molecules
    fig, axs = stylia.create_figure(1, 4, height=0.6)
    plot_targettype_assays(axs.next(), codes, results, tt_colors, "A", normalize=True, show_y=True)
    plot_targettype_molecules(axs.next(), codes, results, "B", normalize=True, show_y=False)
    plot_targettype_assays(axs.next(), codes, results, tt_colors, "C", normalize=False, show_y=False)
    plot_targettype_molecules(axs.next(), codes, results, "D", normalize=False, show_y=False)
    figure_legend(fig, [(lab, region_cmap[lab]) for lab in TT_MOL_REGIONS],
                  loc="upper center", anchor=(0.5, -0.02), ncol=4)
    stylia.save_figure(os.path.join(out_dir, "27_curation_target_type.png"))

    # (b) discard cascade: fraction + count, by assay and by molecule instances
    fig, axs = stylia.create_figure(1, 4, height=0.6)
    plot_discard(axs.next(), codes, results, "step_counts", step_colors, "A", True, "assays", show_y=True)
    plot_discard(axs.next(), codes, results, "step_mols", step_colors, "B", True, "molecule instances (non-unique)", show_y=False)
    plot_discard(axs.next(), codes, results, "step_counts", step_colors, "C", False, "assays", show_y=False)
    plot_discard(axs.next(), codes, results, "step_mols", step_colors, "D", False, "molecule instances (non-unique)", show_y=False)
    figure_legend(fig, [(STEP_LABELS[st], step_colors[i]) for i, st in enumerate(STEP_ORDER)],
                  loc="upper center", anchor=(0.5, -0.02), ncol=5)
    stylia.save_figure(os.path.join(out_dir, "27_curation_discard.png"))

    # (c) retained-data characterization: chem-space retained, dataset sizes, category, reuse.
    #     Three distinct keys placed below their own columns (I / K / L; J needs none).
    rc = pc.reuse_colors()
    fig, axs = stylia.create_figure(1, 4, height=0.6)
    axA = axs.next()
    plot_pct_retained(axA, codes, results, "A", show_y=True)
    inset_legend(axA, loc="lower right")   # of total / of Whole-cell
    plot_box_h(axs.next(), codes, [r["sizes"] for r in results], pc.color("sizes"),
               "Molecules / dataset", "Retained dataset sizes", "B", rng, logx=True, show_y=False)
    axC = axs.next()
    plot_category_h(axC, codes, results, "C", show_y=False)
    inset_legend(axC, loc="lower right")   # Dose Response / Single Point
    axD = axs.next()
    plot_reuse(axD, codes, results, "D", show_y=False)
    inset_legend(axD, loc="lower right")   # 1 / 2 / ≥3 datasets
    stylia.save_figure(os.path.join(out_dir, "27_curation_retained.png"))

    # 4. WHOLE-CELL dataset sizes (before selection): datasets and molecule-instances
    #    distributed across size bins (1x3 + size legend below). A assays per size bin;
    #    B/C molecule-instances per size bin (fraction, then count).
    wc_sizes = [r["wholecell_sizes"] for r in results]
    bin_colors = pc.cyclic(len(SIZE_LABELS))  # ordered size bins (article NPG cyclic)
    fig, axs = stylia.create_figure(1, 3, height=0.6)
    plot_size_stack(axs.next(), codes, wc_sizes, bin_colors, "A", "assays", True, show_y=True)
    plot_size_stack(axs.next(), codes, wc_sizes, bin_colors, "B", "molecules", True, show_y=False)
    plot_size_stack(axs.next(), codes, wc_sizes, bin_colors, "C", "molecules", False, show_y=False)
    figure_legend(fig, list(zip(SIZE_LABELS, bin_colors)),
                  loc="upper center", anchor=(0.5, -0.02), ncol=len(SIZE_LABELS),
                  title="Dataset size (molecules)")
    stylia.save_figure(os.path.join(out_dir, "27_curation_wholecell_sizes.png"))
    print(f"Curation comparison figures (frac_chembl, overlap, target_type, discard, "
          f"retained, wholecell_sizes) -> {out_dir}")


# ============================================================================
# BINARISATION cross-pathogen helpers (per-pathogen figures are in step 22)
# ============================================================================

def plot_ratio_box_grouped(ax, codes, dr_lists, sp_lists, colors, abc, rng, show_y=True,
                           title="Active ratio (DR top · SP bottom)", xlabel="Active ratio",
                           logx=False):
    """Boxplots with jittered points, DR (top) and SP (bottom) per pathogen. logx=True
    -> log x-axis (e.g. molecule counts); otherwise a 0..1 axis (ratios / Jaccard)."""
    pos = np.arange(len(codes))
    for lists, c, dy in [(dr_lists, colors[0], -0.2), (sp_lists, colors[1], 0.2)]:
        ax.boxplot(lists, positions=pos + dy, vert=False, showfliers=False, widths=0.32,
                   patch_artist=False, boxprops={"color": c}, whiskerprops={"color": c},
                   capprops={"color": c}, medianprops={"color": c})
        for i, data in enumerate(lists):
            d = np.asarray([v for v in data if v is not None and not pd.isna(v)], dtype=float)
            if logx:
                d = d[d > 0]
            if len(d) == 0:
                continue
            if len(d) > CAP_POINTS:
                d = rng.choice(d, size=CAP_POINTS, replace=False)
            yj = pos[i] + dy + rng.uniform(-0.12, 0.12, size=len(d))
            ax.scatter(d, yj, s=3, alpha=0.2, color=c, edgecolors="none", zorder=1)
    if logx:
        ax.set_xscale("log")
    else:
        ax.set_xlim(-0.02, 1.02)
    ax.set_yticks(pos); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel=xlabel, ylabel="", title=title, abc=abc)


def plot_stacked_grouped(ax, codes, dr_lists, sp_lists, colors, labels, title, abc,
                         xlabel, show_y=True):
    """Two stacked bars per pathogen (DR top, SP bottom), each split into `labels`
    bins as a fraction summing to 1 within the bar."""
    pos = np.arange(len(codes))
    for lists, dy, lbl in [(dr_lists, -0.2, True), (sp_lists, 0.2, False)]:
        arr = np.array(lists, dtype=float)
        tot = arr.sum(axis=1, keepdims=True); tot[tot == 0] = 1
        frac = arr / tot
        lefts = np.zeros(len(codes))
        for k, color in enumerate(colors):
            ax.barh(pos + dy, frac[:, k], left=lefts, height=0.32, color=color,
                    label=labels[k] if lbl else None)
            lefts += frac[:, k]
    ax.set_yticks(pos); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel=xlabel, ylabel="", title=title, abc=abc)


def plot_grouped_h(ax, codes, a_vals, b_vals, colors, labels, title, abc, xlabel, show_y=True):
    """Two horizontal bars per pathogen (e.g. DR vs SP)."""
    y = np.arange(len(codes)); h = 0.4
    ax.barh(y - h / 2, a_vals, height=h, color=colors[0], label=labels[0])
    ax.barh(y + h / 2, b_vals, height=h, color=colors[1], label=labels[1])
    ax.set_xlim(0, 1.02)
    ax.set_yticks(y); ax.set_yticklabels(codes if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel=xlabel, ylabel="", title=title, abc=abc)


# --- eos1klk cross-pathogen embedding figures -------------------------------

def retained_smiles_map(results):
    """smiles -> set of pathogen codes, restricted to the RETAINED / modelled molecules
    (binarised DR ∪ SP), not the full pre-curation screened set. This is what the eos1klk
    coverage figures highlight, so the picture matches the datasets we actually pool
    (e.g. smansoni goes from 32k screened to 2.3k retained). Maps each pathogen's retained
    InChIKeys to SMILES via its 07_compound_counts, keying the global (SMILES-indexed) coords."""
    sm2path = {}
    for r in results:
        c = r["code"]
        iks = r["tested_iks"]["DR"] | r["tested_iks"]["SP"]
        cc = pd.read_csv(os.path.join(STAGE2PATH, c, "07_compound_counts.csv.gz"),
                         usecols=["inchikey", "smiles"])
        ik2sm = dict(zip(cc["inchikey"], cc["smiles"]))
        for ik in iks:
            s = ik2sm.get(ik)
            if isinstance(s, str):
                sm2path.setdefault(s, set()).add(c)
    return sm2path


def make_eos_pathogens_figure(codes, out_dir, sm2path):
    """Cross-pathogen embedding: colour every RETAINED molecule by the pathogen that
    retains it; molecules retained for >1 pathogen get a single 'shared' colour."""
    coords = load_coords()
    if coords is None:
        return
    idx = coords.index.values
    npath = np.array([len(sm2path.get(s, ())) for s in idx])
    single = np.array([next(iter(sm2path[s])) if (s in sm2path and npath[i] == 1) else ""
                       for i, s in enumerate(idx)])
    pcols = pathogen_colors(codes)
    # per-pathogen groups first (drawn as background), 'shared' LAST so the grey shared
    # molecules are painted on top rather than buried under the single-pathogen colours.
    groups = [(c, single == c, pcols[c]) for c in codes]
    groups.append(("shared (>1 pathogen)", npath > 1, pc.color("neutral")))
    embed_figure(coords, groups, os.path.join(out_dir, "27_eos1klk_all_pathogens.png"),
                 legend_ncol=4)


def make_eos_facets_figure(codes, out_dir, sm2path):
    """t-MAP small multiples: one panel per pathogen, that pathogen's RETAINED / modelled
    molecules highlighted over a faint grey background of the whole eos1klk chemical space.
    Small, low-alpha dots so dense pathogens stay legible."""
    coords = load_coords()
    if coords is None:
        return
    x, y = coords["tmap_x"].values, coords["tmap_y"].values
    idx = coords.index.values
    pcols = pathogen_colors(codes)
    rng = np.random.default_rng(RANDOM_SEED)
    bg = rng.choice(len(idx), size=min(40000, len(idx)), replace=False)  # context sample
    ncols = 5
    nrows = int(np.ceil(len(codes) / ncols))
    fig, axs = stylia.create_figure(nrows, ncols, width=1.0, height=0.21 * nrows)
    flat = np.array(axs.axs).ravel()
    for ax, c in zip(flat, codes):
        mask = np.fromiter((c in sm2path.get(s, ()) for s in idx), dtype=bool, count=len(idx))
        ax.scatter(x[bg], y[bg], s=0.5, alpha=0.08, color=pc.color("neutral"), edgecolors="none", rasterized=True)
        ax.scatter(x[mask], y[mask], s=0.8, alpha=0.35, color=pcols[c], edgecolors="none", rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        stylia.label(ax, xlabel="", ylabel="", title=c)
    for ax in flat[len(codes):]:
        ax.axis("off")
    stylia.save_figure(os.path.join(out_dir, "27_eos1klk_tmap.png"))


def make_binarisation_comparison(results, codes, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    results = sorted(results, key=lambda r: r["n"], reverse=True)
    codes = [r["code"] for r in results]
    rng = np.random.default_rng(RANDOM_SEED)

    drsp = (pc.color("DR"), pc.color("SP"))      # DR / SP, consistent across all figures
    comp_colors = pc.composition_colors()        # all-inactive / mixed / all-active
    hit_colors = pc.cyclic(len(HIT_LABELS))      # distinct colours for hit-count bins

    # 1a. ACTIVE RATIO & composition (1x2 + legends below): active ratio (DR+SP box)
    #     and class composition.
    fig, axs = stylia.create_figure(1, 2, height=0.6)
    axA = axs.next()
    plot_ratio_box_grouped(axA, codes, [r["dr_ratios"] for r in results],
                           [r["sp_ratios"] for r in results], drsp, "A", rng, show_y=True)
    inset_legend(axA, loc="lower right",
                 entries=[("Dose Response", drsp[0]), ("Single Point", drsp[1])])
    axB = axs.next()
    plot_stacked_grouped(axB, codes, [r["composition"]["DR"] for r in results],
                         [r["composition"]["SP"] for r in results], comp_colors, COMPOSITION,
                         "Composition (DR top · SP bottom)", "B", "Fraction of datasets",
                         show_y=False)
    inset_legend(axB, loc="lower right")   # all-inactive/mixed/all-active — title in panel
    stylia.save_figure(os.path.join(out_dir, "27_binarisation_ratio.png"))

    # 1b. HIT REPRODUCIBILITY & DR<->SP agreement (1x3 + legends below): times each active
    #     molecule is hit, contradictions, and DR<->SP concordance.
    fig, axs = stylia.create_figure(1, 3, height=0.6)
    axA = axs.next()
    plot_stacked_grouped(axA, codes, [r["hit_bins"]["DR"] for r in results],
                         [r["hit_bins"]["SP"] for r in results], hit_colors, HIT_LABELS,
                         "Hit reproducibility (DR top · SP bottom)", "A",
                         "Fraction of active molecules", show_y=True)
    inset_legend(axA, loc="lower right")   # 1/2/≥3 assays — 'times hit' is in the title
    axB = axs.next()
    plot_grouped_h(axB, codes, [r["dr_contra"] for r in results],
                   [r["sp_contra"] for r in results], drsp,
                   ("Dose Response", "Single Point"), "Contradictory (tested ≥2×)", "B",
                   "Fraction", show_y=False)
    inset_legend(axB, loc="upper right")   # DR / SP
    axC = axs.next()
    plot_grouped_h(axC, codes, [r["cross_shared"] for r in results],
                   [r["cross_contra"] for r in results], (pc.color("wholecell"), pc.color("related")),
                   ("active in both", "active in one only"), "DR↔SP agreement", "C",
                   "Fraction tested in both", show_y=False)
    inset_legend(axC, loc="upper right")   # active in both / one only — 'DR↔SP' is in the title
    stylia.save_figure(os.path.join(out_dir, "27_binarisation_agreement.png"))

    # 2. ACTIVE chemical-space overlap (square Jaccard heatmaps): which pathogens
    #    share active molecules, DR and SP.
    fig, axs = stylia.create_figure(1, 2, width=1.0, height=0.55)
    plot_overlap_heatmap(axs.next(), codes, [r["active_iks"]["DR"] for r in results],
                         title="Active overlap — DR", abc="A")
    plot_overlap_heatmap(axs.next(), codes, [r["active_iks"]["SP"] for r in results],
                         title="Active overlap — SP", abc="B")
    stylia.save_figure(os.path.join(out_dir, "27_binarisation_active_overlap.png"))

    # (pool structure / quality figures moved out of the binarisation preview: the real
    #  step-23 pool quality is drawn by make_pooled_quality from 23_pool_summary + 24 CV.)

    # cross-pathogen eos1klk embeddings over the RETAINED / modelled molecules:
    #    coloured by pathogen (shared = one colour), and a t-MAP small-multiples facet.
    sm2path = retained_smiles_map(results)
    make_eos_pathogens_figure(codes, out_dir, sm2path)
    make_eos_facets_figure(codes, out_dir, sm2path)
    print(f"Binarisation comparison figures (ratio, agreement, active_overlap) + "
          f"eos1klk (all_pathogens, tmap) -> {out_dir}")


# ============================================================================
# CUTOFF-SENSITIVITY cross-pathogen figure (per-pathogen table is in step 22)
# ============================================================================
# unit -> panel title; order defines panel order. mm/% are the two SP units.
SENS_UNIT_PANELS = [("umol.L-1", "Dose response (µM)"),
                    ("%", "Single point (% effect)"),
                    ("mm", "Single point (zone, mm)")]


def collect_cutoff_sensitivity(codes):
    """Concatenate every pathogen's 22_cutoff_sensitivity.csv, tagged with `pathogen`."""
    frames = []
    for c in codes:
        f = os.path.join(STAGE4PATH, c, "22_cutoff_sensitivity.csv")
        if os.path.exists(f):
            d = pd.read_csv(f)
            d.insert(0, "pathogen", c)
            frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def plot_cutoff_sensitivity_panel(ax, df_unit, codes, abbr, tier_colors, default_tier,
                                  title, abc, show_y=True):
    """One panel per unit: pathogens on y, active ratio on x. Three markers per pathogen
    (low/middle/high candidate cutoff) joined by a faint line to show the achievable
    active-ratio range; the current global default tier is starred."""
    pos = np.arange(len(codes))
    for i, c in enumerate(codes):
        sub = df_unit[df_unit["pathogen"] == c]
        pts = [(t, float(sub.loc[sub["tier"] == t, "active_ratio"].iloc[0]))
               for t in CUTOFF_TIERS if (sub["tier"] == t).any()]
        if not pts:
            continue
        ax.plot([v for _, v in pts], [pos[i]] * len(pts), color=pc.color("neutral"), zorder=1)
        for t, v in pts:
            ax.scatter([v], [pos[i]], color=tier_colors[t], zorder=3)
            if t == default_tier:
                ax.scatter([v], [pos[i]], color=pc.color("ink"), marker="*", zorder=4)
    ax.set_xlim(-0.02, 1.02)
    ax.set_yticks(pos); ax.set_yticklabels([abbr[c] for c in codes] if show_y else [])
    ax.invert_yaxis()
    stylia.label(ax, xlabel="Active ratio", ylabel="", title=title, abc=abc)


def make_cutoff_sensitivity(codes, out_dir):
    """Cross-pathogen cutoff sensitivity: active ratio at each low/middle/high candidate
    cutoff, per (pathogen, unit). DIAGNOSTIC ONLY - it is the evidence for choosing a
    per-pathogen override (config/category_cutoffs_overrides.csv) and changes no label.
    Writes the tidy table and one figure (a panel per unit, current default starred)."""
    sens = collect_cutoff_sensitivity(codes)
    if sens.empty:
        print("  [cutoff-sensitivity] no 22_cutoff_sensitivity.csv found - skipping")
        return
    defaults = load_category_cutoffs(CONFIGPATH)              # {unit: global default cutoff}
    candidates = load_category_cutoff_candidates(CONFIGPATH)  # {unit: {tier: cutoff}}
    sens["is_default"] = [candidates.get(u, {}).get(t) == defaults.get(u)
                          for u, t in zip(sens["unit"], sens["tier"])]
    sens.to_csv(os.path.join(out_dir, "27_cutoff_sensitivity.csv"), index=False)

    def default_tier(unit):
        for t, v in candidates.get(unit, {}).items():
            if unit in defaults and v == defaults[unit]:
                return t
        return None

    panels = [(u, t) for u, t in SENS_UNIT_PANELS if (sens["unit"] == u).any()]
    abbr = pathogen_abbrev(codes)
    tcols = pc.cyclic(len(CUTOFF_TIERS))
    tier_colors = {t: tcols[i] for i, t in enumerate(CUTOFF_TIERS)}
    # order pathogens by DR (µM) middle-cutoff active ratio, descending; missing last
    dr_mid = (sens[(sens["unit"] == "umol.L-1") & (sens["tier"] == "middle")]
              .set_index("pathogen")["active_ratio"])
    order = sorted(codes, key=lambda c: dr_mid.get(c, -1.0), reverse=True)

    fig, axs = stylia.create_figure(1, len(panels))
    for j, (unit, title) in enumerate(panels):
        plot_cutoff_sensitivity_panel(
            axs.next(), sens[sens["unit"] == unit], order, abbr, tier_colors,
            default_tier(unit), title, abc=chr(65 + j), show_y=(j == 0))
    figure_legend(fig, [(t, tier_colors[t]) for t in CUTOFF_TIERS],
                  loc="upper center", anchor=(0.35, -0.02), ncol=len(CUTOFF_TIERS),
                  title="Cutoff tier")
    fig.legend(handles=[Line2D([], [], marker="*", linestyle="none",
                               color=pc.color("ink"), label="current global default")],
               loc="upper center", bbox_to_anchor=(0.8, -0.02), frameon=False)
    stylia.save_figure(os.path.join(out_dir, "27_cutoff_sensitivity.png"))
    print(f"  cutoff-sensitivity -> {os.path.join(out_dir, '27_cutoff_sensitivity.png')} "
          f"(+ 27_cutoff_sensitivity.csv)")


def available_pathogens():
    codes = pd.read_csv(os.path.join(root, "..", "config", "pathogens.csv"))["code"].tolist()
    return [c for c in codes
            if os.path.exists(os.path.join(STAGE4PATH, c, "22_binarisation_summary.csv"))]


# ============================================================================
# POOL-SUMMARY cross-pathogen figures (moved here from the transfer step, step 23)
# ============================================================================

def _collect(codes, fname):
    """Concatenate a per-pathogen CSV across the codes that have it."""
    frames = []
    for c in codes:
        f = os.path.join(STAGE4PATH, c, fname)
        if os.path.exists(f):
            frames.append(pd.read_csv(f))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _read1(code, fname):
    p = os.path.join(STAGE4PATH, code, fname)
    return pd.read_csv(p) if os.path.exists(p) else None


def build_master_table(codes):
    """The pooling backbone: one row per original binarised dataset, traced through
    fold -> pool (23) -> CV (24) -> merge (25) / catch-all (26). Joins the per-pathogen
    CSVs (dataset_id = '<assay_id>_<activity_type>_<unit_slug>', matching the pool-member
    ids). Adds `final_status` (modelled / catch-all / uncovered) and a single `final_auroc`
    (grown_auroc if merged into a modelled pool, else step-24 auroc, else step-26 catch-all
    auroc). Drives every all-pathogen pooling figure + is the tabular summary of all steps."""
    def one(code):
        d22 = _read1(code, "22_binarisation_summary.csv")
        if d22 is None:
            return None
        d22 = d22[d22["category"].isin(["DR", "SP"])].copy()
        at = d22["activity_type"].astype(str).str.replace("/", "-", regex=False)
        slug = d22["unit"].map(lambda u: UNIT_SLUG.get(u, str(u)))
        d22["dataset_id"] = d22["assay_id"].astype(str) + "_" + at + "_" + slug
        d22["pathogen"] = code
        m = d22[["pathogen", "category", "dataset_id", "assay_id", "activity_type",
                 "n_molecules", "active_ratio", "degenerate"]].rename(
            columns={"n_molecules": "ds_n_molecules", "active_ratio": "ds_active_ratio",
                     "degenerate": "ds_degenerate"})
        pm23 = _read1(code, "23_pool_members.csv")
        if pm23 is not None:
            m = m.merge(pm23[["category", "dataset_id", "pool_id"]].rename(columns={"pool_id": "pool_id_23"}),
                        on=["category", "dataset_id"], how="left")
        else:
            m["pool_id_23"] = np.nan
        ps23 = _read1(code, "23_pool_summary.csv")
        if ps23 is not None:
            m = m.merge(ps23[["category", "pool_id", "is_singleton"]].rename(columns={"pool_id": "pool_id_23"}),
                        on=["category", "pool_id_23"], how="left")
        cv24 = _read1(code, "24_cv_summary.csv")
        if cv24 is not None:
            m = m.merge(cv24[["category", "pool_id", "modelled", "auroc", "active_ratio", "reason"]].rename(
                columns={"pool_id": "pool_id_23", "auroc": "auroc_24", "active_ratio": "ar_24",
                         "reason": "reason_24"}),
                on=["category", "pool_id_23"], how="left")
        pm25 = _read1(code, "25_pool_members.csv")
        if pm25 is not None:
            m = m.merge(pm25[["category", "dataset_id", "pool_id"]].rename(columns={"pool_id": "pool_id_25"}),
                        on=["category", "dataset_id"], how="left")
        else:
            m["pool_id_25"] = np.nan
        ps25 = _read1(code, "25_pool_summary.csv")
        if ps25 is not None:
            m = m.merge(ps25[["category", "pool_id", "grown_auroc", "baseline_auroc",
                              "batch_accepted", "ar_after"]].rename(columns={"pool_id": "pool_id_25"}),
                        on=["category", "pool_id_25"], how="left")
        ml25 = _read1(code, "25_merge_log.csv")
        if ml25 is not None:
            m = m.merge(ml25[["category", "leftover", "outcome"]].rename(
                columns={"leftover": "dataset_id", "outcome": "leftover_outcome"}),
                on=["category", "dataset_id"], how="left")
        cv26 = _read1(code, "26_cv_summary.csv")
        if cv26 is not None:
            m = m.merge(cv26[["category", "auroc", "active_ratio", "kept"]].rename(
                columns={"auroc": "catchall_auroc", "active_ratio": "catchall_ar", "kept": "catchall_kept"}),
                on=["category"], how="left")
        for col in ["pool_id_25", "grown_auroc", "baseline_auroc", "batch_accepted", "ar_after",
                    "modelled", "auroc_24", "ar_24", "reason_24", "leftover_outcome",
                    "catchall_auroc", "catchall_ar", "catchall_kept", "is_singleton"]:
            if col not in m.columns:
                m[col] = np.nan
        return m

    frames = [one(c) for c in codes]
    m = pd.concat([f for f in frames if f is not None], ignore_index=True)

    def resolve(r):
        if pd.notna(r["pool_id_25"]):
            return "modelled", r["grown_auroc"], r["ar_after"], r["pool_id_25"]
        if r["catchall_kept"] == True:
            return "catch-all", r["catchall_auroc"], r["catchall_ar"], f'{r["category"]}_catchall'
        if pd.notna(r["pool_id_23"]) and (r["modelled"] is True):
            return "modelled", r["auroc_24"], r["ar_24"], r["pool_id_23"]
        return "uncovered", np.nan, np.nan, np.nan
    res = m.apply(resolve, axis=1, result_type="expand")
    m["final_status"], m["final_auroc"], m["final_active_ratio"], m["final_pool"] = \
        res[0], res[1], res[2], res[3]
    return m


ABBREV_OVERRIDE = {"mtuberculosis": "Mtb"}   # conventional short forms that override the default rule


def pathogen_abbrev(codes):
    """Short axis label per code from the scientific name: genus initial + first two species
    letters (Staphylococcus aureus -> 'Sau'); single-word genera -> first three ('Cam')."""
    out = {}
    for c in codes:
        if c in ABBREV_OVERRIDE:
            out[c] = ABBREV_OVERRIDE[c]
            continue
        parts = str(load_pathogen(c)).split()
        out[c] = (parts[0][0] + parts[1][:2]).capitalize() if len(parts) >= 2 else parts[0][:3].capitalize()
    return out


def plot_pool_sizes(ax, summ, codes, colors, abbr, category, abc):
    """Box of pool sizes (molecules per pool, log scale) per pathogen."""
    sub = summ[summ["category"] == category]
    data = [sub[sub["pathogen"] == c]["n_molecules"].to_numpy() for c in codes]
    bp = ax.boxplot(data, positions=np.arange(len(codes)), widths=0.6, patch_artist=True,
                    showfliers=False)
    for patch, c in zip(bp["boxes"], codes):
        patch.set_facecolor(colors[c])
    for med in bp["medians"]:
        med.set_color("black")
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(codes))); ax.set_xticklabels([abbr[c] for c in codes], rotation=90)
    stylia.label(ax, xlabel="", ylabel="molecules / pool",
                 title=f"Step-23 pool size - {cat_label(category)}", abc=abc)


def plot_chemspace(ax, part, codes, abbr, category, abc):
    """Stacked bar per pathogen: molecule fraction in multi / singleton / not-considered.
    Every pathogen in `codes` is shown; one with no step-23 partition on disk (skipped for
    having < 2 eligible datasets, e.g. campylobacter) is rendered as 100% not-considered."""
    sub = part[part["category"] == category].set_index("pathogen")
    x = np.arange(len(codes))
    f_multi, f_single, f_not = [], [], []
    for c in codes:
        if c in sub.index:
            r = sub.loc[c]
            tot = float(r["n_total"])
            if tot > 0:
                f_multi.append(r["n_multi"] / tot)
                f_single.append(r["n_singleton"] / tot)
                f_not.append(r["n_not_considered"] / tot)
                continue
        # missing partition (pathogen skipped at step 23) or empty universe -> all not-considered
        f_multi.append(0.0); f_single.append(0.0); f_not.append(1.0)
    f_multi, f_single, f_not = np.array(f_multi), np.array(f_single), np.array(f_not)
    ax.bar(x, f_multi, color=pc.color("datasets"), label="multi-dataset pool")
    ax.bar(x, f_single, bottom=f_multi, color=pc.named("periwinkle"), label="singleton")
    ax.bar(x, f_not, bottom=f_multi + f_single, color=pc.color("neutral"), label="not considered")
    ax.set_xticks(x); ax.set_xticklabels([abbr[c] for c in codes], rotation=90)
    ax.set_ylim(0, 1)
    inset_legend(ax, loc="upper right")   # multi / singleton / not considered (labelled bars)
    stylia.label(ax, xlabel="", ylabel="fraction of molecules",
                 title=f"Step-23 pool partition - {cat_label(category)}", abc=abc)


# dot AREA (points^2) encodes pool molecule count: SIZE_MIN..SIZE_MAX mapped linearly from n_molecules
SIZE_MIN, SIZE_MAX = 15.0, 600.0


def _area_for(n, nmax):
    """Marker area proportional to n_molecules (with a small floor so tiny pools stay visible)."""
    return SIZE_MIN + (np.asarray(n, dtype=float) / nmax) * (SIZE_MAX - SIZE_MIN)


def _size_legend_handles(nmax):
    """Proxy handles (gray bubbles) for a shared molecule-count size legend."""
    refs = [n for n in (100, 1000, 10000, 100000) if n <= nmax] or [int(nmax)]
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=np.sqrt(_area_for(n, nmax)),
                      markerfacecolor="lightgray", markeredgecolor="white", markeredgewidth=0.5)
               for n in refs]
    return handles, [f"{n:,}" for n in refs]


def plot_active_ratios(ax, summ, codes, colors, abbr, category, abc, nmax):
    """One dot per CONSIDERED pool — multi-dataset AND singleton (every pool in 23_pool_summary
    is an eligible node) — y = active ratio, area proportional to n_molecules, with a median line
    per pathogen, so pool size and activity are read together."""
    sub = summ[summ["category"] == category]
    cc = [c for c in codes if (sub["pathogen"] == c).any()]
    rng = np.random.default_rng(RANDOM_SEED)
    for xi, c in enumerate(cc):
        g = sub[sub["pathogen"] == c].sort_values("n_molecules", ascending=False)  # big dots first (behind)
        y = g["active_ratio"].to_numpy(dtype=float)
        x = xi + rng.uniform(-0.18, 0.18, len(y))
        ax.scatter(x, y, s=_area_for(g["n_molecules"].to_numpy(), nmax), color=colors[c],
                   alpha=0.55, edgecolor="white", linewidth=0.5, zorder=2)
        ax.hlines(float(np.median(y)), xi - 0.35, xi + 0.35, color="black", zorder=3)  # median = the line
    ax.set_ylim(-0.08, 1.08)                        # headroom so big dots at 0/1 don't clip
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_xticks(np.arange(len(cc))); ax.set_xticklabels([abbr[c] for c in cc], rotation=90)
    stylia.label(ax, xlabel="", ylabel="active ratio",
                 title=f"Step-23 pool active ratio - {cat_label(category)}", abc=abc)


def make_general_plots(codes, out_dir):
    """Per-pathogen step-23 pool figures (pool partition + pool active ratio, both including
    singletons), aggregated from each pathogen's step-23 outputs. Pool size / modellability now
    live in the condensed 27_pool_overview (make_pool_overview)."""
    summ = _collect(codes, "23_pool_summary.csv")
    if summ.empty:
        print("  [general] no 23_pool_summary.csv found - skipping pool-summary figures")
        return
    part = _collect(codes, "23_chemspace_partition.csv")
    colors = pathogen_colors(codes)
    abbr = {c: c for c in codes}   # full pathogen codes on the axis (no abbreviation)

    if not part.empty:
        fig, axs = stylia.create_figure(1, 2)
        plot_chemspace(axs.next(), part, codes, abbr, "DR", "A")
        plot_chemspace(axs.next(), part, codes, abbr, "SP", "B")
        stylia.save_figure(os.path.join(out_dir, "27_pool_partition.png"))

    nmax = max(1.0, float(summ["n_molecules"].max()))   # all considered pools (incl. singletons)
    fig, axs = stylia.create_figure(1, 2)
    plot_active_ratios(axs.next(), summ, codes, colors, abbr, "DR", "A", nmax)
    plot_active_ratios(axs.next(), summ, codes, colors, abbr, "SP", "B", nmax)
    handles, labels = _size_legend_handles(nmax)
    fig.legend(handles, labels, title="molecules / pool", loc="center left",
               bbox_to_anchor=(1.0, 0.5), labelspacing=1.4, borderpad=0.8)
    stylia.save_figure(os.path.join(out_dir, "27_pool_active_ratios.png"))
    print("  [general] wrote 27_pool_partition.png, 27_pool_active_ratios.png")


# step-24 CV outcome bins (reason value -> label, colour). "<100 mol" pools (multi-dataset or
# singleton) are below the symmetric min-molecule floor at step 24 and are rescued downstream by steps 25/26.
CV_OUTCOME = [("modelled", "modelled", "cobalt"),
              ("<100 mol", "below_min_molecules", "silver"),
              ("not CV-able", "not_cv_able", "amber")]


def plot_auroc_box(ax, cv, codes, colors, category, abc):
    """Empty box + one jittered dot per modelled pool (dots coloured by pathogen) of pool CV
    AUROC, per pathogen; a 0.7 guide line marks the modellability floor used downstream."""
    sub = cv[(cv["category"] == category) & (cv["reason"] == "modelled")]
    rng = np.random.default_rng(RANDOM_SEED)
    box_color = pc.color("ink")
    for xi, c in enumerate(codes):
        v = sub[sub["pathogen"] == c]["auroc"].dropna().to_numpy(float)
        if not len(v):
            continue
        ax.boxplot([v], positions=[xi], widths=0.6, showfliers=False, patch_artist=False,
                   boxprops={"color": box_color}, whiskerprops={"color": box_color},
                   capprops={"color": box_color}, medianprops={"color": box_color})
        xj = xi + rng.uniform(-0.18, 0.18, size=len(v))
        ax.scatter(xj, v, s=8, alpha=0.55, color=colors[c], edgecolors="none", zorder=2)
    ax.axhline(0.7, ls="--", lw=0.8, color=pc.color("neutral"))
    ax.set_ylim(0.4, 1.02)
    ax.set_xticks(np.arange(len(codes))); ax.set_xticklabels(codes, rotation=90)
    stylia.label(ax, xlabel="", ylabel="pool CV AUROC",
                 title=f"Modellability - {cat_label(category)}", abc=abc)


def plot_cv_outcome(ax, cv, codes, category, abc, normalize=False):
    """Stacked bar per pathogen of step-23 pools by step-24 CV outcome (modelled /
    <100 mol / not-CV-able). normalize=False -> absolute pool counts (bar height =
    total pools, so this also carries the pool count); True -> fractions summing to 1."""
    sub = cv[cv["category"] == category]
    x = np.arange(len(codes)); bottoms = np.zeros(len(codes))
    for lab, key, cname in CV_OUTCOME:
        raw = np.array([(sub[sub["pathogen"] == c]["reason"] == key).sum() for c in codes], dtype=float)
        vals = raw / np.array([max(len(sub[sub["pathogen"] == c]), 1) for c in codes]) if normalize else raw
        ax.bar(x, vals, bottom=bottoms, color=pc.named(cname), label=lab)
        bottoms += vals
    ax.set_xticks(x); ax.set_xticklabels(codes, rotation=90)
    ylabel = "fraction of pools" if normalize else "pools"
    # headroom above the tallest bar so it never touches the top edge and the in-axis legend fits
    ax.set_ylim(0, 1 if normalize else (float(bottoms.max()) * 1.18 if bottoms.max() > 0 else 1))
    stylia.label(ax, xlabel="", ylabel=ylabel,
                 title=f"Step-24 CV outcome of step-23 pools - {cat_label(category)}", abc=abc)
    inset_legend(ax, loc="upper right")   # legend inside the panel (semi-transparent white bg)


def make_pool_overview(codes, out_dir):
    """Condensed per-pathogen pool figure (3x2, DR left / SP right) spanning step-23 structure
    and step-24 modellability, so pool size, CV AUROC and CV outcome are read together (replaces
    the separate 27_pool_summary + 27_pooled_quality):
        row 1  step-23 pool size (box, log)
        row 2  step-24 pool CV AUROC (empty box + one dot per modelled pool)
        row 3  step-24 CV outcome as absolute pool counts (bar height = pool count; the stack
               composition folds in the old multi/singleton split, so nothing is duplicated)."""
    summ = _collect(codes, "23_pool_summary.csv")
    cv = _collect(codes, "24_cv_summary.csv")
    if summ.empty or cv.empty:
        print("  [general] missing 23_pool_summary / 24_cv_summary - skipping pool_overview")
        return
    colors = pathogen_colors(codes)
    abbr = {c: c for c in codes}
    fig, axs = stylia.create_figure(3, 2, height=1.5)
    plot_pool_sizes(axs.next(), summ, codes, colors, abbr, "DR", "A")
    plot_pool_sizes(axs.next(), summ, codes, colors, abbr, "SP", "B")
    plot_auroc_box(axs.next(), cv, codes, colors, "DR", "C")
    plot_auroc_box(axs.next(), cv, codes, colors, "SP", "D")
    plot_cv_outcome(axs.next(), cv, codes, "DR", "E", normalize=False)  # legend now inside E/F
    plot_cv_outcome(axs.next(), cv, codes, "SP", "F", normalize=False)
    stylia.save_figure(os.path.join(out_dir, "27_pool_overview.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_pool_overview.png')}")


# ============================================================================
# MERGE overview (step 25 — leftover merging into the modelled pools)
# ============================================================================

# 25_merge_log `outcome` value -> (label, colour)
MERGE_OUTCOME = [("merged", "merged", "turquoise"),
                 ("deferred (low fit)", "deferred_low_fit", "crimson"),
                 ("deferred (AUROC drop)", "deferred_auroc_drop", "amber"),
                 ("redundant", "redundant", "silver")]


def plot_auroc_lift(ax, ps, codes, category, abc):
    """Grouped box per pathogen: baseline (pre-merge) vs grown (post-merge) pool CV AUROC.
    Baseline = silver, grown = turquoise; 0.7 guide line. Shows whether merging leftovers into
    a seed pool helped or hurt its cross-validated modellability."""
    base_fill = pc.lighten_named("turquoise", 0.55)   # faded -> full turquoise reads as pre -> post merge
    sub = ps[ps["category"] == category]
    for xi, c in enumerate(codes):
        g = sub[sub["pathogen"] == c]
        for dx, col, fill in [(-0.2, "baseline_auroc", base_fill),
                              (0.2, "grown_auroc", pc.named("turquoise"))]:
            v = g[col].dropna().to_numpy(float)
            if not len(v):
                continue
            bp = ax.boxplot([v], positions=[xi + dx], widths=0.34, showfliers=False,
                            patch_artist=True)
            bp["boxes"][0].set_facecolor(fill)
            bp["medians"][0].set_color("black")
    ax.axhline(0.7, ls="--", lw=0.8, color=pc.color("neutral"))
    ax.set_ylim(0.4, 1.02)
    ax.set_xticks(np.arange(len(codes))); ax.set_xticklabels(codes, rotation=90)
    stylia.label(ax, xlabel="", ylabel="pool CV AUROC",
                 title=f"Step-25 baseline → grown AUROC - {cat_label(category)}", abc=abc)


def plot_leftover_outcome(ax, ml, codes, category, abc):
    """Stacked fraction per pathogen of leftover datasets by step-25 outcome
    (merged into a pool / redundant / deferred because it dropped the pool AUROC)."""
    sub = ml[ml["category"] == category]
    x = np.arange(len(codes)); bottoms = np.zeros(len(codes))
    for lab, key, cname in MERGE_OUTCOME:
        fr = np.array([((sub[sub["pathogen"] == c]["outcome"] == key).sum() /
                        max(len(sub[sub["pathogen"] == c]), 1)) for c in codes])
        ax.bar(x, fr, bottom=bottoms, color=pc.named(cname), label=lab)
        bottoms += fr
    ax.set_xticks(x); ax.set_xticklabels(codes, rotation=90); ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel="fraction of leftovers",
                 title=f"Step-25 leftover outcome - {cat_label(category)}", abc=abc)


def plot_merge_coverage(ax, ps, codes, colors, category, abc):
    """Bar per pathogen: molecule instances absorbed into pools by the step-25 merge (sum of
    n_added over that pathogen's grown pools), coloured by pathogen."""
    sub = ps[ps["category"] == category]
    vals = [float(sub[sub["pathogen"] == c]["n_added"].sum()) for c in codes]
    x = np.arange(len(codes))
    ax.bar(x, vals, color=[colors[c] for c in codes])
    ax.set_xticks(x); ax.set_xticklabels(codes, rotation=90)
    compact_yaxis(ax)
    stylia.label(ax, xlabel="", ylabel="molecules added",
                 title=f"Step-25 molecules absorbed - {cat_label(category)}", abc=abc)


def make_merge_overview(codes, out_dir):
    """Per-pathogen step-25 merge overview (3x2, DR left / SP right) from 25_pool_summary +
    25_merge_log: baseline→grown AUROC lift, leftover-outcome mix, and molecules absorbed."""
    ps = _collect(codes, "25_pool_summary.csv")
    ml = _collect(codes, "25_merge_log.csv")
    if ps.empty and ml.empty:
        print("  [general] no 25_pool_summary / 25_merge_log found - skipping merge_overview")
        return
    colors = pathogen_colors(codes)
    fig, axs = stylia.create_figure(3, 2, height=1.5)
    plot_auroc_lift(axs.next(), ps, codes, "DR", "A")
    plot_auroc_lift(axs.next(), ps, codes, "SP", "B")
    plot_leftover_outcome(axs.next(), ml, codes, "DR", "C")
    plot_leftover_outcome(axs.next(), ml, codes, "SP", "D")
    plot_merge_coverage(axs.next(), ps, codes, colors, "DR", "E")
    plot_merge_coverage(axs.next(), ps, codes, colors, "SP", "F")
    figure_legend(fig, [("baseline (pre-merge)", pc.lighten_named("turquoise", 0.55)),
                        ("grown (post-merge)", pc.named("turquoise"))],
                  loc="upper center", anchor=(0.3, -0.01), ncol=2)
    figure_legend(fig, [(lab, pc.named(cname)) for lab, _key, cname in MERGE_OUTCOME],
                  loc="upper center", anchor=(0.75, -0.01), ncol=3)
    stylia.save_figure(os.path.join(out_dir, "27_merge_overview.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_merge_overview.png')}")


# ============================================================================
# FINAL-DATA overlap (molecule sharing + label conflicts between the FINAL pools)
# ============================================================================

def _final_pool_stats(codes):
    """Load each pathogen's FINAL pools (25 grown + 26 catch-all) ONCE and derive, per category:
      - pairs : one row per pool pair sharing >=1 molecule (n_shared, Jaccard %, discrepancy %);
      - counts: n_pools, n_pairs (all C(n,2)), n_overlap (pairs sharing >=1 molecule);
      - confl : molecule-level n_unique_molecules / n_in_multiple_pools / n_conflicting_label over
                the SAME combined pool set (so C/D matches A/B, catch-all included).
    Discrepancy = % of the shared InChIKeys whose binary label (`bin`) disagrees; a molecule is
    'conflicting' if it appears in >1 pool with >1 distinct bin."""
    pairs, counts, confl = [], [], []
    for code in codes:
        for cat in ("DR", "SP"):
            pools = {}
            for base in ("25_pools", "26_pools"):
                for f in sorted(glob.glob(os.path.join(STAGE4PATH, code, base, cat, "*.csv.gz"))):
                    pid = os.path.basename(f).replace(".csv.gz", "")
                    d = pd.read_csv(f, usecols=["inchikey", "bin"])
                    pools[pid] = d.groupby("inchikey")["bin"].first()
            ids = list(pools)
            ik_npools, ik_bins = defaultdict(int), defaultdict(set)
            for s in pools.values():
                for ik, b in s.items():
                    ik_npools[ik] += 1
                    ik_bins[ik].add(int(b))
            if ik_npools:
                confl.append(dict(pathogen=code, category=cat, n_unique_molecules=len(ik_npools),
                                  n_in_multiple_pools=sum(1 for n in ik_npools.values() if n > 1),
                                  n_conflicting_label=sum(1 for ik, n in ik_npools.items()
                                                          if n > 1 and len(ik_bins[ik]) > 1)))
            n_ov = 0
            for a, b in itertools.combinations(ids, 2):
                sa, sb = pools[a], pools[b]
                ka, kb = set(sa.index), set(sb.index)
                inter = ka & kb
                if not inter:
                    continue
                n_ov += 1
                ba = sa.loc[list(inter)].values
                bb = sb.loc[list(inter)].values
                pairs.append(dict(pathogen=code, category=cat, pool_a=a, pool_b=b,
                                  n_shared=len(inter), jaccard=100.0 * len(inter) / len(ka | kb),
                                  discrepancy=100.0 * float((ba != bb).mean())))
            counts.append(dict(pathogen=code, category=cat, n_pools=len(ids),
                               n_pairs=len(ids) * (len(ids) - 1) // 2, n_overlap=n_ov))
    return pd.DataFrame(pairs), pd.DataFrame(counts), pd.DataFrame(confl)


def plot_pair_box(ax, sub, codes, colors, value_col, title, ylabel, abc, logy=False, counts=None):
    """Empty box + one jittered dot per FINAL-pool pair (dots coloured by pathogen). When `counts`
    is given, annotate each pathogen with 'n_overlap / n_pairs' (how many of its pool pairs share
    >=1 molecule, out of all possible pairs) so the box is read against the total, not in a vacuum."""
    rng = np.random.default_rng(RANDOM_SEED)
    box = pc.color("ink")
    for xi, c in enumerate(codes):
        v = sub.loc[sub["pathogen"] == c, value_col].to_numpy(float)
        v = v[np.isfinite(v)]
        if logy:
            v = v[v > 0]
        if len(v):
            ax.boxplot([v], positions=[xi], widths=0.6, showfliers=False, patch_artist=False,
                       boxprops={"color": box}, whiskerprops={"color": box},
                       capprops={"color": box}, medianprops={"color": box})
            ax.scatter(xi + rng.uniform(-0.18, 0.18, size=len(v)), v, s=8, alpha=0.5,
                       color=colors[c], edgecolors="none", zorder=2)
        if counts is not None:
            cc = counts[counts["pathogen"] == c]
            if len(cc) and int(cc.iloc[0]["n_pairs"]) > 0:
                ax.text(xi, 0.98, f"{int(cc.iloc[0]['n_overlap'])}/{int(cc.iloc[0]['n_pairs'])}",
                        transform=blended_transform_factory(ax.transData, ax.transAxes),
                        ha="center", va="top", rotation=90, fontsize=stylia.FONTSIZE_SMALL,
                        color=pc.color("ink"),
                        bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=0.4))
    if logy:
        ax.set_yscale("log")
    ax.set_xticks(np.arange(len(codes))); ax.set_xticklabels(codes, rotation=90)
    stylia.label(ax, xlabel="", ylabel=ylabel, title=title, abc=abc)


def plot_mol_conflict(ax, lc, codes, category, abc):
    """Stacked bar per pathogen over its FINAL-pool unique molecules: the fraction reused in
    >=2 pools, split into consistent vs conflicting labels. Bar height = reuse rate; the crimson
    (bottom) segment = conflicting molecules / total unique molecules — i.e. the cross-pool
    label conflict normalised over the whole dataset (not just the shared molecules), so it is
    free of the 1-shared-molecule artifact the pairwise discrepancy has. Computed over the same
    final-pool set as A/B (25 grown + 26 catch-all)."""
    sub = lc[lc["category"] == category].set_index("pathogen")
    x = np.arange(len(codes))

    def frac(c, num):
        if c not in sub.index:
            return 0.0
        tot = sub.loc[c, "n_unique_molecules"]
        return float(sub.loc[c, num] / tot) if tot else 0.0

    f_conf = np.array([frac(c, "n_conflicting_label") for c in codes])
    f_cons = np.array([(frac(c, "n_in_multiple_pools") - frac(c, "n_conflicting_label")) for c in codes])
    ax.bar(x, f_conf, color=pc.named("crimson"), label="reused, conflicting")
    ax.bar(x, f_cons, bottom=f_conf, color=pc.named("turquoise"), label="reused, consistent")
    ax.set_xticks(x); ax.set_xticklabels(codes, rotation=90)
    stylia.label(ax, xlabel="", ylabel="fraction of unique molecules",
                 title=f"Cross-pool label conflict - {cat_label(category)}", abc=abc)


def make_final_data_overlap(codes, out_dir):
    """Per-pathogen molecule sharing + label conflicts between the FINAL pools (2x2, DR left /
    SP right):
        row 1  pairwise Jaccard overlap % (log) between pool pairs sharing >=1 molecule;
        row 2  cross-pool label conflict at the MOLECULE level (stacked bar per pathogen: reuse
               rate split into consistent vs conflicting, the conflicting segment = conflicting
               molecules / total unique molecules).
    Writes the tidy pairwise table (incl. per-pair discrepancy) to 27_final_data_overlap.csv."""
    pairs, counts, confl = _final_pool_stats(codes)
    if pairs.empty and confl.empty:
        print("  [general] no final pool files (25_pools/26_pools) - skipping final_data_overlap")
        return
    pairs.to_csv(os.path.join(out_dir, "27_final_data_overlap.csv"), index=False)
    colors = pathogen_colors(codes)
    fig, axs = stylia.create_figure(2, 2, height=0.6)
    plot_pair_box(axs.next(), pairs[pairs.category == "DR"], codes, colors, "jaccard",
                  f"Final-pool overlap - {cat_label('DR')}", "overlap % (Jaccard)", "A",
                  logy=True, counts=counts[counts.category == "DR"])
    plot_pair_box(axs.next(), pairs[pairs.category == "SP"], codes, colors, "jaccard",
                  f"Final-pool overlap - {cat_label('SP')}", "overlap % (Jaccard)", "B",
                  logy=True, counts=counts[counts.category == "SP"])
    if not confl.empty:
        plot_mol_conflict(axs.next(), confl, codes, "DR", "C")
        plot_mol_conflict(axs.next(), confl, codes, "SP", "D")
        figure_legend(fig, [("reused, conflicting", pc.named("crimson")),
                            ("reused, consistent", pc.named("turquoise"))],
                      loc="upper center", anchor=(0.5, -0.01), ncol=2)
    stylia.save_figure(os.path.join(out_dir, "27_final_data_overlap.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_final_data_overlap.png')}")


# ============================================================================
# ACTIVITY RATIO through the pipeline (step 23 -> 26)
# ============================================================================

def plot_stage_ratio(ax, category, d22, p23, c24, p25, c26, abc):
    """Active-ratio distribution at each pipeline stage (all pathogens pooled), one empty box +
    jittered dots per stage, following the funnel:
        22 datasets    - all binarised individual datasets (22_binarisation_summary);
        23 pools       - the eligible pools (pooled + singletons; degenerate / <MIN_DATASET_SIZE
                         datasets are already excluded upstream);
        24 modelled    - the pools that passed cross-validation;
        25 final pools - the grown pools after merging leftovers (ar_after);
        26 catch-all   - the low-data catch-all pools.
    Dots are capped at CAP_POINTS per stage (step 22 has ~28k datasets); the box uses all data."""
    def stage_col(df, value_col, reason=None):
        # a stage whose per-pathogen file was absent for every code yields a column-less
        # empty frame (_collect); treat any frame missing the needed columns as empty
        # rather than letting the ["category"] index raise KeyError.
        if "category" not in df.columns or value_col not in df.columns:
            return pd.Series([], dtype=float)
        sel = df["category"] == category
        if reason is not None:
            if "reason" not in df.columns:
                return pd.Series([], dtype=float)
            sel = sel & (df["reason"] == reason)
        return df.loc[sel, value_col]

    stages = [
        ("22 datasets",    stage_col(d22, "active_ratio")),
        ("23 pools",       stage_col(p23, "active_ratio")),
        ("24 modelled",    stage_col(c24, "active_ratio", reason="modelled")),
        ("25 final pools", stage_col(p25, "ar_after")),
        ("26 catch-all",   stage_col(c26, "active_ratio")),
    ]
    rng = np.random.default_rng(RANDOM_SEED)
    box, dot = pc.color("ink"), pc.color(category)
    for i, (_lab, s) in enumerate(stages):
        v = np.asarray(s.dropna(), dtype=float); v = v[np.isfinite(v)]
        if not len(v):
            continue
        ax.boxplot([v], positions=[i], widths=0.6, showfliers=False, patch_artist=False,
                   boxprops={"color": box}, whiskerprops={"color": box},
                   capprops={"color": box}, medianprops={"color": box})
        d = v if len(v) <= CAP_POINTS else rng.choice(v, size=CAP_POINTS, replace=False)
        ax.scatter(i + rng.uniform(-0.18, 0.18, size=len(d)), d, s=8, alpha=0.35,
                   color=dot, edgecolors="none", zorder=2)
    ax.axhline(0.5, ls="--", lw=0.8, color=pc.color("neutral"))
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels([lab for lab, _ in stages], rotation=30, ha="right")
    stylia.label(ax, xlabel="", ylabel="active ratio",
                 title=f"Active ratio by stage - {cat_label(category)}", abc=abc)


def plot_ratio_beforeafter(ax, p25, c26, codes, category, abc, before_fill, after_fill):
    """Per pathogen, two grouped boxes: seed pool active ratio (ar_before, step 25) vs FINAL pool
    active ratio (step-25 grown ar_after PLUS step-26 catch-all). Catch-all pools have no seed, so
    they only appear in the final box. 0.5 = class balance."""
    sub = p25[p25["category"] == category]
    c26c = c26[c26["category"] == category]
    for xi, c in enumerate(codes):
        before = sub.loc[sub["pathogen"] == c, "ar_before"].dropna().to_numpy(float)
        after = np.concatenate([
            sub.loc[sub["pathogen"] == c, "ar_after"].dropna().to_numpy(float),
            c26c.loc[c26c["pathogen"] == c, "active_ratio"].dropna().to_numpy(float)])
        for dx, v, fill in [(-0.2, before, before_fill), (0.2, after, after_fill)]:
            if not len(v):
                continue
            bp = ax.boxplot([v], positions=[xi + dx], widths=0.34, showfliers=False, patch_artist=True)
            bp["boxes"][0].set_facecolor(fill)
            bp["medians"][0].set_color("black")
    ax.axhline(0.5, ls="--", lw=0.8, color=pc.color("neutral"))
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks(np.arange(len(codes))); ax.set_xticklabels(codes, rotation=90)
    stylia.label(ax, xlabel="", ylabel="active ratio",
                 title=f"Seed → final pool ratio - {cat_label(category)}", abc=abc)


def make_activity_ratio(codes, out_dir):
    """Active ratio through the pooling pipeline (2x2): top row = distribution at each stage
    22->26 (DR/SP, all pathogens pooled); bottom row = per-pathogen seed pool vs final pool ratio
    (final = step-25 grown + step-26 catch-all), DR/SP. Sources: 22_binarisation_summary,
    23_pool_summary, 24_cv_summary, 25_pool_summary, 26_cv_summary."""
    d22 = _collect(codes, "22_binarisation_summary.csv")
    if not d22.empty:
        d22 = d22[d22["category"].isin(["DR", "SP"])]
    p23 = _collect(codes, "23_pool_summary.csv")
    c24 = _collect(codes, "24_cv_summary.csv")
    p25 = _collect(codes, "25_pool_summary.csv")
    c26 = _collect(codes, "26_cv_summary.csv")
    if p23.empty or p25.empty:
        print("  [general] missing pool summaries - skipping activity_ratio")
        return
    before, after = pc.lighten_named("periwinkle", 0.55), pc.named("periwinkle")
    fig, axs = stylia.create_figure(2, 2, height=0.7)
    plot_stage_ratio(axs.next(), "DR", d22, p23, c24, p25, c26, "A")
    plot_stage_ratio(axs.next(), "SP", d22, p23, c24, p25, c26, "B")
    plot_ratio_beforeafter(axs.next(), p25, c26, codes, "DR", "C", before, after)
    plot_ratio_beforeafter(axs.next(), p25, c26, codes, "SP", "D", before, after)
    figure_legend(fig, [("seed pool (before merge)", before),
                        ("final pool (grown + catch-all)", after)],
                  loc="upper center", anchor=(0.5, -0.01), ncol=2)
    stylia.save_figure(os.path.join(out_dir, "27_activity_ratio.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_activity_ratio.png')}")


STAGE_RATIO_SPECS = [   # (stage label, filename, ratio column, optional row filter)
    ("22", "22_binarisation_summary.csv", "active_ratio", None),
    ("23", "23_pool_summary.csv", "active_ratio", None),
    ("24", "24_cv_summary.csv", "active_ratio", lambda d: d[d["reason"] == "modelled"]),
    ("25", "25_pool_summary.csv", "ar_after", None),
    ("26", "26_cv_summary.csv", "active_ratio", None),
]


def _stage_ratio_long(codes):
    """Long-form active ratio per (pathogen, category, stage) across the 22->26 funnel. Reads each
    step's per-pathogen file and tags the pathogen (22_binarisation_summary has no pathogen column,
    so it is tagged on read). One row per dataset/pool."""
    rows = []
    for code in codes:
        for stage, fname, col, filt in STAGE_RATIO_SPECS:
            p = os.path.join(STAGE4PATH, code, fname)
            if not os.path.exists(p):
                continue
            df = pd.read_csv(p)
            if filt is not None:
                df = filt(df)
            if "category" not in df.columns or col not in df.columns:
                continue
            for cat in ("DR", "SP"):
                for v in df.loc[df["category"] == cat, col].dropna():
                    rows.append((code, cat, stage, float(v)))
    return pd.DataFrame(rows, columns=["pathogen", "category", "stage", "ratio"])


def make_activity_ratio_per_pathogen(codes, out_dir):
    """Small multiples (one panel per pathogen) of the active-ratio data flow across the 22->26
    stages, DR vs SP grouped boxes per stage. Per-pathogen breakdown of make_activity_ratio's A/B."""
    long = _stage_ratio_long(codes)
    if long.empty:
        print("  [general] no stage ratios found - skipping activity_ratio_per_pathogen")
        return
    stages = [s for s, _, _, _ in STAGE_RATIO_SPECS]
    drsp = (pc.color("DR"), pc.color("SP"))
    ncols = 3
    nrows = int(np.ceil(len(codes) / ncols))
    fig, axs = stylia.create_figure(nrows, ncols, height=0.28 * nrows)
    flat = np.array(axs.axs).ravel()
    for j, (ax, code) in enumerate(zip(flat, codes)):
        sub = long[long["pathogen"] == code]
        for cat, col, dx in [("DR", drsp[0], -0.18), ("SP", drsp[1], 0.18)]:
            for i, st in enumerate(stages):
                v = sub.loc[(sub["category"] == cat) & (sub["stage"] == st), "ratio"].to_numpy()
                if not len(v):
                    continue
                bp = ax.boxplot([v], positions=[i + dx], widths=0.3, showfliers=False, patch_artist=True)
                bp["boxes"][0].set_facecolor(col)
                bp["medians"][0].set_color("black")
        ax.axhline(0.5, ls="--", lw=0.6, color=pc.color("neutral"))
        ax.set_ylim(-0.02, 1.02); ax.set_xlim(-0.5, len(stages) - 0.5)
        ax.set_xticks(range(len(stages))); ax.set_xticklabels(stages)
        stylia.label(ax, xlabel="", ylabel="active ratio" if j % ncols == 0 else "", title=code)
    for ax in flat[len(codes):]:
        ax.axis("off")
    figure_legend(fig, [("Dose Response", drsp[0]), ("Single Point", drsp[1])],
                  loc="upper center", anchor=(0.5, -0.01), ncol=2)
    stylia.save_figure(os.path.join(out_dir, "27_activity_ratio_per_pathogen.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_activity_ratio_per_pathogen.png')}")


# ============================================================================
# LOW-DATA catch-all overview (moved here from the low-data step, step 26)
# ============================================================================

CATCHALL_MIN_AUROC = 0.70   # matches step 26; evaluation reference line only (NOT enforced)


def make_lowdata_overview(codes, out_dir):
    """Catch-all CV AUROC per deferred category across pathogens (from step-26 26_cv_summary.csv).
    Restricted to the passed `codes` so a subset invocation only aggregates those pathogens."""
    rows = []
    for c in codes:
        f = os.path.join(STAGE4PATH, c, "26_cv_summary.csv")
        if os.path.exists(f):
            rows += pd.read_csv(f).to_dict("records")
    df = pd.DataFrame(rows)
    if df.empty or "auroc" not in df.columns:
        print("  [general] no 26_cv_summary.csv found - skipping low-data overview")
        return
    df = df[df.auroc.notna()].copy()
    if df.empty:
        return
    df["label"] = df.pathogen + " " + df.category
    df = df.sort_values("auroc", ascending=False)
    fig, axs = stylia.create_figure(1, 1)
    ax = axs.next()
    colors = [pc.named("turquoise") if a >= CATCHALL_MIN_AUROC else pc.named("amber")
              for a in df.auroc]
    ax.bar(np.arange(len(df)), df.auroc.to_numpy(), color=colors)
    ax.axhline(CATCHALL_MIN_AUROC, color=pc.color("neutral"), linestyle="--")
    ax.set_xticks(np.arange(len(df))); ax.set_xticklabels(df.label, rotation=90)
    ax.set_ylim(0.5, 1.0)
    stylia.label(ax, xlabel="", ylabel="catch-all CV AUROC",
                 title="Low-data catch-all pools (all kept; teal ≥0.70 / amber <0.70, reference not enforced)")
    stylia.save_figure(os.path.join(out_dir, "27_lowdata_overview.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_lowdata_overview.png')}")


_FATE = [("modelled", "turquoise"), ("catch-all", "periwinkle"), ("uncovered", "silver")]


def make_pipeline_summary(master, codes, out_dir):
    """Cross-step summary from the master table (2x3): the fold->pool->merge consolidation
    funnel, overall coverage, final-pool AUROC and active/inactive ratio, and per-pathogen
    dataset fate + final-pool counts. Shared fate + DR/SP legends sit below the panels.

    Encodes the two non-enforced quality thresholds: (i) step-25's MIN_ASSIGN_BALACC deferrals fall to
    final_status 'uncovered', so panels B (coverage) and E (per-pathogen fate) already carry their effect;
    (ii) the 0.70 AUROC in panel C is drawn as an EVALUATION reference, not a cutoff - modelled and
    catch-all final pools (coloured by fate) may sit below it since steps 24/26 keep on size + CV-ability."""
    abbr = {c: c for c in codes}   # full pathogen codes on the axis (no abbreviation)
    m = master
    order = sorted(codes, key=lambda c: (m[m.pathogen == c].final_status != "uncovered").mean(),
                   reverse=True)
    fcol = {k: pc.named(v) for k, v in _FATE}
    drsp = (pc.color("DR"), pc.color("SP"))
    rng = np.random.default_rng(RANDOM_SEED)
    pools = (m[m.final_status != "uncovered"].dropna(subset=["final_pool", "final_auroc"])
             .drop_duplicates(["pathogen", "category", "final_pool"]))

    fig, axs = stylia.create_figure(2, 3, height=1.1)

    # A. consolidation funnel: datasets -> folded nodes -> step-23 pools -> final pools
    ax = axs.next()
    n_nodes = len(_collect(codes, "23_first_pass.csv"))
    n_pool23 = m.dropna(subset=["pool_id_23"]).drop_duplicates(["pathogen", "pool_id_23"]).shape[0]
    stages = [("Datasets", len(m)), ("Folded\nnodes", n_nodes),
              ("Pools\n(step 23)", n_pool23), ("Final\npools", len(pools))]
    xf = np.arange(len(stages)); vals = [v for _, v in stages]
    ax.bar(xf, vals, color=pc.color("total"))
    ax.set_yscale("log")
    for xi, v in zip(xf, vals):
        ax.text(xi, v, f"{v:,}", ha="center", va="bottom")
    ax.set_xticks(xf); ax.set_xticklabels([s for s, _ in stages])
    compact_yaxis(ax)
    stylia.label(ax, xlabel="", ylabel="count (log)", title="Consolidation: fold → pool → merge", abc="A")

    # B. overall coverage by datasets and by molecules
    ax = axs.next()
    metrics = [("datasets", np.ones(len(m))), ("molecules", m["ds_n_molecules"].to_numpy(float))]
    xb = np.arange(len(metrics)); bottoms = np.zeros(len(metrics))
    for k, _c in _FATE:
        v = np.array([w[m["final_status"].values == k].sum() / w.sum() for _, w in metrics])
        ax.bar(xb, v, bottom=bottoms, color=fcol[k])
        bottoms += v
    ax.set_xticks(xb); ax.set_xticklabels([m0 for m0, _ in metrics]); ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel="fraction", title="Overall coverage", abc="B")

    # C. final-pool AUROC, DR vs SP (x), coloured by fate (modelled vs catch-all). The 0.70 line is an
    #    EVALUATION reference, NOT a cutoff: steps 24 & 26 keep pools on size + CV-ability and only
    #    report AUROC, so final pools - catch-alls especially - may legitimately sit below it.
    ax = axs.next()
    for i, cat in enumerate(["DR", "SP"]):
        v = pools.loc[pools.category == cat, "final_auroc"].to_numpy(float); v = v[np.isfinite(v)]
        if len(v):
            ax.boxplot([v], positions=[i], vert=True, widths=0.5, showfliers=False,
                       boxprops={"color": pc.color("ink")}, medianprops={"color": pc.color("ink")},
                       whiskerprops={"color": pc.color("ink")}, capprops={"color": pc.color("ink")})
        for status in ("modelled", "catch-all"):
            a = pools.loc[(pools.category == cat) & (pools.final_status == status), "final_auroc"].to_numpy(float)
            a = a[np.isfinite(a)]
            if len(a):
                ax.scatter(i + rng.uniform(-0.15, 0.15, size=len(a)), a, s=9, alpha=0.55,
                           color=fcol[status], edgecolors="none")
    ax.axhline(0.7, ls="--", lw=0.8, color=pc.color("neutral"))
    ax.text(1.48, 0.705, "0.70 reference", ha="right", va="bottom",
            fontsize=stylia.FONTSIZE_SMALL, color=pc.color("neutral"))
    ax.set_xticks([0, 1]); ax.set_xticklabels(["DR", "SP"]); ax.set_ylim(0.4, 1.02)
    stylia.label(ax, xlabel="", ylabel="final pool AUROC",
                 title="Final-pool AUROC (not gated)", abc="C")

    # D. final-pool active ratio, DR vs SP (one point per final pool)
    ax = axs.next()
    for i, cat in enumerate(["DR", "SP"]):
        v = pools.loc[pools.category == cat, "final_active_ratio"].to_numpy(float)
        v = v[np.isfinite(v)]
        if len(v):
            ax.boxplot([v], positions=[i], vert=True, widths=0.5, showfliers=False,
                       boxprops={"color": pc.color("ink")}, medianprops={"color": pc.color("ink")},
                       whiskerprops={"color": pc.color("ink")}, capprops={"color": pc.color("ink")})
            ax.scatter(i + rng.uniform(-0.15, 0.15, size=len(v)), v, s=8, alpha=0.4,
                       color=drsp[i], edgecolors="none")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["DR", "SP"]); ax.set_ylim(-0.02, 1.02)
    stylia.label(ax, xlabel="", ylabel="active ratio", title="Final pool active ratio", abc="D")

    # E. per-pathogen dataset fate (fraction)
    ax = axs.next()
    y = np.arange(len(order)); lefts = np.zeros(len(order))
    for k, _c in _FATE:
        fr = np.array([(m[m.pathogen == c].final_status == k).mean() for c in order])
        ax.barh(y, fr, left=lefts, color=fcol[k]); lefts += fr
    ax.set_yticks(y); ax.set_yticklabels([abbr[c] for c in order]); ax.invert_yaxis()
    stylia.label(ax, xlabel="fraction of datasets", ylabel="", title="Dataset fate by pathogen", abc="E")

    # F. final-pool counts per pathogen, DR vs SP
    ax = axs.next()
    yy = np.arange(len(order)); h = 0.4
    for dy, cat, col in [(-h / 2, "DR", drsp[0]), (h / 2, "SP", drsp[1])]:
        n = [pools[(pools.pathogen == c) & (pools.category == cat)].final_pool.nunique() for c in order]
        ax.barh(yy + dy, n, height=h, color=col)
    ax.set_yticks(yy); ax.set_yticklabels([abbr[c] for c in order]); ax.invert_yaxis()
    stylia.label(ax, xlabel="final pools", ylabel="", title="Final pools per pathogen", abc="F")

    figure_legend(fig, [(k, fcol[k]) for k, _ in _FATE], loc="upper center",
                  anchor=(0.3, -0.01), ncol=3)
    figure_legend(fig, [("Dose Response", drsp[0]), ("Single Point", drsp[1])],
                  loc="upper center", anchor=(0.75, -0.01), ncol=2)
    stylia.save_figure(os.path.join(out_dir, "27_pipeline_summary.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_pipeline_summary.png')}")


if __name__ == "__main__":
    # Publication style (also applied at plot_utils import; re-assert so it can't be lost).
    stylia.set_format("print")
    stylia.set_style("article")
    codes = sys.argv[1:] if len(sys.argv) > 1 else available_pathogens()
    if len(codes) < 2:
        print("Need >=2 pathogens for cross-pathogen plots; found:", codes); sys.exit(0)
    out_dir = os.path.join(STAGE4PATH, "general_plots"); os.makedirs(out_dir, exist_ok=True)
    # pooling backbone: dataset -> pool -> AUROC master table (drives the pooling figures)
    master = build_master_table(codes)
    master.to_csv(os.path.join(out_dir, "27_master_table.csv"), index=False)
    print(f"  master table -> {os.path.join(out_dir, '27_master_table.csv')} ({len(master)} datasets)")
    cur = [curation_stats(c) for c in codes]
    binr = [binarisation_stats(c) for c in codes]
    make_curation_comparison(cur, out_dir)
    make_binarisation_comparison(binr, codes, out_dir)
    make_cutoff_sensitivity(codes, out_dir)  # active ratio at low/middle/high per pathogen (diagnostic)
    make_general_plots(codes, out_dir)       # step-23 chemspace partition + active ratios
    make_pool_overview(codes, out_dir)       # condensed pool size + CV AUROC + CV outcome (23+24)
    make_merge_overview(codes, out_dir)      # step-25 merge: AUROC lift, leftover fate, coverage
    make_final_data_overlap(codes, out_dir)  # overlap + label discrepancy between final pools
    make_activity_ratio(codes, out_dir)      # active ratio through the pipeline (22->26)
    make_activity_ratio_per_pathogen(codes, out_dir)  # per-pathogen 22->26 ratio flow (small multiples)
    make_lowdata_overview(codes, out_dir)    # low-data catch-all overview (from step 26 outputs)
    make_pipeline_summary(master, codes, out_dir)   # cross-step summary (from master table)
