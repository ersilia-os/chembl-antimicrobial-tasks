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
    27_eos1klk_all_pathogens.png      eos1klk embedding of retained molecules, coloured by pathogen
    27_eos1klk_tmap.png               t-MAP small multiples (retained molecules), one panel per pathogen
    27_chemspace_partition.png        step-23 molecule-coverage split (multi/singleton/not-considered)
    27_pool_active_ratios.png         step-23 active ratio of multi-dataset pools across pathogens
    27_pool_overview.png              condensed per-pathogen pools (3x2, DR/SP): step-23 pool size,
                                      step-24 CV AUROC, step-24 CV outcome counts
    27_merge_overview.png             step-25 merge per pathogen (3x2, DR/SP): baseline→grown AUROC,
                                      leftover outcome mix, molecules absorbed
    27_lowdata_overview.png           catch-all AUROC per deferred category (from step 26)

Usage:
    python scripts/27_general_plots.py                          # all processed pathogens
    python scripts/27_general_plots.py saureus mtuberculosis    # an explicit subset
"""

import os
import sys
import json
import glob

import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

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
    STEP_LABELS, STEP_ORDER, CATEGORY_ORDER,
    SIZE_LABELS,
    CAT_DISPLAY, COMPOSITION, HIT_LABELS, COMPAT_LABELS, LOC_LABELS,
)
import plot_colors as pc
from plot_base import figure_legend, inset_legend, compact_yaxis
from default import STAGE4PATH, STAGE2PATH, DATABASE_NAME
from pathogen_utils import load_pathogen
from dataset_utils import UNIT_SLUG


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


def plot_chembl_coverage(ax, related, total, abc):
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
                 title=f"Chem-space related to the 15 pathogens\n"
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
        xlabel, title = "Fraction of molecule instances", "Molecule instances by dataset size"
    else:
        xlabel, title = "Molecule instances", "Molecule instances by dataset size"
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
    plot_chembl_coverage(axs.next(), related, total_chembl, abc=None)
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
    plot_discard(axs.next(), codes, results, "step_mols", step_colors, "B", True, "molecule instances", show_y=False)
    plot_discard(axs.next(), codes, results, "step_counts", step_colors, "C", False, "assays", show_y=False)
    plot_discard(axs.next(), codes, results, "step_mols", step_colors, "D", False, "molecule instances", show_y=False)
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


def plot_ratio_vs_size_pooled(ax, results, colors, abc):
    """Active ratio vs dataset size, pooled across all pathogens: faint per-dataset
    scatter + a binned-mean trend per category (log x)."""
    for cat, color, skey, rkey in [("DR", colors[0], "dr_sizes", "dr_ratios"),
                                    ("SP", colors[1], "sp_sizes", "sp_ratios")]:
        sizes = np.concatenate([np.asarray(r[skey], dtype=float) for r in results]) \
            if results else np.array([])
        ratios = np.concatenate([np.asarray(r[rkey], dtype=float) for r in results]) \
            if results else np.array([])
        m = (sizes > 0) & np.isfinite(ratios)
        sizes, ratios = sizes[m], ratios[m]
        if not len(sizes):
            continue
        ax.scatter(sizes, ratios, color=color, alpha=0.08, edgecolors="none", s=6)
        bins = np.logspace(np.log10(sizes.min()), np.log10(sizes.max()), 12)
        centers = np.sqrt(bins[:-1] * bins[1:])
        idx = np.digitize(sizes, bins)
        means = np.array([ratios[idx == k].mean() if (idx == k).any() else np.nan
                          for k in range(1, len(bins))])
        mm = np.isfinite(means)
        ax.plot(centers[mm], means[mm], color=color, marker="o", label=CAT_DISPLAY[cat])
    ax.set_xscale("log")
    ax.set_ylim(-0.02, 1.02)
    stylia.label(ax, xlabel="Molecules / dataset", ylabel="Active ratio",
                 title="Active ratio vs dataset size", abc=abc)


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

    # 1a. ACTIVE RATIO & composition (1x3 + legends below): active ratio (DR+SP box),
    #     class composition, active ratio vs dataset size.
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
        if r["catchall_kept"] is True or r["catchall_kept"] == True:
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


def plot_pool_counts(ax, summ, codes, colors, abbr, category, abc):
    """Stacked bar per pathogen: number of multi-dataset pools vs singletons."""
    sub = summ[summ["category"] == category]
    x = np.arange(len(codes))
    n_multi = [int((~sub[sub["pathogen"] == c]["is_singleton"]).sum()) for c in codes]
    n_single = [int(sub[sub["pathogen"] == c]["is_singleton"].sum()) for c in codes]
    ax.bar(x, n_multi, color=pc.color("datasets"), label="multi-dataset")
    ax.bar(x, n_single, bottom=n_multi, color=pc.named("amber"), label="singleton")
    ax.set_xticks(x); ax.set_xticklabels([abbr[c] for c in codes], rotation=90)
    ax.legend()
    stylia.label(ax, xlabel="", ylabel="pools", title=f"Pool count - {cat_label(category)}", abc=abc)


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
    """Stacked bar per pathogen: molecule fraction in multi / singleton / not-considered."""
    sub = part[part["category"] == category].set_index("pathogen")
    cc = [c for c in codes if c in sub.index]
    x = np.arange(len(cc))
    tot = sub.loc[cc, "n_total"].to_numpy().astype(float)
    tot[tot == 0] = np.nan
    f_multi = sub.loc[cc, "n_multi"].to_numpy() / tot
    f_single = sub.loc[cc, "n_singleton"].to_numpy() / tot
    f_not = sub.loc[cc, "n_not_considered"].to_numpy() / tot
    ax.bar(x, f_multi, color=pc.color("datasets"), label="multi-dataset pool")
    ax.bar(x, f_single, bottom=f_multi, color=pc.named("amber"), label="singleton")
    ax.bar(x, f_not, bottom=f_multi + f_single, color=pc.color("neutral"), label="not considered")
    ax.set_xticks(x); ax.set_xticklabels([abbr[c] for c in cc], rotation=90)
    ax.set_ylim(0, 1)
    ax.legend()
    stylia.label(ax, xlabel="", ylabel="fraction of molecules",
                 title=f"Step-23 chemical-space coverage - {cat_label(category)}", abc=abc)


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
    """One dot per MULTI-dataset pool (y = active ratio, area proportional to n_molecules) with a
    median line per pathogen, so pool size and activity are read together."""
    sub = summ[(~summ["is_singleton"]) & (summ["category"] == category)]
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
                 title=f"Step-23 multi-pool active ratio - {cat_label(category)}", abc=abc)


def make_general_plots(codes, out_dir):
    """Per-pathogen step-23 pool figures (chemical-space partition + multi-pool active ratio),
    aggregated from each pathogen's step-23 outputs. Pool size / modellability now live in the
    condensed 27_pool_overview (make_pool_overview)."""
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
        stylia.save_figure(os.path.join(out_dir, "27_chemspace_partition.png"))

    nmax = max(1.0, float(summ.loc[~summ["is_singleton"], "n_molecules"].max()))
    fig, axs = stylia.create_figure(1, 2)
    plot_active_ratios(axs.next(), summ, codes, colors, abbr, "DR", "A", nmax)
    plot_active_ratios(axs.next(), summ, codes, colors, abbr, "SP", "B", nmax)
    handles, labels = _size_legend_handles(nmax)
    fig.legend(handles, labels, title="molecules / pool", loc="center left",
               bbox_to_anchor=(1.0, 0.5), labelspacing=1.4, borderpad=0.8)
    stylia.save_figure(os.path.join(out_dir, "27_pool_active_ratios.png"))
    print("  [general] wrote 27_chemspace_partition.png, 27_pool_active_ratios.png")


# step-24 CV outcome bins (reason value -> label, colour). "singleton (too small)" pools are
# below the min-molecule floor at step 24 and are the ones rescued downstream by steps 25/26.
CV_OUTCOME = [("modelled", "modelled", "turquoise"),
              ("singleton (too small)", "singleton_below_min_molecules", "silver"),
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
    singleton-too-small / not-CV-able). normalize=False -> absolute pool counts (bar height =
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
    if normalize:
        ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel=ylabel,
                 title=f"Step-24 CV outcome of step-23 pools - {cat_label(category)}", abc=abc)


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
    plot_cv_outcome(axs.next(), cv, codes, "DR", "E", normalize=False)
    plot_cv_outcome(axs.next(), cv, codes, "SP", "F", normalize=False)
    figure_legend(fig, [(lab, pc.named(cname)) for lab, _key, cname in CV_OUTCOME],
                  loc="upper center", anchor=(0.5, -0.01), ncol=3)
    stylia.save_figure(os.path.join(out_dir, "27_pool_overview.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_pool_overview.png')}")


# ============================================================================
# MERGE overview (step 25 — leftover merging into the modelled pools)
# ============================================================================

# 25_merge_log `outcome` value -> (label, colour)
MERGE_OUTCOME = [("merged", "merged", "turquoise"),
                 ("redundant", "redundant", "silver"),
                 ("deferred (AUROC drop)", "deferred_auroc_drop", "amber")]


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
# LOW-DATA catch-all overview (moved here from the low-data step, step 26)
# ============================================================================

CATCHALL_MIN_AUROC = 0.70   # matches step 26; the kept/dropped threshold line


def make_lowdata_overview(out_dir):
    """Catch-all CV AUROC per deferred category across pathogens (from step-26 26_cv_summary.csv)."""
    rows = []
    for f in sorted(glob.glob(os.path.join(STAGE4PATH, "*", "26_cv_summary.csv"))):
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
    colors = [pc.named("turquoise") if k else pc.color("neutral") for k in df.kept]
    ax.bar(np.arange(len(df)), df.auroc.to_numpy(), color=colors)
    ax.axhline(CATCHALL_MIN_AUROC, color=pc.named("amber"), linestyle="--")
    ax.set_xticks(np.arange(len(df))); ax.set_xticklabels(df.label, rotation=90)
    ax.set_ylim(0.5, 1.0)
    stylia.label(ax, xlabel="", ylabel="catch-all CV AUROC",
                 title="Low-data catch-all pools (teal = kept, grey = dropped)")
    stylia.save_figure(os.path.join(out_dir, "27_lowdata_overview.png"))
    print(f"  figure -> {os.path.join(out_dir, '27_lowdata_overview.png')}")


_FATE = [("modelled", "turquoise"), ("catch-all", "periwinkle"), ("uncovered", "silver")]


def make_pipeline_summary(master, codes, out_dir):
    """Cross-step summary from the master table (2x3): the fold->pool->merge consolidation
    funnel, overall coverage, final-pool AUROC and active/inactive ratio, and per-pathogen
    dataset fate + final-pool counts. Shared fate + DR/SP legends sit below the panels."""
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

    # C. final-pool AUROC, DR vs SP
    ax = axs.next()
    for i, cat in enumerate(["DR", "SP"]):
        v = pools.loc[pools.category == cat, "final_auroc"].to_numpy(float); v = v[np.isfinite(v)]
        if len(v):
            ax.boxplot([v], positions=[i], vert=True, widths=0.5, showfliers=False,
                       boxprops={"color": pc.color("ink")}, medianprops={"color": pc.color("ink")},
                       whiskerprops={"color": pc.color("ink")}, capprops={"color": pc.color("ink")})
            ax.scatter(i + rng.uniform(-0.15, 0.15, size=len(v)), v, s=8, alpha=0.4,
                       color=drsp[i], edgecolors="none")
    ax.axhline(0.7, ls="--", lw=0.8, color=pc.color("neutral"))
    ax.set_xticks([0, 1]); ax.set_xticklabels(["DR", "SP"]); ax.set_ylim(0.4, 1.02)
    stylia.label(ax, xlabel="", ylabel="final pool AUROC", title="Modellability (final pools)", abc="C")

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
    make_general_plots(codes, out_dir)       # step-23 chemspace partition + active ratios
    make_pool_overview(codes, out_dir)       # condensed pool size + CV AUROC + CV outcome (23+24)
    make_merge_overview(codes, out_dir)      # step-25 merge: AUROC lift, leftover fate, coverage
    make_lowdata_overview(out_dir)           # low-data catch-all overview (from step 26 outputs)
    make_pipeline_summary(master, codes, out_dir)   # cross-step summary (from master table)
