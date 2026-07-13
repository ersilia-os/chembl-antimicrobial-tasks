"""
Step 23 (pooling rebuild) - model-transfer signal and pool building.

Uses a model-transfer signal (train a RandomForest on dataset A, score its AUROC on
dataset B's disjoint molecules) to build actual pooled datasets. For each
pathogen and category (DR/SP, kept separate), it:

  1. TRANSFER SIGNAL: trains a RandomForest on each eligible dataset A and scores its AUROC
     on every other dataset B's molecules that are NOT shared with A (a disjoint
     generalisation test). "A predicts B well" is the evidence for pooling A and B.
     Directional; pairs with too few disjoint test molecules are left unscored (NaN).

  2. POOLS (the deliverable): the transfer graph at the CHOSEN setting
     (bidirectional, AUROC >= POOL_THRESHOLD) is clustered with GREEDY MODULARITY into a
     few large pools + a singleton tail. Each pool is merged to ONE row per molecule by
     the value-resolution policy below, so no molecule is duplicated inside a pool.

Value-resolution policy (per molecule, inside a pool):
  - average the raw values within the same (unit, direction) and re-binarise against the
    category cutoff (DR <= 10 uM; SP % >= 50; SP mm >= 20). DR averages across activity
    types (MIC/IC50/EC50, all uM) - signed off. Averaging dissolves ties by construction.
  - cross-unit SP (a molecule with both % and mm): "% prevails" (average the %, drop mm).

Scope: pools the ELIGIBLE datasets (non-degenerate, >= MIN_DATASET_SIZE molecules). The
sub-threshold / degenerate datasets stay out - folding them in by chemistry is a later step.

Chosen settings (bidirectional, POOL_THRESHOLD=0.75, greedy modularity, DR cross-type
averaging) were selected empirically on the dev organisms. MIN_DATASET_SIZE and MIN_TEST
are the method's eligibility/disjoint-test floors. Dev-first: defaults to saureus,
campylobacter, mtuberculosis.

Inputs:
    config/category_cutoffs.csv                          (cutoff per unit, via CONFIGPATH)
  output/stage4/<pathogen>/:
    22_binarised/{DR,SP}/<assay>_<type>_<slug>.csv.gz   (inchikey, compound_chembl_id,
                                                          smiles, value, relation, bin)
    22_binarisation_summary.csv                          (unit, direction, cutoff per dataset)

Outputs (per pathogen, output/stage4/<pathogen>/):
    23_transfer/23_transfer_<cat>_<tag>.npz  cached transfer matrices (tag = fold settings, e.g. o90d20)
    23_pools/<cat>/<cat>_####.csv.gz     merged pools (1 row per molecule)
    23_pool_summary.csv                  per-pool sizes / active ratio / n_active_flipped / n_proof_demoted
    23_pool_members.csv                  pool -> member dataset mapping
    23_chemspace_partition.csv           unique-molecule split: multi-pool / singleton / not-considered
    23_first_pass/<cat>/<id>.csv.gz      post-fold "first pass" datasets (1 row per molecule)
    23_first_pass.csv                    first-pass manifest (+ n_active_flipped / n_proof_demoted)
    23_folds.csv                         near-duplicate fold provenance
    23_transfer_heatmap_<cat>.png        pool-ordered transfer AUROC heatmap
    23_pooling_overview.png              2x3 per-pathogen summary: fold-group size hist / pool-size
                                         hist / pool active-ratio boxplot / datasets considered /
                                         chemical space considered / conservative-rebinarisation effect

Cross-pathogen pool figures (pool counts / sizes / active ratios) now live in step 27
(general plots); this step is per-pathogen only.

Usage:
    python scripts/23_transfer_pooling.py                    # 3 dev organisms
    python scripts/23_transfer_pooling.py saureus ecoli      # explicit pathogens
    python scripts/23_transfer_pooling.py --force ecoli      # ignore npz cache
    python scripts/23_transfer_pooling.py --plot-only        # redraw heatmaps + pooling overview from files
"""

import os
import sys
import glob
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import MIN_TEST, MIN_DATASET_SIZE, CONFIGPATH, RANDOM_SEED
from pathogen_utils import load_pathogen
from dataset_utils import UNIT_SLUG, load_category_cutoffs
from pooling_utils import (morgan_matrix, transfer_matrix, greedy_modularity_pools,
                           merge_datasets, load_binarised_datasets, fold_near_duplicates,
                           load_dataset_meta as _load_dataset_meta_shared)
from plot_utils import cat_label, _xticks
import plot_colors as pc
from plot_base import BasePlot, render_figure, axis_legend, compact_yaxis

# Format: print | Style: article - publication house style (shared across stage 4).
stylia.set_format("print")
stylia.set_style("article")

# --- configuration ----------------------------------------------------------
DEV_PATHOGENS = ["saureus", "campylobacter", "mtuberculosis"]
CATEGORIES = ["DR", "SP"]

# MIN_DATASET_SIZE (eligibility floor) is a project-wide constant in src/default.py so step 25
# shares the same value.

# CHOSEN pooling settings (selected empirically on the dev organisms)
POOL_RULE = "bidirectional"
POOL_THRESHOLD = 0.75

# Near-duplicate folding (pre-transfer). An edge joins two datasets with >= NEAR_DUP_OVERLAP
# containment (shared / min size) and < NEAR_DUP_DISCORD shared-label discordance; datasets are
# then folded by CLIQUE COVER (every pair inside a fold group satisfies the edge test directly),
# NOT connected components - so a bridging dataset can no longer transitively fold two datasets
# that fail the test between themselves. Runs over ALL datasets (incl. degenerate/small): transfer
# is blind to near-duplicates (no disjoint test set), so without folding they fragment into
# singletons. NEAR_DUP_MAXPOST caps candidate generation only.
NEAR_DUP_OVERLAP = 0.90
NEAR_DUP_DISCORD = 0.20
NEAR_DUP_MAXPOST = 100

output_base = os.path.join(root, "..", "output")


# --- data loading -----------------------------------------------------------

def load_category_datasets(pathogen_code, category):
    """Load all binarised datasets for one pathogen/category (see pooling_utils)."""
    bdir = os.path.join(output_base, "stage4", pathogen_code, "22_binarised", category)
    return load_binarised_datasets(bdir)


def load_dataset_meta(pathogen_code, category, cutoffs):
    """Map dataset id -> (unit, direction, cutoff) from the binarisation summary."""
    summary = os.path.join(output_base, "stage4", pathogen_code, "22_binarisation_summary.csv")
    return _load_dataset_meta_shared(summary, category, cutoffs, UNIT_SLUG)


def filter_eligible(datasets):
    """Keep non-degenerate datasets (both classes) with >= MIN_DATASET_SIZE molecules."""
    eligible, n_deg, n_small = [], 0, 0
    for d in datasets:
        if d["labels"].min() == d["labels"].max():
            n_deg += 1
        elif len(d["ikset"]) < MIN_DATASET_SIZE:
            n_small += 1
        else:
            eligible.append(d)
    return eligible, n_deg, n_small


# --- transfer matrix (compute or load cache) --------------------------------

def find_cached_transfer(pathogen_code, category):
    """Return the category's transfer npz on disk, or None if there is none.

    Picks WHATEVER npz is present for the category with a deterministic, tag-blind glob
    (sorted; first match if several). The cache is reused as-is with no version/tag/hash
    staleness guard: delete the file (or pass --force) to force recomputation.
    """
    cache_dir = os.path.join(output_base, "stage4", pathogen_code, "23_transfer")
    matches = sorted(glob.glob(os.path.join(cache_dir, f"23_transfer_{category}_*.npz")))
    return matches[0] if matches else None


def compute_or_load_transfer(pathogen_code, category, eligible, fps, force):
    """Reuse the cached transfer matrix if any npz exists for the category; else compute and cache it.

    Read side: take whatever npz is on disk for the category (see find_cached_transfer) - no
    tag/version/hash match. The cache is reused as-is; delete it (or pass --force) to recompute.
    Write side: still writes a fold-tagged filename (o<overlap>d<discord>).
    """
    cache_dir = os.path.join(output_base, "stage4", pathogen_code, "23_transfer")
    os.makedirs(cache_dir, exist_ok=True)

    cached = find_cached_transfer(pathogen_code, category) if not force else None
    if cached is not None:
        z = np.load(cached, allow_pickle=False)
        ids = z["ids"].tolist()
        print(f"    [cache] loaded {cached} ({len(ids)} datasets)")
        return ids, z["AUROC"], z["F1"], z["N"]

    tag = f"o{int(NEAR_DUP_OVERLAP * 100)}d{int(NEAR_DUP_DISCORD * 100)}"
    cache = os.path.join(cache_dir, f"23_transfer_{category}_{tag}.npz")
    triples = [(d["id"], d["iks"], d["labels"]) for d in eligible]
    ids, AUROC, F1, N = transfer_matrix(triples, fps, min_test=MIN_TEST)
    np.savez_compressed(cache, ids=np.array(ids), AUROC=AUROC, F1=F1, N=N)
    print(f"    [cache] wrote {cache}")
    return ids, AUROC, F1, N


# --- pooling (the deliverable) ----------------------------------------------

def merge_pool(members, eligible, meta, ik2info, orig_by_id):
    """Merge a pool to ONE row per molecule, expanding folded nodes back to their ORIGINAL datasets.

    Returns (df, diag): diag = merge transparency counters (n_flipped, n_proof_demoted)."""
    orig = [orig_by_id[oid] for i in members for oid in eligible[i]["members"]]
    return merge_datasets(orig, meta, ik2info, report_diag=True)


def build_and_write_pools(pathogen_code, category, ids, labels, eligible, meta, ik2info, orig_by_id):
    """Cluster (already done -> labels), merge each pool, write files + summary rows.

    `eligible` are folded nodes; each carries `members` (original dataset ids). Pools are merged
    and reported at the ORIGINAL-dataset level. `is_singleton` reflects a single NODE (which may
    itself be a fold of near-duplicates); `n_datasets` counts expanded original datasets.
    """
    by_idx = {d["id"]: i for i, d in enumerate(eligible)}
    comms = defaultdict(list)
    for i, lab in enumerate(labels):
        comms[lab].append(by_idx[ids[i]])
    ranked = sorted(comms.values(),
                    key=lambda mem: len(set().union(*[eligible[i]["ikset"] for i in mem])),
                    reverse=True)

    out_dir = os.path.join(output_base, "stage4", pathogen_code, "23_pools", category)
    os.makedirs(out_dir, exist_ok=True)
    for f in glob.glob(os.path.join(out_dir, "*.csv.gz")):
        os.remove(f)

    summary_rows, member_rows = [], []
    for rank, mem in enumerate(ranked):
        pool_id = f"{category}_{rank:04d}"
        df, diag = merge_pool(mem, eligible, meta, ik2info, orig_by_id)
        df.to_csv(os.path.join(out_dir, pool_id + ".csv.gz"), index=False)
        n = len(df); na = int((df["bin"] == 1).sum())
        orig_ids = [oid for i in mem for oid in eligible[i]["members"]]
        summary_rows.append({"pathogen": pathogen_code, "category": category, "pool_id": pool_id,
                             "rank": rank, "n_nodes": len(mem), "n_datasets": len(orig_ids),
                             "n_molecules": n, "n_active": na, "n_inactive": n - na,
                             "active_ratio": round(na / n, 5) if n else np.nan,
                             "n_active_flipped": diag["n_flipped"],
                             "n_proof_demoted": diag["n_proof_demoted"],
                             "is_singleton": len(mem) == 1})
        member_rows += [{"category": category, "pool_id": pool_id, "dataset_id": oid}
                        for oid in orig_ids]
    return summary_rows, member_rows


def chemspace_partition(datasets, eligible, ids, labels):
    """Split the pathogen/category chemical space (unique molecules) into three buckets.

    Denominator is every unique molecule across ALL datasets for the category (eligible
    plus the small/degenerate ones excluded from pooling). Buckets, priority
    multi > singleton > not-considered:
      - multi:      molecules in a multi-dataset (non-singleton) pool
      - singleton:  molecules only in an eligible singleton pool
      - not:        molecules only in datasets too small/degenerate to be pooled
    """
    by_idx = {d["id"]: i for i, d in enumerate(eligible)}
    comms = defaultdict(list)
    for i, lab in enumerate(labels):
        comms[lab].append(by_idx[ids[i]])
    multi_idx = [i for mem in comms.values() if len(mem) > 1 for i in mem]
    single_idx = [i for mem in comms.values() if len(mem) == 1 for i in mem]

    def _union(idxs):
        return set().union(*[eligible[i]["ikset"] for i in idxs]) if idxs else set()

    multi_set = _union(multi_idx)
    single_set = _union(single_idx) - multi_set
    universe = set().union(*[d["ikset"] for d in datasets])
    n_multi, n_single = len(multi_set), len(single_set)
    return {"n_multi": n_multi, "n_singleton": n_single,
            "n_not_considered": len(universe) - n_multi - n_single, "n_total": len(universe)}


# --- per-pathogen plot ------------------------------------------------------

def make_heatmap_figure(pathogen_code, category, AUROC, labels):
    """Pool-ordered transfer-AUROC heatmap: rows/cols reordered by pool assignment so the
    greedy-modularity pool blocks sit on the diagonal. Drawn at any N (no per-cell ticks)."""
    order = np.argsort(labels, kind="stable")
    M = AUROC[np.ix_(order, order)]
    fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
    ax = axs.next()
    cm = pc.fading("total").cmap.copy()
    cm.set_bad("#eeeeee")
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cm, vmin=0, vmax=1, aspect="equal")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xticks([]); ax.set_yticks([])
    stylia.label(ax, xlabel="test dataset B (pool-ordered)", ylabel="trainer dataset A (pool-ordered)",
                 title=f"Transfer AUROC - {cat_label(category)} ({len(labels)} datasets)")
    path = os.path.join(output_base, "stage4", pathogen_code, f"23_transfer_heatmap_{category}.png")
    stylia.save_figure(path)
    print(f"    heatmap -> {path}")


def replot_heatmaps(pathogen_code):
    """Redraw the per-pathogen transfer heatmaps from the cached npz matrices and
    23_pool_members.csv (for the pool ordering), without recomputing the transfer matrix.
    Pool labels are rebuilt from the members table; datasets not in a multi-dataset pool
    (singletons) each get their own block, matching the on-disk pool structure."""
    out = os.path.join(output_base, "stage4", pathogen_code)
    mem_path = os.path.join(out, "23_pool_members.csv")
    if not os.path.exists(mem_path):
        print(f"  [plot-only] {pathogen_code}: no 23_pool_members.csv - skipping heatmaps")
        return
    mem = pd.read_csv(mem_path)
    for category in CATEGORIES:
        # reuse whatever transfer npz is on disk for this category, as-is (tag-blind, no
        # staleness guard); delete the file to force a recompute via a normal run.
        cached = find_cached_transfer(pathogen_code, category)
        if cached is None:
            continue
        z = np.load(cached, allow_pickle=False)
        ids = [str(i) for i in z["ids"]]
        sub = mem[mem.category == category]
        d2p = dict(zip(sub.dataset_id.astype(str), sub.pool_id.astype(str)))
        pool_int = {pid: i for i, pid in enumerate(sorted(set(d2p.values())))}
        nxt, labels = len(pool_int), []
        for i in ids:
            if i in d2p:
                labels.append(pool_int[d2p[i]])
            else:                                # singleton -> its own block
                labels.append(nxt); nxt += 1
        make_heatmap_figure(pathogen_code, category, z["AUROC"], np.array(labels))


def _pooled_codes():
    """Pathogens that have a 23_pool_summary.csv on disk, sorted (default for --plot-only)."""
    return sorted(os.path.basename(os.path.dirname(f))
                  for f in glob.glob(os.path.join(output_base, "stage4", "*", "23_pool_summary.csv")))


# --- per-pathogen pooling overview (2x3) ------------------------------------
# Bubbles (A/B) = one circle per group, area proportional to molecules; DR/SP coloured.
# Coverage panels (D/E) split datasets and chemical space into multi / singleton /
# not-considered. All colours from the shared semantic map (plot_colors); article style.

_COVER = [("multi-dataset pool", "cobalt"), ("singleton", "periwinkle"),
          ("not considered", "silver")]     # shared legend/colours for the two coverage panels


def _headroom(ax, totals, frac=0.35):
    """Leave head-space above the tallest bar so an in-axis legend never covers data."""
    top = float(np.max(totals)) if len(totals) else 0.0
    if top > 0:
        ax.set_ylim(0, top * (1 + frac))


def _bubble_sizes(mols, smin=8.0, smax=220.0):
    """Marker areas scaled linearly with molecule count (area proportional to molecules,
    so radius tracks sqrt) between smin and smax."""
    m = np.asarray(mols, dtype=float)
    if not len(m):
        return m
    mx = float(np.nanmax(m))
    if mx <= 0:
        return np.full(len(m), smin)
    return smin + (smax - smin) * (m / mx)


_FOLD_LABELS = ["2", "3", "4", "5+"]     # near-dup fold-group size buckets


class FoldSizeHist(BasePlot):
    """How many datasets each near-duplicate fold group combines (grouped bar of fold-group
    size buckets 2 / 3 / 4 / 5+, counted per category). Only groups of >= 2 datasets."""

    def __init__(self, ax, folds, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if (folds is not None and len(folds) and (folds.category == c).any())]
        if not cats:
            self._empty(xlabel="datasets per fold group", ylabel="fold groups",
                        title="Near-duplicate folding")
            return
        x = np.arange(len(_FOLD_LABELS))
        w = 0.8 / len(cats)
        rows = {}
        for i, c in enumerate(cats):
            n = folds.loc[folds.category == c, "n_members"].clip(upper=5).astype(int)
            counts = [int((n == 2).sum()), int((n == 3).sum()), int((n == 4).sum()), int((n >= 5).sum())]
            ax.bar(x + (i - (len(cats) - 1) / 2) * w, counts, width=w, color=pc.color(c))
            rows[c] = counts
        ax.set_xticks(x)
        ax.set_xticklabels(_FOLD_LABELS)
        compact_yaxis(ax)
        self._label(xlabel="datasets per fold group", ylabel="fold groups", title="Near-duplicate folding")
        axis_legend(ax, [(cat_label(c), pc.color(c)) for c in cats],
                    loc="upper right", anchor=(0.98, 0.98), frame=True)
        self.data = pd.DataFrame(rows, index=_FOLD_LABELS)


class PoolSizeHist(BasePlot):
    """Distribution of pool sizes (molecules per pool), DR and SP overlaid on log-spaced bins."""

    def __init__(self, ax, ps, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if (ps["category"] == c).any()]
        allm = ps["n_molecules"].to_numpy(dtype=float)
        allm = allm[allm > 0]
        if not cats or not len(allm):
            self._empty(xlabel="molecules per pool", ylabel="pools", title="Pool sizes")
            return
        bins = np.logspace(np.log10(allm.min()), np.log10(allm.max()), 16)
        for c in cats:
            m = ps.loc[ps.category == c, "n_molecules"].to_numpy(dtype=float)
            m = m[m > 0]
            if len(m):
                ax.hist(m, bins=bins, color=pc.color(c), alpha=0.6)
        ax.set_xscale("log")
        self._label(xlabel="molecules per pool", ylabel="pools", title="Pool sizes")
        axis_legend(ax, [(cat_label(c), pc.color(c)) for c in cats],
                    loc="upper right", anchor=(0.98, 0.98), frame=True)
        self.data = ps[["category", "pool_id", "n_molecules"]].copy()


class ActiveRatioBox(BasePlot):
    """Active ratio of the final pools: a boxplot per category with one dot per pool, dot
    area proportional to the pool's molecule count."""

    def __init__(self, ax, ps, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if (ps["category"] == c).any()]
        if not cats:
            self._empty(xlabel="", ylabel="active ratio", title="Pool active ratio")
            return
        rng = np.random.default_rng(RANDOM_SEED)
        data = [ps.loc[ps.category == c, "active_ratio"].dropna().to_numpy(dtype=float) for c in cats]
        pos = list(range(len(cats)))
        bp = ax.boxplot(data, positions=pos, widths=0.5, showfliers=False, patch_artist=True)
        for patch, c in zip(bp["boxes"], cats):
            patch.set_facecolor("none")               # empty boxes so the pool dots show through
            patch.set_edgecolor(pc.color(c))
        for i, c in enumerate(cats):                   # whiskers/caps come in pairs per box
            for artist in bp["whiskers"][2 * i:2 * i + 2] + bp["caps"][2 * i:2 * i + 2]:
                artist.set_color(pc.color(c))
        for med in bp["medians"]:
            med.set_color(pc.color("ink"))
        for i, c in enumerate(cats):
            sub = ps[ps.category == c]
            r = sub["active_ratio"].to_numpy(dtype=float)
            x = i + rng.uniform(-0.18, 0.18, size=len(r))
            ax.scatter(x, r, s=_bubble_sizes(sub["n_molecules"].to_numpy(dtype=float)),
                       color=pc.color(c), alpha=0.5, edgecolors="none")
        ax.set_xticks(pos)
        ax.set_xticklabels([cat_label(c) for c in cats])
        ax.set_ylim(-0.02, 1.02)
        self._label(xlabel="", ylabel="active ratio", title="Pool active ratio")
        self.data = ps[["category", "pool_id", "active_ratio", "n_molecules"]].copy()


class CoverageDatasets(BasePlot):
    """Original datasets by fate per category: in a multi-dataset pool / singleton pool /
    not considered (degenerate or < MIN_DATASET_SIZE after folding)."""

    def __init__(self, ax, fp, ps, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if (fp["category"] == c).any()]
        if not cats:
            self._empty(ylabel="datasets", title="Datasets considered")
            return
        x = np.arange(len(cats))
        raw = np.array([int(fp[fp.category == c]["n_members"].sum()) for c in cats], dtype=float)
        multi = np.array([int(ps[(ps.category == c) & (~ps.is_singleton)]["n_datasets"].sum())
                          for c in cats], dtype=float)
        singl = np.array([int(ps[(ps.category == c) & (ps.is_singleton)]["n_datasets"].sum())
                          for c in cats], dtype=float)
        notc = raw - multi - singl
        bottom = np.zeros(len(cats))
        for vals, (_, ckey) in zip((multi, singl, notc), _COVER):
            ax.bar(x, vals, bottom=bottom, color=pc.named(ckey))
            bottom += vals
        _xticks(ax, [cat_label(c) for c in cats])
        _headroom(ax, raw)
        compact_yaxis(ax)
        self._label(xlabel="", ylabel="datasets", title="Datasets considered")
        axis_legend(ax, [(l, pc.named(ck)) for l, ck in _COVER], loc="upper right",
                    anchor=(0.98, 0.98), frame=True)
        self.data = pd.DataFrame({"category": cats, "multi": multi, "singleton": singl,
                                  "not_considered": notc})


class CoverageChemspace(BasePlot):
    """Unique molecules (chemical space) by pooling fate, as % per category
    (from 23_chemspace_partition)."""

    def __init__(self, ax, cp, abc=None):
        super().__init__(ax, abc)
        if cp is None or not len(cp):
            self._empty(ylabel="% of chemical space", title="Chemical space considered")
            return
        cats = [c for c in CATEGORIES if (cp["category"] == c).any()]
        if not cats:
            self._empty(ylabel="% of chemical space", title="Chemical space considered")
            return
        x = np.arange(len(cats))
        g = cp.set_index("category")
        tot = np.array([max(int(g.loc[c, "n_total"]), 1) for c in cats], dtype=float)
        multi = np.array([int(g.loc[c, "n_multi"]) for c in cats], dtype=float) / tot * 100
        singl = np.array([int(g.loc[c, "n_singleton"]) for c in cats], dtype=float) / tot * 100
        notc = np.array([int(g.loc[c, "n_not_considered"]) for c in cats], dtype=float) / tot * 100
        bottom = np.zeros(len(cats))
        for vals, (_, ckey) in zip((multi, singl, notc), _COVER):
            ax.bar(x, vals, bottom=bottom, color=pc.named(ckey))
            bottom += vals
        _xticks(ax, [cat_label(c) for c in cats])
        ax.set_ylim(0, 100)
        self._label(xlabel="", ylabel="% of chemical space", title="Chemical space considered")
        axis_legend(ax, [(l, pc.named(ck)) for l, ck in _COVER], loc="upper right",
                    anchor=(0.98, 0.98), frame=True)
        self.data = pd.DataFrame({"category": cats, "multi_pct": multi, "singleton_pct": singl,
                                  "not_considered_pct": notc})


class RebinEffectPlot(BasePlot):
    """Molecule-label accounting across the final pools: active vs the two arms of
    conservative rebinarisation (prove-the-hit and averaging). Counts are per-pool
    occurrences (a molecule in two pools is counted twice)."""

    def __init__(self, ax, ps, cp=None, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if (ps["category"] == c).any()]
        if not cats:
            self._empty(ylabel="molecules", title="Conservative rebinarisation"); return
        x = np.arange(len(cats))
        active = np.array([int(ps[ps.category == c]["n_active"].sum()) for c in cats], dtype=float)
        proof = np.array([int(ps[ps.category == c]["n_proof_demoted"].sum()) for c in cats], dtype=float)
        aver = np.array([int(ps[ps.category == c]["n_active_flipped"].sum()) for c in cats], dtype=float)
        ax.bar(x, active, color=pc.color("active"))
        ax.bar(x, proof, bottom=active, color=pc.named("amber"))
        ax.bar(x, aver, bottom=active + proof, color=pc.named("periwinkle"))
        # annotate each bar with the share of pooled unique molecules affected by the two
        # conservative arms (proof-demoted + averaged-out). Denominator = unique pooled
        # molecules (chemspace multi+singleton); falls back to per-pool label total.
        totals = active + proof + aver
        if cp is not None and len(cp):
            g = cp.set_index("category")
            base = np.array([float(g.loc[c, "n_multi"] + g.loc[c, "n_singleton"])
                             if c in g.index else 0.0 for c in cats])
        else:
            base = totals.copy()
        for xi, aff, b, top in zip(x, proof + aver, base, totals):
            if b > 0:
                ax.text(xi, top, f"{100 * aff / b:.1f}%", ha="center", va="bottom",
                        fontsize=stylia.FONTSIZE_SMALL)
        _xticks(ax, [cat_label(c) for c in cats])
        _headroom(ax, totals, frac=0.22)
        compact_yaxis(ax)
        self._label(xlabel="", ylabel="molecules", title="Conservative rebinarisation")
        axis_legend(ax, [("active", pc.color("active")), ("lost: prove-the-hit", pc.named("amber")),
                         ("lost: averaging", pc.named("periwinkle"))],
                    loc="upper right", anchor=(0.98, 0.98), frame=True)
        self.data = pd.DataFrame({"category": cats, "active": active, "proof_demoted": proof,
                                  "averaged_out": aver, "pct_rebinarised": 100 * (proof + aver) / base
                                  if (base > 0).all() else np.nan})


def make_pooling_overview(pathogen_code):
    """Draw the per-pathogen 2x3 pooling overview from the step-23 CSVs (no recomputation):
    (A) near-duplicate folding, (B) pool sizes, (C) pool active ratio, (D) datasets
    considered, (E) chemical space considered, (F) conservative-rebinarisation effect.
    Skipped if the pool summaries are absent."""
    out = os.path.join(output_base, "stage4", pathogen_code)
    fp_path = os.path.join(out, "23_first_pass.csv")
    ps_path = os.path.join(out, "23_pool_summary.csv")
    if not (os.path.exists(fp_path) and os.path.exists(ps_path)):
        print(f"  [overview] {pathogen_code}: missing 23_first_pass/23_pool_summary - skipping")
        return
    fp = pd.read_csv(fp_path)
    ps = pd.read_csv(ps_path)
    folds = pd.read_csv(os.path.join(out, "23_folds.csv")) if os.path.exists(
        os.path.join(out, "23_folds.csv")) else None
    cp = pd.read_csv(os.path.join(out, "23_chemspace_partition.csv")) if os.path.exists(
        os.path.join(out, "23_chemspace_partition.csv")) else None
    # robust bool parse (read_csv may hand back "True"/"False" strings)
    fp["is_merged"] = fp["is_merged"].astype(str).str.lower().isin(["true", "1"])
    ps["is_singleton"] = ps["is_singleton"].astype(str).str.lower().isin(["true", "1"])
    render_figure(out, "23_pooling_overview", 2, 3, [
        lambda ax: FoldSizeHist(ax, folds, abc="A"),
        lambda ax: PoolSizeHist(ax, ps, abc="B"),
        lambda ax: ActiveRatioBox(ax, ps, abc="C"),
        lambda ax: CoverageDatasets(ax, fp, ps, abc="D"),
        lambda ax: CoverageChemspace(ax, cp, abc="E"),
        lambda ax: RebinEffectPlot(ax, ps, cp=cp, abc="F"),
    ], width=0.85, height=0.85)


# --- first-pass (post-fold) datasets ----------------------------------------

def write_first_pass(pathogen_code, category, nodes, meta, ik2info, orig_by_id):
    """Persist the post-fold 'first pass' datasets: one file per node + manifest rows.

    Each node (a near-dup merge or a pass-through dataset) is merged from its ORIGINAL member
    datasets via `merge_datasets` and written to 23_first_pass/<cat>/<id>.csv.gz (same schema as
    the pools). Merged ids are `<largest-member>_merged<N>`; pass-through keep their dataset id.
    """
    out_dir = os.path.join(output_base, "stage4", pathogen_code, "23_first_pass", category)
    os.makedirs(out_dir, exist_ok=True)
    for f in glob.glob(os.path.join(out_dir, "*.csv.gz")):
        os.remove(f)
    rows = []
    for node in nodes:
        df, diag = merge_datasets([orig_by_id[oid] for oid in node["members"]], meta, ik2info,
                                  report_diag=True)
        df.to_csv(os.path.join(out_dir, node["id"] + ".csv.gz"), index=False)
        n = len(df); na = int((df["bin"] == 1).sum())
        rows.append({"pathogen": pathogen_code, "category": category, "first_pass_id": node["id"],
                     "is_merged": len(node["members"]) > 1, "n_members": len(node["members"]),
                     "members": "|".join(node["members"]), "n_molecules": n,
                     "n_active": na, "active_ratio": round(na / n, 5) if n else np.nan,
                     "n_active_flipped": diag["n_flipped"],
                     "n_proof_demoted": diag["n_proof_demoted"]})
    return rows


# --- per-pathogen driver ----------------------------------------------------

def process(pathogen_code, force):
    pathogen = load_pathogen(pathogen_code)
    print(f"Step 23 - {pathogen_code} ({pathogen})")
    out = os.path.join(output_base, "stage4", pathogen_code)
    cutoffs = load_category_cutoffs(CONFIGPATH)

    all_pool_summary, all_pool_members, partition_rows, all_fold_rows, all_first_pass = [], [], [], [], []
    for category in CATEGORIES:
        datasets, ik2info = load_category_datasets(pathogen_code, category)
        meta_all = load_dataset_meta(pathogen_code, category, cutoffs)
        orig_by_id = {d["id"]: d for d in datasets}

        # Fold near-duplicate datasets (ALL datasets, incl. degenerate/small) BEFORE the transfer.
        nodes, fold_rows = fold_near_duplicates(datasets, meta_all, ik2info, NEAR_DUP_OVERLAP,
                                                NEAR_DUP_DISCORD, NEAR_DUP_MAXPOST)
        for r in fold_rows:
            r["pathogen"] = pathogen_code; r["category"] = category
        all_fold_rows += fold_rows
        print(f"  {category}: {len(datasets)} datasets -> {len(nodes)} after folding "
              f"({len(fold_rows)} near-dup groups, o={NEAR_DUP_OVERLAP}/d={NEAR_DUP_DISCORD})")

        # persist the post-fold 'first pass' datasets (the clean starting group feeding transfer)
        all_first_pass += write_first_pass(pathogen_code, category, nodes, meta_all, ik2info, orig_by_id)

        eligible, n_deg, n_small = filter_eligible(nodes)
        print(f"    {len(nodes)} nodes -> {len(eligible)} eligible "
              f"(dropped {n_deg} degenerate, {n_small} < {MIN_DATASET_SIZE} mols)")
        if len(eligible) < 2:
            print(f"    < 2 eligible datasets - skipping {category}")
            # still record the chemical-space split so a skipped category shows up (as 100%
            # not-considered) in the coverage figure, consistent with the datasets panel
            universe = set().union(*[d["ikset"] for d in datasets]) if datasets else set()
            partition_rows.append({"pathogen": pathogen_code, "category": category,
                                   "n_multi": 0, "n_singleton": 0,
                                   "n_not_considered": len(universe), "n_total": len(universe)})
            continue

        elig_iks = set().union(*[d["ikset"] for d in eligible])
        elig_smiles = {ik: ik2info[ik][0] for ik in elig_iks if ik in ik2info}
        fps = morgan_matrix(elig_smiles)

        ids, AUROC, F1, N = compute_or_load_transfer(pathogen_code, category, eligible, fps, force)

        # pools at the chosen setting
        labels = greedy_modularity_pools(AUROC, POOL_THRESHOLD, rule=POOL_RULE)
        psum, pmem = build_and_write_pools(pathogen_code, category, ids, labels, eligible,
                                           meta_all, ik2info, orig_by_id)
        all_pool_summary += psum; all_pool_members += pmem
        pdf = pd.DataFrame(psum)
        n_multi = int((~pdf["is_singleton"]).sum())
        print(f"    pools ({POOL_RULE}, t={POOL_THRESHOLD}): {len(pdf)} total "
              f"({n_multi} multi-dataset + {len(pdf) - n_multi} singletons); "
              f"largest {int(pdf['n_molecules'].max())} molecules; "
              f"{int(pdf['n_active_flipped'].sum())} molecules averaged out of active")

        make_heatmap_figure(pathogen_code, category, AUROC, labels)
        part = chemspace_partition(datasets, eligible, ids, labels)
        partition_rows.append({"pathogen": pathogen_code, "category": category, **part})

    if all_first_pass:
        pd.DataFrame(all_first_pass).to_csv(os.path.join(out, "23_first_pass.csv"), index=False)
    if all_pool_summary:
        pd.DataFrame(all_pool_summary).to_csv(os.path.join(out, "23_pool_summary.csv"), index=False)
        pd.DataFrame(all_pool_members).to_csv(os.path.join(out, "23_pool_members.csv"), index=False)
        pd.DataFrame(partition_rows).to_csv(os.path.join(out, "23_chemspace_partition.csv"), index=False)
        pd.DataFrame(all_fold_rows).to_csv(os.path.join(out, "23_folds.csv"), index=False)
        print(f"  pools -> {os.path.join(out, '23_pools')} | first_pass -> 23_first_pass/ (+manifest) "
              f"| folds -> 23_folds.csv")
        make_pooling_overview(pathogen_code)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Step 23: transfer signal + pools (per pathogen)")
    ap.add_argument("pathogens", nargs="*", default=None, help="pathogen codes (default: 3 dev organisms)")
    ap.add_argument("--force", action="store_true", help="ignore cached transfer matrices")
    ap.add_argument("--plot-only", action="store_true",
                    help="redraw the per-pathogen transfer heatmaps (from cached npz) and the pooling "
                         "overview (from the 23_ CSVs); no computation. Cross-pathogen pool figures "
                         "live in step 27.")
    args = ap.parse_args()
    if args.plot_only:
        for p in (args.pathogens or _pooled_codes()):
            replot_heatmaps(p)
            make_pooling_overview(p)
    else:
        for p in (args.pathogens if args.pathogens else DEV_PATHOGENS):
            process(p, args.force)
