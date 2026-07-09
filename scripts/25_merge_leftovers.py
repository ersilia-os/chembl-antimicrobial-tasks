"""
Step 25 (pooling rebuild) - merge leftover datasets into the modelled pools.

Steps 23-24 give the source of truth: modelled pools + large single-assay screens, each with a final
RF, a cross-validated AUROC and a Youden probability cutoff. This step folds the LEFTOVERS
(non-modelled singletons, degenerate single-class, and small <20-mol datasets) into the pool whose model
best explains their activity - without disturbing the pools' existing molecules.

Per pathogen x category that has >= 1 modelled pool:
  1. Train a final RF per modelled pool (on its full molecules). modelled_mols = union of their molecules.
  2. Leftovers = every 22_binarised dataset not a member of a modelled pool. DEDUP: keep only each
     leftover's residual molecules (those not already in modelled_mols); drop fully-redundant leftovers.
  3. ASSIGN each leftover to the pool with the highest BALANCED ACCURACY at that pool's Youden cutoff
     (mean of active-recall and inactive-recall on the leftover's molecules; single-class-safe). Merge
     only if that best balanced accuracy >= MIN_ASSIGN_BALACC; below that the pool mis-predicts the
     leftover, so it is DEFERRED (low fit) rather than force-merged. NOTE deferred leftovers are NOT
     revisited by step 26 (which only builds catch-alls for categories with 0 modelled pools), so in a
     category that already has a modelled pool a deferred leftover stays UNCOVERED.
     NOTE: the score is over the leftover's FULL molecule set; if some of those molecules already sit
     inside a candidate pool they are in-sample for that pool's model, so its score can be optimistic.
  4. BATCH-MERGE per pool: append all assigned leftovers' residual molecules (conservatively resolved via
     merge_datasets) to the pool's existing rows, re-CV the grown pool ONCE, and accept the batch only if
     the grown AUROC >= baseline AUROC - AUROC_TOL. Otherwise keep the pool unchanged and flag it; its
     assigned leftovers are deferred.

Low-data pathogens/categories with no modelled pool are skipped (deferred to a later phase). Existing pool
molecules are never re-averaged (source of truth); only new residual molecules are appended.

Pools are NOT disjoint in molecule space: a molecule tested in assays that land in different pools can
belong to several pools and even carry different bins in each (a per-pool task is self-contained). This
is inherent to the pooling design and is not corrected here; 25_label_conflicts.csv reports its extent.

Inputs (output/stage4/<pathogen>/): 23_pool_members.csv, 23_pools/<cat>/<pool_id>.csv.gz, 24_cv_summary.csv,
    23_chemspace_partition.csv, 22_binarised/<cat>/*.csv.gz, 22_binarisation_summary.csv
Outputs (output/stage4/<pathogen>/):
    25_pools/<cat>/<pool_id>.csv.gz   grown (or unchanged) pools
    25_pool_summary.csv               per pool: mol before/after, n_added, AUROC baseline/grown, accepted
    25_pool_members.csv               pool -> member datasets (original + merged leftovers)
    25_merge_log.csv                  per leftover: group, residual, assigned pool, score, accepted/deferred
    25_label_conflicts.csv            per category: unique molecules, # in >1 pool, # with conflicting labels
    25_merge_overview.png             2x2: baseline vs grown AUROC (tolerance band) | coverage
                                      before/after | leftover outcomes | cross-pool label conflicts

Usage:
    python scripts/25_merge_leftovers.py                       # all pathogens with 24_cv_summary
    python scripts/25_merge_leftovers.py saureus mtuberculosis ecoli
"""

import os
import sys
import glob
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH, FP_BITS, MIN_DATASET_SIZE
from pathogen_utils import load_pathogen
from dataset_utils import UNIT_SLUG, load_category_cutoffs
from pooling_utils import morgan_matrix, load_binarised_datasets, load_dataset_meta, merge_datasets
from model_utils import TrainRF, cv_oof_predict, balanced_agreement
from plot_utils import cat_label
import plot_colors as pc
from plot_base import BasePlot, render_figure, inset_legend, compact_yaxis

# Format: print | Style: article - publication house style (shared across stage 4).
stylia.set_format("print")
stylia.set_style("article")

CATEGORIES = ["DR", "SP"]
AUROC_TOL = 0.05           # accept a pool's batch iff grown AUROC >= baseline - AUROC_TOL (trade a little
                           # AUROC for a lot of coverage; the guard still blocks real degradations)
MIN_ASSIGN_BALACC = 0.60   # a leftover is merged only if its best pool reaches this balanced accuracy
                           # ((sens+spec)/2 at the pool's Youden cutoff); below it the pool mis-predicts
                           # the leftover, so it is deferred (low fit) rather than force-merged
N_SPLITS = 5
output_base = os.path.join(root, "..", "output")


def _fp_matrix(iks, fps):
    """(X, kept_iks) for the inchikeys that have a fingerprint."""
    kept = [ik for ik in iks if ik in fps]
    X = np.stack([fps[ik] for ik in kept]) if kept else np.empty((0, FP_BITS), np.uint8)
    return X, kept


def process(pathogen_code):
    out = os.path.join(output_base, "stage4", pathogen_code)
    if not os.path.exists(os.path.join(out, "24_cv_summary.csv")):
        print(f"Step 25 - {pathogen_code}: no 24_cv_summary.csv - skipping")
        return
    print(f"Step 25 - {pathogen_code} ({load_pathogen(pathogen_code)})")
    cutoffs = load_category_cutoffs(CONFIGPATH)
    cv = pd.read_csv(os.path.join(out, "24_cv_summary.csv"))
    members = pd.read_csv(os.path.join(out, "23_pool_members.csv"))
    part = pd.read_csv(os.path.join(out, "23_chemspace_partition.csv"))
    ntot = {r.category: int(r.n_total) for r in part.itertuples(index=False)}

    all_summary, all_members, all_log, coverage, all_conflicts = [], [], [], {}, []
    for category in CATEGORIES:
        cvc = cv[(cv.category == category)]
        modelled = cvc[cvc.modelled]
        if modelled.empty:
            print(f"  {category}: 0 modelled pools - deferred (low-data)")
            continue
        ds, ik2info = load_binarised_datasets(os.path.join(out, "22_binarised", category))
        meta = load_dataset_meta(os.path.join(out, "22_binarisation_summary.csv"), category, cutoffs, UNIT_SLUG)
        fps = morgan_matrix({ik: info[0] for ik, info in ik2info.items()})

        mc = members[members.category == category]
        pid_of = dict(zip(mc.dataset_id, mc.pool_id))
        modelled_ids = set(modelled.pool_id)
        member_of_modelled = {dsid for dsid, pid in pid_of.items() if pid in modelled_ids}

        # --- modelled pools: existing rows, final model, cutoff, baseline AUROC ---
        pool_dir = os.path.join(out, "23_pools", category)
        pool = {}
        modelled_mols = set()
        for r in modelled.itertuples(index=False):
            df = pd.read_csv(os.path.join(pool_dir, r.pool_id + ".csv.gz"))
            X, kept = _fp_matrix(df["inchikey"], fps)
            y = df.set_index("inchikey").loc[kept, "bin"].to_numpy().astype(int)
            model = TrainRF(X, y)
            pool[r.pool_id] = {"df": df, "model": model, "cutoff": float(r.youden_cutoff),
                               "baseline": float(r.auroc), "members": [d for d, p in pid_of.items()
                                                                       if p == r.pool_id]}
            modelled_mols |= set(df["inchikey"])

        # --- leftovers: not a member of a modelled pool; residual = not in modelled_mols ---
        leftovers = []
        for d in ds:
            if d["id"] in member_of_modelled:
                continue
            group = ("singleton_leftover" if d["id"] in pid_of else
                     "degenerate" if d["labels"].min() == d["labels"].max() else
                     "small" if len(d["ikset"]) < MIN_DATASET_SIZE else "other")
            resid = np.array([ik not in modelled_mols for ik in d["iks"]])
            if not resid.any():
                all_log.append({"pathogen": pathogen_code, "category": category, "leftover": d["id"],
                                "group": group, "n_full": len(d["iks"]), "n_residual": 0,
                                "assigned_pool": "", "score": np.nan, "outcome": "redundant"})
                continue
            leftovers.append({"d": d, "group": group, "resid_mask": resid})

        # --- assign each leftover to best pool by balanced accuracy at that pool's cutoff ---
        # predict every leftover molecule once per pool model
        left_iks = sorted({ik for lo in leftovers for ik in lo["d"]["iks"] if ik in fps})
        probs = {}
        if left_iks:
            Xleft = np.stack([fps[ik] for ik in left_iks])
            for pid, P in pool.items():
                probs[pid] = dict(zip(left_iks, P["model"].predict_proba(Xleft)[:, 1]))
        assigned = defaultdict(list)
        for lo in leftovers:
            d = lo["d"]
            iks = [ik for ik in d["iks"] if ik in fps]
            ytrue = np.array([b for ik, b in zip(d["iks"], d["labels"]) if ik in fps])
            best_pid, best_score = None, -1.0
            if len(iks):
                for pid, P in pool.items():
                    preds = np.array([probs[pid][ik] >= P["cutoff"] for ik in iks], dtype=int)
                    s = balanced_agreement(ytrue, preds)
                    if s > best_score:
                        best_score, best_pid = s, pid
            lo["assigned"], lo["score"] = best_pid, best_score
            # Only merge if the best pool actually explains the leftover (balanced accuracy >=
            # MIN_ASSIGN_BALACC at its Youden cutoff). Below that, forcing the merge would inject
            # label noise into the pool, so defer the leftover. (Deferred leftovers are not picked up
            # by step 26 - which handles 0-modelled categories only - so here they stay uncovered.)
            if best_pid is not None and best_score >= MIN_ASSIGN_BALACC:
                assigned[best_pid].append(lo)
            else:
                all_log.append({"pathogen": pathogen_code, "category": category, "leftover": d["id"],
                                "group": lo["group"], "n_full": len(d["iks"]),
                                "n_residual": int(lo["resid_mask"].sum()),
                                "assigned_pool": best_pid or "",
                                "score": round(best_score, 4) if best_pid is not None else np.nan,
                                "outcome": "deferred_low_fit"})

        # --- batch-merge per pool + AUROC-non-decrease guard ---
        os.makedirs(os.path.join(out, "25_pools", category), exist_ok=True)
        for f in glob.glob(os.path.join(out, "25_pools", category, "*.csv.gz")):
            os.remove(f)
        cov_before, cov_after = set(), set()
        ik_npools, ik_bins = defaultdict(int), defaultdict(set)   # cross-pool label-conflict diagnostic
        for pid, P in pool.items():
            base_df = P["df"]
            los = assigned.get(pid, [])
            resid_dsets = []
            for lo in los:
                d = lo["d"]; m = lo["resid_mask"]
                resid_dsets.append({"id": d["id"], "iks": d["iks"][m],
                                    "values": d["values"][m], "relations": d["relations"][m]})
            grown_df, accepted, grown_auroc = base_df, True, P["baseline"]
            if resid_dsets:
                add_df = merge_datasets(resid_dsets, meta, ik2info)
                cand = pd.concat([base_df, add_df], ignore_index=True).drop_duplicates("inchikey")
                X, kept = _fp_matrix(cand["inchikey"], fps)
                y = cand.set_index("inchikey").loc[kept, "bin"].to_numpy().astype(int)
                oof = cv_oof_predict(X, y, n_splits=N_SPLITS)
                grown_auroc = float(roc_auc_score(y, oof))
                accepted = grown_auroc >= P["baseline"] - AUROC_TOL
                if accepted:
                    grown_df = cand
            write_df = grown_df if accepted else base_df
            write_df.to_csv(os.path.join(out, "25_pools", category, pid + ".csv.gz"), index=False)
            cov_before |= set(base_df["inchikey"]); cov_after |= set(write_df["inchikey"])
            for ik, b in zip(write_df["inchikey"], write_df["bin"]):
                ik_npools[ik] += 1; ik_bins[ik].add(int(b))
            nb, na = len(base_df), len(write_df)
            all_summary.append({"pathogen": pathogen_code, "category": category, "pool_id": pid,
                                "n_mol_before": nb, "n_mol_after": na, "n_added": na - nb,
                                "n_leftovers_assigned": len(los),
                                "ar_before": round((base_df.bin == 1).mean(), 4),
                                "ar_after": round((write_df.bin == 1).mean(), 4),
                                "baseline_auroc": round(P["baseline"], 4),
                                "grown_auroc": round(grown_auroc, 4), "batch_accepted": bool(accepted)})
            all_members += [{"category": category, "pool_id": pid, "dataset_id": m} for m in P["members"]]
            for lo in los:
                if accepted:
                    all_members.append({"category": category, "pool_id": pid, "dataset_id": lo["d"]["id"]})
                all_log.append({"pathogen": pathogen_code, "category": category, "leftover": lo["d"]["id"],
                                "group": lo["group"], "n_full": len(lo["d"]["iks"]),
                                "n_residual": int(lo["resid_mask"].sum()), "assigned_pool": pid,
                                "score": round(lo["score"], 4),
                                "outcome": "merged" if accepted else "deferred_auroc_drop"})
        acc = sum(1 for s in all_summary if s["category"] == category and s["batch_accepted"])
        nmerged = sum(1 for x in all_log if x["category"] == category and x["outcome"] == "merged")
        coverage[category] = (len(cov_before), len(cov_after), ntot.get(category, 0))
        n_multi = sum(1 for c in ik_npools.values() if c > 1)
        n_conf = sum(1 for ik, c in ik_npools.items() if c > 1 and len(ik_bins[ik]) > 1)
        all_conflicts.append({"pathogen": pathogen_code, "category": category,
                              "n_unique_molecules": len(ik_npools), "n_in_multiple_pools": n_multi,
                              "n_conflicting_label": n_conf})
        print(f"  {category}: {len(pool)} pools ({acc} accepted merges), {len(leftovers)} leftovers "
              f"-> {nmerged} merged; coverage {cov_pct(coverage[category])}; "
              f"{n_conf} molecules with conflicting labels across pools")

    if all_summary:
        pd.DataFrame(all_summary).to_csv(os.path.join(out, "25_pool_summary.csv"), index=False)
        pd.DataFrame(all_members).drop_duplicates().to_csv(os.path.join(out, "25_pool_members.csv"), index=False)
        pd.DataFrame(all_log).to_csv(os.path.join(out, "25_merge_log.csv"), index=False)
        pd.DataFrame(all_conflicts).to_csv(os.path.join(out, "25_label_conflicts.csv"), index=False)
        make_overview(pathogen_code, pd.DataFrame(all_summary), coverage,
                      pd.DataFrame(all_log), pd.DataFrame(all_conflicts))
        print(f"  -> 25_pools/ | 25_pool_summary.csv | 25_merge_log.csv | 25_label_conflicts.csv")


def cov_pct(cov):
    before, after, total = cov
    return (f"{100*before/total:.1f}% -> {100*after/total:.1f}%" if total else "n/a")


# --- plots (shared stage-4 plot_base framework) -----------------------------
# leftover outcomes (merge_log `outcome` -> label, colour)
_OUTCOME = [("merged", "merged", "turquoise"),
            ("deferred_low_fit", "deferred (low fit)", "crimson"),
            ("deferred_auroc_drop", "deferred (AUROC drop)", "amber"),
            ("redundant", "redundant (no new mols)", "silver")]


class MergeAurocScatter(BasePlot):
    """Baseline vs grown CV AUROC per pool: accepted (turquoise) vs rejected (crimson) batches, marker
    size ~ leftovers assigned. Shaded band = the AUROC_TOL tolerance (accepted despite a small drop);
    below it batches are rejected."""

    def __init__(self, ax, summ, abc=None):
        super().__init__(ax, abc)
        xs = np.linspace(0.5, 1, 100)
        ax.fill_between(xs, xs - AUROC_TOL, xs, color=pc.color("neutral"), alpha=0.15, linewidth=0)
        ax.plot([0.5, 1], [0.5, 1], color=pc.color("neutral"), linestyle="--")
        sz = summ["n_leftovers_assigned"].to_numpy(dtype=float)
        smax = sz.max() if len(sz) and sz.max() > 0 else 1.0
        entries = []
        for acc, cname, lab in [(True, "turquoise", "accepted"), (False, "crimson", "rejected")]:
            s = summ[summ.batch_accepted == acc]
            if len(s):
                ax.scatter(s.baseline_auroc, s.grown_auroc,
                           s=20 + 180 * s.n_leftovers_assigned.to_numpy(dtype=float) / smax,
                           color=pc.named(cname), edgecolors="none", alpha=0.7)
            entries.append((lab, pc.named(cname)))
        ax.set_xlim(0.5, 1); ax.set_ylim(0.5, 1)
        self._label(xlabel="baseline AUROC", ylabel="grown AUROC",
                    title="Pool AUROC after merge (size ~ leftovers)")
        inset_legend(ax, loc="upper left", entries=entries)
        self.data = summ[["pool_id", "baseline_auroc", "grown_auroc", "n_leftovers_assigned",
                          "batch_accepted"]].copy()


class CoverageBars(BasePlot):
    """Fraction of category chemical space covered by the modelled pools, before vs after merging."""

    def __init__(self, ax, coverage, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if c in coverage and coverage[c][2]]
        if not cats:
            self._empty(ylabel="fraction of chem space in modelled pools",
                        title="Coverage before vs after"); return
        x = np.arange(len(cats)); w = 0.35
        before = [coverage[c][0] / coverage[c][2] for c in cats]   # UNION coverage of chem space
        after = [coverage[c][1] / coverage[c][2] for c in cats]
        ax.bar(x - w / 2, before, w, color=pc.color("neutral"))
        ax.bar(x + w / 2, after, w, color=pc.color("datasets"))
        ax.set_ylim(0, 1)
        ax.set_xticks(x); ax.set_xticklabels([cat_label(c) for c in cats])
        self._label(xlabel="", ylabel="fraction of chem space in modelled pools",
                    title="Coverage before vs after")
        inset_legend(ax, loc="upper left",
                     entries=[("before", pc.color("neutral")), ("after", pc.color("datasets"))])
        self.data = pd.DataFrame({"category": cats, "before": before, "after": after})


class LeftoverFate(BasePlot):
    """Leftover outcomes per category (from 25_merge_log): merged / deferred (low fit) /
    deferred (AUROC drop) / redundant (no residual molecules). The dataset-side complement to
    the coverage bars."""

    def __init__(self, ax, log_df, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if len(log_df) and (log_df.category == c).any()]
        if not cats:
            self._empty(ylabel="leftovers", title="Leftover outcomes"); return
        x = np.arange(len(cats)); bottom = np.zeros(len(cats))
        recs = []
        for oc, lab, cname in _OUTCOME:
            vals = np.array([int(((log_df.category == c) & (log_df.outcome == oc)).sum()) for c in cats],
                            dtype=float)
            ax.bar(x, vals, bottom=bottom, color=pc.named(cname))
            bottom += vals
            recs.append({"outcome": lab, **{c: vals[i] for i, c in enumerate(cats)}})
        ax.set_xticks(x); ax.set_xticklabels([cat_label(c) for c in cats])
        compact_yaxis(ax)
        self._label(xlabel="", ylabel="leftovers", title="Leftover outcomes")
        inset_legend(ax, loc="upper right", entries=[(lab, pc.named(cn)) for _, lab, cn in _OUTCOME])
        self.data = pd.DataFrame(recs)


class LabelConflictBars(BasePlot):
    """Unique molecules per category as % of the category chemical space, split by pool membership:
    in one pool / in >1 pool with a consistent label / in >1 pool with conflicting labels. Shown as
    fractions of the total so the conflicting share is seen in context (its % is annotated)."""

    def __init__(self, ax, conf_df, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if conf_df is not None and len(conf_df)
                and (conf_df.category == c).any()]
        if not cats:
            self._empty(ylabel="% of unique molecules", title="Cross-pool label conflicts"); return
        g = conf_df.set_index("category")
        x = np.arange(len(cats))
        tot = np.array([max(int(g.loc[c, "n_unique_molecules"]), 1) for c in cats], dtype=float)
        multi = np.array([int(g.loc[c, "n_in_multiple_pools"]) for c in cats], dtype=float)
        conflicting = np.array([int(g.loc[c, "n_conflicting_label"]) for c in cats], dtype=float)
        segs = [("one pool", (tot - multi) / tot * 100, pc.named("silver")),
                (">1 pool, consistent", (multi - conflicting) / tot * 100, pc.named("periwinkle")),
                (">1 pool, conflicting", conflicting / tot * 100, pc.color("active"))]
        bottom = np.zeros(len(cats))
        for _, vals, col in segs:
            ax.bar(x, vals, bottom=bottom, color=col)
            bottom += vals
        for xi, cf, t in zip(x, conflicting, tot):       # annotate the conflicting share
            ax.text(xi, 101, f"{100 * cf / t:.1f}%", ha="center", va="bottom",
                    fontsize=stylia.FONTSIZE_SMALL, color=pc.color("active"))
        ax.set_ylim(0, 108)
        ax.set_xticks(x); ax.set_xticklabels([cat_label(c) for c in cats])
        self._label(xlabel="", ylabel="% of unique molecules", title="Cross-pool label conflicts")
        inset_legend(ax, loc="center right", entries=[(l, c) for l, _, c in segs])
        self.data = pd.DataFrame({"category": cats, "pct_multi": multi / tot * 100,
                                  "pct_conflicting": conflicting / tot * 100})


def make_overview(pathogen_code, summ, coverage, log_df, conf_df):
    """2x2 via plot_base: A baseline-vs-grown AUROC scatter (tolerance band, size ~ leftovers),
    B chem-space coverage before/after, C leftover outcomes, D cross-pool label conflicts."""
    out = os.path.join(output_base, "stage4", pathogen_code)
    render_figure(out, "25_merge_overview", 2, 2, [
        lambda ax: MergeAurocScatter(ax, summ, abc="A"),
        lambda ax: CoverageBars(ax, coverage, abc="B"),
        lambda ax: LeftoverFate(ax, log_df, abc="C"),
        lambda ax: LabelConflictBars(ax, conf_df, abc="D"),
    ], width=0.7, height=0.7)      # smaller, square-ish panels


def _union_iks(pool_dir, pool_ids):
    """Union of inchikeys across the given pool files in a directory (cheap I/O)."""
    s = set()
    for pid in pool_ids:
        f = os.path.join(pool_dir, pid + ".csv.gz")
        if os.path.exists(f):
            s |= set(pd.read_csv(f, usecols=["inchikey"])["inchikey"])
    return s


def _label_conflicts(pathogen_code, out):
    """Cross-pool label-conflict diagnostic read from the on-disk 25_pools (no models needed):
    per category, unique molecules, how many are in >1 pool, and how many carry conflicting bins."""
    rows = []
    for cat in CATEGORIES:
        ik_npools, ik_bins = defaultdict(int), defaultdict(set)
        for f in glob.glob(os.path.join(out, "25_pools", cat, "*.csv.gz")):
            df = pd.read_csv(f, usecols=["inchikey", "bin"])
            for ik, b in zip(df["inchikey"], df["bin"]):
                ik_npools[ik] += 1; ik_bins[ik].add(int(b))
        if not ik_npools:
            continue
        rows.append({"pathogen": pathogen_code, "category": cat, "n_unique_molecules": len(ik_npools),
                     "n_in_multiple_pools": sum(1 for c in ik_npools.values() if c > 1),
                     "n_conflicting_label": sum(1 for ik, c in ik_npools.items()
                                                if c > 1 and len(ik_bins[ik]) > 1)})
    return rows


def plot_only(pathogen_code):
    """Redraw 25_merge_overview.png from cached outputs (25_pool_summary.csv + the pool
    molecule files). Coverage is recomputed by unioning inchikeys from the 23_pools
    (before) and 25_pools (after) files - no models are retrained."""
    out = os.path.join(output_base, "stage4", pathogen_code)
    sp = os.path.join(out, "25_pool_summary.csv")
    if not os.path.exists(sp):
        print(f"Step 25 (plot-only) - {pathogen_code}: no 25_pool_summary.csv - skipping")
        return
    summ = pd.read_csv(sp)
    part = pd.read_csv(os.path.join(out, "23_chemspace_partition.csv"))
    ntot = {r.category: int(r.n_total) for r in part.itertuples(index=False)}
    coverage = {}
    for cat in CATEGORIES:
        pids = summ[summ.category == cat].pool_id.astype(str).tolist()
        if not pids:
            continue
        before = _union_iks(os.path.join(out, "23_pools", cat), pids)
        after = _union_iks(os.path.join(out, "25_pools", cat), pids)
        coverage[cat] = (len(before), len(after), ntot.get(cat, 0))
    conf_df = pd.DataFrame(_label_conflicts(pathogen_code, out))
    if len(conf_df):
        conf_df.to_csv(os.path.join(out, "25_label_conflicts.csv"), index=False)
    log_path = os.path.join(out, "25_merge_log.csv")
    log_df = pd.read_csv(log_path) if os.path.exists(log_path) else pd.DataFrame(columns=["category", "outcome"])
    make_overview(pathogen_code, summ, coverage, log_df, conf_df)


def _codes_with(artifact):
    return sorted(os.path.basename(os.path.dirname(f))
                  for f in glob.glob(os.path.join(output_base, "stage4", "*", artifact)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Step 25: merge leftovers into modelled pools")
    ap.add_argument("pathogens", nargs="*", help="default: all pathogens with 24_cv_summary.csv")
    ap.add_argument("--plot-only", action="store_true",
                    help="redraw 25_merge_overview.png from cached pool outputs; skip the merge/CV")
    args = ap.parse_args()
    if args.plot_only:
        for p in (args.pathogens or _codes_with("25_pool_summary.csv")):
            plot_only(p)
    else:
        for p in (args.pathogens or _codes_with("24_cv_summary.csv")):
            process(p)
