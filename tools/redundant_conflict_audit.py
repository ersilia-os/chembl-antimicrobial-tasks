"""
redundant_conflict_audit.py  (tools/ — read-only diagnostic; NOT part of the numbered pipeline)

Audits whether step-25 `redundant` leftovers (leftovers whose every molecule is already in a
modelled pool, so they are dropped with `n_residual == 0`) silently lose signal, and whether the
step-25 pool-growth policy itself changes labels. Three sub-analyses, run individually or together:

  conflicts   Raw label comparison: for each redundant leftover, compare its per-assay bin against
              the bin those molecules carry in the pools. Three definitions per molecule --
              consensus (majority of the pools' bins, tie -> active, matching merge_datasets),
              strict (disagrees with EVERY containing pool), any (disagrees with >=1). Plus a
              driver breakdown by activity_type.
  merge-flip  Pipeline-faithful test: merge each redundant leftover into every pool containing its
              molecules and RE-BINARISE via merge_datasets. Reports flip_vs_delivered (would the
              shipped label change), flip_marginal (the leftover's pure effect), and
              delivered_ne_remerge (a step-25 artifact — see below).
  reaverage   Pool-WIDE, de-duplicated version of delivered_ne_remerge: over EVERY molecule in EVERY
              pool, re-merge all listed members and report what fraction of shipped labels would
              flip if step 25 re-averaged (vs its current "base wins / never re-average" policy).
              Molecule set held fixed; cross-pool molecules a naive full re-merge would add are
              reported as `n_would_add`, not flips.

Everything is READ-ONLY — nothing in output/ or config/ is modified. See the write-up
docs/2026-07-15_redundant-conflict-diagnostics.md for the findings and interpretation.

Step-25 artifact (context for merge-flip / reaverage). Step 25 grows a pool by appending an
assigned leftover's RESIDUAL molecules only; for a molecule already in the pool it keeps the base
value and DISCARDS the leftover's reading (`drop_duplicates("inchikey")`,
scripts/25_merge_leftovers.py). Yet the leftover is written as a full member in 25_pool_members.csv.
So a clean re-merge of the listed members can yield a different (re-averaged) label than the one
shipped. `reaverage` measures how often, pool-wide.

Inputs (per pathogen <p>, read-only): output/stage4/<p>/{25_merge_log.csv, 25_pool_members.csv,
  25_pools/<cat>/*.csv.gz, 22_binarised/<cat>/*.csv.gz, 22_binarisation_summary.csv}.
Outputs (docs/, dated 2026-07-15):
  conflicts  -> *_redundant-conflicts_{summary,detail,by_type}.csv + *_redundant-conflicts{,_by_type}.png
  merge-flip -> *_redundant-mergeflip_{summary,detail}.csv + *_redundant-mergeflip.png
  reaverage  -> *_redundant-reaverage_{summary,detail}.csv + *_redundant-reaverage.png

Usage:
  python tools/redundant_conflict_audit.py                       # all sub-analyses, all pathogens
  python tools/redundant_conflict_audit.py conflicts             # one sub-analysis
  python tools/redundant_conflict_audit.py all ecoli saureus     # subset of pathogens
"""

import os
import sys
import glob
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import STAGE4PATH, CONFIGPATH  # noqa: E402
from dataset_utils import UNIT_SLUG, load_category_cutoffs  # noqa: E402
from pooling_utils import load_binarised_datasets, load_dataset_meta, merge_datasets  # noqa: E402
from plot_base import BasePlot, render_figure, inset_legend, compact_yaxis  # noqa: E402
import plot_colors as pc  # noqa: E402

CATEGORIES = ["DR", "SP"]
DATE = "2026-07-15"                     # stable output-name prefix (matches the write-up date)
N_TYPES_FIG = 12                        # activity types shown in the conflict driver figure

docs_dir = os.path.join(root, "..", "docs")
os.makedirs(docs_dir, exist_ok=True)


def _out(name):
    return os.path.join(docs_dir, f"{DATE}_{name}")


# --------------------------------------------------------------------------
# Shared helpers
# --------------------------------------------------------------------------
def _pathogens():
    df = pd.read_csv(os.path.join(CONFIGPATH, "pathogens.csv"))
    return list(zip(df["pathogen"], df["code"]))


def _codes(args_codes, needs_file):
    return args_codes or [c for _, c in _pathogens()
                          if os.path.exists(os.path.join(STAGE4PATH, c, needs_file))]


def _activity_type(dataset_id):
    """activity_type from a dataset_id {assay_id}_{activity_type}_{unit_slug}; assay_id and unit
    slug are single tokens, so activity_type is everything in between (handles underscores)."""
    parts = str(dataset_id).split("_")
    return "_".join(parts[1:-1]) if len(parts) >= 3 else "?"


def _consensus(bins):
    """Majority bin across the pools a molecule is in, tie -> active (1) (merge_datasets' rule)."""
    return 1 if float(np.mean(bins)) >= 0.5 else 0


def _reduce(d, keep):
    """Restrict a dataset dict to molecules in `keep`; None if it has none. Preserves the keys
    merge_datasets consumes (id/iks/values/relations)."""
    mask = np.fromiter((ik in keep for ik in d["iks"]), dtype=bool, count=len(d["iks"]))
    if not mask.any():
        return None
    return {"id": d["id"], "iks": d["iks"][mask], "values": d["values"][mask],
            "relations": d["relations"][mask]}


def _plot_frame(summary, rate_col):
    """Per-pathogen wide frame for the DR/SP rate + direction figures shared by all three analyses."""
    paths = [c for _, c in _pathogens() if c in set(summary["pathogen"])]
    rows = []
    for code in paths:
        s = summary[summary["pathogen"] == code]
        by_cat = {r["category"]: r for _, r in s.iterrows()}
        rows.append({
            "pathogen": code,
            "rate_DR": (by_cat["DR"][rate_col] if "DR" in by_cat else np.nan),
            "rate_SP": (by_cat["SP"][rate_col] if "SP" in by_cat else np.nan),
            "flip_to_active": int(s.get("n_flip_to_active", pd.Series(dtype=int)).sum())
                if "n_flip_to_active" in s else 0,
            "flip_to_inactive": int(s.get("n_flip_to_inactive", pd.Series(dtype=int)).sum())
                if "n_flip_to_inactive" in s else 0,
        })
    return pd.DataFrame(rows)


class _RateBars(BasePlot):
    """Per-pathogen DR/SP grouped rate bars (percent)."""

    def __init__(self, ax, pdf, ylabel, title, abc=None):
        super().__init__(ax, abc)
        x = np.arange(len(pdf)); w = 0.4
        ax.bar(x - w / 2, pdf["rate_DR"].fillna(0) * 100, w, color=pc.color("DR"))
        ax.bar(x + w / 2, pdf["rate_SP"].fillna(0) * 100, w, color=pc.color("SP"))
        ax.set_xticks(x); ax.set_xticklabels(pdf["pathogen"], rotation=45, ha="right")
        self._label(xlabel="", ylabel=ylabel, title=title)
        inset_legend(ax, loc="upper right",
                     entries=[("DR", pc.color("DR")), ("SP", pc.color("SP"))])
        self.data = pdf[["pathogen", "rate_DR", "rate_SP"]].copy()


class _DirBars(BasePlot):
    """Per-pathogen stacked direction counts: -> active over -> inactive."""

    def __init__(self, ax, pdf, ylabel, active_lbl, inactive_lbl, abc=None):
        super().__init__(ax, abc)
        x = np.arange(len(pdf))
        to_i = pdf["flip_to_inactive"].to_numpy(dtype=float)
        to_a = pdf["flip_to_active"].to_numpy(dtype=float)
        ax.bar(x, to_i, color=pc.color("inactive"))
        ax.bar(x, to_a, bottom=to_i, color=pc.color("active"))
        compact_yaxis(ax)
        ax.set_xticks(x); ax.set_xticklabels(pdf["pathogen"], rotation=45, ha="right")
        self._label(xlabel="", ylabel=ylabel, title="Conflict / flip direction")
        inset_legend(ax, loc="upper right",
                     entries=[(active_lbl, pc.color("active")), (inactive_lbl, pc.color("inactive"))])
        self.data = pdf[["pathogen", "flip_to_active", "flip_to_inactive"]].copy()


# --------------------------------------------------------------------------
# Analysis 1: conflicts (raw label comparison + drivers)   [was scripts/29]
# --------------------------------------------------------------------------
def _load_pool_labels(s4, category):
    pool_bins, pool_ids = defaultdict(list), defaultdict(set)
    for f in sorted(glob.glob(os.path.join(s4, "25_pools", category, "*.csv.gz"))):
        pid = os.path.basename(f)[: -len(".csv.gz")]
        d = pd.read_csv(f, usecols=["inchikey", "bin"])
        for ik, b in zip(d["inchikey"], d["bin"]):
            pool_bins[ik].append(int(b)); pool_ids[ik].add(pid)
    return pool_bins, pool_ids


def _trace_leftover(path, pool_bins, pool_ids):
    if not os.path.exists(path):
        return None
    lf = pd.read_csv(path, usecols=["inchikey", "bin"])
    n_matched = n_unmatched = c_cons = c_strict = c_any = 0
    n_active_lost = n_inactive_flip = 0
    pools_touched = set()
    for ik, lb in zip(lf["inchikey"], lf["bin"]):
        lb = int(lb)
        pb = pool_bins.get(ik)
        if not pb:
            n_unmatched += 1; continue
        n_matched += 1
        pools_touched |= pool_ids.get(ik, set())
        cons = _consensus(pb)
        if lb != cons:
            c_cons += 1
            if lb == 1 and cons == 0:
                n_active_lost += 1
            else:
                n_inactive_flip += 1
        if lb not in set(pb):
            c_strict += 1
        if any(lb != x for x in pb):
            c_any += 1
    return {"n_matched": n_matched, "n_unmatched": n_unmatched,
            "n_conflict_consensus": c_cons, "n_conflict_strict": c_strict, "n_conflict_any": c_any,
            "conflict_rate_consensus": round(c_cons / n_matched, 4) if n_matched else np.nan,
            "n_active_lost": n_active_lost, "n_inactive_flipped": n_inactive_flip,
            "n_pools_touched": len(pools_touched)}


def _conflicts_pathogen(code):
    s4 = os.path.join(STAGE4PATH, code)
    mlog_path = os.path.join(s4, "25_merge_log.csv")
    if not os.path.exists(mlog_path):
        return [], []
    mlog = pd.read_csv(mlog_path)
    detail_rows, summary_rows = [], []
    for category in CATEGORIES:
        red = mlog[(mlog["category"] == category) & (mlog["outcome"] == "redundant")]
        if red.empty:
            continue
        pool_bins, pool_ids = _load_pool_labels(s4, category)
        cat_detail, n_missing = [], 0
        for _, r in red.iterrows():
            leftover = str(r["leftover"])
            path = os.path.join(s4, "22_binarised", category, f"{leftover}.csv.gz")
            traced = _trace_leftover(path, pool_bins, pool_ids)
            base = {"pathogen": code, "category": category, "leftover": leftover,
                    "activity_type": _activity_type(leftover), "group": r.get("group"),
                    "n_full": r.get("n_full")}
            if traced is None:
                n_missing += 1; base["file_missing"] = True
            else:
                base.update(traced); base["file_missing"] = False
            cat_detail.append(base)
        detail_rows += cat_detail
        keep = [d for d in cat_detail if not d.get("file_missing")]
        nm = sum(d["n_matched"] for d in keep)
        summary_rows.append({
            "pathogen": code, "category": category, "n_redundant_leftovers": len(cat_detail),
            "n_files_missing": n_missing, "n_matched": nm,
            "n_conflict_consensus": sum(d["n_conflict_consensus"] for d in keep),
            "n_conflict_strict": sum(d["n_conflict_strict"] for d in keep),
            "n_conflict_any": sum(d["n_conflict_any"] for d in keep),
            "pct_conflict_consensus": round(sum(d["n_conflict_consensus"] for d in keep) / nm, 4) if nm else np.nan,
            "pct_conflict_strict": round(sum(d["n_conflict_strict"] for d in keep) / nm, 4) if nm else np.nan,
            "pct_conflict_any": round(sum(d["n_conflict_any"] for d in keep) / nm, 4) if nm else np.nan,
            "n_leftover_active_pool_inactive": sum(d["n_active_lost"] for d in keep),
            "n_leftover_inactive_pool_active": sum(d["n_inactive_flipped"] for d in keep),
            "n_leftovers_with_conflict_consensus": sum(1 for d in keep if d["n_conflict_consensus"] > 0),
            "n_flip_to_active": sum(d["n_active_lost"] for d in keep),      # for the shared direction fig
            "n_flip_to_inactive": sum(d["n_inactive_flipped"] for d in keep),
        })
        print(f"  [conflicts {code}/{category}] {len(cat_detail)} redundant | matched {nm} | "
              f"consensus {summary_rows[-1]['n_conflict_consensus']} "
              f"({summary_rows[-1]['pct_conflict_consensus']})")
    return detail_rows, summary_rows


def _by_type_table(detail):
    d = detail[detail.get("file_missing") != True].copy()  # noqa: E712
    if d.empty:
        return pd.DataFrame()
    g = d.groupby(["activity_type", "category"], dropna=False).agg(
        n_leftovers=("leftover", "count"), n_matched=("n_matched", "sum"),
        n_conflict_consensus=("n_conflict_consensus", "sum"),
        n_conflict_strict=("n_conflict_strict", "sum"), n_conflict_any=("n_conflict_any", "sum"),
        n_active_lost=("n_active_lost", "sum")).reset_index()
    for k in ("consensus", "strict", "any"):
        g[f"pct_conflict_{k}"] = (g[f"n_conflict_{k}"] / g["n_matched"]).round(4)
    return g.sort_values("n_matched", ascending=False).reset_index(drop=True)


class _ConflictRateByType(BasePlot):
    def __init__(self, ax, tab, abc=None):
        super().__init__(ax, abc)
        vol = tab.groupby("activity_type")["n_matched"].sum().sort_values(ascending=False)
        types = list(vol.index[:N_TYPES_FIG])
        y = np.arange(len(types))[::-1]; h = 0.4
        piv = tab.pivot_table(index="activity_type", columns="category",
                              values="pct_conflict_consensus", aggfunc="first")
        gv = lambda t, c: (piv.loc[t, c] * 100 if (c in piv.columns and t in piv.index
                           and pd.notna(piv.loc[t, c])) else 0)
        ax.barh(y + h / 2, [gv(t, "DR") for t in types], h, color=pc.color("DR"))
        ax.barh(y - h / 2, [gv(t, "SP") for t in types], h, color=pc.color("SP"))
        ax.set_yticks(y); ax.set_yticklabels(types)
        self._label(xlabel="consensus conflict rate / %", ylabel="", title="Conflict rate by activity type")
        inset_legend(ax, loc="lower right", entries=[("DR", pc.color("DR")), ("SP", pc.color("SP"))])
        self.data = tab.copy()


class _VolumeByType(BasePlot):
    def __init__(self, ax, tab, abc=None):
        super().__init__(ax, abc)
        vol = tab.groupby("activity_type")["n_matched"].sum().sort_values(ascending=False)
        types = list(vol.index[:N_TYPES_FIG]); y = np.arange(len(types))[::-1]
        ax.barh(y, [vol[t] for t in types], color=pc.named("cobalt"))
        ax.set_yticks(y); ax.set_yticklabels(types)
        self._label(xlabel="matched molecules", ylabel="", title="Volume by activity type")
        self.data = vol.reset_index()


def run_conflicts(codes):
    all_d, all_s = [], []
    for code in codes:
        print(f"[conflicts] {code}")
        d, s = _conflicts_pathogen(code); all_d += d; all_s += s
    detail, summary = pd.DataFrame(all_d), pd.DataFrame(all_s)
    bytype = _by_type_table(detail)
    detail.to_csv(_out("redundant-conflicts_detail.csv"), index=False)
    summary.to_csv(_out("redundant-conflicts_summary.csv"), index=False)
    bytype.to_csv(_out("redundant-conflicts_by_type.csv"), index=False)
    if not summary.empty:
        nm = int(summary["n_matched"].sum())
        for k in ("consensus", "strict", "any"):
            tot = int(summary[f"n_conflict_{k}"].sum())
            print(f"  overall {k}: {tot} ({tot / nm:.1%} of {nm} matched)" if nm else "")
        pdf = _plot_frame(summary, "pct_conflict_consensus")
        render_figure(docs_dir, f"{DATE}_redundant-conflicts", 1, 2, [
            lambda ax: _RateBars(ax, pdf, "consensus conflict rate / %",
                                 "Redundant-vs-pool conflict rate", abc="A"),
            lambda ax: _DirBars(ax, pdf, "conflicting molecules", "active dropped",
                                "inactive flipped", abc="B"),
        ])
        if not bytype.empty:
            render_figure(docs_dir, f"{DATE}_redundant-conflicts_by_type", 1, 2, [
                lambda ax: _ConflictRateByType(ax, bytype, abc="A"),
                lambda ax: _VolumeByType(ax, bytype, abc="B"),
            ])
    print(f"  -> docs/{DATE}_redundant-conflicts_*.csv/.png")


# --------------------------------------------------------------------------
# Analysis 2: merge-flip (merge + re-binarise, every containing pool)   [was scripts/30]
# --------------------------------------------------------------------------
def _mergeflip_pathogen(code):
    s4 = os.path.join(STAGE4PATH, code)
    if not (os.path.exists(os.path.join(s4, "25_merge_log.csv"))
            and os.path.exists(os.path.join(s4, "25_pool_members.csv"))):
        return [], []
    mlog = pd.read_csv(os.path.join(s4, "25_merge_log.csv"))
    pmem = pd.read_csv(os.path.join(s4, "25_pool_members.csv"))
    cutoffs = load_category_cutoffs(CONFIGPATH, code)
    detail_rows, summary_rows = [], []
    for category in CATEGORIES:
        red = mlog[(mlog["category"] == category) & (mlog["outcome"] == "redundant")]
        if red.empty:
            continue
        ds, ik2info = load_binarised_datasets(os.path.join(s4, "22_binarised", category))
        if not ds:
            continue
        by_id = {d["id"]: d for d in ds}
        for d in ds:
            d["_ikset"] = frozenset(d["iks"].tolist())
        meta = load_dataset_meta(os.path.join(s4, "22_binarisation_summary.csv"),
                                 category, cutoffs, UNIT_SLUG)
        pool_bin = {}
        for f in sorted(glob.glob(os.path.join(s4, "25_pools", category, "*.csv.gz"))):
            pid = os.path.basename(f)[: -len(".csv.gz")]
            pdf = pd.read_csv(f, usecols=["inchikey", "bin"])
            pool_bin[pid] = dict(zip(pdf["inchikey"], pdf["bin"].astype(int)))
        pool_members = {pid: [by_id[i] for i in grp["dataset_id"] if i in by_id]
                        for pid, grp in pmem[pmem["category"] == category].groupby("pool_id")}
        pools = [p for p in pool_bin if p in pool_members]

        cat_pairs = c_deliv = c_marg = c_dnr = to_a = to_i = 0
        mol_flip_any, left_flip, n_missing, n_eval = set(), set(), 0, 0
        for _, r in red.iterrows():
            leftover = str(r["leftover"]); dL = by_id.get(leftover)
            if dL is None:
                n_missing += 1; continue
            n_eval += 1; Lmol = dL["_ikset"]
            for pid in pools:
                shared = Lmol & set(pool_bin[pid])
                if not shared:
                    continue
                redL = _reduce(dL, shared)
                redM = [x for x in (_reduce(d, shared) for d in pool_members[pid]
                                    if d["_ikset"] & shared) if x]
                remerge = dict(zip(*(lambda m: (m["inchikey"], m["bin"].astype(int)))(
                    merge_datasets(redM, meta, ik2info))))
                test = dict(zip(*(lambda m: (m["inchikey"], m["bin"].astype(int)))(
                    merge_datasets(redM + [redL], meta, ik2info))))
                nf_d = nf_m = nf_dnr = nf_a = nf_i = ns = 0
                for ik in shared:
                    rb, tb, db = remerge.get(ik), test.get(ik), pool_bin[pid][ik]
                    if rb is None or tb is None:
                        continue
                    ns += 1
                    if db != rb:
                        nf_dnr += 1
                    if rb != tb:
                        nf_m += 1
                    if db != tb:
                        nf_d += 1; mol_flip_any.add((leftover, ik))
                        if db == 0 and tb == 1:
                            nf_a += 1
                        else:
                            nf_i += 1
                if ns == 0:
                    continue
                cat_pairs += ns; c_deliv += nf_d; c_marg += nf_m; c_dnr += nf_dnr
                to_a += nf_a; to_i += nf_i
                if nf_d:
                    left_flip.add(leftover)
                detail_rows.append({"pathogen": code, "category": category, "leftover": leftover,
                                    "activity_type": _activity_type(leftover), "pool": pid,
                                    "n_shared": ns, "n_flip_vs_delivered": nf_d, "n_flip_marginal": nf_m,
                                    "n_delivered_ne_remerge": nf_dnr, "n_flip_to_active": nf_a,
                                    "n_flip_to_inactive": nf_i,
                                    "flip_rate_vs_delivered": round(nf_d / ns, 4) if ns else np.nan})
        summary_rows.append({"pathogen": code, "category": category, "n_redundant_evaluated": n_eval,
                             "n_files_missing": n_missing, "n_molecule_pool_pairs": cat_pairs,
                             "n_flip_vs_delivered": c_deliv,
                             "pct_flip_vs_delivered": round(c_deliv / cat_pairs, 4) if cat_pairs else np.nan,
                             "n_flip_marginal": c_marg,
                             "pct_flip_marginal": round(c_marg / cat_pairs, 4) if cat_pairs else np.nan,
                             "n_delivered_ne_remerge": c_dnr,
                             "pct_delivered_ne_remerge": round(c_dnr / cat_pairs, 4) if cat_pairs else np.nan,
                             "n_flip_to_active": to_a, "n_flip_to_inactive": to_i,
                             "n_molecules_flip_any_pool": len(mol_flip_any),
                             "n_leftovers_with_flip": len(left_flip),
                             "pct_leftovers_with_flip": round(len(left_flip) / n_eval, 4) if n_eval else np.nan})
        print(f"  [merge-flip {code}/{category}] pairs {cat_pairs} | flip_vs_delivered {c_deliv} "
              f"({summary_rows[-1]['pct_flip_vs_delivered']}) | marginal {c_marg} | deliv!=remerge {c_dnr}")
    return detail_rows, summary_rows


def run_mergeflip(codes):
    all_d, all_s = [], []
    for code in codes:
        print(f"[merge-flip] {code}")
        d, s = _mergeflip_pathogen(code); all_d += d; all_s += s
    detail, summary = pd.DataFrame(all_d), pd.DataFrame(all_s)
    detail.to_csv(_out("redundant-mergeflip_detail.csv"), index=False)
    summary.to_csv(_out("redundant-mergeflip_summary.csv"), index=False)
    if not summary.empty:
        pairs = int(summary["n_molecule_pool_pairs"].sum())
        if pairs:
            for lbl, col in [("flip vs delivered", "n_flip_vs_delivered"),
                             ("flip marginal", "n_flip_marginal"),
                             ("delivered!=remerge", "n_delivered_ne_remerge")]:
                v = int(summary[col].sum())
                print(f"  overall {lbl}: {v} ({v / pairs:.1%})")
            print(f"  direction -> active {int(summary['n_flip_to_active'].sum())} | "
                  f"-> inactive {int(summary['n_flip_to_inactive'].sum())}")
        pdf = _plot_frame(summary, "pct_flip_vs_delivered")
        render_figure(docs_dir, f"{DATE}_redundant-mergeflip", 1, 2, [
            lambda ax: _RateBars(ax, pdf, "merge flip rate / %",
                                 "Would-flip rate (merge + re-binarise)", abc="A"),
            lambda ax: _DirBars(ax, pdf, "flipped molecule-pool pairs", "-> active",
                                "-> inactive", abc="B"),
        ])
    print(f"  -> docs/{DATE}_redundant-mergeflip_*.csv/.png")


# --------------------------------------------------------------------------
# Analysis 3: reaverage (pool-wide de-duplicated re-average flip)   [was scripts/31]
# --------------------------------------------------------------------------
def _reaverage_pathogen(code):
    s4 = os.path.join(STAGE4PATH, code)
    if not os.path.exists(os.path.join(s4, "25_pool_members.csv")):
        return [], []
    pmem = pd.read_csv(os.path.join(s4, "25_pool_members.csv"))
    cutoffs = load_category_cutoffs(CONFIGPATH, code)
    detail_rows, summary_rows = [], []
    for category in CATEGORIES:
        pool_files = sorted(glob.glob(os.path.join(s4, "25_pools", category, "*.csv.gz")))
        if not pool_files:
            continue
        ds, ik2info = load_binarised_datasets(os.path.join(s4, "22_binarised", category))
        by_id = {d["id"]: d for d in ds}
        meta = load_dataset_meta(os.path.join(s4, "22_binarisation_summary.csv"),
                                 category, cutoffs, UNIT_SLUG)
        members_by_pool = {pid: [by_id[i] for i in grp["dataset_id"] if i in by_id]
                           for pid, grp in pmem[pmem["category"] == category].groupby("pool_id")}
        cat_mols = cat_flip = to_a = to_i = would_add = 0
        uniq_flip = set()
        for f in pool_files:
            pid = os.path.basename(f)[: -len(".csv.gz")]
            deliv = pd.read_csv(f, usecols=["inchikey", "bin"])
            dbin = dict(zip(deliv["inchikey"], deliv["bin"].astype(int)))
            members = members_by_pool.get(pid, [])
            if not members:
                continue
            rem = merge_datasets(members, meta, ik2info)
            rbin = dict(zip(rem["inchikey"], rem["bin"].astype(int)))
            nm = nf = fa = fi = 0
            for ik, db in dbin.items():
                rb = rbin.get(ik)
                if rb is None:
                    continue
                nm += 1
                if db != rb:
                    nf += 1; uniq_flip.add(ik)
                    if db == 0 and rb == 1:
                        fa += 1
                    else:
                        fi += 1
            n_add = len(set(rbin) - set(dbin))
            cat_mols += nm; cat_flip += nf; to_a += fa; to_i += fi; would_add += n_add
            detail_rows.append({"pathogen": code, "category": category, "pool": pid,
                                "n_molecules": nm, "n_members": len(members), "n_flip": nf,
                                "flip_rate": round(nf / nm, 4) if nm else np.nan,
                                "n_flip_to_active": fa, "n_flip_to_inactive": fi, "n_would_add": n_add})
        summary_rows.append({"pathogen": code, "category": category, "n_pools": len(members_by_pool),
                             "n_pool_molecules": cat_mols, "n_flip": cat_flip,
                             "pct_flip": round(cat_flip / cat_mols, 4) if cat_mols else np.nan,
                             "n_flip_to_active": to_a, "n_flip_to_inactive": to_i,
                             "n_unique_molecules_flip": len(uniq_flip), "n_would_add": would_add})
        print(f"  [reaverage {code}/{category}] {cat_mols} pool-molecules | flips {cat_flip} "
              f"({summary_rows[-1]['pct_flip']}) | would-add {would_add}")
    return detail_rows, summary_rows


def run_reaverage(codes):
    all_d, all_s = [], []
    for code in codes:
        print(f"[reaverage] {code}")
        d, s = _reaverage_pathogen(code); all_d += d; all_s += s
    detail, summary = pd.DataFrame(all_d), pd.DataFrame(all_s)
    detail.to_csv(_out("redundant-reaverage_detail.csv"), index=False)
    summary.to_csv(_out("redundant-reaverage_summary.csv"), index=False)
    if not summary.empty:
        mols = int(summary["n_pool_molecules"].sum()); flips = int(summary["n_flip"].sum())
        if mols:
            print(f"  overall would-flip if re-averaged: {flips} ({flips / mols:.2%} of {mols})")
        print(f"  distinct molecules flipping in >=1 pool: {int(summary['n_unique_molecules_flip'].sum())} | "
              f"would-add (not flips): {int(summary['n_would_add'].sum())}")
        pdf = _plot_frame(summary, "pct_flip")
        render_figure(docs_dir, f"{DATE}_redundant-reaverage", 1, 2, [
            lambda ax: _RateBars(ax, pdf, "re-average flip rate / %",
                                 "Pool labels that flip if step 25 re-averages", abc="A"),
            lambda ax: _DirBars(ax, pdf, "flipped pool molecules", "-> active", "-> inactive", abc="B"),
        ])
    print(f"  -> docs/{DATE}_redundant-reaverage_*.csv/.png")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
RUNNERS = {"conflicts": (run_conflicts, "25_merge_log.csv"),
           "merge-flip": (run_mergeflip, "25_merge_log.csv"),
           "reaverage": (run_reaverage, "25_pool_members.csv")}


def main():
    ap = argparse.ArgumentParser(description="Redundant-leftover vs pool label-conflict audit (read-only)")
    ap.add_argument("which", nargs="?", default="all",
                    choices=["all", "conflicts", "merge-flip", "reaverage"],
                    help="sub-analysis to run (default: all)")
    ap.add_argument("pathogens", nargs="*", help="pathogen codes (default: all on disk)")
    args = ap.parse_args()
    todo = list(RUNNERS) if args.which == "all" else [args.which]
    for name in todo:
        runner, needs = RUNNERS[name]
        print(f"\n===== {name} =====")
        runner(_codes(args.pathogens, needs))


if __name__ == "__main__":
    main()
