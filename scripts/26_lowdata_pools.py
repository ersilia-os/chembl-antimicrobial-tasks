"""
Step 26 (pooling rebuild) - low-data "catch-all" pools.

Some pathogen/categories are too small for steps 23-25 to build any modelled pool: no multi-dataset
transfer cluster forms and no single assay reaches the 200-molecule floor, so there is nothing to merge
leftovers into. For each such **deferred** category (0 modelled pools) this step pools ALL of the
category's datasets into ONE catch-all pool and cross-validates it. The catch-all is KEPT whenever it is
CV-able; its AUROC is EVALUATED and reported but NOT enforced (consistent with step 24, which models
regardless of AUROC), so the user decides downstream how to use a low-AUROC catch-all.

Per pathogen/category with 0 modelled pools (from 24_cv_summary; includes categories with no pools at all):
  1. Merge every 22_binarised dataset in the category into one pool (`merge_datasets`; conservative
     qualifier-aware re-binarise, 1 row/molecule).
  2. If >= CATCHALL_MIN_MOL molecules and >= MIN_PER_CLASS per class: 5-fold CV (RF/Morgan, plain) ->
     AUROC + Youden cutoff, and KEEP the pool. Otherwise not kept (too_small / not_cv_able).
  3. AUROC is reported for evaluation only; `auroc_pass` flags whether it clears CATCHALL_MIN_AUROC.
     CATCHALL_MIN_AUROC is NOT a keep criterion (just a reference line in the figures).

A catch-all merges heterogeneous assays that did not transfer to each other; its cross-validated AUROC
(reported, not enforced) indicates how much to trust it.

Inputs (output/stage4/<pathogen>/): 22_binarised/<cat>/*.csv.gz, 22_binarisation_summary.csv, 24_cv_summary.csv
Outputs:
    output/stage4/<pathogen>/26_pools/<cat>/<cat>_catchall.csv.gz   kept catch-all pool (1 row/molecule)
    output/stage4/<pathogen>/26_cv/<cat>/<cat>_catchall.csv.gz      OOF inchikey,y_true,y_prob (kept only)
    output/stage4/<pathogen>/26_cv_summary.csv                      per deferred category: AUROC, cutoff, kept, reason
    output/stage4/<pathogen>/26_catchall_overview.png              per-pathogen 1x3: catch-all AUROC (kept/
                                                                   dropped) | active vs inactive | datasets merged
    output/stage4/<pathogen>/26_eos1klk_pools_{DR,SP}.png          every FINAL pool (step-25 grown pools +
                                                                   step-26 catch-alls) painted on the eos1klk
                                                                   2D projections; one cyclic colour per pool,
                                                                   a single colour for molecules in >1 pool

The cross-pathogen catch-all overview figure is also drawn by step 27 (general plots).

Usage:
    python scripts/26_lowdata_pools.py                    # all pathogens (auto-detects deferred categories)
    python scripts/26_lowdata_pools.py ngonorrhoeae hpylori
    python scripts/26_lowdata_pools.py --plot-only        # redraw per-pathogen figures from cached 26_cv_summary
"""

import os
import sys
import glob
import argparse
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH
from pathogen_utils import load_pathogen
from dataset_utils import UNIT_SLUG, load_category_cutoffs
from pooling_utils import load_binarised_datasets, load_dataset_meta, merge_datasets, morgan_matrix
from model_utils import cv_oof_predict, youden_threshold
from plot_utils import cat_label, pathogen_coords, EMBEDDINGS
import plot_colors as pc
from plot_base import BasePlot, render_figure, inset_legend, compact_yaxis, figure_legend

# Format: print | Style: article - publication house style (shared across stage 4).
stylia.set_format("print")
stylia.set_style("article")

CATEGORIES = ["DR", "SP"]
CATCHALL_MIN_MOL = 50       # a catch-all pool must have at least this many molecules to be kept
MIN_PER_CLASS = 5           # ... and be CV-able (>= this many actives AND inactives)
CATCHALL_MIN_AUROC = 0.70   # evaluation reference only (NOT enforced): `auroc_pass` flags whether a
                            # catch-all clears it; a kept catch-all below it is still written
N_SPLITS = 5
output_base = os.path.join(root, "..", "output")


def deferred_categories(pathogen_code):
    """Categories that have binarised data but 0 modelled pools (from 24_cv_summary)."""
    cvp = os.path.join(output_base, "stage4", pathogen_code, "24_cv_summary.csv")
    cv = pd.read_csv(cvp) if os.path.exists(cvp) else pd.DataFrame(columns=["category", "modelled"])
    out = []
    for cat in CATEGORIES:
        if not glob.glob(os.path.join(output_base, "stage4", pathogen_code, "22_binarised", cat, "*.csv.gz")):
            continue
        nmod = int(((cv["category"] == cat) & cv["modelled"]).sum()) if len(cv) else 0
        if nmod == 0:
            out.append(cat)
    return out


def process(pathogen_code, cutoffs):
    out = os.path.join(output_base, "stage4", pathogen_code)
    # Clear stale catch-all pools up front (all categories): a category kept as a catch-all in a
    # previous run but now modelled is no longer deferred, so it is not rewritten below and would
    # otherwise leave a stale 26_pools/<cat> that step 27 double-counts.
    for cat in CATEGORIES:
        for f in glob.glob(os.path.join(out, "26_pools", cat, "*.csv.gz")):
            os.remove(f)
    cats = deferred_categories(pathogen_code)
    if not cats:
        return []
    print(f"Step 26 - {pathogen_code} ({load_pathogen(pathogen_code)}): deferred {cats}")
    rows = []
    for cat in cats:
        ds, ik2info = load_binarised_datasets(os.path.join(out, "22_binarised", cat))
        meta = load_dataset_meta(os.path.join(out, "22_binarisation_summary.csv"), cat, cutoffs, UNIT_SLUG)
        m = merge_datasets(ds, meta, ik2info)               # pool EVERYTHING in the category
        n = len(m); na = int((m.bin == 1).sum()); ni = int((m.bin == 0).sum())
        row = {"pathogen": pathogen_code, "category": cat, "n_datasets": len(ds), "n_molecules": n,
               "n_active": na, "n_inactive": ni, "active_ratio": round(na / n, 4) if n else np.nan,
               "auroc": np.nan, "youden_cutoff": np.nan, "sensitivity": np.nan, "specificity": np.nan,
               "auroc_pass": np.nan, "kept": False, "reason": ""}
        if n < CATCHALL_MIN_MOL:
            row["reason"] = "too_small"
        elif na < MIN_PER_CLASS or ni < MIN_PER_CLASS:
            row["reason"] = "not_cv_able"
        else:
            fps = morgan_matrix(dict(zip(m.inchikey, m.smiles)))
            df = m[m.inchikey.isin(fps)]
            # morgan_matrix drops invalid-SMILES molecules; recompute per-class counts on what CV
            # will actually see. A drop can push a class below N_SPLITS, so re-check CV-ability
            # here, otherwise StratifiedKFold(n_splits=N_SPLITS) would raise instead of being
            # flagged. Record the df-based counts (what is modelled).
            y = df.bin.to_numpy().astype(int)
            n = len(df); na = int((y == 1).sum()); ni = int((y == 0).sum())
            row.update(n_molecules=n, n_active=na, n_inactive=ni,
                       active_ratio=round(na / n, 4) if n else np.nan)
            if na < MIN_PER_CLASS or ni < MIN_PER_CLASS:
                row["reason"] = "not_cv_able"
            else:
                X = np.stack([fps[ik] for ik in df.inchikey])
                oof = cv_oof_predict(X, y, n_splits=N_SPLITS)
                auroc = float(roc_auc_score(y, oof))
                cutoff, sens, spec, _ = youden_threshold(y, oof)
                # AUROC is EVALUATED and reported, NOT enforced: keep every CV-able catch-all
                # (consistent with step 24, which models regardless of AUROC). `auroc_pass` records
                # whether it clears the CATCHALL_MIN_AUROC reference; the user decides how to use it.
                row.update(auroc=round(auroc, 4), youden_cutoff=round(cutoff, 4),
                           sensitivity=round(sens, 4), specificity=round(spec, 4),
                           auroc_pass=bool(auroc >= CATCHALL_MIN_AUROC), kept=True, reason="kept")
                pdir = os.path.join(out, "26_pools", cat); os.makedirs(pdir, exist_ok=True)
                m.to_csv(os.path.join(pdir, f"{cat}_catchall.csv.gz"), index=False)
                cdir = os.path.join(out, "26_cv", cat); os.makedirs(cdir, exist_ok=True)
                pd.DataFrame({"inchikey": df.inchikey.to_numpy(), "y_true": y, "y_prob": oof}).to_csv(
                    os.path.join(cdir, f"{cat}_catchall.csv.gz"), index=False)
        rows.append(row)
        print(f"    {cat}: {n} mol ({na}/{ni}) AUROC={row['auroc']} -> {row['reason'].upper()}")
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(out, "26_cv_summary.csv"), index=False)
        make_overview(pathogen_code, df)
    return rows


# --- per-pathogen catch-all figure (shared stage-4 plot_base framework) ------

def _auroc_color(auroc):
    """Colour a catch-all by its EVALUATION AUROC against the (non-enforced) reference:
    teal >= CATCHALL_MIN_AUROC, amber below it, grey if not evaluable (NaN)."""
    if pd.isna(auroc):
        return pc.named("silver")
    return pc.named("turquoise") if auroc >= CATCHALL_MIN_AUROC else pc.named("amber")


class CatchallAuroc(BasePlot):
    """Catch-all CV AUROC: one slot per deferred category; kept catch-alls show an AUROC bar
    (teal >= 0.70 / amber < 0.70), not-kept ones show their reason (too_small / not_cv_able)
    instead of a bar. Dashed line = the 0.70 evaluation reference (NOT enforced)."""

    def __init__(self, ax, df, abc=None):
        super().__init__(ax, abc)
        cats = list(df.category); x = np.arange(len(cats))
        for xi, r in zip(x, df.itertuples(index=False)):
            if pd.notna(r.auroc):
                ax.bar(xi, r.auroc, color=_auroc_color(r.auroc))
            else:
                ax.text(xi, 0.52, str(r.reason).replace("_", "\n"), ha="center", va="bottom",
                        fontsize=stylia.FONTSIZE_SMALL)
        ax.axhline(CATCHALL_MIN_AUROC, color=pc.color("neutral"), linestyle="--")
        ax.set_ylim(0.5, 1.0)
        ax.set_xlim(-0.5, len(cats) - 0.5)      # reserve a slot per category even when one is text-only
        ax.set_xticks(x); ax.set_xticklabels([cat_label(c) for c in cats])
        self._label(xlabel="", ylabel="catch-all CV AUROC", title="Catch-all AUROC (0.70 ref, not enforced)")
        inset_legend(ax, loc="upper right",
                     entries=[("≥ 0.70", pc.named("turquoise")), ("< 0.70", pc.named("amber"))])
        cols = [c for c in ["category", "auroc", "auroc_pass", "kept", "reason"] if c in df.columns]
        self.data = df[cols].copy()


class CatchallComposition(BasePlot):
    """Active vs inactive molecules in each candidate catch-all (before the keep decision)."""

    def __init__(self, ax, df, abc=None):
        super().__init__(ax, abc)
        cats = list(df.category); x = np.arange(len(cats))
        ni = df.n_inactive.to_numpy(dtype=float); na = df.n_active.to_numpy(dtype=float)
        ax.bar(x, ni, color=pc.color("inactive"))
        ax.bar(x, na, bottom=ni, color=pc.color("active"))
        compact_yaxis(ax)
        ax.set_xlim(-0.5, len(cats) - 0.5)
        ax.set_xticks(x); ax.set_xticklabels([cat_label(c) for c in cats])
        self._label(xlabel="", ylabel="molecules", title="Active vs inactive")
        inset_legend(ax, loc="upper right",
                     entries=[("inactive", pc.color("inactive")), ("active", pc.color("active"))])
        self.data = df[["category", "n_active", "n_inactive", "active_ratio"]].copy()


class CatchallDatasets(BasePlot):
    """Datasets merged into each catch-all pool; teal >= 0.70 / amber < 0.70 / grey not evaluable."""

    def __init__(self, ax, df, abc=None):
        super().__init__(ax, abc)
        cats = list(df.category); x = np.arange(len(cats))
        for xi, r in zip(x, df.itertuples(index=False)):
            ax.bar(xi, r.n_datasets, color=_auroc_color(r.auroc))
        compact_yaxis(ax)
        ax.set_xlim(-0.5, len(cats) - 0.5)
        ax.set_xticks(x); ax.set_xticklabels([cat_label(c) for c in cats])
        self._label(xlabel="", ylabel="datasets merged", title="Datasets in catch-all")
        self.data = df[["category", "n_datasets", "kept"]].copy()


def make_overview(pathogen_code, df):
    """Per-pathogen catch-all summary (1x3): AUROC (kept/dropped) | active vs inactive | datasets merged."""
    df = df.copy()
    df["kept"] = df["kept"].astype(str).str.lower().isin(["true", "1"])
    out = os.path.join(output_base, "stage4", pathogen_code)
    render_figure(out, "26_catchall_overview", 1, 3, [
        lambda ax: CatchallAuroc(ax, df, abc="A"),
        lambda ax: CatchallComposition(ax, df, abc="B"),
        lambda ax: CatchallDatasets(ax, df, abc="C"),
    ], width=0.9, height=0.32)


def make_pools_embedding(pathogen_code):
    """Paint every FINAL pool onto the eos1klk 2D projections, DR and SP as separate figures
    (4 projections each -> 8 panels total): 26_eos1klk_pools_DR.png / _SP.png.

    Final pools = the step-25 grown/kept pools (25_pools/<cat>/) plus the step-26 catch-alls
    (26_pools/<cat>/) for the same category. Each pool gets a distinct cyclic colour; molecules
    that belong to more than one pool share a single colour (ink); molecules in no final pool are
    grey. Pools are NOT disjoint in molecule space, hence the dedicated "in >1 pool" colour.
    Reads only pool files + eos1klk coordinates (no model training); skipped if the coordinates
    (step 12) are unavailable, or a category has no final pool."""
    out = os.path.join(output_base, "stage4", pathogen_code)
    cc = pathogen_coords(out)
    if cc is None or not len(cc):
        print(f"  [{pathogen_code}] no eos1klk coordinates - skipping pools embedding")
        return
    iks = cc["inchikey"].values
    for cat in CATEGORIES:
        files = (sorted(glob.glob(os.path.join(out, "25_pools", cat, "*.csv.gz")))
                 + sorted(glob.glob(os.path.join(out, "26_pools", cat, "*.csv.gz"))))
        pool_sets = []
        for f in files:
            s = set(pd.read_csv(f, usecols=["inchikey"])["inchikey"])
            if s:
                pool_sets.append(s)
        if not pool_sets:
            continue
        cnt = Counter()
        for s in pool_sets:
            cnt.update(s)
        membership = np.array([cnt.get(k, 0) for k in iks])
        cols = pc.cyclic(len(pool_sets))
        # draw order: grey background, then each pool's single-membership molecules, then the
        # shared (>1 pool) molecules on top so they are never hidden behind a pool colour.
        layers = [(membership == 0, pc.color("other"))]
        for i, s in enumerate(pool_sets):
            m = np.fromiter((k in s for k in iks), dtype=bool, count=len(iks)) & (membership == 1)
            layers.append((m, cols[i]))
        layers.append((membership >= 2, pc.color("ink")))

        fig, axs = stylia.create_figure(2, 2, width=0.8, height=0.8)
        for proj, xcol, ycol in EMBEDDINGS:
            ax = axs.next()
            for mask, color in layers:
                if mask.any():
                    ax.scatter(cc[xcol].values[mask], cc[ycol].values[mask], s=2, alpha=0.45,
                               color=color, edgecolors="none", rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            stylia.label(ax, xlabel="", ylabel="", title=proj)
        figure_legend(fig, [("not in a pool", pc.color("other")), ("in >1 pool", pc.color("ink"))],
                      loc="upper center", anchor=(0.5, -0.03), ncol=2,
                      title=f"{len(pool_sets)} final {cat_label(cat)} pools (one cyclic colour each)")
        path = os.path.join(out, f"26_eos1klk_pools_{cat}.png")
        stylia.save_figure(path)
        print(f"    embedding -> {path}")


def plot_only(pathogen_code):
    """Redraw 26_catchall_overview.png from the cached 26_cv_summary.csv (no CV re-run)."""
    p = os.path.join(output_base, "stage4", pathogen_code, "26_cv_summary.csv")
    if not os.path.exists(p):
        print(f"Step 26 (plot-only) - {pathogen_code}: no 26_cv_summary.csv - skipping")
        return
    make_overview(pathogen_code, pd.read_csv(p))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Step 26: low-data catch-all pools (per pathogen)")
    ap.add_argument("pathogens", nargs="*", help="default: all pathogens with 22_binarised data")
    ap.add_argument("--plot-only", action="store_true",
                    help="redraw 26_catchall_overview.png + 26_eos1klk_pools_*.png from cached pools/CV; skip the merge/CV")
    args = ap.parse_args()
    if args.plot_only:
        pooldirs = (glob.glob(os.path.join(output_base, "stage4", "*", "25_pools"))
                    + glob.glob(os.path.join(output_base, "stage4", "*", "26_pools")))
        codes = args.pathogens or sorted({os.path.basename(os.path.dirname(d)) for d in pooldirs})
        for p in codes:
            plot_only(p)               # 26_catchall_overview (skips if no 26_cv_summary)
            make_pools_embedding(p)    # 26_eos1klk_pools_{DR,SP}
        sys.exit(0)
    cutoffs = load_category_cutoffs(CONFIGPATH)
    pathogens = args.pathogens or sorted(
        os.path.basename(os.path.dirname(f)) for f in glob.glob(os.path.join(output_base, "stage4", "*", "22_binarised")))
    all_rows = []
    for p in pathogens:
        all_rows += process(p, cutoffs)
        make_pools_embedding(p)        # final-pool eos1klk embedding (step 25 + 26 pools)
    if all_rows:
        kept = sum(1 for r in all_rows if r["kept"])
        print(f"\nTotal deferred categories: {len(all_rows)}  kept: {kept}  dropped: {len(all_rows) - kept}")
    # The cross-pathogen catch-all overview figure is also drawn by step 27 (general plots).
