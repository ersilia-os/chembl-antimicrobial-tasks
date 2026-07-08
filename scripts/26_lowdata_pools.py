"""
Step 26 (pooling rebuild) - low-data "catch-all" pools.

Some pathogen/categories are too small for steps 23-25 to build any modelled pool: no multi-dataset
transfer cluster forms and no single assay reaches the 200-molecule floor, so there is nothing to merge
leftovers into. For each such **deferred** category (0 modelled pools) this step pools ALL of the
category's datasets into ONE catch-all pool, cross-validates it, and keeps it as a real model only if it
clears a quality bar - otherwise the category is documented as insufficient data (no model).

Per pathogen/category with 0 modelled pools (from 24_cv_summary; includes categories with no pools at all):
  1. Merge every 22_binarised dataset in the category into one pool (`merge_datasets`; conservative
     qualifier-aware re-binarise, 1 row/molecule).
  2. If >= CATCHALL_MIN_MOL molecules and >= MIN_PER_CLASS per class: 5-fold CV (RF/Morgan, plain) ->
     AUROC and the Youden probability cutoff.
  3. KEEP iff AUROC >= CATCHALL_MIN_AUROC; else drop (low_auroc / too_small / not_cv_able).

A catch-all merges heterogeneous assays that did not transfer to each other, so it is only trusted when
its cross-validated AUROC actually clears the bar.

Inputs (output/stage4/<pathogen>/): 22_binarised/<cat>/*.csv.gz, 22_binarisation_summary.csv, 24_cv_summary.csv
Outputs:
    output/stage4/<pathogen>/26_pools/<cat>/<cat>_catchall.csv.gz   kept catch-all pool (1 row/molecule)
    output/stage4/<pathogen>/26_cv/<cat>/<cat>_catchall.csv.gz      OOF inchikey,y_true,y_prob (kept only)
    output/stage4/<pathogen>/26_cv_summary.csv                      per deferred category: AUROC, cutoff, kept, reason
    output/stage4/<pathogen>/26_catchall_overview.png              per-pathogen 1x3: catch-all AUROC (kept/
                                                                   dropped) | active vs inactive | datasets merged

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
from plot_utils import cat_label
import plot_colors as pc
from plot_base import BasePlot, render_figure, inset_legend, compact_yaxis

# Format: print | Style: article - publication house style (shared across stage 4).
stylia.set_format("print")
stylia.set_style("article")

CATEGORIES = ["DR", "SP"]
CATCHALL_MIN_MOL = 50       # a catch-all pool must have at least this many molecules
CATCHALL_MIN_AUROC = 0.70   # ... and clear this CV AUROC to be kept as a model
MIN_PER_CLASS = 5           # ... and be CV-able
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
    cats = deferred_categories(pathogen_code)
    if not cats:
        return []
    print(f"Step 26 - {pathogen_code} ({load_pathogen(pathogen_code)}): deferred {cats}")
    out = os.path.join(output_base, "stage4", pathogen_code)
    rows = []
    for cat in cats:
        ds, ik2info = load_binarised_datasets(os.path.join(out, "22_binarised", cat))
        meta = load_dataset_meta(os.path.join(out, "22_binarisation_summary.csv"), cat, cutoffs, UNIT_SLUG)
        m = merge_datasets(ds, meta, ik2info)               # pool EVERYTHING in the category
        n = len(m); na = int((m.bin == 1).sum()); ni = int((m.bin == 0).sum())
        row = {"pathogen": pathogen_code, "category": cat, "n_datasets": len(ds), "n_molecules": n,
               "n_active": na, "n_inactive": ni, "active_ratio": round(na / n, 4) if n else np.nan,
               "auroc": np.nan, "youden_cutoff": np.nan, "sensitivity": np.nan, "specificity": np.nan,
               "kept": False, "reason": ""}
        if n < CATCHALL_MIN_MOL:
            row["reason"] = "too_small"
        elif na < MIN_PER_CLASS or ni < MIN_PER_CLASS:
            row["reason"] = "not_cv_able"
        else:
            fps = morgan_matrix(dict(zip(m.inchikey, m.smiles)))
            df = m[m.inchikey.isin(fps)]
            X = np.stack([fps[ik] for ik in df.inchikey]); y = df.bin.to_numpy().astype(int)
            oof = cv_oof_predict(X, y, n_splits=N_SPLITS)
            auroc = float(roc_auc_score(y, oof))
            cutoff, sens, spec, _ = youden_threshold(y, oof)
            row.update(auroc=round(auroc, 4), youden_cutoff=round(cutoff, 4),
                       sensitivity=round(sens, 4), specificity=round(spec, 4))
            if auroc >= CATCHALL_MIN_AUROC:
                row["kept"] = True; row["reason"] = "kept"
                pdir = os.path.join(out, "26_pools", cat); os.makedirs(pdir, exist_ok=True)
                m.to_csv(os.path.join(pdir, f"{cat}_catchall.csv.gz"), index=False)
                cdir = os.path.join(out, "26_cv", cat); os.makedirs(cdir, exist_ok=True)
                pd.DataFrame({"inchikey": df.inchikey.to_numpy(), "y_true": y, "y_prob": oof}).to_csv(
                    os.path.join(cdir, f"{cat}_catchall.csv.gz"), index=False)
            else:
                row["reason"] = "low_auroc"
        rows.append(row)
        print(f"    {cat}: {n} mol ({na}/{ni}) AUROC={row['auroc']} -> {row['reason'].upper()}")
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(out, "26_cv_summary.csv"), index=False)
        make_overview(pathogen_code, df)
    return rows


# --- per-pathogen catch-all figure (shared stage-4 plot_base framework) ------

def _kept_color(kept):
    return pc.named("turquoise") if kept else pc.named("silver")


class CatchallAuroc(BasePlot):
    """Catch-all CV AUROC per deferred category: teal kept / grey dropped; dashed line = keep bar
    (CATCHALL_MIN_AUROC). Categories with no AUROC (too small / not CV-able) show their reason."""

    def __init__(self, ax, df, abc=None):
        super().__init__(ax, abc)
        cats = list(df.category); x = np.arange(len(cats))
        for xi, r in zip(x, df.itertuples(index=False)):
            if pd.notna(r.auroc):
                ax.bar(xi, r.auroc, color=_kept_color(r.kept))
            else:
                ax.text(xi, 0.52, str(r.reason).replace("_", "\n"), ha="center", va="bottom",
                        fontsize=stylia.FONTSIZE_SMALL)
        ax.axhline(CATCHALL_MIN_AUROC, color=pc.color("neutral"), linestyle="--")
        ax.set_ylim(0.5, 1.0)
        ax.set_xlim(-0.5, len(cats) - 0.5)      # reserve a slot per category even when one is text-only
        ax.set_xticks(x); ax.set_xticklabels([cat_label(c) for c in cats])
        self._label(xlabel="", ylabel="catch-all CV AUROC", title="Catch-all AUROC")
        inset_legend(ax, loc="upper right",
                     entries=[("kept", pc.named("turquoise")), ("dropped", pc.named("silver"))])
        self.data = df[["category", "auroc", "kept", "reason"]].copy()


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
    """Datasets merged into each catch-all pool; teal kept / grey dropped."""

    def __init__(self, ax, df, abc=None):
        super().__init__(ax, abc)
        cats = list(df.category); x = np.arange(len(cats))
        for xi, r in zip(x, df.itertuples(index=False)):
            ax.bar(xi, r.n_datasets, color=_kept_color(r.kept))
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
                    help="redraw 26_catchall_overview.png from cached 26_cv_summary.csv; skip the merge/CV")
    args = ap.parse_args()
    if args.plot_only:
        codes = args.pathogens or sorted(os.path.basename(os.path.dirname(f))
                                         for f in glob.glob(os.path.join(output_base, "stage4", "*", "26_cv_summary.csv")))
        for p in codes:
            plot_only(p)
        sys.exit(0)
    cutoffs = load_category_cutoffs(CONFIGPATH)
    pathogens = args.pathogens or sorted(
        os.path.basename(os.path.dirname(f)) for f in glob.glob(os.path.join(output_base, "stage4", "*", "22_binarised")))
    all_rows = []
    for p in pathogens:
        all_rows += process(p, cutoffs)
    if all_rows:
        kept = sum(1 for r in all_rows if r["kept"])
        print(f"\nTotal deferred categories: {len(all_rows)}  kept: {kept}  dropped: {len(all_rows) - kept}")
    # The cross-pathogen catch-all overview figure is also drawn by step 27 (general plots).
