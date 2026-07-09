"""
Step 24 (pooling rebuild) - cross-validate the step-23 pools and score them.

For every pool built in step 23 (multi-dataset pools AND large single-assay singletons), fit a
RandomForest on Morgan fingerprints under 5-fold stratified cross-validation (80/20 per fold),
collect the out-of-fold (OOF) probability for each molecule, and from that single pooled ROC:
  - report the cross-validated AUROC, and
  - choose the best probability cutoff by the YOUDEN index (max sensitivity + specificity - 1).

Which pools are modelled ("modellable"):
  - ANY pool (multi-dataset OR singleton) with >= MIN_MODELLABLE_MOLECULES molecules AND CV-able
    (both classes, >= MIN_PER_CLASS per class). A single symmetric size floor is applied to multi-dataset
    pools and single-assay singletons alike; pools below it are deferred to step 25 for recovery.
(`chemspace_frac`, from 23_chemspace_partition.csv `n_total`, is reported for context but not gated on.)

Plain CV (no decoys/downsampling); RF(100) on Morgan FPs - the same model family as the transfer.
OOF probabilities are kept so AUROCs and actives/inactives probability distributions can be plotted.

Inputs (output/stage4/<pathogen>/):
    23_pool_summary.csv                      per-pool sizes, is_singleton
    23_pools/<cat>/<pool_id>.csv.gz          merged pools (inchikey, smiles, bin, ...)
    23_chemspace_partition.csv               n_total = category unique chem-space size

Outputs (output/stage4/<pathogen>/):
    24_cv/<cat>/<pool_id>.csv.gz             OOF predictions (inchikey, y_true, y_prob) per modelled pool
    24_cv_summary.csv                        per-pool: modelled flag, AUROC, Youden cutoff, sens/spec
    24_cv_overview.png                       per-pathogen 2x2 (narrow left / wide right): pool
                                             modelling status | per-pool AUROC with active/inactive
                                             pies; DR/SP AUROC boxplot | per-pool OOF violins+dots
    24_roc_curves.png                        per-pathogen grid: one ROC curve per modelled pool with
                                             the area under it filled (DR cobalt / SP tangerine)

Usage:
    python scripts/24_crossvalidate_pools.py                         # all pathogens with step-23 pools
    python scripts/24_crossvalidate_pools.py saureus mtuberculosis
    python scripts/24_crossvalidate_pools.py --plot-only saureus     # redraw the overview from cached CV

Thresholds are module-level constants below (MIN_MODELLABLE_MOLECULES, MIN_PER_CLASS, N_SPLITS);
edit them there rather than via the command line.
"""

import os
import sys
import glob
import math
import argparse

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import RANDOM_SEED
from pathogen_utils import load_pathogen
from pooling_utils import morgan_matrix
from model_utils import cv_oof_predict, youden_threshold
from plot_utils import cat_label, _xticks
import plot_colors as pc
from plot_base import BasePlot, render_figure, axis_legend, compact_yaxis

# Format: print | Style: article - publication house style (shared across stage 4).
stylia.set_format("print")
stylia.set_style("article")

CATEGORIES = ["DR", "SP"]
MIN_PER_CLASS = 5                       # CV-able floor: >= this many actives AND inactives
N_SPLITS = 5                            # 5-fold => 80/20 per fold
MIN_MODELLABLE_MOLECULES = 100          # symmetric size floor: a pool (multi-dataset OR singleton) is
                                        # modelled only if it has >= this many molecules; smaller pools are
                                        # deferred to step 25 for recovery (mirrors the PubChem branch)

output_base = os.path.join(root, "..", "output")


def chemspace_totals(pathogen_code):
    """{category: n_total} from 23_chemspace_partition.csv (union of all datasets' molecules)."""
    f = os.path.join(output_base, "stage4", pathogen_code, "23_chemspace_partition.csv")
    if not os.path.exists(f):
        return {}
    df = pd.read_csv(f)
    return {r.category: int(r.n_total) for r in df.itertuples(index=False)}


def process(pathogen_code):
    out = os.path.join(output_base, "stage4", pathogen_code)
    summ_path = os.path.join(out, "23_pool_summary.csv")
    if not os.path.exists(summ_path):
        print(f"Step 24 - {pathogen_code}: no 23_pool_summary.csv - skipping")
        return None
    print(f"Step 24 - {pathogen_code} ({load_pathogen(pathogen_code)})")
    summ = pd.read_csv(summ_path)
    chemspace = chemspace_totals(pathogen_code)

    rows, oof_by_pool = [], {}          # oof_by_pool[(category, pool_id)] = (y_true, y_prob)
    for category in CATEGORIES:
        sub = summ[summ["category"] == category]
        if sub.empty:
            continue
        total = chemspace.get(category, 0)
        # fingerprints once per category (union of all pool molecules)
        pool_dir = os.path.join(out, "23_pools", category)
        ik2smi = {}
        pool_dfs = {}
        for pid in sub["pool_id"]:
            df = pd.read_csv(os.path.join(pool_dir, pid + ".csv.gz"))
            pool_dfs[pid] = df
            for ik, smi in zip(df["inchikey"], df["smiles"]):
                ik2smi.setdefault(ik, smi)
        fps = morgan_matrix(ik2smi)

        cv_dir = os.path.join(out, "24_cv", category)
        os.makedirs(cv_dir, exist_ok=True)
        for f in glob.glob(os.path.join(cv_dir, "*.csv.gz")):
            os.remove(f)

        for r in sub.itertuples(index=False):
            df = pool_dfs[r.pool_id]
            df = df[df["inchikey"].isin(fps)]                    # keep only fingerprintable molecules
            n_mol = len(df)                                      # count actually modelled/cross-validated
            y = df["bin"].to_numpy().astype(int)
            na, ni = int((y == 1).sum()), int((y == 0).sum())
            frac = round(n_mol / total, 6) if total else np.nan
            include = n_mol >= MIN_MODELLABLE_MOLECULES          # symmetric floor on post-fingerprint count
            cv_able = na >= MIN_PER_CLASS and ni >= MIN_PER_CLASS
            modelled = bool(include and cv_able)
            row = {"pathogen": pathogen_code, "category": category, "pool_id": r.pool_id,
                   "is_singleton": bool(r.is_singleton), "n_molecules": n_mol,
                   "chemspace_frac": frac, "n_active": na, "n_inactive": ni,
                   "active_ratio": round(na / (na + ni), 5) if (na + ni) else np.nan,
                   "modelled": modelled, "auroc": np.nan, "youden_cutoff": np.nan,
                   "sensitivity": np.nan, "specificity": np.nan}
            if not modelled:
                row["reason"] = "below_min_molecules" if not include else "not_cv_able"
                rows.append(row); continue
            X = np.stack([fps[ik] for ik in df["inchikey"]])
            oof = cv_oof_predict(X, y, n_splits=N_SPLITS)
            auroc = float(roc_auc_score(y, oof))
            cutoff, sens, spec, _ = youden_threshold(y, oof)
            pd.DataFrame({"inchikey": df["inchikey"].to_numpy(), "y_true": y, "y_prob": oof}).to_csv(
                os.path.join(cv_dir, r.pool_id + ".csv.gz"), index=False)
            row.update(reason="modelled", auroc=round(auroc, 4), youden_cutoff=round(cutoff, 4),
                       sensitivity=round(sens, 4), specificity=round(spec, 4))
            rows.append(row)
            oof_by_pool[(category, r.pool_id)] = (y, oof)

        nmod = int(sum(1 for x in rows if x["category"] == category and x["modelled"]))
        print(f"  {category}: {len(sub)} pools -> {nmod} modelled "
              f"(molecules >= {MIN_MODELLABLE_MOLECULES})")

    if not rows:
        return None
    summ_df = pd.DataFrame(rows)
    summ_df.to_csv(os.path.join(out, "24_cv_summary.csv"), index=False)
    make_overview(pathogen_code, summ_df, oof_by_pool)
    make_roc_figure(pathogen_code, summ_df, oof_by_pool)
    return rows


# --- plots (shared stage-4 plot_base framework) -----------------------------
# Pool status buckets (reason -> label, colour). Shared by the status panel.
_STATUS = [("modelled", "modelled", "cobalt"),
           ("not_cv_able", "1 class / <5", "amber"),
           ("below_min_molecules", f"<{MIN_MODELLABLE_MOLECULES} mol", "silver")]


def _order_modelled(summ):
    """Modelled pools ordered DR-then-SP, each sorted by AUROC desc. Returns (rows, cat_spans):
    rows = list of dicts with x, category, pool_id, auroc, active_ratio, n_active, n_inactive,
    youden_cutoff; cat_spans = {category: (x0, x1)} for group labels/separators."""
    rows, spans, x = [], {}, 0
    for cat in CATEGORIES:
        sub = summ[(summ.category == cat) & summ.modelled].sort_values("auroc", ascending=False)
        if sub.empty:
            continue
        if rows:
            x += 1                                    # gap between DR and SP groups
        x0 = x
        for r in sub.itertuples(index=False):
            rows.append({"x": x, "category": cat, "pool_id": r.pool_id, "auroc": r.auroc,
                         "active_ratio": r.active_ratio, "n_active": r.n_active,
                         "n_inactive": r.n_inactive, "youden_cutoff": r.youden_cutoff,
                         "is_singleton": bool(r.is_singleton)})
            x += 1
        spans[cat] = (x0, x - 1)
    return rows, spans


def _group_xaxis(ax, spans):
    """Label DR/SP groups on the x-axis and draw a faint separator between them; no per-pool ticks."""
    ax.set_xticks([(a + b) / 2 for a, b in spans.values()])
    ax.set_xticklabels([cat_label(c) for c in spans])
    xs = sorted(b for _, b in spans.values())
    for b in xs[:-1]:
        ax.axvline(b + 0.5, color=pc.color("neutral"), linewidth=0.8, alpha=0.5)


def _pie_point(ax, x, y, frac_active, size=340):
    """Draw a small active/inactive pie glyph at (x, y) using vertex markers (round at any axis
    scale). frac_active in (0,1); active = crimson wedge, inactive = silver wedge."""
    for f0, f1, color in [(0.0, frac_active, pc.color("active")),
                          (frac_active, 1.0, pc.color("inactive"))]:
        if f1 - f0 <= 1e-9:
            continue
        ang = np.linspace(2 * np.pi * f0, 2 * np.pi * f1, 40)
        verts = [(0, 0)] + list(zip(np.cos(ang), np.sin(ang))) + [(0, 0)]
        ax.scatter([x], [y], marker=Path(verts), s=size, facecolor=color, edgecolors="none", zorder=3)


class PoolStatusPlot(BasePlot):
    """Per category: molecules split by pool status (modelled / one-class / small-singleton),
    stacked; each segment annotated with its pool count. Shows what is modelled vs excluded and
    how many molecules that represents."""

    def __init__(self, ax, summ, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if (summ.category == c).any()]
        if not cats:
            self._empty(ylabel="molecules", title="Pool modelling status"); return
        x = np.arange(len(cats))
        bottom = np.zeros(len(cats))
        recs = []
        for reason, label, ckey in _STATUS:
            mols = np.array([int(summ[(summ.category == c) & (summ.reason == reason)]["n_molecules"].sum())
                             for c in cats], dtype=float)
            cnts = [int(((summ.category == c) & (summ.reason == reason)).sum()) for c in cats]
            ax.bar(x, mols, bottom=bottom, color=pc.named(ckey))
            for xi, m, b, n in zip(x, mols, bottom, cnts):
                if m > 0:
                    ax.text(xi, b + m / 2, str(n), ha="center", va="center",
                            fontsize=stylia.FONTSIZE_SMALL)
            bottom += mols
            recs.append({"status": label, **{c: mols[i] for i, c in enumerate(cats)}})
        _xticks(ax, [cat_label(c) for c in cats])
        ax.set_ylim(0, bottom.max() * 1.08 if bottom.max() else 1)   # headroom for the count labels
        compact_yaxis(ax)
        self._label(xlabel="", ylabel="molecules (pool sum)", title="Pool modelling status")
        # legend sits at the top over the shorter bar (right if the rightmost bar is the
        # shorter, else left) so it never overlaps the taller stack's count annotations.
        if int(np.argmin(bottom)) >= len(bottom) / 2:
            loc, anchor = "upper right", (0.98, 0.98)
        else:
            loc, anchor = "upper left", (0.02, 0.98)
        axis_legend(ax, [(l, pc.named(ck)) for _, l, ck in _STATUS],
                    loc=loc, anchor=anchor, frame=True)
        self.data = pd.DataFrame(recs)


def _ring_color(is_singleton):
    """Ring colour marking a single large dataset (amber) vs a real multi-dataset pool (ink)."""
    return pc.named("amber") if is_singleton else pc.color("ink")


class PieAurocPlot(BasePlot):
    """Each modelled pool as an active/inactive pie glyph placed at its CV AUROC (y); the ring
    marks a single large dataset (amber) vs a real multi-dataset pool (ink). Pools ordered
    DR-then-SP by AUROC (same order as panel D). 0.7 = visual guide. Glyphs shrink when there
    are many pools so they do not overlap."""

    def __init__(self, ax, order_rows, spans, abc=None):
        super().__init__(ax, abc)
        if not order_rows:
            self._empty(xlabel="modelled pool", ylabel="CV AUROC", title="AUROC & active ratio"); return
        n = len(order_rows)
        scale = min(1.0, (14.0 / n) ** 2)      # area shrinks ~1/n^2 so radius tracks the x-spacing
        pie_s, ring_s = 340 * scale, 430 * scale
        for d in order_rows:
            ax.scatter([d["x"]], [d["auroc"]], s=ring_s, facecolor="none",
                       edgecolors=_ring_color(bool(d["is_singleton"])), linewidths=1.1, zorder=2)
            _pie_point(ax, d["x"], d["auroc"], d["active_ratio"], size=pie_s)
        ax.axhline(0.7, color=pc.color("neutral"), linestyle="--", zorder=1)
        ax.set_ylim(0.5, 1.02)
        ax.set_xlim(-0.7, max(d["x"] for d in order_rows) + 0.7)
        _group_xaxis(ax, spans)
        self._label(xlabel="", ylabel="CV AUROC",
                    title="AUROC & active ratio per pool")
        handles = [Patch(color=pc.color("active"), label="active"),
                   Patch(color=pc.color("inactive"), label="inactive"),
                   Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none",
                          markeredgecolor=_ring_color(False), label="pool (≥2 datasets)"),
                   Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none",
                          markeredgecolor=_ring_color(True), label="single dataset")]
        leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.02, 0.02), frameon=True)
        leg.get_frame().set_facecolor("white"); leg.get_frame().set_alpha(0.7)
        leg.get_frame().set_edgecolor("none")
        self.data = pd.DataFrame(order_rows)


class AurocBoxplot(BasePlot):
    """Per-pool CV AUROC as a DR and an SP boxplot with one dot per pool; 0.7 = visual guide."""

    def __init__(self, ax, summ, abc=None):
        super().__init__(ax, abc)
        cats = [c for c in CATEGORIES if ((summ.category == c) & summ.modelled).any()]
        if not cats:
            self._empty(ylabel="CV AUROC", title="Pool AUROC"); return
        data = [summ[(summ.category == c) & summ.modelled]["auroc"].dropna().to_numpy() for c in cats]
        pos = list(range(len(cats)))
        bp = ax.boxplot(data, positions=pos, widths=0.5, showfliers=False, patch_artist=True)
        for patch, c in zip(bp["boxes"], cats):
            patch.set_facecolor("none"); patch.set_edgecolor(pc.color(c))
        for i, c in enumerate(cats):
            for art in bp["whiskers"][2 * i:2 * i + 2] + bp["caps"][2 * i:2 * i + 2]:
                art.set_color(pc.color(c))
        for med in bp["medians"]:
            med.set_color(pc.color("ink"))
        rng = np.random.default_rng(RANDOM_SEED)
        for i, c in enumerate(cats):
            a = data[i]
            ax.scatter(i + rng.uniform(-0.15, 0.15, size=len(a)), a, color=pc.color(c),
                       alpha=0.6, edgecolors="none", zorder=3)
        ax.axhline(0.7, color=pc.color("neutral"), linestyle="--", zorder=1)
        ax.set_xticks(pos); ax.set_xticklabels([cat_label(c) for c in cats])
        ax.set_ylim(0.5, 1.0)
        self._label(xlabel="", ylabel="CV AUROC", title="Pool AUROC")
        self.data = pd.DataFrame({"category": np.concatenate([[c] * len(d) for c, d in zip(cats, data)]),
                                  "auroc": np.concatenate(data) if data else np.array([])})


class RocOne(BasePlot):
    """One modelled pool's ROC curve with the area under it filled (AUC), category-coloured;
    dashed diagonal = chance. AUROC value annotated. Used as a small-multiple cell."""

    def __init__(self, ax, d, oof_by_pool, abc=None):
        super().__init__(ax, abc)
        y, p = oof_by_pool[(d["category"], d["pool_id"])]
        fpr, tpr, _ = roc_curve(y, p)
        color = pc.color(d["category"])
        ax.fill_between(fpr, 0, tpr, color=color, alpha=0.3)
        ax.plot(fpr, tpr, color=color, linewidth=1.2)
        ax.plot([0, 1], [0, 1], color=pc.color("neutral"), linestyle="--", linewidth=0.6)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.text(0.94, 0.06, f"{d['auroc']:.2f}", ha="right", va="bottom",
                fontsize=stylia.FONTSIZE_SMALL)
        self._label(xlabel="", ylabel="", title=d["pool_id"])


class PoolViolinPlot(BasePlot):
    """Per modelled pool: split violins of the OOF P(active) for inactives (silver, left) vs actives
    (crimson, right), with that pool's Youden cutoff marked. DR-then-SP-by-AUROC order."""

    def __init__(self, ax, order_rows, oof_by_pool, spans, abc=None):
        super().__init__(ax, abc)
        rows = [d for d in order_rows if (d["category"], d["pool_id"]) in oof_by_pool]
        if not rows:
            self._empty(ylabel="OOF P(active)", title="OOF predictions per pool"); return
        rng = np.random.default_rng(RANDOM_SEED)
        for d in rows:
            y, p = oof_by_pool[(d["category"], d["pool_id"])]
            for side, mask, color in [(-1, y == 0, pc.color("inactive")), (1, y == 1, pc.color("active"))]:
                vals = p[mask]
                if len(vals) < 2:
                    continue
                vp = ax.violinplot([vals], positions=[d["x"] + side * 0.22], widths=0.36,
                                   showextrema=False)
                for body in vp["bodies"]:
                    body.set_facecolor(color); body.set_alpha(0.35); body.set_edgecolor("none")
                # zairachem-style individual dots (subsampled) over the violin
                idx = np.arange(len(vals))
                if len(idx) > 60:
                    idx = rng.choice(idx, 60, replace=False)
                ax.scatter(d["x"] + side * 0.22 + rng.uniform(-0.06, 0.06, size=len(idx)),
                           vals[idx], s=4, color=color, alpha=0.6, edgecolors="none", zorder=3)
            ax.hlines(d["youden_cutoff"], d["x"] - 0.42, d["x"] + 0.42,
                      color=pc.color("ink"), linewidth=1.0, zorder=4)
        ax.set_ylim(-0.06, 1.06)                 # headroom so the violins are never cut at 0/1
        ax.set_xlim(-0.7, max(d["x"] for d in rows) + 0.7)
        _group_xaxis(ax, spans)
        self._label(xlabel="modelled pool (DR then SP, by AUROC)", ylabel="OOF P(active)",
                    title="OOF predictions per pool (line = Youden cutoff)")
        axis_legend(ax, [("inactive", pc.color("inactive")), ("active", pc.color("active"))],
                    loc="upper right", anchor=(0.98, 0.98), frame=True)


def make_overview(pathogen_code, summ, oof_by_pool):
    """2x2, narrow left column + wide right column (width_ratios 1:3):
    A pool modelling status (small) | B per-pool AUROC with active/inactive pies (wide)
    C DR/SP AUROC boxplot + dots (small) | D per-pool OOF violins+dots with Youden cutoffs (wide)."""
    out = os.path.join(output_base, "stage4", pathogen_code)
    order_rows, spans = _order_modelled(summ)

    render_figure(out, "24_cv_overview", 2, 2, [
        lambda ax: PoolStatusPlot(ax, summ, abc="A"),      # status key sits inside panel A
        lambda ax: PieAurocPlot(ax, order_rows, spans, abc="B"),
        lambda ax: AurocBoxplot(ax, summ, abc="C"),
        lambda ax: PoolViolinPlot(ax, order_rows, oof_by_pool, spans, abc="D"),
    ], width_ratios=[1, 3], height=1.0)


def make_roc_figure(pathogen_code, summ, oof_by_pool):
    """Separate figure: one small ROC curve per modelled pool (area under the curve filled,
    DR cobalt / SP tangerine), laid out as a grid (DR-then-SP by AUROC)."""
    out = os.path.join(output_base, "stage4", pathogen_code)
    order_rows, _ = _order_modelled(summ)
    rows = [d for d in order_rows if (d["category"], d["pool_id"]) in oof_by_pool]
    if not rows:
        print(f"  [roc] {pathogen_code}: no modelled pools - skipping")
        return
    ncols = min(6, len(rows))
    nrows = math.ceil(len(rows) / ncols)
    builders = [(lambda ax, d=d: RocOne(ax, d, oof_by_pool)) for d in rows]
    builders += [(lambda ax: ax.axis("off"))] * (nrows * ncols - len(rows))
    render_figure(out, "24_roc_curves", nrows, ncols, builders, height=1.0 * nrows / ncols)


def plot_only(pathogen_code):
    """Redraw 24_cv_overview.png from the cached CV outputs (24_cv_summary.csv + the
    per-pool OOF predictions in 24_cv/<cat>/), without re-running cross-validation."""
    out = os.path.join(output_base, "stage4", pathogen_code)
    summ_path = os.path.join(out, "24_cv_summary.csv")
    if not os.path.exists(summ_path):
        print(f"Step 24 (plot-only) - {pathogen_code}: no 24_cv_summary.csv - skipping")
        return
    summ = pd.read_csv(summ_path)
    oof_by_pool = {}
    for cat in CATEGORIES:
        for f in sorted(glob.glob(os.path.join(out, "24_cv", cat, "*.csv.gz"))):
            pid = os.path.basename(f)[:-len(".csv.gz")]
            d = pd.read_csv(f)
            oof_by_pool[(cat, pid)] = (d["y_true"].to_numpy(), d["y_prob"].to_numpy())
    make_overview(pathogen_code, summ, oof_by_pool)
    make_roc_figure(pathogen_code, summ, oof_by_pool)


def _codes_with(artifact):
    return sorted(os.path.basename(os.path.dirname(f))
                  for f in glob.glob(os.path.join(output_base, "stage4", "*", artifact)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Step 24: cross-validate pools + Youden cutoff")
    ap.add_argument("pathogens", nargs="*", help="default: all pathogens with 23_pool_summary.csv")
    ap.add_argument("--plot-only", action="store_true",
                    help="redraw 24_cv_overview.png from cached CV outputs; skip cross-validation")
    args = ap.parse_args()
    if args.plot_only:
        for p in (args.pathogens or _codes_with("24_cv_summary.csv")):
            plot_only(p)
    else:
        for p in (args.pathogens or _codes_with("23_pool_summary.csv")):
            process(p)
