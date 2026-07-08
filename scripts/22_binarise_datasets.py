"""
Step 22 - Binarise curated datasets and plot the results.

First turns the retained, quantitative datasets from step 21 into binary
active/inactive labels, then renders the per-pathogen binarisation figures.
Runs on ONE pathogen per invocation (like steps 07-20); loop with scripts/run_step.sh.
Cross-pathogen comparison figures are a separate step (23).

Binarisation uses a single fixed cutoff keyed by unit, loaded from
`config/category_cutoffs.csv` (umol.L-1 -> 10, % -> 50, mm -> 20).

The cutoff *value* is keyed by unit; the *direction* (which side is active) is taken
from each dataset's own `direction` column (the curated authority), so the rare SP
assays where the low end is active (direction = -1) are handled correctly.

Binarisation is conservative and qualifier-aware ("prove-the-hit"), via
`dataset_utils.binarize_conservative`: a compound is labelled active (1) only when
its measurement PROVES it lies on the active side of the cutoff -
  direction -1: relation in {'=','<'} and value <= cutoff;
  direction +1: relation in {'=','>'} and value >= cutoff.
Censored records whose true value could straddle the cutoff (e.g. a `MIC > X`
non-inhibitor with X < cutoff) are labelled inactive rather than guessed - they are
kept, not dropped. Step 21 preserves the original qualifier so this is possible; the
`relation` column is carried into the binarised output for transparency.

Recovered qualitative-inactive molecules (step 21) arrive here with value/relation =
NaN; having no value to prove a hit, they are always labelled inactive. Their count
per dataset is reported as `n_qual_inactive` in the summary.

Only ORGANISM datasets with category DR or SP are binarised; UNCATEGORIZED datasets
(no defined cutoff) are skipped. Degenerate datasets (all-active or all-inactive) are
kept and flagged, not dropped.

Plots (via stylia, publication house style: print format + article palette, applied by
src/plot_base.py; colours from the shared semantic map in src/plot_colors.py so DR/SP
and active/inactive stay consistent with step 21), per pathogen
(`output/stage4/<pathogen>/`):
    22_binarisation_overview.png    per-category active-ratio histograms, active ratio
                                    vs dataset size, class composition by datasets and
                                    by molecules.
    22_eos1klk_actives.png, 22_eos1klk_category.png, 22_eos1klk_clusters.png
                                    from the eos1klk (lazychemvis) 2D coordinates:
                                    the four PCA/t-MAP/t-SNE/UMAP projections coloured
                                    by active/inactive, DR/shared/SP/other, and the
                                    largest compatibility-group pools.
    22_eos1klk_*.csv                companion "ready datasets": the plotted points
                                    (inchikey, per-projection x/y, group) for each of
                                    the three embedding figures.

Active-call cut-offs: DR <= 10 uM; SP >= 50 % or >= 20 mm.

Inputs (output/stage4/<pathogen>/):
    21_curation_summary.csv         retained datasets + category/unit/direction
    21_datasets/<assay>_<type>.csv.gz

Outputs (output/stage4/<pathogen>/):
    22_binarised/DR/<assay>_<type>_<slug>.csv.gz   (inchikey, compound_chembl_id,
    22_binarised/SP/<assay>_<type>_<slug>.csv.gz    smiles, value, relation, bin)
    22_binarisation_summary.csv                     one row per binarised dataset
    22_binarisation_overview.png, 22_eos1klk_*.png  per-pathogen figures
    22_eos1klk_*.csv                                companion point CSVs for the eos figures

Usage:
    python scripts/22_binarise_datasets.py <pathogen_code>
"""

import os
import sys
import glob

import numpy as np
import pandas as pd
from tqdm import tqdm

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH, STAGE4PATH
from pathogen_utils import load_pathogen
from dataset_utils import load_category_cutoffs, UNIT_SLUG, binarize_conservative
from plot_utils import (
    binarisation_stats, composition_counts,
    pathogen_coords,
    _xticks,
    CATS, CAT_DISPLAY, COMPOSITION, CUT,
    RANDOM_SEED, CAP_POINTS,
)
import plot_colors as pc
from plot_base import BasePlot, render_figure, EmbeddingFigure, compact_yaxis

# --- binarisation -----------------------------------------------------------

def binarise(pathogen_code):
    pathogen = load_pathogen(pathogen_code)

    OUTPUT = os.path.join(STAGE4PATH, pathogen_code)
    DATASETS_DIR = os.path.join(OUTPUT, "21_datasets")
    BIN_DIR = os.path.join(OUTPUT, "22_binarised")

    print(f"Step 22 - {pathogen_code} ({pathogen})")

    cutoffs = load_category_cutoffs(CONFIGPATH)
    print(f"  cutoffs by unit: {cutoffs}")

    summary = pd.read_csv(os.path.join(OUTPUT, "21_curation_summary.csv"), low_memory=False)
    to_bin = summary[(~summary["discarded"]) & (summary["category"].isin(["DR", "SP"]))].reset_index(drop=True)
    print(f"  datasets to binarise (DR/SP retained): {len(to_bin)} "
          f"(skipping {((~summary['discarded']) & (summary['category'] == 'UNCATEGORIZED')).sum()} UNCATEGORIZED)")

    # Validate the whole batch UP FRONT (unit has a cutoff, direction is ±1) before
    # touching any output, so a config gap can never leave a half-written 22_binarised/
    # with no summary. Both columns come from the step-21 summary.
    missing_units = sorted(set(to_bin["unit"]) - set(cutoffs))
    if missing_units:
        raise KeyError(
            f"no binarisation cutoff for unit(s) {missing_units} "
            f"(add to config/category_cutoffs.csv; configured: {sorted(cutoffs)})")
    bad_dir = to_bin[~to_bin["direction"].isin([-1, 1])]
    if len(bad_dir):
        raise ValueError(
            f"{len(bad_dir)} DR/SP dataset(s) with direction not in {{-1, +1}} "
            f"(e.g. {bad_dir.iloc[0]['assay_id']} {bad_dir.iloc[0]['activity_type']} "
            f"-> {bad_dir.iloc[0]['direction']!r}); step 21 should have discarded these.")

    # Clear prior outputs so the directory always matches the current run (e.g. when an
    # endpoint is reclassified out of DR/SP, its stale file must go). Done only after
    # validation passes.
    for cat in ("DR", "SP"):
        d = os.path.join(BIN_DIR, cat)
        os.makedirs(d, exist_ok=True)
        for f in glob.glob(os.path.join(d, "*.csv.gz")):
            os.remove(f)

    out_cols = ["inchikey", "compound_chembl_id", "smiles", "value", "relation", "bin"]
    rows = []
    for r in tqdm(to_bin.itertuples(index=False), total=len(to_bin), desc="  binarising"):
        assay_id, activity_type, unit, category = r.assay_id, r.activity_type, r.unit, r.category
        cutoff = cutoffs[unit]
        src = os.path.join(DATASETS_DIR, f"{assay_id}_{str(activity_type).replace('/', '-')}.csv.gz")
        df = pd.read_csv(src)
        direction_raw = df["direction"].iloc[0]
        if pd.isna(direction_raw) or int(direction_raw) not in (-1, 1):
            raise ValueError(
                f"dataset {assay_id} {activity_type}: direction must be -1 or +1, got {direction_raw!r}")
        direction = int(direction_raw)

        df = binarize_conservative(df, cutoff, direction)

        n = len(df)
        n_active = int((df["bin"] == 1).sum())
        n_inactive = n - n_active
        n_qual_inactive = int(df["value"].isna().sum())  # recovered qualitative inactives (no value)
        active_ratio = round(n_active / n, 5) if n else float("nan")
        degenerate = n_active == 0 or n_inactive == 0

        slug = UNIT_SLUG.get(unit, str(unit))
        fname = f"{assay_id}_{str(activity_type).replace('/', '-')}_{slug}.csv.gz"
        df[out_cols].to_csv(os.path.join(BIN_DIR, category, fname), index=False)

        rows.append([assay_id, activity_type, unit, category, cutoff, direction,
                     n, n_active, n_inactive, n_qual_inactive, active_ratio, degenerate])

    summary_cols = ["assay_id", "activity_type", "unit", "category", "cutoff", "direction",
                    "n_molecules", "n_active", "n_inactive", "n_qual_inactive", "active_ratio", "degenerate"]
    out = pd.DataFrame(rows, columns=summary_cols).sort_values(
        ["category", "n_molecules"], ascending=[True, False]).reset_index(drop=True)
    out.to_csv(os.path.join(OUTPUT, "22_binarisation_summary.csv"), index=False)

    n_dr = int((out["category"] == "DR").sum())
    n_sp = int((out["category"] == "SP").sum())
    n_deg = int(out["degenerate"].sum())
    print(f"  binarised: {len(out)} (DR {n_dr}, SP {n_sp}) | degenerate: {n_deg} "
          f"({100 * n_deg / len(out):.1f}%)" if len(out) else "  no datasets binarised")
    print(f"  summary -> {os.path.join(OUTPUT, '22_binarisation_summary.csv')}")


# --- per-pathogen panels ----------------------------------------------------
# DR = cobalt, SP = tangerine everywhere (plot_colors); the DR/SP key is a single
# shared legend placed outside the data by process(), so no per-axis legends here.

class RatioHistOne(BasePlot):
    """Distribution of per-dataset active ratio for a single category (own panel, so
    the smaller category is never hidden behind the larger one)."""

    def __init__(self, ax, ratios, color_key, title, abc=None):
        super().__init__(ax, abc)
        r = np.asarray([x for x in ratios if x is not None and not pd.isna(x)], dtype=float)
        if len(r):
            ax.hist(r, bins=np.linspace(0, 1, 21), color=pc.color(color_key))
        ax.set_xlim(0, 1)
        compact_yaxis(ax)
        self._label(xlabel="active ratio", ylabel="datasets", title=title)


class RatioVsSizePlot(BasePlot):
    """Active ratio (y) vs dataset size (x, log) for DR and SP together: faint
    subsampled per-dataset dots + a bold binned-mean trend line per category. Shows
    whether small datasets skew toward the degenerate extremes (0 / 1)."""

    def __init__(self, ax, df, abc=None):
        super().__init__(ax, abc)
        rng = np.random.default_rng(RANDOM_SEED)
        sizes = df["n_molecules"].to_numpy(dtype=float)
        if not len(sizes) or np.nanmax(sizes) <= 0:
            self._empty(xlabel="molecules / dataset", ylabel="active ratio",
                        title="Active ratio vs dataset size")
            return
        bins = np.logspace(np.log10(max(np.nanmin(sizes), 1)), np.log10(np.nanmax(sizes)), 12)
        centers = np.sqrt(bins[:-1] * bins[1:])
        for cat in CATS:
            color = pc.color(cat)
            sub = df[df["category"] == cat]
            if not len(sub):
                continue
            sz = sub["n_molecules"].to_numpy(dtype=float)
            r = sub["active_ratio"].to_numpy(dtype=float)
            # faint dots (subsampled for readability), then a bold binned-mean trend
            idx = np.arange(len(sz))
            if len(idx) > CAP_POINTS:
                idx = rng.choice(idx, size=CAP_POINTS, replace=False)
            ax.scatter(sz[idx], r[idx], color=color, alpha=0.15, edgecolors="none", s=6)
            bi = np.digitize(sz, bins)
            means = np.array([r[bi == k].mean() if (bi == k).any() else np.nan
                              for k in range(1, len(bins))])
            m = np.isfinite(means)
            ax.plot(centers[m], means[m], color=color, marker="o", label=CAT_DISPLAY[cat])
        ax.set_xscale("log")
        ax.set_ylim(-0.02, 1.02)
        ax.legend(loc="upper right")   # DR/SP key inside the panel
        self._label(xlabel="molecules / dataset", ylabel="active ratio",
                    title="Active ratio vs dataset size")


class CompositionPlot(BasePlot):
    """Stacked all-inactive / mixed / all-active per category, weighted either by
    dataset count or by molecule count depending on the supplied `comp`."""

    def __init__(self, ax, comp, ylabel, title, abc=None):
        super().__init__(ax, abc)
        x = np.arange(len(CATS))
        bottoms = np.zeros(len(CATS))
        for i, label in enumerate(COMPOSITION):
            vals = np.array([comp[c][i] for c in CATS], dtype=float)
            ax.bar(x, vals, bottom=bottoms, color=pc.color(label), label=label)
            bottoms += vals
        _xticks(ax, [f"{CAT_DISPLAY[c]}\n({CUT[c]})" for c in CATS])
        ax.legend(loc="upper right")   # composition key inside the panel
        compact_yaxis(ax)
        self._label(xlabel="", ylabel=ylabel, title=title)


def composition_mols(df):
    """Per-category [all-inactive, mixed, all-active] MOLECULE sums (the molecule-weighted
    counterpart of composition_counts, which counts datasets)."""
    out = {}
    for c in CATS:
        sub = df[df["category"] == c]
        out[c] = [int(sub.loc[sub["active_ratio"] == 0, "n_molecules"].sum()),
                  int(sub.loc[(sub["active_ratio"] > 0) & (sub["active_ratio"] < 1), "n_molecules"].sum()),
                  int(sub.loc[sub["active_ratio"] == 1, "n_molecules"].sum())]
    return out


class ActiveMoleculesPlot(BasePlot):
    """Unique molecules tested per category, stacked active (in >=1 assay) vs inactive.
    Ties the binarisation back to how many distinct actives it yields. `counts` =
    {cat: {"active": n, "inactive": n}}."""

    def __init__(self, ax, counts, abc=None):
        super().__init__(ax, abc)
        x = np.arange(len(CATS))
        act = np.array([counts[c]["active"] for c in CATS], dtype=float)
        ina = np.array([counts[c]["inactive"] for c in CATS], dtype=float)
        ax.bar(x, act, color=pc.color("active"), label="active")
        ax.bar(x, ina, bottom=act, color=pc.color("inactive"), label="inactive")
        _xticks(ax, [f"{CAT_DISPLAY[c]}\n({CUT[c]})" for c in CATS])
        ax.legend(loc="upper right")   # active/inactive key inside the panel
        compact_yaxis(ax)
        self._label(xlabel="", ylabel="molecules (unique)", title="Active vs inactive molecules")


# --- per-pathogen processing ------------------------------------------------

def process(pathogen_code):
    pathogen = load_pathogen(pathogen_code)
    out = os.path.join(STAGE4PATH, pathogen_code)
    print(f"Step 22 (plots) - {pathogen_code} ({pathogen})")
    df = pd.read_csv(os.path.join(out, "22_binarisation_summary.csv"))

    stats = binarisation_stats(pathogen_code)

    dr, sp = stats["dr_ratios"], stats["sp_ratios"]
    comp_ds = composition_counts(df)     # by dataset count
    comp_mol = composition_mols(df)      # by molecule count
    # unique molecules tested per category, split active (in >=1 assay) vs inactive
    act_mol = {c: {"active": len(stats["active_iks"][c]),
                   "inactive": len(stats["tested_iks"][c] - stats["active_iks"][c])}
               for c in CATS}

    # overview: filled 2x3, every legend inside its own panel (no shared/external key).
    # Row 1: per-category active-ratio histograms (A, B) + active ratio vs dataset size
    # (C). Row 2: class composition by datasets (D), by molecules (E), and unique
    # active/inactive molecules per category (F).
    render_figure(out, "22_binarisation_overview", 2, 3, [
        lambda ax: RatioHistOne(ax, dr, "DR", "Active ratio — Dose Response", abc="A"),
        lambda ax: RatioHistOne(ax, sp, "SP", "Active ratio — Single Point", abc="B"),
        lambda ax: RatioVsSizePlot(ax, df, abc="C"),
        lambda ax: CompositionPlot(ax, comp_ds, "datasets", "Class composition (datasets)", abc="D"),
        lambda ax: CompositionPlot(ax, comp_mol, "molecules", "Class composition (molecules)", abc="E"),
        lambda ax: ActiveMoleculesPlot(ax, act_mol, abc="F"),
    ])

    print(f"  contradictory -> DR {stats['dr_contra']:.3f} SP {stats['sp_contra']:.3f} "
          f"| cross shared {stats['cross_shared']:.3f} contradictory {stats['cross_contra']:.3f}")

    # eos1klk embedding figures (per pathogen): actives, category, and largest pools,
    # each across the four projections. Skipped if the coordinates are unavailable.
    # Each also emits a companion CSV of the plotted points (coords + group), which live
    # in no other stage-4 CSV.
    cc = pathogen_coords(out)
    if cc is not None and len(cc):
        active = stats["active_iks"]["DR"] | stats["active_iks"]["SP"]
        dr_t, sp_t = stats["tested_iks"]["DR"], stats["tested_iks"]["SP"]
        # 1. active / inactive, binarised (labelled) molecules only
        ccb = cc[cc["inchikey"].isin(dr_t | sp_t)]
        amask = ccb["inchikey"].isin(active).values
        EmbeddingFigure(ccb, [("inactive", ~amask, pc.color("inactive")),
                              ("active", amask, pc.color("active"))]).save(
            out, "22_eos1klk_actives", companion=True)
        # 2. DR-only / shared / SP-only / other (not binarised)
        isdr = cc["inchikey"].isin(dr_t).values
        issp = cc["inchikey"].isin(sp_t).values
        EmbeddingFigure(cc, [("other", ~(isdr | issp), pc.color("other")),
                             ("DR only", isdr & ~issp, pc.color("dr_only")),
                             ("SP only", issp & ~isdr, pc.color("sp_only")),
                             ("shared DR+SP", isdr & issp, pc.color("shared"))]).save(
            out, "22_eos1klk_category", companion=True)
        # 3. the 10 largest compatibility-group pools (across DR+SP), rest grey
        allp = sorted(stats["top_pools"]["DR"] + stats["top_pools"]["SP"],
                      key=lambda t: t[0], reverse=True)[:10]
        pcols = pc.cyclic(len(allp)) if allp else []
        iks = cc["inchikey"].values
        assigned = np.zeros(len(cc), dtype=bool)
        pool_groups = []
        for i, (sz, members) in enumerate(allp):
            m = np.fromiter((k in members for k in iks), dtype=bool, count=len(iks)) & ~assigned
            assigned |= m
            pool_groups.append((f"pool {i + 1} (n={sz})", m, pcols[i]))
        EmbeddingFigure(cc, [("other", ~assigned, pc.color("other"))] + pool_groups).save(
            out, "22_eos1klk_clusters", legend_ncol=4, companion=True)

    return stats


if __name__ == "__main__":
    pathogen_code = sys.argv[1]
    # Publication style (also applied at plot_base import; re-assert here so the import
    # order relative to plot_utils, which sets slide+ersilia, can never win).
    stylia.set_format("print")
    stylia.set_style("article")
    binarise(pathogen_code)
    process(pathogen_code)
