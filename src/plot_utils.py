"""
Shared figures subsystem for the stage-4 plotting scripts (steps 21-27).

Holds the pieces that both the per-pathogen scripts (21, 22) and the cross-pathogen
script (27) rely on, so they never diverge:

  - stylia house-style setup (print + article, applied at import).
  - Plot constants (target-type labels, attrition stages, size bins, category maps,
    eos1klk embedding coordinates, active-call cut-off notes, ...).
  - Generic drawing helpers (mol_region, cat_label, size_bins, _xticks).
  - eos1klk (lazychemvis) embedding helpers (load_coords, pathogen_coords,
    embed_figure, pathogen_colors).
  - Per-pathogen AGGREGATORS that read a pathogen's outputs and RETURN the results
    dict (no figures drawn):
        curation_stats(code)       - what step 21's process() computes.
        binarisation_stats(code)   - what step 22's process() computes.
    plus the cross-assay metric helpers they use (crossassay_metrics,
    composition_counts).

The per-pathogen figure helpers stay in scripts 21/22; the cross-pathogen figure
helpers live in script 27. This module holds only the shared inputs to both.
"""

import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import h5py
from matplotlib.patches import Patch

import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, root)
from default import RANDOM_SEED, STAGE1PATH, STAGE2PATH, STAGE4PATH
from dataset_utils import UNIT_SLUG, dataset_inchikey_sets

# Format: print | Style: article - publication house style, shared by all stage-4 plots
# (steps 21-27). Change with stylia.set_format() / stylia.set_style().
stylia.set_format("print")
stylia.set_style("article")

# --- shared plot constants: curation / attrition ----------------------------
CAP_POINTS = 1000  # max individual points drawn per pathogen in box overlays

# Full-ChEMBL bioactive chemical space (Panel A of the story figure): unique
# compounds with any bioactivity record, the universe step 07 draws pathogens from.
CHEMBL_ACTIVITIES_05 = os.path.join(STAGE1PATH, "05_activities_preprocessed.csv")

# Assay target types (build target_type_curated_extra) and their display labels.
TARGET_TYPES = ["ORGANISM", "SINGLE PROTEIN", "DISCARDED"]
TARGET_TYPE_LABELS = {"ORGANISM": "Whole-cell", "SINGLE PROTEIN": "Protein target",
                      "DISCARDED": "Other"}

# Mutually-exclusive regions for the molecule panel (stacked column summing to 1).
# A molecule's "Other"(DISCARDED) membership is folded into Whole-cell / Protein
# target whenever it also appears there; only molecules seen *exclusively* in
# DISCARDED assays count as "Other".
TT_MOL_REGIONS = ["Whole-cell", "Whole-cell + Protein target", "Protein target", "Other"]


def mol_region(types):
    org = "ORGANISM" in types
    tgt = "SINGLE PROTEIN" in types
    if org and tgt:
        return "Whole-cell + Protein target"
    if org:
        return "Whole-cell"
    if tgt:
        return "Protein target"
    return "Other"

# Ordered attrition stages: (label, discard_step threshold). A pair survives a
# stage iff discard_step == 0 or discard_step > threshold.
# Full curation cascade, one bar per gating step. This is the WATERFALL's stage list;
# the overview figure shows a reduced subset (drops the two near-flat stages Single-unit
# and DR/SP — see OVERVIEW_DROP). "DR/SP" (thr 5) folds the categorisation step (4) and the
# always-empty direction guard (5). Co-ADD and PubChem are separate bars.
STAGES = [
    ("All", -1),
    ("Whole-cell", 1),
    ("Quantitative", 2),
    ("Single-unit", 3),
    ("DR/SP", 5),
    ("Co-ADD", 6),
    ("PubChem", 7),
    ("≥6 mols", 8),
]

# Stages omitted from the OVERVIEW attrition bars (kept in the waterfall). Both are
# near-flat across every pathogen: Single-unit (multi-unit discard <=1 dataset) and
# DR/SP (uncategorised removes few molecules; the direction guard removes none).
OVERVIEW_DROP = {"Single-unit", "DR/SP"}

# Coarse discard categories, keyed by discard_step (see decide() in step 21).
# The two safeguard reasons (multi-unit, no-direction) are pushed last in STEP_ORDER.
STEP_LABELS = {
    0: "retained",
    1: "not Whole-cell",
    2: "qualitative-only",
    3: "multi-unit",
    4: "uncategorized",
    5: "no direction",
    6: "Co-ADD",
    7: "PubChem",
    8: "≤5 molecules",
}
STEP_ORDER = [0, 1, 2, 4, 6, 7, 8, 3, 5]  # headline reasons first, safeguards (3,5) last

CATEGORY_ORDER = ["DR", "SP", "UNCATEGORIZED"]
CATEGORY_DISPLAY = {"DR": "Dose Response", "SP": "Single Point",
                    "UNCATEGORIZED": "Uncategorized"}


def cat_label(c):
    return CATEGORY_DISPLAY.get(c, c)


# Dataset-size bins (molecules per dataset) for the Whole-cell size figure.
SIZE_EDGES = [5, 20, 50, 100, 1000]
SIZE_LABELS = ["≤5", "5–20", "20–50", "50–100", "100–1000", ">1000"]


def size_bins(values):
    """Bin index 0..5 for each value (np.searchsorted, '≤edge' semantics)."""
    return np.searchsorted(SIZE_EDGES, np.asarray(values, dtype=float), side="left")


# --- shared plot constants: binarisation / embeddings -----------------------
CATS = ["DR", "SP"]
CAT_DISPLAY = {"DR": "Dose Response", "SP": "Single Point"}
COMPOSITION = ["all-inactive", "mixed", "all-active"]
HIT_LABELS = ["1 assay", "2 assays", "≥3 assays"]  # times an active molecule is hit
COMPAT_LABELS = ["singleton", "in a pair", "group ≥3"]  # poolable-group size of a dataset

# eos1klk (lazychemvis) precomputed 2D coordinates: 4 projections x (x, y).
EOS_H5 = os.path.join(STAGE2PATH, "general", "12_eos1klk_predictions.h5")
EMBEDDINGS = [("PCA", "pca_x", "pca_y"), ("t-MAP", "tmap_x", "tmap_y"),
              ("t-SNE", "tsne_x", "tsne_y"), ("UMAP", "umap_x", "umap_y")]
_COORDS = None  # cached smiles -> 8-coord DataFrame

# Active-call cut-offs, annotated on every binarised plot.
CUT = {"DR": "≤10 µM", "SP": "≥50% / 20 mm"}


# --- generic drawing helpers ------------------------------------------------

def _xticks(ax, labels):
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")


# --- eos1klk embedding helpers (precomputed 2D coordinates) -----------------

def load_coords():
    """Cached smiles -> 8-coord DataFrame from the consolidated eos1klk HDF5. Returns
    None if the file is absent (embedding plots are then skipped)."""
    global _COORDS
    if _COORDS is not None:
        return _COORDS if len(_COORDS) else None
    if not os.path.exists(EOS_H5):
        print(f"  [eos] {EOS_H5} not found - skipping embedding plots")
        _COORDS = pd.DataFrame()
        return None
    with h5py.File(EOS_H5, "r") as h:
        feats = [f.decode() if isinstance(f, bytes) else f for f in h["Features"][:]]
        smiles = [s.decode() if isinstance(s, bytes) else s for s in h["Inputs"][:]]
        vals = h["Values"][:]
    _COORDS = pd.DataFrame(vals, columns=feats, index=pd.Index(smiles, name="smiles"))
    return _COORDS


def pathogen_coords(out):
    """Per-pathogen molecules (inchikey) joined to their eos1klk coordinates by SMILES."""
    coords = load_coords()
    if coords is None:
        return None
    cc = pd.read_csv(os.path.join(STAGE2PATH, os.path.basename(out), "07_compound_counts.csv.gz"),
                     usecols=["inchikey", "smiles"]).dropna(subset=["smiles"]).drop_duplicates("inchikey")
    cc = cc.join(coords, on="smiles").dropna(subset=[EMBEDDINGS[0][1]])
    return cc


def embed_figure(cc, groups, path, legend_ncol=None):
    """2x2 scatter of the four projections; `groups` = list of (label, bool-mask, colour)
    drawn in order (first = background). Rasterized for the dense point clouds."""
    fig, axs = stylia.create_figure(2, 2, width=0.8, height=0.8)
    for name, xcol, ycol in EMBEDDINGS:
        ax = axs.next()
        for label, mask, color in groups:
            if mask.any():
                ax.scatter(cc[xcol].values[mask], cc[ycol].values[mask], s=1, alpha=0.35,
                           color=color, edgecolors="none", rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        stylia.label(ax, xlabel="", ylabel="", title=name)
    # legend anchored by its top just below the panels (grows downward, never over data)
    fig.legend(handles=[Patch(color=c, label=l) for l, m, c in groups],
               loc="upper center", ncol=legend_ncol or len(groups), bbox_to_anchor=(0.5, -0.02),
               frameon=False)
    stylia.save_figure(path)


def pathogen_colors(codes):
    """Stable per-pathogen colours. Keyed to the canonical config/pathogens.csv order (not
    the caller's list order/length), so a given pathogen keeps the SAME colour across every
    figure and every run — regardless of how the codes were sorted (e.g. by dataset count)."""
    canonical = pd.read_csv(os.path.join(root, "..", "config", "pathogens.csv"))["code"].tolist()
    order = {c: i for i, c in enumerate(canonical)}
    cm = stylia.CyclicColormap("npg")  # article-family cyclic palette (was "ersilia")
    cm.fit(list(range(len(canonical) + 1)))
    palette = cm.transform(list(range(len(canonical))))
    return {c: palette[order[c]] if c in order else palette[i % len(palette)]
            for i, c in enumerate(codes)}


# --- cross-assay metric helpers ---------------------------------------------

def crossassay_metrics(pathogen_code, df):
    """Cross-assay agreement metrics over binarised molecules.

    Reads each binarised dataset file (one row per molecule, bin 0/1) and tracks,
    per InChIKey, the assays it is active/tested in within each category. Returns,
    for each category (DR, SP):

      contradictory = #{active in some assay AND inactive in another}
                      / #{tested in >= 2 assays}
                      (disagreement among molecules tested more than once)

    plus hit_bins, overlap and pool metrics, and a "cross" block comparing
    DR vs SP over molecules tested in BOTH categories
    (active-in-a-category = active in any assay of that category):

      shared        = #{active in both DR and SP} / #{tested in both}
      contradictory = #{active in one format, inactive in the other} / #{tested in both}
      n_both        = #{molecules tested in both categories}
    """
    res = {}
    rng = np.random.default_rng(RANDOM_SEED)
    active_any = {"DR": set(), "SP": set()}
    tested_any = {"DR": set(), "SP": set()}
    for cat in CATS:
        sub = df[df["category"] == cat]
        active_assays, tested_assays = Counter(), Counter()
        has_active, has_inactive = set(), set()
        tested_sets = []  # per-dataset tested-molecule sets (for pairwise overlap)
        group_mols = defaultdict(set)       # compatibility group -> pooled molecule union
        group_active = defaultdict(set)     # compatibility group -> active (any-assay) union
        group_n = defaultdict(int)          # compatibility group -> number of datasets
        for r in sub.itertuples(index=False):
            slug = UNIT_SLUG.get(r.unit, str(r.unit))
            f = os.path.join(STAGE4PATH, pathogen_code, "22_binarised", cat,
                             f"{r.assay_id}_{str(r.activity_type).replace('/', '-')}_{slug}.csv.gz")
            d = pd.read_csv(f, usecols=["inchikey", "bin"])
            actives = d.loc[d["bin"] == 1, "inchikey"].unique()
            inactives = d.loc[d["bin"] == 0, "inchikey"].unique()
            tested_here = set(d["inchikey"].unique())
            tested_sets.append(tested_here)
            sk = "" if pd.isna(getattr(r, "assay_strain_norm", "")) else getattr(r, "assay_strain_norm", "")
            gkey = (r.activity_type, r.unit, r.direction, sk)
            group_mols[gkey] |= tested_here
            group_active[gkey].update(actives)
            group_n[gkey] += 1
            for ik in actives:
                active_assays[ik] += 1
                has_active.add(ik)
            has_inactive.update(inactives)
            for ik in tested_here:
                tested_assays[ik] += 1
        multi = sum(1 for v in tested_assays.values() if v >= 2)
        contradictory = len(has_active & has_inactive)
        hit_bins = [sum(1 for v in active_assays.values() if v == 1),
                    sum(1 for v in active_assays.values() if v == 2),
                    sum(1 for v in active_assays.values() if v >= 3)]
        # pairwise tested-compound overlap (sampled) -> pooling regime signal
        jac, nset = [], len(tested_sets)
        if nset >= 2:
            for _ in range(min(3000, nset * (nset - 1) // 2)):
                i, j = int(rng.integers(0, nset)), int(rng.integers(0, nset))
                if i == j:
                    continue
                a, b = tested_sets[i], tested_sets[j]
                u = len(a | b)
                jac.append((len(a & b) / u) if u else 0.0)
        pool_keys = [k for k in group_mols if group_n[k] >= 2]
        pool_sizes = [len(group_mols[k]) for k in pool_keys]
        pool_ratios = [len(group_active[k]) / len(group_mols[k]) for k in pool_keys]
        res[cat] = {
            "contradictory": (contradictory / multi) if multi else np.nan,
            "active_iks": has_active,
            "hit_bins": hit_bins,
            "overlap_sample": jac,
            "frac_sharing": float(np.mean([1.0 if v > 0 else 0.0 for v in jac])) if jac else np.nan,
            "pool_sizes": pool_sizes,
            "pool_ratios": pool_ratios,
            "tested_iks": set(tested_assays),
        }
        active_any[cat] = has_active
        tested_any[cat] = set(tested_assays)

    both = tested_any["DR"] & tested_any["SP"]
    n_both = len(both)
    both_active = active_any["DR"] & active_any["SP"]
    one_only = ((active_any["DR"] - active_any["SP"]) | (active_any["SP"] - active_any["DR"])) & both
    res["cross"] = {
        "n_both": n_both,
        "shared": (len(both_active & both) / n_both) if n_both else np.nan,
        "contradictory": (len(one_only) / n_both) if n_both else np.nan,
    }
    return res


def composition_counts(df):
    out = {}
    for c in CATS:
        sub = df[df["category"] == c]
        out[c] = [int((sub["active_ratio"] == 0).sum()),
                  int(((sub["active_ratio"] > 0) & (sub["active_ratio"] < 1)).sum()),
                  int((sub["active_ratio"] == 1).sum())]
    return out


# --- per-pathogen aggregators (compute results dict; no figures) ------------

def curation_stats(pathogen_code):
    """Reproduce what step 21's process() computes and returns as its results dict.

    Reads 08_chembl_cleaned_data.csv.gz, 07_compound_counts.csv.gz and
    21_curation_summary.csv. Returns the same dict 21.process() returns, without
    drawing any figures.
    """
    acts = pd.read_csv(os.path.join(STAGE2PATH, pathogen_code, "08_chembl_cleaned_data.csv.gz"), low_memory=False)
    counts = pd.read_csv(os.path.join(STAGE2PATH, pathogen_code, "07_compound_counts.csv.gz"),
                         usecols=["compound_chembl_id", "inchikey"])
    summ = pd.read_csv(os.path.join(STAGE4PATH, pathogen_code, "21_curation_summary.csv"), low_memory=False)

    cid2ik = dict(zip(counts["compound_chembl_id"], counts["inchikey"]))
    acts["inchikey"] = acts["compound_chembl_id"].map(cid2ik)

    key = ["assay_chembl_id", "activity_type"]
    # canonical per-pair dataset membership (numeric ∪ recovered inactives) — the
    # numerator for every molecule count/union so figures match the written datasets.
    ik_dataset = dataset_inchikey_sets(acts, key)
    # full available chemistry (all tested molecules) — kept as the DENOMINATOR only.
    ik_all = acts.groupby(key)["inchikey"].apply(lambda s: set(s.dropna()))

    summ_idx = summ.set_index(["assay_id", "activity_type"])
    total_cs = len(set().union(*ik_all.values)) if len(ik_all) else 0

    # attrition per stage (counts over the dataset membership of survivors)
    stage_labels, stage_assays, stage_cs, stage_mols, stage_strains = [], [], [], [], []
    for label, thr in STAGES:
        survivors = summ_idx[(summ_idx["discard_step"] == 0) | (summ_idx["discard_step"] > thr)]
        keys = list(survivors.index)
        cs, mols = set(), 0
        for k in keys:
            s = ik_dataset.get(k, set())
            cs |= s
            mols += len(s)
        strains = survivors["assay_strain_norm"].dropna()
        strains = strains[strains.astype(str).str.strip() != ""]
        stage_labels.append(label)
        stage_assays.append(len(keys))
        stage_cs.append(len(cs))
        stage_mols.append(mols)
        stage_strains.append(strains.nunique())

    # retained characterization
    retained = summ[~summ["discarded"]].copy()
    ret_keys = list(zip(retained["assay_id"], retained["activity_type"]))
    # written dataset size = measured molecules + recovered qualitative inactives
    sizes = retained["n_mol_dataset"].tolist()
    # molecule reuse over dataset membership (recovered inactives included)
    mol_dataset_counter = Counter()
    for k in ret_keys:
        for ik in ik_dataset.get(k, set()):
            mol_dataset_counter[ik] += 1
    reuse = list(mol_dataset_counter.values())
    category_counts = dict(retained["category"].value_counts())

    # retained chemical space: union of the retained dataset molecule sets (feeds the
    # cross-pathogen overlap heatmap).
    retained_iks = set(mol_dataset_counter)

    # per-category (DR/SP) retained breakdown for the overview H panel: datasets,
    # molecules (sum of dataset sizes = same definition as the "Total molecules" bar),
    # and chemical space (unique InChIKeys per category; DR/SP may overlap via shared
    # molecules). Keyed by category so the panel can show the DR-vs-SP split per measure.
    cat_of = dict(zip(ret_keys, retained["category"]))
    category_datasets, category_mols = Counter(), Counter()
    cat_ik = defaultdict(set)
    for k in ret_keys:
        c = cat_of[k]
        s = ik_dataset.get(k, set())
        category_datasets[c] += 1
        category_mols[c] += len(s)
        cat_ik[c] |= s
    category_chem = {c: len(v) for c, v in cat_ik.items()}

    # cleaned (08) chemical space and its partition by assay target type, for the
    # cross-pathogen story figure. Molecule->assay from 08; assay->target type from
    # the build summary (target_type_curated_extra, constant per assay).
    all_iks = set(acts["inchikey"].dropna())
    assay_tt = summ.drop_duplicates("assay_id").set_index("assay_id")["target_type_curated_extra"]
    acts_tt = acts["assay_chembl_id"].map(assay_tt)
    # per molecule, the set of target types it was tested under -> collapsed region
    _mt = pd.DataFrame({"ik": acts["inchikey"], "tt": acts_tt}).dropna()
    _regions = _mt.groupby("ik")["tt"].agg(set).map(mol_region).value_counts().to_dict()
    tt_combo_counts = {k: int(_regions.get(k, 0)) for k in TT_MOL_REGIONS}
    adf = summ.drop_duplicates("assay_id")
    tt_assay_counts = {T: int((adf["target_type_curated_extra"] == T).sum())
                       for T in TARGET_TYPES}

    # Whole-cell (ORGANISM) dataset sizes BEFORE selection: every quantitative pair,
    # including those later discarded (<=5 molecules, Co-ADD, PubChem, multi-unit). Size
    # = dataset membership (numeric + recovered inactives), matching the written datasets.
    wc = summ[(summ["target_type_curated_extra"] == "ORGANISM") & (summ["n_numeric"] > 0)]
    wholecell_sizes = [s for s in wc["n_mol_dataset"].tolist() if s and not pd.isna(s)]

    # molecule-weighted discard breakdown (dataset membership per pair -> its bucket)
    step_mols = summ.groupby("discard_step")["n_mol_dataset"].sum().to_dict()

    # stats needed for the pct_* fields (mirrors the stats CSV in 21.process)
    org_cs = stage_cs[1] if len(stage_cs) > 1 else 0
    pct_total_cs = round(100 * stage_cs[-1] / total_cs, 2) if total_cs else 0.0
    pct_org_cs = round(100 * stage_cs[-1] / org_cs, 2) if org_cs else 0.0

    return {
        "code": pathogen_code,
        "n_retained": len(retained),
        "retained_iks": retained_iks,
        "all_iks": all_iks,
        "tt_combo_counts": tt_combo_counts,
        "tt_assay_counts": tt_assay_counts,
        "wholecell_sizes": wholecell_sizes,
        "pct_total_cs": float(pct_total_cs) if total_cs else 0.0,
        "pct_org_cs": float(pct_org_cs) if org_cs else 0.0,
        "sizes": [s for s in sizes if s and not pd.isna(s)],
        "reuse": reuse,
        "category_counts": category_counts,
        # per-category (DR/SP) retained metrics for the overview H panel
        "category_datasets": dict(category_datasets),
        "category_mols": dict(category_mols),
        "category_chem": category_chem,
        "step_counts": {int(k): int(v) for k, v in summ["discard_step"].value_counts().items()},
        "step_mols": {int(k): float(v) for k, v in step_mols.items()},
        # per-stage attrition series (sole computation; step 21 process() consumes these
        # instead of re-running the identical loop).
        "stage_labels": stage_labels, "stage_assays": stage_assays,
        "stage_cs": stage_cs, "stage_mols": stage_mols, "stage_strains": stage_strains,
        "total_cs": total_cs, "org_cs": org_cs,
    }


def binarisation_stats(pathogen_code):
    """Reproduce what step 22's process() computes and returns as its results dict.

    Reads 22_binarisation_summary.csv and 21_curation_summary.csv, then reads the
    binarised per-dataset files via crossassay_metrics. Returns the same dict
    22.process() returns.
    """
    df = pd.read_csv(os.path.join(STAGE4PATH, pathogen_code, "22_binarisation_summary.csv"))

    # pooling potential: datasets poolable iff they share (activity_type, unit,
    # direction, strain). Strain comes from the step-21 summary (not carried in the
    # binarisation summary), so merge it on before the file pass that pools molecules.
    strain = pd.read_csv(os.path.join(STAGE4PATH, pathogen_code, "21_curation_summary.csv"),
                         usecols=["assay_id", "activity_type", "assay_strain_norm"])
    dfm = df.merge(strain, on=["assay_id", "activity_type"], how="left")
    cm = crossassay_metrics(pathogen_code, dfm)
    keys = ["activity_type", "unit", "direction", "assay_strain_norm"]
    n_usable, compat_bins, n_pools = {}, {}, {}
    for cat in CATS:
        sub = dfm[dfm["category"] == cat]
        n_usable[cat] = int((~sub["degenerate"]).sum())
        gs = sub.groupby(keys, dropna=False)["assay_id"].transform("size")
        compat_bins[cat] = [int((gs == 1).sum()), int((gs == 2).sum()), int((gs >= 3).sum())]
        grp = sub.groupby(keys, dropna=False).size()
        n_pools[cat] = int((grp >= 2).sum())   # compatible groups that become a merged pool

    def cat_pairs(cat):
        sub = df[(df["category"] == cat) & df["active_ratio"].notna()]
        return sub["n_molecules"].tolist(), sub["active_ratio"].tolist()
    dr_sizes, dr = cat_pairs("DR")
    sp_sizes, sp = cat_pairs("SP")

    return {
        "code": pathogen_code, "n": len(df),
        "dr_ratios": dr, "sp_ratios": sp,
        "dr_sizes": dr_sizes, "sp_sizes": sp_sizes,
        "hit_bins": {c: cm[c]["hit_bins"] for c in CATS},
        "dr_contra": cm["DR"]["contradictory"], "sp_contra": cm["SP"]["contradictory"],
        "cross_shared": cm["cross"]["shared"], "cross_contra": cm["cross"]["contradictory"],
        "active_iks": {c: cm[c]["active_iks"] for c in CATS},
        "composition": composition_counts(df),
        "overlap": {c: cm[c]["overlap_sample"] for c in CATS},
        "frac_sharing": {c: cm[c]["frac_sharing"] for c in CATS},
        "n_usable": n_usable,
        "compat_bins": compat_bins,
        "n_pools": n_pools,
        "pool_sizes": {c: cm[c]["pool_sizes"] for c in CATS},
        "pool_ratios": {c: cm[c]["pool_ratios"] for c in CATS},
        # per-molecule tested/active sets (needed by the per-pathogen eos figures).
        "tested_iks": {c: cm[c]["tested_iks"] for c in CATS},
    }
