"""
Step 21 — Diagnosis plots.

Produces a 3×3 diagnostic figure for a given pathogen and saves it to
output/<pathogen_code>/21_diagnosis.png.

Usage
-----
    python scripts/21_diagnosis.py <pathogen_code>
"""

from collections import Counter, defaultdict
from scipy.stats import gaussian_kde
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib as mpl
import pandas as pd
import numpy as np
import stylia
from stylia import SpectralColormap
import h5py
from tqdm import tqdm
import sys
import os
import matplotlib.pyplot as plt

root = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import load_pathogen, load_expert_cutoffs

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print(f"Step 21: Diagnosis plots for {pathogen_code} ({pathogen})")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)


def load_if_exists(path, **kwargs):
    if os.path.exists(path):
        return pd.read_csv(path, **kwargs)
    return pd.DataFrame()


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

print("Loading data (raw/cleaned activities, assay metadata, final datasets, assay master)...")

raw = pd.read_csv(os.path.join(OUTPUT, "07_chembl_raw_data.csv.gz"), low_memory=False)
cleaned = pd.read_csv(os.path.join(OUTPUT, "08_chembl_cleaned_data.csv.gz"), low_memory=False)
assays_raw_df = pd.read_csv(os.path.join(OUTPUT, "07_assays_raw.csv"))
assays_cleaned_df = pd.read_csv(os.path.join(OUTPUT, "08_assays_cleaned.csv"))
strain_info = pd.read_csv(os.path.join(OUTPUT, "09_assays_parameters_full.csv"))
compound_counts = pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))
final_datasets = pd.read_csv(os.path.join(OUTPUT, "17_final_datasets.csv"))
if "selected" in final_datasets.columns:
    final_datasets["selected"] = final_datasets["selected"].astype(bool)
correlations = pd.read_csv(os.path.join(OUTPUT, "17_dataset_correlations.csv"))

# Load assays master data for rejection analysis
assays_master = pd.read_csv(os.path.join(OUTPUT, "18_assays_master.csv"))

# Additional loads for summary table (step 23 merged)
datasets_12      = load_if_exists(os.path.join(OUTPUT, "12_datasets.csv"))
individual_lm    = load_if_exists(os.path.join(OUTPUT, "13_individual_LM.csv"))
indiv_selected   = load_if_exists(os.path.join(OUTPUT, "14_individual_selected_LM.csv"))
merged_lm        = load_if_exists(os.path.join(OUTPUT, "15_merged_LM.csv"))
merging_analysis = load_if_exists(os.path.join(OUTPUT, "15_merging_analysis.csv"))
merged_selected  = load_if_exists(os.path.join(OUTPUT, "16_merged_selected_LM.csv"))
general_model    = load_if_exists(os.path.join(OUTPUT, "20_general_model.csv"))
expert_cutoffs   = load_expert_cutoffs(CONFIGPATH)

n_selected = int(final_datasets["selected"].sum()) if len(final_datasets) > 0 else 0
any_selected = n_selected > 0
print(f"  Unique compounds in ChEMBL: {len(compound_counts):,}")
print(f"  Assays (raw / cleaned): {len(assays_raw_df):,} / {len(assays_cleaned_df):,}")
print(f"  Final datasets (total / selected): {len(final_datasets):,} / {n_selected:,}")
if any_selected:
    label_counts = final_datasets[final_datasets["selected"]]["label"].value_counts().to_dict()
    print(f"    Selected by label: A={label_counts.get('A', 0)}, B={label_counts.get('B', 0)}, M={label_counts.get('M', 0)}")
else:
    print("  No datasets selected — figure will use 2×2 layout")

pathogen_compounds = set(compound_counts["compound_chembl_id"])

# ---------------------------------------------------------------------------
# Summary counts
# ---------------------------------------------------------------------------

activities_raw = len(raw)
compounds_raw = len(set(raw["compound_chembl_id"]))
del raw

activities_cleaned = len(cleaned)
compounds_cleaned = len(set(cleaned["compound_chembl_id"]))
n_assays_raw = len(assays_raw_df)
n_assays_cleaned = len(assays_cleaned_df)

act_qt = len(cleaned[(cleaned["value"].notna()) & (cleaned["text_flag"] == 0)])
act_mx = len(cleaned[(cleaned["value"].notna()) & (cleaned["text_flag"].isin([-1, 1]))])

# Assay-level fractions (from assays_cleaned)
assay_type_ctr = Counter(assays_cleaned_df["assay_type"])
target_type_ctr = Counter(assays_cleaned_df["target_type"])
organism_ctr = Counter(assays_cleaned_df["assay_organism"])
strain_ctr = Counter(strain_info["assay_strain_curated"])
unit_ctr = Counter(assays_cleaned_df["unit"])

total_assay = sum(assay_type_ctr.values())
assay_type_bars = [
    assay_type_ctr["F"] / total_assay,
    (assay_type_ctr["F"] + assay_type_ctr["B"]) / total_assay,
    1,
][::-1]

total_target = sum(target_type_ctr.values())
target_type_bars = [
    target_type_ctr["ORGANISM"] / total_target,
    (target_type_ctr["ORGANISM"] + target_type_ctr["SINGLE PROTEIN"]) / total_target,
    1,
][::-1]

# Calculate strain bars: [no-strain, other strains, majoritarian strain]
total_strain = sum(strain_ctr.values())
null_strain_count = (
    strain_ctr.get(None, 0) +           # Handle None
    strain_ctr.get('', 0) +             # Handle empty strings
    strain_ctr.get(np.nan, 0)           # Handle pandas NaN (float)
)

if total_strain > 0:
    # Find majoritarian strain (excluding null values)
    non_null_strain_ctr = {k: v for k, v in strain_ctr.items() if k and k != '' and not pd.isna(k)}
    if non_null_strain_ctr:
        majoritarian_strain, majoritarian_count = max(non_null_strain_ctr.items(), key=lambda x: x[1])
        other_strains_count = sum(non_null_strain_ctr.values()) - majoritarian_count
    else:
        majoritarian_strain, majoritarian_count = None, 0
        other_strains_count = 0

    # Create 3-segment bars: [no-strain, other strains, majoritarian strain]
    no_strain_frac = null_strain_count / total_strain
    other_strains_frac = other_strains_count / total_strain
    majoritarian_frac = majoritarian_count / total_strain
    print(f"  Strain distribution: no-strain={no_strain_frac:.1%}, other={other_strains_frac:.1%}, majoritarian='{majoritarian_strain}' ({majoritarian_frac:.1%})")

    strain_bars = [
        majoritarian_frac,
        majoritarian_frac+other_strains_frac,
        1
    ][::-1]  # Reverse to match stacking pattern
else:
    strain_bars = [1, 1, 1]
    majoritarian_strain = "Unknown"

total_unit = sum(unit_ctr.values())
unit_bars = [
    unit_ctr["umol.L-1"] / total_unit,
    (unit_ctr["umol.L-1"] + unit_ctr["%"]) / total_unit,
    1,
][::-1]

# ---------------------------------------------------------------------------
# Compound occurrence rank plot
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Cumulative chemical space coverage per assay
# ---------------------------------------------------------------------------

print(f"Computing cumulative chemical space coverage across {len(assays_cleaned_df):,} assays...")

assay_ids_set = set(assays_cleaned_df["assay_id"])
assay_to_compounds = defaultdict(set)
for assay_id, activity_type, unit, cpd in cleaned[
    ["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"]
].values:
    if assay_id in assay_ids_set:
        assay_to_compounds[(assay_id, activity_type, unit)].add(cpd)

cum_cpds_set = set()
cum_prop = []
cpds_per_assay = []
for assay_id, activity_type, unit in tqdm(
    assays_cleaned_df[["assay_id", "activity_type", "unit"]].values
):
    key = (assay_id, activity_type, unit)
    cum_cpds_set |= assay_to_compounds[key]
    cum_prop.append(len(cum_cpds_set) / len(pathogen_compounds))
    cpds_per_assay.append(len(assay_to_compounds[key]))

print(f"  Coverage after all assays: {cum_prop[-1]:.1%} of the pathogen chemical space" if cum_prop else "  No assays to compute coverage.")

# Count activities per (assay, activity_type, unit) triplet — needed for activity coverage later
triplet_to_n_activities = (
    cleaned[cleaned["assay_chembl_id"].isin(assay_ids_set)]
    .groupby(["assay_chembl_id", "activity_type", "unit"], dropna=False)
    .size()
    .to_dict()
)
del cleaned

# ---------------------------------------------------------------------------
# Final dataset coverage per label
# ---------------------------------------------------------------------------

def _parse_assay_key(s):
    parts = s.split("|")
    return (parts[0], parts[1], np.nan if parts[2] == "" else parts[2])


final_coverage = {label: set() for label in "ABM"}
covered_triplets = {label: set() for label in "ABM"}
all_assay_triplets = set(assay_to_compounds.keys())
merged_dir = os.path.join(OUTPUT, "12_datasets", "M")
if any_selected:
    for label, name, assay_keys in final_datasets[final_datasets["selected"]][
        ["label", "name", "assay_keys"]
    ].values:
        for s in assay_keys.split(";"):
            key = _parse_assay_key(s)
            covered_triplets[label].add(key)
        if label == "M":
            filepath = os.path.join(merged_dir, f"{name}.csv.gz")
            if os.path.exists(filepath):
                df_m = pd.read_csv(filepath)
                final_coverage[label] |= set(df_m[df_m["smiles"] != "decoy"]["compound_chembl_id"])
        else:
            for s in assay_keys.split(";"):
                final_coverage[label] |= assay_to_compounds.get(_parse_assay_key(s), set())

all_covered_cpds = final_coverage["A"] | final_coverage["B"] | final_coverage["M"]
all_covered_triplets = covered_triplets["A"] | covered_triplets["B"] | covered_triplets["M"]
n_total_triplets = len(all_assay_triplets)

activities_in_selected = sum(triplet_to_n_activities.get(k, 0) for k in all_covered_triplets)
n_total_activities_cleaned = sum(triplet_to_n_activities.values())

if any_selected:
    print("Final dataset coverage (selected datasets only):")
    print(
        f"  Compound coverage (fraction of {len(pathogen_compounds):,} unique compounds): "
        f"A={len(final_coverage['A'])/len(pathogen_compounds):.1%}  "
        f"B={len(final_coverage['B'])/len(pathogen_compounds):.1%}  "
        f"M={len(final_coverage['M'])/len(pathogen_compounds):.1%}  "
        f"ALL={len(all_covered_cpds)/len(pathogen_compounds):.1%}"
    )
    print(
        f"  Assay coverage (fraction of {n_total_triplets:,} assays): "
        f"A={len(covered_triplets['A'])/n_total_triplets:.1%}  "
        f"B={len(covered_triplets['B'])/n_total_triplets:.1%}  "
        f"M={len(covered_triplets['M'])/n_total_triplets:.1%}  "
        f"ALL={len(all_covered_triplets)}/{n_total_triplets} ({len(all_covered_triplets)/n_total_triplets:.1%})"
    )
    print(
        f"  Bioactivity endpoint coverage: {activities_in_selected:,} endpoints in selected datasets  "
        f"/ {activities_cleaned:,} cleaned ({activities_in_selected/activities_cleaned:.1%})  "
        f"/ {activities_raw:,} raw ({activities_in_selected/activities_raw:.1%})"
    )

# ---------------------------------------------------------------------------
# Correlation matrices (clustered by compound overlap)
# ---------------------------------------------------------------------------

print(f"Building correlation matrices for {len(set(correlations['name_1']))} datasets (compound overlap + top-100 hit overlap)...")

names = sorted(set(correlations["name_1"]))
co_dict = {(r.name_1, r.name_2): r.compound_overlap for r in correlations.itertuples()}
ho_dict = {(r.name_1, r.name_2): r.hit_overlap_100 for r in correlations.itertuples()}

X_co = np.array([[co_dict[(n1, n2)] for n2 in names] for n1 in names])
X_ho = np.array([[ho_dict[(n1, n2)] for n2 in names] for n1 in names])

if len(names) > 1:
    Z = linkage(squareform(1 - X_co, checks=False), method="average")
    idx = leaves_list(Z)
    X_co = X_co[np.ix_(idx, idx)]
    X_ho = X_ho[np.ix_(idx, idx)]

# ---------------------------------------------------------------------------
# tSNE
# ---------------------------------------------------------------------------

from sklearn.decomposition import PCA
from openTSNE import TSNE as openTSNE

print("Loading ECFPs for tSNE (up to 10,000 pathogen + 30,000 background compounds)...")
h5_path = os.path.join(DATAPATH, "chembl_processed", "06_chembl_ecfps.h5")
with h5py.File(h5_path, "r") as f:
    all_ids_h5 = f["smiles"][:, 3].astype(str)
    all_fps_h5 = f["X_morgan"][:]
ecfps = {cid: fp for cid, fp in zip(all_ids_h5, all_fps_h5)}

reference_ids = compound_counts["compound_chembl_id"].tolist()[:10_000]
n_ref = len(reference_ids)
ids_set = set(reference_ids)
rng = np.random.default_rng(42)
pool = np.array([k for k in ecfps if k not in ids_set])
bg_ids = rng.choice(pool, size=n_ref * 3, replace=False)
all_ids = np.concatenate([reference_ids, bg_ids])

X_tsne = (np.stack([ecfps[i] for i in all_ids]) > 0).astype(np.float32)
print(f"  Input matrix: {X_tsne.shape[0]:,} compounds × {X_tsne.shape[1]} ECFP bits")
print("  PCA (reducing to 16 components before tSNE)...")
X_tsne = PCA(n_components=16, random_state=42, svd_solver="randomized").fit_transform(X_tsne)
print("  tSNE (2D embedding)...")
perp = min(30.0, max(5.0, (X_tsne.shape[0] - 1) / 3.0))
emb = np.asarray(
    openTSNE(
        n_components=2, perplexity=perp, random_state=42,
        n_jobs=16, negative_gradient_method="fft", n_iter=100,
    ).fit(X_tsne)
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def scatter_density(ax, xy, size_min=10, size_max=120, cmap="viridis", alpha=0.9, zorder=3):
    xy = np.asarray(xy)
    x, y = xy[:, 0], xy[:, 1]
    kde = gaussian_kde(np.vstack([x, y]))
    dens = kde(np.vstack([x, y]))
    d0, d1 = dens.min(), dens.max()
    dn = (dens - d0) / (d1 - d0 + 1e-12)
    sizes = size_min + dn * (size_max - size_min)
    return ax.scatter(x, y, c=dens, s=sizes, cmap=cmap, alpha=alpha, edgecolors="none", zorder=zorder)


def left_to_right(y_left):
    return 10 ** (5 * np.asarray(y_left))


def right_to_left(y_right):
    return np.log10(np.maximum(np.asarray(y_right), 1e-100)) / 5


def parse_rejection_categories(comment_series):
    """Parse rejection comments into mutually-exclusive categories.

    Each comment is assigned to exactly one category (the first match in
    priority order). An 'other' catch-all ensures every row is accounted for
    so that the stacked bars always sum to 1.
    """
    named = {
        'selected':              comment_series.str.contains('Retained in final selection', na=False),
        'already_accepted':      comment_series.str.contains('already accepted', na=False),
        'non_organism':          comment_series.str.contains('non-ORGANISM target type', na=False),
        'qualitative_only':      comment_series.str.contains('only qualitative data', na=False),
        'no_activity_data':      comment_series.str.contains('no activity data', na=False),
        'no_cutoff':             comment_series.str.contains('no expert cutoff defined', na=False),
        'too_few_compounds':     comment_series.str.contains('insufficient compounds', na=False),
        'too_few_positives':     comment_series.str.contains('insufficient positives|insufficient actives', na=False, regex=True),
        'ratio_out_of_range':    comment_series.str.contains('active ratio', na=False),
        'middle_cutoff_failure': comment_series.str.contains('middle cutoff', na=False),
        'insufficient_compatible': comment_series.str.contains('insufficient compatible assays', na=False),
        'fractional_contribution': comment_series.str.contains('insufficient_fractional_contribution', na=False),
        'auroc_below':           comment_series.str.contains('below 0.70 threshold', na=False),
        'correlation':           comment_series.str.contains('high correlation', na=False),
    }
    # Build mutually-exclusive masks: each row takes the first category that matches
    assigned = pd.Series(False, index=comment_series.index)
    categories = {}
    for name, mask in named.items():
        categories[name] = mask & ~assigned
        assigned |= categories[name]
    categories['other'] = ~assigned
    return categories


def calculate_rejection_proportions(df):
    """Calculate proportions for each label. Bars are guaranteed to sum to 1."""
    results = {}
    total = len(df)
    if total == 0:
        empty = {cat: 0 for cat in parse_rejection_categories(pd.Series([], dtype=str))}
        return {"A": dict(empty), "B": dict(empty), "M": dict(empty)}

    for label, col in [('A', 'comment_A'), ('B', 'comment_B'), ('M', 'comment_M')]:
        cats = parse_rejection_categories(df[col])
        results[label] = {cat: cats[cat].sum() / total for cat in cats}

    return results


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

layout = "3×3" if any_selected else "2×2"
print(f"Plotting {layout} diagnostic figure ({layout} panels)...")
stylia.set_style("ersilia")
nc = stylia.NamedColors()  

if any_selected:
    fig, axs = stylia.create_figure(3, 3, width=1, height=1)
else:
    fig, axs = stylia.create_figure(2, 2, width=1, height=1)
cmap2 = mpl.colors.LinearSegmentedColormap.from_list("purple_blue", [nc.plum, nc.blue], N=256)
label_to_color = {"A":nc.orange, "B": nc.blue, "M": nc.yellow}

# [0][0] Raw vs cleaned: activities, compounds, assays
ax = axs.next()
ax.bar([0], [activities_raw], color=nc.gray, label="Raw")
ax.bar([0], [activities_cleaned], color=nc.blue, label="Cleaned")
ax.bar([1], [compounds_raw], color=nc.gray)
ax.bar([1], [compounds_cleaned], color=nc.blue)
ax.bar([2], [n_assays_raw], color=nc.gray)
ax.bar([2], [n_assays_cleaned], color=nc.blue)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(["Activities", "Compounds", "Assays"])
ax.legend(loc="upper right")
stylia.label(ax, ylabel="Number", xlabel="", title="Raw vs cleaned data")

# [0][1] Fraction bars: assay type, target type, strain, unit, activity type
ax = axs.next()
ax.bar([0, 0, 0], assay_type_bars, zorder=2, color=[nc.gray, nc.pink,nc.yellow])
ax.text(0, assay_type_bars[-1] / 2, "F", ha="center", va="center")
ax.text(0, (assay_type_bars[-2] + assay_type_bars[-1]) / 2, "B", ha="center", va="center")
ax.bar([1, 1, 1], target_type_bars, zorder=2, color=[nc.gray, nc.pink, nc.yellow])
ax.text(1, target_type_bars[-1] / 2, "ORG", ha="center", va="center")
ax.text(1, (target_type_bars[-2] + target_type_bars[-1]) / 2, "SP", ha="center", va="center")
ax.bar([2, 2, 2], strain_bars, zorder=2, color=[nc.gray, nc.pink,nc.yellow])
# Add text labels for strain information
if 'majoritarian_strain' in locals() and majoritarian_strain:
    ax.text(2, strain_bars[-1] / 2, majoritarian_strain[:6], ha="center", va="center")
# Add label for "other strains" segment
if len(strain_bars) >= 2:
    other_center = (strain_bars[-2] + strain_bars[-1]) / 2
    ax.text(2, other_center, "other", ha="center", va="center")
ax.bar([3, 3, 3], unit_bars, zorder=2, color=[nc.gray, nc.pink,nc.yellow])
ax.text(3, unit_bars[-1] / 2, "uM", ha="center", va="center")
ax.text(3, (unit_bars[-2] + unit_bars[-1]) / 2, "%", ha="center", va="center")
ax.bar(
    [4, 4, 4],
    [1, (act_qt + act_mx) / activities_cleaned, act_qt / activities_cleaned],
    zorder=2, color=[nc.gray, nc.pink,nc.yellow],
)
ax.text(4, (act_qt / activities_cleaned) / 2, "qt.", ha="center", va="center")
ax.set_ylim([0, 1])
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_xticklabels(["Assay\ntype", "Target\ntype", "Strain", "Unit", "Activity\ntype"])
stylia.label(ax, ylabel="Fraction", xlabel="", title="Data classification by:")

# [0][2] tSNE
ax = axs.next()
cmap_tsne = mpl.colors.LinearSegmentedColormap.from_list("purple_yellow", ["#50285A", "#FAD782"], N=256)
ax.scatter(emb[n_ref:][:, 0], emb[n_ref:][:, 1], s=5, c="lightgray")
scatter_density(ax, emb[:n_ref, :2], size_min=0, size_max=20, cmap=cmap_tsne, alpha=1)
ax.set_yticks([])
ax.set_xticks([])
stylia.label(ax, xlabel="tSNE-2", ylabel="tSNE-1", title="tSNE (pathogen compounds colored by density)")

"""
# [1][0] Compound occurrence rank
ax = axs.next()
x_rank = list(range(1, len(count_cpds) + 1))
ax.plot(x_rank, count_cpds, c="#AA96FA")
ax.fill_between(x_rank, count_cpds, 1, color="#AA96FA", alpha=0.6)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([1, 10, 10**2, 10**3, 10**4, 10**5])
for xi in [1, 10, 100, 1000, 10000, 100000]:
    if xi <= len(count_cpds):
        ax.scatter([xi], [count_cpds[xi - 1]], zorder=4, ec="k", s=25, c="#D2D2D2")
stylia.label(ax, xlabel="Compound number", ylabel="Compound occurrences", title="Compound occurrences ranked")
"""

# [1][1] Compounds per assay (yellow) + cumulative coverage (purple), dual axis
ax = axs.next()
x_marks = [i for i in [1, 10, 100, 1000, 10000] if i <= len(cum_prop)]
x_cum = np.arange(1, len(cum_prop) + 1)
x_cpds_arr = np.arange(1, len(cpds_per_assay) + 1)
cum_as_compounds = left_to_right(np.array(cum_prop))
ax.set_xscale("log")
ax.set_xticks([1, 10, 100, 1000, 10000])
ax.set_yscale("log")
ax.set_yticks([1, 10, 100, 1_000, 10_000, 100_000])
ax.set_yticklabels([r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
ax.plot(x_cpds_arr, cpds_per_assay, c=nc.yellow, lw=1.8, zorder=3)
ax.fill_between(x_cpds_arr, cpds_per_assay, 1, color=nc.yellow, alpha=0.6, zorder=2)
for xi in x_marks:
    ax.scatter([xi], [cpds_per_assay[xi - 1]], zorder=4, ec="k", s=25, color=nc.gray)
if len(cpds_per_assay) > 0:
    ymax = max(max(cpds_per_assay), np.max(cum_as_compounds), left_to_right(1.0))
    ax.set_ylim(0.6, ymax * 2)
    secax = ax.secondary_yaxis("right", functions=(right_to_left, left_to_right))
    secax.set_ylabel("Fraction (cum) of the chemical space")
    secax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    secax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.1f"))
    ax.plot(x_cum, cum_as_compounds, c=nc.plum, lw=1.8, zorder=5)
    ax.fill_between(x_cum, cum_as_compounds, 1, color=nc.plum, alpha=0.25, zorder=1)
    for xi in x_marks:
        ax.scatter([xi], [cum_as_compounds[xi - 1]], zorder=6, ec="k", s=25, color=nc.gray)
else:
    ax.text(0.5, 0.5, "No data available", ha="center", va="center",
            transform=ax.transAxes, color="gray")
stylia.label(ax, xlabel="Number of assays", ylabel="Number of compounds", title="Chemical space coverage by assays")

if any_selected:
    # [1][2] Compound overlap heatmap (clustered)
    ax = axs.next()
    if len(names) > 0:
        im = ax.imshow(X_co, vmin=0, vmax=1, cmap=cmap2)
        fig.colorbar(im, ax=ax, fraction=0.045)
    else:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes, color="gray")
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, title="Compound overlap in datasets", ylabel="", xlabel="")

    # [2][0] Hit overlap heatmap (clustered)
    ax = axs.next()
    if len(names) > 0:
        im2 = ax.imshow(X_ho, vmin=0, vmax=1, cmap=cmap2)
        fig.colorbar(im2, ax=ax, fraction=0.045)
    else:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes, color="gray")
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, title="Top-100 hit overlap", xlabel="", ylabel="")

    # [2][1] Chemical space coverage by label
    ax = axs.next()
    all_coverage = final_coverage["A"] | final_coverage["B"] | final_coverage["M"]
    ax.set_xlim([0, 5])
    ax.set_ylim([0, 1])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(["A", "B", "M", "ALL"])
    for xi, lbl, color in [(1, "A", nc.orange), (2, "B", nc.blue), (3, "M", nc.yellow)]:
        frac = len(final_coverage[lbl]) / len(pathogen_compounds)
        ax.bar([xi, xi], [1, frac], zorder=2, color=[nc.gray, color], ec="k")
    frac_all = len(all_coverage) / len(pathogen_compounds)
    ax.bar([4, 4], [1, frac_all], zorder=2, color=[nc.gray, nc.mint], ec="k")
    stylia.label(ax, ylabel="Chemical space percentage", xlabel="", title="Chemical space coverage by label")

    # [2][2] Compounds vs positives per selected dataset
    ax = axs.next()
    selected_df = final_datasets[final_datasets["selected"]].reset_index(drop=True)
    ax.scatter(
        selected_df["cpds"], selected_df["positives"], zorder=2, s=50,
        c=[label_to_color[l] for l in selected_df["label"]], ec="k", alpha=0.7,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([10, 100, 1000, 10000, 100000])
    ax.set_yticks([10, 100, 1000, 10000, 100000])
    stylia.label(ax, xlabel="Number of compounds", ylabel="Number of positives",
                 title=f"Number of datasets: {len(selected_df)}")

    # [2][2] Dataset rejection reasons stacked barplot
    ax = axs.next()

    rejection_props = calculate_rejection_proportions(assays_master)

    # Stack categories (from bottom to top)
    stack_order = [
        'non_organism',
        'qualitative_only',
        'no_activity_data',
        'no_cutoff',
        'insufficient_compatible',
        'fractional_contribution',
        'too_few_compounds',
        'too_few_positives',
        'ratio_out_of_range',
        'middle_cutoff_failure',
        'auroc_below',
        'correlation',
        'already_accepted',
        'other',
        'selected',
    ]

    # Assign one color per category from SpectralColormap fitted to the number of categories
    ccm = SpectralColormap("npg")
    colors = dict(zip(stack_order, ccm.sample(len(stack_order))))

    # Create stacked bars
    bar_labels = ['A', 'B', 'M']
    bar_positions = [0, 1, 2]
    bottoms = [0, 0, 0]

    for category in stack_order:
        heights = []
        for bar_label in bar_labels:
            prop = rejection_props[bar_label].get(category, 0)
            heights.append(prop)

        ax.bar(bar_positions, heights, bottom=bottoms,
               color=colors[category], label=category.replace('_', ' ').title(),
               edgecolor='k', linewidth=0.3, alpha=0.8)

        # Update bottoms for next stack level
        bottoms = [b + h for b, h in zip(bottoms, heights)]

    ax.set_xticks(bar_positions)
    ax.set_xticklabels(bar_labels)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Proportion of assays')
    ax.set_xlabel('Condition label')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    stylia.label(ax, title='Assay rejection reasons by condition')

fig.suptitle(pathogen, size=8, y=1.01)


plt.tight_layout()

outpath = os.path.join(OUTPUT, "21_diagnosis.png")
stylia.save_figure(outpath)
print(f"Saved {layout} diagnostic figure -> {outpath}")

# ---------------------------------------------------------------------------
# Summary table (4-sheet Excel workbook)
# ---------------------------------------------------------------------------

print("\nBuilding summary table (21_summary_table.xlsx)...")

sel_df = final_datasets[final_datasets["selected"]] if len(final_datasets) > 0 else pd.DataFrame()

# ---- Sheet 1: Pipeline Funnel ----

print("  Sheet 1: Pipeline Funnel...")


def _row(stage, n_assays, n_compounds, notes=""):
    return {"Stage": stage, "N assays / datasets": n_assays,
            "N compounds": n_compounds, "Notes": notes}


funnel = []

funnel.append(_row(
    "Raw assays (step 07)", len(assays_raw_df), len(pathogen_compounds),
    "All assays from ChEMBL for this pathogen before any curation"
))
funnel.append(_row(
    "Cleaned assays (step 08)", len(assays_cleaned_df), "",
    "After organism, unit, direction and activity-type curation"
))

if len(datasets_12) > 0:
    eligible = datasets_12[datasets_12["dataset_type"].isin(["quantitative", "mixed"])]
    n_eligible = eligible[["assay_id", "activity_type", "unit"]].drop_duplicates().shape[0]
else:
    n_eligible = 0
funnel.append(_row(
    "Assays with quantitative/mixed datasets (step 12)", n_eligible, "",
    "Have numeric measurements binarizable at ≥1 expert cutoff"
))

if len(individual_lm) > 0:
    n_modelled_ab = individual_lm[["assay_id", "activity_type", "unit"]].drop_duplicates().shape[0]
else:
    n_modelled_ab = 0
funnel.append(_row(
    "Modelled individually A+B (step 13)", n_modelled_ab, "",
    "Met condition A or B criteria; RF model trained (4-fold CV)"
))

funnel.append(_row(
    "Selected individually A+B (step 14)", len(indiv_selected), "",
    "Best-cutoff AUROC > 0.7; one dataset per assay triplet"
))

if len(merging_analysis) > 0:
    n_merge_cands = merging_analysis[["assay_id", "activity_type", "unit"]].drop_duplicates().shape[0]
else:
    n_merge_cands = 0
funnel.append(_row(
    "Sent to merging (step 15)", n_merge_cands, "",
    "Not accepted individually; compatible assays pooled"
))

if len(merged_lm) > 0:
    base_names = merged_lm["name"].apply(lambda n: "_".join(n.split("_")[:2]))
    n_merged_groups = base_names.nunique()
    n_merged_combos = len(merged_lm)
else:
    n_merged_groups = 0
    n_merged_combos = 0
funnel.append(_row(
    "Merged groups modelled (step 15)", n_merged_groups, "",
    f"{n_merged_combos} group-cutoff combinations across {n_merged_groups} unique groups"
))

funnel.append(_row(
    "Merged groups selected (step 16)", len(merged_selected), "",
    "Best cutoff per group with AUROC > 0.7"
))

funnel.append(_row(
    "Before deduplication (step 17)", len(final_datasets), "",
    "All A+B+M candidates combined before correlation-based deduplication"
))

n_final = len(sel_df)
cpds_final = len(all_covered_cpds) if any_selected else 0
pct_coverage = round(100 * cpds_final / len(pathogen_compounds), 1) if len(pathogen_compounds) > 0 else ""
funnel.append(_row(
    "Final selected datasets (step 17)", n_final, cpds_final,
    f"After greedy deduplication (A→B→M priority); ~{pct_coverage}% of pathogen chemical space"
))

if len(general_model) > 0:
    funnel.append(_row(
        "General organism model (step 20)", len(general_model),
        int(general_model["n_compounds"].sum()),
        "One dataset per (activity_type, unit) pair; all ORGANISM assays pooled at middle cutoff"
    ))

sheet1 = pd.DataFrame(funnel)

# ---- Sheet 2: By Condition ----

print("  Sheet 2: By Condition...")


def _cond_stats(label, sub, cpds_col, pos_col, auroc_col):
    if sub is None or len(sub) == 0:
        return {"Condition": label, "N datasets": 0, "N compounds (total)": 0,
                "N actives (total)": 0, "N inactives (total)": 0,
                "Mean active ratio": np.nan, "AUROC mean": np.nan,
                "AUROC median": np.nan, "AUROC std": np.nan,
                "AUROC min": np.nan, "AUROC max": np.nan}
    n_cpds = int(sub[cpds_col].sum()) if cpds_col in sub.columns else 0
    n_pos  = int(sub[pos_col].sum())  if pos_col  in sub.columns else 0
    aurocs = sub[auroc_col].dropna()  if auroc_col in sub.columns else pd.Series(dtype=float)
    return {
        "Condition": label,
        "N datasets": len(sub),
        "N compounds (total)": n_cpds,
        "N actives (total)": n_pos,
        "N inactives (total)": n_cpds - n_pos,
        "Mean active ratio": round(n_pos / n_cpds, 4) if n_cpds > 0 else np.nan,
        "AUROC mean":   round(aurocs.mean(),   3) if len(aurocs) > 0 else np.nan,
        "AUROC median": round(aurocs.median(), 3) if len(aurocs) > 0 else np.nan,
        "AUROC std":    round(aurocs.std(),    3) if len(aurocs) > 0 else np.nan,
        "AUROC min":    round(aurocs.min(),    3) if len(aurocs) > 0 else np.nan,
        "AUROC max":    round(aurocs.max(),    3) if len(aurocs) > 0 else np.nan,
    }


cond_rows = []
for lbl in ["A", "B"]:
    sub = sel_df[sel_df["label"] == lbl] if len(sel_df) > 0 else pd.DataFrame()
    cond_rows.append(_cond_stats(lbl, sub, "cpds", "positives", "auroc"))

m_sub = sel_df[sel_df["label"] == "M"] if len(sel_df) > 0 else pd.DataFrame()
cond_rows.append(_cond_stats("M", m_sub, "cpds", "positives", "auroc"))

if len(general_model) > 0:
    g_row = _cond_stats("G — General (step 20)", general_model, "n_compounds", "n_actives", "auroc")
    g_row["% using middle cutoff"] = 100.0
    cond_rows.append(g_row)

sheet2 = pd.DataFrame(cond_rows)

# ---- Sheet 3: By (activity_type, unit) pair ----

print("  Sheet 3: By Activity-Unit...")

all_pairs = set()
for df in [datasets_12, indiv_selected, merged_lm, general_model]:
    if len(df) > 0 and "activity_type" in df.columns and "unit" in df.columns:
        for at, u in df[["activity_type", "unit"]].drop_duplicates().values:
            all_pairs.add((at, u))


def _mid_cutoff(activity_type, unit, target_type="ORGANISM"):
    key = (activity_type, unit, target_type, pathogen_code)
    cl = expert_cutoffs.get(key)
    if not cl:
        return np.nan
    return cl[1] if len(cl) >= 2 else cl[0]


assay_counts = {}
if len(datasets_12) > 0:
    for (at, u), grp in datasets_12.groupby(["activity_type", "unit"], dropna=False):
        org = grp[grp["target_type_curated_extra"] == "ORGANISM"]["assay_id"].nunique()
        sp  = grp[grp["target_type_curated_extra"] == "SINGLE PROTEIN"]["assay_id"].nunique()
        assay_counts[(at, u)] = (org, sp)


def _selected_for_pair(label_filter, at, u):
    if len(sel_df) == 0:
        return pd.DataFrame()
    unit_mask = sel_df["unit"].eq(u) if isinstance(u, str) else sel_df["unit"].isna()
    sub = sel_df[(sel_df["activity_type"] == at) & unit_mask]
    if label_filter:
        sub = sub[sub["label"] == label_filter]
    return sub


pair_rows = []
for activity_type, unit in sorted(all_pairs, key=lambda x: (x[0], "" if pd.isna(x[1]) else str(x[1]))):
    unit_str = str(unit) if not (isinstance(unit, float) and np.isnan(unit)) else ""
    mid = _mid_cutoff(activity_type, unit)
    org_count, sp_count = assay_counts.get((activity_type, unit), (0, 0))
    sel_a = _selected_for_pair("A", activity_type, unit)
    sel_b = _selected_for_pair("B", activity_type, unit)
    sel_m = _selected_for_pair("M", activity_type, unit)
    if len(general_model) > 0 and "activity_type" in general_model.columns:
        unit_mask_g = general_model["unit"].eq(unit) if isinstance(unit, str) else general_model["unit"].isna()
        g_row_df = general_model[(general_model["activity_type"] == activity_type) & unit_mask_g]
    else:
        g_row_df = pd.DataFrame()
    all_sel = pd.concat([sel_a, sel_b, sel_m], ignore_index=True)
    best_auroc = all_sel["auroc"].max() if len(all_sel) > 0 else np.nan
    best_label = all_sel.loc[all_sel["auroc"].idxmax(), "label"] if len(all_sel) > 0 else ""
    pair_rows.append({
        "activity_type": activity_type,
        "unit": unit_str,
        "middle_cutoff": mid,
        "N ORGANISM assays (total)": org_count,
        "N SINGLE PROTEIN assays (total)": sp_count,
        "N selected A": len(sel_a),
        "N selected B": len(sel_b),
        "N selected M": len(sel_m),
        "N selected G (general)": len(g_row_df),
        "Best AUROC (A/B/M)": round(best_auroc, 3) if not (isinstance(best_auroc, float) and np.isnan(best_auroc)) else np.nan,
        "Best condition": best_label,
        "Compounds in best dataset": int(all_sel.loc[all_sel["auroc"].idxmax(), "cpds"]) if len(all_sel) > 0 else 0,
        "Actives in best dataset": int(all_sel.loc[all_sel["auroc"].idxmax(), "positives"]) if len(all_sel) > 0 else 0,
        "General model AUROC": round(g_row_df["auroc"].iloc[0], 3) if len(g_row_df) > 0 else np.nan,
        "General model compounds": int(g_row_df["n_compounds"].iloc[0]) if len(g_row_df) > 0 else 0,
        "General model actives": int(g_row_df["n_actives"].iloc[0]) if len(g_row_df) > 0 else 0,
    })

sheet3 = pd.DataFrame(pair_rows)

# ---- Sheet 4: Rejection Reasons ----

print("  Sheet 4: Rejection Reasons...")

_CATEGORIES = {
    "selected":               "Retained in final selection",
    "already_accepted":       "already accepted",
    "non_organism":           "non-ORGANISM target type",
    "qualitative_only":       "only qualitative data",
    "no_activity_data":       "no activity data",
    "no_cutoff":              "no expert cutoff defined",
    "too_few_compounds":      "insufficient compounds",
    "too_few_positives":      r"insufficient positives|insufficient actives",
    "ratio_out_of_range":     "active ratio",
    "middle_cutoff_failure":  "middle cutoff",
    "fractional_contribution":"insufficient_fractional_contribution",
    "insufficient_compatible":"insufficient compatible assays",
    "auroc_below":            "below 0.70 threshold",
    "correlation":            "high correlation",
}
_CATEGORY_LABELS = {
    "selected":               "Selected",
    "already_accepted":       "Already accepted in prior condition",
    "non_organism":           "Non-ORGANISM target type",
    "qualitative_only":       "Qualitative data only",
    "no_activity_data":       "No activity data",
    "no_cutoff":              "No expert cutoff defined",
    "too_few_compounds":      "Too few compounds",
    "too_few_positives":      "Too few positives/actives",
    "ratio_out_of_range":     "Active ratio out of range",
    "fractional_contribution":"Insufficient fractional contribution (merging)",
    "middle_cutoff_failure":  "Middle cutoff failure",
    "fractional_contribution":"Insufficient fractional contribution (discarded after rescue pass)",
    "insufficient_compatible":"Insufficient compatible assays for merging",
    "auroc_below":            "AUROC below 0.70 threshold",
    "correlation":            "Removed: high correlation with higher-priority dataset",
    "other":                  "Other / unclassified",
}


def _parse_cats(series):
    assigned = pd.Series(False, index=series.index)
    cats = {}
    for name, pattern in _CATEGORIES.items():
        regex = "|" in pattern
        mask = series.str.contains(pattern, na=False, regex=regex) & ~assigned
        cats[name] = mask
        assigned |= mask
    cats["other"] = ~assigned
    return cats


rejection_rows = []
if len(assays_master) > 0:
    n_total = len(assays_master)
    for cond_label, col in [("A", "comment_A"), ("B", "comment_B"), ("M", "comment_M")]:
        if col not in assays_master.columns:
            continue
        cats = _parse_cats(assays_master[col])
        for cat_key, mask in cats.items():
            n = int(mask.sum())
            if n == 0:
                continue
            rejection_rows.append({
                "Condition": cond_label,
                "Category": _CATEGORY_LABELS.get(cat_key, cat_key),
                "N assays": n,
                "% of all assays": round(100 * n / n_total, 1),
            })

sheet4 = pd.DataFrame(rejection_rows) if rejection_rows else pd.DataFrame(
    columns=["Condition", "Category", "N assays", "% of all assays"]
)

# ---- Write workbook ----

xl_path = os.path.join(OUTPUT, "21_summary_table.xlsx")
with pd.ExcelWriter(xl_path, engine="openpyxl") as writer:
    sheet1.to_excel(writer, sheet_name="1. Pipeline Funnel",    index=False)
    sheet2.to_excel(writer, sheet_name="2. By Condition",       index=False)
    sheet3.to_excel(writer, sheet_name="3. By Activity-Unit",   index=False)
    sheet4.to_excel(writer, sheet_name="4. Rejection Reasons",  index=False)

    for sheet_name, df in [
        ("1. Pipeline Funnel",   sheet1),
        ("2. By Condition",      sheet2),
        ("3. By Activity-Unit",  sheet3),
        ("4. Rejection Reasons", sheet4),
    ]:
        ws = writer.sheets[sheet_name]
        for col_idx, col in enumerate(df.columns, start=1):
            max_len = max(
                len(str(col)),
                df[col].astype(str).str.len().max() if len(df) > 0 else 0,
            )
            ws.column_dimensions[ws.cell(1, col_idx).column_letter].width = min(max_len + 2, 60)

print(f"Saved summary table -> {xl_path}")
print(f"  Sheet 1 (Pipeline Funnel):   {len(sheet1)} rows")
print(f"  Sheet 2 (By Condition):      {len(sheet2)} rows")
print(f"  Sheet 3 (By Activity-Unit):  {len(sheet3)} rows")
print(f"  Sheet 4 (Rejection Reasons): {len(sheet4)} rows")
