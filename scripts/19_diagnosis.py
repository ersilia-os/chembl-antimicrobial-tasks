"""
Step 19 — Diagnosis plots.

Produces a 3×3 diagnostic figure for a given pathogen and saves it to
output/<pathogen_code>/19_diagnosis.png.

Usage
-----
    python scripts/19_diagnosis.py <pathogen_code>
"""

from collections import Counter, defaultdict
from scipy.stats import gaussian_kde
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib as mpl
import pandas as pd
import numpy as np
import stylia
import h5py
from tqdm import tqdm
import sys
import os

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import load_pathogen

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print("Step 19: Diagnosis plots")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

print("Loading data...")

raw = pd.read_csv(os.path.join(OUTPUT, "07_chembl_raw_data.csv.gz"), low_memory=False)
cleaned = pd.read_csv(os.path.join(OUTPUT, "08_chembl_cleaned_data.csv.gz"), low_memory=False)
assays_raw_df = pd.read_csv(os.path.join(OUTPUT, "07_assays_raw.csv"))
assays_cleaned_df = pd.read_csv(os.path.join(OUTPUT, "08_assays_cleaned.csv"))
compound_counts = pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))
final_datasets = pd.read_csv(os.path.join(OUTPUT, "17_final_datasets.csv"))
correlations = pd.read_csv(os.path.join(OUTPUT, "17_dataset_correlations.csv"))

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

total_organism = sum(organism_ctr.values())
assay_organism_bars = [organism_ctr[pathogen] / total_organism, 1][::-1]

total_unit = sum(unit_ctr.values())
unit_bars = [
    unit_ctr["umol.L-1"] / total_unit,
    (unit_ctr["umol.L-1"] + unit_ctr["%"]) / total_unit,
    1,
][::-1]

# ---------------------------------------------------------------------------
# Compound occurrence rank plot
# ---------------------------------------------------------------------------

count_cpds = sorted(Counter(cleaned["compound_chembl_id"]).values(), reverse=True)

# ---------------------------------------------------------------------------
# Cumulative chemical space coverage per assay
# ---------------------------------------------------------------------------

print("Computing cumulative coverage...")

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

del cleaned

# ---------------------------------------------------------------------------
# Final dataset coverage per label
# ---------------------------------------------------------------------------

def _parse_assay_key(s):
    parts = s.split("|")
    return (parts[0], parts[1], np.nan if parts[2] == "" else parts[2])


final_coverage = {label: set() for label in "ABM"}
for label, assay_keys in final_datasets[final_datasets["selected"]][
    ["label", "assay_keys"]
].values:
    for s in assay_keys.split(";"):
        final_coverage[label] |= assay_to_compounds.get(_parse_assay_key(s), set())

# ---------------------------------------------------------------------------
# Correlation matrices (clustered by compound overlap)
# ---------------------------------------------------------------------------

print("Building correlation matrices...")

names = sorted(set(correlations["name_1"]))
co_dict = {(r.name_1, r.name_2): r.compound_overlap for r in correlations.itertuples()}
ho_dict = {(r.name_1, r.name_2): r.hit_overlap_1000 for r in correlations.itertuples()}

X_co = np.array([[co_dict[(n1, n2)] for n2 in names] for n1 in names])
X_ho = np.array([[ho_dict[(n1, n2)] for n2 in names] for n1 in names])

Z = linkage(squareform(1 - X_co, checks=False), method="average")
idx = leaves_list(Z)
X_co = X_co[np.ix_(idx, idx)]
X_ho = X_ho[np.ix_(idx, idx)]

# ---------------------------------------------------------------------------
# tSNE
# ---------------------------------------------------------------------------

from sklearn.decomposition import PCA
from openTSNE import TSNE as openTSNE

print("Loading ECFPs for tSNE...")
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
print("  PCA...")
X_tsne = PCA(n_components=16, random_state=42, svd_solver="randomized").fit_transform(X_tsne)
print("  tSNE...")
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


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

print("Plotting...")

fig, axs = stylia.create_figure(3, 3, width=stylia.TWO_COLUMNS_WIDTH * 1.5, height=stylia.TWO_COLUMNS_WIDTH * 1.35)
cmap2 = mpl.colors.LinearSegmentedColormap.from_list("purple_blue", ["#50285A", "#8DC7FA"], N=256)
label_to_color = {"A": "#FAA08B", "B": "#8DC7FA", "M": "#FAD782"}

# [0][0] Raw vs cleaned: activities, compounds, assays
ax = axs.next()
ax.bar([0, 0], [activities_raw, activities_cleaned], color=["#D2D2D2", "#8DC7FA"], ec="k", zorder=2)
ax.bar([1, 1], [compounds_raw, compounds_cleaned], color=["#D2D2D2", "#8DC7FA"], ec="k", zorder=2)
ax.bar([2, 2], [n_assays_raw, n_assays_cleaned], color=["#D2D2D2", "#8DC7FA"], ec="k", zorder=2)
ax.set_yscale("log")
ax.set_yticks([10**1, 10**2, 10**3, 10**4, 10**5, 10**6])
ax.set_ylim([5, 3e6])
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(["Activities", "Compounds", "Assays"])
stylia.label(ax, ylabel="Number of")

# [0][1] Fraction bars: assay type, target type, organism, unit, activity type
ax = axs.next()
ax.bar([0, 0, 0], assay_type_bars, zorder=2, color=["#D2D2D2", "#DC9FDC", "#FAD782"], ec="k")
ax.text(0, assay_type_bars[-1] / 2, "F", ha="center", va="center")
ax.text(0, (assay_type_bars[-2] + assay_type_bars[-1]) / 2, "B", ha="center", va="center")
ax.bar([1, 1, 1], target_type_bars, zorder=2, color=["#D2D2D2", "#DC9FDC", "#FAD782"], ec="k")
ax.text(1, target_type_bars[-1] / 2, "ORG", ha="center", va="center")
ax.text(1, (target_type_bars[-2] + target_type_bars[-1]) / 2, "SP", ha="center", va="center")
ax.bar([2, 2], assay_organism_bars, zorder=2, color=["#D2D2D2", "#FAD782"], ec="k")
ax.bar([3, 3, 3], unit_bars, zorder=2, color=["#D2D2D2", "#DC9FDC", "#FAD782"], ec="k")
ax.text(3, unit_bars[-1] / 2, "uM", ha="center", va="center")
ax.text(3, (unit_bars[-2] + unit_bars[-1]) / 2, "%", ha="center", va="center")
ax.bar(
    [4, 4, 4],
    [1, (act_qt + act_mx) / activities_cleaned, act_qt / activities_cleaned],
    zorder=2, color=["#D2D2D2", "#DC9FDC", "#FAD782"], ec="k",
)
ax.text(4, (act_qt / activities_cleaned) / 2, "qt.", ha="center", va="center")
ax.set_ylim([0, 1])
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_xticklabels(["Assay\ntype", "Target\ntype", "Organism", "Unit", "Activity\ntype"])
stylia.label(ax, ylabel="Fraction")

# [0][2] tSNE
ax = axs.next()
cmap_tsne = mpl.colors.LinearSegmentedColormap.from_list("purple_yellow", ["#50285A", "#FAD782"], N=256)
ax.scatter(emb[n_ref:][:, 0], emb[n_ref:][:, 1], s=5, c="lightgray")
scatter_density(ax, emb[:n_ref, :2], size_min=0, size_max=20, cmap=cmap_tsne, alpha=1)
ax.set_yticks([])
ax.set_xticks([])
stylia.label(ax, xlabel="tSNE-2", ylabel="tSNE-1")

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
stylia.label(ax, xlabel="Compound number", ylabel="Compound occurrences")

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
ax.plot(x_cpds_arr, cpds_per_assay, c="#FAD782", lw=1.8, zorder=3)
ax.fill_between(x_cpds_arr, cpds_per_assay, 1, color="#FAD782", alpha=0.6, zorder=2)
for xi in x_marks:
    ax.scatter([xi], [cpds_per_assay[xi - 1]], zorder=4, ec="k", s=25, c="#D2D2D2")
ymax = max(max(cpds_per_assay), np.max(cum_as_compounds), left_to_right(1.0))
ax.set_ylim(0.6, ymax * 2)
secax = ax.secondary_yaxis("right", functions=(right_to_left, left_to_right))
secax.set_ylabel("Cumulative fraction of\nthe chemical space")
secax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
secax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.1f"))
ax.plot(x_cum, cum_as_compounds, c="#50285A", lw=1.8, zorder=5)
ax.fill_between(x_cum, cum_as_compounds, 1, color="#50285A", alpha=0.25, zorder=1)
for xi in x_marks:
    ax.scatter([xi], [cum_as_compounds[xi - 1]], zorder=6, ec="k", s=25, c="#D2D2D2")
stylia.label(ax, xlabel="Number of assays", ylabel="Number of compounds")

# [1][2] Compound overlap heatmap (clustered)
ax = axs.next()
im = ax.imshow(X_co, vmin=X_co.min(), vmax=X_co.max(), cmap=cmap2)
ax.set_xticks([])
ax.set_yticks([])
fig.colorbar(im, ax=ax, fraction=0.045).set_label("Compound overlap")
stylia.label(ax, title="Compound overlap")

# [2][0] Hit overlap heatmap (clustered)
ax = axs.next()
im2 = ax.imshow(X_ho, vmin=X_ho.min(), vmax=X_ho.max(), cmap=cmap2)
ax.set_xticks([])
ax.set_yticks([])
fig.colorbar(im2, ax=ax, fraction=0.045).set_label("Hit overlap")
stylia.label(ax, title="Top-1000 hit overlap")

# [2][1] Chemical space coverage by label
ax = axs.next()
all_coverage = final_coverage["A"] | final_coverage["B"] | final_coverage["M"]
ax.set_xlim([0, 5])
ax.set_ylim([0, 1])
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(["A", "B", "M", "ALL"])
for xi, lbl, color in [(1, "A", "#FAA08B"), (2, "B", "#8DC7FA"), (3, "M", "#FAD782")]:
    frac = len(final_coverage[lbl]) / len(pathogen_compounds)
    ax.bar([xi, xi], [1, frac], zorder=2, color=["#D2D2D2", color], ec="k")
frac_all = len(all_coverage) / len(pathogen_compounds)
ax.bar([4, 4], [1, frac_all], zorder=2, color=["#D2D2D2", "#BEE6B4"], ec="k")
stylia.label(ax, ylabel="Chemical space coverage")

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

fig.suptitle(pathogen, size=15, y=1.01)

outpath = os.path.join(OUTPUT, "19_diagnosis.png")
stylia.save_figure(outpath)
print(f"Saved -> {outpath}")
