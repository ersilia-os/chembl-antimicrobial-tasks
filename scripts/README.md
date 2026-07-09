# Scripts

The pipeline is organised into four stages. **Stage 1** (00–06) processes the raw ChEMBL database once
into a clean activity table. **Stage 2** (07–12) prepares per-pathogen assay data and features. **Stages 3
and 4 are two parallel pooling branches** off Stage 2: **Stage 3** (13–20) is conservative and
individual-dataset oriented; **Stage 4** (21–28) is a signal-based branch that pools datasets by ML
transfer signal. Each stage writes to its own bucket: `output/stage1/`, `output/stage2/<code>/` (+
`stage2/general/`), `output/stage3/<code>/`, `output/stage4/<code>/` (+ `stage4/general_plots/`).

Each step's inputs/outputs are documented in its script docstring; this file gives a one-line description
plus the key decisions/thresholds per step. Full design rationale, methodology, and reference tables (the
LLM two-track spec, strain-normalisation design, merge failure-reason table, diagnosis panels, pooling
rationale, …) are in [`docs/pipeline_design.md`](../docs/pipeline_design.md). Project-wide constants live
in `src/default.py`.

## Setup

Create a local ChEMBL PostgreSQL database following [`docs/install_ChEMBL.md`](../docs/install_ChEMBL.md)
(connection defaults in `src/default.py`), then set up the environment:

```sh
conda create -n camt python=3.12 -y
conda activate camt
pip install -r requirements.txt
```

**Environments:** most steps run in `camt`. Step 09 needs a GPU + a local [ollama](https://ollama.com/)
instance; step 12 runs the Ersilia `eos1klk` model in the `ersilia` environment.

---

## Stage 1 — Process ChEMBL (run once, ~4–5 h) → `output/stage1/`

```sh
bash scripts/run_stage1.sh                     # standard
bash scripts/run_stage1.sh --calculate_ecfps   # also compute fingerprints (step 06)
```

**Manual curation pause:** step 01 writes `output/stage1/01_activity_std_units_converted.csv`; a curator
fills `manual_curation_direction` (`1` = higher→more active, `-1` = lower→more active, `0` = unclear) and
saves it as `config/activity_std_units_manual_curation.csv`. `run_stage1.sh` checks for this file
up front and exits before the multi-hour steps if it is missing.

| Step | Script | What it does | Key decisions |
|------|--------|--------------|---------------|
| 00 | `00_export_chembl_activities.py` | Exports raw ChEMBL tables and derives activity/unit/text frequency + assay/target reference tables. | Source is the local ChEMBL PostgreSQL DB. |
| 01 | `01_prepare_manual_files.py` | Harmonises activity-type strings and maps units to UCUM; collapses synonyms. | Uses `ucum_manual.csv` + `synonyms.csv`. **Pauses for manual direction curation.** |
| 02 | `02_get_compound_info.py` | Maps `molregno` → `chembl_id` and computes molecular weight. | — |
| 03 | `03_standardize_compounds.py` | Canonicalises SMILES, strips salts/solvents, extracts parent (ChEMBL structure pipeline). | Pipeline bottleneck (~3–4 h). |
| 04 | `04_merge_activity_and_compounds.py` | Joins activities + assays + targets + standardised compounds into one table. | — |
| 05 | `05_clean_activities.py` | Filters invalid records, flags comments/text as active/inactive, UCUM-converts values, computes pChEMBL. | Missing `relation` treated as `=`; pChEMBL clipped to `[1, 9]`. |
| 06 | `06_calculate_ecfps.py` *(optional)* | ECFP6 fingerprints for all standardised compounds (HDF5). | **radius 3, 2048 bits.** |

---

## Stage 2 — Per-pathogen preparation (07–12) → `output/stage2/<code>/`

```sh
bash scripts/run_stage2.sh                 # all pathogens in config/pathogens.csv
bash scripts/run_stage2.sh mtuberculosis   # a subset
```

`run_stage2.sh` runs the `camt`-env steps **07, 08, 10, 11** per pathogen. **Steps 09 and 12 are
excluded** (different environments) and must be run manually — the script prints a reminder:
- `python scripts/09_curate_assay_parameters.py <code>` — needs **GPU + ollama**; required before Stages 3/4.
- `conda run -n ersilia python scripts/12_calculate_eos1klk.py` — cross-pathogen (`output/stage2/general/`);
  required for the Stage-4 embedding figures.

| Step | Script | What it does | Key decisions |
|------|--------|--------------|---------------|
| 07 | `07_get_pathogen_assays.py` | Selects all ChEMBL records for the pathogen by organism-name match, plus a manual allow/deny list. | `config/assays_to_include|exclude/<code>.csv` (optional). |
| 08 | `08_clean_pathogen_activities.py` | Assigns biological direction per (activity_type, unit) and drops unmodelable records. | Direction from `activity_std_units_manual_curation.csv`; keep only rows with a direction or a text flag. |
| 09 | `09_curate_assay_parameters.py` | LLM extracts curated assay context (target type/name/ID, strain, ATCC, mutations, media). | **GPU + ollama** (`LLM_MODEL`); batches of 6; resumable; target IDs resolved by exact match, never guessed. |
| 10 | `10_calculate_assay_clusters.py` | Clusters each assay's compounds (BitBirch on ECFP4) for a diversity measure. | **Tanimoto 0.3 / 0.6 / 0.85.** Independent of 09 (parallelisable). |
| 11 | `11_get_assay_overlap.py` | Pairwise compound overlap between assays. | Only assays with **≥ 50 compounds**. |
| 12 | `12_calculate_eos1klk.py` | Runs the Ersilia `eos1klk` model once over the cross-pathogen SMILES union → embeddings. | Pooled & predicted once (~556k unique vs ~897k summed); `ersilia` env; `BATCH_SIZE=1000`. |

---

## Stage 3 — Conservative pooling (13–20) → `output/stage3/<code>/`

The conservative, individual-dataset branch: expert-cutoff binarisation → per-assay baseline models →
merge small assays → deduplicate → master table + diagnosis. Runs per pathogen; a **parallel** branch to
Stage 4 (not a prerequisite), so it's **opt-in** — not part of `run_all.sh`. Needs step 09 output.

```sh
bash scripts/run_stage3.sh                 # all pathogens (13→20 each)
bash scripts/run_stage3.sh mtuberculosis   # a subset
```

| Step | Script | What it does | Key decisions |
|------|--------|--------------|---------------|
| 13 | `13_prepare_assay_datasets.py` | Binarises each assay with expert cutoffs into quantitative/qualitative/mixed datasets. | Cutoffs from `expert_cutoffs.csv`; COADD assays excluded; `%` endpoints use cutoffs `50;75;90`. |
| 14 | `14_lightmodel_individual.py` | RandomForest (100 trees, 4-fold CV, Morgan) per dataset → AUROC. | Condition A (≥1000 cpds, ≥50 pos, ratio<0.5) vs B (≥100 pos, ratio≥0.5, decoys to 10%); reference set = 10 000 pathogen-naive compounds (seed 42). |
| 15 | `15_select_datasets_individual.py` | Picks the best cutoff per assay. | **AUROC > 0.70** (`AUROC_MIN_THRESHOLD`); prefer reference cutoff, switch to best only if +0.10 (`AUROC_IMPROVEMENT_THRESHOLD`). |
| 16 | `16_lightmodel_merged.py` | Merges context-compatible assays to rescue small ones, then re-models. | Group key = activity_type/unit/direction/assay_type/target_type[/strain]; ≥2 assays, ≥5% fractional contribution each, post-merge ≥1000/≥100 cpds & >50 pos. |
| 17 | `17_select_datasets_merged.py` | Picks the best cutoff per merged group. | Same AUROC/mid-cutoff rule as step 15; single-cutoff groups use the only cutoff. |
| 18 | `18_evaluate_correlations.py` | Similarity (Spearman + hit-overlap + compound-overlap) between ORGANISM models on the reference set, then greedy dedup. | Discard iff mean model-correlation **> 0.50 AND** compound_overlap **> 0.50** (`CORRELATION_THRESHOLD`). |
| 19 | `19_prepare_assay_master.py` | Master per-assay table merging all step metadata + pipeline status + final selection. | — |
| 20 | `20_diagnosis.py` *(optional)* | 3×3 diagnostic figure + 4-sheet summary workbook. | 2×2 layout if nothing selected yet. |

---

## Stage 4 — Signal-based pooling (21–28) → `output/stage4/<code>/`

The signal-based branch: rebuild fewer/larger curated datasets → conservative binarisation → pool by ML
transfer signal → cross-validate → fold in leftovers and low-data catch-alls, then a terminal
general/reporting block. Independent of Stage 3 (reads only Stage-2 outputs). All figures share one
publication house style (`src/plot_base.py`) and semantic palette (`src/plot_colors.py`).

```sh
bash scripts/run_stage4.sh                 # all pathogens
bash scripts/run_stage4.sh mtuberculosis   # a subset
```

`run_stage4.sh` runs the per-pathogen data chain **21 → 22 → 23 → 24 → 25** for each code, then the
**general/reporting block** once over all pathogens on disk: **26** (low-data catch-alls; loops pathogens
internally), **27** (ALL cross-pathogen figures), **28** (dataset trace/flow audit). Per-pathogen figures
stay with their step (transfer heatmaps → 23, cv overview → 24, merge overview → 25); every cross-pathogen
figure is consolidated in **27**. Embedding figures need step 12 output (skipped gracefully if absent).

**Refreshing figures without recomputing:** steps 23/24/25 accept `--plot-only` (redraw their per-pathogen
figures from cached `.npz` / CV predictions / pool files); **27** is plot-only by nature — re-run it any
time to refresh every cross-pathogen figure from the on-disk summaries.

| Step | Script | What it does | Key decisions |
|------|--------|--------------|---------------|
| 21 | `21_build_curated_datasets.py` | Rebuilds one curated dataset per `(assay_id, activity_type)`: InChIKey dedup (most-active), ORGANISM + quantitative only. | Discard pairs with **≤ 5 molecules** (`MIN_MOLECULES`); exclude COADD + PubChem-superseded; recover text-only qualitative *inactives*; category (DR/SP/UNCATEGORIZED) from `activity_type_categories.csv`. |
| 22 | `22_binarise_datasets.py` | Conservative "prove-the-hit" binarisation of DR/SP datasets; also writes `22_cutoff_sensitivity.csv` (active ratio at each low/middle/high candidate per category/unit — diagnostic only). | Cutoff per unit from `category_cutoffs.csv` (`µM→10`, `%→50`, `mm→20`), **overridable per pathogen** via `category_cutoffs_overrides.csv` (tier picked from `category_cutoff_candidates.csv`); the overrides file holds per-pathogen `(pathogen_code, unit, tier)` entries — currently `pfalciparum` and `spneumoniae` use the DR **low** tier (**5 µM**, signed off), and pathogens not listed fall back to the global default. Changing or adding an override is a threshold decision requiring sign-off. Direction from the dataset; UNCATEGORIZED skipped; degenerate kept + flagged. |
| 23 | `23_transfer_pooling.py` | Trains RF on each dataset A, scores AUROC on B's disjoint molecules → pools by greedy-modularity communities, merged to 1 row/molecule. Per-pathogen figures: transfer heatmap + `23_pooling_overview.png` (folding / pool vs singleton / pool active ratio / conservative-rebinarisation effect); cross-pathogen pool figures in step 27. | Near-dup fold (containment **≥0.90** / discordance **<0.20**); **bidirectional, `POOL_THRESHOLD=0.75`**; `MIN_DATASET_SIZE=20`; `MIN_TEST=5`; `.npz` cache reused manually (`--force` to recompute). |
| 24 | `24_crossvalidate_pools.py` | 5-fold CV (RF/Morgan) per pool → AUROC + Youden cutoff. | Model a pool (multi-dataset **or** singleton, symmetric) iff **≥100 molecules** (`MIN_MODELLABLE_MOLECULES`) and CV-able (≥5/class); smaller pools deferred to step 25. |
| 25 | `25_merge_leftovers.py` | Folds leftover datasets into the best-matching modelled pool, then batch re-CVs. | Assign by balanced accuracy at the pool's Youden cutoff, but merge only if that balanced accuracy ≥ **`MIN_ASSIGN_BALACC=0.60`** (else deferred, low fit); accept a batch iff grown AUROC ≥ baseline − **`AUROC_TOL=0.05`**. |
| 26 | `26_lowdata_pools.py` | For categories with 0 modelled pools, merges everything into one catch-all and CVs it. Per-pathogen figures `26_catchall_overview.png` (AUROC / active-inactive / datasets merged) plus final-pool embedding figures `26_eos1klk_pools_{DR,SP}.png`; the cross-pathogen overview is also in step 27. | Keep iff **≥ 50 molecules** (`CATCHALL_MIN_MOL`) and CV-able (≥5/class); else `too_small`/`not_cv_able`. AUROC is **evaluated, not enforced** — `auroc_pass` flags whether it clears **`CATCHALL_MIN_AUROC=0.70`** (a reference, not a keep gate). |
| 27 | `27_general_plots.py` | **All** cross-pathogen figures in one terminal step: curation + binarisation comparisons, the cutoff-sensitivity figure/table (`27_cutoff_sensitivity.*`, aggregating each pathogen's `22_cutoff_sensitivity.csv`), pool summaries (from step 23), the step-24 pool/CV overview, the step-25 merge overview, a final-pool overlap figure, an activity-ratio figure, the low-data catch-all overview (from step 26), a cross-step pipeline summary (`27_pipeline_summary.png`), and `27_master_table.csv`. | Runs once across all pathogens with ≥2 available → `output/stage4/general_plots/`. Plot-only by nature (no model training). |
| 28 | `28_trace_sample_datasets.py` | Read-only audit: samples 10 datasets per pathogen and traces each one's full journey (stages 1/2/4) to its final endpoint (modelled pool / merged-leftover / catch-all / discarded), emitting a single self-contained interactive HTML panel (traces + population aggregates embedded; no CSVs written). | **Sampling only, no data touched.** Seeded (`RANDOM_SEED`), stratified across fate, activity-type family, size (`≤5`/`6–20`/`21–100`/`>100` mol) and active/inactive ratio (`<0.1`/`0.1–0.5`/`>0.5`). These bucket edges are illustrative sampling strata, **not** pipeline cutoffs. |

---

## Configuration files (`config/`)

| File | Step(s) | Purpose |
|------|---------|---------|
| `pathogens.csv` | all | Pathogen codes and names. |
| `ucum_manual.csv` | 01, 05 | Unit → UCUM mapping and value-conversion formulas. |
| `synonyms.csv` | 01, 05 | Collapses activity-type variants to canonical names. |
| `activity_comments_manual_curation.csv`, `standard_text_manual_curation.csv` | 05 | Map comment / text values to active / inactive / unknown. |
| `activity_std_units_manual_curation.csv` | 08 | Biological direction per (activity_type, unit). |
| `assays_to_include/<code>.csv`, `assays_to_exclude/<code>.csv` | 07 | Manual assay allow / deny lists (optional). |
| `expert_cutoffs.csv` | 13 | Stage-3 binarisation thresholds per (activity_type, unit, target_type, pathogen). |
| `pubchem_aids/<code>.csv` | 16, 19, 21 | Assays with a PubChem counterpart (modelled/handled separately). |
| `strains/strain_equivalences_*.csv` | 16, 19, 21 | Strain alias map for the normalised grouping key (built by `tools/build_strain_equivalences.py`). |
| `activity_type_categories.csv` | 21 | Explicit `(activity_type, unit) → DR/SP` map; anything unlisted → UNCATEGORIZED. |
| `category_cutoffs.csv` | 22 | Global default cutoff per (category, unit): **DR = 10 µmol/L**, **SP-% = 50 %**, **SP-zone = 20 mm**. |
| `category_cutoff_candidates.csv` | 22 | The low/middle/high candidate values a per-pathogen override may pick per unit (`µM 5;10;20`, `% 50;75;90`, `mm 10;20;30`). For µM (direction −1) a **lower** value is more stringent; for %/mm a **higher** value is. |
| `category_cutoffs_overrides.csv` | 22 | Per-pathogen `(pathogen_code, unit, tier)` cutoff overrides. Currently `pfalciparum` and `spneumoniae` use the DR **low** tier (**5 µM**, signed off); pathogens not listed fall back to the global default. Changing or adding an override is a threshold decision — do so only after reviewing `27_cutoff_sensitivity.png` and with sign-off. |

## Accessing the datasets

Two families of ML-ready datasets are produced per pathogen: the **Stage-3 conservative** datasets
(`output/stage3/<code>/`: `18_final_datasets.csv`, `13_datasets/`, `19_assays_master.csv`) and the
**Stage-4 signal-based** datasets (`output/stage4/<code>/`: `22_binarised/`, and the pooled
`23_`/`25_` outputs). See [`docs/dataset_access_guide.md`](../docs/dataset_access_guide.md) for
programmatic access.
