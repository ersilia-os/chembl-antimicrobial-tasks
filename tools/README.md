# tools/

On-demand helper / reporting utilities. These are **not** part of the numbered `scripts/` pipeline or
`run_all.sh` — run them manually. All are **read-only** over `config/` and `output/` **except**
`build_strain_equivalences.py`, which writes a draft config file for review.

They read from the current stage layout: **Stage 2** `output/stage2/<code>/` (07–11) + `stage2/general/`
(12); **Stage 3** `output/stage3/<code>/` (13–20); **Stage 4** `output/stage4/<code>/` (21 curation ·
22 binarised · 23 transfer pools · 24 CV · 25 merged pools · 26 low-data) + `stage4/general_plots/` (27).

## Wired into the pipeline

### `build_strain_equivalences.py`
Builds a **draft strain-equivalence table** for one pathogen from the **BacDive API** (authoritative
culture-collection cross-references), for human review — the alternative to mining ChEMBL co-occurrence
(which over-merges). Its output is **consumed at runtime** by `strain_norm_key` / `load_strain_equivalences`
(the assay-merge grouping key, used by steps 16, 19, 21).
- Run: `python tools/build_strain_equivalences.py <pathogen_code>`
- Reads: `output/stage2/<code>/09_assays_parameters_full.csv` + BacDive.
- Writes: `config/strains/strain_equivalences_<code>.csv` (draft; the manual seed
  `config/strains/strain_equivalences_manual.csv` takes precedence).

## Validation / QC

### `validate_datasets.py`
Integrity validator for Stage-4 steps **21** (build curated) & **22** (binarise): confirms no rows are
lost, no assay cross-contamination, correct InChIKey dedup (keeps the most-active record), and that
binarisation matches the stated direction/cutoff rule. Collects **all** violations (never crashes on the
first), prints per-pathogen worked examples, and exits non-zero on any failure (CI-gateable).
- Run: `python tools/validate_datasets.py [codes...]` (default: all pathogens).
- Reads: `stage2/<code>/` (07/08) + `stage4/<code>/` (21/22).
- Writes: `docs/<date>_step21-22-validation.md` (PASS/FAIL report).

## Result summaries

### `pool_summary_report.py`
Cross-pathogen summary of the **final Stage-4 pools**: pool counts, chemical-space coverage, activity
ratios, and per-pool CV AUROC (the grown-pool AUROC when the step-25 merge batch was accepted, else the
step-24 baseline; step-26 catch-alls folded in).
- Run: `python tools/pool_summary_report.py`
- Reads: `stage4/<code>/` `25_pool_summary.csv`, `25_pools/`, `23_chemspace_partition.csv`, `26_cv_summary.csv`.
- Writes: `docs/2026-07-03_final-pool-summary.csv` (per pool) + `-bycategory.csv`.

### `pool_quality_audit.py`
Per-pathogen audit of **final Stage-4 pools that do not reach the AUROC threshold** (`AUROC_MIN_THRESHOLD`
= 0.70 from `src/default.py`). For each pathogen it counts sub-threshold pools (step-25 shipped pools +
step-26 catch-alls dropped with `low_auroc`) and reports their chemical-space coverage two ways —
**gross** (union of their molecules) and **exclusive** (those in no ≥-threshold pool anywhere in the
pathogen) — each against the category and the total-pathogen chem-space. Regenerates in place.
- Run: `python tools/pool_quality_audit.py`
- Reads: `stage4/<code>/` `{25_pool_summary,24_cv_summary,26_cv_summary}.csv`, `25_pools/`, `23_pools/`,
  `26_pools/`, `22_binarised/`.
- Writes: `docs/pool_quality_audit.md` + `output/stage4/general_plots/pool_quality_audit.csv`.

### `generate_pipeline_report.py`
Cross-pathogen status snapshot of the **Stage-3 conservative (A/B/M) branch**: one row per pathogen of how
molecules and assays flow (assessed → SINGLE PROTEIN / ORGANISM / discarded → final selected A/B/M
datasets). Intended as a before/after snapshot for tracking pipeline changes. Does **not** cover the
Stage-4 signal-based branch.
- Run: `python tools/generate_pipeline_report.py [output_csv]` (default `docs/260622_report.csv`).
- Reads: `stage2/<code>/` (07/08) + `stage3/<code>/` (`13_datasets`, `18_final_datasets`, `19_assays_master`).

## Diagnostics & threshold analyses (one-off; re-runnable)

### `near_dup_report.py`
Quantifies how much step-23 near-duplicate folding reduces the dataset count and sweeps the
(overlap, discordance) cutoffs — the analysis behind the **0.90 / 0.20** choice.
- Run: `python tools/near_dup_report.py`
- Reads: `stage4/<code>/22_binarised/`. Writes: `docs/2026-07-02_neardup-reduction.csv` (+ `.png`).

### `singleton_inclusion_report.py`
Sweeps the singleton chemical-space-fraction threshold (1 % / 2 % / 5 %) to inform step-24's
singleton-promotion cutoff (the alternative that was rejected for the flat **≥ 200-molecule** floor).
- Run: `python tools/singleton_inclusion_report.py`
- Reads: `stage4/<code>/` `23_pool_summary.csv`, `23_chemspace_partition.csv`. Writes:
  `docs/2026-07-02_singleton-inclusion.csv`.

### `evaluate_censoring_interference.py`
Quantifies how often censored-measurement handling flips a binary label near the activity cutoff
(ambiguity straddle + face-value change) — the evidence that motivated the **conservative "prove-the-hit"**
binarisation in step 22.
- Run: `python tools/evaluate_censoring_interference.py`
- Reads: `stage2/<code>/` (07/08) + `stage4/<code>/` (`21_curation_summary`, `22_binarised`). Writes:
  `docs/<date>_censoring-interference.md` + `.csv`.

### `pooling_conflict_report.py`
Per Stage-4 transfer pool (step 23): how many pools, how many "conflicts before handling" (molecules
tested in ≥ 2 of a pool's member datasets with disagreeing raw labels), and how they resolve after the
merge — read directly from the pool files (does not re-implement the resolution).
- Run: `python tools/pooling_conflict_report.py [codes...]`
- Reads: `stage4/<code>/` `23_pool_members.csv`, `23_pools/`, `22_binarised/`. Writes:
  `docs/<date>_step23-conflict-perpool.csv` + `-summary.csv`.

### `redundant_conflict_audit.py`
Whether step-25 `redundant` leftovers (dropped because all their molecules are already in a modelled
pool) silently lose signal, and whether the step-25 pool-growth policy itself changes labels. Three
sub-analyses: **conflicts** (raw leftover-vs-pool label comparison — consensus/strict/any + a driver
breakdown by activity_type), **merge-flip** (merge each leftover into every containing pool and
re-binarise via `merge_datasets`; reports flip_vs_delivered / flip_marginal / delivered_ne_remerge),
**reaverage** (pool-wide, de-duplicated: what fraction of shipped labels would flip if step 25
re-averaged its members instead of "base wins"). Findings: raw conflict 38 % (MIC-driven) but only
~17 % would change a shipped label (~10 % marginal), and re-averaging pool-wide changes just 0.04 %.
- Run: `python tools/redundant_conflict_audit.py [all|conflicts|merge-flip|reaverage] [codes...]`
- Reads: `stage4/<code>/` `25_merge_log.csv`, `25_pool_members.csv`, `25_pools/`, `22_binarised/`,
  `22_binarisation_summary.csv`. Writes: `docs/2026-07-15_redundant-{conflicts,mergeflip,reaverage}_*.csv`
  + `.png`. Write-up: `docs/2026-07-15_redundant-conflict-diagnostics.md`.
