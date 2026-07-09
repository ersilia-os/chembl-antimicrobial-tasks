# Dataset flow audit — tracing sampled datasets through pipeline stages 1, 2 and 4

**Date:** 2026-07-08
**Script:** `scripts/28_trace_sample_datasets.py` (read-only; no data mutated, no models re-run)
**Output:** `output/stage4/general_plots/28_dataset_trace_panel.html` (self-contained interactive panel / Artifact). The 150 sampled traces and the population aggregates are computed in-memory and embedded in the panel; no intermediate CSVs are persisted.

> This is an AI-generated audit report. It **presents** what the pipeline outputs show; it does **not** decide what they mean or what (if anything) to change. Please review before it is committed.

## Why this audit
A single mislabelling case prompted it: assay `CHEMBL5163887` (a plain *Acinetobacter baumannii* MIC growth assay) was assigned `assay_organism_curated = Homo sapiens` by the step-09 LLM (`gpt-oss:20b`). To see whether that is a one-off or a pattern — and more generally to make the pipeline's per-dataset behaviour legible — we sampled 10 datasets per pathogen (15 pathogens = **150 datasets**) spanning activity types, sizes and active/inactive ratios, and traced each one's full journey.

## What a "dataset" is, and what we traced
The modelling unit is one `(assay_id, activity_type, unit)` triplet, identified as e.g. `CHEMBL5163887_MIC_uM`. We traced each through **stages 1, 2 and 4**; stage 3 (the parallel conservative pooling branch) is deliberately excluded.

The journey, stage by stage, and the file each fact comes from:

| Node | Question | Source |
|---|---|---|
| ChEMBL raw | what did ChEMBL record? | `data/chembl_activities/assays.csv`, `stage1/00_assay_descriptions.csv` |
| Stage-2 clean (07→08) | what survived unit/relation cleaning? | `stage2/<p>/07_assays_raw.csv`, `08_assays_cleaned.csv` |
| Curate (09) | what did the LLM assign (organism, target, strain)? | `stage2/<p>/09_assays_parameters_full.csv` |
| Curate (21) | how many molecules; retained or discarded (why)? | `stage4/<p>/21_curation_summary.csv` |
| Binarise (22) | how many actives / inactives at what cutoff? | `stage4/<p>/22_binarisation_summary.csv` |
| Pool (23) | which provisional transfer pool? | `stage4/<p>/23_pool_members.csv` |
| Model / merge (24/25/26) | modelled pool, merged leftover, or catch-all? | `24_cv_summary.csv`, `25_pool_members.csv`, `25_merge_log.csv`, `26_cv_summary.csv` |
| Final | the single resolved endpoint | derived |

**Fate branches per category.** Each category (DR / SP) is handled independently: a category is either *modelled* (steps 23→24→25) or *fully deferred* to the step-26 catch-all. E.g. for `abaumannii`, DR follows the modelled path while SP goes to a catch-all. The resolver handles all three observed file-coverage shapes (normal 23→25; deferred→26; fully catch-all, as in `campylobacter`) without assuming a file exists.

## How the sample was chosen
A **seeded** (`RANDOM_SEED = 42`), **coverage-greedy** picker selects 10 datasets per pathogen that jointly span: fate class, activity-type family (potency-MIC / IC-EC-K / percent / zone / qualitative), size bucket (`≤5`, `6–20`, `21–100`, `>100` molecules) and active/inactive ratio bucket (`<0.1`, `0.1–0.5`, `>0.5`), and preferentially surfaces anomalies (organism mismatch, degenerate binarisation, multi-unit). Re-running reproduces the identical sample.

> The size and ratio bucket edges are **illustrative sampling strata used only to choose examples** — they are not cutoffs applied to the data, and are exposed as editable constants at the top of the script. They are not pipeline thresholds and require no sign-off.

The 150-row sample is intentionally anomaly-biased (it exists to show variety), so **rates within the sample are not population rates**. Population-level counts below come from all 113,331 datasets, not the sample.

## Worked examples (from the sample)

**1 — the motivating mislabel, `CHEMBL5163887_MIC_uM` (abaumannii).**
Raw organism *A. baumannii* → **curated *Homo sapiens*** (`organism_mismatch = True`), yet `target_type_curated = ORGANISM` and `target_chembl_id_curated = CHEMBL614425` (*A. baumannii*) are correct. Binarised 9 active / 26 inactive at cutoff 10 µM → provisional pool `DR_0020` (singleton) → re-homed as a leftover into modelled pool `DR_0001` (AUROC 0.935). **The wrong organism label did not stop it entering a model**, because downstream logic keys off `target_chembl_id_curated`, not the organism field.

**2 — retained and modelled, `CHEMBL4195109_MIC50_uM` (abaumannii).** 10 molecules, all active (degenerate) → original member of modelled pool `DR_0002` (AUROC 0.948). Illustrates that a degenerate (all-one-class) dataset is kept, flagged, and still contributes to a pool.

**3 — discarded despite huge size, `CHEMBL4296188 / INHIBITION%` (abaumannii).** 21,506 molecules, discarded at curation with reason **"Co-ADD"** — discarding is driven by provenance/target rules, not only by size.

**4 — UNCATEGORIZED, now formally discarded at step 21.** Datasets whose `(activity_type, unit)` is not in the DR/SP inclusion list (e.g. `ED50 / ug.mg-1`) used to be retained-but-never-binarised (a silent dead end). This has been fixed: they are now discarded at curation with reason *"uncategorized ((activity_type, unit) not in DR/SP)"* (580 datasets; 0 retained UNCATEGORIZED remain).

## Cross-cutting observations (population-level, all 113,331 datasets)

*Curation-stage fate* (from `21_curation_summary.csv`, all pathogens; current run):

| Outcome at curation (step 21) | Count |
|---|---|
| discarded — ≤5 molecules | 47,650 |
| discarded — no quantitative data (qualitative-only) | 16,619 |
| discarded — not ORGANISM (single protein) | 12,817 |
| discarded — not ORGANISM (other/discarded target) | 3,064 |
| discarded — uncategorized (activity_type/unit not DR/SP) | 580 |
| discarded — Co-ADD / PubChem-superseded / multi-unit | 29 / 11 / 4 |
| **retained → binarisation + pooling** | **32,557** |

Only ~29% of datasets survive curation; the rest are dominated by the ≤5-molecule and qualitative-only discards. The retained 32,557 proceed to binarisation (step 22) and pooling — see the panel for their downstream fate.

## Severity of the organism mislabel

**Question posed:** does the raw-vs-curated organism mismatch affect (a) which pathogen a dataset is assigned to, or (b) whether it is retained?

**(a) Pathogen assignment — no effect (severity: none).**
Pathogen membership is decided at **step 07**, which runs *before* the step-09 LLM, using only the **raw** ChEMBL fields (`target_organism`/`assay_organism` substring match) plus the manual `assays_to_include|exclude` lists. `assay_organism_curated` does not exist yet and is only *carried as metadata* downstream (step-21 `ENRICH_COLS`) — it is never a filter or assignment key, and each pathogen is a separate output folder that nothing after step 07 re-partitions.
Empirically, of **3,936** raw-vs-curated mismatches, **3,923 (99.7%)** are still justified by the raw field for their folder; the remaining **13** had *empty* raw organism fields and the LLM filled in the **correct** pathogen. **Zero** datasets were pulled toward a wrong pathogen.

**(b) Retention — the organism label is not a causal input.**
Discard-as-"not ORGANISM" is decided by `extra_curation_target_type(target_type, target_type_curated)` (`src/pathogen_utils.py:418`), which reads the **raw ChEMBL target type** and the **LLM target type** — *never* `assay_organism_curated`. The rules are raw-anchored: raw `ORGANISM` is always kept as ORGANISM (the LLM cannot demote a whole-organism assay), and raw `SINGLE PROTEIN/PROTEIN*` is always SINGLE PROTEIN. The LLM's call only changes the outcome when the raw type is `UNCHECKED`/`NON-MOLECULAR`.

So the organism mislabel cannot itself cause a discard. It is a **symptom** of the same step-09 pass that also emits `target_type_curated`; harm, if any, lives in that separate field. Narrowing to where it could actually matter:

| Filter (cumulative) | Count |
|---|---|
| Datasets with a raw-vs-curated organism mismatch (raw present) | 3,887 |
| …that are a *genuine* organism error (raw names the pathogen, curated drops it; excludes strain-suffix stripping) | 2,246 |
| …and were discarded as "not ORGANISM" | 898 |
| …**and** raw ChEMBL target type was LLM-decidable (`UNCHECKED`/`NON-MOLECULAR`) — the true risk set | **373** |

Of the 2,057 mismatched-and-discarded datasets, **1,004 were already `SINGLE PROTEIN` in raw ChEMBL** (correct discard, organism irrelevant) and most of the rest were raw `NUCLEIC-ACID / SUBCELLULAR / PROTEIN COMPLEX / ADMET` (also non-organism by raw). The realistic upper bound on datasets whose retention hinged on the step-09 LLM call *and* that carry a genuine organism error is **373 (0.3% of 113,331)** — and only a subset of those are true errors, since `UNCHECKED` assays are frequently legitimately non-organism.

### What the 373 actually are — and whether any are real errors
We inspected all 373. They are **not** whole-cell antibacterial screens wrongly demoted; they are genuinely single-target assays:

- **Activity types** are entirely target/binding readouts — IC50 (153), INHIBITION (77), ACTIVITY (57), KI (49), KD (25), EC50 (8). **Zero MIC-family** (no whole-cell potency).
- **Target type:** the LLM called 370/373 `SINGLE PROTEIN`; all 373 had raw type `UNCHECKED`.
- **The organism relabel usually names the real assay/expression host, not an error.** 284 were relabelled to *E. coli* — heterologous-expression assays where the pathogen's *protein* is expressed in a surrogate host, e.g. *E. faecium* ClpP in *E. coli* BL21(DE3); *K. pneumoniae* NDM1 β-lactamase in *E. coli* SNO3; *C. albicans* CDR1/MDR1 transporter in *S. cerevisiae*. 17 relabelled to *Homo sapiens* are toxin/immunomodulation assays in human cells (*E. coli* Stx2 in HeLa; *E. coli* LPS-stimulated IL6 in human blood). In these, `SINGLE PROTEIN` + exclusion from whole-organism pools is the **correct** call, and the "E. coli"/"S. cerevisiae"/"Homo sapiens" label is arguably *more* accurate than ChEMBL's raw pathogen label.
- **A description-based single-target heuristic** fires for **361/373 (97%)**.
- **No identifiable real error.** A whole-cell/growth-phenotype keyword scan surfaced only 3 candidates, and all three remain *target-directed* (inhibition of EthR / MetAP1c / lipase in infection or surrogate-growth models) rather than general antibacterial screens. The remaining ambiguous handful (cytochrome-bd oxidase in *M. bovis* BCG; fbpD/desA2 RNA-seq at 1×MIC; PPIs such as myoA–MTIP) are also single-target or mechanism-of-action, not activity datasets.

**Verdict:** the organism mislabel is **cosmetic for both pathogen assignment and retention**, and even its highest-risk slice (the 373) resolves to **~0 real data loss** — where the LLM had discretion (raw `UNCHECKED`) it classified single-target correctly. The field remains a useful *canary* for step-09 quality (≈2,200 genuine organism errors exist) but changes no downstream decision on its own.

These are observations, not conclusions. Open questions for you to decide:
1. The only items that are a *definitional* judgment (not a bug): the ~3 intracellular **target-modulation** assays (EthR / MetAP1c / lipase in infection models) — do target-directed infection-model assays belong in whole-organism pools or not?
2. Should `assay_organism_curated` be corrected/deprecated given it drives nothing (canary value vs. misleading-metadata risk)?
3. Are the Co-ADD / PubChem-supersession discards behaving as intended on large assays like example 3?

## How to explore
Open `output/stage4/general_plots/28_dataset_trace_panel.html` (or the published Artifact) and filter by pathogen, fate, or activity type; toggle "anomalies only" to isolate mismatches/degenerate/multi-unit. The panel embeds all 150 sampled traces and the population aggregates (fate/discard counts and the organism-severity funnel, including the 373 true-risk count); re-run `scripts/28_trace_sample_datasets.py` to regenerate it.
