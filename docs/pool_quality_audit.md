# Stage-4 pool quality audit — sub-threshold final pools and their chem-space coverage

**Generated:** 2026-07-09 by `tools/pool_quality_audit.py` (regenerate any time).
**AUROC threshold:** `AUROC_MIN_THRESHOLD = 0.7` (from `src/default.py`). A pool *reaches* the threshold when its shipped CV AUROC ≥ 0.7; it is *sub-threshold* below that.

Across **15 pathogens**, **3** have at least one sub-threshold final pool (5 such pools in total).

## What this measures

Final pools are the shipped datasets from **step 25** (leftover-merged modelled pools) and **step 26** (low-data catch-alls). For each pathogen, sub-threshold pools are counted and their molecule coverage reported in two views:

- **Gross** — union of molecules in the sub-threshold pools, *ignoring* whether the same molecules are also covered by a ≥-threshold pool. *(What passes through a weak model.)*
- **Exclusive** — the subset in **no** ≥-threshold pool anywhere in the pathogen. *(What would lose all modelled coverage if the weak pools were dropped.)*

Each is given against two denominators: the pathogen-**category** chemical space and the **total-pathogen** chemical space (DR ∪ SP unique InChIKeys, computed from `22_binarised/`). Step-25 shipped AUROC = `grown_auroc` if the batch was accepted else `baseline_auroc`. Step-26 catch-alls dropped with reason `low_auroc` count as sub-threshold (their molecules = the category's whole binarised union); catch-alls dropped as `too_small`/`not_cv_able` were never scored and are listed separately, not in the coverage math.

## Per-pathogen summary

`excl %` / `gross %` are of the **total pathogen** chem-space. `modelled source` flags pathogens where step 25 did not run.

| pathogen | #pools <thr | gross mol | gross % path | excl mol | excl % path | pathogen chem-space | modelled source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| pfalciparum | 3 | 1182 | 0.425% | 450 | 0.162% | 277794 | 25 |
| efaecium | 1 | 135 | 2.038% | 85 | 1.283% | 6625 | 25 |
| hpylori | 1 | 103 | 12.827% | 100 | 12.453% | 803 | 25 |
| abaumannii | 0 | 0 | 0.0% | 0 | 0.0% | 7102 | 25 |
| calbicans | 0 | 0 | 0.0% | 0 | 0.0% | 29278 | 25 |
| campylobacter | 0 | 0 | 0.0% | 0 | 0.0% | 96 | none (catch-alls only) |
| ecoli | 0 | 0 | 0.0% | 0 | 0.0% | 52171 | 25 |
| enterobacter | 0 | 0 | 0.0% | 0 | 0.0% | 3962 | 25 |
| kpneumoniae | 0 | 0 | 0.0% | 0 | 0.0% | 16891 | 25 |
| mtuberculosis | 0 | 0 | 0.0% | 0 | 0.0% | 121564 | 25 |
| ngonorrhoeae | 0 | 0 | 0.0% | 0 | 0.0% | 473 | 24 (step 25 not run) |
| paeruginosa | 0 | 0 | 0.0% | 0 | 0.0% | 32892 | 25 |
| saureus | 0 | 0 | 0.0% | 0 | 0.0% | 64470 | 25 |
| smansoni | 0 | 0 | 0.0% | 0 | 0.0% | 2256 | 25 |
| spneumoniae | 0 | 0 | 0.0% | 0 | 0.0% | 10341 | 25 |

## Per-category detail (pathogens with sub-threshold pools)

| pathogen | category | #pools | gross mol | gross % cat | gross % path | excl mol | excl % cat | excl % path |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| efaecium | SP | 1 | 135 | 100.0% | 2.038% | 85 | 62.963% | 1.283% |
| hpylori | SP | 1 | 103 | 100.0% | 12.827% | 100 | 97.087% | 12.453% |
| pfalciparum | DR | 3 | 1182 | 2.847% | 0.425% | 450 | 1.084% | 0.162% |

## Sub-threshold pools

| pathogen | category | pool | source | AUROC | molecules |
| --- | --- | --- | --- | --- | --- |
| efaecium | SP | SP_catchall | step 26 | 0.662 | 135 |
| hpylori | SP | SP_catchall | step 26 | 0.586 | 103 |
| pfalciparum | DR | DR_0008 | step 25 | 0.6664 | 818 |
| pfalciparum | DR | DR_0011 | step 25 | 0.6873 | 759 |
| pfalciparum | DR | DR_0013 | step 25 | 0.6471 | 640 |

## Catch-all categories never modelled (no AUROC — context, not in coverage math)

These categories had a catch-all attempt dropped before CV for insufficient data, so their molecules have no model at all — a different gap from low-quality models.

| pathogen | category | reason | molecules |
| --- | --- | --- | --- |
| campylobacter | SP | too_small | 14 |
| ngonorrhoeae | SP | too_small | 43 |

---

*Sources: `output/stage4/<code>/{25_pool_summary,24_cv_summary,26_cv_summary}.csv`, the pool molecule files under `25_pools/`, `23_pools/`, `26_pools/`, and the category chem-space from `22_binarised/`. Threshold from `src/default.py`. Machine-readable table: `output/stage4/general_plots/pool_quality_audit.csv`.*
