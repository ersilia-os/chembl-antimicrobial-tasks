# Redundant-leftover label-conflict diagnostics

**Date:** 2026-07-15  **Tool:** [`tools/redundant_conflict_audit.py`](../tools/redundant_conflict_audit.py) (read-only)

## Question

Step 25 (`25_merge_leftovers.py`) drops a leftover dataset as **`redundant`** when *every* one of
its molecules is already present in a modelled pool (`n_residual == 0`). This is treated as a
no-loss outcome. But it is only lossless if the leftover's own binary label **agrees** with the
label those molecules already carry in the pool. This audit measures the disagreement three ways,
and separately checks whether step 25's pool-growth policy itself alters shipped labels.

Population: **2,455** redundant leftovers across the 11 pathogens with modelled pools
(ecoli 765, saureus 622, mtuberculosis 340, pfalciparum 151, paeruginosa 148, calbicans 136,
abaumannii 116, spneumoniae 62, kpneumoniae 55, efaecium 38, enterobacter 3).

The tool has three sub-analyses; run `python tools/redundant_conflict_audit.py [all|conflicts|merge-flip|reaverage] [codes...]`.

---

## 1. `conflicts` — raw label comparison

For each redundant leftover, compare its per-assay bin against the bin its molecules carry in the
pools. Three per-molecule definitions (pools are non-disjoint, so a molecule can sit in several):

- **consensus** (headline): leftover ≠ majority of the pools' bins (**tie → active**, matching
  `merge_datasets`).
- **strict**: leftover disagrees with *every* containing pool ("represented nowhere").
- **any**: leftover disagrees with ≥1 pool.

**Result** (29,811 matched molecule-instances):

| definition | conflicts | rate |
|---|---|---|
| consensus | 11,459 | **38.4 %** |
| strict | 8,264 | 27.7 % |
| any | 16,230 | 54.4 % |

- 87.1 % of redundant leftovers have ≥1 consensus conflict; 4,666 are "actives dropped"
  (leftover active, pool consensus inactive).

**Driver (by activity_type).** The conflict is dominated by **numeric MIC-vs-MIC** disagreement,
*not* qualitative assays:

| activity_type | cat | matched mols | consensus conflict |
|---|---|---|---|
| MIC | DR | 21,269 (71 %) | 41.4 % |
| IC50 | DR | 2,226 | 19.8 % |
| MIC90 | DR | 1,797 | 32.3 % |
| INHIBITION | SP | 1,006 | 49.8 % |
| ACTIVITY (qualitative) | DR+SP | 332 | ~22–62 % |

So the disagreement is the same compound getting MIC values across screens that straddle the 10 µM
DR cutoff — cutoff sensitivity, not a qualitative-vs-quantitative artifact.

Outputs: `2026-07-15_redundant-conflicts_{summary,detail,by_type}.csv`,
`2026-07-15_redundant-conflicts{,_by_type}.png`.

---

## 2. `merge-flip` — merge + re-binarise (pipeline-faithful)

The raw comparison overstates loss: a pool molecule usually carries many readings, so one extra
assay rarely swings the averaged, re-binarised label. This analysis merges each redundant leftover
into **every** pool containing its molecules and re-binarises through `merge_datasets`. Three
quantities per molecule × pool (66,902 pairs):

| quantity | meaning | result |
|---|---|---|
| **flip_vs_delivered** | would the shipped label change | 11,584 (**17.3 %**) |
| **flip_marginal** | the leftover's *pure* effect (member set fixed) | 6,466 (9.7 %) |
| **delivered_ne_remerge** | step-25 artifact (see §3) | 5,388 (8.1 %) |

- Flips are overwhelmingly **→ inactive** (10,408 vs 1,176 → active).
- 84.7 % of redundant leftovers flip ≥1 shipped label somewhere.
- So of the 38 % raw conflict, only ~17 % would actually change a shipped label; ~10 % is the
  dropped assay's genuine marginal effect.

Outputs: `2026-07-15_redundant-mergeflip_{summary,detail}.csv`, `…_redundant-mergeflip.png`.

### The step-25 no-re-average artifact

Step 25 grows a pool by appending an assigned leftover's **residual** molecules only; for a molecule
already in the pool it keeps the base value and **discards** the leftover's reading
(`pd.concat([base, add]).drop_duplicates("inchikey")`, `25_merge_leftovers.py`). Yet the leftover is
recorded as a full member in `25_pool_members.csv`. So a clean re-merge of the listed members can
give a different label than the one shipped — the `delivered_ne_remerge` term.

Worked example (abaumannii, pool `DR_0001`, molecule `WZPBZJONDBGPKJ-VEHQQRBSSA-N`):

| member | MIC | role |
|---|---|---|
| `CHEMBL1638694_MIC_uM` | 4.593 µM (`<`) | base member (step 23) |
| `CHEMBL1028858_MIC_uM` | 18.372 µM | merged leftover |

Delivered: value 4.593 µM → **active** (only base kept). Clean re-merge: mean 11.483 µM →
**inactive**. The leftover's reading was silently dropped.

---

## 3. `reaverage` — pool-wide impact of re-averaging

`delivered_ne_remerge` above is counted per leftover × pool, so a promiscuous molecule is counted
many times. This analysis de-duplicates: over **every** molecule in **every** pool, re-merge all
listed members and count how many shipped labels would flip if step 25 re-averaged (molecule set
held fixed; cross-pool molecules a naive full re-merge would add are reported as `n_would_add`, not
flips).

**Result:** **267 / 697,471 molecule × pool = 0.04 %** (232 distinct molecules); ≤ 0.3 % per
pathogen; direction again mostly → inactive.

**Reconciliation with §2's 8.1 %.** The two are consistent — the difference is multiplicity. For
abaumannii DR, step 2's 337 `(leftover, pool, molecule)` instances collapse to **16** unique
`(pool, molecule)` differences; each differing molecule is promiscuous (~42 redundant leftovers
touch it). So the step-25 no-re-average policy is, pool-wide, **effectively inert** (0.04 %).

Outputs: `2026-07-15_redundant-reaverage_{summary,detail}.csv`, `…_redundant-reaverage.png`.

---

## Bottom line

- Raw disagreement between redundant leftovers and pools is high (38 %), driven by MIC-vs-MIC
  cutoff sensitivity.
- Most of it is absorbed by merge averaging: only ~17 % would change a shipped label, ~10 % is the
  dropped assay's marginal effect — almost always pushing a label **toward inactive**.
- The step-25 "base wins / never re-average" policy changes essentially nothing pool-wide (0.04 %).

Open decisions are tracked in [`items-to-improve.md`](items-to-improve.md). No pipeline behaviour
was changed by this audit.
