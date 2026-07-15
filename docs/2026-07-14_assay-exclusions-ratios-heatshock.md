# Assay-level exclusions: ratio/selectivity metrics and heat-shock readouts

**Date:** 2026-07-14
**Scope:** `config/assays_to_exclude/<code>.csv` — the manual deny list enforced at steps 16 (merging) and 21 (curation).

## Context

The deny list previously held one entry each for `ecoli` (a PubChem growth-inhibition
**counterscreen**, `CHEMBL1794494`) and `mtuberculosis` (the Cornell **heat-shock** assay,
`CHEMBL4649971`). We swept the descriptions of all ~108,893 curation-feeding assays across the
15 pathogens (joining each `output/stage3/<code>/19_assays_master.csv` to
`output/stage1/00_assay_descriptions.csv`) for the same two flavours: counter-screens and
non-potency "weird" readouts.

## Findings and decisions

### 1. Selectivity / ratio metrics that survive as DR datasets (the real leak)

These are ORGANISM-typed and **were not caught** by any existing gate — they reach the curated
DR datasets even though their measured value is a **ratio**, not an antimicrobial potency.
Excluded (signed off):

| Pathogen | Assay | n mol | Reason |
|---|---|---|---|
| pfalciparum | `CHEMBL761397` | 7 | "Selectivity for FM3A cells over *P. falciparum* as mean EC50 ratio" — host/parasite selectivity ratio; wrong magnitude **and** direction for antiplasmodial potency. |
| pfalciparum | `CHEMBL763544` | 9 | "Mean EC50 FM3A cells / mean EC50 *P. falciparum*" — same selectivity ratio. |

MIC fold-change ratios (`ecoli` CHEMBL999344, CHEMBL1691264; `mtuberculosis` CHEMBL4775336)
and the `calbicans` fluconazole chemo-sensitization assay (CHEMBL4390556) were **reviewed and
left in** — not excluded this round.

### 2. Heat-shock / stress-response readouts (documentary; already discarded)

Seventeen ORGANISM-typed assays measuring thermotolerance, heat-shock-protein expression
(microarray / MALDI-TOF / RT-qPCR) or σ32 heat-shock promoter reporter activity — the same
class as the existing `mtuberculosis` exclusion. **All are already discarded at step 21** (single
-molecule qualitative, or ≤5 molecules), so excluding them changes no output; they are listed to
record intent explicitly.

- `ecoli`: CHEMBL933980, CHEMBL934966, CHEMBL934967, CHEMBL934968, CHEMBL934969, CHEMBL934970, CHEMBL934978
- `paeruginosa`: CHEMBL1052272, CHEMBL1052273, CHEMBL1052974
- `calbicans`: CHEMBL1665535, CHEMBL1665548, CHEMBL1665550, CHEMBL1665552, CHEMBL1665555
- `pfalciparum`: CHEMBL1695391
- `mtuberculosis`: CHEMBL4031583

### 3. Counterscreens

The only ORGANISM-level assay labelled "Counterscreen" in the entire set is the pre-existing
`ecoli` exclusion. No new counterscreens found.

## Not pursued

Broad classes matched but **not** excluded because they are dominated by legitimate potency
assays or accepted proxies: gene-expression/microarray profiling (908 hits, incl. cell-free
transcription-translation inhibition), reporter/promoter readouts (503, incl. pLDH/luciferase
antimalarial readouts), host-cell cytotoxicity (210, incl. intracellular-macrophage Mtb potency).
Excluding these wholesale would be a curation-policy change, not a targeted counter-screen removal.
