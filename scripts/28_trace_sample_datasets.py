"""
28_trace_sample_datasets.py  (Stage 4, cross-pathogen)

Read-only AUDIT utility. Samples a small, deliberately-diverse set of datasets
per pathogen and traces each one's full journey through pipeline stages 1, 2 and
4 (stage 3, the parallel conservative branch, is intentionally excluded), from
the raw ChEMBL record through curation, binarisation and pooling, down to where
it finally ends up (modelled pool / merged-leftover / catch-all / dropped /
discarded).

It NEVER re-runs any model and NEVER mutates data: it only reads the summary
files already produced by steps 00/07/08/09 and 21-27 and joins them.

A "dataset" here is the modelling unit: one (assay_id, activity_type, unit)
triplet, whose id is  f"{assay_id}_{activity_type}_{UNIT_SLUG[unit]}"  (the same
`dataset_id` used by the pooling scripts, e.g. CHEMBL5163887_MIC_uM). This maps
1:1 onto the rows of 21_curation_summary.csv (retained datasets carry exactly one
unit; multi-unit pairs are discarded at curation step 3).

Sampling is STRATIFIED and seeded (RANDOM_SEED): a coverage-greedy picker chooses
datasets that jointly span fate, activity-type family, dataset size and
active/inactive ratio, and preferentially surfaces anomalies (raw-vs-curated
organism mismatch, degenerate binarisation, uncategorised, multi-unit). The
bucket edges below are ILLUSTRATIVE sampling strata used only to pick examples —
they are NOT scientific cutoffs applied to the data, and can be edited freely.

Inputs  (per pathogen <p>, all read-only):
  output/stage1/00_assay_descriptions.csv          raw assay description text
  data/chembl_activities/assays.csv                raw ChEMBL assay row
  output/stage2/<p>/07_assays_raw.csv              pre-clean per-assay summary
  output/stage2/<p>/08_assays_cleaned.csv          post-clean per-assay summary
  output/stage2/<p>/09_assays_parameters_full.csv  LLM curation (raw vs curated)
  output/stage4/<p>/21_curation_summary.csv        one row per dataset (incl. discarded)
  output/stage4/<p>/22_binarisation_summary.csv    binary labels (DR/SP only)
  output/stage4/<p>/23_pool_members.csv            provisional step-23 pool membership
  output/stage4/<p>/24_cv_summary.csv              per-pool CV / modelled flag
  output/stage4/<p>/25_pool_members.csv            FINAL pool membership (post leftover-merge)
  output/stage4/<p>/25_pool_summary.csv            grown pool metrics
  output/stage4/<p>/25_merge_log.csv               leftover -> outcome / assigned_pool
  output/stage4/<p>/26_cv_summary.csv              catch-all for fully-deferred categories

Output (output/stage4/general_plots/):
  28_dataset_trace_panel.html    self-contained interactive panel (also publishable as an Artifact).
                                 The 150 sampled traces and the population-level aggregates (fate/discard
                                 counts + the organism-mislabel severity funnel) are computed in-memory and
                                 embedded in the panel; no intermediate CSVs are written.
"""

import os
import sys
import json
import random

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (  # noqa: E402
    RANDOM_SEED,
    STAGE1PATH,
    STAGE2PATH,
    STAGE4PATH,
    DATAPATH,
    CONFIGPATH,
)
from dataset_utils import UNIT_SLUG  # noqa: E402

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
general_dir = os.path.join(STAGE4PATH, "general_plots")
os.makedirs(general_dir, exist_ok=True)

PANEL_HTML = os.path.join(general_dir, "28_dataset_trace_panel.html")

# raw ChEMBL target types where the LLM's curated call decides ORGANISM vs not
# (for every other raw type, extra_curation_target_type is raw-anchored and the
# LLM cannot change the outcome -- see src/pathogen_utils.py:extra_curation_target_type).
LLM_DECIDABLE_RAW_TYPES = {"UNCHECKED", "NON-MOLECULAR"}

# Description signatures of a genuinely single-target / non-whole-organism assay
# (used ONLY to annotate the risk CSV; not a pipeline decision).
SINGLE_TARGET_HINT = (
    r"expressed in|transfected|recombinant|purified|BL21|SNO3|toxin|beta.?lactam|"
    r"transporter|enzyme|protease|kinase|gyrase|reductase|synthase|dehydrogenase|"
    r"lactamase|oxidase|interaction|binding affinity|displacement|lysate|"
    r"gene expression|mutant"
)

# --------------------------------------------------------------------------
# Sampling configuration -- illustrative strata ONLY (safe to edit; not applied
# to the data, only used to choose which examples to show).
# --------------------------------------------------------------------------
N_PER_PATHOGEN = 10

# Size buckets on n_mol_dataset (number of modelled molecules).
SIZE_EDGES = [5, 20, 100]                       # -> "<=5", "6-20", "21-100", ">100"
SIZE_LABELS = ["<=5", "6-20", "21-100", ">100"]

# Active/inactive ratio buckets (binarised datasets only).
RATIO_EDGES = [0.10, 0.50]                       # -> "<0.1", "0.1-0.5", ">0.5"
RATIO_LABELS = ["<0.1", "0.1-0.5", ">0.5"]

# Activity-type families (coarse grouping for diversity; unknown -> "other").
TYPE_FAMILY = {
    "MIC": "potency(MIC)", "MIC50": "potency(MIC)", "MIC80": "potency(MIC)",
    "MIC90": "potency(MIC)", "MBC": "potency(MIC)",
    "IC50": "potency(IC/EC/K)", "IC90": "potency(IC/EC/K)", "EC50": "potency(IC/EC/K)",
    "EC90": "potency(IC/EC/K)", "GI": "potency(IC/EC/K)", "GI50": "potency(IC/EC/K)",
    "KI": "potency(IC/EC/K)", "KD": "potency(IC/EC/K)",
    "INHIBITION": "percent", "SURVIVAL": "percent", "DRUGUPTAKE": "percent",
    "IZ": "zone", "GROWTH": "percent",
    "ACTIVITY": "qualitative",
}


def type_family(activity_type):
    return TYPE_FAMILY.get(str(activity_type).upper(), "other")


def size_bucket(n):
    if pd.isna(n):
        return "n/a"
    n = float(n)
    for edge, label in zip(SIZE_EDGES, SIZE_LABELS):
        if n <= edge:
            return label
    return SIZE_LABELS[-1]


def ratio_bucket(r):
    if pd.isna(r):
        return "n/a"
    r = float(r)
    for edge, label in zip(RATIO_EDGES, RATIO_LABELS):
        if r < edge:
            return label
    return RATIO_LABELS[-1]


# --------------------------------------------------------------------------
# Small IO helpers (all read-only, tolerate missing files)
# --------------------------------------------------------------------------
def _read(path, **kw):
    if os.path.exists(path):
        try:
            return pd.read_csv(path, **kw)
        except Exception as exc:  # pragma: no cover - defensive
            print(f"  ! could not read {path}: {exc}")
    return None


def _pathogens():
    df = pd.read_csv(os.path.join(CONFIGPATH, "pathogens.csv"))
    return list(zip(df["pathogen"], df["code"]))


def dataset_id_of(assay_id, activity_type, unit):
    slug = UNIT_SLUG.get(unit, str(unit))
    return f"{assay_id}_{str(activity_type).replace('/', '-')}_{slug}"


# --------------------------------------------------------------------------
# Load everything for one pathogen and build a per-dataset trace table.
# --------------------------------------------------------------------------
def build_pathogen_trace(pathogen, code):
    s2 = os.path.join(STAGE2PATH, code)
    s4 = os.path.join(STAGE4PATH, code)

    curation = _read(os.path.join(s4, "21_curation_summary.csv"))
    if curation is None or curation.empty:
        print(f"  [{code}] no 21_curation_summary.csv -- skipped")
        return pd.DataFrame()

    binr = _read(os.path.join(s4, "22_binarisation_summary.csv"))
    params = _read(os.path.join(s2, "09_assays_parameters_full.csv"))
    raw07 = _read(os.path.join(s2, "07_assays_raw.csv"))
    clean08 = _read(os.path.join(s2, "08_assays_cleaned.csv"))
    desc = _read(os.path.join(STAGE1PATH, "00_assay_descriptions.csv"))

    m24 = _read(os.path.join(s4, "23_pool_members.csv"))
    m26 = _read(os.path.join(s4, "25_pool_members.csv"))
    cv25 = _read(os.path.join(s4, "24_cv_summary.csv"))
    sum26 = _read(os.path.join(s4, "25_pool_summary.csv"))
    mlog = _read(os.path.join(s4, "25_merge_log.csv"))
    cv27 = _read(os.path.join(s4, "26_cv_summary.csv"))

    # ---- lookups --------------------------------------------------------
    map24 = dict(zip(m24["dataset_id"], m24["pool_id"])) if m24 is not None else {}
    map26 = dict(zip(m26["dataset_id"], m26["pool_id"])) if m26 is not None else {}

    merge_by_leftover = {}
    if mlog is not None:
        for _, r in mlog.iterrows():
            merge_by_leftover[r["leftover"]] = r

    # pool_id -> auroc / modelled (prefer step-25 grown AUROC, fall back to 24)
    pool_auroc, pool_modelled, pool_size = {}, {}, {}
    if cv25 is not None:
        for _, r in cv25.iterrows():
            pool_auroc[r["pool_id"]] = r.get("auroc")
            pool_modelled[r["pool_id"]] = bool(r.get("modelled"))
            pool_size[r["pool_id"]] = r.get("n_molecules")
    if sum26 is not None:
        for _, r in sum26.iterrows():
            if pd.notna(r.get("grown_auroc")):
                pool_auroc[r["pool_id"]] = r.get("grown_auroc")
            pool_size[r["pool_id"]] = r.get("n_mol_after")

    # categories handled by the step-26 catch-all (fully deferred categories)
    cat27 = {}
    if cv27 is not None:
        for _, r in cv27.iterrows():
            cat27[str(r["category"])] = r

    # per-dataset provenance lookups (keyed on assay_id, activity_type)
    def _index(df):
        if df is None:
            return {}
        return {(str(r["assay_id"]), str(r["activity_type"])): r
                for _, r in df.iterrows()}

    params_ix = _index(params)
    raw07_ix = _index(raw07)
    clean08_ix = _index(clean08)
    binr_ix = _index(binr)
    desc_map = dict(zip(desc["chembl_id"], desc["description"])) if desc is not None else {}

    rows = []
    for _, cur in curation.iterrows():
        assay_id = str(cur["assay_id"])
        atype = str(cur["activity_type"])
        key = (assay_id, atype)
        category = str(cur.get("category"))
        discarded = bool(cur.get("discarded"))
        discard_step = int(cur.get("discard_step", 0) or 0)

        b = binr_ix.get(key)
        unit = (b["unit"] if b is not None else cur.get("unit"))
        did = dataset_id_of(assay_id, atype, unit) if b is not None else ""

        p = params_ix.get(key)
        raw_org = (p["assay_organism"] if p is not None else cur.get("assay_organism"))
        cur_org = (p["assay_organism_curated"] if p is not None
                   else cur.get("assay_organism_curated"))
        organism_mismatch = (
            pd.notna(raw_org) and pd.notna(cur_org)
            and str(raw_org).strip().lower() != str(cur_org).strip().lower()
        )

        r07 = raw07_ix.get(key)
        r08 = clean08_ix.get(key)

        rec = {
            "pathogen": code,
            "pathogen_name": pathogen,
            "assay_id": assay_id,
            "activity_type": atype,
            "unit": unit,
            "dataset_id": did,
            "category": category,
            "direction": cur.get("direction"),
            "doc_chembl_id": cur.get("doc_chembl_id"),
            "bao_label": cur.get("bao_label"),
            "description": desc_map.get(assay_id, ""),
            # stage 1/2 provenance + curation (where mislabelling lives)
            "assay_type": (p["assay_type"] if p is not None else None),
            "assay_confidence_score": (p["assay_confidence_score"] if p is not None else None),
            "assay_organism_raw": raw_org,
            "assay_organism_curated": cur_org,
            "organism_mismatch": organism_mismatch,
            "target_type_raw": cur.get("target_type"),
            "target_type_curated": cur.get("target_type_curated"),
            "target_type_curated_extra": cur.get("target_type_curated_extra"),
            "target_chembl_id_curated": cur.get("target_chembl_id_curated"),
            "target_name_curated": cur.get("target_name_curated"),
            "assay_strain_curated": cur.get("assay_strain_curated"),
            "atcc_id": cur.get("atcc_id"),
            "source_label": cur.get("source_label"),
            # stage-2 counts (pre/post clean)
            "n_activities_raw07": (r07["activities"] if r07 is not None else None),
            "n_activities_clean08": (r08["activities"] if r08 is not None else None),
            # stage-4 curation counts (21)
            "n_activities": cur.get("n_activities"),
            "n_numeric": cur.get("n_numeric"),
            "n_text": cur.get("n_text"),
            "n_mol_all": cur.get("n_mol_all"),
            "n_mol_numeric": cur.get("n_mol_numeric"),
            "n_recovered_inactive": cur.get("n_recovered_inactive"),
            "n_text_only_active_discarded": cur.get("n_text_only_active_discarded"),
            "n_mol_dataset": cur.get("n_mol_dataset"),
            "flag_multi_unit": bool(cur.get("flag_multi_unit")),
            "flag_uncategorized": bool(cur.get("flag_uncategorized")),
            # curation decision
            "discarded": discarded,
            "discard_step": discard_step,
            "discard_reason": cur.get("discard_reason"),
        }

        # binarisation (22)
        if b is not None:
            rec.update({
                "bin_cutoff": b.get("cutoff"),
                "bin_n_molecules": b.get("n_molecules"),
                "bin_n_active": b.get("n_active"),
                "bin_n_inactive": b.get("n_inactive"),
                "bin_n_qual_inactive": b.get("n_qual_inactive"),
                "bin_active_ratio": b.get("active_ratio"),
                "bin_degenerate": bool(b.get("degenerate")),
            })
        else:
            rec.update({k: None for k in (
                "bin_cutoff", "bin_n_molecules", "bin_n_active", "bin_n_inactive",
                "bin_n_qual_inactive", "bin_active_ratio")})
            rec["bin_degenerate"] = None

        # pooling / fate resolution
        rec.update(resolve_fate(rec, did, category, discarded, discard_step, b is not None,
                                 map24, map26, merge_by_leftover, pool_modelled,
                                 pool_auroc, pool_size, cat27))
        rows.append(rec)

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Fate resolution -- the definitive "where does it end up" trace.
# Branches per category: a category is either modelled (24/25/26) or fully
# deferred to the step-26 catch-all.
# --------------------------------------------------------------------------
def resolve_fate(rec, did, category, discarded, discard_step, binarised,
                 map24, map26, merge_by_leftover, pool_modelled, pool_auroc,
                 pool_size, cat27):
    out = {
        "pool24": map24.get(did, ""),
        "pool_final": "",
        "pool_auroc": None,
        "pool_n_molecules": None,
        "leftover_outcome": "",
        "fate_class": "",
        "final_endpoint": "",
    }

    if discarded:
        out["fate_class"] = f"discarded_step{discard_step}"
        out["final_endpoint"] = f"Discarded @ step {discard_step}: {rec.get('discard_reason')}"
        return out

    if not binarised:
        out["fate_class"] = "uncategorized_unbinarised"
        out["final_endpoint"] = f"Retained but not binarised (category {category})"
        return out

    # deferred category -> step-26 catch-all
    if category in cat27:
        c = cat27[category]
        kept = bool(c.get("kept"))
        out["pool_final"] = f"{category}_catchall"
        out["pool_auroc"] = c.get("auroc")
        out["pool_n_molecules"] = c.get("n_molecules")
        out["fate_class"] = "catchall_kept" if kept else "catchall_dropped"
        state = "kept" if kept else f"dropped ({c.get('reason')})"
        out["final_endpoint"] = (
            f"Catch-all pool {category}_catchall {state}"
            + (f", AUROC {float(c.get('auroc')):.3f}" if pd.notna(c.get('auroc')) else "")
        )
        return out

    # modelled category -> 24/25/26 path
    lo = merge_by_leftover.get(did)
    if lo is not None:
        outcome = str(lo.get("outcome"))
        out["leftover_outcome"] = outcome
        if outcome == "merged":
            pool = lo.get("assigned_pool")
            out["pool_final"] = pool
            out["pool_auroc"] = pool_auroc.get(pool)
            out["pool_n_molecules"] = pool_size.get(pool)
            out["fate_class"] = "leftover_merged"
            au = f", AUROC {float(pool_auroc.get(pool)):.3f}" if pd.notna(pool_auroc.get(pool)) else ""
            out["final_endpoint"] = f"Merged as leftover into modelled pool {pool}{au}"
        else:
            out["fate_class"] = "leftover_dropped"
            out["final_endpoint"] = f"Leftover not merged ({outcome})"
        return out

    if did in map26:
        pool = map26[did]
        out["pool_final"] = pool
        out["pool_auroc"] = pool_auroc.get(pool)
        out["pool_n_molecules"] = pool_size.get(pool)
        out["fate_class"] = "modelled_pool"
        au = f", AUROC {float(pool_auroc.get(pool)):.3f}" if pd.notna(pool_auroc.get(pool)) else ""
        out["final_endpoint"] = f"Original member of modelled pool {pool}{au}"
        return out

    if did in map24:
        pool = map24[did]
        if pool_modelled.get(pool):
            out["pool_final"] = pool
            out["pool_auroc"] = pool_auroc.get(pool)
            out["pool_n_molecules"] = pool_size.get(pool)
            out["fate_class"] = "modelled_pool"
            au = f", AUROC {float(pool_auroc.get(pool)):.3f}" if pd.notna(pool_auroc.get(pool)) else ""
            out["final_endpoint"] = f"Member of modelled pool {pool}{au}"
        else:
            out["fate_class"] = "leftover_dropped"
            out["final_endpoint"] = f"In non-modelled pool {pool}; not re-homed"
        return out

    out["fate_class"] = "leftover_dropped"
    out["final_endpoint"] = "Binarised but never pooled (ineligible / not re-homed)"
    return out


# --------------------------------------------------------------------------
# Stratified, seeded, coverage-greedy sampler.
# --------------------------------------------------------------------------
def _cells(row):
    """The diversity cells a candidate covers, across the sampling axes."""
    cells = {
        ("fate", row["fate_class"]),
        ("type", type_family(row["activity_type"])),
        ("size", size_bucket(row["n_mol_dataset"])),
        ("ratio", ratio_bucket(row["bin_active_ratio"])),
        ("cat", row["category"]),
    }
    if bool(row.get("organism_mismatch")):
        cells.add(("anom", "organism_mismatch"))
    if row.get("bin_degenerate") is True:
        cells.add(("anom", "degenerate"))
    if bool(row.get("flag_multi_unit")):
        cells.add(("anom", "multi_unit"))
    return cells


def sample_pathogen(df, seed):
    if df.empty:
        return df
    rng = random.Random(seed)
    candidates = df.to_dict("records")
    n = min(N_PER_PATHOGEN, len(candidates))

    covered = set()
    chosen_idx = []
    remaining = list(range(len(candidates)))

    def priority(i):
        # anomalies and discarded first when tied, then deterministic jitter
        r = candidates[i]
        anom = bool(r.get("organism_mismatch")) or (r.get("bin_degenerate") is True)
        disc = str(r["fate_class"]).startswith("discarded")
        return (anom, disc, rng.random())

    while len(chosen_idx) < n and remaining:
        best, best_gain, best_pri = None, -1, None
        for i in remaining:
            gain = len(_cells(candidates[i]) - covered)
            pri = priority(i)
            if gain > best_gain or (gain == best_gain and pri > (best_pri or (False, False, -1))):
                best, best_gain, best_pri = i, gain, pri
        chosen_idx.append(best)
        covered |= _cells(candidates[best])
        remaining.remove(best)
        if best_gain == 0:
            # no more diversity to gain: fill the rest with a seeded shuffle
            rng.shuffle(remaining)
            for i in remaining[: n - len(chosen_idx)]:
                chosen_idx.append(i)
            break

    return df.iloc[chosen_idx].reset_index(drop=True)


# --------------------------------------------------------------------------
# HTML panel (self-contained fragment: no <html>/<head>/<body> wrappers so it
# both opens in a browser and publishes directly as a Claude Artifact).
# --------------------------------------------------------------------------
def build_panel(sample_df, pop_stats):
    records = json.loads(sample_df.to_json(orient="records"))
    payload = json.dumps({"records": records, "stats": pop_stats}, ensure_ascii=False)
    # data is inlined; the template below reads window.__TRACE__
    template = _PANEL_TEMPLATE.replace("/*__DATA__*/", payload)
    return template


_PANEL_TEMPLATE = r"""
<style>
  :root{
    --bg:#ffffff;--fg:#1a1a2e;--muted:#6b7280;--card:#f7f7fb;--border:#e3e3ee;
    --crimson:#d1495b;--silver:#adb5bd;--cobalt:#3d5a80;--tangerine:#ee9b00;
    --amber:#e07a1f;--green:#2a9d8f;--chip:#ececf5;
  }
  @media (prefers-color-scheme: dark){
    :root{--bg:#14141c;--fg:#e8e8f0;--muted:#9aa0ad;--card:#1d1d29;--border:#2c2c3a;--chip:#25253a;}
  }
  :root[data-theme="dark"]{--bg:#14141c;--fg:#e8e8f0;--muted:#9aa0ad;--card:#1d1d29;--border:#2c2c3a;--chip:#25253a;}
  :root[data-theme="light"]{--bg:#ffffff;--fg:#1a1a2e;--muted:#6b7280;--card:#f7f7fb;--border:#e3e3ee;--chip:#ececf5;}
  *{box-sizing:border-box}
  .wrap{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;color:var(--fg);background:var(--bg);padding:18px;max-width:1200px;margin:0 auto;}
  h1{font-size:20px;margin:0 0 2px}
  .sub{color:var(--muted);font-size:13px;margin:0 0 14px}
  .bar{display:flex;flex-wrap:wrap;gap:14px;margin:10px 0 16px;padding:12px;background:var(--card);border:1px solid var(--border);border-radius:10px}
  .stat{font-size:12px;color:var(--muted)}
  .stat b{display:block;font-size:18px;color:var(--fg)}
  .controls{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:14px;align-items:center}
  select,input[type=search]{padding:6px 8px;border:1px solid var(--border);border-radius:8px;background:var(--bg);color:var(--fg);font-size:13px}
  label.tgl{font-size:12px;color:var(--muted);display:flex;gap:5px;align-items:center}
  .card{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:12px 14px;margin-bottom:10px}
  .card h3{margin:0 0 3px;font-size:14px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
  .meta{font-size:12px;color:var(--muted);margin-bottom:9px}
  .desc{font-size:12px;margin:6px 0 10px;line-height:1.4}
  .journey{display:flex;flex-wrap:wrap;gap:6px;overflow-x:auto}
  .node{flex:1 1 110px;min-width:110px;border:1px solid var(--border);border-radius:8px;padding:7px 9px;background:var(--bg)}
  .node .t{font-size:10px;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}
  .node .v{font-size:12px;margin-top:3px;word-break:break-word}
  .node.term{border-width:2px}
  .badges{margin-top:8px;display:flex;flex-wrap:wrap;gap:5px}
  .chip{font-size:11px;padding:2px 8px;border-radius:20px;background:var(--chip);color:var(--fg)}
  .chip.warn{background:var(--amber);color:#fff}
  .chip.bad{background:var(--crimson);color:#fff}
  .fate{font-weight:600}
  .fate.discarded_step1,.fate.discarded_step2,.fate.discarded_step3,
  .fate.discarded_step4,.fate.discarded_step5,.fate.discarded_step6,
  .fate.discarded_step7,.fate.discarded_step8,.fate.discarded_step9,
  .fate.leftover_dropped,.fate.catchall_dropped{color:var(--amber)}
  .fate.modelled_pool,.fate.leftover_merged,.fate.catchall_kept{color:var(--green)}
  .fate.uncategorized_unbinarised{color:var(--muted)}
  .empty{color:var(--muted);font-size:13px;padding:20px;text-align:center}
</style>
<div class="wrap">
  <h1>Dataset flow audit &mdash; 10 datasets &times; pathogen</h1>
  <p class="sub">Each dataset traced through pipeline stages 1, 2 and 4. Read-only; sampled with a fixed seed.</p>
  <div class="bar" id="statbar"></div>
  <div class="controls">
    <select id="fPath"></select>
    <select id="fFate"></select>
    <select id="fType"></select>
    <label class="tgl"><input type="checkbox" id="fAnom"> anomalies only</label>
    <input type="search" id="fText" placeholder="search assay id / text&hellip;">
  </div>
  <div id="cards"></div>
</div>
<script>
const DATA = /*__DATA__*/;
const R = DATA.records;
const fmt = (v)=> (v===null||v===undefined||v==="")?"&mdash;":v;
const nz = (v)=> (v===null||v===undefined||Number.isNaN(v))?"&mdash;":v;

function statbar(){
  const s=DATA.stats;
  const el=document.getElementById('statbar');
  const items=[["datasets sampled",s.n_sampled],["pathogens",s.n_pathogens],
    ["organism mismatches (all)",s.n_organism_mismatch],
    ["organism true-risk (all)",s.org_true_risk_llm_decidable],
    ["discarded (sampled)",s.n_discarded_sampled],["reach a modelled pool (sampled)",s.n_modelled_sampled]];
  el.innerHTML=items.map(([k,v])=>`<div class="stat"><b>${v}</b>${k}</div>`).join("");
}

function opts(id,vals,label){
  const el=document.getElementById(id);
  el.innerHTML=`<option value="">all ${label}</option>`+vals.map(v=>`<option>${v}</option>`).join("");
}

function node(t,v,term){return `<div class="node ${term?'term':''}"><div class="t">${t}</div><div class="v">${v}</div></div>`;}

function card(r){
  const badges=[];
  if(r.organism_mismatch) badges.push(`<span class="chip bad">organism mismatch</span>`);
  if(r.bin_degenerate===true) badges.push(`<span class="chip warn">degenerate</span>`);
  if(r.flag_multi_unit) badges.push(`<span class="chip warn">multi-unit</span>`);
  if(r.flag_uncategorized) badges.push(`<span class="chip">uncategorized</span>`);
  badges.push(`<span class="chip">${fmt(r.category)}</span>`);

  const bin = (r.bin_n_active===null||r.bin_n_active===undefined)?"skipped"
      :`${r.bin_n_active} act / ${r.bin_n_inactive} inact`;
  const j=[
    node("ChEMBL raw org", fmt(r.assay_organism_raw)),
    node("Stage-2 clean (08)", nz(r.n_activities_clean08)+" acts"),
    node("Curate org (09)", fmt(r.assay_organism_curated)),
    node("Curate (21)", nz(r.n_mol_dataset)+" mol"+(r.discarded?" &rarr; DISCARDED":"")),
    node("Binarise (22)", bin+(r.bin_cutoff?` @${r.bin_cutoff}`:"")),
    node("Pool (24)", fmt(r.pool24)),
    node("Model/merge", fmt(r.pool_final)+(r.pool_auroc!=null?` (AUROC ${(+r.pool_auroc).toFixed(3)})`:"")),
    node(`Final &mdash; <span class="fate ${r.fate_class}">${r.fate_class}</span>`, r.final_endpoint, true),
  ].join("");

  return `<div class="card">
    <h3>${r.dataset_id||r.assay_id} <span style="color:var(--muted);font-weight:400">&middot; ${r.pathogen}</span></h3>
    <div class="meta">${fmt(r.activity_type)} / ${fmt(r.unit)} &middot; dir ${fmt(r.direction)} &middot; conf ${fmt(r.assay_confidence_score)} &middot; ${fmt(r.doc_chembl_id)}</div>
    <div class="desc">${fmt(r.description)}</div>
    <div class="journey">${j}</div>
    <div class="badges">${badges.join("")}</div>
  </div>`;
}

function render(){
  const p=document.getElementById('fPath').value,
        f=document.getElementById('fFate').value,
        t=document.getElementById('fType').value,
        a=document.getElementById('fAnom').checked,
        q=document.getElementById('fText').value.toLowerCase();
  const rows=R.filter(r=>{
    if(p&&r.pathogen!==p)return false;
    if(f&&r.fate_class!==f)return false;
    if(t&&r.activity_type!==t)return false;
    if(a&&!(r.organism_mismatch||r.bin_degenerate===true||r.flag_multi_unit))return false;
    if(q){const hay=((r.assay_id||"")+" "+(r.dataset_id||"")+" "+(r.description||"")).toLowerCase();
          if(!hay.includes(q))return false;}
    return true;
  });
  const c=document.getElementById('cards');
  c.innerHTML=rows.length?rows.map(card).join(""):`<div class="empty">No datasets match these filters.</div>`;
}

statbar();
opts('fPath',[...new Set(R.map(r=>r.pathogen))].sort(),'pathogens');
opts('fFate',[...new Set(R.map(r=>r.fate_class))].sort(),'fates');
opts('fType',[...new Set(R.map(r=>r.activity_type))].sort(),'activity types');
['fPath','fFate','fType','fAnom','fText'].forEach(id=>{
  const ev=(id==='fText')?'input':'change';
  document.getElementById(id).addEventListener(ev,render);
});
render();
</script>
"""


# --------------------------------------------------------------------------
# Organism-mislabel severity: funnel from "any mismatch" down to the datasets
# whose RETENTION could actually hinge on the step-09 LLM call, and write those
# out for human review. Establishes (empirically) that the curated organism is
# never itself a filter -- retention is decided by target_type, not organism.
# --------------------------------------------------------------------------
def organism_severity(full):
    if full.empty:
        return {}, pd.DataFrame()

    raw = full["assay_organism_raw"].astype(str).str.strip()
    cur = full["assay_organism_curated"].astype(str).str.strip()
    name = full["pathogen_name"].astype(str)
    rawtt = full["target_type_raw"].astype(str).str.upper().str.strip()
    disc_notorg = full["discard_reason"].astype(str).str.contains("not ORGANISM", na=False)

    raw_present = raw.ne("nan") & raw.ne("")
    # substring checks are row-wise because the pattern (pathogen name) varies per row
    raw_is_path = pd.Series([str(n).lower() in str(r).lower() for r, n in zip(raw, name)], index=full.index)
    cur_keeps_path = pd.Series([str(n).lower() in str(c).lower() for c, n in zip(cur, name)], index=full.index)

    mism_rawpresent = raw_present & raw.str.lower().ne(cur.str.lower())
    # a *genuine* organism error excludes strain-suffix stripping (raw "E. coli K-12"
    # -> curated "Escherichia coli" still keeps the pathogen name -> not an error)
    genuine = mism_rawpresent & raw_is_path & ~cur_keeps_path
    llm_decidable = rawtt.isin(LLM_DECIDABLE_RAW_TYPES)
    true_risk = genuine & disc_notorg & llm_decidable

    counts = {
        "org_mismatch_raw_present": int(mism_rawpresent.sum()),
        "org_genuine_error": int(genuine.sum()),
        "org_genuine_error_discarded_notORG": int((genuine & disc_notorg).sum()),
        "org_true_risk_llm_decidable": int(true_risk.sum()),
    }

    risk = full[true_risk].copy()
    risk["single_target_hint"] = (
        risk["description"].astype(str).str.contains(SINGLE_TARGET_HINT, case=False, regex=True)
    )
    cols = ["pathogen", "assay_id", "activity_type", "unit", "bao_label",
            "assay_organism_raw", "assay_organism_curated",
            "target_type_raw", "target_type_curated", "target_type_curated_extra",
            "n_mol_dataset", "discard_reason", "single_target_hint", "description"]
    risk = risk[[c for c in cols if c in risk.columns]]
    return counts, risk


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    pathogens = _pathogens()
    all_traces, samples = [], []

    for idx, (name, code) in enumerate(pathogens):
        print(f"[{code}] tracing...")
        df = build_pathogen_trace(name, code)
        if df.empty:
            continue
        all_traces.append(df)
        sample = sample_pathogen(df, seed=RANDOM_SEED * 1000 + idx)
        samples.append(sample)
        print(f"  {len(df)} datasets -> sampled {len(sample)}")

    full = pd.concat(all_traces, ignore_index=True) if all_traces else pd.DataFrame()
    sample_df = pd.concat(samples, ignore_index=True) if samples else pd.DataFrame()
    print(f"\ntraced {len(full)} datasets; sampled {len(sample_df)} for the panel")

    # population-level summary (over ALL datasets, not just the sample); computed
    # in-memory and embedded in the panel -- no CSVs are persisted.
    def _counts(series):
        return {str(k): int(v) for k, v in series.value_counts(dropna=False).items()}

    org_mask = full["target_type_curated_extra"].astype(str).eq("ORGANISM") if not full.empty else None
    pop = {
        "n_datasets_total": int(len(full)),
        "n_sampled": int(len(sample_df)),
        "n_pathogens": int(sample_df["pathogen"].nunique()) if not sample_df.empty else 0,
        "n_organism_mismatch": int(full["organism_mismatch"].sum()) if not full.empty else 0,
        "n_organism_mismatch_organism_rows":
            int(full.loc[org_mask, "organism_mismatch"].sum()) if org_mask is not None else 0,
        "n_degenerate_sampled": int((sample_df["bin_degenerate"] == True).sum()) if not sample_df.empty else 0,
        "n_discarded_sampled": int(sample_df["fate_class"].astype(str).str.startswith("discarded").sum()) if not sample_df.empty else 0,
        "n_modelled_sampled": int(sample_df["fate_class"].isin(["modelled_pool", "leftover_merged", "catchall_kept"]).sum()) if not sample_df.empty else 0,
        "fate_class_all": _counts(full["fate_class"]) if not full.empty else {},
        "discard_reason_all": _counts(full.loc[full["discarded"], "discard_reason"]) if not full.empty else {},
    }

    # organism-mislabel severity funnel (counts only; embedded in the panel stats)
    org_counts, _risk_df = organism_severity(full)
    pop.update(org_counts)
    print(f"  organism true-risk datasets: {pop.get('org_true_risk_llm_decidable', 0)}")

    panel = build_panel(sample_df, pop)
    with open(PANEL_HTML, "w") as fh:
        fh.write(panel)
    print(f"wrote {PANEL_HTML}")


if __name__ == "__main__":
    main()
