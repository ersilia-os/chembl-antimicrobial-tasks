"""
Pooling engine for the fresh ORG dataset-creation strategy (post-step-12 redesign).

This module is STRICTLY ADDITIVE: it only reads existing `output/<code>/` and `config/`
files and never modifies them. It implements one shared engine used by three dataset-shaping
strategies (see scripts/22_build_pooled_org_datasets.py):

    S1  two pooled binary datasets per pathogen (one Dose-Response, one Single-Point)
    S2  endpoint-family binary datasets (a handful)
    S3  adaptive concordance-clustering (data-driven number of binary datasets)

Core ideas
----------
- Work on ORGANISM assays only, numeric-valued records only.
- Split by assay format using `direction`: Dose-Response (direction == -1, potency endpoints)
  vs Single-Point (direction == +1, % / zone readouts).
- Deduplicate molecules by standardized InChIKey; within an assay keep the most-potent row.
- Binarize each (assay, activity_type, unit) endpoint at its reference expert cutoff, then pool
  binary labels across assays (a molecule is active if active under any contributing endpoint).
- Resolve contradictions with a concordance graph + a model-in-the-loop check: an assay that
  disagrees with the reference on shared molecules is dropped only if removing it does not reduce
  the pooled dataset's cross-validated AUROC.

Reused utilities: dataset_utils.{get_cut_value, adjust_relation, disambiguate_compounds,
binarize_with_expert_cutoff}, pathogen_utils.{load_expert_cutoffs, reference_cutoff},
model_utils.{load_ecfp_subset_by_chembl_id, add_decoys, downsample_negatives, KFoldTrain}.
"""

import os
import pickle
import numpy as np
import pandas as pd

from dataset_utils import (
    get_cut_value,
    adjust_relation,
    disambiguate_compounds,
    binarize_with_expert_cutoff,
)
from pathogen_utils import load_expert_cutoffs, reference_cutoff
from model_utils import (
    load_ecfp_subset_by_chembl_id,
    add_decoys,
    downsample_negatives,
    KFoldTrain,
)
from default import DATAPATH, DECOY_RATIO

SEED = 42  # matches the hardcoded seed used across model_utils (no project-wide RANDOM_SEED exists)
NA = "<NA>"

# ---- default engine parameters (overridable from the CLI; flagged in the plan) ----
DEFAULTS = dict(
    min_shared=5,            # K: min shared molecules for a concordance comparison to count
    concordance_threshold=0.70,   # below this two endpoints are "discordant"
    auroc_tolerance=0.02,    # max AUROC drop tolerated when dropping a discordant assay
    min_actives_to_model=10, # need at least this many actives to attempt CV AUROC
    max_model_checks=200,    # cap on in-loop model evaluations per pathogen (cost guard)
    max_datasets_per_bucket=6,  # S3: cap; smallest clusters merged into the residual dataset
    n_estimators_inloop=50,  # lighter RF for the in-loop modellability check
    max_model_rows=20000,    # subsample cap for the in-loop AUROC proxy (speed)
)

ECFP_H5 = os.path.join(DATAPATH, "chembl_processed", "06_chembl_ecfps.h5")


# ===========================================================================
# Loading & endpoint preparation
# ===========================================================================

def _norm(series):
    return series.astype("string").fillna(NA).str.strip()


def _key(assay, act, unit):
    a = assay if isinstance(assay, str) else str(assay)
    u = unit if (isinstance(unit, str) and unit.strip()) else NA
    return (a.strip(), str(act).strip(), u)


def assign_family(activity_type, unit, bucket):
    """Coarse endpoint family for strategy S2."""
    at = str(activity_type).upper()
    if bucket == "DR":
        if at.startswith("MIC"):
            return "DR_growth_MIC"          # MIC, MIC50, MIC80, MIC90
        return "DR_other_potency"           # IC50, EC50, GI50, ED50, KI, KD, ...
    # SP
    if unit == "mm":
        return "SP_zone"                    # inhibition zone diameter
    return "SP_percent"                     # % inhibition / activity / growth


def load_endpoints(pathogen_code, output_root, config_path, cache=True):
    """Build per-(assay, activity_type, unit) binarized ORG endpoints.

    Strategy/reference-independent, so it is cached per pathogen under
    `output/<code>/22_pooled/_cache_endpoints.pkl` (additive) to make repeated runs fast.

    Returns
    -------
    endpoints : dict[key -> DataFrame[inchikey, compound_chembl_id, smiles, value, score, bin]]
    meta : dict[key -> dict(bucket, family, activity_type, unit, n_cpds, has_cutoff)]
    org_inchikeys : set   (all unique InChIKeys assessed in ORG assays, numeric+qual, for coverage)
    """
    pdir = os.path.join(output_root, pathogen_code)
    cache_path = os.path.join(pdir, "22_pooled", "_cache_endpoints.pkl")
    if cache and os.path.exists(cache_path):
        with open(cache_path, "rb") as fh:
            return pickle.load(fh)
    expert = load_expert_cutoffs(config_path)

    # InChIKey + smiles map
    counts = pd.read_csv(os.path.join(pdir, "07_compound_counts.csv.gz"),
                         usecols=["compound_chembl_id", "inchikey", "smiles"])
    cid2ik = dict(zip(counts["compound_chembl_id"], counts["inchikey"]))
    cid2smi = dict(zip(counts["compound_chembl_id"], counts["smiles"]))

    # ORG (assay, activity_type, unit) keys from the step-12 classification (read-only)
    info = pd.read_csv(os.path.join(pdir, "12_assay_data_info.csv"), low_memory=False)
    org = info[info["target_type_curated_extra"] == "ORGANISM"]
    org_keys = {_key(a, act, u) for a, act, u in
                zip(org["assay_id"], org["activity_type"], org["unit"])}

    cleaned = pd.read_csv(
        os.path.join(pdir, "08_chembl_cleaned_data.csv.gz"),
        usecols=["assay_chembl_id", "compound_chembl_id", "activity_type", "unit",
                 "value", "relation", "direction", "pchembl_calculated"],
        low_memory=False,
    )
    cleaned["key"] = [_key(a, act, u) for a, act, u in
                      zip(cleaned["assay_chembl_id"], cleaned["activity_type"], cleaned["unit"])]

    # all ORG InChIKeys (for coverage denominator), before numeric filtering
    org_rows_all = cleaned[cleaned["key"].isin(org_keys)]
    org_inchikeys = {cid2ik.get(c) for c in org_rows_all["compound_chembl_id"].unique()}
    org_inchikeys.discard(None)
    org_inchikeys.discard(np.nan)

    # numeric-valued ORG records with a usable direction
    df = org_rows_all[org_rows_all["value"].notna() & org_rows_all["direction"].isin([-1.0, 1.0])]

    endpoints, meta = {}, {}
    for key, g in df.groupby("key"):
        assay, act, unit = key
        direction = int(g["direction"].iloc[0])
        bucket = "DR" if direction == -1 else "SP"

        # censored-value fix + most-potent measurement per compound (reused utilities)
        cut = get_cut_value(g, direction)
        g2 = adjust_relation(g, direction, cut)
        g2 = disambiguate_compounds(g2, direction)

        # map to InChIKey; if two compound ids collapse to one key, keep the most active
        g2 = g2.copy()
        g2["inchikey"] = g2["compound_chembl_id"].map(cid2ik)
        g2 = g2[g2["inchikey"].notna()]
        if g2.empty:
            continue
        g2 = disambiguate_compounds(g2.rename(columns={"compound_chembl_id": "_cid_tmp",
                                                       "inchikey": "compound_chembl_id"}),
                                    direction).rename(columns={"compound_chembl_id": "inchikey",
                                                               "_cid_tmp": "compound_chembl_id"})

        unit_key = unit if unit != NA else np.nan
        cutoffs = expert.get((act, unit_key, "ORGANISM", pathogen_code))
        refc = reference_cutoff(cutoffs, unit_key) if cutoffs else None
        has_cutoff = refc is not None
        if has_cutoff:
            g2 = binarize_with_expert_cutoff(g2, refc, direction)
        else:
            g2["bin"] = np.nan  # no expert cutoff -> cannot binarize; excluded from binary pools

        g2["smiles"] = g2["compound_chembl_id"].map(cid2smi)
        g2["score"] = g2["pchembl_calculated"] if bucket == "DR" else g2["value"]

        endpoints[key] = g2[["inchikey", "compound_chembl_id", "smiles", "value", "score", "bin"]].reset_index(drop=True)
        meta[key] = dict(bucket=bucket, family=assign_family(act, unit, bucket),
                         activity_type=act, unit=unit, n_cpds=g2["inchikey"].nunique(),
                         has_cutoff=has_cutoff)

    result = (endpoints, meta, org_inchikeys)
    if cache:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, "wb") as fh:
            pickle.dump(result, fh)
    return result


# ===========================================================================
# Concordance graph
# ===========================================================================

def load_overlap_pairs(pathogen_code, output_root):
    """Candidate overlapping assay-endpoint pairs from 11_assays_overlap.csv (read-only)."""
    path = os.path.join(output_root, pathogen_code, "11_assays_overlap.csv")
    if not os.path.exists(path) or os.path.getsize(path) < 5:
        return []
    o = pd.read_csv(path)
    pairs = []
    for r in o.itertuples(index=False):
        k1 = _key(r.assay_id_1, r.activity_type_1, r.unit_1)
        k2 = _key(r.assay_id_2, r.activity_type_2, r.unit_2)
        if k1 != k2:
            pairs.append((k1, k2))
    return pairs


def concordance(ep_a, ep_b, min_shared):
    """Fraction of shared InChIKeys with the same binary label; None if < min_shared or unbinarised."""
    a = ep_a.dropna(subset=["bin"])[["inchikey", "bin"]]
    b = ep_b.dropna(subset=["bin"])[["inchikey", "bin"]]
    m = a.merge(b, on="inchikey", suffixes=("_a", "_b"))
    if len(m) < min_shared:
        return None, len(m)
    return float((m["bin_a"] == m["bin_b"]).mean()), len(m)


def build_concordance(endpoints, meta, pairs, params):
    """Return dict[(k1,k2)] -> (concordance, n_shared) for same-bucket overlapping pairs."""
    out = {}
    seen = set()
    for k1, k2 in pairs:
        if k1 not in endpoints or k2 not in endpoints:
            continue
        if meta[k1]["bucket"] != meta[k2]["bucket"]:
            continue
        pair = tuple(sorted([k1, k2]))
        if pair in seen:
            continue
        seen.add(pair)
        c, n = concordance(endpoints[k1], endpoints[k2], params["min_shared"])
        if c is not None:
            out[pair] = (c, n)
    return out


# ===========================================================================
# Pooling, modellability, conflict pruning
# ===========================================================================

# Features in 06_chembl_ecfps.h5 are keyed by compound_chembl_id, but pooled rows are keyed by
# InChIKey. We carry a representative compound_chembl_id per InChIKey to fetch fingerprints.

def pool_endpoints_with_cid(endpoints, keys):
    frames = []
    for k in keys:
        if k not in endpoints:
            continue
        f = endpoints[k].dropna(subset=["bin"])[["inchikey", "compound_chembl_id", "smiles", "bin"]]
        if len(f):
            frames.append(f)
    if not frames:
        return pd.DataFrame(columns=["inchikey", "compound_chembl_id", "smiles", "bin"])
    data = pd.concat(frames, ignore_index=True)
    data = (data.sort_values("bin", ascending=False)
                .groupby("inchikey", as_index=False)
                .agg(bin=("bin", "max"),
                     compound_chembl_id=("compound_chembl_id", "first"),
                     smiles=("smiles", "first")))
    return data


def cv_auroc(pool_df, ecfps, decoys_pool, params):
    """Relative modellability proxy: 4-fold CV AUROC on the pooled binary labels.

    Self-contained on the assay-measured molecules (no synthetic decoys) — this is used only to
    compare pools (does dropping a discordant assay change AUROC?), so absolute realism is not
    required. Returns None when there are too few actives or inactives to cross-validate.
    """
    rows = pool_df[pool_df["compound_chembl_id"].isin(ecfps)]
    n_pos, n_neg = int(rows["bin"].sum()), int((rows["bin"] == 0).sum())
    if n_pos < max(params["min_actives_to_model"], 4) or n_neg < 4:
        return None
    # subsample large pools (relative proxy only) — stratified, preserving the class ratio
    cap = params.get("max_model_rows", 20000)
    if len(rows) > cap:
        rng = np.random.RandomState(SEED)
        pos = rows[rows["bin"] == 1]
        neg = rows[rows["bin"] == 0]
        frac = cap / len(rows)
        n_pos_keep = max(4, int(round(len(pos) * frac)))
        n_neg_keep = max(4, int(round(len(neg) * frac)))
        pos = pos.iloc[rng.choice(len(pos), size=min(n_pos_keep, len(pos)), replace=False)]
        neg = neg.iloc[rng.choice(len(neg), size=min(n_neg_keep, len(neg)), replace=False)]
        rows = pd.concat([pos, neg])
    X = np.array([ecfps[c] for c in rows["compound_chembl_id"]])
    Y = rows["bin"].to_numpy().astype(int)
    mean_auroc, _ = KFoldTrain(X, Y, n_splits=4,
                               n_estimators=params["n_estimators_inloop"], random_state=SEED)
    return mean_auroc


def select_reference(group_keys, meta, conc, reference):
    """Pick the reference endpoint within a group: 'medoid' or 'largest'."""
    if reference == "largest":
        return max(group_keys, key=lambda k: meta[k]["n_cpds"])
    # medoid: highest shared-weighted mean concordance with overlapping neighbours
    best, best_score = None, -1.0
    for k in group_keys:
        num = den = 0.0
        for (a, b), (c, n) in conc.items():
            if k == a or k == b:
                num += c * n
                den += n
        score = (num / den) if den else -1.0
        if score > best_score or (score == best_score and (best is None or meta[k]["n_cpds"] > meta[best]["n_cpds"])):
            best, best_score = k, score
    return best if best is not None else max(group_keys, key=lambda k: meta[k]["n_cpds"])


def _neighbor_concordance(group_keys, conc):
    """Shared-weighted mean concordance of each endpoint with its overlap neighbours."""
    kset = set(group_keys)
    num = {k: 0.0 for k in group_keys}
    den = {k: 0.0 for k in group_keys}
    for (a, b), (c, n) in conc.items():
        if a in kset and b in kset:
            num[a] += c * n; den[a] += n
            num[b] += c * n; den[b] += n
    return {k: (num[k] / den[k]) if den[k] else None for k in group_keys}


def prune_discordant(group_keys, endpoints, meta, conc, reference, ecfps, decoys_pool, params, log):
    """Resolve each discordant assay pair by dropping the 'loser', guarded by modellability.

    For every overlapping pair whose concordance is below threshold, the loser is chosen by the
    reference policy — `largest`: drop the smaller assay; `medoid`: drop the one less concordant
    with its other neighbours (the consensus outlier). The loser is dropped only if removing it
    does not lower the pooled CV AUROC by more than the tolerance; otherwise both are kept.

    Returns (kept_keys, dropped_keys).
    """
    group_keys = list(group_keys)
    if len(group_keys) < 2:
        return group_keys, []

    nbr = _neighbor_concordance(group_keys, conc)
    kset = set(group_keys)
    discordant_pairs = sorted(
        [(p, c, n) for p, (c, n) in conc.items()
         if c < params["concordance_threshold"] and p[0] in kset and p[1] in kset],
        key=lambda t: t[1])

    kept = set(group_keys)
    dropped = []
    base = cv_auroc(pool_endpoints_with_cid(endpoints, list(kept)), ecfps, decoys_pool, params)
    checks = 0
    for (a, b), c, n in discordant_pairs:
        if a not in kept or b not in kept:
            continue  # one already dropped resolving an earlier pair
        # choose the loser per policy
        if reference == "largest":
            loser = a if meta[a]["n_cpds"] <= meta[b]["n_cpds"] else b
        else:  # medoid: drop the consensus outlier (lower neighbour concordance)
            ca, cb = nbr.get(a) or -1, nbr.get(b) or -1
            loser = a if ca <= cb else b
        if checks >= params["max_model_checks"]:
            log.append("      model-check cap reached; keeping remaining discordant pairs")
            break
        trial = list(kept - {loser})
        auroc_wo = cv_auroc(pool_endpoints_with_cid(endpoints, trial), ecfps, decoys_pool, params)
        checks += 1
        if base is None or auroc_wo is None or auroc_wo >= base - params["auroc_tolerance"]:
            kept.discard(loser)
            dropped.append(loser)
            base = auroc_wo if auroc_wo is not None else base
            log.append(f"      drop {loser[0]} ({meta[loser]['activity_type']}) — conc={c:.2f} "
                       f"n={n} vs {(b if loser==a else a)[0]}; AUROC={base}")
        else:
            log.append(f"      keep both {a[0]}/{b[0]} conc={c:.2f} n={n} "
                       f"(dropping {loser[0]} would lower AUROC {auroc_wo:.3f}<{base:.3f})")
    return list(kept), dropped


# ===========================================================================
# Strategy builders
# ===========================================================================

def _dataset_record(name, bucket, keys, endpoints, meta, ecfps, decoys_pool, params, dropped=None):
    pool = pool_endpoints_with_cid(endpoints, keys)
    auroc = cv_auroc(pool, ecfps, decoys_pool, params)
    return dict(
        name=name, bucket=bucket,
        n_assays=len(keys), n_dropped=len(dropped or []),
        n_mols=len(pool), n_actives=int(pool["bin"].sum()) if len(pool) else 0,
        auroc=auroc,
        data=pool, keys=list(keys), dropped=list(dropped or []),
    )


def build_S1(endpoints, meta, conc, reference, ecfps, decoys_pool, params, log):
    datasets = []
    for bucket in ["DR", "SP"]:
        keys = [k for k, m in meta.items() if m["bucket"] == bucket and m["has_cutoff"]]
        if not keys:
            continue
        log.append(f"  [{bucket}] {len(keys)} endpoints -> pruning")
        kept, dropped = prune_discordant(keys, endpoints, meta, conc, reference,
                                         ecfps, decoys_pool, params, log)
        datasets.append(_dataset_record(f"ORG_{bucket}", bucket, kept, endpoints, meta,
                                        ecfps, decoys_pool, params, dropped))
    return datasets


def build_S2(endpoints, meta, conc, reference, ecfps, decoys_pool, params, log):
    datasets = []
    families = sorted({m["family"] for m in meta.values() if m["has_cutoff"]})
    for fam in families:
        keys = [k for k, m in meta.items() if m["family"] == fam and m["has_cutoff"]]
        if not keys:
            continue
        bucket = meta[keys[0]]["bucket"]
        log.append(f"  [{fam}] {len(keys)} endpoints -> pruning")
        kept, dropped = prune_discordant(keys, endpoints, meta, conc, reference,
                                         ecfps, decoys_pool, params, log)
        datasets.append(_dataset_record(fam, bucket, kept, endpoints, meta,
                                        ecfps, decoys_pool, params, dropped))
    return datasets


def _connected_components(keys, conc, params):
    """Union-find over high-concordance edges; returns list of components (lists of keys)."""
    parent = {k: k for k in keys}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    kset = set(keys)
    for (a, b), (c, n) in conc.items():
        if a in kset and b in kset and c >= params["concordance_threshold"]:
            union(a, b)
    comps = {}
    for k in keys:
        comps.setdefault(find(k), []).append(k)
    return list(comps.values())


def build_S3(endpoints, meta, conc, reference, ecfps, decoys_pool, params, log):
    """Adaptive clustering: high-concordance connected components become datasets;
    isolated/non-overlapping endpoints are pooled into a per-bucket residual dataset."""
    datasets = []
    for bucket in ["DR", "SP"]:
        keys = [k for k, m in meta.items() if m["bucket"] == bucket and m["has_cutoff"]]
        if not keys:
            continue
        comps = _connected_components(keys, conc, params)
        multi = sorted([c for c in comps if len(c) > 1], key=len, reverse=True)
        singletons = [c[0] for c in comps if len(c) == 1]

        # cap: keep the largest clusters, merge the rest into the residual pool
        cap = params["max_datasets_per_bucket"] - 1  # reserve one slot for the residual
        if len(multi) > cap:
            overflow = multi[cap:]
            multi = multi[:cap]
            singletons += [k for comp in overflow for k in comp]
        log.append(f"  [{bucket}] {len(keys)} endpoints -> {len(multi)} clusters + "
                   f"{len(singletons)} pooled into residual")
        for i, comp in enumerate(multi):
            datasets.append(_dataset_record(f"ORG_{bucket}_cluster{i}", bucket, comp, endpoints,
                                            meta, ecfps, decoys_pool, params))
        if singletons:
            datasets.append(_dataset_record(f"ORG_{bucket}_residual", bucket, singletons, endpoints,
                                            meta, ecfps, decoys_pool, params))
    return datasets


STRATEGIES = {"S1": build_S1, "S2": build_S2, "S3": build_S3}


def load_features_and_decoys(endpoints, pathogen_code, output_root, cache=True):
    """Load ECFP fingerprints for all compounds in the endpoints, plus a decoy pool.

    The per-pathogen ECFP subset is cached (additive) under
    `output/<code>/22_pooled/_cache_ecfps.npz`, since the source HDF5 is 1.4 GB and scanning it
    for every run is the dominant cost.
    """
    cids = set()
    for ep in endpoints.values():
        cids.update(ep["compound_chembl_id"].tolist())

    cache_path = os.path.join(output_root, pathogen_code, "22_pooled", "_cache_ecfps.npz")
    if cache and os.path.exists(cache_path):
        z = np.load(cache_path, allow_pickle=False)
        ecfps = {cid: fp for cid, fp in zip(z["ids"].astype(str), z["fps"])}
    else:
        ecfps = load_ecfp_subset_by_chembl_id(ECFP_H5, cids)
        if cache and ecfps:
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            ids = np.array(list(ecfps.keys()))
            fps = np.array([ecfps[i] for i in ids])
            np.savez_compressed(cache_path, ids=ids, fps=fps)

    decoys_pool = list(ecfps.keys())  # unused by the current cv_auroc proxy; kept for signature
    return ecfps, decoys_pool
