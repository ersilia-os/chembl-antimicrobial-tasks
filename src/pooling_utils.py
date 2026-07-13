"""
Pooling utilities (dataset-pooling rebuild, steps 23+).

*Transfer* signal: train a RandomForest on one dataset A and use it to predict
another dataset B's molecules that are NOT shared with A (a true generalisation
test). A high AUROC ("A predicts B well") is evidence that A and B can be pooled.
The signal is directional; for a pair with too few disjoint test molecules (or a
single-class test set) the cell is left undefined (NaN).

This clean module depends only on the current src tree (`default`).

Functions
---------
morgan_matrix(inchikey_to_smiles, radius, n_bits)
    Morgan fingerprints (uint8 bit vectors) per InChIKey, from SMILES.
transfer_matrix(datasets, fps, min_test)
    Directional train-A -> predict-B (disjoint) AUROC / F1 / count matrices.
"""

import os
import sys
import glob
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from scipy.sparse import csr_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, f1_score

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, root)
from default import FP_RADIUS, FP_BITS, MIN_TEST, RANDOM_SEED


# --- fingerprints -----------------------------------------------------------

def _fingerprint(smiles, radius, n_bits):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    bv = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(bv, arr)
    return arr


def morgan_matrix(inchikey_to_smiles, radius=FP_RADIUS, n_bits=FP_BITS):
    """Return {inchikey: uint8[n_bits]} Morgan fingerprints (invalid SMILES skipped)."""
    fps = {}
    for ik, smi in inchikey_to_smiles.items():
        fp = _fingerprint(smi, radius, n_bits)
        if fp is not None:
            fps[ik] = fp
    return fps


# --- transfer matrix --------------------------------------------------------
# Worker globals for the parallel transfer matrix (inherited via fork on Linux, so
# the large universe array is shared copy-on-write rather than pickled per task).
_TX = _TD = _TL = None
_TMIN = MIN_TEST


def _init_transfer(Xuni, didx, dlab, min_test):
    global _TX, _TD, _TL, _TMIN
    _TX, _TD, _TL, _TMIN = Xuni, didx, dlab, min_test


def _transfer_row(ai):
    """Compute one trainer A's row of the transfer matrix (worker process)."""
    Xuni, didx, dlab, min_test = _TX, _TD, _TL, _TMIN
    n, U = len(didx), Xuni.shape[0]
    auroc = np.full(n, np.nan, dtype=np.float32)
    f1 = np.full(n, np.nan, dtype=np.float32)
    nn = np.zeros(n, dtype=np.int32)
    idxA, yA = didx[ai], dlab[ai]
    if idxA.size == 0 or yA.min() == yA.max():
        return ai, auroc, f1, nn
    rf = RandomForestClassifier(n_estimators=100, max_features="sqrt", n_jobs=1,
                                random_state=RANDOM_SEED).fit(Xuni[idxA], yA)
    prob = rf.predict_proba(Xuni)[:, 1]
    A_mask = np.zeros(U, dtype=bool)
    A_mask[idxA] = True
    for bj in range(n):
        if bj == ai:
            continue
        idxB, yB = didx[bj], dlab[bj]
        if idxB.size == 0:
            continue
        keep = ~A_mask[idxB]                       # B's molecules NOT in A (disjoint test)
        k = int(keep.sum())
        if k < min_test:
            continue
        yt = yB[keep]
        if yt.min() == yt.max():                   # single-class test set -> AUROC undefined
            continue
        pb = prob[idxB[keep]]
        auroc[bj] = roc_auc_score(yt, pb)
        f1[bj] = f1_score(yt, (pb >= 0.5).astype(int), zero_division=0)
        nn[bj] = k
    return ai, auroc, f1, nn


def transfer_matrix(datasets, fps, min_test=MIN_TEST, log=print):
    """Directional model-transfer matrix on disjoint molecules.

    Parameters
    ----------
    datasets : list of (id, inchikeys: np.ndarray, labels: np.ndarray)
        Non-degenerate datasets only (each has both classes).
    fps : dict {inchikey: fingerprint}
    min_test : int
        Minimum disjoint test molecules (with both classes) to score a cell.

    Returns
    -------
    ids : list[str]
    AUROC, F1 : (n, n) float32   M[A, B] = transfer from trainer A to test B (NaN if unscored)
    N : (n, n) int32             number of disjoint test molecules scored
    """
    ids = [d[0] for d in datasets]
    n = len(datasets)
    uni = sorted({ik for _, iks, _ in datasets for ik in iks if ik in fps})
    pos = {ik: i for i, ik in enumerate(uni)}
    Xuni = np.vstack([fps[ik] for ik in uni]) if uni else np.empty((0, FP_BITS), dtype=np.uint8)

    didx, dlab = [], []
    for _, iks, ys in datasets:
        ii = [pos[ik] for ik in iks if ik in pos]
        ll = [y for ik, y in zip(iks, ys) if ik in pos]
        didx.append(np.array(ii, dtype=int))
        dlab.append(np.array(ll, dtype=int))

    AUROC = np.full((n, n), np.nan, dtype=np.float32)
    F1 = np.full((n, n), np.nan, dtype=np.float32)
    N = np.zeros((n, n), dtype=np.int32)
    if n == 0:
        return ids, AUROC, F1, N

    # Parallelise across trainers A (embarrassingly parallel). Workers inherit the
    # universe array via fork (copy-on-write); each worker RF uses n_jobs=1 to avoid
    # oversubscribing cores.
    n_workers = max(1, min(12, cpu_count() - 2))
    _init_transfer(Xuni, didx, dlab, min_test)     # set globals; workers inherit via fork
    done = 0
    with Pool(processes=n_workers) as pool:
        for ai, a_row, f_row, n_row in pool.imap_unordered(_transfer_row, range(n), chunksize=4):
            AUROC[ai] = a_row
            F1[ai] = f_row
            N[ai] = n_row
            done += 1
            if done % 200 == 0:
                log(f"    transfer: {done}/{n} trainers ({n_workers} workers)")
    return ids, AUROC, F1, N


# --- pools ------------------------------------------------------------------

def _transfer_adjacency(AUROC, threshold, rule):
    with np.errstate(invalid="ignore"):
        if rule == "bidirectional":
            M = np.minimum(AUROC, AUROC.T)          # both directions required
        elif rule == "unidirectional":
            M = np.fmax(AUROC, AUROC.T)             # either scored direction counts
        else:
            raise ValueError(f"unknown rule: {rule!r}")
        adj = M >= threshold                        # NaN >= t -> False
    np.fill_diagonal(adj, False)
    return adj


def greedy_modularity_pools(AUROC, threshold, rule="bidirectional"):
    """Pools = greedy-modularity communities of the transfer graph at `threshold`.

    Unlike connected components, modularity does not percolate a graph into one blob via
    transitive weak links; it partitions it into communities and stops at the modularity
    peak (pool count is emergent, no internal threshold). Isolated datasets become
    singleton communities. Community ids are assigned largest-first (0 = biggest pool).

    Returns
    -------
    labels : np.ndarray (n,)   community id per dataset
    """
    n = AUROC.shape[0]
    adj = _transfer_adjacency(AUROC, threshold, rule)
    G = nx.from_scipy_sparse_array(csr_matrix(np.triu(adj)))
    labels = np.empty(n, dtype=int)
    for cid, comm in enumerate(sorted(greedy_modularity_communities(G), key=len, reverse=True)):
        for node in comm:
            labels[node] = cid
    return labels


# --- pool building (shared by the pooling steps 23-26) ----------------------

def load_binarised_datasets(bdir):
    """Load every binarised dataset in one pathogen/category directory.

    Returns
    -------
    datasets : list of dicts {id, iks, labels, values, relations, ikset}
    ik2info  : {inchikey: (smiles, compound_chembl_id)}   first occurrence wins
    """
    datasets, ik2info = [], {}
    for f in sorted(glob.glob(os.path.join(bdir, "*.csv.gz"))):
        dsid = os.path.basename(f)[:-len(".csv.gz")]
        df = pd.read_csv(f).drop_duplicates("inchikey")
        iks = df["inchikey"].to_numpy()
        for ik, smi, cid in zip(iks, df["smiles"].to_numpy(), df["compound_chembl_id"].to_numpy()):
            ik2info.setdefault(ik, (smi, cid))
        rels = (df["relation"].astype(str).to_numpy() if "relation" in df.columns
                else np.array(["="] * len(iks)))
        datasets.append({"id": dsid, "iks": iks,
                         "labels": df["bin"].to_numpy().astype(int),
                         "values": df["value"].to_numpy(dtype=float),
                         "relations": rels,
                         "ikset": frozenset(iks)})
    return datasets, ik2info


def load_dataset_meta(summary_csv, category, cutoffs, unit_slug):
    """Map dataset id -> (unit, direction, cutoff) from a 22_binarisation_summary.csv.

    The cutoff is taken from the value STEP 22 recorded for each dataset (`r.cutoff`), which has
    already resolved any per-pathogen override. This is the single source of truth: pool
    re-binarisation then uses exactly the cutoff the binarised file used, so an override is never
    silently dropped by a caller that loaded `cutoffs` without the pathogen code. The passed
    `cutoffs` {unit: cutoff} dict is only a fallback for a dataset whose recorded cutoff is missing.

    `unit_slug` is passed in (dataset_utils.UNIT_SLUG) to keep this module decoupled.
    """
    summ = pd.read_csv(summary_csv)
    summ = summ[summ["category"] == category]
    meta = {}
    for r in summ.itertuples(index=False):
        stem = f"{r.assay_id}_{str(r.activity_type).replace('/', '-')}_{unit_slug.get(r.unit, str(r.unit))}"
        cutoff = r.cutoff if pd.notna(r.cutoff) else cutoffs.get(r.unit)
        meta[stem] = (r.unit, int(r.direction), cutoff)
    return meta


def merge_datasets(datasets, meta, ik2info, report_diag=False):
    """Merge a set of datasets to ONE row per molecule via the value-resolution policy.

    For each molecule, group its observations by (unit, direction, cutoff); pick the group
    with highest priority ('%' first, then most observations - implements "% prevails" for
    the rare cross-unit SP molecule and averaging-within-unit everywhere else); average the
    values in that group and re-binarise CONSERVATIVELY against the unit's cutoff/direction.

    Conservative (qualifier-aware) re-binarisation mirrors step 22's `binarize_conservative`:
    a molecule is active (1) only if the AVERAGED value lies on the active side of the cutoff
    AND at least one contributing measurement PROVES a hit (an uncensored/right-way qualifier
    on the active side). Censored records (e.g. "< cutoff") and value-less rows never count as
    active. This makes a single-dataset merge reproduce the stored step-22 bins exactly.

    `datasets` is a list of dicts each carrying `id`, `iks`, `values`, `relations`.

    If `report_diag` is True, also return a dict of two transparency counters (labels are
    never changed by this):
      - `n_flipped`       : molecules with `proven and not on_side` - a proving measurement
                            exists but averaging dragged the mean to the inactive side.
      - `n_proof_demoted` : molecules with `on_side and not proven` - the averaged value is on
                            the active side but no measurement proves the hit, so the
                            conservative rule keeps them inactive (a plain threshold would not).
    Returns `(df, diag)` when set, else just `df`.
    """
    permol = defaultdict(lambda: defaultdict(list))       # ik -> (unit, direction, cutoff) -> [(value, relation)]
    for d in datasets:
        unit, direction, cutoff = meta[d["id"]]
        rels = d.get("relations")
        if rels is None:
            rels = ["="] * len(d["iks"])
        for ik, v, rel in zip(d["iks"], d["values"], rels):
            permol[ik][(unit, direction, cutoff)].append((v, str(rel)))

    rows, n_flipped, n_proof_demoted = [], 0, 0
    for ik, groups in permol.items():
        # priority: '%' unit first, then the group with the most measurements
        key = max(groups, key=lambda k: (1 if k[0] == "%" else 0, len(groups[k])))
        unit, direction, cutoff = key
        obs = groups[key]
        vals = [v for v, _ in obs if v == v]              # drop NaN (value-less rows)
        avg = float(np.mean(vals)) if vals else float("nan")
        if direction < 0:
            on_side = (avg == avg) and (avg <= cutoff)
            proven = any((v == v) and (r in ("=", "<")) and (v <= cutoff) for v, r in obs)
        else:
            on_side = (avg == avg) and (avg >= cutoff)
            proven = any((v == v) and (r in ("=", ">")) and (v >= cutoff) for v, r in obs)
        b = int(on_side and proven)
        if proven and not on_side:                        # proving measurement averaged out
            n_flipped += 1
        elif on_side and not proven:                      # on active side but no proof (censored)
            n_proof_demoted += 1
        smi, cid = ik2info[ik]
        rows.append((ik, cid, smi, avg, unit, b))
    df = pd.DataFrame(rows, columns=["inchikey", "compound_chembl_id", "smiles",
                                     "value", "unit", "bin"])
    if report_diag:
        return df, {"n_flipped": n_flipped, "n_proof_demoted": n_proof_demoted}
    return df


# --- near-duplicate folding (step-23 pre-transfer) --------------------------

def _neardup_edges(datasets, overlap_thr, discord_thr, maxpost):
    """Edges between near-duplicate datasets, via a frequency-capped inverted index.

    Two datasets get an edge iff containment overlap |A n B| / min(|A|,|B|) >= overlap_thr AND
    the fraction of their SHARED molecules with disagreeing stored bins < discord_thr.

    Candidate pairs are generated from molecules present in <= maxpost datasets (ubiquitous
    molecules are skipped for candidate generation only - a genuine near-duplicate shares its
    rarer molecules too). Exact overlap / discordance are then computed on the FULL molecule sets.
    Returns list of (i, j, overlap, discordance) with i < j.
    """
    n = len(datasets)
    ikset = [d["ikset"] for d in datasets]
    bins = [dict(zip(d["iks"], d["labels"])) for d in datasets]
    freq = Counter()
    for s in ikset:
        for ik in s:
            freq[ik] += 1
    index = defaultdict(list)
    for i, s in enumerate(ikset):
        for ik in s:
            if freq[ik] <= maxpost:
                index[ik].append(i)

    edges = []
    for i in range(n):
        co = Counter()
        for ik in ikset[i]:
            if freq[ik] <= maxpost:
                for j in index[ik]:
                    if j > i:
                        co[j] += 1
        si = len(ikset[i])
        for j in co:                                   # only pairs sharing >= 1 non-ubiquitous mol
            inter = ikset[i] & ikset[j]
            ni = len(inter)
            if ni == 0:
                continue
            ov = ni / min(si, len(ikset[j]))
            if ov < overlap_thr:
                continue
            disc = sum(1 for ik in inter if bins[i][ik] != bins[j][ik]) / ni
            if disc < discord_thr:
                edges.append((i, j, ov, disc))
    return edges


def _greedy_clique_cover(n, edges):
    """Partition {0..n-1} into cliques of the near-dup graph (deterministic greedy peel).

    Repeatedly peel a maximal clique from the remaining nodes: seed with the highest-degree
    remaining node (ties -> lowest index), then grow by the candidate adjacent to the whole
    current clique with the most remaining connections (ties -> lowest index). Because every
    returned group is a clique, ALL within-group pairs satisfy the 90/20 edge test directly -
    unlike connected components, no pair is folded in transitively via a bridging dataset.

    Returns a list of index-lists (nodes with no edges come back as singletons).
    """
    nbr = {i: set() for i in range(n)}
    for a, b, *_ in edges:
        nbr[a].add(b)
        nbr[b].add(a)
    remaining = set(range(n))
    groups = []
    while remaining:
        seed = max(remaining, key=lambda x: (len(nbr[x] & remaining), -x))
        clique = {seed}
        cands = nbr[seed] & remaining
        while cands:
            v = max(cands, key=lambda x: (len(cands & nbr[x]), -x))
            clique.add(v)
            cands = cands & nbr[v]
        groups.append(sorted(clique))
        remaining -= clique
    return groups


def fold_near_duplicates(datasets, meta, ik2info, overlap_thr=0.90, discord_thr=0.20, maxpost=100):
    """Fold near-duplicate datasets (cliques of the near-dup graph) into single nodes.

    Runs over ALL datasets (eligible + degenerate + small). Datasets are grouped into CLIQUES
    of the 90/20 near-dup graph (see `_greedy_clique_cover`), so every pair inside a fold group
    satisfies the overlap/discordance test directly - a middle dataset can no longer transitively
    fold two datasets that fail the test between themselves. Each group with >1 member is merged
    via `merge_datasets` into one node; every node carries `members` (its original dataset ids) so
    downstream pooling can expand back to the source datasets.

    Returns
    -------
    nodes : list of dicts {id, iks, labels, values, ikset, members}
    fold_rows : list of provenance dicts (only for folded groups of >= 2 members)
    """
    n = len(datasets)
    edges = _neardup_edges(datasets, overlap_thr, discord_thr, maxpost)
    groups = _greedy_clique_cover(n, edges)
    edge_by_pair = {(e[0], e[1]): (e[2], e[3]) for e in edges}

    nodes, fold_rows = [], []
    for idxs in groups:
        if len(idxs) == 1:
            d = datasets[idxs[0]]
            node = dict(d)
            node["members"] = [d["id"]]
            nodes.append(node)
            continue
        orig = [datasets[i] for i in idxs]
        df = merge_datasets(orig, meta, ik2info)
        rep = max(orig, key=lambda d: len(d["ikset"]))["id"]
        iks = df["inchikey"].to_numpy()
        node = {"id": f"{rep}_merged{len(idxs)}", "iks": iks,
                "labels": df["bin"].to_numpy().astype(int),
                "values": df["value"].to_numpy(dtype=float),
                "relations": np.array(["="] * len(iks)),   # resolved node; unused downstream
                "ikset": frozenset(iks),
                "members": [datasets[i]["id"] for i in idxs]}
        nodes.append(node)
        ovs = [edge_by_pair[(a, b)][0] for a in idxs for b in idxs if (a, b) in edge_by_pair]
        dcs = [edge_by_pair[(a, b)][1] for a in idxs for b in idxs if (a, b) in edge_by_pair]
        fold_rows.append({"fold_id": node["id"], "n_members": len(idxs),
                          "members": "|".join(node["members"]),
                          "n_molecules": len(node["ikset"]),
                          "min_edge_overlap": round(min(ovs), 4) if ovs else float("nan"),
                          "max_edge_discordance": round(max(dcs), 4) if dcs else float("nan")})
    return nodes, fold_rows
