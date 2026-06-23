"""
Build a DRAFT strain-equivalence table for one pathogen from the BacDive API
(authoritative culture-collection cross-references) — for human review.

What it does:
  1. Taxon-search BacDive for the pathogen's species -> every BacDive record for it.
  2. Read each record's strain designation + ALL culture-collection numbers.
  3. Group records by (species, normalized strain designation) -> one canonical strain,
     unioning every collection number and the designation as equivalent members.
  4. Match each canonical group back to the data's normalized strain keys (by designation
     OR collection number) and report how many assays it would unify; write a draft CSV.
  5. List the data strain keys that did NOT match any BacDive record (manual-review tail).

This is the authoritative alternative to mining ChEMBL co-occurrence (which over-merges).
It performs NO transitive closure beyond what BacDive states for a single strain.

Usage:
    python tools/build_strain_equivalences.py mtuberculosis

Output (loaded at runtime by strain_norm_key via load_strain_equivalences):
    config/strains/strain_equivalences_<pathogen>.csv
Manual curation seed (takes precedence over BacDive; can split/reassign groups):
    config/strains/strain_equivalences_manual.csv
"""
import os, re, sys, json, time, urllib.parse, urllib.request
from collections import Counter, defaultdict

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from pathogen_utils import load_pathogen, normalize_strain, strain_norm_key, prefix_key  # noqa: E402
import pandas as pd  # noqa: E402

API = "https://api.bacdive.dsmz.de/v2"
OUT = os.path.join(root, "..", "output")
STRAINS = os.path.join(root, "..", "config", "strains")          # per-pathogen equivalence files
MANUAL = os.path.join(STRAINS, "strain_equivalences_manual.csv")  # curation seed (all pathogens)
COLS = ["pathogen", "canonical_id", "canonical_label", "species", "type_strain",
        "member_value", "member_key", "match_type", "data_assays", "source", "verified_by", "notes"]


def get_json(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "ersilia-strain-curation"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def taxon_ids(genus, species):
    """All BacDive IDs for a (genus[, species]), following pagination."""
    path = f"{API}/taxon/{urllib.parse.quote(genus)}"
    if species:
        path += f"/{urllib.parse.quote(species)}"
    ids, url = [], path
    while url:
        res = get_json(url)
        ids.extend(res.get("results", []) or [])
        url = res.get("next")
        time.sleep(0.2)
    return ids


def all_collection_nos(rec):
    found = []
    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if k == "culture collection no.":
                    found.extend(str(v).split(","))
                walk(v)
        elif isinstance(o, list):
            for i in o:
                walk(i)
    walk(rec)
    return sorted({x.strip() for x in found if x.strip()})


def main(pathogen):
    df = pd.read_csv(os.path.join(OUT, pathogen, "09_assays_parameters_full.csv"), low_memory=False).drop_duplicates("assay_id")
    # data strain keys + assay counts (what we are trying to unify)
    keys = strain_norm_key(df["assay_strain_curated"], df["atcc_id"])
    data_counts = Counter(k for k in keys if k)
    # prefix-preserving counts (for 'prefixed' split rules): count assays by the
    # original effective value WITH its collection acronym kept.
    prefix_counts = Counter()
    for cur, atc in zip(df["assay_strain_curated"].fillna("").astype(str),
                        df["atcc_id"].fillna("").astype(str)):
        eff = cur.strip() or atc.strip()
        pk = prefix_key(eff)
        if pk:
            prefix_counts[pk] += 1

    binomial = load_pathogen(pathogen).split()
    genus, species = binomial[0], (binomial[1] if len(binomial) > 1 else "")
    print(f"{pathogen}: taxon search BacDive for {genus} {species} ...")
    all_ids = sorted(set(taxon_ids(genus, species)))
    print(f"  {len(all_ids)} BacDive records; fetching...")

    # fetch records in batches of 100
    records = {}
    for i in range(0, len(all_ids), 100):
        batch = ";".join(str(x) for x in all_ids[i:i + 100])
        res = get_json(f"{API}/fetch/{batch}")
        records.update(res.get("results", {}))
        time.sleep(0.2)

    # group records by (species, normalized designation)
    groups = defaultdict(lambda: {"ids": set(), "numbers": set(), "designations": set(), "species": "", "type": set()})
    for rid, rec in records.items():
        nt = rec.get("Name and taxonomic classification", {})
        species = nt.get("species", "") or ""
        design = nt.get("strain designation", "") or ""
        gkey = (species, normalize_strain(design)) if design else (species, f"BACDIVE:{rid}")
        g = groups[gkey]
        g["ids"].add(int(rec.get("General", {}).get("BacDive-ID", rid)))
        g["numbers"].update(all_collection_nos(rec))
        if design:
            g["designations"].add(design)
        g["species"] = species
        ts = nt.get("type strain", "")
        if ts:
            g["type"].add(ts)

    # --- Curated manual rows FIRST: they take precedence over BacDive, so a manual
    # row can split or reassign a BacDive grouping (any key it claims is removed from
    # BacDive below). match_type "exact" matches one key; "contains" lumps every data
    # key whose normalized form contains the (normalized) member_value substring.
    # "exact" matches one base key; "contains" lumps every base key containing the
    # substring; "prefixed" SPLITS a same-number strain out by its collection acronym
    # (e.g. LMG 17978 ≠ ATCC 17978) — it claims a prefix-preserving key, not a base key.
    manual_out, contains_canon, manual_canon, manual_keys, prefixed_keys = [], set(), set(), set(), set()
    if os.path.exists(MANUAL):
        man = pd.read_csv(MANUAL)
        for _, r in man[man["pathogen"] == pathogen].iterrows():
            mt = (str(r.get("match_type", "")).strip() or "exact").lower()
            base = dict(pathogen=pathogen, canonical_id=r["canonical_id"], canonical_label=r["canonical_label"],
                        species=load_pathogen(pathogen), type_strain="", source="manual",
                        verified_by=r.get("verified_by", ""), notes=r.get("notes", ""))
            if mt == "contains":
                sub = normalize_strain(r["member_value"])
                matched = sorted(k for k in data_counts if sub and sub in k)
                if not matched:
                    continue
                manual_out.append({**base, "member_value": r["member_value"], "member_key": sub,
                                   "match_type": "contains", "data_assays": sum(data_counts[k] for k in matched),
                                   "notes": f"{base['notes']} [contains '{sub}': {len(matched)} keys, e.g. {', '.join(matched[:6])}]"})
                contains_canon.add(r["canonical_id"]); manual_canon.add(r["canonical_id"]); manual_keys |= set(matched)
            elif mt == "prefixed":
                pk = prefix_key(r["member_value"])
                if prefix_counts.get(pk, 0) <= 0:
                    continue
                manual_out.append({**base, "member_value": r["member_value"], "member_key": pk,
                                   "match_type": "prefixed", "data_assays": prefix_counts.get(pk, 0)})
                manual_canon.add(r["canonical_id"]); prefixed_keys.add(pk)
                # NB: do NOT claim the base key — ATCC/bare forms stay with their own group
            else:
                mk = normalize_strain(r["member_value"])
                if data_counts.get(mk, 0) <= 0:
                    continue
                manual_out.append({**base, "member_value": r["member_value"], "member_key": mk,
                                   "match_type": "exact", "data_assays": data_counts.get(mk, 0)})
                manual_canon.add(r["canonical_id"]); manual_keys.add(mk)
    manual_rows = pd.DataFrame(manual_out, columns=COLS)
    matched_keys = set(manual_keys)

    # assays pulled out of their base-key group by a 'prefixed' split rule — subtract
    # these from the base-key member counts so data_assays stays accurate.
    split_by_base = Counter()
    if prefixed_keys:
        for cur, atc in zip(df["assay_strain_curated"].fillna("").astype(str),
                            df["atcc_id"].fillna("").astype(str)):
            eff = cur.strip() or atc.strip()
            if eff and prefix_key(eff) in prefixed_keys:
                split_by_base[normalize_strain(eff)] += 1

    # --- BacDive groups: emit members present in our data, skipping any base key
    # claimed by a manual row, and any exact spelling split off by a 'prefixed' rule
    # (so e.g. the LMG 17978 member is removed from the ATCC 17978 group).
    rows = []
    for (sp, dkey), g in groups.items():
        canonical_id = "BacDive:" + ";".join(str(x) for x in sorted(g["ids"]))
        label = sorted(g["designations"], key=len)[-1] if g["designations"] else sp
        members = set(g["numbers"]) | set(g["designations"])
        member_keys = {m: normalize_strain(m) for m in members}
        covered = {mk for mk in member_keys.values() if mk in data_counts and mk not in manual_keys}
        if not covered:
            continue
        matched_keys |= covered
        for m in sorted(members):
            mk = member_keys[m]
            if mk not in data_counts or mk in manual_keys or prefix_key(m) in prefixed_keys:
                continue
            n = data_counts.get(mk, 0) - split_by_base.get(mk, 0)
            if n <= 0:  # every assay of this base key was split off by a prefixed rule
                continue
            rows.append({
                "pathogen": pathogen, "canonical_id": canonical_id, "canonical_label": label,
                "species": sp, "type_strain": ";".join(sorted(g["type"])),
                "member_value": m, "member_key": mk, "match_type": "exact",
                "data_assays": n, "source": "bacdive", "verified_by": "", "notes": "",
            })
    bacdive_rows = pd.DataFrame(rows, columns=COLS)

    out = pd.concat([manual_rows, bacdive_rows], ignore_index=True)
    out = out.drop_duplicates(subset=["canonical_id", "member_key"])
    # drop inert singleton groups (one key unifies nothing) — but keep 'contains'
    # rules (they lump many keys) and any explicitly curated manual group.
    grp_n = out.groupby("canonical_id")["member_key"].transform("nunique")
    out = out[(grp_n >= 2) | (out["canonical_id"].isin(contains_canon | manual_canon))]

    # Canonical label: ALWAYS the ATCC ID when the group has one (dominant ATCC member
    # by assays); else, for BacDive groups, the dominant in-data spelling (BacDive's
    # raw designation is often a historical name like "Seattle 1945"); else keep the
    # curator's manual label.
    for cid, grp in out.groupby("canonical_id"):
        atcc = grp[grp["member_value"].astype(str).str.upper().str.startswith("ATCC")]
        if len(atcc):
            out.loc[out["canonical_id"] == cid, "canonical_label"] = atcc.loc[atcc["data_assays"].idxmax(), "member_value"]
        elif (grp["source"] == "bacdive").all():
            out.loc[out["canonical_id"] == cid, "canonical_label"] = grp.loc[grp["data_assays"].idxmax(), "member_value"]

    # Order so each group's rows are contiguous, groups ranked by total size, and
    # the biggest member first within a group — easy to scan/validate top-down.
    out["group_assays"] = out.groupby("canonical_id")["data_assays"].transform("sum")
    out = out.sort_values(
        ["group_assays", "canonical_label", "canonical_id", "data_assays", "member_value"],
        ascending=[False, True, True, False, True],
    )
    out = out[["group_assays"] + [c for c in COLS]]
    os.makedirs(STRAINS, exist_ok=True)
    dest = os.path.join(STRAINS, f"strain_equivalences_{pathogen}.csv")
    out.to_csv(dest, index=False)

    # summary
    grp_assays = out.groupby(["canonical_id", "canonical_label"])["data_assays"].sum().sort_values(ascending=False)
    print(f"\n  {len(grp_assays)} canonical strains matched to data. Top groups (assays unified):")
    for (cid, label), a in grp_assays.head(15).items():
        keys_in = out[out.canonical_id == cid]["member_value"].tolist()
        print(f"    {a:5d}  {label:14} {cid:24}  data: {', '.join(keys_in)[:60]}")

    # data strain keys NOT covered by BacDive (manual-review tail), ranked by assays
    unmatched = [(k, c) for k, c in data_counts.items() if k not in matched_keys]
    unmatched.sort(key=lambda kc: -kc[1])
    covered_assays = sum(data_counts[k] for k in matched_keys)
    total_assays = sum(data_counts.values())
    print(f"\n  Coverage: {covered_assays}/{total_assays} strain-keyed assays "
          f"({100*covered_assays/max(total_assays,1):.0f}%) matched to BacDive.")
    print(f"  Top unmatched keys (need manual rows or another source):")
    for k, c in unmatched[:15]:
        print(f"    {c:5d}  {k}")
    print(f"\nDraft saved -> {dest}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "mtuberculosis")
