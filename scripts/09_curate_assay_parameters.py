# Enhanced Assay Information Extraction (Step 09)
# Extracts enhanced target information from ChEMBL assays using OLLAMA
# This script needs to be run with a GPU machine available

from pydantic import BaseModel
from typing import List
from tqdm import tqdm
import pandas as pd
import numpy as np
import ollama
import json
import csv
import sys
import os
import re

import subprocess
print("=== Python GPU Check ===\n")
try:
    result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
    print(result.stdout)
except Exception as e:
    print(f"nvidia-smi failed: {e}\n")

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH, LLM_MODEL

print("Step 09 - Enhanced Assay Information Extraction")
print("This script needs to run on a GPU-enabled machine.")

# Organism-specific ChEMBL ID mapping for all 16 supported pathogens
ORGANISM_CHEMBL_IDS = {
    "Mycobacterium tuberculosis": "CHEMBL360",
    "Acinetobacter baumannii": "CHEMBL614425",
    "Escherichia coli": "CHEMBL354",
    "Pseudomonas aeruginosa": "CHEMBL348",
    "Klebsiella pneumoniae": "CHEMBL350",
    "Streptococcus pneumoniae": "CHEMBL347",
    "Enterococcus faecium": "CHEMBL357",
    "Candida albicans": "CHEMBL366",
    "Plasmodium falciparum": "CHEMBL364",
    "Campylobacter": "CHEMBL612492",
    "Helicobacter pylori": "CHEMBL612600",
    "Neisseria gonorrhoeae": "CHEMBL614430",
    "Enterobacter": "CHEMBL614439",
    "Staphylococcus aureus": "CHEMBL352",
    "Schistosoma mansoni": "CHEMBL612893"
}

# Get pathogen code from command line
if len(sys.argv) != 2:
    print("Usage: python 09_extract_assay_details.py <pathogen_code>")
    sys.exit(1)

pathogen_code = sys.argv[1]

# Load pathogen info
df = pd.read_csv(os.path.join(CONFIGPATH, 'pathogens.csv'))
row = df.loc[df["code"].eq(pathogen_code)]
if row.empty:
    raise SystemExit(f"Unknown pathogen code: {pathogen_code}")
pathogen = row.iloc[0]["pathogen"]

# Validate pathogen is supported
if pathogen not in ORGANISM_CHEMBL_IDS:
    print(f"Warning: {pathogen} not found in organism ChEMBL ID mapping")
    print(f"Available organisms: {list(ORGANISM_CHEMBL_IDS.keys())}")
    organism_chembl_id = ""
else:
    organism_chembl_id = ORGANISM_CHEMBL_IDS[pathogen]

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

print(f"Processing pathogen: {pathogen} ({pathogen_code})")

# Load all ChEMBL source data once for efficient lookup
print("Loading ChEMBL source data...")

# Primary assay metadata
assays_data = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "assays.csv"), low_memory=False)
assay_info = assays_data.set_index("chembl_id")[["doc_id", "description", "assay_strain", "assay_organism",
                                                  "assay_category", "assay_test_type", "assay_cell_type", "tid"]].to_dict("index")

# Detailed assay descriptions
try:
    assay_descriptions = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "assay_descriptions.csv"))
    descriptions = assay_descriptions.set_index("chembl_id")["description"].to_dict()
except:
    descriptions = {}
    print("  Warning: assay_descriptions.csv not found, skipping detailed descriptions")

# Publication information
docs_data = pd.read_csv(os.path.join(DATAPATH, "chembl_activities", "docs.csv"), low_memory=False)
doc_info = docs_data.set_index("doc_id")[["title", "abstract", "journal", "pubmed_id", "doi"]].to_dict("index")

# Target information — loaded from pre-joined file produced by step 00
_tds_path = os.path.join(DATAPATH, "chembl_processed", "00_target_dictionary_synonyms.csv")
target_dictionary = pd.read_csv(_tds_path, low_memory=False)
target_dictionary = target_dictionary[
    target_dictionary["organism"].str.contains(pathogen, case=False, na=False)
]
target_info = target_dictionary.set_index("tid")[["pref_name", "organism", "target_type", "chembl_id"]].to_dict("index")
target_synonym_map = (
    target_dictionary.set_index("tid")["synonyms"]
    .apply(lambda s: s.split(";") if isinstance(s, str) and s else [])
    .to_dict()
)

# Pre-built index and prefix list for resolve_target_chembl_id (O(1) lookup)
def _build_lookup_index(df):
    index = {}
    for _, row in df.iterrows():
        cid = str(row["chembl_id"])
        for name in [row.get("pref_name", "")] + str(row.get("synonyms", "") or "").split(";"):
            key = str(name or "").strip().lower()
            if key and key not in index:
                index[key] = cid
    return index

def _get_organism_prefixes(df):
    orgs = df["organism"].dropna().unique().tolist()
    return sorted(set(orgs), key=len, reverse=True)

def _strip_organism_prefix(name, prefixes):
    lower = name.lower()
    for org in prefixes:
        if lower.startswith(org.lower()):
            rest = name[len(org):].strip()
            if rest:
                return rest
    return name

_target_lookup = _build_lookup_index(target_dictionary)
_org_prefixes = _get_organism_prefixes(target_dictionary)

print(f"  Loaded {len(assay_info)} assays, {len(descriptions)} detailed descriptions")
print(f"  Loaded {len(doc_info)} documents, {len(target_info)} targets")

# Assay type mapping
assay_type_map = {"F": "Functional", "B": "Binding", "T": "Toxicity", "A": "ADME", "P": "Physicochemical", "U": "Uncategorized"}

AMBIGUOUS_TARGET_TYPES = {"UNCHECKED", "NON-MOLECULAR", "NON-PROTEIN TARGET"}
KNOWN_PROTEIN_TYPES = {"SINGLE PROTEIN", "PROTEIN COMPLEX", "PROTEIN FAMILY"}


class AssayDetails(BaseModel):
    """Pydantic schema for structured extraction of enhanced assay information"""
    target_type_curated: str           # More specific than ChEMBL's basic classification
    target_name_curated: str            # Extracted from descriptions/abstracts
    target_chembl_id_curated: str      # Cross-referenced with descriptions
    assay_organism_curated: str               # Validated organism name (what organism the assay is actually performed in/on)
    assay_strain_curated: str           # Standardized strain identifier
    atcc_id: str                        # ATCC catalog numbers (e.g., "ATCC 27294")
    mutations: List[str]                # Specific mutations (e.g., ["S450L", "H526Y"])
    known_drug_resistances: List[str]   # Drugs with explicitly stated resistance (e.g., ["rifampicin"])
    culture_media: str                  # Growth medium (e.g., "Middlebrook 7H9")


class AssayDetailsBatch(BaseModel):
    """Batch processing schema for multiple assays"""
    results: List[AssayDetails]


class AssayDetailsTextOnly(BaseModel):
    """Fields for assays where target identity is already known from structured data."""
    assay_organism_curated: str
    assay_strain_curated: str
    atcc_id: str
    mutations: List[str]
    known_drug_resistances: List[str]
    culture_media: str


class AssayDetailsTextOnlyBatch(BaseModel):
    results: List[AssayDetailsTextOnly]


def classify_assay(row):
    """Return 'text_only' or 'full_inference' based on how much the LLM needs to infer."""
    tt = str(row.get("target_type", "")).strip()
    target_chembl_id = str(row.get("target_chembl_id", "")).strip()
    bao = str(row.get("bao_label", "")).strip().lower()
    assay_org = str(row.get("assay_organism", "")).strip()

    if tt in KNOWN_PROTEIN_TYPES and target_chembl_id != "CHEMBL612545":
        return "text_only"
    if tt == "ORGANISM" and pathogen.lower() in assay_org.lower():
        return "text_only"
    if tt in AMBIGUOUS_TARGET_TYPES and "organism-based" in bao and pathogen.lower() in assay_org.lower():
        return "text_only"
    return "full_inference"


def prefill_target_fields(row):
    """Return pre-filled target fields for Group A assays (target already known)."""
    tt = str(row.get("target_type", "")).strip()
    if tt in KNOWN_PROTEIN_TYPES:
        assay_id = row.get("assay_id", "")
        ai = assay_info.get(assay_id, {})
        tid = ai.get("tid") if isinstance(ai, dict) else None
        ti = target_info.get(tid, {}) if tid and not pd.isna(tid) else {}
        return {
            "target_type_curated": tt,
            "target_chembl_id_curated": str(row.get("target_chembl_id", "")),
            "target_name_curated": ti.get("pref_name", "") if isinstance(ti, dict) else "",
            "assay_organism_curated": str(row.get("assay_organism", "")),
        }
    else:  # ORGANISM or UNCHECKED+organism-based
        return {
            "target_type_curated": "ORGANISM",
            "target_chembl_id_curated": organism_chembl_id,
            "target_name_curated": pathogen,
            "assay_organism_curated": str(row.get("assay_organism", "")),
        }


def build_text_extraction_prompt(contexts):
    """Build simplified prompt for Group A assays where target identity is already known."""
    assay_contexts = []
    for i, context in enumerate(contexts):
        context_str = f"\n== ASSAY {i+1} ==\n"
        context_str += f"Assay ID: {context['assay_id']}\n"
        context_str += f"Assay Organism: {context['assay_organism']}\n"
        context_str += f"Assay Strain: {context['assay_strain']}\n"
        context_str += f"Assay Description: {context['assay_description']}\n"
        context_str += f"Detailed Description: {context['detailed_description']}\n"
        context_str += f"Target Organism (ChEMBL): {context['target_organism']}\n"
        context_str += f"Document Abstract: {context['doc_abstract'][:1000]}{'...' if len(context['doc_abstract']) > 1000 else ''}\n"
        assay_contexts.append(context_str)

    all_contexts = "\n".join(assay_contexts)

    prompt = f"""
You are an information extraction assistant analyzing antimicrobial bioactivity assay records.
Below are {len(contexts)} assay records. The target identity is already known for these assays.

For each assay, extract ONLY the following fields from the text:

- assay_organism_curated: The organism where the assay is actually performed. Validate and correct the "Assay Organism" field using description and abstract. If the description mentions cytotoxicity in human cells, use "Homo sapiens". If the assay organism field is empty or clearly wrong, extract from text. Otherwise confirm or correct the provided value.
- assay_strain_curated: Biological strain name only (e.g., "H37Rv", "K12", "PAO1"). Do NOT include ATCC, DSM, or other catalog identifiers — those belong in atcc_id. If the Assay Strain field is provided but contains mixed content (e.g. "H37Rv ATCC 27294", "ATCC 25618 / H37Rv"), parse it: extract the biological strain name only into assay_strain_curated and the ATCC number into atcc_id. If no strain name is found in the text, also check Target Organism (ChEMBL): the format "Organism (strain ATCC XXXXX / StrainName)" contains the strain name after the slash.
- atcc_id: ATCC catalog number in format "ATCC XXXXX". Extract from text or Target Organism field.
- mutations: Explicit mutations only as ["S450L", "H526Y"]. Format: one-letter amino acid + position + one-letter amino acid. Convert long form (e.g., Ser450Leu) to short form if unambiguous.
- known_drug_resistances: Drugs for which resistance is explicitly stated (e.g., ["rifampicin", "isoniazid"]). Do NOT infer from mutations. Only include explicitly stated resistances.
- culture_media: Growth medium (e.g., "Middlebrook 7H9", "Mueller-Hinton").

RULES:
1. Use empty strings "" for missing values, empty lists [] for missing arrays.
2. Replace commas in strings with semicolons.
3. Only extract information explicitly stated in the text. Do not infer.
4. Never put an ATCC number in assay_strain_curated — those belong in atcc_id.

Return JSON with "results" array containing one object per assay, in order.

{all_contexts}
"""
    return prompt


def process_text_only_batch(assay_triplets_batch, pathogen_name, organism_chembl_id, _is_fallback=False):
    """Process Group A assays: pre-fill target fields, use LLM only for text extraction."""
    empty_text_result = {
        "assay_organism_curated": "",
        "assay_strain_curated": "",
        "atcc_id": "",
        "mutations": [],
        "known_drug_resistances": [],
        "culture_media": ""
    }

    contexts = [build_extraction_context(triplet) for triplet in assay_triplets_batch]
    prompt = build_text_extraction_prompt(contexts)
    schema = AssayDetailsTextOnlyBatch.model_json_schema()
    MAX_RETRIES = 3

    for attempt in range(MAX_RETRIES):
        try:
            response = ollama.chat(
                messages=[{"role": "user", "content": prompt}],
                model=LLM_MODEL,
                options={"temperature": 0, "num_ctx": 12288},
                keep_alive="1h",
                format=schema,
            )
            result = json.loads(response.message["content"])
            if len(result.get("results", [])) == len(assay_triplets_batch):
                merged = []
                for triplet, text_result in zip(assay_triplets_batch, result["results"]):
                    prefilled = prefill_target_fields(triplet)
                    # LLM assay_organism_curated takes priority; fall back to prefilled if LLM returned empty
                    if not text_result.get('assay_organism_curated'):
                        text_result['assay_organism_curated'] = prefilled.get('assay_organism_curated', '')
                    merged.append({**prefilled, **text_result})
                return merged
        except Exception:
            pass

    # All retries failed
    if _is_fallback or len(assay_triplets_batch) == 1:
        print(f"  Warning: Failed to process text-only batch of {len(assay_triplets_batch)} assay(s) after {MAX_RETRIES} attempts; filling with empty values")
        results = []
        for triplet in assay_triplets_batch:
            prefilled = prefill_target_fields(triplet)
            results.append({**prefilled, **empty_text_result})
        return results

    results = []
    for triplet in assay_triplets_batch:
        single = process_text_only_batch([triplet], pathogen_name, organism_chembl_id, _is_fallback=True)
        results.append(single[0])
    return results


def build_extraction_context(assay_triplet):
    """Build comprehensive context from all available text sources"""
    assay_id = assay_triplet['assay_id']

    # Structured metadata
    ai = assay_info.get(assay_id, {})

    # Natural language descriptions
    detailed_desc = descriptions.get(assay_id, "")

    # Publication context
    doc_id = ai.get("doc_id") if isinstance(ai, dict) else None
    di = doc_info.get(doc_id, {}) if doc_id and not pd.isna(doc_id) else {}

    # Target information
    tid = ai.get("tid") if isinstance(ai, dict) else None
    ti = target_info.get(tid, {}) if tid and not pd.isna(tid) else {}

    return {
        "assay_id": assay_id,
        "activity_type": assay_triplet['activity_type'],
        "unit": str(assay_triplet['unit']) if pd.notna(assay_triplet['unit']) else "",
        "assay_type": assay_type_map.get(assay_triplet.get('assay_type', ''), assay_triplet.get('assay_type', '')),
        "assay_organism": assay_triplet.get('assay_organism', ''),
        "assay_strain": ai.get("assay_strain", "") if isinstance(ai, dict) else "",
        "assay_description": ai.get("description", "") if isinstance(ai, dict) else "",
        "detailed_description": detailed_desc,
        "assay_category": ai.get("assay_category", "") if isinstance(ai, dict) else "",
        "assay_test_type": ai.get("assay_test_type", "") if isinstance(ai, dict) else "",
        "assay_cell_type": ai.get("assay_cell_type", "") if isinstance(ai, dict) else "",
        "target_name": ti.get("pref_name", "") if isinstance(ti, dict) else "",
        "target_type_chembl": ti.get("target_type", "") if isinstance(ti, dict) else "",
        "target_organism": ti.get("organism", "") if isinstance(ti, dict) else "",
        "target_chembl_id": ti.get("chembl_id", "") if isinstance(ti, dict) else "",
        "doc_title": di.get("title", "") if isinstance(di, dict) else "",
        "doc_abstract": str(di.get("abstract", "")) if isinstance(di, dict) and di.get("abstract") is not None else "",
        "doc_journal": di.get("journal", "") if isinstance(di, dict) else "",
        "doc_pubmed_id": str(di.get("pubmed_id", "")) if isinstance(di, dict) else "",
        "doc_doi": di.get("doi", "") if isinstance(di, dict) else ""
    }


def build_pathogen_target_list(pathogen_name):
    """Build compact list of known SINGLE PROTEIN targets for the pathogen.

    target_dictionary is already filtered to the pathogen at load time, so
    we only need to restrict to SINGLE PROTEIN target type here.
    """
    pathogen_targets = target_dictionary[
        target_dictionary['target_type'] == 'SINGLE PROTEIN'
    ].copy()

    if len(pathogen_targets) == 0:
        return ""

    lines = [f"\nKnown {pathogen_name} SINGLE PROTEIN targets (use ONLY these ChEMBL IDs for SINGLE PROTEIN assignments):"]
    for _, row in pathogen_targets.iterrows():
        syns = target_synonym_map.get(row['tid'], [])
        syn_str = f" (also known as: {', '.join(syns)})" if syns else ""
        lines.append(f"  {row['chembl_id']}: {row['pref_name']}{syn_str}")
    return "\n".join(lines)


def build_organism_context(pathogen_name, organism_chembl_id):
    """Build organism-specific context for the current pathogen"""
    if not organism_chembl_id:
        organism_chembl_id = "UNKNOWN"

    # Determine organism type for context
    if pathogen_name in ["Candida albicans"]:
        organism_type = "fungus"
        activity_examples = ["antifungal activity", "growth inhibition"]
    elif pathogen_name in ["Plasmodium falciparum", "Schistosoma mansoni"]:
        organism_type = "parasite"
        activity_examples = ["antiparasitic activity", "growth inhibition"]
    else:
        organism_type = "bacteria"
        activity_examples = ["antibacterial activity", "growth inhibition", "MIC determination"]

    return f"""
Key Target ChEMBL IDs for {pathogen_name} and common patterns:
- {organism_chembl_id}: {pathogen_name} (ORGANISM)
- CHEMBL612545: Unchecked/Unknown target

Organism-specific inference patterns for {pathogen_name} ({organism_type}):
- If Target ChEMBL ID is CHEMBL612545 (Unchecked) AND Assay Organism is "{pathogen_name}"
  AND assay appears to be organism-level ({', '.join(activity_examples)}, phenotypic, whole cell)
  → target_type_curated: ORGANISM, target_chembl_id_curated: {organism_chembl_id}
- If Target ChEMBL ID is CHEMBL612545 (Unchecked) AND description mentions specific protein names
  → target_type_curated: SINGLE PROTEIN, extract protein name into target_name_curated (ChEMBL ID will be resolved automatically)
- If Target Type (ChEMBL) is "NON-MOLECULAR" or "NON-PROTEIN TARGET"
  AND Assay Organism is "{pathogen_name}"
  AND assay appears to be organism-level (phenotypic, whole cell, growth-based, {', '.join(activity_examples)})
  → target_type_curated: ORGANISM, target_chembl_id_curated: {organism_chembl_id}
- If Target Type (ChEMBL) is "NON-MOLECULAR" or "NON-PROTEIN TARGET"
  AND description mentions specific protein names
  → target_type_curated: SINGLE PROTEIN, extract protein name into target_name_curated (ChEMBL ID will be resolved automatically)
- If assay mentions cytotoxicity, cell viability, or human cell lines
  → assay_organism_curated should be the cell line organism (e.g., Homo sapiens), NOT {pathogen_name}
"""


def build_inference_examples(pathogen_name, organism_chembl_id):
    """Build organism-agnostic inference examples"""
    if not organism_chembl_id:
        organism_chembl_id = "UNKNOWN"

    # Determine activity type based on organism
    if pathogen_name in ["Candida albicans"]:
        activity_prefix = "anti-fungal"
        protein_example = "ergosterol biosynthesis enzyme"
    elif pathogen_name in ["Plasmodium falciparum"]:
        activity_prefix = "anti-malarial"
        protein_example = "dihydrofolate reductase"
    elif pathogen_name in ["Schistosoma mansoni"]:
        activity_prefix = "anti-schistosomal"
        protein_example = "thioredoxin glutathione reductase"
    else:
        activity_prefix = "antibacterial"
        protein_example = "DNA gyrase"

    pathogen_short = pathogen_name.split()[-1].lower()  # e.g., "tuberculosis", "coli", "aureus"

    return f"""
INFERENCE EXAMPLES for {pathogen_name}:

- Assay Organism: "{pathogen_name}", Target: CHEMBL612545 (Unchecked), Description: "growth inhibition"
  → target_type_curated: ORGANISM, target_name_curated: {pathogen_name},
    target_chembl_id_curated: {organism_chembl_id}, assay_organism_curated: {pathogen_name}

- Assay Organism: "{pathogen_name}", Target: CHEMBL612545, Description: "{protein_example} inhibition"
  → target_type_curated: SINGLE PROTEIN, target_name_curated: {protein_example},
    target_chembl_id_curated: "",
    assay_organism_curated: {pathogen_name}
  NOTE: extract the protein name as-is (including abbreviations like "InhA", "rpoB", "DPRE1").
  The ChEMBL ID will be resolved automatically from target_name_curated — leave it empty ("").

- Description: "cytotoxicity in HepG2 cells", compounds tested for anti-{pathogen_short} activity
  → assay_organism_curated: Homo sapiens (NOT {pathogen_name}!)

- Assay Strain field: "StrainName ATCC XXXXX" or "ATCC XXXXX / StrainName"
  → assay_strain_curated: StrainName, atcc_id: ATCC XXXXX
  NOTE: always split combined Assay Strain strings; never copy the raw value as-is.

- Assay text mentions both strain and ATCC: "{pathogen_name} strain H37Rv, ATCC 25618"
  → assay_strain_curated: H37Rv, atcc_id: ATCC 25618

- Target Organism (ChEMBL) shows: "{pathogen_name} (strain ATCC 25618 / H37Rv)"
  (no strain mentioned in assay text)
  → assay_strain_curated: H37Rv, atcc_id: ATCC 25618
  NOTE: extract the strain name (H37Rv) from after the slash in Target Organism.

- Only ATCC number is mentioned, no strain name anywhere:
  → assay_strain_curated: "", atcc_id: ATCC XXXXX
  NEVER put an ATCC number in assay_strain_curated even if it is the only identifier available.

- Assay Organism: "{pathogen_name}", Target Type: NON-MOLECULAR, Target ChEMBL ID: CHEMBL3879801,
  Description: "growth inhibition / MIC / agar dilution method against {pathogen_name}"
  → target_type_curated: ORGANISM, target_name_curated: {pathogen_name},
    target_chembl_id_curated: {organism_chembl_id}, assay_organism_curated: {pathogen_name}
"""


def resolve_target_chembl_id(target_name: str) -> str:
    """Resolve a protein name to its ChEMBL ID using multi-strategy lookup.

    Strategy order (returns on first hit):
    1. Exact case-insensitive match
    2. Strip organism prefix → exact match
    3. Slash split on full name → exact match each part
    4. Strip organism prefix + slash split → exact match each part
    5. GyrA/B-style reconstruction (short second part after slash)
    6. Token-level slash split — for each whitespace token containing "/",
       split and exact-match each sub-part (e.g. "DNA GyrA/B heterotetramer")
    """
    if not target_name:
        return ""

    def hit(q):
        return _target_lookup.get(q.strip().lower(), "")

    # 1. Exact
    r = hit(target_name)
    if r:
        return r

    # 2. Strip organism prefix
    stripped = _strip_organism_prefix(target_name, _org_prefixes)
    if stripped != target_name:
        r = hit(stripped)
        if r:
            return r

    # 3. Slash split on original name
    if "/" in target_name:
        for part in target_name.split("/"):
            r = hit(part.strip())
            if r:
                return r

    # 4. Slash split after organism strip
    if stripped != target_name and "/" in stripped:
        for part in stripped.split("/"):
            r = hit(part.strip())
            if r:
                return r

    # 5. GyrA/B-style reconstruction (second part ≤ 3 chars → expand using prefix of first)
    if "/" in stripped:
        parts = [p.strip() for p in stripped.split("/")]
        if len(parts) == 2 and len(parts[1]) <= 3:
            base = re.sub(r"[A-Za-z0-9]+$", "", parts[0])
            if base:
                r = hit(base + parts[1])
                if r:
                    return r

    # 6. Token-level slash split — only fires on tokens containing "/", safer than word-token matching
    for tok in stripped.split():
        if "/" in tok:
            for part in tok.split("/"):
                if part.strip():
                    r = hit(part.strip())
                    if r:
                        return r
            tok_parts = [p.strip() for p in tok.split("/")]
            if len(tok_parts) == 2 and len(tok_parts[1]) <= 3:
                base = re.sub(r"[A-Za-z0-9]+$", "", tok_parts[0])
                if base:
                    r = hit(base + tok_parts[1])
                    if r:
                        return r

    return ""


def build_multi_assay_prompt(contexts, pathogen_name, organism_chembl_id):
    """Build optimized prompt for batch processing multiple assays with organism-agnostic target inference"""

    # Build organism-specific target dictionary context for inference
    target_context = f"""
REFERENCE INFORMATION FOR TARGET INFERENCE:

{build_organism_context(pathogen_name, organism_chembl_id)}
"""

    assay_contexts = []
    for i, context in enumerate(contexts):
        context_str = f"\n== ASSAY {i+1} ==\n"
        context_str += f"Assay ID: {context['assay_id']}\n"
        context_str += f"Activity Type: {context['activity_type']}\n"
        context_str += f"Unit: {context['unit']}\n"
        context_str += f"Assay Type: {context['assay_type']}\n"
        context_str += f"Assay Organism: {context['assay_organism']}\n"
        context_str += f"Assay Strain: {context['assay_strain']}\n"
        context_str += f"Assay Description: {context['assay_description']}\n"
        context_str += f"Detailed Description: {context['detailed_description']}\n"
        context_str += f"Assay Category: {context['assay_category']}\n"
        context_str += f"Assay Test Type: {context['assay_test_type']}\n"
        context_str += f"Assay Cell Type: {context['assay_cell_type']}\n"
        context_str += f"Target Name (ChEMBL): {context['target_name']}\n"
        context_str += f"Target Type (ChEMBL): {context['target_type_chembl']}\n"
        context_str += f"Target Organism (ChEMBL): {context['target_organism']}\n"
        context_str += f"Target ChEMBL ID: {context['target_chembl_id']}\n"
        context_str += f"Document Title: {context['doc_title']}\n"
        context_str += f"Document Abstract: {context['doc_abstract'][:1000]}{'...' if len(context['doc_abstract']) > 1000 else ''}\n"
        context_str += f"Document Journal: {context['doc_journal']}\n"
        context_str += f"Document PubMed ID: {context['doc_pubmed_id']}\n"
        context_str += f"Document DOI: {context['doc_doi']}\n"
        assay_contexts.append(context_str)

    all_contexts = "\n".join(assay_contexts)

    prompt = f"""
You are an information extraction assistant specialized in analyzing antimicrobial bioactivity data from ChEMBL.
Below are {len(contexts)} assay records that need enhanced information extraction with intelligent target inference.

{target_context}

For each assay, extract and infer the following fields:

- target_type_curated: Enhanced target classification. For UNCHECKED targets (CHEMBL612545),
  NON-MOLECULAR targets, and NON-PROTEIN TARGET targets, infer from context:
  * If assay organism is pathogen AND assay is phenotypic/growth-based → ORGANISM
  * If description mentions specific proteins → SINGLE PROTEIN
  * Otherwise use: SINGLE PROTEIN, PROTEIN COMPLEX, PROTEIN FAMILY, ORGANISM, CELL-LINE, or DISCARDED
  * Set to DISCARDED if: (a) assay is explicitly not a bioactivity assay (e.g., transcriptomics), OR
    (b) both an explicit target name/ID AND an explicit organism/cell-line are missing from annotations.

- target_name_curated: Target name from text OR inferred target name:
  * For ORGANISM assays → use organism name (e.g., "Mycobacterium tuberculosis"), except where we are assessing cytotoxicity in human cells for instance, where the organism name would be Homo sapiens
  * For SINGLE PROTEIN → extract protein name from descriptions

- target_chembl_id_curated: Target ChEMBL ID from text OR inferred ID:
  * For {pathogen_name} organism assays → {organism_chembl_id}
  * For SINGLE PROTEIN assays → leave this field empty (""). The ChEMBL ID will be resolved automatically from target_name_curated.

- assay_organism_curated: CRITICAL - The organism where the assay is actually performed:
  * For pathogen growth assays → pathogen organism (e.g., "{pathogen_name}")
  * For cytotoxicity/human cell assays → "Homo sapiens" (NOT the pathogen!)
  * For animal model assays → animal organism
  * BE VERY CAREFUL: compounds active against {pathogen_name} tested in human cells = Homo sapiens

- assay_strain_curated: Biological strain name only (e.g., "H37Rv", "K12", "PAO1"). Do NOT include ATCC, DSM, or other catalog identifiers — those belong in atcc_id. If the Assay Strain field is provided but contains mixed content (e.g. "H37Rv ATCC 27294", "ATCC 25618 / H37Rv"), parse it: extract the biological strain name only into assay_strain_curated and the ATCC number into atcc_id. If no strain name is found in the text, also check Target Organism (ChEMBL): the format "Organism (strain ATCC XXXXX / StrainName)" contains the strain name after the slash. If no strain name can be found anywhere, leave assay_strain_curated empty ("").
- atcc_id: ATCC catalog number in format "ATCC XXXXX"
- mutations: Explicit mutations only as ["S450L", "H526Y"]. Format: one-letter amino acid + position + one-letter amino acid. If text uses long form (e.g., Ser450Leu) and conversion is unambiguous, convert to short form.
- known_drug_resistances: Drugs for which resistance is explicitly stated (e.g., ["rifampicin", "isoniazid"]).
  Do NOT infer resistance from mutations. Only include resistances explicitly stated in the text.
- culture_media: Growth medium (e.g., "Middlebrook 7H9", "Mueller-Hinton")

{build_inference_examples(pathogen_name, organism_chembl_id)}

CRITICAL RULES:
1. For UNCHECKED, NON-MOLECULAR, and NON-PROTEIN TARGET targets, ALWAYS try to infer the correct target type and ChEMBL ID
2. Be very careful about assay_organism_curated - distinguish between assay organism and pathogen being studied
3. Use empty strings "" for missing values, empty lists [] for missing arrays
4. Replace commas in strings with semicolons
5. Mutations in single-letter format (S450L not Ser450Leu)
6. If assay_type is Binding and Target Type (ChEMBL) is UNCHECKED, target_type_curated can only be DISCARDED or SINGLE PROTEIN — never ORGANISM.
7. For all fields except target_type_curated, only extract information explicitly stated in the text. Do not infer strain, media, ATCC, or mutations from context.

Return JSON with "results" array containing one AssayDetails object per assay, in order.

{all_contexts}
"""

    return prompt


def process_assay_batch(assay_triplets_batch, pathogen_name, organism_chembl_id, batch_size=6, _is_fallback=False):
    """Process multiple assays in single OLLAMA call for efficiency"""

    empty_result = {
        "target_type_curated": "",
        "target_name_curated": "",
        "target_chembl_id_curated": "",
        "assay_organism_curated": "",
        "assay_strain_curated": "",
        "atcc_id": "",
        "mutations": [],
        "known_drug_resistances": [],
        "culture_media": ""
    }

    # Build comprehensive contexts for all assays in batch
    contexts = [build_extraction_context(triplet) for triplet in assay_triplets_batch]

    # Build multi-assay prompt with organism-specific context
    prompt = build_multi_assay_prompt(contexts, pathogen_name, organism_chembl_id)

    schema = AssayDetailsBatch.model_json_schema()
    MAX_RETRIES = 3

    for attempt in range(MAX_RETRIES):
        try:
            response = ollama.chat(
                messages=[{"role": "user", "content": prompt}],
                model=LLM_MODEL,
                options={"temperature": 0, "num_ctx": 12288},
                keep_alive="1h",
                format=schema,
            )

            result = json.loads(response.message["content"])

            if len(result.get("results", [])) == len(assay_triplets_batch):
                # Resolve SINGLE PROTEIN ChEMBL IDs in Python
                for item in result["results"]:
                    if item.get("target_type_curated") == "SINGLE PROTEIN":
                        if not item.get("target_chembl_id_curated"):
                            item["target_chembl_id_curated"] = resolve_target_chembl_id(
                                item.get("target_name_curated", "")
                            )
                return result["results"]
            # Wrong count — retry

        except Exception:
            pass  # retry

    # All retries failed — fall back to individual processing (only for non-singleton batches)
    if _is_fallback or len(assay_triplets_batch) == 1:
        print(f"  Warning: Failed to process batch of {len(assay_triplets_batch)} assay(s) after {MAX_RETRIES} attempts; filling with empty values")
        return [empty_result] * len(assay_triplets_batch)

    results = []
    for triplet in assay_triplets_batch:
        single = process_assay_batch([triplet], pathogen_name, organism_chembl_id, 1, _is_fallback=True)
        results.append(single[0])
    return results


def get_processed_triplets(output_file):
    """Get set of already processed (assay_id, activity_type, unit) triplets for resume support"""
    try:
        processed_df = pd.read_csv(output_file)
        return set(zip(processed_df['assay_id'], processed_df['activity_type'], processed_df['unit']))
    except:
        return set()


def save_batch_results(assay_triplets_batch, extracted_info_batch, output_file):
    """Save batch results to output file"""

    # Prepare output columns - all step 08 columns plus curated columns adjacent to originals
    output_cols = [
        "assay_id", "assay_type",
        "assay_organism", "assay_organism_curated",
        "assay_tax_id",
        "assay_strain", "assay_strain_curated",
        "assay_confidence_score",
        "doc_chembl_id",
        "target_type", "target_type_curated",
        "target_chembl_id", "target_chembl_id_curated",
        "target_organism", "target_name_curated",
        "target_tax_id", "tid",
        "bao_label", "source_label",
        "activity_type", "unit", "activities", "nan_values", "cpds", "act_flag", "inact_flag",
        "frac_cs", "direction",
        "atcc_id", "mutations", "known_drug_resistances", "culture_media"
    ]

    # Write header if file doesn't exist
    file_exists = os.path.exists(output_file)

    with open(output_file, "a", newline="") as f:
        writer = csv.writer(f)

        if not file_exists:
            writer.writerow(output_cols)

        # Write each triplet with its extracted info
        for triplet, extracted in zip(assay_triplets_batch, extracted_info_batch):
            # Fallback: propagate assay_organism to assay_organism_curated when LLM returned empty
            if not extracted.get('assay_organism_curated'):
                extracted['assay_organism_curated'] = str(triplet.get('assay_organism', '') or '')

            row = [
                triplet['assay_id'], triplet['assay_type'],
                triplet['assay_organism'], extracted.get('assay_organism_curated', ''),
                triplet.get('assay_tax_id', ''),
                triplet.get('assay_strain', ''), extracted.get('assay_strain_curated', ''),
                triplet.get('assay_confidence_score', ''),
                triplet['doc_chembl_id'],
                triplet['target_type'], extracted.get('target_type_curated', ''),
                triplet['target_chembl_id'], extracted.get('target_chembl_id_curated', ''),
                triplet['target_organism'], extracted.get('target_name_curated', ''),
                triplet.get('target_tax_id', ''), triplet.get('tid', ''),
                triplet['bao_label'], triplet['source_label'],
                triplet['activity_type'], triplet['unit'], triplet['activities'],
                triplet['nan_values'], triplet['cpds'], triplet['act_flag'],
                triplet['inact_flag'], triplet['frac_cs'], triplet['direction'],
                extracted.get('atcc_id', ''),
                ';'.join(extracted.get('mutations', [])),
                ';'.join(extracted.get('known_drug_resistances', [])),
                extracted.get('culture_media', '')
            ]
            writer.writerow(row)


def main():
    """Main processing function with batching and resume support"""

    # Load assay triplets from step 08
    input_file = os.path.join(OUTPUT, "08_assays_cleaned.csv")
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        print("Please run step 08 first.")
        sys.exit(1)

    assays_df = pd.read_csv(input_file)
    print(f"Loaded {len(assays_df)} assay triplets from step 08")

    # Resume support - skip already processed triplets
    output_file = os.path.join(OUTPUT, "09_assays_parameters.csv")
    processed_triplets = get_processed_triplets(output_file)

    if processed_triplets:
        print(f"Resuming - {len(processed_triplets)} triplets already processed, skipping")
        mask = assays_df.apply(
            lambda r: (r['assay_id'], r['activity_type'], r['unit']) in processed_triplets, axis=1
        )
        assays_df = assays_df[~mask]
        print(f"Processing remaining {len(assays_df)} triplets")

    if len(assays_df) == 0:
        print("All triplets already processed!")
        return expand_to_triplets(assays_df, output_file)

    # Classify assays into groups
    assays_df['_group'] = assays_df.apply(classify_assay, axis=1)
    group_a = assays_df[assays_df['_group'] == 'text_only']
    group_b = assays_df[assays_df['_group'] == 'full_inference']
    print(f"Group A (text extraction only, target pre-filled): {len(group_a)} assays")
    print(f"Group B (full inference): {len(group_b)} assays")

    batch_size = 6  # Optimized for context window

    # Process Group A with simplified prompt (no target list, no inference examples)
    if len(group_a) > 0:
        print("Processing Group A assays...")
        for i in tqdm(range(0, len(group_a), batch_size), desc="Group A batches"):
            batch = group_a.iloc[i:i+batch_size]
            batch_dicts = batch.to_dict('records')
            extracted_info = process_text_only_batch(batch_dicts, pathogen, organism_chembl_id)
            save_batch_results(batch_dicts, extracted_info, output_file)

    # Process Group B with full inference prompt
    if len(group_b) > 0:
        print("Processing Group B assays...")
        for i in tqdm(range(0, len(group_b), batch_size), desc="Group B batches"):
            batch = group_b.iloc[i:i+batch_size]
            batch_dicts = batch.to_dict('records')
            extracted_info = process_assay_batch(batch_dicts, pathogen, organism_chembl_id, batch_size)
            save_batch_results(batch_dicts, extracted_info, output_file)

    print(f"Enhanced assay information saved to: {output_file}")

    # Expand assay-level results to all triplets
    return expand_to_triplets(assays_df, output_file)


def expand_to_triplets(assays_df, enhanced_file):
    """Expand assay-level results to all triplets and create final output"""

    print("Expanding assay-level results to all triplets...")

    # Load enhanced results
    enhanced_df = pd.read_csv(enhanced_file)

    # Create lookup dictionary for enhanced info by assay_id
    enhanced_lookup = {}
    for _, row in enhanced_df.iterrows():
        enhanced_lookup[row['assay_id']] = {
            'target_type_curated': row['target_type_curated'],
            'target_name_curated': row['target_name_curated'],
            'target_chembl_id_curated': row['target_chembl_id_curated'],
            'assay_organism_curated': row['assay_organism_curated'],
            'assay_strain_curated': row['assay_strain_curated'],
            'atcc_id': row['atcc_id'],
            'mutations': row['mutations'],
            'known_drug_resistances': row['known_drug_resistances'],
            'culture_media': row['culture_media']
        }

    # Add enhanced columns to all triplets and reorder for easier analysis
    for col in ['target_type_curated', 'target_name_curated', 'target_chembl_id_curated',
                'assay_organism_curated', 'assay_strain_curated', 'atcc_id', 'mutations', 'known_drug_resistances', 'culture_media']:
        assays_df[col] = assays_df['assay_id'].map(lambda x: enhanced_lookup.get(x, {}).get(col, ''))

    # Reorder columns to place _curated next to originals for easier comparison
    ordered_cols = [
        'assay_id', 'assay_type',
        'assay_organism', 'assay_organism_curated',
        'assay_tax_id',
        'assay_strain', 'assay_strain_curated',
        'assay_confidence_score',
        'doc_chembl_id',
        'target_type', 'target_type_curated',
        'target_chembl_id', 'target_chembl_id_curated',
        'target_organism', 'target_name_curated',
        'target_tax_id', 'tid',
        'bao_label', 'source_label',
        'activity_type', 'unit',
        'activities', 'nan_values', 'cpds', 'act_flag', 'inact_flag', 'frac_cs', 'direction',
        'atcc_id', 'mutations', 'known_drug_resistances', 'culture_media'
    ]

    assays_df = assays_df[ordered_cols]

    # Save final result
    final_output = os.path.join(OUTPUT, "09_assays_parameters_full.csv")
    assays_df.to_csv(final_output, index=False)

    print(f"Final enhanced dataset with all triplets saved to: {final_output}")
    print(f"Total triplets: {len(assays_df)}")

    # Clean up intermediary file
    if os.path.exists(enhanced_file):
        os.remove(enhanced_file)
        print(f"Removed intermediary file: {os.path.basename(enhanced_file)}")

    # Summary statistics
    enhanced_count = (assays_df['target_type_curated'] != '').sum()
    strain_count = (assays_df['assay_strain_curated'] != '').sum()
    atcc_count = (assays_df['atcc_id'] != '').sum()
    mutations_count = (assays_df['mutations'] != '').sum()
    resistance_count = (assays_df['known_drug_resistances'] != '').sum()

    print(f"\nEnhancement Summary:")
    print(f"  Triplets with enhanced target info: {enhanced_count}/{len(assays_df)} ({100*enhanced_count/len(assays_df):.1f}%)")
    print(f"  Triplets with strain info: {strain_count}/{len(assays_df)} ({100*strain_count/len(assays_df):.1f}%)")
    print(f"  Triplets with ATCC ID: {atcc_count}/{len(assays_df)} ({100*atcc_count/len(assays_df):.1f}%)")
    print(f"  Triplets with mutations: {mutations_count}/{len(assays_df)} ({100*mutations_count/len(assays_df):.1f}%)")
    print(f"  Triplets with known drug resistances: {resistance_count}/{len(assays_df)} ({100*resistance_count/len(assays_df):.1f}%)")

    return final_output


if __name__ == "__main__":
    main()