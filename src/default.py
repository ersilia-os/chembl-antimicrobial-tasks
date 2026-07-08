import os

abspath = os.path.dirname(os.path.abspath(__file__))

# Database defaults
DATABASE_NAME = "chembl_36"
CHEMBL_USR = "chembl_user"
CHEMBL_PWD = "aaa"

# Path defaults
CONFIGPATH = os.path.join(abspath, "..", "config")
DATAPATH = os.path.join(abspath, "..", "data")          # raw ChEMBL DB dumps (data/chembl_activities)
OUTPUTPATH = os.path.join(abspath, "..", "output")
STAGE1PATH = os.path.join(OUTPUTPATH, "stage1")         # Stage-1 processed results (moved from data/chembl_processed)
STAGE2PATH = os.path.join(OUTPUTPATH, "stage2")         # Stage-2 results: per-pathogen (07-11) + general/ (eos1klk)
STAGE3PATH = os.path.join(OUTPUTPATH, "stage3")         # Stage-3 results: conservative pooling branch (13-20), per-pathogen
STAGE4PATH = os.path.join(OUTPUTPATH, "stage4")         # Stage-4 results: signal-based pooling branch (21-28), per-pathogen

# LLM model for assay parameter curation (step 09)
LLM_MODEL = "gpt-oss:20b"

# Project-wide random seed (train/test splits, sampling, plot subsampling)
RANDOM_SEED = 42

# Curation constants (steps 22–24)
MIN_MOLECULES = 5              # datasets with <= MIN_MOLECULES unique InChIKeys are discarded
COADD_SOURCE_LABEL = "COADD"   # source label marking single-shot CO-ADD screens (excluded)
KEPT_TARGET_TYPE = "ORGANISM"  # target type retained after step-1 filtering

# Modeling thresholds (steps 13–17, 22)
DECOY_RATIO = 0.1                  # target active fraction when adding decoys
MIN_CPDS_CONDITION_A = 1000        # min compounds for condition A
MIN_POSITIVES_CONDITION_A = 50     # min actives for condition A
MIN_POSITIVES_CONDITION_B = 100    # min actives for condition B
AUROC_MIN_THRESHOLD = 0.70         # minimum AUROC to accept a modeled dataset
AUROC_IMPROVEMENT_THRESHOLD = 0.1  # AUROC gain required to prefer best cutoff over mid cutoff
CORRELATION_THRESHOLD = 0.5        # pairwise correlation threshold for dataset deduplication

# Pooling / fingerprint constants (steps 23+)
FP_RADIUS = 3          # Morgan fingerprint radius (ECFP6-equivalent)
FP_BITS = 2048         # Morgan fingerprint length (bits)
MIN_TEST = 5           # min disjoint test molecules (with both classes) to score a transfer cell
MIN_DATASET_SIZE = 20  # eligibility floor: pool only datasets with >= this many unique InChIKeys (step 23);
                       # also the "small" cutoff used to label leftovers in step 25

# Strain-name normalization (step 15; reported in 18, 21)
# Culture-collection / catalog acronyms folded out of curated strain names so the
# designator alone identifies the strain (e.g. "ATCC 25922" == "NCTC 25922" == "25922").
# WHO and FDA are included by decision: their panel designators (e.g. "WHO X") are
# preserved — only the acronym is stripped. See scripts/README.md (step 15).
STRAIN_CATALOG_PREFIXES = [
    "ATCC", "DSMZ", "DSM", "NCTC", "NCIMB", "NCIB", "MTCC", "CIP", "CGMCC", "CPCC",
    "CCUG", "KCTC", "JCM", "NRRL", "NBRC", "NRBC", "CECT", "BCRC", "CGSC", "NCCP",
    "NCPF", "IFO", "KCCM", "LMG", "CCM", "CMCC", "CICC", "GIM", "WHO", "FDA",
]