import os

abspath = os.path.dirname(os.path.abspath(__file__))

# Database defaults
DATABASE_NAME = "chembl_36"
CHEMBL_USR = "chembl_user"
CHEMBL_PWD = "aaa"

# Path defaults
CONFIGPATH = os.path.join(abspath, "..", "config")
DATAPATH = os.path.join(abspath, "..", "data")
TMPDIR = os.path.join(abspath, "..", "tmp")

# LLM model for assay parameter curation (step 09)
LLM_MODEL = "gpt-oss:20b"

MIN_ASSAY_SIZE = 0
MIN_POSITIVES = 10

# Modeling thresholds (steps 13–17, 22)
DECOY_RATIO = 0.1                  # target active fraction when adding decoys
MIN_CPDS_CONDITION_A = 1000        # min compounds for condition A
MIN_POSITIVES_CONDITION_A = 50     # min actives for condition A
MIN_POSITIVES_CONDITION_B = 100    # min actives for condition B
AUROC_MIN_THRESHOLD = 0.70         # minimum AUROC to accept a modeled dataset
AUROC_IMPROVEMENT_THRESHOLD = 0.1  # AUROC gain required to prefer best cutoff over mid cutoff
CORRELATION_THRESHOLD = 0.5        # pairwise correlation threshold for dataset deduplication