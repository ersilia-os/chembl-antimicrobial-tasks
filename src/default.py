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